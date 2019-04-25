'''
Created on May 20, 2014

@author: dgrewal
'''
import os
import logging
import warnings
from bamutils import Classifier
from datetime import datetime
from string import Template
import pybamapi
import numpy
import resource

mutationSeq_version = "4.3.8"
MUSEQ_VERSION = mutationSeq_version

class PreProcess(Classifier):
    def __init__(self,args):

        args.samples, self.tumour = self.update_samples(args)
        args.single=True

        super(PreProcess, self).__init__(args) 
        self.__update_args(self.tumour)

    def update_samples(self, args):
        samples = []
        for sample in args.samples:
            s = sample.split(':')
            if s[0] == 'tumour':
                tumour = s[1]
                continue

            samples.append(sample)

        return samples, tumour

   
    def __update_args(self, tumour):
        
        if self.args.deep:
            rmdups = False
        else:
            rmdups = True
            

        logging.info("initializing a PairedBam")
        self.bam = pybamapi.PairedBam(tumour=tumour,
                                      normal=self.samples.get("normal"),
                                      reference=self.samples.get(
                                          "reference"),
                                      coverage=self.coverage,
                                      rmdups=rmdups,
                                      mapq_threshold=self.mapq_threshold,
                                      baseq_threshold=self.baseq_threshold)

        
    def get_features(self):
        tuples = self.bam.get_tuples(self.target_positions)
        for tt,it in tuples: 
            self._update_coverage_info(it)
            chromosome_id = it[-1]
            position = it[0]


            refbase = self.bam.get_reference_base(chromosome_id, position,
                                                  index=True)

            nonrefbases = [x for x in range(4) if x != refbase]

            # ignore tumour tuples with no/few variants in the bam file
            if not self.no_filter:
                if it[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
                        it[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
                        it[nonrefbases[2] + 1][0] < self.args.tumour_variant:
                    continue

          
            ## get corresponding reference tuple
            rt = self.bam.get_reference_tuple(chromosome_id, position)
            
            #MUT-238 If the ref base is 4(N) ignore the position
            if rt[0] >= 4: 
                logging.error(str(position)+' position references base N and has been ignored')
                continue

            # swapped tum and norm when calling outstr, so it returns the values on wrong index. swap them later
            # we're using counts from normal to decide ref and alt, then grabbing counts for those in tumour
            # normally its the other way around.
            outstr = self._make_outstr(it, rt[0], tt)
            normoutstr = outstr[-2]
            tumoutstr = outstr[-1]
            outstr[-2] = tumoutstr
            outstr[-1] = normoutstr

            if outstr[-1][5] == '0/1':
                ## calculate features     
                feature_set = self.features_module.Features(it, rt,'n')
                temp_feature = feature_set.get_features()
            
                self.features_buffer.append(temp_feature)
                self.outstr_buffer.append(outstr)
            
            ## check the buffer size and flush
            if len(self.features_buffer) >= self.buffer_size:
                yield self._flush()
        
        yield self._flush()


    def print_results(self, probabilities_outstrs):
        """
        collects the outstr and probabilities and writes to
        output file in vcf format
        """
        # open the output vcf file to write
        out = open(self.args.out, 'w')
        out_path = str(self.args.out)

        format_str = 'RC:AC:NI:ND:DP:GT:PL'

        # print the vcf header to the output
        header = self._meta_data()
        if header is not None:
            print >> out, header.strip()

        # print the results
        for probabilities, outstrs in probabilities_outstrs:

            logging.info("printing results to: " + out_path)
            for i in xrange(len(probabilities)):
                outstr = outstrs[i]
                prob = probabilities[i]

                filter_flag = outstr[5]

                # do not print positions with p < threshold if --all option is
                # not set
                if not self.args.all and prob < self.args.threshold:
                    continue

                # set the filter_flag
                if filter_flag is None:
                    if prob >= self.args.threshold:
                        filter_flag = "PASS"
                    else:
                        filter_flag = "FAIL"

                info_str = ["PR=", "%.2f" % prob, ";TC=", outstr[6]]

                tum_str = ':'.join([str(v) for v in outstr[7]])
                norm_str = ':'.join([str(v) for v in outstr[8]])

                info_str = ''.join(info_str)

                # calculate phred quality
                phred_quality = self.get_phred_score(prob, typ='quality')

                # alternative base
                altbase = self.base[outstr[4]]

                # make sure it is all strings
                outstr = [outstr[0], outstr[1], outstr[2],
                          self.base[outstr[3]], altbase,
                          "%.2f" % phred_quality, filter_flag, info_str,
                          format_str]

                outstr.append(tum_str)
                outstr.append(norm_str)

                outstr = [str(outstrval) for outstrval in outstr]

                print >> out, "\t".join(outstr)

        out.close()


    def _meta_data(self):
        tumour = self.tumour
        normal = self.samples.get("normal")
        reference = self.samples.get("reference")
        model = self.model

        contigs = self.bam.get_reference_chromosome_lengths()
        contigs = ['##contig=<ID=%s,length=%s>' %(k,v) for k,v in contigs.iteritems()]
        contigs = '\n'.join(contigs)

        if tumour is None:
            tumour = "N/A"

        elif normal is None:
            normal = "N/A"

        cfg_file = self.args.config
        if not cfg_file or not os.path.exists(cfg_file):
            cfg_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'metadata.config')
        try:
            cfg_file = open(cfg_file, 'r')
            header = ""

            for hdrline in cfg_file:
                hdrline = Template(hdrline).substitute(
                        DATETIME=datetime.now().strftime("%Y%m%d"),
                        VERSION=MUSEQ_VERSION,
                        REFERENCE=reference,
                        TUMOUR=tumour,
                        NORMAL=normal,
                        MODEL=model,
                        THRESHOLD=self.args.threshold,
                        INDLTHRESHOLD=self.args.indl_threshold,
                        MAPQTHRESHOLD=self.bam.mapq_threshold,
                        BASEQTHRESHOLD=self.bam.baseq_threshold,
                        COVERAGE=self.bam.coverage,
                        RMDUPS=self.bam.rmdups,
                        CONTIG=contigs)

                # add format section headings
                if hdrline.startswith('#CHROM'):
                    hdrline = hdrline.strip('\n')
                    hdrline += '\tTUMOUR\tNORMAL'
                    hdrline += '\n'

                header += hdrline
            cfg_file.close()

            return header
        except AttributeError as exc:
            exc = exc.strerror if hasattr(exc, 'strerror') else str(exc)
            logging.warning(
                "warning: failed to load metadata file due to error: %s",
                exc)
            return

        except KeyError as exc:
            exc = exc.strerror if hasattr(exc, 'strerror') else str(exc)
            logging.warning(
                "warning: failed to load metadata file due to error: %s",
                exc)
            return

