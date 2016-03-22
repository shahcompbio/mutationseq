'''
Created on May 20, 2014

@author: dgrewal
'''
import logging
import warnings
from bamutils import Classifier
import pybamapi
import numpy
import resource

mutationSeq_version = "4.3.8"

class PreProcess(Classifier):
    def __init__(self,args):

        args.samples, tumour = self.update_samples(args)
        args.single=True

        super(PreProcess, self).__init__(args) 
        self.__update_args(tumour)

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

            outstr = self._make_outstr(it, rt[0],tt)
            

            if outstr[-1][7] == '0/1':
                ## calculate features     
                feature_set = self.features_module.Features(it, rt,'n')
                temp_feature = feature_set.get_features()
            
                self.features_buffer.append(temp_feature)
                self.outstr_buffer.append(outstr)
            
            ## check the buffer size and flush
            if len(self.features_buffer) >= self.buffer_size:
                yield self._flush()
        
        yield self._flush()

