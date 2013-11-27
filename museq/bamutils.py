# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
from __future__ import division
import logging
import sys
import numpy
import pybamapi
import resource
import re
import newfeatures, newfeatures_single, newfeatures_deep, newfeatures_deep_single
from math import log10
from sklearn.ensemble import RandomForestClassifier
from string import Template
from datetime import datetime

mutationSeq_version = "4.0.0"
 
class BamHelper(object):
    def __init__(self, args):
        logging.info("initializig BamHelper")
        self.samples = {}
        self.args = args
        self.base = ['A', 'C', 'G', 'T', 'N']
        self.outstr_buffer = []    
        self.features_buffer = []
        self.__get_buffer_size()
        
        ## parse the positional argument to get tumour/normal bam, reference fasta and model
        for s in self.args.samples:
            self.samples[s.split(':')[0]] = s.split(':')[1]
        
        ## check if there is a referece in the input
        if not self.samples.get("reference"):
            logging.error("error: bad input: reference must be specified")
            raise Exception("no referece file in the input.")
        
        self.ref = self.samples.get("reference")
        
        ## check if the model is specified correctly
        self.model = self.samples.get("model")
        if not self.model:
            logging.error("error: bad input: model must be specified in the input")
            raise Exception("no model")
        
        ## check if there is any bam files in the input
        if not self.samples.get("normal") and not self.samples.get("tumour"):   
            logging.error("error: bad input: no bam files specified in the input")
            raise Exception("no bam file")
        
        ## check if it is single mode but there are two bam files instead of one in the input
        if self.samples.get("normal") and self.samples.get("tumour") and self.args.single:   
            logging.error("error: bad input: single mode but two bam files specified in the input")
            raise Exception("single mode but two bam files specified")
            
        ## check if it is not single mode but there is only one bam file in the input
        if (not self.samples.get("normal") or not self.samples.get("tumour")) and not self.args.single:   
            logging.error("error: bad input: one bam file specified in the input but it does not seem to be the single mode.")
            raise Exception("one bam file specified but not the single mode")
         
        ## single mode 
        if self.args.single:
            if self.args.deep:
                self.features = newfeatures_deep_single
                rmdups = False
            
            else:
                self.features = newfeatures_single
                rmdups = True

            if not self.samples.get("tumour"):
                self.type = 'n'

                logging.info("initializig a normal Bam")
                self.bam = pybamapi.Bam(bam=self.samples.get("normal"), reference=self.ref, coverage=self.args.coverage, rmdups=rmdups)
        
            else:
                self.type = 't'
                
                logging.info("initializig a tumour Bam")
                self.bam = pybamapi.Bam(bam=self.samples.get("tumour"), reference=self.ref, coverage=self.args.coverage, rmdups=rmdups)
        
        ## paired mode
        else:
            if self.args.deep:
                self.features = newfeatures_deep
                rmdups = False
                
            else:
                self.features = newfeatures
                rmdups = True
        
            logging.info("initializig a PairedBam")
            self.bam  = pybamapi.PairedBam(tumour=self.samples.get("tumour"), normal=self.samples.get("normal"), 
                                                reference=self.samples.get("reference"), coverage=self.args.coverage, rmdups=rmdups)
        
        ## check if the version of the input model matches that of the feature set  
        try:
            logging.info("loading model")
            self.npz = numpy.load(self.samples.get("model"))
        
        except:
            logging.error("error: failed to load model")
            raise Exception("failed to load model")
        
        if self.npz["arr_0"] != self.features.version:
            logging.error("mismatched feature set versions:"+ str(self.npz["arr_0"]), "and", str(self.features.version))
            raise Exception("mismatched model")
                            
    def __get_buffer_size(self):
        s = re.split('(\d+)', self.args.buffer_size)
        d = int(s[1])
        l = s[2].upper()
        if d < 0 or l not in ('G', 'M'):
            ## relax the restriction on the memory usage
            self.buffer_size = float('inf')
        
        elif d < 100 and l == 'M':
            logging.warning("warning: buffer size is ingnored. It should be >= 100M.")
            self.buffer_size = float('inf')

        elif l == 'G':
            ## every million output and feature strings together
            ## takes about 2G of memory. Also, ~50M memory is required 
            ## to initialize and about another ~50M is required to fit the model.
            self.buffer_size = d * 200000 
        
        elif l == 'M':
            self.buffer_size = (d / 1024) * 200000
        
        
    def __parse_positions(self, positions_list, pch=':'):
        chromosome = positions_list.split(pch)[0]
        try:
            ##check for "chr" in the input interval
            chromosome = chromosome.split('r')[1] 
        
        except:
            pass
        
        try:
            position = positions_list.split(pch)[1]
            start = int(position.split('-')[0])
            
            try:
                stop = int(position.split('-')[1])

            except:
                stop = start
            return [chromosome, start, stop]

        except:
            return [chromosome, None, None]
        
    def get_positions(self, pch=':'):
        logging.info("getting positions")
        target_positions = []
        
        if self.args.interval:
            temp_tp = self.__parse_positions(self.args.interval, pch)
            target_positions.append(temp_tp) 
        
        elif self.args.positions_file:
            logging.info("parsing the position_file")
            try:
                positions_file = open(self.args.positions_file, 'r')
                for l in positions_file.readlines():
                    temp_tp = self.__parse_positions(l.strip(), pch)
                    target_positions.append(temp_tp)
                positions_file.close()
            
            except:
                logging.error("error: failed to load the positions file " + self.args.positions_file)
                raise Exception("failed to load the positions file")

        else:
            ## get all the common chromosome names
            tcn = self.bam.get_tumour_refnames.keys() # chromosome names in tumour bam
            ncn = self.bam.get_normal_refnames.keys() # chromosome names in normal bam
            chromosome_names = set(tcn).intersection(set(ncn))

            for cn in chromosome_names:
                temp_tp = [cn, None, None]
                target_positions.append(temp_tp)
                
        return target_positions
    
    def __make_outstr(self, tt, refbase, nt=None):
        ## flag insertions and deletions 
        if tt[-4] > 0 or tt[-2] > 0:
            filter_flag = "INDL"
        else:
            filter_flag = None
        
        ## find alternative base   
        if refbase == tt[6]:
            altbase = tt[7]
        else:
            altbase = tt[6]
        
        ## get tri-nucleotide context
        chromosome_id = tt[-1]
        position = tt[0]
        tc = self.bam.get_trinucleotide_context(chromosome_id, position)

        ## generate informatio for the INFO column in the output
        TR = tt[refbase + 1][0] # tumour to reference base count 
        TA = tt[altbase + 1][0] # tumour to alternative base count 
       
        if nt is not None:
            NR = nt[refbase + 1][0] # normal to reference base count 
            NA = nt[altbase + 1][0] # normal to alternative base count 

        else:
            NR = "N/A"
            NA = "N/A"
        
        ## format the out string based on the mode and the type of the bam file
        if self.args.single and self.type == 'n':
            info = [NR, NA, TR, TA, tc, tt[-4], tt[-2]]
        
        else:
            info = [TR, TA, NR, NA, tc, tt[-4], tt[-2]]
            
        info = map(str, info) 
        
        ## reserved to be filled later
        qual   = None # phred_quality
        out_id = "."  # database ID   
        
        ## get chromosome name of the given chromosome ID
        chromosome_name = self.bam.get_chromosome_name(chromosome_id)
        
        outstr = [chromosome_name, position, out_id, refbase, altbase, qual, filter_flag, info]
        
        return outstr  
    
    def __flush(self):
        logging.info("flushing memory. Usage was: " + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024) + "M")
        
        ## a numpy array required as an input to the random forest predictor
        features = numpy.array(self.features_buffer)
        outstrs  = self.outstr_buffer
        
        ## empty the buffers
        self.features_buffer = []
        self.outstr_buffer = []
        
        ## should return a list
        if features is None:
            return [], []
            
        return features, outstrs
    
    def get_features(self, tuples):
        logging.info("getting features")
        
        if self.args.single:
            return self.__get_features_single(tuples)
       
        else:
            return self.__get_features_paired(tuples)
    
    def __get_features_single(self, tuples):
        for it in tuples: 
            chromosome_id = it[-1]
            position = it[0]
           
            refbase = self.bam.get_reference_base(chromosome_id, position, index=True)
            nonrefbases = [x for x in range(4) if x != refbase]
            
            ## ignore tumour tuples with no/few variants in the tumour or too many variants in the normal
            if not self.args.no_filter:
                if  it[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
                    it[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
                    it[nonrefbases[2] + 1][0] < self.args.tumour_variant:
                        continue
            
            ## get corresponding reference tuple
            rt = self.bam.get_reference_tuple(chromosome_id, position)
            
            ## calculate features and buffer it     
            feature_set = self.features.Features(it, rt)
            temp_feature = feature_set.get_features()
            self.features_buffer.append(temp_feature)
        
            ## generate output string and buffer it
            outstr = self.__make_outstr(it, rt[0], nt=None)
            self.outstr_buffer.append(outstr)
            
            ## check the buffer size and flush
            if len(self.features_buffer) >= self.buffer_size:
                yield self.__flush()
        
        yield self.__flush()
        
    def __get_features_paired(self, tuples):
        for tt, nt in tuples: 
            chromosome_id = tt[-1]
            position = tt[0]
            
            refbase = self.bam.get_reference_base(chromosome_id, position, index=True)
            nonrefbases = [x for x in range(4) if x != refbase]
            
            ## ignore tumour tuples with no/few variants in the tumour or too many variants in the normal
            if not self.args.no_filter:
                if  tt[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
                    tt[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
                    tt[nonrefbases[2] + 1][0] < self.args.tumour_variant or \
                    (nt[5][0] - nt[refbase + 1][0]) / nt[5][0] > (self.args.normal_variant / 100):
                        continue
            
            ## get corresponding reference tuple
            rt = self.bam.get_reference_tuple(chromosome_id, position)
            
            ## calculate features and buffer it     
            feature_set = self.features.Features(tt, nt, rt)
            temp_feature = feature_set.get_features()
            self.features_buffer.append(temp_feature)
        
            ## generate output string and buffer it
            outstr = self.__make_outstr(tt, rt[0], nt)
            self.outstr_buffer.append(outstr)
            
            ## check the buffer size and flush
            if len(self.features_buffer) >= self.buffer_size:
                yield self.__flush()
        
        yield self.__flush()

    def __fit_model(self):
        train  = self.npz["arr_1"]
        labels = self.npz["arr_2"]
        
        logging.info("running random forest")
        model  = RandomForestClassifier(random_state=0, n_estimators=1000, n_jobs=1, compute_importances=True)
        
        logging.info("fitting model")
        model.fit(train, labels)
        
        return model   
        
    def predict(self, features_outstrs):
        model = self.__fit_model()

        logging.info("predicting probabilities ")
        for features, outstrs in features_outstrs:
            if len(features) == 0:
                continue
            
            probabilities = model.predict_proba(features)
           
           ## return only probabilities of being somatic
            probabilities = [x[1] for x in probabilities] 
            
            yield probabilities, outstrs
    
    def __meta_data(self):
        tumour = self.samples.get("tumour")
        normal = self.samples.get("normal")        
        reference = self.samples.get("reference")
        
        if tumour is None:
            tumour = "N/A"
        
        elif normal is None:
            normal = "N/A"
            
        try:
            cfg_file = open(self.args.config, 'r')
            header = ""
            
            for l in cfg_file:
                l = Template(l).substitute(DATETIME=datetime.now().strftime("%Y%m%d"),
                                           VERSION=mutationSeq_version,
                                           REFERENCE=reference,
                                           TUMOUR=tumour,
                                           NORMAL=normal,
                                           THRESHOLD=self.args.threshold
                                           )
                header += l
            cfg_file.close()
            return header
        
        except:
            logging.warning("warning: failed to load metadata file.")
            return
            
    def print_results(self, probabilities_outstrs):
        ## open the output vcf file to write        
        if self.args.out is None:
            logging.warning("warning: --out is not specified, standard output is used to write the results")
            out = sys.stdout
            out_path = "stdout"
        
        else:
            out = open(self.args.out, 'w')
            out_path = str(self.args.out)
        
        ## print the vcf header to the output
        header = self.__meta_data() 
        if header is not None:
            print >> out, header.strip()
        
        ## print the results
        any_result = False
        for probabilities, outstrs in probabilities_outstrs:
            if len(probabilities) == 0:
                continue
            
            logging.info("printing results to: " + out_path)
            for i in xrange(len(probabilities)):
                outstr = outstrs[i]
                p = probabilities[i]
                
                ## do not print positions with p < threshold if --all option is not set
                if not self.args.all and p < self.args.threshold:
                    continue 

                any_result = True

                ## set the filter_flag to INDL, PASS, or FAIL
                filter_flag = outstr[-2]
                
                if filter_flag is None:
                    if p >= self.args.threshold:
                        filter_flag = "PASS"
                    
                    else:
                        filter_flag = "FAIL"
                
                info_str = "PR=" + "%.2f" % p + ";TR=" + outstr[-1][0] + \
                            ";TA=" + outstr[-1][1] + ";NR=" + outstr[-1][2] + \
                            ";NA=" + outstr[-1][3] + ";TC=" + outstr[-1][4] + \
                            ";NI=" + outstr[-1][5] + ";ND=" + outstr[-1][6]
                
                ## calculate phred quality
                try:
                    phred_quality = -10 * log10(1 - p)
                
                except:
                    phred_quality = 99
                
                ## make sure it is all strings
                outstr = map(str, [outstr[0], outstr[1], outstr[2], self.base[outstr[3]], 
                                   self.base[outstr[4]], "%.2f" % phred_quality, filter_flag, info_str])
                
                print >> out, "\t".join(outstr)
            
        if not any_result:
            print "**no somatic mutation call**"

        out.close()

    def get_feature_names(self):
        tmp_obj = self.features.Features()
        names = tmp_obj.get_feature_names()

        return names
        
    def export_features(self, features):
        logging.info("exporting features to: " + self.args.export_features)
        version = self.features.version
        names = self.get_feature_names()
        
        with open(self.args.export_features, 'w') as export_file:
            print >> export_file, "##features_version:" + version            
            print >> export_file, "\t".join(names)

            for fs in features:
                for f in fs:
                    print >> export_file, f
#                print >> export_file, fs



