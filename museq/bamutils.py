# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
##TODO: needs a logger
from __future__ import division
import sys
import numpy
import newfeatures, newfeatures_single, newfeatures_deep
from math import log10
from collections import namedtuple
from sklearn.ensemble import RandomForestClassifier
from string import Template
from datetime import datetime
import resource 

mutationSeq_version = "4.0.0"
#DEBUG = False

class BamHelper:
    def __init__(self, bam, args):
        self.samples = {}
        for s in args.samples:
            self.samples[s.split(':')[0]] = s.split(':')[1]
        
        ## check the input
        if not self.samples.get("reference"):
            print "bad input: reference must be specified"
            sys.exit(1)
        
        if not self.samples.get("normal") and not self.samples.get("tumour"):   
            print "bad input: no bam files specified"
            sys.exit(1) 
            
        ## check if it is single sample analysis as well as type of the sample
        if not self.samples.get("normal"):
            single_flag = True
            single_type = "t"
            
        elif not self.samples.get("tumour"):
            single_flag = True
            single_type = "n"
            
        else:
            single_flag = False
            single_type = None
            
        ## set the right feature set
        if single_flag:
            self.features = newfeatures_single

        elif args.deep:
            self.features = newfeatures_deep
        
        else:
            self.features = newfeatures
           
        ## flags for single_sample and deepseq analysis
        Flags = namedtuple("Flags", "deep, single, type")         
        self.flags = Flags._make([args.deep, single_flag, single_type])
       
        self.outstr_buffer   = []    
        self.features_buffer = []
        
        self.bam  = bam
        self.args = args
        self.base = ['A', 'C', 'G', 'T', 'N']
       
        ## set the buffer size to limit the memory usage       
        if self.args.buffer_size <= 0:
            self.buffer_size = float('inf')
        else:
            self.buffer_size = self.args.buffer_size

    def __parse_positions(self, positions_list):
        chromosome = positions_list.split(':')[0]
        try:
            ##check for "chr" in the input interval
            chromosome = chromosome.split('r')[1] 
        except:
            pass
        
        try:
            position = positions_list.split(':')[1]
            start = int(position.split('-')[0])
            try:
                stop = int(position.split('-')[1])
            except:
                stop = start
            return [chromosome, start, stop]
        except:
            return [chromosome, None, None]
        
    def get_positions(self):
        target_positions = []
        
        if self.args.interval:
            temp_tp = self.__parse_positions(self.args.interval)
            target_positions.append(temp_tp) 
        
        elif self.args.positions_file:
            try:
                positions_file = open(self.args.positions_file, 'r')
                for l in positions_file.readlines():
                    temp_tp = self.__parse_positions(l.strip())
                    target_positions.append(temp_tp)
                positions_file.close()
            
            except:
                print "Failed to load the positions file from " + self.args.positions_file
                sys.exit(1)

        else:
            ## get all the common chromosome names
            tcn = self.bam.get_tumour_refnames.keys() # chromosome names in tumour bam
            ncn = self.bam.get_normal_refnames.keys() # chromosome names in normal bam
            chromosome_names = set(tcn).intersection(set(ncn))

            for cn in chromosome_names:
                temp_tp = [cn, None, None]
                target_positions.append(temp_tp)
                
        return target_positions
    
    def __make_outstr(self, tt, nt, refbase):
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
        tc = self.bam.get_trinucleotide_context(tt[-1], tt[0])

        ## generate informatio for the INFO column in the output
        TR = tt[refbase + 1][0] # tumour to reference base count 
        NR = nt[refbase + 1][0] # normal to reference base count 
        TA = tt[altbase + 1][0] # tumour to alternative base count 
        NA = nt[altbase + 1][0] # normal to alternative base count 
        info = [TR, TA, NR, NA, tc, tt[-4], tt[-2]]
        info = map(str, info) 
        
        ## reserved to be filled later
        qual   = None # phred_quality
        out_id = "."  # database ID   
        
        ## get chromosome name of the given chromosome ID
        chromosome_name = self.bam.get_tumour_chromosome_name(tt[-1])
        
        outstr = [chromosome_name, tt[0], out_id, refbase, altbase, qual, filter_flag, info]
        
        return outstr  
    
    def __flush(self):
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
        for tt, nt in tuples:          
            ## ignore tumour tuples with no/few variants in the tumour and too many variants in the normal 
            refbase = self.bam.get_reference_base(tt[-1], tt[0], index=True)
            if tt[5][0] - tt[refbase + 1][0] < self.args.tumour_variant or (nt[5][0] - nt[refbase + 1][0]) / nt[5][0] > (self.args.normal_variant / 100):
                continue
            
            ## get corresponding reference tuple
            chromosome_id = tt[-1]
            position = tt[0]
            rt = self.bam.get_reference_tuple(chromosome_id, position)
            
            ## calculate features and buffer it     
            feature_set = self.features.Features(tt, nt, rt)
            temp_feature = feature_set.get_features()
            self.features_buffer.append(temp_feature)
        
            ## generate output string and buffer it
            outstr = self.__make_outstr(tt, nt, rt[0])
            self.outstr_buffer.append(outstr)
            
            ## check the buffer size and flush
            if len(self.features_buffer) >= self.buffer_size:
                yield self.__flush()
        
        yield self.__flush()
            
    def __fit_model(self):
        try:
            npz = numpy.load(self.samples["model"])
        except:
            print "Failed to load model"
            print sys.exc_info()[0]
            sys.exit(1)

        train  = npz["arr_1"]
        labels = npz["arr_2"]
        model  = RandomForestClassifier(random_state=0, n_estimators=1000, n_jobs=1, compute_importances=True)
        model.fit(train, labels)
        return model   
        
    def predict(self, features_outstrs):
        model = self.__fit_model()

        for features, outstrs in features_outstrs:
            if len(features) == 0:
                yield [], []
                
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
        
        if normal is None:
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
            print "Failed to load metadata file"
            return
            
    def print_results(self, probabilities_outstrs):
        ## open the output vcf file to write        
        if self.args.out is None:
            print "--out is not specified, standard output is used to write the results"
            out = sys.stdout
        else:
            out = open(self.args.out, 'w')
        
        ## print the vcf header to the output
        header = self.__meta_data() 
        if header is not None:
            print >> out, header.strip()
        

        ## print the results
        for probabilities, outstrs in probabilities_outstrs:
            if len(probabilities) == 0:
                continue
                
            for i in xrange(len(probabilities)):
                outstr = outstrs[i]
                p = probabilities[i]
                
                ## check if the filter_flag is INDL or None
                if outstr[-2] is None:
                    if p > self.args.threshold:
                        filter_flag = "PASS"
                    else:
                        filter_flag = "FAIL"
                else:
                    filter_flag = outstr[-2]
                    
                info_str = "PR=" + "%.2f" % p + ";TR=" + outstr[-1][0] + \
                            ";TA=" + outstr[-1][1] + ";NR=" + outstr[-1][2] + \
                            ";NA=" + outstr[-1][3] + ",TC=" + outstr[-1][4] + \
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
                
                ##TODO: remove this line
                print "total mem usage after flush: " + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

        out.close()
        
    def export_features(self, features):
#        names = 
        
        with open(self.args.export, 'w') as export_file:
            for f in features:
                print >> export_file, f


#==============================================================================
# old get_featurs
#==============================================================================
#    def get_features(self, tumour_tuples, normal_tuples):
#        tuples_buffer   = deque()        
#        features_buffer = []
#        
#        ## TODO: remove this, write tuples in a file
#        if DEBUG:        
#            tuple_file = open("tuples.f", 'w')
#            print >> tuple_file, "Tumour tuples"
#            
#        for tt in tumour_tuples:
##            print tt
#            ## ignore tumour tuples with no/few variants compared to reference            
#            refbase = self.bam.get_reference_base(tt[-1], tt[0], index=True)
#            if (tt[5][0] - tt[refbase + 1][0]) / tt[5][0] < 0.1:
#                continue
#
#            ##TODO: remove this            
#            if DEBUG:
#                print >> tuple_file, tt
#                
#            ## buffer tumour tuple to campare against normal tuples
#            tuples_buffer.append(tt)
#            
#        ## return if all tuples were filterd
#        if len(tuples_buffer) == 0:
#            return []
#        
#        tt = tuples_buffer.popleft()
#        
#        ##TODO: remove this
#        if DEBUG:
#            print >> tuple_file, "Normal tuples"
#            
#        for nt in normal_tuples:
##            print nt
#            ## find positions where tuples for both tumour and normal exist
#            while tt[0] < nt[0]:
#                if len(tuples_buffer) == 0:
#                    break
#                tt = tuples_buffer.popleft()
#
#            if tt[0] != nt[0]:
#                continue
#            
#            ##TODO: remvoe this:
#            if DEBUG:
#                print >> tuple_file, nt
#                
#            ## extract reference tuples 
#            print tt
#            rt = self.bam.get_reference_tuple(tt[-1], tt[0])
#
#            ## calculate features      
#            feature_set = self.features.Features(tt, nt, rt)
#            tf = feature_set.get_features()
#            features_buffer.append(tf)
#            
#            ## generate output string and buffer it
#            outstr = self.__make_outstr(tt, rt, nt)
#            self.outstr_buffer.append(outstr)
#            
#            if len(tuples_buffer) == 0:
#                break
#        
#        ##TODO: remove this
#        if DEBUG:
#            tuple_file.close()
#        
#        ## make a numpy array required as an input to the random forest predictor
#        features_buffer = numpy.array(features_buffer)
#        
#        ## make sure a list is returned         
#        if features_buffer is None:
#            return []
#        else:
#            return features_buffer



