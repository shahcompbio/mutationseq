# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
from __future__ import division
import logging
import sys
import os
import numpy
import pybamapi
import resource
import re
import features, features_single, features_deep, features_deep_single
import matplotlib
matplotlib.use('Agg',warn=False)
import matplotlib.pyplot as plt
from math import log10
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
from sklearn import cross_validation
from sklearn.metrics import roc_curve, auc
from string import Template
from datetime import datetime
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages

mutationSeq_version = "4.2.0"

#==============================================================================
# Classifier class 
#==============================================================================
class Classifier(object):
    def __init__(self, args):
        self.samples = {}
        self.args = args
        self.base = ['A', 'C', 'G', 'T', 'N']
        self.outstr_buffer = []    
        self.features_buffer = []
        self.__get_buffer_size()
        
        ## parse the positional argument to get tumour/normal bam, reference fasta and model
        for s in self.args.samples:
            self.samples[s.split(':')[0]] = s.split(':')[1]
        
        ## check if there is a reference in the input
        if not self.samples.get("reference"):
            logging.error("error: bad input: reference must be specified")
            raise Exception("no reference file in the input.")
        
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
                self.features_module = features_deep_single
                rmdups = False
            
            else:
                self.features_module = features_single
                rmdups = True

            if not self.samples.get("tumour"):
                self.type = 'n'

                logging.info("initializing a normal Bam")
                self.bam = pybamapi.Bam(bam=self.samples.get("normal"), reference=self.ref, coverage=self.args.coverage, rmdups=rmdups)
        
            else:
                self.type = 't'
                
                logging.info("initializing a tumour Bam")
                self.bam = pybamapi.Bam(bam=self.samples.get("tumour"), reference=self.ref, coverage=self.args.coverage, rmdups=rmdups)
        
        ## paired mode
        else:
            if self.args.deep:
                self.features_module = features_deep
                rmdups = False
                
            else:
                self.features_module = features
                rmdups = True
        
            logging.info("initializing a PairedBam")
            self.bam  = pybamapi.PairedBam(tumour=self.samples.get("tumour"), normal=self.samples.get("normal"), 
                                                reference=self.samples.get("reference"), coverage=self.args.coverage, rmdups=rmdups)
        
        ## check if the version of the input model matches that of the feature set  
#         try:
#             logging.info("loading model")
#             self.npz = numpy.load(self.samples.get("model"))
#         
#         except:
#             logging.error("error: failed to load model")
#             raise Exception("failed to load model")
#         
#         if self.npz["arr_0"] != self.features_module.version:
#             logging.error("mismatched feature set versions:"+ str(self.npz["arr_0"]), "and", str(self.features_module.version))
#             raise Exception("mismatched model")
        
        if not self.bam.is_matched_reference():
            logging.error("mismatched reference, sounds like the input reference is not the same as the reference used for alignment")
            raise Exception("mismatched reference")
                            
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
#        try:
#            ##check for "chr" in the input interval
#            chromosome = chromosome.split('r')[1] 
#        
#        except:
#            pass
        
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
            tcn = self.bam.get_refnames().keys() # chromosome names in tumour bam
            ncn = self.bam.get_refnames().keys() # chromosome names in normal bam
            chromosome_names = set(tcn).intersection(set(ncn))

            for cn in chromosome_names:
                temp_tp = [cn, None, None]
                target_positions.append(temp_tp)
                
        return target_positions
    
    def __make_outstr(self, tt, refbase, nt=None):
        t_coverage = tt[5][0]
        n_coverage = 0
        
        ## tumour information to print
        if t_coverage == 0:
            ## alternative base
            altbase = "N/A"
            
            ## tumour counts
            TR = 0 # tumour to reference base count
            TA = 0 # tumour to alternative base count
            
            ## indel
            insertion = 0
            deletion  = 0
            
            ## filter flag
            filter_flag = "NOCV"
        
        else: 
            ## alternative base
            major = tt[6]
            minor = tt[7]
            if refbase == major:
                altbase = minor
        
            else:
                altbase = major
                  
            ## tumour counts
            TR = tt[refbase + 1][0] # tumour to reference base count 
            TA = tt[altbase + 1][0] # tumour to alternative base count

            ## indel 
            insertion = tt[-4]
            deletion  = tt[-2]
                        
            ## filter flag
            if deletion > 0 or insertion > 0:
                filter_flag = "INDL"

            else:
                filter_flag = None
            
        ## normal information to print         
        if nt is not None:
            n_coverage = nt[5][0]
        
        if n_coverage == 0:
            ## normal counts
            NR = 0 # normal to reference base count
            NA = 0 # normal to alternative base count
            
        else:
            NR = nt[refbase + 1][0] # normal to reference base count 
            
            ## if it is zero coverage in tumour bam then altbase is "N/A" and so NA should be 0
            if t_coverage == 0:
                NA = 0 
            
            else:
                NA = nt[altbase + 1][0] # normal to alternative base count     
        
        ## tri_nucleotide context
        chromosome_id = tt[-1]
        position = tt[0]
        tc = self.bam.get_trinucleotide_context(chromosome_id, position)

        ## generate information for the INFO column in the output vcf               
        if self.args.single and self.type == 'n':
            info = [NR, NA, TR, TA, tc, insertion, deletion]
            
        else:
            info = [TR, TA, NR, NA, tc, insertion, deletion]
            
        info = map(str, info) 

        ## reserved for database ID, to be filled later
        out_id = "."  
        
        ## get chromosome name of the given chromosome ID
        chromosome_name = self.bam.get_chromosome_name(chromosome_id)
        
        outstr = [chromosome_name, position, out_id, refbase, altbase, filter_flag, info]
        
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
        
            ## ignore tumour tuples with no/few variants in the bam file
            if not self.args.no_filter:
                if  it[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
                    it[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
                    it[nonrefbases[2] + 1][0] < self.args.tumour_variant:
                        continue
            
            ## get corresponding reference tuple
            rt = self.bam.get_reference_tuple(chromosome_id, position)
            
            #MUT-238 If the ref base is 4(N) ignore the position
            if rt[0] >= 4: 
                logging.error(str(position)+' position references base N and has been ignored')
                continue
            
            ## calculate features     
            feature_set = self.features_module.Features(it, rt,self.type)
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
            
            #MUT-238 If the ref base is 4(N) ignore the position
            if rt[0] >= 4: 
                logging.error(str(position)+' position references base N and has been ignored')
                continue

            ## calculate features and buffer it     
            feature_set = self.features_module.Features(tt, nt, rt)
            temp_feature = feature_set.get_features()
            self.features_buffer.append(temp_feature)
        
            ## generate output string and buffer it
            outstr = self.__make_outstr(tt, rt[0], nt)
            self.outstr_buffer.append(outstr)
            
            ## check the buffer size and flush
            if len(self.features_buffer) >= self.buffer_size:
                yield self.__flush()
        
        yield self.__flush()

#     def __fit_model(self):
#         train  = self.npz["arr_1"]
#         labels = self.npz["arr_2"]
#         
#         logging.info("running random forest")
#         model = RandomForestClassifier(random_state=0, n_estimators=1000, n_jobs=1, compute_importances=True)
#         
#         logging.info("fitting model")
#         model.fit(train, labels)

#         return model   
        
    def __load_model(self):
        try:
            logging.info("loading model")
            return joblib.load(self.model)
        
        except:
            logging.error("error: failed to load model")
            raise Exception("failed to load model")
        
    def predict(self, features_outstrs):
        #model = self.__fit_model()
        model = self.__load_model()
        
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
        model = self.samples.get("model")
        
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
                                           MODEL=model,
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
                
                ## set p = 0 for positions with coverage == 0
                filter_flag = outstr[-2]
                if filter_flag == "NOCV":
                    p = 0
                
                ## do not print positions with p < threshold if --all option is not set
                if not self.args.all and p < self.args.threshold:
                    continue 

                any_result = True

                ## set the filter_flag
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
                if p == 0:
                    phred_quality = 0.0
                    
                elif p == 1:
                    phred_quality = 99
                
                else:
                    phred_quality = -10 * log10(1 - p)
                
                ## alternative base
                altbase = outstr[4]
                if altbase != "N/A":
                    altbase = self.base[altbase]
                
                ## make sure it is all strings
                outstr = map(str, [outstr[0], outstr[1], outstr[2], self.base[outstr[3]], 
                                   altbase, "%.2f" % phred_quality, filter_flag, info_str])
                
                print >> out, "\t".join(outstr)
            
        if not any_result:
            print "**no somatic mutation calls**"

        out.close()

    def get_feature_names(self):
        tmp_obj = self.features_module.Features()
        names = tmp_obj.get_feature_names()

        return names
        
    def export_features(self, features):
        version = self.features_module.version
        names = self.get_feature_names()
        
        with open(self.args.export_features, 'w') as export_file:
            print >> export_file, "##features_version:" + version            
            print >> export_file, "\t".join(names)

            for fs in features:
                for f in fs:
                    for i in range(len(f)):
                        for j in range(len(f[i])):
                            print >> export_file, str(f[i][j]) + "\t",
            
            print >> export_file
#                print >> export_file, fs

#==============================================================================
# Trainer class
#==============================================================================
class Trainer(object):
    def __init__(self, args):
        self.args = args
        self.filename = os.path.basename(self.args.out)

        if self.args.single:
            self.feature_module = features_single
    
        elif self.args.deep:
            self.feature_module = features_deep

        else:
            self.feature_module = features
            
        self.version = self.feature_module.version
        
    def __isvalid_label(self, labels):
        for l in labels:
            if l not in ("SOMATIC", "WILDTYPE", "GERMLINE", "HET", "HET_ONE", "HET_GERMLINE", "HOM", "HOM_ONE", "HOM_GERMLINE","CLEAN"):
                return False
                
        return True
        
    def __parse_infiles(self, infiles):
        self.data = defaultdict(list)
       
        positive_labels = [x.upper() for x in self.args.labels.split(',')]
        if not self.__isvalid_label(positive_labels):
            logging.error("unknown labels specified in the input")
            raise Exception("unknown labels specified in the input")
        
        for case in infiles:
            tfile = None
            nfile = None
            rfile = None

            if self.args.deep:
                contamination = (float(10000), float(10000), float(70), float(0))
                
            else:
                contamination = (float(30), float(30), float(70), float(0))
    
            for line in case:
                l = line.strip().split()
                if len(l) < 3:
                    continue
    
                ## parse the line
                if l[0] == "#":
                    if self.args.single:
                        if l[1] == "tumour":
                            tfile = nfile = l[2]
                            self.type = 't'
                        elif l[1] == "normal":
                            tfile = nfile = l[2]
                            self.type = 'n'
                        #if l[1] == "tumour" or l[1] == "normal":
                        #   tfile = nfile = l[2]
                    
                    else:
                        if l[1] == "tumour":
                            tfile = l[2]
    
                        elif l[1] == "normal":
                            nfile = l[2]
    
                    if l[1] == "reference":
                        rfile = l[2]

                    elif l[1] == "contamination":
                        contamination = (float(l[2]), float(l[3]), float(l[4]), float(1))
                    
                    continue

                ## check if required bam files/reference are specified in the training data
                if not all((tfile, nfile, rfile)):
                    logging.warning("'%s' does not contain the required paths to bam/reference files" % os.path.basename(case.name))
                    continue
            
                chromosome = l[0]
                position = int(l[1])
                label_name = l[2]
                if label_name in positive_labels:
                    label = 1

                else:
                    label = -1
                    
                self.data[(tfile, nfile, rfile)].append((chromosome, position, label, contamination,label_name))
    
    def __get_features(self):
        features_buffer = []
        labels_buffer = []
        keys_buffer = []
        file_stream_w = open(self.args.out+'feature_db.txt','w')

        for tfile, nfile, rfile in self.data.keys():            
            logging.info(tfile)
            if not self.args.single:
                logging.info(nfile)
            
            t_bam = pybamapi.Bam(bam=tfile, reference=rfile, coverage=1)
            if not self.args.single:
                n_bam = pybamapi.Bam(bam=nfile, reference=rfile, coverage=1)
            
            for chromosome, position, label, c, label_name in self.data[(tfile, nfile, rfile)]:
                chromosome_id = t_bam.get_chromosome_id(chromosome)
                tt = t_bam.get_tuple(chromosome, position)
                if not self.args.single:
                    nt = n_bam.get_tuple(chromosome, position)            
            
                rt = t_bam.get_reference_tuple(chromosome_id, position)            
                
                #MUT-238 If the ref base is 4(N) ignore the position
                if rt[0] >= 4: 
                    logging.error(str(position)+' position references base N and has been ignored')
                    continue

                ## check for None tuples
                if self.args.single:
                    if not all([tt, rt]):
                        logging.warning(" ".join(["None tuple", tfile, rfile]))
                        continue

                elif not all([tt, nt, rt]):
                    logging.warning(" ".join(["None tuple", tfile, nfile, rfile]))
                    continue
                
                ## calculate features
                if self.args.single:
                    feature_set = self.feature_module.Features(tt, rt, self.type)
                    
                else:
                    feature_set = self.feature_module.Features(tt, nt, rt)

                temp_features = feature_set.get_features()   
                file_stream_w.write(rfile+';'+nfile+';'+tfile+';'+chromosome+';'+str(position)+'\t'+ str(temp_features)+'\t'+label_name+'\n' )
                
                features_buffer.append(temp_features)
                labels_buffer.append(label)
                keys_buffer.append((rfile, nfile, tfile, chromosome, position, label))
        
        file_stream_w.close()        
        self.features = numpy.array(features_buffer)
        self.labels = numpy.array(labels_buffer)
        self.keys = numpy.array(keys_buffer)
        
    def generate(self, infiles=None):
        if infiles is None:
            infiles = self.args.infiles
        
        logging.info("parsing infiles")
        self.__parse_infiles(infiles)

        logging.info("getting features")        
        self.__get_features()
        
        
#     def load(self):
#         try:
#             npz = numpy.load(self.args.model)
#         
#         except:
#             logging.error("failed to load the model: " + self.args.model)
#             raise Exception("failed to load the model")
#         
#         self.version = npz["arr_0"]
#         self.features = npz["arr_1"]
#         self.labels = npz["arr_2"]
    
    def load(self):
        self.model = joblib.load(self.args.model)
        
    def fit(self):
        self.model = RandomForestClassifier(random_state=0, n_estimators=3000, n_jobs=1, compute_importances=True) 
        self.model.fit(self.features, self.labels)
        
    def save(self):
#         numpy.savez(self.args.out, self.version, self.features, self.labels)
        joblib.dump(self.model, self.args.out, compress=9)
    
    def get_feature_importance(self):
        return self.model.feature_importances_
        
    def get_feature_names(self):
        return self.feature_module.Features().get_feature_names()
 
    def print_feature_importance(self):
        with open(self.args.out + "_importance.txt", 'w') as importance_file:
            feature_importance = self.get_feature_importance()
            feature_names = self.get_feature_names()
             
            for importance, feature_name in sorted(zip(feature_importance, feature_names)):
                print >> importance_file, feature_name, importance
        
    def validate(self):
        self.generate(infiles=self.args.validate)
        
        ## NOTE:copied from Jeff/Fatemeh old code
        logging.info("predicting probabilities")
        probs = self.model.predict_proba(self.features)
        voted = probs[:,1]
        fpr, tpr, thresholds = roc_curve(self.labels, voted)
        roc_auc = auc(fpr, tpr)
        logging.info("AUC:%f" % roc_auc)
        
        fd = open(self.args.out + '_result.txt', 'w')      
        for f, k, vf in zip(voted, self.keys, self.features):
            vf = " ".join(map(str, vf))        
            k = " ".join(map(str, k))
            print >> fd, k + " " + str(f) + " " + vf
        
        fd.close()  
        
        logging.info("plotting ROC curves")
        plt.plot(fpr,tpr,'k--', label='ROC curve (area = %0.3f)' %float(roc_auc))
        plt.title('ROC curve (area = %0.3f)' %float(roc_auc))
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.savefig(self.args.out + "_roc.png")
        #pylab.plot(fpr, tpr, 'k--', label='ROC curve (area = %0.3f)' %float(roc_auc))
        #pylab.plot([0, 1], [0, 1], 'k--')
        #pylab.xlim([0.0, 1.0])
        #pylab.ylim([0.0, 1.0])
        #pylab.xlabel('False Positive Rate')
        #pylab.ylabel('True Positive Rate')
        #pylab.savefig(args.out + "_roc.png")
    
    def cross_validate(self):
        ## NOTE:copied from Jeff/Fatemeh old code
        cv = cross_validation.StratifiedKFold(self.labels, n_folds=3)
        for i, (train, test) in enumerate(cv):
            cv_model = self.model.fit(self.features[train], self.labels[train])
            logging.info("predicting probabilities")
            probs = cv_model.predict_proba(self.features[test])
            
            voted = probs[:,1]
            fpr, tpr, thresholds = roc_curve(self.labels[test], voted)
            roc_auc = auc(fpr, tpr)
            logging.info("AUC_%d:%f" % (i, roc_auc))
            
            fd = open(self.args.out + '_result_' + str(i) + '.txt', 'w')  
            if self.args.model is None:
                for f, k in zip(voted, self.keys[test]):
                    print >> fd, k[0] + ' ' + k[1] + ' ' + k[2] + ' ' + k[3] + ' ' + k[4] + ' ' + k[5] + ' ' + str(f)
                
            else:
                ## if the model is loaded from input arguments, then the keys are not known
                for f in voted:
                    print >> fd, str(f)
                
            fd.close()
            
            logging.info("plotting ROC curves")
            plt.plot(fpr, tpr, 'k--', lw=1, label="Fold %i (AUC=%0.3f)" % (i + 1, float(roc_auc)))
    
        plt.legend(loc="lower right", numpoints=1,)
        plt.savefig(self.args.out + "_rocxval.png")
    #=================
    #Boxplot
    #=================
    def generate_boxplot(self):
        
        logging.info('Starting the plotting process')
            
        self.tot_features_name = self.__get_feature_names()
        
        self.features_vals_dict = self.__generate_feature_dict()
        self.top_features_name, self.top_features_map = self.__get_top_features()

        labels = self.__get_label_names(self.args.infiles)
        self.features_list_label= [self.__get_features_from_feature_db(self.features_vals_dict, label) for label in labels]
        
        
        self.boxplot_plot(labels)

    def __get_features_from_feature_db(self,feature_dict,label):
        features = []
        for key,value in feature_dict.iteritems():
            key = key.strip().split(';')
            label_feature_dict = value[1]
            if not label == label_feature_dict:
                continue
            features.append(eval(value[0]))
        return features

    def __get_label_names(self, reference_files):
        labels = set()
        for files in reference_files:
            file_stream = open(files.name)
            for line in file_stream:
                if line[0] == '#':
                    continue
                line = line.strip().split()
                labels.add(line[2])
        return list(labels)
        
    def __get_feature_names(self):
        if self.args.single:
            feature_set = features_single.Features()
        else:
            feature_set = features.Features()
        feature_names = feature_set.get_feature_names()
        return feature_names

        
    def __get_dict_bylabels(self,inputs,labels):
        out_dict = defaultdict(list)
        tfile = None
        nfile = None
        rfile = None
        
        for labelname in labels:
            for infilename in inputs:
                tfile = None
                nfile = None
                rfile = None
                file_stream = open(infilename.name,'r')
                for line in file_stream:
                    if line[0] == '#':
                        line = line.strip().split()
                        if self.args.single:
                            if line[1] == "tumour" or line[1] == "normal":
                                tfile = nfile = line[2]
                        else:
                            if line[1]=='tumour':
                                tfile = line[2]
                            elif line[1]=='normal':
                                nfile = line[2]
                        if line[1]=='reference':
                            rfile = line[2]
                    else:
                        l = line.strip().split()
                        if l[2] == labelname:
                            chromosome = l[0]
                            pos = l[1]
                            out_dict[labelname].append( (tfile,nfile,rfile,chromosome,pos) )
                file_stream.close()
        return out_dict
        
    def __generate_feature_dict(self):
        features_vals_dict = {}
        file_stream = open(self.args.out+'feature_db.txt','r')
        for line in file_stream:
            l = line.strip().split('\t')
            key = l[0].split()[0]
            value = l[1]
            label = l[2]
            features_vals_dict[key] = (value,label)
        return features_vals_dict
    
    def __get_top_features(self):
        top_features_names = []
        file_stream = open(self.args.out + "_importance.txt",'r')
        for line in file_stream:
            line = line.strip().split()
            top_features_names.append(line[0])
        top_features_names.reverse()
        
        top_features_names_dict =  {v:i for i, v in enumerate(self.tot_features_name) for x in top_features_names if v == x}
        
        return top_features_names,top_features_names_dict
    
    def __get_features_from_dict(self,label,f_dict,ref_dict):
        features = []
        for tfile,nfile,rfile,chromosome,pos in ref_dict.get(label):      
            key = ';'.join([rfile,nfile,tfile,chromosome,pos])
            try:
                features.append(eval(f_dict[key]))
            except KeyError:
                logging.error('error: cannot find key "%s"\n' % str(key))
        return features 
    
    def boxplot_plot(self,labels):
        pdfout = PdfPages(self.args.out+'_plots_train.pdf')
        for i,v in enumerate(self.top_features_name):
            index = self.top_features_map[v]
            if not self.tot_features_name[index] == v:
                logging.error('the feature name and feature values don\'t match')
                
            fvalue = []
            
            for feature in self.features_list_label:
                fv = []
                for p in feature:
                    fv.append(p[index])
                fvalue.append(fv)
            
            fvalue_count = [len(x) for x in fvalue]

            fig = plt.figure()
            fig.text(0.90, 0.97, 'importance:'+str(i+1), rotation='horizontal',horizontalalignment='center', verticalalignment='bottom',fontsize = 8)
            
            xlabel_description = 'Features (Count of positions) (Outliers above the plot, Outliers below the plot)'
            bplot = plt.boxplot(fvalue)
            
            fliers = []
            for i in xrange(len(bplot['boxes'])):
                fliers_above = len(bplot['fliers'][i*2]._y)
                fliers_below = len(bplot['fliers'][i*2+1]._y)
                fliers.append(str(fliers_above)+','+str(fliers_below))
                
            #xlabel_names = ['%s(%s)\n(%s)' % (labels[i], y,fliers[i]) for i,y in enumerate(fvalue_count)]
             
            try:   
                xlabel_names = ['%s(%s)\n(%s)' % (labels[i], y,fliers[i]) for i,y in enumerate(fvalue_count)]
            except:
                print '\nlabels',labels
                print '\ny',y 
                print '\nfliers',fliers
                print '\ni',i 
                print '\nfvalue_count',fvalue_count
            
            label_ylim_upper = None
            label_upper_cap =None
            label_lower_cap = None
            label_ylim_lower = None   
            for i in xrange(len(bplot['boxes'])):
                try:
                    uppercap = bplot['caps'][i*2]._y[0]
                    highestflier = max(bplot['fliers'][i*2]._y)
                    if highestflier > uppercap*100:
                        if not uppercap == 0:
                            if uppercap>label_upper_cap:
                                label_upper_cap = uppercap
                            if highestflier>label_ylim_upper:
                                label_ylim_upper = highestflier
                except:
                    #if unable to set the axis, continue without changing them
                    pass
             
                try:
                    lowercap = bplot['caps'][i*2+1]._y[0]
                    lowestflier = min(bplot['fliers'][i*2+1]._y)
                    if lowestflier > lowercap/100:
                        if not lowercap == 0:
                            if lowercap < label_lower_cap:
                                label_lower_cap = lowercap
                            if lowestflier < label_ylim_lower:
                                label_ylim_lower = lowestflier  
                except:
                    #if unable to set the axis, continue without changing them
                    pass
            
            if label_ylim_upper and label_upper_cap is not None:
                if label_ylim_upper > label_upper_cap*100:
                    label_ylim_upper = int(label_upper_cap*100)
            
            if not label_ylim_lower and label_lower_cap is not None:
                if label_ylim_lower < label_lower_cap/100:
                    label_ylim_upper = int(label_lower_cap/100)
            if label_ylim_upper is not None:
                prev_val = plt.ylim()
                if not label_ylim_upper == 0:
                    plt.ylim(prev_val[0],label_ylim_upper)
                    fig.text(0.02, 0.03,'*Upper limit (y-axis) rescaled, some datapoints are not shown', horizontalalignment='left',verticalalignment='top',fontsize=6)
            if label_ylim_lower is not None:
                prev_val = plt.ylim()
                if not label_ylim_lower == 0:
                    plt.ylim(label_ylim_lower,prev_val[1])
                    fig.text(0.02, 0.02,'*Lower limit (y-axis) rescaled, some datapoints are not shown', horizontalalignment='left',verticalalignment='top',fontsize=6)
            
            #just to show boundaries(if overlaps with axis)
            ylims = plt.ylim()
            if ylims[0] == 0:
                plt.ylim((ylims[0]-0.05*ylims[1]), (1.10*ylims[1]))
            else:
                plt.ylim( (0.90*ylims[0]), (1.10*ylims[1]) )   
                               
            plt.title(v.replace('_',' '))
            plt.ylabel('Distribution', fontsize=8)
            plt.xlabel(xlabel_description, fontsize=8)
            plt.xticks(range(1, len(xlabel_names) + 1), xlabel_names, rotation=45, fontsize=8)
            plt.tight_layout()
            pdfout.savefig(fig)
        pdfout.close()
