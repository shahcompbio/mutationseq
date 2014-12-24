'''
Created on May 20, 2014

@author: dgrewal
'''
import logging
import pybamapi
import features_single
import re
from scipy.stats import binom
from math import log10
from sklearn.externals import joblib
import sys
import numpy
import resource
from string import Template
from datetime import datetime

mutationSeq_version = "4.3.3"

class PreProcess(object):
    def __init__(self,args):
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
 
        self.features_module = features_single
        rmdups = True

        logging.info("initializing the NormalBam")
        self.bam = pybamapi.Bam(bam=self.samples.get("normal"), reference=self.ref, coverage=self.args.coverage,
                                  rmdups=rmdups, qual_threshold=self.args.quality_threshold)
        
        logging.info("initializing the TumourBam")
        self.t_bam = pybamapi.Bam(bam=self.samples.get("tumour"), reference=self.ref, coverage=self.args.coverage,
                                  rmdups=rmdups, qual_threshold=self.args.quality_threshold)
        
        
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
    
    def __filter_positions(self,temp_tp, pch = ':'):
        if self.args.interval:
            int_pos = self.__parse_positions(self.args.interval, pch)          
            if int_pos[0] == temp_tp[0]:
                #if only chr is provided
                if None in int_pos:
                    return temp_tp
                start = max(int_pos[1],temp_tp[1])
                stop = min(int_pos[2],temp_tp[2])
                if start > stop:
                    return None
                return [int_pos[0],start,stop]
                
                
        
    def get_positions(self, pch=':'):
        target_positions = []
        
        logging.info("parsing the position_file")
        try:
            positions_file = open(self.args.positions_file, 'r')
            for l in positions_file.readlines():
                temp_tp = self.__parse_positions(l.strip(), pch)
                if self.args.interval:
                    int_pos = self.__parse_positions(self.args.interval,pch)
                temp_tp = self.__filter_positions(temp_tp)
                if temp_tp is not None:
                    target_positions.append(temp_tp)
            positions_file.close()
        
        except:
            logging.error("error: failed to load the positions file " + self.args.positions_file)
            raise Exception("failed to load the positions file")
        
        return target_positions

    def __get_genotype(self,nonref_count,count_all):
        aa = 0.01
        ab = 0.50
        bb = 0.99
        prob = [aa, ab, bb]
        
        ## binomial pmf for three different probabilities
        binom_val = [binom.pmf(nonref_count, count_all, p) for p in prob]
        
        if sum(binom_val)==0:
            binom_val = [val/(sum(binom_val)+1e-150) for val in binom_val]
        else:
            binom_val = [val/sum(binom_val) for val in binom_val]
        
        pr_aa = round(binom_val[0],4)
        pr_ab = round(binom_val[1],4)
        pr_bb = round(binom_val[2],4)
        
        if pr_aa == max(pr_aa,pr_ab,pr_bb):
            gt = '0/0'
        if pr_ab == max(pr_aa,pr_ab,pr_bb):
            gt = '0/1'
        if pr_bb == max(pr_aa,pr_ab,pr_bb):
            gt = '1/1'
        
        if pr_aa == 0:
            pr_aa = 255
        elif pr_aa == 1:
            pr_aa = 0
        else:
            pr_aa = -10*log10(pr_aa) 
        
        if pr_ab == 0:
            pr_ab = 255
        elif pr_ab == 1:
            pr_ab = 0
        else:
            pr_ab = -10*log10(pr_ab) 
            
        if pr_bb == 0:
            pr_bb = 255
        elif pr_bb == 1:
            pr_bb = 0
        else:
            pr_bb = -10*log10(pr_bb) 
            
        return int(pr_aa),int(pr_ab),int(pr_bb),gt
    
    def __get_allele_counts(self,chromosome,position,chr_id):
        # get the tumour tuple
        tt = self.t_bam.get_tuple(chromosome, position)
        rt = self.t_bam.get_reference_tuple(chr_id, position)
        refbase = rt[0]
        
        #if tumour tuple is None
        if tt == None:
            return 0,0
        
         
        cov = tt[5][0]
        if cov == 0:
            altbase = 'N/A'
            TR = 0
            TA = 0
        else:
            major = tt[6]
            minor = tt[7]
            if refbase == major:
                altbase = minor
            else:
                altbase = major
                
            TR = tt[refbase + 1][0]
            TA = tt[altbase + 1][0]
                
        return TR,TA
 
    
    def __make_outstr(self, nt, refbase):
        n_coverage = nt[5][0]
        
        ## tumour information to print
        if n_coverage == 0:
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
            major = nt[6]
            minor = nt[7]
            if refbase == major:
                altbase = minor
        
            else:
                altbase = major
                  
            ## tumour counts
            NR = nt[refbase + 1][0] # tumour to reference base count 
            NA = nt[altbase + 1][0] # tumour to alternative base count

            ## indel 
            insertion = nt[-4]
            deletion  = nt[-2]
                        
            ## filter flag
            if deletion > 0 or insertion > 0:
                filter_flag = "INDL"

            else:
                filter_flag = None
                    
        ## tri_nucleotide context
        chromosome_id = nt[-1]
        position = nt[0]
        tc = self.bam.get_trinucleotide_context(chromosome_id, position)

        ## generate information for the INFO column in the output vcf               
        pr_aa,pr_ab,pr_bb,gt = self.__get_genotype(NA,nt[5][0])
        
        chromosome = self.bam.get_chromosome_name(chromosome_id)
        TR,TA = self.__get_allele_counts(chromosome, position, chromosome_id)

        ## generate information for the INFO column in the output vcf               
        info = [TR, TA, NR, NA, tc, insertion, deletion, gt, pr_aa, pr_ab, pr_bb]
            
        info = map(str, info) 

        ## reserved for database ID, to be filled later
        out_id = "."  
        
        ## get chromosome name of the given chromosome ID
        chromosome_name = self.bam.get_chromosome_name(chromosome_id)
        
        ##remove chr from chromosome
        chromosome_name = chromosome_name.replace('chr','')
        
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
        for it in tuples: 
            chromosome_id = it[-1]
            position = it[0]
          
            ## get corresponding reference tuple
            rt = self.bam.get_reference_tuple(chromosome_id, position)
            
            #MUT-238 If the ref base is 4(N) ignore the position
            if rt[0] >= 4: 
                logging.error(str(position)+' position references base N and has been ignored')
                continue
            
            ## calculate features     
            feature_set = self.features_module.Features(it, rt,'n')
            temp_feature = feature_set.get_features()
            
            ## generate output string and buffer it
            outstr = self.__make_outstr(it, rt[0])
            
            ##add the features to buffer if genotype is 0/1
            if outstr[-1][7] == '0/1':
                self.features_buffer.append(temp_feature)
                self.outstr_buffer.append(outstr)
            
            ## check the buffer size and flush
            if len(self.features_buffer) >= self.buffer_size:
                yield self.__flush()
        
        yield self.__flush()
        
        
    def __load_model(self):
        try:
            logging.info("loading model")
            return joblib.load(self.model)
        
        except:
            logging.error("error: failed to load model")
            raise Exception("failed to load model")
        
    def __verify_model_features(self,model):
        features = self.features_module.Features()
        if not model.version == features.version:
            return False
        if not model.name == features.name:
            return False
        return True
        
        
    def predict(self, features_outstrs):
        #model = self.__fit_model()
        model = self.__load_model()
        
        #verify the model against the features
        if not self.__verify_model_features(model):
            logging.error('The features and the model do not match')
            raise Exception('mismatched model')
        
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
                
                ## do not print positions with p < threshold
                if p < self.args.threshold:
                    continue
                
                ## only print heterozygous positions
                if outstr[-1][7] != '0/1':
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
                            ";NI=" + outstr[-1][5] + ";ND=" + outstr[-1][6] + \
                            ";GT=" + outstr[-1][7] +";PL="+outstr[-1][8]+\
                           ','+outstr[-1][9]+','+outstr[-1][10]
                
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
