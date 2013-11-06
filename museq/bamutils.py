# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
import newpybam as np # new pybam
import sys
import numpy
import newfeatures
from math import log10
from collections import deque
from collections import namedtuple
from sklearn.ensemble import RandomForestClassifier
from string import Template
from datetime import datetime

mutationSeq_version = "4.0.0"
DEBUG = False

class Bam:
    def __init__(self, **kwargs):       
        self.t_bam = kwargs.get("tumour")
        self.n_bam = kwargs.get("normal")
        self.ref   = kwargs.get("reference") 
        self.rmdup = kwargs.get("rmdup")
        self.base = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
        
        ## check if duplicates need to be removed
        if  self.rmdup is None:
            self.rmdup = True
            
        ## make a pileup for normal bam
        if  self.n_bam is not None:
            self.n_pileup = np.pileup()
            self.n_pileup.open(self.n_bam)
    
        ## make a pileup for tumour bam
        if  self.t_bam is not None:
            self.t_pileup = np.pileup() 
            self.t_pileup.open(self.t_bam)
        
        if  self.ref is not None:
            self.__load_reference()
        
    def __load_reference(self):
        ## make a fasta object to laod reference
        self.fasta = np.fasta() 
        self.fasta.open(self.ref) 

    def get_normal_refnames(self):
        """ get the list of chromosome names and their corresponding IDs as a dictionary"""
        
        return dict(self.n_pileup.refnames)
    
    def get_tumour_refnames(self):
        """ get the list of chromosome names and their corresponding IDs as a dictionary"""
    
        return dict(self.t_pileup.refnames)
    
    def get_normal_chromosome_id(self, chromosome_name):
        rn = self.get_normal_refnames()
        return rn.get(chromosome_name)
    
    def get_tumour_chromosome_id(self, chromosome_name):
        rn = self.get_tumour_refnames()
        return rn.get(chromosome_name)
    
    def get_normal_chromosome_name(self, chromosome_id):
        rn = self.get_normal_refnames()
        for k in rn.keys():
            if rn[k] == chromosome_id:
                return k
        return None
        
    def get_tumour_chromosome_name(self, chromosome_id):
        rn = self.get_tumour_refnames()
        for k in rn.keys():
            if rn[k] == chromosome_id:
                return k
        return None
        
    ##TODO: change the interface to chromosome name instead of chromosome_id for FASTA object
    ## Maybe it is not a good idea since the tuples have chromosome_id and converting them to chromosome names
    ## require extra function call to get_..._chromosome_name() which is a performance over head
    def get_reference_base(self, chromosome_id, position, index=False):
        b = self.fasta.get_base(chromosome_id, int(position))
        
        if index:
            b = self.base.get(b)
            if b is None:
                b = 4
                
        return b
        
    def get_reference_sequence(self, chromosome_id, position, windowLength=500):
        return self.fasta.get_sequence(chromosome_id, position, windowLength)
        
    def get_reference_sequence_bybase(self, chromosome_id, position, windowLength=500):
        return self.fasta.get_sequence_base(chromosome_id, position, windowLength)
    
    def get_trinucleotide_context(self, chromosome_id, position):
        tc = []
        for i in range(-1,2):
            tc.append(self.get_reference_base(chromosome_id, position+i))
        return ''.join(tc)

    def get_reference_tuple(self, chromosome_id, position):
        temp_tuple = self.fasta.get_tuple(chromosome_id, int(position))
        refbase = temp_tuple[0]

        # index refbase
        refbase = self.base[refbase]

        ## replace the base with its index (required for feature extraction)
        temp_tuple = (refbase, temp_tuple[1], temp_tuple[2], temp_tuple[3], temp_tuple[4])
        return temp_tuple
            
    def get_normal_tuple(self, chromosome, position):
        self.n_pileup.set_region(chromosome, position, position)
        return self.n_pileup.get_tuple()
    
    def get_tumour_tuple(self, chromosome, position):
        self.t_pileup.set_region(chromosome, position, position)
        return self.t_pileup.get_tuple()
    
    def get_normal_tuples(self, target_positions):
        for tp in target_positions: 
            if tp[1] is None:
                self.n_pileup.set_region(tp[0])
            else:
                self.n_pileup.set_region(tp[0], tp[1], tp[2])

            while True:
                tt = self.n_pileup.get_tuple() 
                if tt is None:
                    break
                else:
                    yield tt
            
    def get_tumour_tuples(self, target_positions):
        for tp in target_positions: 
            if tp[1] is None:
                self.t_pileup.set_region(tp[0])
            else:
                self.t_pileup.set_region(tp[0], tp[1], tp[2])

            while True:
                tt = self.t_pileup.get_tuple() 
                if tt is None:
                    break
                else:
                    yield tt

    def get_pair_tuples(self, target_positions):
        for tp in target_positions:
            if tp[1] is None:
                self.t_pileup.set_region(tp[0])
                self.n_pileup.set_region(tp[0])
            else:
                self.t_pileup.set_region(tp[0], tp[1], tp[2])
                self.n_pileup.set_region(tp[0], tp[1], tp[2])
                
            while True:
                tt = self.t_pileup.get_tuple()
                nt = self.n_pileup.get_tuple()
                
                ## check for none tuples
                if not all((tt, nt)):
                    break

                ## check if the position is the same for both tumour/normal
                while tt[0] != nt[0]:
                    if tt[0] < nt[0]:
                        tt = self.t_pileup.get_tuple()
                        if tt is None:
                            break
                    else:
                        nt = self.n_pileup.get_tuple()
                        if nt is None:
                            break
                        
                if not all((tt, nt)):
                    break
                yield (tt, nt)
            
class BamUtils:
    def __init__(self, bam, args):
        self.samples = {}
        for s in args.samples:
            self.samples[s.split(':')[0]] = s.split(':')[1]
        
        ## check if it is single sample analysis as well as type of the sample
        Flags = namedtuple("Flags", "single, type")        
        if "normal" not in self.samples:
            single_flag = True
            single_type = "t"
            
        elif "tumour" not in self.samples:
            single_flag = True
            single_type = "n"
        
        else:
            single_flag = False
            single_type = None

        self.outstr_buffer = []        
        self.flags = Flags._make([single_flag, single_type])
        self.base = ['A', 'C', 'G', 'T', 'N']
        self.bam  = bam
        self.args = args

    def __parse_positions(self, positions_list):
        chromosome = positions_list.split(':')[0]
        try:
            ##check for "chr"
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
        
        if self.args.interval is not None:
            temp_tp = self.__parse_positions(self.args.interval)
            target_positions.append(temp_tp) 
        
        elif self.args.positions_file is not None:
            try:
                positions_file = open(self.args.positions_file, 'r')
                for l in positions_file.readlines():
                    temp_tp = self.__parse_positions(l.strip())
                    target_positions.append(temp_tp)
                positions_file.close()
            
            except:
                print >> sys.stderr, "\tFailed to load the positions file from " + self.args.positions_file
                sys.exit(1)

        else:
            ## get all the common chromosome names
            tcn = self.bam.get_tumour_refnames.keys() # chromosome names of tumour bam
            ncn = self.bam.get_normal_refnames.keys() # chromosome names of normal bam
            chromosome_names = set(tcn).intersection(set(ncn))

            for cn in chromosome_names:
                temp_tp = [cn, None, None]
                target_positions.append(temp_tp)
                
        return target_positions
    
    def __make_outstr(self, tt, rt, nt):
        ## flag insertions and deletions 
        if tt[-4] > 0 or tt[-2] > 0:
            filter_flag = "INDL"
        else:
            filter_flag = None
        
        ## find alternative base   
        refbase = rt[0]
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
    
    def get_features(self, tumour_tuples, normal_tuples):
        tuples_buffer   = deque()        
        features_buffer = []
        
        ## TODO: remove this, write tuples in a file
        if DEBUG:
            tuple_file = open("tuples.f", 'w')
        for tt in tumour_tuples:
            ## ignore tumour tuples with no/few variants compared to reference            
            refbase = self.bam.get_reference_base(tt[-1], tt[0], index=True)
            if tt[5][0] - tt[refbase + 1][0] < 3:
                continue
            
            ## buffer tumour tuple to campare against normal tuples
            tuples_buffer.append(tt)
            
            ##TODO: remove this
            if DEBUG:
                print >> tuple_file, tt

        ## return if all tuples were filterd
        if len(tuples_buffer) == 0:
            return []
        
        tt = tuples_buffer.popleft()
        for nt in normal_tuples:
            ## find positions where tuples for both tumour and normal exist
            while tt[0] < nt[0]:
                if len(tuples_buffer) == 0:
                    break
                tt = tuples_buffer.popleft()

            if tt[0] != nt[0]:
                continue

            ## extract reference tuples            
            rt = self.bam.get_reference_tuple(tt[-1], tt[0])

            ## calculate features      
            feature_set = newfeatures.Features(tt, nt, rt)
            tf = feature_set.get_features()
            features_buffer.append(tf)
            
            ## generate output string and buffer it
            outstr = self.__make_outstr(tt, rt, nt)
            self.outstr_buffer.append(outstr)
            
            if len(tuples_buffer) == 0:
                break
        
        ##TODO: remove this
        if DEBUG:
            tuple_file.close()
        
        ## make a numpy array required as an input to the random forest predictor
        features_buffer = numpy.array(features_buffer)
        
        ## make sure a list is returned         
        if features_buffer is None:
            return []
        else:
            return features_buffer
        
    def __fit_model(self):
        try:
            npz = numpy.load(self.samples["model"])
        except:
            print >> sys.stderr, "\tFailed to load model"
            print >> sys.stderr, sys.exc_info()[0]
            sys.exit(1)

        train  = npz["arr_1"]
        labels = npz["arr_2"]
        model  = RandomForestClassifier(random_state=0, n_estimators=1000, n_jobs=1, compute_importances=True)
        model.fit(train, labels)
        return model   
        
    def predict(self, features):
        model = self.__fit_model()
        probabilities = model.predict_proba(features)
        
        ## return only probabilities of being somatic
        probabilities = [x[1] for x in probabilities] 
        return probabilities
    
    def __meta_data(self):
        reference = self.samples.get("reference")
        tumour   = self.samples.get("tumour")
        normal   = self.samples.get("normal")
        
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
            #print "Failed to load metadata file"
            return
            
    def print_results(self, out, probabilities):
        header = self.__meta_data() 
        if header is not None:
            print >> out, header.strip()
        
        if probabilities is None:
            print >> out, "***No nominated positions***"
            return
        
        for i in xrange(len(probabilities)):
            outstr = self.outstr_buffer[i]
            p = probabilities[i]
            
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
            
            try:
                phred_quality = -10 * log10(1 - p)
            except:
                phred_quality = 99
            
            outstr = map(str, [outstr[0], outstr[1], outstr[2], self.base[outstr[3]], 
                               self.base[outstr[4]], "%.2f" % phred_quality, filter_flag, info_str])
            
            print >> out, "\t".join(outstr)
    








