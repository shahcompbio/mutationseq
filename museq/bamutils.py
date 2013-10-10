# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
#import BamClass_newPybam
import newpybam as np # new pybam
import features
import Nfeatures
import features_single
import features_deep
import sys
import numpy
from math import log, log10
from collections import deque
from collections import namedtuple
from sklearn.ensemble import RandomForestClassifier
#from string import Template
#from datetime import datetime

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
        self.fasta = np.fasta() # make a fasta object to laod reference
        self.fasta.open(self.ref) 

    def get_reference_base(self, chromosomeId, position, index=False):
        b = self.fasta.get_base(chromosomeId, int(position))
        
        if index:
            return self.base[b]
        return b
        
    def get_trinucleotide_context(self, chromosomeId, position):
        tc = []
        for i in range(-1,2):
            tc.append(self.get_reference_base(chromosomeId, position+i))
        return ''.join(tc)

    def get_reference_tuple(self, chromosomeId, position):
        temp_tuple = self.fasta.get_tuple(chromosomeId, int(position))
        refBase = temp_tuple[0]

        # index refBase
        refBase = self.base[refBase]

        ## replace the base with its index (required for feature extraction in the BamUtils calss)
        temp_tuple = (refBase, temp_tuple[1], temp_tuple[2], temp_tuple[3], temp_tuple[4])
        return temp_tuple
            
    def get_normal_refnames(self):
        return self.n_pileup.refnames
    
    def get_tumour_refnames(self):
        return self.t_pileup.refnames
        
    def get_reference_sequence(self, chromosomeId, position, windowLength=500):
        return self.fasta.get_sequence(chromosomeId, position, windowLength)
        
    def get_reference_sequence_bybase(self, chromosomeId, position, windowLength=500):
        return self.fasta.get_sequence_base(chromosomeId, position, windowLength)
    
    def get_normal_tuples(self, target_positions):
        for tp in target_positions: 
            if tp.start is None:
                self.n_pileup.set_region(tp.chromosome)
            else:
                self.n_pileup.set_region(tp.chromosome, tp.start, tp.stop)

            while True:
                tt = self.n_pileup.get_tuple() 
                if tt is None:
                    break
                else:
                    yield tt
            
    def get_tumour_tuples(self, target_positions):
        for tp in target_positions: 
            if tp.start is None:
                self.t_pileup.set_region(tp.chromosome)
            else:
                self.t_pileup.set_region(tp.chromosome, tp.start, tp.stop)

            while True:
                tt = self.t_pileup.get_tuple() 
                if tt is None:
                    break
                else:
                    yield tt
            
class BamUtils:
    def __init__(self, bam, args):
        self.samples = {}
        for s in args.samples:
            self.samples[s.split(':')[0]] = s.split(':')[1]
        
        ## check if it is single sample analysis as well as type of the sample
        if "normal" not in self.samples:
            single_flag = True
            single_type = "t"
            
        elif "tumour" not in self.samples:
            single_flag = True
            single_type = "n"
        
        else:
            single_flag = False
            single_type = None
        
        Flags = namedtuple("Flags", "single, type")
        self.flags = Flags._make([single_flag, single_type])
        self.bam  = bam
        self.args = args
        self.base = ['A', 'C', 'G', 'T', 'N']
        self.outstr_buffer = []
        
        ## get chromosome names and thier corresponding IDs
        self.chromosome_names = self.bam.get_tumour_refnames()
        
    def __parse_positions(self, positions_list):
        chromosome = positions_list.split(':')[0]
        try:
            chromosome = chromosome.split('r')[1] #check if "chr" is used
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
        
    ## NOTE: COPIED FROM JEFF, MIGHT BE WRONG
    def __xentropy(self, tumour_counts, normal_counts):
        total_tc = tumour_counts[4]
        total_nc = normal_counts[4]
        ent = 0 # entropy
        
        for i in xrange(4):
            base_probability_tumour = tumour_counts[i] / total_tc
            base_probability_normal = normal_counts[i] / total_nc            
            if base_probability_tumour != 0:
                if base_probability_normal == 0:
                    ent -= -7 * base_probability_tumour
                else:
                    ent -= log(base_probability_normal) * base_probability_tumour
        return ent
    
    def __extract_feature(self, tumour_tuple, reference_tuple, normal_tuple=None):
        if self.args.normalized:
            feature_set = Nfeatures.feature_set
            coverage_features = Nfeatures.coverage_features
            coverage_data = (float(30), float(30), int(self.args.purity), float(0))
    #        extra_features = (("xentropy", 0), ("SENTINEL", 0)) # can be used for args.verbose
    #        version = Nfeatures.version
            
        elif self.args.deep: 
            feature_set = features_deep.feature_set
            coverage_features = features_deep.coverage_features
            coverage_data = (float(10000), float(10000), int(self.args.purity), float(0)) 
    #        extra_features = (("xentropy", 0), ("SENTINEL", 0))
    #        version = features_deep.version 
        
        elif self.flags.single:
            feature_set = features_single.feature_set
            coverage_features = features_single.coverage_features
            coverage_data = (float(30), float(30), int(self.args.purity), float(0))
    #        extra_features = (("xentropy", 0), ("SENTINEL", 0)) 
    #        version = features_single.version
            
        else:
            feature_set = features.feature_set
            coverage_features = features.coverage_features
            coverage_data = (float(30), float(30), int(self.args.purity), float(0))
    #        extra_features = (("xentropy", 0), ("SENTINEL", 0))
    #        version = features.version
       
        features_buffer = []
        for _, featureFunc in feature_set:
            if self.flags.single:
                features_buffer.append(featureFunc(tumour_tuple, reference_tuple))
            else:
                features_buffer.append(featureFunc(tumour_tuple, normal_tuple, reference_tuple))  
    
        for _, featureFunc in coverage_features:
            if self.flags.single:
                features_buffer.append(featureFunc(tumour_tuple, coverage_data))
            else:
                features_buffer.append(featureFunc(tumour_tuple, normal_tuple, coverage_data))
    
        if not self.flags.single:
            tumour_counts = (tumour_tuple.A_tuple[0], tumour_tuple.C_tuple[0], 
                             tumour_tuple.G_tuple[0], tumour_tuple.T_tuple[0], tumour_tuple.all_tuple[0])
                         
            normal_counts = (normal_tuple.A_tuple[0], normal_tuple.C_tuple[0], 
                             normal_tuple.G_tuple[0], normal_tuple.T_tuple[0], normal_tuple.all_tuple[0])

            features_buffer.append(self.__xentropy(tumour_counts, normal_counts))
    
        return features_buffer
    
    def __make_outstr(self, tt, rt, nt):
        ## flag insertions and deletions 
        if tt.insertion > 0 or tt.deletion > 0:
            filter_flag = "INDL"
        else:
            filter_flag = None
        
        ## find alternative base   
        ref_base = rt[0]
        if ref_base == tt.major:
            alt_base = tt.minor
        else:
            alt_base = tt.major
        
        ## get tri-nucleotide context
        tc = self.bam.get_trinucleotide_context(tt.chromosomeId, tt.position)

        ## generate informatio for the INFO column in the output
        TR = tt[ref_base + 1][0] # tumour to reference base count 
        NR = nt[ref_base + 1][0] # normal to reference base count 
        TA = tt[alt_base + 1][0] # tumour to alternative base count 
        NA = nt[alt_base + 1][0] # normal to alternative base count 
        info = [TR, TA, NR, NA, tc, tt.insertion, tt.deletion]
        info = map(str, info) # change the int to string
        
        ## reserved to be filled later
        qual   = None # phred_quality
        out_id = "."  # ID   
        
        ## get chromosome name of the given chromosome ID
        chromosome_name = self.chromosome_names[tt.chromosomeId][1]
        
        OutStr = namedtuple("OutStr", "chrom, pos, id, ref, alt, qual, filter, info")
        outstr = OutStr._make([chromosome_name, tt.position, out_id, ref_base, alt_base, qual, filter_flag, info])
        
        return outstr             
        
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
        
    def get_positions(self):
        target_positions = []
        Pos = namedtuple("Position", "chromosome, start, stop")
        
        if self.args.interval is not None:
            temp_tp = self.__parse_positions(self.args.interval)
            temp_tp = Pos._make(temp_tp)
            target_positions.append(temp_tp) 
        
        elif self.args.positions_file is not None:
            try:
                positions_file = open(self.args.positions_file, 'r')
                for l in positions_file.readlines():
                    temp_tp = self.__parse_positions(l.strip())
                    temp_tp = Pos._make(temp_tp)
                    target_positions.append(temp_tp)
                positions_file.close()
            
            except:
                print >> sys.stderr, "\tFailed to load the positions file from " + self.args.positions_file
                sys.exit(1)

        else:
            for cname in self.chromosome_names: 
                ## return only chromosome ID's
                temp_tp = [cname[0], None, None]
                temp_tp = Pos._make(temp_tp)
                target_positions.append(temp_tp)
                
        return target_positions
    
    def get_features(self, tumour_tuples, normal_tuples):
        tuples_buffer   = deque()        
        features_buffer = []
        Tuple = namedtuple("Tuple", "position, A_tuple, C_tuple, G_tuple, T_tuple, all_tuple,\
                           major, minor, ambiguous_reads, insertion, entropy, deletion, chromosomeId")
                   
        ## get tumour tuples           
        for tt in tumour_tuples:
            if tt is None:
                continue
            tt = Tuple._make(tt)

            ## ignore tumour tuples with no/few variants compared to reference            
            refBase = self.bam.get_reference_base(tt.chromosomeId, tt.position, index=True)
            if tt.all_tuple[0] - tt[refBase + 1][0] < 3:
                continue
            
            ## buffer tumour tuple to campare against normal tuples
            tuples_buffer.append(tt)

        ## return if there are no tuples for tumour
        if len(tuples_buffer) == 0:
            return []
        
        tt = tuples_buffer.popleft()
        
        ## get normal tuples
        for nt in normal_tuples:
            if nt is None:
                continue
            nt = Tuple._make(nt)
            
            ## find positions where tuples for both tumour and normal exist
            while tt.position < nt.position:
                if len(tuples_buffer) == 0:
                    break
                tt = tuples_buffer.popleft()

            if tt.position != nt.position:
                continue
            
            ## extract reference tuples            
            rt = self.bam.get_reference_tuple(tt.chromosomeId, tt.position)

            ## calculate features            
            temp_features = self.__extract_feature(tt, rt, nt)

            ## check for Inf/NaN values            
            if not (numpy.isnan(sum(temp_features)) or numpy.isinf(sum(temp_features))):
                features_buffer.append(temp_features)
                
                ## generate output string and buffer it
                outstr = self.__make_outstr(tt, rt, nt)
                self.outstr_buffer.append(outstr)

            if len(tuples_buffer) == 0:
                break
        
        ## make a numpy array required as an input to the random forest predictor
        features_buffer = numpy.array(features_buffer)
        
        ## make sure a list is returned         
        if features_buffer is None:
            return []
        else:
            return features_buffer
        
    def predict(self, features):
        model = self.__fit_model()
        probabilities = model.predict_proba(features)
        
        ## return only probabilities of being somatic
        probabilities = [x[1] for x in probabilities] 
        return probabilities
    
    def print_results(self, out, probabilities=None):
        if probabilities is None:
            #for i in xrange(pos.start, pos.stop):
            #    print >> out, chrom + "\t" + str(i) + "\t" + "N/A\t" * 6            
            return
        
        for p in probabilities:
            #self.outstr_buffer.reverse()
            outstr = self.outstr_buffer.pop()
            
            if outstr.filter is None:
                if p > self.args.threshold:
                    filter_flag = "PASS"
                else:
                    filter_flag = "FAIL"
            else:
                filter_flag = outstr.filter
                
            info_str = "PR=" + "%.2f" % p + ";TR=" + outstr.info[0] + \
                        ";TA=" + outstr.info[1] + ";NR=" + outstr.info[2] + \
                        ";NA=" + outstr.info[3] + ",TC=" + outstr.info[4] + \
                        ";NI=" + outstr.info[5] + ";ND=" + outstr.info[6]
            
            try:
                phred_quality = -10 * log10(1 - p)
            except:
                phred_quality = 99
            
            outstr = map(str, [outstr.chrom, outstr.pos, outstr.id, self.base[outstr.ref], 
                               self.base[outstr.alt], "%.2f" % phred_quality, filter_flag, info_str])
            
            print >> out, "\t".join(outstr)
    
    def get_features_names():
        pass
     
    def export_features():
#==============================================================================
#         if args.export is not None or args.features_only:
#             print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " exporting features"
#             for p in xrange(u_pos-l_pos):
#                 print >> expfile, chrom + "\t" + str(l_pos + p) + "\t" + ("\t").join(map(str,batch[p]))
#==============================================================================
        pass
                
    def __print_meta_data(self):
#        try:
#            cfg_file = open(self.args.config, 'r')
#            tmp_file = ""
#            for l in cfg_file:
#                l = Template(l).substitute(DATETIME=datetime.now().strftime("%Y%m%d"),
#                                           VERSION=mutationSeq_version,
#                                           REFERENCE="N\A", #samples["reference"],
#                                           TUMOUR=samples["tumour"],
#                                           NORMAL=samples["normal"],
#                                           THRESHOLD="N/A"#args.threshold
#                                           )
#                tmp_file += l
#            cfg_file.close()
#            print >> out, tmp_file,
#        except:
#            warn("Failed to load metadata file")
        pass

#    def __removeNanInf(self, features_buffer):
#        ind = -1        
#        for f in features_buffer:
#            ind += 1
#            if numpy.isnan(numpy.sum(f)) or numpy.isinf(numpy.sum(f)):
#                features_buffer.pop(ind)
#                self.outstr_buffer.pop(ind)
#                
#        return features_buffer