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
        self.refBase = None # to store reference nucleotide for a position
        
        if  self.n_bam is not None:
            self.n_pileup = np.pileup() # make a pileup for normal bam
            self.n_pileup.open(self.n_bam)
    
        if  self.t_bam is not None:
            self.t_pileup = np.pileup() # make a pileup for tumour bam
            self.t_pileup.open(self.t_bam)
        
        if  self.ref is not None:
            self.__loadReference()
        
    def __loadReference(self):
        self.fasta = np.fasta() # make a fasta object to laod reference
        self.fasta.open(self.ref) 

    def getReferenceTuple(self, chromosomeId, position):
        """get reference tuple of chromosomeId:position"""
        
        temp_tuple = self.fasta.get(chromosomeId, int(position))
        self.refBase = temp_tuple[0]

        # index ACGTN
        if  self.refBase == "A":
            self.refBase = 0
        elif self.refBase == "C":
            self.refBase = 1
        elif self.refBase == "G":
            self.refBase = 2
        elif self.refBase == "T":
            self.refBase = 3
        else:
            self.refBase = 4

        ## replace the base with its index
        temp_tuple = (self.refBase, temp_tuple[1], temp_tuple[2], temp_tuple[3], temp_tuple[4])
        return temp_tuple
            
    def getReferenceBase(self, chromosomeId, position):
        """get reference nucleotide of chromosomId:position"""

        if  self.refBase is None:
            self.refBase = self.getReferenceTuple(chromosomeId, int(position))[0]

        return self.refBase
    
    def getNormalTuple(self, chromosome, position=None):
        if position is None:
            self.n_pileup.jump(chromosome)
            while True:
                temp_tuple = self.n_pileup.next()
                if temp_tuple is None:
                    break
                else:
                    yield temp_tuple
        else:
            self.n_pileup.jump(chromosome, int(position))
            yield self.n_pileup.next()
            
    def getTumourTuple(self, chromosome, position=None):
        if position is None:
            self.t_pileup.jump(chromosome)
            while True:
                temp_tuple = self.t_pileup.next()
                if temp_tuple is None:
                    break
                else:
                    yield temp_tuple
        else:
            self.t_pileup.jump(chromosome, int(position))
            yield self.t_pileup.next()

    def getNormalRefNames(self):
        return self.n_pileup.refnames
    
    def getTumourRefNames(self):
        return self.t_pileup.refnames
        
    def getReferenceSequenceByBase(self, chromosomeId, position, windowLength=500):
        print "getReferenceSequenceByBase"
        return self.fasta.getSequenceByBase(chromosomeId, position, windowLength)
        
    def getReferenceSequence(self, chromosomeId, position, windowLength=500):
        print "getReferenceSequence"
        return self.fasta.getSequence(chromosomeId, position, windowLength)

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
        self.outstr_buffer = []
        
    def __parsePositions(self, positions_list):
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
        
    def __getTuples(self, bam_type, target_positions):
        """ iterator over tuples of target positions of a bam file"""
        
        if bam_type == "tumour":
            func = self.bam.getTumourTuple
        else:
            func = self.bam.getNormalTuple
            
        for tp in target_positions: 
            if tp.start is None:
                return func(tp.chromosome)    
            else:
                for position in xrange(tp.start, tp.stop):
                    temp_tuple = func(tp.chromosome, position)
                    return temp_tuple
    
    ## NOTE: COPIED FROM JEFF, MIGHT BE WRONG
    def __xEntropy(self, tumour_counts, normal_counts):
        total_tc = tumour_counts[4]
        total_nc = normal_counts[4]
        ent = 0 # entropy
        
        for i in xrange(4):
            base_probability_tumour = tumour_counts[i] / total_tc
            base_probability_normal = normal_counts[i] / total_nc            
            if base_probability_tumour != 0:
                if base_probability_normal == 0:
                    ent -= -7 * base_probability_tumour ## WHAT IS IT?`
                else:
                    ent -= log(base_probability_normal) * base_probability_tumour
        return ent
    
    def __extractFeature(self, tumour_tuple, reference_tuple, normal_tuple=None):
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

            features_buffer.append(self.__xEntropy(tumour_counts, normal_counts))
    
        return features_buffer
    
    def __makeOutStr(self, tt, rt, nt):
        ## flag insertions and deletions 
        if tt.insertion > 0 or tt.deletion > 0:
            filter_flag = "INDL"
        else:
            filter_flag = "PASS"
        
        ## find alternative base   
        ref_base = rt[0]
        if ref_base == tt.major:
            alt_base = tt.minor
        else:
            alt_base = tt.major
         
        ## generate informatio for the INFO column in the output
        TR = tt[ref_base + 1][0] # tumour to reference base count 
        NR = nt[ref_base + 1][0] # normal to reference base count 
        TA = tt[alt_base + 1][0] # tumour to alternative base count 
        NA = nt[alt_base + 1][0] # normal to alternative base count 
        info = [TR, TA, NR, NA, tt.insertion, tt.deletion]
        info = map(str, info) # change the int to string
        
        ## reserved to be filled later
        qual   = None # phred_quality
        out_id = "."  # ID        
        
        OutStr = namedtuple("OutStr", "chrom, pos, id, ref, alt, qual, filter, info")
        outstr = OutStr._make([tt.chromosomeId, tt.position, out_id, ref_base, alt_base, qual, filter_flag, info])
        
        return outstr             
        
    def __removeNanInf(self, features_buffer):
        ind = -1        
        for f in features_buffer:
            ind += 1
            if numpy.isnan(numpy.sum(f)) or numpy.isinf(numpy.sum(f)):
                features_buffer.pop(ind)
                self.outstr_buffer.pop(ind)
                
        return features_buffer
    
    def __fitModel(self):
        print "fitModel ..."
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
        
    def getPositions(self):
        target_positions = []
        Pos = namedtuple("Position", "chromosome, start, stop")
        
        if self.args.interval is not None:
            temp_tp = self.__parsePositions(self.args.interval)
            temp_tp = Pos._make(temp_tp)
            target_positions.append(temp_tp) 
        
        elif self.args.positions_file is not None:
            try:
                positions_file = open(self.args.positions_file, 'r')
                for l in positions_file.readlines():
                    temp_tp = self.__parsePositions(l.strip())
                    temp_tp = Pos._make(temp_tp)
                    target_positions.append(temp_tp)
                positions_file.close()
            
            except:
                print >> sys.stderr, "\tFailed to load the positions file from " + self.args.positions_file
                sys.exit(1)

        else:
            if self.flags.single:
                ## read all the chromosome ID's
                chrom_refname = self.bam.getTumourRefNames() 
            else:
                ## read all the chromosome ID's
                chrom_refname = set(self.bam.getTumourRefNames()) & set(self.bam.getNormalRefNames())    
        
            for crf in chrom_refname:
                ## return only chromosome ID's
                temp_tp = [crf, None, None]
                temp_tp = Pos._make(temp_tp)
                target_positions.append(temp_tp)
                
        return target_positions
        
    def getTumourTuples(self, target_positions):
        """ return iterator over tuples of tumour bam file"""
        
        return self.__getTuples("tumour", target_positions)
        
    def getNormalTuples(self, target_positions):
        """ return iterator over tuples of normal bam file"""
        
        return self.__getTuples("normal", target_positions)
    
    def getFeatures(self, tumour_tuples, normal_tuples):
        tuples_buffer   = deque()        
        features_buffer = []
        Tuple = namedtuple("Tuple", "position, A_tuple, C_tuple, G_tuple, T_tuple, all_tuple,\
                           major, minor, ambiguous_reads, insertion, entropy, deletion, chromosomeId")
                   
        ## read all tumour tuples           
        for tt in tumour_tuples:
            tt = Tuple._make(tt)
            tuples_buffer.append(tt)

        ## return if there are no tuples for tumour
        if len(tuples_buffer) == 0:
            return None
        
        tt = tuples_buffer.popleft()
        
        ## read all normal tuples
        for nt in normal_tuples:
            nt = Tuple._make(nt)
            
            ## find positions where tuples for both tumour and normal exist
            while tt.position < nt.position:
                if len(tuples_buffer) == 0:
                    break
                tt = tuples_buffer.popleft()

            if tt.position != nt.position:
                continue
            
            ## calculate features
            rt = self.bam.getReferenceTuple(tt.chromosomeId, tt.position)
            temp_features = self.__extractFeature(tt, rt, nt)
            features_buffer.append(temp_features)
            
            ## generate output string and buffer it
            outstr = self.__makeOutStr(tt, rt, nt)
            self.outstr_buffer.append(outstr)
            
            if len(tuples_buffer) == 0:
                break
        
        ## remove potential NaN and/or Inf values 
        features_buffer = self.__removeNanInf(features_buffer)
        
        ## make a numpy array required as an input to the random forest predictor
        features_buffer = numpy.array(features_buffer)
        
        ## make sure a list is returned         
        if features_buffer is None:
            features_buffer = 0
        
        return features_buffer

    def predict(self, features):
        model = self.__fitModel()
        probabilities = model.predict_proba(features)
        
        ## return only probabilities
        probabilities = [x[0] for x in probabilities] 
        return probabilities
    
    def printResults(self, out, probabilities=None):
        if probabilities is not None:
            for p in probabilities:
                ## outsrt: [chrom, pos, id, ref, alt, qual, filter, info]
                self.outstr_buffer.reverse()
                outstr = self.outstr_buffer.pop()
                
                ## outstr.info: [info_TR, info_TA, info_NR, info_NA, tt.insertion, tt.deletion]
                info_str = "PR=" + "%.2f" % p + ";TR=" + outstr.info[0] + \
                            ";TA=" + outstr.info[1] + ";NR=" + outstr.info[2] + \
                            ";NA=" + outstr.info[3] + ";NI=" + outstr.info[4] + ";ND=" + outstr.info[5]
                
                ## caculate phred quality
                try:
                    phred_quality = -10 * log10(1 - p)
                except:
                    phred_quality = 99
                
                outstr = map(str, [outstr.chrom, outstr.pos, outstr.id, outstr.ref, outstr.alt,
                          "%.2f" % phred_quality, outstr.filter, info_str])
                
                print >> out, "\t".join(outstr)
        
         ## print N/A if there are no candidates
#        else:
#            for i in xrange(pos.start, pos.stop):
#                print >> out, chrom + "\t" + str(i) + "\t" + "N/A\t" * 6 

    def getFeatresNames():
        pass
     
    def exportFeatures():
#==============================================================================
#         if args.export is not None or args.features_only:
#             print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " exporting features"
#             for p in xrange(u_pos-l_pos):
#                 print >> expfile, chrom + "\t" + str(l_pos + p) + "\t" + ("\t").join(map(str,batch[p]))
#==============================================================================
        pass
                
    def __printMetaData(self):
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
