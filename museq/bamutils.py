# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
#import BamClass_newPybam
import features
import Nfeatures
import features_single
import features_deep
import sys
from collections import namedtuple
from collections import deque
#from string import Template
#from datetime import datetime

class BamUtils:
    def __init__(self, bam, args):
        ## check whether it is single sample analysis and whether it is a normal("n") or tumour("t") sample
        samples = []
        for s in args.samples:
            samples.append(s.split(':')[0])

        if "normal" not in samples:
            single_flag = True
            single_type = "t"
            
        elif "tumour" not in samples:
            single_flag = True
            single_type = "n"

        else:
            single_flag = False
            single_type = None
        
        Flags = namedtuple("Flags", "single, type")
        self.flags = Flags._make([single_flag, single_type])
        self.bam  = bam
        self.args = args
        self.tuples_buffer = []
        
    def __parsePositions(self, positions_list):
        chromosome = positions_list.split('\t')[0]
        try:
            chromosome = chromosome.split('r')[1] #check if "chr" is used
        except:
            pass
        try:
            position = positions_list.split('\t')[1]
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
                position = 1
                while True:
                    temp_tuple = func(tp.chromosome, position)
                    if temp_tuple is None:
                        break
                    else:
                        position +=1 
                        yield temp_tuple
            else:
                for position in xrange(tp.start, tp.stop):
                    temp_tuple = func(tp.chromosome, position)
                    yield temp_tuple
    
    def __getAlt(ref_base, major, minor):
        if ref_base == major:
            print ref_base, major, minor
            return minor
        else:
            print ref_base, major, minor
            return major

    def getPositions(self):
        chromIds  = None
        target_positions = []
        Pos = namedtuple(["Position", "chromosome, start, stop"])
        
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
                # read all the chromosome ID's
                chromIds = self.bam.getTumourChromosomeIds() 
            else:
                # read all the chromosome ID's
                chromIds = set(self.bam.getTumourChromosomeIds()) & set(self.bam.getNormalChromosomeIds())    
        
            for cid in chromIds:
                # return only chromosome ID's
                temp_tp = [cid, None, None]
                temp_tp = Pos._make(temp_tp)
                target_positions.append(temp_tp) 
                
        return target_positions
        
    def getNormalTuples(self, target_positions):
        """ return iterator over tuples of normal bam file"""
        
        return self.__getTuples("normal", target_positions)
        
    def getTumourTuples(self, target_positions):
        """ return iterator over tuples of tumour bam file"""
        
        return self.__getTuples("tumour", target_positions)
        
    
    def getFeatures(self, tumour_tuples, normal_tuples): 
        tuples_buffer   = deque()        
        features_buffer = []
        Tuple = namedtuple("Tuple", "position, A_tuple, C_tuple, G_tuple, T_tuple, all_tuple,\
                           major, minor, ambiguous_reads, insertion, entropy, deletion, ChromosomeId")
                   
        ## read all the tumour tuples           
        for tt in tumour_tuples:
            tt = Tuple._make(tt)
            tuples_buffer.append(tt)

        if len(tuples_buffer) == 0:
            return None
        
        tt = tuples_buffer.popleft()
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
            rt = self.bam.getReferenceTuple(tt.ChromosomeId, tt.position)
            if self.flags.single:
                temp_features = self.__extractFeature(tt, rt)
            else:
                temp_features = self.__extractFeature(tt, rt, nt)
            features_buffer.append(temp_features)
            
            if len(tuples_buffer) == 0:
               break

        return features_buffer

    def __extractFeature(self, tumour_tuple, reference_tuple, normal_tuple=None):
        if self.args.normalized:
            feature_set = Nfeatures.feature_set
            coverage_features = Nfeatures.coverage_features
    #        extra_features = (("xentropy", 0), ("SENTINEL", 0)) # can be used for args.verbose
    #        version = Nfeatures.version
            coverage_data = (float(30), float(30), int(self.args.purity), float(0))
            
        elif self.args.deep: 
            feature_set = features_deep.feature_set
            coverage_features = features_deep.coverage_features
    #        extra_features = (("xentropy", 0), ("SENTINEL", 0))
    #        version = features_deep.version 
            coverage_data = (float(10000), float(10000), int(self.args.purity), float(0)) 
        
        elif self.flags.single:
            feature_set = features_single.feature_set
            coverage_features = features_single.coverage_features
    #        extra_features = (("xentropy", 0), ("SENTINEL", 0)) 
    #        version = features_single.version
            coverage_data = (float(30), float(30), int(self.args.purity), float(0))
            
        else:
            feature_set = features.feature_set
            coverage_features = features.coverage_features
    #        extra_features = (("xentropy", 0), ("SENTINEL", 0))
    #        version = features.version
            coverage_data = (float(30), float(30), int(self.args.purity), float(0))
       
        features_buffer = []
        for _, feature in feature_set:
            if self.flags.single:
                features_buffer.append(feature(tumour_tuple, reference_tuple))
            else:
                features_buffer.append(feature(tumour_tuple, normal_tuple, reference_tuple))  
    
        for _, feature in coverage_features:
            if self.flags.single:
                features_buffer.append(feature(tumour_tuple, coverage_data))
            else:
                features_buffer.append(feature(tumour_tuple, normal_tuple, coverage_data))
    
        tumour_counts = (tumour_data[1][1], tumour_data[2][1], tumour_data[3][1],
                    tumour_data[4][1], tumour_data[5][1])
        if not self.flags.single:
            normal_counts = (normal_data[1][1], normal_data[2][1], normal_data[3][1],
                        normal_data[4][1], normal_data[5][1])
            features_buffer.append(n.xentropy(n_counts, t_counts))
    
        return features_tmp
        
#    def getFeatures(self, positions):
#                ## remove NaN/Inf values from the features extracted
#        features_buffer = self.__removeNanInf(features_buffer)
        pass
    
    def getFeatresNames():
        pass
    
    def __removeNanInf():
        pass
    
    def predict():
        pass
    
#    def __printMetaData(self):
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
    
    def printResult(out, outstrs):
        OutStr = namedtuple("OutStr", "chrom, pos, ref, alt, filter, TR, TA, NR, NA, NI, ND")
        
        for outstr in outstrs:
            outstr = map(str, outstr)
            outstr = OutStr._make(outstr)
            info_str = "PR=0" + ";TR=" + outstr.TR + ";TA=" + outstr.TA +\
            ";NR=" + outstr.NR + ";NA=" + outstr.NA + ";NI=" + outstr.NI + ";ND=" + outstr.ND
            phred_qual = 0
            print >> out, outstr.chrom + "\t" + outstr.pos + "\t" + "." + "\t" + outstr.ref + "\t" +\
                  outstr.alt + "\t" + "%.2f" % phred_qual + "\t" + outstr.filter + "\t" + info_str

    
#    def makeOutStr():
#         ref_base = bam.getRefBase(tt.position)
#         alt = self.__getAlt(ref_base, tt.major, tt.minor)
#         if tt.insertion > 0 or tt.deletion > 0:
#             filter_flag = "INDL"
#         else:
#             filter_flag = "PASS"
#         info_TR = tt[ref_base + 1][0]
#         info_NR = nt[ref_base + 1][0]
#         info_TA = tt[alt + 1][0]
#         info_NA = nt[alt + 1][0]
#         outstr = OutStr._make([chrom, tt.position, bases[ref_base], bases[alt],
#                     filter_flag, info_TR, info_TA, info_NR, info_NA, tt.insertion, tt.deletion])
#         cacheOutStr(outstr_buffer, outstr)                    

    