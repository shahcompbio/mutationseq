"""
created on Mon Nov 18 10:26:34 2013

@author: jtaghiyar
"""

from __future__ import division
import numpy
from math import log
from scipy.stats import binom
from scipy.special import gammaln
import logging


class Features:
    def __init__(self, input_tuple=None, reference_tuple=None, input_bg=None, reference_bg=None, sample_type = None, indexes = None, purity=70):
        base_index = {'A':1,'C':2,'G':3,'T':4}

        self.name = "TCGA Benchmark 4 feature set with coverage info"
        self.version = "deep_single_0.1"

        self.it = input_tuple
        self.rt = reference_tuple
        
        self.it_bg = input_bg
        self.rt_bg = reference_bg
        
        self.type = sample_type
        
        if indexes is None or indexes[0] == 'N/A' or indexes[1] == 'N/A':
            self.man_ref = None
            self.man_nonref = None
        else:
            self.man_ref = base_index.get(indexes[0])
            self.man_nonref = base_index.get(indexes[1])
            
        ## check for zero coverage or None tuple
        if self.it is None or self.it[5][0] == 0:
            self.it = (None, [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1]*6, 1, 1, 1, 1, 1, 1, 1, 1, [1], None)
        
        if self.rt is None:
            self.rt = (0, 0, 0, 0, 0)

        if self.it_bg is None:
            self.it_bg = [[1, [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1]*6, 1, 1, 1, 1, 1, 1, None]]*10

        if self.rt_bg is None:
            self.rt_bg = [[1,(1, 0, 0, 0, 0)]]*10
            
        ## reference base index + 1 = index of the same base in the tumour/normal bam tuple
        self.b  = self.rt[0] + 1  
        
        ## coverage data        
        self.cd = (float(30), float(30), int(purity), float(0))
        
        ## to avoid division by zero
        self.ep = 1e-5 

        ## get the the count of nonref bases from both tumour/normal tuples. Get only for the variant that is common between both and has the maximum depth.
        self.nonref_index = self.__get_nonref_index(self.b,self.it)
        self.__get_nonref_count()
        self.__snvmix()

        self.__background_var_freq_mean()
        
#=============================================================================
# features
#=============================================================================
        self.feature_set = (
        ("tumour_indels", self.it[9] / self.it[5][0]),
        #("tumour_depth", self.it[5][0]),
        ("tumour_entropy", self.it[10]),
        
        
        ("tumour_mapq_ratio", self.it[5][2] / self.it[5][0]),
        ("tumour_quality_ratio", self.it[5][1] / self.it[5][0]),
        ("tumour_distance_ratio", self.it[5][3] / self.it[5][0]),
        ("tumour_direction_ratio", self.it[5][4] / self.it[5][0]),
        #("tumour_tail_distance", self.it[5][3] / self.it[5][0]),
                
                        
        ("tumour_ref_depth", self.it[self.b][0] / self.it[5][0]),
        ("tumour_ref_quality", self.it[self.b][1] / (self.it[self.b][0]+ self.ep)),
        ("tumour_ref_mapq", self.it[self.b][2] / (self.it[self.b][0]+ self.ep)), 
        ("tumour_ref_direction_total", self.it[self.b][4] / self.it[5][0]),
        ("tumour_ref_direction", self.it[self.b][4] / (self.it[self.b][0] + self.ep)),
        
        
        ("region_entropy", self.rt[4]),
        ("region_gc_content", self.rt[3]),
        ("homopolymer_f", self.rt[1]),
        ("homopolymer_b", self.rt[2]),
        
        ("tumour_variant_quality_ratio", ((self.it[5][1] - self.it[self.b][1]) / (self.it[5][0] - self.it[self.b][0] + self.ep))),
        #("tumour_variant_quality", (self.it[5][1] - self.it[self.b][1])),
        ("tumour_variant_direction_ratio", (self.it[5][4] - self.it[self.b][4]) / (self.it[5][0] - self.it[self.b][0] + self.ep)),
        ("tumour_variant_distance", (self.it[5][3] - self.it[self.b][3]) / (self.it[5][0] - self.it[self.b][0] + self.ep)),
        #("tumour_variant_depth_ratio", ((self.it[5][0] - self.it[self.b][0]) / self.it[5][0])),
        ("tumour_variant_mapq_mean", (self.it[5][2] - self.it[self.b][2]) / (self.it[5][0] - self.it[self.b][0] + self.ep)),
        
        ("probability_hom", self.hom),
        ("probability_het", self.het),
        ("probability_hom_alt", self.hom_alt),
        ("background_var_freq_mean",self.it_bg_avg_var_freq),
        ("p_value",self.p_val)
        )
        
        self.coverage_features = (
        ("tumour_depth_coverage", self.it[5][0] / self.cd[1]),
        ("tumour_contamination", self.cd[2] / 100),
        ("whole_genome", self.cd[3])
        )

    def __get_nonref_index(self,ref,it):
        nonrefbases = [x for x in range(4) if x != ref - 1]
        max_nonrefbase_depth = 0

        nonref_index = nonrefbases[1]  
        if nonref_index == ref:
            nonref_index = nonrefbases[2]

        for nb in nonrefbases:
            index = nb + 1
            nonrefbase_depth = it[index][0]
            if nonrefbase_depth != 0 and nonrefbase_depth > max_nonrefbase_depth:
                max_nonrefbase_depth = it[index][0]
                nonref_index = index
                
        return nonref_index
    
    def __get_nonref_count(self):
        self.it_nonref_count = self.it[self.nonref_index][0]

    ##TODO: eventually this function should be removed    
    def __isvalid(self, x):
        if numpy.isnan(x) or numpy.isinf(x):
            
            ##TODO: remove this line
            print "NaN"
            return False
        return True
    
    def get_features(self):
        features = []
        for _, f in self.feature_set:
            if self.__isvalid(f):
                features.append(f)
            else:
                features.append(0)
       
        for _,f in self.coverage_features:
            if self.__isvalid(f):
                features.append(f)
            else:
                features.append(0)
        
        return features
       
    def get_feature_names(self):
        feature_names = []
        for n, _ in self.feature_set:
            feature_names.append(n)

        for n, _ in self.coverage_features:
            feature_names.append(n)
        
        return feature_names

    def __snvmix(self):
        ## normal allele frequency
        aa = 0.01
        ab = 0.50
        bb = 0.99
        n_prob = [aa, ab, bb]
        
        ## tumour allele frequency
        aa = 0.01
        ab = 0.30
        bb = 0.90
        t_prob = [aa, ab, bb]
        
        ## binomial pmf for three different probabilities
        if self.type == 't':
            binom_val = [binom.pmf(self.it_nonref_count, self.it[5][0], p) for p in t_prob]
        else:
            binom_val = [binom.pmf(self.it_nonref_count, self.it[5][0], p) for p in n_prob]
        binom_val = [val/(sum(binom_val)+self.ep) for val in binom_val]
        
        self.hom = binom_val[0]
        self.het = binom_val[1]
        self.hom_alt = binom_val[2]
                
    def __background_var_freq_mean(self):
        it_bg_avg_var_freq = []
        
        for it in self.it_bg:
            rt = [val for val in self.rt_bg if val[0] == it[0]]
            
            if rt == []:
                raise Exception('Error: reference tuple for position '+ str(it[0])+ ' couldn\'t be extracted')
            
            ref = rt[0][1][0] + 1
            nonref_index = self.__get_nonref_index(ref,it)
                
            it_bg_avg_var_freq.append(it[nonref_index][0]/ (it[nonref_index][0]+it[ref][0]) )
            
        try:
            self.it_bg_avg_var_freq = sum(it_bg_avg_var_freq)/len(it_bg_avg_var_freq)
        except Exception,e:
            logging.error('Zero division error due to error in retreiving tuples (See MUT-255) for position: '+str(self.tt[0]))
            self.it_bg_avg_var_freq = float('nan')
 
        if self.man_ref:
            if self.b != self.man_ref:
                logging.error('Error: ref base in manifest and tuples don\'t match')
                raise Exception('ref base in manifest and tuples don\'t match')
            
            nref_count = self.it[self.man_nonref][0]
            depth = self.it[self.man_nonref][0]+self.it[self.b][0]
        else:
            nref_count = self.it_nonref_count
            depth = self.it[self.nonref_index][0]+self.it[self.b][0]

        self.p_val = binom.sf(nref_count-1,depth,self.it_bg_avg_var_freq)

