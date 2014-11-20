# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 10:47:18 2013

@author: jtaghiyar
"""

from __future__ import division
import numpy
from scipy.stats import binom


        
class Features:
    def __init__(self, input_tuple=None, reference_tuple=None, sample_type = None, purity=70):
        self.it = input_tuple
        self.rt = reference_tuple
        self.type = sample_type
        self.name = "TCGA Benchmark 4 feature set with coverage info"
        self.version = "single_4.0.2"
        
        if self.it is None or self.it[5][0] == 0:
            self.it = (None, [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1]*6, 1, 1, 1, 1, 1, 1, 1, 1, [1], None)

        if self.rt is None:
            self.rt = (0, 0, 0, 0, 0)
            
        ## reference base index + 1 = index of the same base in the tumour/normal bam tuple
        self.b  = self.rt[0] + 1  
        
        ## coverage data        
        self.cd = (float(30), float(30), int(purity), float(0))
        
        ## to avoid division by zero
        self.ep = 1e-5 
        
        self.__get_nonref_index()
        self.__get_nonref_count()
        self.__snvmix()

        
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
        ("tumour_variant_depth_ratio", ((self.it[5][0] - self.it[self.b][0]) / self.it[5][0])),
        ("tumour_variant_mapq_mean", (self.it[5][2] - self.it[self.b][2]) / (self.it[5][0] - self.it[self.b][0] + self.ep)),
        
        ("probability_hom", self.hom),
        ("probability_het", self.het),
        ("probability_hom_alt", self.hom_alt),
        )
        
        self.coverage_features = (
        ("tumour_depth_coverage", self.it[5][0] / self.cd[1]),
        ("tumour_contamination", self.cd[2] / 100),
        ("whole_genome", self.cd[3])
        )
        
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


    def __get_nonref_index(self):
        nonrefbases = [x for x in range(4) if x != self.b - 1]
        max_nonrefbase_depth = 0
        nonref_index = nonrefbases[1]  
        if nonref_index == self.b:
            nonref_index = nonrefbases[2]
        
        for nb in nonrefbases:
            index = nb + 1
            tumour_nonrefbase_depth = self.it[index][0]
            if tumour_nonrefbase_depth != 0 and tumour_nonrefbase_depth > max_nonrefbase_depth:
                max_nonrefbase_depth = self.it[index][0]
                nonref_index = index
                
        self.nonref_index = nonref_index
    
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

