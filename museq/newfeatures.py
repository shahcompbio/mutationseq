# -*- codingutf-8 -*-
"""
Created on Wed Oct 16 11:31:56 2013

@author: jtaghiyar
"""
import numpy
from math import log

class Features:
    def __init__(self, tumour_tuple, normal_tuple, reference_tuple, purity=70):
        self.name = "TCGA Benchmark 4 featureset with coverage info"
        self.version = "5"
        self.tt = tumour_tuple
        self.nt = normal_tuple
        self.rt = reference_tuple
        self.cd = (float(30), float(30), int(purity), float(0))
        
        self.feature_set = (
        ("tumour_indels", self.tt[9] / self.tt[5][0]),
        ("normal_indels", self.nt[9] / self.nt[5][0]),
        
        ("tumour_ref_depth", self.tt[self.rt[0] + 1][0] / self.tt[5][0]),
        ("normal_ref_depth", self.nt[self.rt[0] + 1][0] / self.nt[5][0]),
        ("normal_mapq_ratio", self.nt[5][2] / self.nt[5][0]),
        ("tumour_mapq_ratio", self.tt[5][2] / self.tt[5][0]),
        ("normal_ref_quality", self.nt[self.rt[0] + 1][1] / self.nt[5][0]),
        ("tumour_ref_quality", self.tt[self.rt[0] + 1][1] / self.tt[5][0]),
        ("normal_quality_ratio", self.nt[5][1] / self.nt[5][0]),
        ("tumour_quality_ratio", self.tt[5][1] / self.tt[5][0]),
        
        ("normal_tumour_quality", (self.tt[5][1] / self.tt[5][0]) / ((self.nt[5][1] / self.nt[5][0]) + 0.00001)),
        
        ("normal_tumour_mapq", (self.tt[5][2] / self.tt[5][0]) / ((self.nt[5][2] / self.nt[5][0]) + 0.00001)),
        ("normal_tumour_ref_depth", ((self.tt[self.rt[0] + 1][1] / self.tt[5][0]) + 0.00001) / ((self.nt[self.rt[0] + 1][1] / self.nt[5][0]) + 0.00001)),
        ("normal_tumour_direction", (self.tt[5][4] / self.tt[5][0]) / ((self.nt[5][4] / self.nt[5][0]) + 0.00001)),
        ("normal_tumour_ref_direction", (self.tt[self.rt[0] + 1][4] / (self.tt[self.rt[0] + 1][0] + 0.00001)) / ((self.nt[self.rt[0] + 1][4] / (self.nt[self.rt[0] + 1][0] + 0.00001)) + 0.00001)),
        ("normal_distance_ratio", self.nt[5][3] / self.nt[5][0]),
        ("tumour_distance_ratio", self.tt[5][3] / self.tt[5][0]),
        
        ("normal_direction_ratio", self.nt[5][4] / self.nt[5][0]),
        ("tumour_direction_ratio", self.tt[5][4] / self.tt[5][0]),
        
        ("normal_ref_direction_total", self.nt[self.rt[0] + 1][4] / self.nt[5][0]),
        ("tumour_ref_direction_total", self.tt[self.rt[0] + 1][4] / self.tt[5][0]),
        
        ("normal_ref_direction", self.nt[self.rt[0] + 1][4] / (self.nt[self.rt[0] + 1][0] + 0.00001)),
        ("tumour_ref_direction", self.tt[self.rt[0] + 1][4] / (self.tt[self.rt[0] + 1][0] + 0.00001)),
        
        
        ("region_entropy", self.rt[4]),
        ("region_gc_content", self.rt[3]),
        ("homopolymer_f", self.rt[1]),
        ("homopolymer_b", self.rt[2]),
        
        ("tumour_variant_quality_ratio", ((self.tt[5][1] - self.tt[self.rt[0] + 1][1]) / (self.tt[5][0] - self.tt[self.rt[0] + 1][0] + 0.00001))),
        ("normal_variant_quality_ratio", ((self.nt[5][1] - self.nt[self.rt[0] + 1][1]) / (self.nt[5][0] - self.nt[self.rt[0] + 1][0] + 0.00001))),
        ("tumour_variant_quality", (self.tt[5][1] - self.tt[self.rt[0] + 1][1])),
        ("normal_variant_quality", (self.nt[5][1] - self.nt[self.rt[0] + 1][1])),
        
        
        #("normal_direction", self.nt[5][4]),
        #("tumour_direction", self.tt[5][4]),
        
        ("normal_variant_direction_ratio", (self.nt[5][4] - self.nt[self.rt[0] + 1][4]) / (self.nt[5][0] - self.nt[self.rt[0] + 1][0] + 0.00001)),
        ("tumour_variant_direction_ratio", (self.tt[5][4] - self.tt[self.rt[0] + 1][4]) / (self.tt[5][0] - self.tt[self.rt[0] + 1][0] + 0.00001)),
        
        #("normal_variant_direction", self.nt[5][4] - self.nt[self.rt[0] + 1][4]),
        #("tumour_variant_direction", self.tt[5][4] - self.tt[self.rt[0] + 1][4]),
        
        #("normal_tumour_variant_ratio", (self.tt[5][4] - self.tt[self.rt[0] + 1][4]) / (self.nt[5][4] - self.nt[self.rt[0] + 1][4] + 0.00001)),
        
        #("normal_tumour_variant_quality", (self.tt[5][1] - self.tt[self.rt[0] + 1][1]) / (self.nt[5][1] - self.nt[self.rt[0] + 1][1] + 0.00001)),
        
        #("normal_tumour_variant_distance", (self.tt[5][3] - self.tt[self.rt[0] + 1][3]) / (self.nt[5][3] - self.nt[self.rt[0] + 1][3] + 0.00001)),
        
        #("normal_quality", self.nt[5][1]),
        #("tumour_quality", self.tt[5][1]),
        
        ("normal_tail_distance", self.nt[5][3] / self.nt[5][0]),
        ("tumour_tail_distance", self.tt[5][3] / self.tt[5][0]),
        
        
        #("normal_tumour_variant_mapq", (self.tt[5][2] - self.tt[self.rt[0] + 1][2]) / (self.nt[5][2] - self.nt[self.rt[0] + 1][2] + 0.00001)),
        ("normal_variant_distance", (self.nt[5][3] - self.nt[self.rt[0] + 1][3]) / (self.nt[5][0] - self.nt[self.rt[0] + 1][0] + 0.00001)),
        ("tumour_variant_distance", (self.tt[5][3] - self.tt[self.rt[0] + 1][3]) / (self.tt[5][0] - self.tt[self.rt[0] + 1][0] + 0.00001)),
        
        #("normal_minor_allele", (self.rt[0] != self.tt[7] and self.nt[self.tt[6] + 1][0] > 0) or (self.rt[0] != self.tt[6] and self.nt[self.tt[7] + 1][0] > 0)),
        
        
        #==============================================================================
        ## Classic features
        #==============================================================================
        ("normal_depth", self.nt[5][0]),
        ("tumour_depth", self.tt[5][0]),
        
        #("normal_variant_depth", (self.nt[5][0] - self.nt[self.rt[0] + 1][0])),
        #("tumour_variant_depth", (self.tt[5][0] - self.tt[self.rt[0] + 1][0])),
        ("normal_variant_depth_ratio", ((self.nt[5][0] - self.nt[self.rt[0] + 1][0]) / self.nt[5][0])),
        ("tumour_variant_depth_ratio", ((self.tt[5][0] - self.tt[self.rt[0] + 1][0]) / self.tt[5][0])),
        
        ("normal_tumour_depth", (self.tt[5][0] / self.nt[5][0])),
        #("normal_tumour_variant_depth", (self.tt[5][0] - self.tt[self.rt[0] + 1][0]) / (self.nt[5][0] - self.nt[self.rt[0] + 1][0] + 0.00001)),
        ("normal_tumour_variant_depth_ratio", ((self.tt[5][0] - self.tt[self.rt[0] + 1][0]) / self.tt[5][0]) / (((self.nt[5][0] - self.nt[self.rt[0] + 1][0]) / self.nt[5][0]) + 0.00001)),
        
        ("tumour_entropy", self.tt[10]),
        ("normal_entropy", self.nt[10]),
        ("normal_tumour_entropy", self.nt[10] / (self.tt[10] + 0.00000001)),
        
        ("normal_variant_mapq_mean", (self.nt[5][2] - self.nt[self.rt[0] + 1][2]) / (self.nt[5][0] - self.nt[self.rt[0] + 1][0] + 0.00001)),
        ("tumour_variant_mapq_mean", (self.tt[5][2] - self.tt[self.rt[0] + 1][2]) / (self.tt[5][0] - self.tt[self.rt[0] + 1][0] + 0.00001)),
        
        ("normal_tumour_variant_direction_ratio", ((self.tt[5][4] - self.tt[self.rt[0] + 1][4]) / (self.tt[5][0] - self.tt[self.rt[0] + 1][0] + 0.00001)) / ((self.nt[5][4] - self.nt[self.rt[0] + 1][4]) / (self.nt[5][0] - self.nt[self.rt[0] + 1][0] + 0.00001) + 0.00001)),
        ("normal_tumour_variant_mapq_ratio", ((self.tt[5][2] - self.tt[self.rt[0] + 1][2]) / (self.tt[5][0] - self.tt[self.rt[0] + 1][0] + 0.00001)) / ((self.nt[5][2] - self.nt[self.rt[0] + 1][2]) / (self.nt[5][0] - self.nt[self.rt[0] + 1][0] + 0.00001) + 0.00001)),
        ("normal_tumour_variant_quality_ratio", ((self.tt[5][1] - self.tt[self.rt[0] + 1][1]) / (self.tt[5][0] - self.tt[self.rt[0] + 1][0] + 0.00001)) / ((self.nt[5][1] - self.nt[self.rt[0] + 1][1]) / (self.nt[5][0] - self.nt[self.rt[0] + 1][0] + 0.00001) + 0.00001)),
        ("normal_tumour_variant_distance_ratio", ((self.tt[5][3] - self.tt[self.rt[0] + 1][3]) / (self.tt[5][0] - self.tt[self.rt[0] + 1][0] + 0.00001)) / ((self.nt[5][3] - self.nt[self.rt[0] + 1][3]) / (self.nt[5][0] - self.nt[self.rt[0] + 1][0] + 0.00001) + 0.00001)),
        
        
        ("normal_tumour_direction_ratio", (self.tt[5][4] / self.tt[5][0]) / ((self.nt[5][4] / self.nt[5][0]) + 0.00001)),
        ("normal_tumour_mapq_ratio", (self.tt[5][2] / self.tt[5][0]) / ((self.nt[5][2] / self.nt[5][0]) + 0.00001)),
        ("normal_tumour_distance_ratio", (self.tt[5][3] / self.tt[5][0]) / ((self.nt[5][3] / self.nt[5][0]) + 0.00001)),
        ("normal_tumour_quality_ratio", (self.tt[5][1] / self.tt[5][0]) / ((self.nt[5][1] / self.nt[5][0]) + 0.00001)),
        
        
        #==============================================================================
        ## Unscaled features
        #==============================================================================
        #("tumour_mapq", self.tt[5][2]),
        #("normal_mapq", self.nt[5][2]),
        #("normal_variant_mapq", (self.nt[5][2] - self.nt[self.rt[0] + 1][2])),
        #("tumour_variant_mapq", (self.tt[5][2] - self.tt[self.rt[0] + 1][2])),
        )
        
        self.coverage_features = (
        ("normal_depth_coverage", self.nt[5][0] / self.cd[0]),
        ("tumour_depth_coverage", self.tt[5][0] / self.cd[1]),
        ("tumour_contamination", self.cd[2] / 100),
        ("whole_genome", self.cd[3])
        )
    
    def __isvalid(self, x):
        if numpy.isnan(x) or numpy.isinf(x):
            return False
        return True
    
    def __xentropy(self):
        tumour_counts = (self.tt[1][0], self.tt[2][0], self.tt[3][0], self.tt[4][0], self.tt[5][0])
        normal_counts = (self.nt[1][0], self.nt[2][0], self.nt[3][0], self.nt[4][0], self.nt[5][0])
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
        
        features.append(self.__xentropy())
        return features
       
    def get_feature_names(self):
        feature_names = []
        for n, _ in self.feature_set:
            feature_names.append(n)

        for n, _ in self.coverage_features:
            feature_names.append(n)
        
        return feature_names








    
    
        
