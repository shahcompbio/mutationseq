# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 10:26:34 2013

@author: jtaghiyar
"""

from __future__ import division
import numpy
from math import log
from scipy.stats import binom
from scipy.special import gammaln
import logging

class Features:
    def __init__(self, tumour_tuple=None, normal_tuple=None, reference_tuple=None, tumour_bg=None, normal_bg=None, reference_bg=None, purity=70):
        self.name = "TCGA Benchmark 4 featureset with coverage info"
        self.version = "5_deep"

        self.tt = tumour_tuple
        self.nt = normal_tuple
        self.rt = reference_tuple
        
        self.tt_bg = tumour_bg
        self.nt_bg = normal_bg
        self.rt_bg = reference_bg
        
        ## check for zero coverage or None tuple
        if self.tt is None or self.tt[5][0] == 0:
            self.tt = (None, [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1]*6, 1, 1, 1, 1, 1, 1, None)

        if self.nt is None or self.nt[5][0] == 0: 
            self.nt = (None, [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1]*6, 1, 1, 1, 1, 1, 1, None)
        
        if self.rt is None:
            self.rt = (0, 0, 0, 0, 0)
            
        ## reference base index + 1 = index of the same base in the tumour/normal bam tuple
        self.b  = self.rt[0] + 1  
        
        ## coverage data        
        self.cd = (float(30), float(30), int(purity), float(0))
        
        ## to avoid division by zero
        self.ep = 1e-5 

        ## get the the count of nonref bases from both tumour/normal tuples. Get only for the variant that is common between both and has the maximum depth.
        self.nonref_index = self.__get_tumour_nonref_index(self.b,self.tt)
        self.__get_tumour_nonref_count()
        self.__get_normal_nonref_count()
        self.__jointsnvmix()
        self.__background_var_freq_mean()
        
#=============================================================================
# features
#=============================================================================
        self.feature_set = (
        ("homopolymer_b", self.rt[2]),
        ("homopolymer_f", self.rt[1]),
        ("normal_direction_ratio", self.nt[5][4] / self.nt[5][0]),
        ("normal_distance_ratio", self.nt[5][3] / self.nt[5][0]),
        ("normal_entropy", self.nt[10]),
        ("normal_indels", (self.nt[9]+self.nt[11]) / self.nt[5][0]),
        ("normal_mapq_ratio", self.nt[5][2] / self.nt[5][0]),
        ("normal_quality_ratio", self.nt[5][1] / self.nt[5][0]),
        ("normal_ref_depth", self.nt[self.b][0] / (self.nt[self.b][0] + self.ep)),
        ("normal_ref_direction", self.nt[self.b][4] / (self.nt[self.b][0] + self.ep)),
        ("normal_ref_direction_total", self.nt[self.b][4] / self.nt[5][0]),
        ("normal_ref_quality", self.nt[self.b][1] / (self.nt[self.b][0] + self.ep)),#self.nt[5][0]),
        ("normal_tumour_depth", (self.tt[5][0] / self.nt[5][0])),
        ("normal_tumour_direction", (self.tt[5][4] / self.tt[5][0]) / ((self.nt[5][4] / self.nt[5][0]) + self.ep)),
        ("normal_tumour_entropy", self.nt[10] / (self.tt[10] + self.ep)),       
        ("normal_tumour_mapq", (self.tt[5][2] / self.tt[5][0]) / ((self.nt[5][2] / self.nt[5][0]) + self.ep)),
        ("normal_tumour_quality", (self.tt[5][1] / self.tt[5][0]) / ((self.nt[5][1] / self.nt[5][0]) + self.ep)),
        ("normal_tumour_ref_depth", ((self.tt[self.b][0] / self.tt[5][0]) + self.ep) / ((self.nt[self.b][0] / self.nt[5][0]) + self.ep)),
        ("normal_tumour_ref_direction", (self.tt[self.b][4] / (self.tt[self.b][0] + self.ep)) / ((self.nt[self.b][4] / (self.nt[self.b][0] + self.ep)) + self.ep)),
        ("normal_variant_allele_frequency", self.nt_nonref_count / self.nt[5][0]),
        ("normal_variant_depth_ratio", ((self.nt[5][0] - self.nt[self.b][0]) / self.nt[5][0])),
        ("normal_variant_direction_ratio", self.nt[self.nonref_index][4] / (self.nt_nonref_count + self.ep)),
        ("normal_variant_distance", self.nt[self.nonref_index][3] / (self.nt_nonref_count + self.ep)),
        ("normal_variant_mapq_mean", self.nt[self.nonref_index][2] / (self.nt_nonref_count + self.ep)),
        ("normal_variant_quality_ratio", ((self.nt[5][1] - self.nt[self.b][1]) / (self.nt[5][0] - self.nt[self.b][0] + self.ep))),
        ("region_entropy", self.rt[4]),
        ("region_gc_content", self.rt[3]),
        ("tumour_direction_ratio", self.tt[5][4] / self.tt[5][0]),
        ("tumour_distance_ratio", self.tt[5][3] / self.tt[5][0]),
        ("tumour_entropy", self.tt[10]),
        ("tumour_indels", (self.tt[9]+self.tt[11]) / self.tt[5][0]),
        ("tumour_mapq_ratio", self.tt[5][2] / self.tt[5][0]),
        ("tumour_quality_ratio", self.tt[5][1] / self.tt[5][0]),
        ("tumour_ref_depth", self.tt[self.b][0] / (self.tt[self.b][0] + self.ep)),
        ("tumour_ref_direction", self.tt[self.b][4] / (self.tt[self.b][0] + self.ep)),
        ("tumour_ref_direction_total", self.tt[self.b][4] / self.tt[5][0]),
        ("tumour_ref_quality", self.tt[self.b][1] / (self.tt[self.b][0] + self.ep)),#self.tt[5][0]),
        ("tumour_variant_allele_frequency", self.tt_nonref_count / self.tt[5][0]),
        ("tumour_variant_depth_ratio", ((self.tt[5][0] - self.tt[self.b][0]) / self.tt[5][0])),        
        ("tumour_variant_direction_ratio", self.tt[self.nonref_index][4] / (self.tt_nonref_count + self.ep)),
        ("tumour_variant_distance", self.tt[self.nonref_index][3] / (self.tt_nonref_count + self.ep)),
        ("tumour_variant_mapq_mean", self.tt[self.nonref_index][2] / (self.tt_nonref_count + self.ep)),
        ("tumour_variant_quality_ratio", ((self.tt[5][1] - self.tt[self.b][1]) / (self.tt[5][0] - self.tt[self.b][0] + self.ep))),

        ("normal_direction_ratio_mean", self.nt[self.nonref_index][4]/(self.nt[self.nonref_index][0]+self.ep) / (self.nt[5][0]+self.ep)),
        ("normal_mapq_ratio_mean", self.nt[self.nonref_index][2]/(self.nt[self.nonref_index][0]+self.ep) / (self.nt[5][0]+self.ep) ),
        ("normal_quality_ratio_mean", self.nt[self.nonref_index][1]/(self.nt[self.nonref_index][0]+self.ep) / (self.nt[5][0]+self.ep) ),
        ("normal_distance_ratio_mean", self.nt[self.nonref_index][3]/(self.nt[self.nonref_index][0]+self.ep) / (self.nt[5][0]+self.ep)),
        #("normal_nonref_depth_mean", self.nt[self.nonref_index][0]/(self.nt[self.nonref_index][0]+self.ep) / (self.nt[5][0]+self.ep)),
        
        ("tumour_direction_ratio_mean", self.tt[self.nonref_index][4]/(self.tt[self.nonref_index][0]+self.ep) / (self.tt[5][0]+self.ep)),
        ("tumour_mapq_ratio_mean", self.tt[self.nonref_index][2]/(self.tt[self.nonref_index][0]+self.ep) / (self.tt[5][0]+self.ep)),
        ("tumour_quality_ratio_mean", self.tt[self.nonref_index][1]/(self.tt[self.nonref_index][0]+self.ep) / (self.tt[5][0]+self.ep)),
        ("tumour_distance_ratio_mean", self.tt[self.nonref_index][3]/(self.tt[self.nonref_index][0]+self.ep) / (self.tt[5][0]+self.ep)),
        #("tumour_nonref_depth_mean", self.tt[self.nonref_index][0]/(self.tt[self.nonref_index][0]+self.ep) / (self.tt[5][0]+self.ep)),
                
        #("normal_tumour_depth_mean", ((self.nt[self.nonref_index][0]/(self.nt[self.nonref_index][0]+self.ep))/self.nt[5][0]) / ((self.tt[self.nonref_index][0]/(self.tt[self.nonref_index][0]+self.ep))/self.tt[5][0]) ),
        ("normal_tumour_direction_mean", ((self.nt[self.nonref_index][4]/(self.nt[self.nonref_index][0]+self.ep))/self.nt[5][0]) / ((self.tt[self.nonref_index][4]/(self.tt[self.nonref_index][0]+self.ep))/self.tt[5][0] +self.ep) ),
        ("normal_tumour_mapq_mean",((self.nt[self.nonref_index][2]/(self.nt[self.nonref_index][0]+self.ep))/self.nt[5][0]) / ((self.tt[self.nonref_index][2]/(self.tt[self.nonref_index][0]+self.ep))/self.tt[5][0]+self.ep) ),
        ("normal_tumour_quality_mean", ((self.nt[self.nonref_index][1]/(self.nt[self.nonref_index][0]+self.ep))/self.nt[5][0]) / ((self.tt[self.nonref_index][1]/(self.tt[self.nonref_index][0]+self.ep))/self.tt[5][0]+self.ep) ),
        ("normal_tumour_distance_mean", ((self.nt[self.nonref_index][3]/(self.nt[self.nonref_index][0]+self.ep))/self.nt[5][0]) / ((self.tt[self.nonref_index][3]/(self.tt[self.nonref_index][0]+self.ep))/self.tt[5][0]+self.ep) ),
       
        ## jointsnvmix
        ("tumour_normal_jointsnvmix_somatic", self.p_somatic),
        ("tumour_normal_jointsnvmix_germline", self.p_germline),
        ("tumour_normal_jointsnvmix_wildtype", self.p_wildtype),
        ("tumour_normal_jointsnvmix_loh", self.p_loh),
        ("tumour_normal_jointsnvmix_error", self.p_error),
        
        ("normal_reads_indel",self.nt[-2]),
        ("tumour_reads_indel",self.tt[-2]),
        
        ("tumour_background_average_variant_mean",self.tt_bg_avg_var_freq),
        ("tumour_p_value",self.tum_p_val),
        
        ("normal_background_average_variant_mean",self.nt_bg_avg_var_freq),
        ("normal_p_value",self.norm_p_val),
        
       # ("normal_tumour_ref_quality", ((self.tt[self.b][1] / (self.tt[self.b][0]+ self.ep)) + self.ep) / ((self.nt[self.b][1] / (self.nt[self.b][0]+ self.ep)) + self.ep)),
       # ("normal_tumour_ref_distance", ((self.tt[self.b][3] / (self.tt[self.b][0]+ self.ep)) + self.ep) / ((self.nt[self.b][3] / (self.nt[self.b][0]+ self.ep)) + self.ep)),        
       #("normal_minor_allele", lambda t, n, r: (r[0] != t[7] and n[t[6] + 1][0] > 0) or (r[0] != t[6] and n[t[7] + 1][0] > 0)),
       # ("normal_tumour_ref_mapq", ((self.tt[self.b][2] / (self.tt[self.b][0]+ self.ep)) + self.ep) / ((self.nt[self.b][2] / (self.nt[self.b][0]+ self.ep)) + self.ep)),
       # ("normal_tumour_distance", (self.tt[5][3] / self.tt[5][0]) / ((self.nt[5][3] / self.nt[5][0]) + self.ep)),        
       # ("tumour_variant_quality", (self.tt[5][1] - self.tt[self.b][1])/self.tt[5][0]),
       # ("normal_variant_quality", (self.nt[5][1] - self.nt[self.b][1])/self.nt[5][0]),       
       # ("normal_tumour_variant_depth_ratio", ((self.tt[5][0] - self.tt[self.b][0]) / self.tt[5][0]) / (((self.nt[5][0] - self.nt[self.b][0]) / self.nt[5][0]) + self.ep)),
       # ("normal_tumour_variant_mapq_ratio", ((self.tt[5][2] - self.tt[self.b][2]) / (self.tt[5][0] - self.tt[self.b][0] + self.ep)) / ((self.nt[5][2] - self.nt[self.b][2]) / (self.nt[5][0] - self.nt[self.b][0] + self.ep) + self.ep)),
       # ("normal_tumour_variant_quality_ratio", ((self.tt[5][1] - self.tt[self.b][1]) / (self.tt[5][0] - self.tt[self.b][0] + self.ep)) / ((self.nt[5][1] - self.nt[self.b][1]) / (self.nt[5][0] - self.nt[self.b][0] + self.ep) + self.ep)),
       # ("normal_tumour_variant_distance_ratio", ((self.tt[5][3] - self.tt[self.b][3]) / (self.tt[5][0] - self.tt[self.b][0] + self.ep)) / ((self.nt[5][3] - self.nt[self.b][3]) / (self.nt[5][0] - self.nt[self.b][0] + self.ep) + self.ep)),
        )
        
        self.coverage_features = (
        ("tumour_contamination", self.cd[2] / 100),
        ("whole_genome", self.cd[3])
        )

    def __get_tumour_nonref_index(self,ref,tt):
        nonrefbases = [x for x in range(4) if x != ref - 1]
        max_nonrefbase_depth = 0
        nonref_index = nonrefbases[1]  
        
        for nb in nonrefbases:
            index = nb + 1
            tumour_nonrefbase_depth = tt[index][0]
            if tumour_nonrefbase_depth != 0 and tumour_nonrefbase_depth > max_nonrefbase_depth:
                max_nonrefbase_depth = tt[index][0]
                nonref_index = index
                
        return nonref_index
    
    def __get_tumour_nonref_count(self):
        self.tt_nonref_count = self.tt[self.nonref_index][0]

    def __get_normal_nonref_count(self):
        self.nt_nonref_count = self.nt[self.nonref_index][0]

    ##TODO: eventually this function should be removed    
    def __isvalid(self, x):
        if numpy.isnan(x) or numpy.isinf(x):
            
            ##TODO: remove this line
            print "NaN"
            return False
        return True
    
    def __xentropy(self):
        tumour_base_qualities = (self.tt[1][1], self.tt[2][1], self.tt[3][1], self.tt[4][1], self.tt[5][1])
        normal_base_qualities = (self.nt[1][1], self.nt[2][1], self.nt[3][1], self.nt[4][1], self.nt[5][1])
        total_tbq = tumour_base_qualities[4]
        total_nbq = normal_base_qualities[4]
        ent = 0 # entropy
        
        if total_tbq == 0 or total_nbq == 0:
            return ent
            
        for i in xrange(4):
            tumour_base_probability = tumour_base_qualities[i] / total_tbq
            normal_base_probability = normal_base_qualities[i] / total_nbq            
            if tumour_base_probability != 0:
                if normal_base_probability == 0:
                    ent -= -7 * tumour_base_probability
                else:
                    ent -= log(normal_base_probability) * tumour_base_probability
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
        feature_names.append('xentropy')
        
        return feature_names

    def __jointsnvmix(self):
        ## a_priori
        aa_aa = 1e3
        aa_ab = 1e1
        aa_bb = 1e1
        ab_aa = 1e1
        ab_ab = 1e2
        ab_bb = 1e1
        bb_aa = 1
        bb_ab = 1
        bb_bb = 1e2
        a_priori_list = numpy.array([aa_aa, aa_ab, aa_bb, ab_aa, ab_ab, ab_bb, bb_aa, bb_ab, bb_bb])
        a_priori_mat  = a_priori_list.reshape(3,3)
        
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
        t_binom = [binom.pmf(self.tt_nonref_count, self.tt[5][0], p) for p in t_prob]
        n_binom = [binom.pmf(self.nt_nonref_count, self.nt[5][0], p) for p in n_prob]
        
        ## transpose(t_binom) x n_binome = 3x3 binom_mat matrix
        binom_mat = [[i*j for i in n_binom] for j in t_binom]
        binom_mat = numpy.array(binom_mat)
        
        ## element by element product of a_priori matrix to the m_binom matrix
        p_mat = numpy.multiply(a_priori_mat, binom_mat)
        
        ## normalize 
        weight = sum(sum(p_mat))
        p_mat = p_mat / (weight + self.ep)
    
        ## p_SOMATIC = p_aa,ab + p_aa,bb
        self.p_somatic = p_mat[1][0] + p_mat[2][0]
        
        ## p_WILDTYPE = p_aa,aa
        self.p_wildtype = p_mat[0][0]
        
        ## p_GERMLINE = p_ab,ab + p_bb,bb
        self.p_germline = p_mat[1][1] + p_mat[2][2]
        
        ## p_LOH = p_ab,aa + p_ab,bb
        self.p_loh = p_mat[0][1] + p_mat[2][1]
        
        ## p_ERROR = p_bb,aa + p_bb,ab
        self.p_error = p_mat[0][2] + p_mat[1][2]
        
    def __background_var_freq_mean(self):
        tt_bg_avg_var_freq = []
        nt_bg_avg_var_freq = []
        
        for tt in self.tt_bg:
            rt = [val for val in self.rt_bg if val[0] == tt[0]]
            if rt == []:
                raise Exception('Error: reference tuple for position '+ str(tt[0])+ ' couldn\'t be extracted')
            
            ref = rt[0][1][0] + 1
            nonref_index = self.__get_tumour_nonref_index(ref,tt)
            tt_bg_avg_var_freq.append(tt[nonref_index][0]/tt[5][0])
            
        for nt in self.nt_bg:
            tt = [val for val in self.tt_bg if val[0] == nt[0]]
            if tt == []:
                logging.warning('Skipping the position:'+ str(nt[0])+' as tumour tuple not present')
                continue
            
            rt = [val for val in self.rt_bg if val[0] == tt[0][0]]
            if rt == []:
                raise Exception('Error: reference tuple for position '+ str(tt[0])+ ' couldn\'t be extracted')
            

            
            ref = rt[0][1][0] +1
            nonref_index = self.__get_tumour_nonref_index(ref,tt[0])
            nt_bg_avg_var_freq.append(nt[nonref_index][0]/nt[5][0])
            
        try:
            self.tt_bg_avg_var_freq = sum(tt_bg_avg_var_freq)/len(tt_bg_avg_var_freq)
            self.nt_bg_avg_var_freq = sum(nt_bg_avg_var_freq)/len(nt_bg_avg_var_freq)
        except Exception,e:
            print e
            
        self.tum_p_val = binom.sf(self.tt_nonref_count,self.tt[5][0],self.tt_bg_avg_var_freq)
        self.norm_p_val = binom.sf(self.nt_nonref_count,self.nt[5][0],self.nt_bg_avg_var_freq)
            
        
        
        