# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 11:14:37 2013

@author: jtaghiyar
"""
import newpybam as np # new pybam

class BamApi:
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        