# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:10:07 2013
 
"""
import newpybam as np # new pybam

class Bam:
    def __init__(self, **kwargs):       
        self.t_bam = kwargs.get("tumour_bam")
        self.n_bam = kwargs.get("normal_bam")
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

    def getReferenceTuple(self, chromosome, position):
        """get reference tuple of chromosome:position"""
        
        temp_tuple = self.fasta.get(chromosome, position)
        self.refBase = temp_tuple[0]
        return temp_tuple
            
    def getReferenceBase(self, chromosome, position):
        """get reference nucleotide of chromosom:position"""

        if  self.refBase is None:
            self.refBase = self.getReferenceTuple(chromosome, position)[0]

        # index ACGTN
        if  self.refBase == "A":
            return 0
        elif self.refBase == "C":
            return 1
        elif self.refBase == "G":
            return 2
        elif self.refBase == "T":
            return 3
        else:
            return 4
    
    def getNormalTuple(self, chromosome, position):
        self.n_pileup.jump(chromosome, position)
        return self.n_pileup.next()
            
    def getTumourTuple(self, chromosome, position):
        self.t_pileup.jump(chromosome, position)
        return self.t_pileup.next()

    def getNormalChromosomeIds(self):
        return self.n_pileup.refnames
    
    def getTumourChromosomeIds(self):
        return self.t_pileup.refnames
        