# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:10:07 2013
 
"""
import newpybam as np # new pybam

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

    def getReferenceTuple(self, chromosome, position):
        """get reference tuple of chromosome:position"""
        
        temp_tuple = self.fasta.get(str(chromosome), int(position))
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
            
    def getReferenceBase(self, chromosome, position):
        """get reference nucleotide of chromosom:position"""

        if  self.refBase is None:
            self.refBase = self.getReferenceTuple(chromosome, int(position))[0]

        return self.refBase
    
    def getNormalTuple(self, chromosome, position):
        self.n_pileup.jump(chromosome, int(position))
        return self.n_pileup.next()
            
    def getTumourTuple(self, chromosome, position):
        self.t_pileup.jump(chromosome, int(position))
        return self.t_pileup.next()

    def getNormalChromosomeIds(self):
        return self.n_pileup.refnames
    
    def getTumourChromosomeIds(self):
        return self.t_pileup.refnames
        
