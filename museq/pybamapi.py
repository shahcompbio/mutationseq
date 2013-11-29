# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 11:14:37 2013

@author: jtaghiyar
"""
import newpybam as np # new pybam
import cStringIO

class Bam(object):
    def __init__(self, **kwargs):       
        self.bam = kwargs.get("bam")
        self.ref = kwargs.get("reference") 
        rmdups   = kwargs.get("rmdups") 
        coverage = kwargs.get("coverage")
        self.base = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
        
        ## check if duplicates need to be removed
        if rmdups is None:
            rmdups = True
            
        ## set the default for the coverage
        if coverage is None:
            coverage = 4
            
        ## make a pileup for bam file
        if self.bam is not None:
            self.pileup = np.pileup(coverage, rmdups)
            self.pileup.open(self.bam)
    
        if  self.ref is not None:
            self.__load_reference()
        
    def __load_reference(self):
        ## make a fasta object to laod reference
        self.fasta = np.fasta() 
        self.fasta.open(self.ref) 

    def get_reference_base(self, chromosome_id, position, index=False):
        if position < 1:
            return 'N'
            
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
        refbase = self.base.get(refbase)
        if refbase is None:
            refbase = 4

        ## replace the base with its index (required for feature extraction)
        temp_tuple = (refbase, temp_tuple[1], temp_tuple[2], temp_tuple[3], temp_tuple[4])
        return temp_tuple
    
    def get_refnames(self):
        """ get the list of chromosome names and their corresponding IDs as a dictionary """
        
        return dict(self.pileup.refnames)
    
    def get_chromosome_id(self, chromosome_name):
        rn = self.get_refnames()
        return rn.get(chromosome_name)
    
    def get_chromosome_name(self, chromosome_id):
        rn = self.get_refnames()
        for k in rn.keys():
            if rn[k] == chromosome_id:
                return k
        return None
        
    def get_samheader(self):
        return self.pileup.samheader
        
    def get_chromosome_lengths(self):
        """ parse the sam file header to get the chromosome lengths """
        
        chromosome_lengths = {}
        sam_header  = self.pileup.samheader
        
        ## make a file_like object
        sam_file = cStringIO.StringIO(sam_header)
        
        for l in sam_file:
            l = l.strip().split()
            if l[0] == "@SQ" :
                SN = l[1]
                chrom_name = SN.split(':')[1]
                LN = l[2]
                chrom_len  = LN.split(':')[1]
                chromosome_lengths[chrom_name] = int(chrom_len)

        sam_file.close()
        return chromosome_lengths
        
    def get_tuple(self, chromosome, position):
        self.pileup.set_region(chromosome, position, position)
        return self.pileup.get_tuple()
    
    def get_tuples(self, target_positions):
        for tp in target_positions: 
            if tp[1] is None:
                self.pileup.set_region(tp[0])
            else:
                self.pileup.set_region(tp[0], tp[1], tp[2])

            while True:
                t = self.pileup.get_tuple() 
                if t is None:
                    break
                else:
                    yield t
            

class PairedBam(object):
    def __init__(self, tumour, normal, reference, rmdups, coverage):
        self.t_bam = Bam(bam=tumour, reference=reference, coverage=coverage, rmdups=rmdups)
        self.n_bam = Bam(bam=normal, reference=reference, coverage=coverage, rmdups=rmdups)
    
    def get_chromosome_name(self, chromosome_id):
        return self.t_bam.get_chromosome_name(chromosome_id)
     
    def get_trinucleotide_context(self, chromosome_id, position):
        return self.t_bam.get_trinucleotide_context(chromosome_id, position)
        
    def get_reference_tuple(self, chromosome_id, position):
        return self.t_bam.get_reference_tuple(chromosome_id, position)
    
    def get_reference_base(self, chromosome_id, position, index=False):
        return self.t_bam.get_reference_base(chromosome_id, position, index)
    
    def get_refnames(self):
        return self.t_bam.get_refnames()
        
    def get_tuples(self, target_positions):
        for tp in target_positions: 
            if tp[1] is None:
                self.t_bam.pileup.set_region(tp[0])
                self.n_bam.pileup.set_region(tp[0])
            else:
                self.t_bam.pileup.set_region(tp[0], tp[1], tp[2])
                self.n_bam.pileup.set_region(tp[0], tp[1], tp[2])
                
            while True:
                tt = self.t_bam.pileup.get_tuple() 
                nt = self.n_bam.pileup.get_tuple() 
            
                ## check for none tuples
                if not all((tt, nt)):
                    break
                
                ## check if the position is the same for both tumour/normal 
                while tt[0] != nt[0]:
                    if tt[0] < nt[0]:
                        tt = self.t_bam.pileup.get_tuple()
                        if tt is None:
                            break
                    else:
                        nt = self.n_bam.pileup.get_tuple()
                        if nt is None:
                            break
                        
                if not all((tt, nt)):
                    break
                yield (tt, nt)    

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        