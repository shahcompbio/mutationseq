# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:13:27 2013

@author: jtaghiyar
"""

import bamutils
from collections import namedtuple

class AmpDepth:
    def __init__(self, bam_file):
        self.bam = bamutils.Bam(tumour = bam_file)
        
    def __getTuples(self, chromosome, start, stop):
        Pos = namedtuple("Position", "chromosome, start, stop")
        target_position = Pos._make([chromosome, start, stop])
        print target_position
        return self.bam.getTumourTuples([target_position]) #returns a generator
        
    def getAverageDepth(self, chromosome, start, stop):
        count = 0
        total_depth = 0
        for t in self.__getTuples(chromosome, start, stop):
            if t is not None:
                total_depth += t[5][0]
                count += 1
        
        if count != 0:
            return total_depth/count
        else:
            return 0
        
#def main():
#    print "reading bam ..."    
#    bam_file = "/share/lustre/projects/breast_xeno_evolution/samples/SA530X2/illumina_wgss/A30648/bwa_aligned/SA530X2_A30648_5_lanes_hg19_dupsFlagged.bam"
#    ad = AmpDepth(bam_file)
#    
#    print "getting average depth ..."
#    for n in ad.bam.getTumourRefNames():
#        print ad.getAverageDepth(n, 1000000, 1000001)
#    
#if __name__ == '__main__':
#    main()