# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
import logging
import classifyui

mutationSeq_version="4.0.0"
args = classifyui.args 

if args.verbose:
    level = logging.DEBUG    

else:
    level = logging.WARNING
    
logging.basicConfig(filename = args.log_file, 
                    format   = '%(asctime)s %(message)s', 
                    #datefmt = '%m/%d/%Y %I:%M:%S %p', 
                    level = level)

logging.warning("<<< mutationSeq_" + mutationSeq_version + " started >>>")
logging.info("importing required modules")
import bamutils

logging.info(args)
#==============================================================================
# main body
#==============================================================================
bam_helper = bamutils.BamHelper(args)
target_positions = bam_helper.get_positions()

logging.info("generating tuple iterator")
tuples = bam_helper.bam.get_tuples(target_positions)

logging.info("generating features iterator")
features = bam_helper.get_features(tuples)

if args.export_features is not None:
    bam_helper.export_features(features)
    
probabilities = bam_helper.predict(features)
bam_helper.print_results(probabilities)

logging.info("successfully completed.\n")
