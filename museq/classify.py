# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
import logging
import classifyui

mutationSeq_version="4.3.0"
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
logging.info("initializing a Classifier")
classifier = bamutils.Classifier(args)

logging.info("getting positions")
target_positions = classifier.get_positions()

logging.info("generating tuple iterator")
tuples = classifier.bam.get_tuples(target_positions)

logging.info("generating features iterator")
features = classifier.get_features(tuples)

if args.export_features is not None:
    logging.info("exporting features")
    classifier.export_features(features)

    
probabilities = classifier.predict(features)
classifier.print_results(probabilities)

logging.warning("successfully completed.\n")
