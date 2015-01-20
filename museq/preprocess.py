'''
Created on May 20, 2014

@author: dgrewal
'''
import logging
import preprocessui

mutationSeq_version="4.3.4"
args = preprocessui.args 

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
import preprocessutils

logging.info(args)
#==============================================================================
# main body
#==============================================================================
logging.info("initializing a Classifier")
classifier = preprocessutils.PreProcess(args)

logging.info("getting positions")
target_positions = classifier.get_positions()

logging.info("generating tuple iterator")
tuples = classifier.bam.get_tuples(target_positions)

logging.info("generating features iterator")
features = classifier.get_features(tuples)
    
probabilities = classifier.predict(features)
classifier.print_results(probabilities)

logging.warning("successfully completed.\n")
