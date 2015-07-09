'''
Created on May 20, 2014

@author: dgrewal
'''
import logging
import classifyui

mutationSeq_version="4.3.6"
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
import preprocessutils

logging.info(args)
#==============================================================================
# main body
#==============================================================================
logging.info("initializing a Classifier")
classifier = preprocessutils.PreProcess(args)

logging.info("getting positions")
classifier.get_positions()

logging.info("generating features iterator")
features = classifier.get_features()
    
probabilities = classifier.predict(features)
classifier.print_results(probabilities)

logging.warning("successfully completed.\n")
