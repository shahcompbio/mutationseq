# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
import logging
import bamutils
import classifyui
import resource

mutationSeq_version="4.0.0"
args = classifyui.args                    

if args.verbose:
    level = logging.DEBUG
else:
    level = logging.INFO
    
logging.basicConfig(filename = args.log_file, 
                    format   = '%(asctime)s %(message)s', 
                    #datefmt = '%m/%d/%Y %I:%M:%S %p', 
                    level    = level)

logging.info("mutationSeq_" + mutationSeq_version + " started >>>>>>")
logging.info(args)

#==============================================================================
# main body
#==============================================================================
bam_helper = bamutils.BamHelper(args)

logging.info("getting positions ...")
target_positions = bam_helper.get_positions()

logging.info("getting tuples ...")
tuples = bam_helper.bam.get_pair_tuples(target_positions)

logging.info("getting features ...")
features = bam_helper.get_features(tuples)

##TODO: remove this line
logging.info("total mem usage: " + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

if args.export is not None:
    logging.info("exporting features to:" + str(args.export))
    bam_helper.export_features(features)
    
logging.info("fitting model and predict probabilities ...")
probabilities = bam_helper.predict(features)

        
logging.info("printting results to:" + str(args.out))
bam_helper.print_results(probabilities)

logging.info("successfully completed.\n")