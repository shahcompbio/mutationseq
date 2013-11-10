# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
import logging
import bamutils
import pybamapi
import classifyui
import resource

mutationSeq_version="4.0.0"
args = classifyui.args                    

## parse out tumour/normal/reference/model from the positional arguments
samples = {}
for s in args.samples:
    samples[s.split(':')[0]] = s.split(':')[1]

## min coverage of a position to get tuples for
coverage = args.coverage
   
if args.verbose:
    level = logging.DEBUG
else:
    level = logging.INFO
    
logging.basicConfig(filename = args.log_file, 
                    format   = '%(asctime)s %(message)s', 
                    #datefmt = '%m/%d/%Y %I:%M:%S %p', 
                    level    = level)

logging.info(args)

#==============================================================================
# main body
#==============================================================================
bam = pybamapi.BamApi(tumour=samples.get("tumour"), normal=samples.get("normal"), 
                      reference=samples.get("reference"), coverage=coverage, rmdups=None)
                      
bam_helper = bamutils.BamHelper(bam, args)

logging.info("getting positions ...")
target_positions = bam_helper.get_positions()

logging.info("getting tuples ...")
#tumour_tuples = bam.get_tumour_tuples(target_positions)
#normal_tuples = bam.get_normal_tuples(target_positions)

tuples = bam.get_pair_tuples(target_positions)

logging.info("getting features ...")
#features = bam_helper.get_features(tumour_tuples, normal_tuples)

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