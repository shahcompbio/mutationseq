# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
import logging
import bamutils
import pybamapi
import classifyui

mutationSeq_version="4.0.0"
args = classifyui.args                    
samples = {}
for s in args.samples:
    samples[s.split(':')[0]] = s.split(':')[1]
    
if args.verbose:
    level = logging.DEBUG
else:
    level = logging.INFO
    
logging.basicConfig(filename = "mutationSeq_run.log", 
                    format   = '%(asctime)s %(message)s', 
                    #datefmt = '%m/%d/%Y %I:%M:%S %p', 
                    level    = level)

info_str = " mutationSeq_" + mutationSeq_version + " started\n" + \
            " tumour:"     + samples.get("tumour") + \
            " normal:"     + samples.get("normal") + \
            " reference:"  + samples.get("reference") + \
            " --interval " + str(args.interval)
            
logging.info(info_str)

#==============================================================================
# main body
#==============================================================================
bam = pybamapi.BamApi(tumour=samples.get("tumour"), normal=samples.get("normal"), reference=samples.get("reference"))
bam_helper = bamutils.BamHelper(bam, args)

logging.info("getting positions ...")
target_positions = bam_helper.get_positions()

logging.info("getting tuples ...")
tuples = bam.get_pair_tuples(target_positions)

logging.info("getting features ...")
features = bam_helper.get_features(tuples)

if args.export is not None:
    logging.info("exporting features to:" + str(args.export))
    bam_helper.export_features(features)
    
logging.info("fitting model and predict probabilities ...")
if len(features) != 0:
    probabilities = bam_helper.predict(features)
else:
    probabilities = None
        
logging.info("printting results to:" + str(args.out))
bam_helper.print_results(probabilities)

logging.info("successfully completed.\n")