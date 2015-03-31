# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 10:43:37 2013

@author: jtaghiyar
"""
import logging
import trainui

mutationSeq_version="4.3.5"
args = trainui.args 

if args.verbose:
    level = logging.DEBUG    

else:
    level = logging.WARNING
    
logging.basicConfig(filename=args.log_file, format='%(asctime)s %(message)s', 
                    #datefmt='%m/%d/%Y %I:%M:%S %p', 
                    level=level)

logging.warning("<<< mutationSeq_training_" + mutationSeq_version + " started >>>")
logging.info("importing required modules")
import bamutils

logging.info(args)
#==============================================================================
# main body
#==============================================================================
logging.info("initializing a trainer")
model = bamutils.Trainer(args)

## check if should generate a new model or load a given model
if args.model is None:
    logging.info("generating a model")
    model.generate()
    
    logging.info("fitting a model")
    model.fit()
    
    logging.info("saving the trained model")
    model.save()

else:
    logging.info("loading the model")
    model.load()

logging.info("writing the sorted list of features' importance")
model.print_feature_importance()

## validate the model and save the results of validation into args.out/
if args.validate:
    logging.info("validating the new model")
    model.validate()
    
else:
    logging.info("cross validating the new model")
    model.cross_validate()

logging.info("generating boxplots for the new model")
model.generate_boxplot() 
    
logging.warning("model training successfully completed.\n")
