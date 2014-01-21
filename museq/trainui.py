# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 10:41:52 2013

@author: jtaghiyar
"""
import argparse
import sys

mutationSeq_version="4.1.0"

#==============================================================================
# make a UI 
#==============================================================================
parser = argparse.ArgumentParser(description='''train a model''')

## positional arguments
parser.add_argument("infiles", 
                    metavar='FILE', nargs='*', default=[sys.stdin], type=argparse.FileType('r'), 
                    help= '''A list of space delimited samples: chromosome position label, 
                    the file has a header of normal, tumour and reference file names''')

## optional arguments 
parser.add_argument("-d", "--deep", 
                    default=False, action="store_true", 
                    help="If you want to test on deep data you need to change contamination rate")
                    
parser.add_argument("--labels", 
                    default="SOMATIC",
                    help="Labels in the training file list to be considered as positive labels. Valid labels are SOMATIC, WILDTYPE, GERMLINE, HET, HOM, CLEAN")

parser.add_argument("-l", "--log_file",
                    default="mutationSeq_training.log",
                    help ='''specify name or path of the log file''')

parser.add_argument("-m", "--model",
                    default=None,
                    help="specify an existing model. Usually used for validation")

parser.add_argument("-n", "--normalized", 
                    default=False, action="store_true", 
                    help="If you want to train with normalized features")
                    
parser.add_argument("-o", "--out", 
                    default=None, 
                    help="save output to file")

parser.add_argument("-s", "--single", 
                    default=False, action="store_true", 
                    help="single sample training")

parser.add_argument("--validate", 
                    metavar='FILE', nargs='*', type=argparse.FileType('r'), default=None, 
                    help="To validate the same format file with known label")
                    
parser.add_argument("--version", 
                    action="version", version=mutationSeq_version)
                    
parser.add_argument("-v", "--verbose", 
                    action="store_true", default=False,
                    help="verbose")

args = parser.parse_args()
