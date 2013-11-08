# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
import sys
import bamutils
import argparse
import logging
import resource

mutationSeq_version="4.0.0"

#==============================================================================
# make a UI 
#==============================================================================
parser = argparse.ArgumentParser(prog='mutationSeq', 
                                 description = '''mutationSeq: a feature-based classifier 
                                 for somatic mutation detection in
                                 tumour-normal paired sequencing data''')
## positional arguments                    
parser.add_argument("samples", 
                    nargs='*', 
                    help='''A list of colon delimited sample names; normal:normal.bam
                    tumour:tumour.bam model:model.npz reference:reference.fasta''')

parser.add_argument("-l", "--log_file",
                    default="mutationSeq_run.log",
                    help='''specify name or path of the log file''')                    

parser.add_argument("--coverage",
                    default=4,
                    type=int,
                    help='''specify the depth of the coverage to be considered''')
                    
parser.add_argument("--normal_variant",
                    default=25,
                    type=int,
                    help='''specify the max variant percentage in the normal bam file''')

parser.add_argument("--tumour_variant",
                    default=2,
                    type=int,
                    help='''specify the min number of variants in the tumour bam file''')

parser.add_argument("-a", "--all", 
                    default=None, choices=["no", "yes"], 
                    help= '''force to print out even if the position(s) does not satisfy 
                    the initial criteria for Somatic call''')
                    
parser.add_argument("-e" , "--export", 
                    default=None, 
                    help='''save exported feature vector to the specified path''')
                    
parser.add_argument("-u", "--features_only", 
                    default=False, action="store_true", 
                    help='''if true, only extracted features are exported''')
                    
parser.add_argument("-d", "--deep", 
                    default=False, action="store_true", 
                    help='''for deepseq data''')
                    
parser.add_argument("-n", "--normalized", 
                    default=False, action="store_true",
                    help='''If you want to test with normalized features 
                    (the number of features are also different from non-deep)''')
                    
parser.add_argument("-p", "--purity", 
                    default=70, 
                    type=int,
                    help='''pass sample purity to features''')
                    
parser.add_argument("-v", "--verbose", 
                    action="store_true", default=False,
                    help='''verbose''')
                    
parser.add_argument("--version", 
                    action="version", version=mutationSeq_version)
                    
parser.add_argument("-t", "--threshold", 
                    default=0.5, type=float,
                    help='''set threshold for positive call''') 

## mandatory options                   
mandatory_options = parser.add_argument_group("required arguments")
mandatory_options.add_argument("-o", "--out", 
                               default=None, 
                               #required=True, 
                               help='''specify the path/to/out.vcf to save output to a file''')  
                               
mandatory_options.add_argument("-c", "--config", 
                               default=None, 
                               #required=True,
                               help='''specify the path/to/metadata.config file used to add 
                               meta information to the output file''')
                    
## mutually exclusive options
exgroup = parser.add_mutually_exclusive_group()
exgroup.add_argument("-f", "--positions_file", 
                     default=None, 
                     help='''input a file containing a list of positions each of which in
                     a separate line, e.g. chr1:12345\nchr2:23456''')                  
                     
exgroup.add_argument("-i", "--interval",
                     default=None,
                     help='''specify an interval "chr[:start-stop]"''')

args = parser.parse_args()

#==============================================================================
# check the input 
#==============================================================================
logging.basicConfig(filename=args.log_file, 
                    format='%(asctime)s %(message)s', 
                    #datefmt='%m/%d/%Y %I:%M:%S %p', 
                    level=logging.DEBUG)
                    
samples = {}
for s in args.samples:
    samples[s.split(':')[0]] = s.split(':')[1]

coverage = args.coverage

logging.info(args)

if len(args.samples) < 3:
    logging.error("""bad input, usage: 'classify.py normal:<normal.bam> 
    tumour:<tumour.bam> reference:<ref.fasta> model:<model.npz> [--options]'""")
    sys.exit(1)

#if not args.features_only:
if args.out is None:
    logging.warning("--out is not specified, standard output is used to write the results")
    out = sys.stdout
else:
    out = open(args.out, 'w')

#==============================================================================
# beginning of the main body
#==============================================================================
bam = bamutils.Bam(tumour=samples.get("tumour"), normal=samples.get("normal"),
                   reference=samples.get("reference"), coverage=coverage, rmdups=None)

bam_helper = bamutils.BamUtils(bam, args)

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
    logging.info("exporting features ...")
    with open(args.export, 'w') as ff:
        for i in features:
            print >> ff, i
    
logging.info("fitting model and predict probabilities ...")
if len(features) != 0:
    probabilities = bam_helper.predict(features)
else:
    probabilities = None
        
info_str = "printting results to:" + str(args.out) 
logging.info(info_str)
bam_helper.print_results(out, probabilities)

#==============================================================================
# household cleaning
#==============================================================================
try:
    out.close()
except:
    pass
logging.info("successfully completed.\n")
