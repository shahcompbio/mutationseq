import argparse
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
parser.add_argument("--version", 
                    action="version", 
                    version=mutationSeq_version)
 
## mandatory options                   
mandatory_options = parser.add_argument_group("required arguments")
mandatory_options.add_argument("-o", "--out", 
                               default=None, 
                               #required=True, 
                               help="specify the path/to/out.vcf to save output to a file")                               
mandatory_options.add_argument("-c", "--config", 
                               default=None, 
                               #required=True,
                               help="specify the path/to/metadata.config file used to add \
                               meta information to the output file")
                    
## mutually exclusive options
exgroup = parser.add_mutually_exclusive_group()
exgroup.add_argument("-f", "--positions_file", 
                     default=None, 
                     help="input a file containing a list of positions each of which in a \
                     separate line, e.g. chr1:12345\nchr2:23456")                     
exgroup.add_argument("-i", "--interval",
                     default=None,
                     help="classify given chromosome[:start-end] range")

args = parser.parse_args()

#==============================================================================
# check the input 
#==============================================================================
import sys
from warnings import warn

if len(args.samples) < 3:
    print >> sys.stderr, "bad input, usage: 'classify.py normal:<normal.bam> \
    tumour:<tumour.bam> reference:<ref.fasta> model:<model.npz> [--options]'"
    sys.exit(1)

#if not args.features_only:
if args.out is None:
    warn("--out is not specified, standard output is used to write the results")
    out = sys.stdout
else:
    out = open(args.out, 'w')

samples = {}
for s in args.samples:
    samples[s.split(':')[0]] = s.split(':')[1]

#==============================================================================
# import required modules here to save time when only checking version or help
#==============================================================================
from datetime import datetime
import BamClass_newPybam
import bamutils

print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " mutationSeq_" + mutationSeq_version + " started"
#==============================================================================
# beginning of the main body
#==============================================================================
bam = BamClass_newPybam.Bam(tumour=samples["tumour"], normal=samples["normal"], reference=samples["reference"])
bam_helper = bamutils.BamUtils(bam, args)

## get target positions "chromosome:start-stop"
target_positions = bam_helper.getPositions()

## get an iterator over tuples of the target positions of tumour/normal bam files
tumour_tuples = bam_helper.getTumourTuple(target_positions)
normal_tuples = bam_helper.getNormalTuple(target_positions)

## calculate features for the candidate positions
features = bam_helper.getFeatures(tumour_tuples, normal_tuples)

## predict the probabilities
if len(features) != 0:
    probabilities = bam_helper.predict(features)
else:
    probabilities = None
        
## print the output string to out
bam_helper.printResults(out, probabilities)

#==============================================================================
# end of the main body
#==============================================================================
try:
#    expfile.close()
    out.close()
except:
    pass
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " done."