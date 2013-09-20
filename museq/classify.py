import argparse
mutationSeq_version="4.0.0"

#==============================================================================
# making a UI 
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
for sample in args.samples:
    samples[sample.split(':')[0]] = sample.split(':')[1]

#==============================================================================
# import required modules here to save time when only checking version or help
#==============================================================================
from datetime import datetime
from collections import deque, defaultdict, namedtuple
from string import Template
import BamClass_newPybam
import bamutils
#import features
#import Nfeatures
#import features_single
#import features_deep

#==============================================================================
# start of the main body
#==============================================================================
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " mutationSeq_" + mutationSeq_version + " started"

bam = BamClass_newPybam.Bam(tumour=samples["tumour"], normal=samples["normal"], reference=samples["reference"])
bam_helper = bamutils.BamUtils(bam, args)

## get target positions "chromosome:start-stop"
target_positions = bam_helper.getPositions()

## get nominated tuples from the tumour/normal bams for target positions 
tumour_tuples = bam_helper.getTumourTuple(target_positions)
normal_tuples = bam_helper.getNormalTuple(target_positions)

## calculate features for the candidate positions
features = bam_helper.getFeatures(tumour_tuples, normal_tuples)

if len(features) != 0:
    ## compute the probabilities
    probabilities = bam_helper.predict(features)
    
    ## generate the output string 
    outstr = bam_helper.makeOutStr()
    
    ## print the output string to out
    bam_helper.printResults(out, outstr)

else: 
    ## print N/A if there are no candidates
    for i in xrange(pos.start, pos.stop):
            print >> out, chrom + "\t" + str(i) + "\t" + "N/A\t" * 6 
#==============================================================================
# end of the main body
#==============================================================================
try:
#    expfile.close()
    out.close()
except:
    pass
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " done."
        
               
        
        
        
        

        

#==============================================================================
#       export features
#==============================================================================
        if args.export is not None or args.features_only:
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " exporting features"
            for p in xrange(u_pos-l_pos):
                print >> expfile, chrom + "\t" + str(l_pos + p) + "\t" + ("\t").join(map(str,batch[p]))
                
#==============================================================================
#       remove nan/inf values
#==============================================================================
        if not args.features_only:
            batch = numpy.array(batch)
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " removing potential nan/inf values"
            batch = removeNanInf(coords, batch, strings, info_strs)
            #print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        
#==============================================================================
#       filter and print the results to out
#==============================================================================
            if len(batch) != 0:            
                print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " predicting probabilities"
                results = model.predict_proba(batch)   
                print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " filtering and printing results"      
                filterAndPrintResult(coords, results, strings, info_strs)

            else: 
                results = None
                print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " filtering and printing results"                          
                for i in xrange(l_pos, u_pos):
                    print >> out, chrom + "\t" + str(i) + "\t" + "N/A\t" * 6                
#==============================================================================
# end of the main body
#==============================================================================
try:
#    expfile.close()
    out.close()
except:
    pass
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " done."

