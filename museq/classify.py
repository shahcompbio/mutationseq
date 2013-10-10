import sys
import bamutils
import argparse
import logging

mutationSeq_version="4.0.0"

logging.basicConfig(filename="mutationSeq_run.log", 
                    format='%(asctime)s %(message)s', 
                    #datefmt='%m/%d/%Y %I:%M:%S %p', 
                    level=logging.DEBUG)
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

samples = {}
for s in args.samples:
    samples[s.split(':')[0]] = s.split(':')[1]

info_str = " mutationSeq_" + mutationSeq_version + " started\n" + \
            " tumour:" + samples["tumour"] + \
            " normal:" + samples["normal"] + \
            " reference:" + samples["reference"] + \
            " --interval " + args.interval
logging.info(info_str)
#==============================================================================
# beginning of the main body
#==============================================================================
bam = bamutils.Bam(tumour=samples["tumour"], normal=samples["normal"], reference=samples["reference"])
bam_helper = bamutils.BamUtils(bam, args)

## get target positions "chromosome:start-stop"
logging.info("getting positions ...")
target_positions = bam_helper.get_positions()

## get an iterator over tuples of the target positions of tumour/normal bam files
logging.info("getting tuples ...")
tumour_tuples = bam.get_tumour_tuples(target_positions)
normal_tuples = bam.get_normal_tuples(target_positions)

## calculate features for the candidate positions
try:
    logging.info("extracting features ...")
    features = bam_helper.get_features(tumour_tuples, normal_tuples)

except KeyboardInterrupt:
    logging.info("pressed ctrl+c ...\n")
    sys.exit(1)

## predict the probabilities
logging.info("fitting model and predict probabilities ...")

if len(features) != 0:
    probabilities = bam_helper.predict(features)
else:
    probabilities = None
        
## print the output string to out
info_str = "printting results to:" + str(args.out) 
logging.info(info_str)

bam_helper.print_results(out, probabilities)

#==============================================================================
# end of the main body
#==============================================================================
try:
    #expfile.close()
    out.close()
except:
    pass
logging.info("successfully completed.\n")
