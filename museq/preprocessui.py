'''
Created on May 20, 2014

@author: dgrewal
'''
import argparse

mutationSeq_version="4.3.3"

#==============================================================================
# make a UI 
#==============================================================================
parser = argparse.ArgumentParser(prog='preprocess', 
                                 description = '''generate the allele counts for the 
                                 tumour at heterozygous positions in the normal.''')
## positional arguments                    
parser.add_argument("samples", 
                    nargs='*', 
                    help='''A list of colon delimited sample names; normal:normal.bam
                    tumour:tumour.bam model:model.npz reference:reference.fasta''')

## optional arguments
parser.add_argument("-b", "--buffer_size",
                    default="2G",
                    help='''specify max amount of memory usage''')

parser.add_argument("--coverage", 
                    default=0,
                    type=int,
                    help='''specify min depth (coverage) to be considered''')

parser.add_argument("-c", "--config", 
                    default=None, 
                    required=True,
                    help='''specify the path/to/metadata.config file used to add 
                          meta information to the output file''')

parser.add_argument("-i", "--interval",
                     default=None,
                     help='''specify an interval "chr[:start-stop]"''')

parser.add_argument("-t", "--threshold", 
                    default=0.5, type=float,
                    help='''set threshold for positive call''') 

parser.add_argument("-q", "--quality_threshold", 
                    default=0, 
                    type=int,
                    help='''set threshold for the mapping quality''')

parser.add_argument("-l", "--log_file",
                    default="mutationSeq_run.log",
                    help='''specify name or path of the log file''')

parser.add_argument("-v", "--verbose", 
                    action="store_true", default=False,
                    help='''verbose''')
                    
parser.add_argument("--version", 
                    action="version", version=mutationSeq_version)

## mandatory options                   
mandatory_options = parser.add_argument_group("required arguments")
mandatory_options.add_argument("-o", "--out", 
                               default=None, 
                               required=True, 
                               help='''specify the path/to/out.vcf to save output to a file''')  
                               
mandatory_options.add_argument("-f", "--positions_file", 
                     default=None,
                     required = True, 
                     help='''input a file containing a list of positions each of which in
                     a separate line, e.g. chr1:12345\nchr2:23456''')                  

args = parser.parse_args()
