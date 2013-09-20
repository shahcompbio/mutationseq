import argparse
from datetime import datetime
import sys
from warnings import warn

mutationSeq_version="3.4.1"

#==============================================================================
# making a UI 
#==============================================================================
parser = argparse.ArgumentParser(prog='mutationSeq', description = '''mutationSeq:
                                 a feature-based classifier for somatic mutation detection in
                                 tumour-normal paired sequencing data''')
parser.add_argument("samples", nargs='*', help='''
                    A list of colon delimited sample names; normal:normal.bam
                    tumour:tumour.bam model:model.npz reference:reference.fasta''')

parser.add_argument("--version", action="version", version=mutationSeq_version)
mandatory_options = parser.add_argument_group("required arguments")
mandatory_options.add_argument("-o", "--out", default=None, #required=True, 
                               help="specify the path/to/out.vcf to save output to a file")
mandatory_options.add_argument("-c", "--config", default=None, #required=True,
                               help="specify the path/to/metadata.config file used to add \
                               meta information to the output file")
                               
exgroup = parser.add_mutually_exclusive_group()
exgroup.add_argument("-f", "--positions_file", default=None, 
                     help="input a file containing a list of positions each of which in a \
                     separate line, e.g. chr1:12345\nchr2:23456")
exgroup.add_argument("-i", "--interval", default=None, 
                     help="classify given chromosome[:start-end] range")

args = parser.parse_args()

#==============================================================================
# check the input 
#==============================================================================
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

#==============================================================================
# import required modules here to save time when only checking version or help
#==============================================================================
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " importing required modules"
#import pybam
#import old_pybam
#import numpy
from collections import deque, defaultdict, namedtuple
#from sklearn.ensemble import RandomForestClassifier
#from math import log10
from string import Template
import BamClass_newPybam

#==============================================================================
# helper functions
#==============================================================================
def parseTargetPos(poslist):
    target = poslist.split('\t')[0]
    try:
        target = target.split('r')[1] #check if "chr" is used
    except:
        pass
    try:
        pos = poslist.split('\t')[1]
        l_pos = int(pos.split('-')[0])
        try:
            u_pos = int(pos.split('-')[1])
        except:
            u_pos = int(l_pos)
        return [target, l_pos, u_pos]
    except:
        return [target, None, None]

def getAlt(ref_base, major, minor):
    if ref_base == major:
        print ref_base, major, minor
        return minor
    else:
        print ref_base, major, minor
        return major           

def printResult(outstrs):
    for outstr in outstrs:
        outstr = map(str, outstr)
        outstr = OutStr._make(outstr)
        info_str = "PR=0" + ";TR=" + outstr.TR + ";TA=" + outstr.TA +\
        ";NR=" + outstr.NR + ";NA=" + outstr.NA + ";NI=" + outstr.NI + ";ND=" + outstr.ND
        phred_qual = 0
        print >> out, outstr.chrom + "\t" + outstr.pos + "\t" + "." + "\t" + outstr.ref + "\t" +\
              outstr.alt + "\t" + "%.2f" % phred_qual + "\t" + outstr.filter + "\t" + info_str

def makeNamedTuple(tuple):
    tuple = Tuple._make(tuple)
    return tuple

def cacheTuple(tuples_buffer, tuple):
    tuples_buffer.append(tuple)

def cacheOutStr(outstr_buffer, outstr):
    outstr_buffer.append(outstr)

def printMetaData():
    try:
        cfg_file = open(args.config, 'r')
        tmp_file = ""
        for l in cfg_file:
            l = Template(l).substitute(DATETIME=datetime.now().strftime("%Y%m%d"),
                                       VERSION=mutationSeq_version,
                                       REFERENCE="N\A", #samples["reference"],
                                       TUMOUR=samples["tumour"],
                                       NORMAL=samples["normal"],
                                       THRESHOLD="N/A"#args.threshold
                                       )
            tmp_file += l
        cfg_file.close()
        print >> out, tmp_file,
    except:
        warn("Failed to load metadata file")
#==============================================================================
# start of the main body
#==============================================================================
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " mutationSeq_" + mutationSeq_version + " started"

bases = ('A', 'C', 'G', 'T', 'N')
samples = {}
for sample in args.samples:
    samples[sample.split(':')[0]] = sample.split(':')[1]

## check whether it is single sample analysis and if it is a normal("n") or tumour("t") sample
#if "normal" not in samples:
#    single_flag = 1
#    single_type = "t"
#elif "tumour" not in samples:
#    single_flag = 1
#    single_type = "n"
#else:
#    single_flag = 0 
#    single_type = None
#    
#if args.deep:
#    deep_flag = 1
#else:
#    deep_flag = 0
#
#Flags = namedtuple("Flags", "deep, single, type")
#flags = Flags._make([deep_flag, single_flag, single_type])

Tuple = namedtuple("Tuple", "position, A_tuple, C_tuple, G_tuple, T_tuple, all_tuple,\
                   major, minor, ambiguous_reads, insertion, entropy, deletion, ID")
    
OutStr = namedtuple("OutStr", "chrom, pos, ref, alt, filter, TR, TA, NR, NA, NI, ND")

#==============================================================================
# parse the positions or the file of list of positions to get targets
#==============================================================================
targets = None
l_pos = None
target_positions = defaultdict(list)

if args.interval:
    tmp_tp = parseTargetPos(args.interval)
    target_positions[tmp_tp[0]].append([tmp_tp[1], tmp_tp[2]])
elif args.positions_file:
    try:
        pos_file = open(args.positions_file, 'r')
        for l in pos_file.readlines():
            tmp_tp = parseTargetPos(l.strip())
            target_positions[tmp_tp[0]].append([tmp_tp[1], tmp_tp[2]])
        pos_file.close()
    except:
        print >> sys.stderr, "\tFailed to load the positions file from " + args.positions_file
        sys.exit(1)

else:
    print >> sys.stderr, "\tNo position/interval was specified"
    sys.exit(1)

#==============================================================================
# run for each chromosome/position
#==============================================================================
if args.config is not None:
    printMetaData()
bam = BamClass_newPybam.Bam(tumour_bam=samples["tumour"],
                            normal_bam=samples["normal"],
                            reference =samples["reference"])

for chrom in target_positions.keys(): # each key is a chromosome
    bam.setChromosome(chrom)
    
    for pos in xrange(len(target_positions[chrom])):
        l_pos = target_positions[chrom][pos][0]
        u_pos = target_positions[chrom][pos][1]

        if l_pos is None:
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + \
            " reading chromosome " + chrom 
        elif l_pos > u_pos:
            print >> sys.stderr, "bad input: lower bound is greater than or equal to upper bound in the interval"
            continue
        else:
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + \
            " reading chromosome " + chrom + " at positions " + str(l_pos) + " to " + str(u_pos)
    
        tuples_buffer = deque()
        outstr_buffer = []
         
        tumour_tuples = bam.getTumourTuplesIter(l_pos, u_pos)         
        for tt in tumour_tuples:
            tt = makeNamedTuple(tt)
            #if tt.insertion > 2 or tt.deletion > 2: # check the insertions/deletions
            #    continue
            cacheTuple(tuples_buffer, tt)

        if len(tuples_buffer) == 0:
            continue
        tt = tuples_buffer.popleft()

        normal_tuples = bam.getNormalTuplesIter(l_pos, u_pos)
        for nt in normal_tuples:
            nt = makeNamedTuple(nt)
            ## find positions where tuples for both tumour and normal exist
            while tt.position < nt.position:
                if len(tuples_buffer) == 0:
                    break
                tt = tuples_buffer.popleft()

            if tt.position != nt.position:
                continue
            
            ## generate the output strings
            ref_base = bam.getRefBase(tt.position)
            alt = getAlt(ref_base, tt.major, tt.minor)
            if tt.insertion > 0 or tt.deletion > 0:
                filter_flag = "INDL"
            else:
                filter_flag = "PASS"
            info_TR = tt[ref_base + 1][0]
            info_NR = nt[ref_base + 1][0]
            info_TA = tt[alt + 1][0]
            info_NA = nt[alt + 1][0]
            outstr = OutStr._make([chrom, tt.position, bases[ref_base], bases[alt],
                        filter_flag, info_TR, info_TA, info_NR, info_NA, tt.insertion, tt.deletion])
            cacheOutStr(outstr_buffer, outstr)
            if len(tuples_buffer) == 0:
                break

        ## print results 
        printResult(outstr_buffer)

#==============================================================================
# end of the main body
#==============================================================================
try:
#    expfile.close()
    out.close()
except:
    pass
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " done."

