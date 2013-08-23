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
#parser.add_argument("-s", "--single", default=None, help="specify the bam file for single sample analysis")
parser.add_argument("-a", "--all", default=None, choices=["no", "yes"], 
                    help= "force to print out even if the position(s) does not satisfy the initial criteria for Somatic calls")
parser.add_argument("-e" , "--export", default=None, help="save exported feature vector to the specified path")
parser.add_argument("-u", "--features_only", default=False, action="store_true", help="if true, only extracted features are exported")
parser.add_argument("-d", "--deep", default=False, action="store_true", help="for deepseq data")
parser.add_argument("-n", "--normalized", default=False, action="store_true",
                    help="If you want to test with normalized features(the number of features are also ifferent from non-deep)")
parser.add_argument("-p", "--purity", default=70, help="pass sample purity to features")
parser.add_argument("-v", "--verbose", action="store_true", default=False)
parser.add_argument("--version", action="version", version=mutationSeq_version)
parser.add_argument("-t", "--threshold", default=0.5, help="set threshold for positive call", type=float)
mandatory_options = parser.add_argument_group("required arguments")
mandatory_options.add_argument("-c", "--config", default=None, #required=True, 
                    help="specify the path/to/metadata.config file used to add meta information to the output file")
mandatory_options.add_argument("-o", "--out", default=None, #required=True, 
                               help="specify the path/to/out.vcf to save output to a file")
exgroup = parser.add_mutually_exclusive_group()
exgroup.add_argument("-f", "--positions_file", default=None, 
                   help="input a file containing a list of positions each of which in a separate line, e.g. chr1:12345\nchr2:23456")
exgroup.add_argument("-i", "--interval", default=None, help="classify given chromosome[:start-end] range")

args = parser.parse_args()

#==============================================================================
# check the input 
#==============================================================================
if len(args.samples) < 3:
    print >> sys.stderr, "bad input, usage: 'classify.py normal:<normal.bam> \
    tumour:<tumour.bam> reference:<ref.fasta> model:<model.npz> [--options]'"
    sys.exit(1)

if not args.features_only:
    if args.out is None:
        warn("--out is not specified, standard output is used to write the results")
        out = sys.stdout
    else:
        out = open(args.out, 'w')
    
    if args.config is None:
        warn("--config is not specified, no meta information used for the output VCF file")

else: 
    if args.export is None:
        warn("-u/--features_only is used without -e/--export to specify a path to save the features")
        expfile = sys.stdout
    else:    
        try:
            expfile = open(args.export, 'w')       
        except:
            print >> sys.stderr, "\tFailed to open the export file for writing: " + args.export
            expfile = sys.stdout
#==============================================================================
# import required modules here to save time when only checking version or help
#==============================================================================
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " importing required modules"
import os
import pybam
import numpy
#import scipy
#import pickle
import resource
import features
import Nfeatures
import features_single
import features_deep
from collections import deque, defaultdict, namedtuple
from sklearn.ensemble import RandomForestClassifier
from math import log10
from string import Template
#print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

#==============================================================================
# helper functions
#==============================================================================
def parseTargetPos(poslist):
    target = poslist.split(':')[0]
    try:
        target = target.split('r')[1] #check if "chr" is used
    except:
        pass
    try:
        pos = poslist.split(':')[1]
        l_pos = int(pos.split('-')[0])
        try:
            u_pos = int(pos.split('-')[1])
        except:
            u_pos = int(l_pos) + 1
        return [target, l_pos, u_pos]
    except:
        return [target, None, None]

def printMetaData():
    try:
        cfg_file = open(args.config, 'r')
        tmp_file = ""
        for l in cfg_file:
            l = Template(l).substitute(DATETIME=datetime.now().strftime("%Y%m%d"),
                                       VERSION=mutationSeq_version,
                                       REFERENCE=samples["reference"],
                                       TUMOUR=samples["tumour"],
                                       NORMAL=samples["normal"],
				       THRESHOLD=args.threshold)
            tmp_file += l
        cfg_file.close()
        print >> out, tmp_file,
    except:
        warn("Failed to load metadata file")

def getFeatureNames():
    feature_names = []
    if args.normalized:
        feature_set = Nfeatures.feature_set
        coverage_features = Nfeatures.coverage_features
        
    elif args.deep: 
        feature_set = features_deep.feature_set
        coverage_features = features_deep.coverage_features
    
    elif flags.single:
        feature_set = features_single.feature_set
        coverage_features = features_single.coverage_features
        
    else:
        feature_set = features.feature_set
        coverage_features = features.coverage_features

    for i in xrange(len(feature_set)):
        feature_names.append(feature_set[i][0])

    for i in xrange(len(coverage_features)):
        feature_names.append(coverage_features[i][0])

    return feature_names
       
def nominateMutPos(tumour_data):
    ref_data = f.vector(int(tumour_data[0])) # tumour_data[0] is position
    return tumour_data[5][0] - tumour_data[ref_data[0] + 1][0] > 2 
            
def getPositions(tumour_data):
    position = tumour_data[0]
    ref_data = f.vector(int(position))
    try:
        pre_refBase = f.vector(int(position) - 1)[0]
    except:
        pre_refBase = 4
    try:
        nxt_refBase = f.vector(int(position) + 1)[0]
    except:
        nxt_refBase = 4
    tri_nucleotide = bases[pre_refBase] + bases[ref_data[0]] + bases[nxt_refBase]
    positions.append((position, tumour_data, ref_data, tri_nucleotide))

def getAlt(base, ref_nuc, major, minor):
    if bases[ref_nuc] == bases[major]:
        return minor
    else:
        return major   
        
def extractFeature(tumour_data, ref_data, *arguments):
    if len(arguments) != 0:
        normal_data = arguments[0]
        
    if args.normalized:
        feature_set = Nfeatures.feature_set
        coverage_features = Nfeatures.coverage_features
#        extra_features = (("xentropy", 0), ("SENTINEL", 0)) # can be used for args.verbose
#        version = Nfeatures.version
        c = (float(30), float(30), int(args.purity), float(0))
        
    elif args.deep: 
        feature_set = features_deep.feature_set
        coverage_features = features_deep.coverage_features
#        extra_features = (("xentropy", 0), ("SENTINEL", 0))
#        version = features_deep.version 
        c = (float(10000), float(10000), int(args.purity), float(0)) 
    
    elif flags.single:
        feature_set = features_single.feature_set
        coverage_features = features_single.coverage_features
#        extra_features = (("xentropy", 0), ("SENTINEL", 0)) 
#        version = features_single.version
        c = (float(30), float(30), int(args.purity), float(0))
        
    else:
        feature_set = features.feature_set
        coverage_features = features.coverage_features
#        extra_features = (("xentropy", 0), ("SENTINEL", 0))
#        version = features.version
        c = (float(30), float(30), int(args.purity), float(0))
   
    features_tmp = []
    for _, feature in feature_set:
        if flags.single:
            features_tmp.append(feature(tumour_data, ref_data))
        else:
            features_tmp.append(feature(tumour_data, normal_data, ref_data))  

    coverage_data = c
    for _, feature in coverage_features:
        if flags.single:
            features_tmp.append(feature(tumour_data, coverage_data))
        else:
            features_tmp.append(feature(tumour_data, normal_data, coverage_data))

    t_counts = (tumour_data[1][1], tumour_data[2][1], tumour_data[3][1],
                tumour_data[4][1], tumour_data[5][1])
    if not flags.single:
        n_counts = (normal_data[1][1], normal_data[2][1], normal_data[3][1],
                    normal_data[4][1], normal_data[5][1])
        features_tmp.append(n.xentropy(n_counts, t_counts))

    return features_tmp   
    
def removeNanInf(coords, batch, strings, info_strs):
    b = []
    i = -1
    for l in batch:
        i = i + 1
        if numpy.isnan(l).any() or numpy.isinf(l).any():
            print >> sys.stderr, "\tnan/inf value removed"
            coords.pop(i)
            strings.pop(i)
            info_strs.pop(i)
            continue
        b.append(l)
    batch = numpy.array(b)
    return batch

def getInfoStr(alt, tumour_data, ref_base, *arguments):
    if len(arguments) != 0:
        normal_data = arguments[0]

    if alt != ref_base:
        if flags.single:
            info_strs.append((alt, int(tumour_data[ref_base + 1][0]), int(tumour_data[alt + 1][0])))
        else:
            info_strs.append((alt, int(tumour_data[ref_base + 1][0]), int(tumour_data[alt + 1][0]),
                              int(normal_data[ref_base + 1][0]), int(normal_data[alt + 1][0])))
    
    else: # take care of the non-somatic positions
        if flags.single:
            info_strs.append((alt, int(tumour_data[ref_base + 1][0]), 0))
        else:
            info_strs.append((alt, int(tumour_data[ref_base + 1][0]), 0, 
                              int(normal_data[ref_base + 1][0]), 0))
                              
def getOutStr(**kwargs):
    if flags.type == "n":
        rr = ";NR="
        aa = ";NA="
    else:
        rr = ";TR="
        aa = ";TA="
        
    if flags.single:
        info_str = "PR=" + "%.3f" % kwargs["PR"] + rr + str(kwargs["RR"]) + aa + \
        str(kwargs["AA"])+ ";TC=" + str(kwargs["TC"])
    else:
        info_str = "PR=" + "%.3f" % kwargs["PR"] + ";TR=" + str(kwargs["TR"]) + ";TA=" + \
        str(kwargs["TA"])+ ";NR=" + str(kwargs["NR"]) + ";NA=" + str(kwargs["NA"])+ ";TC=" + str(kwargs["TC"])
        
    out_str = str(kwargs["CHROM"]) + "\t" + str(kwargs["POS"]) + "\t" + kwargs["ID"] + "\t" + \
    kwargs["REF"] + "\t" + kwargs["ALT"] + "\t" + "%.2f" % kwargs["QUAL"] + "\t" + kwargs["FILTER"] + \
    "\t" + info_str

    return out_str
            
def filterAndPrintResult(coords, results, strings, info_strs):
    for coord, result, string, info in zip(coords, results, strings, info_strs):
        if result[1] >= args.threshold and coord[-1] <= 3:
            try:
                phred_qual = -10 * log10(1 - result[1])
            except:
                phred_qual = 99
            
#            info_str = "PR=" + "%.3f" % result[1] + ";TR=" + str(info[1]) + ";TA=" + str(info[2]) + \
#            ";NR=" + str(info[3]) + ";NA=" + str(info[4]) + ";TC=" + string[2]
#            print >> out, str(string[0]) + "\t" + str(coord[0]) + "\t" + "." + "\t" + string[1] + "\t" \
#            + bases[info[0]] + "\t"+ "%.2f" % phred_qual + "\t" + "PASS" + "\t" + info_str
            if flags.single:
                print >> out, getOutStr(PR=result[1], RR=info[1], AA=info[2],
                                        TC=string[2], CHROM=string[0], POS=coord[0],
                                        ID=".", REF=string[1], ALT=bases[info[0]],
                                        QUAL = phred_qual, FILTER="PASS")
            else:
                print >> out, getOutStr(PR=result[1], TR=info[1], TA=info[2], NR=info[3],
                                        NA=info[4], TC=string[2], CHROM=string[0], POS=coord[0],
                                        ID=".", REF=string[1], ALT=bases[info[0]],
                                        QUAL = phred_qual, FILTER="PASS")
        
        elif args.all == "yes":
            phred_qual = 0
            if flags.single:
                print >> out, getOutStr(PR=result[1], RR=info[1], AA=info[2],
                                        TC=string[2], CHROM=string[0], POS=coord[0],
                                        ID=".", REF=string[1], ALT=bases[info[0]],
                                        QUAL = phred_qual, FILTER="FAIL")
            else:
                print >> out, getOutStr(PR=result[1], TR=info[1], TA=info[2], NR=info[3],
                                        NA=info[4], TC=string[2], CHROM=string[0], POS=coord[0],
                                        ID=".", REF=string[1], ALT=bases[info[0]],
                                        QUAL = phred_qual, FILTER="FAIL")
       
#            info_str = "PR=" + "%.3f" % result[1] + ";TR=" + str(info[1]) + ";TA=" + str(info[2]) + \
#            ";NR=" + str(info[3]) + ";NA=" + str(info[4]) + ";TC=" + string[2]
#            print >> out, str(string[0]) + "\t" + str(coord[0]) + "\t" + "." + "\t" + string[1] + "\t" \
#            + bases[info[0]] + "\t"+ "%.2f" % phred_qual + "\t" + "FAIL" + "\t" + info_str
     
#==============================================================================
# start of the main body
#==============================================================================
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " mutationSeq_" + mutationSeq_version + " started"

bases = ('A', 'C', 'G', 'T', 'N')
samples = {}
for sample in args.samples:
    samples[sample.split(':')[0]] = sample.split(':')[1]

## check whether it is single sample analysis and if it is a normal("n") or tumour("t") sample
if "normal" not in samples:
    single_flag = 1
    single_type = "t"
elif "tumour" not in samples:
    single_flag = 1
    single_type = "n"
else:
    single_flag = 0 
    single_type = None
    
if args.deep:
    deep_flag = 1
else:
    deep_flag = 0

Flags = namedtuple("Flags", "deep, single, type")
flags = Flags._make([deep_flag, single_flag, single_type])    
#==============================================================================
# fit a model
#==============================================================================
#print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " pickle started"
try:
    npz = numpy.load(samples["model"])
#    model = pickle.load(open(samples["model"], 'rb'))
except:
    print >> sys.stderr, "\tFailed to load model"
    print >> sys.stderr, sys.exc_info()[0]
    sys.exit(1)
#print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + "done"
train = npz["arr_1"]
labels = npz["arr_2"]
model = RandomForestClassifier(random_state=0, n_estimators=1000, n_jobs=1, compute_importances=True)
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " fitting model"
model.fit(train, labels)
#print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

# if npz["arr_0"] != version:
#     print >> sys.stderr, "\tmismatched feature set versions:",
#     print >> sys.stderr, "\t" + str(npz["arr_0"]), "and", str(version)
#     out.close()
#     sys.exit(1)
#if args.verbose:
#    print >> sys.stderr, "\tunsorted feature names"
#    for importance, _feature in zip(model.feature_importances_, feature_set + coverage_features + extra_features):
#        print >> sys.stderr, _feature[0], importance
#
#    print >> sys.stderr, "\tsorted by importance:"
#    for importance, _feature in sorted(zip(model.feature_importances_, feature_set + coverage_features + extra_features)):
#        print >> sys.stderr, _feature[0], importance

#==============================================================================
# read in bam files
#==============================================================================
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " reading bam files"
if flags.single:

    if flags.type == "n":
        type_name = "normal" 
    else:
        type_name = "tumour"
        
    try:
        t = pybam.Bam(samples[type_name]) # single bam variable is treated the same as tumour bam variable
    except:
        print >> sys.stderr, "\tFailed to load " + type_name + " bam file"
        sys.exit(1)

else:
    try:
        n = pybam.Bam(samples["normal"])
    except:
        print >> sys.stderr, "\tFailed to load normal bam file"
        sys.exit(1)

    try:
        t = pybam.Bam(samples["tumour"])
    except:
        print >> sys.stderr, "\tFailed to load tumour bam file"
        sys.exit(1)

try:
    f = pybam.Fasta(samples["reference"])
except:
    print >> sys.stderr, "\tFailed to load reference"
    sys.exit(1)
#print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

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
elif flags.single:
    targets = set(t.targets)
else:
    targets = set(n.targets) & set(t.targets)    

if len(target_positions.keys()) == 0:
    for targ in targets:
        target_positions[targ].append([None, None])

if args.config and not args.features_only: # add VCF format meta-information
    printMetaData()

if args.export is not None or args.features_only: # add feature names to each column of the exported file
    print >> expfile, "CHR" + "\t" + "POS" +  "\t" + ('\t').join(getFeatureNames())
#==============================================================================
# run for each chromosome/position
#==============================================================================
for chrom in target_positions.keys(): # each key is a chromosome
    for pos in xrange(len(target_positions[chrom])):
        l_pos = target_positions[chrom][pos][0]
        u_pos = target_positions[chrom][pos][1]
        if l_pos is None:
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + \
            " reading chromosome " + chrom 
        elif l_pos >= u_pos:
            print >> sys.stderr, "bad input: lower bound is greater than or equal to upper bound in the interval"
            continue
        else:
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + \
            " reading chromosome " + chrom + " at positions " + str(l_pos) + " to " + str(u_pos)
    
        batch = []
        coords = []
        strings = []
        info_strs = []
        positions = deque([])
        
        print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " loading reference for chromosome " + chrom
        try:
            f.load(chrom)
        except:
            message = "\treference does not have chromosome " + chrom
            warn(message)
            continue
        #print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " reading tumour data"
        if l_pos is None:
            g = t.vector(chrom, flags.deep)        
        else:
            g = t.vector(chrom, flags.deep, l_pos, u_pos)
            if args.all is None:
                args.all = "yes"
            
        #print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " nominating mutation positions in tumour"
        if args.all == "yes":
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " mapping"
            map(getPositions, g)
        else:
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " filtering"
            g = filter(nominateMutPos, g)
            #print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " mapping"
            map(getPositions, g)   
        #print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        
        if flags.single:
            while True:
                if len(positions) == 0:
                    break
                position, tumour_data, ref_data, tc = positions.popleft()

                ## feature extraction                
                features_tmp = extractFeature(tumour_data, ref_data)
                batch.append(features_tmp)
                coords.append((position, ref_data[0], tumour_data[6], tumour_data[tumour_data[6] + 1][0], 
                               tumour_data[7],tumour_data[tumour_data[7] + 1][0], tumour_data[11]))
                strings.append((chrom, bases[ref_data[0]], tc))
    
                ## find the ALT and generate the info_strs used for INFO column in the output
                alt = getAlt(bases, ref_data[0], tumour_data[6], tumour_data[7])    
                getInfoStr(alt, tumour_data, ref_data[0])
                
        else:
            if len(positions) == 0:
                continue
            position, tumour_data, ref_data, tc = positions.popleft()
        
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " reading normal"  
            if l_pos is None:
                g = n.vector(chrom, flags.deep)
            else:
                g = n.vector(chrom, flags.deep, l_pos, u_pos)
                
            #print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            skip = False
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " extracting features"      
            for normal_data in g:
                while (position < normal_data[0]):
                    if len(positions) == 0:
                        skip = True
                        break
                    position, tumour_data, ref_data, tc = positions.popleft()
                
                if skip:
                    break
                if normal_data[0] != position:
                    continue
                ## feature extraction                
                features_tmp = extractFeature(tumour_data, ref_data, normal_data)
                batch.append(features_tmp)
                coords.append((position, ref_data[0], normal_data[6], normal_data[normal_data[6] + 1][0],
                               normal_data[7], normal_data[normal_data[7] + 1][0], normal_data[11],
                               tumour_data[6], tumour_data[tumour_data[6] + 1][0], tumour_data[7],
                               tumour_data[tumour_data[7] + 1][0], tumour_data[11]))
                strings.append((chrom, bases[ref_data[0]], tc))       
        
                ## find the ALT and generate the info_strs used for INFO column in the output
                alt = getAlt(bases, ref_data[0], tumour_data[6], tumour_data[7])    
                getInfoStr(alt, tumour_data, ref_data[0], normal_data)
                        
                if len(positions) == 0:
                    break
        #print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss    
                
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
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " predicting probabilities"
            results = model.predict_proba(batch)        
            print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " filtering and printing results"      
            filterAndPrintResult(coords, results, strings, info_strs)
        #print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        
    print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " done printing"        
    #print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss    
#        print >> sys.stderr, "***No position has been nominated (does not satisfy initial criteria for Somatic calls )"
#==============================================================================
# end of the main body
#==============================================================================
try:
    expfile.close()
    out.close()
except:
    pass
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " done."
#print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss    
