import argparse
from datetime import datetime
import sys
from warnings import warn

mutationSeq_version="3.1.1"

#==============================================================================
# making a UI 
#==============================================================================
parser = argparse.ArgumentParser(description='''classify a dataset''')
parser.add_argument("samples", nargs='*', help='''
                    A list of colon delimited sample names; normal:normal.bam
                    tumour:tumour.bam model:model.npz reference:reference.fasta''')
parser.add_argument("-o", "--out", default=None, help="save output to file")
parser.add_argument("-t", "--threshold", default=0.5, help="set threshold for positive call", type=float)
parser.add_argument("-i", "--interval", default=None, help="classify given chromosome[:start-end] range")
parser.add_argument("-a", "--all", default=None, 
                    help= "force to print out even if the position(s) does not satisfy the initial criteria for Somatic calls")
parser.add_argument("-f", "--positions_file", default=None, 
                    help="input a file containing a list of positions each of which in a separate line, e.g. chr1:12345\nchr2:23456")
parser.add_argument("--export", default=None, help="save exported feature vector")
parser.add_argument("-n", "--normalized", default=False, action="store_true",
                    help="If you want to test with normalized features(the number of features are also ifferent from non-deep)")
parser.add_argument("-v", "--verbose", action="store_true", default=False)
parser.add_argument("--version", action="version", version=mutationSeq_version)
parser.add_argument("-p", "--purity", default=70, help="pass sample purity to features")
parser.add_argument("-c", "--config", default=None,
                    help="Specify the path/to/metadata.config file used to add meta information to the output file")
args = parser.parse_args()

#==============================================================================
# check the input 
#==============================================================================
if len(args.samples) != 4:
    print >> sys.stderr, "bad input, should follow: 'classify.py normal:<normal.bam> \
    tumour:<tumour.bam> reference:<ref.fasta> model:<model.npz> [--options]'"
    sys.exit(1)
    
if args.out is None:
    warn("--out is not specified, standard output is used to write the results")
    out = sys.stdout
else:
    out = open(args.out, 'w')

    
if args.config is None:
    warn("--config is not specified, no meta information used for the output VCF file")

if args.all is not None and args.all not in ("yes", "no"):
    print >> sys.stderr, "bad input for -a/--all, please specify only 'yes' or 'no'"
    sys.exit(1)

#==============================================================================
# import required modules here to save time when only checking version or help
#==============================================================================
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " importing required modules"
import pybam
import numpy
#import scipy
#import pickle
import features
import Nfeatures
from collections import deque, defaultdict
from sklearn.ensemble import RandomForestClassifier
from math import log10
from string import Template

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
        
def removeNanInf(coords, batch, strings, info_strs):
    b = []
    i = -1
    for coord, l, string in zip(coords, batch, strings):
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

def filterAndPrintResult(coords, batch, strings, info_strs):
    for coord, result, string, info in zip(coords, model.predict_proba(batch), strings, info_strs):
        if result[1] >= args.threshold:
            try:
                phred_qual = -10 * log10(1 - result[1])
            except:
                phred_qual = 99
            info_str = "PR=" + "%.3f" % result[1] + ";TR=" + str(info[1]) + ";TA=" + str(info[2]) + \
            ";NR=" + str(info[3]) + ";NA=" + str(info[4]) + ";TC=" + string[6]
            print >> out, str(string[0]) + "\t" + str(coord[0]) + "\t" + "." + "\t" + string[1] + "\t" \
            + bases[info[0]] + "\t"+ "%.2f" % phred_qual + "\t" + "PASS" + "\t" + info_str
        elif args.all == "yes":
            phred_qual = 0
            info_str = "PR=" + "%.3f" % result[1] + ";TR=" + str(info[1]) + ";TA=" + str(info[2]) + \
            ";NR=" + str(info[3]) + ";NA=" + str(info[4]) + ";TC=" + string[6]
            print >> out, str(string[0]) + "\t" + str(coord[0]) + "\t" + "." + "\t" + string[1] + "\t" \
            + bases[info[0]] + "\t"+ "%.2f" % phred_qual + "\t" + "FAIL" + "\t" + info_str

def extractFeature(tumour_data, normal_data, ref_data):
    if not args.normalized:
        feature_set = features.feature_set
        coverage_features = features.coverage_features
#        version = features.version
#        c = (float(30), float(30), float(70), float(0))
    else:
        feature_set = Nfeatures.feature_set
        coverage_features = Nfeatures.coverage_features
#        version = Nfeatures.version
#        c = (float(10000), float(10000), float(70), float(0))

    features_tmp = []
    for _, feature in feature_set:
        features_tmp.append(feature(tumour_data, normal_data, ref_data))

    coverage_data = (30, 30, int(args.purity), 1)
    for _, feature in coverage_features:
        features_tmp.append(feature(tumour_data, normal_data, coverage_data))
    n_counts = (normal_data[1][1], normal_data[2][1], normal_data[3][1],
                normal_data[4][1], normal_data[5][1])
    t_counts = (tumour_data[1][1], tumour_data[2][1], tumour_data[3][1],
                tumour_data[4][1], tumour_data[5][1])
    features_tmp.append(n.xentropy(n_counts, t_counts))
    return features_tmp    

#==============================================================================
# start of the main body
#==============================================================================
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " mutationSeq_" + mutationSeq_version + " started"
deep_flag = 0
bases = ('A', 'C', 'G', 'T', 'N')
samples = {}
for sample in args.samples:
    samples[sample.split(':')[0]] = sample.split(':')[1]

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
model = RandomForestClassifier(random_state=0, n_estimators=1000, n_jobs=-1, compute_importances=True)
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " model fitting started"
model.fit(train, labels)
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " done model fitting"

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
try:
    n = pybam.Bam(samples["normal"])
except:
    print >> sys.stderr, "\tFailed to load normal"
    sys.exit(1)
try:
    t = pybam.Bam(samples["tumour"])
except:
    print >> sys.stderr, "\tFailed to load tumour"
    sys.exit(1)
try:
    f = pybam.Fasta(samples["reference"])
except:
    print >> sys.stderr, "\tFailed to load reference"
    sys.exit(1)

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
    targets = set(n.targets) & set(t.targets)
    
#==============================================================================
# add VCF format meta-information from metadata.cfg file to the output file
#==============================================================================
if args.config:
    try:
        cfg_file = open(args.config, 'r')
        tmp_file = ""
        for l in cfg_file:
            l = Template(l).substitute(DATETIME=datetime.now().strftime("%Y%m%d"),
                                       REFERENCE=samples["reference"],
                                       TUMOUR=samples["tumour"],
                                       NORMAL=samples["normal"],
				       THRESHOLD=args.threshold)
            tmp_file += l
        cfg_file.close()
        print >> out, tmp_file,
    except:
        warn("Failed to load metadata file")
        
#==============================================================================
# run for each chromosome/position
#==============================================================================
for chrom in target_positions.keys(): # each key is a chromosomes
    for pos in range(len(target_positions[chrom])):
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
        
        try:
            f.load(chrom)
        except:
            message = "\treference does not have chromosome " + chrom
            warn(message)
            continue
        print >> sys.stderr, "\treading tumour data"
        if l_pos is None:
            g = t.vector(chrom, deep_flag)        
        else:
            g = t.vector(chrom, deep_flag, l_pos, u_pos)
            args.all = "yes"
    
        print >> sys.stderr, "\tnominating mutation positions in tumour"
        for tumour_data in g:
            position = tumour_data[0]
            ref_data = f.vector(int(position))
            
            ## get tri-nucleotide context, seems not really an efficient way, but tested it is indeed readlly fast
            if tumour_data[5][0] - tumour_data[ref_data[0] + 1][0] > 2 or args.all == "yes":
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
        
        if len(positions) == 0:
            continue
        position, tumour_data, ref_data, tc = positions.popleft()
    
        print >> sys.stderr, "\treading normal"  
        if l_pos is None:
            g = n.vector(chrom, deep_flag)
        else:
            g = n.vector(chrom, deep_flag, l_pos, u_pos)
    
        skip = False
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
            
#==============================================================================
#             feature extraction
#==============================================================================
            features_tmp = extractFeature(tumour_data, normal_data, ref_data)
            batch.append(features_tmp)
            coords.append((position, ref_data[0], normal_data[6], normal_data[normal_data[6] + 1][0],
                           normal_data[7], normal_data[normal_data[7] + 1][0], normal_data[11],
                           tumour_data[6], tumour_data[tumour_data[6] + 1][0], tumour_data[7],
                           tumour_data[tumour_data[7] + 1][0], tumour_data[11]))
            strings.append((chrom, bases[ref_data[0]], bases[normal_data[6]], bases[normal_data[7]], \
            bases[tumour_data[6]], bases[tumour_data[7]], tc))       
    
            ## find the ALT 
            if bases[ref_data[0]]==bases[tumour_data[6]]:
                alt = tumour_data[7]
            else:
                alt = tumour_data[6]        
    
            ## generate the values of info fields in the vcf output
            if alt != ref_data[0]:
                info_strs.append((alt, int(tumour_data[ref_data[0] + 1][0]), int(tumour_data[alt + 1][0]), \
                int(normal_data[ref_data[0] + 1][0]), int(normal_data[alt + 1][0])))
            else: # take care of the non-somatic positions
                info_strs.append((alt, int(tumour_data[ref_data[0] + 1][0]), 0, int(normal_data[ref_data[0] + 1][0]), 0))
    
            if len(positions) == 0:
                break
        batch = numpy.array(batch)
    
#==============================================================================
#       remove nan/inf values
#==============================================================================
        print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " removing potential nan/inf values"
        batch = removeNanInf(coords, batch, strings, info_strs)
    
#==============================================================================
#       filter and print the results to out
#==============================================================================
        print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " filtering and printing results"      
        filterAndPrintResult(coords, batch, strings, info_strs)
    print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " done printing"        
#        print >> sys.stderr, "***No position has been nominated (does not satisfy initial criteria for Somatic calls )"
#==============================================================================
# end of the main body
#==============================================================================
print >> sys.stderr, datetime.now().strftime("%H:%M:%S") + " done."
    
