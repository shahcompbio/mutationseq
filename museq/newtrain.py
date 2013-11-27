# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 12:21:18 2013

@author: jtaghiyar
"""
##TODO: needs a logger
import pybamapi
import numpy
import sys
import argparse
import matplotlib
matplotlib.use("Agg")
import newfeatures
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn import cross_validation
#from sklearn import tree
#from math import log
from sklearn.metrics import roc_curve, auc
from collections import defaultdict

mutationSeq_version="4.0.0"

#==============================================================================
# make a UI
#==============================================================================
parser = argparse.ArgumentParser(description='''train a model''')

parser.add_argument("infiles", 
                    metavar='FILE', nargs='*', default=[sys.stdin], type=argparse.FileType('r'), 
                    help= "A list of space delimited samples: chromosom position label, \
                    the file has a header of normal, tumour and reference file names")

parser.add_argument("--version", 
                    action="version", version=mutationSeq_version)
                    
parser.add_argument("--normalized", 
                    default=False, action="store_true", 
                    help="If you want to train with normalized features")
                    
parser.add_argument("--deep", 
                    default=False, action="store_true", 
                    help="If you want to test on deep data you need to change contamination rate")
                    
parser.add_argument("-o", "--out", 
                    default=None, 
                    help="save output to file")
                    
parser.add_argument("--validate", 
                    metavar='FILE', nargs='*', type=argparse.FileType('r'), default=None, 
                    help="To validat the same format file with known lable")
                    
parser.add_argument("--label", 
                    default="SOMATIC",
                    help="Label in training file list")

parser.add_argument("--export", 
                    default=None, 
                    help="save output to file")
                
parser.add_argument("--model",
                    default=None,
                    help="specify an existing model. Usually used for validation")

args = parser.parse_args()

#==============================================================================
# helper functions
#==============================================================================
def extract_labels(infiles):
    data=defaultdict(list)
    
    for case in infiles:
        tfile = None
        nfile = None
        rfile = None
        contamination = (float(30), float(30), float(70), float(0))

        for line in case:
            l = line.strip().split()
            if len(l) < 3:
                continue

            ## parse the line
            if l[0] == "#":
                if l[1] == "tumour":
                    tfile = l[2]
                elif l[1] == "normal":
                    nfile = l[2]
                elif l[1] == "reference":
                    rfile = l[2]
                elif l[1] == "contamination":
                    contamination = (float(l[2]), float(l[3]), float(l[4]), float(1))
                continue
            
            ## ignore cases where no reference is given
            if not all([tfile, nfile, rfile]):
                continue
        
            chromosome = l[0]
            position   = int(l[1])
           
            if l[2] == args.label:
                label = 1
            else:
                label = -1
                
            data[(tfile, nfile, rfile)].append((chromosome, position, label, contamination))
    
    return data
    
def extract_features(data):	
    features_buffer = []
    labels_buffer   = []
    keys_buffer    = []

    for tfile, nfile, rfile in data.keys():
        print "tumour:", tfile
        print "normal:", nfile
        t_bam = pybamapi.Bam(bam=tfile, reference=rfile, coverage=1)
        n_bam = pybamapi.Bam(bam=nfile, reference=rfile, coverage=1)
        
        for chromosome, position, label, c in data[(tfile, nfile, rfile)]:
            chromosome_id = t_bam.get_chromosome_id(chromosome)
            tt = t_bam.get_tuple(chromosome, position)
            nt = n_bam.get_tuple(chromosome, position)            
            rt = t_bam.get_reference_tuple(chromosome_id, position)            
            
            if not all([tt, nt, rt]):
                print "None tuple"
                continue
            
            ## calculate features            
            feature_version = newfeatures.version
            feature_set = newfeatures.Features(tt, nt, rt)
            temp_features = feature_set.get_features()   
            
            features_buffer.append(temp_features)
            labels_buffer.append(label)
            keys_buffer.append((rfile, nfile, tfile, chromosome, position, label))
            
    features_buffer = numpy.array(features_buffer)
    labels_buffer = numpy.array(labels_buffer)
    keys_buffer = numpy.array(keys_buffer)
    
    return feature_version, features_buffer, labels_buffer, keys_buffer
		                
#==============================================================================
# beginnig of the main body
#==============================================================================
if args.deep:
    c = (float(10000), float(10000), float(70), float(0))

if args.model is None:
    data=extract_labels(args.infiles)
    feature_version, m_features, labels, keys = extract_features(data)
    numpy.savez(args.out, feature_version, m_features, labels)
    
else:
    npz = numpy.load(args.model)
    m_features = npz["arr_1"]
    labels = npz["arr_2"]

model = RandomForestClassifier(random_state=0, n_estimators=3000, n_jobs=-1, compute_importances=True)     
model.fit(m_features, labels)

##==============================================================================
## extra stuff
##==============================================================================
## save important features
importance_file = open(args.out + "_importance.txt", 'w')

for importance, _feature in sorted(zip(model.feature_importances_, feature_set + coverage_features + extra_features)):
    print >> importance_file, _feature[0], importance
importance_file.close()


## validation
if args.validate:
    validation = extract_labels(args.validate)
    _, val_features, val_labels, val_keys = extract_features(validation)
    
    ## predict probabilities
    probs = model.predict_proba(val_features)
    voted = probs[:,1]
    fpr, tpr, thresholds = roc_curve(val_labels, voted)
    roc_auc = auc(fpr, tpr)
    print roc_auc
    fd = open(args.out + '_result.txt', 'w')
    
    for f, k, vf in zip(voted, val_keys, val_features):
        vf = " ".join(map(str, vf))        
        k = " ".join(map(str, k))
        print >> fd, k+ " "+ str(f)+ " "+ vf
    
    fd.close()  
    plt.plot(fpr,tpr,'k--', label='ROC curve (area = %0.3f)' %float(roc_auc))
    plt.title('ROC curve (area = %0.3f)' %float(roc_auc))
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.savefig(args.out + "_roc.png")
    #pylab.plot(fpr, tpr, 'k--', label='ROC curve (area = %0.3f)' %float(roc_auc))
    #pylab.plot([0, 1], [0, 1], 'k--')
    #pylab.xlim([0.0, 1.0])
    #pylab.ylim([0.0, 1.0])
    #pylab.xlabel('False Positive Rate')
    #pylab.ylabel('True Positive Rate')
    #pylab.savefig(args.out + "_roc.png")

## cross validation
else:
    cv = cross_validation.StratifiedKFold(labels, n_folds=3)
    for i, (train, test) in enumerate(cv):
        _model = model.fit(features[train], labels[train])

        probs = _model.predict_proba(features[test])
        voted = probs[:,1]
        fpr, tpr, thresholds = roc_curve(labels[test], voted)
        roc_auc = auc(fpr, tpr)
        fd = open(args.out + '_result_' + str(i), 'w')
        
        for f, k in zip(voted, keys[test]):
            print >> fd, k[0]+' '+k[1]+' '+k[2]+' '+k[3]+' '+k[4]+' '+k[5]+' '+str(f)

        fd.close()
        plt.plot(fpr, tpr, 'k--', lw=1, label="Fold %i (AUC=%0.3f)" % (i + 1, float(roc_auc)))

    plt.legend(loc="lower right", numpoints=1,)
    plt.savefig(args.out + "_rocxval.png")

#scores = []
#candidates = model.estimators_
#for _ in candidates:
#    scores.append(0)
#
#for f, label in zip(features, labels):
#    for estimator, index in zip(candidates, range(len(candidates))):
#        if estimator.predict(f)[0] == label:
#            scores[index] += 1 

#winner = candidates[numpy.argmax(numpy.array(scores))]
#tree.export_graphviz(winner, out_file=open(args.out + "_tree.dot", 'w'),)
