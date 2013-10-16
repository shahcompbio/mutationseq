# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 12:21:18 2013

@author: jtaghiyar
"""

import bamutils
import numpy
import features
import Nfeatures
import sys
import argparse
#import pylab
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
from sklearn import cross_validation
#from sklearn.metrics import roc_curve, auc, confusion_matrix
from sklearn import tree
from math import log

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
                    
parser.add_argument("--C", 
                    default=False, action="store_true", 
                    help="If you want to test on deep data you need to change containation rate")
                    
parser.add_argument("--out", 
                    default=None, 
                    help="save output to file")
                    
parser.add_argument("--validate", 
                    metavar='FILE', nargs='*', type=argparse.FileType('r'), default=None, 
                    help="To validat the same format file with known lable")
                    
parser.add_argument("--label", 
                    default="SOMATIC",
                    help="Label in training file list")

args = parser.parse_args()

#==============================================================================
# helper functions
#==============================================================================
def extract_data(infiles):
    data={}	
    for case in infiles:
        #print >> sys.stderr, case
        rfile = None
        nfile = None
        tfile = None

        ## contamination
        c = (float(30), float(30), float(70), float(0))

        f = case
        for line in f:
            l = line.strip().split()
            if len(l) < 3:
                continue

            if line[0] == '#':
                if l[1] == "tumour":
                    tfile = l[2]
                elif l[1] == "normal":
                    nfile = l[2]
                elif l[1] == "reference":
                    rfile = l[2]
                elif l[1] == "contamination":
                    c = (float(l[2]), float(l[3]), float(l[4]), float(1))
                continue
       
            if not rfile:
                continue
        
            if rfile not in data:
                data[rfile] = {}
                
            chromosome = l[0]
            position = l[1]
            
            if chromosome not in data[rfile]:
                data[rfile][chromosome] = {}
                
            if (nfile, tfile) not in data[rfile][chromosome]:
                data[rfile][chromosome][(nfile, tfile)] = []

            if l[2] == args.label:
                label = 1
            else:
                label = -1

            data[rfile][chromosome][(nfile, tfile)].append((int(position), label, c))
            
    return data
 
def xentropy(tumour_counts, normal_counts):
        total_tc = tumour_counts[4]
        total_nc = normal_counts[4]
        ent = 0 # entropy
        
        for i in xrange(4):
            base_probability_tumour = tumour_counts[i] / total_tc
            base_probability_normal = normal_counts[i] / total_nc            
            if base_probability_tumour != 0:
                if base_probability_normal == 0:
                    ent -= -7 * base_probability_tumour
                else:
                    ent -= log(base_probability_normal) * base_probability_tumour
        return ent
           
def extract_features(data):	
    features_buffer = []
    labels_buffer = []
    keys_buffer = []

    for rfile in data.keys():
        for chromosome in data[rfile].keys():
            for nfile, tfile in data[rfile][chromosome].keys():
                print tfile
                print nfile
                bam = bamutils.Bam(tumour=tfile, normal=nfile, reference=rfile)
                for (position, label, c) in data[rfile][chromosome][(nfile, tfile)]:
                    temp_features = []
                    
                    chromosome_id = bam.get_tumour_chromosome_id(chromosome)
                    rt = bam.get_reference_tuple(chromosome_id, position)
                    tt = bam.get_normal_tuple(chromosome, position)
                    nt = bam.get_tumour_tuple(chromosome, position)
                    
                    if tt is None or nt is None or rt is None:
                        print "None tuple"
                        continue
                    
                    for _, feature_func in feature_set:
                        temp_features.append(feature_func(tt, nt, rt))
                    for _, feature in coverage_features:
                        temp_features.append(feature(tt, nt, c))
                    t_counts = (tt[1][0], tt[2][0], tt[3][0], tt[4][0], tt[5][0])
                    n_counts = (nt[1][0], nt[2][0], nt[3][0], nt[4][0], nt[5][0])
                    temp_features.append(xentropy(n_counts, t_counts))
                    features_buffer.append(temp_features)
                    labels_buffer.append(label)
                    keys_buffer.append((rfile, nfile, tfile, chromosome, position, label))
                    
    keys_buffer=numpy.array(keys_buffer)
    features_buffer = numpy.array(features_buffer)
    labels_buffer = numpy.array(labels_buffer)
    return features_buffer, labels_buffer, keys_buffer
		                
#==============================================================================
# parse input arguments                  
#==============================================================================
if not args.normalized:
    feature_set = features.feature_set
    coverage_features = features.coverage_features
    extra_features = (("xentropy", 0), ("SENTINEL", 0))
    version = features.version
    c = (float(30), float(30), float(70), float(0))
    if args.C:
        c = (float(100000), float(100000), float(70), float(0))
        print "Warning: you are trying to test deep data with unnormalized features"
else:
    feature_set = Nfeatures.feature_set
    coverage_features = Nfeatures.coverage_features
    extra_features = (("xentropy", 0), ("SENTINEL", 0))
    version = Nfeatures.version
    c = (float(30), float(30), float(70), float(0))
    if args.C:
        c = (float(10000), float(10000), float(70), float(0))
        
#==============================================================================
# beginnig of the main body
#==============================================================================
data=extract_data(args.infiles)
feature, labels, keys = extract_features(data)
numpy.savez(args.out, version, feature, labels)

model = RandomForestClassifier(random_state=0, n_estimators=3000, n_jobs=-1, compute_importances=True)
model.fit(feature, labels)

##==============================================================================
## extra stuff
##==============================================================================
#print "saving the model with pickle ..."
#saveObject(model, "model2.npz")
#print "saving done!"
##features_file=open("features_name.txt", 'w')
##for _feature in feature_set + coverage_features + extra_features:
#    #print >> _feature[0]
#    
##features_file.close()
#
#
#importance_file = open(args.out + "_importance.txt", 'w')
#
##for importance, _feature in (zip(model.feature_importances_, feature_set + coverage_features + extra_features)):
##    print >> importance_file, 
#for importance, _feature in sorted(zip(model.feature_importances_,
#                                       feature_set + coverage_features + extra_features)):
#    print >> importance_file, _feature[0], importance
#importance_file.close()
#
#if args.validate:
#    validation=extract_data(args.validate)
#    val_features, val_labels, val_keys=extract_features(validation)
#    probs = model.predict_proba(val_features)
#    voted = probs[:,1]
#    fpr, tpr, thresholds = roc_curve(val_labels, voted)
#    roc_auc = auc(fpr, tpr)
#    print roc_auc
#    fd=open(args.out+'_result', 'w')
#    for f,k in zip(voted,val_keys):
#        print >> fd, k[0]+' '+k[1]+' '+k[2]+' '+k[3]+' '+k[4]+' '+k[5]+' '+str(f)
#    fd.close()  
#    plt.plot(fpr,tpr,'k--', label='ROC curve (area = %0.3f)' %float(roc_auc))
#    plt.title('ROC curve (area = %0.3f)' %float(roc_auc))
#    plt.xlabel('False Positive Rate')
#    plt.ylabel('True Positive Rate')
#    plt.savefig(args.out + "_roc.png")
#    #pylab.plot(fpr, tpr, 'k--', label='ROC curve (area = %0.3f)' %float(roc_auc))
#    #pylab.plot([0, 1], [0, 1], 'k--')
#    #pylab.xlim([0.0, 1.0])
#    #pylab.ylim([0.0, 1.0])
#    #pylab.xlabel('False Positive Rate')
#    #pylab.ylabel('True Positive Rate')
#    #pylab.savefig(args.out + "_roc.png")
#
#else:
#    cv = cross_validation.StratifiedKFold(labels, n_folds=3)
#    for i, (train, test) in enumerate(cv):
#        _model = model.fit(feature[train], labels[train])
#
#        probs = _model.predict_proba(feature[test])
#        voted = probs[:,1]
#        fpr, tpr, thresholds = roc_curve(labels[test], voted)
#        roc_auc = auc(fpr, tpr)
#        fd=open(args.out+'_result_'+str(i), 'w')
#        for f,k in zip(voted,keys[test]):
#            print >> fd, k[0]+' '+k[1]+' '+k[2]+' '+k[3]+' '+k[4]+' '+k[5]+' '+str(f)
#        fd.close()
#        plt.plot(fpr, tpr, 'k--', lw=1, label="Fold %i (AUC=%0.3f)" % (i + 1, float(roc_auc)))
#
#    
#    plt.legend(loc="lower right", numpoints=1,)
#    plt.savefig(args.out + "_rocxval.png")
#
#scores = []
#candidates = model.estimators_
#for _ in candidates:
#    scores.append(0)
#
#for f, label in zip(feature, labels):
#    for estimator, index in zip(candidates, range(len(candidates))):
#        if estimator.predict(f)[0] == label:
#            scores[index] += 1 
#
#winner = candidates[numpy.argmax(numpy.array(scores))]
#tree.export_graphviz(winner, out_file=open(args.out + "_tree.dot", 'w'),)
