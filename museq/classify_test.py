import unittest
import pybamapi
import bamutils
import sys
import numpy
import scipy
import sklearn
import matplotlib

import logging
import os
import resource
import re
import features, features_single, features_deep, features_deep_single
import matplotlib.pyplot as plt
from math import log10
from sklearn.ensemble import RandomForestClassifier
from sklearn import cross_validation
from sklearn.metrics import roc_curve, auc
from string import Template
from datetime import datetime
from collections import defaultdict


class base_class():
    args = None
    def create_dict_tuples_paired(self,filename):
        tuples_dict_labels={}
        input_tt=[]
        input_nt=[]
        input_rt=[]
        input_position=[]
    
        file_stream = open(filename)
        for line in file_stream:
            line=line.strip().split(':')
            if line[0]=='pos':
                input_position=line[1]
            if line[0]=='tt':
                input_tt=line[1]
            if line[0]=='nt':
                input_nt=line[1]
            if line[0]=='rt':
                input_rt=line[1]
            if len(input_position)!=0 and len(input_tt)!=0 and len(input_nt)!=0 and len(input_rt)!=0:
                tuples_dict_labels[input_position]=[input_tt,input_nt,input_rt]
                input_position=[] 
                input_tt=[]
                input_nt=[]
                input_rt=[] 
        return tuples_dict_labels;
    
    def create_dict_tuples_single(self,filename):
        tuples_dict_labels={}
        input_it=[]
        input_rt=[]
        input_position=[]
    
        file_stream = open(filename)
        for line in file_stream:
            line=line.strip().split(':')
            if line[0]=='pos':
                input_position=line[1]
            if line[0]=='it':
                input_it=line[1]
            if line[0]=='rt':
                input_rt=line[1]
            if len(input_position)!=0 and len(input_it)!=0 and len(input_rt)!=0:
                tuples_dict_labels[input_position]=[input_it,input_rt]
                input_position=[] 
                input_it=[]
                input_rt=[] 
        return tuples_dict_labels
    
    def create_dict_vcf(self,filename):
        #create a dict from known inputs
        label_dict = {}
        file_stream = open(filename)
        for line in file_stream.readlines():
            line = line.strip().split()
            if line[0][0] == "#":
                continue
            else:
                label_dict[line[0]+":"+line[1]] = [line[3],line[4],line[5],line[6],line[7]]
        return label_dict
    
    def run_classifier(self,args):
        classifier = bamutils.Classifier(args)
        target_positions = classifier.get_positions()
        tuples = classifier.bam.get_tuples(target_positions)
        features = classifier.get_features(tuples)
        if args.export_features is not None:
            self.classifier.export_features(features)
        probabilities = classifier.predict(features)
        classifier.print_results(probabilities)
    
    def get_tuples(self,args):
        classifier = bamutils.Classifier(self.args)
        target_positions = classifier.get_positions()
        tuples = classifier.bam.get_tuples(target_positions)
        return tuples,classifier
    
    def compare_vcf_dicts(self,vcf_list1,vcf_list2):
        self.assertEqual(len(vcf_list1),len(vcf_list2),
                                     'The result generated doesn\'t match')
        for i in range(4):
            if i==2:
                self.assertAlmostEqual(float(vcf_list1[2]), float(vcf_list2[2]), 
                                       None,'The result generated doesn\'t match',0.5 )
            else:
                self.assertEqual(vcf_list1[i], vcf_list2[i], 
                                         'The result generated doesn\'t match')
                    #the info field
            vcf_list1_info=vcf_list1[4].strip().replace('=',';').split(';')
            vcf_list2_info=vcf_list2[4].strip().replace('=',';').split(';')
            for i in [1,3,5,7,9,11,13,15]:
                if i==1:
                    self.assertAlmostEqual(float(vcf_list2_info[i]), float(vcf_list1_info[i]), 
                                           None,'The result generated doesn\'t match', 0.5) 
                if i == 11:
                    self.assertEqual(vcf_list1_info[i], vcf_list2_info[i], 
                                           'The result generated doesn\'t match')
                else:
                    if vcf_list1_info[i]=="N/A" or vcf_list2_info[i] == "N/A":
                        continue
                    self.assertAlmostEqual(float(vcf_list1_info[i]), float(vcf_list2_info[i]), 
                                           None,'The result generated doesn\'t match', 0.5)

    def get_features_paired(self,tt,nt,classifier):
        chromosome_id = tt[-1]
        position = tt[0]
        refbase = classifier.bam.get_reference_base(chromosome_id, position, index=True)
        nonrefbases = [x for x in range(4) if x != refbase]
            
        ## ignore tumour tuples with no/few variants in the tumour or too many variants in the normal
        if not self.args.no_filter:
            if  tt[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
                tt[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
                tt[nonrefbases[2] + 1][0] < self.args.tumour_variant or \
                (nt[5][0] - nt[refbase + 1][0]) / nt[5][0] > (self.args.normal_variant / 100):
                    return None,refbase
           
        ## get corresponding reference tuple
        rt = classifier.bam.get_reference_tuple(chromosome_id, position)
        feature_set = classifier.features_module.Features(tt, nt, rt)
        return feature_set,refbase

        
    def get_features_single(self,it,classifier):
        chromosome_id = it[-1]
        position = it[0]
        refbase = classifier.bam.get_reference_base(chromosome_id, position, index=True)
        nonrefbases = [x for x in range(4) if x != refbase]
            
        ## ignore tumour tuples with no/few variants in the tumour or too many variants in the normal
        if not self.args.no_filter:
            if  it[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
                it[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
                it[nonrefbases[2] + 1][0] < self.args.tumour_variant:
                    return None,refbase
           
        ## get corresponding reference tuple
        rt = classifier.bam.get_reference_tuple(chromosome_id, position)
        feature_set = classifier.features_module.Features(it, rt)
        return feature_set,refbase

#============================================
# Check the versions of all the dependencies
#============================================
class verify_dependencies(unittest.TestCase):
        
    def test_verify_python(self):
        #since the code works on 2.7.x, so we need to check if the major version
        #and the minor version. Don't worry about micro version
        python_version = sys.version_info
        self.assertEqual(python_version[0], 2, 'the python version doesn\'t'
                        ' match the required version' )
        self.assertEqual(python_version[1], 7, 'The python version doesn\'t'
                         ' match the required version')
    
    def test_verify_numpy(self):
        numpy_version = numpy.version.full_version
        numpy_version = numpy_version.strip().split('.')
        self.assertEqual(numpy_version[0], str(1), 'the numpy version doesn\'t'
                         ' match the required version')
        self.assertEqual(numpy_version[1], str(7), 'The numpy version doesn\'t'
                        'match the required version')
        self.assertEqual(numpy_version[2], str(1), 'The numpy version doesn\'t'
                         ' match the required version')
        
    def test_verify_scipy(self):
        scipy_version = scipy.version.full_version
        scipy_version = scipy_version.strip().split('.')
        self.assertEqual(scipy_version[0], str(0), 'the scipy version doesn\'t'
                         ' match the required version')
        self.assertEqual(scipy_version[1], str(12),'the scipy version doesn\'t'
                         ' match the required version')
        self.assertEqual(scipy_version[2], str(0), 'the scipy version doesn\'t'
                         ' match the required version')

    def test_verify_sklearn(self):
        sklearn_version = sklearn.__version__
        sklearn_version = sklearn_version.strip().split('.')
        self.assertEqual(sklearn_version[0], str(0), 'the SciKit Learn version'
                         ' doesn\'t match the required version')
        self.assertEqual(sklearn_version[1], str(14), 'the SciKit Learn version'
                         ' doesn\'t match the required version')
        self.assertEqual(sklearn_version[2], str(1), 'the SciKit Learn version'
                         ' doesn\'t match the required version')
    
    def test_verify_matplotlib(self):
        matplotlib_version = matplotlib.__version__
        matplotlib_version = matplotlib_version.strip().split('.')
        self.assertEqual(matplotlib_version[0], str(1), 'the matplotlib version'
                         ' doesn\'t match the required version')
        self.assertEqual(matplotlib_version[1], str(2), 'the matplotlib version'
                         ' doesn\'t match the required version')
        self.assertEqual(matplotlib_version[2], str(1), 'the matplotlib version'
                         ' doesn\'t match the required version')
    
    def test_verify_imports(self):
        self.assertNotEqual(logging, None, 'logging couldn\'t be imported')
        self.assertNotEqual(pybamapi, None, 'pybamapi couldn\'t be imported')
        self.assertNotEqual(os, None, 'os couldn\'t be imported')
        self.assertNotEqual(resource, None, 'resource couldn\'t be imported')
        self.assertNotEqual(re, None, 're couldn\'t be imported')
        self.assertNotEqual(log10, None, 'log10 couldn\'t be imported')
        self.assertNotEqual(RandomForestClassifier, None, 'RandomForestClassifier'
                            'couldn\'t be imported')
        self.assertNotEqual(features, None, 'features couldn\'t be imported')
        self.assertNotEqual(features_single, None, 'features_single couldn\'t'
                            ' be imported')
        self.assertNotEqual(features_deep, None, 'features_deep couldn\'t '
                            'be imported')
        self.assertNotEqual(features_deep_single, None, 'features_deep_single'
                            ' couldn\'t be imported')
        self.assertNotEqual(plt, None, 'plt couldn\'t be imported')
        self.assertNotEqual(cross_validation, None, 'cross_validation '
                            'couldn\'t be imported')
        self.assertNotEqual(roc_curve, None, 'roc_curve couldn\'t be '
                            'imported')
        self.assertNotEqual(auc, None, 'auc couldn\'t be imported')
        self.assertNotEqual(Template, None, 'Template couldn\'t be imported')
        self.assertNotEqual(datetime, None, 'datetime couldn\'t be imported')
        self.assertNotEqual(defaultdict, None, 'defaultdict couldn\'t be imported')        
    
#===================================================
#Test for error in reading single alignment(MUT-140)
#===================================================
class single_position_error(unittest.TestCase,base_class):

    def setUp(self):
        self.args.interval = None
        self.args.positions_file = "unit_test/positions_missed.pos"
        self.classifier = bamutils.Classifier(self.args)
        self.target_positions = self.classifier.get_positions()
    
    def test_error_reading_single_position(self):
        #get the alignment for the single alignment(n) and 
        #then for the range(n-1,n), compare
        temp_target_positions=[]
        for target_position in self.target_positions:
            if target_position[1]==None and target_position[2]==None:
                break
            if target_position[1] == target_position[2]:
                temp_target_positions.append(target_position)
                tuples = self.classifier.bam.get_tuples(temp_target_positions)
                if self.args.single:
                    for it in tuples:
                        temp_target_positions = []
                        target_position[1] = target_position[1] -1
                        temp_target_positions.append(target_position)
                        newtuples = self.classifier.bam.get_tuples(temp_target_positions)
                        for new_it in newtuples:
                            if it[0] == new_it[0]:
                                self.assertEqual(new_it, it, 'Error reading this'
                                             ' sequence from the bam file('+str(it[0])+')'
                                             ', please use range instead')
                else:
                    for tt,nt in tuples:
                        temp_target_positions = []
                        target_position[1] = target_position[1] -1
                        temp_target_positions.append(target_position)
                        newtuples = self.classifier.bam.get_tuples(temp_target_positions)
                        for new_tt,new_nt in newtuples:
                            if nt[0] == new_nt[0] and tt[0] == new_tt[0]:
                                self.assertEqual(new_tt, tt, 'Error reading this'
                                             ' sequence from the bam file('+str(tt[0])+')'
                                             ', please use range instead')
                                self.assertEqual(new_nt, nt, 'Error reading this'
                                             ' sequence from the bam file('+str(tt[0])+')'
                                             ', please use range instead')
     
#==========================================================
# Check the validity of the sequences returned. 
#==========================================================
class verify_tuples_positions(unittest.TestCase,base_class):

    def test_tuples(self):
        self.args.interval = None
        self.args.positions_file = "unit_test/verify_tuples_positions"
        tuples,classifier = self.get_tuples(self.args)
        
        if self.args.single:
            #create dict from known tuples(from IGV) 
            tuples_dict_single = self.create_dict_tuples_single(
                            'unit_test/verify_tuples_outputs_single')
            #Iterate over tuples and compare with the dict
            #MUT-1: Check validity of sequences
            #MUT-8: Counts will be verified
            for it in tuples:
                chromosome_id = it[-1]
                position = it[0]
                rt = classifier.bam.get_reference_tuple(chromosome_id, position)
                resulting_tuple = [str(it),str(rt)]
                input_tuple=[]
                if tuples_dict_single.has_key(str(position)):
                    input_tuple=tuples_dict_single.get(str(position))
                    self.assertEqual(resulting_tuple, input_tuple,
                                      'The tuple read by Museq is incorrect')
 
        else:
            tuples_dict_paired = self.create_dict_tuples_paired('unit_test/verify_tuples_outputs_paired')
            for tt,nt in tuples:
                chromosome_id = tt[-1]
                position = tt[0]
                rt = classifier.bam.get_reference_tuple(chromosome_id, position)
                resulting_tuple = [str(tt),str(nt),str(rt)]
                input_tuple=[]
                if tuples_dict_paired.has_key(str(position)):
                    input_tuple=tuples_dict_paired.get(str(position))
                    self.assertEqual(resulting_tuple, input_tuple, 
                                     'The tuple read by Museq is incorrect')

    def test_tuples_returned(self):
        self.args.interval = None
        self.args.positions_file = None
        file_stream = open("unit_test/test_tuples_range")
        for line in file_stream:
            self.args.interval = line
            tuples,_ = self.get_tuples(self.args)
            interval_list = self.args.interval.replace(':','-').split('-')
            #Ensure the tuples are from the specified range
            #MUT-9: Ensures the tuples are from specified position only
            #MUT-13: Infinite loop would fail since the read from first position is out of range 
            if self.args.single:
                for it in tuples:
                    if len(interval_list)==3:
                        self.assertGreaterEqual(it[0], int(interval_list[-2]),
                                        'The tuples('+str(it[0])+') read by museq aren\'t in range')
                        self.assertLessEqual(it[0],int(interval_list[-1])+1,
                                    'The tuples('+str(it[0])+') read by museq aren\'t in range' )
                    elif len(interval_list)==2:
                        self.assertEqual(it[0],int(interval_list[-1]),
                                     'The tuples('+str(it[0])+') read by museq aren\'t in range' )
            else: 
                for tt,_ in tuples:
                    if len(interval_list)==3:
                        self.assertGreaterEqual(tt[0], int(interval_list[-2]),
                                        'The tuples read('+str(tt[0])+') by museq aren\'t in range')
                        self.assertLessEqual(tt[0],int(interval_list[-1])+1,
                                    'The tuples('+str(tt[0])+') read by museq aren\'t in range' )
                    elif len(interval_list)==2:
                        self.assertEqual(tt[0],int(interval_list[-1]),
                                     'The tuples('+str(tt[0])+') read by museq aren\'t in range' )

    def test_circular_reads(self):
        self.args.interval = None
        self.args.positions_file = None
        #MUT-10: ensure the tuples are always in increasing order and
        #there are no circular reads or reads out of order
        file_stream = open("unit_test/test_tuples_range")
        for line in file_stream:
            self.args.interval = line
            tuples,_ = self.get_tuples(self.args)
            prev_pos = 0
            if self.args.single:
                for it in tuples:
                    current_pos = it[0]
                    self.assertGreater(current_pos, prev_pos,
                                    'The positions should increase at every iteration')
                    prev_pos = current_pos
            else:
                for tt,_ in tuples:
                    current_pos = tt[0]
                    self.assertGreater(current_pos, prev_pos,
                                    'The positions should increase at every iteration')
                    prev_pos = current_pos

#check invalid positions and ensure no tuples are returned
class verify_boundary_invalid_reads(unittest.TestCase,base_class):

    def test_empty_positions(self):
        self.args.interval = None
        self.args.positions_file = "unit_test/emptyreads_positions"
        self.args.out = "unit_test/empty_reads.vcf"
        
        tuples,_ = self.get_tuples(self.args)
        if self.args.single:
            for it in tuples:
                #if generator is not empty throw an error
                self.assertEqual(True, False,
                                 'These positions('+str(it[0])+') shouldn\'t return tuples ')
        else:
            for tt,_ in tuples:
                self.assertEqual(True, False, 
                                 'These positions('+str(tt[0])+') shouldn\'t return tuples ')
        try:
            os.remove(self.args.out)
        except:
            logging.error('couldn\'t delete: '+self.args.out)


        
    def test_invalid_reads(self):
        #Ensure that positions beyond the chromosome's end and invalid positions don't return tuples
        self.args.interval = None
        self.args.positions_file = "unit_test/invalidreads_positions"
        self.args.out = "unit_test/invalid_reads.vcf" 
        tuples,_ = self.get_tuples(self.args)
        if self.args.single:
            try:
                for it in tuples:
                    #Should not start reading from beginning of the chromosome(MUT-13)
                    self.assertEqual(True, False,'These positions('+it[0]+')'
                                     ' shouldn\'t return tuples ')
            except Exception as ex:
                self.assertEqual(type(ex).__name__,'RuntimeError',
                                  'These positions should return Runtime Error only')
        else:
            try:
                for tt,_ in tuples:
                    self.assertEqual(True, False,'These positions('+tt[0]+')'
                                     ' shouldn\'t return tuples')
            except Exception as ex:
                self.assertEqual(type(ex).__name__,'RuntimeError',
                                 'These positions should return Runtime Error only')
        try:
            os.remove(self.args.out)
        except:
            logging.error('couldn\'t delete: '+self.args.out)
                
    def test_boundary_positions(self):
        #MUT-2 : Checks the boundary positions(loaded from file)
        self.args.interval = None
        self.args.positions_file = "unit_test/verify_boundary_positions"
        tuples,classifier = self.get_tuples(self.args)
        
        if self.args.single:
            #create dict from known tuples(from IGV) 
            tuples_dict_single = self.create_dict_tuples_single(
                                'unit_test/verify_boundary_outputs_single')
            #Iterate over tuples and compare
            for it in tuples:
                chromosome_id = it[-1]
                position = it[0]
                rt = classifier.bam.get_reference_tuple(chromosome_id, position)
                resulting_tuple = [str(it),str(rt)]
                input_tuple=[]
                if tuples_dict_single.has_key(str(position)):
                    input_tuple=tuples_dict_single.get(str(position))
                    self.assertEqual(resulting_tuple, input_tuple, 
                                     'The tuple read by Museq is incorrect')
 
        else:
            tuples_dict_paired = self.create_dict_tuples_paired(
                                'unit_test/verify_boundary_outputs_paired')
            for tt,nt in tuples:
                chromosome_id = tt[-1]
                position = tt[0]
                rt = classifier.bam.get_reference_tuple(chromosome_id, position)
                resulting_tuple = [str(tt),str(nt),str(rt)]
                input_tuple=[]
                if tuples_dict_paired.has_key(str(position)):
                    input_tuple=tuples_dict_paired.get(str(position))
                    self.assertEqual(resulting_tuple, input_tuple,
                                      'The tuple read by Museq is incorrect')
        
        
#=====================================================
#Check the output against the known results
#=====================================================
class verify_outputs(unittest.TestCase,base_class):
    
    def setUp(self):
        self.args.positions_file= "unit_test/correct_positions"
        self.args.out = "unit_test/correct_results.vcf"
              
    def test_verify_with_labels(self):
        #compares the results from museq against the known output
        #Checks if the probabilities match.
        self.run_classifier(self.args)
        if self.args.single:
            label_dict = self.create_dict_vcf("unit_test/correct_labels_single")
            result_dict = self.create_dict_vcf("unit_test/correct_results.vcf")
        
            for key,value in result_dict.iteritems():
                generated_result = value
                if label_dict.has_key(key):
                    expected_result = label_dict.get(key)
                    self.compare_vcf_dicts(generated_result, expected_result)
        
        else:
            label_dict = self.create_dict_vcf("unit_test/correct_labels_paired")
            result_dict = self.create_dict_vcf("unit_test/correct_results.vcf")
        
            for key,value in result_dict.iteritems():
                generated_result = value
                if label_dict.has_key(key):
                    expected_result = label_dict.get(key)
                    self.compare_vcf_dicts(generated_result, expected_result)
    
    def tearDown(self):
        try:
            os.remove(self.args.out)
        except:
            logging.error('couldn\'t delete: '+self.args.out)

               
#==================================================
#Ensure memory used stays under tolerance for -b
#==================================================
class check_memory_footprint(unittest.TestCase,base_class):
    
    def setUp(self):
        self.args.buffer_size = "500M"
        self.args.positions_file = None
        self.args.interval = "1:1-200000"
        
        #not correct: run it in a separate thread
        #python threads are lightweight and slow(very slow when fitting)
        #This test should be the first to run
    def test_memory_footprint(self):
        self.run_classifier(self.args)
        
        memory_used = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        tolerance = 0
        s = re.split('(\d+)', self.args.buffer_size)
        if s[-1] == 'M':
            tolerance = int(s[1])*1024*2
        elif s[-1] == 'G':
            tolerance = int(s[1])*1024*1024*1.5
        self.assertLess(memory_used, tolerance, 'Using more memory'
                        ' than the tolerance limit')

#=================================
#check for  Reference file
#=================================

class verify_reference_file(unittest.TestCase,base_class):
    
    def setUp(self):
        self.args.positions_file = None
        self.args.interval = None
           
    def test_reference_file(self):
        #Ensure that refbase is always less than 5
        #MUT-20 : Test fails if invalid pos is present in ref file
        tuples,classifier = self.get_tuples(self.args)
        if self.args.single:
            for it in tuples:
                chromosome_id = it[-1]
                position = it[0]
                refbase = classifier.bam.get_reference_base(chromosome_id,
                                                             position, index=True)
                chromosome_name = classifier.bam.get_chromosome_name(chromosome_id)
                #Museq bails if it gets anything other than 0-3
                self.assertLess(refbase, 5, 'The reference base should be'
                                ' A,C,G,T or N at '+""+str(chromosome_name)+
                                ":"+str(position)+"\n")
        else:
            for tt,_ in tuples: 
                chromosome_id = tt[-1]
                position = tt[0]
                refbase = classifier.bam.get_reference_base(chromosome_id,
                                                             position, index=True)
                chromosome_name = classifier.bam.get_chromosome_name(chromosome_id)
                #Museq bails if it gets anything other than 0-3
                self.assertLess(refbase, 5, 'The reference base should be A,C,G,T or N'
                                ' at '+""+str(chromosome_name)+":"+str(position)+"\n")

    #made changes in bamutils for this code. line:222,321,289
    def test_invalid_refbase(self):
        #MUT-20:Checks if invalid refbase raises any errors
        #try inserting invalid values for refbase and see if museq throws error
        #trinucleotide just overwrites
        #the filter doesnt work for those wrong positions.
        #could compare outputs instead too
        #paired mode fails in the filter,single doesnt
        self.args.invalid = True
        self.args.out = "unit_test/output_ref"
        try:
            self.run_classifier(self.args)
            self.assertEqual(True,False,'MuSeq should fail when invalid'
                             ' positions are inserted in refbase')
        except TypeError:
            pass
#===========================================
#Ensure all features are setup correctly
#===========================================
class verify_features(unittest.TestCase,base_class):
    
    def setUp(self):
        self.args.positions_file = None
        self.args.interval = "1:1000-5000"
        self.tuples,self.classifier = self.get_tuples(self.args)
        
    def __check_filter_single(self,feature,refbase,it):
        nonrefbases = []
        for i in range(4):
            if i == refbase: 
                continue
            nonrefbases.append(i)
        #check if the feature
        if it[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
        it[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
        it[nonrefbases[2] + 1][0] < self.args.tumour_variant:
            self.assertEqual(feature,None,'This feature should be filtered')
        else:
            self.assertNotEqual(feature,None,'This feature should'
                                    ' not be filtered')
    
    def __check_features_single(self,features):
        self.assertEqual(len(features.coverage_features),3, 
                                 'Incorrect coverage_features')
        self.assertEqual(len(features.feature_set),23,
                                  'Incorrect features')
        #contamination is less than 1
        self.assertLess(features.coverage_features[1][1], 1.0,
                                 'invalid coverage feature')
        self.assertEqual(features.coverage_features[2][1], 0.0,
                                  'invalid coverage feature') 
                        
    def __check_features_single_deep(self,features):
        pass
    
    def __check_filter_paired(self,feature,refbase,tt,nt):
        nonrefbases = []        
        if feature!=None:
            for i in range(4):
                if i == refbase:
                    continue
                nonrefbases.append(i)
        normal_val = (nt[5][0] - nt[refbase + 1][0]) / nt[5][0] > (self.args.normal_variant / 100)
        if not normal_val:
            if tt[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
                tt[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
                tt[nonrefbases[2] + 1][0] < self.args.tumour_variant:
                    self.assertEqual(feature,None,
                                    'This feature should be filtered')
            else:
                self.assertNotEqual(feature,None,
                                   'This feature should not be filtered')
        else:
            self.assertEqual(feature,None,
                          'This feature should be filtered')           
                 
    def __check_features_paired(self,feature):
        self.assertEqual(len(feature.coverage_features),4,
                                      'Incorrect coverage_features')
        self.assertEqual(len(feature.feature_set),60,
                                      'Incorrect features')
                        
        #contamination is less than 1
        self.assertLess(feature.coverage_features[2][1], 1.0,
                                     'invalid coverage feature')
        self.assertEqual(feature.coverage_features[3][1], 0.0,
                                      'invalid coverage feature')
        #self.assertLessEqual(feature.feature_set[0][1],1.0,
        #                     'tumour_indels should be less than 1')
        #self.assertLessEqual(feature.feature_set[1][1],1.0,
        #                     'normal_indels should be less than 1')
        #self.assertLessEqual(feature.feature_set[2][1],1.0,
        #                     'tumour_ref_depth should be less than 1')
        self.assertLessEqual(feature.feature_set[3][1],1.0,
                             'normal_ref_depth should be less than 1')
        #Add more tests?? according to the range of the feature values
            
    
    
    
    def __check_features_paired_deep(self,feature):
        self.assertEqual(len(feature.coverage_features),2,
                                      'Incorrect coverage_features')
        self.assertEqual(len(feature.feature_set),41,
                                      'Incorrect features')
    
                
    def test_tumour_normal_filtering(self):
        self.assertLess(self.args.tumour_variant, 10,
                         'Tumour Variant set to a very high value')
        self.assertGreater(self.args.normal_variant, 10,
                            'Normal Variant set to a very low value')
        
        if self.args.single:
            for it in self.tuples:
                feature,refbase = self.get_features_single(it,self.classifier)
                if feature !=None:
                    self.assertEqual(it, feature.it, 
                                     'the position in tuple and feature are different')
                    self.__check_filter_single(feature, refbase, it)
                    
                    if self.args.deep:
                        self.__check_features_single_deep(feature)
                    else:
                        self.__check_features_single(feature)  
                       
        else:
            for tt,nt in self.tuples:
                feature,refbase = self.get_features_paired(tt, nt, self.classifier)
                if feature!=None:
                    self.assertEqual(tt, feature.tt, 
                                     'the position in tuple and feature are different')
                    self.assertEqual(nt, feature.nt,
                                     'the position in tuple and feature are different')
                    self.__check_filter_paired(feature, refbase, tt, nt)
                    
                    if self.args.deep:
                        self.__check_features_paired_deep(feature)
                    else:
                        self.__check_features_paired(feature)  
       

#====================
#check rmdups flag
#====================
class verify_flags(unittest.TestCase,base_class):

    def setUp(self):
        self.args.positions_file = None
        self.args.interval = None
        self.classifier = bamutils.Classifier(self.args)
    
    def test_rmdups(self):
        rmdups = self.classifier.rmdups
        if self.args.deep:
            self.assertEqual(rmdups,False,'rmdups not set for deepseq analysis')
        else:
            self.assertEqual(rmdups,True,'rmdups set for non deepseq analysis')

class verify_individual_functions(unittest.TestCase,base_class):
    def setUp(self):
        self.args.interval = "1"
        self.args.positions_file = None
        
    def test_chromosome_ids(self):
        tuples,classifier = self.get_tuples(self.args)
        if self.args.single:
            for it in tuples:
                chromosome_id = it[-1]
                self.assertLessEqual(chromosome_id,24,'Chromosome_id should be less than 24')
                chromosome_name = classifier.bam.get_chromosome_name(chromosome_id)
                self.assertRegexpMatches(chromosome_name, "[1-22XYMT]", "Error Chromosome name invalid")
        else:
            for tt,_ in tuples:
                chromosome_id = tt[-1]
                self.assertLessEqual(chromosome_id,24,'Chromosome_id should be less than 24')
                chromosome_name = classifier.bam.get_chromosome_name(chromosome_id)
                self.assertRegexpMatches(chromosome_name, "[1-22XY]|[MT]", "Error Chromosome name invalid")
    
    def test_trinucleotide_context(self):
        #only need to run it once.
        if self.args.single:
            classifier = bamutils.Classifier(self.args)
            for chromosome_id in xrange(25):
                for position in xrange(1000):
                    #Throws runtimeerror if unable to get base
                    try:
                        tc = classifier.bam.get_trinucleotide_context(chromosome_id, position)
                        self.assertRegexpMatches(tc, '[ACGTN]|[ACGTN]|[ACGTN]', 'Invalid Trinucleotide_context')
                    except RuntimeError:
                        pass
    def test_reference_base(self):
        if self.args.single:
            classifier = bamutils.Classifier(self.args)
            for chromosome_id in xrange(25):
                for position in xrange(100):
                    #Throws runtimeerror if unable to get base
                    try:
                        refbase = classifier.bam.get_reference_base(chromosome_id, position, index=True)
                        self.assertRegexpMatches(str(refbase), '[0-4]|[ACGTN]', 'Invalid Trinucleotide_context')
                    except RuntimeError:
                        pass
            for chromosome_id in xrange(25):
                for position in xrange(100):
                    #Throws runtimeerror if unable to get base
                    try:
                        refbase = classifier.bam.get_reference_base(chromosome_id, position, index=False)
                        self.assertRegexpMatches(refbase, '[0-4]|[ACGTN]', 'Invalid Trinucleotide_context')
                    except RuntimeError:
                        pass

#========================================
#Ensure chr only fails with runtime error
#========================================
class verify_position_with_chr(unittest.TestCase,base_class):
    
    def setUp(self):
        self.args.positions_file = "unit_test/position_chr"
        self.args.interval = None

    def test_position_with_chr(self):
        try:
            self.run_classifier(self.args)
            self.assertEqual(True, False, 'chr should fail in classifier')
        except Exception as ex:
            self.assertEqual(type(ex).__name__,'RuntimeError', 
                             'Chr not working properly')


def suite_all_tests():   
    suite = unittest.TestSuite()
    loader = unittest.TestLoader()
    
    memory_footprint = loader.loadTestsFromTestCase(check_memory_footprint)
    dependencies = loader.loadTestsFromTestCase(verify_dependencies)
    single_pos_error = loader.loadTestsFromTestCase(single_position_error)
    boundary_invalid_reads = loader.loadTestsFromTestCase(verify_boundary_invalid_reads)
    verify_outputs_test = loader.loadTestsFromTestCase(verify_outputs)
    #ref_file = loader.loadTestsFromTestCase(verify_reference_file)
    features_test = loader.loadTestsFromTestCase(verify_features)
    flags = loader.loadTestsFromTestCase(verify_flags)
    chr_pos = loader.loadTestsFromTestCase(verify_position_with_chr)
    tuples = loader.loadTestsFromTestCase(verify_tuples_positions)
    individual_function = loader.loadTestsFromTestCase(verify_individual_functions)
    
    suite.addTests(memory_footprint)
    suite.addTests(dependencies)
    suite.addTests(single_pos_error)
    suite.addTests(boundary_invalid_reads)
    suite.addTests(verify_outputs_test)
#    suite.addTests(ref_file)
    suite.addTests(features_test)
    suite.addTests(flags)
    suite.addTests(chr_pos)
    suite.addTests(tuples)
    suite.addTests(individual_function)

    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
