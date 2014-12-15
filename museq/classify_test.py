import unittest

import bamutils
import sys
import logging
import os
import resource
import re
from classify_test_api import base_class,initargs

mutationSeq_version="4.3.1"

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
        import numpy
        numpy_version = numpy.version.full_version
        numpy_version = numpy_version.strip().split('.')
        numpy_version = map(int,numpy_version)
        
        self.assertGreaterEqual(numpy_version[0], 1, 'the numpy version doesn\'t'
                         ' match the required version')
        
        if numpy_version == 1:
            self.assertGreaterEqual(numpy_version[1], 7, 'The numpy version doesn\'t'
                                    'match the required version')

        
    def test_verify_scipy(self):
        import scipy
        scipy_version = scipy.version.full_version
        scipy_version = scipy_version.strip().split('.')
        
        scipy_version = map(int,scipy_version)
        
        self.assertGreaterEqual(scipy_version[0], 0, 'the numpy version doesn\'t'
                         ' match the required version')
        
        if scipy_version == 0:
            self.assertGreaterEqual(scipy_version[1], 12, 'The numpy version doesn\'t'
                                    'match the required version')



    def test_verify_sklearn(self):
        import sklearn
        sklearn_version = sklearn.__version__
        sklearn_version = sklearn_version.strip().split('.')
        self.assertEqual(sklearn_version[0], str(0), 'the SciKit Learn version'
                         ' doesn\'t match the required version')
        self.assertEqual(sklearn_version[1], str(14), 'the SciKit Learn version'
                         ' doesn\'t match the required version')
        self.assertEqual(sklearn_version[2], str(1), 'the SciKit Learn version'
                         ' doesn\'t match the required version')
    
    def test_verify_matplotlib(self):
        import matplotlib
        matplotlib_version = matplotlib.__version__
        matplotlib_version = matplotlib_version.strip().split('.')
        matplotlib_version = map(int,matplotlib_version)
        
        self.assertGreaterEqual(matplotlib_version[0], 1, 'the numpy version doesn\'t'
                         ' match the required version')
        
        if matplotlib_version == 1:
            self.assertGreaterEqual(matplotlib_version[1], 2, 'The numpy version doesn\'t'
                                    'match the required version')
    
    def test_verify_imports(self):
        self.assertNotEqual(logging, None, 'logging couldn\'t be imported')
        
        import pybamapi
        self.assertNotEqual(pybamapi, None, 'pybamapi couldn\'t be imported')
        
        
        self.assertNotEqual(os, None, 'os couldn\'t be imported')
        self.assertNotEqual(resource, None, 'resource couldn\'t be imported')
        self.assertNotEqual(re, None, 're couldn\'t be imported')
        
        from math import log10
        self.assertNotEqual(log10, None, 'log10 couldn\'t be imported')

        from sklearn.ensemble import RandomForestClassifier
        self.assertNotEqual(RandomForestClassifier, None, 'RandomForestClassifier'
                            'couldn\'t be imported')
        
        import features, features_single, features_deep, features_deep_single
        self.assertNotEqual(features, None, 'features couldn\'t be imported')
        self.assertNotEqual(features_single, None, 'features_single couldn\'t'
                            ' be imported')
        self.assertNotEqual(features_deep, None, 'features_deep couldn\'t '
                            'be imported')
        self.assertNotEqual(features_deep_single, None, 'features_deep_single'
                            ' couldn\'t be imported')
        
        import matplotlib.pyplot as plt
        self.assertNotEqual(plt, None, 'plt couldn\'t be imported')
        
        from sklearn import cross_validation
        self.assertNotEqual(cross_validation, None, 'cross_validation '
                            'couldn\'t be imported')
        
        from sklearn.metrics import roc_curve, auc
        self.assertNotEqual(roc_curve, None, 'roc_curve couldn\'t be '
                            'imported')
        self.assertNotEqual(auc, None, 'auc couldn\'t be imported')
        
        from string import Template
        self.assertNotEqual(Template, None, 'Template couldn\'t be imported')
        
        from datetime import datetime
        self.assertNotEqual(datetime, None, 'datetime couldn\'t be imported')
        
        from collections import defaultdict
        self.assertNotEqual(defaultdict, None, 'defaultdict couldn\'t be imported')        
    
#===================================================
#Test for error in reading single position(MUT-140)
#===================================================
class single_position_error(unittest.TestCase,base_class):
    def setUp(self):
        self.args_paired = initargs()
        args_single = initargs()
        args_single.set_single()
        self.args_single = args_single
            
    def test_error_single_pos_single(self):
        self.args_single.interval = None
        self.args_single.positions_file = "unit_test/positions_missed.pos"
        self.classifier = bamutils.Classifier(self.args_single)
        self.target_positions = self.classifier.get_positions()
        temp_target_positions=[]
        for target_position in self.target_positions:
            if target_position[1]==None and target_position[2]==None:
                break
            if target_position[1] == target_position[2]:
                temp_target_positions.append(target_position)
                tuples = self.classifier.bam.get_tuples(temp_target_positions)
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

    
    def test_error_single_pos_paired(self):
        self.args_paired.interval = None
        self.args_paired.positions_file = "unit_test/positions_missed.pos"
        self.classifier = bamutils.Classifier(self.args_paired)
        self.target_positions = self.classifier.get_positions()
        #get the alignment for the single alignment(n) and 
        #then for the range(n-1,n), compare
        temp_target_positions=[]
        for target_position in self.target_positions:
            if target_position[1]==None and target_position[2]==None:
                break
            if target_position[1] == target_position[2]:
                temp_target_positions.append(target_position)
                tuples = self.classifier.bam.get_tuples(temp_target_positions)
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
    def setUp(self):
        self.args_paired = initargs()
        args_single = initargs()
        args_single.set_single()
        self.args_single = args_single
        
    def test_tuples_single(self):
        self.args_single.interval = None
        self.args_single.positions_file = "unit_test/verify_tuples_positions"
        tuples,classifier = self.get_tuples(self.args_single)
        if self.args_single.single:
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
 
    def test_tuples_paired(self):
        self.args_paired.interval = None
        self.args_paired.positions_file = "unit_test/verify_tuples_positions"
        tuples,classifier = self.get_tuples(self.args_paired)
        if not self.args_paired.single:
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
            

    def test_tuples_returned_single(self):
        self.args_single.interval = None
        self.args_single.positions_file = None
        file_stream = open("unit_test/test_tuples_range")
        for line in file_stream:
            self.args_single.interval = line
            tuples,_ = self.get_tuples(self.args_single)
            interval_list = self.args_single.interval.replace(':','-').split('-')
            #Ensure the tuples are from the specified range
            #MUT-9: Ensures the tuples are from specified position only
            #MUT-13: Infinite loop would fail since the read from first position is out of range 
            if self.args_single.single:
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
                self.assertEqual(True, False, 'Single mode flag not set')
                

    def test_tuples_returned_paired(self):
        self.args_paired.interval = None
        self.args_paired.positions_file = None
        file_stream = open("unit_test/test_tuples_range")
        for line in file_stream:
            self.args_paired.interval = line
            tuples,_ = self.get_tuples(self.args_paired)
            interval_list = self.args_paired.interval.replace(':','-').split('-')
            if not self.args_paired.single:
                for tt,_ in tuples:
                    if len(interval_list)==3:
                        self.assertGreaterEqual(tt[0], int(interval_list[-2]),
                                        'The tuples read('+str(tt[0])+') by museq aren\'t in range')
                        self.assertLessEqual(tt[0],int(interval_list[-1])+1,
                                    'The tuples('+str(tt[0])+') read by museq aren\'t in range' )
                    elif len(interval_list)==2:
                        self.assertEqual(tt[0],int(interval_list[-1]),
                                     'The tuples('+str(tt[0])+') read by museq aren\'t in range' )  
            else:
                self.assertEqual(True, False, 'Single mode flag set')
        
    def test_circular_reads_single(self):
        self.args_single.interval = None
        self.args_single.positions_file = None
        #MUT-10: ensure the tuples are always in increasing order and
        #there are no circular reads or reads out of order
        file_stream = open("unit_test/test_tuples_range")
        for line in file_stream:
            self.args_single.interval = line
            tuples,_ = self.get_tuples(self.args_single)
            prev_pos = 0
            if self.args_single.single:
                for it in tuples:
                    current_pos = it[0]
                    self.assertGreater(current_pos, prev_pos,
                                    'The positions should increase at every iteration')
                    prev_pos = current_pos
            else:
                self.assertEqual(True, False, 'Single Flag not set')
                    
    def test_circular_reads_paired(self):
        self.args_paired.interval = None
        self.args_paired.positions_file = None
        file_stream = open("unit_test/test_tuples_range")
        for line in file_stream:
            self.args_paired.interval = line
            tuples,_ = self.get_tuples(self.args_paired)
            prev_pos = 0
            if not self.args_paired.single:
                for tt,_ in tuples:
                    current_pos = tt[0]
                    self.assertGreater(current_pos, prev_pos,
                                    'The positions should increase at every iteration')
                    prev_pos = current_pos
            else:
                self.assertEqual(True, False, 'Single flag set')
                    

#check invalid positions and ensure no tuples are returned
class verify_boundary_invalid_reads(unittest.TestCase,base_class):
    def setUp(self):
        self.args_paired = initargs()
        args_single = initargs()
        args_single.set_single()
        self.args_single = args_single
        self.delete_files = []

    def test_empty_positions_single(self):
        self.args_single.interval = None
        self.args_single.positions_file = "unit_test/emptyreads_positions"
        self.args_single.out = "unit_test/empty_reads.vcf"
        
        tuples,_ = self.get_tuples(self.args_single)
        if self.args_single.single:
            for it in tuples:
                #if generator is not empty throw an error
                self.assertEqual(True, False,
                                 'These positions('+str(it[0])+') shouldn\'t return tuples ')
        else:
            self.assertEqual(True, False, 'Single flag not set')


    def test_empty_positions_paired(self):
        self.args_paired.interval = None
        self.args_paired.positions_file = "unit_test/emptyreads_positions"
        self.args_paired.out = "unit_test/empty_reads.vcf"
        
        tuples,_ = self.get_tuples(self.args_paired)
        if not self.args_paired.single:
            for tt,_ in tuples:
                self.assertEqual(True, False, 
                                 'These positions('+str(tt[0])+') shouldn\'t return tuples ')
        else:
            self.assertEqual(True, False, 'Single flag set')
        self.delete_files.append(self.args_paired.out)
        
    def test_invalid_reads_single(self):
        #Ensure that positions beyond the chromosome's end and invalid positions don't return tuples
        self.args_single.interval = None
        self.args_single.positions_file = "unit_test/invalidreads_positions"
        self.args_single.out = "unit_test/invalid_reads.vcf" 
        tuples,_ = self.get_tuples(self.args_single)
        if self.args_single.single:
            try:
                for it in tuples:
                    #Should not start reading from beginning of the chromosome(MUT-13)
                    self.assertEqual(True, False,'These positions('+it[0]+')'
                                     ' shouldn\'t return tuples ')
            except Exception as ex:
                self.assertEqual(type(ex).__name__,'RuntimeError',
                                  'These positions should return Runtime Error only')
        else:
            self.assertEqual(True, False, 'Single flag not set')

            
    def test_invalid_reads_paired(self):
        #Ensure that positions beyond the chromosome's end and invalid positions don't return tuples
        self.args_paired.interval = None
        self.args_paired.positions_file = "unit_test/invalidreads_positions"
        self.args_paired.out = "unit_test/invalid_reads.vcf" 
        tuples,_ = self.get_tuples(self.args_paired)
        if  self.args_paired.single:
            self.assertEqual(True, False, 'Single flag set')
                   
        else:
            try:
                for tt,_ in tuples:
                    self.assertEqual(True, False,'These positions('+tt[0]+')'
                                     ' shouldn\'t return tuples')
            except Exception as ex:
                self.assertEqual(type(ex).__name__,'RuntimeError',
                                 'These positions should return Runtime Error only')
        self.delete_files.append(self.args_paired.out)
                
    def test_boundary_positions_single(self):
        #MUT-2 : Checks the boundary positions(loaded from file)
        self.args_single.interval = None
        self.args_single.positions_file = "unit_test/verify_boundary_positions"
        tuples,classifier = self.get_tuples(self.args_single)
        
        if self.args_single.single:
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
            self.assertEqual(True, False, 'Single flag not set')
           
           
    def test_boundary_positions_paired(self):
        #MUT-2 : Checks the boundary positions(loaded from file)
        self.args_paired.interval = None
        self.args_paired.positions_file = "unit_test/verify_boundary_positions"
        tuples,classifier = self.get_tuples(self.args_paired)
        
        if self.args_paired.single:
            self.assertEqual(True, False, 'single flag set')
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
    def tearDown(self):
        for files in self.delete_files:
            try:
                os.remove(files)
            except:
                logging.error('couldn\'t remove file:'+files)
        
#=====================================================
#Check the output against the known results
#=====================================================
class verify_outputs(unittest.TestCase,base_class):
    
    def setUp(self):
        self.args_paired = initargs()
        args_single = initargs()
        args_single.set_single()
        self.args_single = args_single
              
    def test_verify_with_labels_single(self):
        #compares the results from museq against the known output
        #Checks if the probabilities match.
        self.args_single.positions_file= "unit_test/correct_positions"
        self.args_single.out = "unit_test/correct_results.vcf"
        self.run_classifier(self.args_single)
        if self.args_single.single:
            label_dict = self.create_dict_vcf("unit_test/correct_labels_single")
            result_dict = self.create_dict_vcf("unit_test/correct_results.vcf")
        
            for key,value in result_dict.iteritems():
                generated_result = value
                if label_dict.has_key(key):
                    expected_result = label_dict.get(key)
                    self.compare_vcf_dicts(generated_result, expected_result)
        
        else:
            self.assertEqual(True, False, 'single flag not set')
            
    def test_verify_with_labels(self):
        #compares the results from museq against the known output
        #Checks if the probabilities match.
        self.args_paired.positions_file= "unit_test/correct_positions"
        self.args_paired.out = "unit_test/correct_results.vcf"
        self.run_classifier(self.args_paired)
        if self.args_paired.single:
            self.assertEqual(True,False,'Single flag set')
            
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
            os.remove(self.args_single.out)
        except:
            logging.error('couldn\'t delete: '+self.args.out)

               
#==================================================
#Ensure memory used stays under tolerance for -b
#==================================================
class check_memory_footprint(unittest.TestCase,base_class):
    
    def setUp(self):
        self.args_paired = initargs()
        args_single = initargs()
        args_single.set_single()
        self.args_single = args_single
        
        #not correct: run it in a separate thread
        #python threads are lightweight and slow(very slow when fitting)
        #This test should be the first to run
    def test_memory_footprint_paired(self):
        self.args_paired.buffer_size = "500M"
        self.args_paired.positions_file = None
        self.args_paired.interval = "1:1-200000"

        self.run_classifier(self.args_paired)
        
        memory_used = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        tolerance = 0
        s = re.split('(\d+)', self.args_paired.buffer_size)
        if s[-1] == 'M':
            tolerance = int(s[1])*1024*2
        elif s[-1] == 'G':
            tolerance = int(s[1])*1024*1024*1.5
        self.assertLess(memory_used, tolerance, 'Using more memory'
                        ' than the tolerance limit')
    
    def test_memory_footprint_single(self):
        self.args_single.buffer_size = "500M"
        self.args_single.positions_file = None
        self.args_single.interval = "1:1-200000"

        self.run_classifier(self.args_single)
        
        memory_used = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        tolerance = 0
        s = re.split('(\d+)', self.args_paired.buffer_size)
        if s[-1] == 'M':
            tolerance = int(s[1])*1024*2
        elif s[-1] == 'G':
            tolerance = int(s[1])*1024*1024*1.5
        self.assertLess(memory_used, tolerance, 'Using more memory'
                        ' than the tolerance limit')

#=================================
#check for  Reference file
#=================================


#===========================================
#Ensure all features are setup correctly
#===========================================
class verify_features(unittest.TestCase,base_class):
    
    def setUp(self):
        self.args_paired = initargs()
        args_single = initargs()
        args_single.set_single()
        self.args_single = args_single
        
        self.args_single.positions_file = None
        self.args_single.interval = "1:1000-5000"
        self.tuples_single,self.classifier_single = self.get_tuples(self.args_single)
        
        self.args_paired.positions_file = None
        self.args_paired.interval = "1:1000-5000"
        self.tuples_paired,self.classifier_paired = self.get_tuples(self.args_paired)
        
    def __check_filter_single(self,feature,refbase,it):
        nonrefbases = []
        for i in range(4):
            if i == refbase: 
                continue
            nonrefbases.append(i)
        #check if the feature
        if it[nonrefbases[0] + 1][0] < self.args_single.tumour_variant and \
        it[nonrefbases[1] + 1][0] < self.args_single.tumour_variant and \
        it[nonrefbases[2] + 1][0] < self.args_single.tumour_variant:
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
        normal_val = (nt[5][0] - nt[refbase + 1][0]) / nt[5][0] > (self.args_paired.normal_variant / 100)
        if not normal_val:
            if tt[nonrefbases[0] + 1][0] < self.args_paired.tumour_variant and \
                tt[nonrefbases[1] + 1][0] < self.args_paired.tumour_variant and \
                tt[nonrefbases[2] + 1][0] < self.args_paired.tumour_variant:
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
        self.assertLess(self.args_single.tumour_variant, 10,
                         'Tumour Variant set to a very high value')
        self.assertGreater(self.args_single.normal_variant, 10,
                            'Normal Variant set to a very low value')

    def test_features_single(self):        
        if self.args_single.single:
            for it in self.tuples_single:
                feature,refbase = self.get_features_single(it,self.classifier_single,self.args_single)
                if feature !=None:
                    self.assertEqual(it, feature.it, 
                                     'the position in tuple and feature are different')
                    self.__check_filter_single(feature, refbase, it)
        else:
            self.assertEqual(True, False, 'Single Flag not set')

    def test_features_single_deep(self):     
        self.assertEqual(True, True, 'nothing for now')
#         if self.args.deep:
#             self.__check_features_single_deep(feature)
#         else:
#             self.__check_features_single(feature)  
                       
    def test_features_paired(self):
        if self.args_paired.single:
            self.assertEqual(True, False, 'Single Flag set')
        else:
            for tt,nt in self.tuples_paired:
                feature,refbase = self.get_features_paired(tt, nt, self.classifier_paired,self.args_paired)
                if feature!=None:
                    self.assertEqual(tt, feature.tt, 
                                     'the position in tuple and feature are different')
                    self.assertEqual(nt, feature.nt,
                                     'the position in tuple and feature are different')
                    self.__check_filter_paired(feature, refbase, tt, nt)

    def test_features_paired_deep(self):     
        self.assertEqual(True, True, 'nothing for now')           
#         if self.args.deep:
#             self.__check_features_paired_deep(feature)
#         else:
#             self.__check_features_paired(feature)  
       

#====================
#check rmdups flag
#====================
class verify_flags(unittest.TestCase,base_class):

    def setUp(self):
        self.args = initargs()
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
        self.args_paired = initargs()
        args_single = initargs()
        args_single.set_single()
        self.args_single = args_single

        
    def test_chromosome_ids_single(self):
        self.args_single.interval = "1"
        self.args_single.positions_file = None
        tuples,classifier = self.get_tuples(self.args_single)
        if self.args_single.single:
            for it in tuples:
                chromosome_id = it[-1]
                self.assertLessEqual(chromosome_id,24,'Chromosome_id should be less than 24')
                chromosome_name = classifier.bam.get_chromosome_name(chromosome_id)
                self.assertRegexpMatches(chromosome_name, "[1-22XYMT]", "Error Chromosome name invalid")
        else:
            self.assertEqual(True, False, 'Single flag not set')        
    
    def test_chromosome_ids_paired(self):
        self.args_paired.interval = "1"
        self.args_paired.positions_file = None
        tuples,classifier = self.get_tuples(self.args_paired)
        if self.args_paired.single:
            self.assertEqual(True, False, 'Single Flag set')
        else:
            for tt,_ in tuples:
                chromosome_id = tt[-1]
                self.assertLessEqual(chromosome_id,24,'Chromosome_id should be less than 24')
                chromosome_name = classifier.bam.get_chromosome_name(chromosome_id)
                self.assertRegexpMatches(chromosome_name, "[1-22XY]|[MT]", "Error Chromosome name invalid")
    
    def test_trinucleotide_context_single(self):
        #only need to run it once.
        if self.args_single.single:
            classifier = bamutils.Classifier(self.args_single)
            for chromosome_id in xrange(25):
                for position in xrange(1000):
                    #Throws runtimeerror if unable to get base
                    try:
                        tc = classifier.bam.get_trinucleotide_context(chromosome_id, position)
                        self.assertRegexpMatches(tc, '[ACGTN]|[ACGTN]|[ACGTN]', 'Invalid Trinucleotide_context')
                    except RuntimeError:
                        pass
        else:
            self.assertEqual(True, False, 'Single flag not set')
            
    def test_trinucleotide_context_paired(self):
        #only need to run it once.
        if self.args_paired.single:
            self.assertEqual(True, False, 'Single flag set')
        else:
            classifier = bamutils.Classifier(self.args_paired)
            for chromosome_id in xrange(25):
                for position in xrange(1000):
                    #Throws runtimeerror if unable to get base
                    try:
                        tc = classifier.bam.get_trinucleotide_context(chromosome_id, position)
                        self.assertRegexpMatches(tc, '[ACGTN]|[ACGTN]|[ACGTN]', 'Invalid Trinucleotide_context')
                    except RuntimeError:
                        pass

                    
    def test_reference_base_single(self):
        if self.args_single.single:
            classifier = bamutils.Classifier(self.args_single)
            for chromosome_id in xrange(25):
                for position in xrange(100):
                    #Throws runtimeerror if unable to get base
                    try:
                        refbase = classifier.bam.get_reference_base(chromosome_id, position, index=True)
                        self.assertRegexpMatches(str(refbase), '[0-4]|[ACGTN]', 'Invalid Trinucleotide_context')
                    except RuntimeError:
                        pass

    def test_reference_base_paired(self):
        if self.args_paired.single:
            self.assertEqual(True, False, 'Single flag set')
        else:
            classifier = bamutils.Classifier(self.args_paired)
            for chromosome_id in xrange(25):
                for position in xrange(100):
                    #Throws runtimeerror if unable to get base
                    try:
                        refbase = classifier.bam.get_reference_base(chromosome_id, position, index=True)
                        self.assertRegexpMatches(str(refbase), '[0-4]|[ACGTN]', 'Invalid Trinucleotide_context')
                    except RuntimeError:
                        pass

#========================================
#tests for deep seuencing mode - paired
#========================================

class verify_individual_functions(unittest.TestCase,base_class):
    def setUp(self):
        pass
    
    def test_tuples(self):
        """
        ensure correctness of the tuples
        MUT-255
        test will fail due to a bug in samtools api
        """
        pass
    
    def tearDown(self):
        pass


    
