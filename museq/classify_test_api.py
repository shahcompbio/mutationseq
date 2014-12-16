'''
Created on Jan 28, 2014

@author: dgrewal
'''
import classify_test
import time
import logging
import bamutils

mutationSeq_version="4.3.1"

#==============================================================================
# Get the arguments
#==============================================================================

class initargs():  
    def __init__(self):
        self.samples = ["tumour:unit_test/DAH147_A12957_3_lanes_dupsFlagged_1Mbp.bam",
                     "normal:unit_test/DAH147N_A12969_3_lanes_hg18_dupsFlagged_1Mbp.bam", 
                     "reference:unit_test/human_all.fasta", "model:model_v4.1.1.npz"]
        
        #self.samples = args.samples[1:]+["model:model_v4.1.1.npz"]
        self.all = False
        self.buffer_size = "2G"
        self.config = "metadata.config"
        self.coverage = 4
        self.deep = False
        self.export_features = None
        self.features_only = False
        self.interval = None
        self.log_file = "./unit_test/UnitTest_run.log"
        self.no_filter = False
        self.normal_variant = 25
        self.out = "./unit_test/output_unittest.vcf"
        self.positions_file = None
        self.purity = 70
        self.mapq_threshold = 0
        self.baseq_threshold = 0
        self.single = False
        self.threshold = 0.5
        self.tumour_variant =2
        self.verbose = True
        self.invalid = False
        self.indl_threshold = 0.5
    
    def set_single(self):
        self.single = True
        self.samples = ["tumour:unit_test/DAH147_A12957_3_lanes_dupsFlagged_1Mbp.bam", 
                     "reference:unit_test/human_all.fasta", "model:model_single_v4.0.2.npz"]
        #self.samples = args.samples[1:]+["model:model_single_v4.0.2.npz"]
        #self.samples=[val for val in self.samples if val.strip().split(':')[0] != 'normal']

    def set_paired(self):
        self.single = False
        self.samples = ["tumour:unit_test/DAH147_A12957_3_lanes_dupsFlagged_1Mbp.bam", 
                        "normal:unit_test/DAH147N_A12969_3_lanes_hg18_dupsFlagged_1Mbp.bam", 
                        "reference:unit_test/human_all.fasta", "model:model_single_v4.0.2.npz"]
        
    
#==============================================================================
# utility functions
#==============================================================================   
class base_class():
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
        classifier = bamutils.Classifier(args)
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

    def get_features_paired(self,tt,nt,classifier,args):
        chromosome_id = tt[-1]
        position = tt[0]
        refbase = classifier.bam.get_reference_base(chromosome_id, position, index=True)
        nonrefbases = [x for x in range(4) if x != refbase]
            
        ## ignore tumour tuples with no/few variants in the tumour or too many variants in the normal
        if not args.no_filter:
            if  tt[nonrefbases[0] + 1][0] < args.tumour_variant and \
                tt[nonrefbases[1] + 1][0] < args.tumour_variant and \
                tt[nonrefbases[2] + 1][0] < args.tumour_variant or \
                (nt[5][0] - nt[refbase + 1][0]) / nt[5][0] > (args.normal_variant / 100):
                    return None,refbase
           
        ## get corresponding reference tuple
        rt = classifier.bam.get_reference_tuple(chromosome_id, position)
        feature_set = classifier.features_module.Features(tt, nt, rt)
        return feature_set,refbase

        
    def get_features_single(self,it,classifier,args):
        chromosome_id = it[-1]
        position = it[0]
        refbase = classifier.bam.get_reference_base(chromosome_id, position, index=True)
        nonrefbases = [x for x in range(4) if x != refbase]
            
        ## ignore tumour tuples with no/few variants in the tumour or too many variants in the normal
        if not args.no_filter:
            if  it[nonrefbases[0] + 1][0] < args.tumour_variant and \
                it[nonrefbases[1] + 1][0] < args.tumour_variant and \
                it[nonrefbases[2] + 1][0] < args.tumour_variant:
                    return None,refbase
           
        ## get corresponding reference tuple
        rt = classifier.bam.get_reference_tuple(chromosome_id, position)
        feature_set = classifier.features_module.Features(it, rt)
        return feature_set,refbase

