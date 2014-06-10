'''
Created on Jan 28, 2014

@author: dgrewal
'''
import classify_test
import time
import logging

mutationSeq_version="4.2.1"

#==============================================================================
# Wrapper around the test class. Runs all the tests
#==============================================================================

class init_args():  
    def __init__(self,args):
        #self.samples = ["tumour:unit_test/DAH147_A12957_3_lanes_dupsFlagged_1Mbp.bam",
        #             "normal:unit_test/DAH147N_A12969_3_lanes_hg18_dupsFlagged_1Mbp.bam", 
        #             "reference:unit_test/human_all_1Mbp.fasta", "model:model_v4.1.1.npz"]
        
        self.samples = args.samples[1:]+["model:model_v4.1.1.npz"]
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
        self.quality_threshold = 0
        self.single = False
        self.threshold = 0.5
        self.tumour_variant =2
        self.verbose = True
        self.invalid = False
    
    def set_single(self,args):
        self.single = True
        #self.samples = ["tumour:unit_test/DAH147_A12957_3_lanes_dupsFlagged_1Mbp.bam", 
        #             "reference:unit_test/human_all_1Mbp.fasta", "model:model_single_v4.0.2.npz"]
        self.samples = args.samples[1:]+["model:model_single_v4.0.2.npz"]
        self.samples=[val for val in self.samples if val.strip().split(':')[0] != 'normal']
    
    
def run_tests(args_orig):
    log_file = "./unit_test/UnitTest_run.log"
    logging.basicConfig(filename = log_file, 
                        format   = '%(asctime)s %(message)s', 
                        #datefmt = '%m/%d/%Y %I:%M:%S %p', 
                        level = logging.DEBUG)
    print "*"*80
    print "Running Single Sample"
    print "*"*80
    args = init_args(args_orig)
    args.set_single(args_orig)
    classify_test.base_class.args = args
    classify_test.suite_all_tests()
    del args
    
    time.sleep(1)
    print "*"*80
    print "Running Paired Mode"
    print "*"*80
    args = init_args(args_orig)
    classify_test.base_class.args = args
    classify_test.suite_all_tests()