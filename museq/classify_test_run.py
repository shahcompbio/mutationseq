'''
Created on Jun 10, 2014

@author: dgrewal
'''
import classify_test
import unittest
import logging
mutationSeq_version="4.3.7"

def suite_all_tests():   
    suite = unittest.TestSuite()
    loader = unittest.TestLoader()

    memory_footprint = loader.loadTestsFromTestCase(classify_test.check_memory_footprint)
    dependencies = loader.loadTestsFromTestCase(classify_test.verify_dependencies)
    single_pos_error = loader.loadTestsFromTestCase(classify_test.single_position_error)
    boundary_invalid_reads = loader.loadTestsFromTestCase(classify_test.verify_boundary_invalid_reads)
    verify_outputs_test = loader.loadTestsFromTestCase(classify_test.verify_outputs)
    #ref_file = loader.loadTestsFromTestCase(classify_test.verify_reference_file)
    features_test = loader.loadTestsFromTestCase(classify_test.verify_features)
    flags = loader.loadTestsFromTestCase(classify_test.verify_flags)
    #chr_pos = loader.loadTestsFromTestCase(classify_test.verify_position_with_chr)
    tuples = loader.loadTestsFromTestCase(classify_test.verify_tuples_positions)
    individual_function = loader.loadTestsFromTestCase(classify_test.verify_individual_functions)
    get_positions = loader.loadTestsFromTestCase(classify_test.verify_get_positions)

    suite.addTests(memory_footprint)
    suite.addTests(dependencies)
    suite.addTests(single_pos_error)
    suite.addTests(boundary_invalid_reads)
    suite.addTests(verify_outputs_test)
    suite.addTests(features_test)
    suite.addTests(flags)
    suite.addTests(tuples)
    suite.addTests(individual_function)
    suite.addTests(get_positions)
    #suite.addTests(ref_file)
    #suite.addTests(chr_pos)
 
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
    


if __name__ == '__main__':
    log_file = "./unit_test/UnitTest_run.log"
    logging.basicConfig(filename = log_file, 
                        format   = '%(asctime)s %(message)s', 
                        #datefmt = '%m/%d/%Y %I:%M:%S %p', 
                        level = logging.DEBUG)

    suite_all_tests()
