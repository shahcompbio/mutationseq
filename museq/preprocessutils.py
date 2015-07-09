'''
Created on May 20, 2014

@author: dgrewal
'''
import logging
import warnings
from bamutils import Classifier
import pybamapi
import numpy
import resource

mutationSeq_version = "4.3.6"

class PreProcess(Classifier):
    def __init__(self,args):
        self.args = args 
        self.__check_args()
        self.__update_args()
        super(PreProcess, self).__init__(self.args) 
   
    def __update_args(self):
        samples = []
        for sample in self.args.samples:
            s = sample.split(':')
            if s[0] == 'tumour':
                tumour = s[1]
                continue
            if s[0] == 'reference':
                reference = s[1]
                
            samples.append(sample)
        
        
        if self.args.deep:
            rmdups = False
        else:
            rmdups = True
            
#         self.t_bam = pybamapi.Bam(bam=tumour, reference=reference, coverage=self.args.coverage,
#                                   rmdups=rmdups, qual_threshold=self.args.quality_threshold)
        
        self.t_bam = pybamapi.Bam(bam=tumour, reference=reference, coverage=self.args.coverage,
                                  rmdups=rmdups, mapq_threshold=self.args.mapq_threshold,
                                  baseq_threshold=self.args.baseq_threshold)
        
        self.args.samples = samples
        
        
    def __check_args(self):
        if not self.args.single:
            logging.warning('paired mode specified in preprocessing mode, setting single to True')
            warnings.warn('paired mode specified in preprocessing mode, setting single to True')
            self.args.single = True
        
        #if not self.args.all:
            #logging.warning('preprocessing mode selected, Turning filtering off')
            #warnings.warn('preprocessing mode selected, Turning filtering off')
            #self.args.all = True
            #self.args.mapq_threshold = 0
            #self.args.baseq_threshold = 0
            #self.args.no_filter = True
            #self.args.coverage = 0

        if not self.args.positions_file:
            logging.warning('positions_file not specified in preprocessing mode')
            raise Exception('positions_file not specified in preprocessing mode')
            
        
        pass
    
    def __get_allele_counts(self,chromosome,position,chr_id):
        # get the tumour tuple
        tt = self.t_bam.get_tuple(chromosome, position)
        rt = self.t_bam.get_reference_tuple(chr_id, position)
        refbase = rt[0]
        
        #if tumour tuple is None
        if tt == None:
            return 0,0
        
         
        cov = tt[5][0]
        if cov == 0:
            altbase = 'N/A'
            TR = 0
            TA = 0
        else:
            major = tt[6]
            minor = tt[7]
            if refbase == major:
                altbase = minor
            else:
                altbase = major
                
            TR = tt[refbase + 1][0]
            TA = tt[altbase + 1][0]
                
        return TR,TA
    
 #    def __flush(self):
 #        logging.info("flushing memory. Usage was: %s M",
 #                     str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 #                         / 1024))
 #
 #        # a numpy array required as an input to the random forest predictor
 #        features = numpy.array(self.features_buffer)
 #        outstrs = self.outstr_buffer
 #
 #        # empty the buffers
 #        self.features_buffer = []
 #        self.outstr_buffer = []
 #
 #        # should return a list
 #        if features is None:
 #            return [], []
 #        # gc.collect()
 #        return features, outstrs
        
    def get_features(self):
        tuples = self.bam.get_tuples(self.target_positions)
        for it in tuples: 
            self._update_coverage_info(it)
            chromosome_id = it[-1]
            position = it[0]
          
            ## get corresponding reference tuple
            rt = self.bam.get_reference_tuple(chromosome_id, position)
            
            #MUT-238 If the ref base is 4(N) ignore the position
            if rt[0] >= 4: 
                logging.error(str(position)+' position references base N and has been ignored')
                continue
            
            ## calculate features     
            feature_set = self.features_module.Features(it, rt,'n')
            temp_feature = feature_set.get_features()
            
            ## generate output string and buffer it
            chromosome = self.t_bam.get_chromosome_name(chromosome_id)
            tt = self.t_bam.get_tuple(chromosome, position)
            
            outstr = self._make_outstr(it, rt[0],tt)
            
            ##add the features to buffer if genotype is 0/1
            if outstr[-1][7] == '0/1':
                self.features_buffer.append(temp_feature)
                self.outstr_buffer.append(outstr)
            
            ## check the buffer size and flush
            if len(self.features_buffer) >= self.buffer_size:
                yield self._flush()
        
        yield self._flush()
