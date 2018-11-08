# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:22:08 2013

@author: jtaghiyar
"""
from __future__ import division
import logging
import sys
import os
import numpy
import pybamapi
import resource
import re
import features
import features_single
import features_deep
import features_deep_single
import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
from math import log10
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
from sklearn import cross_validation
from sklearn.metrics import roc_curve, auc
from string import Template
from datetime import datetime
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import binom
from sklearn.linear_model import ElasticNetCV
from intervaltree import IntervalTree
#from intervaltree.bio import GenomeIntervalTree

mutationSeq_version = "4.3.8"
MUSEQ_VERSION = mutationSeq_version


"""
==============================================================================
Classifier class
==============================================================================
"""

class Classifier(object):

    def __init__(self, args):
        self.samples = {}
        self.args = args
        self.base = ['A', 'C', 'G', 'T', 'N']
        self.outstr_buffer = []
        self.features_buffer = []
        self.__get_buffer_size()
        self.coverage_info = [0] * 4

        self.threshold = self.args.threshold
        self.coverage = self.args.coverage
        self.mapq_threshold = self.args.mapq_threshold
        self.baseq_threshold = self.args.baseq_threshold
        self.no_filter = self.args.no_filter

        # if the all option is specified then override other options and run
        # on all positions
        if self.args.all:
            self.mapq_threshold = 0
            self.baseq_threshold = 0
            self.no_filter = True
            self.threshold = 0.0
            self.coverage = 0

        # parse the positional argument to get tumour/normal bam, reference
        # fasta and model
        for s in self.args.samples:
            self.samples[s.split(':')[0]] = s.split(':')[1]

        # check if there is a reference in the input
        if not self.samples.get("reference"):
            logging.error("error: bad input: reference must be specified")
            raise Exception("no reference file in the input.")

        self.ref = self.samples.get("reference")

        dirname = os.path.dirname(os.path.realpath(__file__))
        # check if the model is specified correctly
        model = self.samples.get("model")
        if not model:
            self.model = "model_single_v4.0.2.pickle" if self.args.single else "model_v4.1.2.pickle" 
            self.model = os.path.join(dirname, self.model)
        if model and not os.path.exists(model):
            self.model = os.path.join(dirname, model)

        if not self.model:
            logging.error("error: bad input: model must be specified in the \
                          input")
            raise Exception("no model")

        # check if there are any bam files in the input
        if not self.samples.get("normal") and not self.samples.get("tumour"):
            logging.error("error: bad input: no bam files specified in the \
                            input")
            raise Exception("no bam file")

        # check if it is single mode but there are two bam files instead of
        # one in the input
        if (self.samples.get("normal") and self.samples.get("tumour") and
                self.args.single):
            logging.error("error: bad input: single mode but two bam files\
                        specified in the input")
            raise Exception("single mode but two bam files specified")

        # check if it is not single mode but there is only one bam file in the
        # input
        if (not self.samples.get("normal") or not self.samples.get("tumour")) \
                and not self.args.single:
            logging.error("error: bad input: one bam file specified in the\
                        input but it does not seem to be the single mode.")
            raise Exception("one bam file specified but not the single mode")

        rmdups = False if self.args.count_duplicate_reads else True

        # single mode
        if self.args.single:
            if self.args.deep:
                self.features_module = features_deep_single

            else:
                self.features_module = features_single

            if not self.samples.get("tumour"):
                self.type = 'n'

                logging.info("initializing a normal Bam")
                self.bam = pybamapi.Bam(bam=self.samples.get("normal"),
                                        reference=self.ref,
                                        coverage=self.coverage,
                                        rmdups=rmdups,
                                        mapq_threshold=self.mapq_threshold,
                                        baseq_threshold=self.baseq_threshold)

            else:
                self.type = 't'

                logging.info("initializing a tumour Bam")
                self.bam = pybamapi.Bam(bam=self.samples.get("tumour"),
                                        reference=self.ref,
                                        coverage=self.coverage,
                                        rmdups=rmdups,
                                        mapq_threshold=self.mapq_threshold,
                                        baseq_threshold=self.baseq_threshold)

        # paired mode
        else:
            if self.args.deep:
                self.features_module = features_deep
            else:
                self.features_module = features

            logging.info("initializing a PairedBam")
            self.bam = pybamapi.PairedBam(tumour=self.samples.get("tumour"),
                                          normal=self.samples.get("normal"),
                                          reference=self.samples.get(
                                              "reference"),
                                          coverage=self.coverage,
                                          rmdups=rmdups,
                                          mapq_threshold=self.mapq_threshold,
                                          baseq_threshold=self.baseq_threshold)

        if not self.bam.is_matched_reference():
            logging.error("mismatched reference, sounds like the input\
                        reference is not the same as the reference used\
                        for alignment")
            raise Exception("mismatched reference")

        if self.args.deep:
            self.manifest = self.__parse_manifest()

    def __parse_manifest(self):
        if not self.args.manifest:
            raise Exception('Manifest file is required in deep mode')

        man_stream = open(self.args.manifest)
        gtree = defaultdict(lambda: IntervalTree())

        for line in man_stream:
            line = line.strip().split()
            if line[0] == 'chrom':
                continue
            if len(line) > 3:
                raise Exception('Please ensure that the '
                                'manifest is formatted properly')
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])

            #empty intervals (points) are not allowed.
            assert end>start, 'amplicon end should be greater than start.'\
                              ' Please check the manifest file'

            gtree[chrom].addi(start, end)
        return gtree

    def __get_flanking_regions(self, chromosome, start,stop):
        if start == stop:
            vals = self.manifest[chromosome][start]
        else:
            vals = self.manifest[chromosome][start:stop]
        vals = sorted(vals)
            
        if len(vals) > 1:
            raise Exception('The position %s:%s-%s falls in more than one '  
                            'region in manifest file' %(chromosome,start,stop))
                
        elif vals == []:
            logging.info('Skipping position %s:%s-%s as it doesn\'t fall in any '
                         'amplicon region' %(chromosome,start,stop))
            return None,None
        else:
            start = vals[0].begin
            end = vals[0].end
            return start, end

    def __get_buffer_size(self):
        s = re.split('(\d+)', self.args.buffer_size)
        d = int(s[1])
        l = s[2].upper()
        if d < 0 or l not in ('G', 'M'):
            # relax the restriction on the memory usage
            self.buffer_size = float('inf')

        elif d < 100 and l == 'M':
            logging.warning("warning: buffer size is ingnored. \
                            It should be >= 100M.")
            self.buffer_size = float('inf')

        elif l == 'G':
            # every million output and feature strings together
            # takes about 2G of memory. Also, ~50M memory is required
            # to initialize and about another ~50M is required to fit the
            # model.
            self.buffer_size = d * 200000

        elif l == 'M':
            self.buffer_size = (d / 1024) * 200000

       
    def __parse_position_file(self, pch=':'):
        positions = []
        posfile = open(self.args.positions_file)
        for line in posfile:
            temp_tp = self.__parse_position(line)
            if temp_tp is not None:
                positions.extend([temp_tp])
        return positions
    
    def __parse_position(self, pos, pch=':'):
        chromosome = pos.split(pch)[0]
        try:
            pos_range = pos.split(pch)[1]
            start = int(pos_range.split('-')[0])
            
            try:
                stop = int(pos_range.split('-')[1])

            except:
                stop = start
            return [chromosome, start, stop]

        except:
            return [chromosome, None, None]
    
    def __filter_positions(self, target_positions, pch = ':'):
        int_pos = self.__parse_position(self.args.interval, pch)
        output = []
        for temp_tp in target_positions:
            if int_pos[0] == temp_tp[0]:
                #if only chr is provided
                if None in int_pos:
                    output.append(temp_tp)
                    continue
                start = max(int_pos[1],temp_tp[1])
                stop = min(int_pos[2],temp_tp[2])
                if start > stop:
                    continue
                output.append([int_pos[0],start,stop])
        return output
        
    def get_positions_deep_filter(self, target_positions):
        '''
        gets the overlapping amplicons for each positions and 
        adds the intersection to the target list
        '''
        output = []
        for val in target_positions:
            if None in val:
                output.append(val)
                continue
            
            if val[1] == val[2]:
                regions = self.manifest[val[0]][val[1]]
            else:
                regions = self.manifest[val[0]][val[1]:val[2]]

            for region in regions:
                start = max(region.begin,val[1])
                stop = min(region.end,val[2])
                if start > stop:
                    ##TODO: exception might make more sense here.
                    continue
                output.append([val[0],start,stop])
        return output

        
    def get_positions_deep(self, target_positions):
        '''
        If target positions isn't empty then filter them based on the manifest
        else return all positions from the manifest
        '''
        if target_positions:
            target_positions = self.get_positions_deep_filter(target_positions)
        else:
            #converting the tree to list    
            [target_positions.extend([[val[0], i.begin, i.end] for i in val[1]])
             for val in self.manifest.items()]
        return target_positions
    
    def get_positions(self, pch=':'):
        target_positions = []
        
        if self.args.positions_file:
            target_positions = self.__parse_position_file()
            if self.args.interval:
                target_positions = self.__filter_positions(target_positions)
            
        if self.args.deep:
            target_positions = self.get_positions_deep(target_positions)
            if self.args.interval:
                target_positions = self.__filter_positions(target_positions)

        if target_positions == [] and self.args.interval and not self.args.positions_file and not self.args.deep:
            target_positions.append(self.__parse_position(self.args.interval))

                
        #if no pos found- then run over whole genome
        if target_positions == [] and not self.args.interval and not self.args.positions_file and not self.args.deep:
            # get all the common chromosome names
            # chromosome names in tumour bam
            tcn = self.bam.get_refnames().keys()
            # chromosome names in normal bam
            ncn = self.bam.get_refnames().keys()
            chr_names = set(tcn).intersection(set(ncn))
        
            for cn in chr_names:
                temp_tp = [cn, None, None]
                target_positions.append(temp_tp)
            
        #print target_positions
        #print self.args.interval
        #print self.args.positions_file    
        self.target_positions = target_positions

    @staticmethod
    def get_phred_score(prob, typ):
        """
        convert probability to phred scaled likelihood
        when returning likelihood, higher value has lower probability
        when quality: higher value for high prob
        """
        if typ == 'likelihood':
            if prob == 0:
                prob = 255
            elif prob == 1:
                prob = 0
            else:
                prob = -10 * log10(prob)
        elif typ == 'quality':
            if prob == 0:
                prob = 0
            elif prob == 1:
                prob = 100
            else:
                prob = -10 * log10(1 - prob)

        return prob

    def __get_genotype(self, nonref_count, count_all, typ):
        """
        uses snvmix priors and binomial test to calculate
        genotype and their likelihoods
        """
        if typ == 'n':
            # snvmix priors - normal allele frequency
            aa_freq = 0.01
            ab_freq = 0.50
            bb_freq = 0.99
        elif typ == 't':
            # tumour allele frequency
            aa_freq = 0.01
            ab_freq = 0.30
            bb_freq = 0.90
        else:
            raise ValueError('Invalid input parameter.')

        prior = [aa_freq, ab_freq, bb_freq]

        # binomial pmf for three different probabilities
        bnom_pmf = [binom.pmf(nonref_count, count_all, p) for p in prior]

        # ensure binom results sum to 1 (get probabilities)
        bnom_sum = 1e-150 if sum(bnom_pmf) == 0 else sum(bnom_pmf)
        bnom_pmf = [round(val / bnom_sum, 4) for val in bnom_pmf]

        pr_aa, pr_ab, pr_bb = bnom_pmf

        if pr_aa == max(pr_aa, pr_ab, pr_bb):
            gtyp = '0/0'
        elif pr_ab == max(pr_aa, pr_ab, pr_bb):
            gtyp = '0/1'
        elif pr_bb == max(pr_aa, pr_ab, pr_bb):
            gtyp = '1/1'

        pr_aa = self.get_phred_score(pr_aa, typ='likelihood')
        pr_ab = self.get_phred_score(pr_ab, typ='likelihood')
        pr_bb = self.get_phred_score(pr_bb, typ='likelihood')

        prob_lkl = [int(pr_aa), int(pr_ab), int(pr_bb)]
        prob_lkl = ','.join([str(v) for v in prob_lkl])

        return gtyp, prob_lkl


    def _make_outstr(self, tumtup, refbase, normtup=None):
        """
        generate output string with all the relavant values
        for vcf output
        """
        # Assuming that coverage cannot be 0 since we skip those tuples in
        # pybam

        # alternative base, tt[6] is major and tt[7] is minor
        altbase = tumtup[7] if refbase == tumtup[6] else tumtup[6]

        # calculate ratio of reads with an indl
        ratio = sum(tumtup[-2]) / tumtup[5][0]
        filter_flag = "INDL" if ratio > self.args.indl_threshold else None

        tumdata = [tumtup[refbase + 1][0], tumtup[altbase + 1][0],
                   tumtup[9], tumtup[11], tumtup[5][0]]

        in_typ = self.type if self.args.single else 't'
        tumdata.extend(self.__get_genotype(tumdata[1], tumtup[5][0], in_typ))

        normdata = []
        if normtup:
            normdata.extend([normtup[refbase + 1][0], normtup[altbase + 1][0],
                             normtup[9], normtup[11], normtup[5][0]])
            normdata.extend(self.__get_genotype(normdata[1],
                                                normtup[5][0], 'n'))

        # reserved for database ID, to be filled later
        out_id = "."

        # tri_nucleotide context
        chromosome_id = tumtup[-1]
        position = tumtup[0]
        trinuc = self.bam.get_trinucleotide_context(chromosome_id, position)

        # get chromosome name of the given chromosome ID
        chromosome_name = self.bam.get_chromosome_name(chromosome_id)

        outstr = [chromosome_name, position, out_id, refbase, altbase,
                  filter_flag, trinuc, tumdata, normdata]

        return outstr


    def _flush(self):
        logging.info("flushing memory. Usage was: %s M",
                     str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                         / 1024))

        # a numpy array required as an input to the random forest predictor
        features = numpy.array(self.features_buffer)
        outstrs = self.outstr_buffer

        # empty the buffers
        self.features_buffer = []
        self.outstr_buffer = []

        # should return a list
        if features is None:
            return [], []
        # gc.collect()
        return features, outstrs

    def _update_coverage_info(self, tt=None, nt=None):
        if tt is not None:
            self.coverage_info[0] += tt[-3]
            self.coverage_info[1] += 1
        if nt is not None:
            self.coverage_info[2] += nt[-3]
            self.coverage_info[3] += 1

    def get_features(self):
        logging.info("getting features")

        if self.args.single:
            if self.args.deep:
                return self.__get_features_single_deep()
            else:
                return self.__get_features_single()
        else:
            if self.args.deep:
                return self.__get_features_paired_deep()
            else:
                return self.__get_features_paired()

    def __get_features_single(self):
        tuples = self.bam.get_tuples(self.target_positions)
        for it in tuples:
            self._update_coverage_info(it)

            chromosome_id = it[-1]
            position = it[0]

            refbase = self.bam.get_reference_base(chromosome_id, position,
                                                  index=True)

            nonrefbases = [x for x in range(4) if x != refbase]

            # ignore tumour tuples with no/few variants in the bam file
            if not self.no_filter:
                if it[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
                        it[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
                        it[nonrefbases[2] + 1][0] < self.args.tumour_variant:
                    continue

            # get corresponding reference tuple
            rt = self.bam.get_reference_tuple(chromosome_id, position)

            # MUT-238 If the ref base is 4(N) ignore the position
            if rt[0] >= 4:
                chromosome_name = self.bam.get_chromosome_name(chromosome_id)
                logging.error("%s:%s position references base N and has been " \
                             "ignored" % (chromosome_name ,str(position)))
                continue

            # calculate features
            feature_set = self.features_module.Features(it, rt, self.type)
            temp_feature = feature_set.get_features()
            self.features_buffer.append(temp_feature)

            # generate output string and buffer it
            outstr = self._make_outstr(it, rt[0], normtup=None)
            self.outstr_buffer.append(outstr)

            # check the buffer size and flush
            if len(self.features_buffer) >= self.buffer_size:
                yield self._flush()

        yield self._flush()

    def __get_features_paired(self):
        tuples = self.bam.get_tuples(self.target_positions)
        for tt, nt in tuples:
            self._update_coverage_info(tt, nt)
            chromosome_id = tt[-1]
            position = tt[0]

            refbase = self.bam.get_reference_base(chromosome_id, position,
                                                  index=True)
            nonrefbases = [x for x in range(4) if x != refbase]

            # ignore tumour tuples with no/few variants in the tumour or too
            # many variants in the normal
            if not self.no_filter:
                if tt[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
                        tt[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
                        tt[nonrefbases[2] + 1][0] < self.args.tumour_variant or \
                        (nt[5][0] - nt[refbase + 1][0]) / nt[5][0] >\
                        (self.args.normal_variant / 100):
                    continue

            # get corresponding reference tuple
            rt = self.bam.get_reference_tuple(chromosome_id, position)

            # MUT-238 If the ref base is 4(N) ignore the position
            if rt[0] >= 4:
                chromosome_name = self.bam.get_chromosome_name(chromosome_id)
                logging.error("%s:%s position references base N and has been " \
                             "ignored" %  (chromosome_name ,str(position)))
                continue

            # calculate features and buffer it
            feature_set = self.features_module.Features(tt, nt, rt)
            temp_feature = feature_set.get_features()
            self.features_buffer.append(temp_feature)

            # generate output string and buffer it
            outstr = self._make_outstr(tt, rt[0], nt)
            self.outstr_buffer.append(outstr)

            # check the buffer size and flush
            if len(self.features_buffer) >= self.buffer_size:
                yield self._flush()

        yield self._flush()

    def __get_tuples(self, chromosome, start, end):
        tt_int= []
        rt_int = []
        nt_int = []
        tuples = self.bam.get_tuples([[chromosome,start,end]])
        if self.args.single:
            for tt in tuples:
                tt_int.append(tt)
                rt = [tt[0], self.bam.get_reference_tuple(tt[-1],tt[0])]
                rt_int.append(rt)            
        else:
            for tt,nt in tuples:
                tt_int.append(tt)
                nt_int.append(nt)
                rt = [tt[0], self.bam.get_reference_tuple(tt[-1],tt[0])]
                rt_int.append(rt)
        return tt_int,rt_int,nt_int
    
    
    def __get_features_paired_deep(self):
        """
        get the target_positions and then get the flanking region for that position
        get tuples for the flanking region 
        classify the positions
        target positions cannot have full chr as manifest is required
        any range in target positions shouldn't fall in 2 manifest regions-see get_positions
        """
        for val in self.target_positions:
            chromosome = val[0]
            start,end = self.__get_flanking_regions(chromosome, val[1], val[2])
            tt_int, rt_int, nt_int = self.__get_tuples(chromosome, start, end)
            
            for tt,rt,nt in zip(tt_int, rt_int, nt_int):
                if tt[0] < val[1] or tt[0] > val[2]:
                    continue
                if tt[0] != rt[0]:
                    raise Exception()
                rt = rt[1]
                self._update_coverage_info(tt, nt)
                position = tt[0]
                chromosome_id = tt[-1]
                refbase = self.bam.get_reference_base(chromosome_id,
                                                      position,
                                                      index=True)
                nonrefbases = [x for x in range(4) if x != refbase]

                # ignore tumour tuples with no/few variants in the tumour or too
                # many variants in the normal
                if not self.no_filter:
                    if tt[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
                       tt[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
                       tt[nonrefbases[2] + 1][0] < self.args.tumour_variant or \
                       (nt[5][0] - nt[refbase + 1][0]) / nt[5][0] >\
                       (self.args.normal_variant / 100):
                        continue

                # MUT-238 If the ref base is 4(N) ignore the position
                if rt[0] >= 4:
                    chromosome_name = self.bam.get_chromosome_name(chromosome_id)
                    logging.error("%s:%s position references base N and has been " \
                                 "ignored" % (chromosome_name ,str(position)))
                    continue

                tt_bg = [tval for tval in tt_int if not tval == tt]
                nt_bg = [nval for nval in nt_int if not nval == nt]
            
                feature_set = self.features_module.Features(tt, nt, rt, tt_bg,
                                                            nt_bg, rt_int)
                temp_feature = feature_set.get_features()
                self.features_buffer.append(temp_feature)

                # generate output string and buffer it
                outstr = self._make_outstr(tt, rt[0], nt)
                self.outstr_buffer.append(outstr)

                # check the buffer size and flush
                if len(self.features_buffer) >= self.buffer_size:
                    yield self._flush()

        yield self._flush()

    def __get_features_single_deep(self):
        for val in self.target_positions:
            chromosome = val[0]
            start,end = self.__get_flanking_regions(chromosome, val[1], val[2])
            it_int, rt_int, _ = self.__get_tuples(chromosome, start, end)
            
            for it,rt in zip(it_int, rt_int):
                if it[0] < val[1] or it[0] > val[2]:
                    continue
                if it[0] != rt[0]:
                    raise Exception()
                rt = rt[1]
                self._update_coverage_info(it)
                position = it[0]
                chromosome_id = it[-1]
                refbase = self.bam.get_reference_base(chromosome_id,
                                                      position,
                                                      index=True)
                nonrefbases = [x for x in range(4) if x != refbase]
        
                # ignore tumour tuples with no/few variants in the tumour or too
                # many variants in the normal
                if not self.no_filter:
                    if it[nonrefbases[0] + 1][0] < self.args.tumour_variant and \
                       it[nonrefbases[1] + 1][0] < self.args.tumour_variant and \
                       it[nonrefbases[2] + 1][0] < self.args.tumour_variant:
                        continue

                # MUT-238 If the ref base is 4(N) ignore the position
                if rt[0] >= 4:
                    chromosome_name = self.bam.get_chromosome_name(chromosome_id)
                    logging.error("%s:%s position references base N and has been " \
                                 "ignored" % (chromosome_name ,str(position)))
                    continue

                it_bg = [tval for tval in it_int if not tval == it]

                feature_set = self.features_module.Features(it, rt, it_bg, rt_int,
                                                            self.type)

                temp_feature = feature_set.get_features()
                self.features_buffer.append(temp_feature)

                # generate output string and buffer it
                outstr = self._make_outstr(it, rt[0], normtup=None)
                self.outstr_buffer.append(outstr)

                # check the buffer size and flush
                if len(self.features_buffer) >= self.buffer_size:
                    yield self._flush()
        yield self._flush()  

    def __load_model(self):
        try:
            logging.info("loading model")
            return joblib.load(self.model)

        except Exception as e:
            if hasattr(e, 'strerror'):
                logging.error("error loading model: "+e.strerror)
            else:
                logging.error("error loading model: "+str(e))
            raise Exception("failed to load model")

    def __verify_model_features(self, model):
        features = self.features_module.Features()
        if not model.version == features.version:
            return False
        if not model.name == features.name:
            return False
        return True

    def predict(self, features_outstrs):
        # model = self.__fit_model()
        model = self.__load_model()

        # verify the model against the features
        if not self.__verify_model_features(model):
            logging.error('The features and the model do not match')
            raise Exception('mismatched model')

        logging.info("predicting probabilities ")
        for features, outstrs in features_outstrs:
            if len(features) == 0:
                continue

            if hasattr(model, 'coefs'):
                zeros = [i for i, val in enumerate(model.coefs.tolist())
                         if val == 0]

                features = [[featureval for i, featureval in enumerate(feature)
                             if i not in zeros] for feature in features]

                probabilities = model.predict_proba(features)
            else:
                probabilities = model.predict_proba(features)

            # return only probabilities of being somatic
            probabilities = [x[1] for x in probabilities]

            yield probabilities, outstrs

    def __update_header(self):
        tumour_cov_mean = self.coverage_info[0]/(self.coverage_info[1] + 1e-5)
        normal_cov_mean = self.coverage_info[2]/(self.coverage_info[3] + 1e-5)

        # Doing a byte by byte substitution to avoid rewriting file.
        try:
            header = []
            outfile = open(self.args.out, 'r+')
            for line in outfile:
                if not line[0] == '#':
                    break
                header.append(line)

            len_header = len(''.join(header))

            # the length of the mean should be exactly 10 chars
            tumour_cov_mean = str("%.10f" % tumour_cov_mean)[:10]
            normal_cov_mean = str("%.10f" % normal_cov_mean)[:10]

            for i, val in enumerate(header):
                if '$TUMOURCOV' in val:
                    header[i] = header[i].replace('$TUMOURCOV',
                                                  tumour_cov_mean)
                if '$NORMALCOV' in val:
                    header[i] = header[i].replace('$NORMALCOV',
                                                  normal_cov_mean)

            # if the new header is not the same length as old header, it will
            # overwrite the record in next line/leave
            # some text from old header intact
            if not len(''.join(header)) == len_header:
                raise Exception('The new header is not the same length as the\
                                older ')

            outfile.seek(0)
            for val in header:
                outfile.write(val)
            outfile.flush()

        except Exception as e:
            logging.error('Error writing coverage information to the header\n')
            logging.error(str(e) + '\n')
            logging.info(
                'The mean coverage for Tumour: ' +
                str(tumour_cov_mean))
            logging.info(
                'The mean coverage for normal: ' +
                str(normal_cov_mean))

    def _meta_data(self):
        tumour = self.samples.get("tumour")
        normal = self.samples.get("normal")
        reference = self.samples.get("reference")
        model = self.model

        contigs = self.bam.get_reference_chromosome_lengths()
        contigs = ['##contig=<ID=%s,length=%s>' %(k,v) for k,v in contigs.iteritems()]
        contigs = '\n'.join(contigs)

        if tumour is None:
            tumour = "N/A"

        elif normal is None:
            normal = "N/A"

        cfg_file = self.args.config
        if not cfg_file or not os.path.exists(cfg_file):
            cfg_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'metadata.config')
        try:
            cfg_file = open(cfg_file, 'r')
            header = ""

            for hdrline in cfg_file:
                hdrline = Template(hdrline).substitute(
                        DATETIME=datetime.now().strftime("%Y%m%d"),
                        VERSION=MUSEQ_VERSION,
                        REFERENCE=reference,
                        TUMOUR=tumour,
                        NORMAL=normal,
                        MODEL=model,
                        THRESHOLD=self.args.threshold,
                        INDLTHRESHOLD=self.args.indl_threshold,
                        MAPQTHRESHOLD=self.bam.mapq_threshold,
                        BASEQTHRESHOLD=self.bam.baseq_threshold,
                        COVERAGE=self.bam.coverage,
                        RMDUPS=self.bam.rmdups,
                        CONTIG=contigs)

                # add format section headings
                if hdrline.startswith('#CHROM'):
                    hdrline = hdrline.strip('\n')
                    if self.args.single:
                        hdrline += '\tTUMOUR' if self.type == 't' else '\tNORMAL'
                    else:
                        hdrline += '\tTUMOUR\tNORMAL'
                    hdrline += '\n'

                header += hdrline
            cfg_file.close()

            return header
        except AttributeError as exc:
            exc = exc.strerror if hasattr(exc, 'strerror') else str(exc)
            logging.warning(
                "warning: failed to load metadata file due to error: %s",
                exc)
            return

        except KeyError as exc:
            exc = exc.strerror if hasattr(exc, 'strerror') else str(exc)
            logging.warning(
                "warning: failed to load metadata file due to error: %s",
                exc)
            return


    def print_results(self, probabilities_outstrs):
        """
        collects the outstr and probabilities and writes to
        output file in vcf format
        """
        # open the output vcf file to write
        out = open(self.args.out, 'w')
        out_path = str(self.args.out)

        format_str = 'RC:AC:NI:ND:DP:GT:PL'

        # print the vcf header to the output
        header = self._meta_data()
        if header is not None:
            print >> out, header.strip()

        # print the results
        for probabilities, outstrs in probabilities_outstrs:

            logging.info("printing results to: " + out_path)
            for i in xrange(len(probabilities)):
                outstr = outstrs[i]
                prob = probabilities[i]

                filter_flag = outstr[5]

                # do not print positions with p < threshold if --all option is
                # not set
                if not self.args.all and prob < self.args.threshold:
                    continue

                # set the filter_flag
                if filter_flag is None:
                    if prob >= self.args.threshold:
                        filter_flag = "PASS"
                    else:
                        filter_flag = "FAIL"

                info_str = ["PR=", "%.2f" % prob, ";TC=", outstr[6]]

                if not self.args.single:
                    tum_str = ':'.join([str(v) for v in outstr[7]])
                    norm_str = ':'.join([str(v) for v in outstr[8]])
                else:
                    samp_str = ':'.join([str(v) for v in outstr[7]])

                info_str = ''.join(info_str)

                # calculate phred quality
                phred_quality = self.get_phred_score(prob, typ='quality')

                # alternative base
                altbase = self.base[outstr[4]]

                # make sure it is all strings
                outstr = [outstr[0], outstr[1], outstr[2],
                          self.base[outstr[3]], altbase,
                          "%.2f" % phred_quality, filter_flag, info_str,
                          format_str]
                if self.args.single:
                    outstr.append(samp_str)
                else:
                    outstr.append(tum_str)
                    outstr.append(norm_str)

                outstr = [str(outstrval) for outstrval in outstr]

                print >> out, "\t".join(outstr)

        out.close()


    def get_feature_names(self):
        tmp_obj = self.features_module.Features()
        names = tmp_obj.get_feature_names()

        return names

    def export_features(self, features):
        version = self.features_module.Features().version
        names = ['chromosome', 'position'] + self.get_feature_names()

        if not hasattr(self, 'model_obj'):
            self.model_obj = self.__load_model()

        with open(self.args.export_features, 'w') as export_file:
            print >> export_file, "##features_version:" + version
            print >> export_file, "\t".join(names)


            for features, outstrs in features:
                if len(features) == 0:
                    continue

                if hasattr(self.model_obj, 'coefs'):
                    zeros = [i for i, val in enumerate(self.model_obj.coefs.tolist())
                             if val == 0]

                    features = [[featureval for i, featureval in enumerate(feature)
                             if i not in zeros] for feature in features]

                for feature,outstr in zip(features,outstrs):
                    self.features_buffer.append(feature)
                    self.outstr_buffer.append(outstr)
                    chromosome = str(outstr[0])
                    position = str(outstr[1])
                    feature = map(str, feature.tolist())
                    output_string = '\t'.join([chromosome, position] +feature)+'\n'
                    export_file.write(output_string)
                    if len(self.outstr_buffer)>self.buffer_size:
                        yield self._flush()
        yield self._flush()

    def print_features(self, features):
        '''
        since we aren't predicting, just exhaust the generator.
        '''
        for _, _ in features:
            pass

"""
==============================================================================
Trainer class
==============================================================================
"""


class Trainer(object):

    def __init__(self, args):
        self.args = args
        self.filename = os.path.basename(self.args.out)

        if self.args.single:
            if self.args.deep:
                self.feature_module = features_deep_single
                self.rmdups = False
            else:
                self.feature_module = features_single
                self.rmdups = True

        elif self.args.deep:
            self.feature_module = features_deep
            self.rmdups = False

        else:
            self.feature_module = features
            self.rmdups = True

        if self.args.deep:
            self.mapq_threshold = 10
            self.baseq_threshold = 10
        else:
            self.mapq_threshold = 0
            self.baseq_threshold = 0

    def __parse_manifest(self, manifest_file):
        manifest = defaultdict(list)
        if not self.args.deep:
            logging.error('manifest is only required in deep mode')
            raise Exception('Manifest file is only required in deep mode')

        if not manifest_file:
            return None

        man_stream = open(manifest_file)
        for line in man_stream:
            line = line.strip().split()
            if line[0] == 'chrom':
                continue
            chrom = line[0]
            pos = line[1]
            ref = line[2]
            alt = line[3]
            start = line[4]
            end = line[5]
            manifest[(chrom, pos)].append((ref, alt, start, end))
        return manifest

    def __get_flanking_regions(self, chromosome, position, manifest):
        if manifest:
            vals = manifest.get((chromosome, position))
            if vals is not None:
                start = vals[2]
                end = vals[3]
                return start, end

        # If manifest file is not provided or region not available in manifest
        start = position - 25
        end = position + 26
        return start, end

    def __isvalid_label(self, labels):
        for l in labels:
            if l not in ("SOMATIC", "WILDTYPE", "GERMLINE", "HET", "HET_ONE",
                         "HET_GERMLINE", "HOM", "HOM_ONE", "HOM_GERMLINE",
                         "CLEAN"):
                return False

        return True

    def __parse_infiles(self, infiles):
        self.data = defaultdict(list)

        positive_labels = [x.upper() for x in self.args.labels.split(',')]
        if not self.__isvalid_label(positive_labels):
            logging.error("unknown labels specified in the input")
            raise Exception("unknown labels specified in the input")

        for case in infiles:
            tfile = None
            nfile = None
            rfile = None
            manfile = None

            if self.args.deep:
                contamination = (float(10000), float(10000), float(70),
                                 float(0))
            else:
                contamination = (float(30), float(30), float(70), float(0))

            for line in case:
                l = line.strip().split()
                if len(l) < 3:
                    continue

                # parse the line
                if l[0] == "#":
                    if self.args.single:
                        if l[1] == "tumour":
                            tfile = nfile = l[2]
                            self.type = 't'
                        elif l[1] == "normal":
                            tfile = nfile = l[2]
                            self.type = 'n'
                        # if l[1] == "tumour" or l[1] == "normal":
                        #   tfile = nfile = l[2]

                    else:
                        if l[1] == "tumour":
                            tfile = l[2]

                        elif l[1] == "normal":
                            nfile = l[2]

                    if l[1] == "reference":
                        rfile = l[2]

                    if l[1] == "manifest":
                        manfile = l[2]

                    elif l[1] == "contamination":
                        contamination = (float(l[2]), float(l[3]), float(l[4]),
                                         float(1))

                    continue

                # check if required bam files/reference are specified in the
                # training data
                if not all((tfile, nfile, rfile)):
                    logging.warning("%s' does not contain the required paths "\
                                    "to bam/reference files"
                                    % os.path.basename(case.name))
                    continue

                chromosome = l[0]
                position = int(l[1])
                label_name = l[2]
                if label_name in positive_labels:
                    label = 1

                else:
                    label = -1

                if self.args.deep:
                    if manfile:
                        manifest = self.__parse_manifest(manfile)
                        start, end = self.__get_flanking_regions(chromosome,
                                                                 position,
                                                                 manifest)

                        indexes = manifest.get((chromosome, position))[:2]
                    else:
                        start = position - 25
                        end = position + 26
                        indexes = None
                    self.data[(tfile, nfile, rfile)].append((chromosome,
                                                             position, label,
                                                             contamination,
                                                             label_name,
                                                             indexes, start,
                                                             end))
                else:
                    self.data[(tfile, nfile, rfile)].append((chromosome,
                                                             position, label,
                                                             contamination,
                                                             label_name))

    def __get_features(self):
        if self.args.deep:
            features = self.__get_features_deep()
        else:
            features = self.__get_features_nondeep()

    def __get_features_nondeep(self):
        features_buffer = []
        labels_buffer = []
        keys_buffer = []
        file_stream_w = open(self.args.out + '_feature_db_train.txt', 'w')

        for tfile, nfile, rfile in self.data.keys():
            logging.info(tfile)
            if not self.args.single:
                logging.info(nfile)

            t_bam = pybamapi.Bam(bam=tfile, reference=rfile, coverage=1,
                                 rmdups=self.rmdups,
                                 mapq_threshold=self.mapq_threshold,
                                 baseq_threshold=self.baseq_threshold)

            if not self.args.single:
                n_bam = pybamapi.Bam(bam=nfile, reference=rfile, coverage=1,
                                     rmdups=self.rmdups,
                                     mapq_threshold=self.mapq_threshold,
                                     baseq_threshold=self.baseq_threshold)

            for chromosome, position, label, c, label_name in\
                    self.data[(tfile, nfile, rfile)]:
                # unpack pos to get flanking regions in deep mode
                chromosome_id = t_bam.get_chromosome_id(chromosome)
                tt = t_bam.get_tuple(chromosome, position)
                if not self.args.single:
                    nt = n_bam.get_tuple(chromosome, position)

                rt = t_bam.get_reference_tuple(chromosome_id, position)

                # MUT-238 If the ref base is 4(N) ignore the position
                if rt[0] >= 4:
                    chromosome_name = self.bam.get_chromosome_name(chromosome_id)
                    logging.error("%s:%s position references base N and has been " \
                                 "ignored" % (chromosome_name ,str(position)))
                    continue

                # check for None tuples
                if self.args.single:
                    if not all([tt, rt]):
                        logging.warning(" ".join(["None tuple", tfile, rfile]))
                        continue

                elif not all([tt, nt, rt]):
                    logging.warning(
                        " ".join(["None tuple", tfile, nfile, rfile]))
                    continue

                # calculate features
                if self.args.single:
                    feature_set = self.feature_module.Features(tt, rt,
                                                               self.type)
                else:
                    feature_set = self.feature_module.Features(tt, nt, rt)

                    
                temp_features = feature_set.get_features()

                outstr = ';'.join([rfile, nfile, tfile, chromosome,
                                  '\t'.join([str(position),
                                            str(temp_features),
                                            label_name])]) + '\n'
                file_stream_w.write(outstr)

                features_buffer.append(temp_features)
                labels_buffer.append(label)
                keys_buffer.append((rfile, nfile, tfile, chromosome, position,
                                    label))

        file_stream_w.close()
        self.features = numpy.array(features_buffer)
        self.labels = numpy.array(labels_buffer)
        self.keys = numpy.array(keys_buffer)
        self.feature_set_name = feature_set.name
        self.feature_set_version = feature_set.version

    def __get_features_deep(self):
        features_buffer = []
        labels_buffer = []
        keys_buffer = []
        file_stream_w = open(self.args.out + '_feature_db_train.txt', 'w')

        for tfile, nfile, rfile in self.data.keys():
            logging.info(tfile)
            if not self.args.single:
                logging.info(nfile)

            t_bam = pybamapi.Bam(bam=tfile, reference=rfile, coverage=1,
                                 rmdups=self.rmdups,
                                 mapq_threshold=self.mapq_threshold,
                                 baseq_threshold=self.baseq_threshold)

            if not self.args.single:
                n_bam = pybamapi.Bam(bam=nfile, reference=rfile, coverage=1,
                                     rmdups=self.rmdups,
                                     mapq_threshold=self.mapq_threshold,
                                     baseq_threshold=self.baseq_threshold)

            for chromosome, pos, label, _, label_name, indexes, start, end in self.data[
                    (tfile, nfile, rfile)]:
                # unpack pog to get flanking regions in deep mode
                chr_id = t_bam.get_chromosome_id(chromosome)

                # get all tuples in flanking region and remove the
                # pos(foreground)
                tt_bg = t_bam.get_tuples([[chromosome, start, end]])
                tt_bg = [tt for tt in tt_bg if tt[0] != pos]

                nt_bg = n_bam.get_tuples([[chromosome, start, end]])
                nt_bg = [nt for nt in nt_bg if nt[0] != pos]

                rt_bg = []
                for posval in range(start, end + 1):
                    if posval == pos:
                        continue
                    rt_bg.append([posval, t_bam.get_reference_tuple(chr_id,
                                                                    posval)])

                chromosome_id = t_bam.get_chromosome_id(chromosome)
                tt = t_bam.get_tuple(chromosome, pos)
                if not self.args.single:
                    nt = n_bam.get_tuple(chromosome, pos)

                rt = t_bam.get_reference_tuple(chromosome_id, pos)

                # MUT-238 If the ref base is 4(N) ignore the pos
                if rt[0] >= 4:
                    chromosome_name = self.bam.get_chromosome_name(chromosome_id)
                    logging.error("%s:%s position references base N and has been " \
                                 "ignored" % (chromosome_name ,str(pos)))
                    continue

                # check for None tuples
                if self.args.single:
                    if not all([tt, rt]):
                        logging.warning(" ".join(["None tuple", tfile, rfile]))
                        continue

                elif not all([tt, nt, rt]):
                    logging.warning(" ".join(["None tuple", tfile, nfile,
                                              rfile]))
                    continue

                # calculate features
                feature_set = self.feature_module.Features(tt, nt, rt, tt_bg,
                                                           nt_bg, rt_bg,
                                                           indexes)

                temp_features = feature_set.get_features()
                
                file_stream_w.write(rfile + ';' + nfile + ';' + tfile + ';' +
                                    chromosome + ';' + str(pos) + '\t' +
                                    str(temp_features) + '\t' + label_name +
                                    '\n')

                features_buffer.append(temp_features)
                labels_buffer.append(label)
                keys_buffer.append((rfile, nfile, tfile, chromosome,
                                    pos, label))

        file_stream_w.close()
        self.features = numpy.array(features_buffer)
        self.labels = numpy.array(labels_buffer)
        self.keys = numpy.array(keys_buffer)
        self.feature_set_name = feature_set.name
        self.feature_set_version = feature_set.version

    def generate(self, infiles=None):
        if infiles is None:
            infiles = self.args.infiles

        logging.info("parsing infiles")
        self.__parse_infiles(infiles)

        logging.info("getting features")
        self.__get_features()


#     def load(self):
#         try:
#             npz = numpy.load(self.args.model)
#
#         except:
#             logging.error("failed to load the model: " + self.args.model)
#             raise Exception("failed to load the model")
#
#         self.version = npz["arr_0"]
#         self.features = npz["arr_1"]
#         self.labels = npz["arr_2"]

    def load(self):
        self.model = joblib.load(self.args.model)

    def learn_coefs(self):
        alphas = [0.0001, 0.0003, 0.0005, 0.0007, 0.0009, 0.002, 0.004,
                  0.01, 0.05]
        l1_ratio = [.1, 0.25, .3]
        clf = ElasticNetCV(l1_ratio=l1_ratio, alphas=alphas)
        clf.fit(self.features, self.labels)
        logging.info('l1_ratio ' + str(clf.l1_ratio_))
        logging.info('alpha ' + str(clf.alpha_))
        

        coefs = clf.coef_
        coefs[-6] = 1
        logging.info('coefficients ' + str(coefs))
        return coefs

    def fit(self):
        self.model = RandomForestClassifier(random_state=0, n_estimators=3000,
                                            n_jobs=1)

        if self.args.enet:
            coefs = self.learn_coefs()
            
            self.zeros = [i for i,val in enumerate(coefs.tolist()) if val == 0]
            self.features = [[feature for i, feature in enumerate(features)
                              if i not in self.zeros]
                             for features in self.features]
            self.model.coefs = coefs
            self.features = numpy.array(self.features)

        self.model.fit(self.features, self.labels)
        self.model.name = self.feature_set_name
        self.model.version = self.feature_set_version

    def save(self):
        # numpy.savez(self.args.out, self.version, self.features, self.labels)
        joblib.dump(self.model, self.args.out, compress=9)

    def get_feature_importance(self):
        return self.model.feature_importances_

    def get_feature_names(self):
        return self.feature_module.Features().get_feature_names()

    def print_feature_importance(self):
        with open(self.args.out + "_importance.txt", 'w') as importance_file:
            feature_importance = self.get_feature_importance()
            feature_names = self.get_feature_names()

            
            if self.args.enet:
                feature_names = [val for i,val in enumerate(feature_names) if i not in self.zeros]

            for importance, feature_name in sorted(
                    zip(feature_importance, feature_names)):
                print >> importance_file, feature_name, importance

    def validate(self):
        self.generate(infiles=self.args.validate)

        # NOTE:copied from Jeff/Fatemeh old code
        logging.info("predicting probabilities")
        probs = self.model.predict_proba(self.features)
        voted = probs[:, 1]
        fpr, tpr, _ = roc_curve(self.labels, voted)
        roc_auc = auc(fpr, tpr)
        logging.info("AUC:%f" % roc_auc)

        fd = open(self.args.out + '_result.txt', 'w')
        for f, k, vf in zip(voted, self.keys, self.features):
            vf = " ".join(map(str, vf))
            k = " ".join(map(str, k))
            print >> fd, k + " " + str(f) + " " + vf

        fd.close()

        logging.info("plotting ROC curves")
        plt.plot(fpr, tpr, 'k--', label='ROC curve (area = %0.3f)' %
                 float(roc_auc))

        plt.title('ROC curve (area = %0.3f)' % float(roc_auc))
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.savefig(self.args.out + "_roc.png")

    def cross_validate(self):
        # NOTE:copied from Jeff/Fatemeh old code
        cv = cross_validation.StratifiedKFold(self.labels, n_folds=3)
        for i, (train, test) in enumerate(cv):
            cv_model = self.model.fit(self.features[train], self.labels[train])
            logging.info("predicting probabilities")
            probs = cv_model.predict_proba(self.features[test])

            voted = probs[:, 1]
            fpr, tpr, thresholds = roc_curve(self.labels[test], voted)
            roc_auc = auc(fpr, tpr)
            logging.info("AUC_%d:%f" % (i, roc_auc))

            fd = open(self.args.out + '_result_' + str(i) + '.txt', 'w')
            if self.args.model is None:
                for f, k in zip(voted, self.keys[test]):
                    print >> fd, k[0] + ' ' + k[1] + ' ' + k[2] + \
                        ' ' + k[3] + ' ' + k[4] + ' ' + k[5] + ' ' + str(f)

            else:
                # if the model is loaded from input arguments, then the keys
                # are not known
                for f in voted:
                    print >> fd, str(f)

            fd.close()

            logging.info("plotting ROC curves")
            plt.plot(
                fpr, tpr, 'k--', lw=1, label="Fold %i (AUC=%0.3f)" %
                (i + 1, float(roc_auc)))

        plt.legend(loc="lower right", numpoints=1,)
        plt.savefig(self.args.out + "_rocxval.png")

    # =================
    # Boxplot
    # =================

    def generate_boxplot(self):

        logging.info('Starting the plotting process')

        self.tot_features_name = self.__get_feature_names()

        self.features_vals_dict = self.__generate_feature_dict()
        self.top_features_name, self.top_features_map = self.__get_top_features()

        labels = self.__get_label_names(self.args.infiles)
        self.features_list_label = [
            self.__get_features_from_feature_db(
                self.features_vals_dict,
                label) for label in labels]

        self.boxplot_plot(labels)

    def __get_features_from_feature_db(self, feature_dict, label):
        features = []
        for key, value in feature_dict.iteritems():
            key = key.strip().split(';')
            label_feature_dict = value[1]
            if not label == label_feature_dict:
                continue
            features.append(eval(value[0]))
        return features

    def __get_label_names(self, reference_files):
        labels = set()
        for files in reference_files:
            file_stream = open(files.name)
            for line in file_stream:
                if line[0] == '#':
                    continue
                line = line.strip().split()
                labels.add(line[2])
        return list(labels)

    def __get_feature_names(self):
        if self.args.single:
            feature_set = features_single.Features()
        else:
            feature_set = features.Features()
        feature_names = feature_set.get_feature_names()
        return feature_names

    def __get_dict_bylabels(self, inputs, labels):
        out_dict = defaultdict(list)
        tfile = None
        nfile = None
        rfile = None

        for labelname in labels:
            for infilename in inputs:
                tfile = None
                nfile = None
                rfile = None
                file_stream = open(infilename.name, 'r')
                for line in file_stream:
                    if line[0] == '#':
                        line = line.strip().split()
                        if self.args.single:
                            if line[1] == "tumour" or line[1] == "normal":
                                tfile = nfile = line[2]
                        else:
                            if line[1] == 'tumour':
                                tfile = line[2]
                            elif line[1] == 'normal':
                                nfile = line[2]
                        if line[1] == 'reference':
                            rfile = line[2]
                    else:
                        l = line.strip().split()
                        if l[2] == labelname:
                            chromosome = l[0]
                            pos = l[1]
                            out_dict[labelname].append(
                                (tfile, nfile, rfile, chromosome, pos))
                file_stream.close()
        return out_dict

    def __generate_feature_dict(self):
        features_vals_dict = {}
        file_stream = open(self.args.out + '_feature_db_train.txt', 'r')
        for line in file_stream:
            l = line.strip().split('\t')
            key = l[0].split()[0]
            value = l[1]
            label = l[2]
            features_vals_dict[key] = (value, label)
        return features_vals_dict

    def __get_top_features(self):
        top_features_names = []
        file_stream = open(self.args.out + "_importance.txt", 'r')
        for line in file_stream:
            line = line.strip().split()
            top_features_names.append(line[0])
        top_features_names.reverse()

        top_features_names_dict = {
            v: i for i,
            v in enumerate(
                self.tot_features_name) for x in top_features_names if v == x}

        return top_features_names, top_features_names_dict

    def __get_features_from_dict(self, label, f_dict, ref_dict):
        features = []
        for tfile, nfile, rfile, chromosome, pos in ref_dict.get(label):
            key = ';'.join([rfile, nfile, tfile, chromosome, pos])
            try:
                features.append(eval(f_dict[key]))
            except KeyError:
                logging.error('error: cannot find key "%s"\n' % str(key))
        return features

    def boxplot_plot(self, labels):
        pdfout = PdfPages(self.args.out + '_plots_train.pdf')
        for i, v in enumerate(self.top_features_name):
            index = self.top_features_map[v]
            if not self.tot_features_name[index] == v:
                logging.error(
                    'the feature name and feature values don\'t match')

            fvalue = []

            for feature in self.features_list_label:
                fv = []
                for p in feature:
                    fv.append(p[index])
                fvalue.append(fv)

            fvalue_count = [len(x) for x in fvalue]

            fig = plt.figure()
            fig.text(0.90, 0.97, 'importance:' + str(i + 1),
                     rotation='horizontal', horizontalalignment='center',
                     verticalalignment='bottom', fontsize=8)

            xlabel_description = 'Features (Count of positions) (Outliers \
                                 above the plot, Outliers below the plot)'
            bplot = plt.boxplot(fvalue)

            fliers = []
            for i in xrange(len(bplot['boxes'])):
                fliers_above = len(bplot['fliers'][i * 2]._y)
                fliers_below = len(bplot['fliers'][i * 2 + 1]._y)
                fliers.append(str(fliers_above) + ',' + str(fliers_below))

            try:
                xlabel_names = ['%s(%s)\n(%s)' % (labels[i], y, fliers[i])
                                for i, y in enumerate(fvalue_count)]
            except:
                print '\nlabels', labels
                print '\ny', y
                print '\nfliers', fliers
                print '\ni', i
                print '\nfvalue_count', fvalue_count
                raise Exception()

            label_ylim_upper = None
            label_upper_cap = None
            label_lower_cap = None
            label_ylim_lower = None
            for i in xrange(len(bplot['boxes'])):
                try:
                    uppercap = bplot['caps'][i * 2]._y[0]
                    highestflier = max(bplot['fliers'][i * 2]._y)
                    if highestflier > uppercap * 100:
                        if not uppercap == 0:
                            if uppercap > label_upper_cap:
                                label_upper_cap = uppercap
                            if highestflier > label_ylim_upper:
                                label_ylim_upper = highestflier
                except:
                    # if unable to set the axis, continue without changing them
                    pass

                try:
                    lowercap = bplot['caps'][i * 2 + 1]._y[0]
                    lowestflier = min(bplot['fliers'][i * 2 + 1]._y)
                    if lowestflier > lowercap / 100:
                        if not lowercap == 0:
                            if lowercap < label_lower_cap:
                                label_lower_cap = lowercap
                            if lowestflier < label_ylim_lower:
                                label_ylim_lower = lowestflier
                except:
                    # if unable to set the axis, continue without changing them
                    pass

            if label_ylim_upper and label_upper_cap is not None:
                if label_ylim_upper > label_upper_cap * 100:
                    label_ylim_upper = int(label_upper_cap * 100)

            if not label_ylim_lower and label_lower_cap is not None:
                if label_ylim_lower < label_lower_cap / 100:
                    label_ylim_upper = int(label_lower_cap / 100)
            if label_ylim_upper is not None:
                prev_val = plt.ylim()
                if not label_ylim_upper == 0:
                    plt.ylim(prev_val[0], label_ylim_upper)
                    fig.text(0.02, 0.03, '*Upper limit (y-axis) rescaled, some\
                            datapoints are not shown',
                             horizontalalignment='left',
                             verticalalignment='top',
                             fontsize=6)
            if label_ylim_lower is not None:
                prev_val = plt.ylim()
                if not label_ylim_lower == 0:
                    plt.ylim(label_ylim_lower, prev_val[1])
                    fig.text(0.02, 0.02, '*Lower limit (y-axis) rescaled,\
                            some datapoints are not shown',
                             horizontalalignment='left',
                             verticalalignment='top',
                             fontsize=6)

            # just to show boundaries(if overlaps with axis)
            ylims = plt.ylim()
            if ylims[0] == 0:
                plt.ylim((ylims[0] - 0.05 * ylims[1]), (1.10 * ylims[1]))
            else:
                plt.ylim((0.90 * ylims[0]), (1.10 * ylims[1]))

            plt.title(v.replace('_', ' '))
            plt.ylabel('Distribution', fontsize=8)
            plt.xlabel(xlabel_description, fontsize=8)
            plt.xticks(range(1, len(xlabel_names) + 1), xlabel_names,
                       rotation=45, fontsize=8)
            plt.tight_layout()
            pdfout.savefig(fig)
        pdfout.close()
