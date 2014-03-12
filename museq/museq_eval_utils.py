'''
Created on Feb 11, 2014

@author: dgrewal
'''

import numpy
from sklearn import metrics
import logging
from collections import defaultdict
import pybamapi
import features
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import os

#====================================================
#Probability and ROC plots
#====================================================
class museq_plots(object):
    
    def __init__(self,args):
        self.args = args
            
        self.generated_reffiles = []
        
        #if the callfiles are vcf then generate their corresponding callfiles
        for reference_files in self.args.reference_files:
            for reference_file in reference_files.strip().split(','):
                if self.__is_vcf(reference_file):
                    outputname = self.args.out+''+reference_file.strip().split('/')[-1]+"_generated_referencefile"           
                    self.__generate_ref_file(outputname,reference_file)    
                    self.generated_reffiles.append(outputname)
                else:
                    self.generated_reffiles = self.args.reference_files
               
        #parse the text files generated into dictionary
        ref_dict = self.__parse_ref_files(self.generated_reffiles)
        
        self.y_true, self.y_score, self.vr, pos_accessed = self.__parse_input_vcf_files(ref_dict)
        
        self.__get_pos_neg_lists_counts()
        
        self.__get_pos_present_ref_absent_input(pos_accessed,ref_dict)
        
            
    def get_ref_files(self):
        return self.generated_reffiles
    
    # get the positions that are present in reference_file but not in input
    def __get_pos_present_ref_absent_input(self,pos_accessed,ref_dict):     
        file_stream = open(self.args.out+'present_in_ref_absent_in_input','w')
        keys = ref_dict.keys()
        listvals = [key for key in keys if key not in pos_accessed]
        for val in listvals:
            file_stream.write(str(val)+'\n')
        file_stream.close()  
         
    def __is_vcf(self,filename):
        file_stream = open(filename)
        for line in file_stream:
            line = line.strip().split('=')
            if line[0]=='##fileformat' and line[1].split('v')[0]=='VCF':
                return True
            else:
                return False
    
    def __generate_ref_file(self,outputname,reffile):
        self.neg_pos_call = []
        file_stream = open(reffile)
        output = []
        tfile = None
        rfile = None
        nfile = None
        for line in file_stream:
            if line[0] == '#':
                line = line.strip().split('=')
                if line[0].replace('#','').replace('=','') == 'tumour':
                    tfile = line[1]
                elif line[0].replace('#','').replace('=','') == 'normal':
                    nfile = line[1]
                elif line[0].replace('#','').replace('=','') == 'reference':
                    rfile = line[1] 
                continue  
            else:
                line = line.strip().split('\t')
                chromosome_id = line[0]
                position = line[1]
                info = line[7].strip().split(';')
                probability = info[0].split('=')[1]
                
                if float(probability) >= self.args.filter_threshold:
                    output.append(""+chromosome_id+" "+position+"  TRUE \n")
                else:
                    output.append(""+chromosome_id+" "+position+"  FALSE \n")
        
        output_file_stream = open(outputname,'w')    
        output_file_stream.write("# normal " + nfile +"\n")
        output_file_stream.write("# tumour " + tfile +"\n")
        output_file_stream.write("# reference " + rfile+"\n")
        for line in output:
            output_file_stream.write(line)
        output_file_stream.close()
        file_stream.close()
            
            
    def __get_pos_neg_lists_counts(self):
        self.positives_count = 0
        self.tp_count = 0
        self.fp_count = 0
        self.negatives_count = 0
        self.tn_count = 0
        self.fn_count= 0
        self.positives_list = []
        self.negatives_list = []
        self.positives_vr = []
        self.negatives_vr = []
        
        for (y_true,y_score,vr) in zip(self.y_true,self.y_score,self.vr):
            if y_true ==1:
                self.positives_count +=1
                self.positives_list.append(y_score)
                self.positives_vr.append(vr)
                if y_score >= self.args.filter_threshold:
                    self.tp_count +=1
                else:
                    self.fn_count +=1
            elif y_true ==0:
                self.negatives_count +=1
                self.negatives_list.append(y_score)
                self.negatives_vr.append(vr)
                if y_score<self.args.filter_threshold:
                    self.tn_count +=1
                else:
                    self.fp_count +=1
        try:
            logging.info("# of Positive labels = "+ str(self.positives_count))
            logging.info("# of Negative labels = "+ str(self.negatives_count))
            logging.info("TPR (TH= 0.9) = "+ 
                     str(float(self.tp_count)/(self.tp_count+self.fn_count)))
            logging.info("TNR (TH= 0.9) = "+ 
                     str(float(self.tn_count)/(self.tn_count+self.fp_count)))
        except ZeroDivisionError:
            logging.error("Either Positives count or Negatives count is zero")
    
    def __parse_ref_files(self,reffiles):   
        ref_dict = defaultdict(list)
        for cases in reffiles:
            for case in cases.strip().split(','):
                tfile = None
                nfile = None
                rfile = None
            
                stream_case = open(case)
                for line in stream_case:
                    l = line.strip().split()
                    if len(l) < 3:
                        logging.error('The call file should have 3 space delimited columns')
                        continue
                    
                    if l[0] == "#":
                        if l[1] == "tumour":
                            tfile = l[2]
    
                        elif l[1] == "normal":
                            nfile = l[2]
    
                        if l[1] == "reference":
                            rfile = l[2]
                        continue

                    ## check if required bam files/reference are specified in the training data
                    if not all((tfile, nfile, rfile)):
                        logging.warning("'%s' does not contain the required paths to bam/reference files" % case)
                        continue
            
                    chromosome = l[0]
                    position = l[1]
                    label = l[2]              
                    
                    ref_dict[(tfile, nfile, rfile,chromosome,position)].append(label)
                stream_case.close()
        return ref_dict    
            
        
    def __parse_input_vcf_files(self,ref_dict):
        self.input_vcf_dict = defaultdict(list)
        y_true = []
        y_score = []
        vr = []
        ta_tr_counter = 0 
        poslist_stream = open(self.args.out+'posList'+".txt",'w')
        tplist_stream = open(self.args.out+
                             'posList_TP_th%f'%self.args.filter_threshold+".txt",'w')
        fnlist_stream = open(self.args.out+
                             'posList_FN_th%f'%self.args.filter_threshold+".txt",'w')
        absent_reference_stream = open(self.args.out+"present_input_absent_reference.txt", 'a')
        pos_accessed = []
        
        for case in self.args.input_files:
            tfile = None
            nfile = None
            rfile = None
            stream_case = open(case)
            for line in stream_case:
                ## parse the line
                if line[0] == "#":
                    l = line.strip().split('=')
                    if l[0].replace('#','') == "tumour":
                        tfile = l[1]
    
                    elif l[0].replace('#','') == "normal":
                        nfile = l[1]
    
                    if l[0].replace('#','') == "reference":
                        rfile = l[1]
                    continue

                ## check if required bam files/reference are specified in the training data
                if not all((tfile, nfile, rfile)):
                    logging.error("'%s' does not contain the required paths to bam/reference files" % case)
                    continue
                
                l = line.strip().split()
                chromosome = l[0]
                pos = l[1]
                info = l[7].split(';')
                try:
                    pr = float(info[0].split('=')[1])
                    tr = float(info[1].split('=')[1])
                    ta = float(info[2].split('=')[1])
                except:
                    logging.error("invalid vcf file")
                    
                if ref_dict.has_key((tfile,nfile,rfile,chromosome,pos)):
                    ref = ref_dict.get((tfile,nfile,rfile,chromosome,pos))
                    if len(ref)>1:
                        logging.error('reference dictionary should have only one label per key')
                    ref = ref[0]
                    #flag the positions that have been accessed
                    pos_accessed.append((tfile,nfile,rfile,chromosome,pos))
                else:
                    absent_reference_stream.write(tfile+" "+nfile+" "+rfile+" "+chromosome+" "+pos+" Input_file: "+case+"\n")
                    continue
                        
                vr.append(float(ta)/(ta+tr))
                     
                if ref == 'FALSE':
                    y_true.append(0)
                elif ref == 'TRUE':
                    y_true.append(1)
                else:
                    logging.error('The call value is neither True nor False')
                y_score.append(pr)
                poslist_stream.write(tfile +'\n'+nfile+'\n'+rfile+'\n'+' '+chromosome+' '+
                                     pos+' '+ref+' '+str(pr)+'\n')
                if ref == 'TRUE':
                    if pr > self.args.filter_threshold:
                        tplist_stream.write(tfile +'\n'+nfile+'\n'+rfile+'\n'+' '+chromosome+' '+
                                            pos+' '+ref+' '+str(pr)+'\n')
                    if pr < self.args.filter_threshold:
                        fnlist_stream.write(tfile +'\n'+nfile+'\n'+rfile+'\n'+' '+chromosome+' '+
                                            pos+' '+ref+' '+str(pr)+'\n')  
            stream_case.close()            
        poslist_stream.close()
        tplist_stream.close()
        fnlist_stream.close()
        absent_reference_stream.close()               
        
        logging.info('Count data where TA/TR is less than or equal 0.04: ' 
                     + str(ta_tr_counter)) 
        #Ensure everything is valid
        if len(y_true) != len(y_score):
            logging.error('invalid values returned for y_true and y_score')
        logging.info('number of inputs : '+ str(len(y_true))+'\n')
        if (not 0 in y_true) or (not 1 in y_true):
            logging.error('Error: Missing a mix of either TRUE or FALSE in input')
        
        y_true = numpy.array(y_true)
        y_score = numpy.array(y_score)
        
        return y_true,y_score,vr,pos_accessed 
                
         
    def __write_info_log(self,fpr,tpr,thresholds,roc_auc):
        roc_info_stream = open(self.args.out+
                               self.args.model_name+'_roc_info.txt','w')
        roc_info_stream.write('fpr'+' '+'tpr'+' '+'threshold'+'\n')
        
        for (fpr_val,tpr_val,thresholds_val) in zip(fpr,tpr,thresholds):
            roc_info_stream.write(str(fpr_val)+' '+str(tpr_val)+' '+
                                  str(thresholds_val)+'\n')
        roc_info_stream.close()
        
        logging.info('fpr length'+ str(len(fpr)))
        logging.info('tpr length'+ str(len(tpr)))
        try:
            logging.info("Area under the ROC curve :"+str(roc_auc)+'\n')
        except ValueError:
            logging.error("ERROR: fpr or tpr Array contains NaN or infinity")
        
    def __roc_plot(self,fig,subplot_pos):
        fpr,tpr,thresholds = metrics.roc_curve(self.y_true,self.y_score)
        
        try:
            roc_auc = metrics.auc(fpr,tpr)
            self.__write_info_log(fpr,tpr,thresholds,roc_auc)
        except ValueError:
            logging.error('ERROR: fpr or tpr array contains NaN or infinity')
            return
        
        ax=fig.add_subplot(subplot_pos[0],subplot_pos[1],subplot_pos[2])
        ax.plot(fpr[0:50],tpr[0:50],'k--', 
                    label='ROC curve (area = %0.3f)' %float(roc_auc))
        
        ax.set_title('ROC curve (area = %0.3f)' %float(roc_auc),fontsize = 8)
        ax.set_xlabel('False Positive Rate',fontsize = 8)
        ax.set_ylabel('True Positive Rate',fontsize = 8)
        ax.yaxis.set_tick_params(labelsize=8)
        ax.xaxis.set_tick_params(labelsize=8)
                         
    def __probability_plot(self,probabilities,fig,subplot_pos):
        
        probability_values = list(set(probabilities))
        probability_values.sort()
        
        probability_counts = [probabilities.count(p) for p in probability_values]
        
        ax=fig.add_subplot(subplot_pos[0],subplot_pos[1],subplot_pos[2])
        plot_title = ""+self.args.model_name+" calls, "+str(len(probabilities))+" positions"
        
        ax.plot(probability_values,probability_counts)
        ax.set_xlabel('probability',fontsize = 8)
        ax.set_ylabel('frequency',fontsize = 8)
        ax.set_title(plot_title,fontsize = 8)
        ax.yaxis.set_tick_params(labelsize=8)
        ax.xaxis.set_tick_params(labelsize=8)
    
    def __scatter_plot(self,fig,subplot_pos):
        ax=fig.add_subplot(subplot_pos[0],subplot_pos[1],subplot_pos[2])
        ax.scatter(self.positives_vr,self.positives_list,color = 'b', 
                       label = 'Positive Labels',alpha = 0.2)
        
        ax.scatter(self.negatives_vr,self.negatives_list,color = 'r', 
                       label = 'Negative Labels',alpha = 0.2)
        
        ax.set_xlabel('Variant Ratio',fontsize = 8)
        ax.set_ylabel('Probability',fontsize = 8)
        ax.yaxis.set_tick_params(labelsize=8)
        ax.xaxis.set_tick_params(labelsize=8)
    
    def __feature_importance_bar_plot(self,fig,subplot_pos):
        features=[]
        importance=[]
        file_stream = open(self.args.ranked_features,'r')
        for line in file_stream:
            line = line.strip().split()
            features.append(line[0].replace('_',' '))
            importance.append(float(line[1]) )
            
        ax = fig.add_subplot(subplot_pos[0],subplot_pos[1],subplot_pos[2])
        ax.bar([x for x in xrange(len(importance))],importance)
        ax.set_xlabel('Features',fontsize=8)
        ax.set_xticks(numpy.arange(len(importance)))
        ax.set_xticklabels(features,rotation = 90,fontsize=6)
        ax.set_ylabel('Importance',fontsize = 8)
        ax.yaxis.set_tick_params(labelsize=8)
        ax.xaxis.set_tick_params(labelsize=6)
        
    def generate_plots(self):
        fig = pyplot.figure(figsize=(15, 15))
        self.__probability_plot(self.positives_list,fig,(2,3,1))
        self.__probability_plot(self.negatives_list,fig,(2,3,2))
        self.__scatter_plot(fig,(2,3,3))
        self.__roc_plot(fig,(2,3,4))
        self.__feature_importance_bar_plot(fig,(2,3,5))
        fig.set_dpi(150)
        pyplot.tight_layout()
        return fig
    
#====================================
# Box Plots
#====================================
class box_plots():
    def __init__(self,args,reffiles):
        self.args = args
        
        if not self.args.separate_plots:
            fileargs = ','.join(reffiles)
            reffiles=[]
            reffiles.append(fileargs)
            
        model_impfile = self.args.ranked_features.strip().split('/')[-1]
        if not self.args.model_name in model_impfile:
            logging.error('The importance file doesn\'t match with model provided')
        
        if self.args.boxplot_labels != None:
            if len(reffiles) != len(self.args.boxplot_labels):
                logging.error('The length of reference files and labels do not match')
        else:
            self.args.boxplot_labels = [('input '+str(x) ) for x in range(len(reffiles))]
        
        self.tot_features_name = self.__get_feature_names(reffiles)
        
        self.input_files = []
        for input_args in reffiles:
            self.input_files.append(input_args.split(','))
            
        inputs = []
        for input_args in reffiles:
            for files in input_args.split(','):
                inputs.append(files)
                
        #generate files filtered by labels
        self.finalargs = self.__update_pos_files_bylabels(reffiles,inputs)
            
        #filter by correctness
        if self.args.input_files:
            self.vcf_dict = self.__get_vcf_dict()
            self.posneg_files,self.posneg_labels = self.__get_posneg_files(inputs)
                
        #generate feature_db 
        if self.args.feature_db == None:
            outputname = self.args.out+"feature_db"
            self.__generate_feature_db(outputname,reffiles)
            self.args.feature_db = outputname
        
        self.features_vals_dict = self.__get_features_dict()
        self.top_features_name, self.top_features_map = self.__get_top_features()
        self.features_list_label = [self.__get_features_from_dict(infile, self.features_vals_dict) for infile in self.finalargs]
        self.features_list_normal = [self.__get_features_from_dict(infile, self.features_vals_dict) for infile in self.input_files]
        if self.args.input_files:
            self.features_list_posneg = [self.__get_features_from_dict(infile, self.features_vals_dict) for infile in self.posneg_files]        

    def __update_pos_files_bylabels(self,reffiles,inputs):
        outputargs = []
        finalargs = []
                
        self.labels = self.__get_labels(inputs)
        for label in self.labels:
            for input_file in inputs:
                output = []
                file_stream = open(input_file)
                for line in file_stream:
                    if line[0]=='#':
                        output.append(line)
                        continue
                    l=line.strip().split()
                    if l[2]==label:
                        output.append(line)
                file_stream.close()
                outfile_stream = open(self.args.out+''+input_file.strip().split('/')[-1]+'_'+label,'w')
                for l in output:
                    outfile_stream.write(l)
                outfile_stream.close()
                outputargs.append(self.args.out+''+input_file.strip().split('/')[-1]+'_'+label)
            finalargs.append(outputargs)
            outputargs = []
        return finalargs
    
    def __get_feature_names(self,reffiles):
        tfile= None
        nfile=None
        rfile=None
        feature_names = None
        for infiles in reffiles:
            for infile in infiles.strip().split(','):
                infile_stream = open(infile)
                for line in infile_stream:
                    line=line.strip().split()
                    if line[0]=='#':
                        if line[1]=='tumour':
                            tfile = line[2]
                        if line[1]=='normal':
                            nfile = line[2]
                        if line[1]=='reference':
                            rfile = line[2]
                    else:
                        chromosome = line[0]
                        position = line[1]
                        t_bam = pybamapi.Bam(bam=tfile, reference=rfile, coverage=1)
                        n_bam = pybamapi.Bam(bam=nfile, reference=rfile, coverage=1)
                        tt = t_bam.get_tuple(chromosome, int(position) )
                        nt = n_bam.get_tuple(chromosome, int(position) )
                        chromosome_id = tt[-1]
                        rt = t_bam.get_reference_tuple(chromosome_id, position)
                        feature_set = features.Features(tt, nt, rt)
                        feature_names = feature_set.get_feature_names()
                        break
        return feature_names
                        
    def __get_posneg_files(self,inputs):
        tp_files=[]
        fp_files=[]
        tn_files=[]
        fn_files=[]
        for files in inputs:
            tfile = None
            rfile = None
            nfile = None
            tp=[]
            fp=[]
            tn=[]
            fn=[]
            tp_files.append(self.args.out+''+files.strip().split('/')[-1]+'_tp')
            fp_files.append(self.args.out+''+files.strip().split('/')[-1]+'_fp')
            tn_files.append(self.args.out+''+files.strip().split('/')[-1]+'_tn')
            fn_files.append(self.args.out+''+files.strip().split('/')[-1]+'_fn')
            
            file_stream = open(files)
            for line in file_stream:
                if line[0] == '#':
                    line = line.strip().split()
                    if line[1]=='tumour':
                        tfile = line[2]
                    if line[1]=='normal':
                        nfile = line[2]
                    if line[1] == 'reference':
                        rfile = line[2]
                else:
                    if not all((tfile,nfile,rfile)):
                        logging.error('Couldn\'t retreive paths for bams ')
                    line = line.strip().split()
                    chromosome = line[0]
                    position = line[1]
                    label = None
                    if (line[2] =='SOMATIC' or line[2] == 'TRUE'):
                        label = 'TRUE'
                    else:
                        label = 'FALSE'
                    if self.vcf_dict.has_key((tfile,nfile,rfile,chromosome,position)):
                        pr = self.vcf_dict[(tfile,nfile,rfile,chromosome,position)]
                        
                        if label == 'TRUE':
                            if pr[0]>=self.args.filter_threshold:
                                tp.append(chromosome+' '+position+' TP')
                            else:
                                fn.append(chromosome+' '+position+' FN')
                        else:
                            if pr[0]<self.args.filter_threshold:
                                tn.append(chromosome+' '+position+' TN')
                            else:
                                fp.append(chromosome+' '+position+' FP')
                        
            tp_file_stream = open(self.args.out+''+files.strip().split('/')[-1]+'_tp','w')
            tp_file_stream.write('# tumour '+tfile+'\n')
            tp_file_stream.write('# normal '+nfile+'\n')
            tp_file_stream.write('# reference '+rfile+'\n')
            for line in tp:
                tp_file_stream.write(line+'\n')
                    
            fp_file_stream = open(self.args.out+''+files.strip().split('/')[-1]+'_fp','w')
            fp_file_stream.write('# tumour '+tfile+'\n')
            fp_file_stream.write('# normal '+nfile+'\n')
            fp_file_stream.write('# reference '+rfile+'\n')
            for line in fp:
                fp_file_stream.write(line+'\n')
                    
            tn_file_stream = open(self.args.out+''+files.strip().split('/')[-1]+'_tn','w')
            tn_file_stream.write('# tumour '+tfile+'\n')
            tn_file_stream.write('# normal '+nfile+'\n')
            tn_file_stream.write('# reference '+rfile+'\n')
            for line in tn:
                tn_file_stream.write(line+'\n')
                   
            fn_file_stream = open(self.args.out+''+files.strip().split('/')[-1]+'_fn','w')  
            fn_file_stream.write('# tumour '+tfile+'\n')
            fn_file_stream.write('# normal '+nfile+'\n')
            fn_file_stream.write('# reference '+rfile+'\n')
            for line in fn:
                fn_file_stream.write(line+'\n')
                    
        return [tp_files,fn_files,tn_files,fp_files],['True Positives','False Negatives','True Negatives','False Positives']  
                     
    def __get_vcf_dict(self):
        output = defaultdict(list)
        for vcffile in self.args.input_files:
            tfile = None
            nfile = None
            rfile = None
            file_stream = open(vcffile)
            for line in file_stream:
                ## parse the line
                if line[0] == "#":
                    l = line.strip().split('=')
                    if l[0].replace('#','') == "tumour":
                        tfile = l[1]
    
                    elif l[0].replace('#','') == "normal":
                        nfile = l[1]
    
                    if l[0].replace('#','') == "reference":
                        rfile = l[1]
                    continue

                ## check if required bam files/reference are specified in the training data
                if not all((tfile, nfile, rfile)):
                    logging.warning("'%s' does not contain the required paths to bam/reference files" % vcffile)
                    continue
                
                l = line.strip().split()
                chromosome = l[0]
                pos = l[1]
                info = l[7].split(';')
                pr = float(info[0].split('=')[1])
                output[(tfile,nfile,rfile,chromosome,pos)].append(pr)
        return output                      
            
    def __get_labels(self,inputs):
        labels =set()
        for files in inputs:
            file_stream = open(files)
            for line in file_stream:
                if line[0] == '#':
                    continue
                line = line.strip().split()
                labels.add(line[2])
        return list(labels)
                    
    def remove_temp_files(self):
        try:
            for files in self.posneg_files:
                for filenames in files:
                    os.remove(filenames)
            for files in self.finalargs:
                for filenames in files:
                    os.remove(filenames)
        except:
            logging.error('couldn\'t delete temp files')
        
    def __generate_feature_db(self,outputname,reffiles):
        missing_positions = self.__get_missing_positions(reffiles)
        self.__write_feature_db(missing_positions,outputname)   
        
    def __get_missing_positions(self,reffiles):
        missing_positions = []
        for infile in reffiles:
            infile = infile.strip().split(',')
            for files in infile:
                file_stream = open(files,'r')
                for line in file_stream:
                    l = line.strip().split()
                    if len(l) < 3:
                        logging.error('The call file should have 3 space delimited columns')
                        continue
                    ## parse the line
                    if l[0] == "#":
                        if l[1] == "tumour":
                            tfile = l[2]
    
                        elif l[1] == "normal":
                            nfile = l[2]
    
                        if l[1] == "reference":
                            rfile = l[2]
                        continue

                    ## check if required bam files/reference are specified in the training data
                    if not all((tfile, nfile, rfile)):
                        logging.error("'%s' does not contain the required paths to bam/reference files" % infile)
                        continue
            
                    chromosome = l[0]
                    pos = l[1]
                    missing_positions.append((chromosome,pos,tfile,nfile,rfile))
        return missing_positions
    
    def __write_feature_db(self,missing_positions,outputname):
        feature, keys = self.__extract_features(missing_positions)
        file_stream = open(outputname,'w')
        for i,_ in enumerate(keys):
            
            k = ';'.join(keys[i])
            f = feature[i].tolist()
            file_stream.write(k+'\t'+str(f)+'\n')
        file_stream.close()
        
    def __get_features_dict(self):
        features_vals_dict = {}
        file_stream = open(self.args.feature_db,'r')
        for line in file_stream:
            key = line.strip().split('\t')[0].split()[0]
            value = line.strip().split('\t')[1]
            features_vals_dict[key] = value
        return features_vals_dict    
        
    def __get_top_features(self):
        top_features_names = []
        file_stream = open(self.args.ranked_features,'r')
        for line in file_stream:
            line = line.strip().split()
            top_features_names.append(line[0])
        top_features_names.reverse()
        if self.args.top_features:
            top_features_names = top_features_names[:self.args.top_features]
        
        top_features_names_dict =  {v:i for i, v in enumerate(self.tot_features_name) for x in top_features_names if v == x}
        
        return top_features_names,top_features_names_dict
        
    def __get_features_from_dict(self,infile,f_dict):
        features = []
        for files in infile:
            file_stream = open(files,'r')
            for line in file_stream:
                l = line.strip().split()
                if len(l) < 3:
                    logging.error('The call file should have 3 space delimited columns')
                    continue
                ## parse the line
                if l[0] == "#":
                    if l[1] == "tumour":
                        tfile = l[2]
    
                    elif l[1] == "normal":
                        nfile = l[2]
    
                    if l[1] == "reference":
                        rfile = l[2]
                    continue

                ## check if required bam files/reference are specified in the training data
                if not all((tfile, nfile, rfile)):
                    logging.warning("'%s' does not contain the required paths to bam/reference files" % infile)
                    continue
            
                chromosome = l[0]
                pos = l[1]
                key = ';'.join([rfile,nfile,tfile,chromosome,pos])
                try:
                    features.append(eval(f_dict[key]))
                except KeyError:
                    logging.error('error: cannot find key "%s"\n' % str(key))
        return features 
       
       
    def __extract_features(self,missing_positions):
        data = defaultdict(list)
        contamination = (float(30), float(30), float(70), float(0))
        
        for chromosome,pos,tfile,nfile,rfile in missing_positions:
            data[(tfile,nfile,rfile)].append((chromosome,int(pos),contamination))
        
        features_buffer = []
        keys_buffer    = []
    
        for tfile, nfile, rfile in data.keys():
            logging.info("reading from tumour:"+ tfile)
            logging.info("reading from normal:"+ nfile)
            t_bam = pybamapi.Bam(bam=tfile, reference=rfile, coverage=1)
            n_bam = pybamapi.Bam(bam=nfile, reference=rfile, coverage=1)
        
            for chromosome, position, _ in data[(tfile, nfile, rfile)]:
                tt = t_bam.get_tuple(chromosome, position)
                nt = n_bam.get_tuple(chromosome, position)
                if tt == None or nt == None:
                    logging.warning('Could not retreive the tuples from bams')
                    continue
                chromosome_id = tt[-1]
                rt = t_bam.get_reference_tuple(chromosome_id, position)
                if not all([tt, nt, rt]):
                    logging.error('ERROR: Could not retreive all three tuples')
                    continue
            
                ## calculate features
                feature_set = features.Features(tt, nt, rt)
                temp_features = feature_set.get_features()
            
                features_buffer.append(temp_features)
                keys_buffer.append((rfile, nfile, tfile, chromosome, position))
        
        features_buffer = numpy.array(features_buffer)
        keys_buffer = numpy.array(keys_buffer)
        return features_buffer, keys_buffer


    def boxplot_plot(self):
        plots = []
        for i,v in enumerate(self.top_features_name):
            index = self.top_features_map[v]
            if not self.tot_features_name[index] == v:
                logging.error('the feature name and feature values don\'t match')
                
            xlabel_description = 'Features (Count of positions) (Outliers above the plot, Outliers below the plot)'
            
            f = pyplot.figure(figsize=(15, 15))
            f.set_dpi(150)
            gs1 = gridspec.GridSpec(1, 3)           
            f.text(0.90, 0.98, 'importance:'+str(i+1), rotation='horizontal',horizontalalignment='center', verticalalignment='bottom',fontsize = 8)
            
            #normal
            fnvalue = []
            for feature in self.features_list_normal:
                fnv = []
                for p in feature:
                    fnv.append(p[index])
                fnvalue.append(fnv)
            fnvalue_count = [len(x) for x in fnvalue]
            ax1 = f.add_subplot(gs1[0])
            bplot = ax1.boxplot(fnvalue)
            if ax1.axis()[3] >40000:
                ax1.set_yscale('symlog')
            #get the outliers
            ax1_fliers = []
            for i in xrange(len(bplot['boxes'])):
                fliers_above = len(bplot['fliers'][i*2]._y)
                fliers_below = len(bplot['fliers'][i*2+1]._y)
                ax1_fliers.append( str(fliers_above)+','+str(fliers_below) )
            
            #if some value in vector is [],then no boxplot.so we need to adjust the value to appropriate plot    
                if not len(ax1_fliers) == len(fnvalue):
                    indices = [i for i, j in enumerate(fnvalue) if j == [] ]
                    for val in indices:
                        ax1_fliers.insert(val, [])
                                    
            normal_xlabel_names = ['%s(%s)\n(%s)' % (self.args.boxplot_labels[i], y,ax1_fliers[i]) for i,y in enumerate(fnvalue_count)]
                         
            normal_ylim_upper = None
            normal_upper_cap =None
            normal_lower_cap = None
            normal_ylim_lower = None   
            for i in xrange(len(bplot['boxes'])):
                try:
                    uppercap = bplot['caps'][i*2]._y[0]
                    highestflier = max(bplot['fliers'][i*2]._y)
                    if highestflier > uppercap*100:
                        if not uppercap == 0:
                            if uppercap>normal_upper_cap:
                                normal_upper_cap = uppercap
                            if highestflier>normal_ylim_upper:
                                normal_ylim_upper = highestflier
                except:
                    #if unable to set the axis, continue without changing them
                    pass
             
                try:
                    lowercap = bplot['caps'][i*2+1]._y[0]
                    lowestflier = min(bplot['fliers'][i*2+1]._y)
                    if lowestflier > lowercap/100:
                        if not lowercap == 0:
                            if lowercap < normal_lower_cap:
                                normal_lower_cap = lowercap
                            if lowestflier < normal_ylim_lower:
                                normal_ylim_lower = lowestflier  
                except:
                    #if unable to set the axis, continue without changing them
                    pass
            
            if normal_ylim_upper and normal_upper_cap is not None:
                if normal_ylim_upper > normal_upper_cap*100:
                    normal_ylim_upper = int(normal_upper_cap*100)
            
            if normal_ylim_lower and normal_lower_cap is not None:
                if normal_ylim_lower < normal_lower_cap/100:
                    normal_ylim_upper = int(normal_lower_cap/100)
         
            ax1.tick_params(axis='x', labelsize=8)
            ax1.set_ylabel('Distribution',fontsize = 8)
            ax1.set_xticklabels( normal_xlabel_names, rotation=45, fontsize=8)
            ax1.set_title('Feature Distributions for whole data',fontsize = 8)
            ax1.yaxis.set_tick_params(labelsize=8)
                                      
            #labels
            flvalue = []
            for feature in self.features_list_label:
                flv = []
                for p in feature:
                    flv.append(p[index])
                flvalue.append(flv)
                
            flvalue_count = [len(x) for x in flvalue]
            
            ax2 = f.add_subplot(gs1[1],sharey=ax1)
            bplot = ax2.boxplot(flvalue)
            #get the outliers
            ax2_fliers = []
            for i in xrange(len(bplot['boxes'])):
                fliers_above = len(bplot['fliers'][i*2].get_data()[1])
                fliers_below = len(bplot['fliers'][i*2+1].get_data()[1])
                ax2_fliers.append(str(fliers_above)+','+str(fliers_below))
            #if some value in vector is [],then no boxplot.so we need to adjust the value to appropriate plot    
                if not len(ax2_fliers) == len(flvalue):
                    indices = [i for i, j in enumerate(flvalue) if j == [] ]
                    for val in indices:
                        ax2_fliers.insert(val, [])
                                    
            label_xlabel_names = ['%s(%s)\n(%s)' % (self.labels[i], y,ax2_fliers[i]) for i,y in enumerate(flvalue_count)]
                
            label_ylim_upper = None
            label_upper_cap =None
            label_lower_cap = None
            label_ylim_lower = None   
            for i in xrange(len(bplot['boxes'])):
                try:
                    uppercap = bplot['caps'][i*2]._y[0]
                    highestflier = max(bplot['fliers'][i*2]._y)
                    if highestflier > uppercap*100:
                        if not uppercap == 0:
                            if uppercap>label_upper_cap:
                                label_upper_cap = uppercap
                            if highestflier>label_ylim_upper:
                                label_ylim_upper = highestflier
                except:
                    #if unable to set the axis, continue without changing them
                    pass
             
                try:
                    lowercap = bplot['caps'][i*2+1]._y[0]
                    lowestflier = min(bplot['fliers'][i*2+1]._y)
                    if lowestflier > lowercap/100:
                        if not lowercap == 0:
                            if lowercap < label_lower_cap:
                                label_lower_cap = lowercap
                            if lowestflier < label_ylim_lower:
                                label_ylim_lower = lowestflier  
                except:
                    #if unable to set the axis, continue without changing them
                    pass
            
            if label_ylim_upper and label_upper_cap is not None:
                if label_ylim_upper > label_upper_cap*100:
                    label_ylim_upper = int(label_upper_cap*100)
            
            if not label_ylim_lower and label_lower_cap is not None:
                if label_ylim_lower < label_lower_cap/100:
                    label_ylim_upper = int(label_lower_cap/100)
                              
            ax2.set_xticklabels( label_xlabel_names, rotation=45, fontsize=8)
            ax2.set_title('Feature Distributions by labels',fontsize = 8)
            ax2.yaxis.set_tick_params(labelsize=8)
            
            posneg_ylim_upper = None
            posneg_upper_cap =None
            posneg_lower_cap = None
            posneg_ylim_lower = None
            ax3_fliers = []
            #posneg
            if self.args.input_files:
                fpnvalue = []
                for feature in self.features_list_posneg:
                    fpnv = []
                    for p in feature:
                        fpnv.append(p[index])
                    fpnvalue.append(fpnv)
                fpnvalue_count = [len(x) for x in fpnvalue]
                
                ax3 = f.add_subplot(gs1[2],sharey=ax1)
                bplot = ax3.boxplot(fpnvalue)
                
                #get the outliers
                for i in xrange(len(bplot['boxes'])):
                    fliers_above = len(bplot['fliers'][i*2].get_data()[1])
                    fliers_below = len(bplot['fliers'][i*2+1].get_data()[1])
                    ax3_fliers.append(str(fliers_above)+','+str(fliers_below))
                
                #if some value in vector is [],then no boxplot.so we need to adjust the value to appropriate plot    
                if not len(ax3_fliers) == len(fpnvalue):
                    indices = [i for i, j in enumerate(fpnvalue) if j == [] ]
                    for val in indices:
                        ax3_fliers.insert(val, [])
                    
                posneg_xlabel_names = ['%s(%s)\n(%s)' % (self.posneg_labels[i], y,ax3_fliers[i]) for i,y in enumerate(fpnvalue_count)]
                
                for i in xrange(len(bplot['boxes'])):
                    try:
                        uppercap = bplot['caps'][i*2]._y[0]
                        highestflier = max(bplot['fliers'][i*2]._y)
                        if highestflier > uppercap*100:
                            if not uppercap == 0:
                                if uppercap>posneg_upper_cap:
                                    posneg_upper_cap = uppercap
                                if highestflier>posneg_ylim_upper:
                                    posneg_ylim_upper = highestflier
                    except:
                        #if unable to set the axis, continue without changing them
                        pass
                    
                    try:
                        lowercap = bplot['caps'][i*2+1]._y[0]
                        lowestflier = min(bplot['fliers'][i*2+1]._y)
                        if lowestflier > lowercap/100:
                            if not lowercap == 0:
                                if lowercap < posneg_lower_cap:
                                    posneg_lower_cap = lowercap
                                if lowestflier < posneg_ylim_lower:
                                    posneg_ylim_lower = lowestflier  
                    except:
                        #if unable to set the axis, continue without changing them
                        pass
            
                if posneg_ylim_upper and posneg_upper_cap is not None:
                    if posneg_ylim_upper > posneg_upper_cap*100:
                        posneg_ylim_upper = int(posneg_upper_cap*100)
               
                if posneg_ylim_lower and posneg_lower_cap is not None: 
                    if posneg_ylim_lower < posneg_lower_cap/100:
                        posneg_ylim_upper = int(posneg_lower_cap/100)
                    
                ax3.set_xticklabels( posneg_xlabel_names, rotation=70, fontsize=8)
                ax3.set_title('Feature Distributions by T/F P/N',fontsize = 8 )
                ax3.yaxis.set_tick_params(labelsize=8)
            
            if not normal_ylim_lower==None or label_ylim_lower==None or posneg_ylim_lower==None:
                prev_val = pyplot.ylim()
                lower_val = min(normal_ylim_lower, label_ylim_lower, posneg_ylim_lower)
                if not lower_val == 0:
                    pyplot.ylim(lower_val,prev_val[1])
                    f.text(0.02, 0.01,'*Lower limit (y-axis) rescaled, some datapoints are not shown', horizontalalignment='left',verticalalignment='top',fontsize=8)
                
            if normal_ylim_upper or label_ylim_upper or posneg_ylim_upper is not None:
                prev_val = pyplot.ylim()
                upper_val = max(normal_ylim_upper, label_ylim_upper, posneg_ylim_upper)
                if not upper_val == 0:
                    pyplot.ylim(prev_val[0], upper_val)
                    f.text(0.02, 0.02,'*Upper limit (y-axis) rescaled, some datapoints are not shown', horizontalalignment='left',verticalalignment='top',fontsize=8)
            
            #just +-1 to show boundaries(if overlaps with axis)
            ylims = pyplot.ylim() 
            pyplot.ylim(ylims[0]-1, ylims[1]+1)    
            
            f.text(0.90, 0.97,'Threshold ='+str(self.args.filter_threshold), horizontalalignment='center',verticalalignment='top',fontsize=8)
            f.text(0.50, 0.02, xlabel_description, rotation='horizontal',horizontalalignment='center', verticalalignment='top',fontsize = 10)
            f.text(0.50, 0.97, v+' ('+self.args.model_name+')', rotation='horizontal',horizontalalignment='center', verticalalignment='bottom',fontsize = 10)
            
            gs1.tight_layout(f, rect=[0, 0, 1, .96], h_pad=0.5)
            plots.append(f)
        return plots    