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
import features_deep
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import os
from scipy.integrate import trapz,simps

mutationSeq_version="4.3.2"
#====================================================
#Probability and ROC plots
#====================================================
class museq_plots(object):
    
    def __init__(self,args):
        self.args = args
            
        self.generated_reffiles = []
        
        #if the reference files are vcf then generate their corresponding reference files
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
        
        #y_true list with 0 for False and 1 for True. y_score contains corresponding probability, vr = ta/ta+tr
        self.y_true, self.y_score, self.vr, pos_accessed = self.__parse_input_vcf_files(ref_dict)
        
        #lists with T/F P/N positions and corresponding counts
        self.__get_pos_neg_lists_counts()
        
        #get the positions present in reference but missing in input vcf
        self.__get_pos_present_ref_absent_input(pos_accessed,ref_dict)
        
            
    def get_ref_files(self):
        return self.generated_reffiles
    
    # get the positions that are present in reference_file but not in input
    def __get_pos_present_ref_absent_input(self,pos_accessed,ref_dict):     
        file_stream = open(self.args.out+'present_in_ref_absent_in_input.txt','w')
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
                    print line
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
                     
                if ref in self.args.positive_labels:
                    y_true.append(1)
                else:
                    y_true.append(0)
                    
                y_score.append(pr)
                poslist_stream.write(tfile +'\n'+nfile+'\n'+rfile+'\n'+' '+chromosome+' '+
                                     pos+' '+ref+' '+str(pr)+'\n')
                if ref in self.args.positive_labels:
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
        
        #Get the index and then area for vals 0-0.1
        prev_val = 0
        index = 0
        roc_auc = 0
        for i,val in enumerate(fpr.tolist()):
            if val > 0.1 and prev_val < 0.1:
                index = i
                break
            prev_val = val
        fpr_10 = fpr.tolist()[:index+1]
        tpr_10 = tpr.tolist()[:index+1]
        
        try:
            roc_auc = metrics.auc(fpr,tpr)
            roc_auc_10 = metrics.auc(fpr_10,tpr_10)
            self.__write_info_log(fpr,tpr,thresholds,roc_auc)
        except Exception, e:
            print roc_auc
            print e
            logging.error('ERROR: fpr or tpr array contains NaN or infinity')
            return
        
        ax=fig.add_subplot(subplot_pos[0],subplot_pos[1],subplot_pos[2])
        ax.plot(fpr[0:50],tpr[0:50],'k--', 
                    label='ROC curve (area = %0.3f)' %float(roc_auc))
        
        ax.annotate(str(round(roc_auc_10,4))+'\n('+str(round(val,4))+')',(fpr[index+1],tpr[index+1]),
                xytext=(-20, 30), textcoords='offset points',
                arrowprops=dict(facecolor='gray', shrink=0.1),
                horizontalalignment='right', verticalalignment='bottom',fontsize = '4' )
        
        ax.set_title('ROC curve (area = %0.3f)' %float(roc_auc),fontsize = 8)
        ax.set_xlabel('False Positive Rate',fontsize = 8)
        ax.set_ylabel('True Positive Rate',fontsize = 8)
        ax.yaxis.set_tick_params(labelsize=8)
        ax.xaxis.set_tick_params(labelsize=8)
                         
    def __probability_plot(self,probabilities,fig,title,subplot_pos):
        file_stream = open(self.args.out+title+'_prob_plot_auc','w')
        probability_values = list(set(probabilities))
        probability_values.sort()
        
        probability_counts = [probabilities.count(p) for p in probability_values]
        
        ax=fig.add_subplot(subplot_pos[0],subplot_pos[1],subplot_pos[2])
        plot_title = ""+self.args.model_name+" "+title+" calls, "+str(len(probabilities))+" positions"
        ax.plot(probability_values,probability_counts)

        for i,vals in enumerate(probability_values):
            index = probability_values.index(vals)
            if index == 0:
                continue
            diffx = ax.get_xticks()[1]-ax.get_xticks()[0]
            auc = simps(probability_counts[:index],dx=diffx)   
            file_stream.write(str(vals)+'\t'+str(auc)+'\n')
            if len(probability_values)>1000:
                if i%100:
                    continue
            if len(probability_values)>100:
                if i%20:
                    continue
            elif len(probability_values)>20:
                if i%5:
                    continue
            elif len(probability_values)>5:
                if i%3:
                    continue
                
            index = probability_values.index(vals)
            if index == 0:
                continue
            diffx = ax.get_xticks()[1]-ax.get_xticks()[0]
            auc = simps(probability_counts[:index],dx=diffx)   
            auc=round(auc,3)         
            ax.annotate(str(auc)+'\n('+str(vals)+')',(vals,probability_counts[index]),
                xytext=(-20, 30), textcoords='offset points',
                arrowprops=dict(facecolor='gray', shrink=0.1),
                horizontalalignment='right', verticalalignment='bottom',fontsize = '4' )
                
        ax.set_xlabel('probability',fontsize = 8)
        ax.set_ylabel('frequency',fontsize = 8)
        ax.set_title(plot_title,fontsize = 8)
        ax.yaxis.set_tick_params(labelsize=8)
        ax.xaxis.set_tick_params(labelsize=8)
        file_stream.close()
        
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
        self.__probability_plot(self.positives_list,fig,'Positive',(2,3,1))
        self.__probability_plot(self.negatives_list,fig,'Negative',(2,3,2))
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
        
        #if not set then generate single boxplot
        if not self.args.separate_plots:
            fileargs = ','.join(reffiles)
            reffiles=[]
            reffiles.append(fileargs)
        
        #ensure that the name of importance file and model matches.    
        model_impfile = self.args.ranked_features.strip().split('/')[-1]
        if not self.args.model_name in model_impfile:
            logging.error('The importance file doesn\'t match with model provided')
        
        #length of labels and boxplots should be same.
        if self.args.boxplot_labels != None:
            if len(reffiles) != len(self.args.boxplot_labels):
                logging.error('The length of reference files and labels do not match')
        else:
            self.args.boxplot_labels = [('input '+str(x) ) for x in range(len(reffiles))]
        
        #get names of features from features file
        self.tot_features_name = self.__get_feature_names()
        
        #generate feature_db 
        if self.args.feature_db == None:
            outputname = self.args.out+"feature_db_test.txt"
            self.__generate_feature_db(outputname,reffiles)
            self.args.feature_db = outputname
        
        #read feature_db into dict
        self.features_vals_dict = self.__get_features_dict()
        self.top_features_name, self.top_features_map = self.__get_top_features()
        
        #read ref files into a dict and generate list of features.
        self.input_ref_dict = self.__get_input_dict(reffiles)            
        inputs = []
        for input_args in reffiles:
            for files in input_args.split(','):
                inputs.append(files)
        self.features_list = [self.__get_features_from_dict(infile, self.features_vals_dict, self.input_ref_dict) for infile in self.input_ref_dict]
        
        #get all labels in ref data and generate list of all features
        self.labels = self.__get_label_names(reffiles)
        self.features_list_label= [self.__get_features_from_feature_db(self.features_vals_dict, label) for label in self.labels]
       
        #get T/P F/N and generate list 
        self.features_list_posneg,self.posneg_labels = self.__get_features_posneg()


    def __parse_manifest(self,manifest_file):
        manifest = defaultdict(list)
        if not self.args.deep:
            logging.error('manifest is only required in deep mode')
            raise Exception('Manifest file is only required in deep mode')
        
        if not manifest_file:
            return None
        
        man_stream = open(manifest_file)
        for line in man_stream:
            line = line.strip().split()
            if line[0]=='chrom':
                continue
            chrom = line[0]
            pos = line[1]
            ref = line[2]
            alt = line[3]
            start = line[4]
            end = line[5]
            manifest[(chrom,pos)].append((ref,alt,start,end))
        return manifest
        
    def __get_label_names(self, reference_files):
        labels = set()
        for files in reference_files:
            for filename in files.split(','):
                file_stream = open(filename)
                for line in file_stream:
                    if line[0] == '#':
                        continue
                    line = line.strip().split()
                    labels.add(line[2])
        return list(labels)
    
                  
    def __get_input_dict(self,reffiles):
        out_dict = defaultdict(list)
        tfile = None
        rfile = None
        nfile = None
        for files in reffiles:
            files = files.strip().split(',')
            for filename in files:
                tfile = None
                rfile = None
                nfile = None
                file_stream = open(filename)
                for line in file_stream:
                    if line[0] =='#':
                        l = line.strip().split()
                        if l[1]=='tumour':
                            tfile = l[2]
                        if l[1]=='normal':
                            nfile = l[2]
                        if l[1] == 'reference':
                            rfile = l[2]
                    else:
                        if not all( (tfile,nfile,rfile) ):
                            logging.error('could not read bam paths in ref file')
                        l = line.strip().split()
                        chromosome = l[0]
                        pos = l[1]
                        out_dict[str(files)].append( (tfile,nfile,rfile,chromosome,pos) )
                file_stream.close()
        return out_dict
    
    def __get_feature_names(self):
        if self.args.deep:
            feature_set = features_deep.Features()
        else:
            feature_set = features.Features()
        feature_names = feature_set.get_feature_names()
        return feature_names
                                             
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
                                
    def __generate_feature_db(self,outputname,reffiles):
        missing_positions = self.__get_missing_positions(reffiles)
        self.__write_feature_db(missing_positions,outputname)   
        
    def __get_missing_positions(self,reffiles):
        missing_positions = []
        manfile = None
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

                        if l[1] == 'manifest':
                            manfile = l[2]
                        continue

                    ## check if required bam files/reference are specified in the training data
                    if not all((tfile, nfile, rfile)):
                        logging.error("'%s' does not contain the required paths to bam/reference files" % infile)
                        continue
            
                    chromosome = l[0]
                    pos = l[1]
                    label = l[2]
                    missing_positions.append((chromosome,pos,tfile,nfile,rfile,manfile,label))
        return missing_positions
    
    def __write_feature_db(self,missing_positions,outputname):
        feature, keys = self.__extract_features(missing_positions)
        file_stream = open(outputname,'w')
        for i,_ in enumerate(keys):
                     
            k = ';'.join(keys[i])
            f = feature[i].tolist()[0]
            label = feature[i].tolist()[1]
            file_stream.write(k+'\t'+str(f)+'\t'+ str(label)+'\n')
        file_stream.close()
        
    def __get_features_dict(self):
        features_vals_dict = {}
        file_stream = open(self.args.feature_db,'r')
        for line in file_stream:
            l = line.strip().split('\t')
            key = l[0].split()[0]
            value = l[1]
            label = l[2]
            features_vals_dict[key] = (value,label)
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
    
    #new function that reads all stuff from feature_db
    def __get_features_from_feature_db(self,feature_dict,label):
        features = []
        for key,value in feature_dict.iteritems():
            key = key.strip().split(';')
            label_feature_dict = value[1]
            if not label == label_feature_dict:
                continue
            features.append(eval(value[0]))
        return features
    
    #retrieve and filter features from feature_db_dict        
    def __get_features_posneg(self):
        vcf_dict = self.__get_vcf_dict()
        true_positives = []
        true_negatives = []
        false_positives = []
        false_negatives = []
        for key,value in self.features_vals_dict.items():
            rfile,nfile,tfile,chromosome,pos = key.strip().split(';')
            label = value[1]
            features = value[0]
            pr = vcf_dict[(tfile,nfile,rfile,chromosome,pos)]
            
            if label in self.args.positive_labels:
                if pr[0]>=self.args.filter_threshold:
                    true_positives.append( eval(features) )
                else:
                    false_negatives.append( eval(features) )
            else:
                if pr[0]<self.args.filter_threshold:
                    true_negatives.append( eval(features) )
                else:
                    false_positives.append( eval(features) )
        
        return [true_positives,true_negatives,false_positives,false_negatives],['True Positives','True Negatives','False Positives','False Negatives']
    
        
    def __get_features_from_dict(self,label,f_dict,ref_dict):
        features = []
        for tfile,nfile,rfile,chromosome,pos in ref_dict.get(label):      
            key = ';'.join([rfile,nfile,tfile,chromosome,pos])
            try:
                features.append(eval(f_dict[key][0]))
            except KeyError:
                logging.error('error: cannot find key "%s"\n' % str(key))
        return features 

    def __get_flanking_regions(self,chromosome,position,manifest):
        if manifest != None:
            vals = manifest.get((chromosome,position))
            if vals != None:
                start = vals[2]
                end = vals[3]
                return start,end

        #If manifest file is not provided or region not available in manifest
        start = position - 25
        end = position + 26
        return start,end
       
    def __extract_features(self,missing_positions):
        data = defaultdict(list)
        contamination = (float(30), float(30), float(70), float(0))
        
        for chromosome,pos,tfile,nfile,rfile,manfile,label in missing_positions:
            data[(tfile,nfile,rfile,manfile)].append((chromosome,int(pos),label))
        
        features_buffer = []
        keys_buffer    = []
    
        for tfile, nfile, rfile,manfile in data.keys():
            logging.info("reading from tumour:"+ tfile)
            logging.info("reading from normal:"+ nfile)
            t_bam = pybamapi.Bam(bam=tfile, reference=rfile, coverage=1)
            n_bam = pybamapi.Bam(bam=nfile, reference=rfile, coverage=1)
        
            for chromosome, position, label in data[(tfile, nfile, rfile, manfile)]:
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
                
                if self.args.deep:
                    chromosome = t_bam.get_chromosome_name(chromosome_id)
                    self.manifest = self.__parse_manifest(manfile)
                    start,end = self.__get_flanking_regions(chromosome, position, self.manifest)
                    bam  = pybamapi.PairedBam(tumour=tfile,
                                              normal=nfile, 
                                              reference=rfile,
                                              coverage=self.args.coverage, rmdups=False,
                                              mapq_threshold=self.args.mapq_threshold,
                                              baseq_threshold = self.args.baseq_threshold)
            
                    tt_bg = bam.t_bam.get_tuples([[chromosome,start,end]])
                    tt_bg = [tval for tval in tt_bg if tval[0]!=position]
                
                    nt_bg = bam.n_bam.get_tuples([[chromosome,start,end]])
                    nt_bg = [nval for nval in nt_bg if nval[0]!=position]
                    
                    rt_bg = []
                    for posval in range(start,end+1):
                        if posval == position:
                            continue
                        rt_bg.append([posval,bam.t_bam.get_reference_tuple(chromosome_id,posval)])
                    
                    if self.manifest:
                        indexes = self.manifest.get((chromosome,position))
                    else:
                        indexes = None

                ## calculate features
                if self.args.deep:
                    feature_set = features_deep.Features(tt,nt,rt,tt_bg,nt_bg,rt_bg,indexes)
                else:
                    feature_set = features.Features(tt,nt,rt)
                temp_features = feature_set.get_features()
            
                features_buffer.append((temp_features,label))
                keys_buffer.append((rfile, nfile, tfile, chromosome, position))
        
        features_buffer = numpy.array(features_buffer)
        keys_buffer = numpy.array(keys_buffer)
        return features_buffer, keys_buffer

    def __get_flier_counts(self, bplot_obj,fval):
        #get the outliers
            ax_fliers = []
            for i in xrange(len(bplot_obj['boxes'])):
                fliers_above = len(bplot_obj['fliers'][i*2]._y)
                fliers_below = len(bplot_obj['fliers'][i*2+1]._y)
                ax_fliers.append( str(fliers_above)+','+str(fliers_below) )
            
            #if some value in vector is [],then no boxplot.so we need to adjust the value to appropriate plot    
            if not len(ax_fliers) == len(fval):
                indices = [i for i, j in enumerate(fval) if j == [] ]
                for val in indices:
                    ax_fliers.insert(val, [])
                                    
            #get the upper most and lower most flier positions and caps for boxplot         
            ylim_upper = None
            upper_cap =None
            lower_cap = None
            ylim_lower = None   
            for i in xrange(len(bplot_obj['boxes'])):
                try:
                    uppercap = bplot_obj['caps'][i*2]._y[0]
                    highestflier = max(bplot_obj['fliers'][i*2]._y)
                    if highestflier > uppercap*100:
                        if not uppercap == 0:
                            if uppercap>upper_cap:
                                upper_cap = uppercap
                            if highestflier>ylim_upper:
                                ylim_upper = highestflier
                except:
                    #if unable to set the axis, continue without changing them
                    pass
             
                try:
                    lowercap = bplot_obj['caps'][i*2+1]._y[0]
                    lowestflier = min(bplot_obj['fliers'][i*2+1]._y)
                    if lowestflier > lowercap/100:
                        if not lowercap == 0:
                            if lowercap < lower_cap:
                                lower_cap = lowercap
                            if lowestflier < ylim_lower:
                                ylim_lower = lowestflier  
                except:
                    #if unable to set the axis, continue without changing them
                    pass
            
            if ylim_upper and upper_cap is not None:
                if ylim_upper > upper_cap*100:
                    ylim_upper = int(upper_cap*100)
            
            if ylim_lower and lower_cap is not None:
                if ylim_lower < lower_cap/100:
                    ylim_upper = int(lower_cap/100)
            return ylim_lower,ylim_upper,ax_fliers
        
    def __plot_fvals(self,fval,axes,labels):
        fval_count = [len(x) for x in fval]
        bplot = axes.boxplot(fval)
        if axes.axis()[3] >40000:
            axes.set_yscale('symlog')
        ylim_lower,ylim_upper,ax_fliers = self.__get_flier_counts(bplot, fval)
        normal_xlabel_names = ['%s(%s)\n(%s)' % (labels[i], y,ax_fliers[i]) for i,y in enumerate(fval_count)]
         
        axes.tick_params(axis='x', labelsize=8)
        axes.set_ylabel('Distribution',fontsize = 8)
        axes.set_xticklabels( normal_xlabel_names, rotation=45, fontsize=8)
        axes.set_title('Feature Distributions for whole data',fontsize = 8)
        axes.yaxis.set_tick_params(labelsize=8)
        return ylim_lower, ylim_upper
   
    def __rescale_plots(self,normal_ylim_lower,normal_ylim_upper,label_ylim_lower,label_ylim_upper,posneg_ylim_lower,posneg_ylim_upper,f):
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
            
        #just to show boundaries(if overlaps with axis)
        ylims = pyplot.ylim()
        if ylims[0] == 0:
            pyplot.ylim((ylims[0]-0.05*ylims[1]), (1.10*ylims[1]))
        else:
            pyplot.ylim( (0.90*ylims[0]), (1.10*ylims[1]) )  
     
    def boxplot_plot(self):
        plots = []
        for i,v in enumerate(self.top_features_name):
            index = self.top_features_map[v]
            if not self.tot_features_name[index] == v:
                logging.error('the feature name and feature values don\'t match')
            
            f = pyplot.figure(figsize=(15, 15))
            f.set_dpi(150)
            gs1 = gridspec.GridSpec(1, 3)           
   
            #normal
            fnvalue = []
            for feature in self.features_list:
                fnv = []
                for p in feature:
                    fnv.append(p[index])
                fnvalue.append(fnv)
            ax1 = f.add_subplot(gs1[0])
            normal_ylim_lower,normal_ylim_upper = self.__plot_fvals(fnvalue,ax1,self.args.boxplot_labels)
                                         
            #labels
            flvalue = []
            for feature in self.features_list_label:
                flv = []
                for p in feature:
                    flv.append(p[index])
                flvalue.append(flv)
            ax2 = f.add_subplot(gs1[1],sharey=ax1)
            label_ylim_lower, label_ylim_upper = self.__plot_fvals(flvalue, ax2,self.labels)

            posneg_ylim_upper = None
            posneg_ylim_lower = None
            #posneg
            if self.args.input_files:
                fpnvalue = []
                for feature in self.features_list_posneg:
                    fpnv = []
                    for p in feature:
                        fpnv.append(p[index])
                    fpnvalue.append(fpnv)
                ax3 = f.add_subplot(gs1[2],sharey=ax1)
                posneg_ylim_lower,posneg_ylim_upper = self.__plot_fvals(fpnvalue, ax3,self.posneg_labels)

            if self.args.rescale:
                self.__rescale_plots(normal_ylim_lower, normal_ylim_upper, label_ylim_lower, label_ylim_upper, posneg_ylim_lower, posneg_ylim_upper, f)
            
            xlabel_description = 'Features (Count of positions) (Outliers above the plot, Outliers below the plot)'
            f.text(0.90, 0.98, 'importance:'+str(i+1), rotation='horizontal',horizontalalignment='center', verticalalignment='bottom',fontsize = 8)
            f.text(0.90, 0.97,'Threshold ='+str(self.args.filter_threshold), horizontalalignment='center',verticalalignment='top',fontsize=8)
            f.text(0.50, 0.02, xlabel_description, rotation='horizontal',horizontalalignment='center', verticalalignment='top',fontsize = 10)
            f.text(0.50, 0.97, v+' ('+self.args.model_name+')', rotation='horizontal',horizontalalignment='center', verticalalignment='bottom',fontsize = 10)
            
            gs1.tight_layout(f, rect=[0, 0, 1, .96], h_pad=0.5)
            plots.append(f)
        return plots    
