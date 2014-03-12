'''
Created on Feb 21, 2014

@author: dgrewal
'''
import museq_eval_utils
import museq_eval_ui
import logging
import copy
import os
from matplotlib.backends.backend_pdf import PdfPages


mutationSeq_version="4.1.0"
args = museq_eval_ui.args 

if args.verbose:
    level = logging.DEBUG    

else:
    level = logging.WARNING
    
logging.basicConfig(filename = args.log_file, 
                    format   = '%(asctime)s %(message)s', 
                    #datefmt = '%m/%d/%Y %I:%M:%S %p', 
                    level = level)

logging.warning("<<< mutationSeq_" + mutationSeq_version + " Classifier started >>>")
logging.info("importing required modules")
import bamutils

logging.info(args)

def run_classifier(arguments,reffiles):
    output_vcf = []
    if not args.out==None and len(args.out) == len(reffiles):
        logging.error('There should be an outputfile for each positions file specified')
    arguments.interval = None
    
    for i in xrange( len(reffiles) ):   
        reference_file = reffiles[i]
        
        if not args.out:
            arguments.out = os.getcwd()+'/'+reference_file.strip().split('/')[-1]+'.vcf'
        else:
            arguments.out = args.out[i]
        output_vcf.append(arguments.out)
        #parse the pos file
        file_stream = open(reference_file,'r')
        tfile =None
        nfile = None
        rfile = None
        output = []
        for line in file_stream:
            l=line.strip().split()
            if line[0] == '#':
                if l[1] == 'tumour':
                    tfile = l[2]
                elif l[1] == 'normal':
                    nfile = l[2]
                if l[1] == 'reference':
                    rfile = l[2]
            else:
                output.append(l[0]+':'+l[1]+'\n')
        file_stream.close()
        
        #create a positions file for classifier
        file_stream_w = open(arguments.out+'.tmp','w')
        for line in output:
            file_stream_w.write(line)
        file_stream_w.close()
        
        #update arguments 
        if not all((tfile,nfile,rfile)):
            logging.error('Invalid input (one of paths is missing)')
        
        arguments.positions_file = arguments.out +'.tmp'
        arguments.samples = ['tumour:'+tfile, 'normal:'+nfile, 'reference:'+rfile, 'model:'+args.model]
        
        logging.info("initializing a Classifier")
        classifier = bamutils.Classifier(arguments)

        logging.info("getting positions")
        target_positions = classifier.get_positions()

        logging.info("generating tuple iterator")
        tuples = classifier.bam.get_tuples(target_positions)

        logging.info("generating features iterator")
        features = classifier.get_features(tuples)

        if args.export_features is not None:
            logging.info("exporting features")
            classifier.export_features(features)

        probabilities = classifier.predict(features)
        classifier.print_results(probabilities)

        logging.warning("successfully completed.\n")
        os.remove(arguments.out+'.tmp')
        
    return output_vcf
    
def run_museqeval(arguments,features_only):
    if args.out:
        pdfout = PdfPages(arguments.out+''+arguments.model_name+'_test.pdf')
    else:
        pdfout = PdfPages('boxplots_test.pdf')
        
    if not features_only:
        plots = museq_eval_utils.museq_plots(args)
        fig = plots.generate_plots()
        pdfout.savefig(fig)

    boxplot = museq_eval_utils.box_plots(args,arguments.reference_files)
    plots = boxplot.boxplot_plot()
    for plot in plots:
        pdfout.savefig(plot)
    pdfout.close()
    
    boxplot.remove_temp_files()
    
#============================================
#Run code according to the arguments provided
#============================================
    
if args.input_files == None and args.plot_features_only == False:
    reffiles = []
    for refs in args.reference_files:
        for ref in refs.strip().split(','):
            reffiles.append(ref)
    arguments = copy.deepcopy(args)
    
    vcf_outputs = run_classifier(arguments,reffiles)
    
    args.input_files = vcf_outputs
    if args.out:
        #get path from out
        args.out = args.out[0].strip().split('/')
        if len(args.out) == 1:
            args.out = os.getcwd()+'/'
        else:
            args.out = '/'.join(args.out[:-1])+'/'
    else:
        args.out= os.getcwd()  +'/'  
    args.model_name = args.model.strip().split('/')[-1]
    
    run_museqeval(args,args.plot_features_only)
else:
    if args.out:
        args.out = args.out[0].strip().split('/')
        if len(args.out) == 1:
            args.out = os.getcwd()+'/'
        else:
            args.out = '/'.join(args.out[:-1])+'/'
    else:
        args.out = os.getcwd()+'/'
            
    args.model_name = args.model.strip().split('/')[-1]
    run_museqeval(args,args.plot_features_only)
    
    
    
