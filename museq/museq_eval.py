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


mutationSeq_version="4.3.2"
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
    output_folder = arguments.out
        
    for i in xrange( len(reffiles) ):   
        reference_file = reffiles[i]
        #parse the pos file
        file_stream = open(reference_file,'r')
        tfile =None
        nfile = None
        rfile = None
        manfile = None
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
                if l[1] == 'manifest':
                    manfile = l[2]
            else:
                output.append(l[0]+':'+l[1]+'\n')
        file_stream.close()
        #update arguments 
        if not all((tfile,nfile,rfile)):
            logging.error('Invalid input (one of paths is missing)')    
            
        arguments.out = output_folder+reffiles[i].strip().split('/')[-1]+'.vcf'
        #create a positions file for classifier
        file_stream_w = open(arguments.out + '.tmp','w')
        for line in output:
            file_stream_w.write(line)
        file_stream_w.close()
        
        output_vcf.append(arguments.out )
        
        arguments.interval = None
        arguments.positions_file = arguments.out + '.tmp'
        arguments.samples = ['tumour:'+tfile, 'normal:'+nfile, 'reference:'+rfile, 'model:'+arguments.model]
        arguments.manifest = manfile
        
        logging.info("initializing a Classifier")
        classifier = bamutils.Classifier(arguments)

        logging.info("getting positions")
        target_positions = classifier.get_positions()

        logging.info("generating tuple iterator")
        tuples = classifier.bam.get_tuples(target_positions)

        logging.info("generating features iterator")
        features = classifier.get_features(tuples)

        if arguments.export_features is not None:
            logging.info("exporting features")
            classifier.export_features(features)

        probabilities = classifier.predict(features)
        classifier.print_results(probabilities)

        logging.warning("successfully completed.\n")
        
        #remove the positions file
        os.remove(arguments.out + '.tmp')
        
    return output_vcf
    
def run_museqeval(arguments,features_only):
    if arguments.out:
        pdfout = PdfPages(arguments.out+''+arguments.model_name+'_test.pdf')
    else:
        pdfout = PdfPages('boxplots_test.pdf')
        
    if not features_only:
        plots = museq_eval_utils.museq_plots(arguments)
        fig = plots.generate_plots()
        pdfout.savefig(fig)

    boxplot = museq_eval_utils.box_plots(arguments,arguments.reference_files)
    plots = boxplot.boxplot_plot()
    for plot in plots:
        pdfout.savefig(plot)
    pdfout.close()
    
    
#============================================
#Run code according to the arguments provided
#============================================
    
if args.input_files == None and args.plot_features_only == False:
    reffiles = []
    #create a list of all ref files, irrespective of comma
    for refs in args.reference_files:
        for ref in refs.strip().split(','):
            reffiles.append(ref)          
    #check for trailing /, o.w. files are saved in parent dir
    if not args.out[-1] == '/':
        args.out = args.out + '/'
        
    arguments = copy.deepcopy(args)
    vcf_outputs = run_classifier(arguments,reffiles)
    
    args.input_files = vcf_outputs
    args.model_name = args.model.strip().split('/')[-1]
    
    run_museqeval(args,args.plot_features_only)
else:
            
    args.model_name = args.model.strip().split('/')[-1]
    run_museqeval(args,args.plot_features_only)
    
    
    
