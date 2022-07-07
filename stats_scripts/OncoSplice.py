#!/usr/bin/env python

#Copyright 2017 Cincinnati Children's Hospital Medical Center, Research Foundation
#Author Meenakshi Venkatasubramanian - altanalyze@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

""" OncoSplice Module
https://github.com/venkatmi/oncosplice

This workflow identifies alternative splicing events (MultiPath-PSI - DEFAULT)
that can be clustered into coherent subtypes in a disease context. In cancer,
a patient may possess multiple disease mutations/rearangments, hence the software
considers overlapping subtypes across patients and identifies/excludes potential
confounding broad splicing signatures that represent important cellular impacts
or batch effects (signature depletion).

Steps applied in this workflow:
1 - Run splice-ICGS (Feature Selection)
2 - Block identification (k analysis)
3 - NMF Analysis (Initial subtype identification)
4 - Filter Event Annotation
5 - Meta data analysis (differential expression)
6 - Expand clusters (SVM sample classification)
7 - Mutation enrichment (MAF or VCF - optional)
8 - Correlation depletion (excluded biological confounding signatures)
"""

test_mode = False

import sys, string, os, shutil
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import traceback
import sys, string, os
import RNASeq
import numpy as np
from stats_scripts import RNASeq_blockIdentification
from stats_scripts import NMF_Analysis; reload(NMF_Analysis)
from stats_scripts import filterEventAnnotation
from stats_scripts import metaDataAnalysis
from stats_scripts import ExpandSampleClusters; reload(ExpandSampleClusters)
from import_scripts import sampleIndexSelection
from stats_scripts import CorrelationDepletion
import UI
import multiprocessing as mlp
import export
upd_guides=[]
import operator
from collections import OrderedDict
from collections import defaultdict
from stats_scripts import Kmeans
from stats_scripts import MutationEnrichment_adj as ME
from visualization_scripts import Orderedheatmap
import warnings
warnings.filterwarnings('ignore')

def filterPSIValues(filename,percentCutoff=0.75,filterStatus=True):
    """ Filter the PSI file to only include events in which >75% of the samples have PSI values """
    
    firstRow=True
    header = True
    rows=0
    filtered=0
    ### Filtered Export file
    new_file = filename[:-4]+'-'+str(int(100*percentCutoff))+'p.txt'
   
    if filterStatus:
        ea = export.ExportFile(new_file)

    for line in open(filename,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            header = False
            eventindex=t.index('EventAnnotation')  ### This field is the last annotation column before sample PSI values 
            t = [t[1]]+t[eventindex+1:]
            header_length = len(t)-1
            minimum_values_present = int(float((header_length)-1.0)*percentCutoff)
            not_detected = header_length-minimum_values_present
            new_line = line
            if filterStatus:
                ea.write(new_line)
        else:
            if filterStatus: ### If sufficient PSI detected samples, write to new file
                t = [t[1]]+t[eventindex+1:]
                missing_values_at_the_end = (header_length+1)-len(t)
                missing = missing_values_at_the_end+t.count('')
                if missing<not_detected:
                    new_line = line
                    ea.write(new_line)
                    filtered+=1
        rows+=1

    if filterStatus:
        ea.close()
        return new_file,header_length
    else:
        ### Just return the number of samples in the file
        return header_length
    
def header_list(EventAnnot):
    head=0
    header=[]
    with open(EventAnnot, 'rU') as fin:
        for line in fin:
            if head==0:
                line = line.rstrip(os.linesep)
                line=string.split(line,'\t')
                startpos=line.index('EventAnnotation')
                #header.append('UID')
                for i in range(startpos+1,len(line)):          
                        header.append(line[i])
                head=1
            else:break
    return header

def FindTopUniqueEvents(dPSI_results_fn,psi,dPSI_results_dir):
    """ This function reports signatures with at least 100 splicing events.
    If none are present the function indicates that these results should be ommitted """
    
    head=0
    guidekeys=[]
    exportnam=os.path.join(dPSI_results_dir,"SplicingeventCount1.txt")
    export_class=open(exportnam,"a")

    tempkeys={}
    global upd_guides
    global train
    omitcluster=0
    
    unique_clusters={}

    for line in open(dPSI_results_fn,'rU').xreadlines():
        if head==0:
            line1=line.rstrip('\r\n')
            q= string.split(line1,'\t')
            head=1
            try:
                uid=q.index('UID')
                adjp=q.index('rawp')
                dpsi=q.index('dPSI')
                Clusterid=q.index('UpdatedClusterID')
                cutoff=0.1
                continue
            except Exception:
                uid=q.index('GeneID')
                adjp=q.index('rawp')
                dpsi=q.index('LogFold')
                Clusterid=q.index('GeneID')
                cutoff=0.58 ### 1.5 fold
        else:
            line1=line.rstrip('\r\n')
            q= string.split(line1,'\t')
            if abs(float(q[dpsi]))>cutoff and float(q[adjp])<0.05: ### Adjusted p-value cutoff of 0.01
                try:
                    tempkeys[q[Clusterid]].append([q[uid],float(q[adjp]),q[adjp+1]])
                except KeyError:
                    tempkeys[q[Clusterid]]=[[q[uid],float(q[adjp]),q[adjp+1]],]
                    
    for i in tempkeys:
        if len(tempkeys[i])>1:
            tempkeys[i].sort(key=operator.itemgetter(1),reverse=False)
            try:
                unique_clusters[0].append(tempkeys[i][0])
            except KeyError:
                unique_clusters[0]=[tempkeys[i][0],]
          
        else:
            try:
                unique_clusters[0].append(tempkeys[i][0])
            except KeyError:
                unique_clusters[0]=[tempkeys[i][0],]
    
    try:
        if len(unique_clusters[0])>1:
            unique_clusters[0].sort(key=operator.itemgetter(1))
            if len(unique_clusters[0])>100: ### ORIGINALY SET TO 120
                guidekeys=unique_clusters[0][0:150]
                for i in range(0,len(guidekeys)):
                    upd_guides.append(guidekeys[i][0])
            else:
                omitcluster=1
        else:
            omitcluster=1
        export_class.write(psi+"\t"+str(len(unique_clusters[0]))+"\n")
    except Exception:
        omitcluster=1
    
    return omitcluster

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def MergeResults(root_dir):
    files_to_merge=[]
    
    for filename in os.listdir(root_dir): 
        if ("_Results" in filename or "Kmeans" in filename)  and "._" not in filename and "ordered" not in filename:
            files_to_merge.append(root_dir+'/'+filename)
            
    import mergeFiles
    exportnam=os.path.join(root_dir,"MergedResult.txt")
    join_option='Intersection'
    uniqueOnly=False
    outputfile = mergeFiles.remoteCombineAllLists(files_to_merge,'',includeDatasetName=False,ExportPath=exportnam)
    return exportnam

def CompleteWorkflow(full_PSI_InputFile,EventAnnot,rho_cutoff,strategy,seq,gsp,forceBroadClusters,AnalysisRound,platform='PSI'):
    """ This function is used perform a single-iteration of the OncoSplice workflow (called from main),
    including the unsupervised splicing analysis (splice-ICGS) and signature depletion """
    
    ### Filter the EventAnnotation PSI file with non-depleted events from the prior round
    filtered_EventAnnot_dir=filterEventAnnotation.FilterFile(full_PSI_InputFile,EventAnnot,AnalysisRound)
    
    try:
        print "Running splice-ICGS for feature selection - Round"+str(AnalysisRound)

        ### Reset the below variables which can be altered in prior rounds
        gsp.setGeneSelection('')
        gsp.setGeneSet('None Selected')
        gsp.setPathwaySelect([])
        species = gsp.Species()
        if forceBroadClusters == True:
            ### Find Broad clusters with at least 25% of all samples
            originalSamplesDiffering = gsp.SamplesDiffering()
            gsp.setSamplesDiffering(int(SampleNumber*0.25))
            
        print 'Number varying samples to identify:',gsp.SamplesDiffering()
        
        #if AnalysisRound > 1:
        graphic_links3,events = RNASeq.singleCellRNASeqWorkflow(species, 'exons', full_PSI_InputFile,mlp,exp_threshold=0, rpkm_threshold=0, parameters=gsp)
        if forceBroadClusters == True:
            gsp.setSamplesDiffering(originalSamplesDiffering)

        #if AnalysisRound > 1:
        dPSI_results_fn=graphic_links3[-1][-1]
        dPSI_results_fn=dPSI_results_fn[:-4]+'.txt'
        print 'using:',[dPSI_results_fn]
        #dPSI_results_fn = '/Users/saljh8/Documents/GitHub/oncosplice/demo/test1/ExpressionInput/amplify//DataPlots/Clustering-exp.input-Guide3-hierarchical_euclidean_correlation.txt'
       
        print "Running block identification for k analyses - Round"+str(AnalysisRound)
        ### Parameters are fixed as they are distinct 
        RNASeq_blockIdentification.correlateClusteredGenesParameters(dPSI_results_fn,rho_cutoff=0.4,hits_cutoff=4,hits_to_report=50,ReDefinedClusterBlocks=True,filter=True) 
        dPSI_results_fn_block=dPSI_results_fn[:-4]+'-BlockIDs.txt'
        NMFinput, k = NMF_Analysis.FilterGuideGeneFile(dPSI_results_fn,dPSI_results_fn_block,full_PSI_InputFile,AnalysisRound,platform)

    except ZeroDivisionError:
        print 'UNKNOWN ERROR!!!!! Setting k=0' 
        print traceback.format_exc()
        k=0
    
    print "Round =", AnalysisRound,'and k =', k
    if AnalysisRound == 1:
        if force_broad_round1:
            k = 2
        else:
            NMFinput,k = NMF_Analysis.FilterGuideGeneFile(dPSI_results_fn,dPSI_results_fn,full_PSI_InputFile,AnalysisRound,platform) ### Just use the Guide 3 file alone
    if k < 2:
        NMFinput,k = NMF_Analysis.FilterGuideGeneFile(dPSI_results_fn,dPSI_results_fn,full_PSI_InputFile,AnalysisRound,platform) ### Just use the Guide 3 file alone
        #k = 2
        NMFinput,Rank=NMF_Analysis.FilterGuideGeneFile(GuideFile,GuideFile_block,processedInputExpFile,iteration,platform,uniqueIDs,symbolIDs)
        
    print "Round =", AnalysisRound,'and k =', k
    if k>1:
        ### ADJUST THE k - MUST UPDATE!!!!
        if AnalysisRound == 1:
            if k < 2:
                k = 30
        else:
            if k > 2:
                k = 30
        print "Round =", AnalysisRound,'and k =', k
        
        try:
            flag,full_PSI_InputFile = performNMF(species, NMFinput, full_PSI_InputFile, filtered_EventAnnot_dir, k,AnalysisRound, strategy,expFile=full_PSI_InputFile)
        except:
            print traceback.format_exc()
            k+=1
            print 'Adjusted k =',k
            try:
                flag,full_PSI_InputFile = performNMF(species, NMFinput, full_PSI_InputFile, filtered_EventAnnot_dir, k, AnalysisRound, strategy,expFile=full_PSI_InputFile)
                print traceback.format_exc()
            except:
                k = 30
                print 'Adjusted k = 30'
                try:
                    flag,full_PSI_InputFile = performNMF(species, NMFinput, full_PSI_InputFile, filtered_EventAnnot_dir, k,AnalysisRound, strategy,expFile=full_PSI_InputFile)
                    print traceback.format_exc()
                except:
                    flag = True
                    pass ### will force k-means below
    
    if k<2:
        if k==1:
            try:
                print "Running K-means analyses instead of NMF - Round"+str(AnalysisRound)
                header=[]
                header=Kmeans.header_file(dPSI_results_fn_block)
                Kmeans.KmeansAnalysis(dPSI_results_fn_block,header,full_PSI_InputFile,AnalysisRound)
                if AnalysisRound == 1:
                    flag=True
                else:
                    flag=False
            except Exception:
                print 'WARNING!!!!! DID NOT RUN K-MEANS!!!!!'
                print traceback.format_exc()
                AnalysisRound = True
        else:
            flag=False
     
    return flag,full_PSI_InputFile,filtered_EventAnnot_dir

def performNMF(species, NMFinput, full_PSI_InputFile, filtered_EventAnnot_dir, k, AnalysisRound, strategy,expFile=None):
    """ Run NMF and determine the number of valid clusters based on the magnitude of detected differential splicing """
    
    use_adjusted_p=True
    platform = 'PSI'
           
    print "Running NMF analyses for dimension reduction using "+str(k)+" k - Round"+str(AnalysisRound)
    NMFResult,BinarizedOutput,metaData,Annotation=NMF_Analysis.NMFAnalysis(expFile,NMFinput,k,AnalysisRound,strategy) ### This is where we get the correct version
    print "Running metaData Analyses for finding differential splicing events"
    print [Annotation]
    rootdir,CovariateQuery=metaDataAnalysis.remoteAnalysis(species,filtered_EventAnnot_dir,metaData,'PSI',0.1,use_adjusted_p,0.05,meta_description_file=Annotation)
    counter=1
    dPSI_results_dir=rootdir+CovariateQuery
    global upd_guides
    upd_guides=[]
    name=[]
    group=[]
    grplst=[]
    for filename in os.listdir(dPSI_results_dir):
        if filename.startswith("PSI."):
            dPSI_results_fn=os.path.join(dPSI_results_dir, filename)
            dPSI_comparison_alt_name=string.replace(filename,"PSI.","")
            omitcluster=FindTopUniqueEvents(dPSI_results_fn,dPSI_comparison_alt_name,dPSI_results_dir)
            if omitcluster==0: ### Hence, clustering succeeded and did not fail in this dPSI comparison
                group.append(counter)
                name.append(string.replace(filename,"PSI.",""))
                counter+=1
                
    print counter, 'robust splicing subtypes identified in round',AnalysisRound
    if counter>0: #counter>2 --- changed to 0 to force NMF
        root_dir = export.findParentDir(full_PSI_InputFile)
        output_dir = root_dir+'OncoInputs'
        if os.path.exists(output_dir)==False:
            export.createExportFolder(output_dir)

        output_file = output_dir+'/SVMInput-Round'+str(AnalysisRound)+'.txt'
        ExpandSampleClusters.filterRows(full_PSI_InputFile,output_file,filterDB=upd_guides,logData=False)
        header=ExpandSampleClusters.header_file(output_file)
        print "Running SVM prediction for improved subtypes - Round"+str(AnalysisRound)

        train,centroid_heatmap_input=ExpandSampleClusters.TrainDataGeneration(output_file,BinarizedOutput,name)
        grplst.append(group)
        ExpandSampleClusters.Classify(header,train,output_file,grplst,name,AnalysisRound,platform,output_dir,root_dir)
        header=CorrelationDepletion.header_file(NMFResult)
        
        output_file=output_dir+'/DepletionInput-Round'+str(AnalysisRound)+".txt"
        sampleIndexSelection.filterFile(full_PSI_InputFile,output_file,header)
        print "Running Correlation Depletion - Round"+str(AnalysisRound)
        commonkeys,count=CorrelationDepletion.FindCorrelations(NMFResult,output_file,name)
        Depleted=CorrelationDepletion.DepleteSplicingevents(commonkeys,output_file,count,full_PSI_InputFile)
        full_PSI_InputFile=Depleted
    
        flag=True ### Indicates that K-means was not run - hence, another round of splice-ICGS should be performed

    return flag,full_PSI_InputFile

def formatMetaData(filename):
    """ Export metadata annotations from a matrix that consist of mutations, numerical values
    categorical values and true/false values to a simple two column format for OncoSplice"""
    
    export_path = filename[:-4]+'_reorganized'+'.txt'
    eo = export.ExportFile(export_path)
    metadata=[]
    quantitative_values={}
    quantitative_sample_values={}
    row_count=1
    for line in open(filename, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if row_count==1:
            headers = t[1:]
        elif row_count==2:
            if t[0] == 'FORMAT':
                dataFormat = t[1:]
        else:
            sampleID = t[0]
            index=0
            for value in t[1:]:
                if dataFormat[index] == 'TRUE-FALSE':
                    if value == 'TRUE' or value == 'y' or value == 'yes' or value == 'true' or value == 'YES' or value == 'True':
                        ### The sample is labeled as the column header (e.g., Alive)
                        metadata.append([sampleID,headers[index]])
                if dataFormat[index] == 'QUANTITATIVE':
                    if len(headers[index])>0 and value != '':
                        ### Prior to export, define low, median and high clusters
                        try:
                            value = float(value)
                            try:
                                quantitative_values[headers[index]].append(value)
                            except:
                                quantitative_values[headers[index]] = [value]
                            try:
                                quantitative_sample_values[headers[index]].append([sampleID,value])
                            except:
                                quantitative_sample_values[headers[index]] = [[sampleID,value]]
                        except:
                            pass ### Invalid non-numeric
                if dataFormat[index] == 'VERBOSE':
                    if len(headers[index])>0 and value != '':
                        metadata.append([sampleID,headers[index]+'('+value+')'])
                if dataFormat[index] == 'MUTATION':
                    if len(headers[index])>0 and value != '':
                        metadata.append([sampleID,headers[index]])
                        if 'p.' in value:
                            value = string.replace(value,'p.','')
                        if '; ' in value:
                            value = string.split(value,';')[0]
                        if len(value)>1:
                            metadata.append([sampleID,headers[index]+'-'+value[:-1]])
                index+=1
        row_count+=1
    
    for annotation in quantitative_values:
        values = quantitative_values[annotation]
        one_third = len(values)/3
        bottom = values[:one_third]
        middle = values[one_third:-1*one_third]
        top = values[-1*one_third:]
        for (sampleID, value) in quantitative_sample_values[annotation]:
            if value in bottom:
                metadata.append([sampleID,annotation+'-low'])
            elif value in middle:
                metadata.append([sampleID,annotation+'-mid'])
            elif value in top:
                metadata.append([sampleID,annotation+'-high'])
            else:
                print value,'value is out-of-range!!!'; sys.exit()
                
    ### Write these metadata annotations out to a two column file
    for (sampleID,annotation) in metadata:
        if len(sampleID)>1:
            eo.write(sampleID+'\t'+annotation+'\n') 
    eo.close()
    return export_path
    
def checkmetaDataFormat(filename):
    ### Allow more complex metadata file formats with a matrix of diverse annotations
    for line in open(filename, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
    if len(t)>3:
        try:
            filename = formatMetaData(filename)
        except:
            pass
    return filename  

def getPSIHeader(full_PSI_InputFile):
    firstRow=True
    PSI_sample_headers=[]
    for line in open(full_PSI_InputFile,'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstRow:
            PSI_sample_headers = values
            firstRow=False
        else:
            break
    return PSI_sample_headers
    
if __name__ == '__main__':
    """ OncoSplice user-input and default parameters for iterative analysis of PSI splicing files """
    
    # python stats_scripts/OncoSplice.py --i /Users/saljh8/Desktop/dataAnalysis/SalomonisLab/TCGA/COAD/Non-Broad/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt --EnrichmentOnly yes --metadata /Users/saljh8/Desktop/dataAnalysis/SalomonisLab/TCGA/COAD/final-metadata.txt

    import getopt
    """ Below are the default variables for OncoSplice """
    seq = "bulk"
    rho_cutoff = 0.3
    strategy = "conservative"
    filters = True
    mode = "iterative"
    Mutationref = ""
    percentCutoff = 0.75
    
    ### Splice-ICGS Options
    species = "Hs"
    row_method = 'hopach'
    column_method = 'hopach'
    row_metric = 'correlation'
    column_metric = 'euclidean'
    platform = "PSI"
    excludeCellCycle = False
    restrictBy = 'protein_coding'
    restrictBy = None
    featurestoEvaluate = 'Genes'
    ExpressionCutoff = 0
    CountsCutoff = 0
    FoldDiff = 1.2
    SamplesDiffering = 4
    removeOutliers = False
    forceBroadClusters = False
    EnrichmentOnly = False
    metaDataMatrixFormat = False
    force_broad_round1 = False
    
    """ Below are the user-defined variables for OncoSplice """
    flag=True
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['species=','platform=','EventAnnotation=','rho=',
                                            'strategy=','filter=','mode=','Mutationref=','Assoc=','row_method',
                                            'column_method=','ExpressionCutoff=','normalization=','CountsCutoff=',
                                            'FoldDiff=','SamplesDiffering=','removeOutliers=','percentCutoff=',
                                            'forceBroadClusters=','i=','metadata=','EnrichmentOnly=', 'm=',
                                            'subtypeDiscovery=','row_metric=','column_metric=','SamplesDiffering-',
                                            'removeOutliers=','FoldDiff=','ExpressionCutoff=','normalization=',
                                            'metaDataMatrixFormat=','Expand=','force_broad_round1='])
        print options
        for opt, arg in options:
            if opt == '--EventAnnotation' or opt == '--i':EventAnnot=arg
            elif opt == '--strategy' or opt == '--subtypeDiscovery':
                ### How stringent should OncoSplice be to EXCLUDE clusters with evidence of sample overlap
                ### Options: stringent (DEFAULT), conservative
                strategy=arg 
            elif opt == '--filter':filters=arg
            elif opt == '--mode':mode=arg
            elif opt == '--metaDataMatrixFormat' or opt == '--Expand':
                if string.lower(arg) == 'true' or string.lower(arg) == 'yes':
                    metaDataMatrixFormat = True
                else:
                    metaDataMatrixFormat = False
            elif opt == '--force_broad_round1':
                if string.lower(arg) == 'yes' or string.lower(arg) == 'yes':
                    ### Rather than forcing 30 or 2 clusters in Round1 - uses the blockID based k
                    force_broad_round1 = True
            elif opt == '--Mutationref' or opt == '--metadata' or opt == '--m':
                Mutationref=arg
            elif opt == '--Assoc':
                arg = string.lower(arg)
                if 'yes' in arg or 'true' in arg:
                    EnrichmentOnly=False
                else:
                    EnrichmentOnly = True
            elif  opt == '--EnrichmentOnly':
                ### Skip the full OncoSplice workflow and just perform mutation/metadata enrichment
                ### True, False (Default), yes, no
                arg = string.lower(arg)
                if 'yes' in arg or 'true' in arg:
                    EnrichmentOnly=True
                else:
                    EnrichmentOnly = False                
            elif opt == '--percentCutoff':
                ### % cutoff for # of samples with detected PSI values for inclusion
                ### Integer: 0-100
                percentCutoff=float(percentCutoff)/100.00
            elif opt == '--species':EventAnnot=arg
            elif opt == '--platform':platform=arg
            
            ### The below options are specifically for splice-ICGS
            
            elif opt == '--forceBroadClusters':
                if 'rue' in string.lower(arg) or 'yes' in string.lower(arg):
                    forceBroadClusters = True
            elif opt == '--row_method':
                row_method=arg
                if row_method == 'None': row_method = None
            elif opt == '--column_method':
                column_method=arg
                if column_method == 'None': column_method = None
            elif opt == '--row_metric': row_metric=arg
            elif opt == '--column_metric': column_metric=arg
            elif opt == '--ExpressionCutoff': ExpressionCutoff=arg
            elif opt == '--normalization': normalization=arg
            elif opt == '--rho': rho_cutoff=float(arg)
            elif opt == '--CountsCutoff':CountsCutoff=int(float(arg))
            elif opt == '--FoldDiff':FoldDiff=float(arg)
            elif opt == '--SamplesDiffering':SamplesDiffering=int(float(arg))
            elif opt == '--removeOutliers':
                removeOutliers=arg
                if removeOutliers=='yes' or removeOutliers=='True':
                    removeOutliers = True
    
    print "Subtype discovery stringency:",strategy
    root_dir = export.findParentDir(EventAnnot)

    if EnrichmentOnly==False:
        
        print 'PSI input files:',EventAnnot
        print 'Using a rho-cutoff of:',rho_cutoff
    
        if filters==True: ### Filter based on a default percentage of samples with detected PSI values
            EventAnn,SampleNumber=filterPSIValues(EventAnnot,percentCutoff=percentCutoff,filterStatus=True)
        else:
            SampleNumber=filterPSIValues(EventAnnot,percentCutoff=percentCutoff,filterStatus=False)
        output_dir = root_dir+'ExpressionInput'
    
        export.createExportFolder(output_dir)
        full_PSI_InputFile=output_dir+"/exp.input.txt"
        header=header_list(EventAnnot)
        sampleIndexSelection.filterFile(EventAnnot,full_PSI_InputFile,header)
        
        ### Set Splice-ICGS defaults
        gsp = UI.GeneSelectionParameters(species,platform,platform)
        gsp.setNormalize('median')
        gsp.setGeneSelection('')
        gsp.setGeneSet('None Selected')
        gsp.setPathwaySelect([])
        gsp.setJustShowTheseIDs('')
        dynamicCorrelation = False
        gsp.setSampleDiscoveryParameters(ExpressionCutoff,CountsCutoff,FoldDiff,SamplesDiffering,dynamicCorrelation,
            removeOutliers,featurestoEvaluate,restrictBy,excludeCellCycle,column_metric,column_method,rho_cutoff)
        
        AnalysisRound=1
        if mode == "single":
            """ Perform a single round of Splice-ICGS (RNASeq.py module) """
            flag,full_PSI_InputFile,EventAnnot=CompleteWorkflow(full_PSI_InputFile,EventAnnot,rho_cutoff,strategy,seq,gsp,forceBroadClusters,AnalysisRound)
      
        else:
            """ Optionally iterate through multiple rounds of Splice-ICGS (RNASeq.py module) """
            while flag:
                if AnalysisRound != 1:
                    forceBroadClusters = False
                flag,full_PSI_InputFile,EventAnnot=CompleteWorkflow(full_PSI_InputFile,EventAnnot,rho_cutoff,strategy,seq,gsp,forceBroadClusters,AnalysisRound)
                AnalysisRound+=1
                if AnalysisRound>3:
                    flag = False
                if flag == False:
                    break
        
        output_dir = root_dir+'SVMOutputs'
        Combinedres=MergeResults(output_dir)
    else:
        output_dir = root_dir+'SVMOutputs'
        Combinedres=os.path.join(output_dir,"MergedResult.txt")
    
    ### Peform enrichment anlayses for each subtype against user-supplied clusters
    mutlabels={}
    if Mutationref!="":
        PSI_sample_headers = getPSIHeader(EventAnnot)
        print "Running Mutation Enrichment Analyses"
        Expand="yes"
        Mutationref = checkmetaDataFormat(Mutationref)
        mutdict=defaultdict(list)
        #print Mutationref
        print metaDataMatrixFormat
        header=ME.returnSamplesInMetaData(Mutationref,metaDataMatrixFormat=metaDataMatrixFormat)
        #print len(header)
        mutdict=ME.findsiggenepermut(Mutationref,valid_filenames=PSI_sample_headers)
        mutlabels=ME.Enrichment(Combinedres,mutdict,Mutationref,metaDataMatrixFormat,header)

    print "Generating the final consolidated results"
    Orderedheatmap.Classify(Combinedres,mutlabels,root_dir)
    Orderedheatmap.Classify(Combinedres,mutlabels,root_dir,False) 
    
    print "successfully completed"
