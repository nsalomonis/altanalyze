### Based on code from AltAnalyze's LineageProfiler (http://altanalyze.org)
#Author Nathan Salomonis - nsalomonis@gmail.com

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


"""
This script iterates the LineageProfiler algorithm (correlation based classification method) to identify sample types relative to one
of two references given one or more gene models. The main function is runLineageProfiler.

The program performs the following actions:
    1) Import a tab-delimited reference expression file with three columns (ID, biological group 1, group 2) and a header row (biological group names)
    2) Import a tab-delimited expression file with gene IDs (column 1), sample names (row 1) and normalized expression values (e.g., delta CT values)
    3) (optional - import existing models) Import a tab-delimited file with comma delimited gene-models for analysis
    4) (optional - find new models) Identify all possible combinations of gene models for a supplied model size variable (e.g., --s 7)
    5) Iterate through any supplied or identified gene models to obtain predictions for novel or known sample types
    6) Export prediction results for all analyzed models to the folder CellClassification.
    7) (optional) Print the top 20 scores and models for all possible model combinations of size --s
    
"""
    
import sys, string, shutil
import math
import os.path
import copy
import time
import getopt
try: import scipy
except Exception: pass
import traceback
import warnings
import random

try: import unique ### Not required (used in AltAnalyze)
except Exception: None
try: import export ### Not required (used in AltAnalyze)
except Exception: None

command_args = string.join(sys.argv,' ')
if len(sys.argv[1:])>1 and '-' in command_args and '--GUI' not in command_args:
    runningCommandLine = True
else:
    runningCommandLine = False
    
#from stats_scripts import salstat_stats; reload(salstat_stats)
try:
    from scipy import stats
    use_scipy = True
except Exception:
    use_scipy = False ### scipy is not required but is used as a faster implementation of Fisher Exact Test when present

def filepath(filename):
    try: fn = unique.filepath(filename)
    except Exception: fn = filename
    return fn

def exportFile(filename):
    try: export_data = export.ExportFile(filename)
    except Exception: export_data = open(filename,'w')
    return export_data

def makeUnique(item):
    db1={}; list1=[]; k=0
    for i in item:
        try: db1[i]=[]
        except TypeError: db1[tuple(i)]=[]; k=1
    for i in db1:
        if k==0: list1.append(i)
        else: list1.append(list(i))
    list1.sort()
    return list1

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def verifyFile(filename):
    status = True
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = True;break
    except Exception: status = False
    return status

def getHeader(filename):
    header=[]
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        header = string.split(data,'\t')[1:];break
    return header

def getDataMatrix(filename):
    header=[]
    matrix={}
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if len(header)==0:
            header = string.split(data,'\t')[1:]
        else:
            t = string.split(data,'\t')
            values = map(float,t[1:])
            matrix[t[0]]=values
    return header,matrix

def returnLargeGlobalVars():
    ### Prints all large global variables retained in memory (taking up space)
    all = [var for var in globals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all:
        try:
            if len(globals()[var])>1:
                print var, len(globals()[var])
        except Exception: null=[]
        
def clearObjectsFromMemory(db_to_clear):
    db_keys={}
    for key in db_to_clear: db_keys[key]=[]
    for key in db_keys:
        try: del db_to_clear[key]
        except Exception: 
            try:
                for i in key: del i ### For lists of tuples
            except Exception: del key ### For plain lists

def int_check(value):
    val_float = float(value)
    val_int = int(value)
    if val_float == val_int:
        integer_check = 'yes'
    if val_float != val_int:
        integer_check = 'no'
    return integer_check

def IQR(array):
    k1 = 75
    k2 = 25
    array.sort()
    n = len(array)
    value1 = float((n*k1)/100)
    value2 = float((n*k2)/100)
    if int_check(value1) == 'no':
        k1_val = int(value1) + 1
    if int_check(value1) == 'yes':
        k1_val = int(value1)
    if int_check(value2) == 'no':
        k2_val = int(value2) + 1
    if int_check(value2) == 'yes':
        k2_val = int(value2)
    try: median_val = scipy.median(array)
    except Exception: median_val = Median(array)
    upper75th = array[k1_val]
    lower25th = array[k2_val]
    
    int_qrt_range = upper75th - lower25th
    T1 = lower25th-(1.5*int_qrt_range)
    T2 = upper75th+(1.5*int_qrt_range)
    return lower25th,median_val,upper75th,int_qrt_range,T1,T2

class IQRData:
    def __init__(self,maxz,minz,medz,iq1,iq3):
        self.maxz = maxz; self.minz = minz
        self.medz = medz; self.iq1 = iq1
        self.iq3 = iq3
    def Max(self): return self.maxz
    def Min(self): return self.minz
    def Medium(self): return self.medz
    def IQ1(self): return self.iq1
    def IQ3(self): return self.iq3
    def SummaryValues(self):
        vals = string.join([str(self.IQ1()),str(self.Min()),str(self.Medium()),str(self.Max()),str(self.IQ3())],'\t')
        return vals
    
def importGeneModels(filename):
    x=0
    geneModels={}; thresholds=None
    #filename = None ### Override file import with default reference data (hard-coded)
    if filename != None:
        fn=filepath(filename)
        fileRead = open(fn,'rU').xreadlines()
    else:
        fileRead = defaultGeneModels()
        
    for line in fileRead:
        try:
            data = cleanUpLine(line)
            t = string.split(data,'\t')
        except Exception:
            t = line
        genes = t[0]
        genes = string.replace(genes,"'",'')
        genes = string.replace(genes,' ',',')
        genes = string.split(genes,',')
        if t>1:
            try: thresholds = map(float,t[1:])
            except Exception: thresholds = None
        try:
            if len(thresholds)==0: thresholds = None
        except Exception: pass
        models=[]
        for gene in genes:
            if len(gene)>0:
                models.append(gene)
        if len(models)>0:
            geneModels[tuple(models)] = thresholds
    return geneModels
                
######### Below code deals is specific to this module #########
def runLineageProfiler(species,array_type,exp_input,exp_output,codingtype,compendium_platform,modelSize=None,customMarkers=False,geneModels=False,permute=False,useMulti=False,fl=None):
    """ This code differs from LineageProfiler.py in that it is able to iterate through the LineageProfiler functions with distinct geneModels
    that are either supplied by the user or discovered from all possible combinations. """
    #global inputType
    global exp_output_file; exp_output_file = exp_output; global targetPlatform
    global tissues; global sample_headers; global collapse_override
    global analysis_type; global coding_type; coding_type = codingtype
    global tissue_to_gene; tissue_to_gene = {}; global platform; global cutoff
    global customMarkerFile; global delim; global keyed_by; global pearson_list
    global Permute; Permute=permute; global useMultiRef; useMultiRef = useMulti
    pearson_list={}
    cellHarmony=True
    #global tissue_specific_db
    
    try: returnCentroids = fl.ReturnCentroids()
    except Exception: returnCentroids = False
    
    print 'platform:',array_type
    
    collapse_override = True
    
    if 'ICGS' in customMarkers or 'MarkerGene' in customMarkers:
        """ When performing cellHarmony, build an ICGS expression reference with log2 TPM values rather than fold """
        print 'Converting ICGS folds to ICGS expression values as a reference first...'
        try: customMarkers = convertICGSClustersToExpression(customMarkers,exp_input,returnCentroids=returnCentroids,species=species)
        except:
            print "Using the supplied reference file only (not importing raw expression)...Proceeding without differential expression analsyes..."
            pass

    customMarkerFile = customMarkers
    if geneModels == False or geneModels == None: geneModels = []
    else:
        geneModels = importGeneModels(geneModels)

    exp_input = string.replace(exp_input,'\\','/')
    exp_input = string.replace(exp_input,'//','/')
    exp_output = string.replace(exp_output,'\\','/')
    exp_output = string.replace(exp_output,'//','/')
    
    delim = "/"
    
    print '\nRunning cellHarmony analysis on',string.split(exp_input,delim)[-1][:-4]
    
    global correlate_by_order; correlate_by_order = 'no'
    global rho_threshold; rho_threshold = -1
    global correlate_to_tissue_specific; correlate_to_tissue_specific = 'no'
    cutoff = 0.01
    global value_type
    platform = array_type
    
    if 'stats.' in exp_input:
        value_type = 'calls'
    else:
        value_type = 'expression'
    
    tissue_specific_db={}; sample_headers=[]; tissues=[]
    if len(array_type)==2:
        ### When a user-supplied expression is provided (no ExpressionOutput files provided - importGeneIDTranslations)
        array_type, vendor = array_type
        if ':' in vendor:
            vendor = string.split(vendor,':')[1]
        if array_type == 'RNASeq':
            vendor = 'Symbol'
        platform = array_type
    else: vendor = 'Symbol'

    if 'RawSplice' in exp_input or 'FullDatasets' in exp_input or coding_type == 'AltExon':
        analysis_type = 'AltExon'
        if platform != compendium_platform: ### If the input IDs are not Affymetrix Exon 1.0 ST probesets, then translate to the appropriate system
            translate_to_genearray = 'no'
            targetPlatform = compendium_platform
            translation_db = importExonIDTranslations(array_type,species,translate_to_genearray)
            keyed_by = 'translation'
        else: translation_db=[]; keyed_by = 'primaryID'; targetPlatform = compendium_platform
    else:
        try:
            ### Get arrayID to Ensembl associations
            if vendor != 'Not needed':
                ### When no ExpressionOutput files provided (user supplied matrix)
                translation_db = importVendorToEnsemblTranslations(species,vendor,exp_input)
            else:
                translation_db = importGeneIDTranslations(exp_output)
            keyed_by = 'translation'
            targetPlatform = compendium_platform
            analysis_type = 'geneLevel' 
        except Exception:
            translation_db=[]; keyed_by = 'primaryID'; targetPlatform = compendium_platform; analysis_type = 'geneLevel'

    targetPlatform = compendium_platform ### Overides above
    try: importTissueSpecificProfiles(species,tissue_specific_db)
    except Exception:
        try:
            try:
                targetPlatform = 'exon'
                importTissueSpecificProfiles(species,tissue_specific_db)
            except Exception:
                try:
                    targetPlatform = 'gene'
                    importTissueSpecificProfiles(species,tissue_specific_db)
                except Exception: 
                    targetPlatform = "3'array"
                    importTissueSpecificProfiles(species,tissue_specific_db)
        except Exception:
            print 'No compatible compendiums present...'
            print traceback.format_exc() 
            forceError

    gene_expression_db, sample_headers = importGeneExpressionValuesSimple(exp_input,translation_db)
    
    ### Make sure each sample ID is unique
    samples_added=[]; i=1
    ### Organize sample expression, with the same gene order as the tissue expression set
    for s in sample_headers:
        if s in samples_added:
            samples_added.append(s+"."+str(i)) ### Ensure only unique sample IDs exist
            i+=1
        else:
            samples_added.append(s)
    sample_headers = samples_added

    pruneTissueSpecific=False
    for gene in tissue_specific_db:
        if gene not in gene_expression_db:
            pruneTissueSpecific = True
            break
        
    if pruneTissueSpecific:
        tissue_specific_db2={}
        for gene in gene_expression_db:
            if gene in tissue_specific_db:
                tissue_specific_db2[gene] = tissue_specific_db[gene]
            elif gene in translation_db:
                altGeneID = translation_db[gene]
                if altGeneID in tissue_specific_db:
                    tissue_specific_db2[gene] = tissue_specific_db[altGeneID]
                    
        tissue_specific_db = tissue_specific_db2
    
    all_marker_genes=[]            
    for gene in tissue_specific_db:
        all_marker_genes.append(gene)
    #print [modelSize]
    
    if len(geneModels)>0:
        allPossibleClassifiers = geneModels
    elif modelSize == None or modelSize == 'optimize' or modelSize == 'no':
        allPossibleClassifiers={}
        allPossibleClassifiers[tuple(all_marker_genes)]=None
    else:
        ### A specific model size has been specified (e.g., find all 10-gene models)
        allPossibleClassifiers = getRandomSets(all_marker_genes,modelSize)
                
    num=1
    all_models=[]
    if len(allPossibleClassifiers)<16:
        print 'Using:'
        for model in allPossibleClassifiers:
            print 'model',num, 'with',len(model),'genes'  #model
            num+=1
            all_models+=model
    
    #all_models = unique.unique(all_models)
    #print len(all_models);sys.exit()
    
    ### This is the main analysis function
    print 'Number of reference samples to compare to:',len(tissues)
    if len(tissues)<16:
        print tissues

    if modelSize != 'optimize':
        hit_list, hits, fails, prognostic_class_db,sample_diff_z, evaluate_size, prognostic_class1_db, prognostic_class2_db = iterateLineageProfiler(exp_input,
                            tissue_specific_db, allPossibleClassifiers,translation_db,compendium_platform,modelSize,species,gene_expression_db, sample_headers)
    else:
        summary_hit_list=[]
        try: evaluate_size = len(allPossibleClassifiers[0])
        except Exception:
            ### Occurs when custom models loaded
            for i in allPossibleClassifiers:
                evaluate_size = len(i); break
        hit_list, hits, fails, prognostic_class_db,sample_diff_z, evaluate_size, prognostic_class1_db, prognostic_class2_db = iterateLineageProfiler(exp_input,
                            tissue_specific_db, allPossibleClassifiers,translation_db,compendium_platform,None,species,gene_expression_db, sample_headers)
        while evaluate_size > 4:
            hit_list.sort()
            top_model = hit_list[-1][-1]
            top_model_score = hit_list[-1][0]
            """
            try: ### Used for evaluation only - gives the same top models
                second_model = hit_list[-2][-1]
                second_model_score = hit_list[-2][0]
                if second_model_score == top_model_score:
                    top_model = second_model_score ### Try this
                    print 'selecting secondary'
            except Exception: None
            """
            allPossibleClassifiers = [hit_list[-1][-1]]
            
            hit_list, hits, fails, prognostic_class_db,sample_diff_z, evaluate_size, prognostic_class1_db, prognostic_class2_db = iterateLineageProfiler(exp_input,
                            tissue_specific_db, allPossibleClassifiers,translation_db,compendium_platform,modelSize,species,gene_expression_db, sample_headers)
            summary_hit_list+=hit_list
        hit_list = summary_hit_list
    
    root_dir = string.join(string.split(exp_output_file,'/')[:-1],'/')+'/'
    dataset_name = string.replace(string.split(exp_input,'/')[-1][:-4],'exp.','')
    output_classification_file = root_dir+'CellClassification/'+dataset_name+'-CellClassification.txt'
    try: os.mkdir(root_dir+'CellClassification')
    except Exception: None
    export_summary = exportFile(output_classification_file)
    models = []
    for i in allPossibleClassifiers:
        i = string.replace(str(i),"'",'')[1:-1]
        models.append(i)

    ### If multiple class-headers with a common phenotype (different source), combine into a single report
    class_list=[]
    processed=0
    for h in tissues:
        if ':' in h and collapse_override==False:
            try:
                phenotype, source = string.split(h,':')
                processed+=1
                if phenotype not in class_list:
                    class_list.append(phenotype)
            except Exception: pass
    
    
    if len(class_list)==2 and len(tissues) == processed and collapse_override==False and cellHarmony==False: ### Ensures all reference headers have : in them
        tissue_list = class_list
        collapse = True
    else:
        tissue_list = tissues
        collapse = False
    
    print ''
    class_headers = map(lambda x: x+' Predicted Hits',tissue_list)
    headers = string.join(['Samples']+class_headers+['Composite Classification Score','Combined Correlation DiffScore','Predicted Class','Max-Rho'],'\t')+'\n'
    export_summary.write(headers)
    
    sorted_results=[] ### sort the results
    try: numberOfModels = len(allPossibleClassifiers)
    except Exception: numberOfModels = 1

    accuracy=[]
    ar=[]
    non=[]
    no_intermediate_accuracy=[]
    verboseReport = False
    for sample in prognostic_class_db:
        if len(tissues)==2:
            class1_score = prognostic_class1_db[sample]
            class2_score = prognostic_class2_db[sample]
        zscore_distribution = map(lambda x: str(x[0]), sample_diff_z[sample])
        pearson_max_values = map(lambda x: str(x[1]), sample_diff_z[sample])
        dist_list=[]
        for i in zscore_distribution:
            try: dist_list.append(float(i))
            except Exception: None ### Occurs for 'NA's
        #try: median_score = scipy.median(dist_list)
        #except Exception: median_score = Median(dist_list)
        try: sum_score = sum(dist_list) #scipy.median
        except Exception: sum_score = sum(dist_list)

        correlations=[]
        for i in pearson_max_values:
            try: correlations.append(float(i))
            except Exception: None ### Occurs for 'NA's
        median_correlation = scipy.median(correlations)
        if median_correlation<0.8 and verboseReport:
            print 'Sample: %s has a low median model Pearson correlation coefficient (%s)' % (sample,str(median_correlation))
        if verboseReport==False:
            print '.',
        class_db = prognostic_class_db[sample]
        class_scores=[]; class_scores_str=[]; class_scores_refs=[]; collapsed_pheno_scores={}
        
        for tissue in tissues:
            if collapse and collapse_override==False:
                phenotype,source = string.split(tissue,':')
                try: collapsed_pheno_scores[phenotype]+=class_db[tissue]
                except Exception: collapsed_pheno_scores[phenotype]=class_db[tissue]
            else:
                class_scores_str.append(str(class_db[tissue]))
                class_scores.append(class_db[tissue])
                class_scores_refs.append((class_db[tissue],tissue))
        if collapse and collapse_override == False:
            for phenotype in tissue_list:
                class_scores.append(collapsed_pheno_scores[phenotype]) ### Collapse the scores and report in the original phenotype order
                class_scores_str.append(str(collapsed_pheno_scores[phenotype]))
                class_scores_refs.append((collapsed_pheno_scores[phenotype],phenotype))
        
        """
        for tissue in tissues:
            class_scores_str.append(str(class_db[tissue]))
            class_scores.append(class_db[tissue])
            class_scores_refs.append((class_db[tissue],tissue))
        """
        overall_prog_score = str(max(class_scores)-min(class_scores))
        if len(tissues)==2 and cellHarmony==False:
            class_scores_str = [str(class1_score),str(class2_score)] ### range of positive and negative scores for a two-class test
            if class1_score == 0 and class2_score == 0:
                call = 'Intermediate Risk '+ tissues[0]
            elif class1_score == numberOfModels:
                call = 'High Risk '+ tissues[0]
            elif class2_score == numberOfModels:
                call = 'Low Risk '+ tissues[0]
            elif class1_score == 0:
                call = 'Intermediate Risk '+ tissues[0]
            elif class2_score == 0:
                call = 'Intermediate Risk '+ tissues[0]
            else:
                call = 'Itermediate Risk '+ tissues[0]
            overall_prog_score = str(class1_score-class2_score)
        else:
            class_scores_refs.sort()
            call=class_scores_refs[-1][1] ### This is the reference with the max reported score
            if call == tissue_list[-1]: ### Usually the classifier of interest should be listed first in the reference file, not second
                overall_prog_score = str(float(overall_prog_score)*-1)
                sum_score = sum_score*-1
        values = [sample]+class_scores_str+[overall_prog_score,str(sum_score),call]
        values = string.join(values+zscore_distribution[:-1]+[str(max(correlations))],'\t')+'\n'
        if ':' in sample:
            sample = string.split(sample,':')[0]
        if ':' in call:
            call = string.split(call,':')[0]
        if call==sample:
            accuracy.append(float(1))
            if float(overall_prog_score) > 10 or float(overall_prog_score) < -10:
                no_intermediate_accuracy.append(float(1))
                if 'non' in call: non.append(float(1))
                else: ar.append(float(1))
        else:
            accuracy.append(float(0))
            if float(overall_prog_score) > 10 or float(overall_prog_score) < -10:
                no_intermediate_accuracy.append(float(0))
                if 'non' in call: non.append(float(0))
                else: ar.append(float(0))
        sorted_results.append([float(overall_prog_score),sum_score,values])
        sample_diff_z[sample] = dist_list

    if verboseReport:
        print len(no_intermediate_accuracy)
        print no_intermediate_accuracy
        print 'Overall Acuracy:',Average(accuracy)*100
        print 'Sensititivity:', sum(ar), len(ar)
        print 'Specificity:', sum(non), len(non)
        print str(Average(accuracy)*100)+'\t'+str(Average(ar)*100)+'\t'+str(Average(non)*100)+'\t'+str(Average(no_intermediate_accuracy)*100)+'\t'+str(sum(ar))+'\t'+str(len(ar))+'\t'+str(sum(non))+'\t'+str(len(non))
    else:
        print '\nClassification analysis completed...'
    sorted_results.sort()
    sorted_results.reverse()
    for i in sorted_results:
        export_summary.write(i[-1])

    export_summary.close()
    print '\nResults file written to:',root_dir+'CellClassification/'+dataset_name+'-CellClassification.txt','\n'

    hit_list.sort(); hit_list.reverse()
    top_hit_list=[]
    top_hit_db={}
    hits_db={}; fails_db={}
    
    ### Only look at the max correlation for each sample
    max_pearson_list=[]
    for sample in pearson_list:
        pearson_list[sample].sort()
        for rho in pearson_list[sample][-2:]: ### get the top two correlations
            max_pearson_list.append(rho)
    
    avg_pearson_rho = Average(max_pearson_list)
    maxPearson = max(max_pearson_list)

    try:
        for i in sample_diff_z:
            zscore_distribution = sample_diff_z[i]
            maxz = max(zscore_distribution); minz = min(zscore_distribution)
            sample_diff_z[i] = string.join(map(str,zscore_distribution),'\t')   
            try:
                lower25th,medz,upper75th,int_qrt_range,T1,T2 = IQR(zscore_distribution)
                if float(maxz)>float(T2): maxz = T2
                if float(minz) < float(T1): minz = T1
                #iqr = IQRData(maxz,minz,medz,lower25th,upper75th)
                #sample_diff_z[i] = iqr
            except Exception:
                pass
        for i in hits:
            try: hits_db[i]+=1
            except Exception: hits_db[i]=1
        for i in fails:
            try: fails_db[i]+=1
            except Exception: fails_db[i]=1
        for i in fails_db:
            if i not in hits:
                try:
                    #print i+'\t'+'0\t'+str(fails_db[i])+'\t'+ sample_diff_z[i]
                    None
                except Exception:
                    #print i
                    None    
    except Exception:
        pass
    
    exportModelScores = True
    if modelSize != False:
        #print 'Returning all model overall scores'
        hits=[]
        for i in hits_db:
            hits.append([hits_db[i],i])
        hits.sort()
        hits.reverse()
        for i in hits:
            if i[1] in fails_db: fail = fails_db[i[1]]
            else: fail = 0
            try:
                #print i[1]+'\t'+str(i[0])+'\t'+str(fail)+'\t'+sample_diff_z[i[1]]
                None
            except Exception:
                #print i[1]
                None
        if modelSize == 'optimize': threshold = 80
        else: threshold = 0
        #print 'threshold:',threshold
        for i in hit_list:
            if i[0]>threshold:
                top_hit_list.append(i[-1])
                top_hit_db[tuple(i[-1])]=i[0]
        
        if len(geneModels) > 0 and exportModelScores==False:
            for i in hit_list:
                #print i[:5],i[-1],i[-2] ### print all
                pass
        else:
            """
            print 'Returning all over 90'
            for i in hit_list:
                if i[0]>85:
                    print i[:5],i[-1],i[-2] ### print all
                
            sys.exit()"""
            #print 'Top hits'
            output_model_file = root_dir+'CellClassification/'+dataset_name+'-ModelScores.txt'
            export_summary = exportFile(output_model_file)
            print 'Exporting top-scoring models to:',output_model_file
            title = 'Classification-Rate\tClass1-Hits\tClass1-Total\tClass2-Hits\tClass2-Total\tModel\tModel-Gene-Number\n'
            export_summary.write(title)
            
            for i in hit_list: #hit_list[:500]
                overall_scores=[]
                for x in i[:5]: overall_scores.append(str(x))
                model = string.replace(str(i[-1])[1:-1],"'",'')
                values = string.join(overall_scores+[model]+[str(i[-2])],'\t')+'\n'
                export_summary.write(values)
            export_summary.close()
            """
            try:
                if hit_list[0][0] == hit_list[20][0]:
                    for i in  hit_list[20:]:
                        if hit_list[0][0] == i[0]:
                            print i[:5],i[-1],i[-2]
                        else: sys.exit()
            except Exception: None ### Occurs if less than 20 entries here
            """

    print 'Average Pearson correlation coefficient:', avg_pearson_rho
    if avg_pearson_rho<0.9 and verboseReport:
        print '\n\nWARNING!!!!!!!!!'
        print '\tThe average Pearson correlation coefficient for all example models is less than 0.9.'
        print '\tYour data may not be comparable to the provided reference (quality control may be needed).\n\n'
    elif verboseReport:
        print 'No unusual warning.\n'
        
    reference_exp_file = customMarkers
    query_exp_file = exp_input
    classification_file = output_classification_file
    harmonizeClassifiedSamples(species,reference_exp_file, query_exp_file, classification_file,fl=fl)
    
    return top_hit_db


def iterateLineageProfiler(exp_input,tissue_specific_db,allPossibleClassifiers,translation_db,compendium_platform,
                modelSize,species,gene_expression_db,sampleHeaders):

    classifyBasedOnRho=True
    hit_list=[]
    ### Iterate through LineageProfiler for all gene models (allPossibleClassifiers)
    times = 1; k=1000; l=1000; hits=[]; fails=[]; f=0; s=0; sample_diff_z={}; prognostic_class1_db={}; prognostic_class2_db={}
    prognostic_class_db={}
    begin_time = time.time()
    
    try: evaluate_size=len(allPossibleClassifiers[0]) ### Number of reference markers to evaluate
    except Exception:
        for i in allPossibleClassifiers: evaluate_size = len(i); break
    if modelSize=='optimize':
        evaluate_size -= 1
        allPossibleClassifiers = getRandomSets(allPossibleClassifiers[0],evaluate_size)

    ### Determine if we should collapse the entries or not based on common phenotype references
    class_list=[]; processed=0; alternate_class={}
    for h in tissues:
        if ':' in h and 'ENS' not in h:
            try:
                phenotype, source = string.split(h,':'); processed+=1
                if phenotype not in class_list: class_list.append(phenotype)
            except Exception: pass

    try:
        alternate_class[class_list[0]] = class_list[1]
        alternate_class[class_list[1]] = class_list[0]
    except Exception: pass
    
    if len(class_list)==2 and len(tissues) == processed and collapse_override==False: ### Ensures all reference headers have : in them
        tissue_list = class_list; collapse = True
    else: tissue_list = tissues; collapse = False
    
    cellHarmony=True
    mean_percent_positive=[]
    for classifiers in allPossibleClassifiers:
        try: thresholds = allPossibleClassifiers[classifiers]
        except Exception: thresholds = None
        tissue_to_gene={}; expession_subset=[]; sample_headers=[]; classifier_specific_db={}
        for gene in classifiers:
            try: classifier_specific_db[gene] = tissue_specific_db[gene]
            except Exception: None
        #print len(gene_expression_db), len(classifier_specific_db), len(expession_subset), len(translation_db)
        expession_subset = filterGeneExpressionValues(gene_expression_db,classifier_specific_db,translation_db,expession_subset)

        ### If the incorrect gene system was indicated re-run with generic parameters
        if len(expession_subset)==0:
            translation_db=[]; keyed_by = 'primaryID'; targetPlatform = compendium_platform; analysis_type = 'geneLevel'
            tissue_specific_db={}
            importTissueSpecificProfiles(species,tissue_specific_db)
            expession_subset = filterGeneExpressionValues(gene_expression_db,tissue_specific_db,translation_db,expession_subset)
        if len(sample_diff_z)==0: ### Do this for the first model examine only
            for h in sampleHeaders:
                sample_diff_z[h]=[] ### Create this before any data is added, since some models will exclude data for some samples (missing dCT values)
        if len(expession_subset)!=len(classifiers): f+=1
        #if modelSize=='optimize': print len(expession_subset), len(classifiers);sys.exit()
        if (len(expession_subset) != len(classifiers)) and  modelSize=='optimize':
            print "Please provide a reference set of equal length or smaller to the input analysis set"; kill
        #print len(expession_subset), len(classifiers);sys.exit()
        if len(expession_subset)==len(classifiers): ### Sometimes a gene or two are missing from one set
            s+=1
            #print classifiers,'\t',
            zscore_output_dir,tissue_scores,sampleHeaders = analyzeTissueSpecificExpressionPatterns(tissue_specific_db,expession_subset,sampleHeaders)
            #except Exception: print len(classifier_specific_db), classifiers; error
            headers = list(tissue_scores['headers']); del tissue_scores['headers']
            if times == k:
                end_time = time.time()
                print int(end_time-begin_time),'seconds'
                k+=l
            times+=1; index=0; positive=0; positive_score_diff=0
            sample_number = (len(headers)-1)
            population1_denom=0; population1_pos=0; population2_pos=0; population2_denom=0
            diff_positive=[]; diff_negative=[]
            while index < sample_number:
                ### The scores are now a tuple of (Z-score,original_pearson_rho)
                scores = map(lambda x: tissue_scores[x][index], tissue_scores)
                zscores_only = map(lambda x: tissue_scores[x][index][0], tissue_scores)
                scores_copy = list(scores); scores_copy.sort()
                max_pearson_model = scores_copy[-1][1]  ### e.g., tumor rho (this is the one we want to QC on later)
                min_pearson_model = scores_copy[-2][1] ### e.g., non-tumor rho
                diff_rho = max_pearson_model - min_pearson_model
                diff_z = (scores_copy[-1][0]-scores_copy[-2][0])*100 ### Diff between the top two scores (z-scores are the first item)
                if classifyBasedOnRho == True:
                    diff_z = diff_rho*10
                positive_class=None
                j=0
                for tissue in tissue_scores:
                    if ':' in tissue and 'ENS' not in tissue:
                        group_name = string.split(tissue,':')[0]
                    else:
                        group_name = tissue
                    if scores[j][0] == max(zscores_only):
                        hit_score = 1; positive_class = tissue
                    else: hit_score = 0
                    if len(tissues)>2 or cellHarmony==True:
                        if group_name+':' in headers[index+1] and hit_score==1:
                            g = string.split(headers[index+1],':')[0]+':'
                            if g in group_name+':': ### reciprocol of above
                                positive+=1
                    try:
                        class_db = prognostic_class_db[headers[index+1]]
                        try: class_db[tissue]+=hit_score
                        except Exception: class_db[tissue]=hit_score
                    except Exception:
                        class_db={}
                        class_db[tissue]=hit_score
                        prognostic_class_db[headers[index+1]] = class_db
                    j+=1
                if collapse and collapse_override==False:
                    phenotype, source = string.split(positive_class,':')
                    baseline_positive_class = alternate_class[phenotype]+':'+source
                    denom_score,denom_rho = tissue_scores[baseline_positive_class][index]
                    old_diff = diff_z
                    diff_z = scores_copy[-1][0]-denom_score ### Diff between the top two scores of the SAME SOURCE
                    diff_rho = (max_pearson_model - denom_rho)
                    if classifyBasedOnRho == True:
                        diff_z = diff_rho*10
                    #print headers[index+1], scores_copy[-1]
                if len(tissues)==2:
                    if ':' in headers[index+1]:
                        pheno = string.split(headers[index+1],':')[0]
                    else:
                        pheno = None
                    diff_z = tissue_scores[tissues[0]][index][0]-tissue_scores[tissues[-1]][index][0] ### z-scores are the first item and pearson-rho is the second
                    diff_rho = (tissue_scores[tissues[0]][index][1]-tissue_scores[tissues[-1]][index][1])
                    if classifyBasedOnRho == True:
                        diff_z = diff_rho*10
                        
                    if thresholds == None:
                        threshold1 = 0; threshold2 = 0
                    else:
                        threshold1, threshold2 = thresholds ### emperically derived cutoffs provided by the user for each model (e.g., mean+2SD of diff_rho)
                        
                    if headers[index+1] not in prognostic_class1_db:
                        prognostic_class1_db[headers[index+1]]=0 ### Create a default value for each sample
                    if headers[index+1] not in prognostic_class2_db:
                        prognostic_class2_db[headers[index+1]]=0 ### Create a default value for each sample
                    if diff_z>threshold1:
                        prognostic_class1_db[headers[index+1]]+=1
                    elif diff_z<threshold2:
                        prognostic_class2_db[headers[index+1]]+=1
                    if diff_z>0 and (tissues[0] == pheno):
                        positive+=1; positive_score_diff+=abs(diff_z)
                        population1_pos+=1; diff_positive.append(abs(diff_z))
                        hits.append(headers[index+1]) ### see which are correctly classified
                    elif diff_z<0 and (tissues[-1] == pheno):
                        positive+=1; positive_score_diff+=abs(diff_z)
                        population2_pos+=1; diff_positive.append(abs(diff_z))
                        hits.append(headers[index+1]) ### see which are correctly classified
                    elif diff_z>0 and (tissues[-1] == pheno): ### Incorrectly classified
                        diff_negative.append(abs(diff_z))
                        fails.append(headers[index+1])
                    elif diff_z<0 and (tissues[0] == pheno): ### Incorrectly classified
                        #print headers[index+1]
                        diff_negative.append(abs(diff_z))
                        fails.append(headers[index+1])
                    if (tissues[0] == pheno):
                        population1_denom+=1
                    else:
                        population2_denom+=1
                sample_diff_z[headers[index+1]].append((diff_z,max([max_pearson_model,min_pearson_model])))  ### Added pearson max here
                index+=1
            try: percent_positive = (float(positive)/float(index))*100
            except ZeroDivisionError:
                print 'WARNING!!!! No matching genes. Make sure your gene IDs are the same ID type as your reference.'
                forceNoMatchingGeneError
            mean_percent_positive.append(percent_positive)
            if len(tissues)==2 and cellHarmony==False:
                try:
                    pos = float(population1_pos)/float(population1_denom)
                    neg = float(population2_pos)/float(population2_denom)
                    #percent_positive = (pos+neg)/2
                except Exception: pos = 0; neg = 0
                hit_list.append([percent_positive,population1_pos, population1_denom,population2_pos,population2_denom,[Average(diff_positive),Average(diff_negative)],positive_score_diff,len(classifiers),classifiers])
            else:
                hit_list.append([percent_positive,len(classifiers),classifiers])
            for sample in sample_diff_z:
                if len(sample_diff_z[sample]) != (times-1): ### Occurs when there is missing data for a sample from the analyzed model
                    sample_diff_z[sample].append(('NA','NA')) ### add a null result
    #print Average(mean_percent_positive), '\tAverage'
    return hit_list, hits, fails, prognostic_class_db, sample_diff_z, evaluate_size, prognostic_class1_db, prognostic_class2_db

def factorial(n):
    ### Code from http://docs.python.org/lib/module-doctest.html
    if not n >= 0:
        raise ValueError("n must be >= 0")
    if math.floor(n) != n:
        raise ValueError("n must be exact integer")
    if n+1 == n:  # catch a value like 1e300
        raise OverflowError("n too large")
    result = 1
    factor = 2
    while factor <= n:
        result *= factor
        factor += 1
    return result

def choose(n,x):
    """Equation represents the number of ways in which x objects can be selected from a total of n objects without regard to order."""
    #(n x) = n!/(x!(n-x)!)
    f = factorial
    result = f(n)/(f(x)*f(n-x))
    return result

def getRandomSets(a,size):
    #a = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    #size = 4
    select_set={'ENSG00000140678':'ITGAX','ENSG00000105835':'NAMPT','ENSG00000027697':'IFNGR1','ENSG00000120129':'DUSP1','ENSG00000003402':'CFLAR','ENSG00000113269':'RNF130'}
    select_set={}
    
    select_set2={'ENSG00000163602': 'RYBP'}
    negative_select = {'ENSG00000105352':'CEACAM4'}
    negative_select={}

    possible_sets = choose(len(a),size)
    print 'Possible',size,'gene combinations to test',possible_sets
    permute_ls = []; done = 0; permute_db={}
    while done == 0:
        b = list(tuple(a)); random.shuffle(b)
        bx_set={}
        i = 0
        while i < len(b):
            try:
                bx = b[i:i+size]; bx.sort()
                if len(bx)==size: permute_db[tuple(bx)]=None
                else: break
            except Exception: break
            i+=1
        if len(permute_db) == possible_sets:
            done=1; break
    for i in permute_db:
        add=0; required=0; exclude=0
        for l in i:
            if len(select_set)>0:
                if l in select_set: add+=1
                #if l in select_set2: required+=1
            #if l in negative_select: exclude+=1
            else: add = 1000
        if add>2 and exclude==0:# and required==1:
            permute_ls.append(i)
    #print len(permute_ls)
    return permute_ls

def importVendorToEnsemblTranslations(species,vendor,exp_input):
    translation_db={}
    """
    ### Faster method but possibly not as good
    uid_db = simpleUIDImport(exp_input)
    import gene_associations
    ### Use the same annotation method that is used to create the ExpressionOutput annotations
    array_to_ens = gene_associations.filterGeneToUID(species,'Ensembl',vendor,associated_IDs)
    for arrayid in array_to_ens:
        ensembl_list = array_to_ens[arrayid]
        try: translation_db[arrayid] = ensembl_list[0] ### This first Ensembl is ranked as the most likely valid based on various metrics in getArrayAnnotationsFromGOElite
        except Exception: None
    """
    translation_db={}
    
    ### Use the same annotation method that is used to create the ExpressionOutput annotations
    use_go = 'yes'
    conventional_array_db={}
    import gene_associations
    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-'+vendor))
    
    for EnsGeneID in gene_to_symbol:
        vendorID = gene_to_symbol[EnsGeneID][0]
        try:
            translation_db[vendorID] = EnsGeneID
            translation_db[EnsGeneID] = vendorID
        except Exception: None
    
    return translation_db

def importTissueSpecificProfiles(species,tissue_specific_db):
    if analysis_type == 'AltExon':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_'+targetPlatform +'_tissue-specific_AltExon_protein_coding.txt'
    else:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_'+targetPlatform +'_tissue-specific_'+coding_type+'.txt'
    if customMarkerFile != False:
        filename = customMarkerFile

    if value_type == 'calls':
        filename = string.replace(filename,'.txt','_stats.txt')
    fn=filepath(filename); x=0
    tissues_added={}
    dataType = 'matrix'
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        
        if x==0:
            #print 'Importing the tissue compedium database:',string.split(filename,delim)[-1][:-4]
            headers = t; x=1; index=0
            for i in headers:
                if 'UID' == i or 'uid' == i: ens_index = index; uid_index = index
                if analysis_type == 'AltExon': ens_index = ens_index ### Assigned above when analyzing probesets
                elif 'Ensembl' in i: ens_index = index
                if 'marker-in' in i:
                    tissue_index = index+1; marker_in = index
                    dataType = 'MarkerFinder'
                if 'row_clusters-flat' in i:
                    tissue_index = 2 ### start reading fold data at column 3
                    dataType = 'clustering'
                index+=1
            try:
                for i in t[tissue_index:]: tissues.append(i)
            except Exception:
                for i in t[1:]: tissues.append(i)
            if keyed_by == 'primaryID':
                try: ens_index = uid_index
                except Exception: None
        else:
            if 'column_clusters-flat' == t[0]: ### clustered heatmap - skip this row
                continue
            if dataType=='clustering':
                gene = t[0]
                tissue_exp = map(float, t[2:])
                tissue_specific_db[gene]=x,tissue_exp
            else:  
                try:
                    gene = t[0]
                    tissue_exp = map(float, t[1:])
                    tissue_specific_db[gene]=x,tissue_exp ### Use this to only grab relevant gene expression profiles from the input dataset
                except Exception:
                    gene = string.split(t[ens_index],'|')[0] ### Only consider the first listed gene - this gene is the best option based on ExpressionBuilder rankings
                    #if 'Pluripotent Stem Cells' in t[marker_in] or 'Heart' in t[marker_in]:
                    #if t[marker_in] not in tissues_added: ### Only add the first instance of a gene for that tissue - used more for testing to quickly run the analysis
                    tissue_exp = map(float, t[tissue_index:])
                    if value_type == 'calls':
                        tissue_exp = produceDetectionCalls(tissue_exp,platform) ### 0 or 1 calls
                    
                    tissue_specific_db[gene]=x,tissue_exp ### Use this to only grab relevant gene expression profiles from the input dataset
                    tissues_added[t[marker_in]]=[]
            x+=1
    #print len(tissue_specific_db), 'genes in the reference database'

    if correlate_to_tissue_specific == 'yes':
        try: importTissueCorrelations(filename)
        except Exception:
            null=[]
            #print '\nNo tissue-specific correlations file present. Skipping analysis.'; kill
    return tissue_specific_db

def importTissueCorrelations(filename):
    filename = string.replace(filename,'specific','specific_correlations')
    fn=filepath(filename); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x==0: x=1 ### Ignore header line
        else:
            uid,symbol,rho,tissue = string.split(data,'\t')
            if float(rho)>rho_threshold: ### Variable used for testing different thresholds internally
                try: tissue_to_gene[tissue].append(uid)
                except Exception: tissue_to_gene[tissue] = [uid]

def simpleUIDImport(filename):
    """Import the UIDs in the gene expression file"""
    uid_db={}
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        uid_db[string.split(data,'\t')[0]]=[]
    return uid_db

def importGeneExpressionValuesSimple(filename,translation_db,translate=False):
    ### Import gene-level expression raw values           
    fn=filepath(filename); x=0; genes_added={}; gene_expression_db={}
    dataset_name = string.split(filename,delim)[-1][:-4]
    #print 'importing:',dataset_name
    sampleIndex=1
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            if 'row_clusters-flat' in t:
                sampleIndex = 2 ### start reading fold data at column 3
                dataType = 'clustering'
            if '#' not in data[0]:
                for i in t[sampleIndex:]: sample_headers.append(i)
                x=1
        else:
            gene = t[0]
            if 'column_clusters-flat' == t[0]: ### clustered heatmap - skip this row
                continue
            #if '-' not in gene and ':E' in gene: print gene;sys.exit()
            if analysis_type == 'AltExon':
                try: ens_gene,exon = string.split(gene,'-')[:2]
                except Exception: exon = gene
                gene = exon
            if keyed_by == 'translation': ### alternative value is 'primaryID'
                """if gene == 'ENSMUSG00000025915-E19.3':
                    for i in translation_db: print [i], len(translation_db); break
                    print gene, [translation_db[gene]];sys.exit()"""
                if translate:
                    ### This caused two many problems so we removed it with the translate default above
                    try: gene = translation_db[gene] ### Ensembl annotations
                    except Exception: pass
                
            try: genes_added[gene]+=1
            except Exception: genes_added[gene]=1
            try: exp_vals = map(float, t[sampleIndex:])
            except Exception:
                ### If a non-numeric value in the list
                exp_vals=[]
                for i in t[sampleIndex:]:
                    try: exp_vals.append(float(i))
                    except Exception: exp_vals.append(i)
            gene_expression_db[gene] = exp_vals
    #print len(gene_expression_db), 'matching genes in the dataset and tissue compendium database'
    
    return gene_expression_db, sample_headers

def filterGeneExpressionValues(all_gene_expression_db,tissue_specific_db,translation_db,expession_subset):
    ### Filter all imported gene expression values
    gene_expression_db={}; genes_added={}
    for gene in all_gene_expression_db:
        exp_vals = all_gene_expression_db[gene]
        if gene in tissue_specific_db:
            index,tissue_exp=tissue_specific_db[gene]
            try: genes_added[gene]+=1
            except Exception: genes_added[gene]=1
            if value_type == 'calls': ### Hence, this is a DABG or RNA-Seq expression
                exp_vals = produceDetectionCalls(exp_vals,targetPlatform) ### 0 or 1 calls
            gene_expression_db[gene] = [index,exp_vals]
    #print len(gene_expression_db), 'matching genes in the dataset and tissue compendium database'
    
    for gene in genes_added:
        if genes_added[gene]>1: del gene_expression_db[gene] ### delete entries that are present in the input set multiple times (not trustworthy)
        else: expession_subset.append(gene_expression_db[gene]) ### These contain the rank order and expression
    #print len(expession_subset);sys.exit()
    expession_subset.sort() ### This order now matches that of 
    gene_expression_db=[]
    return expession_subset

def importGeneExpressionValues(filename,tissue_specific_db,translation_db,expession_subset):
    ### Import gene-level expression raw values           
    fn=filepath(filename); x=0; genes_added={}; gene_expression_db={}
    dataset_name = string.split(filename,delim)[-1][:-4]
    #print 'importing:',dataset_name
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        
        if x==0:
            if '#' not in data[0]:
                for i in t[1:]: sample_headers.append(i)
                x=1
        else:
            gene = t[0]
            #if '-' not in gene and ':E' in gene: print gene;sys.exit()
            if analysis_type == 'AltExon':
                try: ens_gene,exon = string.split(gene,'-')[:2]
                except Exception: exon = gene
                gene = exon
            if keyed_by == 'translation': ### alternative value is 'primaryID'
                """if gene == 'ENSMUSG00000025915-E19.3':
                    for i in translation_db: print [i], len(translation_db); break
                    print gene, [translation_db[gene]];sys.exit()"""
                try: gene = translation_db[gene] ### Ensembl annotations
                except Exception: gene = 'null'
            
            if gene in tissue_specific_db:
                index,tissue_exp=tissue_specific_db[gene]
                try: genes_added[gene]+=1
                except Exception: genes_added[gene]=1
                try: exp_vals = map(float, t[1:])
                except Exception:
                    ### If a non-numeric value in the list
                    exp_vals=[]
                    for i in t[1:]:
                        try: exp_vals.append(float(i))
                        except Exception: exp_vals.append(i)
                if value_type == 'calls': ### Hence, this is a DABG or RNA-Seq expression
                    exp_vals = produceDetectionCalls(exp_vals,targetPlatform) ### 0 or 1 calls
                gene_expression_db[gene] = [index,exp_vals]
    #print len(gene_expression_db), 'matching genes in the dataset and tissue compendium database'
    
    for gene in genes_added:
        if genes_added[gene]>1: del gene_expression_db[gene] ### delete entries that are present in the input set multiple times (not trustworthy)
        else: expession_subset.append(gene_expression_db[gene]) ### These contain the rank order and expression
    #print len(expession_subset);sys.exit()
    expession_subset.sort() ### This order now matches that of 
    gene_expression_db=[]
    return expession_subset, sample_headers

def produceDetectionCalls(values,Platform):
    # Platform can be the compendium platform (targetPlatform) or analyzed data platform (platform or array_type)
    new=[]
    for value in values:
        if Platform == 'RNASeq':
            if value>1:
                new.append(1) ### expressed
            else:
                new.append(0)
        else:
            if value<cutoff: new.append(1)
            else: new.append(0)
    return new

def importGeneIDTranslations(filename):
    ### Import ExpressionOutput/DATASET file to obtain Ensembl associations (typically for Affymetrix 3' arrays)
    fn=filepath(filename); x=0; translation_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            headers = t; x=1; index=0
            for i in headers:
                if 'Ensembl' in i: ens_index = index; break
                index+=1
        else:
            uid = t[0]
            ens_geneids = t[ens_index]
            ens_geneid = string.split(ens_geneids,'|')[0] ### In v.2.0.5, the first ID is the best protein coding candidate
            if len(ens_geneid)>0:
                translation_db[uid] = ens_geneid
    return translation_db

def remoteImportExonIDTranslations(array_type,species,translate_to_genearray,targetplatform):
    global targetPlatform; targetPlatform = targetplatform
    translation_db = importExonIDTranslations(array_type,species,translate_to_genearray)
    return translation_db
    
def importExonIDTranslations(array_type,species,translate_to_genearray):
    gene_translation_db={}; gene_translation_db2={}
    if targetPlatform == 'gene' and translate_to_genearray == 'no':
        ### Get gene array to exon array probeset associations
        gene_translation_db = importExonIDTranslations('gene',species,'yes')
        for geneid in gene_translation_db:
            exonid = gene_translation_db[geneid]
            gene_translation_db2[exonid] = geneid
            #print exonid, geneid
        translation_db = gene_translation_db2
    else:

        filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_'+array_type+'-exon_probesets.txt'
        ### Import exon array to target platform translations (built for DomainGraph visualization)
        fn=filepath(filename); x=0; translation_db={}
        print 'Importing the translation file',string.split(fn,delim)[-1][:-4]
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0:  x=1
            else:
                platform_id,exon_id = t
                if targetPlatform == 'gene' and translate_to_genearray == 'no':
                    try:
                        translation_db[platform_id] = gene_translation_db[exon_id] ### return RNA-Seq to gene array probeset ID
                        #print platform_id, exon_id, gene_translation_db[exon_id];sys.exit()
                    except Exception: null=[]
                else:
                    translation_db[platform_id] = exon_id
        del gene_translation_db; del gene_translation_db2
    return translation_db

def analyzeTissueSpecificExpressionPatterns(tissue_specific_db,expession_subset,sampleHeaders):
    tissue_specific_sorted = []; genes_present={}; tissue_exp_db={}; gene_order_db={}; gene_order=[]
    gene_list=[]
    for (index,vals) in expession_subset: genes_present[index]=[]
    for gene in tissue_specific_db:
        gene_list.append(gene)
        tissue_specific_sorted.append(tissue_specific_db[gene])
        gene_order_db[tissue_specific_db[gene][0]] = gene ### index order (this index was created before filtering)
    tissue_specific_sorted.sort()

    new_index=0
    for (index,tissue_exp) in tissue_specific_sorted:
        try:
            null=genes_present[index]
            i=0
            gene_order.append([new_index,gene_order_db[index]]); new_index+=1
            for f in tissue_exp:
                ### The order of the tissue specific expression profiles is based on the import gene order
                try: tissue_exp_db[tissues[i]].append(f)
                except Exception: tissue_exp_db[tissues[i]] = [f]
                i+=1
            
        except Exception: null=[] ### Gene is not present in the input dataset

    sample_exp_db={}; 
    for (index,exp_vals) in expession_subset:
        i=0
        for f in exp_vals:
            try: s = sampleHeaders[i]
            except Exception: print i, len(sampleHeaders); print index, len(exp_vals), f;sys.exit()
            ### The order of the tissue specific expression profiles is based on the import gene order
            try: sample_exp_db[s].append(f)
            except Exception: sample_exp_db[s] = [f]
            i+=1

    if correlate_by_order == 'yes':
        ### Rather than correlate to the absolute expression order, correlate to the order of expression (lowest to highest)
        sample_exp_db = replaceExpressionWithOrder(sample_exp_db)
        tissue_exp_db = replaceExpressionWithOrder(tissue_exp_db)

    global tissue_comparison_scores; tissue_comparison_scores={}
    
    if correlate_to_tissue_specific == 'yes':
        ### Create a gene_index that reflects the current position of each gene
        gene_index={}
        for (i,gene) in gene_order: gene_index[gene] = i
        ### Create a tissue to gene-index from the gene_index 
        tissue_to_index={}
        for tissue in tissue_to_gene:
            for gene in tissue_to_gene[tissue]:
                if gene in gene_index: ### Some are not in both tissue and sample datasets
                    index = gene_index[gene] ### Store by index, since the tissue and expression lists are sorted by index
                    try: tissue_to_index[tissue].append(index)
                    except Exception: tissue_to_index[tissue] = [index]
            tissue_to_index[tissue].sort()
        sample_exp_db,tissue_exp_db = returnTissueSpecificExpressionProfiles(sample_exp_db,tissue_exp_db,tissue_to_index)
    
    distributionNull = True
    if Permute:
        import copy
        sample_exp_db_original = copy.deepcopy(sample_exp_db)
        tissue_exp_db_original = copy.deepcopy(tissue_exp_db)
        group_list=[]; group_db={}
        for sample in sample_exp_db:
            group = string.split(sample,':')[0]
            try: group_db[group].append(sample)
            except Exception: group_db[group] = [sample]
        
        if distributionNull:
            group_lengths=[]
            for group in group_db:
                group_lengths.append(len(group_db[group]))
            group_db={}
            for sample in sample_exp_db:
                group = 'null1'
                try: group_db[group].append(sample)
                except Exception: group_db[group] = [sample]
            group_db['null2'] = group_db['null1']
            
            choice = random.sample
            tissue_groups = ['null1','null2']
        else:
            choice = random.choice
            tissue_groups = tuple(tissues)

        permute_groups=[]
        groups=[]
        gn=0
        for group in group_db:
            samples = group_db[group]
            permute_db={}; x=0
            while x<200:
                if distributionNull:
                    size = group_lengths[gn]
                    psamples = choice(samples,size)
                else: psamples = [choice(samples) for _ in xrange(len(samples))] ### works for random.sample or choice (with replacement)
                permute_db[tuple(psamples)]=None
                x+=1
            permute_groups.append(permute_db)
            groups.append(group); gn+=1 ### for group sizes
        
        groups.sort()
        permute_group1 = permute_groups[0]
        permute_group2 = permute_groups[1]
        
        permute_group1_list=[]
        permute_group2_list=[]
        for psamples in permute_group1:
            permute_group1_list.append(psamples)
        for psamples in permute_group2:
            permute_group2_list.append(psamples)

        i=0; diff_list=[]
        group_zdiff_means={}
        sample_diff_zscores=[]
        for psamples1 in permute_group1_list:
            psamples2 = permute_group2_list[i] #this is the second group to compare to
            x=0; permute_sample_exp_db={}
            for sample in psamples1:
                if distributionNull:
                    nsample = 'null1:'+string.split(sample,':')[1] ### reassign group ID
                    new_sampleID=nsample+str(x)
                else: new_sampleID=sample+str(x)
                try: permute_sample_exp_db[new_sampleID]=sample_exp_db[sample]
                except Exception: print 'Error:', sample, new_sampleID, sample_exp_db[sample];sys.exit()
                x+=1
            for sample in psamples2:
                if distributionNull:
                    nsample = 'null2:'+string.split(sample,':')[1] ### reassign group ID
                    new_sampleID=nsample+str(x)
                else: new_sampleID=sample+str(x)
                permute_sample_exp_db[new_sampleID]=sample_exp_db[sample]
                x+=1
            i+=1
            
            new_tissue_exp_db={}
            ### Create a new reference from the permuted data
            for sample in permute_sample_exp_db:
                group = string.split(sample,':')[0]
                try: new_tissue_exp_db[group].append(permute_sample_exp_db[sample])
                except Exception: new_tissue_exp_db[group] = [permute_sample_exp_db[sample]]
            for group in new_tissue_exp_db:
                k = new_tissue_exp_db[group]
                new_tissue_exp_db[group] = [Average(value) for value in zip(*k)] ### create new reference from all same group sample values
            
            PearsonCorrelationAnalysis(permute_sample_exp_db,new_tissue_exp_db)
            
            zscore_output_dir,tissue_scores = exportCorrelationResults()
            tissue_comparison_scores={}
            
            headers = list(tissue_scores['headers']); del tissue_scores['headers']
            index=0; positive=0; positive_score_diff=0
            sample_number = (len(headers)-1)
            diff_z_list=[]
            population1_denom=0; population1_pos=0; population2_pos=0; population2_denom=0
                        
            group_diff_z_scores={} ### Keep track of the differences between the z-scores between the two groups
            while index < sample_number:
                j=0
                #ref1 = tissue_groups[0]+':'; ref2 = tissue_groups[-1]+':'
                sample = headers[index+1]
                diff_z = tissue_scores[tissue_groups[0]][index]-tissue_scores[tissue_groups[-1]][index]
                diff_list.append([diff_z,sample])
                group = string.split(sample,':')[0]
                try: group_diff_z_scores[group].append(diff_z)
                except Exception: group_diff_z_scores[group] = [diff_z]
                sample_diff_zscores.append(diff_z)
                index+=1
                
            for group in group_diff_z_scores:
                avg_group_zdiff = Average(group_diff_z_scores[group])
                try: group_zdiff_means[group].append(avg_group_zdiff)
                except Exception: group_zdiff_means[group] = [avg_group_zdiff]
                
        diff_list.sort()
        
        all_group_zdiffs=[]
        for group in group_zdiff_means:
            all_group_zdiffs += group_zdiff_means[group]
        all_group_zdiffs.sort()

        #print sample_diff_zscores;sys.exit()

        #for i in diff_list: print i
        #sys.exit()
        
        i=1
        groups.reverse()
        group1,group2 = groups[:2]
        group1+=':'; group2+=':'
        scores=[]
        #print max(diff_list), min(diff_list);sys.exit()
        while i < len(diff_list):
            g1_hits=0; g2_hits=0
            list1 = diff_list[:i]
            list2 = diff_list[i:]
            for (z,s) in list1:
                if group1 in s: g1_hits+=1
            for (z,s) in list2:
                if group2 in s: g2_hits+=1
            sensitivity = float(g1_hits)/len(list1)
            specificity = float(g2_hits)/len(list2)
            accuracy = sensitivity+specificity
            #accuracy = g1_hits+g2_hits
            
            #print g1_hits, len(list1)
            #print g2_hits, len(list2)
            #print sensitivity, specificity;sys.exit()
            z_cutoff = Average([list1[-1][0],list2[0][0]])
            scores.append([accuracy,z_cutoff])
            i+=1
        
        scores.sort(); scores.reverse()
        print scores[0][0],'\t',scores[0][1]
        
        sample_exp_db = sample_exp_db_original
        tissue_exp_db = tissue_exp_db_original
        
    PearsonCorrelationAnalysis(sample_exp_db,tissue_exp_db)
    sample_exp_db=[]; tissue_exp_db=[]
    zscore_output_dir,tissue_scores = exportCorrelationResults()    
    return zscore_output_dir, tissue_scores, sampleHeaders

def returnTissueSpecificExpressionProfiles(sample_exp_db,tissue_exp_db,tissue_to_index):
    tissue_exp_db_abreviated={}
    sample_exp_db_abreviated={} ### This db is designed differently than the non-tissue specific (keyed by known tissues)

    ### Build the tissue specific expression profiles    
    for tissue in tissue_exp_db:
        tissue_exp_db_abreviated[tissue] = []
        for index in tissue_to_index[tissue]:
            tissue_exp_db_abreviated[tissue].append(tissue_exp_db[tissue][index]) ### populate with just marker expression profiles

    ### Build the sample specific expression profiles
    for sample in sample_exp_db:
        sample_tissue_exp_db={}
        sample_exp_db[sample]
        for tissue in tissue_to_index:
            sample_tissue_exp_db[tissue] = []
            for index in tissue_to_index[tissue]:
                sample_tissue_exp_db[tissue].append(sample_exp_db[sample][index])
        sample_exp_db_abreviated[sample] = sample_tissue_exp_db
    return sample_exp_db_abreviated, tissue_exp_db_abreviated

def replaceExpressionWithOrder(sample_exp_db):
    for sample in sample_exp_db:
        sample_exp_sorted=[]; i=0
        for exp_val in sample_exp_db[sample]: sample_exp_sorted.append([exp_val,i]); i+=1
        sample_exp_sorted.sort(); sample_exp_resort = []; order = 0
        for (exp_val,i) in sample_exp_sorted: sample_exp_resort.append([i,order]); order+=1
        sample_exp_resort.sort(); sample_exp_sorted=[] ### Order lowest expression to highest
        for (i,o) in sample_exp_resort: sample_exp_sorted.append(o) ### The expression order replaces the expression, in the original order
        sample_exp_db[sample] = sample_exp_sorted ### Replace exp with order
    return sample_exp_db

def PearsonCorrelationAnalysis(sample_exp_db,tissue_exp_db):
    #print "Beginning LineageProfiler analysis"
    k=0
    original_increment = int(len(tissue_exp_db)/15.00); increment = original_increment
    p = 1 ### Default value if not calculated

    for tissue in tissue_exp_db:
        #print k,"of",len(tissue_exp_db),"classifier tissue/cell-types"
        if k == increment: increment+=original_increment; #print '*',
        k+=1
        tissue_expression_list = tissue_exp_db[tissue]
        for sample in sample_exp_db:
            if correlate_to_tissue_specific == 'yes':
                ### Keyed by tissue specific sample profiles
                sample_expression_list = sample_exp_db[sample][tissue] ### dictionary as the value for sample_exp_db[sample]
                #print tissue, sample_expression_list
                #print tissue_expression_list; sys.exit()
            else: sample_expression_list = sample_exp_db[sample]
            try:
                ### p-value is likely useful to report (not supreemly accurate but likely sufficient)
                if len(tissue_expression_list) != len(sample_expression_list):
                    print len(tissue_expression_list), len(sample_expression_list)
                    print "Program Error!!! The length of the input expression list does not match the reference expression list (could indicate duplicate sample names)"
                    pass
                rho,p = stats.pearsonr(tissue_expression_list,sample_expression_list)
                #print rho,p
                #rho,p = stats.spearmanr(tissue_expression_list,sample_expression_list)
                try: pearson_list[sample].append(rho)
                except Exception: pearson_list[sample] = [rho]
                try: tissue_comparison_scores[tissue].append([rho,p,sample])
                except Exception: tissue_comparison_scores[tissue] = [[rho,p,sample]]
            except Exception:
                ### simple pure python implementation - no scipy required (not as fast though and no p-value)
                try:
                    rho = pearson(tissue_expression_list,sample_expression_list); p=0
                    try: pearson_list[sample].append(rho)
                    except Exception: pearson_list[sample] = [rho]
                    try: tissue_comparison_scores[tissue].append([rho,p,sample])
                    except Exception: tissue_comparison_scores[tissue] = [[rho,p,sample]]
                    pearson_list.append(rho)
                except Exception: None ### Occurs when an invalid string is present - ignore and move onto the next model
            """
            from stats_scripts import salstat_stats
            tst = salstat_stats.TwoSampleTests(tissue_expression_list,sample_expression_list)
            pp,pr = tst.PearsonsCorrelation()
            sp,sr = tst.SpearmansCorrelation()
            print tissue, sample
            if rho>.5: print [rho, pr, sr],[pp,sp];sys.exit()
            if rho<.5: print [rho, pr, sr],[pp,sp];sys.exit()
            """
    sample_exp_db=[]; tissue_exp_db=[]
    #print 'Correlation analysis finished'
    
def pearson(array1,array2):
    item = 0; sum_a = 0; sum_b = 0; sum_c = 0    
    while item < len(array1):
        a = (array1[item] - Average(array1))*(array2[item] - Average(array2))
        b = math.pow((array1[item] - Average(array1)),2)
        c = math.pow((array2[item] - Average(array2)),2)        
        sum_a = sum_a + a
        sum_b = sum_b + b
        sum_c = sum_c + c
        item = item + 1
    r = sum_a/math.sqrt(sum_b*sum_c)
    return r


def Median(array):
    array.sort()
    len_float = float(len(array))
    len_int = int(len(array))
    if (len_float/2) == (len_int/2):
        try: median_val = avg([array[(len_int/2)-1],array[(len_int/2)]])
        except IndexError: median_val = ''
    else:
        try: median_val = array[len_int/2]
        except IndexError: median_val = ''
    return median_val

def Average(array):
    try: return sum(array)/len(array)
    except Exception: return 0

def adjustPValues():
    """ Can be applied to calculate an FDR p-value on the p-value reported by scipy.
        Currently this method is not employed since the p-values are not sufficiently
        stringent or appropriate for this type of analysis """
        
    from stats_scripts import statistics
    all_sample_data={}
    for tissue in tissue_comparison_scores:
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            all_sample_data[sample] = db = {} ### populate this dictionary and create sub-dictionaries
        break
    
    for tissue in tissue_comparison_scores:
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            gs = statistics.GroupStats('','',p)
            all_sample_data[sample][tissue] = gs 
    for sample in all_sample_data:
        statistics.adjustPermuteStats(all_sample_data[sample])
     
    for tissue in tissue_comparison_scores:
        scores = []
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            p = all_sample_data[sample][tissue].AdjP()
            scores.append([r,p,sample])
        tissue_comparison_scores[tissue] = scores

def stdev(array):
    sum_dev = 0
    try: x_bar = scipy.average(array)
    except Exception: x_bar=Average(array)
    n = float(len(array))
    for x in array:
        x = float(x)
        sq_deviation = math.pow((x-x_bar),2)
        sum_dev += sq_deviation

    try:
        s_sqr = (1.0/(n-1.0))*sum_dev #s squared is the variance
        s = math.sqrt(s_sqr)
    except Exception:
        s = 'null'
    return s

def replacePearsonPvalueWithZscore():
    adjust_rho=False
    all_sample_data={}
    for tissue in tissue_comparison_scores:
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            all_sample_data[sample] = [] ### populate this dictionary and create sub-dictionaries
        break

    for tissue in tissue_comparison_scores:
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            if adjust_rho:
                try: r = 0.5*math.log(((1+r)/(1-r)))
                except Exception: print 'Error1:',tissue, sample, r, p; sys.exit()    
            all_sample_data[sample].append(r)
            #print tissue, sample, r

    #sample_stats={}
    all_dataset_rho_values=[]
    ### Get average and standard deviation for all sample rho's
    for sample in all_sample_data:
        try:
            all_sample_data[sample].sort() ### Sort, since when collapsing references, only the top two matter
            all_dataset_rho_values+=all_sample_data[sample][-2:]
            #try: avg=scipy.average(all_sample_data[sample])
            #except Exception: avg=Average(all_sample_data[sample])
        except Exception:
            all_dataset_rho_values+=all_sample_data[sample]
            #try: avg=scipy.average(all_sample_data[sample])
            #except Exception: avg=Average(all_sample_data[sample])
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides import warnings
            st_dev=stdev(all_sample_data[sample])
        #sample_stats[sample]=avg,st_dev
        
    try: global_rho_avg = scipy.average(all_dataset_rho_values)
    except Exception: global_rho_avg=Average(all_sample_data[sample])
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides import warnings
        global_rho_stdev = stdev(all_dataset_rho_values)
        
    ### Replace the p-value for each rho
    for tissue in tissue_comparison_scores:
        scores = []
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            pearson_rho = r
            if adjust_rho:
                try: r = 0.5*math.log(((1+r)/(1-r)))
                except Exception: print tissue, sample, r, p; sys.exit()
            #u,s=sample_stats[sample]
            #z = (r-u)/s
            z = (r-global_rho_avg)/global_rho_stdev ### Instead of doing this for the sample background, do it relative to all analyzed samples
            #z_alt = (r-global_rho_avg)/global_rho_stdev
            scores.append([pearson_rho, r,z,sample])
            #print sample, r, global_rho_avg, global_rho_stdev, z
        tissue_comparison_scores[tissue] = scores

def exportCorrelationResults():
    corr_output_file = string.replace(exp_output_file,'DATASET','LineageCorrelations')
    corr_output_file = string.replace(corr_output_file,'.txt','-'+coding_type+'.txt')
    if analysis_type == 'AltExon':
        corr_output_file = string.replace(corr_output_file,coding_type,'AltExon')
    filename =  string.split(corr_output_file,delim)[-1][:-4]
    #score_data = exportFile(corr_output_file)

    zscore_output_dir = string.replace(corr_output_file,'.txt','-zscores.txt')
    #probability_data = exportFile(zscore_output_dir)
    #adjustPValues()
    #sample_pearson_db = copy.deepcopy(tissue_comparison_scores) ### prior to pearson normalization
    replacePearsonPvalueWithZscore()
    
    ### Make title row
    headers=['Sample_name']
    for tissue in tissue_comparison_scores:
        for (pearson_rho, r,z,sample) in tissue_comparison_scores[tissue]: headers.append(sample)
        break
    title_row = string.join(headers,'\t')+'\n'
    #score_data.write(title_row)
    #if use_scipy: probability_data.write(title_row)
    ### Export correlation data
    tissue_scores = {}; tissue_probabilities={}; tissue_score_list = [] ### store and rank tissues according to max(score)
    tissue_scores2={}
    for tissue in tissue_comparison_scores:
        correlations=[]
        probabilities=[]
        for (pearson_rho, r,z,sample) in tissue_comparison_scores[tissue]:
            correlations.append(pearson_rho) ### un-adjusted correlation
            probabilities.append((z,pearson_rho))
        tissue_score_list.append((max(correlations),tissue))
        tissue_scores[tissue] = probabilities ### These are actually z-scores
        tissue_scores2[tissue] = string.join(map(str,[tissue]+correlations),'\t')+'\n' ### export line
        if use_scipy:
            tissue_probabilities[tissue] = string.join(map(str,[tissue]+probabilities),'\t')+'\n'
        
    tissue_score_list.sort()
    tissue_score_list.reverse()
    #for (score,tissue) in tissue_score_list: score_data.write(tissue_scores2[tissue])
        #if use_scipy: probability_data.write(tissue_probabilities[tissue])
    #score_data.close()
    #if use_scipy: probability_data.close()
    #print filename,'exported...'
    tissue_scores['headers'] = headers
    return zscore_output_dir, tissue_scores

def visualizeLineageZscores(zscore_output_dir,grouped_lineage_zscore_dir,graphic_links):
    from visualization_scripts import clustering
    ### Perform hierarchical clustering on the LineageProfiler Zscores
    graphic_links = clustering.runHCOnly(zscore_output_dir,graphic_links)   
    return graphic_links

def crossValidation(filename,setsToOutput=10,outputName=None):
    """ Function for performing 2-level cross-validation. This entails randomly dividing the input set into 2/3rds
    of the samples from each indicated biological group represented (equal prortion of disease and control each time)
    calculating a mean reference for each gene in the 2/3rds set and then exporting the remaining 1/3rd of samples """
    
    print 'Importing data to permute for 2-level cross-validation'
    fn = filepath(filename)
    dct_db={}; target_db={}; sample_list=[]; genes=[]
    
    ### Get the known group types (only works on two groups right now)
    firstLine = True
    group_list=[]
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            headers = t[1:]
            firstLine = False
            for h in headers:
                if ':' in h:
                    group = string.split(h,':')[0]
                    if group not in group_list:
                        group_list.append(group)
                  
    try:
        import collections
        group_samples=collections.OrderedDict()
    except Exception: import ordereddict as collections
    group_samples=collections.OrderedDict()
    
    firstLine = True
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            firstLine = False
            headers = t[1:]
            group1=[]; group2=[]
            for h in headers:
                group = string.split(h,':')[0]
                try: group_samples[group].append(headers.index(h))
                except Exception: group_samples[group] = [headers.index(h)]
        else:
            genes.append(t[0])
            for sample in headers:
                i = headers.index(sample)
                try:
                    dct_db[sample].append(t[i+1])
                except Exception:
                    dct_db[sample] = [t[i+1]]
            
    permute_set = 1
    
    #for i in dct_db: print i, dct_db[i]
    #sys.exit()
    
    inputTwoThirdsFiles=[]; inputOneThirdFiles=[]; referenceFiles=[]
    while permute_set < (setsToOutput+1):
        if setsToOutput == 1:
            output_file = string.replace(filename,'.txt','_training.txt')
        else:
            output_file = string.replace(filename,'.txt','-'+str(permute_set)+'_2-3rds'+'.txt')
        inputTwoThirdsFiles.append(output_file)
        export_obj = export.ExportFile(output_file)
        
        group_samples_two_thirds={}
        if setsToOutput == 1:
            ref_file = string.replace(filename,'.txt','-reference.txt')
        else:
            ref_file = string.replace(filename,'.txt','-'+str(permute_set)+'-reference'+'.txt')
        referenceFiles.append(ref_file)
        ref_obj = export.ExportFile(ref_file)
        for group in group_list:
            group_names = group_samples[group]
            two_thirds = int((2.0/3.0)*len(group_names)) ### e.g., AR
            group_random = random.sample(group_names,two_thirds)
            group_samples_two_thirds[group] = group_random
            
        permute_db={}; group_permute_db={}
        
        for group in group_list:
            for i in group_samples_two_thirds[group]:
                s = headers[i]
                values = dct_db[s]
                permute_db[s] = values
                try:
                    db = group_permute_db[group]
                    db[s] = values
                except Exception:
                    db = {}
                    db[s] = values
                    group_permute_db[group] = db
                    
        ### Get the 1/3rd remain samples for export to a separate file
        if setsToOutput == 1:
            outputR_file = string.replace(filename,'.txt','_test.txt')
        else:
            outputR_file = string.replace(filename,'.txt','-'+str(permute_set)+'_1-3rds'+'.txt')
        inputOneThirdFiles.append(outputR_file)

        exportR_obj = export.ExportFile(outputR_file)
        remaining_headers=[]; remaining_sample_db={}
        for s in dct_db:
            if s not in permute_db:
                remaining_sample_db[s] = dct_db[s]
                remaining_headers.append(s)
        exportR_obj.write(string.join(['UID']+remaining_headers,'\t')+'\n')

        for gene in genes:
            i = genes.index(gene)
            values = [gene]
            for s in remaining_headers:
                values.append(remaining_sample_db[s][i])
            values = string.join(values,'\t')+'\n'
            exportR_obj.write(values)
        exportR_obj.close()
        
        ### Export samples and references fro the 2/3rds set
        updated_headers = []
        for s in permute_db:
            updated_headers.append(s)
        export_obj.write(string.join(['UID']+updated_headers,'\t')+'\n')
        
        group_headers={}
        for group in group_list:
            for s in group_permute_db[group]:
                try: group_headers[group].append(s)
                except Exception: group_headers[group] = [s]
            
        for gene in genes:
            i = genes.index(gene)
            values = [gene]
            for s in updated_headers:
                values.append(permute_db[s][i])
            values = string.join(values,'\t')+'\n'
            export_obj.write(values)
        export_obj.close()
        
        ref_obj.write(string.join(['UID']+group_list,'\t')+'\n')
        for gene in genes:
            i = genes.index(gene)
            group_avgs=[]
            for group in group_list:
                group_values = []
                gdb = group_permute_db[group]
                for s in gdb:
                    try: group_values.append(float(gdb[s][i]))
                    except Exception: pass ### Exclude columns with NA from mean calculation
                group_avg = str(Average(group_values))
                group_avgs.append(group_avg)
            values = string.join([gene]+group_avgs,'\t')+'\n'
            ref_obj.write(values)
        ref_obj.close()
        permute_set+=1
    
    return inputTwoThirdsFiles, inputOneThirdFiles, referenceFiles
                 
def crossValidationAnalysis(species,platform,exp_input,exp_output,codingtype,compendium_platform,
                        modelSize,geneModels,permute, useMulti, finalNumberSetsToOutput):

    inputTwoThirdsFiles, inputOneThirdFiles, referenceFiles = crossValidation(exp_input,setsToOutput=1,outputName=None) ### This is the set to validate
    #inputTwoThirdsFiles = '/Users/saljh8/Desktop/Sarwal-New/UCRM_bx-transposed-train8.txt'
    #inputOneThirdFiles = '/Users/saljh8/Desktop/Sarwal-New/UCRM_bx-transposed-test_most.txt'
    #referenceFiles = ['/Users/saljh8/Desktop/Sarwal-New/manual-ref.txt']
    inputTwoThirdsFile = inputTwoThirdsFiles[0]; inputOneThirdFile = inputOneThirdFiles[0]; reference = referenceFiles[0]
    
    exp_input = string.replace(exp_input,'.txt','_training.txt')

    inputTwoThirdsFiles, inputOneThirdFiles, referenceFiles = crossValidation(exp_input,setsToOutput=finalNumberSetsToOutput) ### This is the set to validate
    #referenceFiles = ['/Users/saljh8/Desktop/Sarwal-New/manual-ref.txt']*len(referenceFiles)
    index = 0
    comparison_db={} ### store all of the 1/3rd validation results
    all_models={}
    
    ### Iterate through each 2/3rds set -> discover the best models over 80% -> evalute in 1/3rd set -> store those above 80%
    for exp_in in inputTwoThirdsFiles:
        exp_out = exp_in
        customMarkers = referenceFiles[index]
        top_hit_db = runLineageProfiler(species,platform,exp_in,exp_out,codingtype,compendium_platform,modelSize='optimize',
                        customMarkers=customMarkers,geneModels=geneModels,permute=permute,useMulti=useMulti)
        model_file = exportModels(exp_in,top_hit_db)

        exp_in = string.replace(exp_in,'_2-3rds','_1-3rds')
        
        try:
            top_hit_db = runLineageProfiler(species,platform,exp_in,exp_out,codingtype,compendium_platform,
                        customMarkers=customMarkers,geneModels=model_file,permute=permute)
        except Exception: top_hit_db={}
        #comparison_db[exp_in] = top_hit_db
        all_models.update(top_hit_db)
        index+=1

    ### Create a combined models file
    exp_in = string.replace(exp_input,'.txt','-combinedModels.txt')
    model_file = exportModels(exp_in,all_models)
    index = 0
    ### Re-analyze all of the 2/3rd and 1/3rd files with these
    for exp_in in inputTwoThirdsFiles:
        customMarkers = referenceFiles[index]
        top_hit_db = runLineageProfiler(species,platform,exp_in,exp_out,codingtype,compendium_platform,modelSize=modelSize,
                        customMarkers=customMarkers,geneModels=model_file,permute=permute)
        comparison_db[exp_in] = top_hit_db
        exp_in = string.replace(exp_in,'_2-3rds','_1-3rds')
        top_hit_db = runLineageProfiler(species,platform,exp_in,exp_out,codingtype,compendium_platform,modelSize=modelSize,
                        customMarkers=customMarkers,geneModels=model_file,permute=permute)
        comparison_db[exp_in] = top_hit_db
        index+=1
        
    score_db={}
    for one_third_input_name in comparison_db:
        top_hit_db = comparison_db[one_third_input_name]
        for model in top_hit_db:
            score = top_hit_db[model]
            try: score_db[model].append(float(score))
            except Exception: score_db[model] = [float(score)]
            
    score_ranked_list=[]
    for model in score_db:
        avg_score = Average(score_db[model])
        min_score = min(score_db[model])
        score = avg_score*min_score
        score_ranked_list.append([score,avg_score,min_score,model])
    score_ranked_list.sort()
    score_ranked_list.reverse()
    for s in score_ranked_list:
        print s
    
    
def exportModels(exp_in,top_hit_db):
    model_file = string.replace(exp_in,'.txt','_topModels.txt')
    export_obj = export.ExportFile(model_file)
    for m in top_hit_db:
        if float(top_hit_db[m])>70:
            m = list(m)
            m = string.join(m,',')+'\n'
            export_obj.write(m)
    export_obj.close()
    return model_file

def modelScores(results_dir):
    score_db={}
    files = unique.read_directory(results_dir+'/')
    for file in files:
        firstLine = True
        if '-ModelScore' in file:
            filename = results_dir+'/'+file
            fn=filepath(filename)
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                if firstLine == True: firstLine = False
                else:
                    t = string.split(data,'\t')
                    model=t[-2]
                    score = t[0]
                    sensitivity = float(t[1])/float(t[2])*100
                    specificity = float(t[3])/float(t[4])*100
                    score = sensitivity+specificity
                    try: score_db[model].append(float(score))
                    except Exception: score_db[model] = [float(score)]
    score_ranked_list=[]
    for model in score_db:
        avg_score = Average(score_db[model])
        min_score = min(score_db[model])
        stdev_min = stdev(score_db[model])
        l = (avg_score*min_score)
        score_ranked_list.append([l,avg_score,min_score,stdev_min,model,score_db[model]])
    score_ranked_list.sort()
    score_ranked_list.reverse()
    x = 100; y=1
    
    output_dir = results_dir+'/overview.txt'
    oe = export.ExportFile(output_dir)
    for s in score_ranked_list:
        y+=1
        s = map(lambda x: str(x), s)
        #if y==x: sys.exit()
        #print s
        oe.write(string.join(s,'\t')+'\n')
    oe.close()
    #sys.exit()

def allPairwiseSampleCorrelation(fn):
    import numpy
    firstLine = True
    group_list=[]
    control_set={}
    exp_set={}
    all_samples_names=[]
    all_sample_values=[]
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            firstLine=False
        else:
            values = map(float,t[1:])
            group = string.split(t[0],':')[0]
            if 'no' in group:
                #control_set[t[0]] = []
                s = 0 ### control
            else:
                #exp_set[t[0]] = []
                s = 1 ### experiment
            all_samples_names.append((s,t[0]))
            all_sample_values.append(values)
    all_sample_values = numpy.array(all_sample_values)

    ### calculate all pairwise comparisons and store exp and control correlations in two separate dictionaries
    D1 = numpy.ma.corrcoef(all_sample_values)
    i=0
    for score_ls in D1:
        scores = []
        s,sample = all_samples_names[i]
        if s==0: control_set[sample] = score_ls
        if s==1: exp_set[sample] = score_ls
        i+=1
        k=0
        #for rho in score_ls: sample2 = all_samples_names[k]; k+=1
        
    def subtract(values):
        return values[0]-values[1]
    
    from stats_scripts import statistics
    results=[]
    n = len(all_samples_names)
    e = len(exp_set)
    c = len(control_set)
    for cs in control_set:
        ccor = control_set[cs] ### all sample correlations for the control
        for es in exp_set:
            ecor = exp_set[es] ### all sample correlations for the exp
            diff_corr = [subtract(value) for value in zip(*[ecor,ccor])]
            k=0
            score=0; sensitivity=0; specificity=0
            for diff in diff_corr:
                s,sample = all_samples_names[k]; k+=1
                if s==1 and diff>0: score+=1; sensitivity+=1
                elif s==0 and diff<0: score+=1; specificity+=1
            sens = float(sensitivity)/e
            spec = float(specificity)/c
            accuracy = float(score)/n
            avg_accuracy = statistics.avg([sens,spec])
            results.append([avg_accuracy,accuracy,sens,spec,es,cs])
    results.sort()
    results.reverse()
    for i in results[:20]:
        #if i[1]>.70: print i
        print i
        
def harmonizeClassifiedSamples(species,reference_exp_file, query_exp_file, classification_file,fl=None):
    """
    The goal of this function is to take LineageProfilerIterative classified samples to a reference matrix,
    combine the reference matrix and the query matrix at the gene symbol level, retain the original reference
    column and row orders, then re-order the query samples within the context of the reference samples the
    correlations were derived from. In doing so, the harmonization occurs by attempting to map the expression
    values within a common numerical format to minimize dataset specific effects (folds, non-log expression).
    The function culminates in a heatmap showing the genes and re-ordered samples in thier listed orders.
    Outlier samples with low correlations will ultimately need to be represented outside of the reference
    sample continuum. The function attempts to propagate sample group labels from the reference and query sets
    (if avaialble, in the sample headers or in a groups file), and indicate original gene and column clusters
    if available. We call this approach cellHarmony.
    """
    """
    Alternative description: To represent new cells/samples within the continuum of established scRNA-Seq profiles
    we have developed a robust and fast data harmonization function. This data harmonization approach applies a
    k-nearest neighbor classification approach to order submitted cells along the continuum of reference gene expression
     profiles, without alternating the order of the reference cells and cluster groups in original reference matrix.
    The established group labels are in effect propagated to new queried samples. The final data is represented in an
    intuitive queriable heatmap and dimensionality reduction scatter plot (t-SNE, PCA). To harmonize the input exp. data,
    the input file formats are emperically determined (log, non-log, log-fold) and standardized to log2 expression.
    Submitted FASTQ files to AltAnalyze will be optionally pseudoaligned and normalized to maximize compatibility with the
    existing reference datasets. Alternative classification algoriths and options for submission of new reference sets will
    be supported. Issues may arrise in the expression harmonization step, however, profiles are already biased towards genes
    with cell population restricted expression, largely abrogating this problem.

    """
    ### Set the user or default options for cellHarmony Differential Expression Analyses
    try: pearsonThreshold=fl.PearsonThreshold()
    except: pearsonThreshold = 0.1

    try: peformDiffExpAnalysis=fl.PeformDiffExpAnalysis()
    except: peformDiffExpAnalysis = True
    
    try: use_adjusted_pval=fl.UseAdjPvalue()
    except: use_adjusted_pval = False
    
    try: pvalThreshold=float(fl.PvalThreshold())
    except: pvalThreshold = 0.05
    
    try: FoldCutoff = fl.FoldCutoff()
    except: FoldCutoff = 1.5
            
    ### Output the alignment results and perform the differential expression analysis
    output_file,query_output_file,folds_file,DEGs_combined = importAndCombineExpressionFiles(species,reference_exp_file,
            query_exp_file,classification_file,pearsonThreshold=pearsonThreshold,peformDiffExpAnalysis=peformDiffExpAnalysis,
            pvalThreshold=pvalThreshold,fold_cutoff=FoldCutoff,use_adjusted_pval=use_adjusted_pval)
    
    ### Output the cellHarmony heatmaps
    from visualization_scripts import clustering
    row_method = None; row_metric = 'cosine'; column_method = None; column_metric = 'euclidean'; color_gradient = 'yellow_black_blue'
    transpose = False; Normalize='median'
    
    if runningCommandLine:
        display = False
    else:
        display = True
        display = False
        
    graphics = clustering.runHCexplicit(query_output_file, [], row_method, row_metric, column_method,
            column_metric, color_gradient, transpose, Normalize=Normalize, contrast=5, display=display)

    graphics = clustering.runHCexplicit(output_file, [], row_method, row_metric, column_method,
            column_metric, color_gradient, transpose, Normalize=Normalize, contrast=5, display=display)
    
    zscore = True
    graphics=[]
    transpose='no'
    
    ### Build-UMAP plot
    import UI
    try:
        UI.performPCA(output_file, 'no', 'UMAP', False, None, plotType='2D',
            display=False, geneSetName=None, species=species, zscore=True, reimportModelScores=False,
            separateGenePlots=False, returnImageLoc=True)
    except:
        try:
            print 'UMAP error encountered (dependency not met), trying t-SNE'
            UI.performPCA(output_file, 'no', 't-SNE', False, None, plotType='2D',
                display=False, geneSetName=None, species=species, zscore=True, reimportModelScores=False,
                separateGenePlots=False, returnImageLoc=True)
        except: pass
    
    try: fl.setSpecies(species); fl.setVendor("3'array")
    except:
        import UI
        fl = UI.ExpressionFileLocationData(folds_file,'','','')
        fl.setSpecies(species); fl.setVendor("3'array")
    fl.setOutputDir(export.findParentDir(folds_file))
        
    ### Run MarkerFinder
    if len(DEGs_combined)>0:
        fl.setRPKMThreshold(0.00)
        fl.setCorrelationDirection('up')
        compendiumType = 'protein_coding'
        genesToReport = 100000
        correlateAll = True
        import markerFinder
        markerFinder.analyzeData(folds_file,species,platform,compendiumType,geneToReport=genesToReport,correlateAll=correlateAll,AdditionalParameters=fl,logTransform=False)
        export.deleteFolder(fl.OutputDir()+'/MarkerFinder-positive')
        os.rename(fl.OutputDir()+'/MarkerFinder',fl.OutputDir()+'/MarkerFinder-positive')
        fl.setCorrelationDirection('down')
        markerFinder.analyzeData(folds_file,species,platform,compendiumType,geneToReport=genesToReport,correlateAll=correlateAll,AdditionalParameters=fl,logTransform=False)
        export.deleteFolder(fl.OutputDir()+'/MarkerFinder-negative')
        os.rename(fl.OutputDir()+'/MarkerFinder',fl.OutputDir()+'/MarkerFinder-negative')
        exportCombinedMarkerFinderResults(folds_file,fl.OutputDir()+'/MarkerFinder-positive',fl.OutputDir()+'/MarkerFinder-negative',DEGs_combined,species)

def exportCombinedMarkerFinderResults(folds_file,positive_dir,negative_dir,DEGs_combined,species):
    positive_dir+= '/AllGenes_correlations-ReplicateBased.txt'
    negative_dir+= '/AllGenes_correlations-ReplicateBased.txt'
    import collections
    marker_db=collections.OrderedDict()
    
    def importMarkerFinderResults(input_file,direction):
        header_row=True
        for line in open(input_file,'rU').xreadlines():
            data = cleanUpLine(line)
            data = string.replace(data,'"','')
            geneID,	symbol,	rho, p, cell_state = string.split(data,'\t')
            if header_row:
                header_row=False
            else:
                rho = float(rho)
                if geneID in marker_db:
                    ### Positive marker genes
                    alt_cell_state,alt_direction,alt_rho = marker_db[geneID]
                    if rho>alt_rho:
                        marker_db[geneID]=cell_state, direction, rho
                elif geneID in DEGs_combined:
                    marker_db[geneID]=cell_state, direction, rho
    importMarkerFinderResults(positive_dir,'positive')
    importMarkerFinderResults(negative_dir,'negative')

    ordered_markers={}
    for gene in marker_db:
        cell_state, direction, rho = marker_db[gene]
        try: ordered_markers[cell_state].append([direction, rho, gene])
        except Exception: ordered_markers[cell_state]= [[direction, rho, gene]]
    for cell_state in ordered_markers:
        ordered_markers[cell_state].sort()
        ordered_markers[cell_state].reverse()
        
    headers,matrix = getDataMatrix(folds_file)
    export_file = folds_file[:-4]+'-MarkerFinder.txt'
    export_object = export.ExportFile(export_file)
    headers2=[]
    for h in headers:
        headers2.append(h+':'+h)
    export_object.write(string.join(['UID']+headers2,'\t')+'\n')
    for h in headers:
        if h in ordered_markers:
            for (direction, rho, gene) in ordered_markers[h]:
                export_object.write(string.join([h+':'+gene]+map(str,matrix[gene]),'\t')+'\n')
    export_object.close()
    
    from visualization_scripts import clustering
    row_method = None; row_metric = 'cosine'; column_method = None; column_metric = 'euclidean'; color_gradient = 'yellow_black_blue'
    transpose = False; Normalize=False;  display = False
    import UI
    gsp = UI.GeneSelectionParameters(species,"RNASeq","RNASeq")
    gsp.setPathwaySelect('None Selected')
    gsp.setGeneSelection('')
    gsp.setOntologyID('')
    gsp.setGeneSet('None Selected')
    gsp.setJustShowTheseIDs('')                
    gsp.setTranspose(False)
    gsp.setNormalize(False)
    gsp.setGeneSelection('')
    gsp.setClusterGOElite('GeneOntology')
    
    graphics = clustering.runHCexplicit(export_file, [], row_method, row_metric,
                column_method, column_metric, color_gradient, gsp, Normalize=Normalize,
                contrast=3, display=display)

class ClassificationData:
    def __init__(self,sample,score,assigned_class):
        self.sample = sample; self.score = score; self.assigned_class = assigned_class
    def Sample(self): return self.sample
    def Score(self): return self.score
    def AssignedClass(self): return self.assigned_class
      
def importAndCombineExpressionFiles(species,reference_exp_file,query_exp_file,classification_file,platform='RNASeq',pearsonThreshold=0.1,
            peformDiffExpAnalysis=True, pvalThreshold=0.05,fold_cutoff=1.5, use_adjusted_pval=False):
    """Harmonize the numerical types and feature IDs of the input files, then combine """
     
    import collections
    from visualization_scripts import clustering
    original_reference_exp_file = reference_exp_file
    reference_exp_file_alt = string.replace(reference_exp_file,'-centroid.txt','.txt')
    status = verifyFile(reference_exp_file_alt)
    if status==True:
        reference_exp_file = reference_exp_file_alt
        includeAllReferences=True
    else:
        includeAllReferences=False
    ref_filename = export.findFilename(reference_exp_file)
    query_filename = export.findFilename(query_exp_file)
    root_dir = string.replace(export.findParentDir(query_exp_file),'ExpressionInput','')+'/cellHarmony/'
    output_dir = ref_filename[:-4]+'__'+query_filename ### combine the filenames and query path
    output_dir = root_dir+'/exp.'+string.replace(output_dir,'exp.','')
    output_dir =string.replace(output_dir,'-OutliersRemoved','')
    groups_dir = string.replace(output_dir,'exp.','groups.')
    ref_exp_db,ref_headers,ref_col_clusters,cluster_format_reference = importExpressionFile(reference_exp_file)
    cluster_results = clustering.remoteImportData(query_exp_file,geneFilter=ref_exp_db)
    if len(cluster_results[0])>0: filterIDs = ref_exp_db
    else: filterIDs = False
    query_exp_db,query_headers,query_col_clusters,cluster_format_query = importExpressionFile(query_exp_file,ignoreClusters=True,filterIDs=False)
    ref_filename_clean = string.replace(ref_filename,'exp.','')
    ref_filename_clean = string.replace(ref_filename_clean,'.txt','')
    ref_filename_clean = 'cellHarmony-Reference'
    query_filename_clean = string.replace(query_filename,'exp.','')
    query_filename_clean = string.replace(query_filename_clean,'.txt','')
    
    ### Define Reference groups and import prior reference file
    groups_file = groups_dir[:-4]+'-AllCells.txt'
    groups_file = string.replace(groups_file,'OutliersRemoved-','')
    if includeAllReferences:
        groups_file = string.replace(groups_file,'-centroid__','__')
    comps_file  = string.replace(groups_file,'groups.','comps.')
    groups_reference_file = string.replace(groups_file,'-AllCells.txt','-Reference.txt')
    try:
        import ExpressionBuilder ### If using centroids for classification, merge with original reference groups
        sample_group_ref_db = ExpressionBuilder.simplerGroupImport(groups_reference_file)
        groups_to_reference_cells = {}
        if '-centroid' in original_reference_exp_file:
            for cell in sample_group_ref_db:
                group = sample_group_ref_db[cell]
                try: groups_to_reference_cells[group].append(cell)
                except: groups_to_reference_cells[group]=[cell]
    except:
        sample_group_ref_db={}
        groups_to_reference_cells={}

    column_clusters=[]
    final_clusters=[]
    numeric_cluster=True
    comps=[]
    
    ### If using all cells as a reference, rather than centroid, then len(ref_col_clusters)==len(sample_group_ref_db)
    if len(sample_group_ref_db)>len(ref_col_clusters):
        for sample in sample_group_ref_db:
            column_clusters.append(sample_group_ref_db[sample]) ### this is an ordered dictionary
            try: final_clusters.append([ref_filename_clean,int(float(sample_group_ref_db[sample])),sample])  ### used later for creating the apples-to-apples group file
            except Exception:
                final_clusters.append([ref_filename_clean,sample_group_ref_db[sample],sample])
                numeric_cluster=False
    elif len(ref_col_clusters)>0:
        for sample in ref_col_clusters:
            column_clusters.append(ref_col_clusters[sample]) ### this is an ordered dictionary
            try: final_clusters.append([ref_filename_clean,int(float(ref_col_clusters[sample])),sample])  ### used later for creating the apples-to-apples group file
            except Exception:
                final_clusters.append([ref_filename_clean,ref_col_clusters[sample],sample])
                numeric_cluster=False
        
    """ Store an alternative reference for each header """
    ### In case a ":" is in one header but not in another, create an alternative reference set
    alt_ref_headers = []
    for sample in ref_headers:
        if ':' in sample:
            alt_sample_id = string.split(sample,':')[1]
            alt_ref_headers.append(alt_sample_id)
            if sample in ref_col_clusters: ### If clusters assigned in original file
                cluster = ref_col_clusters[sample]
                ref_col_clusters[alt_sample_id]=cluster  ### append to the original
    
    """ Store alternative query sample names with groups added """
    original_sampleID_translation={}
    for sample in query_headers:
        if ':' in sample:
            original_sampleID = string.split(sample,':')[1]
            original_sampleID_translation[original_sampleID] = sample
            
    """ Import the Classification data """
    input_dir = output_dir
    sample_classes={}
    new_headers=[]
    query_header_proppegated_clusters={}
    firstLine = True
    exclude=[]
    for line in open(classification_file,'rU').xreadlines():
        data = line.rstrip()
        data = string.replace(data,'"','')
        values = string.split(data,'\t')
        if firstLine:
            firstLine=False
            header_row = values
            score_index = header_row.index('Max-Rho')
            class_index = header_row.index('Predicted Class')
        else:
            sample = values[0]
            score = float(values[score_index])
            
            assigned_class = values[class_index]
            if sample in original_sampleID_translation:
                sample = original_sampleID_translation[sample]
            if len(groups_to_reference_cells)>0: ### Hence, centroids - replace centroids with individual reference cells
                assigned_class = groups_to_reference_cells[assigned_class][-1]
            cd = ClassificationData(sample,score,assigned_class)
            #print sample, assigned_class;sys.exit()
            try: sample_classes[assigned_class].append([score,cd])
            except Exception: sample_classes[assigned_class] = [[score,cd]]
            query_header_proppegated_clusters[sample]=assigned_class
            if score<pearsonThreshold: ### Minimum allowed correlation threshold
                exclude.append(sample)
    
    print len(exclude), 'cells excluded due to correlation below the indicated threshold'
    """ Assign a cluster label to the query sample if applicable """
    query_clusters=[]
    classified_samples={}
    if numeric_cluster: prefix = 'c'
    else: prefix = ''

    selected_query_headers=[]
    for sample in query_headers:
        original_sample_id=sample
        if ':' in sample:
            sample_alt = string.split(sample,':')[1]
        if sample in query_header_proppegated_clusters or sample_alt in query_header_proppegated_clusters:
            try: ref_sample = query_header_proppegated_clusters[sample]
            except: ref_sample = query_header_proppegated_clusters[sample_alt]
            if ref_sample in ref_col_clusters:
                if ':' in sample:
                    sample = string.split(sample,':')[1]
                cluster = ref_col_clusters[ref_sample]
                #column_clusters.append(cluster) ### Assign a predicted cluster label to the query sample
                try: final_clusters.append([query_filename_clean,int(float(cluster)),sample]) ### used later for creating the apples-to-apples group file
                except Exception: final_clusters.append([query_filename_clean,cluster,sample])
                try: ranked_cluster = int(cluster)
                except Exception: ranked_cluster = cluster
                comps.append([ranked_cluster,prefix+cluster+'_'+query_filename_clean,prefix+cluster+'_'+ref_filename_clean])
                classified_samples[sample]=cluster
                selected_query_headers.append(original_sample_id)
                
    for assigned_class in sample_classes:
        ### Ordered by score
        sample_classes[assigned_class].sort()
        sample_classes[assigned_class].reverse()
        
    """ Integrate the queried samples in with the reference continuum """
    ### Loop through the reference headers and add those in the query where they arise
    new_query_headers = []
    for sample in ref_headers:
        if sample not in new_headers:
            new_headers.append(sample)
            if sample in sample_classes:
                sample_data = sample_classes[sample]
                for (score,cd) in sample_data:
                    if cd.Sample() not in new_headers and cd.Sample() not in alt_ref_headers: ### In case some of the same reference samples were added to the query
                        if cd.Sample() not in new_headers:
                            new_headers.append(cd.Sample())
                            new_query_headers.append(cd.Sample())

    for sample in query_headers:
        if ':' in sample:
            sample_alt = string.split(sample,':')[1]
        try: cluster = classified_samples[sample]
        except: cluster = classified_samples[sample_alt]
        column_clusters.append(cluster)
        
    """ Combine the two datasets, before re-ordering """
    ### In case the UIDs in the two datasets are in different formats (we assume either Ensembl or Symbol)
    ### we allow conversion between these systems
    for refUID in ref_exp_db: break
    for queryUID in query_exp_db: break
    
    if ('ENS' in refUID or 'ENS' in queryUID) and '.' not in refUID:
        ### Covert one or both to gene symbol
        try:
            import gene_associations
            gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
        except Exception: gene_to_symbol={}
        if 'ENS' in refUID:
            ref_exp_db = convertFromEnsemblToSymbol(ref_exp_db,gene_to_symbol)
        if 'ENS' in queryUID:
            query_exp_db = convertFromEnsemblToSymbol(query_exp_db,gene_to_symbol)
    
    ### Write out the combined data
    export_object = export.ExportFile(output_dir)
    if len(column_clusters)>0:
        ### Output with gene and column clusters
        export_object.write(string.join(['UID','row_clusters-flat']+ref_headers+query_headers,'\t')+'\n') #
        export_object.write(string.join(['column_clusters-flat','']+column_clusters,'\t')+'\n')
        new_headers = ['row_clusters-flat']+new_headers
        new_query_headers = ['row_clusters-flat']+new_query_headers
    else:
        export_object.write(string.join(['UID']+ref_headers+query_headers,'\t')+'\n')
    
    for uid in ref_exp_db:
        if uid in query_exp_db:
            if cluster_format_reference: ### If the reference is an AltAnalyze heatmap format text file
                export_object.write(string.join([uid]+ref_exp_db[uid]+query_exp_db[uid],'\t')+'\n')
            else:
                export_object.write(string.join([uid,'1']+ref_exp_db[uid]+query_exp_db[uid],'\t')+'\n')
    export_object.close()
    
    """ Write out a groups file that lists samples per file for optional visualization """
    export_object = export.ExportFile(groups_dir[:-4]+'-ReOrdered.txt')
    for sample in ref_headers:
        if ':' in sample:
            sample = string.split(sample,':')[1]
        export_object.write(sample+'\t1\t'+ref_filename_clean+'\n')
    for sample in query_headers:
        if ':' in sample:
            sample = string.split(sample,':')[1]
        export_object.write(sample+'\t2\t'+query_filename_clean+'\n')
    export_object.close()
    
    """ Re-order the samples based on the classification analysis """
    ### The ref_headers has the original reference sample order used to guide the query samples
    from import_scripts import sampleIndexSelection
    input_file=output_dir
    output_file = input_file[:-4]+'-ReOrdered.txt'
    query_output_file = input_file[:-4]+'-ReOrdered-Query.txt'
    filter_names = new_headers
    
    """ Remove the correlation outliers """
    group_export_object = export.ExportFile(root_dir+'/groups.cellHarmony.txt')
    filter_names2=[]
    new_query_headers2=[]
    added=[]
    order=0
    ordering=[]
    for sample_name in filter_names:
        if sample_name not in added:
            added.append(sample_name)
            if sample_name not in exclude:
                filter_names2.append(sample_name)
                if sample_name in new_query_headers:
                    new_query_headers2.append(sample_name)
                    if ':' in sample_name:
                        sample_name = string.split(sample_name,':')[1]
                    try:
                        cluster = classified_samples[sample_name]
                        group_export_object.write(sample_name+'\t'+cluster+'\t'+prefix+cluster+'\n')
                    except Exception: pass
    group_export_object.close()
    
    sampleIndexSelection.filterFile(input_file,output_file,filter_names2)
    sampleIndexSelection.filterFile(input_file,query_output_file,new_query_headers2)
    
    """Export a combined groups and comps file to perform apples-to-apples comparisons"""
    folds_file=''
    all_DEGs=[]
    if len(final_clusters)>0 and peformDiffExpAnalysis:
        try:
            final_clusters.sort()
            new_query_headers=[]
            expression_file = string.replace(groups_file,'groups.','exp.')
            export_object = export.ExportFile(groups_file)
            group_counter=0
            added_groups={}
            group_numbers={}
            
            try: ### If using centroids for classification, merge with original reference groups
                import ExpressionBuilder 
                sample_group_ref_db = ExpressionBuilder.simplerGroupImport(groups_reference_file)
            except:
                pass
            
            ### The below tuple order controls which clusters and groups are listed in which orders
            groups_to_samples={}
            cluster_lookup={}
            for (source,cluster,sample) in final_clusters:
                cluster_name = prefix+str(cluster)+'_'+source
                cluster_lookup[cluster_name]=prefix+str(cluster)
                if cluster_name in added_groups:
                    added_groups[cluster_name]+=1
                    group_number = group_counter
                else:
                    added_groups[cluster_name]=1
                    group_counter += 1
                    group_number = group_counter
                    group_numbers[cluster_name]=group_number
                if ':' in sample:
                    sample = string.split(sample,':')[1]
                export_object.write(sample+'\t'+str(group_number)+'\t'+cluster_name+'\n')
                new_query_headers.append(sample)
                ### For summary plot export of cell cluster frequencies
                try: groups_to_samples[cluster_name].append(sample)
                except: groups_to_samples[cluster_name] = [sample]
            export_object.close()
            expression_file_reordered = expression_file[:-4]+'-reordered.txt'
            status = verifyFile(expression_file)
            if status==False: ### If not AllCells expression file, use the merged file (assuming all genes are in both the query and reference)
                export_object = export.ExportFile(expression_file)
                export_object.write(string.join(['UID']+ref_headers+query_headers,'\t')+'\n') #
                for uid in ref_exp_db:
                    if uid in query_exp_db:
                        export_object.write(string.join([uid]+ref_exp_db[uid]+query_exp_db[uid],'\t')+'\n')
                export_object.close()
                
            sampleIndexSelection.filterFile(expression_file,expression_file_reordered,new_query_headers)
            shutil.move(expression_file_reordered,expression_file) ### replace our temporary file
            
            export_object = export.ExportFile(comps_file)
            freq_object = export.ExportFile(root_dir+'/cell-frequency-stats.txt')
            comps = unique.unique(comps)
            comps.sort()
            
            ### Export the comps and frequency results
            freq_object.write('Cluster-name\t'+ref_filename_clean+'\t'+query_filename_clean+'\tFisher exact test p-value\t'+'% '+ref_filename_clean+'\t'+'% '+query_filename_clean+'\n')
            added=[]
            g2_all = len(ref_headers)
            g1_all = len(query_headers)
            for (null,group1,group2) in comps:
                g = cluster_lookup[group1]
                g1_len = len(groups_to_samples[group1])
                g2_len = len(groups_to_samples[group2])
                g1_frq = (g1_len*100.00)/g1_all
                g2_frq = (g2_len*100.00)/g2_all
                oddsratio, pvalue = stats.fisher_exact([[g2_len, g1_len], [g2_all-g2_len, g1_all-g1_len]])
                #print g2_len, g1_len, g2_all-g2_len, g1_all-g1_len, g2_all, g1_all, pvalue
                added.append(group1); added.append(group2)
                freq_object.write(g+'\t'+str(g2_len)+'\t'+str(g1_len)+'\t'+str(pvalue)+'\t'+str(g2_frq)[:4]+'\t'+str(g1_frq)[:4]+'\n')
                if added_groups[group1]>3 and added_groups[group2]>3:
                    ### Atleast 4 cells present in the two compared groups
                    g1 = str(group_numbers[group1])
                    g2 = str(group_numbers[group2])
                    export_object.write(g1+'\t'+g2+'\n')
            for group in groups_to_samples:
                g = cluster_lookup[group]
                if group not in added:
                    g2_len = len(groups_to_samples[group])
                    g2_frq = (g2_len*100.00)/g2_all
                    freq_object.write(g+'\t'+str(g2_len)+'\t'+str(0)+'\t'+str(pvalue)+'\t'+str(g2_frq)[:4]+'\t'+str(0)+'\n')
            export_object.close()
            freq_object.close()
            try:
                index1=1;index2=2; x_axis='Number of cells'; y_axis = 'Reference clusters'; title='Assigned Cell Frequencies'
                clustering.barchart(root_dir+'/cell-frequency-stats.txt',index1,index2,x_axis,y_axis,title)
            except Exception:
                pass
    
            from stats_scripts import metaDataAnalysis
            """
            strictCutoff = True
            pvalThreshold=0.05
            use_adjusted_pval = False
            if platform == 'RNASeq':
                log_fold_cutoff=0.585
                output_dir = root_dir+'/DEGs-LogFold_0.585_rawp'
                if strictCutoff:
                    use_adjusted_pval = True
                    log_fold_cutoff = 1
                    pvalThreshold = 0.05
                    output_dir = root_dir+'/DEGs-LogFold_1_adjp'
                
            else:
                log_fold_cutoff=0.1
                output_dir = root_dir+'Events-LogFold_0.1_rawp'
            """
            if use_adjusted_pval:
                pval_type = '_adjp_'+str(pvalThreshold)
            else:
                pval_type = '_rawp_'+str(pvalThreshold)
            use_custom_output_dir = 'DifferentialExpression_Fold_'+str(fold_cutoff)[:4]+pval_type
            output_dir = root_dir+use_custom_output_dir
            
            log_fold_cutoff = math.log(float(fold_cutoff),2)
            
            metaDataAnalysis.remoteAnalysis(species,expression_file,groups_file,platform=platform,use_custom_output_dir=use_custom_output_dir,
                                log_fold_cutoff=log_fold_cutoff,use_adjusted_pval=use_adjusted_pval,pvalThreshold=pvalThreshold)
            gene_summary_file = output_dir+'/gene_summary.txt'
            gene_summary_file2 = string.replace(gene_summary_file,use_custom_output_dir,'')
            export.copyFile(gene_summary_file,gene_summary_file2)
            index1=2;index2=3; x_axis='Number of DEGs'; y_axis = 'Reference clusters'; title='cellHarmony Differentially Expressed Genes'
            clustering.barchart(gene_summary_file2,index1,index2,x_axis,y_axis,title,color1='IndianRed',color2='SkyBlue')
            
            all_DEGs = aggregateRegulatedGenes(output_dir)
            display_genes = string.join(list(all_DEGs),' ')
            ICGS_DEGs_combined = ref_exp_db.keys()
            for gene in all_DEGs:
                if gene not in ICGS_DEGs_combined:
                    ICGS_DEGs_combined.append(gene) ### Add these genes at the end
            ICGS_DEGs_combined.reverse()
            all_DEGs2 = string.join(ICGS_DEGs_combined,' ')
            import UI
            vendor = 'Ensembl'
            gsp = UI.GeneSelectionParameters(species,platform,vendor)
            gsp.setGeneSet('None Selected')
            gsp.setPathwaySelect('')
            gsp.setGeneSelection(all_DEGs2)
            gsp.setJustShowTheseIDs(display_genes)
            gsp.setNormalize('median')
            transpose = gsp
            column_method = None #'hopach'
            column_metric = 'cosine'
            row_method = None #'hopach'
            row_metric = 'correlation'
            graphic_links=[]
            color_gradient = 'yellow_black_blue'
            from visualization_scripts import clustering
            if runningCommandLine:
                display = False
            else:
                display = True
                display = False
                
            outputDEGheatmap=False
            if outputDEGheatmap:
                graphic_links = clustering.runHCexplicit(expression_file, graphic_links, row_method, row_metric, column_method, column_metric, color_gradient, transpose, display=display, Normalize=True)
        except Exception:
            print traceback.format_exc()
            print '!!!!! NO merged expression file availble for differential expression analyis (apples-to-apples).'
            
    try:  
        folds_file = expression_file[:-4]+'-folds.txt'
        comparisons = sampleIndexSelection.getComparisons(groups_file)
        filter_names,group_index_db = sampleIndexSelection.getFilters(groups_file,calculateCentroids=True)
        sampleIndexSelection.filterFile(expression_file,folds_file,(filter_names,group_index_db),force=False,calculateCentroids=True,comparisons=comparisons)
        
        export_object = export.ExportFile(string.replace(folds_file,'exp.','groups.'))
        header = getHeader(folds_file)
        i=0
        for h in header:
            i+=1
            export_object.write(h+'\t'+str(i)+'\t'+h+'\n')
        export_object.close()
        export_object = export.ExportFile(string.replace(folds_file,'exp.','comps.'))
        export_object.close()
    except:
        pass
    return output_file, query_output_file, folds_file, all_DEGs

def aggregateRegulatedGenes(folder):
    import UI
    geneIDs=[]
    files = UI.read_directory(folder)
    for file in files:
        comp_geneIDs=[]
        if '.txt' in file and ('PSI.' in file or 'GE.' in file):
            ls=[]
            fn = folder+'/'+file
            firstLine = True
            for line in open(fn,'rU').xreadlines():
                data = line.rstrip()
                t = string.split(data,'\t')
                if firstLine:
                    if 'Event-Direction' in t:
                        uid_index = t.index('Event-Direction')
                        fold_index = t.index('dPSI')
                    else:
                        uid_index = t.index('Symbol')
                        fold_index = t.index('LogFold')
                    firstLine= False
                    continue
                uid = t[uid_index]
                fold = t[fold_index]
                comp_geneIDs.append((fold,uid))
        comp_geneIDs.sort() ### Sort by fold
        comp_geneIDs.reverse()
        for (fold,uid) in comp_geneIDs:
            if uid not in geneIDs: ### Remove this to inlcude redundant genes
                if len(uid)>0:
                    geneIDs.append(uid)
    return geneIDs
                    
def convertFromEnsemblToSymbol(exp_db,gene_to_symbol):
    ### covert primary ID to symbol
    import collections
    exp_db_symbol=collections.OrderedDict()
    for UID in exp_db:
        if UID in gene_to_symbol:
            symbol = gene_to_symbol[UID][0]
            exp_db_symbol[symbol] = exp_db[UID]
    return exp_db_symbol

def checkForGroupsFile(filename,headers):
    new_headers = headers
    import collections
    groups_db=collections.OrderedDict()
    if ('exp.' in filename or 'filteredExp.' in filename):
        filename = string.replace(filename,'-steady-state.txt','.txt')
        try:
            import ExpressionBuilder
            sample_group_db = ExpressionBuilder.simplerGroupImport(filename)
            new_headers = []
            for v in headers:
                if v in sample_group_db:
                    groups_db[v]=sample_group_db[v]
                    v = sample_group_db[v]+':'+v
                new_headers.append(v)
        except Exception:
            #print traceback.format_exc()
            pass
    return headers, new_headers, groups_db

def importExpressionFile(input_file,ignoreClusters=False,filterIDs=False):
    """ Import the expression value and harmonize (cellHarmony) to log, non-fold values """
    import collections
    expression_db=collections.OrderedDict()
    column_cluster_index=collections.OrderedDict()
    row_cluster_index=collections.OrderedDict()
    ### Check the file format (log, non-log, fold) - assume log is log2
    import ExpressionBuilder
    expressionDataFormat,increment,convertNonLogToLog = ExpressionBuilder.checkExpressionFileFormat(input_file,reportNegatives=True,filterIDs=filterIDs)
    inputFormat = 'NumericMatrix'
    firstLine = True
    cluster_format_file = True
    if filterIDs!=False:
        try:
            ### Conserve memory by only importing reference matching IDs
            import gene_associations; from import_scripts import OBO_import
            gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
            symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
        except Exception:
            gene_to_symbol={}
            symbol_to_gene={}
    row_count=0
    for line in open(input_file,'rU').xreadlines():
        data = line.rstrip()
        data = string.replace(data,'"','')
        row_count+=1
        if '.csv' in input_file:
            values = string.split(data,',')
        else:
            values = string.split(data,'\t')
        if firstLine:
            firstLine=False
            if 'row_clusters-flat' in values[1]:
                numStart = 2 ### Start numeric import at column 3
                inputFormat = 'Clustering'
                header_row = values[numStart:]
            else:
                numStart = 1 ### Start numeric import at column 2
                header_row = values[numStart:]
                original_header, header_row, column_cluster_index_alt=checkForGroupsFile(input_file,header_row) ### in case there is a groups file
        else:
            ### store the column to cluster relationships (reference only - propegated to query)
            ### We can simply retain this column for the output
            if 'column_clusters-flat' in values:
                clusters = values[numStart:]
                i=0
                for header in header_row:
                    cluster = clusters[i]
                    column_cluster_index[header]=cluster 
                    i+=1
                
                ### Replace the cluster with the sample name if there is only one sample per cluster (centroid name)
                if len(unique.unique(clusters)) == len(column_cluster_index):
                    for header in column_cluster_index:
                        column_cluster_index[header]=header
                continue
            elif row_count==2:
                header_row=original_header
                column_cluster_index = column_cluster_index_alt
                cluster_format_file = False
            uid = values[0]
            if filterIDs !=False:
                if uid not in filterIDs:
                    if uid in gene_to_symbol:
                        alt_id = gene_to_symbol[uid][0]
                        if alt_id not in filterIDs:
                            continue ### Skip additional processing of this line
                    elif uid in symbol_to_gene:
                        alt_id = symbol_to_gene[uid][0]
                        if alt_id not in filterIDs:
                            continue ### Skip additional processing of this line
                    else:
                        continue ### Skip additional processing of this line
            if inputFormat == 'Clustering':
                ### store the row cluster ID
                row_cluster_index[uid]=values[1] 
            numericVals = map(float, values[numStart:])
            if increment<-1 and convertNonLogToLog == False:
                ### Indicates the data is really fold changes (if the increment is a small negative,
                ### it could be due to quantile normalization and the values are really raw expression values)
                row_min = min(numericVals)
                numericVals = map(lambda x: x-row_min, numericVals) ### move into positive log expression space
            elif increment<-1 and convertNonLogToLog: ### shouldn't be encountered since fold values should all be logged
                row_min = min(numericVals)
                numericVals = map(lambda x: math.log(x-row_min+1,2), numericVals) ### move into positive log expression space
            elif convertNonLogToLog:
                try: numericVals = map(lambda x: math.log(x+increment,2), numericVals) ### log2 and increment
                except Exception:
                    print 'increment',increment
                    print numericVals[0:10],numericVals[-10:]
                    print traceback.format_exc()
                    kill
            numericVals = map(str,numericVals) ### we are saving to a file
            if inputFormat == 'Clustering' and ignoreClusters==False:
                expression_db[uid] = [values[1]]+numericVals
            else:
                expression_db[uid] = numericVals
    print len(expression_db),'IDs imported'

    return expression_db, header_row, column_cluster_index, cluster_format_file
            
def createMetaICGSResults(ICGS_files,outputDir,CenterMethod='median',species='Hs',PearsonThreshold=0.9):
    """This function combines ICGS or MarkerFinder results from multiple outputs"""
    # import all ICGS or MarkerFinder variable genes
    gene_db = simpleICGSGeneImport(ICGS_files)

    # import raw expression values for each ICGS or MarkerFinder
    files_to_merge = []

    for heatmap_file in ICGS_files:
        ### returnCentroids = False to return all cells
        cellHarmonyReferenceFile = convertICGSClustersToExpression(heatmap_file,heatmap_file,returnCentroids=True,
                    CenterMethod=CenterMethod,geneOverride=gene_db,combineFullDatasets=False,species=species)
        files_to_merge.append(cellHarmonyReferenceFile)
    
    join_option = 'Intersection'
    ID_option = False ### Require unique-matches only when True
    output_merge_dir = outputDir
    
    # output a combined centroid file, combining all input ICGS or MarkerFinder results with raw expression values
    from import_scripts import mergeFiles
    outputfile = mergeFiles.joinFiles(files_to_merge, join_option, ID_option, output_merge_dir)

    # collapse similar medoids into centroids to remove redundant references
    query_output_file = collapseSimilarMedoids(outputfile,cutoff=PearsonThreshold)
    
    # re-cluster this merged file with HOPACH to produce the final combined medoid reference
    from visualization_scripts import clustering
    row_method = 'hopach'; row_metric = 'correlation'; column_method = 'hopach'; column_metric = 'cosine'; color_gradient = 'yellow_black_blue'
    transpose = False; Normalize=False
    graphics = clustering.runHCexplicit(query_output_file, [], row_method, row_metric,
                column_method, column_metric, color_gradient, transpose, Normalize=Normalize,
                contrast=3, display=False)
    
    revised_cellHarmony_reference = graphics[-1][-1][:-4]+'.txt'
    final_output_dir = outputDir+'/CellHarmonyReference/ICGS-merged-reference.txt'
    exportMergedReference(revised_cellHarmony_reference,final_output_dir)
    
    return final_output_dir
    
def exportMergedReference(input,output):
    row_count=0
    eo = export.ExportFile(output)
    eo.write(string.join([],'\t')+'')
    for line in open(input,'rU').xreadlines():
        row_count+=1
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if row_count==1:    
            clusters = len(t)-1
            eo.write(string.join(['UID','row_clusters-flat']+map(str,range(1,clusters)),'\t')+'\n')
        elif row_count==2:
            clusters = len(t)-1
            eo.write(string.join(['column_clusters-flat','']+map(str,range(1,clusters)),'\t')+'\n')
        else:
            eo.write(line)
    eo.close()
    
def collapseSimilarMedoids(outputfile,cutoff=0.9):
    print 'correlation cutoff =',cutoff
    from visualization_scripts import clustering
    import numpy
    from stats_scripts import statistics
    matrix, column_header, row_header, dataset_name, group_db, priorColumnClusters, priorRowClusters = clustering.remoteImportData(outputfile)
    # Filter for rows with median values
    filteredMatrix = []
    filtered_row_header = []
    i=0
    for values in matrix:
        if (max(values)-min(values))==0:
            pass
        else:
            filteredMatrix.append(values)
            filtered_row_header.append(row_header[i])
        i+=1
    matrix = numpy.array(filteredMatrix)
    row_header=filtered_row_header
    Tmatrix = zip(*matrix) ### transpose this back to normal
    D1 = numpy.corrcoef(Tmatrix)
    
    i=0
    combine={}
    combined_list=[]
    for score_ls in D1:
        k=0
        for v in score_ls:
            if v!=1:
                if v>cutoff: ### Pearson cutoff
                    c1=column_header[i]
                    c2=column_header[k]
                    clusters = [c1,c2]
                    clusters.sort()
                    a,b = clusters
                    if a in combine:
                        if b not in combine[a]:
                            if b not in combine:
                                combine[a].append(b)
                                combined_list.append(b)
                    else:
                        if b in combine:
                            if a not in combine[b]:
                                combine[b].append(a)
                                combined_list.append(a)
                        else:
                            present=False
                            for x in combine:
                                if b in combine[x]:
                                    present = True
                            if present==False:
                                try: combine[a].append(b)
                                except Exception: combine[a] = [a,b]
                                combined_list.append(a)
                                combined_list.append(b)
            k+=1
        i+=1
        
    # Add each cluster not in correlated set to obtain all final clusters
    for cluster in column_header:
        if cluster not in combined_list:
            combine[cluster]=[cluster]
    
    group_index_db={} ### Order is not important yet
    for cluster in combine:
        ci_list=[]
        for c in combine[cluster]:
            ci = column_header.index(c) ### get the index of the cluster name
            ci_list.append(ci)
        combined_cluster_name = string.join(combine[cluster],'|')
        group_index_db[combined_cluster_name] = ci_list

    root_dir = export.findParentDir(outputfile)
    collapsed_dir = root_dir+'/CellHarmonyReference/collapsedMedoids.txt'
    eo = export.ExportFile(collapsed_dir)
    
    eo.write(string.join(['UID']+map(str,group_index_db),'\t')+'\n')
    i=0
    for values in matrix:
        avg_matrix=[]
        for cluster in group_index_db:
            try: avg_matrix.append(str(statistics.avg(map(lambda x: matrix[i][x], group_index_db[cluster]))))
            except Exception: ### Only one value
                try: avg_matrix.append(str(map(lambda x: matrix[i][x], group_index_db[cluster])[0]))
                except Exception:
                    print matrix[i]
                    print group_index_db[cluster];sys.exit()
        eo.write(string.join([row_header[i]]+avg_matrix,'\t')+'\n')
        i+=1

    return collapsed_dir
    
def convertICGSClustersToExpression(heatmap_file,query_exp_file,returnCentroids=False,
                CenterMethod='median',geneOverride=None,combineFullDatasets=True,species='Hs'):
    """This function will import an ICGS row normalized heatmap and return raw
    expression values substituted for the values. """

    from visualization_scripts import clustering
    graphic_links=[]
    filename = export.findFilename(heatmap_file)
    ICGS_dir = export.findParentDir(heatmap_file)
    if 'DataPlots' in ICGS_dir: ### For MarkerFinder input
        ### Go one more level up
        ICGS_dir = export.findParentDir(ICGS_dir[:-1])
    root_dir = export.findParentDir(ICGS_dir[:-1])
    files = unique.read_directory(root_dir+'/ExpressionInput')

    exp_dir_prefix = string.split(string.replace(filename,'Clustering-',''),'-')[0]

    ### Look for the best expression file match
    specific_matches = []
    steady_state_files = []
    exp_files = []
    filteredExp_files = []
    for file in files:
        if 'exp.' in file and '~exp.' not in file:
            exp_files.append([os.path.getsize(root_dir+'/ExpressionInput/'+file),file])
            if 'steady-state' in file:
                steady_state_files.append([os.path.getsize(root_dir+'/ExpressionInput/'+file),file])
                expression_files_ss
            if exp_dir_prefix in file:
                specific_matches.append([os.path.getsize(root_dir+'/ExpressionInput/'+file),file])
        if 'filteredExp.' in file: 
            filteredExp_files.append([os.path.getsize(root_dir+'/ExpressionInput/'+file),file])

    specific_matches.sort()
    steady_state_files.sort()
    exp_files.sort()
    filteredExp_files.sort()
    
    if len(specific_matches)>0:
        expdir = root_dir+'/ExpressionInput/'+specific_matches[-1][1]
    elif len(steady_state_files)>0:
        expdir = root_dir+'/ExpressionInput/'+steady_state_files[-1][1]
    else:
        expdir = root_dir+'/ExpressionInput/'+exp_files[-1][1]
    try: filtered_expdir = root_dir+'/ExpressionInput/'+filteredExp_files[-1][1]
    except: filtered_expdir = ''

    print 'Selected the full expression file:',expdir
    if '-Guide' in filename:
        guide = string.split(filename,'-Guide')[1][:1]
        cellHarmonyReferenceFile = root_dir+'/CellHarmonyReference/Guide'+guide+'-cellHarmony-reference.txt'
    elif 'Marker' in filename:
        cellHarmonyReferenceFile = root_dir+'/CellHarmonyReference/MarkerFinder-cellHarmony-reference.txt'
    else:
        cellHarmonyReferenceFile = root_dir+'/CellHarmonyReference/ICGS-cellHarmony-reference.txt'
    eo = export.ExportFile(cellHarmonyReferenceFile)
    
    ### Import the heatmap and expression files
    matrix, column_header, row_header, dataset_name, group_db, priorColumnClusters, priorRowClusters = clustering.remoteImportData(heatmap_file,reverseOrder=False)
    gene_cluster_db={}
    if len(priorRowClusters)==0:
        clusters=[]
        cluster_number=0
        new_row_header = []
        for uid in row_header:
            if ':' in uid:
                cluster,uid = string.split(uid,':')
                if cluster not in clusters:
                    clusters.append(cluster)
                    cluster_number+=1
                if ' ' in uid:
                    uid = string.split(uid,' ')[0] ### Typically source ID and then symbol
                gene_cluster_db[uid] = str(cluster_number)
            new_row_header.append(uid)
        clusters=[]
        cluster_number=0
        new_column_header = []
        for uid in column_header:
            if ':' in uid:
                cluster,uid = string.split(uid,':')
                if cluster not in clusters:
                    clusters.append(cluster)
                    cluster_number+=1
                priorColumnClusters.append(str(cluster_number))
            new_column_header.append(uid)
        row_header = new_row_header
        column_header = new_column_header
    else:
        ### Store the row clusters in a dictionary rather than a list since matching IDs in the expression file may change (e.g., different ID systems in the ICGS results)
        index=0
        for gene in row_header:
            gene_cluster_db[gene] = priorRowClusters[index]
            index+=1

    ### Correct fileheader if group prefix present
    try: column_header = map(lambda x: string.split(x,':')[1],column_header)
    except Exception: pass ### If no ":" in cell library names
    
    updated_column_header=[]
    for i in column_header:
        if ':' in i:
            i = string.split(i,':')[1]
        updated_column_header.append(i)
    column_header = updated_column_header
    
    ### Record the index for each sample name in the ICGS result order in the original expression file (exp.*)
    priorColumnClusters = map(str,priorColumnClusters)
        
    if geneOverride != None:
        row_header = geneOverride.keys() ### Replace with a custom or expanded list of gene IDs to use

    ### Expand the possible set of usable IDs

    matrix_exp, column_header_exp, row_header_exp, dataset_name, group_db_exp = clustering.importData(expdir,geneFilter=row_header)
    percent_found = (len(row_header_exp)*1.00)/len(row_header)
    if percent_found<0.5:
        print "...Incompatible primary ID (Symbol), converting to Ensembl"
        import gene_associations; from import_scripts import OBO_import
        gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
        symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
        ### Assume the original IDs were Ensembl IDs and get IDs from the filtered Expression File
        filteredIDs = simpleICGSGeneImport([filtered_expdir])
        row_header2=[]
        index=0
        for symbol in row_header:
            if symbol in symbol_to_gene:
                for gene in symbol_to_gene[symbol]:
                    if 'LRG' not in gene and gene in filteredIDs:
                        row_header2.append(gene)
                    gene_cluster_db[gene]=priorRowClusters[index]
            index+=1
        row_header = row_header2
        matrix_exp, column_header_exp, row_header_exp, dataset_name, group_db_exp = clustering.importData(expdir,geneFilter=row_header)

    ### Correct fileheader if group prefix present
    try: column_header_exp = map(lambda x: string.split(x,':')[1],column_header_exp)
    except Exception: pass ### If no ":" in cell library names
    
    try: sample_index_list = map(lambda x: column_header_exp.index(x), column_header)
    except Exception: ### Hence, there are missing cells in the expression file
        sample_index_list=[]
        for cell in column_header:
            if cell in column_header_exp:
                sample_index_list.append(column_header_exp.index(cell))
        print "WARNING... only",len(sample_index_list), "out of",len(column_header), "cells found in the expression file."
        
    if geneOverride != None:
        eo.write(string.join(['UID']+column_header,'\t')+'\n')
    else:
        eo.write(string.join(['UID','row_clusters-flat']+column_header,'\t')+'\n')
        eo.write(string.join(['column_clusters-flat','']+priorColumnClusters,'\t')+'\n')
    index=0
    
    import collections
    reference_matrix = collections.OrderedDict() ### store the reordered data for later medioid calculation
    for uid in row_header:
        if uid in row_header_exp:
            ### Get the heatmap ordered gene entries
            exp_row_index = row_header_exp.index(uid)
            ### Re-order the samples according to the heatmap
            try: reordered_values = map(lambda x: str(matrix_exp[exp_row_index][x]), sample_index_list) ### simple and fast way to reorganize the samples
            except Exception:
                ### For PSI files with missing values at the end of each line, often
                if len(column_header_exp) != len(matrix_exp[exp_row_index]):
                    diff = len(column_header_exp)-len(matrix_exp[exp_row_index])
                    matrix_exp[exp_row_index]+=diff*['']
                reordered_values = map(lambda x: str(matrix_exp[exp_row_index][x]), sample_index_list)
            if geneOverride != None:
                eo.write(string.join([uid]+reordered_values,'\t')+'\n')
                reference_matrix[uid,0]=map(float,reordered_values)
            else:
                eo.write(string.join([uid,str(gene_cluster_db[uid])]+reordered_values,'\t')+'\n')
                reference_matrix[uid,str(gene_cluster_db[uid])]=map(float,reordered_values)
        index+=1
    eo.close()
    
    print 'New ICGS CellHarmony reference saved to:',cellHarmonyReferenceFile
    
    ### Create the groups for median calculation
    group_index_db = collections.OrderedDict()
    index=0
    for cluster in priorColumnClusters:
        try: group_index_db[cluster].append(index) #group_index_db['cluster-'+cluster].append(index)
        except Exception: group_index_db[cluster]=[index] #group_index_db['cluster-'+cluster]=[index]
        index+=1
    
    cellHarmonyReferenceFileMediod = cellHarmonyReferenceFile[:-4]+'-centroid.txt'
    eo = export.ExportFile(cellHarmonyReferenceFileMediod)
    
    if geneOverride != None:
        eo.write(string.join(['UID']+map(str,group_index_db),'\t')+'\n')
    else:
        eo.write(string.join(['UID','row_clusters-flat']+map(str,group_index_db),'\t')+'\n')
        eo.write(string.join(['column_clusters-flat','']+map(lambda x: string.replace(x,'cluster-',''),group_index_db),'\t')+'\n')
    from stats_scripts import statistics
    for uid in reference_matrix:
        median_matrix=[]
        for cluster in group_index_db:
            if CenterMethod == 'mean' or CenterMethod == 'centroid':
                try: median_matrix.append(str(statistics.avg(map(lambda x: reference_matrix[uid][x], group_index_db[cluster]))))
                except Exception: ### Only one value
                    median_matrix.append(str(map(lambda x: reference_matrix[uid][x], group_index_db[cluster])))
            else: ### Median
                try: median_matrix.append(str(statistics.median(map(lambda x: reference_matrix[uid][x], group_index_db[cluster]))))
                except Exception: ### Only one value
                    median_matrix.append(str(map(lambda x: reference_matrix[uid][x], group_index_db[cluster])))
        if geneOverride != None:
            eo.write(string.join([uid[0]]+median_matrix,'\t')+'\n')
        else:
            eo.write(string.join([uid[0],uid[1]]+median_matrix,'\t')+'\n')
    eo.close()
    
    if combineFullDatasets:
        """ Merge the query and the reference expression files """
        query_exp_matrix, query_header, query_row_header, null, group_db_exp = clustering.importData(query_exp_file)
        reference_exp_matrix, reference_header, ref_row_header, null, group_db_exp = clustering.importData(expdir)

        """ Simple combine of the two expression files without normalization """
        try:
            #ref_filename = export.findFilename(expdir)
            ref_filename = export.findFilename(cellHarmonyReferenceFile)
            query_filename = export.findFilename(query_exp_file)
            root_dir = export.findParentDir(query_exp_file)
            output_dir = ref_filename[:-4]+'__'+query_filename[:-4]+'-AllCells.txt' ### combine the filenames and query path
            output_dir = string.replace(output_dir,'OutliersRemoved-','')
            output_dir = string.replace(output_dir,'-centroid__','__')
            ### Matches the outputdir name from the function importAndCombineExpressionFiles
            output_dir = string.replace(root_dir,'/ExpressionInput','')+'/cellHarmony/exp.'+string.replace(output_dir,'exp.','')
            group_ref_export_dir = string.replace(output_dir,'exp.','groups.')
            group_ref_export_dir = string.replace(group_ref_export_dir,'-AllCells.','-Reference.')
        
            eo = export.ExportFile(output_dir)
            eo.write(string.join(['UID']+reference_header+query_header,'\t')+'\n')
            for uid in ref_row_header:
                ri=ref_row_header.index(uid)
                if uid in query_row_header:
                    qi=query_row_header.index(uid)
                    eo.write(string.join([uid]+map(str,reference_exp_matrix[ri])+map(str,query_exp_matrix[qi]),'\t')+'\n')
            eo.close()
            print 'Combined query and reference expression file exported to:'
            eo = export.ExportFile(group_ref_export_dir)
            
            index=0
            for cell in column_header:
                eo.write(cell+'\t'+str(priorColumnClusters[index])+'\t'+str(priorColumnClusters[index])+'\n')
                index+=1
            eo.close()
        except:
            print 'Failed to join the query and reference expression files due to:'
            print traceback.format_exc()
    if returnCentroids == True or returnCentroids == 'yes' or returnCentroids == 'centroid':
        print 'Using centroids rather than individual cells for alignment.'
        return cellHarmonyReferenceFileMediod
    else:
        print 'Using invidual cells rather than cell centroids for alignment.'
        return cellHarmonyReferenceFile

def simpleExpressionFileImport(filename):
    fn=filepath(filename)
    expression_db = {}
    firstRow = True
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            header = t[1:]
            firstRow = False
        else:
            expression_db[t[0]] = t[1:]
    return expression_db, header

def simpleICGSGeneImport(files):
    """ Import the gene IDs from different ICGS or MarkerFinder results prior
    to combining to derive combined ICGS results and making combined medoid file"""
    try:
        import collections
        gene_db=collections.OrderedDict()
    except Exception:
        import ordereddict as collections
        gene_db=collections.OrderedDict()

    for file in files:
        fn=filepath(file)
        firstRow = True
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if firstRow:
                header = t[1:]
                firstRow = False
            else:
                gene = t[0]
                if '-flat' not in gene:
                    if gene not in gene_db:
                        if ':' in gene: ### If MarkerFinder as the input
                            geneID = string.split(gene,':')[1]
                            gene = string.split(geneID,' ')[0]
                        gene_db[gene]={}
    return gene_db

def compareICGSpopulationFrequency(folder):
    ea = export.ExportFile(folder+'/overlap.tsv')
    files = unique.read_directory(folder)
    cluster_counts={}
    datasets = []
    clusters=[]
    from collections import Counter
    for file in files:
        if '.txt' in file:
            datasets.append(file[:-4])
            fn = folder+'/'+file
            LineCount = 0
            for line in open(fn,'rU').xreadlines():
                data = line.rstrip()
                t = string.split(data,'\t')
                LineCount += 1
                if LineCount<3:
                    if LineCount==2:
                        """ Examine the second row only for cluster counts """
                        dataset_cluster_counts = Counter(t[2:])
                        for cluster in dataset_cluster_counts:
                            try: cluster_counts[cluster].append(str(dataset_cluster_counts[cluster]))
                            except Exception: cluster_counts[cluster] = [str(dataset_cluster_counts[cluster])]
                            if cluster not in clusters:
                                clusters.append(cluster)
                else:
                    break
    clusters.sort()
    print len(cluster_counts)
    print len(datasets)
    ea.write(string.join(['Cluster']+datasets,'\t')+'\n')
    for cluster in clusters:
        ea.write(string.join([cluster]+cluster_counts[cluster],'\t')+'\n')
    ea.close()
                
if __name__ == '__main__':
    icgs_dir = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/ICGS/Clustering-exp.NaturePanorma-Ly6G-Guide3-Augment-F2r-hierarchical_cosine_correlation.txt'
    exp_dir = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/Ly6g/ExpressionInput/exp.Guide3-cellHarmony-revised.txt'
    #convertICGSClustersToExpression(icgs_dir,exp_dir);sys.exit()
    """compareICGSpopulationFrequency('/Users/saljh8/Desktop/dataAnalysis/Collaborative/Jose/NewTranscriptome/cellHarmonyResults/');sys.exit()
    """
    reference_exp_file = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/cellHarmony-evaluation/Grimes/exp.WT.txt'
    query_exp_file = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/cellHarmony-evaluation/Grimes/exp.AML.txt'
    classification_file= '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/cellHarmony-evaluation/Grimes/CellClassification/AML-CellClassification.txt'
    harmonizeClassifiedSamples('Hs',reference_exp_file,query_exp_file,classification_file);sys.exit()
    
    #modelScores('/Users/saljh8/Desktop/dataAnalysis/LineageProfiler/Training/CellClassification');sys.exit()
    #allPairwiseSampleCorrelation('/Users/saljh8/Desktop/Sarwal-New/UCRM_bx.txt');sys.exit()

    try:
        import multiprocessing as mlp
        mlp.freeze_support()
    except Exception:
        mpl = None

    ################  Default Variables ################
    species = 'Hs'
    platform = "exon"
    vendor = 'Affymetrix'
    compendium_platform = "exon"
    codingtype = 'protein_coding'
    platform = vendor, platform
    exp_output = None
    geneModels = False
    modelSize = None
    permute = False
    useMulti = False
    finalNumberSetsToOutput = 10
    cross_validation = False

    """ This script iterates the LineageProfiler algorithm (correlation based classification method) to identify sample types relative
    two one of two references given one or more gene models. The program '
    """
    #python LineageProfilerIterate.py --i "/Users/nsalomonis/Desktop/dataAnalysis/qPCR/ExpressionInput/exp.ABI_Pediatric.txt" --r "/Users/nsalomonis/Desktop/dataAnalysis/qPCR/ExpressionOutput/MarkerFinder/MarkerFinder-ABI_Pediatric.txt" --m "/Users/nsalomonis/Desktop/dataAnalysis/qPCR/ExpressionInput/7GeneModels.txt"
    #python LineageProfilerIterate.py --i "/Users/nsalomonis/Desktop/dataAnalysis/qPCR/deltaCT/LabMeeting/ExpressionInput/exp.ABI_PediatricSNS.txt" --r "/Users/nsalomonis/Desktop/dataAnalysis/qPCR/ExpressionOutput/MarkerFinder/MarkerFinder-ABI_PediatricSNS.txt" --s 4

    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a tab-delimited input expression file in the command-line"
        print 'Example: python LineageProfilerIterate.py --i "/Users/me/qPCR.txt" --r "/Users/me/reference.txt" --m "/Users/me/models.txt"'
    else:
        try:
            options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','species=','o=','platform=','codingtype=',
                                                    'compendium_platform=','r=','m=','v=','s=','permute=','useMulti=',
                                                    'cross_validation=','setsToOutput='])
        except Exception,e:
            print ""
        for opt, arg in options:
            if opt == '--i': exp_input=arg
            elif opt == '--o': exp_output=arg
            elif opt == '--platform': platform=arg
            elif opt == '--codingtype': codingtype=arg
            elif opt == '--compendium_platform': compendium_platform=arg
            elif opt == '--r': customMarkers=arg
            elif opt == '--m': geneModels=arg
            elif opt == '--v': vendor=arg
            elif opt == '--permute': permute=True
            elif opt == '--useMulti': useMulti=True
            elif opt == '--cross_validation': cross_validation = True
            elif opt == '--setsToOutput': finalNumberSetsToOutput = int(arg)
            elif opt == '--s':
                try: modelSize = int(arg)
                except Exception:
                    modelSize = arg
                    if modelSize != 'optimize':
                        print 'Please specify a modelSize (e.g., 7-gene model search) as a single integer (e.g., 7)'
                        sys.exit()
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
            
        if exp_output == None: exp_output = exp_input
        
        if cross_validation == True:
            ### Generate 2/3rd and 1/3rd sets for testing and validation
            crossValidationAnalysis(species,platform,exp_input,exp_output,codingtype,compendium_platform,modelSize,
                        geneModels,permute,useMulti,finalNumberSetsToOutput)
            sys.exit()
            
        runLineageProfiler(species,platform,exp_input,exp_output,codingtype,compendium_platform,modelSize=modelSize,
                           customMarkers=customMarkers,geneModels=geneModels,permute=permute,useMulti=useMulti)