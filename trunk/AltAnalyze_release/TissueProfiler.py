###TissueProfiler
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
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

import sys, string
import statistics
import math
import os.path
import unique
import copy
import time
import export
import salstat_stats; reload(salstat_stats)

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

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
            
######### Below code deals is specific to this module #########
def runTissueProfiler(species,array_type,exp_input,exp_output,codingtype,compendium_platform):
    global exp_output_file; exp_output_file = exp_output; global targetPlatform
    global tissue_specific_db; global expession_subset; global tissues; global sample_headers
    global analysis_type; global coding_type; coding_type = codingtype
    global tissue_to_gene; tissue_to_gene = {}
    
    global correlate_by_order; correlate_by_order = 'no'
    global rho_threshold; rho_threshold = -1
    global correlate_to_tissue_specific; correlate_to_tissue_specific = 'no'
    
    tissue_specific_db={}; expession_subset=[]; sample_headers=[]; tissues=[]
    if 'RawSplice' in exp_input or 'FullDatasets' in exp_input:
        analysis_type = 'exonLevel'
        if array_type != 'exon': ### If the input IDs are not Affymetrix Exon 1.0 ST probesets, then translate to the appropriate system
            translate_to_genearray = 'no'
            targetPlatform = compendium_platform
            translation_db = importExonIDTranslations(array_type,species,translate_to_genearray)
            keyed_by = 'translation'
        else: translation_db=[]; keyed_by = 'primaryID'; targetPlatform = 'exon'
    elif array_type == "3'array" or array_type == 'AltMouse':
        ### Get arrayID to Ensembl associations
        translation_db = importGeneIDTranslations(exp_output)
        keyed_by = 'translation'
        targetPlatform = 'exon'
        analysis_type = 'geneLevel'
    else: translation_db=[]; keyed_by = 'primaryID'; targetPlatform = 'exon'; analysis_type = 'geneLevel'
    try: importTissueSpecificProfiles(species)
    except Exception:
        try:
            targetPlatform = 'gene'
            importTissueSpecificProfiles(species)
        except Exception: 
            targetPlatform = "3'array"
            importTissueSpecificProfiles(species)
            
    importGeneExpressionValues(exp_input,tissue_specific_db,keyed_by,translation_db)
    analyzeTissueSpecificExpressionPatterns()

def importTissueSpecificProfiles(species):
    if analysis_type == 'exonLevel':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_'+targetPlatform +'_tissue-specific_AltExon_protein_coding.txt'
    else:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_'+targetPlatform +'_tissue-specific_'+coding_type+'.txt'
    print 'Importing the tissue compedium database:',filename
    fn=filepath(filename); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            headers = t; x=1; index=0
            for i in headers:
                if 'UID' == i: ens_index = index
                if analysis_type == 'exonLevel': ens_index = ens_index ### Assigned above when analyzing probesets
                elif 'Ensembl' in i: ens_index = index
                if 'marker-in' in i: tissue_index = index+1; marker_in = index
                index+=1
            for i in t[tissue_index:]: tissues.append(i)
        else:
            gene = string.split(t[ens_index],'|')[0] ### Only consider the first listed gene - this gene is the best option based on ExpressionBuilder rankings
            #if 'Pluripotent Stem Cells' in t[marker_in] or 'Heart' in t[marker_in]:
            tissue_exp = map(float, t[tissue_index:])
            tissue_specific_db[gene]=x,tissue_exp ### Use this to only grab relevant gene expression profiles from the input dataset
            x+=1
    print len(tissue_specific_db), 'genes in the tissue compendium database'

    if correlate_to_tissue_specific == 'yes':
        try: importTissueCorrelations(filename)
        except Exception:
            null=[]
            #print '\nNo tissue-specific correlations file present. Skipping analysis.'; kill
        
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
                
def importGeneExpressionValues(filename,tissue_specific_db,keyed_by,translation_db):
    ### Import gene-level expression raw values           
    fn=filepath(filename); x=0; genes_added={}; gene_expression_db={}
    print 'importing:',filename
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        
        if x==0:
            if '#' not in data:
                for i in t[1:]: sample_headers.append(i)
                x=1
        else:
            gene = t[0]
            #if '-' not in gene and ':E' in gene: print gene;sys.exit()
            if analysis_type == 'exonLevel':
                try: ens_gene,exon = string.split(gene,'-')[:2]
                except Exception: exon = gene
                gene = exon
            if keyed_by == 'translation': ### alternative value is 'primaryID'
                """if gene == 'ENSMUSG00000025915-E19.3':
                    for i in translation_db: print [i], len(translation_db); break
                    print gene, [translation_db[gene]];sys.exit()"""
                try: gene = translation_db[gene] ### Ensembl annotations
                except Exception: gene = 'null'
            try:
                index,tissue_exp=tissue_specific_db[gene]
                try: genes_added[gene]+=1
                except Exception: genes_added[gene]=1
                exp_vals = map(float, t[1:])
                gene_expression_db[gene] = [index,exp_vals]
            except Exception: null=[] ### gene is not a tissue specific marker

    print len(gene_expression_db), 'matching genes in the dataset and tissue compendium database'
    
    for gene in genes_added:
        if genes_added[gene]>1: del gene_expression_db[gene] ### delete entries that are present in the input set multiple times (not trustworthy)
        else: expession_subset.append(gene_expression_db[gene]) ### These contain the rank order and expression
    #print len(expession_subset);sys.exit()
    expession_subset.sort() ### This order now matches that of 
    gene_expression_db=[]
    
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
        gene_translation_db = gene_translation_db2
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_'+array_type+'-exon_probesets.txt'
    ### Import exon array to target platform translations (built for DomainGraph visualization)
    fn=filepath(filename); x=0; translation_db={}
    print 'Importing the translation file',fn
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
    #print len(translation_db)
    return translation_db

def analyzeTissueSpecificExpressionPatterns():
    tissue_specific_sorted = []; genes_present={}; tissue_exp_db={}; gene_order_db={}; gene_order=[]
    for (index,vals) in expession_subset: genes_present[index]=[]
    for gene in tissue_specific_db:
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

    ### Organize sample expression, with the same gene order as the tissue expression set
    sample_exp_db={}
    for (index,exp_vals) in expession_subset:
        i=0
        for f in exp_vals:
            ### The order of the tissue specific expression profiles is based on the import gene order
            try: sample_exp_db[sample_headers[i]].append(f)
            except Exception: sample_exp_db[sample_headers[i]] = [f]
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
        
    PearsonCorrelationAnalysis(sample_exp_db,tissue_exp_db)
    sample_exp_db=[]; tissue_exp_db=[]
    exportCorrelationResults()

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
    print "Analyzing"; k=1
    for tissue in tissue_exp_db:
        print k,"of",len(tissue_exp_db),"classifier tissue/cell-types"; k+=1
        tissue_expression_list = tissue_exp_db[tissue]
        for sample in sample_exp_db:
            if correlate_to_tissue_specific == 'yes':
                ### Keyed by tissue specific sample profiles
                sample_expression_list = sample_exp_db[sample][tissue] ### dictionary as the value for sample_exp_db[sample]
                #print tissue, sample_expression_list
                #print tissue_expression_list; sys.exit()
            else: sample_expression_list = sample_exp_db[sample]
            rho = pearson(tissue_expression_list,sample_expression_list)
            #tst = salstat_stats.TwoSampleTests(tissue_expression_list,sample_expression_list)
            #pp,pr = tst.PearsonsCorrelation()
            #sp,sr = tst.SpearmansCorrelation()
            #print tissue, sample
            #if rho>.5: print [rho, pr, sr],[pp,sp];sys.exit()
            #if rho<.5: print [rho, pr, sr],[pp,sp];sys.exit()
            try: tissue_comparison_scores[tissue].append([rho,sample])
            except Exception: tissue_comparison_scores[tissue] = [[rho,sample]]
    sample_exp_db=[]; tissue_exp_db=[]
    
def pearson(array1,array2):
    item = 0; sum_a = 0; sum_b = 0; sum_c = 0    
    while item < len(array1):
        a = (array1[item] - avg(array1))*(array2[item] - avg(array2))
        b = math.pow((array1[item] - avg(array1)),2)
        c = math.pow((array2[item] - avg(array2)),2)        
        sum_a = sum_a + a
        sum_b = sum_b + b
        sum_c = sum_c + c
        item = item + 1
    r = sum_a/math.sqrt(sum_b*sum_c)
    return r

def avg(array):
    return sum(array)/len(array)

def exportCorrelationResults():
    corr_output_file = string.replace(exp_output_file,'DATASET','TissueCorrelations')
    corr_output_file = string.replace(corr_output_file,'.txt','-'+coding_type+'.txt')
    if analysis_type == 'exonLevel':
        corr_output_file = string.replace(corr_output_file,coding_type,'AltExon')
    data = export.ExportFile(corr_output_file)
    ### Make title row
    headers=['Sample_name']
    for tissue in tissue_comparison_scores:
        for (r,sample) in tissue_comparison_scores[tissue]: headers.append(sample)
        break
    title_row = string.join(headers,'\t')+'\n'; data.write(title_row)
    ### Export correlation data
    for tissue in tissue_comparison_scores:
        values=[tissue]
        for (r,sample) in tissue_comparison_scores[tissue]: values.append(str(r))
        values = string.join(values,'\t')+'\n'; data.write(values)
    data.close()
    print corr_output_file,'exported...'

if __name__ == '__main__':
    species = 'Hs'
    array_type = "RNASeq"
    compendium_platform = "3'array"
    compendium_platform = "gene"
    compendium_platform = "exon"
    #array_type = "junction"
    codingtype = 'ncRNA'
    codingtype = 'protein_coding'
    
    exp_input = 'C:/Users/Nathan Salomonis/Desktop/temp2/ExpressionInput/exp.meta-steady-state.txt'
    exp_output = 'C:/Users/Nathan Salomonis/Desktop/temp2/ExpressionOutput/DATASET-Hs-Exon-meta.txt'

    exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/ExonArray/CP-hESC/ExpressionInput/exp.hESC_differentiation-steady-state.txt'
    exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/ExonArray/CP-hESC/ExpressionOutput/DATASET-hESC_differentiation.txt'

    exp_input = 'C:/Users/Nathan Salomonis/Desktop/Eugene-Custom array/ExpressionInput/exp.test-steady-state.txt'
    exp_output = 'C:/Users/Nathan Salomonis/Desktop/Eugene-Custom array/ExpressionOutput/DATASET-test.txt'

    exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/JunctionArray/Lilach/ExpressionInput/exp.test-steady-state.txt'  
    exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/JunctionArray/Lilach/ExpressionOutput/DATASET-test.txt'
    
    exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/hESC-NP/TopHat-hESC_differentiation/AltResults/RawSpliceData/Hs/splicing-index/test.txt'
    exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/hESC-NP/TopHat-hESC_differentiation/ExpressionOutput/DATASET-test.txt'
    
    exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/ExonArray/CP-hESC/AltExpression/FullDatasets/ExonArray/Hs/hESC_differentiation.txt'
    exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/ExonArray/CP-hESC/ExpressionOutput/DATASET-hESC_differentiation.txt'

    exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/Combined-GSE14588_RAW/junction/AltExpression/FullDatasets/JunctionArray/Hs/test.txt'
    exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/Combined-GSE14588_RAW/junction/ExpressionOutput/DATASET-test.txt'
                
    exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/Combined-GSE14588_RAW/junction/ExpressionInput/exp.test-steady-state.txt'
    exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/Combined-GSE14588_RAW/junction/ExpressionOutput/DATASET-test.txt'

    exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/Combined-GSE14588_RAW/junction/AltResults/RawSpliceData/Hs/splicing-index/Hs_Junction_d14_vs_d7.p5_average.txt'
    exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/Combined-GSE14588_RAW/junction/ExpressionOutput/DATASET-test.txt'

    exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/hESC-NP/TopHat-hESC_differentiation/AltExpression/FullDatasets/RNASeq/Hs/test.txt'
    exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/hESC-NP/TopHat-hESC_differentiation/ExpressionOutput/DATASET-test.txt'

    exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/hESC-NP/TopHat-hESC_differentiation/ExpressionInput/exp.test-steady-state.txt'
    exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/hESC-NP/TopHat-hESC_differentiation/ExpressionOutput/DATASET-test.txt'
    
    #exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/r4_Bruneau_TopHat/ExpressionInput/exp.test-steady-state.txt'
    #exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/r4_Bruneau_TopHat/ExpressionOutput/DATASET-test.txt'  

    #exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/r4_Bruneau_TopHat/AltExpression/FullDatasets/RNASeq/Mm/test.txt'
    #exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/r4_Bruneau_TopHat/ExpressionOutput/DATASET-test.txt'  

    #exp_input = 'C:/Users/Nathan Salomonis/Desktop/Mouse-Compendium-AVERAGE/MoEx-AVERAGE-meta.txt'
    #exp_output = 'C:/Users/Nathan Salomonis/Desktop/Mouse-Compendium-AVERAGE/DATASET-meta.txt'  

    #exp_input = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/Amy-RNASeq/ExpressionInput/exp.miR1_KO-steady-state.txt'
    #exp_output = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/Amy-RNASeq/ExpressionOutput/DATASET-miR1_KO.txt'
    
    runTissueProfiler(species,array_type,exp_input,exp_output,codingtype,compendium_platform)