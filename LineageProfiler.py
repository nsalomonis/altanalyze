###LineageProfiler
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
from stats_scripts import statistics
import math
import os.path
import unique
import copy
import time
import export
import traceback
import warnings
#from stats_scripts import salstat_stats; reload(salstat_stats)
try:
    from scipy import stats
    use_scipy = True
    import numpy
except Exception:
    use_scipy = False ### scipy is not required but is used as a faster implementation of Fisher Exact Test when present

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
def runLineageProfiler(species,array_type,exp_input,exp_output,codingtype,compendium_platform,customMarkers=False):
    """
    print species
    print array_type
    print export.findFilename(exp_input)
    print export.findFilename(exp_output)
    print codingtype
    print compendium_platform
    print customMarkers
    """
    
    global exp_output_file; exp_output_file = exp_output; global targetPlatform
    global tissue_specific_db; global expession_subset; global tissues; global sample_headers
    global analysis_type; global coding_type; coding_type = codingtype
    global tissue_to_gene; tissue_to_gene = {}; global platform; global cutoff
    global customMarkerFile; global keyed_by; global compendiumPlatform
    customMarkerFile = customMarkers; compendiumPlatform = compendium_platform
    
    global correlate_by_order; correlate_by_order = 'no'
    global rho_threshold; rho_threshold = -1
    global correlate_to_tissue_specific; correlate_to_tissue_specific = 'no'
    platform = array_type
    cutoff = 0.01
    global value_type
    global missingValuesPresent
    
    if 'stats.' in exp_input:
        value_type = 'calls'
    else:
        value_type = 'expression'
    
    tissue_specific_db={}; expession_subset=[]; sample_headers=[]; tissues=[]
    if len(array_type)==2:
        ### When a user-supplied expression is provided (no ExpressionOutput files provided - importGeneIDTranslations)
        vendor, array_type = array_type
        platform = array_type
    else: vendor = 'Not needed'
    
    if 'other:' in vendor:
        vendor = string.replace(vendor,'other:','')
        array_type = "3'array"
    if 'RawSplice' in exp_input or 'FullDatasets' in exp_input or coding_type == 'AltExon':
        analysis_type = 'AltExon'
        if platform != compendium_platform: ### If the input IDs are not Affymetrix Exon 1.0 ST probesets, then translate to the appropriate system
            translate_to_genearray = 'no'
            targetPlatform = compendium_platform
            translation_db = importExonIDTranslations(array_type,species,translate_to_genearray)
            keyed_by = 'translation'
        else: translation_db=[]; keyed_by = 'primaryID'; targetPlatform = compendium_platform
    elif array_type == "3'array" or array_type == 'AltMouse':
        ### Get arrayID to Ensembl associations
        if vendor != 'Not needed':
            ### When no ExpressionOutput files provided (user supplied matrix)
            translation_db = importVendorToEnsemblTranslations(species,vendor,exp_input)
        else:
            try: translation_db = importGeneIDTranslations(exp_output)
            except: translation_db = importVendorToEnsemblTranslations(species,'Symbol',exp_input)
                
        keyed_by = 'translation'
        targetPlatform = compendium_platform
        analysis_type = 'geneLevel'
    else:
        translation_db=[]; keyed_by = 'primaryID'; targetPlatform = compendium_platform; analysis_type = 'geneLevel'

    if compendium_platform == "3'array" and array_type != "3'array":
        keyed_by = 'ensembl' ### ensembl is not indicated anywhere but avoides key by primaryID and translation -> works for RNASeq
        
    targetPlatform = compendium_platform ### Overides above
    
    """ Determine if a PSI file with missing values """
    if vendor == 'PSI':
        missingValuesPresent = True
    else:
        missingValuesPresent = importTissueSpecificProfiles(species,checkForMissingValues=True)
        
    try: importTissueSpecificProfiles(species)
    except Exception:
        try:
            try:
                targetPlatform = 'exon'
                importTissueSpecificProfiles(species)
            except Exception:
                try:
                    targetPlatform = 'gene'
                    importTissueSpecificProfiles(species)
                except Exception: 
                    targetPlatform = "3'array"
                    importTissueSpecificProfiles(species)
        except Exception,e:
            print traceback.format_exc()
            print 'No compatible compendiums present...'
            forceTissueSpecificProfileError
            
    try: importGeneExpressionValues(exp_input,tissue_specific_db,translation_db,species=species)
    except:
        print "Changing platform to 3'array"
        array_type = "3'array"
        exp_input = string.replace(exp_input,'-steady-state.txt','.txt')
        importGeneExpressionValues(exp_input,tissue_specific_db,translation_db,species=species)
    ### If the incorrect gene system was indicated re-run with generic parameters

    if len(expession_subset)==0 and (array_type == "3'array" or array_type == 'AltMouse' or array_type == 'Other'):
        translation_db=[]; keyed_by = 'primaryID'; targetPlatform = compendium_platform; analysis_type = 'geneLevel'
        tissue_specific_db={}
        try: importTissueSpecificProfiles(species)
        except Exception:
            try: targetPlatform = 'exon'; importTissueSpecificProfiles(species)
            except Exception:
                try: targetPlatform = 'gene'; importTissueSpecificProfiles(species)
                except Exception: targetPlatform = "3'array"; importTissueSpecificProfiles(species)
        importGeneExpressionValues(exp_input,tissue_specific_db,translation_db,species=species)
    zscore_output_dir = analyzeTissueSpecificExpressionPatterns(expInput=exp_input)

    return zscore_output_dir

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
    from import_scripts import BuildAffymetrixAssociations
    
    ### Use the same annotation method that is used to create the ExpressionOutput annotations
    use_go = 'yes'
    conventional_array_db={}
    conventional_array_db = BuildAffymetrixAssociations.getUIDAnnotationsFromGOElite(conventional_array_db,species,vendor,use_go)
    for arrayid in conventional_array_db:
        ca = conventional_array_db[arrayid]
        ens = ca.Ensembl()
        try: translation_db[arrayid] = ens[0] ### This first Ensembl is ranked as the most likely valid based on various metrics in getArrayAnnotationsFromGOElite
        except Exception: None
    
    return translation_db

def importTissueSpecificProfiles(species,checkForMissingValues=False):
    if analysis_type == 'AltExon':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_'+targetPlatform +'_tissue-specific_AltExon_protein_coding.txt'
    else:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_'+targetPlatform +'_tissue-specific_'+coding_type+'.txt'
    if customMarkerFile != False and customMarkerFile != None:
        if len(customMarkerFile)>0:
            filename = customMarkerFile

    #filename = 'AltDatabase/ensembl/'+species+'/random.txt'
    #print 'Target platform used for analysis:',species, targetPlatform, coding_type
    if value_type == 'calls':
        filename = string.replace(filename,'.txt','_stats.txt')
    fn=filepath(filename); x=0
    tissue_index = 1
    tissues_added={}
    missing_values_present = False
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            print 'Importing the tissue compedium database:',export.findFilename(filename)
            headers = t; x=1; index=0
            for i in headers:
                if 'UID' == i: ens_index = index; uid_index = index
                if analysis_type == 'AltExon': ens_index = ens_index ### Assigned above when analyzing probesets
                elif 'Ensembl' in i: ens_index = index
                if 'marker-in' in i: tissue_index = index+1; marker_in = index
                index+=1
            try:
                for i in t[tissue_index:]: tissues.append(i)
            except Exception:
                for i in t[1:]: tissues.append(i)
            if keyed_by == 'primaryID':
                try: ens_index = uid_index
                except Exception: None
        else:
            try:
                gene = t[0]
                try: gene = string.split(gene,'|')[0] ### Only consider the first listed gene - this gene is the best option based on ExpressionBuilder rankings
                except Exception: pass
                tissue_exp = map(float, t[1:])
                tissue_specific_db[gene]=x,tissue_exp ### Use this to only grab relevant gene expression profiles from the input dataset
            except Exception:
                try: gene = string.split(gene,'|')[0] ### Only consider the first listed gene - this gene is the best option based on ExpressionBuilder rankings
                except Exception: pass
                #if 'Pluripotent Stem Cells' in t[marker_in] or 'Heart' in t[marker_in]:
                #if t[marker_in] not in tissues_added: ### Only add the first instance of a gene for that tissue - used more for testing to quickly run the analysis
                tissue_exp = t[tissue_index:]
                if '' in tissue_exp:
                    missing_values_present = True
                    ### If missing values present (PSI values)
                    tissue_exp = ['0.000101' if i=='' else i for i in tissue_exp]
                tissue_exp = map(float,tissue_exp)
                if value_type == 'calls':
                    tissue_exp = produceDetectionCalls(tissue_exp,platform) ### 0 or 1 calls
                tissue_specific_db[gene]=x,tissue_exp ### Use this to only grab relevant gene expression profiles from the input dataset
                try: tissues_added[t[marker_in]]=[] ### Not needed currently
                except Exception: pass
            x+=1
    print len(tissue_specific_db), 'genes in the tissue compendium database'

    if correlate_to_tissue_specific == 'yes':
        try: importTissueCorrelations(filename)
        except Exception:
            null=[]
            #print '\nNo tissue-specific correlations file present. Skipping analysis.'; kill
    if checkForMissingValues:
        return missing_values_present
        
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
        
def importGeneExpressionValues(filename,tissue_specific_db,translation_db,useLog=False,previouslyRun=False,species=None):
    ### Import gene-level expression raw values           
    fn=filepath(filename); x=0; genes_added={}; gene_expression_db={}
    dataset_name = export.findFilename(filename)
    max_val=0
    print 'importing:',dataset_name
    
    try:
        import gene_associations, OBO_import
        gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
        symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
    except Exception: symbol_to_gene={}
    
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        
        if x==0:
            if '#' not in data:
                for i in t[1:]: sample_headers.append(i)
                x=1
        else:
            gene = t[0]
            try: gene = string.split(t[0],'|')[0]
            except Exception: pass
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
                except Exception: pass
            try: gene = symbol_to_gene[gene][0] ### If RNASeq is the selected platform and Symbol is the uid
            except Exception: pass
            if gene in tissue_specific_db:
                index,tissue_exp=tissue_specific_db[gene]
                try: genes_added[gene]+=1
                except Exception: genes_added[gene]=1
                proceed=True
                try:
                    exp_vals = t[1:]
                    if '' in exp_vals:
                        ### If missing values present (PSI values)
                        exp_vals = ['0.000101' if i=='' else i for i in exp_vals]
                        useLog = False
                    exp_vals = map(float, exp_vals)
                    if platform == 'RNASeq':
                        if max(exp_vals)>max_val: max_val = max(exp_vals)
                        #if max(exp_vals)<3: proceed=False
                        if useLog==False:
                            exp_vals = map(lambda x: math.log(x+1,2),exp_vals)
                    if value_type == 'calls': ### Hence, this is a DABG or RNA-Seq expression
                        exp_vals = produceDetectionCalls(exp_vals,targetPlatform) ### 0 or 1 calls
                    if proceed:
                        gene_expression_db[gene] = [index,exp_vals]
                except Exception:
                    print 'Non-numeric values detected:'
                    x = 5
                    print t[:x]
                    while x < t:
                        t[x:x+5]
                        x+=5
                    print 'Formatting error encountered in:',dataset_name; forceError
            """else:
                for gene in tissue_specific_db:
                    if 'Ndufa9:ENSMUSG00000000399:I2.1-E3.1' in gene:
                        print gene, 'dog';sys.exit()
                print gene;kill"""
        
    print len(gene_expression_db), 'matching genes in the dataset and tissue compendium database'
    
    for gene in genes_added:
        if genes_added[gene]>1:
            del gene_expression_db[gene] ### delete entries that are present in the input set multiple times (not trustworthy)
        else: expession_subset.append(gene_expression_db[gene]) ### These contain the rank order and expression
    #print len(expession_subset);sys.exit()
    expession_subset.sort() ### This order now matches that of 
    gene_expression_db=[]
    
    if max_val<20 and platform == 'RNASeq' and previouslyRun==False: ### Only allow to happen once
        importGeneExpressionValues(filename,tissue_specific_db,translation_db,useLog=True,previouslyRun=True,species=species)

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
        print 'Importing the translation file',export.findFilename(fn)
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

def analyzeTissueSpecificExpressionPatterns(expInput=None):
    tissue_specific_sorted = []; genes_present={}; tissue_exp_db={}; gene_order_db={}; gene_order=[]
    for (index,vals) in expession_subset: genes_present[index]=[]
    for gene in tissue_specific_db:
        tissue_specific_sorted.append(tissue_specific_db[gene])
        # tissue_specific_db[gene][1]
        #print tissue_specific_db[gene][1].count(0.000101);sys.exit()
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
            
        except Exception:
            #print gene;sys.exit()
            null=[] ### Gene is not present in the input dataset

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
    zscore_output_dir = exportCorrelationResults(expInput)
    return zscore_output_dir

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
    print "Beginning LineageProfiler analysis"; k=0
    original_increment = int(len(tissue_exp_db)/15.00); increment = original_increment
    p = 1 ### Default value if not calculated
    for tissue in tissue_exp_db:
        #print k,"of",len(tissue_exp_db),"classifier tissue/cell-types"
        if k == increment: increment+=original_increment; print '*',
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
                if missingValuesPresent:
                    ### For PSI values

                    tissue_expression_list = numpy.ma.masked_values(tissue_expression_list,0.000101)
                    #tissue_expression_list = numpy.ma.array([numpy.nan if i==0.000101 else i for i in tissue_expression_list])
                    sample_expression_list = numpy.ma.masked_values(sample_expression_list,0.000101)
                    #tissue_expression_list = numpy.ma.array([numpy.nan if i==0.000101 else i for i in tissue_expression_list])
                    
                    updated_tissue_expression_list=[]
                    updated_sample_expression_list=[]
                    i=0
                    
                    coefr=numpy.ma.corrcoef(tissue_expression_list,sample_expression_list)
                    rho = coefr[0][1]
                    """
                    if sample == 'Cmp.21':
                        #print rho
                        #print tissue_expression_list[:10]
                        #print string.join(map(str,tissue_expression_list[:20]),'\t')
                        #print sample_expression_list[:10]
                        #print string.join(map(str,sample_expression_list[:20]),'\t')
                        #coefr=numpy.ma.corrcoef(numpy.array(tissue_expression_list[:10]),numpy.array(sample_expression_list[:10]))
                        print tissue, sample, rho, len(tissue_expression_list), len(sample_expression_list)
                        """
                else:
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides import warnings
                        rho,p = stats.pearsonr(tissue_expression_list,sample_expression_list)
            except Exception:
                #print traceback.format_exc(); sys.exit()
                ### simple pure python implementation - no scipy required (not as fast though and no p-value)
                rho = pearson(tissue_expression_list,sample_expression_list)
            #tst = salstat_stats.TwoSampleTests(tissue_expression_list,sample_expression_list)
            #pp,pr = tst.PearsonsCorrelation()
            #sp,sr = tst.SpearmansCorrelation()
            #print tissue, sample
            #if rho>.5: print [rho, pr, sr],[pp,sp];sys.exit()
            #if rho<.5: print [rho, pr, sr],[pp,sp];sys.exit()
            try: tissue_comparison_scores[tissue].append([rho,p,sample])
            except Exception: tissue_comparison_scores[tissue] = [[rho,p,sample]]
    sample_exp_db=[]; tissue_exp_db=[]
    print 'Correlation analysis finished'
    
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

def adjustPValues():
    """ Can be applied to calculate an FDR p-value on the p-value reported by scipy.
        Currently this method is not employed since the p-values are not sufficiently
        stringent or appropriate for this type of analysis """
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

def replacePearsonPvalueWithZscore():
    all_sample_data={}
    for tissue in tissue_comparison_scores:
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            all_sample_data[sample] = [] ### populate this dictionary and create sub-dictionaries
        break

    for tissue in tissue_comparison_scores:
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            all_sample_data[sample].append(r)

    sample_stats={}
    all_dataset_rho_values=[]
    ### Get average and standard deviation for all sample rho's
    for sample in all_sample_data:
        all_dataset_rho_values+=all_sample_data[sample]
        avg=statistics.avg(all_sample_data[sample])
        stdev=statistics.stdev(all_sample_data[sample])
        sample_stats[sample]=avg,stdev
    
    global_rho_avg = statistics.avg(all_dataset_rho_values)
    global_rho_stdev = statistics.stdev(all_dataset_rho_values)
    
    ### Replace the p-value for each rho
    for tissue in tissue_comparison_scores:
        scores = []
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            #u,s=sample_stats[sample]
            #z = (r-u)/s
            z = (r-global_rho_avg)/global_rho_stdev ### Instead of doing this for the sample background, do it relative to all analyzed samples
            scores.append([r,z,sample])
        tissue_comparison_scores[tissue] = scores

def exportCorrelationResults(exp_input):
    input_file = export.findFilename(exp_input)
    if '.txt' in exp_output_file:
        corr_output_file = string.replace(exp_output_file,'DATASET','LineageCorrelations')
    else: ### Occurs when processing a non-standard AltAnalyze file
        corr_output_file = exp_output_file+'/'+input_file
    corr_output_file = string.replace(corr_output_file,'.txt','-'+coding_type+'-'+compendiumPlatform+'.txt')
    if analysis_type == 'AltExon':
        corr_output_file = string.replace(corr_output_file,coding_type,'AltExon')
    filename = export.findFilename(corr_output_file)
    score_data = export.ExportFile(corr_output_file)
    if use_scipy:
        zscore_output_dir = string.replace(corr_output_file,'.txt','-zscores.txt')
        probability_data = export.ExportFile(zscore_output_dir)
        #adjustPValues()
        replacePearsonPvalueWithZscore()
    ### Make title row
    headers=['Sample_name']
    for tissue in tissue_comparison_scores:
        for (r,p,sample) in tissue_comparison_scores[tissue]: headers.append(sample)
        break
    title_row = string.join(headers,'\t')+'\n'
    score_data.write(title_row)
    if use_scipy:
        probability_data.write(title_row)
    ### Export correlation data
    tissue_scores = {}; tissue_probabilities={}; tissue_score_list = [] ### store and rank tissues according to max(score)
    for tissue in tissue_comparison_scores:
        scores=[]
        probabilities=[]
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            scores.append(r)
            probabilities.append(p)
        tissue_score_list.append((max(scores),tissue))
        tissue_scores[tissue] = string.join(map(str,[tissue]+scores),'\t')+'\n' ### export line
        if use_scipy:
            tissue_probabilities[tissue] = string.join(map(str,[tissue]+probabilities),'\t')+'\n'
        
    tissue_score_list.sort()
    tissue_score_list.reverse()
    for (score,tissue) in tissue_score_list:
        score_data.write(tissue_scores[tissue])
        if use_scipy:
            probability_data.write(tissue_probabilities[tissue])
    score_data.close()
    if use_scipy:
        probability_data.close()
    print filename,'exported...'
    return zscore_output_dir

def visualizeLineageZscores(zscore_output_dir,grouped_lineage_zscore_dir,graphic_links):
    from visualization_scripts import clustering
    ### Perform hierarchical clustering on the LineageProfiler Zscores
    graphic_links = clustering.runHCOnly(zscore_output_dir,graphic_links)   
    return graphic_links
    
if __name__ == '__main__':
    species = 'Hs'
    array_type = "3'array"
    vendor = 'Affymetrix'
    vendor = 'other:Symbol'
    vendor = 'other:Ensembl'
    #vendor = 'RNASeq'
    array_type = "exon"
    #array_type = "3'array"
    #array_type = "RNASeq"
    compendium_platform = "3'array"
    compendium_platform = "exon"
    #compendium_platform = "gene"
    #array_type = "junction"
    codingtype = 'ncRNA'
    codingtype = 'protein_coding'
    #codingtype = 'AltExon'
    array_type = vendor, array_type

    exp_input = "/Users/saljh8/Documents/1-conferences/GE/LineageMarkerAnalysis/Synapse-ICGS-EB-Ensembl.txt"
    exp_output = "/Users/saljh8/Documents/1-conferences/GE/LineageMarkerAnalysis/temp.txt"
    #customMarkers = "/Users/nsalomonis/Desktop/dataAnalysis/qPCR/PAM50/AltAnalyze/ExpressionOutput/MarkerFinder/AVERAGE-training.txt"
    customMarkers = False

    runLineageProfiler(species,array_type,exp_input,exp_output,codingtype,compendium_platform,customMarkers)