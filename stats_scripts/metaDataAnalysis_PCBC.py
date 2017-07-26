import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import re
import unique
import export
import math
from stats_scripts import statistics
import traceback

"""
    This script takes existing formatted metadata (C4 PCBC approved fields) and filters them to determine unique and non-unique
    donors for a specific covariate and derives comparison and group relationships to extract from an existing expression file
"""

def BatchCheck(sample_id,nonpreferential_batchs,preferential_samples,platform):
    priority=0
    for batch in nonpreferential_batchs:
        if batch in sample_id:
            priority=1
    if platform == 'RNASeq' or platform == 'exon':
        if sample_id not in preferential_samples: priority = 1
        elif sample_id in preferential_samples: priority = 0
    return priority
    
class MetaData:
    def __init__(self,cellLine,donor_id,sex,diffState,quality,coi,vector,genes,lab,public,priority):
        self.donor_id = donor_id; self.sex = sex; self.diffState = diffState; self.quality = quality; self.coi = coi
        self.cellLine = cellLine; self.vector = vector; self.genes = genes ; self.lab = lab; self.public = public
        self.priority = priority
    def CellLine(self): return self.cellLine
    def Gender(self): return self.sex
    def DiffState(self): return self.diffState
    def Quality(self): return self.quality
    def COI(self): return self.coi
    def Vector(self): return self.vector
    def Genes(self): return self.genes
    def Lab(self): return self.lab
    def Public(self): return self.public
    def Priority(self): return self.priority
    def __repr__(self): return self.CellLine()
    
def prepareComparisonData(input_file,diffStateQuery,CovariateQuery,uniqueDonors,genderRestricted,platform=None,compDiffState=None,restrictCovariateTerm=None):
    removeConfoundedCovariates=True
    firstLine = True
    notation_db={}; donor_sex_db={}
    failed_QC = ['FAIL','bad','EXCLUDE']
    nonpreferential_batchs = ['.144.7.','.144.6.','.219.2','.219.5','H9'] ### used when removing non-unique donors
    preferential_samples = ['SC11-004A.133.1.7', 'SC11-008B.149.3.14', 'SC11-010A.149.5.16', 'SC11-010B.149.5.18', 'SC11-012A.149.5.19', 'SC11-012B.149.5.20', 'SC11-013A.149.5.21', 'SC11-013B.149.5.22', 'SC12-002A.154.1.4', 'SC12-005A.144.7.19', 'SC12-005B.144.7.21', 'SC12-007.181.7.1', 'SC13-043.420.12.3', 'SC13-044.219.2.9', 'SC13-045.219.5.10', 'SC14-066.558.12.18', 'SC14-067.558.12.19', 'SC14-069.569.12.25', 'IPS18-4-1.102.2.2', 'IPS18-4-2.102.2.4', 'SC11-005B.119.5.9', 'SC11-006A.119.1.4', 'SC11-006B.119.1.10', 'SC11-007A.119.3.5', 'SC11-007B.119.3.1', 'SC11-014A.133.1.13', 'SC11-014B.133.2.4', 'SC11-015A.133.1.14', 'SC11-015B.133.2.5', 'SC11-016A.133.1.8', 'SC11-016B.133.2.15', 'SC11-017A.144.6.16', 'SC11-017B.154.1.2', 'SC11-018A.144.6.18', 'SC11-018B.154.1.3', 'SC12-006A.144.7.22', 'SC12-006B.144.7.23', 'SC12-019.181.7.2', 'SC12-020.181.7.3', 'SC12-022A.172.5.8', 'SC12-024.219.2.7', 'SC12-025A.172.5.9', 'SC12-028.181.7.5', 'SC12-029.181.7.6', 'SC12-030.181.7.4', 'SC12-031.181.7.7', 'SC12-034.182.1.7', 'SC12-035.569.12.16', 'SC12-036.182.2.20', 'SC12-037.182.1.12', 'SC12-038.420.12.1', 'SC13-049.219.2.8']
    preferential_samples += ['SC11-010BEB.144.6.6', 'SC13-045EB.219.6.11', 'SC12-005EB.585.2.13', 'SC11-013EB.558.12.7', 'SC13-043BEB.419.12.13', 'SC14-067EB.558.12.1', 'SC13-044EB.219.6.10', 'SC11-012BEB.144.6.7', 'SC14-066EB.585.2.14'] # variable XIST 'SC11-009A.133.3.5', 'SC11-009A.149.5.15'
    preferential_samples += ['H9.119.3.7']
    unique_covariate_samples={}
    covariate_samples={}
    sample_metadata={}
    cellLineIDs={}
    uniqueDonor_db={}
    donor_db={}
    allCovariateSamples={}
    blood = ['CD34+ cells','mononuclear']
    integrating = ['lentivirus','retrovirus']
    nonintegrating = ['plasmid','RNA','Sendai Virus']

    for line in open(input_file,'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstLine:
            if 'Cell_Line_Type' in data:
                values = string.replace(data,'_',' ')
                values = string.replace(values,'"','')
                values = string.split(values,'\t')
            headers = values
            covariateIndex=None
            mergeBlood = False
            index=0
            for h in headers:
                if 'CellType' in h or 'Diffname short' in h: cellTypeIndex = index ### diff-state
                #if 'uid' in h: uidIndex = index
                if 'Sample' == h or 'Decorated Name' == h: uidIndex = index
                if 'CellLine' in h or 'C4 Cell Line ID' in h: cellLineIndex = index
                if 'Pass QC' in h: qualityIndex = index
                if 'Cell Type of Origin' in h: coiIndex = index
                if 'Reprogramming Vector Type' in h: vectorIndex = index
                if 'Reprogramming Gene Combination' in h: geneIndex = index
                if 'Gender' in h: sex_index = index
                if 'originating lab' in h or 'Originating Lab' in h: labIndex = index
                if 'High Confidence Donor ID (HCDID)' in h or 'Donor ID' == h:
                    unique_donor_index = index
                if 'ublic' in h: publicIndex = index
                if 'C4 Karyotype Result' in h: karyotypeIndex = index
                if 'Donor Life Stage' in h: donorStageIndex = index
                if 'Other Conditions During Reprogramming' in h: otherConditionIndex = index
                if 'Small Molecules' in h: smallMoleculeIndex = index
                if CovariateQuery == h: covariateIndex = index
                if 'UID' == h: uidIndexAlt = index 
                index+=1
            firstLine = False
            if CovariateQuery == 'Cell Type of Origin Combined' and covariateIndex==None:
                covariateIndex = headers.index('Cell Type of Origin')
                mergeBlood = True
            if CovariateQuery == 'originating lab' and covariateIndex==None:
                covariateIndex = headers.index('Originating Lab ID')
        else:
            try: sample_id = values[uidIndex]
            except Exception:
                uidIndex = uidIndexAlt
                sample_id = values[uidIndex]
            cellLine = values[cellLineIndex]
            try: donor_id = values[unique_donor_index]
            except Exception:
                #print values
                continue
                #print len(values), unique_donor_index;sys.exit()
            if donor_id == '': donor_id = cellLine
            sex = values[sex_index]
            diffState = values[cellTypeIndex]
            try: quality = values[qualityIndex]
            except Exception: quality = ''
            COI = values[coiIndex]
            vector = values[vectorIndex]
            genes = values[geneIndex]
            lab = values[labIndex]
            public = values[publicIndex]
            karyotype = values[karyotypeIndex]
            priority = BatchCheck(sample_id,nonpreferential_batchs,preferential_samples,platform)
            md = MetaData(cellLine,donor_id,sex,diffState,quality,COI,vector,genes,lab,public,priority)
            sample_metadata[sample_id] = md ### store all relevant metadata for future access (e.g., comparing omics types)
            try:  covariateType =  values[covariateIndex]
            except Exception: covariateType = None
            if CovariateQuery == 'integrating vectors': ### Build this covariate type interactively
                if values[vectorIndex] in integrating:
                    covariateType = 'integrating'
                elif values[vectorIndex] in nonintegrating:
                    covariateType = 'non-integrating'
                else:
                    covariateType = 'NA'
                #print covariateType, values[vectorIndex]
            #print covariateType, CovariateQuery, diffState, diffStateQuery, public, sample_id;sys.exit()
            if covariateType!=None:
                try: covariateType =  values[covariateIndex] ### e.g., COI
                except Exception: pass
                if mergeBlood:
                    if covariateType in blood: covariateType = 'blood'
                if covariateType == 'N/A': covariateType = 'NA'
                elif covariateType == '': covariateType = 'NA'
                elif '/' in covariateType: covariateType = string.replace(covariateType,'/','-')
                if (diffState==diffStateQuery or diffState==compDiffState or diffStateQuery=='NA') and (public == 'yes' or public == 'TRUE') and '_PE' not in sample_id and karyotype != 'abnormal':
                    if quality not in failed_QC:
                        proceed=True
                        if genderRestricted!=None:
                            if genderRestricted == sex: proceed = True
                            else: proceed = False
                        if restrictCovariateTerm != None:
                            if restrictCovariateTerm == covariateType:
                                proceed = True
                                if diffState==compDiffState:
                                    covariateType = compDiffState +' '+covariateType ### Make this a unique name
                                else:
                                    covariateType = diffStateQuery +' '+covariateType ### Make this a unique name
                            else:
                                proceed = False
                        if proceed:
                            try:
                                donor_db = unique_covariate_samples[covariateType]
                                try: donor_db[donor_id].append((priority,sample_id))
                                except Exception: donor_db[donor_id] = [(priority,sample_id)]
                            except Exception:
                                ### get unique donor samples sorted by priority
                                donor_db={donor_id:[(priority,sample_id)]}
                                unique_covariate_samples[covariateType] = donor_db
                            try: ### do the same for non-unique data
                                covariate_samples[covariateType].append(sample_id)
                            except Exception:
                                covariate_samples[covariateType] = [sample_id]
            else:
                if len(cellLine)>0 and diffState==diffStateQuery and public == 'yes' and '_PE' not in sample_id and quality not in failed_QC:
                    proceed=True
                    if genderRestricted!=None:
                        if genderRestricted == sex: proceed = True
                        else: proceed = False
                    if proceed:
                        try: cellLineIDs[cellLine].append((priority,sample_id))
                        except Exception: cellLineIDs[cellLine] = [(priority,sample_id)]
                        try: uniqueDonor_db[donor_id].append(cellLine)
                        except Exception: uniqueDonor_db[donor_id] = [cellLine]
                        try: donor_db[donor_id].append((priority,cellLine))
                        except Exception: donor_db[donor_id] = [(priority,cellLine)]
            ### Now, do an exhaustive search to exclude covariateType that are confounded by another covariate such that either could attribute to the effect
            ### Only do this relative to the CovariateQuery being compared
            if removeConfoundedCovariates:
                index=0
                for value in values:
                    header = headers[index]
                    if len(value)>0:
                        try: allCovariateSamples[header+':'+value].append(sample_id)
                        except Exception: allCovariateSamples[header+':'+value] = [sample_id]
                    index+=1

    if len(covariate_samples)==0:
        cellLinesPerUniqueDonor={}
        for cellLine in cellLineIDs:
            cellLineIDs[cellLine].sort()
            cellLineIDs[cellLine] = cellLineIDs[cellLine][0][1]
        for donor_id in uniqueDonor_db:
            for cellLine in uniqueDonor_db[donor_id]:
                cellLinesPerUniqueDonor[cellLine]=donor_id
        return cellLineIDs,sample_metadata,cellLinesPerUniqueDonor, donor_db

    if uniqueDonors:
        ### Select a single representative sample for each donor
        allUniqueDonors={}
        for covariateType in unique_covariate_samples:
            unique_donors=[]
            donor_db = unique_covariate_samples[covariateType]
            #print covariateType, donor_db
            for donor_id in donor_db:
                donor_db[donor_id].sort() ### should now be consistent between different types of covariate analyses in terms of ranking
                #print donor_db[donor_id]
                uniqueDonorSample = donor_db[donor_id][0][1]
                unique_donors.append(uniqueDonorSample)
                allUniqueDonors[uniqueDonorSample]=None
            covariate_samples[covariateType] = unique_donors

    ### Check to see which covariates completely conflict
    conflictingCovariates={}
    for variable in allCovariateSamples:
        samples=[]
        if uniqueDonors:
            for sample in allCovariateSamples[variable]:
                if sample in allUniqueDonors: samples.append(sample)
        else:
            samples = allCovariateSamples[variable]
        samples.sort()
        try: conflictingCovariates[tuple(samples)].append(variable)
        except Exception: conflictingCovariates[tuple(samples)] = [variable]
    conflicting=[]
    for samples in conflictingCovariates:
        if len(conflictingCovariates[samples])>1:
            if variable not in conflicting:
                if 'Run' not in variable and 'uid' not in variable and 'Sample' not in variable:
                    if 'other conditions during reprogramming:N/A' not in variable:
                        conflicting.append(variable)
    if len(conflicting)>0:
        print 'There were conflicting covariates present:',conflicting

    covariates_to_consider=[]
    for covariateType in covariate_samples:
        if len(covariate_samples[covariateType])>1: ### Thus at least two samples to compare
            covariates_to_consider.append(covariateType)
    comps_db={}
    for covariateType in covariates_to_consider:
        covariatePairs=[]
        for covariateType2 in covariates_to_consider:
            if covariateType!=covariateType2:
                covariatePairs = [covariateType,covariateType2]
                covariatePairs.sort()
                comps_db[tuple(covariatePairs)]=None
    
    groups_db={}
    for covariateType in covariates_to_consider:
        groups_db[covariateType] = covariate_samples[covariateType]
        ### print out the associated unique donor samples
        for i in covariate_samples[covariateType]: print covariateType,i
    #sys.exit()
    return sample_metadata,groups_db,comps_db
    
def performDifferentialExpressionAnalysis(species,platform,input_file,sample_metadata,groups_db,comps_db,CovariateQuery,uniqueDonors):
    ### filter the expression file for the samples of interest and immediately calculate comparison statistics
    firstLine = True
    group_index_db={}
    pval_summary_db={}
    group_avg_exp_db={}
    export_object_db={}
    rootdir=export.findParentDir(input_file)
    #print rootdir;sys.exit()
    
    try:
        gene_to_symbol,system_code = getAnnotations(species,platform)
        from import_scripts import OBO_import
        symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
    except ZeroDivisionError: gene_to_symbol={}; system_code=''
    
    for groups in comps_db:
        group1, group2 = groups
        pval_summary_db[groups] = {} ### setup this data structure for later
        filename='ExpressionProfiles/'+string.join(groups,'_vs_')+'.txt'
        eo = export.ExportFile(rootdir+filename) ### create and store an export object for the comparison (for raw expression)
        export_object_db[groups] = eo
        
    for line in open(input_file,'rU').xreadlines():
        if '.bed' in line:
            line = string.replace(line,'.bed','')
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstLine:
            header = values
            for group in groups_db:
                samplesToEvaluate = groups_db[group]
                try: sample_index_list = map(lambda x: values.index(x), samplesToEvaluate)
                except Exception: ### For datasets with mising samples (removed due to other QC issues)
                    sample_index_list=[]
                    filteredSamples=[]
                    for x in samplesToEvaluate:
                        try:
                            sample_index_list.append(values.index(x))
                            filteredSamples.append(x)
                        except Exception: pass
                    samplesToEvaluate = filteredSamples
                groups_db[group] = samplesToEvaluate  
                #print group, sample_index_list
                group_index_db[group] = sample_index_list
                group_avg_exp_db[group] = {}
                
            ### Write out headers for grouped expression values
            for (group1,group2) in comps_db:
                eo = export_object_db[(group1,group2)]
                g1_headers = groups_db[group1]
                g2_headers = groups_db[group2]
                g1_headers = map(lambda x: group1+':'+x,g1_headers)
                g2_headers = map(lambda x: group2+':'+x,g2_headers)
                eo.write(string.join(['GeneID']+g1_headers+g2_headers,'\t')+'\n')
            firstLine = False
        else:
            geneID = values[0]
            if ',' in geneID:
                geneID = string.replace(geneID,',','_')
            if 'ENS' in geneID and '.' in geneID and ':' not in geneID:
                geneID = string.split(geneID,'.')[0] ### for cufflinks
            elif platform == 'RNASeq':
                try: geneID = symbol_to_gene[geneID][0]
                except Exception: pass
            group_expression_values={}
            original_group={}
            for group in group_index_db:
                sample_index_list = group_index_db[group]
                if platform != 'exon':
                    filtered_values = map(lambda x: float(values[x]), sample_index_list) ### simple and fast way to reorganize the samples
                else: ### for splice-event comparisons
                    if len(header) != len(values):
                        diff = len(header)-len(values)
                        values+=diff*['']
                    initial_filtered=[] ### the blanks can cause problems here so we loop through each entry and catch exceptions
                    unfiltered=[]
                    for x in sample_index_list:
                        initial_filtered.append(values[x])
                    filtered_values=[]
                    for x in initial_filtered:
                        if x != '':
                            filtered_values.append(float(x))
                        unfiltered.append(x)
                    #if geneID == 'ENSG00000105321:E3.2-E4.2 ENSG00000105321:E2.3-E4.2' and 'inner cell mass' in group:
                    #print filtered_values;sys.exit()
                if platform == 'exon':
                    original_group[group]=unfiltered
                else:
                    original_group[group]=filtered_values
                if platform == 'RNASeq':# or platform == 'miRSeq':
                    filtered_values = map(lambda x: math.log(x+1,2),filtered_values) ### increment and log2 adjusted
                group_expression_values[group] = filtered_values
            for groups in comps_db:
                group1,group2 = groups
                data_list1 = group_expression_values[group1]
                data_list2 = group_expression_values[group2]
                if len(data_list1)>1 and len(data_list2)>1: ### For splicing data
                    p = statistics.runComparisonStatistic(data_list1,data_list2,probability_statistic) ### this is going to just be a oneway anova first
                    avg1 = statistics.avg(data_list1)
                    avg2 = statistics.avg(data_list2)
                    log_fold = avg1-avg2
                    if platform == 'RNASeq':# or platform == 'miRSeq':
                        max_avg = math.pow(2,max([avg1,avg2]))-1
                    else: max_avg = 10000
                    #if platform == 'miRSeq': max_avg = 10000
                    valid = True
                    if max_avg<minRPKM:
                        log_fold = 'Insufficient Expression'
                    gs = statistics.GroupStats(log_fold,None,p)
                    gs.setAdditionalStats(data_list1,data_list2) ### Assuming equal variance
                    pval_db = pval_summary_db[groups] ### for calculated adjusted statistics
                    pval_db[geneID] = gs ### store the statistics here
                    proceed = True
                    if len(restricted_gene_denominator)>0:
                        if geneID not in restricted_gene_denominator:
                            proceed = False
                    if uniqueDonors == False and proceed:
                        ### store a new instance
                        gsg = statistics.GroupStats(log_fold,None,p)
                        gsg.setAdditionalStats(data_list1,data_list2) ### Assuming equal variance
                        global_adjp_db[CovariateQuery,groups,geneID] = gsg ### for global adjustment across comparisons
                    #if 'lentivirus' in group1 and geneID == 'ENSG00000184470':
                    #print groups,log_fold, p, avg1,avg2, max_avg, original_group[group1],original_group[group2]

                    if geneID == 'ENSG00000140416:I1.2-E6.4 ENSG00000140416:E3.6-E6.4' and 'male' in groups and 'female' in groups:
                        print groups,log_fold, p, avg1,avg2, max_avg, original_group[group1],original_group[group2],data_list1,data_list2
                    group_avg_exp_db[group1][geneID] = avg1 ### store the group expression values
                    group_avg_exp_db[group2][geneID] = avg2 ### store the group expression values
                    #if geneID == 'ENSG00000213973': print log_fold
                    if 'Insufficient Expression2' != log_fold:
                        if geneID == 'hsa-mir-512-1_hsa-miR-512-3p' and 'CD34+ cells' in groups and 'mononuclear' in groups:
                            ls1 = map(str,original_group[group1])
                            ls2 = map(str,original_group[group2])
                            #print groups, log_fold, logfold_threshold, p, ls1, ls2, avg1, avg2
                            #print data_list1, data_list2
                            if abs(log_fold)>logfold_threshold and p<pval_threshold:
                                print 'yes'
                            else: print 'no'
                            #sys.exit()
                            pass
                        #if abs(log_fold)>logfold_threshold:
                        eo = export_object_db[groups]
                        ls1 = map(str,original_group[group1])
                        ls2 = map(str,original_group[group2])
                        eo.write(string.join([geneID]+ls1+ls2,'\t')+'\n')
                            
    for groups in export_object_db:
        export_object_db[groups].close()
    
    ### Calculate adjusted p-values for all pairwise comparisons
    for groups in pval_summary_db:
        group1,group2 = groups
        if uniqueDonors:
            filename=CovariateQuery+'/GE.'+string.join(groups,'_vs_')+'-UniqueDonors.txt'
        else:
            filename=CovariateQuery+'/GE.'+string.join(groups,'_vs_')+'.txt'
        eo = export.ExportFile(rootdir+'/'+filename)
        do = export.ExportFile(rootdir+'/Downregulated/'+filename)
        uo = export.ExportFile(rootdir+'/Upregulated/'+filename)
        so = export.ExportFile(rootdir+'PValues/'+CovariateQuery+'-'+string.join(groups,'_vs_')+'.txt')
        header = 'GeneID\tSystemCode\tLogFold\trawp\tadjp\tSymbol\tavg-%s\tavg-%s\n' % (group1,group2)
        eo.write(header)
        do.write(header)
        uo.write(header)
        so.write('Gene\tPval\n')
        pval_db = pval_summary_db[groups]
        if 'moderated' in probability_statistic:
            try: statistics.moderateTestStats(pval_db,probability_statistic) ### Moderates the original reported test p-value prior to adjusting
            except Exception: print 'Moderated test failed... using student t-test instead'
        statistics.adjustPermuteStats(pval_db) ### sets the adjusted p-values for objects
        for geneID in pval_db:
            gs = pval_db[geneID]
            group1_avg = str(group_avg_exp_db[group1][geneID])
            group2_avg = str(group_avg_exp_db[group2][geneID])
            if use_adjusted_p:
                pval = float(gs.AdjP())
            else:
                pval = gs.Pval()
            if platform == 'miRSeq':
                symbol=[]
                altID = string.replace(geneID,'hsa-mir-','MIR')
                altID = string.replace(altID,'hsa-miR-','MIR')
                altID = string.replace(altID,'3p','')
                altID = string.replace(altID,'5p','')
                altID = string.upper(string.replace(altID,'hsa-let-','LET'))
                altID = string.replace(altID,'-','')
                altIDs = string.split(altID,'_')
                altIDs+=string.split(geneID,'_')
                altIDs = unique.unique(altIDs)
                for id in altIDs:
                    if id in gene_to_symbol:
                        symbol.append(gene_to_symbol[id][0])
                        symbol.append(id)
                symbol = string.join(symbol,'|')
            elif geneID in gene_to_symbol:
                symbols = unique.unique(gene_to_symbol[geneID])
                symbol = string.join(symbols,'|')
            elif 'ENS' in geneID and ':' in geneID:
                ens_gene = string.split(geneID,':')[0]
                try: symbol = gene_to_symbol[ens_gene][0]
                except Exception: symbol=''
            else:
                symbol = ''
            proceed = True
            ### Remove genes not a predetermined list (optional)
            if len(restricted_gene_denominator)>0:
                if geneID not in restricted_gene_denominator:
                    if symbol not in restricted_gene_denominator:
                        proceed = False

            if geneID == 'ENSG00000105321:E3.2-E4.2 ENSG00000105321:E2.3-E4.2' and 'fibroblast' in groups and 'inner cell mass' in groups:
                print groups, gs.LogFold(), logfold_threshold, pval, pval_threshold, gs.Pval(), gs.AdjP(), symbol#;sys.exit()
            if 'Insufficient Expression' != gs.LogFold() and proceed:
                #if geneID == 'hsa-mir-512-1_hsa-miR-512-3p' and 'CD34+ cells' in groups and 'mononuclear' in groups:
                #print groups, gs.LogFold()
                if abs(gs.LogFold())>logfold_threshold and pval<pval_threshold:
                    values = string.join([geneID,system_code,str(gs.LogFold()),str(gs.Pval()),str(gs.AdjP()),symbol,group1_avg,group2_avg],'\t')+'\n'
                    eo.write(values)
                    try: chr = gene_location_db[ug.GeneID()][0]
                    except Exception: chr = ''
                    proceed = True
                    if 'Gender' in filename:
                        if 'Y' in chr: proceed = False
                    if proceed:
                        if gs.LogFold()>0:
                            uo.write(values)
                        if gs.LogFold()<0:
                            do.write(values)
            so.write(geneID+'\t'+str(gs.Pval())+'\n')
        eo.close()
        do.close()
        uo.close()
        so.close()
        
def getAnnotations(species,platform):
    import gene_associations
    if platform == 'RNASeq' or platform == 'exon' or platform == 'miRSeq':
        gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
        system_code = 'En'
        if platform == 'miRSeq':
            from import_scripts import OBO_import
            gene_to_symbol = OBO_import.swapKeyValues(gene_to_symbol)
    if platform == 'methylation':
        gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
        gene_to_symbol = importMethylationAnnotations(species,gene_to_symbol)
        system_code = 'Ilm'
    return gene_to_symbol, system_code

def importMethylationAnnotations(species,gene_to_symbol):
    filename = 'AltDatabase/ucsc/'+species+'/illumina_genes.txt'
    from import_scripts import OBO_import
    symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
    firstLine=True
    probe_gene_db={}
    for line in open(OBO_import.filepath(filename),'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstLine:
            geneIndex = values.index('UCSC_RefGene_Name')
            locationIndex = values.index('UCSC_RefGene_Group')
            firstLine = False
        else:
            probeID = values[0]
            try: genes = string.split(values[geneIndex],';')
            except Exception: genes=[]
            try:
                locations = unique.unique(string.split(values[locationIndex],';'))
                locations = string.join(locations,';')
            except Exception:
                locations = ''
            for symbol in genes:
                if len(symbol)>0:
                    if symbol in symbol_to_gene:
                        for geneID in symbol_to_gene[symbol]:
                            try: probe_gene_db[probeID].append(geneID)
                            except Exception: probe_gene_db[probeID] = [geneID]
                        probe_gene_db[probeID].append(symbol)
                        probe_gene_db[probeID].append(locations)
    return probe_gene_db
   
def getDatasetSamples(expression_file,sample_metadata,cellLines):
    ### Required as samples may exist in the metadata but were excluded due to QC
    for line in open(expression_file,'rU').xreadlines():
        if '.bed' in line:
            line = string.replace(line,'.bed','')
        data = line.rstrip()
        headers = string.split(data,'\t')[1:]
        break
    supported_cellLines={}
    for s in headers:
        if s in sample_metadata:
            metadata = sample_metadata[s]
            if len(metadata.CellLine())>0:
                try:
                    if s == cellLines[metadata.CellLine()]: ### Make sure the correct samples is being matched
                        supported_cellLines[metadata.CellLine()]=None
                except Exception: pass
    return supported_cellLines

def importExpressionData(species,platform,expression_file,cell_line_db,common_lines):
    ### Imports the omics data, filters/orders samples, transforms (if necessary), keys by a common geneID
    filtered_exp_data={}
    try: gene_to_symbol,system_code = getAnnotations(species,platform)
    except ZeroDivisionError: gene_to_symbol={}; system_code=''
    
    firstLine=True
    for line in open(expression_file,'rU').xreadlines():
        if '.bed' in line:
            line = string.replace(line,'.bed','')
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstLine:
            samplesToEvaluate = map(lambda x: cell_line_db[x], common_lines)
            sample_index_list = map(lambda x: values.index(x), samplesToEvaluate)
            header = values
            #print len(samplesToEvaluate),platform, samplesToEvaluate
            firstLine = False
        else:
            try: filtered_values = map(lambda x: float(values[x]), sample_index_list) ### simple and fast way to reorganize the samples
            except Exception: ### for splice-event comparisons
                if len(header) != len(values):
                    diff = len(header)-len(values)
                    values+=diff*['']
                initial_filtered=[] ### the blanks can cause problems here so we loop through each entry and catch exceptions
                initial_filtered = map(lambda x: values[x], sample_index_list)
                filtered_values=[]
                for x in initial_filtered:
                    if x != '': filtered_values.append(float(x))
            if platform == 'RNASeq':# or platform == 'miRSeq':
                filtered_values = map(lambda x: math.log(x+1,2),filtered_values) ### increment and log2 adjusted
            uid = values[0]
            geneIDs = []
            if platform == 'miRSeq':
                altID = string.replace(uid,'hsa-mir-','MIR')
                altID = string.replace(altID,'hsa-miR-','MIR')
                altID = string.replace(altID,'3p','')
                altID = string.replace(altID,'5p','')
                altID = string.upper(string.replace(altID,'hsa-let-','MIRLET'))
                altID = string.replace(altID,'-','')
                altIDs = string.split(altID,'_')
                altIDs+=string.split(uid,'_')
                altIDs = unique.unique(altIDs)
                for id in altIDs:
                    if id in gene_to_symbol:
                        geneIDs.append((gene_to_symbol[id][0],uid))
            original_uid = uid
            if platform == 'methylation' and ':' in uid:
                uid=string.split(uid,':')[1]
            if 'ENS' in uid and '.' in uid and ':' not in uid:
                uid = string.split(uid,'.')[0] ### for cufflinks
            if 'ENS' in uid:
                if uid in gene_to_symbol:
                    symbol = gene_to_symbol[uid][0]
                else:
                    symbol = ''
                geneIDs = [(uid,symbol)]
            elif uid in gene_to_symbol:
                if 'uid'== 'cg21028156':
                    print gene_to_symbol[uid]
                for g in gene_to_symbol[uid]:
                    if 'ENS' in g:
                        if platform == 'methylation' and ':' in original_uid:
                            uid = original_uid
                        geneIDs.append((g,uid))
            for (geneID,uid) in geneIDs:
                try: filtered_exp_data[geneID].append((uid,filtered_values))
                except Exception: filtered_exp_data[geneID] = [(uid,filtered_values)]
    print len(filtered_exp_data)
    return filtered_exp_data, samplesToEvaluate

def combineAndCompareMatrices(input_file,filtered_exp_data1,filtered_exp_data2,platform1,platform2,samplesToEvaluate):
    ### Get the matching genes and identify anti-correlated rows (same sample order) to export to a merged file
    rootdir=export.findParentDir(input_file)
    compared_uid_pairs={}
    count=0
    import warnings
    from  scipy import stats
    p1_samples = map(lambda x: platform1+':'+x, samplesToEvaluate)
    p2_samples = map(lambda x: platform2+':'+x, samplesToEvaluate)
    exportFile = rootdir+'/MergedOmicsTables/'+platform1+'-'+platform2+'.txt'
    if uniqueDonors:
        exportFile = rootdir+'/MergedOmicsTables/'+platform1+'-'+platform2+'-UniqueDonors.txt'
    co = export.ExportFile(exportFile)
    co.write(string.join(['GeneID',platform1+'-UID','Pearson-rho']+p1_samples+[platform2+'-UID']+p2_samples,'\t')+'\n')
    correlated_geneIDs={}
    for geneID in filtered_exp_data1:
        if geneID in filtered_exp_data2:
            rows1 = filtered_exp_data1[geneID]
            rows2 = filtered_exp_data2[geneID]
            for (uid1,row1) in rows1:
                """
                if platform1 == 'RNASeq' and platform2 == 'methylation':
                    try: row1 = map(lambda x: math.pow(2,x)-1,row1)
                    except Exception: print uid1,row1;sys.exit()
                """
                for (uid2,row2) in rows2:
                    try: null=compared_uid_pairs[(uid1,uid2)] ### already compared
                    except Exception:
                        with warnings.catch_warnings():
                            warnings.filterwarnings("ignore") ### hides import warnings
                            try: rho,p = stats.pearsonr(row1,row2)
                            except Exception: print 'The rows are not of equal length, likely due to missing values in that row:',uid1,uid2;sys.exit()
                            compared_uid_pairs[(uid1,uid2)]=None
                            if rho < -0.5:
                                values = [geneID,uid1,rho]+row1+[uid2]+row2
                                values = string.join(map(str, values),'\t')+'\n'
                                correlated_geneIDs[geneID]=None
                                co.write(values)
                                count+=1
    co.close()
    print 'Writing out %d entries to %s:' % (count,exportFile)
    correlated_geneIDs_ls=[]
    for i in correlated_geneIDs:
        correlated_geneIDs_ls.append(i)
    print len(correlated_geneIDs_ls)
    

def getFiles(sub_dir,directories=True):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    for entry in dir_list:
        if directories:
            if '.' not in entry: dir_list2.append(entry)
        else:
            if '.' in entry: dir_list2.append(entry)
    return dir_list2

class GeneData:
    def __init__(self,geneID, systemCode, logFold, rawp, adjp, symbol, avg1, avg2):
        self.geneID = geneID; self.systemCode = systemCode; self.logFold = logFold; self.rawp = rawp; self.adjp = adjp
        self.symbol = symbol; self.avg1 = avg1; self.avg2 = avg2
    def GeneID(self): return self.geneID
    def LogFold(self): return self.logFold
    def Rawp(self): return self.rawp
    def Adjp(self): return self.adjp
    def Symbol(self): return self.symbol
    def Avg1(self): return self.avg1
    def Avg2(self): return self.avg2
    def SystemCode(self): return self.systemCode    
    def __repr__(self): return self.GeneID()
    
def importResultsSummary(filepath,comparison,gene_associations_db):
    #'GeneID\tSystemCode\tLogFold\trawp\tadjp\tSymbol\tavg-%s\tavg-%s\n
    firstLine=True
    for line in open(filepath,'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstLine:
            firstLine=False
        else:
            geneID, systemCode, logFold, rawp, adjp, symbol, avg1, avg2 = values
            gd = GeneData(geneID, systemCode, logFold, rawp, adjp, symbol, avg1, avg2)
            if comparison in gene_associations_db:
                gene_db = gene_associations_db[comparison]
                gene_db[geneID]=gd
            else:
                gene_db = {}
                gene_db[geneID]=gd
                gene_associations_db[comparison]=gene_db
    return gene_associations_db

def compareGOEliteEnrichmentProfiles(expressionDir,eliteDir):
    up_elite_miRNAs={}
    down_elite_miRNAs={}
    folders = getFiles(expressionDir)
    for folder in folders:
        subdirs = getFiles(expressionDir+'/'+folder)
        for sub_dir in subdirs:
            subdir = expressionDir+'/'+folder + '/'+ sub_dir
            elite_dirs = getFiles(subdir) ### Are there any files to analyze?
            if 'GO-Elite_results' in elite_dirs:
                elitedir = subdir + '/GO-Elite_results/pruned-results_z-score_elite.txt'
                if 'Down' in folder:
                    try: down_elite_miRNAs = getMIRAssociations(elitedir,sub_dir,down_elite_miRNAs)
                    except Exception: pass
                else:
                    try: up_elite_miRNAs = getMIRAssociations(elitedir,sub_dir,up_elite_miRNAs)
                    except Exception: pass

    if '.txt' not in eliteDir:
        if 'CombinedResults' not in eliteDir:
            eliteDir += '/CombinedResults/allTopGenes.txt'
        else: eliteDir += '/allTopGenes.txt'
    firstLine=True
    for line in open(eliteDir,'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstLine:
            firstLine=False
        else:
            comparison, gene, symbol, up_rawp, ng_adjp, up_logfold,ng_logofold, ng_avg1, ng_avg2 = values[:9]
            comparison = string.replace(comparison,'.txt','')
            miRNAs = string.split(gene,'_')+string.split(symbol,'|')
            log_fold = float(up_logfold)
            if log_fold>0:
                if comparison in down_elite_miRNAs:
                    for miRNA in miRNAs:
                        miRNA = string.lower(miRNA)
                        if miRNA in down_elite_miRNAs[comparison]:
                            if 'lab' not in comparison and 'CD34+ cells_vs_mononuclear' not in comparison:
                                print miRNA, comparison, 'down'
            else:
                if comparison in up_elite_miRNAs:
                    for miRNA in miRNAs:
                        miRNA = string.lower(miRNA)
                        if miRNA in up_elite_miRNAs[comparison]:
                            if 'lab' not in comparison and 'CD34+ cells_vs_mononuclear' not in comparison:
                                print miRNA, comparison, 'up'
                
def importRestrictedSetOfGenesToQuery(filepath):
    ### Applied to predetermined expressed genes matching some criterion (e.g., FPKM > 5 and 20% expression in EBs)
    restricted_gene_denominator_db={}
    firstLine=True
    for line in open(filepath,'rU').xreadlines():
        data = line.rstrip()
        gene = string.split(data,'\t')[0]
        if firstLine:
            firstLine=False
        else:
            restricted_gene_denominator_db[gene]=[]
    return restricted_gene_denominator_db

def getMIRAssociations(filepath,diffstate_comparison,elite_miRNAs):
        
    firstLine=True
    for line in open(filepath,'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstLine:
            firstLine=False
        else:
            try:
                #if 'Combined' in values[0]: ### Restrict comparison to these
                regulated_geneset_name = values[0] # GE.group1_vs_group2-microRNATargets.txt
                regulated_geneset_name = string.split(regulated_geneset_name,'-')[0]
                #regulated_geneset_name = string.replace(regulated_geneset_name,'-UniqueDonors','')
                #regulated_geneset_name = string.replace(regulated_geneset_name,'-Combined','')
                miRNA = string.lower(values[2])
                miRNA = string.replace(miRNA,'*','')
                miRNA2 = string.replace(miRNA,'hsa-mir-','MIR')
                miRNA2 = string.replace(miRNA2,'hsa-let-','LET')
                miRNA3 = string.replace(miRNA2,'-5p','')
                miRNA3 = string.replace(miRNA3,'-3p','')
                miRNAs = [miRNA,miRNA2,miRNA3]
                for miRNA in miRNAs:
                    try: elite_miRNAs[diffstate_comparison+':'+regulated_geneset_name].append(string.lower(miRNA))
                    except Exception: elite_miRNAs[diffstate_comparison+':'+regulated_geneset_name] = [string.lower(miRNA)]   
            except Exception:
                pass

    return elite_miRNAs

def runGOEliteAnalysis(species,resultsDirectory):
    mod = 'Ensembl'
    pathway_permutations = 'FisherExactTest'
    filter_method = 'z-score'
    z_threshold = 1.96
    p_val_threshold = 0.05
    change_threshold = 2
    resources_to_analyze = ['microRNATargets','pictar','miRanda','mirbase','RNAhybrid','TargetScan','microRNATargets_All']
    returnPathways = 'no'
    root = None
    import GO_Elite
    print '\nBeginning to run GO-Elite analysis on all results'
    
    folders = getFiles(resultsDirectory)
    for folder in folders:
        subdirs = getFiles(resultsDirectory+'/'+folder)
        for subdir in subdirs:
            subdir = resultsDirectory+'/'+folder + '/'+ subdir
            file_dirs = subdir,None,subdir
            input_files = getFiles(subdir,directories=False) ### Are there any files to analyze?
            if len(input_files)>0:
                variables = species,mod,pathway_permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,returnPathways,file_dirs,root
                try: GO_Elite.remoteAnalysis(variables,'non-UI')
                except Exception: 'GO-Elite failed for:',subdir

def identifyCommonGenes(resultsDirectory):
    """ Compares results from parallel statistical analyses for unique and non-unique genetic donor workflows """
    uniqueDonorGenes = {}
    nonUniqueDonorGenes={}
    folders = getFiles(resultsDirectory)
    for folder in folders:
        files = getFiles(resultsDirectory+'/'+folder,directories=False)
        for file in files:
            if '.txt' in file and 'GE.'== file[:3]:
                filepath = resultsDirectory+'/'+folder+'/'+file
                comparison = folder+':'+string.replace(file,'-UniqueDonors.txt','.txt')
                if 'UniqueDonors.txt' in filepath:
                    uniqueDonorGenes = importResultsSummary(filepath,comparison,uniqueDonorGenes)
                else:
                    nonUniqueDonorGenes = importResultsSummary(filepath,comparison,nonUniqueDonorGenes)

    #nonUniqueDonorGenes = uniqueDonorGenes
    from build_scripts import EnsemblImport
    try: gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,platform,'key_by_array')
    except Exception: gene_location_db={}
    
    includeGlobalAdjustedPvals = False
    if len(global_adjp_db)>0: ### When all comparisons are run together
        #global_adjp_db[CovariateQuery,uniqueDonors,groups,geneID] = gs
        if 'moderated' in probability_statistic:
            try: statistics.moderateTestStats(global_adjp_db,probability_statistic) ### Moderates the original reported test p-value prior to adjusting
            except Exception: print 'Moderated test failed... using student t-test instead'
        statistics.adjustPermuteStats(global_adjp_db) ### sets the adjusted p-values for objects
        includeGlobalAdjustedPvals = True

    output_dir = resultsDirectory+'/CombinedResults/allTopGenes.txt'
    eo = export.ExportFile(output_dir)
    header = 'Comparison\tGeneID\tSymbol\tUniqueDonor-rawp\tNonUnique-adjp\tUniqueDonor-LogFold\tNonUnique-LogFold\tNonUnique-Avg1\tNonUnique-Avg2'
    if includeGlobalAdjustedPvals:
        header+='\tGlobalAdjustedP'
    eo.write(header+'\n')
    topComparisonAssociations={}
    for comparison in uniqueDonorGenes:
        if comparison in nonUniqueDonorGenes:
            CovariateQuery,groups = string.split(comparison[:-4],':')
            groups = tuple(string.split(groups[3:],'_vs_'))
            comparison_dir = string.replace(comparison,':','/')[:-4]
            do = export.ExportFile(resultsDirectory+'/Downregulated/'+comparison_dir+'-Combined.txt')
            uo = export.ExportFile(resultsDirectory+'/Upregulated/'+comparison_dir+'-Combined.txt')
            header = 'GeneID\tSy\tFoldChange\trawp\n'
            uo.write(header)
            do.write(header)
            unique_gene_db = uniqueDonorGenes[comparison]
            nonunique_gene_db = nonUniqueDonorGenes[comparison]
            for gene in unique_gene_db: ### loop through the gene dictionary
                if gene in nonunique_gene_db: ### common genes between unique and non-unique donors
                    ug = unique_gene_db[gene]
                    ng = nonunique_gene_db[gene]
                    values = [comparison,gene, ug.Symbol(),ug.Rawp(),ng.Adjp(),ug.LogFold(),ng.LogFold(),ng.Avg1(),ng.Avg2()]
                    if includeGlobalAdjustedPvals:
                        try:
                            gs = global_adjp_db[CovariateQuery,groups,gene]
                            ng_adjp = float(gs.AdjP())
                            values+=[str(ng_adjp)]
                        
                            if platform == 'miRSeq' or platform == 'exon' and use_adjusted_p == False:
                                ng_adjp = float(ug.Rawp())
                        except Exception:
                            if platform == 'miRSeq' or platform == 'exon' and use_adjusted_p == False:
                                ng_adjp = float(ug.Rawp())
                    else:
                        ng_adjp = float(ug.Rawp())
                    values = string.join(values,'\t')+'\n'
                    eo.write(values)
                    if ng_adjp<pval_threshold:
                        try: topComparisonAssociations[gene].append((float(ug.Rawp()),values))
                        except Exception: topComparisonAssociations[gene] = [(float(ug.Rawp()),values)]
                        values = [ug.GeneID(), ug.SystemCode(), ug.LogFold(), ug.Rawp()]
                        values = string.join(values,'\t')+'\n'
                        try: chr = gene_location_db[ug.GeneID()][0]
                        except Exception: chr = ''
                        proceed = True
                        if 'Gender' in comparison:
                            if 'Y' in chr: proceed = False
                        if proceed:
                            if float(ug.LogFold())>0:
                                uo.write(values)
                            else:
                                do.write(values)    
            do.close()
            uo.close()
                        
    eo.close()
    print 'Matching Unique-Donor and NonUnique Donor results written to:',output_dir
    
    ### Write out the comparison for each gene with the most significant result (best associations)
    output_dir = resultsDirectory+'/CombinedResults/eliteTopGenes.txt'
    eo = export.ExportFile(output_dir)
    eo.write('Comparison\tGeneID\tSymbol\tUniqueDonor-rawp\tNonUnique-adjp\tUniqueDonor-LogFold\tNonUnique-LogFold\tNonUnique-Avg1\tNonUnique-Avg2\n')
    for gene in topComparisonAssociations:
        topComparisonAssociations[gene].sort()
        eo.write(topComparisonAssociations[gene][0][1])
    eo.close()
    print 'The most significant comparisons for each gene reported to:',output_dir

def getCommonCellLines(cellLines1,cellLines2,exp_cellLines1,exp_cellLines2,uniqueDonor_db,uniqueDonors,donor_db):
    common_lines = list(cellLines1.viewkeys() & cellLines2.viewkeys() & exp_cellLines1.viewkeys() & exp_cellLines2.viewkeys())
    common_lines.sort()
    exclude = ['SC11-009','SC11-008'] ### selected the line with the greatest XIST for this donor
    if uniqueDonors:
        common_lines2=[]; donor_added=[]
        for donorID in donor_db:
            donor_db[donorID].sort()
            for (priority,cellLine) in donor_db[donorID]: ### Prioritize based on sample QC and donor preference
                if cellLine in common_lines and cellLine not in exclude:
                    if donorID not in donor_added and 'H9' not in cellLine:
                        common_lines2.append(cellLine)
                        donor_added.append(donorID)    
        common_lines = common_lines2
    print 'Common Approved Lines:',common_lines
    return common_lines

def downloadSynapseFile(synid,output_dir):
    import synapseclient
    import os,sys,string,shutil
    syn = synapseclient.Synapse()
    syn.login()
    matrix = syn.get(synid, downloadLocation=output_dir, ifcollision="keep.local")
    return matrix.path

def buildAdditionalMirTargetGeneSets():
    miR_association_file = 'AltDatabase/ensembl/Hs/Hs_microRNA-Ensembl.txt'
    output_dir = 'AltDatabase/EnsMart72/goelite/Hs/gene-mapp/Ensembl-'
    eo = export.ExportFile(output_dir+'microRNATargets_All.txt')
    header = 'GeneID\tSystemCode\tmiRNA\n'
    eo.write(header)
    miRNA_source_db={}
    for line in open(unique.filepath(miR_association_file),'rU').xreadlines():
        data = line.rstrip()
        miRNA,ensembl,source = string.split(data,'\t')
        output = ensembl+'\t\t'+miRNA+'\n'
        eo.write(output)
        sources = string.split(source,'|')
        for source in sources:
            try: miRNA_source_db[source].append(output)
            except KeyError:  miRNA_source_db[source]=[output]
    eo.close()
    
    for source in miRNA_source_db:
        eo = export.ExportFile(output_dir+source+'.txt')
        eo.write(header)
        for line in miRNA_source_db[source]:
            eo.write(line)
        eo.close()

def returnSynFileLocations(file1,file2,output_dir):
    if 'syn' in file1:
        try: file1 = downloadSynapseFile(file1,output_dir)
        except Exception:
            print 'Is the destination file %s already open?' % file1;sys.exit()
    if 'syn' in file2:
        try: file2 = downloadSynapseFile(file2,output_dir)
        except Exception:
            print 'Is the destination file %s already open?' % file2;sys.exit()
    return file1,file2

def synapseStore(file_dirs,root,parent_syn,executed_urls,used):
    for file in file_dirs:
        file_dir = root+'/'+file
        file = synapseclient.File(file_dir, parent=parent_syn)
        file = syn.store(file,executed=executed_urls,used=used)
        
def synapseStoreFolder(dir_path,parent_syn):
    data_folder = synapseclient.Folder(dir_path, parent=parent_syn)
    data_folder = syn.store(data_folder)
    sub_parent = data_folder.id
    return sub_parent

def synapseDirectoryUpload(expressionDir, parent_syn, executed_urls, used):
    root = string.split(expressionDir,'/')[-1]
    sub_parent = synapseStoreFolder(root,parent_syn)
    folders = getFiles(expressionDir)
    for folder in folders:
        ### Create the folder in Synapse
        dir_path_level1 = expressionDir+'/'+folder
        sub_parent1 = synapseStoreFolder(folder,sub_parent)
        f2 = getFiles(dir_path_level1)
        files = getFiles(dir_path_level1,False)
        synapseStore(files,dir_path_level1,sub_parent1,executed_urls,used)
        for folder in f2:
            dir_path_level2 = dir_path_level1+'/'+folder
            sub_parent2 = synapseStoreFolder(folder,sub_parent1)
            f3 = getFiles(dir_path_level2)
            files = getFiles(dir_path_level2,False)
            synapseStore(files,dir_path_level2,sub_parent2,executed_urls,used)
            for folder in f3:
                dir_path_level3 = dir_path_level2+'/'+folder
                sub_parent3 = synapseStoreFolder(folder,sub_parent2)
                files = getFiles(dir_path_level3,False)
                ### These are the GO-Elite result files (not folders)
                synapseStore(files,dir_path_level3,sub_parent3,executed_urls,used)

def exportGeneSetsFromCombined(filename):
    firstLine=True
    synapse_format = True
    simple_format = False
    reverse = True
    comparison_to_gene={}
    rootdir = export.findParentDir(filename)
    file = export.findFilename(filename)
    for line in open(filename,'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstLine:
            firstLine=False
        else:
            comparison, gene, symbol, up_rawp, ng_adjp, up_logfold,ng_logofold, ng_avg1, ng_avg2 = values[:9]
            comparison = string.replace(comparison,'GE.','')
            prefix = string.split(comparison,':')[0]
            state = string.split(prefix,'-')[0]
            if state=='NA': state = 'All'
            comparison = file[:-4]+'-'+string.split(comparison,':')[1]
            comparison = string.replace(comparison,'allTopGenes-','')[:-4]
            c1,c2 = string.split(comparison,'_vs_')
            c1 = re.sub('[^0-9a-zA-Z]+', '', c1)
            c2 = re.sub('[^0-9a-zA-Z]+', '', c2)
            comparison = c1+'_vs_'+c2
            #comparison = string.replace(comparison,':','-')
            log_fold = float(up_logfold)
            #print c1,c2,reverse,log_fold
            if c1>c2:
                comparison = c2+'_vs_'+c1
                log_fold = log_fold *-1
            if log_fold<0:
                if synapse_format:
                    comparison+='__down'
                else:
                    comparison = c2+'_vs_'+c1 ### reverse the regulation direction
            else:
                if synapse_format:
                    comparison+='__up'
            if synapse_format:
                if len(symbol) == 0: symbol = gene
                gene = symbol
                #comparison = string.replace(comparison,'_vs_','_')
                comparison = string.replace(comparison,'NA','NotApplicable')
                #comparison = string.replace(comparison,'MESO-5','MESO-EARLY')
                comparison = string.replace(comparison,' ','_')
                #if state == 'MESO': state = 'MESOEARLY'
            if 'ENS' in gene:
                SystemCode = 'En'
                if ':' in gene:
                    gene = string.split(gene,':')[0]
                genes = [gene]
            elif 'cg' in gene:
                SystemCode = 'En'
                genes = []
                for i in string.split(symbol,'|'):
                    if 'ENS' in i:
                        genes.append(i)
            else:
                SystemCode = 'Sy'
                genes = [gene]
            for g in genes:
                if synapse_format:
                    if simple_format:
                        try: comparison_to_gene[comparison].append([g,log_fold])
                        except Exception: comparison_to_gene[comparison] = [[g,log_fold]]
                    else:
                        try: comparison_to_gene[state+'__'+comparison].append([g,log_fold])
                        except Exception: comparison_to_gene[state+'__'+comparison] = [[g,log_fold]]
                elif simple_format:
                    try: comparison_to_gene[comparison].append([g,log_fold])
                    except Exception: comparison_to_gene[comparison] = [[g,log_fold]]
                else:
                    try: comparison_to_gene[state+'-'+comparison].append([g,log_fold])
                    except Exception: comparison_to_gene[state+'-'+comparison] = [[g,log_fold]]
                
    aro = export.ExportFile(rootdir+'/Regulated/combined.txt')
    aro.write('Gene\tLogFold\tComparison\n')
    for comparison in comparison_to_gene:
        ro = export.ExportFile(rootdir+'/Regulated/'+comparison+'.txt')
        ro.write('Gene\tSystemCode\n')
        for (gene,logfold) in comparison_to_gene[comparison]:
            ro.write(gene+'\t'+SystemCode+'\n')
            aro.write(gene+'\t'+str(logfold)+'\t'+string.replace(comparison,'.txt','')+'\n')
        ro.close()
    aro.close()

if __name__ == '__main__':
    ################  Comand-line arguments ################
    #buildAdditionalMirTargetGeneSets();sys.exit()
    filename = '/Users/saljh8/Desktop/PCBC_MetaData_Comparisons/eXpress/CombinedResults/allTopGenes.txt' #DiffStateComps Reprogramming
    #exportGeneSetsFromCombined(filename);sys.exit()
    
    platform='RNASeq'
    species='Hs'
    probability_statistic = 'moderated t-test'
    #probability_statistic = 'unpaired t-test'
    minRPKM=-1000
    logfold_threshold=math.log(1.3,2)
    pval_threshold=0.05
    use_adjusted_p = False
    expression_files=[]
    platforms=[]
    metadata_files=[]
    gender_restricted = None
    runGOElite=False
    compareEnrichmentProfiles=False
    restrictCovariateTerm=None
    compDiffState=None
    runAgain=False
    output_dir=None
    include_only=''
    diffStateQuery = 'NA'
    used=[]
    executed_urls=[]
    restricted_gene_denominator={}
    global_adjp_db={}
    covariate_set = ['Cell Type of Origin Combined', 'Cell Type of Origin', 'Cell Line Type','Reprogramming Vector Type']
    covariate_set+= ['Reprogramming Gene Combination','Gender','originating lab','Donor Life Stage','Culture_Conditions','C4 Karyotype Result','Small Molecules','Other Conditions During Reprogramming','XIST Level']
    #covariate_set = ['Donor Life Stage','Other Conditions During Reprogramming','C4 Karyotype Result','Small Molecules','Other Conditions During Reprogramming']
    #covariate_set = ['Cell Line Type']
    diffState_set = ['SC','ECTO','DE','MESO-5','EB','MESO-15','MESO-30']
    import getopt
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print 'Supply the argument --i location'
        ### python metaDataAnalysis.py --i /Volumes/SEQ-DATA\ 1/PCBC/RNASeq/July2014/MetaData/RNASeq_MetaData_July2014.txt --key "CellType" --value "Cell Type of Origin"
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['m=','i=','d=','c=','u=','p=','s=','g=','e=','ce=','rc=','cd=','o=','md=','in=','target=','parent=','urls=','used='])
        for opt, arg in options:
            if opt == '--m': metadata_files.append(arg)
            if opt == '--o':
                if output_dir==None:
                    output_dir = arg
                else:
                    output_dir = [output_dir,arg]
            if opt == '--i': expression_files.append(arg)
            if opt == '--e': runGOElite=True
            if opt == '--ce': compareEnrichmentProfiles = True
            if opt == '--d': diffStateQuery=arg
            if opt == '--c': CovariateQuery=arg
            if opt == '--p': platforms.append(arg)
            if opt == '--g': gender_restricted=arg
            if opt == '--s': species=arg
            if opt == '--rc': restrictCovariateTerm=arg
            if opt == '--cd': compDiffState=arg
            if opt == '--md': mirDataDir=arg
            if opt == '--in': include_only=arg
            if opt == '--target': target_dir=arg
            if opt == '--parent': parent_syn=arg
            if opt == '--urls': executed_urls.append(arg) ### options are: all, junction, exon, reference
            if opt == '--used': used.append(arg)  
            if opt == '--u':
                if string.lower(arg) == 'yes' or string.lower(arg) == 'true':
                    uniqueDonors=True
                    use_adjusted_p = False
                else:
                    uniqueDonors = False
                    if string.lower(arg) == 'both':
                        runAgain = True
    if len(used)>0:
        ###Upload existing results folder to Synapse
        import synapseclient
        import os,sys,string,shutil,getopt
        syn = synapseclient.Synapse()
        syn.login()
        synapseDirectoryUpload(target_dir, parent_syn, executed_urls, used)
        sys.exit()
    elif compareEnrichmentProfiles:
        #print expression_files
        compareGOEliteEnrichmentProfiles(expression_files[0],expression_files[1])
    elif runGOElite:
        runGOEliteAnalysis(species,expression_files[0])
    elif (len(expression_files)==1 and '.txt' in expression_files[0]) or (len(expression_files)==1 and 'syn' in expression_files[0]):
        ### Perform a covariate based analysis on the lone input expression file
        metadata_file = metadata_files[0]
        if 'syn' in metadata_file:
            try: metadata_file = downloadSynapseFile(metadata_file,output_dir)
            except Exception:
                print 'Is the destination file %s already open?' % metadata_file;sys.exit()
        expression_file = expression_files[0]
        if 'syn' in expression_file:
            try:expression_file = downloadSynapseFile(expression_file,output_dir)
            except Exception:
                print 'Is the destination file %s already open?' % expression_file;sys.exit()
        if 'syn' in include_only:
            try:include_only = downloadSynapseFile(include_only,output_dir)
            except Exception:
                print 'Is the destination file %s already open?' % include_only;sys.exit()

        if len(platforms)>0: platform = platforms[0]
        if platform == 'exon' or platform == 'methylation':
            logfold_threshold=math.log(1.1892,2) ### equivalent to a 0.25 dPSI or 0.25 beta differences
        if platform == 'exon':
            logfold_threshold=math.log(1,2)
            use_adjusted_p = False ### Too many drop-outs with lowdepth seq that adjusting will inherently exclude any significant changes
        if platform == 'methylation':
            use_adjusted_p = True
        if platform == 'miRSeq':
            use_adjusted_p = False
            logfold_threshold=math.log(1,2)
        if CovariateQuery != 'all': covariate_set = [CovariateQuery]
        if diffStateQuery != 'all': diffState_set = [diffStateQuery]
        print 'Filtering on adjusted p-value:',use_adjusted_p

        from build_scripts import EnsemblImport
        try: gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,platform,'key_by_array')
        except Exception: gene_location_db={}
        
        if len(include_only)>0:
            restricted_gene_denominator = importRestrictedSetOfGenesToQuery(include_only)
        
        for CovariateQuery in covariate_set:
          for diffStateQuery in diffState_set:
            print 'Analyzing the covariate:',CovariateQuery, 'and diffState:',diffStateQuery, 'unique donor analysis:',uniqueDonors
            if 'XIST' in CovariateQuery: gender_restricted='female'
            
            genderRestricted = gender_restricted
            try:
                sample_metadata,groups_db,comps_db = prepareComparisonData(metadata_file,diffStateQuery,CovariateQuery,uniqueDonors,genderRestricted,platform=platform,compDiffState=compDiffState,restrictCovariateTerm=restrictCovariateTerm)
                performDifferentialExpressionAnalysis(species,platform,expression_file,sample_metadata,groups_db,comps_db,diffStateQuery+'-'+CovariateQuery,uniqueDonors)
            except Exception:
                print traceback.format_exc()
            if runAgain:
                uniqueDonors=True
                use_adjusted_p = False
                print 'Analyzing the covariate:',CovariateQuery, 'and diffState:',diffStateQuery, 'unique donor analysis:',uniqueDonors
                try:
                    sample_metadata,groups_db,comps_db = prepareComparisonData(metadata_file,diffStateQuery,CovariateQuery,uniqueDonors,genderRestricted,platform=platform,compDiffState=compDiffState,restrictCovariateTerm=restrictCovariateTerm)
                    performDifferentialExpressionAnalysis(species,platform,expression_file,sample_metadata,groups_db,comps_db,diffStateQuery+'-'+CovariateQuery,uniqueDonors)
                except Exception: pass
                uniqueDonors=False; use_adjusted_p = True
                if platform == 'miRSeq' or platform == 'exon': use_adjusted_p = False
            
        if runAgain:
            root_exp_dir=export.findParentDir(expression_file)
            identifyCommonGenes(root_exp_dir)
            runGOEliteAnalysis(species,root_exp_dir)
            try: compareGOEliteEnrichmentProfiles(root_exp_dir,mirDataDir)
            except Exception: pass
    elif '.txt' not in expression_files[0] and 'syn' not in expression_files[0]:
        ### Compare unique and non-unique results to get overlapps (single high-confidence result file)
        ### The parent directory that contain results from the above analysis serve as input
        resultsDirectory = expression_files[0]
        identifyCommonGenes(resultsDirectory)
    else:
        ### Perform a correlation analysis between two omics technologies
        
        expression_file1,expression_file2 = expression_files
        platform1,platform2 = platforms
        metadata_file1,metadata_file2 = metadata_files
        output_dir1, output_dir2 = output_dir

        print expression_files, output_dir1
        print metadata_files, output_dir2
        expression_file1,metadata_file1 = returnSynFileLocations(expression_file1,metadata_file1,output_dir1)
        metadata_file2,expression_file2 = returnSynFileLocations(metadata_file2,expression_file2,output_dir2)
        print expression_file1, expression_file2
        print metadata_file1, metadata_file2
        
        cellLines1,sample_metadata1,uniqueDonor_db,donor_db = prepareComparisonData(metadata_file1,diffStateQuery,None,uniqueDonors,gender_restricted,platform=platform1)
        cellLines2,sample_metadata2,uniqueDonor_db,donor_db = prepareComparisonData(metadata_file2,diffStateQuery,None,uniqueDonors,gender_restricted,platform=platform2)
        exp_cellLines1 = getDatasetSamples(expression_file1,sample_metadata1,cellLines1)
        exp_cellLines2 = getDatasetSamples(expression_file2,sample_metadata2,cellLines2)
        common_lines = getCommonCellLines(cellLines1,cellLines2,exp_cellLines1,exp_cellLines2,uniqueDonor_db,uniqueDonors,donor_db)
        filtered_exp_data1, samplesToEvaluate = importExpressionData(species,platform1,expression_file1,cellLines1,common_lines)
        filtered_exp_data2, samplesToEvaluate = importExpressionData(species,platform2,expression_file2,cellLines2,common_lines)
        combineAndCompareMatrices(expression_file1,filtered_exp_data1,filtered_exp_data2,platform1,platform2,samplesToEvaluate) ### export results to two different files and a combined file
