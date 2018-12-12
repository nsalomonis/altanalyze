import sys,string,os,re
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
#import altanalyze.unique as unique

import unique
import export
import math
from stats_scripts import statistics
import traceback
import collections
import junctionGraph

"""
    This script takes existing formatted metadata (C4 PCBC approved fields) and filters them to determine unique and non-unique
    donors for a specific covariate and derives comparison and group relationships to extract from an existing expression file
"""

def BatchCheck(sample_id,nonpreferential_batchs,preferential_samples,platform):
    priority=0
    for batch in nonpreferential_batchs:
        if batch in sample_id:
            priority=1
    if platform == 'RNASeq' or platform == 'PSI':
        if sample_id not in preferential_samples: priority = 1
        elif sample_id in preferential_samples: priority = 0
    return priority
        
class MetaDataQuery:
    def __init__(self, field_name, field_type, field_values):
        self.field_name = field_name
        self.field_type = field_type
        self.field_values = field_values
    def FieldName(self): return self.field_name
    def FieldType(self): return self.field_type
    def FieldValues(self): ### can be a list of values
        if self.field_type == 'NumFilter':
            cutoff = float(self.field_values[1:])
            if self.field_values[0]=='>':
                direction = '>'
            else:
                direction = '<'
            return cutoff, direction
        else:
            ls_results = string.split(self.field_values,',')
            return ls_results
    def __repr__(self): return self.FieldName(), self.FieldName(), self.FieldType(),FieldValues()
    
def importMetaDataDescriptions(field_description_file):
    ###Import the looks ups to indicate which fields to consider for metadata filtering
    metadata_filters={}
    for line in open(field_description_file,'rU').xreadlines():
        data = line.rstrip()
        data = string.replace(data,'"','')
        values = string.split(data,'\t')
        if len(values)>1: ### Hence this field should be considered for filtering
            field_name = values[0]
            field_type = values[1]
            field_values = values[2]
            mdq = MetaDataQuery(field_name,field_type,field_values)
            metadata_filters[field_name] = mdq
    return metadata_filters

def prepareComparisonData(metadata_file,metadata_filters,groups_db,comps_db):
    """Import the metadata and include/exclude fields/samples based on the user inputs"""
    firstLine = True
    samplesToRetain=[]
    samplesToRemove=[]
    uniqueDB={}
    indexesToFilterOn={}
    covariateSamples={}
    for line in open(metadata_file,'rU').xreadlines():
        data = line.rstrip()
        data = string.replace(data,'"','')
        values = string.split(data,'\t')
        if len(values)==1: continue
        if firstLine:
            headers = values
            covariateIndex=None
            index=0
            for h in headers:
                if h in metadata_filters:
                    md = metadata_filters[h]
                    if md.FieldType()=='UID':
                        uid_index = index
                    else:
                        indexesToFilterOn[index] = md
                index+=1
            firstLine = False
        else:
            try:  covariateType =  values[covariateIndex]
            except Exception: covariateType = None
            try: sampleID = values[uid_index] ### Must always be present
            except: continue ### Typically caused by rows with no values or variable blank values
            if '.bed' in sampleID:
                sampleID = string.replace(sampleID,'.bed','')
                
            for index in indexesToFilterOn:
                sample_value = values[index]
                md = indexesToFilterOn[index]
                if md.FieldType() == 'Unique' and 'TRUE' in md.FieldValues():
                    uniqueID_sample_name = values[index]
                    ### For example, unique donor ID (remove duplicate samples)
                    try: uniqueDB[uniqueID_sample_name].append(sampleID) 
                    except Exception: uniqueDB[uniqueID_sample_name] = [sampleID]
                elif md.FieldType() == 'Restrict':
                    if sample_value in md.FieldValues():
                        samplesToRetain.append(sampleID)
                    else:
                        #if 'M' not in md.FieldValues() and 'Training' not in md.FieldValues() and 'Bone marrow' not in md.FieldValues():
                        #print md.FieldValues(), md.FieldType(), sample_value, sampleID;sys.exit()
                        samplesToRemove.append(sampleID)
                        #if sampleID == 'SRS874639':
                        #print md.FieldValues(), md.FieldType(), [sample_value], sampleID;sys.exit()
                elif md.FieldType() == 'Exclude': 
                    if sample_value in md.FieldValues():
                        samplesToRemove.append(sampleID)
                elif md.FieldType() == 'Covariate':
                    """ This is the field we are deriving our groups and comps from """
                    ### If the covariateSamples key already present
                    if sample_value in md.FieldValues():
                        ### Hence, this sample is TRUE to the field name (e.g., diseaseX)
                        sample_type = md.FieldName()
                    else:
                        sample_type = 'Others'
                    if md.FieldName() in covariateSamples:
                        groups = covariateSamples[md.FieldName()]
                    else:
                        groups={}
                    if sample_type in groups:
                        groups[sample_type].append(sampleID)
                    else:
                        groups[sample_type] = [sampleID]
                    covariateSamples[md.FieldName()]=groups
                elif md.FieldType() == 'NumFilter':
                    cutoff,direction = md.FieldValues()
                    try:
                        if direction==">": 
                            if float(sample_value)>cutoff:
                                samplesToRetain.append(sampleID)
                            else:
                                samplesToRemove.append(sampleID)
                        if direction=="<": 
                            if float(sample_value)<cutoff:
                                samplesToRetain.append(sampleID)
                            else:
                                samplesToRemove.append(sampleID)
                    except Exception: ### Sample value not a float
                        samplesToRemove.append(sampleID)
    if len(samplesToRetain)==0 and len(samplesToRemove) ==0:
        for sample_type in groups:
            samplesToRetain+=groups[sample_type]
    #print len(list(set(samplesToRetain)))
    #print len(list(set(samplesToRemove)));sys.exit()
    for unique_donor_id in uniqueDB:
        if len(uniqueDB[unique_donor_id])>1:
            ### This method needs to be updated with a prioritization schema for duplicates, based on metadata
            for sample in uniqueDB[unique_donor_id][1:]:
                samplesToRemove.append(sample)
    samplesToRetainFinal=[]
    for sample in samplesToRetain:
        if sample not in samplesToRemove:
            samplesToRetainFinal.append(sample)
            
    for field in covariateSamples:
        groups = covariateSamples[field]
        for group_id in groups:
            updated_samples = []
            #print field, group_id, len(groups[group_id]), 
            for sample in groups[group_id]:
                if sample in samplesToRetainFinal:
                    updated_samples.append(sample)
            #print len(updated_samples)
            groups[group_id] = updated_samples
        if len(covariateSamples[field])>0:
            groups_db[field] = covariateSamples[field]
    
    ### Determine comparisons
    for field in groups_db:
        group_names=[]
        comps=[]
        groups = groups_db[field]
        for group_id in groups:
            if len(groups[group_id])>1: ### Ensure sufficient samples to compare
                if group_id not in group_names:
                    group_names.append(group_id)
        if len(group_names) == 2:
            if 'Others'in group_names:
                comps.append(tuple(group_names))
            if 'Others'not in group_names:
                comps.append(tuple(group_names))
        else:
            for group1 in group_names:
                for group2 in group_names:
                    if group1!=group2:
                        if (group2,group1) not in comps:
                            comps.append((group1,group2))
        comps_db[field]=comps
                        
    print len(list(set(samplesToRetainFinal))), 'samples considered for analysis.', len(list(set(samplesToRemove))), 'removed.'
    return groups_db,comps_db

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def importGroupsComps(groups_file):
    initial_groups={}
    group_id_lookup={}
    groups_db={}
    comps_db={}
    for line in open(groups_file,'rU').xreadlines():
        data = cleanUpLine(line)
        sample,ID,group_name = string.split(data,'\t')
        sample = string.replace(sample,'.bed','')
        group_id_lookup[ID]=group_name
        try: initial_groups[group_name].append(sample)
        except Exception: initial_groups[group_name] = [sample]
    groups_db[1] = initial_groups
    
    comps=[]
    comps_file = string.replace(groups_file,'groups.','comps.')
    for line in open(comps_file,'rU').xreadlines():
        data = cleanUpLine(line)
        g1,g2 = string.split(data,'\t')
        try: 
            g1_name = group_id_lookup[g1]
            g2_name = group_id_lookup[g2]
            comps.append((g1_name,g2_name))
        except Exception:
            """ The group ID in the comps file does not exist in the groups file """
            print g1, "or", g2, "groups do not exist in the groups file, but do in the comps."
    groups_db[1] = initial_groups
    comps_db[1] = comps
    return groups_db,comps_db
        
class PSIData:
    def __init__(self, clusterID, altexons, event_annotation, protein_predictions, coordinates):
        self.clusterID = clusterID
        self.altexons = altexons
        self.event_annotation = event_annotation
        self.protein_predictions = protein_predictions
        self.coordinates = coordinates
    def setUpdatedClusterID(self,updated_clusterID):
        self.updatedClusterID = updated_clusterID ### update this object
    def ClusterID(self): return self.clusterID
    def UpdatedClusterID(self): return self.updatedClusterID
    def AltExons(self): return self.altexons
    def EventAnnotation(self): return self.event_annotation
    def ProteinPredictions(self,dPSI):
        if dPSI<0:
            pp = string.replace(self.protein_predictions,'(+)','(--)')
            pp = string.replace(pp,'(-)','(+)')
            pp = string.replace(pp,'(--)','(-)')
            return pp
        else:
            return self.protein_predictions
    def Coordinates(self): return self.coordinates
    def Inclusion(self):
        """ Determines if the first junction is the inclusion or the second """
        junction1,junction2=string.split(self.Coordinates(),'|')
        j1s,j1e = string.split(string.split(junction1,':')[1],'-')
        j2s,j2e = string.split(string.split(junction2,':')[1],'-')
        j1_dist = abs(int(j1s)-int(j1e))
        j2_dist = abs(int(j2s)-int(j2e))
        if j1_dist>j2_dist:
            self.junction_type = 'exclusion'
        else:
            self.junction_type = 'inclusion'
        return self.junction_type
    def InclusionJunction(self):
        if self.junction_type == 'inclusion':
            return 'True'
        else:
            return 'False'
    def RelativeInclusion(self,dPSI):
        try: junction_type = self.junction_type
        except Exception: junction_type = self.Inclusion()
        if junction_type == 'exclusion' and dPSI<0:
            return 'inclusion'
        elif junction_type == 'inclusion' and dPSI<0:
            return 'exclusion'
        else:
            return junction_type

def performDifferentialExpressionAnalysis(species,platform,input_file,groups_db,comps_db,CovariateQuery,splicingEventTypes={}):
    ### filter the expression file for the samples of interest and immediately calculate comparison statistics

    firstLine = True
    group_index_db={}
    pval_summary_db={}
    group_avg_exp_db={}
    export_object_db={}
    psi_annotations={}
    rootdir=export.findParentDir(input_file)
    #print rootdir;sys.exit()
    if platform != 'PSI':
        try:
            import ExpressionBuilder
            expressionDataFormat,increment,convertNonLogToLog = ExpressionBuilder.checkExpressionFileFormat(input_file)
        except:
            increment = 0
            convertNonLogToLog=True
    else:
        convertNonLogToLog=True
        
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

    compared_ids={}
    row_count=0
    header_db={}
    headers_compared={}
    for line in open(input_file,'rU').xreadlines():
        row_count+=1
        if '.bed' in line:
            line = string.replace(line,'.bed','')
        data = line.rstrip()
        data = string.replace(data,'"','')
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
                group_index_db[group] = sample_index_list
                group_avg_exp_db[group] = {}
                
            ### Write out headers for grouped expression values
            for (group1,group2) in comps_db:
                eo = export_object_db[(group1,group2)]
                g1_headers = groups_db[group1]
                g2_headers = groups_db[group2]
                g1_headers = map(lambda x: group1+':'+x,g1_headers)
                g2_headers = map(lambda x: group2+':'+x,g2_headers)
                eo.write(string.join(['UID']+g1_headers+g2_headers,'\t')+'\n')
                #print (group1,group2),len(g1_headers),len(g2_headers)
                header_db[(group1,group2)] = g1_headers, g2_headers
            index=0
            uid_header=0
            for i in header:
                if i=='UID': uid_header=index
                if i=='AltExons': alt_exon_header=index
                if i=='ProteinPredictions': protein_prediction=index
                if i=='ClusterID': cluster_id_index=index
                if i=='Coordinates': coordinates_index=index
                if i=='EventAnnotation': event_annotation_index=index
                index+=1
            firstLine = False
        else:
            uid = values[uid_header]
            if platform == 'RNASeq':
                try: uid = symbol_to_gene[uid][0]
                except Exception: pass
            elif platform == 'PSI':
                clusterID = values[cluster_id_index]
                altexons = values[alt_exon_header]
                protein_predictions = values[protein_prediction]
                coordinates = values[coordinates_index]
                try: event_annotation = values[event_annotation_index]
                except:
                    continue ### Occurs with a filtered PSI value and no event annotation and anything follwoing
                ps = PSIData(clusterID, altexons, event_annotation, protein_predictions, coordinates)
                psi_annotations[uid]=ps
            group_expression_values={}
            original_group={}
            for group in group_index_db:
                sample_index_list = group_index_db[group]
                if platform != 'PSI':
                    try: filtered_values = map(lambda x: float(values[x]), sample_index_list) ### simple and fast way to reorganize the samples
                    except ValueError:
                        ### Strings rather than values present - skip this row
                        continue
                else: ### for splice-event comparisons where there a missing values at the end of the row
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
                    #if uid == 'ENSG00000105321:E3.2-E4.2 ENSG00000105321:E2.3-E4.2' and 'inner cell mass' in group:
                    #print filtered_values;sys.exit()
                if platform == 'PSI':
                    original_group[group]=unfiltered
                else:
                    original_group[group]=filtered_values
                if platform == 'RNASeq':# or platform == 'miRSeq':
                    if convertNonLogToLog:
                        filtered_values = map(lambda x: math.log(x+increment,2),filtered_values) ### increment and log2 adjusted
                group_expression_values[group] = filtered_values
            for groups in comps_db:
                group1,group2 = groups
                g1_headers,g2_headers=header_db[groups]
                header_samples = (tuple(g1_headers),tuple(g2_headers))
                if header_samples not in headers_compared:
                    headers_compared[header_samples]=[]
                    print len(g1_headers),len(g2_headers),groups
                try:
                    data_list1 = group_expression_values[group1]
                except KeyError:
                    ### This error is linked to the above Strings rather than values present error
                    continue
                data_list2 = group_expression_values[group2]
                combined = data_list1+data_list2
                if g1_headers==0 or g2_headers==0:
                    continue ### no samples matching the criterion
                try: diff = max(combined)-min(combined) ### Is there any variance?
                except Exception:
                    ### No PSI values for this splice-event and these groups
                    continue
                if diff==0:
                    continue
                elif (100.00*len(data_list1))/len(g1_headers)>=PercentExp and (100.00*len(data_list2))/len(g2_headers)>=PercentExp:
                    if len(data_list1)<minSampleNumber or len(data_list2)<minSampleNumber:
                        continue ### Don't analyze
                    compared_ids[uid]=[]
                    pass
                else:
                    continue ### Don't analyze
                if (len(data_list1)>1 and len(data_list2)>1) or (len(g1_headers)==1 or len(g2_headers)==1): ### For splicing data
                    ### Just a One-Way ANOVA at first - moderation happens later !!!!
                    try:
                        p = statistics.runComparisonStatistic(data_list1,data_list2,probability_statistic)
                    except Exception: p = 1
                    avg1 = statistics.avg(data_list1)
                    avg2 = statistics.avg(data_list2)
                    log_fold = avg1-avg2
                    if platform == 'RNASeq':# or platform == 'miRSeq':
                        max_avg = math.pow(2,max([avg1,avg2]))-1
                    else: max_avg = 10000

                    valid = True
                    try:
                        if max_avg<minRPKM:
                            log_fold = 'Insufficient Expression'
                    except Exception: pass
                    gs = statistics.GroupStats(log_fold,None,p)
                    gs.setAdditionalStats(data_list1,data_list2) ### Assuming equal variance
                    pval_db = pval_summary_db[groups] ### for calculated adjusted statistics
                    pval_db[uid] = gs ### store the statistics here
                    if len(restricted_gene_denominator)>0:
                        if uid not in restricted_gene_denominator:
                            proceed = False

                    ### store a new instance
                    gsg = statistics.GroupStats(log_fold,None,p)
                    gsg.setAdditionalStats(data_list1,data_list2) ### Assuming equal variance
                    global_adjp_db[CovariateQuery,groups,uid] = gsg ### for global adjustment across comparisons

                    group_avg_exp_db[group1][uid] = avg1 ### store the group expression values
                    group_avg_exp_db[group2][uid] = avg2 ### store the group expression values
                    if 'Insufficient Expression2' != log_fold:
                        #if abs(log_fold)>logfold_threshold:
                        eo = export_object_db[groups]
                        ls1 = map(str,original_group[group1])
                        ls2 = map(str,original_group[group2])
                        eo.write(string.join([uid]+ls1+ls2,'\t')+'\n')
                            
    for groups in export_object_db:
        export_object_db[groups].close()
    
    ### Calculate adjusted p-values for all pairwise comparisons
    global_count=0
    ensembls_found=False
    try: to = export.ExportFile(rootdir+'/top50/MultiPath-PSI.txt')
    except: pass
    for groups in pval_summary_db:
        psi_results = collections.OrderedDict()
        group1,group2 = groups
        group_comp_name = string.join(groups,'_vs_')
        if platform == 'PSI':
            filename=CovariateQuery+'/PSI.'+group_comp_name+'.txt'
        else:
            filename=CovariateQuery+'/GE.'+group_comp_name+'.txt'
        eo = export.ExportFile(rootdir+'/'+filename)
        #do = export.ExportFile(rootdir+'/Downregulated/'+filename)
        #uo = export.ExportFile(rootdir+'/Upregulated/'+filename)
        so = export.ExportFile(rootdir+'PValues/'+CovariateQuery+'-'+string.join(groups,'_vs_')+'.txt')
        header = 'GeneID\tSystemCode\tLogFold\trawp\tadjp\tSymbol\tavg-%s\tavg-%s\n' % (group1,group2)
        if platform == 'PSI':
            header = 'UID\tInclusion-Junction\tEvent-Direction\tClusterID\tUpdatedClusterID\tAltExons\tEventAnnotation\tCoordinates\tProteinPredictions\tdPSI\trawp\tadjp\tavg-%s\tavg-%s\n' % (group1,group2)
        eo.write(header)
        #do.write(header)
        #uo.write(header)
        so.write('Gene\tPval\n')
        pval_db = pval_summary_db[groups]
        if 'moderated' in probability_statistic: # and pval_threshold<1:
            try: statistics.moderateTestStats(pval_db,probability_statistic) ### Moderates the original reported test p-value prior to adjusting
            except Exception:
                print 'Moderated test failed... using student t-test instead',group1,group2
                #print traceback.format_exc()
        else:
            'Skipping moderated t-test...'
        statistics.adjustPermuteStats(pval_db) ### sets the adjusted p-values for objects
        
        pval_sort=[]
        for uid in pval_db:
            gs = pval_db[uid]
            pval_sort.append((gs.Pval(),uid))
        pval_sort.sort()

        ranked_count = 0
        for (pval,uid) in pval_sort:
            if ranked_count<11 and global_count<51:
                try: to.write(uid+'\n')
                except: pass
            gs = pval_db[uid]
            ranked_count+=1
            global_count+=1
            group1_avg = str(group_avg_exp_db[group1][uid])
            group2_avg = str(group_avg_exp_db[group2][uid])
            if use_adjusted_p:
                pval = float(gs.AdjP())
            else:
                pval = gs.Pval()
            if platform == 'miRSeq':
                symbol=[]
                altIDs = getAltID(uid)
                for id in altIDs:
                    if id in gene_to_symbol:
                        symbol.append(gene_to_symbol[id][0])
                        symbol.append(id)
                symbol = string.join(symbol,'|')
            elif uid in gene_to_symbol:
                symbols = unique.unique(gene_to_symbol[uid])
                symbol = string.join(symbols,'|')
            elif 'ENS' in uid and ':' in uid:
                ens_gene = string.split(uid,':')[0]
                try: symbol = gene_to_symbol[ens_gene][0]
                except Exception: symbol=''
            else:
                symbol = ''
            proceed = True
            ### Remove genes not a predetermined list (optional)
            if len(restricted_gene_denominator)>0:
                if uid not in restricted_gene_denominator:
                    if symbol not in restricted_gene_denominator:
                        proceed = False
            if 'Insufficient Expression' != gs.LogFold() and proceed:
                #pval_threshold=1; logfold_threshold=0
                if abs(gs.LogFold())>logfold_threshold and pval<=pval_threshold:
                    #print uid, groups, abs(gs.LogFold()), logfold_threshold;sys.exit()
                    if platform == 'PSI':
                        psi_results[uid]=gs,group1_avg,group2_avg
                        ### Write these out once all entries for that gene have been parsed
                    else:
                        if gs.LogFold()>0: fold = 'upregulated'
                        else: fold = 'downregulated'
                        try: splicingEventTypes[group_comp_name].append(fold)
                        except Exception: splicingEventTypes[group_comp_name] = [fold]
                        
                        if 'ENS' in uid and ensembls_found==False:
                            ensembls_found = True
                        else:
                            if 'ENS' not in uid:
                                system_code = 'Sy'; symbol = uid
                        if 'ENS' in uid:
                            system_code = 'En'
                        values = string.join([uid,system_code,str(gs.LogFold()),str(gs.Pval()),str(gs.AdjP()),symbol,group1_avg,group2_avg],'\t')+'\n'
                        eo.write(values)
                    proceed = True
                    if proceed:
                        if gs.LogFold()>0:
                            #uo.write(values)
                            pass
                        if gs.LogFold()<0:
                            #do.write(values)
                            pass
            so.write(uid+'\t'+str(gs.Pval())+'\n')
            
        if platform == 'PSI':
            ### Although we could have easily written the results above, we need to update the cluster ID here
            ### based on the significant results for this comparison only.
            gene_db={}
            for uid in psi_results:
                ps=psi_annotations[uid]
                geneID = string.split(uid,':')[1]
                if geneID in gene_db:
                    event_coordinates = gene_db[geneID]
                    event_coordinates[uid]=ps.Coordinates()
                else:
                    event_coordinates={}
                    event_coordinates[uid]=ps.Coordinates()
                    gene_db[geneID] = event_coordinates
            event_cluster_db=junctionGraph.createFeaturesFromEvents(gene_db)
            for uid in psi_results:
                updated_cluster_id=event_cluster_db[uid]
                ps=psi_annotations[uid]
                ps.setUpdatedClusterID(updated_cluster_id)
                gs,group1_avg,group2_avg=psi_results[uid]
                event_direction = ps.RelativeInclusion(gs.LogFold())
                incl_junction = ps.InclusionJunction()
                values = [uid,incl_junction,event_direction,ps.ClusterID(),ps.UpdatedClusterID(),ps.AltExons(),ps.EventAnnotation(),ps.Coordinates(),
                              ps.ProteinPredictions(gs.LogFold()),str(gs.LogFold()),str(gs.Pval()),str(gs.AdjP()),group1_avg,group2_avg]
                values = string.join(values,'\t')+'\n'
                eo.write(values)
                
                try: splicingEventTypes[group_comp_name].append((ps.ClusterID(),event_direction,ps.EventAnnotation()))
                except Exception: splicingEventTypes[group_comp_name] = [(ps.ClusterID(),event_direction,ps.EventAnnotation())]

        eo.close()
        #do.close()
        #uo.close()
        so.close()
        to.close()
    print len(compared_ids),'unique IDs compared (denominator).'
    return rootdir, splicingEventTypes

def outputGeneExpressionSummaries(rootdir,DEGs):
    eo = export.ExportFile(rootdir+'/gene_summary.txt')
    header = ['ComparisonName','RegulatedGenes','Upregulated','Downregulated']
    eo.write(string.join(header,'\t')+'\n')
    for comparison in DEGs:
        folds={}
        genes = DEGs[comparison]
        for fold in genes:
            try: folds[fold]+=1
            except Exception: folds[fold]=1
        try: up = folds['upregulated']
        except: up = 0
        try: down = folds['downregulated']
        except: down = 0
        eo.write(string.join([comparison,str(len(genes)),str(up),str(down)],'\t')+'\n')
    eo.close()
    
def outputSplicingSummaries(rootdir,splicingEventTypes):
    events_file = rootdir+'/event_summary.txt'
    eo = export.ExportFile(events_file)
    header = ['ComparisonName','UniqueJunctionClusters','InclusionEvents','ExclusionEvents']
    header+= ["alt-3'_exclusion","alt-3'_inclusion","alt-5'_exclusion","alt-5'_inclusion"]
    header+= ["alt-C-term_exclusion","alt-C-term_inclusion","altPromoter_exclusion","altPromoter_inclusion"]
    header+= ["cassette-exon_exclusion","cassette-exon_inclusion","intron-retention_exclusion"]
    header+= ["intron-retention_inclusion","trans-splicing_exclusion","trans-splicing_inclusion"]
    header = string.join(header,'\t')+'\n'
    eo.write(header)
    for comparison in splicingEventTypes:
        uniqueEvents={}
        for (clusterID,inclusionType,eventAnnotation) in splicingEventTypes[comparison]:
            try: uniqueEvents[clusterID].append((eventAnnotation,inclusionType))
            except Exception: uniqueEvents[clusterID] =[(eventAnnotation,inclusionType)]
        eventAnnotations={}
        inclusionEvents={}
        for clusterID in uniqueEvents:
            events = uniqueEvents[clusterID]
            events = list(set(events)) ### get unique (multiple distinct events can occur together)
            events.sort()
            ### Try to ignore events for the same clustering without an eventAnnotation
            for (eventAnnotation,inclusionType) in events:
                if len(eventAnnotation)>1:
                    if '|' in eventAnnotation: ### Occurs rarely - alt-3'|cassette-exon
                        annotations = string.split(eventAnnotation,'|')
                    else:
                        annotations = [eventAnnotation]
                    for annotation in annotations:
                        try: eventAnnotations[eventAnnotation+'_'+inclusionType]+=1
                        except Exception: eventAnnotations[eventAnnotation+'_'+inclusionType]=1
                try: inclusionEvents[inclusionType]+=1
                except Exception: inclusionEvents[inclusionType]=1

        try: a3e = str(eventAnnotations["alt-3'_exclusion"])
        except Exception: a3e = "0"
        try: a3i = str(eventAnnotations["alt-3'_inclusion"])
        except Exception: a3i = "0"
        try: a5e = str(eventAnnotations["alt-5'_exclusion"])
        except Exception: a5e = "0"
        try: a5i = str(eventAnnotations["alt-5'_inclusion"])
        except Exception: a5i = "0"
        try: aCe = str(eventAnnotations["alt-C-term_exclusion"])
        except Exception: aCe = "0"
        try: aCi = str(eventAnnotations["alt-C-term_inclusion"])
        except Exception: aCi = "0"
        try: ae = str(eventAnnotations["altPromoter_exclusion"])
        except Exception: ae = "0"
        try: ai = str(eventAnnotations["altPromoter_inclusion"])
        except Exception: ai = "0"
        try: ce = str(eventAnnotations["cassette-exon_exclusion"])
        except Exception: ce = "0"
        try: ci = str(eventAnnotations["cassette-exon_inclusion"])
        except Exception: ci = "0"
        try: ie = str(eventAnnotations["intron-retention_exclusion"])
        except Exception: ie = "0"
        try: ii = str(eventAnnotations["intron-retention_inclusion"])
        except Exception: ii = "0"
        try: te = str(eventAnnotations["trans-splicing_exclusion"])
        except Exception: te = "0"
        try: ti = str(eventAnnotations["trans-splicing_inclusion"])
        except Exception: ti = "0"
        try: incEvents = str(inclusionEvents['inclusion'])
        except Exception: incEvents = "0"
        try: exclEvents = str(inclusionEvents['exclusion'])
        except Exception: exclEvents = "0"
        values = [comparison,str(len(uniqueEvents)),incEvents,
                  exclEvents,a3e,a3i,a5e,a5i,
                  aCe,aCi,ae,ai,ce,ci,ie,ii,te,ti]
        values = string.join(values,'\t')+'\n'
        eo.write(values)
    eo.close()
    graphics = []
    try:
        from visualization_scripts import clustering
        parent_dir = export.findParentDir(events_file)
        parent_dir = export.findParentDir(parent_dir[:-1])
        OutputFile1 = parent_dir+'/SpliceEvent-Types.png'
        clustering.stackedbarchart(events_file,display=False,output=OutputFile1)
        OutputFile2 = parent_dir+'/SignificantEvents.png'
        index1=2;index2=3; x_axis='Number of Alternative Events'; y_axis = 'Comparisons'; title='MultiPath-PSI Alternative Splicing Events'
        clustering.barchart(events_file,index1,index2,x_axis,y_axis,title,display=False,color1='Orange',color2='SkyBlue',output=OutputFile2)
        graphics.append(['SpliceEvent-Types',OutputFile1])
        graphics.append(['Significant MultiPath-PSI Events',OutputFile2])
    except Exception:
        print traceback.format_exc()
    return graphics
            
def getAltID(uid):
    altID = string.replace(uid,'hsa-mir-','MIR')
    altID = string.replace(altID,'hsa-miR-','MIR')
    altID = string.replace(altID,'3p','')
    altID = string.replace(altID,'5p','')
    altID = string.upper(string.replace(altID,'hsa-let-','LET'))
    altID = string.replace(altID,'-','')
    altIDs = string.split(altID,'_')
    altIDs+=string.split(uid,'_')
    altIDs = unique.unique(altIDs)
    return altIDs

def getAnnotations(species,platform):
    import gene_associations
    if platform == 'RNASeq' or platform == 'PSI' or platform == 'miRSeq':
        try: gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
        except Exception: gene_to_symbol={}
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
                        
                            if platform == 'miRSeq' or platform == 'PSI' and use_adjusted_p == False:
                                ng_adjp = float(ug.Rawp())
                        except Exception:
                            if platform == 'miRSeq' or platform == 'PSI' and use_adjusted_p == False:
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

def exportUpDownGenes(results_dir):
    files = os.listdir(results_dir)
    for file in files:
        filename = results_dir+'/'+file
        output_dir = results_dir+'/regulated/'+file
        firstLine=True
        if '.txt' in filename and 'GE.' in filename:
            ou = export.ExportFile(output_dir[:-4]+'-up.txt')
            od = export.ExportFile(output_dir[:-4]+'-down.txt')
            for line in open(filename,'rU').xreadlines():
                data = line.rstrip()
                values = string.split(data,'\t')
                if firstLine:
                    firstLine=False
                    ou.write(line)
                    od.write(line)
                    lfi = values.index('LogFold')
                else:
                    if float(values[lfi]) >0:
                        ou.write(line)
                    else:
                        od.write(line)
            ou.close()
            od.close()
                    
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

def remoteAnalysis(species,expression_file,groups_file,platform='PSI',log_fold_cutoff=0.1,
                use_adjusted_pval=True,pvalThreshold=0.05,use_custom_output_dir=''):
    global pval_threshold
    global PercentExp
    global restricted_gene_denominator
    global global_adjp_db
    global probability_statistic
    global use_adjusted_p
    global logfold_threshold
    global minSampleNumber
    minSampleNumber = 2
    logfold_threshold = log_fold_cutoff
    use_adjusted_p = use_adjusted_pval
    global_adjp_db={}
    restricted_gene_denominator={}
    probability_statistic = 'moderated t-test'
    PercentExp = 0.5
    if platform == 'PSI':
        CovariateQuery = 'Events'
    else:
        CovariateQuery = 'DEGs'
    pval_threshold = pvalThreshold
    metadata_files = [groups_file]
    meta_description_file = None
    
    if platform == 'PSI' or platform == 'methylation':
        if log_fold_cutoff==None:
            logfold_threshold=0.1 ### equivalent to a 0.15 dPSI or 0.15 beta differences
        else:
            logfold_threshold= log_fold_cutoff
    if platform == 'PSI':
        print 'Using a dPSI of:',logfold_threshold
    if platform == 'methylation':
        use_adjusted_p = True
    if platform == 'miRSeq':
        use_adjusted_p = False
        logfold_threshold=math.log(1,2)
    print 'Filtering on adjusted p-value:',use_adjusted_p

    if platform != 'PSI':
        from build_scripts import EnsemblImport
        try: gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,platform,'key_by_array')
        except Exception: gene_location_db={}
    
    if meta_description_file !=None:
        if '.txt' in meta_description_file:
            meta_description_files = [meta_description_file]
        else:
            meta_description_files=[]
            files = os.listdir(meta_description_file)
            for file in files:
                if '.txt' in file:
                    meta_description_files.append(meta_description_file+'/'+file)
                
    splicingEventTypes={}
    all_groups_db={}
    all_comps_db={}

    if 'groups.' in metadata_files[0]:
        all_groups_db, all_comps_db = importGroupsComps(metadata_files[0])
    else:
        for meta_description_file in meta_description_files:
            metadata_filters = importMetaDataDescriptions(meta_description_file)
            all_groups_db, all_comps_db = prepareComparisonData(metadata_files[0],metadata_filters,all_groups_db, all_comps_db)

    for i in all_groups_db:
        #print i
        for k in all_groups_db[i]: print '  ',k,'\t',len(all_groups_db[i][k])
    #print all_comps_db
    
    if platform == 'PSI':
        result_type = 'dPSI'
    else:
        result_type = 'LogFold'
    
    if len(use_custom_output_dir)>0:
        CovariateQuery = use_custom_output_dir
    elif use_adjusted_p:
        CovariateQuery += '-'+result_type+'_'+str(logfold_threshold)[:4]+'_adjp'
    else:
        CovariateQuery += '-'+result_type+'_'+str(logfold_threshold)+'_rawp'

    for specificCovariate in all_groups_db:
        comps_db = all_comps_db[specificCovariate]
        groups_db = all_groups_db[specificCovariate]
        rootdir,splicingEventTypes = performDifferentialExpressionAnalysis(species,platform,expression_file,
                        groups_db,comps_db,CovariateQuery,splicingEventTypes)

    if platform == 'PSI':
        graphics = outputSplicingSummaries(rootdir+'/'+CovariateQuery,splicingEventTypes)
    else:
        graphics=[]
        outputGeneExpressionSummaries(rootdir+'/'+CovariateQuery,splicingEventTypes)
        #except: pass
    return graphics

def compareDomainComposition(folder):
    ### Compare domain composition

    import collections
    event_db = collections.OrderedDict()
    background_db = collections.OrderedDict()
    groups_list=['']
    import UI
    files = UI.read_directory(folder)
    for file in files:
        if '.txt' in file and 'PSI.' in file:
            db={}
            event_db[file[:-4]]=db
            all_unique_events = {}
            background_db[file[:-4]]=all_unique_events
            groups_list.append(file[:-4])
            fn = folder+'/'+file
            firstLine = True
            output_file = fn[:-4]+'-ProteinImpact.txt'
            output_file = string.replace(output_file,'PSI.','Impact.')
            eo = export.ExportFile(output_file)
            eo.write('GeneSymbol\tSy\tImpact\tImpactDescription\tUID\n')
            for line in open(fn,'rU').xreadlines():
                data = line.rstrip()
                t = string.split(data,'\t')
                #(+)alt-coding, (+)AA:232(ENSP00000230685)->296(ENSP00000009530)|(-)MHC_II-assoc_invar_chain-IPR022339, (+)DISULFID, (+)DOMAIN-Thyroglobulin type-1, (+)HELIX, (+)MHC_II-assoc_invar_chain-IPR022339, (+)SITE-Breakpoint for translocation to form a CD74-ROS1 fusion protein, (+)STRAND, (+)TOPO_DOM-Extracellular, (+)TURN, (+)Thyroglobulin_1-IPR000716
                if firstLine:
                    protein_index = t.index('ProteinPredictions')
                    clusterID_index = t.index('UpdatedClusterID')
                    event_index = t.index('EventAnnotation')
                    firstLine= False
                    continue
                uid = t[0]
                symbol = string.split(uid,':')[0]
                clusterID=t[clusterID_index]
                all_unique_events[clusterID]=None
                if len(t[protein_index])>0:
                    protein_prediction = t[protein_index]
                    protein_predictions = string.split(protein_prediction,',')
                    if 'Promoter' in t[event_index]:
                        continue
                    for annotation in protein_predictions:
                        if 'AA:' in annotation: 
                            if annotation[0]==' ':
                                direction = annotation[2]
                                diff = string.split(annotation[7:],'|')[0]
                            else:
                                direction = annotation[1]
                                diff = string.split(annotation[6:],'|')[0]
                            
                            try: diff1,diff2 = string.split(diff,'->')
                            except:
                                diff1,diff2 = string.split(diff,'+>')
                            try: diff1 = int(string.split(diff1,'(')[0])
                            except: print diff,diff1,diff2,[annotation];sys.exit()
                            diff2 = int(string.split(diff2,'(')[0])
                            coding_diff = (diff1)/float(diff2)
                            if coding_diff < 0.333:
                                type = 'Truncation'
                            elif coding_diff < 0.75:
                                type = 'Protein Length'
                            else:
                                type = None
                            if '|' in annotation:
                                domain_difference = True
                            else:
                                domain_difference = False
                            if type==None and domain_difference == False:
                                pass
                            else:
                                if direction == '+':
                                    direction = 'decreased' ### seems counter intuitive but we are specificall looking at decreased protein length
                                    """ specifically, if a (+)200->1000, this means in the experiment the protein is longer but truncation is decreased """
                                elif direction == '-':
                                    direction = 'increased'
                                if direction != '~':
                                    if type == 'Truncation':
                                        score = 2
                                    elif domain_difference:
                                        score = 1
                                    elif type == 'Protein Length':
                                        score = 0.5
                                    if direction == 'decreased':
                                        score = score*-1
                                    if type != None:
                                        updated_type = type
                                        if domain_difference:
                                            updated_type += '|domain-impact'
                                    else:
                                        updated_type = ''
                                        if domain_difference:
                                            updated_type += 'domain-impact'
                                    eo.write(string.join([symbol,'Sy',str(score),updated_type,uid],'\t')+'\n')
                                    db[clusterID]=type,domain_difference,direction
                            continue

    eo = export.ExportFile(folder+'/Protein-Level-Impact-Summary.txt')
    file_list = []
    for file in event_db:
        file_list.append(file)
        domain_difference_db={}
        impact_type_db={}
        impacted_events = event_db[file]
        for clusterID in impacted_events:
            type,domain_difference,direction = impacted_events[clusterID]
            ### Create the entries in the dictionary so these can be consistently populated below
            if domain_difference!=False:
                domain_difference_db[direction,domain_difference]=[]
            if type !=None:
                impact_type_db[direction,type]=[]

    for file in event_db:
        impacted_events = event_db[file]
        unique_events = len(background_db[file])
        file_domain_difference_db={}
        file_impact_type_db={}
        for clusterID in impacted_events:
            type,domain_difference,direction = impacted_events[clusterID]
            try:
                if domain_difference!=False:
                    file_domain_difference_db[direction,domain_difference]+=1
                if type!=None:
                    file_impact_type_db[direction,type]+=1
            except:
                if domain_difference!=False:
                    file_domain_difference_db[direction,domain_difference]=1
                if type!=None:
                    file_impact_type_db[direction,type]=1

        for (direction,type) in impact_type_db:
            try: impact_type_db[direction,type].append(str(file_impact_type_db[direction,type]/float(unique_events)))
            except: impact_type_db[direction,type].append('0')
        for (direction,domain_difference) in domain_difference_db: 
            try: domain_difference_db[direction,domain_difference].append(str(file_domain_difference_db[direction,domain_difference]/float(unique_events)))
            except: domain_difference_db[direction,domain_difference].append('0')
            
    eo.write(string.join(['UID']+file_list,'\t')+'\n')
    for (direction,domain_difference) in domain_difference_db:
            out = [direction+' Domain Impact']+domain_difference_db[direction,domain_difference]
            out = string.join(map(str,out),'\t')+'\n'

            eo.write(out)
    for (direction,type) in impact_type_db:
            out = [direction+' '+type]+impact_type_db[direction,type]
            out = string.join(map(str,out),'\t')+'\n'
            #print file, out, unique_events;sys.exit()
            eo.write(out)
    eo.close()
    sys.exit()
    
if __name__ == '__main__':
    species = 'Hs';
    expression_file = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Kumar/July-26-2017/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt'
    groups_file = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Kumar/July-26-2017/groups.KD.txt'
    computed_results_dir = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Leucegene/July-2017/PSI/SpliceICGS.R1.Depleted.12.27.17/all-depleted-and-KD'
    #exportUpDownGenes('/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/cellHarmony-evaluation/HCA-alignment/DEGs');sys.exit()
    #remoteAnalysis(species,expression_file,groups_file,platform='PSI',log_fold_cutoff=0.1,use_adjusted_pval=True,pvalThreshold=0.05);sys.exit()
    #compareDomainComposition(computed_results_dir)
    ################  Comand-line arguments ################
    #buildAdditionalMirTargetGeneSets();sys.exit()
    filename = '/Users/saljh8/Desktop/PCBC_MetaData_Comparisons/eXpress/CombinedResults/allTopGenes.txt' #DiffStateComps Reprogramming
    #exportGeneSetsFromCombined(filename);sys.exit()
    
    platform='RNASeq'
    species='Hs'
    probability_statistic = 'moderated t-test'
    #probability_statistic = 'Kolmogorov Smirnov'
    #probability_statistic = 'unpaired t-test'
    minRPKM=-1000
    PercentExp = 75
    minSampleNumber = 2
    logfold_threshold=math.log(2,2)
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
    meta_description_file=None
    sampleSetQuery = 'NA'
    used=[]
    executed_urls=[]
    restricted_gene_denominator={}
    global_adjp_db={}
    splicingEventTypes={}
    CovariateQuery = None
    log_fold_cutoff=None
    print sys.argv[1:]
    import getopt
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print 'Supply the argument --i location'
        ### python metaDataAnalysis.py --i /Volumes/SEQ-DATA\ 1/PCBC/RNASeq/July2014/MetaData/RNASeq_MetaData_July2014.txt --key "CellType" --value "Cell Type of Origin"
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['m=','i=','d=','c=','u=','p=','s=','f=',
            'g=','e=','ce=','rc=','cd=','o=','md=','in=','target=','parent=','urls=','used=','mf=',
            'adjp=','dPSI=','pval=','percentExp='])
        for opt, arg in options:
            if opt == '--m': metadata_files.append(arg)
            if opt == '--o':
                if output_dir==None:
                    output_dir = arg
                else:
                    output_dir = [output_dir,arg]
            if opt == '--i': expression_files.append(arg)
            if opt == '--e': runGOElite=True
            if opt == '--f':
                try: logfold_threshold = math.log(float(arg),2)
                except Exception: logfold_threshold = 0
            if opt == '--ce': compareEnrichmentProfiles = True
            if opt == '--d': sampleSetQuery=arg
            if opt == '--c': CovariateQuery=arg
            if opt == '--p': platforms.append(arg)
            if opt == '--g': gender_restricted=arg
            if opt == '--s': species=arg
            if opt == '--rc': restrictCovariateTerm=arg
            if opt == '--cd': compDiffState=arg
            if opt == '--md': mirDataDir=arg
            if opt == '--in': include_only=arg
            if opt == '--pval': pval_threshold = float(arg) 
            if opt == '--target': target_dir=arg
            if opt == '--mf': meta_description_file=arg
            if opt == '--parent': parent_syn=arg
            if opt == '--urls': executed_urls.append(arg) ### options are: all, junction, exon, reference
            if opt == '--used': used.append(arg)
            if opt == '--percentExp': PercentExp = int(arg)
            if opt == '--adjp':
                if string.lower(arg) == 'yes' or string.lower(arg) == 'true':
                    use_adjusted_p = True
            if opt == '--dPSI':
                log_fold_cutoff = float(arg)
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
        if platform == 'PSI' or platform == 'methylation':
            if log_fold_cutoff==None:
                logfold_threshold=0.1 ### equivalent to a 0.15 dPSI or 0.15 beta differences
            else:
                logfold_threshold= log_fold_cutoff
        if platform == 'PSI':
            print 'Using a dPSI of:',logfold_threshold
        if platform == 'methylation':
            use_adjusted_p = True
        if platform == 'miRSeq':
            use_adjusted_p = False
            logfold_threshold=math.log(1,2)
        print 'Filtering on adjusted p-value:',use_adjusted_p

        if platform == 'PSI' and CovariateQuery == None:
            CovariateQuery = 'Events'
        else:
            CovariateQuery = 'DEGs'
        
        if platform != 'PSI':
            from build_scripts import EnsemblImport
            try: gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,platform,'key_by_array')
            except Exception: gene_location_db={}
        
        if len(include_only)>0:
            restricted_gene_denominator = importRestrictedSetOfGenesToQuery(include_only)
        
        if meta_description_file !=None:
            if '.txt' in meta_description_file:
                meta_description_files = [meta_description_file]
            else:
                meta_description_files=[]
                files = os.listdir(meta_description_file)
                for file in files:
                    if '.txt' in file:
                        meta_description_files.append(meta_description_file+'/'+file)
                    
        splicingEventTypes={}
        all_groups_db={}
        all_comps_db={}

        if 'groups.' in metadata_files[0]:
            all_groups_db, all_comps_db = importGroupsComps(metadata_files[0])
        else:
            for meta_description_file in meta_description_files:
                metadata_filters = importMetaDataDescriptions(meta_description_file)
                all_groups_db, all_comps_db = prepareComparisonData(metadata_files[0],metadata_filters,all_groups_db, all_comps_db)

        for i in all_groups_db:
            print i
            for k in all_groups_db[i]: print '  ',k,'\t',len(all_groups_db[i][k])
        print all_comps_db
        
        if platform == 'PSI':
            result_type = 'dPSI'
        else:
            result_type = 'LogFold'
        if use_adjusted_p:
            CovariateQuery += '-'+result_type+'_'+str(logfold_threshold)[:4]+'_adjp'
        else:
            CovariateQuery += '-'+result_type+'_'+str(logfold_threshold)+'_rawp'

        for specificCovariate in all_groups_db:
            comps_db = all_comps_db[specificCovariate]
            groups_db = all_groups_db[specificCovariate]
            rootdir,splicingEventTypes = performDifferentialExpressionAnalysis(species,platform,expression_file,groups_db,comps_db,CovariateQuery,splicingEventTypes)

        if platform == 'PSI':
            graphics = outputSplicingSummaries(rootdir+'/'+CovariateQuery,splicingEventTypes)
        else:
            outputGeneExpressionSummaries(rootdir+'/'+CovariateQuery,splicingEventTypes)
            #except: pass
            
        sys.exit()
        for CovariateQuery in covariate_set:
          for sampleSetQuery in Sample_set:
            print 'Analyzing the covariate:',CovariateQuery, 'and diffState:',sampleSetQuery, 'unique donor analysis:',uniqueDonors
            if 'XIST' in CovariateQuery: gender_restricted='female'
            
            genderRestricted = gender_restricted
            try:
                sample_metadata,groups_db,comps_db = prepareComparisonData(metadata_file,metadata_filters,sampleSetQuery,CovariateQuery,uniqueDonors,genderRestricted,platform=platform,compDiffState=compDiffState,restrictCovariateTerm=restrictCovariateTerm)
                performDifferentialExpressionAnalysis(species,platform,expression_file,sample_metadata,groups_db,comps_db,sampleSetQuery+'-'+CovariateQuery,uniqueDonors)
            except Exception:
                print traceback.format_exc()
            if runAgain:
                uniqueDonors=True
                use_adjusted_p = False
                print 'Analyzing the covariate:',CovariateQuery, 'and diffState:',sampleSetQuery, 'unique donor analysis:',uniqueDonors
                try:
                    sample_metadata,groups_db,comps_db = prepareComparisonData(metadata_file,sampleSetQuery,CovariateQuery,uniqueDonors,genderRestricted,platform=platform,compDiffState=compDiffState,restrictCovariateTerm=restrictCovariateTerm)
                    performDifferentialExpressionAnalysis(species,platform,expression_file,sample_metadata,groups_db,comps_db,sampleSetQuery+'-'+CovariateQuery,uniqueDonors)
                except Exception: pass
                uniqueDonors=False; use_adjusted_p = True
                if platform == 'miRSeq' or platform == 'PSI': use_adjusted_p = False
            
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
        
        cellLines1,sample_metadata1,uniqueDonor_db,donor_db = prepareComparisonData(metadata_file1,sampleSetQuery,None,uniqueDonors,gender_restricted,platform=platform1)
        cellLines2,sample_metadata2,uniqueDonor_db,donor_db = prepareComparisonData(metadata_file2,sampleSetQuery,None,uniqueDonors,gender_restricted,platform=platform2)
        exp_cellLines1 = getDatasetSamples(expression_file1,sample_metadata1,cellLines1)
        exp_cellLines2 = getDatasetSamples(expression_file2,sample_metadata2,cellLines2)
        common_lines = getCommonCellLines(cellLines1,cellLines2,exp_cellLines1,exp_cellLines2,uniqueDonor_db,uniqueDonors,donor_db)
        filtered_exp_data1, samplesToEvaluate = importExpressionData(species,platform1,expression_file1,cellLines1,common_lines)
        filtered_exp_data2, samplesToEvaluate = importExpressionData(species,platform2,expression_file2,cellLines2,common_lines)
        combineAndCompareMatrices(expression_file1,filtered_exp_data1,filtered_exp_data2,platform1,platform2,samplesToEvaluate) ### export results to two different files and a combined file
