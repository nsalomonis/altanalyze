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
import os.path
import unique
import export
import ExpressionBuilder
import copy
import traceback

try:
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
        from stats_scripts import statistics
        import math
        from scipy import stats
        use_scipy = True
except Exception:
    use_scipy = False ### scipy is not required but is used as a faster implementation of Fisher Exact Test when present
    
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt"or entry[-4:] == ".tab" or entry[-4:] == ".csv" or '.fa' in entry: dir_list2.append(entry)
    return dir_list2

class GrabFiles:
    def setdirectory(self,value):
        self.data = value
    def display(self):
        print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        file_dirs = getDirectoryFiles(self.data,str(search_term))
        if len(file_dirs)<1: print search_term,'not found',self.data
        return file_dirs
    
def getDirectoryFiles(import_dir, search_term):
    exact_file = ''; exact_file_dirs=[]
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir: exact_file_dirs.append(affy_data_dir)
    return exact_file_dirs

########## End generic file import ##########

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def createCorrelationTemplate(tissues):
    i=0; tissue_template_db={}
    null_template = [0]*len(tissues)
    for tissue in tissues:
        tissue_template = list(null_template) ### creates a modifiable copy
        tissue_template[i] = correlationDirection ### make -1 reverse profile (anti-correlated)
        #if 'Atrium' in tissue:
        tissue_template_db[tissue] = tissue_template
        i+=1
    return tissue_template_db

def findHousekeepingGenes(uid,data_list1):
    ### Pearson won't work for no variance, so have to use an alternative approach
    stdev = statistics.stdev(data_list1)
    housekeeping.append([stdev,uid]) #stdev/min(data_list1)
    
def expressedIndexes(values):
    filtered_indexes=[]
    i=0
    for value in values:
        if value!=None:
            filtered_indexes.append(i)
        i+=1
    return filtered_indexes
    
def advancedPearsonCorrelationAnalysis(uid,data_list1,tissue_template_db):
    expIndexes = expressedIndexes(data_list1)
    Queried[uid]=[]
    if (float(len(expIndexes))/len(data_list1))>0.0: ### Atleast 50% of samples evaluated express the gene
        data_list = map(lambda i: data_list1[i],expIndexes) ### Only expressed values (non-None)
        max_diff = max(data_list)-statistics.avg(data_list)
        if max_diff>-1000 and max(data_list)>-1000:
            if correlateAllGenes:
                min_rho = -1
                min_p = 1
            else:
                min_rho = 0.3
                min_p = 0.05
            for tissue in tissue_template_db:
                tissue_template = tissue_template_db[tissue]
                c1 = tissue_template.count(1)
                filtered_template = map(lambda i: tissue_template[i],expIndexes)
                c2 = filtered_template.count(1)
                if len(data_list)!= len(filtered_template): kill
                if c1 == c2 or c1 != c2: ### If number of 1's in list1 matches list2
                    rho,p = rhoCalculation(data_list,filtered_template)
                    if tissue == 'Housekeeping':
                        print tissue, p, uid;sys.exit()
                    #if rho>min_rho:
                    if p<min_p and rho>min_rho:
                        Added[uid]=[]
                        try: tissue_scores[tissue].append([(rho,p),uid])
                        except Exception: tissue_scores[tissue] = [[(rho,p),uid]]

def PearsonCorrelationAnalysis(uid,data_list1,tissue_template_db):
    if correlateAllGenes:
        min_rho = -2
        min_p = 1
    else:
        min_rho = 0.3
        min_p = 0.05
    for tissue in tissue_template_db:
        tissue_template = tissue_template_db[tissue]
        rho,p = rhoCalculation(data_list1,tissue_template)
        #print rho, uid, tissue
        if tissue == 'Housekeeping':
            print tissue, rho, uid;sys.exit()
        #if rho>min_rho:
        if p<min_p and rho>min_rho:
            try: tissue_scores[tissue].append([(rho,p),uid])
            except Exception: tissue_scores[tissue] = [[(rho,p),uid]]

def rhoCalculation(data_list1,tissue_template):
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides import warnings
            try:
                rho,p = stats.pearsonr(data_list1,tissue_template)
                return rho,p
            except Exception:
                #data_list_alt = [0 if x==None else x for x in data_list1]
                #rho,p = stats.pearsonr(data_list1,tissue_template)
                kill
    except Exception:
        rho = pearson(data_list1,tissue_template)
        return rho

def simpleScipyPearson(query_lists,reference_list):
    """ Get the top correlated values referenced by index of the query_lists (e.g., data matrix) """
    i=0
    rho_results=[]
    for query_list in query_lists:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides import warnings
            rho,p = stats.pearsonr(query_list,reference_list)
            #if query_list == reference_list: print query_list,reference_list, rho;sys.exit()
        if str(rho)!='nan':
            rho_results.append([(float(rho),float(p)),i])
        i+=1
    rho_results.sort()
    rho_results.reverse()
    return rho_results
    
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
    try: r = sum_a/math.sqrt(sum_b*sum_c)
    except Exception: r =0
    return r

def avg(array):
    return sum(array)/len(array)
            
def getArrayData(filename,filtered_probeset):
    fn=filepath(filename); x=0; k=0; expression_data={}; annotations={}; expression_data_str={}
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)  #remove endline
        t = string.split(data,'\t')
        if x == 0:
            i=0
            for h in t:
                if 'Definition' in h: lti = i ### last tissue index
                if 'Description' in h: lti = i ### last tissue index
                i+=1
            x=1
            try:
                header = t[1:lti]
                annotation_header = t[lti:]
            except Exception:
                lti = len(t)
                header=[]
                annotation_header=[]
        else:
            probeset = t[0]
            try:
                if len(filtered_probeset)>0: ### Otherwise, just get all annotations
                    null = filtered_probeset[probeset]
                    exp_values = map(float, t[1:lti])
                    if log_transform:
                        try: exp_values =  map(lambda x: math.log(x,2), exp_values)
                        except Exception:
                            exp_values =  map(lambda x: math.log(x+1,2), exp_values)
                    expression_data_str[probeset] = map(str,exp_values)
                    expression_data[probeset] = exp_values
                try: annotations[probeset] = t[lti:]
                except KeyError:
                    annotations[probeset] = []
            except KeyError:
                null=[]
    
    sum_tissue_exp = {}
    for probeset in expression_data:
        i=0
        for fold in expression_data[probeset]:
            try: sum_tissue_exp[i].append(fold)
            except Exception: sum_tissue_exp[i] = [fold]
            i+=1

    expression_relative={}
    for probeset in expression_data:
        i=0
        for fold in expression_data[probeset]:
            ratio = str(fold/max(sum_tissue_exp[i]))
            try: expression_relative[probeset].append(ratio)
            except Exception: expression_relative[probeset] = [ratio]
            i+=1
    return expression_data_str, annotations, header, annotation_header

def importMarkerProfiles(filename,fl):
    x=0
    ### Import correlated marker results
    fn=filepath(filename)
    marker_list = []
    condition_list=[]
    marker_db={}
    probeset_symbol_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            x=1
        else:
            uid,symbol,rho,p,condition=t
            probeset_symbol_db[uid]=symbol
            try: marker_db[condition].append(uid)
            except Exception: marker_db[condition] = [uid]
            if condition not in condition_list:
                condition_list.append(condition)
    marker_condition_db={}
    
    try: condition_list = getOrderedGroups(fl.DatasetFile()) ### Order this set with the same way as the samples on the opposite axis
    except Exception,e:
        #print 'failed',e
        condition_list.sort()
        condition_list.reverse() ### This makes this ordered top-down for clustering
        
    for condition in condition_list:
        if condition in marker_db:
            for uid in marker_db[condition]:
                if uid not in marker_list:
                    marker_list.append(uid) ### ranked and unique marker list
                    marker_condition_db[uid] = condition
    
    exportMarkersForGOElite(filename,marker_db,fl) ### Export these lists for GO-Elite
    
    return marker_list, probeset_symbol_db, marker_condition_db

def exportMarkersForGOElite(filename,gene_db,fl):
    if fl.Vendor() == 'Affymetrix': system = 'X'
    elif fl.Vendor() == 'Agilent': system = 'Ag'
    elif fl.Vendor() == 'Illumina': system = 'Il'
    elif 'other:' in fl.Vendor():
        system = string.replace(fl.Vendor(),'other:','')
        if system == 'Symbol': system = 'Sy'
        else: system = 'Ma'
    else: system = 'Sy'
    root_dir = fl.OutputDir()
    if 'ReplicateBased' in filename: suffix = '-ReplicateBased'
    if 'MeanBased' in filename: suffix = '-MeanBased'
    
    for markerSet in gene_db:
        header = string.join(['Gene','System','Hit'],'\t')+'\n'
        filename = root_dir+'GO-Elite/MarkerFinder/'+markerSet+suffix+'.txt'
        export_obj = export.ExportFile(filename)
        export_obj.write(header)
        for gene in gene_db[markerSet]:
            if 'ENS' in gene:
                system = 'En'
            try: system = system
            except Exception: system = 'Swiss'
            values = string.join([gene,system,'1'],'\t')+'\n'
            export_obj.write(values)
        export_obj.close()
            
def reorderInputFile(custom_path,marker_list,marker_condition_db):
    x=0
    ### Import correlated marker results
    fn=filepath(custom_path)
    exp_db={}
    probeset_symbol_db={}
    #print custom_path;sys.exit()
    #print fn
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            header = line
            x=1
        else:
            uid = string.join(string.split(t[0],' ')[0:-1],' ') ### IDs can have spaces
            if '///' in uid:
                uid = string.split(uid,' ')[0] ### Affy ID with multiple annotations
            exp_db[uid] = line
                
    ### Over-write read in file
    export_obj = export.ExportFile(custom_path)
    export_obj.write(header)
    marker_list.reverse() ### Reverse the order of the MarkerFinder results
    for uid in marker_list:
        condition = marker_condition_db[uid]
        new_uid = condition+':'+uid
        if uid in exp_db:
            export_obj.write(condition+':'+exp_db[uid])
        elif new_uid in exp_db:
            export_obj.write(exp_db[new_uid])
        else:
            """
            print [uid], len(exp_db)
            for i in exp_db:
                print [i];break
            print 'Error encountered with the ID:',uid, 'not in exp_db'; kill
            """
            pass
    export_obj.close()
    
def getOrderedGroups(filename):
    group_list=[]
    filename = string.replace(filename,'///','/')
    filename = string.replace(filename,'//','/')
    if 'ExpressionOutput' in filename:
        filename = string.replace(filename,'ExpressionOutput','ExpressionInput')
        filename = string.replace(filename,'-steady-state','')
        filename = string.replace(filename,'DATASET-','exp.')
        groups_dir = string.replace(filename,'exp.','groups.')
        group_names_db = ExpressionBuilder.simpleGroupImport(groups_dir)[3]
        
    for i in group_names_db: group_list.append(i)
    group_list.reverse()
    return group_list
        
def generateMarkerHeatMaps(fl,platform,convertNonLogToLog=False,graphics=[],Species=None):
    from visualization_scripts import clustering
    """ From the generated marker sets, output the replicate input data """
    marker_root_dir = fl.OutputDir()+'/'+'ExpressionOutput/MarkerFinder'
    #print 1,fl.DatasetFile()
    #print 2, fl.Vendor()
    for marker_dir in read_directory(marker_root_dir):
        if 'MarkerGenes' in marker_dir and 'correlation' in marker_dir:
            marker_dir = marker_root_dir+'/'+marker_dir
            marker_list, probeset_symbol_db, marker_condition_db = importMarkerProfiles(marker_dir,fl)
            custom_path = string.replace(marker_dir,'MarkerGenes','Clustering/MarkerGenes')
            """
            print fl.DatasetFile()
            print len(marker_list), marker_list[:3]
            print len(probeset_symbol_db)
            print custom_path
            print convertNonLogToLog
            print Species
            print platform
            print len(probeset_symbol_db)
            print custom_path
            print fl.Vendor()
            print convertNonLogToLog
            """
            ExpressionBuilder.exportGeometricFolds(fl.DatasetFile(),platform,marker_list,probeset_symbol_db,exportOutliers=False,exportRelative=False,customPath=custom_path,convertNonLogToLog=convertNonLogToLog)
            reorderInputFile(custom_path,marker_list, marker_condition_db)
            row_method = None; row_metric = 'cosine'; column_method = None; column_metric = 'euclidean'; color_gradient = 'yellow_black_blue'; transpose = False
            import UI
            gsp = UI.GeneSelectionParameters(Species,platform,fl.Vendor())
            gsp.setPathwaySelect('None Selected')
            gsp.setGeneSelection('')
            gsp.setOntologyID('')
            gsp.setGeneSet('None Selected')
            gsp.setJustShowTheseIDs('')                
            gsp.setTranspose(False)
            gsp.setNormalize('median')
            gsp.setGeneSelection('')
            gsp.setClusterGOElite('GeneOntology')
            #gsp.setClusterGOElite('BioMarkers')
            """
            print custom_path
            print graphics
            print row_method
            print column_method
            print column_metric
            """
            reload(clustering)
            try:
                graphics = clustering.runHCexplicit(custom_path, graphics, row_method, row_metric,
                                    column_method, column_metric, color_gradient, gsp, contrast=4, display=False)
            except Exception:
                print traceback.format_exc()
                print 'Error occured in generated MarkerGene clusters... see ExpressionOutput/MarkerFinder files.'
    return graphics

def reorderMultiLevelExpressionFile(input_file):
    ### Takes an input file and re-orders it based on the order in the groups file... needed for multi-level expression file with replicates
    from import_scripts import sampleIndexSelection
    output_file = input_file[:-4]+'-output.txt'
    filter_file = string.replace(input_file,'-steady-state','')
    filter_file = string.replace(filter_file,'exp.','groups.')
    filter_file = string.replace(filter_file,'stats.','groups.')
    filter_file = string.replace(filter_file,'topSplice.','groups.')
    filter_file = string.replace(filter_file,'filter.','groups.')
    filter_names = sampleIndexSelection.getFilters(filter_file)
    sampleIndexSelection.filterFile(input_file,output_file,filter_names)
    c1 = verifyFileLength(input_file)
    c2 = verifyFileLength(output_file)
    if c1==c2:
        os.remove(input_file)
        export.copyFile(output_file, input_file)

def verifyFileLength(filename):
    count = 0
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            if line[0]!='#':
                count+=1
    except Exception: null=[]
    return count
    
def analyzeData(filename,Species,Platform,codingType,geneToReport=60,correlateAll=True,AdditionalParameters=None,logTransform=False,binarize=False):
    global genesToReport; genesToReport = geneToReport
    global correlateAllGenes; correlateAllGenes = correlateAll
    global all_genes_ranked; all_genes_ranked={}
    global RPKM_threshold; global correlationDirection
    global Added; Added={}; global Queried; Queried={}
    
    """
    print 4,Platform, codingType, geneToReport, correlateAll, logTransform,
    try:
        #print AdditionalParameters.CorrelationDirection()
        print AdditionalParameters.RPKMThreshold()
    except Exception:
        print 'nope'
    """
    global AvgExpDir
    if len(filename) == 2:
        filename, AvgExpDir = filename  #### Used when there are replicate samples: avg_exp_dir is non-replicate
        if AvgExpDir==None:
            AvgExpDir = string.replace(filename,'-steady-state','')
            AvgExpDir = string.replace(AvgExpDir,'exp.','AVERAGE-')
            AvgExpDir = string.replace(AvgExpDir,'ExpressionInput','ExpressionOutput')
    if 'ExpressionOutput' in filename:
        use_replicates = False
    else:
        use_replicates = True
    
    import RNASeq
    try: Platform = RNASeq.checkExpressionFileFormat(filename,Platform)
    except Exception: Platform = "3'array"

    try: RPKM_threshold = AdditionalParameters.RPKMThreshold() ### Used for exclusion of non-expressed genes
    except Exception:
        pass
    if Platform == 'RNASeq':
        try: RPKM_threshold = AdditionalParameters.RPKMThreshold() ### Used for exclusion of non-expressed genes
        except Exception: RPKM_threshold = 1; logTransform = True

    correlationDirection = 1.00 ### Correlate to a positive or negative idealized pattern
    try:
        if AdditionalParameters.CorrelationDirection() != 'up' and AdditionalParameters.CorrelationDirection() != 'positive':
            correlationDirection = -1.00
    except Exception: pass
    #print correlationDirection
    fn=filepath(filename); x=0; t2=['ID']; cluster_db={}; cluster_list = []; global coding_type; coding_type = codingType
    global cluster_comps; cluster_comps = []; global compare_clusters; compare_clusters = 'no'
    global housekeeping; housekeeping=[]; global analyze_housekeeping; analyze_housekeeping = 'no'
    global species; global platform; species = Species; platform = Platform; global log_transform
    log_transform=logTransform
    
    #if 'topSplice.' not in fn and 'steady' not in fn and 'AVERAGE' not in fn and 'DATASET' not in fn: reorderMultiLevelExpressionFile(fn)
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x == 0:
            updated_names = ['ID']
            correlations = 'single'
            tissue_template_db,group_sample_db = getReplicateData(fn,t[1:])
            if '~' in data: correlations = 'multiple'
            elif t[1] in group_sample_db:
                if '~' in group_sample_db[t[1]]:
                    correlations = 'multiple'
                    for i in t[1:]:
                        updated_names.append(group_sample_db[i])
                    t = updated_names
                    
            i=0
            for h in t:
                if correlations == 'multiple':
                    if '~' in h:
                        cluster, group_name = string.split(h,'~')
                        cluster = int(cluster)
                        try: cluster_db[cluster].append(i)
                        except Exception: cluster_db[cluster] = [i]
                i+=1
            if correlations == 'multiple':
                compare_clusters = 'yes'
                ### If there are multiple sample group clusters
                for cluster in cluster_db: cluster_list.append(cluster)
                cluster_list.sort()
                if cluster_list[0]==0: ### All clusters should then be compared to this first cluster if the first is named 0
                    cluster_comp_db={}; cluster_comp_db[0] = cluster_db[0]
                    cluster_comps.append(cluster_comp_db) ### Add the first by itself (cluster 0 samples only comapred among themselves initially)
                    for cluster in cluster_list[1:]:
                        cluster_comp_db={}; cluster_comp_db[0] = cluster_db[0]
                        cluster_comp_db[cluster] = cluster_db[cluster]
                        cluster_comps.append(cluster_comp_db) ### Each non-zero cluster compared to cluster 0 as pairwise group combinations
                else:
                    for cluster in cluster_list:
                        cluster_comp_db={}
                        cluster_comp_db[cluster] = cluster_db[cluster]
                        cluster_comps.append(cluster_comp_db)
            x = 1
            break

    iteration=1
    if compare_clusters == 'yes':
        tissue_specific_IDs_combined={}; correlations_combined={}
        for cluster_comp_db in cluster_comps:
            ###Interate through each comparison
            print 'Iteration',iteration,'of',len(cluster_comps)
            tissue_specific_IDs,interim_correlations,annotation_headers,tissues = identifyMarkers(filename,cluster_comp_db,binarize=binarize)
            iteration+=1
            for tissue in tissue_specific_IDs:
                if tissue not in tissue_specific_IDs_combined: ### Combine the tissue results from all of the cluster group analyses, not over-writing the existing 
                    tissue_specific_IDs_combined[tissue] = tissue_specific_IDs[tissue]
                    correlations_combined[tissue] = interim_correlations[tissue]

        tissue_specific_IDs={}; interim_correlations={}
        for tissue in tissue_specific_IDs_combined:
            for probeset in tissue_specific_IDs_combined[tissue]:
                try: tissue_specific_IDs[probeset].append(tissue)
                except Exception: tissue_specific_IDs[probeset] = [tissue]
            for (probeset,symbol,(rho,p)) in correlations_combined[tissue]:
                try: interim_correlations[tissue].append([probeset,symbol,(rho,p)])
                except Exception: interim_correlations[tissue] = [[probeset,symbol,(rho,p)]]
        analyze_housekeeping = 'yes'; compare_clusters = 'no'
        original_tissue_headers2 = original_tissue_headers ### The last function will overwrite the group~ replacement
        #identifyMarkers(filename,[]) ### Used to get housekeeping genes for all conditions
    else:
        tissue_specific_IDs,interim_correlations,annotation_headers,tissues = identifyMarkers(filename,[],binarize=binarize)
        original_tissue_headers2 = original_tissue_headers
    ### Add a housekeeping set (genes that demonstrate expression with low variance
    housekeeping.sort(); ranked_list=[]; ranked_lookup=[]; tissue = 'Housekeeping'
    for (stdev,(probeset,symbol)) in housekeeping:
        if probeset not in tissue_specific_IDs: ### Shouldn't be if it is a housekeeping gene
            if symbol not in ranked_list:
                ranked_list.append(symbol); ranked_lookup.append([probeset,symbol,(stdev,0)])

    ### Replicates code in identifyMarkers - but only applied to housekeeping genes to add those in addition to the existing ones in tissue_specific_IDs
    for (probeset,symbol,(stdev,p)) in ranked_lookup[:genesToReport]:
        try: tissue_specific_IDs[probeset].append(tissue)
        except Exception: tissue_specific_IDs[probeset] = [tissue]
        try: interim_correlations[tissue].append([probeset,symbol,(stdev,p)])
        except Exception: interim_correlations[tissue] = [[probeset,symbol,(stdev,p)]]
    ### If no mean file provided
    #print [use_replicates, filename, tissue]
    if use_replicates:
        try: filename = AvgExpDir
        except Exception: pass ### For AltExon queries
    try:
        expression_relative,annotations,tissue_headers, annot_header = getArrayData(filename,tissue_specific_IDs)
        if use_replicates:
            original_tissue_headers2, annotation_headers = tissue_headers, annot_header

        tissue_specific_IDs2 = copy.deepcopy(tissue_specific_IDs)
        for probeset in tissue_specific_IDs2:
            if probeset in annotations:
                annotations[probeset]+=[string.join(list(tissue_specific_IDs[probeset]),'|')] ### Save as string
        title_row = ['UID']+annotation_headers+['marker-in']+original_tissue_headers2
        export_dir = exportMarkerGeneProfiles(filename,annotations,expression_relative,title_row)
    except Exception,e:
        #print traceback.format_exc()
        pass

    exportCorrelations(filename,interim_correlations)
    if correlateAllGenes:
        exportAllGeneCorrelations(filename,all_genes_ranked)
    try: return export_dir
    except Exception: pass
    
def getReplicateData(expr_input,t):
    groups_dir = string.replace(expr_input,'exp.','groups.')
    groups_dir = string.replace(groups_dir,'stats.','groups.')
    groups_dir = string.replace(groups_dir,'topSplice.','groups.')
    groups_dir = string.replace(groups_dir,'filter.','groups.')
    groups_dir = string.replace(groups_dir,'-steady-state','') ### groups is for the non-steady-state file
    #groups_dir = string.replace(groups_dir,'-average.txt','.txt') ### groups is for the non-steady-state file
    if 'groups.' not in groups_dir and 'AltResults' in groups_dir:
            parent_dir = string.split(expr_input,'AltResults')[0]
            file = export.findFilename(expr_input)
            file = string.replace(file,'AltExonConfirmed-','groups.')
            file = string.replace(file,'AltExon-','groups.')
            groups_dir = parent_dir+'ExpressionInput/'+file
            
    group_index_db={}
    splitHeaders=False
    for i in t:
        if '~' not in i:splitHeaders=True

    ### use comps in the future to visualize group comparison changes
    sample_list,group_sample_db,group_db,group_name_sample_db,comp_groups,comps_name_db = ExpressionBuilder.simpleGroupImport(groups_dir,splitHeaders=splitHeaders, ignoreComps=True)
    sample_list = t ### This is the actual order in the input expression files

    for x in t:
        try: group_name = group_db[x]
        except Exception:
            try:
                y = string.split(x,':')[-1] ### for an alternative exon file with the name wt:sample1.bed
                group_name = group_db[y]
            except Exception: pass
        try:
            group_name = group_db[x]
            sample_index = t.index(x)
            try: group_index_db[group_name].append(sample_index)
            except Exception: group_index_db[group_name] = [sample_index] ### dictionary of group to input file sample indexes
        except Exception: pass

    sample_template_db = createReplicateCorrelationTemplate(sample_list,group_index_db)
    return sample_template_db,group_sample_db
    
def createReplicateCorrelationTemplate(samples,group_index_db):
    ### This create multiple binary indicates, but for replicates as opposed to an individual mean of multiple groups
    sample_template_db={}
    null_template = [0.00]*len(samples)
    for group in group_index_db:
        sample_template = list(null_template) ### creates a modifiable copy
        group_indeces = group_index_db[group]
        for index in group_indeces:
            sample_template[index] = correlationDirection ### make -1 to inverse in silico pattern (anti-correlated)
        sample_template_db[group] = sample_template
    return sample_template_db
            
def selectiveFloats(values):
    float_values=[]
    for i in values:
        try:
            float_values.append(float(i))
        except Exception: float_values.append(None)
    return float_values

def binaryExp(value):
    if value>1:
        return 2
    else:
        return 0
    
def identifyMarkers(filename,cluster_comps,binarize=False):
    """ This function is the real workhorse of markerFinder, which coordinates the correlation analyses and data import """
    
    global tissue_scores; tissue_scores={}; print_interval=2000; print_limit=2000
    fn=filepath(filename); x=0; k=0; probeset_db={}; tissues_with_lowest_expression={}; index_sets = []
    global original_tissue_headers
    global use_replicates
    
    count=0

    import gene_associations; from import_scripts import OBO_import
    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
              
    try: coding_db = ExpressionBuilder.importTranscriptBiotypeAnnotations(species)
    except Exception: coding_db = {}
    if 'ExpressionOutput' in filename:
        use_replicates = False
    else:
        use_replicates = True
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)  #remove endline
        t = string.split(data,'\t')
        if data[0] == '#':
            x = 0
        elif x == 0:
            i=0
            for h in t:
                if 'Definition' in h: lti = i ### last tissue index
                if 'Description' in h: lti = i ### last tissue index
                if 'Select Protein Classes' in h: ct = i
                i+=1
            try: original_tissue_headers = t[1:lti]
            except Exception:
                ### Occurs when analyzing a simple expression file with no annotations
                original_tissue_headers = t[1:]
            if len(cluster_comps) == 0: ### No group clusters to separately analyze present
                tissues = list(original_tissue_headers)
            else:
                if len(cluster_comps)>1: ### 2 groups clusters to compare
                    indexes = cluster_comps[0]
                    tissues = t[indexes[0]:indexes[-1]+1]
                    index_sets = [[indexes[0],indexes[-1]+1]]
                    for cluster in cluster_comps:
                        if cluster>0:
                            indexes = cluster_comps[cluster]
                            tissues += t[indexes[0]:indexes[-1]+1]
                            index_sets.append([indexes[0],indexes[-1]+1])
                            #print tissues
                            print 'being analyzed now'
                else: ### First reference set of tissues looked at
                    for cluster in cluster_comps: ### There is only one here!
                        indexes = cluster_comps[cluster]
                        tissues = t[indexes[0]:indexes[-1]+1]
                        index_sets = [[indexes[0],indexes[-1]+1]]
                        #print tissues;sys.exit()
                        print 'being analyzed only in round 1'
                original_tissue_headers2=[]
                for tissue in original_tissue_headers: ### This is the original full header for all clusters
                    try: cluster,tissue = string.split(tissue,'~')
                    except Exception: pass
                    original_tissue_headers2.append(tissue)
                original_tissue_headers = original_tissue_headers2
            try: annotation_headers = t[lti:]
            except Exception: annotation_headers = []

            if len(cluster_comps) > 0:
                tissues2=[]
                for tissue in tissues:
                    if '~' in tissue:
                        cluster, tissue = string.split(tissue,'~')
                    tissues2.append(tissue)
                tissues = tissues2
                #print tissues, cluster_comps;sys.exit()
            if use_replicates:
                if len(cluster_comps)>0:
                    tissue_template_db,group_sample_db = getReplicateData(fn,tissues)
                else:
                    tissue_template_db,group_sample_db = getReplicateData(fn,original_tissue_headers)
                try: annotation_db = getArrayData(AvgExpDir,[])[1]
                except Exception: pass
            else:
                tissue_template_db = createCorrelationTemplate(tissues)
                
            x = 1
        else: #elif x<500:
            probeset = t[0]; proceed = 'no'; symbol=''; geneID=''
            try:
                lti = len(tissues)+1
                try: description,symbol=annotation_db[probeset][:2] ### See above annotation_db download
                except Exception: symbol = probeset; description = ''
                try:
                    probeset,geneID = string.split(probeset,':')
                    if 'ENS' in probeset:
                        geneID, probeset = probeset,geneID
                        probeset=geneID+':'+probeset
                except Exception:
                    if 'ENS' in probeset:
                        geneID = probeset
                try: symbol = gene_to_symbol[geneID][0]; description = ''
                except Exception: pass
                except Exception: symbol = probeset; description = ''
                try: coding_class = coding_db[probeset][-1]
                except Exception:
                    try:
                        geneID = symbol_to_gene[probeset][0]
                        symbol = probeset
                        coding_class = coding_db[geneID][-1]
                    except Exception:
                        coding_class = 'protein_coding'
            except Exception: pass
            if symbol =='' or symbol == probeset:
                try: coding_class = t[ct]; symbol = t[lti+1]; description = t[lti]
                except Exception: coding_class = 'protein_coding'
            if coding_type == 'protein_coding':
                if coding_type in coding_class:
                    if 'MT-' not in symbol and '.' not in symbol:
                        proceed = 'yes'
            elif coding_type == 'AltExon':
                proceed = 'yes'
            else:
                if 'protein_coding' not in coding_class and 'pseudogene' not in coding_class and len(description)>0:
                    if 'MT-' not in symbol and '.' not in symbol:
                        proceed = 'yes'
            proceed = 'yes' ### Force it to anlayze all genes
            count+=1
            #if coding_class != 'protein_coding':
            #print coding_class, coding_type, proceed, probeset, symbol, species, len(gene_to_symbol),coding_db[probeset];sys.exit()
            #proceed = 'yes'
            if len(coding_class) == 0 or proceed == 'yes':
                if compare_clusters == 'yes':
                    exp_values=[]
                    for (i1,i2) in index_sets:
                        try: exp_values += map(float, t[i1:i2])
                        except Exception: exp_values+=selectiveFloats(t[i1:i2])
                    #print len(exp_values), len(tissues)
                    #print exp_values
                    #print tissues; kill
                else:
                    try: exp_values = map(float, t[1:lti]) ### map allows you to apply the function to all elements in the object
                    except Exception: exp_values=selectiveFloats(t[1:lti])
                    if log_transform:
                        try: exp_values =  map(lambda x: math.log(x,2), exp_values)
                        except Exception:
                            exp_values =  map(lambda x: math.log(x+1,2), exp_values)
                    if binarize:
                        exp_values =  map(lambda x: binaryExp(x), exp_values)
                if analyze_housekeeping == 'yes': ### Only grab these when analyzing all tissues
                    findHousekeepingGenes((probeset,symbol),exp_values)
                elif platform == 'RNASeq': ### Exclude low expression (RPKM) genes
                    if max(exp_values)>RPKM_threshold:
                        PearsonCorrelationAnalysis((probeset,symbol),exp_values,tissue_template_db)
                    else:
                        pass
                        #print max(exp_values), RPKM_threshold;sys.exit()
                else:
                    if 'exp.' in filename:
                        try: PearsonCorrelationAnalysis((probeset,symbol),exp_values,tissue_template_db)
                        except Exception: ### For missing values
                            advancedPearsonCorrelationAnalysis((probeset,symbol),exp_values,tissue_template_db)
                    else:
                        advancedPearsonCorrelationAnalysis((probeset,symbol),exp_values,tissue_template_db)
                x+=1
                if x == print_limit:
                    #if print_limit == 2000: break
                    #print print_limit,'genes analyzed'
                    print '*',
                    print_limit+=print_interval
    
    #print len(Added),len(Queried),len(tissue_scores),count;sys.exit()
    tissue_specific_IDs={}; interim_correlations={}
    gene_specific_rho_values = {}
    tissue_list=[]
    for tissue in tissue_scores:
        tissue_scores[tissue].sort()
        tissue_scores[tissue].reverse()
        ranked_list=[]; ranked_lookup=[]
        if tissue not in tissue_list: tissue_list.append(tissue) ### Keep track of the tissue order
        for ((rho,p),(probeset,symbol)) in tissue_scores[tissue]:
            
            ### Get a matrix of all genes to correlations
            try: gene_specific_rho_values[symbol].append(rho)
            except Exception: gene_specific_rho_values[symbol] = [rho]
            
            if symbol == '': symbol = probeset
            #print tissue, tissue_scores[tissue];sys.exit()
            if symbol not in ranked_list:
                ranked_list.append(symbol); ranked_lookup.append([probeset,symbol,(rho,p)])
        for (probeset,symbol,(rho,p)) in ranked_lookup[:genesToReport]:  ### Here is where we would compare rho values between tissues with the same probesets
            if rho>0.1 and p<0.1:
                if compare_clusters == 'yes':
                    try: tissue_specific_IDs[tissue].append(probeset)
                    except Exception: tissue_specific_IDs[tissue] = [probeset]
                else:
                    try: tissue_specific_IDs[probeset].append(tissue)
                    except Exception: tissue_specific_IDs[probeset] = [tissue]
                try: interim_correlations[tissue].append([probeset,symbol,(rho,p)])
                except Exception: interim_correlations[tissue] = [[probeset,symbol,(rho,p)]]    
        if correlateAllGenes:
            for tissue in tissue_scores:
                for ((rho,p),(probeset,symbol)) in tissue_scores[tissue]:
                    try: all_genes_ranked[probeset,symbol].append([(rho,p),tissue])
                    except Exception:all_genes_ranked[probeset,symbol] = [[(rho,p),tissue]]
    """
    
    for ID in all_genes_ranked:
        ag = all_genes_ranked[ID]
        ag.sort()
        all_genes_ranked[ID] = ag[-1] ### topcorrelated
    """
    #"""
    data = export.ExportFile(string.replace(filename[:-4]+'-all-correlations.txt','exp.','MarkerFinder.'))
    data.write(string.join(tissue_list,'\t')+'\n')
    for gene in gene_specific_rho_values:
        data.write(string.join([gene]+map(str,gene_specific_rho_values[gene]),'\t')+'\n')
    #sys.exit()
    #"""
    #print len(tissue_specific_IDs);sys.exit()
    return tissue_specific_IDs,interim_correlations,annotation_headers,tissues

def exportMarkerGeneProfiles(original_filename,annotations,expression_relative,title_row):
    destination_dir = 'AltDatabase/ensembl/'+species+'/' ### Original default
    destination_dir = export.findParentDir(original_filename)
    if 'AltResults' in original_filename: dataset_type = '_AltExon'
    elif 'FullDatasets' in original_filename: dataset_type = '_AltExon'
    else: dataset_type = ''
    #filename = species+'_'+platform+'_tissue-specific'+dataset_type+'_'+coding_type+'.txt'
    filename = 'MarkerFinder/MarkerGenes.txt'
    try:
        if use_replicates:
            filename = string.replace(filename,'.txt','-ReplicateBased.txt')
        else:
            filename = string.replace(filename,'.txt','-MeanBased.txt')
    except Exception: None
    filename = destination_dir+filename
    filename = string.replace(filename,'ExpressionInput','ExpressionOutput')
    data = export.ExportFile(filename)
    title_row = string.join(title_row,'\t')
    data.write(title_row+'\n')
    for probeset in expression_relative:
        values = string.join([probeset]+annotations[probeset]+expression_relative[probeset],'\t')+'\n'
        data.write(values)
    data.close()
    print '\nexported:',filepath(filename)
    return filepath(filename)

def exportAllGeneCorrelations(filename,allGenesRanked):
    destination_dir = export.findParentDir(filename)
    filename = destination_dir+'MarkerFinder/AllGenes_correlations.txt'
    filename = string.replace(filename,'ExpressionInput','ExpressionOutput')
    try:
        if use_replicates:
            filename = string.replace(filename,'.txt','-ReplicateBased.txt')
        else:
            filename = string.replace(filename,'.txt','-MeanBased.txt')
    except Exception: pass
    data = export.ExportFile(filename)
    title_row = string.join(['UID','Symbol','Pearson rho','Pearson p-value','Cell State'],'\t')
    data.write(title_row+'\n')
    rho_sorted=[]
    for (probeset,symbol) in allGenesRanked:
        try: (rho,p),tissue = allGenesRanked[(probeset,symbol)]
        except Exception:
            ### Applies to tiered analysis
            allGenesRanked[(probeset,symbol)].sort()
            (rho,p),tissue = allGenesRanked[(probeset,symbol)][-1]
        values = string.join([probeset,symbol,str(rho),str(p),tissue],'\t')+'\n'
        rho_sorted.append([(tissue,1.0/rho),values])
    rho_sorted.sort()
    for (x,values) in rho_sorted:
        data.write(values)
    data.close()
    
def exportCorrelations(original_filename,interim_correlations):
    destination_dir = 'AltDatabase/ensembl/'+species+'/'
    destination_dir = export.findParentDir(original_filename)
    if 'AltResults' in original_filename: dataset_type = '_AltExon'
    elif 'FullDatasets' in original_filename: dataset_type = '_AltExon'
    else: dataset_type = ''
    filename = species+'_'+platform+'_tissue-specific_correlations'+dataset_type+'_'+coding_type+'.txt'
    filename = destination_dir+filename
    filename = destination_dir+'MarkerFinder/MarkerGenes_correlations.txt'
    filename = string.replace(filename,'ExpressionInput','ExpressionOutput')
    try:
        if use_replicates:
            filename = string.replace(filename,'.txt','-ReplicateBased.txt')
        else:
            filename = string.replace(filename,'.txt','-MeanBased.txt')
    except Exception: pass
    data = export.ExportFile(filename)
    title_row = string.join(['UID','Symbol','Pearson rho','Pearson p-value','Cell State'],'\t')
    data.write(title_row+'\n')
    for tissue in interim_correlations:
        for key in interim_correlations[tissue]:
            probeset,symbol,rho_p = key
            rho,p = rho_p
            values = string.join([probeset,symbol,str(rho),str(p),tissue],'\t')+'\n'
            data.write(values)
    data.close()
    #print 'exported:',filepath(filename)
############### Second set of methods for extracting out average expression columns from initial RMA data ##########

def verifyFile(filename):
    status = 'not found'
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = 'found';break
    except Exception: status = 'not found'
    return status

def getAverageExpressionValues(filename,platform):
    
    """ This function imports two file sets: (A) The original raw input expression files and groups and (B) the DATASET file with annotations.
    It outputs a new file with annotations and average stats for all groups (some group avgs can be missing from the DATASET file)."""

    ### Get the original expression input file location
    if 'ExpressionInput' in filename:
        exp_input_dir = filename
    else:
        exp_input_dir = string.replace(filename,'ExpressionOutput','ExpressionInput')
        exp_input_dir = string.replace(exp_input_dir,'DATASET-','exp.')
        if verifyFile(exp_input_dir) == 'not found':
            exp_input_dir = string.replace(exp_input_dir,'exp.','') ### file may not have exp.
            if verifyFile(exp_input_dir) == 'not found':
                exp_input_dir = string.replace(exp_input_dir,'ExpressionInput/','') ### file may be in a root dir
        if platform != "3'array":
            exp_input_dir = string.replace(exp_input_dir,'.txt','-steady-state.txt')
    ### Find the DATASET file if this is not a DATASET file dir (can be complicated)
    if 'ExpressionInput' in filename:
        filename = string.replace(filename,'ExpressionInput','ExpressionOutput')
        parent_dir = export.findParentDir(filename)
        file = export.findFilename(filename)
        if 'exp.' in file:
            file = string.replace(file,'exp.','DATASET-')
        else:
            file = 'DATASET-'+file
        filename = parent_dir+'/'+file
    elif 'ExpressionOutput' not in filename:
        parent_dir = export.findParentDir(filename)
        file = export.findFilename(filename)
        if 'exp.' in file:
            file = string.replace(file,'exp.','DATASET-')
        else:
            file = 'DATASET-'+file
        filename = parent_dir+'/ExpressionOutput/'+file
    filename = string.replace(filename,'-steady-state.txt','.txt')
    
    ### Import the DATASET file annotations    
    fn=filepath(filename); x=0; k=0; expression_data={}; annotation_columns={}
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)  #remove endline
        t = string.split(data,'\t')
        if x == 0:
            if 'Constitutive_exons_used' in data:
                array_type = 'exon'
                annotation_columns[1]=['Description']
                annotation_columns[2]=['Symbol']
                annotation_columns[4]=['Constitutive_exons_used']
                annotation_columns[7]=['Select Cellular Compartments']
                annotation_columns[8]=['Select Protein Classes']
                annotation_columns[9]=['Chromosome']
                annotation_columns[10]=['Strand']
                annotation_columns[11]=['Genomic Gene Corrdinates']
                annotation_columns[12]=['GO-Biological Process']
            else:
                array_type = "3'array"
                annotation_columns[1]=['Description']
                annotation_columns[2]=['Symbol']
                annotation_columns[3]=['Ensembl_id']
                annotation_columns[9]=['Pathway_info']
                annotation_columns[11]=['Select Cellular Compartments']
                annotation_columns[12]=['Select Protein Classes']
                
            i=0; columns_to_save={}; title_row=[t[0]]; annotation_row=[]
            for h in t:
                if 'avg-' in h:
                    columns_to_save[i]=[]
                    #title_row.append(string.replace(h,'avg-',''))
                if i in annotation_columns: annotation_row.append(h)
                i+=1
            x=1
            if array_type == "3'array":
                annotation_row2 = [annotation_row[1]]+[annotation_row[0]]+annotation_row[2:]### Switch Description and Symbol columns
                annotation_row = annotation_row2
            title_row = string.join(title_row,'\t')+'\t'
            annotation_headers = annotation_row
            annotation_row = string.join(annotation_row,'\t')+'\n'
            title_row+=annotation_row

        else:
            uid = t[0]; exp_vals=[]; annotations=[]
            i=0
            for val in t:
                try:
                    null=columns_to_save[i]
                    exp_vals.append(val)
                except Exception: null=[]
                try:
                    null = annotation_columns[i]
                    annotations.append(val)
                except Exception: null=[]
                i+=1
            if array_type == "3'array":
                annotations2 = [annotations[1]]+[annotations[0]]+annotations[2:]### Switch Description and Symbol columns
                annotations = annotations2
            expression_data[uid]=annotations # exp_vals+annotations
                
    ### This function actually takes the average of all values - not biased by groups indicated in comps.
    filename = string.replace(filename,'DATASET-','AVERAGE-')
    #filename = exportSimple(filename,expression_data,title_row)
    importAndAverageExport(exp_input_dir,platform,annotationDB=expression_data,annotationHeader=annotation_headers,customExportPath=filename)
    return filename

def getAverageExonExpression(species,platform,input_exp_file):

    ### Determine probesets with good evidence of expression
    global max_exp_exon_db; max_exp_exon_db={}; global expressed_exon_db; expressed_exon_db={}
    global alternative_exon_db; alternative_exon_db={}; global alternative_annotations; alternative_annotations={}
    global exon_db
    importRawExonExpData(input_exp_file)
    importDABGData(input_exp_file)
    importAltExonData(input_exp_file)
    
    includeOnlyKnownAltExons = 'yes'

    if includeOnlyKnownAltExons == 'yes':
        ### Optionally only include exons with existing alternative splicing annotations    
        import AltAnalyze
        if platform == 'exon' or platform == 'gene': probeset_type = 'full'
        else: probeset_type = 'all'
        avg_all_for_ss = 'no'; root_dir = ''; filter_by_known_AE = 'no' ### Exclude exons that are not known to be alternatively expressed from cDNA databases
        exon_db,constitutive_probeset_db = AltAnalyze.importSplicingAnnotations(platform,species,probeset_type,avg_all_for_ss,root_dir)
        del constitutive_probeset_db
        protein_coding_db = ExpressionBuilder.importTranscriptBiotypeAnnotations(species)
        
        delete_exon_entries={}; x=0
        for probeset in exon_db:
            try: null = alternative_exon_db[probeset]
            except Exception: delete_exon_entries[probeset]=[]
            ### If a splicing annotation exists
            ed = exon_db[probeset]
            as_call = ed.SplicingCall()
            gene = ed.GeneID()
            compartment, custom_class = protein_coding_db[gene]
            if 'protein_coding' not in custom_class: delete_exon_entries[probeset]=[]
            if filter_by_known_AE == 'yes':
                if as_call == 0: delete_exon_entries[probeset]=[]
            x+=1
        ### Delete where not expressed, alternatively expressed or known to be an alternative exon
        print len(exon_db)-len(delete_exon_entries), 'out of', len(exon_db), 'with known alternative exon annotations analyzed'
        for probeset in delete_exon_entries: del exon_db[probeset] ### Clear objects from memory

    ### Get gene-level annotations
    gene_annotation_file = "AltDatabase/ensembl/"+species+"/"+species+"_Ensembl-annotations.txt"
    from build_scripts import ExonAnalyze_module; global annotate_db
    annotate_db = ExonAnalyze_module.import_annotations(gene_annotation_file,platform)
    
    importRawSpliceData(input_exp_file)
    AltAnalyze.clearObjectsFromMemory(exon_db)

def importAltExonData(filename):
    file = export.findFilename(filename)
    if 'AltResults' in filename:
        root_dir = string.split(filename,'AltResults')[0]
    elif 'AltExpression' in filename:
        root_dir = string.split(filename,'AltExpression')[0]
    alt_exon_dir = root_dir+'AltResults/AlternativeOutput/'
    dir_list = read_directory(alt_exon_dir)
    #print [file];kill
    methods_assessed=0
    for alt_exon_file in dir_list:
        if file[:-4] in alt_exon_file and '-exon-inclusion-results.txt' in alt_exon_file:
            alt_exon_file = alt_exon_dir+alt_exon_file
            fn = filepath(alt_exon_file); x=0
            print 'AltExonResults:',fn
            methods_assessed+=1
            for line in open(fn,'rU').xreadlines():             
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if x == 0: x=1
                else:
                    probeset = t[6]; isoform_description=t[14]; domain_description=t[15]; indirect_domain=t[17]; x+=1
                    #if probeset == '2681775': print '2681775 in AltResults file'
                    exon_annotation = t[4]; block_exon_annotation = t[26]; ens_exon = t[24]; genomic_location = t[-1]
                    splicing_annotation = t[27]
                    if len(block_exon_annotation)>0: exon_annotation = block_exon_annotation
                    if probeset in expressed_exon_db: ### If the probeset is expressed and in the alternative exon file
                        try: alternative_exon_db[probeset] += 1
                        except Exception: alternative_exon_db[probeset] = 1
                        alternative_annotations[probeset] = genomic_location,ens_exon,exon_annotation,splicing_annotation,isoform_description,domain_description,indirect_domain
    print len(alternative_exon_db), 'expressed exons considered alternative out of',x

    if methods_assessed > 1: ### Likely both FIRMA and splicing-index
        ### only include those that are present with both methods
        single_algorithm  = {}
        for probeset in alternative_exon_db:
            if alternative_exon_db[probeset]==1:
                single_algorithm[probeset] = []
        for probeset in single_algorithm:
            del alternative_exon_db[probeset]
            del alternative_annotations[probeset]
            
def importDABGData(filename,filterIndex=False,DABGFile=True):
    try:
        file = export.findFilename(filename)
        if 'AltResults' in filename:
            root_dir = string.split(filename,'AltResults')[0]
        elif 'AltExpression' in filename:
            root_dir = string.split(filename,'AltExpression')[0]
        dabg_file = root_dir+'ExpressionInput/stats.'+file  
    except Exception:
        dabg_file = filename ### when supplying this file directly

    expressed_exons={}
    fn = filepath(dabg_file); x=0
    print 'DABGInput:',fn
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x == 0:
            if '#' != data[0]:
                header = t
                x=1
        else:
            if DABGFile:
                dabg_values = map(float, t[1:])
                min_dabg = min(dabg_values)
                if min_dabg<0.01:
                    proceed = True
                    try:
                        if t[0] not in max_exp_exon_db: proceed = False
                    except Exception: pass
                    #if t[0] == '2681775': print '2681775 in DABG file'
                    if filterIndex and proceed:
                        filtered_index=[]
                        i=0
                        for dabg in dabg_values:
                            if dabg<0.01: ### Expressed samples for that gene 
                                filtered_index.append(i)
                            i+=1
                        if len(filtered_index)>5:
                            expressed_exon_db[t[0]] = filtered_index
                    elif proceed:
                        expressed_exon_db[t[0]] = []
            else:
                #{4: 29460, 5: 150826, 6: 249487, 7: 278714, 8: 244304, 9: 187167, 10: 135514, 11: 84828, 12: 39731, 13: 10834, 14: 500, 15: 34}
                max_val = max(map(float, t[1:]))
                if max_val > 10:
                    expressed_exons[t[0]] = []

    if 'steady-state' in filename: type = 'genes'
    else: type = 'exons'
    print len(expressed_exon_db), type, 'meeting expression and dabg thesholds'
    return expressed_exons

def filterAltExonResults(fn,filterDB=None):
    ### take a stats file as input
    fn = string.replace(fn,'-steady-state.txt','-splicing-index-exon-inclusion-results.txt')
    fn = string.replace(fn,'-steady-state-average.txt','-splicing-index-exon-inclusion-results.txt')
    fn = string.replace(fn, 'ExpressionInput/stats.','AltResults/AlternativeOutput/')
    firstRow=True
    regulatedProbesets = {}
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            header = t
            ea = header.index('exon annotations')
            fp = header.index('functional_prediction')
            ps = header.index('probeset')
            up = header.index('uniprot-ens_feature_predictions')
            sy = header.index('symbol')
            de = header.index('description')
            ex = header.index('exons')
            pl = header.index('probeset location')
            ee = header.index('ensembl exons')
            eo = header.index('ens_overlapping_domains')
            en = header.index('Ensembl')
            firstRow=False
            
        else:
            if filterDB!=None:
                probeset = t[ps]
                #UID	Definition	Symbol	Ensembl	Genomic Location	ExonExternalID	Exon ID	Alternative Exon Annotation	Isoform Associations	Inferred Domain Modification	Direct Domain Modification
                if probeset in filterDB:
                    markersIn = filterDB[probeset]
                    markersIn = string.join(markersIn,'|')
                    values = string.join([probeset,t[de],t[sy],t[en],t[pl],t[ee],t[ex],t[ea],t[fp],t[up],t[eo],markersIn],'\t')
                    regulatedProbesets[probeset] = values
            elif len(t[ea])>2: # and len(t[fp])>1 and t[up]>1:
                if 'cassette' in t[ea] or 'alt-3' in t[ea] or 'alt-5' in t[ea] or 'alt-C-term' in t[ea] or 'exon' in t[ea] or 'intron' in t[ea] or 'altFive' in t[ea] or 'altThree' in t[ea]:
                    if 'altPromoter' not in t[ea] and 'alt-N-term' not in t[ea]:
                        regulatedProbesets[t[ps]]=t[up]
                #regulatedProbesets[t[ps]]=t[up]
            if filterDB==None: pass #regulatedProbesets[t[ps]]=t[up]
            evaluatedGenes[t[en]]=[]
        
    """
    print len(regulatedProbesets)          
    if filterDB==None:
        import AltAnalyze
        if platform == 'exon' or platform == 'gene': probeset_type = 'full'
        else: probeset_type = 'all'
        avg_all_for_ss = 'no'; root_dir = ''; filter_by_known_AE = 'no' ### Exclude exons that are not known to be alternatively expressed from cDNA databases
        exon_db,constitutive_probeset_db = AltAnalyze.importSplicingAnnotations('exon',species,probeset_type,avg_all_for_ss,root_dir)
        protein_coding_db = ExpressionBuilder.importTranscriptBiotypeAnnotations(species)
        for probeset in exon_db:
            if probeset not in constitutive_probeset_db:
                ### If a splicing annotation exists
                ed = exon_db[probeset]
                as_call = ed.SplicingCall()
                gene = ed.GeneID()
                compartment, custom_class = protein_coding_db[gene]
                if 'protein_coding' in custom_class: #as_call != 0: #
                    regulatedProbesets[probeset]=''
                evaluatedGenes[gene]=[]
        print len(regulatedProbesets), 'AltExon probesets retained...'
    """
    return regulatedProbesets

def averageNIValues(fn,dabg_gene_dir,regulatedProbesets):
    export_data = export.ExportFile(fn)
    fn = string.replace(fn, '-average','')
    groups_dir = string.replace(dabg_gene_dir,'exp.','groups.')
    groups_dir = string.replace(groups_dir,'stats.','groups.')
    groups_dir = string.replace(groups_dir,'-steady-state','')
    groups_dir = string.replace(groups_dir,'-average','')
    firstRow=True
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            headers=[]
            for h in t:
                try: h = string.split(h,':')[-1]
                except Exception: pass
                headers.append(h)
            
            group_index_db={}
            ### use comps in the future to visualize group comparison changes
            sample_list,group_sample_db,group_db,group_name_sample_db,comp_groups,comps_name_db = ExpressionBuilder.simpleGroupImport(groups_dir)
            for x in sample_list:
                group_name = group_db[x]
                sample_index = headers.index(x)
                try: group_index_db[group_name].append(sample_index)
                except Exception: group_index_db[group_name] = [sample_index] ### dictionary of group to input file sample indexes
            groups = map(str, group_index_db) ### store group names
            new_sample_list = map(lambda x: group_db[x], sample_list) ### lookup index of each sample in the ordered group sample list  
            groups.sort()
            export_data.write(string.join(t[:3]+groups,'\t')+'\n')
            firstRow=False

        else:
            geneID = t[0]
            probeset = t[2]
            if probeset in regulatedProbesets:
                avg_z=[]
                for group_name in groups:
                    group_values = map(lambda x: float(t[x]), group_index_db[group_name]) ### simple and fast way to reorganize the samples
                    avg = statistics.avg(group_values)
                    avg_z.append(str(avg))
                values = string.join(t[:3]+avg_z,'\t')+'\n'
                export_data.write(values)
    export_data.close()
    
def calculateSplicingIndexForExpressedConditions(species,dabg_gene_dir,regulatedProbesets,genes_to_exclude):
    protein_coding_db = ExpressionBuilder.importTranscriptBiotypeAnnotations(species)
    fn = string.replace(dabg_gene_dir, 'ExpressionInput/stats.','AltResults/RawSpliceData/'+species+'/splicing-index/')
    fn = string.replace(fn, '-steady-state','')
    
    if '-average.txt' in fn:
        averageNIValues(fn,dabg_gene_dir,regulatedProbesets)
    firstRow=True
    onlyIncludeOldExons=False
    splicedExons = {}
    splicedExons1={}
    splicedExons2={}
    splicedExons3={}
    splicedExons4={}
    splicedExons5={}
    splicedExons6={}
    splicedExons7={}
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            header = t[3:]
            ### Test to below to make sure we have the right headers - it works - only heart for a heart specific gene
            #exp_indexes = expressed_exon_db['ENSG00000134571']
            #filtered_headers = map(lambda i: header[i],exp_indexes)
            #print filtered_headers;sys.exit()
            headers=[]
            for h in header:
                h = string.split(h,':')[-1]
                headers.append(h)
            headers = string.join(headers,'\t')
            firstRow = False
        else:
            geneID = t[0]
            probeset = t[2]
            compartment, custom_class = protein_coding_db[geneID]
            if 'protein_coding' in custom_class and onlyIncludeOldExons == False: ### Only consider protein coding genes
                if probeset in regulatedProbesets and geneID in expressed_exon_db and probeset in expressed_exon_db: ### If the probeset is associated with a protein modifying splicing event and has gene expressed conditions
                    exp_indexes = expressed_exon_db[geneID]
                    NI_values = map(float, t[3:])
                    exp_NI_values = map(lambda i: NI_values[i],exp_indexes)
                    #print len(exp_NI_values), len(exp_indexes), len(NI_values)
                    si = max(exp_NI_values)-statistics.avg(exp_NI_values)
                    if si>2 and max(NI_values)>-1: ### Ensures that the splicing event is large and highly expressed
                        i=0; mod_NI_values=[]
                        for ni in NI_values: ### Replaced non-expressed samples with empty values
                            if i in exp_indexes:
                                mod_NI_values.append(str(ni))
                            else:
                                mod_NI_values.append('')
                            i+=1
                        mod_NI_values = string.join(mod_NI_values,'\t')
                        if geneID not in genes_to_exclude: ### gene markers
                            splicedExons[probeset,geneID] = mod_NI_values
                            splicedExons6[probeset] = []
                    splicedExons7[probeset] = []
            if onlyIncludeOldExons:
                if probeset in prior_altExons and geneID in expressed_exon_db:
                    exp_indexes = expressed_exon_db[geneID]
                    i=0; NI_values = map(float, t[3:]); mod_NI_values=[]
                    for ni in NI_values: ### Replaced non-expressed samples with empty values
                        if i in exp_indexes:
                            mod_NI_values.append(str(ni))
                        else:
                            mod_NI_values.append('')
                        i+=1
                    mod_NI_values = string.join(mod_NI_values,'\t')
                    splicedExons[probeset,geneID] = mod_NI_values
            if probeset in regulatedProbesets:
                splicedExons1[probeset]=[]
            if probeset in regulatedProbesets and 'protein_coding' in custom_class:
                splicedExons2[probeset]=[]
            if probeset in regulatedProbesets and probeset in expressed_exon_db:
                splicedExons3[probeset]=[]
            if probeset in regulatedProbesets and geneID in expressed_exon_db:
                splicedExons4[probeset]=[]
            if probeset in regulatedProbesets and geneID not in genes_to_exclude:
                splicedExons5[probeset]=[]


    print len(splicedExons),'top alternative exons...'
    output_file = string.replace(dabg_gene_dir,'stats.','topSplice.')
    eo = export.ExportFile(output_file)
    eo.write('Probeset:GeneID\t'+headers+'\n')
    for (probeset,geneID) in splicedExons:
        mod_NI_values = splicedExons[(probeset,geneID)]
        eo.write(probeset+':'+geneID+'\t'+mod_NI_values+'\n')
    eo.close()
    
    expressed_uids1 = splicedExons1.viewkeys() & prior_altExons.viewkeys()
    print len(expressed_uids1)
    
    expressed_uids1 = splicedExons1.viewkeys() & prior_altExons.viewkeys()
    print len(expressed_uids1)

    expressed_uids1 = splicedExons2.viewkeys() & prior_altExons.viewkeys()
    print len(expressed_uids1)
  
    expressed_uids1 = splicedExons3.viewkeys() & prior_altExons.viewkeys()
    print len(expressed_uids1)
    
    expressed_uids1 = splicedExons4.viewkeys() & prior_altExons.viewkeys()
    print len(expressed_uids1)
    
    expressed_uids1 = splicedExons5.viewkeys() & prior_altExons.viewkeys()
    print len(expressed_uids1)
    
    expressed_uids1 = splicedExons6.viewkeys() & prior_altExons.viewkeys()
    print len(expressed_uids1)
    
    expressed_uids1 = splicedExons7.viewkeys() & prior_altExons.viewkeys()
    print len(expressed_uids1)

    return output_file, splicedExons
        
def filterRNASeqSpliceEvents(Species,Platform,fl,psi_file_dir):
    global species; import AltAnalyze
    global platform; from import_scripts import sampleIndexSelection
    global max_exp_exon_db
    global evaluatedGenes; evaluatedGenes={'Gene_ID':[]}
    species = Species
    platform = Platform
    output_file = psi_file_dir
    compendiumType = 'AltExon'
    #analyzeData(output_file,species,'altSplice',compendiumType,geneToReport=200,AdditionalParameters=fl,logTransform=False,correlateAll=False)
    output_file = '/Volumes/SEQ-DATA/IlluminaBodyMap/exons/SplicingMarkers/exp.BodyMap.txt'
    compendiumType = 'AltExon'
    analyzeData(output_file,species,'RNASeq',compendiumType,geneToReport=200,AdditionalParameters=fl,logTransform=False,correlateAll=False)

def filterDetectionPvalues(Species,Platform,fl,dabg_gene_dir):
    global species; import AltAnalyze
    global platform; from import_scripts import sampleIndexSelection
    global max_exp_exon_db; global prior_altExons
    global evaluatedGenes; evaluatedGenes={'Gene_ID':[]}
    global prior_altExons
    species = Species
    platform = Platform
    averageProfiles = False
    includeMarkerGeneExons = False
    alsoConsiderExpression = False
    filterForHighExpressionExons = True
    
    input_exp_file = string.replace(dabg_gene_dir,'stats.','exp.')
    input_exp_file = string.replace(input_exp_file,'-steady-state','')

    global expressed_exon_db; expressed_exon_db={}
    #"""
    prior_altExons = importMarkerFinderExons(dabg_gene_dir,species=species,type='MarkerFinder-Exon')
    #reorderMultiLevelExpressionFile(dabg_gene_dir) ### Re-sort stats file
    regulatedProbesets = filterAltExonResults(dabg_gene_dir) ### AltResults/AlternativeOutput SI file
    if filterForHighExpressionExons:
        high_expression_probesets = importDABGData(input_exp_file,filterIndex=True,DABGFile=False) ### Get highly expressed exons by RMA
        regulatedProbesets2={}
        print len(regulatedProbesets)
        for i in high_expression_probesets:
            if i in regulatedProbesets:
                regulatedProbesets2[i]=regulatedProbesets[i]
        regulatedProbesets=regulatedProbesets2
    print len(regulatedProbesets)
        
    dabg_exon_dir = string.replace(dabg_gene_dir,'-steady-state','')

    if averageProfiles:
        regulatedProbesets['probeset_id']=[] ### Will force to get the header row
        filtered_exon_exp = AltAnalyze.importGenericFiltered(input_exp_file,regulatedProbesets)
        input_exp_file = filterAverageExonProbesets(input_exp_file,filtered_exon_exp,regulatedProbesets,exportType='expression')
        filtered_gene_dabg = AltAnalyze.importGenericFiltered(dabg_gene_dir,evaluatedGenes)
        dabg_gene_dir = filterAverageExonProbesets(dabg_gene_dir,filtered_gene_dabg,evaluatedGenes,exportType='expression')
        filtered_exon_exp = AltAnalyze.importGenericFiltered(dabg_exon_dir,regulatedProbesets)
        dabg_exon_dir = filterAverageExonProbesets(dabg_exon_dir,filtered_exon_exp,regulatedProbesets,exportType='expression')
        reorderMultiLevelExpressionFile(input_exp_file)
        reorderMultiLevelExpressionFile(dabg_gene_dir)
        reorderMultiLevelExpressionFile(dabg_exon_dir)
        input_exp_file = string.replace(input_exp_file,'.txt','-average.txt')
        dabg_gene_dir = string.replace(dabg_gene_dir,'.txt','-average.txt')
        dabg_exon_dir = string.replace(dabg_exon_dir,'.txt','-average.txt')
    
    expressed_uids1 = regulatedProbesets.viewkeys() & prior_altExons.viewkeys()

    importDABGData(dabg_gene_dir,filterIndex=True) ### Get expressed conditions by DABG p-value for genes
    max_exp_exon_db = regulatedProbesets
    importDABGData(dabg_exon_dir,filterIndex=False) ### Get expressed conditions by DABG p-value for selected exons
    
    genes_to_exclude = importMarkerFinderExons(dabg_gene_dir,species=species,type='MarkerFinder-Gene')

    output_file,splicedExons = calculateSplicingIndexForExpressedConditions(species,dabg_gene_dir,regulatedProbesets,genes_to_exclude)
   
    output_file = string.replace(dabg_gene_dir,'stats.','topSplice.')
    array_type = 'altSplice'; compendiumType = 'AltExon'
 
    analyzeData(output_file,species,array_type,compendiumType,geneToReport=200,AdditionalParameters=fl,logTransform=False,correlateAll=False)
    
    addSuffixToMarkerFinderFile(dabg_gene_dir,'AltExon')
  
    #"""  
    splicedExons = importMarkerFinderExons(dabg_gene_dir,type='AltExon')
    
    if alsoConsiderExpression:
        exp_file = string.replace(dabg_gene_dir,'stats.','exp.')
        exp_file = string.replace(exp_file,'-steady-state','')
        output_file = string.replace(exp_file,'exp.','filter.')
        sampleIndexSelection.filterRows(exp_file,output_file,filterDB=splicedExons)
        analyzeData(output_file,species,platform,'AltExon',geneToReport=400,AdditionalParameters=fl,logTransform=False,correlateAll=False)
        
        addSuffixToMarkerFinderFile(dabg_gene_dir,'Filtered')
        filteredExons = importMarkerFinderExons(dabg_gene_dir,type='Filtered')
        splicedExons = intersectMarkerFinderExons(splicedExons,filteredExons)

    regulatedProbesets = filterAltExonResults(dabg_gene_dir,filterDB=splicedExons)
    ### Make a probeset-level average expression file
    
    splicedExons['probeset_id']=[] ### Will force to get the header row

    if includeMarkerGeneExons:
        markerGenes = importMarkerFinderExons(dabg_gene_dir,species=species,type='MarkerFinder-Gene')
        splicedExons, regulatedProbesets = getConstitutiveMarkerExons(species, platform, markerGenes, splicedExons, regulatedProbesets)
    filtered_exon_exp = AltAnalyze.importGenericFiltered(input_exp_file,splicedExons)
    filterAverageExonProbesets(input_exp_file,filtered_exon_exp,regulatedProbesets,exportType='MarkerFinder');sys.exit()

def getConstitutiveMarkerExons(species, platform, markerGenes, splicedExons, regulatedProbesets):
    #Identify a single consitutive probeset for prior identified MarkerGenes
    genesAdded={}; count=0
    import AltAnalyze
    if platform == 'exon' or platform == 'gene': probeset_type = 'full'
    else: probeset_type = 'all'
    avg_all_for_ss = 'no'; root_dir = ''; filter_by_known_AE = 'no' ### Exclude exons that are not known to be alternatively expressed from cDNA databases
    exon_db,constitutive_probeset_db = AltAnalyze.importSplicingAnnotations(platform,species,probeset_type,avg_all_for_ss,root_dir)
    filler = ['']*7
    import gene_associations
    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    
    for probeset in constitutive_probeset_db:
        geneID = constitutive_probeset_db[probeset]
        if geneID in markerGenes and geneID not in genesAdded:
            count+=1
            splicedExons[probeset]=markerGenes[geneID]
            tissues = string.join(markerGenes[geneID],'|')
            try: symbol = gene_to_symbol[geneID][0]; description = ''
            except Exception: symbol = ''; description = ''
            regulatedProbesets[probeset]=string.join([probeset,description,symbol,geneID]+filler+[tissues],'\t')
            genesAdded[geneID]=[]
    print count, 'constitutive marker gene exons added'
    return splicedExons, regulatedProbesets

def intersectMarkerFinderExons(splicedExons,filteredExons):
    splicedExons_filtered={}
    for exon in splicedExons:
        tissues = splicedExons[exon]
        filtered_tissues=[]
        if exon in filteredExons:
            tissues2 = filteredExons[exon]
            for tissue in tissues:
                if tissue in tissues2:
                    filtered_tissues.append(tissue)
        if len(filtered_tissues)>0:
            splicedExons_filtered[exon] = filtered_tissues
    print len(splicedExons_filtered), 'exons found analyzing both NI values and raw expression values'
    return splicedExons_filtered

def addSuffixToMarkerFinderFile(expr_input,suffix):
    markerfinder_results = string.split(expr_input,'ExpressionInput')[0]+'ExpressionOutput/MarkerFinder/MarkerGenes_correlations-ReplicateBased.txt'
    new_file = string.replace(markerfinder_results,'.txt','-'+suffix+'.txt')
    export.copyFile(markerfinder_results, new_file)
    
def filterAverageExonProbesets(expr_input,filtered_exon_exp,regulatedProbesets,exportType='expression'):
    if exportType == 'expression':
        export_path = string.replace(expr_input,'.txt','-average.txt')
    else:
        export_path = string.split(expr_input,'ExpressionInput')[0]+'ExpressionOutput/MarkerFinder/MarkerGenes-ReplicateBased-AltExon.txt'
    export_data = export.ExportFile(export_path)
    groups_dir = string.replace(expr_input,'exp.','groups.')
    groups_dir = string.replace(groups_dir,'stats.','groups.')
    groups_dir = string.replace(groups_dir,'-steady-state','')
    
    group_index_db={}
    ### use comps in the future to visualize group comparison changes
    sample_list,group_sample_db,group_db,group_name_sample_db,comp_groups,comps_name_db = ExpressionBuilder.simpleGroupImport(groups_dir)
    for x in sample_list:
        group_name = group_db[x]
        try: sample_index = filtered_exon_exp['probeset_id'].index(x)
        except Exception: sample_index = filtered_exon_exp['Gene_ID'].index(x)
        try: group_index_db[group_name].append(sample_index)
        except Exception: group_index_db[group_name] = [sample_index] ### dictionary of group to input file sample indexes
    groups = map(str, group_index_db) ### store group names
    new_sample_list = map(lambda x: group_db[x], sample_list) ### lookup index of each sample in the ordered group sample list
    
    if exportType == 'expression':
        headers = 'probeset_id\t'
    else:
        headers = string.join(['UID','Definition','Symbol','Ensembl','Genomic Location','ExonExternalID','Exon ID','Alternative Exon Annotation','Isoform Associations','Inferred Domain Modification','Direct Domain Modification', 'marker-in'],'\t')+'\t'
    headers += string.join(groups,'\t')
    export_data.write(headers+'\n')

    
    for uid in filtered_exon_exp:
        if uid != 'probeset_id' and uid != 'Gene_ID':
            values = filtered_exon_exp[uid]
            if uid in regulatedProbesets:
                annotations = regulatedProbesets[uid]
                if platform == 'RNASeq':
                    ### Convert to log2 RPKM values - or counts
                    values = map(lambda x: math.log(float(x),2), values)
                else:
                    values = map(float,values)
                avg_z=[]
                for group_name in group_index_db:
                    group_values = map(lambda x: values[x], group_index_db[group_name]) ### simple and fast way to reorganize the samples
                    avg = statistics.avg(group_values)
                    avg_z.append(str(avg))
                if exportType == 'expression':
                    values = string.join([uid]+avg_z,'\t')+'\n'
                else:
                    values = annotations+'\t'+string.join(avg_z,'\t')+'\n'
                export_data.write(values)
    export_data.close()
    return export_path
    
def importMarkerFinderExons(dabg_gene_dir,species=None,type=None):
    if type=='AltExon':
        fn = string.split(dabg_gene_dir,'ExpressionInput')[0]+'ExpressionOutput/MarkerFinder/MarkerGenes_correlations-ReplicateBased-AltExon.txt'
    elif type=='Filtered':
        fn = string.split(dabg_gene_dir,'ExpressionInput')[0]+'ExpressionOutput/MarkerFinder/MarkerGenes_correlations-ReplicateBased-Filtered.txt'
    elif type == 'MarkerFinder-Gene':
        coding_type = 'protein_coding'
        fn = 'AltDatabase/ensembl/'+species+'/'+species+'_'+platform +'_tissue-specific_'+coding_type+'.txt'
        fn = filepath(fn)
    elif type == 'MarkerFinder-Exon':
        coding_type = 'protein_coding'
        fn = 'AltDatabase/ensembl/'+species+'/'+species+'_'+platform +'_tissue-specific_AltExon_'+coding_type+'.txt'
        fn = filepath(fn)
        
    print 'Importing',fn
    firstRow=True
    splicedExons = {}
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            header = t
            firstRow=False
        else:
            try: splicedExons[t[0]].append(t[-1])
            except Exception: splicedExons[t[0]] = [t[-1]]
    print 'Stored', len(splicedExons), 'markers for further evaluation.'
    return splicedExons
            
def importRawExonExpData(filename):
    file = export.findFilename(filename)
    if 'AltResults' in filename:
        root_dir = string.split(filename,'AltResults')[0]
    elif 'AltExpression' in filename:
        root_dir = string.split(filename,'AltExpression')[0]
    exp_file = root_dir+'ExpressionInput/exp.'+file
    fn = filepath(exp_file); x=0
    log_exp_threshold = math.log(70,2)
    print 'RawExpInput:',fn
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x == 0:
            if '#' != data[0]: x=1
        else:
            max_exp = max(map(float, t[1:]))
            if max_exp>log_exp_threshold:
                max_exp_exon_db[t[0]] = []
                #if t[0] == '2681775': print '2681775 in expression file'
    print len(max_exp_exon_db), 'exons with maximal expression greater than 50'
    
def importRawSpliceData(filename):
    """Import the RawSplice file normalized intensity exon data (or Raw Expression - optional) and then export after filtering, averaging and annotating for each group"""
    fn=filepath(filename); x=0; group_db={}
    print 'RawExonInput:',fn
    
    if 'RawSplice' in fn:    
        output_file = string.replace(fn,'RawSplice','AVERAGESplice')
    elif 'FullDatasets' in fn:    
        output_file = string.replace(fn,'FullDatasets','AVERAGE-FullDatasets')
    else: print 'WARNING! The text "RawSplice" must be in the input filename to perform this analysis'; sys.exit()
    export_data = export.ExportFile(output_file) ### Write out the averaged data once read and processed
    
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x == 0:
            i=0; group_list=[]; group_list2=[]; headers = ['ExonID']
            for h in t:
                try:
                    ### Store the index for each group numerical ID
                    group_name,sample = string.split(h,':')
                    try: group_db[group_name].append(i)
                    except KeyError: group_db[group_name] = [i]
                except Exception: null=[]
                i+=1
                
            for group_name in group_db:
                group_list.append(group_db[group_name]) ### list of group indexes in order
                group_list2.append([group_db[group_name],group_name]) 
            group_list.sort(); group_list2.sort()

            for values in group_list2: headers.append(values[-1]) ### Store the original order of group names
            headers = string.join(headers+['Definition','Symbol','Ensembl','Genomic Location','ExonExternalID','Exon ID','Alternative Exon Annotation','Isoform Associations','Inferred Domain Modification','Direct Domain Modification'],'\t')+'\n'
            export_data.write(headers)
            x+=1
        else:
            avg_NI=[]
            if '-' in t[0]: exonid = string.join(string.split(t[0],'-')[1:],'-') ### For the NI data
            else: exonid = t[0] ### When importing the FullDataset organized expression data
            try:
                    #if exonid == '2681775': print '2681775 is in the RawSplice file'
                    ed = exon_db[exonid]; gene = ed.GeneID()
                    #if exonid == '2681775': print '2681775 is in the filtered exon_db'
                    y = annotate_db[gene]; symbol = y.Symbol(); description = y.Description()
                    genomic_location,ens_exon,exon_annotation,splicing_annotation,isoform_description,domain_description,indirect_domain = alternative_annotations[exonid]
                    for indexes in group_list:
                        avg_NI.append(avg(map(float, t[indexes[0]:indexes[-1]+1])))
                    values = string.join([exonid]+map(str, avg_NI)+[description,symbol,gene,genomic_location,ens_exon,exon_annotation,splicing_annotation,isoform_description,domain_description,indirect_domain],'\t')+'\n'
                    export_data.write(values)
            except Exception: null=[]

def getExprValsForNICorrelations(array_type,altexon_correlation_file,rawsplice_file):
    """Import the FullDatasets expression file to replace NI values from the built tissue correlation file"""

    ### Import the AltExon correlation file
    fn = filepath(altexon_correlation_file); marker_probesets={}; x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            headers = t; x=1; index=0
            for i in headers:
                if 'marker-in' in i: tissue_index = index+1
                index+=1
            annotation_headers = headers[:tissue_index]
        else:
            marker_probeset = t[0]
            marker_probesets[marker_probeset] = t[:tissue_index] ### Store annotations

    print len(marker_probesets), 'marker probests imported'
    ### Import the corresponding full-dataset expression file (must average expression for replicates)          
    if array_type == 'exon': array_type_dir = 'ExonArray'
    elif array_type == 'gene': array_type_dir = 'GeneArray'
    elif array_type == 'junction': array_type_dir = 'JunctionArray'
    else: array_type_dir = array_type
            
    file = export.findFilename(rawsplice_file); x=0; group_db={}
    root_dir = string.split(rawsplice_file,'AltResults')[0]

    output_file = string.replace(altexon_correlation_file,'.txt','-exp.txt')
    export_data = export.ExportFile(output_file) ### Write out the averaged data once read and processed        
    exp_file = root_dir+'AltExpression/FullDatasets/'+array_type_dir+'/'+species+'/'+file
    fn = filepath(exp_file)
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x == 0:
            i=0; group_list=[]; group_list2=[]; headers = []
            for h in t:
                try:
                    ### Store the index for each group numerical ID
                    group_name,sample = string.split(h,':')
                    if '~' in group_name:
                        group_name = string.split(group_name,'~')[1]
                    try: group_db[group_name].append(i)
                    except KeyError: group_db[group_name] = [i]
                except Exception: null=[]
                i+=1
                
            for group_name in group_db:
                group_list.append(group_db[group_name]) ### list of group indexes in order
                group_list2.append([group_db[group_name],group_name]) 
            group_list.sort(); group_list2.sort()

            for values in group_list2: headers.append(values[-1]) ### Store the original order of group names
            headers = string.join(annotation_headers+headers,'\t')+'\n'
            export_data.write(headers)
            x+=1
        else:
            exonid = t[0]; avg_exp=[]
            try:
                annotations = marker_probesets[exonid]
                for indexes in group_list:
                    avg_exp.append(avg(map(float, t[indexes[0]:indexes[-1]+1])))
                    values = string.join(annotations+map(str, avg_exp),'\t')+'\n'
                export_data.write(values)
            except Exception: null=[]
            
def avg(ls):
    return sum(ls)/len(ls)

def exportSimple(filename,expression_data,title_row):
    filename = string.replace(filename,'DATASET-','AVERAGE-')
    data = export.ExportFile(filename)
    data.write(title_row)
    for uid in expression_data:
        values = string.join([uid]+expression_data[uid],'\t')+'\n'
        data.write(values)
    data.close()
    #print 'exported...'
    #print filename
    return filename

def returnCommonProfiles(species):
    ###Looks at exon and gene array AltExon predictions to see which are in common

    targetPlatforms = ['exon','gene']; tissue_to_gene={}; rho_threshold = 0; p_threshold = 0.2
    import TissueProfiler
    gene_translation_db = TissueProfiler.remoteImportExonIDTranslations('gene',species,'no','exon')
    for targetPlatform in targetPlatforms:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_'+targetPlatform +'_tissue-specific_correlations_AltExon_protein_coding.txt'
        #filename = 'AltDatabase/ensembl/'+species+'/'+species+'_'+targetPlatform +'_tissue-specific_correlations_'+coding_type+'.txt'
        
        fn=filepath(filename); x=0
        #print filename
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            if x==0: x=1 ### Ignore header line
            else:
                uid,symbol,(rho,p),tissue = string.split(data,'\t')
                if targetPlatform=='gene':
                    try: uid = gene_translation_db[uid] ### translate from gene to exon array probesets
                    except Exception: uid = ''
                if float(rho)>rho_threshold and float(p)<p_threshold and len(uid)>0 : ### Variable used for testing different thresholds internally
                    try: tissue_to_gene[tissue,uid,symbol]+=1
                    except Exception: tissue_to_gene[tissue,uid,symbol] = 1
                    
    count=0
    for (tissue,uid,symbol) in tissue_to_gene:
        if tissue_to_gene[(tissue,uid,symbol)]>1:
            count+=1
    #print count

def importAndAverageStatsData(expr_input,compendium_filename,platform):
    """ This function takes LineageProfiler z-scores and organizes the samples into groups
    takes the mean results for each group and looks for changes in lineage associations """
    
    groups_dir = string.replace(expr_input,'exp.','groups.')
    groups_dir = string.replace(groups_dir,'stats.','groups.')
    groups_dir = string.replace(groups_dir,'-steady-state.txt','.txt') ### groups is for the non-steady-state file
    export_path = string.replace(compendium_filename,'.txt','_stats.txt')
    export_data = export.ExportFile(export_path)
    print 'Export LineageProfiler database dabg to:',export_path

    ### Import compendium file
    fn=filepath(compendium_filename); row_number=0; compendium_annotation_db={}; x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            print 'Importing the tissue compedium database:',compendium_filename
            headers = t; x=1; index=0
            for i in headers:
                if 'UID' == i: ens_index = index
                if 'AltExon' in compendium_filename: ens_index = ens_index ### Assigned above when analyzing probesets
                elif 'Ensembl' in i: ens_index = index
                if 'marker-in' in i: tissue_index = index+1; marker_in = index
                index+=1
            new_headers = t[:tissue_index]
        else:
            uid = string.split(t[ens_index],'|')[0]
            annotations = t[:tissue_index]
            compendium_annotation_db[uid] = annotations
            
    fn=filepath(expr_input); row_number=0; exp_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if data[0]=='#': row_number = 0
        elif row_number==0:
            group_index_db={}
            ### use comps in the future to visualize group comparison changes
            sample_list,group_sample_db,group_db,group_name_sample_db,comp_groups,comps_name_db = ExpressionBuilder.simpleGroupImport(groups_dir)
            for x in sample_list:
                group_name = group_db[x]
                try: sample_index = t[1:].index(x)
                except Exception, e:
                    print x
                    print t[1:]
                    print e; sys.exit()
                try: group_index_db[group_name].append(sample_index)
                except Exception: group_index_db[group_name] = [sample_index] ### dictionary of group to input file sample indexes
            groups = map(str, group_index_db) ### store group names
            new_sample_list = map(lambda x: group_db[x], sample_list) ### lookup index of each sample in the ordered group sample list
            title = string.join(new_headers+groups,'\t')+'\n' ### output the new sample order (group file order)
            export_data.write(title)
            row_number=1
        else:
            uid = t[0]
            if uid in compendium_annotation_db:
                if platform == 'RNASeq':
                    ### Convert to log2 RPKM values - or counts
                    values = map(lambda x: math.log(float(x),2), t[1:])
                else:
                    try: values = map(float,t[1:])
                    except Exception: values = logTransformWithNAs(t[1:])
                avg_z=[]
                for group_name in group_index_db:
                    group_values = map(lambda x: values[x], group_index_db[group_name]) ### simple and fast way to reorganize the samples
                    avg = statistics.avg(group_values)
                    avg_z.append(str(avg))
                export_data.write(string.join(compendium_annotation_db[uid]+avg_z,'\t')+'\n')
    export_data.close()
    return export_path

def floatWithNAs(values):
    values2=[]
    for x in values:
        try: values2.append(float(x)) #values2.append(math.log(float(x),2))
        except Exception:
            #values2.append(0.00001)
            values2.append('')
    return values2

def importAndAverageExport(expr_input,platform,annotationDB=None,annotationHeader=None,customExportPath=None):
    """ More simple custom used function to convert a exp. or stats. file from sample values to group means """
            
    export_path = string.replace(expr_input,'.txt','-AVERAGE.txt')
    if customExportPath != None: ### Override the above
        export_path = customExportPath
        
    groups_dir = string.replace(expr_input,'exp.','groups.')
    groups_dir = string.replace(groups_dir,'stats.','groups.')
    groups_dir = string.replace(groups_dir,'-steady-state.txt','.txt') ### groups is for the non-steady-state file
    if 'AltExon' in expr_input:
        groups_dir = string.replace(expr_input,'AltExonConfirmed-','groups.')
        groups_dir = string.replace(groups_dir,'AltExon-','groups.')
        groups_dir = string.replace(groups_dir,'AltResults/Clustering','ExpressionInput')
    export_data = export.ExportFile(export_path)
    
    if annotationDB == None:
        annotationDB = {}
        annotationHeader = []
    
    ### CRITICAL!!!! ordereddict is needed to run clustered markerFinder downstream analyses!!!
    import collections

    count=0
    fn=filepath(expr_input); row_number=0; exp_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if data[0]=='#' and row_number==0: row_number = 0
        elif row_number==0:
            if ':' in data:
                tab_data = [] ### remove the group prefixes
                for h in t:
                    h = string.split(h,':')[-1]
                    tab_data.append(h)
                t = tab_data
            try:
                ### OrderedDict used to return the keys in the orders added for markerFinder
                group_index_db=collections.OrderedDict()
            except Exception:
                try:
                    import ordereddict
                    group_index_db = ordereddict.OrderedDict()
                except Exception:
                    group_index_db={}
            ### use comps in the future to visualize group comparison changes
            sample_list,group_sample_db,group_db,group_name_sample_db,comp_groups,comps_name_db = ExpressionBuilder.simpleGroupImport(groups_dir)
            for x in sample_list:
                group_name = group_db[x]
                try: sample_index = t[1:].index(x)
                except Exception, e:
                    print [x]
                    print t[1:]
                    print expr_input
                    print e; sys.exit()
                try: group_index_db[group_name].append(sample_index)
                except Exception: group_index_db[group_name] = [sample_index] ### dictionary of group to input file sample indexes
            groups = map(str, group_index_db) ### store group names
            new_sample_list = map(lambda x: group_db[x], sample_list) ### lookup index of each sample in the ordered group sample list
            title = string.join([t[0]]+groups+annotationHeader,'\t')+'\n' ### output the new sample order (group file order)
            export_data.write(title)
            row_number=1
        else:
            uid = t[0]
            try: values = map(float,t[1:])
            except Exception: values = floatWithNAs(t[1:])
            avg_z=[]
            for group_name in group_index_db:
                group_values = map(lambda x: values[x], group_index_db[group_name]) ### simple and fast way to reorganize the samples
                group_values = [x for x in group_values if x != ''] ### Remove NAs from the group
                five_percent_len = int(len(group_values)*0.05)
                if len(group_values)>five_percent_len:
                    avg = statistics.avg(group_values) #stdev
                else:
                    avg = ''
                avg_z.append(str(avg))
            if uid in annotationDB:
                annotations = annotationDB[uid] ### If provided as an option to the function
            else:
                annotations=[]
            export_data.write(string.join([uid]+avg_z+annotations,'\t')+'\n'); count+=1
    export_data.close()
    print 'Export',count,'rows to:',export_path
    
    return export_path

if __name__ == '__main__':
    Species='Mm'
    filename = ('/Users/saljh8/Desktop/dataAnalysis/Damien/Revised/Guide1/PopulationComparisons/4mo/ExpressionInput/exp.restricted.txt','/Users/saljh8/Desktop/dataAnalysis/Damien/Revised/Guide1/PopulationComparisons/4mo/ExpressionOutput/AVERAGE-restricted.txt')
    analyzeData(filename,Species,"RNASeq","protein_coding",geneToReport=60,correlateAll=True,AdditionalParameters=None,logTransform=True)
    sys.exit()
    averageNIValues('/Users/saljh8/Desktop/LineageProfiler/AltResults/RawSpliceData/Hs/splicing-index/meta-average.txt','/Users/saljh8/Desktop/LineageProfiler/ExpressionInput/stats.meta-steady-state.txt',{});sys.exit()
    dabg_file_dir = '/Users/saljh8/Desktop/LineageProfiler/ExpressionInput/stats.meta-steady-state.txt'
    filterDetectionPvalues(species, dabg_file_dir);sys.exit()
    Platform = 'RNASeq'
    codingType = 'AltExon'
    expr_input = '/Users/nsalomonis/Desktop/Mm_Gene_Meta/ExpressionInput/exp.meta.txt'
    #expr_input = '/home/socr/c/users2/salomoni/conklin/nsalomonis/normalization/Hs_Exon-TissueAtlas/ExpressionInput/stats.meta-steady-state.txt'
    expr_input = '/Volumes/My Passport/dataAnalysis/CardiacRNASeq/BedFiles/ExpressionInput/exp.CardiacRNASeq-steady-state.txt'
    expr_input = '/Volumes/My Passport/dataAnalysis/CardiacRNASeq/BedFiles/AltResults/Clustering/AltExonConfirmed-CardiacRNASeq.txt'
    compendium_filename = 'Hs_exon_tissue-specific_AltExon_protein_coding.txt'
    compendium_filename = '/Users/nsalomonis/Desktop/AltAnalyze/AltDatabase/EnsMart65/ensembl/Hs/Hs_exon_tissue-specific_AltExon_protein_coding.txt'
    compendium_filename = '/home/socr/c/users2/salomoni/AltAnalyze_v.2.0.7-Py/AltDatabase/ensembl/Mm/Mm_gene_tissue-specific_protein_coding.txt'
    #compendium_filename = '/home/socr/c/users2/salomoni/AltAnalyze_v.2.0.7-Py/AltDatabase/ensembl/Hs/Hs_exon_tissue-specific_ncRNA.txt'
    compendium_filename = '/home/socr/c/users2/salomoni/AltAnalyze_v.2.0.7-Py/AltDatabase/ensembl/Mm/Mm_gene_tissue-specific_AltExon_protein_coding.txt'
    importAndAverageExport(expr_input,Platform); sys.exit()
    importAndAverageStatsData(expr_input,compendium_filename,Platform); sys.exit()
    
    returnCommonProfiles(species);sys.exit()
    #codingType = 'ncRNA'
    """
    input_exp_file = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/ExonArray/CP-hESC/AltResults/RawSpliceData/Hs/splicing-index/hESC_differentiation.txt'
    #input_exp_file = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/ExonArray/CP-hESC/AltExpression/FullDatasets/ExonArray/Hs/hESC_differentiation.txt'
    input_exp_file = '/home/socr/c/users2/salomoni/other/boxer/normalization/Hs_Exon-TissueAtlas/AltResults/RawSpliceData/Mm/splicing-index/meta.txt'
    #input_exp_file = '/home/socr/c/users2/salomoni/other/boxer/normalization/Hs_Exon-TissueAtlas/AltExpression/FullDatasets/ExonArray/Hs/meta.txt'
    getAverageExonExpression(Species,Platform,input_exp_file)#;sys.exit()
    
    dataset_file = 'DATASET-Hs_meta-Exon_101111.txt'
    
    #getAverageExpressionValues(dataset_file); sys.exit()
    group_exp_file = 'AVERAGE-Hs_meta-Exon_101111.txt'
    if 'Raw' in input_exp_file:
        group_exp_file = string.replace(input_exp_file,'Raw','AVERAGE')
    else:
        group_exp_file = string.replace(input_exp_file,'FullDatasets','AVERAGE-FullDatasets')
    altexon_correlation_file = analyzeData(group_exp_file,Species,Platform,codingType)
    """
    altexon_correlation_file = 'temp.txt'
    input_exp_file = '/home/socr/c/users2/salomoni/conklin/nsalomonis/normalization/Hs_Exon-TissueAtlas/AltResults/RawSpliceData/Hs/meta.txt'
    getExprValsForNICorrelations(Platform,altexon_correlation_file,input_exp_file)    