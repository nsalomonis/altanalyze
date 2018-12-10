###ExpressionBuilder
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
import bisect
import unique
from stats_scripts import statistics
import math
import reorder_arrays
try: from build_scripts import ExonArray
except: pass
import export
import copy
import time
import traceback
import UI
from import_scripts import BuildAffymetrixAssociations; reload(BuildAffymetrixAssociations)


try:
    from scipy import average as Average
except Exception:
    try: from statistics import avg as Average
    except: pass

use_Tkinter = 'no'
try:
    from Tkinter import *
    use_Tkinter = 'yes'
except ImportError: use_Tkinter = 'yes'; print "\nPmw or Tkinter not found... Tkinter print out not available";
debug_mode = 'no'

### Method specific global variables most easily initialized here
cluster_id=0
cluster_name='clu_0'

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

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
    for key in db_keys: del db_to_clear[key]
    
################# Begin Analysis from parsing files

def checkArrayHeaders(expr_input_dir,expr_group_dir):
    array_names, array_linker_db = getArrayHeaders(expr_input_dir)
    expr_group_list,expr_group_db = importArrayGroups(expr_group_dir,array_linker_db)

def getArrayHeaders(expr_input_dir):
    ### This method is used to check to see if the array headers in the groups and expression files match
    fn=filepath(expr_input_dir); x = 0
    for line in open(fn,'rU').xreadlines():             
      data = cleanUpLine(line)
      headers = string.split(data,'\t')
      if data[0] != '#':
        ### differentiate data from column headers
        if x == 1: break ### Exit out of loop, since we only want the array names
        if x == 0: ### only grab headers if it's the first row
            array_names = []; array_linker_db = {}; d = 0
            for entry in headers[1:]: entry = string.replace(entry,'"',''); array_names.append(entry)
            for array in array_names: array = string.replace(array,'\r',''); array_linker_db[array] = d; d +=1
            x = 1
    return array_names, array_linker_db

def checkExpressionFileFormat(expFile,reportNegatives=False,filterIDs=False):
    """ Determine if the data is log, non-log and increment value for log calculation """
    firstLine=True; convert=False
    inputMax=0; inputMin=10000; increment=0
    expressed_values={}
    startIndex = 1
    for line in open(expFile,'rU').xreadlines():
        line = cleanUpLine(line)
        key = string.split(line,'\t')[0]
        t = string.split(line,'\t')
        if firstLine:
            headers = line
            if 'row_clusters-flat' == t[1]:
                startIndex = 2
            firstLine = False
        else:
            if 'column_clusters-flat' in t:
                continue ### skip this row if analyzing a clustered heatmap file
            try: uid, coordinates = string.split(key,'=')
            except Exception: uid = key
            if filterIDs!=False:
                if uid not in filterIDs:
                    continue
            if '' in t[1:]:
                values = [0 if x=='' else x for x in t[startIndex:]]
            elif 'NA' in t[1:]:
                values = [0 if x=='NA' else x for x in t[startIndex:]]
            else:
                values = t[1:]
            try: values = map(lambda x: float(x), values)
            except Exception:
                print values
                print traceback.format_exc()
            
            if max(values)>inputMax: inputMax = max(values)
            if min(values)<inputMin: inputMin = min(values)
            
    if inputMax>100: ### Thus, not log values
        expressionDataFormat = 'non-log'
        if inputMin<=1: #if inputMin<=1:
            increment = inputMin+1
        convert = True
    else:
        expressionDataFormat = "log"
    #print expressionDataFormat,increment,convert
    if reportNegatives == False:
        return expressionDataFormat,increment,convert
    else:
        ### Report if negative values are present
        increment = inputMin
        if convert: ### Should rarely be the case, as this would indicate that a non-log folds are present in the file
            increment = increment+1
        return expressionDataFormat,increment,convert

def calculate_expression_measures(expr_input_dir,expr_group_dir,experiment_name,comp_group_dir,probeset_db,annotate_db):
    print "Processing the expression file:",expr_input_dir
    
    try: expressionDataFormat,increment,convertNonLogToLog = checkExpressionFileFormat(expr_input_dir)
    except Exception:
        expressionDataFormat = expression_data_format; increment = 0
        if expressionDataFormat == 'non-log': convertNonLogToLog=True
        else: convertNonLogToLog = False
    
    #print convertNonLogToLog, expressionDataFormat, increment
    global array_fold_headers; global summary_filtering_stats; global raw_data_comp_headers; global array_folds
    fn1=filepath(expr_input_dir)
    x = 0; y = 0; d = 0
    blanksPresent=False
    array_folds={}
    for line in open(fn1,'rU').xreadlines():             
      data = cleanUpLine(line)
      if data[0] != '#' and data[0] != '!':
        fold_data = string.split(data,'\t')
        try: arrayid = fold_data[0]
        except Exception: arrayid = 'UID'
        if len(arrayid)>0:
            if arrayid[0]== ' ':
                try: arrayid = arrayid[1:] ### Cufflinks issue
                except Exception: arrayid = ' ' ### can be the first row UID column as blank
            if 'ENSG' in arrayid and '.' in arrayid:
                arrayid = string.split(arrayid,'.')[0]
        else:
            arrayid = 'UID'
        #if 'counts.' in expr_input_dir: arrayid,coordinates = string.split(arrayid,'=') ### needed for exon-level analyses only
        ### differentiate data from column headers
        if x == 1:
            fold_data = fold_data[1:]; fold_data2=[]
            for fold in fold_data:
                fold = string.replace(fold,'"','')
                try:
                    fold = float(fold); fold_data2.append(fold)
                except Exception:
                    fold_data2.append('')
                    blanksPresent = True
                    """
                    print_out = 'WARNING!!! The ID'+arrayid+ 'has an invalid expression value:'+[fold]+'\n. Correct and re-run'
                    try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); sys.exit()
                    except NameError: print print_out; sys.exit()
                    """
            if expressionDataFormat == 'non-log' and (convertNonLogToLog or array_type == 'RNASeq'):
                fold_data3=[] ###Convert numeric expression to log fold (previous to version 2.05 added 1)
                for fold in fold_data2:
                    try:
                        log_fold = math.log((float(fold)+increment),2) ### changed from - log_fold = math.log((float(fold)+1),2) - version 2.05
                        fold_data3.append(log_fold)
                    except ValueError:  ###Not an ideal situation: Value is negative - Convert to zero
                        if float(fold)<=0: log_fold = math.log(1.01,2); fold_data3.append(log_fold)
                        else:
                            fold_data3.append('')
                            blanksPresent = True
                            """
                            print_out = 'WARNING!!! The ID'+arrayid+ 'has an invalid expression value:'+fold+'\n. Correct and re-run'
                            try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); sys.exit()
                            except NameError: print print_out; sys.exit()
                            """
                fold_data2 = fold_data3
            if (array_type == "AltMouse"):
                if arrayid in probeset_db: array_folds[arrayid] = fold_data2; y = y+1
            else: array_folds[arrayid] = fold_data2; y = y+1
        else: #only grab headers if it's the first row
            array_names = []; array_linker_db = {}
            for entry in fold_data[1:]:
                entry = string.replace(entry,'"','')
                if len(entry)>0: array_names.append(entry)
            for array in array_names: #use this to have an orignal index order of arrays
                array = string.replace(array,'\r','') ###This occured once... not sure why
                array_linker_db[array] = d; d +=1
                #add this aftwards since these will also be used as index values
            x = 1
    print len(array_folds),"IDs imported...beginning to calculate statistics for all group comparisons"
    expr_group_list,expr_group_db = importArrayGroups(expr_group_dir,array_linker_db)
    comp_group_list, comp_group_list2 = importComparisonGroups(comp_group_dir)
    
    if 'RPKM' in norm and 'counts.' in expr_input_dir: normalization_method = 'RPKM-counts' ### process as counts if analyzing the counts file
    else: normalization_method = norm
    
    if expressionDataFormat == 'non-log': logvalues=False
    else: logvalues=True
    if convertNonLogToLog: logvalues = True
    try:
        array_folds, array_fold_headers, summary_filtering_stats,raw_data_comp_headers = reorder_arrays.reorder(array_folds,array_names,expr_group_list,
                                comp_group_list,probeset_db,include_raw_data,array_type,normalization_method,fl,logvalues=logvalues,blanksPresent=blanksPresent)
    except Exception:
        print traceback.format_exc(),'\n'
        print_out = 'AltAnalyze encountered an error with the format of the expression file.\nIf the data was designated as log intensities and it is not, then re-run as non-log.'
        try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); root.destroy(); force_exit ### Forces the error log to pop-up
        except NameError: print print_out; sys.exit()

    ### Integrate maximum counts for each gene for the purpose of filtering (RNASeq data only)
    if array_type == 'RNASeq' and 'counts.' not in expr_input_dir: addMaxReadCounts(expr_input_dir)
    
    ### Export these results to a DATASET statistics and annotation results file
    if 'counts.' not in expr_input_dir:
        if array_type == 'RNASeq' and norm == 'RPKM':
            filterRNASeq(count_statistics_db)
            ### Export count summary in GenMAPP format
            if include_raw_data == 'yes': headers = removeRawCountData(array_fold_headers)
            else: headers = array_fold_headers
            exportDataForGenMAPP(headers,'counts')
            
        exportAnalyzedData(comp_group_list2,expr_group_db)
    
        ### Export formatted results for input as an expression dataset into GenMAPP or PathVisio
        if data_type == 'expression':
            if include_raw_data == 'yes': headers = removeRawData(array_fold_headers)
            else: headers = array_fold_headers
            exportDataForGenMAPP(headers,'expression')
            
        try: clearObjectsFromMemory(summary_filtering_stats); clearObjectsFromMemory(array_folds)
        except Exception: null=[]
        try: clearObjectsFromMemory(summary_filtering_stats); summary_filtering_stats=[]
        except Exception: null=[]
    
    else:
        ### When performing an RNASeq analysis on RPKM data, we first perform these analyses on the raw counts to remove fold changes for low expressing genes
        """count_statistics_db={}; count_statistics_headers=[]
        for key in array_folds:
            count_statistics_db[key] = array_folds[key]
        for name in array_fold_headers: count_statistics_headers.append(name)"""

        try: clearObjectsFromMemory(summary_filtering_stats)
        except Exception: null=[]
        try: clearObjectsFromMemory(summary_filtering_stats); summary_filtering_stats=[]
        except Exception: null=[]
    
        return array_folds, array_fold_headers

def filterRNASeq(counts_db):
    ### Parse through the raw count data summary statistics and annotate any comparisons considered NOT EXPRESSED by read count filtering as not expressed (on top of RPKM filtering)
    reassigned = 0; re = 0
    for gene in counts_db:
        i=0 ### keep track of the index (same as RPKM index)
        for val in counts_db[gene]:
            if val =='Insufficient Expression':
                #print val, i, array_folds[gene][i];kill
                if array_folds[gene][i] != 'Insufficient Expression': reassigned = gene, array_folds[gene][i]
                array_folds[gene][i] = 'Insufficient Expression' ### Re-assign the fold changes to this non-numeric value
                re+=1
            i+=1
    #print reassigned, re
    
    
def addMaxReadCounts(filename):
    import RNASeq
    max_count_db,array_names = RNASeq.importGeneCounts(filename,'max')
    for gene in summary_filtering_stats:
        gs = summary_filtering_stats[gene]
        gs.setMaxCount(max_count_db[gene]) ### Shouldn't cause an error, but we want to get an exception if it does (something is wrong with the analysis)        

def simplerGroupImport(group_dir):
    if 'exp.' in group_dir or 'filteredExp.' in group_dir:
        group_dir = string.replace(group_dir,'exp.','groups.')
        group_dir = string.replace(group_dir,'filteredExp.','groups.')
    import collections
    try: sample_group_db = collections.OrderedDict()
    except Exception:
        try:
            import ordereddict
            sample_group_db = ordereddict.OrderedDict()
        except Exception:
            sample_group_db={}
    fn = filepath(group_dir)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        try: sample_filename,group_number,group_name = string.split(data,'\t')
        except Exception:
            #print 'Non-Standard Groups file or missing relationships'
            print string.split(data,'\t')[:10], 'more than 3 columns present in groups file'
            kill
        sample_group_db[sample_filename] = group_name
    return sample_group_db

def simpleGroupImport(group_dir,splitHeaders=False, ignoreComps=False):
    
    """ Used for calculating fold changes prior to clustering for individual samples (genomtric folds) """
    import collections
    try:
        ### OrderedDict used to return the keys in the orders added for markerFinder
        group_sample_db=collections.OrderedDict()
        group_name_db=collections.OrderedDict()
        group_name_sample_db=collections.OrderedDict()
        group_db=collections.OrderedDict()
    except Exception:
        try:
            import ordereddict
            group_sample_db = ordereddict.OrderedDict()
            group_name_db=ordereddict.OrderedDict()
            group_name_sample_db=ordereddict.OrderedDict()
            group_db=ordereddict.OrderedDict()
        except Exception:
            group_sample_db={}
            group_name_db={}
            group_name_sample_db={}
            group_db={}
    sample_list=[]
    group_dir = verifyExpressionFile(group_dir)
    fn = filepath(group_dir)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        try: sample_filename,group_number,group_name = string.split(data,'\t')
        except Exception:
            print traceback.format_exc()
            print "\nWARNING!!! Impropper groups file format detected. Terminating AltAnalyze. The groups file must have only three columns (sampleName, groupNumber, groupName).\n"
            forceGroupsError
        if splitHeaders:
            if '~' in sample_filename: sample_filename = string.split(sample_filename,'~')[-1]
        group_sample_db[sample_filename] = group_name+':'+sample_filename
        try: group_name_sample_db[group_name].append(group_name+':'+sample_filename)
        except Exception: group_name_sample_db[group_name] = [group_name+':'+sample_filename]
        sample_list.append(sample_filename)
        group_db[sample_filename] = group_name
        
        group_name_db[group_number]=group_name ### used by simpleCompsImport
        
    ### Get the comparisons indicated by the user
    if ignoreComps==False: ### Not required for some analyses
        comps_name_db,comp_groups = simpleCompsImport(group_dir,group_name_db)
    else:
        comps_name_db={}; comp_groups=[]
    return sample_list,group_sample_db,group_db,group_name_sample_db,comp_groups,comps_name_db
                
def simpleCompsImport(group_dir,group_name_db):
    """ Used for calculating fold changes prior to clustering for individual samples (genomtric folds) """
    comps_dir = string.replace(group_dir,'groups.','comps.')
    comps_name_db={}
    comp_groups=[]
    comps_dir = verifyExpressionFile(comps_dir)
    fn = filepath(comps_dir)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        try:
            exp_group_num,con_group_num = string.split(data,'\t')
            exp_group_name = group_name_db[exp_group_num]
            con_group_name = group_name_db[con_group_num]
            try: comps_name_db[con_group_name].append(exp_group_name)
            except Exception:
                #comps_name_db[con_group_name] = [exp_group_name] ### If we don't want to include the control samples
                comps_name_db[con_group_name] = [con_group_name] ### Add the control group versus itself the first time
                comps_name_db[con_group_name].append(exp_group_name)
            ### Keep track of the order of the groups for ordering the cluster inputs
            if con_group_name not in comp_groups:
                comp_groups.append(con_group_name)
            if exp_group_name not in comp_groups:
                comp_groups.append(exp_group_name)
        except Exception: pass ### Occurs if there are dummy lines in the file (returns with no values)
    return comps_name_db,comp_groups

def importArrayGroups(expr_group_dir,array_linker_db):
    new_index_order = 0
    import collections
    updated_groups = collections.OrderedDict()
    expr_group_list=[]
    expr_group_db = {} ### use when writing out data
    fn=filepath(expr_group_dir)
    try:
        try: 
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                length = string.join(t,'') ### Some lines can be blank
                if len(length)>2:
                    array_header,group,group_name = t
                    group = int(group)
                    ### compare new to original index order of arrays
                    try:
                        original_index_order = array_linker_db[array_header]
                    except:
                        if array_header+'.bed' in array_linker_db:
                            new_header = array_header+'.bed'
                            original_index_order = array_linker_db[new_header]
                            updated_groups[new_header]=group,group_name
                        elif array_header[:-4] in array_linker_db:
                            new_header = array_header[:-4]
                            original_index_order = array_linker_db[new_header]
                            updated_groups[new_header]=group,group_name      
                        else:
                            print_out = 'WARNING!!! At least one sample-ID listed in the "groups." file (e.g.,'+array_header+')'+'\n is not in the sample "exp." file. See the new file "arrays." with all "exp." header names\nand correct "groups."' 
                            try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!')
                            except Exception: print print_out
                            exportArrayHeaders(expr_group_dir,array_linker_db)
                            try: root.destroy(); sys.exit()
                            except Exception: sys.exit()
                    entry = new_index_order, original_index_order, group, group_name
                    expr_group_list.append(entry)
                    new_index_order += 1 ### add this aftwards since these will also be used as index values
                    expr_group_db[str(group)] = group_name
            expr_group_list.sort() ### sorting put's this in the original array order
        except ValueError:
            print_out = 'The group number "'+group+'" is not a valid integer. Correct before proceeding.'
            try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); root.destroy(); sys.exit()
            except Exception: print print_out; sys.exit()
        
    except Exception,e:
        print traceback.format_exc(),'\n'
        exportArrayHeaders(expr_group_dir,array_linker_db)
        print_out = 'No groups or comps files found for'+expr_group_dir+'... exiting program.'
        try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); root.destroy(); sys.exit()
        except Exception: print print_out; sys.exit()
    if len(updated_groups)>0:
        exportUpdatedGroups(expr_group_dir,updated_groups)
    return expr_group_list,expr_group_db

def exportUpdatedGroups(expr_group_dir,updated_groups):
    eo = export.ExportFile(expr_group_dir)
    for sample in updated_groups:
        eo.write(sample+'\t'+str(updated_groups[sample][0])+'\t'+updated_groups[sample][1]+'\n')
    eo.close()
    print 'The groups file has been updated with bed file sample names'

def exportArrayHeaders(expr_group_dir,array_linker_db):
    new_file = string.replace(expr_group_dir,'groups.','arrays.')
    new_file = string.replace(new_file,'exp.','arrays.')
    new_file = string.replace(new_file,'counts.','arrays.')
    if 'arrays.' not in new_file: new_file = 'arrays.' + new_file ### Can occur if the file does not have 'exp.' in it
    fn=filepath(new_file); data = open(fn,'w')
    for array in array_linker_db: data.write(array+'\n')
    data.close()
    
def importComparisonGroups(comp_group_dir):
    comp_group_list=[]; comp_group_list2=[]
    try:
        fn=filepath(comp_group_dir)
        for line in open(fn,'rU').xreadlines():            
            data = cleanUpLine(line)
            groups = string.split(data,'\t')
            groups2 = groups[0],groups[1] #as a list these would be unhashable
            comp_group_list.append(groups)
            comp_group_list2.append(groups2)
    except Exception: null=[] ### Occcurs when no file present
    return comp_group_list, comp_group_list2

def importMicrornaAssociations(species,report):
    filename = 'AltDatabase/Ensembl/'+species+'/'+species+'_microRNA-Ensembl.txt'
    fn=filepath(filename); ensembl_microRNA_db={}
    for line in open(fn,'rU').xreadlines():            
        data = cleanUpLine(line)
        miR,ens_geneid,sources = string.split(data,'\t')
        miR_annot = miR+'('+sources+')'
        try: ensembl_microRNA_db[ens_geneid].append(miR_annot)
        except KeyError: ensembl_microRNA_db[ens_geneid] = [miR_annot]

    ###Optionally filter out miRs with evidence from just one algorithm (options are 'any' and 'muliple'
    for gene in ensembl_microRNA_db:
        miRs = ensembl_microRNA_db[gene]; miRs.sort()
        if report == 'multiple':
            miRs2=[]
            for mir in miRs:
                if '|' in mir: miRs2.append(mir)
            miRs=miRs2
        miRs = string.join(miRs,', ')
        ensembl_microRNA_db[gene] = miRs
        
    return ensembl_microRNA_db

def importSystemCodes():
    filename = 'Config/source_data.txt'
    fn=filepath(filename); x=0; systems={}
    for line in open(fn,'rU').readlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        system_name=t[0];system_code=t[1]
        if x==0: x=1
        else: systems[system_name] = system_code
    return systems
            
def exportDataForGenMAPP(headers,input_type):
    ###Export summary columns for GenMAPP analysis
    systems = importSystemCodes()
    GenMAPP_file = expression_dataset_output_dir + 'GenMAPP-'+experiment_name+'.txt'
    if 'counts' in input_type:
        GenMAPP_file = string.replace(GenMAPP_file,'GenMAPP-','COUNTS-')
    try: genmapp = export.createExportFile(GenMAPP_file,expression_dataset_output_dir[:-1])
    except RuntimeError:
        export.isFileOpen(GenMAPP_file,expression_dataset_output_dir[:-1])
        genmapp = export.createExportFile(GenMAPP_file,expression_dataset_output_dir[:-1])
    if array_type == "3'array" and 'Ensembl' not in vendor:
        if vendor == 'Affymetrix': system_code = 'X'
        elif vendor == 'Illumina': system_code = 'Il'
        elif vendor == 'Agilent': system_code = 'Ag'
        elif vendor == 'Codelink': system_code = 'Co'
        else:
            ### This is another system selected by the user
            system = string.replace(vendor,'other:','')
            try: system_code = systems[system]
            except Exception: system_code = 'Sy'
    elif array_type != 'AltMouse': system_code = 'En'
    else:
        try: system_code = systems[vendor]
        except Exception: system_code = 'X'
    genmapp_title = ['GeneID','SystemCode'] + headers
    genmapp_title = string.join(genmapp_title,'\t')+'\t'+'ANOVA-rawp'+'\t'+'ANOVA-adjp'+'\t'+'largest fold'+'\n'
    genmapp.write(genmapp_title)
    
    for probeset in array_folds:
        if 'ENS' in probeset and (' ' in probeset or '_' in probeset or ':' in probeset or '-' in probeset) and len(probeset)>9:
            system_code = 'En'
            ensembl_gene = 'ENS'+string.split(probeset,'ENS')[1]
            if ' ' in ensembl_gene:
                ensembl_gene = string.split(ensembl_gene,' ')[0]
            if '_' in ensembl_gene:
                ensembl_gene = string.split(ensembl_gene,'_')[0]
            if ':' in ensembl_gene:
                ensembl_gene = string.split(ensembl_gene,':')[0]
            if '-' in ensembl_gene:
                ensembl_gene = string.split(ensembl_gene,'-')[0]
            data_val = ensembl_gene+'\t'+system_code
        elif ('ENS' in probeset or 'ENF' in probeset) and system_code == 'Sy' and len(probeset)>9:
            system_code = 'En'
            data_val = probeset+'\t'+system_code
        else:
            data_val = probeset+'\t'+system_code
        for value in array_folds[probeset]: data_val += '\t'+ str(value)
        gs = summary_filtering_stats[probeset]
        data_val += '\t'+ str(gs.Pval()) +'\t'+ str(gs.AdjP()) +'\t'+ str(gs.LogFold()) +'\n'
        genmapp.write(data_val)
    genmapp.close()
    exportGOEliteInput(headers,system_code)
    print 'Exported GO-Elite input files...'

def buildCriterion(ge_fold_cutoffs, ge_pvalue_cutoffs, ge_ptype, main_output_folder, operation, UseDownRegulatedLabel=False, genesToExclude={}):
    global array_folds; global m_cutoff; global p_cutoff; global expression_dataset_output_dir
    global ptype_to_use; global use_downregulated_label; use_downregulated_label = UseDownRegulatedLabel
    m_cutoff = math.log(float(ge_fold_cutoffs),2); p_cutoff = ge_pvalue_cutoffs; ptype_to_use = ge_ptype
    expression_dataset_output_dir = string.replace(main_output_folder,'GO-Elite','ExpressionOutput/')
    dir_list = read_directory(expression_dataset_output_dir[:-1])
    if operation == 'summary': filetype = 'DATASET-'
    else: filetype = 'GenMAPP-'
    for filename in dir_list:
        if filetype in filename:
            fn=filepath(expression_dataset_output_dir+filename)
            array_folds = {}; x=0
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line); t = string.split(data,'\t')
                if x==0: x=1; headers = t[1:-2]
                else: 
                    values = t[1:-2]; probeset = t[0]; system_code = t[1]
                    if probeset not in genesToExclude: ### E.g., sex-associated or pseudogenes
                        array_folds[probeset] = values
            if operation == 'summary':
                exportGeneRegulationSummary(filename,headers,system_code)
            else:
                input_files_exported = exportGOEliteInput(headers,system_code)
    array_folds=[]

def excludeGenesImport(filename):
    fn=filepath(filename)
    exclude_genes = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        uid = string.split(data,'\t')[0]
        exclude_genes[uid] = None
    return exclude_genes

def importCountSummary():
    ### Copied code from buildCriterion
    count_summary_db={}
    indexed_headers={}
    filetype = 'COUNTS-'
    dir_list = read_directory(expression_dataset_output_dir[:-1])
    for filename in dir_list:
        if filetype in filename: 
            fn=filepath(expression_dataset_output_dir+filename)
            count_summary_db = {}; x=0
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line); t = string.split(data,'\t')
                if x==0:
                    x=1; i=0
                    for header in t:
                        indexed_headers[header]=i
                        i+=1
                else: 
                    values = t[1:-2]; probeset = t[0]; system_code = t[1]
                    count_summary_db[probeset] = values
            return count_summary_db, indexed_headers
        
def exportGOEliteInput(headers,system_code):
    ### Filter statistics based on user-defined thresholds as input for GO-Elite analysis
    criterion_db={}; denominator_geneids={}; index = 0; ttest=[]
    for column in headers:
        if 'ANOVA' in ptype_to_use and ptype_to_use in column: ttest.append(index) ### Not currently implemented
        elif ptype_to_use in column and 'ANOVA' not in column: ttest.append(index)
        lfi = 2 ### relative logfold index position
        if ptype_to_use == 'adjp': lfi = 3
        index+=1
    
    ### Had to introduce the below code to see if any p-values for a criterion are < 1 (otherwise, include them for GO-Elite)
    exclude_p1={}
    for probeset in array_folds:
        index = 0; af = array_folds[probeset]
        for value in array_folds[probeset]:
            if index in ttest:
                criterion_name = headers[index][5:]
                if criterion_name not in exclude_p1:
                    try: p_value = float(value)
                    except Exception: p_value = 1 ### Occurs when a p-value is annotated as 'Insufficient Expression'
                    if p_value < 1:
                        exclude_p1[criterion_name] = True # Hence, at least one gene has a p<1
            index+=1

    for probeset in array_folds:
        index = 0; af = array_folds[probeset]
        for value in array_folds[probeset]:
            denominator_geneids[probeset]=[]
            if index in ttest:
                criterion_name = headers[index][5:]
                if use_downregulated_label==False:
                    rcn = string.split(criterion_name,'_vs_'); rcn.reverse() ### re-label all downregulated as up (reverse the numerator/denominator)
                    reverse_criterion_names = string.join(rcn,'_vs_')
                    regulation_call = '-upregulated'
                else:
                    reverse_criterion_names = criterion_name
                    regulation_call = '-downregulated'
                try: log_fold = float(af[index-lfi])
                except Exception: log_fold = 0 ### Occurs when a fold change is annotated as 'Insufficient Expression'
                try: p_value = float(value)
                except Exception: p_value = 1 ### Occurs when a p-value is annotated as 'Insufficient Expression'
                try: excl_p1 = exclude_p1[criterion_name] ### You can have adjusted p-values that are equal to 1
                except Exception: excl_p1 = False #Make True to exclude ALL non-adjp < sig value entries
                #print [log_fold, m_cutoff, p_value, p_cutoff];sys.exit()
                if abs(log_fold)>m_cutoff and (p_value<p_cutoff or (p_value==1 and excl_p1==False)):
                    #if p_value == 1: print log_fold, probeset,[value]; sys.exit()
                    try: criterion_db[criterion_name].append((probeset,log_fold,p_value))
                    except KeyError: criterion_db[criterion_name] = [(probeset,log_fold,p_value)]
                    if log_fold>0:
                        try: criterion_db[criterion_name+'-upregulated'].append((probeset,log_fold,p_value))
                        except KeyError: criterion_db[criterion_name+'-upregulated'] = [(probeset,log_fold,p_value)]
                    else:
                        if use_downregulated_label==False:
                            log_fold = abs(log_fold)
                        try: criterion_db[reverse_criterion_names+regulation_call].append((probeset,log_fold,p_value))
                        except KeyError: criterion_db[reverse_criterion_names+regulation_call] = [(probeset,log_fold,p_value)]                    
            index += 1

    ### Format these statistical filtering parameters as a string to include in the file as a record
    if m_cutoff<0: fold_cutoff = -1/math.pow(2,m_cutoff)     
    else: fold_cutoff = math.pow(2,m_cutoff)        
    stat_filters = ' (Regulation criterion: fold > '+str(fold_cutoff)+' and '+ptype_to_use+ ' p-value < '+str(p_cutoff)+')'
    stat_filters_filename = '-fold'+str(fold_cutoff)+'_'+ptype_to_use+str(p_cutoff)
    
    ### Format these lists to export as tab-delimited text files
    if len(criterion_db)>0:
        ### Export denominator gene IDs
        input_files_exported = 'yes'
        expression_dir = string.replace(expression_dataset_output_dir,'ExpressionOutput/','')
        goelite_file = expression_dir +'GO-Elite/denominator/GE.denominator.txt'
        goelite = export.createExportFile(goelite_file,expression_dir+'GO-Elite/denominator')        
        goelite_title = ['GeneID','SystemCode']
        goelite_title = string.join(goelite_title,'\t')+'\n'; goelite.write(goelite_title)
        for probeset in denominator_geneids:
            try:
                if 'ENS' in probeset and (' ' in probeset or '_' in probeset or ':' in probeset or '-' in probeset) and len(probeset)>9:
                    system_code = 'En'
                    ensembl_gene = 'ENS'+string.split(probeset,'ENS')[1]
                    if ' ' in ensembl_gene:
                        ensembl_gene = string.split(ensembl_gene,' ')[0]
                    if '_' in ensembl_gene:
                            ensembl_gene = string.split(ensembl_gene,'_')[0]
                    if ':' in ensembl_gene:
                        ensembl_gene = string.split(ensembl_gene,':')[0]
                    if '-' in ensembl_gene:
                        ensembl_gene = string.split(ensembl_gene,'-')[0]
                    probeset = ensembl_gene
                elif ':' in probeset:
                    probeset = string.split(probeset,':')[0]
                    system_code = 'Sy'
            except Exception:
                pass
            if ('ENS' in probeset or 'ENF' in probeset) and system_code == 'Sy' and len(probeset)>9:
                system_code = 'En'
            values = string.join([probeset,system_code],'\t')+'\n'; goelite.write(values)
        goelite.close()

        ### Export criterion gene IDs and minimal data     
        for criterion_name in criterion_db:
            if criterion_name[-1] == ' ': criterion_file_name = criterion_name[:-1]
            else: criterion_file_name = criterion_name
            if 'upregulated' in criterion_name: elitedir = 'upregulated'
            elif 'downregulated' in criterion_name: elitedir = 'downregulated'
            else: elitedir = 'regulated'
            goelite_file = expression_dir + 'GO-Elite/'+elitedir+'/GE.'+criterion_file_name+stat_filters_filename+'.txt'
            goelite = export.ExportFile(goelite_file)        
            goelite_title = ['GeneID'+stat_filters,'SystemCode',criterion_name+'-log_fold',criterion_name+'-p_value']
            goelite_title = string.join(goelite_title,'\t')+'\n'; goelite.write(goelite_title)
            for (probeset,log_fold,p_value) in criterion_db[criterion_name]:
                try:
                    if 'ENS' in probeset and (' ' in probeset or '_' in probeset or ':' in probeset or '-' in probeset):
                        system_code = 'En'
                        ensembl_gene = 'ENS'+string.split(probeset,'ENS')[1]
                        if ' ' in ensembl_gene:
                            ensembl_gene = string.split(ensembl_gene,' ')[0]
                        if '_' in ensembl_gene:
                            ensembl_gene = string.split(ensembl_gene,'_')[0]
                        if ':' in ensembl_gene:
                            ensembl_gene = string.split(ensembl_gene,':')[0]
                        if '-' in ensembl_gene:
                            ensembl_gene = string.split(ensembl_gene,'-')[0]
                        probeset = ensembl_gene
                    elif ':' in probeset:
                        probeset = string.split(probeset,':')[0]
                        system_code = 'Sy'
                except Exception:
                    pass
                values = string.join([probeset,system_code,str(log_fold),str(p_value)],'\t')+'\n'
                goelite.write(values)
            goelite.close()
    else: input_files_exported = 'no'
    return input_files_exported

def exportGeneRegulationSummary(filename,headers,system_code):
    """
    1) Exports summary results description - Performs a series of targetted queries to report the number
    of coding and non-coding genes expressed along with various regulation and annotation parameters.
    
    2) Exports a global regulated expression table - Values are log2 geometric folds relative to baseline
    of the entire row (all samples) for any criterion met (see ptype_to_use, m_cutoff, p_cutoff). Optionally
    cluster these results downstream and perform QC analyses."""
    
    criterion_db={}; detected_exp_db={}; denominator_geneids={}; index = 0; ttest=[]; avg_columns=[]; all_criterion=[]; all_groups=[]
    search_miR = 'miR-1('
    coding_types = ['protein_coding','ncRNA']
    
    for column in headers:
        if 'ANOVA' in ptype_to_use and ptype_to_use in column: ttest.append(index) ### Not currently implemented
        elif ptype_to_use in column and 'ANOVA' not in column: ttest.append(index)
        lfi = 2 ### relative logfold index position
        if ptype_to_use == 'adjp': lfi = 3
        index+=1
        if 'Protein Classes' in column: pc = index-1
        if 'microRNA' in column: mi = index-1
        if 'avg-' in column: avg_columns.append(index-1)
        if 'Symbol' in column: sy = index-1

    try: count_summary_db,indexed_headers = importCountSummary()
    except Exception: count_summary_db={}

    ### Had to introduce the below code to see if any p-values for a criterion are < 1 (otherwise, include them for GO-Elite)
    exclude_p1={}
    for probeset in array_folds:
        index = 0; af = array_folds[probeset]
        for value in array_folds[probeset]:
            if index in ttest:
                criterion_name = headers[index][5:]
                if criterion_name not in exclude_p1:
                    try: p_value = float(value)
                    except Exception: p_value = 1 ### Occurs when a p-value is annotated as 'Insufficient Expression'
                    if p_value < 1:
                        exclude_p1[criterion_name] = True # Hence, at least one gene has a p<1
            index+=1
            
    genes_to_import={}; probeset_symbol={}
    for probeset in array_folds:
        index = 0; af = array_folds[probeset]
        probeset_symbol[probeset] = af[sy]
        for value in array_folds[probeset]:
            denominator_geneids[probeset]=[]
            if index in avg_columns:
                column_name = headers[index]
                group_name = column_name[4:]
                try: protein_class = af[pc]
                except Exception: protein_class = 'NULL'
                proceed = False
                if array_type == 'RNASeq':
                    if norm == 'RPKM':
                        try: ### Counts file should be present but if not, still proceed
                            i2 = indexed_headers[column_name]
                            if float(af[index])>gene_rpkm_threshold and count_summary_db[probeset][i2]>gene_exp_threshold:
                            #if float(af[index])>5 and count_summary_db[probeset][i2]>50:
                                proceed = True
                        except Exception:
                            proceed = True
                        exp_info = probeset, af[index],count_summary_db[probeset][i2] ### keep track of the expression info
                    else:
                        if float(af[index])>expr_threshold:
                            proceed = True
                        exp_info = probeset, expr_threshold,expr_threshold
                    if proceed==True:
                        if group_name not in all_groups: all_groups.append(group_name)
                        if 'protein_coding' in protein_class:
                            try: detected_exp_db[group_name,'protein_coding'].append(exp_info)
                            except KeyError: detected_exp_db[group_name,'protein_coding']=[exp_info]
                        else:
                            try: detected_exp_db[group_name,'ncRNA'].append(exp_info)
                            except KeyError: detected_exp_db[group_name,'ncRNA']=[exp_info]
            if index in ttest:
                criterion_name = headers[index][5:]
                try: log_fold = float(af[index-lfi])
                except Exception: log_fold = 0 ### Occurs when a fold change is annotated as 'Insufficient Expression'
                try: p_value = float(value)
                except Exception: p_value = 1 ### Occurs when a p-value is annotated as 'Insufficient Expression'
                try: excl_p1 = exclude_p1[criterion_name] ### You can have adjusted p-values that are equal to 1
                except Exception: excl_p1 = False #Make True to exclude ALL non-adjp < sig value entries
                try: protein_class = af[pc]
                except Exception: protein_class = 'NULL'
                if abs(log_fold)>m_cutoff and (p_value<p_cutoff or (p_value==1 and excl_p1==False)):
                    if criterion_name not in all_criterion: all_criterion.append(criterion_name)
                    try: criterion_db[criterion_name]+=1
                    except KeyError: criterion_db[criterion_name] = 1
                    genes_to_import[probeset]=[] ### All, regulated genes (any criterion)
                    
                    if 'protein_coding' in protein_class:
                        if log_fold>0:
                            try: criterion_db[criterion_name,'upregulated','protein_coding']+=1
                            except KeyError: criterion_db[criterion_name,'upregulated','protein_coding'] = 1
                            try:
                                if 'miR-1(' in af[mi]:
                                    try: criterion_db[criterion_name,'upregulated','protein_coding',search_miR[:-1]]+=1
                                    except KeyError: criterion_db[criterion_name,'upregulated','protein_coding',search_miR[:-1]] = 1
                            except Exception: None ### occurs when mi not present
                        else:
                            try: criterion_db[criterion_name,'downregulated','protein_coding']+=1
                            except KeyError: criterion_db[criterion_name,'downregulated','protein_coding'] = 1
                            try:
                                if 'miR-1(' in af[mi]:
                                    try: criterion_db[criterion_name,'downregulated','protein_coding',search_miR[:-1]]+=1
                                    except KeyError: criterion_db[criterion_name,'downregulated','protein_coding',search_miR[:-1]] = 1
                            except Exception: None ### occurs when mi not present
                    else:
                        if protein_class == 'NULL':
                            class_name = 'unclassified'
                        else:
                            class_name = 'ncRNA'
                        if log_fold>0:
                            try: criterion_db[criterion_name,'upregulated',class_name]+=1
                            except KeyError: criterion_db[criterion_name,'upregulated',class_name] = 1
                            try:
                                if 'miR-1(' in af[mi]:
                                    try: criterion_db[criterion_name,'upregulated',class_name,search_miR[:-1]]+=1
                                    except KeyError: criterion_db[criterion_name,'upregulated',class_name,search_miR[:-1]] = 1
                            except Exception: None ### occurs when mi not present
                        else:
                            try: criterion_db[criterion_name,'downregulated',class_name]+=1
                            except KeyError: criterion_db[criterion_name,'downregulated',class_name] = 1
                            try:
                                if 'miR-1(' in af[mi]:
                                    try: criterion_db[criterion_name,'downregulated',class_name,search_miR[:-1]]+=1
                                    except KeyError: criterion_db[criterion_name,'downregulated',class_name,search_miR[:-1]] = 1
                            except Exception: None ### occurs when mi not present
            index += 1
            
    if len(criterion_db)>0:
        try: exportGeometricFolds(expression_dataset_output_dir+filename,array_type,genes_to_import,probeset_symbol)
        except Exception,e:
            print 'Failed to export geometric folds due to:'
            print e ### Don't exit the analysis just report the problem
            print traceback.format_exc()
            None
        
        ### Export lists of expressed genes
        all_expressed={}
        for (group_name,coding_type) in detected_exp_db:
            eo = export.ExportFile(expression_dataset_output_dir+'/ExpressedGenes/'+group_name+'-'+coding_type+'.txt')
            eo.write('GeneID\tRPKM\tCounts\n')
            for (gene,rpkm,counts) in detected_exp_db[(group_name,coding_type)]:
                eo.write(gene+'\t'+str(rpkm)+'\t'+str(counts)+'\n')
                all_expressed[gene]=[]
        try: eo.close()
        except Exception: pass
        
        filename = string.replace(filename,'DATASET-','SUMMARY-')
        filename = string.replace(filename,'GenMAPP-','SUMMARY-')
        summary_path = expression_dataset_output_dir +filename
        export_data = export.ExportFile(summary_path)
        print 'Export summary gene expression results to:',filename

        ### Output Number of Expressed Genes        
        title = ['Biological group']
        for group_name in all_groups: title.append(group_name)
        title = string.join(title,'\t')+'\n'; export_data.write(title)

        if array_type == 'RNASeq':
            ### Only really informative for RNA-Seq data right now, since DABG gene-level stats are not calculated (too time-intensive for this one statistic)
            for coding_type in coding_types:
                if coding_type == 'protein_coding': values = ['Expressed protein-coding genes']
                else: values = ['Expressed ncRNAs']
                for group in all_groups:
                    for group_name in detected_exp_db:
                        if group in group_name and coding_type in group_name:
                            values.append(str(len(detected_exp_db[group_name])))
                values = string.join(values,'\t')+'\n'; export_data.write(values)
            export_data.write('\n')

        if m_cutoff<0: fold_cutoff = -1/math.pow(2,m_cutoff)     
        else: fold_cutoff = math.pow(2,m_cutoff)        
        ### Export criterion gene IDs and minimal data
        export_data.write('Regulation criterion: fold > '+str(fold_cutoff)+' and '+ptype_to_use+ ' p-value < '+str(p_cutoff)+'\n\n')
        for criterion in all_criterion:
            title = [criterion,'up','down','up-'+search_miR[:-1],'down-'+search_miR[:-1]]
            title = string.join(title,'\t')+'\n'; export_data.write(title)
            for coding_type in coding_types:
                values = ['Regulated '+coding_type+' genes']
                for criterion_name in criterion_db:
                    if len(criterion_name)==3:
                        if criterion in criterion_name and ('upregulated',coding_type) == criterion_name[1:]:
                            values.append(str(criterion_db[criterion_name]))
                if len(values)==1: values.append('0')
                for criterion_name in criterion_db:
                    if len(criterion_name)==3:
                        if criterion in criterion_name and ('downregulated',coding_type) == criterion_name[1:]:
                            values.append(str(criterion_db[criterion_name]))
                if len(values)==2: values.append('0')
                for criterion_name in criterion_db:
                    if len(criterion_name)==4:
                        if criterion in criterion_name and ('upregulated',coding_type) == criterion_name[1:-1]:
                            values.append(str(criterion_db[criterion_name]))
                if len(values)==3: values.append('0')
                for criterion_name in criterion_db:
                    if len(criterion_name)==4:
                        if criterion in criterion_name and ('downregulated',coding_type) == criterion_name[1:-1]:
                            values.append(str(criterion_db[criterion_name]))
                if len(values)==4: values.append('0')
                #print values;sys.exit()
                values = string.join(values,'\t')+'\n'; export_data.write(values)
            export_data.write('\n')
        export_data.close()

def exportGeometricFolds(filename,platform,genes_to_import,probeset_symbol,exportOutliers=True,exportRelative=True,customPath=None,convertNonLogToLog=False):
    #print expression_data_format
    #print platform
    #print len(genes_to_import)
    #print exportOutliers
    #print exportRelative
    #print customPath
    """ Import sample and gene expression values from input file, filter, calculate geometric folds
    and export. Use for clustering and QC."""
    #print '\n',filename
    filename = string.replace(filename,'///','/')
    filename = string.replace(filename,'//','/')
    status = 'yes'
    convertGeneToSymbol = True

    if 'ExpressionOutput' in filename:
        filename = string.replace(filename,'-steady-state.txt','.txt')
        export_path = string.replace(filename,'ExpressionOutput','ExpressionOutput/Clustering')
        export_path = string.replace(export_path,'DATASET-','SampleLogFolds-') ### compared to all-sample mean
        export_path2 = string.replace(export_path,'SampleLogFolds-','OutlierLogFolds-') ### compared to all-sample mean
        export_path3 = string.replace(export_path,'SampleLogFolds-','RelativeSampleLogFolds-') ### compared to control groups
        filename = string.replace(filename,'ExpressionOutput','ExpressionInput')
        filename = string.replace(filename,'DATASET-','exp.')
        groups_dir = string.replace(filename,'exp.','groups.')
        if platform != "3'array" and platform != "AltMouse":
            ### This is the extension for gene-level results for exon sensitive platfomrs
            filename1 = string.replace(filename,'.txt','-steady-state.txt')
            status = verifyFile(filename1)
            if status == 'yes':
                filename = filename1
        status = verifyFile(filename)

        if status != 'yes':
            filename = string.replace(filename,'exp.','')
            filename = string.replace(filename,'ExpressionInput','')
            status = verifyFile(filename)
            
    if customPath!=None:
        ### If an alternative output path is desired
        export_path = customPath

    print len(genes_to_import), 'genes with data to export...'
    
    try:
        if expression_data_format == 'non-log' and platform != 'RNASeq': convertNonLogToLog = True
    except Exception: pass

    if status != 'yes':
        print "Clustering expression file not exported due to missing file:"
        print filename
    if status == 'yes':
        export_data = export.ExportFile(export_path)
        if exportOutliers: export_outliers = export.ExportFile(export_path2)
        if exportRelative: export_relative = export.ExportFile(export_path3)
        print 'Export inputs for clustering to:',export_path

        expressionDataFormat,increment,convertNonLogToLog = checkExpressionFileFormat(filename)
        #print expressionDataFormat,increment,convertNonLogToLog
        fn=filepath(filename); row_number=0; exp_db={}; relative_headers_exported = False
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')

            if data[0]=='#' and row_number==0: row_number = 0
            elif row_number==0:
                sample_list,group_sample_db,group_db,group_name_sample_db,comp_groups,comps_name_db = simpleGroupImport(groups_dir)

                try: sample_index_list = map(lambda x: t[1:].index(x), sample_list) ### lookup index of each sample in the ordered group sample list
                except Exception:
                    missing=[]
                    for x in sample_list:
                        if x not in t[1:]: missing.append(x)
                    print 'missing:',missing
                    print t
                    print sample_list
                    print filename, groups_dir
                    print 'Unknown Error!!! Skipping cluster input file build (check column and row formats for conflicts)'; forceExit
                new_sample_list = map(lambda x: group_sample_db[x], sample_list) ### lookup index of each sample in the ordered group sample list
                title = string.join([t[0]]+new_sample_list,'\t')+'\n' ### output the new sample order (group file order)
                export_data.write(title)
                if exportOutliers: export_outliers.write(title)
                if exportRelative:
                    ### Used for the relative fold calculation
                    group_index_db={}
                    for x in sample_list:
                        group_name = group_db[x]
                        sample_index = t[1:].index(x)
                        try: group_index_db[group_name].append(sample_index)
                        except Exception: group_index_db[group_name] = [sample_index] ### dictionary of group to input file sample indexes
                row_number=1
            else:
                gene = t[0]

                if platform == 'RNASeq':
                    ### Convert to log2 RPKM values - or counts
                    try: values = map(lambda x: math.log(float(x)+increment,2), t[1:])
                    except Exception:
                        if convertNonLogToLog:
                            values = logTransformWithNAs(t[1:],increment)
                        else:
                            values = TransformWithNAs(t[1:])
                else:
                    try:
                        if convertNonLogToLog:
                            values = map(lambda x: math.log(float(x)+increment,2), t[1:])
                        else:
                            values = map(float,t[1:])
                    except Exception:
                        if convertNonLogToLog:
                            values = logTransformWithNAs(t[1:],increment)
                        else:
                            values = TransformWithNAs(t[1:])
                
                ### Calculate log-fold values relative to the mean of all sample expression values

                values = map(lambda x: values[x], sample_index_list) ### simple and fast way to reorganize the samples

                try: avg = statistics.avg(values)
                except Exception:
                    values2=[]
                    for v in values:
                        try: values2.append(float(v))
                        except Exception: pass
                    try: avg = statistics.avg(values2)
                    except Exception:
                        if len(values2)>0: avg = values2[0]
                        else: avg = 0

                try: log_folds = map(lambda x: (x-avg), values)
                except Exception: 
                    log_folds=[]
                    for x in values:
                        try: log_folds.append(x-avg)
                        except Exception: log_folds.append('')

                if gene in genes_to_import:
                    ### Genes regulated in any user-indicated comparison according to the fold and pvalue cutoffs provided
                    log_folds = map(lambda x: str(x), log_folds)
                    try:
                        """
                        if convertGeneToSymbol:
                            if gene == probeset_symbol[gene]:
                                gene2 = gene
                                convertGeneToSymbol = False
                            else:
                                gene2 = gene+' '+probeset_symbol[gene]
                        else:
                            gene2 = gene
                        """
                        gene2 = gene+' '+probeset_symbol[gene]
                    except Exception: gene2 = gene
                    #print [gene2,gene];sys.exit()
                    if len(t[1:])!=len(log_folds):
                        log_folds = t[1:] ### If NAs - output the original values
                    export_data.write(string.join([gene2]+log_folds,'\t')+'\n')

                    if exportRelative:
                        ### Calculate log-fold values relative to the mean of each valid group comparison
                        control_group_avg={}; comps_exp_db={}
                        for group_name in comps_name_db: ### control group names
                            con_group_values = map(lambda x: values[x], group_index_db[group_name]) ### simple and fast way to reorganize the samples
                            try: control_group_avg[group_name] = statistics.avg(con_group_values) ### store the mean value of each control group
                            except Exception:
                                con_group_values2=[]
                                for val in con_group_values:
                                    try: con_group_values2.append(float(val))
                                    except Exception: pass
                                    try: control_group_avg[group_name] = statistics.avg(con_group_values)
                                    except Exception:
                                        if len(con_group_values)>0:
                                            control_group_avg[group_name] = con_group_values[0]
                                        else: control_group_avg[group_name] = 0.0
                            for exp_group in comps_name_db[group_name]:
                                try: comps_exp_db[exp_group].append(group_name) ### Create a reversed version of the comps_name_db, list experimental as the key
                                except Exception: comps_exp_db[exp_group] = [group_name]
                                
                        relative_log_folds=[] ### append all new log folds to this list
                        relative_column_names=[]
                        for group_name in comp_groups:
                            if group_name in comps_exp_db: ### Hence, the group has a designated control (controls may not) - could have the control group be a control for the control samples
                                group_values = map(lambda x: values[x], group_index_db[group_name]) ### simple and fast way to reorganize the samples
                                for control_group_name in comps_exp_db[group_name]:
                                    con_avg = control_group_avg[control_group_name]
                                    
                                    try:
                                        relative_log_folds += map(lambda x: str(x-con_avg), group_values) ### calculate log-folds and convert to strings
                                    except Exception:
                                        relative_log_folds=[]
                                        for x in group_values:
                                            try: relative_log_folds.append(str(x-con_avg))
                                            except Exception: relative_log_folds.append('')
                            
                                    if relative_headers_exported == False:
                                        exp_sample_names = group_name_sample_db[group_name]
                                        relative_column_names += map(lambda x: (x+' vs '+control_group_name), exp_sample_names) ### add column names indicating the comparison
        
                        if relative_headers_exported == False:
                            title = string.join(['UID']+relative_column_names,'\t')+'\n' ### Export column headers for the relative fold changes
                            export_relative.write(title)
                            relative_headers_exported = True
                        if len(t[1:])!=len(relative_log_folds):
                            relative_log_folds = t[1:] ### If NAs - output the original values
                        export_relative.write(string.join([gene2]+relative_log_folds,'\t')+'\n')
                            
                elif exportOutliers:
                    ### When a gene is regulated and not significant, export to the outlier set
                    try: gene2 = gene+' '+probeset_symbol[gene]
                    except Exception: gene2 = gene
                    ### These are defaults we may allow the user to control later
                    log_folds = [0 if x=='' else x for x in log_folds] ### using list comprehension, replace '' with 0
                    if max([max(log_folds),abs(min(log_folds))])>1:
                        proceed = True
                        if platform == 'RNASeq':
                            if max(values)<0.1: proceed = False
                        if proceed == True:
                            log_folds = map(lambda x: str(x), log_folds)
                            if len(t[1:])!=len(log_folds):
                                log_folds = t[1:] ### If NAs - output the original values
                            export_outliers.write(string.join([gene2]+log_folds,'\t')+'\n')
                            
                row_number+=1 ### Keep track of the first gene as to write out column headers for the relative outputs
                
        export_data.close()
        if exportOutliers: export_outliers.close()
        if exportRelative: export_relative.close()

def logTransformWithNAs(values,increment):
    values2=[]
    for x in values:
        try: values2.append(math.log(float(x)+increment,2))
        except Exception:
            values2.append('')
    return values2

def TransformWithNAs(values):
    values2=[]
    for x in values:
        try: values2.append(float(x))
        except Exception:
            values2.append('')
    return values2

def importAndOrganizeLineageOutputs(expr_input,filename,platform):
    """ This function takes LineageProfiler z-scores and organizes the samples into groups
    takes the mean results for each group and looks for changes in lineage associations """
    
    groups_dir = string.replace(expr_input,'exp.','groups.')
    groups_dir = string.replace(groups_dir,'-steady-state.txt','.txt') ### groups is for the non-steady-state file
    export_path = string.replace(filename,'ExpressionOutput','ExpressionOutput/Clustering')
    export_path = string.replace(export_path,'.txt','-groups.txt')
    export_data = export.ExportFile(export_path)
    export_pathF = string.replace(export_path,'.txt','_filtered.txt')
    export_dataF = export.ExportFile(export_pathF)
    print 'Export inputs for clustering to:',export_path
        
    fn=filepath(filename); row_number=0; exp_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if data[0]=='#': row_number = 0
        elif row_number==0:
            group_index_db={}
            ### use comps in the future to visualize group comparison changes
            sample_list,group_sample_db,group_db,group_name_sample_db,comp_groups,comps_name_db = simpleGroupImport(groups_dir)
            for x in sample_list:
                group_name = group_db[x]
                sample_index = t[1:].index(x)
                try: group_index_db[group_name].append(sample_index)
                except Exception: group_index_db[group_name] = [sample_index] ### dictionary of group to input file sample indexes
            groups = map(str, group_index_db) ### store group names
            new_sample_list = map(lambda x: group_db[x], sample_list) ### lookup index of each sample in the ordered group sample list
            title = string.join([t[0]]+groups,'\t')+'\n' ### output the new sample order (group file order)
            export_data.write(title)
            export_dataF.write(title)
            row_number=1
        else:
            tissue = t[0]
            if platform == 'RNASeq' and 'LineageCorrelations' not in filename:
                ### Convert to log2 RPKM values - or counts
                values = map(lambda x: math.log(float(x),2), t[1:])
            else:
                values = map(float,t[1:])
            avg_z=[]; avg_z_float=[]
            for group_name in group_index_db:
                group_values = map(lambda x: values[x], group_index_db[group_name]) ### simple and fast way to reorganize the samples
                avg = statistics.avg(group_values)
                avg_z.append(str(avg))
                avg_z_float.append(avg)
            export_data.write(string.join([tissue]+avg_z,'\t')+'\n')
            if max(avg_z_float)>1:
                export_dataF.write(string.join([tissue]+avg_z,'\t')+'\n')
    export_data.close(); export_dataF.close()
    return export_path,export_pathF

def removeRawData(array_fold_headers):
    ### Prior to exporting data for GenMAPP, remove raw data columns
    columns_with_stats=[]; i=0; stat_headers = ['avg', 'log_fold', 'fold', 'rawp', 'adjp']; filtered_headers=[]
    for header in array_fold_headers:
        broken_header = string.split(header,'-')
        ### Only keep those headers and indexes with recognized ExpressionBuilder inserted prefixes
        if broken_header[0] in stat_headers: columns_with_stats.append(i); filtered_headers.append(header)
        i+=1

    for probeset in array_folds:
        filtered_list=[]
        for i in columns_with_stats: filtered_list.append(array_folds[probeset][i])
        array_folds[probeset] = filtered_list ### Re-assign values of the db
    return filtered_headers

def removeRawCountData(array_fold_headers):
    ### Prior to exporting data for GenMAPP, remove raw data columns
    columns_with_stats=[]; i=0; stat_headers = ['avg', 'log_fold', 'fold', 'rawp', 'adjp']; filtered_headers=[]
    for header in array_fold_headers:
        broken_header = string.split(header,'-')
        ### Only keep those headers and indexes with recognized ExpressionBuilder inserted prefixes
        if broken_header[0] in stat_headers: columns_with_stats.append(i); filtered_headers.append(header)
        i+=1

    for probeset in count_statistics_db:
        filtered_list=[]
        for i in columns_with_stats: filtered_list.append(count_statistics_db[probeset][i])
        count_statistics_db[probeset] = filtered_list ### Re-assign values of the db
    return filtered_headers

def exportAnalyzedData(comp_group_list2,expr_group_db):
    report = 'multiple'; report = 'single'
    try: ensembl_microRNA_db = importMicrornaAssociations(species,report)
    except IOError: ensembl_microRNA_db={}
    if array_type != "AltMouse" and array_type != "3'array":
        try:
            from build_scripts import EnsemblImport
            gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,array_type,'key_by_array')
        except Exception: gene_location_db={} 
    if data_type == 'expression':
        new_file = expression_dataset_output_dir + 'DATASET-'+experiment_name+'.txt'
        try: data = export.createExportFile(new_file,expression_dataset_output_dir[:-1])
        except RuntimeError:
            export.isFileOpen(new_file,expression_dataset_output_dir[:-1])
            data = export.createExportFile(new_file,expression_dataset_output_dir[:-1])
        try: custom_annotation_dbase = importCustomAnnotations()
        except Exception: custom_annotation_dbase={}
        x=0;y=0;z=0
        for arrayid in array_folds:
            if arrayid in annotate_db and arrayid in probeset_db: x = 1
            if arrayid in annotate_db: y = 1
            if arrayid in conventional_array_db: z = 1
            break
        Vendor = vendor ### Need to rename as re-assigning will cause a global conflict error
        for arrayid in array_folds:
            if 'ENS' in arrayid and Vendor == 'Symbol': 
                Vendor = 'Ensembl'
                break
        if array_type != "AltMouse" and (array_type != "3'array" or 'Ensembl' in Vendor):
            #annotate_db[gene] = symbol, definition,rna_processing
            #probeset_db[gene] = transcluster_string, exon_id_string
            title = ['Ensembl_gene','Definition','Symbol','Transcript_cluster_ids','Constitutive_exons_used','Constitutive_IDs_used','Putative microRNA binding sites','Select Cellular Compartments','Select Protein Classes','Chromosome','Strand','Genomic Gene Corrdinates','GO-Biological Process','GO-Molecular Function','GO-Cellular Component','WikiPathways']
            title = string.join(title,'\t')
        elif arrayCode == 3: ### Code indicates this array probes only for small RNAs
            title = ['Probeset ID','Sequence Type','Transcript ID','Species Scientific Name','Genomic Location']
            title = string.join(title,'\t')
        elif x == 1:
            title = "Probesets" +'\t'+ 'Definition' +'\t'+ 'Symbol' +'\t'+ 'affygene' +'\t'+ 'exons' +'\t'+ 'probe_type_call' +'\t'+ 'ensembl'
        elif y==1: title = "Probesets" +'\t'+ 'Symbol' +'\t'+ 'Definition'
        elif array_type == "3'array":
             title = ['Probesets','Symbol','Definition','Ensembl_id','Entrez_id','Unigene_id','GO-Process','GO-Function','GO-Component','Pathway_info','Putative microRNA binding sites','Select Cellular Compartments','Select Protein Classes']
             title = string.join(title,'\t')
        else: title = "Probesets"
        for entry in array_fold_headers: title = title + '\t' + entry
        title += '\t'+ 'ANOVA-rawp' +'\t'+ 'ANOVA-adjp' +'\t'+'largest fold'
        if array_type == 'RNASeq': title += '\t'+ 'maximum sample read count'
        data.write(title+'\n')
        for arrayid in array_folds:
            if arrayCode == 3:
                ca = conventional_array_db[arrayid]
                definition = ca.Description()
                symbol = ca.Symbol()
                data_val = [arrayid,ca.Description(),ca.Symbol(),ca.Species(),ca.Coordinates()]
                data_val = string.join(data_val,'\t')
            elif array_type != 'AltMouse' and (array_type != "3'array" or 'Ensembl' in Vendor):
                try: definition = annotate_db[arrayid][0]; symbol = annotate_db[arrayid][1]; rna_processing = annotate_db[arrayid][2]
                except Exception: definition=''; symbol=''; rna_processing=''
                report = 'all'
                try: miRs = ensembl_microRNA_db[arrayid]
                except KeyError: miRs = ''
                try:
                    trans_cluster = probeset_db[arrayid][0]
                    exon_ids = probeset_db[arrayid][1]
                    probesets = probeset_db[arrayid][2]
                except Exception:
                    trans_cluster='';exon_ids='';probesets=''
                try: compartment,custom_class = custom_annotation_dbase[arrayid]
                except KeyError: compartment=''; custom_class=''
                try: chr,strand,start,end = gene_location_db[arrayid]
                except Exception: chr=''; strand=''; strand=''; start=''; end=''
                try: pi = conventional_array_db[arrayid]; process = pi.Process(); function=pi.Function(); component=pi.Component(); pathway = pi.Pathway()
                except Exception: process=''; function=''; component=''; pathway=''
                data_val = [arrayid,symbol,definition,trans_cluster,exon_ids,probesets,miRs,compartment,custom_class,chr,strand,start+'-'+end,process,function,component,pathway]
                data_val = string.join(data_val,'\t')
            elif arrayid in annotate_db and arrayid in probeset_db: ### This is for the AltMouse array
                symbol = annotate_db[arrayid][0]
                definition = annotate_db[arrayid][1]
                affygene = probeset_db[arrayid][0][0:-1]     #probeset_db[probeset] = affygene,exons,probe_type_call,ensembl
                exons = probeset_db[arrayid][1]
                probe_type_call = probeset_db[arrayid][2]
                ensembl = probeset_db[arrayid][3]
                data_val = arrayid +'\t'+ definition +'\t'+ symbol +'\t'+ affygene +'\t'+ exons +'\t'+ probe_type_call +'\t'+ ensembl
            elif arrayid in annotate_db:
                definition = annotate_db[arrayid][0]
                symbol = annotate_db[arrayid][1]
                data_val = arrayid +'\t'+ definition +'\t'+ symbol
            elif array_type == "3'array" and 'Ensembl' not in Vendor:
                try:
                    ca = conventional_array_db[arrayid]
                    definition = ca.Description()
                    symbol = ca.Symbol()
                    ens = ca.EnsemblString()
                    entrez = ca.EntrezString()
                    unigene = ca.UnigeneString()
                    pathway_info = ca.PathwayInfo()
                    component = ca.GOComponentNames(); process = ca.GOProcessNames(); function = ca.GOFunctionNames()
                    compartment=''; custom_class=''; miRs=''
                    if len(ens)>0:
                        if ens[0]=='|': ens = ens[1:]
                    store=[]
                    for ens_gene in ca.Ensembl(): ### Add Custom Annotation layer
                        try: compartment,custom_class = custom_annotation_dbase[ens_gene]
                        except KeyError: null=[]
                        try: miRs = ensembl_microRNA_db[ens_gene]
                        except KeyError: null=[]
                        if 'protein_coding' in custom_class and len(store)==0: ### Use the first instance only
                            store = miRs,compartment,custom_class+'('+ens_gene+')'
                    if len(store)>0: ### pick the Ensembl with protein coding annotation to represent (as opposed to aligning annotated pseudo genes)
                        miRs,compartment,custom_class = store
                except KeyError:
                    definition=''; symbol=''; ens=''; entrez=''; unigene=''; pathway_info=''
                    process=''; function=''; component=''; compartment='' ;custom_class=''; miRs=''
                data_val = [arrayid,symbol,definition,ens,entrez,unigene,process,function,component,pathway_info,miRs,compartment,custom_class]
                data_val = string.join(data_val,'\t')
            else:
                data_val = arrayid
            for value in array_folds[arrayid]:
                data_val = data_val + '\t' + str(value)
            gs = summary_filtering_stats[arrayid]
            #if arrayid == '1623863_a_at': print [gs.LogFold()]
            data_val += '\t'+ str(gs.Pval()) +'\t'+ str(gs.AdjP()) +'\t'+ str(gs.LogFold())
            if array_type == 'RNASeq': data_val+= '\t'+ gs.MaxCount()
            data_val = string.replace(data_val,'\n','')
            data.write(data_val+'\n')
        data.close()
        print "Full Dataset with statistics:",'DATASET-'+experiment_name+'.txt', 'written'
        gene_location_db=[]
        ensembl_microRNA_db=[]
        custom_annotation_dbase=[]
        
def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def convert_to_list(database):
    db1=[]; db2=[]; temp_list=[]
    for key in database:
        list = database[key]
        #print key,list,dog  #32 [(2, 1.1480585565447154), (3, 0.72959188370731742), (0, 0.0), (1, -0.60729064216260165)]
        list.sort()
        temp_list=[]
        temp_list.append(key)
        for entry in list:
            avg_fold = entry[1]
            temp_list.append(avg_fold)
        #print temp_list, dog  #[32, 0.0, -0.60729064216260165, 1.1480585565447154, 0.72959188370731742]
        db1.append(temp_list)
    db1.sort()
    return db1

def import_annotations(filename):
    fn=filepath(filename)
    annotation_dbase = {}
    for line in open(fn,'rU').xreadlines():
        try:
            data = cleanUpLine(line)
            try: probeset,definition,symbol,rna_processing = string.split(data,'\t')
            except ValueError:
                probeset,definition,symbol = string.split(data,'\t')
                rna_processing  = ''
            annotation_dbase[probeset] = definition, symbol,rna_processing
        except ValueError: continue
    return annotation_dbase

def importCustomAnnotations():
    ### Combine non-coding Ensembl gene annotations with UniProt functional annotations
    try: custom_annotation_dbase = importTranscriptBiotypeAnnotations(species)
    except Exception: custom_annotation_dbase = {}
    try: housekeeping_genes=BuildAffymetrixAssociations.getHousekeepingGenes(species)
    except Exception: housekeeping_genes=[]
    print len(custom_annotation_dbase),'Ensembl Biotypes and', len(housekeeping_genes),'housekeeping genes.'
    
    for ens_gene in housekeeping_genes:
        if ens_gene not in custom_annotation_dbase: custom_annotation_dbase[ens_gene] = '','housekeeping'
        else:
            compartment,custom_class = custom_annotation_dbase[ens_gene]
            custom_class+='|housekeeping'
            custom_annotation_dbase[ens_gene] = compartment,custom_class

    filename = 'AltDatabase/uniprot/'+species+'/custom_annotations.txt'
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        ens_gene,compartment,custom_class = t[:3]
        if ens_gene in custom_annotation_dbase:
            biotype = custom_annotation_dbase[ens_gene][1]
            if len(custom_class)>0: custom_class+='|'+biotype
            else: custom_class=biotype
        custom_annotation_dbase[ens_gene] = compartment,custom_class
    return custom_annotation_dbase

def importTranscriptBiotypeAnnotations(species):
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-biotypes.txt'
    fn=filepath(filename); biotype_db = {}; custom_annotation_dbase={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        gene,transcript,biotype = string.split(data,'\t')
        ### Determine if only one annotation is associated with each gene
        try: biotype_db[gene][biotype]=[]
        except Exception: biotype_db[gene] = db = {biotype:[]}
        
    for gene in biotype_db:
        db = biotype_db[gene]
        if len(db)==1 and 'protein_coding' not in db:
            for biotype in db: ### Non-coding gene annotation
                custom_annotation_dbase[gene] = '',biotype
        elif 'protein_coding' in db:
            custom_annotation_dbase[gene] = '','protein_coding'
        elif 'transcribed_unprocessed_pseudogene' in db:
            custom_annotation_dbase[gene] = '','transcribed_unprocessed_pseudogene'
        else:
            ls=[] ### otherwise include all gene annotations
            for i in db: ls.append(i)
            ls = string.join(ls,'|')
            custom_annotation_dbase[gene] = '',ls
            
    return custom_annotation_dbase

def importAltMerge(import_type):
    ### Import Probeset annotations
    try:
        ensembl_db={}; fn=filepath('AltDatabase/Mm/AltMouse/AltMouse-Ensembl.txt')
        for line in open(fn,'rU').xreadlines():             
            data = cleanUpLine(line)
            affygene,ensembl = string.split(data,'\t')
            ensembl_db[affygene]=ensembl
        #print len(ensembl_db),'Ensembl-AltMouse relationships imported'
    except TypeError: null=[]
        
    ### Import Probeset annotations
    probeset_annotation_file = "AltDatabase/"+species+'/'+array_type+'/'+ "MASTER-probeset-transcript.txt"
    probeset_db = {}; constitutive_db = {}; fn=filepath(probeset_annotation_file); replacements=0
    for line in open(fn,'rU').xreadlines():             
        probeset_data = cleanUpLine(line)
        probeset,affygene,exons,transcript_num,transcripts,probe_type_call,ensembl,block_exon_ids,block_structure,comparison_info = string.split(probeset_data,'\t')
        if probeset == "Probeset": continue
        else:
            if affygene[:-1] in ensembl_db: ensembl = ensembl_db[affygene[:-1]]; replacements+=1
            if import_type == 'full': ### Mimics the structure of ExonArrayEnsemblRules.reimportEnsemblProbesets() dictionary probe_association_db
                probe_data = affygene,affygene,exons,'','core'
                probeset_db[probeset] = probe_data
            else: probeset_db[probeset] = affygene,exons,probe_type_call,ensembl
            if probe_type_call == 'gene':
                try: constitutive_db[affygene].append(probeset)
                except KeyError: constitutive_db[affygene] = [probeset]
    return probeset_db, constitutive_db

def parse_custom_annotations(filename):
    custom_array_db = {}
    x=0
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        array_data = data
        array_id,probeset,other = string.split(array_data,'\t')  #remove endline
        custom_array_db[array_id] = probeset,other
    print len(custom_array_db), "custom array entries process"
    return custom_array_db

def remoteLineageProfiler(params,expr_input_dir,ArrayType,Species,Vendor,customMarkers=False,specificPlatform=False,visualizeNetworks=True):
    global species
    global array_type
    global vendor
    global remoteAnalysis
    global fl
    remoteAnalysis = True
    species = Species
    array_type = ArrayType
    vendor = Vendor
    fl = params
    graphics_links = []
    if 'ExpressionInput' in expr_input_dir:
        output_dir = string.replace(expr_input_dir,'ExpressionInput', 'ExpressionOutput')
        root_dir = export.findParentDir(output_dir)
    else:
        root_dir = export.findParentDir(expr_input_dir)+'ExpressionOutput/'
    try:
        ### If this directory exists, create a global variable for it
        dir_list = read_directory(root_dir[:-1])
        global expression_dataset_output_dir
        global experiment_name
        experiment_name = string.replace(export.findFilename(expr_input_dir)[:-4],'exp.','')
        expression_dataset_output_dir = root_dir
    except Exception:
        #print traceback.format_exc()
        None
    
    graphic_links = performLineageProfiler(expr_input_dir,graphics_links,customMarkers,specificPlatform=specificPlatform,visualizeNetworks=visualizeNetworks)
    return graphic_links
    
def performLineageProfiler(expr_input_dir,graphic_links,customMarkers=False,specificPlatform=False,visualizeNetworks=True):
    try:
        from visualization_scripts import WikiPathways_webservice
        import LineageProfiler; reload(LineageProfiler)
        start_time = time.time()
        try:
            compendium_type = fl.CompendiumType()
            compendium_platform = fl.CompendiumPlatform()
        except Exception:
            compendium_type = 'protein_coding'
            compendium_platform = 'exon'
        #print 'Compendium platform selected:',compendium_platform
        print 'Biological data type to examine:',compendium_type
        
        try: ### Works when expression_dataset_output_dir is defined
            exp_output = expression_dataset_output_dir + 'DATASET-'+experiment_name+'.txt'
            #try: array_type_data = vendor, array_type
            array_type_data = array_type
        except Exception: ### Otherwise, user directly supplied file is used
            array_type_data = vendor, array_type
            exp_output = export.findParentDir(expr_input_dir)+'/LineageCorrelations-'+export.findFilename(expr_input_dir)
        
        if specificPlatform == False:
            compendium_platform = 'exon'
        status = False
        compareToAll=False
        """
        print species
        print array_type_data
        print expr_input_dir
        print exp_output
        print compendium_type
        print compendium_platform
        print customMarkers
        """
        if 'steady-state.txt' in expr_input_dir:
            status = verifyFile(expr_input_dir)
            if status != 'yes':
                expr_input_dir = string.replace(expr_input_dir,'-steady-state.txt','.txt')
                array_type_data = "3'array"
        try:
            zscore_output_dir1 = LineageProfiler.runLineageProfiler(species,array_type_data,
                    expr_input_dir,exp_output,compendium_type,compendium_platform,customMarkers); status = True
            #zscore_output_dir1 = None
        except Exception:
            print traceback.format_exc(),'\n'
            zscore_output_dir1 = None
        if compareToAll:
            try:
                zscore_output_dir2 = LineageProfiler.runLineageProfiler(species,array_type_data,expr_input_dir, exp_output,compendium_type,'gene',customMarkers); status = True
                #zscore_output_dir2 = None
            except Exception: zscore_output_dir2 = None
            try:
                zscore_output_dir3 = LineageProfiler.runLineageProfiler(species,array_type_data,expr_input_dir, exp_output,compendium_type,"3'array",customMarkers); status = True
                #zscore_output_dir3 = None
            except Exception: zscore_output_dir3 = None

            zscore_output_dirs=[zscore_output_dir1,zscore_output_dir2,zscore_output_dir3]
        else:
            zscore_output_dirs=[zscore_output_dir1]
        ### Create a combined zscore_output_dir output using all predictions
        zscore_output_dir = combineLPResultFiles(zscore_output_dirs)
        if status == False:
            #print traceback.format_exc(),'\n'
            if species != 'Mm' and species != 'Hs':
                print 'LineageProfiler analysis failed (possibly unsupported species).'
        else:
            time_diff = str(round(time.time()-start_time,1))
            print 'LineageProfiler analysis completed in %s seconds' % time_diff
            try: ### If the file a groups file exists in the expected directory structure -> export a groups z-score file
                export_path, export_filtered_path = importAndOrganizeLineageOutputs(expr_input_dir,zscore_output_dir,array_type)
            except Exception:
                export_path = zscore_output_dir ### keeps the sample z-score file as input
                export_filtered_path = None
            
            ### Output a heat map of the sample Z-score correlations
            graphic_links = LineageProfiler.visualizeLineageZscores(zscore_output_dir,export_path,graphic_links)
            if export_filtered_path != None:
                try: LineageProfiler.visualizeLineageZscores(export_filtered_path,export_path,graphic_links) ### Just output a heatmap of filtered grouped terms
                except Exception: pass
                
            ### Color the TissueMap from WikiPathways using their webservice
            if customMarkers==False and visualizeNetworks:
                print 'Coloring LineageMap profiles using WikiPathways webservice...'
                graphic_links = WikiPathways_webservice.viewLineageProfilerResults(export_path,graphic_links)
    except Exception:
        print traceback.format_exc(),'\n'
        ### Analysis may not be supported for species or data is incompatible
        try:
            if remoteAnalysis:
                if species != 'Mm' and species != 'Hs':
                    print 'LineageProfiler analysis failed (possibly unsupported species).'
                #print traceback.format_exc(),'\n'
        except Exception:
            pass
    return graphic_links

def combineLPResultFiles(input_files):
    combined_sample_cell_db={}
    celltypes=[]
    for fn in input_files:
        if fn != None:
            firstLine=True
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t=string.split(data,'\t')
                if firstLine:
                    headers = t[1:]
                    firstLine=False
                else:
                    cell_type = t[0]
                    if cell_type not in celltypes:
                        celltypes.append(cell_type)
                    zscores = t[1:]
                    for sample in headers:
                        z = float(zscores[headers.index(sample)])
                        try:
                            cell_zscores = combined_sample_cell_db[sample]
                            try: cell_zscores[cell_type].append(z)
                            except Exception: cell_zscores[cell_type]=[z]
                        except Exception:
                            combined_sample_cell_db[sample] = {cell_type:[z]}
    try:
        headers.sort()
        celltypes.sort()
        for i in input_files:
            if i!=None:
                output_file = string.join(string.split(i,'-')[:-2],'-')+'-zscores.txt'
                o = export.ExportFile(output_file)
                o.write(string.join(['LineagePredictions']+headers,'\t')+'\n')
                break
        for cell_type in celltypes:
            values = [cell_type]
            for sample in headers:
                cell_zscores = combined_sample_cell_db[sample]
                #try: cell_zscores[cell_type].sort()
                #except Exception: cell_zscores[cell_type] = [0]
                selectedZ=str(cell_zscores[cell_type][0])
                #if 'Breast' in sample and cell_type=='Breast': print cell_zscores['Breast'],sample, selectedZ;sys.exit()
                #selectedZ=str(statistics.avg(cell_zscores[cell_type]))
                values.append(selectedZ)
            o.write(string.join(values,'\t')+'\n')
        o.close()
    except Exception: pass
    
    try: returnRowHeaderForMaxEntry(output_file,10)
    except Exception: pass
    return output_file

def visualizeQCPlots(expr_input_dir):
    original_expr_input_dir = expr_input_dir
    expr_input_dir = string.replace(expr_input_dir,'-steady-state','')  ### We want the full exon/probeset-level expression file
    try:
        from visualization_scripts import clustering
        from visualization_scripts import QC

        print 'Building quality control graphs...'
        if array_type == 'RNASeq':
            counts_input_dir = string.replace(expr_input_dir,'exp.','counts.')
            graphic_links = QC.outputRNASeqQC(counts_input_dir)
        else:
            graphic_links = QC.outputArrayQC(expr_input_dir)
        
        print 'Building hierarchical cluster graphs...'
        paths = getSampleLogFoldFilenames(expr_input_dir)
        graphic_links = clustering.outputClusters(paths,graphic_links, Normalize='median',Species=species)
        try: graphic_links = clustering.runPCAonly(original_expr_input_dir,graphic_links,False,plotType='2D',display=False)
        except Exception: pass
    except Exception:
        print 'Unable to generate QC plots:'
        print traceback.format_exc()
        try: graphic_links = graphic_links
        except Exception: graphic_links=None ### Matplotlib likely not installed - or other unknown issue
    return graphic_links

def getSampleLogFoldFilenames(filename):
    if '/' in filename: delim = '/'
    else: delim = '\\'
    if 'ExpressionInput' in filename:
        export_path = string.replace(filename,'ExpressionInput','ExpressionOutput/Clustering')
        path1 = string.replace(export_path,'exp.','SampleLogFolds-')
        path2 = string.replace(export_path,'exp.','OutlierLogFolds-')
        path3 = string.replace(export_path,'exp.','RelativeSampleLogFolds-')
        paths = [path1,path2,path3]
    else:
        paths = string.split(filename,delim)
        path1 = string.join(paths[:-1],delim)+'/ExpressionOutput/Clustering/SampleLogFolds-'+paths[-1]
        path2 = string.replace(path1,'SampleLogFolds-','OutlierLogFolds-')
        path3 = string.replace(path1,'SampleLogFolds-','RelativeSampleLogFolds-')
        paths = [path1,path2,path3]
    return paths

def importGeneAnnotations(species):
    ### Used for internal testing
    gene_annotation_file = "AltDatabase/ensembl/"+species+"/"+species+"_Ensembl-annotations_simple.txt"
    fn=filepath(gene_annotation_file)
    count = 0; annotate_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if count == 0: count = 1
        else:
            gene, description, symbol = string.split(data,'\t')
            annotate_db[gene] = symbol,description,''
    return annotate_db

def remoteExpressionBuilder(Species,Array_type,dabg_p,expression_threshold,
                            avg_all_for_ss,Expression_data_format,Vendor,
                            constitutive_source,data_source,Include_raw_data,
                            perform_alt_analysis,GE_fold_cutoffs,GE_pvalue_cutoffs,
                            GE_ptype,exp_file_location_db,Root):
    start_time = time.time()
    global root; root = Root
    #def remoteExpressionBuilder():
    global species; global array_type ; species = Species; array_type = Array_type; global altanalyze_files; global vendor; vendor = Vendor
    global filter_by_dabg; filter_by_dabg = 'yes' ### shouldn't matter, since the program should just continue on without it
    global expression_data_format; global expression_dataset_output_dir; global root_dir; global data_type
    global conventional_array_db; global custom_array_db; global constitutive_db; global include_raw_data; global experiment_name
    global annotate_db; global probeset_db; global process_custom; global m_cutoff; global p_cutoff; global ptype_to_use; global norm
    global arrayCode; arrayCode = 0; global probability_statistic; global fl; global use_downregulated_label; use_downregulated_label = True

    global count_statistics_db; global count_statistics_headers; count_statistics_db = {}
    include_raw_data = Include_raw_data; expression_data_format = Expression_data_format
    data_type = 'expression' ###Default, otherwise is 'dabg'
    d = "core"; e = "extendend"; f = "full"; exons_to_grab = d ### Currently, not used by the program... intended as an option for ExonArrayAffymetrixRules full annotation (deprecated)
  
    ### Original options and defaults
    """
    dabg_p = 0.75; data_type = 'expression' ###used for expression analysis when dealing with AltMouse arrays
    a = "3'array"; b = "exon"; c = "AltMouse"; e = "custom"; array_type = c
    l = 'log'; n = 'non-log'; expression_data_format = l
    w = 'Agilent'; x = 'Affymetrix'; y = 'Ensembl'; z = 'default'; data_source = y; constitutive_source = z; vendor = x
    hs = 'Hs'; mm = 'Mm'; dr = 'Dr'; rn = 'Rn'; species = mm
    include_raw_data = 'yes'  
    expression_threshold = 70 ### Based on suggestion from BMC Genomics. 2006 Dec 27;7:325. PMID: 17192196, for hu-exon 1.0 st array
    avg_all_for_ss = 'no'  ###Default is 'no' since we don't want all probes averaged for the exon arrays
    """

    ct = 'count'; avg = 'average'; filter_method = avg
    filter_by_dabg = 'yes'

    m_cutoff = m_cutoff = math.log(float(GE_fold_cutoffs),2); p_cutoff = float(GE_pvalue_cutoffs); ptype_to_use = GE_ptype
  
    if array_type=="3'array":
        platform_description = "gene-level"
    else:
        platform_description = array_type
    print "Beginning to process the",species,platform_description,'dataset'
  
    process_custom = 'no'  
    if array_type == "custom": ### Keep this code for now, even though not currently used
        import_dir = '/AltDatabase/affymetrix/custom'
        dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
        for affy_data in dir_list:    #loop through each file in the directory to output results
            affy_data_dir = 'AltDatabase/affymetrix/custom/'+affy_data
            custom_array_db = parse_custom_annotations(affy_data_dir)
            array_type = a; process_custom = 'yes'
    if array_type == "AltMouse":
        print "Processing AltMouse splicing data"
        original_probeset_db,constitutive_db = importAltMerge('basic')
        probe_annotation_file = "AltDatabase/"+species+'/'+ array_type+'/'+array_type+"_annotations.txt"
        original_annotate_db = import_annotations(probe_annotation_file)
        conventional_array_db = []
    elif array_type == "3'array" and 'Ensembl' not in vendor: ### If user supplied IDs are from Ensembl - doesn't matter the vendor
        original_vendor = vendor
        if 'other:' in vendor:
            vendor = string.replace(vendor,'other:','')
        process_go='yes';extract_go_names='yes';extract_pathway_names='yes'
        probeset_db = []; annotate_db = []
        constitutive_db = ""; conventional_array_db = {}
        affy_data_dir = 'AltDatabase/affymetrix'
        if vendor == 'Affymetrix':
            try: conventional_array_db, arrayCode = BuildAffymetrixAssociations.importAffymetrixAnnotations(affy_data_dir,species,process_go,extract_go_names,extract_pathway_names)
            except Exception: print 'Error in processing CSV data. Getting this data from GO-Elite annotations instead.'
        if vendor == 'Affymetrix' and len(conventional_array_db)>0: use_go='no'
        else: use_go = 'yes'
        try:
            print "Adding additional gene, GO and WikiPathways annotations"
            conventional_array_db = BuildAffymetrixAssociations.getUIDAnnotationsFromGOElite(conventional_array_db,species,vendor,use_go)
        except Exception: print "Additional annotation import failed"
        print len(conventional_array_db), "Array IDs with annotations from",vendor,"annotation files imported."
        vendor = original_vendor
    elif array_type != "AltMouse":
        probeset_db = []; annotate_db = []; constitutive_db = []; conventional_array_db = []
        ### The below function gathers GO annotations from the GO-Elite database (not Affymetrix as the module name implies)
        conventional_array_db = BuildAffymetrixAssociations.getEnsemblAnnotationsFromGOElite(species)
    if 'Ensembl' in vendor:
        annotate_db = importGeneAnnotations(species) ### populate annotate_db - mimicking export structure of exon array

    global expr_threshold; global dabg_pval; global gene_exp_threshold; global gene_rpkm_threshold; dabg_pval = dabg_p
    altanalyze_files = []; datasets_with_all_necessary_files=0
    for dataset in exp_file_location_db:
        experiment_name = string.replace(dataset,'exp.',''); experiment_name = string.replace(experiment_name,'.txt','')
        fl = exp_file_location_db[dataset]
        expr_input_dir = fl.ExpFile()
        stats_input_dir = fl.StatsFile()
        expr_group_dir = fl.GroupsFile()
        comp_group_dir = fl.CompsFile()
        try: batch_effects = fl.BatchEffectRemoval()
        except Exception: batch_effects = 'NA'
        try: norm = fl.FeatureNormalization()
        except Exception: norm = 'NA'
      
        try: probability_statistic = fl.ProbabilityStatistic()
        except Exception: probability_statistic = 'unpaired t-test'
        try: gene_exp_threshold = fl.GeneExpThreshold()
        except Exception: gene_exp_threshold = 0
        try: gene_rpkm_threshold = fl.RPKMThreshold()
        except Exception: gene_rpkm_threshold = 0
        if expression_data_format == 'log':
            try: expr_threshold = math.log(float(expression_threshold),2)
            except Exception: expr_threshold = 0 ### Applies to RNASeq datasets
        else:
            try: expr_threshold = float(expression_threshold)
            except Exception: expr_threshold = 0
        
        residuals_input_dir = string.replace(expr_input_dir,'exp.','residuals.')
        root_dir = fl.RootDir()
        datasets_with_all_necessary_files +=1
        checkArrayHeaders(expr_input_dir,expr_group_dir)
        expression_dataset_output_dir = root_dir+"ExpressionOutput/"
        if batch_effects == 'yes':
            try:
                from stats_scripts import combat
                combat.runPyCombat(fl)
            except Exception:
                print_out = 'Batch effect removal analysis (py-combat) failed due to an uknown error:'
                print traceback.format_exc()
                try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); root.destroy(); sys.exit()
                except Exception: print print_out; sys.exit()
        
        if array_type != "3'array": #array_type != 'AltMouse' and 
            try: probeset_db,annotate_db,comparison_filename_list = ExonArray.getAnnotations(fl,array_type,dabg_p,expression_threshold,data_source,vendor,constitutive_source,species,avg_all_for_ss,filter_by_dabg,perform_alt_analysis,expression_data_format)
            except Exception, e:
                print traceback.format_exc()
                print_out = 'Error ecountered for the '+species+', '+array_type+' dataset. Check to ensure that:\n(1) the correct platform and species were selected and\n(2) some expression values are present in ExpressionInput/exp.YourDataset.txt'
                try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); root.destroy(); sys.exit()
                except Exception: print print_out; sys.exit()
            if array_type != 'AltMouse': expr_input_dir = expr_input_dir[:-4]+'-steady-state.txt'
            else: probeset_db = original_probeset_db; annotate_db = original_annotate_db
            for file in comparison_filename_list: altanalyze_files.append(file)
            residual_file_status = ExonArray.verifyFile(residuals_input_dir)
            ### Separate residual file into comparison files for AltAnalyze (if running FIRMA)
            if residual_file_status == 'found': ExonArray.processResiduals(fl,Array_type,Species,perform_alt_analysis)
        """
        from build_scripts import ExonArrayEnsemblRules
        source_biotype = array_type, root_dir
        probeset_db,annotate_db,constitutive_gene_db,splicing_analysis_db = ExonArrayEnsemblRules.getAnnotations('no',constitutive_source,source_biotype,species)
        expr_input_dir = expr_input_dir[:-4]+'-steady-state.txt'
        """
        if norm == 'RPKM' and array_type == 'RNASeq':
            ### Separately analyze steady-state counts first, to replace fold changes
            counts_expr_dir = string.replace(expr_input_dir,'exp.','counts.')
            if 'counts.' not in counts_expr_dir: counts_expr_dir = 'counts.'+counts_expr_dir ### Occurs if 'exp.' not in the filename
            count_statistics_db, count_statistics_headers = calculate_expression_measures(counts_expr_dir,expr_group_dir,experiment_name,comp_group_dir,probeset_db,annotate_db)
        
        calculate_expression_measures(expr_input_dir,expr_group_dir,experiment_name,comp_group_dir,probeset_db,annotate_db)
        buildCriterion(GE_fold_cutoffs, p_cutoff, ptype_to_use, root_dir+'/ExpressionOutput/','summary') ###Outputs a summary of the dataset and all comparisons to ExpressionOutput/summary.txt
        #except Exception: null=[]
        graphic_links = None
        if fl.ProducePlots() == 'yes':
            graphic_links = visualizeQCPlots(expr_input_dir)
        if fl.PerformLineageProfiler() == 'yes':
            if graphic_links==None: graphic_links = []
            graphic_links = performLineageProfiler(expr_input_dir,graphic_links) ### Correlate gene-level expression values with known cells and tissues
        if graphic_links != None:
            fl.setGraphicLinks(graphic_links) ### Uses Matplotlib to export QC and clustering plots
        
    annotate_db={}; probeset_db={}; constitutive_db={}; array_fold_db={}; raw_data_comps={}; conventional_array_db=[]
    clearObjectsFromMemory(conventional_array_db); conventional_array_db=[]
    try: clearObjectsFromMemory(summary_filtering_stats); summary_filtering_stats=[]
    except Exception: null=[]
    try: clearObjectsFromMemory(array_folds); array_folds=[]
    except Exception: null=[]
    try: clearObjectsFromMemory(count_statistics_db); count_statistics_db=[]
    except Exception: null=[]
  
    #print 'after deleted'; returnLargeGlobalVars()
    
    ### Generate the NI file if possible
    try: calculateNormalizedIntensities(root_dir,species,array_type,avg_all_for_SS=avg_all_for_ss)
    except Exception: pass
    
    if datasets_with_all_necessary_files == 0:
        ###Thus no files were found with valid inputs for all file types
        print 'WARNING....No propperly named datasets were found. ExpressionBuilder requires that there are at least 3 files with the prefixes "exp.", "groups." and "comps.", with the following dataset name being identical with all three files.'
        print "...check these file names before running again."
        inp = sys.stdin.readline(); sys.exit()
    altanalyze_files = unique.unique(altanalyze_files) ###currently not used, since declaring altanalyze_files a global is problematic (not available from ExonArray... could add though)
    if array_type != "3'array" and perform_alt_analysis != 'expression':
        from stats_scripts import FilterDabg; reload(FilterDabg)
        altanalyze_output = FilterDabg.remoteRun(fl,species,array_type,expression_threshold,filter_method,dabg_p,expression_data_format,altanalyze_files,avg_all_for_ss)
        return 'continue',altanalyze_output
    else:
        end_time = time.time(); time_diff = int(end_time-start_time)
        return 'stop'

def verifyFile(filename):
    fn=filepath(filename)
    try:
        for line in open(fn,'rU').xreadlines(): found = 'yes'; break
    except Exception: found = 'no'
    return  found

def verifyExpressionFile(filename):
    """ Unlike the above, corrects the expression file path if not found """
    fn=filepath(filename)
    try:
        for line in open(fn,'rU').xreadlines(): break
    except Exception:
        fn = string.replace(fn,'ExpressionInput/','') ### This file is in the parent path presumably (otherwise can't find it)
    return fn

def exportSignatures(db,directory,species):
    import gene_associations
    export_file = export.findParentDir(directory[:-2]+'.txt')
    if 'AltExon' in directory:
        export_file+='signatures/AltExon-signatures.txt'
    else:
        export_file+='signatures/GeneExpression-signatures.txt'
    export_object = export.ExportFile(export_file)
    if 'AltExon' in directory:
        export_object.write('symbol\tentrez\tensembl\tsource\turl\tname\tAltAnalyze-ExonID\tASPIRE|Splicing-Index LogFold\tGenomicLocation\n') ### Header line
    else:
        export_object.write('symbol\tentrez\tsource\turl\tname\tLogFold\tBH-adjp-value\n') ### Header line
    url = 'http://code.google.com/p/altanalyze/wiki/PCBC_C4_compendium'
    source = 'AltAnalyze'
    sig = 0
    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    for filename in gene_conversion_db:
        db,input_data_db = gene_conversion_db[filename]
        filename = string.replace(filename,'.txt','')
        filename = string.replace(filename,'GE.','')
        filename = string.replace(filename,'-upregulated','')
        for ensembl in db:
            sig+=1
            try: symbol = gene_to_symbol[ensembl][0]
            except Exception: continue
            entrezgenes = db[ensembl][0]
            entrezgenes = string.split(entrezgenes,'|')
            statistics = input_data_db[ensembl][0][1:]
            #print [statistics];sys.exit()
            for entrez in entrezgenes:
                ### output format: symbol, entrez, source, url, name
                values = string.join([symbol,entrez,ensembl,source,url,filename]+statistics,'\t')+'\n'
                export_object.write(values)
    export_object.close()
    print 'Exported',sig,'signatures to:'
    print export_file
    
def runGOElite(species,directory):
    """ Separate pipeline for automating GO-Elite when re-generating criterion - Currently used outside of any pipelines """
    mod = 'Ensembl'
    pathway_permutations = 'FisherExactTest'
    filter_method = 'z-score'
    z_threshold = 1.96
    p_val_threshold = 0.05
    change_threshold = 2
    resources_to_analyze = 'local'
    returnPathways = 'yes'
    root = None
    import GO_Elite
    directory = string.replace(directory,'ExpressionOutput','')
    results_dir = directory
    print '\nBeginning to run GO-Elite analysis on all results'
    
    elite_input_dirs = ['regulated']#,'upregulated','downregulated','MarkerFinder'] ### 'AltExon' Run GO-Elite multiple times to ensure heatmaps are useful and to better organize results
    for elite_dir in elite_input_dirs:
        if elite_dir == 'AltExon': returnPathways = 'no'
        else: returnPathways = 'yes'
        file_dirs = results_dir+'GO-Elite/'+elite_dir,results_dir+'GO-Elite/denominator',results_dir+'GO-Elite/'+elite_dir
        input_dir = results_dir+'GO-Elite/'+elite_dir
        try: input_files = read_directory(input_dir) ### Are there any files to analyze?
        except Exception: input_files = []
        if len(input_files)>0:
            variables = species,mod,pathway_permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,returnPathways,file_dirs,root
            try: GO_Elite.remoteAnalysis(variables,'non-UI',Multi=mlp)
            except Exception: pass

def filterDatasetFile(main_output_folder):
    global array_folds; global expression_dataset_output_dir
    expression_dataset_output_dir = string.replace(main_output_folder,'GO-Elite','ExpressionOutput/')
    dir_list = read_directory(expression_dataset_output_dir[:-1])
    for filename in dir_list:
        if 'DATASET-' in filename: 
            dataset_fn=filepath(expression_dataset_output_dir+filename)
            x=0
            for line in open(dataset_fn,'rU').xreadlines():
                data = cleanUpLine(line); t = string.split(data,'\t')
                if x==0: headers = t; break

    ### Get a list of column header titles to include
    fn=filepath('Config/DATASET-headers.txt')
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line); t = string.split(data,'\t')
        columns_to_include = t
        
        
    ### Filter statistics based on user-defined thresholds as input for GO-Elite analysis
    criterion_db={}; denominator_geneids={}; index = 0; indexes=[]; avg_indexes=[]
    for column in headers:
        if 'avg-' in column: avg_indexes.append(index)
        elif 'log_fold-' in column: indexes.append(index)
        elif 'fold-' in column: indexes.append(index)
        elif 'rawp-' in column: indexes.append(index)
        elif 'adjp-' in column: indexes.append(index)
        elif 'ANOVA' in column: indexes.append(index)
        elif column in columns_to_include: indexes.append(index)
        index+=1
        
    ###Export out the filtered file
    dataset_fn_filtered = string.replace(dataset_fn,'.txt','-abreviated.txt')
    dataset_fn_filtered = string.replace(dataset_fn_filtered,'ExpressionOutput','ExpressionOutput/filtered')
    data_filter = export.ExportFile(dataset_fn_filtered)
    firstLine=True
    for line in open(dataset_fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            firstLine=False
            if 'maximum sample read count' in t: ### For RNA-Seq
                h=[]
                for l in t:
                    if 'avg-' in l: h.append(string.replace(l,'avg-','avg_RPKM-'))
                    else: h.append(l)
                t=h
        values = map(lambda x: t[x], indexes)
        avg_values = map(lambda x: t[x], avg_indexes) ### Specifically grab only the average-expression values (not anything else)
        values = string.join(values+avg_values,'\t')+'\n'
        data_filter.write(values)
    data_filter.close()

def importProbesetRegions(species,platform):
    filename = 'AltDatabase/'+species+'/'+platform+'/'+species+'_Ensembl_probesets.txt'
    fn=filepath(filename)
    region_db = {}
    firstRow=True
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow: firstRow = False
        else:
            probeset = t[0]
            gene = t[2]
            region = gene+':'+string.replace(t[12],'-','.')
            region_db[region] = probeset
    return region_db

def buildAltExonClusterInputs(input_folder,species,platform,dataType='AltExonConfirmed'):
    alternative_exon_db={}
    dir_list = read_directory(input_folder)
    
    if platform == 'junction':
        region_db = importProbesetRegions(species,platform)
        
    for filename in dir_list:
        if '.txt' in filename:
            proceed = True
            if platform == 'RNASeq':
                if 'ASPIRE' in filename or 'egress' in filename: proceed = False ### Need to create a special analysis just for reciprocal junctions
            if proceed: ### Don't include splicing-index results for RNA-Seq
                fn=filepath(input_folder+filename)
                x=0
                for line in open(fn,'rU').xreadlines():
                    data = cleanUpLine(line)
                    t = string.split(data,'\t')
                    if x==0: x=1
                    else:
                        if platform == 'junction' and 'GO-Elite' in input_folder:
                            altExon = t[2]
                            if altExon in region_db:
                                altExon = region_db[altExon]
                        elif 'splicing-index' in filename:
                            altExon = t[-1]
                        elif 'ASPIRE' in filename or 'egress' in filename:
                            altExon = t[-1]
                            altExon = string.split(altExon,' | ')
                        else:
                            altExon = t[2]
                        #if float(t[3])<0.001:
                        alternative_exon_db[altExon]=None
    print len(alternative_exon_db), 'alternative exon IDs imported'
    import gene_associations
    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    
    input_folder = string.split(input_folder,'GO-Elite')[0]+'AltResults/RawSpliceData/'+species+'/splicing-index/'
    dir_list = read_directory(input_folder)
    exported_IDs=0
    added={}
    for filename in dir_list:
        if '.txt' in filename and ('_average' not in filename) or len(dir_list)==1: ### We only want the group-comparison SI file
            export_dir = string.split(input_folder,'RawSpliceData')[0]+'Clustering/'+dataType+'-'+filename
            export_data = export.ExportFile(export_dir)
            fn=filepath(input_folder+filename)
            x=0
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if x==0:
                    headers = t[2:]  ### first two columsn are gene and ExonID
                    export_data.write(string.join(headers,'\t')+'\n') ### write header row
                    x=1
                else:
                    if platform != 'RNASeq':
                        if 'ENS' in t[0]: #ENSG human
                            gene = t[0]; altExon = t[2]
                        else:
                            altExon = t[0]
                    else:
                        gene = t[0]; altExon = t[2]
                    if ';' in altExon:
                        altExon1, altExon2 = string.split(altExon,';')
                        altExons = [altExon1,gene+':'+altExon2]
                    else:
                        altExons = [altExon]
                    for altExon in altExons:
                        if altExon in alternative_exon_db and altExon not in added:
                            added[altExon]=[]
                            #values = map(lambda x: float(x)*-1, t[3:]) #reverse the fold for cluster visualization
                            values = map(lambda x: float(x), t[3:])
                            #print len(headers),len(values);sys.exit()
                            avg = statistics.avg(values)
                            log_folds = map(lambda x: str(x-avg), values) ### calculate log folds and make these strings
                            values = string.join([altExon]+log_folds,'\t')+'\n' ### [t[3]]+ before log_folds?
                            if gene in gene_to_symbol: symbol = gene_to_symbol[gene][0]+"  "
                            else: symbol = ''
                            export_data.write(symbol+values)
                            exported_IDs+=1
    
    print exported_IDs, 'exported ID values for clustering'
    export_data.close()
    return export_dir, exported_IDs
                    
def exportHeatmap(filename,useHOPACH=True, color_gradient='red_black_sky',normalize=False,columnMethod='average',size=0,graphics=[]):
    from visualization_scripts import clustering
    row_method = 'weighted'; row_metric = 'cosine'; column_method = 'average'; column_metric = 'euclidean'; transpose = False
    try:
        if columnMethod !=None:
            column_method = columnMethod
        if size < 7000:
            graphics = clustering.runHCexplicit(filename, graphics, row_method, row_metric, column_method, column_metric, color_gradient, transpose, display=False, Normalize=normalize)
    except Exception:
        print 'Clustering failed for:',filename
    return graphics

def meanCenterPSI(filename):
    firstLine=True
    output_file = filename[:-4]+'-cluster.txt'
    export_obj = export.ExportFile(output_file)
    for line in open(filename,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            export_obj.write(line)
            firstLine = False
        else:
            try:
                avg = averageWithNulls(t[1:])
                values = map(lambda x: replaceNulls(x,avg), t[1:])
                export_obj.write(string.join([t[0]]+values,'\t')+'\n')  
            except Exception: pass
    return output_file

def filterJunctionExpression(filename,minPercentPresent=None):
    output_file = filename[:-4]+'-filter.txt'
    export_obj = export.ExportFile(output_file)
    filtered_lines = []; filtered_size=0; filtered_size_stringent=0
    firstLine = True
    size=0; imported_num=0
    for line in open(filename,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            export_obj.write(line)
            firstLine = False
        else:
            imported_num+=1
            if '_' in t[0] or '-ENS' in t[0] or '-' not in t[0]:
                pass
            else:
                try: a,b = string.split(t[0],'|')
                except Exception:
                    try: n,a,b = string.split(t[0],' ')
                    except Exception: continue
                if '-' in a and '-' in b:
                    vals = map(lambda x: countNulls(x),t[1:])
                    percent_present = sum(vals)/(len(vals)*1.00)
                    if percent_present>0.1:
                        filtered_lines.append([percent_present,line])
                        size+=1
                        if percent_present>0.5:
                            filtered_size+=1
                            if percent_present>0.85:
                                filtered_size_stringent+=1
                            
    percent_imported =  ((1.00*size)/imported_num)
    if minPercentPresent!=None:
        size=0
        filtered_lines.sort(); filtered_lines.reverse()
        for (percent_present,line) in filtered_lines:
            if percent_present>minPercentPresent:
                    export_obj.write(line)
                    size+=1
    elif size < 8000 and percent_imported<0.5: ### The filtered is at least 30% the size of the imported
        for (percent_present,line) in filtered_lines:
            export_obj.write(line)
    else:
        size=0
        to_import = int(imported_num*0.3)
        #print "Filtering down to", to_import, "entries."
        filtered_lines.sort(); filtered_lines.reverse()
        for (percent_present,line) in filtered_lines:
            if filtered_size_stringent>8000:
                if percent_present>0.95:
                    export_obj.write(line)
                    size+=1
            elif filtered_size>8000:
                if percent_present>0.85:
                    export_obj.write(line)
                    size+=1
            else:   
                if percent_present>0.5:
                    export_obj.write(line)
                    size+=1
                
    print 'Filtered RI-PSI entries to',size
    export_obj.close()
    return output_file,size

def countNulls(val):
    if float(val) == 0:
        return 0
    else: return 1

def calculateNormalizedIntensities(root_dir, species, array_type, avg_all_for_SS = 'no', probeset_type = 'core', analysis_type = 'processed', expFile = False):
    """ Since it is a time-consuming step that is needed for visualizing SI values for exons, it is faster
    to efficiently calculate NI values in advance. """
    if array_type == 'gene': platform = 'GeneArray'
    elif  array_type == 'exon': platform = 'ExonArray'
    elif  array_type == 'junction': platform = 'JunctionArray'
    else: platform = array_type
    #alt_exon_exp_dir = root_dir+'/AltExpression/FullDatasets/'+platform+'/'+species ### This file doesn't exist if only one pairwise comparison group
    alt_exon_exp_dir = root_dir+'/ExpressionInput/' ### This file doesn't exist if only one pairwise comparison group
    dir_list = read_directory(alt_exon_exp_dir)
    for file in dir_list:
        if '.txt' in file and 'exp.' in file and 'steady' not in file:
            selected_file = file[4:]
            fn=filepath(alt_exon_exp_dir+'/'+file)

    #sample_group_db = simplerGroupImport(string.replace(fn,'exp.','groups.'))
    
    if analysis_type == 'raw':
        ### Create a filtered exon and junction expression file outside the typical workflow
        import RNASeq
        exp_threshold=5; rpkm_threshold=5
        expressed_uids_rpkm = RNASeq.getMaxCounts(expFile,rpkm_threshold)
        expressed_uids_counts = RNASeq.getMaxCounts(string.replace(expFile,'exp.','counts.'),exp_threshold)
        expressed_uids = expressed_uids_rpkm.viewkeys() & expressed_uids_counts.viewkeys() ### common
        fn = root_dir+'/AltExpression/FilteredDataset/'+platform+'/'+species+'/'+ export.findFilename(expFile)### This file doesn't exist if only one pairwise comparison group
        expressed_uids_rpkm = RNASeq.getMaxCounts(expFile,rpkm_threshold,filterExport=expressed_uids,filterExportDir=fn)
        
    ### Get the gene values
    gene_db={};
    if analysis_type == 'raw':
        ### Get these directly from the steady-state file
        exp_dir=string.replace(expFile,'.txt','-steady-state.txt')
        firstLine = True
        low_diff_exp_genes={}
        for line in open(exp_dir,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if firstLine:
                firstLine = False
                ge_header = t
            else:
                values = map(lambda x: math.log(float(x),2), t[1:])
                gene_db[t[0]]=values

    elif array_type == 'RNASeq':
        firstLine=True
        fn_steady = fn[:-4]+'-steady-state.txt'
        for line in open(fn_steady,'rU').xreadlines():
            if firstLine: firstLine = False
            elif ':' in line[:50]: pass
            else:
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                try:
                    values = map(lambda x: float(x), t[1:])
                    values = map(lambda x: math.log(x,2),values)
                    gene_db[t[0]]=values
                except Exception: pass
    else:
        import AltAnalyze
        exon_db, constitutive_probeset_db = AltAnalyze.importSplicingAnnotations(array_type,species,probeset_type,avg_all_for_SS,root_dir)
        gene_db={}; firstLine=True
        for line in open(fn,'rU').xreadlines():
            if firstLine: firstLine = False
            else:
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if t[0] in constitutive_probeset_db:
                    gene = constitutive_probeset_db[t[0]]
                    values = map(lambda x: float(x), t[1:])
                    try: gene_db[gene].append(values) ### combine these
                    except Exception: gene_db[gene] = [values]
        for gene in gene_db:
            #print gene, gene_db[gene]
            constitutive_values = zip(*gene_db[gene]) ### combine values from multiple lists into lists of lists
            #print constitutive_values
            values = map(lambda x: Average(x), constitutive_values)
            #print values; sys.exit()
            gene_db[gene] = values
         
    if len(gene_db)>0:
        alt_exon_output_dir = root_dir+'/AltResults/RawSpliceData/'+species+'/splicing-index/'+selected_file
        if analysis_type == 'raw':
            alt_exon_output_dir = string.replace(alt_exon_output_dir,'RawSpliceData','RawSpliceDataTemp')
        export_obj = export.ExportFile(alt_exon_output_dir)
        print 'Exporting exon-level Normalized Intensity file to:',alt_exon_output_dir
        print len(gene_db),'genes with data imported'
        
        firstLine=True
        exon_entries=0; saved_entries=0
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if firstLine:
                firstLine = False
                values = string.join(['Gene','ExonID','probesetID']+t[1:],'\t')+'\n'
                export_obj.write(values)
            else:
                if (':' in t[0]) or array_type != 'RNASeq':
                    feature_values = map(lambda x: float(x), t[1:])
                    if analysis_type == 'raw' or array_type == 'RNASeq':
                        feature_values = map(lambda x: math.log(x,2), feature_values)
                    if array_type == 'RNASeq':
                        gene = string.split(t[0],':')[0]
                    else:
                        try: gene = exon_db[t[0]].GeneID()
                        except Exception: gene = None
                    if gene != None:
                        exon_entries+=1
                        if gene in gene_db:
                            gene_values = gene_db[gene]       
                            if ('-' not in t[0]) or analysis_type == 'raw':
                                NI_values = [logratio(value) for value in zip(*[feature_values,gene_values])]
                                NI_values = map(lambda x: str(x), NI_values)
                                NI_values = string.join([gene,'NA',t[0]]+NI_values,'\t')+'\n'
                                export_obj.write(NI_values)
                                saved_entries+=1
        export_obj.close()
        print exon_entries, 'found with',saved_entries,'entries normalized.'
    return alt_exon_output_dir

def compareRawJunctionExpression(root_dir,platform,species,critical_exon_db,expFile,min_events=0,med_events=0):
    expFile = exportSorted(expFile, 0) ### Sort the expression file
    print expFile
    from  scipy import stats
    exported=[]
    retained_introns=[]
    features_examined={}
    junction_locations = {}
    critical_exon_gene_db={}
    critical_junction_pair_db={}
    psi_db={}
    #min_events = 4; med_events = 13
    #min_events = 0; med_events = 0
    
    print 'Begining reciprocal junction/intron retention unbiased analyses (min_events=%d, med_events=%d)' % (min_events,med_events)
    for altexon in critical_exon_db:
        gene = string.split(altexon,':')[0]
        inclusion_list,exclusion_list = critical_exon_db[altexon]
        try:
            critical_exon_gene_db[gene].append(critical_exon_db[altexon])
        except Exception:
            critical_exon_gene_db[gene] = [critical_exon_db[altexon]]
        
        for i in inclusion_list:
            for e in exclusion_list:
                try: critical_junction_pair_db[i,e].append(altexon)
                except Exception: critical_junction_pair_db[i,e] = [altexon]
            
    export_dir = root_dir+'AltResults/AlternativeOutput/'+species+'_'+platform+'_top_alt_junctions-PSI.txt'
    export_data = export.ExportFile(export_dir)
    
    clust_export_dir = root_dir+'AltResults/AlternativeOutput/'+species+'_'+platform+'_top_alt_junctions-PSI-clust.txt'
    clust_export_data = export.ExportFile(clust_export_dir)

    ### Included a new file with cluster and junction information
    clusterID_export_dir=root_dir+'AltResults/AlternativeOutput/'+species+'_'+platform+'_top_alt_junctions-PSI-ClusterIDs.txt'
    clusterID_export_data = export.ExportFile(clusterID_export_dir)

    def ratio(values):
        #print values
        #return float(values[0]+1)/(values[1]+1)
        try:return float(values[0])/(values[1])
        except Exception: return 0.0
    
    def dPSI(incl_exp,total_exp,indexes,min_events):
        indexes.sort()
        incl_filt = map(lambda x: incl_exp[x],indexes) ### Only consider expressed cells
        total_filt = map(lambda x: total_exp[x],indexes)
        dpsi_values = [ratio(value) for value in zip(*[incl_filt,total_filt])]
        dpsi_values.sort()
        if len(dpsi_values)==1:
            return dpsi_values[0]
        elif min_events == 0:
            return dpsi_values[-1]-dpsi_values[0]
        else: return dpsi_values[-1*min_events]-dpsi_values[min_events]
        #return max(dpsi_values)-min(dpsi_values) ### percent change in isoform expression
    
    def getIndexes(events,min_exp_thresh):
        i=0
        indexes=[]
        for e in events:
            if e>min_exp_thresh: indexes.append(i) ### minimum expression value (changed from 5 to 10 8/5/2016)
            i+=1
        return indexes
    
    def diff(values,demon,num):
        return values[num]-values[demon]
    
    def diffCompare(incl_exp,excl_exp,incl,excl):
        if max(incl_exp)>max(excl_exp): denom = 1;num = 0
        else: denom = 0; num=1
        diff_values = [diff(value,denom,num) for value in zip(*[incl_exp,excl_exp])]
        if min(diff_values)==0 and ('-' not in incl or '-' not in excl): ### subtract out the overlap if all reads in the inclusion exon are confounded by the junction 1nt reads
            if denom == 1:
                feature_exp_db[incl] = diff_values
            else:
                feature_exp_db[excl] = diff_values
            return True
        else: return False
            
    def junctionComparisonMethod(incl,excl):
        useAllJunctionsForExcl=True
        #min_exp = 4; med_exp = 9 #min_exp = 19; med_exp = 19
        min_exp=9;med_exp=19
        incl_exp = feature_exp_db[incl]
        excl_exp = feature_exp_db[excl]
        if useAllJunctionsForExcl:
            gene = string.split(excl,':')[0]
            try: gene_exp = gene_junction_denom[gene]
            except Exception:
                try:
                    gene = string.split(incl,':')[0]
                    gene_exp = gene_junction_denom[gene]
                except Exception: gene_exp=[]
        rho,p = stats.pearsonr(incl_exp,excl_exp)
        ### See if the data for one feature is confounded by another
        ### The example is exon-inclusion with 1nt overlap getting the associated reads from BedTools
        status = diffCompare(incl_exp,excl_exp,incl,excl)
        if status: ### update the difference generated values
            excl_exp = feature_exp_db[excl]
            rho,p = stats.pearsonr(incl_exp,excl_exp)
        #if 'ENSMUSG00000009350:E14.2_87617106-E15.1' in incl: print feature_exp_db[incl],'a'
        num_incl_events = sum(i > min_exp for i in incl_exp)
        num_excl_events = sum(i > min_exp for i in excl_exp)
    
        combined = [sum(value) for value in zip(*[incl_exp,excl_exp])]
        total_number_junctions = sum(i > min_exp for i in combined)
        ### Calculate a delta PSI value for each expressed cell
        combined = [sum(value) for value in zip(*[incl_exp,excl_exp])]
        exp_indexes = getIndexes(combined,med_exp)
        try: dpsi = dPSI(incl_exp,combined,exp_indexes,min_events)
        except Exception: dpsi = 1
        
        dpsi_values = [ratio(value) for value in zip(*[incl_exp,combined])]
        dpsi_values = nullReplace(dpsi_values,combined,min_exp)
        psi_db[incl,excl] = dpsi_values ### optionally output these at the end
        #if incl == 'ENSMUSG00000019505:E3.2_62365772-E3.6_62366457' and excl == 'ENSMUSG00000029162:E7.2-E7.4':
        #dpsi2 = dPSI(excl_exp,combined,exp_indexes)
        #if 'ENSMUSG00000009350:E14.2_87617106-E15.1' in incl: print feature_exp_db[incl],'b'
        if '-' in incl and '-' in excl:
            try: max_ratio = expressionRatio(incl_exp,excl_exp,num_incl_events,num_excl_events) # Checks to see if the predominant isoform is expressed at significantly higher levels
            except Exception: max_ratio = 0
            try:
                max_gene_ratio = max([ratio(value) for value in zip(*[incl_exp,gene_exp])])
            except Exception: max_gene_ratio = 0
        else: max_ratio = 1; max_gene_ratio = 0
        #if 'ENSMUSG00000009350:E14.2_87617106-E15.1' in incl: print feature_exp_db[incl],'c'
        #if num_incl_events > 15 and num_excl_events > 15 and total_number_junctions>30 and max_ratio>0.5: ### ensures high expression of the minor isoform
        #if ((num_incl_events > 15 and num_excl_events > 7) or (num_incl_events > 7 and num_excl_events > 15)) and total_number_junctions>20 and max_ratio>0.3:
        #if 'ENSG00000100650' in incl: print incl,excl,max_ratio,num_incl_events,num_excl_events,dpsi,rho
        if ((num_incl_events > med_events and num_excl_events > min_events) or (num_incl_events > min_events and num_excl_events > med_events)) and total_number_junctions>(min_events*2) and max_ratio>0.1: ### ensures high expression of the minor isoform
            #print rho
            if dpsi > 0.15:
                return max_ratio,num_incl_events,num_excl_events,dpsi,rho,max_gene_ratio,True
            else:
                return max_ratio,num_incl_events,num_excl_events,dpsi,rho,max_gene_ratio,False
        else:
            return max_ratio,num_incl_events,num_excl_events,dpsi,rho,max_gene_ratio,False
    
    def intronComparisonMethod(query):
        upstream_exon = string.replace(query,'I','E')
        if '_' in query:
            intronic_region = string.split(query,'_')[0] ### Unclear why this is needed, but sometimes only an _ intronic region exists for the intron
            if intronic_region not in junction_locations: ### Ensures we are not using this instead of a real intron region (e.g., I2.1_1234 and I2.1)
                upstream_exon = string.replace(intronic_region,'I','E')
        try:
            pos1,pos2 = junction_locations[upstream_exon][0]
            upstream_exon_len = abs(pos2-pos1)
        except Exception: upstream_exon_len = None
        downstream_exon = getDownstreamExon(upstream_exon)
        try:
            pos1,pos2 = junction_locations[downstream_exon][0]
            downstream_exon_len = abs(pos2-pos1)
        except Exception: downstream_exon_len=None
        pos1,pos2 = junction_locations[query][0]
        intron_len = abs(pos2-pos1)
        try:
            up_rpk = max(feature_exp_db[upstream_exon])/float(upstream_exon_len)
            if downstream_exon_len!=None:
                down_rpk = max(feature_exp_db[downstream_exon])/float(downstream_exon_len)
                adjacent_rpk = max([up_rpk,down_rpk]) ### get the most conservative estimate
            adjacent_rpk = up_rpk
            intron_rel_exp = (max(feature_exp_db[query])/float(intron_len))/adjacent_rpk
            return intron_rel_exp
        except Exception:
            return 0
    
    def compareJunctionExpression(gene):
        regulated_junctions={}
        inclusion_max_psi={}
        try: symbol,description = gene_annotations[gene]
        except Exception: symbol='';description=''
        if gene in critical_exon_gene_db:
            critical_exon_list = critical_exon_gene_db[gene]
            critical_exon_list = unique.unique(critical_exon_list)
            for (inclusion_list,exclusion_list) in critical_exon_list:
                inclusion_list = unique.unique(inclusion_list)
                exclusion_list = unique.unique(exclusion_list)
                for incl in inclusion_list:
                    for excl in exclusion_list:
                      if excl != incl:
                        pair = [incl,excl]; pair.sort()
                        features_examined[incl]=[]
                        try:
                            max_ratio,num_incl_events,num_excl_events,dpsi,rho,max_all_psi,proceed = junctionComparisonMethod(incl,excl)
                            inclusion_max_psi[incl] = max_all_psi
                            #if 'ENSMUSG00000027680:E21.1-E22.1' in incl: print incl, excl,max_ratio,num_incl_events,num_excl_events,dpsi,rho,proceed
                            if proceed:
                                """if max_ratio<0.1:
                                    if num_excl_events > num_incl_events:
                                        print max_ratio
                                        print max(incl_exp)
                                        print statistics.median(excl_exp)
                                        print incl_exp
                                        print excl_exp;sys.exit()"""
                                #if 'ENSMUSG00000009350:E14.2_87617106-E15.1' in incl: print feature_exp_db[incl]
                                altexons = unique.unique(critical_junction_pair_db[incl,excl])
                                altexons = string.join(altexons,'|')
                                if num_excl_events > num_incl_events:
                                    #print max_ratio, '\t',gene
                                    regulated_junctions[incl]=rho,excl,num_incl_events,num_excl_events,'incl',dpsi,rho,max_ratio,altexons
                                else:
                                    #print max_ratio, '\t',gene
                                    regulated_junctions[excl]=rho,incl,num_excl_events,num_incl_events,'excl',dpsi,rho,max_ratio,altexons
                                #if 'ENSMUSG00000009350:E14.2_87617106-E15.1' in incl: print feature_exp_db[incl]
                                #if rho < -0.3:
                                #print incl, excl, rho
                                #print incl_exp
                                #print excl_exp
                                #print sum(i > 5 for i in incl_exp)
                                #print sum(i > 5 for i in excl_exp)
                        except Exception: pass
                #lower_index = int(len(rpkms)*0.2)
                #upper_index = len(rpkms)-lower_index
            for reg_junction in regulated_junctions:
                ref_junction = regulated_junctions[reg_junction][1]
                altexons = regulated_junctions[reg_junction][-1]
                max_ratio = regulated_junctions[reg_junction][-2]
                rho = regulated_junctions[reg_junction][-3]
                dpsi = regulated_junctions[reg_junction][-4]
                #if 'ENSMUSG00000009350:E14.2_87617106-E15.1' in incl: print feature_exp_db[incl]
                ### Perform a sanity check for junction comparisons
                reg_pos,reg_loc = junction_locations[reg_junction]
                ref_pos,ref_loc = junction_locations[ref_junction]
                positions = reg_pos+ref_pos
                positions.sort()
                values = map(lambda x: str(x+1),feature_exp_db[reg_junction])
                try: values = psi_db[reg_junction,ref_junction]
                except Exception: values = psi_db[ref_junction,reg_junction]; reg_junction,ref_junction = ref_junction,reg_junction
                try: max_incl_psi = str(inclusion_max_psi[reg_junction])
                except Exception: max_incl_psi = '0'
                if reg_pos == positions[:2] or reg_pos == positions[-2:]:
                    #print 'Impropper junctions excluded',reg_junction,ref_junction,positions,reg_pos
                    pass
                elif ('-' not in reg_junction or '-' not in ref_junction) and platform != 'junction': ### Possible retained intron
                    if '-' not in reg_junction: query = reg_junction
                    elif '-' not in ref_junction: query = ref_junction
                    if 'I' in query:
                        intron_rel_exp = intronComparisonMethod(query)
                        if intron_rel_exp>0.15:
                            """
                            print upstream_exon, query, intron_len, upstream_exon_len,(max(feature_exp_db[query])/float(intron_len)),(max(feature_exp_db[upstream_exon])/float(upstream_exon_len))
                            print max(feature_exp_db[query])
                            print junction_locations[query]
                            print max(feature_exp_db[upstream_exon])
                            print junction_locations[upstream_exon]
                            print intron_rel_exp
                            """
                            if platform == 'junction':
                                try: reg_junction = probeset_junction_db[reg_junction] + ' ' +reg_junction
                                except Exception: pass
                                try: ref_junction = probeset_junction_db[ref_junction] + ' ' +ref_junction
                                except Exception: pass
                            export_data.write(string.join([symbol,description,reg_junction,ref_junction,altexons,str(max_ratio),str(dpsi),str(rho),max_incl_psi,reg_loc+'|'+ref_loc,'intron-retained']+values,'\t')+'\n')
                            avg = averageWithNulls(values)
                            values = map(lambda x: replaceNulls(x,avg), values)
                            clust_export_data.write(string.join([symbol+':'+reg_junction+'|'+ref_junction]+values,'\t')+'\n')
                            retained_introns.append(reg_junction)
                            retained_introns.append(ref_junction)
                else:
                    try: avg = averageWithNulls(values)
                    except Exception: continue
                    if platform == 'junction':
                        try: reg_junction = probeset_junction_db[reg_junction] + ' ' +reg_junction
                        except Exception: pass
                        try: ref_junction = probeset_junction_db[ref_junction] + ' ' +ref_junction
                        except Exception: pass
                    export_data.write(string.join([symbol,description,reg_junction,ref_junction,altexons,str(max_ratio),str(dpsi),str(rho),max_incl_psi,reg_loc+'|'+ref_loc,'junctions']+values,'\t')+'\n')
                    values = map(lambda x: replaceNulls(x,avg), values)
                    clust_export_data.write(string.join([symbol+':'+reg_junction+'|'+ref_junction]+values,'\t')+'\n')
                    exported.append(reg_junction)
                    exported.append(ref_junction)
                    
        ### Predict novel undetected events from above
        junctions=[]; other=[]; features=[]; feature_pos=[]; loc_db={}; coord_db={}; junctions_to_compare=[]; regulated_junctions={};inclusion_max_psi={}
        #max_ratio,proceed = junctionComparisonMethod(incl,excl)
        for feature in feature_exp_db:
            ### The below code is for overlapping junctions not found from the above analysis (could include exons and junctions)
            if feature not in exported and feature not in retained_introns:
                if '-' in feature and platform != 'junction':
                    junctions.append(feature)
                else:
                    other.append(feature)
                features.append(feature)
                try:
                    pos1,pos2 = junction_locations[feature][0]
                    feature_pos.append(pos1); feature_pos.append(pos2)
                    try:loc_db[pos1].append(feature)
                    except Exception: loc_db[pos1] = [feature]
                    try:loc_db[pos2].append(feature)
                    except Exception: loc_db[pos2] = [feature]
                    try: coord_db[pos1,pos2].append(feature)
                    except Exception: coord_db[pos1,pos2] = [feature]
                except Exception: pass ### occurs for junction arrays if the probeset ID is not in the database
        feature_pos = unique.unique(feature_pos); feature_pos.sort() ### These are the unique positions sorted
        overlapping_features=[]
        additional_possible_retained_introns=[] ### catch some funky intron retention events
        ### Get initial list of overlapping features
        for feature in features: ### e.g., some junction
            try:
                pos1,pos2 = junction_locations[feature][0] ### coordinates of that junction
                i1 = feature_pos.index(pos1) ### index position of the junctions
                i2 = feature_pos.index(pos2)
                #print feature, i1, i2, pos1, pos2
                if (i1-i2) != 1:
                    overlapping_features.append(feature)
                if 'I' in feature and '_' in feature:
                    possible_intron = string.split(feature,'_')[0]
                    if possible_intron not in features:
                        additional_possible_retained_introns.append((i1,i2,feature))
            except Exception:
                pass
        for feature in overlapping_features:
                #if feature not in features_examined: ### Remove this to allow for other reasonable junction or junction intron pairs that were not significant above
                ### get overlapping feature pairs
                pos1,pos2 = junction_locations[feature][0]
                i1 = feature_pos.index(pos1)
                i2 = feature_pos.index(pos2)
                ### Include a search for funky intron retention events (needed due to some weird intron retention issue)
                for (in1,in2,f2) in additional_possible_retained_introns:
                    if i1<=in1 and i2>=in2 and feature!=f2:
                        junctions_to_compare.append([feature,f2])
                for i in range(i1+1,i2):
                    overlapping = loc_db[feature_pos[i]]
                    for o in overlapping:
                        if o not in features_examined and '-' in o and '-' in feature and platform != 'junction':
                            #print feature, o
                            #print junction_locations[feature][0]
                            #print junction_locations[o][0]
                            junctions_to_compare.append([feature,o])
        #print 'junctions_to_compare:',junctions_to_compare
        ### Since this is the same coordinates, should be finding intron retention pairs
        for coord in coord_db:
            features = unique.unique(coord_db[coord])
            if len(features)>1:
                for f in features:
                    for g in features:
                        if g!=f:
                            fs = [g,f]; fs.sort()
                            if g not in exported and f not in exported:
                                if g not in retained_introns and f not in retained_introns:
                                    junctions_to_compare.append(fs)

        unique.unique(junctions_to_compare)
        for (incl,excl) in junctions_to_compare: #Not really incl, excl, just features
            max_ratio,num_incl_events,num_excl_events,dpsi,rho,max_all_psi,proceed = junctionComparisonMethod(incl,excl)
            inclusion_max_psi[incl] = max_all_psi
            #if 'ENSG00000100650' in incl: print incl,excl, max_ratio, proceed, rho, num_incl_events, num_excl_events, 'k'
            if proceed:
                altexons = ''
                if num_excl_events > num_incl_events:
                    #print max_ratio, '\t',gene
                    regulated_junctions[incl]=rho,excl,num_incl_events,num_excl_events,'incl',dpsi,rho,max_ratio,altexons
                else:
                    #print max_ratio, '\t',gene
                    regulated_junctions[excl]=rho,incl,num_excl_events,num_incl_events,'excl',dpsi,rho,max_ratio,altexons
                        
        for reg_junction in regulated_junctions:
            ref_junction = regulated_junctions[reg_junction][1]
            altexons = regulated_junctions[reg_junction][-1]
            max_ratio = regulated_junctions[reg_junction][-2]
            rho = regulated_junctions[reg_junction][-3]
            dpsi = regulated_junctions[reg_junction][-4]
            ### Perform a sanity check for junction comparisons
            reg_pos,reg_loc = junction_locations[reg_junction]
            ref_pos,ref_loc = junction_locations[ref_junction]
            positions = reg_pos+ref_pos
            positions.sort()
            values = map(lambda x: str(x+1),feature_exp_db[reg_junction])
            try: values = psi_db[reg_junction,ref_junction]
            except Exception: values = psi_db[ref_junction,reg_junction]; reg_junction,ref_junction = ref_junction,reg_junction
            try: max_incl_psi = str(inclusion_max_psi[reg_junction])
            except Exception: max_incl_psi = '0'
            if ('-' not in reg_junction or '-' not in ref_junction) and platform != 'junction': ### Possible retained intron
                if '-' not in reg_junction: query = reg_junction
                elif '-' not in ref_junction: query = ref_junction
                intron_rel_exp = intronComparisonMethod(query)
                if intron_rel_exp>0.15:
                    if platform == 'junction':
                        try: reg_junction = probeset_junction_db[reg_junction] + ' ' +reg_junction
                        except Exception: pass
                        try: ref_junction = probeset_junction_db[ref_junction] + ' ' +ref_junction
                        except Exception: pass
                    export_data.write(string.join([symbol,description,reg_junction,ref_junction,altexons,str(max_ratio),str(dpsi),str(rho),max_incl_psi,reg_loc+'|'+ref_loc,'exon-retained']+values,'\t')+'\n')
                    avg = averageWithNulls(values)
                    values = map(lambda x: replaceNulls(x,avg), values)
                    clust_export_data.write(string.join([symbol+':'+reg_junction+'|'+ref_junction]+values,'\t')+'\n')
                    retained_introns.append(reg_junction)
                    retained_introns.append(ref_junction)
                    #print query
            else:
                if platform == 'junction':
                    try: reg_junction = probeset_junction_db[reg_junction] + ' ' +reg_junction
                    except Exception: pass
                    try: ref_junction = probeset_junction_db[ref_junction] + ' ' +ref_junction
                    except Exception: pass
                export_data.write(string.join([symbol,description,reg_junction,ref_junction,altexons,str(max_ratio),str(dpsi),str(rho),'0',reg_loc+'|'+ref_loc,'others']+values,'\t')+'\n')
                avg = averageWithNulls(values)
                values = map(lambda x: replaceNulls(x,avg), values)
                clust_export_data.write(string.join([symbol+':'+reg_junction+'|'+ref_junction]+values,'\t')+'\n')
                exported.append(reg_junction)
                exported.append(ref_junction)
                #print ref_junction,reg_junction
        
    if platform == 'junction' or platform == 'AltMouse':
        probeset_gene_db={}
        from build_scripts import ExonArrayEnsemblRules
        if platform == 'junction':
            export_exon_filename = 'AltDatabase/'+species+'/'+platform+'/'+species+'_Ensembl_probesets.txt'
        if platform == 'AltMouse':
            export_exon_filename = 'AltDatabase/'+species+'/'+platform+'/'+species+'_Ensembl_probesets.txt'
        ensembl_probeset_db = ExonArrayEnsemblRules.reimportEnsemblProbesetsForSeqExtraction(export_exon_filename,'only-junctions',{})
        for gene in ensembl_probeset_db:
            for uid,pos,location in ensembl_probeset_db[gene]:
                junction_locations[uid] = pos,location
                probeset_gene_db[uid]=gene
                
    def filterByLocalJunctionExp(gene,features):
        try: symbol,description = gene_annotations[gene]
        except Exception: symbol='';description=''
        global cluster_id
        global cluster_name
        begin_time = time.time()
        count=0
        junctions_to_compare={}
        overlapping_junctions_exp={}
        ovelapping_pos={}
        existing=[]
        overlapping_junctions_test={}
        for feature in feature_exp_db:
            feature_exp = feature_exp_db[feature]
            if '-' in feature:
                pos1,pos2 = junction_locations[feature][0]
                for f2 in feature_exp_db:
                    flag=False
                    if '-' in f2:
                        if f2!=feature:
                            f2_exp = feature_exp_db[f2]
                            alt_pos1,alt_pos2 = junction_locations[f2][0]
                            positions = [pos1,pos2,alt_pos1,alt_pos2]
                            positions.sort()
                            diff = positions.index(pos2)-positions.index(pos1)
                            if diff!=1:
                                #try: junctions_to_compare[feature].append(f2)
                                #except Exception: junctions_to_compare[feature] = [f2]
                                try:
                                    overlapping_junctions_exp[feature].append([f2_exp,f2])               
                                except Exception:
                                    overlapping_junctions_exp[feature] = [[f2_exp,f2]]
                                flag=True
                            else:
                                diff = positions.index(alt_pos2)-positions.index(alt_pos1)
                                if diff!=1:
                                    try:
                                        overlapping_junctions_exp[feature].append([f2_exp,f2])
                                       
                                    except Exception:
                                        overlapping_junctions_exp[feature] = [[f2_exp,f2]]
                                    #overlapping_junctions_test[count]=[feature,]
                                    flag=True
                                    
                            if flag==True:  
                                if feature not in existing and f2 not in existing:
                                    count=count+1
                                    overlapping_junctions_test[count]=[feature,]
                                    overlapping_junctions_test[count].append(f2)
                                    existing.append(feature)
                                    existing.append(f2)
                                    
                                if feature in existing and f2 not in existing:
                                    for i in overlapping_junctions_test:
                                        if feature in overlapping_junctions_test[i]:
                                            overlapping_junctions_test[i].append(f2)
                                            existing.append(f2)
                                if f2 in existing and feature not in existing:
                                    for i in overlapping_junctions_test:
                                        if f2 in overlapping_junctions_test[i]:
                                            overlapping_junctions_test[i].append(feature)
                                            existing.append(feature)
                                
                                if feature in existing and f2 in existing:
                                    for i in overlapping_junctions_test:
                                        if feature in overlapping_junctions_test[i]:
                                            loc1=i
                                        if f2 in overlapping_junctions_test[i]:
                                            loc2=i
                                    if loc1!=loc2:
                                        for jun in overlapping_junctions_test[loc2]:
                                            if jun not in overlapping_junctions_test[loc1]:
                                                overlapping_junctions_test[loc1].append(jun)
                                        del overlapping_junctions_test[loc2]
                                          
        cluster_junction={}
        #Finding clusters and corresponding junctions
        for count in overlapping_junctions_test:
            for feature in overlapping_junctions_test[count]:
                   cluster_junction[feature]=cluster_name
                   clusterID_export_data.write(cluster_name+"\t"+feature+"\t"+junction_locations[feature][1]+"\n")
            cluster_id+=1
            cluster_name=string.join('clu_'+str(cluster_id))
            cluster_name=cluster_name.replace(" ","")

        #duration = time.time() - begin_time

        #print duration, 'seconds'
        expressed_junctions=[]
        for feature in overlapping_junctions_exp:
            counts = map(lambda x: x[0], overlapping_junctions_exp[feature])
            combined = [sum(value) for value in zip(*counts)]
            
            #if feature == 'ENSG00000002586:E1.5-E4.1':
            #print feature
            #print combined
            #print overlapping_junctions_exp[feature][0];sys.exit()
            #dpsi_values = [ratio(value) for value in zip(*[overlapping_junctions_exp[feature][0],combined])]
            #print feature
            #print overlapping_junctions[feature]
            #print overlapping_junctions_exp[feature]
            #print combined;sys.exit()
            exclusion_id = feature+'|exclusion'
            feature_exp_db[exclusion_id] = combined
            #rho : replace with clusterid
            max_ratio,num_incl_events,num_excl_events,dpsi,rho,max_all_psi,proceed = junctionComparisonMethod(feature,exclusion_id)
            if proceed:
                fe1,fe2 = string.split(feature,'-')
                if '_' in fe1 and '_' in fe2: pass
                else:
                    #"""
                    top_excl_junction=[]
                    for (exp_ls,f2) in overlapping_junctions_exp[feature]:
                        top_excl_junction.append([statistics.avg(exp_ls),f2])
                    top_excl_junction.sort()
                    #print top_excl_junction[-8:]
                    #print statistics.avg(feature_exp_db[feature])
                    top_excl_junction = top_excl_junction[-1][-1]
                    
                    t1,t2 = string.split(top_excl_junction,'-')
                    altexons = []
                    if t1!=fe1: altexons.append(fe1)
                    if t2!=fe2: altexons.append(gene+':'+fe2)
                    altexons = string.join(altexons,'|')
                    reg_pos,reg_loc = junction_locations[feature]
                    ref_pos,ref_loc = junction_locations[top_excl_junction]
                    #print [feature, dpsi,rho]
                    #top_excl_junctions = map(lambda x: x[-1], top_excl_junction[-5:])
                    #print top_excl_junctions;sys.exit()
                    #for i in top_excl_junctions: max_ratio,num_incl_events,num_excl_events,dpsi,rho,max_all_psi,proceed = junctionComparisonMethod(feature,i); print i, dpsi,rho
                    values = psi_db[feature,exclusion_id]
                    max_incl_psi = str(getMax(values))
                    # adding cluster information to PSI file
                    combined_ID = symbol+':'+feature+"|"+top_excl_junction
                    export_data.write(string.join([symbol,description,feature,top_excl_junction,altexons,str(max_ratio),str(dpsi),cluster_junction[feature],combined_ID,reg_loc+'|'+ref_loc,'junctions']+values,'\t')+'\n')
                    avg = averageWithNulls(values)
                    values_imputed = map(lambda x: replaceNulls(x,avg), values)
                    clust_export_data.write(string.join([symbol+':'+feature+'|'+top_excl_junction]+values_imputed,'\t')+'\n')
                    exported.append(feature)
                    exported.append(top_excl_junction)
        #sys.exit()
        
    gene_annotations = getGeneAnnotations(species)
    firstLine = True
    feature_exp_db={}
    gene_junction_denom={} ### Determine the max junction counts per gene per sample
    regulated_junctions = {}
    genes_examined=0; gene_increment=1000
    prior_gene = None
    gene = None
    for line in open(expFile,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            firstLine = False
            ge_header = t
            if len(t)<4:
                ### If only two samples reduce the below numbers to the minimum number of observations
                min_events=0;med_events=0
            additional_headers = string.join(['Symbol','Description','Examined-Junction','Background-Major-Junction','AltExons',"PME","dPSI",'ClusterID','UID','Coordinates','feature']+t[1:],'\t')+'\n'
            export_data.write(additional_headers)
            clust_export_data.write(line)
        else:
            uid = t[0]
            if '=' in uid:
               
                try: uid,location = string.split(uid,'=')
                except Exception: print t[0];sys.exit()
                pos1,pos2 = string.split(string.split(location,':')[-1],'-')
                pos = [int(pos1),int(pos2)]
                pos.sort()
                junction_locations[uid] = pos,location ### use the to report position and verify compared junctions
            if gene == 'ENSG00000100650': ### For testing
                proceed = True
            else: proceed = True
            if platform == 'RNASeq':
                gene = string.split(uid,':')[0]
                
            else:
                if uid in probeset_gene_db:
                    gene = probeset_gene_db[uid]
                else: proceed = False
            if proceed:
                counts = map(lambda x: float(x), t[1:])
                if platform == 'junction' or platform == 'AltMouse':
                    counts = map(lambda x: int(math.pow(2,x)), counts)  #log transform these instead, to make like junction counts
                if '-' in uid or uid in junction_locations:
                    #try: gene_junction_denom[gene].append(counts)
                    #except Exception: gene_junction_denom[gene] = [counts]
                    pass
                if genes_examined==gene_increment:
                    gene_increment+=1000
                    print '*',
                if gene != prior_gene and prior_gene !=None:
                    genes_examined+=1
                    #if len(gene_junction_denom)>0:
                    if prior_gene == '!ENSG00000198001': ### For testing
                        
                        filterByLocalJunctionExp(prior_gene,feature_exp_db)
                        #try: gene_junction_denom[prior_gene] = [max(value) for value in zip(*gene_junction_denom[prior_gene])] # sum the junction counts for all junctions across the gene
                        #except Exception: pass
                    if platform == 'RNASeq':
                        
                        filterByLocalJunctionExp(prior_gene,feature_exp_db)
                    else:
                        compareJunctionExpression(prior_gene)
                    feature_exp_db={}
                    gene_junction_denom={}
                if max(counts)>4:
                    feature_exp_db[uid] = counts
                prior_gene = gene
                    
    #compareJunctionExpression(gene)
    export_data.close()
    clust_export_data.close()
    clusterID_export_data.close()
    
    graphic_links=[]
    if (len(exported)/2)<7000:
        if (len(exported)/2)<4000:
            graphic_links = exportHeatmap(clust_export_dir,useHOPACH=False,color_gradient='yellow_black_blue',normalize=True,columnMethod='hopach',size=len(exported)/2)
    else:
        clust_export_dir,size = filterJunctionExpression(clust_export_dir)
        if size<4000:
            try: graphic_links = exportHeatmap(clust_export_dir,useHOPACH=False,color_gradient='yellow_black_blue',normalize=True,columnMethod='hopach',size=len(exported)/2,filter=True)
            except Exception: graphic_links=[]
    print len(exported)/2,'junctions exported' #,len(retained_introns)/2, 'retained introns exported...'
    return graphic_links, clust_export_dir
        
def getGeneAnnotations(species):
    gene_annotations={}
    fn = filepath('AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations_simple.txt')
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        ensembl,description,symbol = string.split(data,'\t')
        gene_annotations[ensembl] = symbol,description
    return gene_annotations

def getDownstreamExon(upstream_exon):
    gene,exon = string.split(upstream_exon,':')
    downstream_exon = gene+':E'+str(int(string.split(exon,'.')[0][1:])+1)+'.1'
    return downstream_exon

def expressionRatio(incl_exp,excl_exp,num_incl_events,num_excl_events):
    ### Selects the minor isoform and looks at its relative expression
    if num_incl_events>num_excl_events:
        max_ratio = max(excl_exp)/max(incl_exp)
    else:
        max_ratio = max(incl_exp)/max(excl_exp)
    return max_ratio

def unbiasedComparisonSpliceProfiles(root_dir,species,platform,expFile=None,min_events=-1,med_events=-1): # 4 9
    """ This is prototype code to identify critical splicing events (SI-exon-level) from single cell data prior to group assignment """
    begin_time = time.time()
    if platform == 'RNASeq': avg_all_for_SS = 'yes'
    else: avg_all_for_SS = 'no'
    agglomerate_inclusion_probesets = 'no'
    probeset_type = 'core'
    try:from build_scripts import JunctionArray
    except:
        try: import JunctionArray
        except: pass
    try: import AltAnalyze
    except: pass
    buildFromDeNovoJunctionsOnly=True
    if buildFromDeNovoJunctionsOnly and platform=='RNASeq':
        alt_junction_db={}
    else:
        exon_db, constitutive_probeset_db = AltAnalyze.importSplicingAnnotations(platform,species,probeset_type,avg_all_for_SS,root_dir)
        alt_junction_db,critical_exon_db,exon_dbase,exon_inclusion_db,exon_db = JunctionArray.getPutativeSpliceEvents(species,platform,exon_db,agglomerate_inclusion_probesets,root_dir)
    #print 'Number of Genes with Examined Splice Events:',len(alt_junction_db)
    
    if platform == 'junction':
        global probeset_junction_db; probeset_junction_db={}
        
    #alt_junction_db = {'ENSG00000100650':alt_junction_db['ENSG00000100650']}
    critical_exon_db={}
    for affygene in alt_junction_db:
        for event in alt_junction_db[affygene]:
            for critical_exon in event.CriticalExonList():
                critical_exon = affygene+':'+critical_exon
                try:
                    #print event.InclusionJunction(), event.ExclusionJunction();sys.exit()
                    inclusion_list,exclusion_list = critical_exon_db[critical_exon]
                    if '-' in event.InclusionProbeset() or (platform == 'junction' and '-' in event.InclusionJunction()):
                        inclusion_list.append(event.InclusionProbeset())
                    exclusion_list.append(event.ExclusionProbeset())
                    if platform == 'junction':
                        probeset_junction_db[event.InclusionProbeset()] = event.InclusionJunction()
                        probeset_junction_db[event.ExclusionProbeset()] = event.ExclusionJunction()
                except Exception:
                    if '-' in event.InclusionProbeset() or (platform == 'junction' and '-' in event.InclusionJunction()):
                        inclusion_list = [event.InclusionProbeset()]
                    else: inclusion_list=[]
                    exclusion_list = [event.ExclusionProbeset()]
                    #inclusion_list.append(critical_exon)
                    inclusion_list = unique.unique(inclusion_list)
                    exclusion_list = unique.unique(exclusion_list)
                    if len(inclusion_list)>0 and len(exclusion_list)>0:
                        critical_exon_db[critical_exon] = inclusion_list,exclusion_list
                    elif 'I' in critical_exon and '_' not in critical_exon and '.1' in critical_exon:
                        critical_exon_db[critical_exon] = [critical_exon],exclusion_list
                    #if affygene == 'ENSMUSG00000004952':
                        #if '.1' not in critical_exon: print critical_exon,inclusion_list,exclusion_list

    if expFile != None:
        graphic_links, cluster_input = compareRawJunctionExpression(root_dir,platform,species,critical_exon_db,expFile,min_events=min_events,med_events=med_events)
        print 'finished in',int(time.time()-begin_time), 'seconds'
        return graphic_links, cluster_input
    ### Determine the location of the gene expression file
    input_folder = root_dir+'AltResults/RawSpliceDataTemp/'+species+'/splicing-index/'
    dir_list = read_directory(input_folder) ### get all of the RawSplice files
    for filename in dir_list:
        if '.txt' in filename and ('_average' not in filename):
            dataset_name = filename
            input_dir = input_folder + dataset_name

    exportSorted(input_dir, 2)
    for filename in dir_list:
        if '.txt' in filename and ('_average' not in filename):
            dataset_name = filename
            input_dir = input_folder + dataset_name

    import RNASeq
    biological_categories = RNASeq.importBiologicalRelationships(species)
    genes = biological_categories['protein_coding']
    #genes = biological_categories['BioMarker']
    genes.update(biological_categories['transcription regulator'])
    genes.update(biological_categories['splicing regulator'])
    genes.update(biological_categories['kinase'])
    genes.update(biological_categories['GPCR'])
    
    ### Import gene expression summaries to exclude high differential genes
    fn=filepath(root_dir+'/ExpressionInput/exp.'+dataset_name[:-4]+'-steady-state.txt')
    firstLine = True
    low_diff_exp_genes={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            firstLine = False
            ge_header = t
        else:
            rpkms = map(lambda x: float(x), t[1:])
            rpkms.sort()
            lower_index = int(len(rpkms)*0.2); upper_index = len(rpkms)-lower_index
            gene = t[0]
            #if (max(rpkms)/min(rpkms))<5: ### Max allowed differential expression
            #if (rpkms[upper_index]/rpkms[lower_index])<5:
            #if gene == 'ENSMUSG00000078812': print statistics.avg(rpkms)
            if gene in genes and statistics.avg(rpkms)>5:
                low_diff_exp_genes[gene]=rpkms
    #print low_diff_exp_genes['ENSMUSG00000078812']
    
    print len(low_diff_exp_genes), 'genes with less than 5-fold differential expression'
            
    import gene_associations; from  scipy import stats
    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    
    exported_IDs=0; ids_exmained=0; genes_examined=0; gene_increment=1000
    added={}; prior_gene = False; gene_comparison_db={}
    print 'Begining to quickly find alternative exons...'
    for filename in dir_list:
        if '.txt' in filename and ('_average' not in filename) or len(dir_list)==1: ### We only want the group-comparison SI file
            export_dir = string.split(input_folder,'RawSpliceData')[0]+'Unbiased/'+filename
            export_data = export.ExportFile(export_dir)
            fn=filepath(input_folder+filename)
            x=0
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if x==0:
                    headers = t[2:]  ### first two columsn are gene and ExonID
                    #print headers
                    #print ge_header;sys.exit()
                    export_data.write(string.join(headers,'\t')+'\n') ### write header row
                    x=1
                else:
                    if platform != 'RNASeq':
                        if 'ENS' in t[0]: #ENSG human
                            gene = t[0]; altExon = t[2]
                        else:
                            altExon = t[0]
                    else:
                        gene = t[0]; altExon = t[2]
                    if genes_examined==gene_increment:
                        gene_increment+=1000
                        print '*',
                    if gene != prior_gene:
                        genes_examined+=1
                        hits={}
                        for gene in gene_comparison_db:
                            feature_db = gene_comparison_db[gene]
                            if gene == 'ENSMUSG00000078812':
                                for i in feature_db: print i
                                
                            for altexon in feature_db:
                                #print altexon
                                if altexon in critical_exon_db:
                                    inclusion_list,exclusion_list = critical_exon_db[altexon]
                                    #print inclusion_list, 'incl'
                                    altexon_SI = feature_db[altexon]
                                    for incl in inclusion_list:
                                        if incl in feature_db:
                                            incl_SI = feature_db[incl]
                                            with warnings.catch_warnings():
                                                warnings.filterwarnings("ignore") ### hides import warnings
                                                try: rho,p = stats.pearsonr(altexon_SI,incl_SI)
                                                except Exception: rho = 1
                                            #print rho, altexon,incl
                                            #print string.join(map(str,altexon_SI),'\t')
                                            #print string.join(map(str,incl_SI),'\t');sys.exit()
                                            if rho>0.4:
                                                hits[altexon]=[]
                                            if gene == 'ENSMUSG00000078812':print '***', incl
                                    #print inclusion_list, 'excl'
                                    for excl in inclusion_list:
                                        if excl in feature_db:
                                            excl_SI = feature_db[excl]
                                            with warnings.catch_warnings():
                                                warnings.filterwarnings("ignore") ### hides import warnings
                                                rho,p = stats.pearsonr(altexon_SI,excl_SI)
                                            if rho<-0.4:
                                                hits[altexon]=[]
                                            if gene == 'ENSMUSG00000078812': print '***', excl
                            if gene == 'ENSMUSG00000078812': print hits
                        for altExon in hits:
                            added[altExon]=[]
                            log_folds = feature_db[altExon]
                            log_folds = map(str, log_folds)
                            values = string.join([altExon]+log_folds,'\t')+'\n' ### [t[3]]+ before log_folds?
                            if gene in gene_to_symbol: symbol = gene_to_symbol[gene][0]+"  "
                            else: symbol = ''
                            export_data.write(symbol+values)
                            exported_IDs+=1
                        gene_comparison_db={}
                        #if exported_IDs> 1: sys.exit()
                    prior_gene = gene
                    
                    if ';' in altExon:
                        altExon1, altExon2 = string.split(altExon,';')
                        altExons = [altExon1,gene+':'+altExon2]
                    else:
                        altExons = [altExon]
                    for altExon in altExons:
                        if altExon not in added and gene in low_diff_exp_genes: #altExon in alternative_exon_db and
                            #if altExon == 'ENSMUSG00000022841:E7.2':
                            ids_exmained+=1
                            #values = map(lambda x: float(x)*-1, t[3:]) #reverse the fold for cluster visualization
                            values = map(lambda x: float(x), t[3:])
                            #print len(headers),len(values);sys.exit()
                            avg = statistics.avg(values)
                            log_folds = map(lambda x: x-avg, values) ### calculate log folds and make these strings
                            i=0; si_exp_list = [] ### store the pairs of SI and gene expression for each sample
                            rpkms = list(low_diff_exp_genes[gene])
                            for si in log_folds: si_exp_list.append([rpkms[i],si]); i+=1
                            rpkms.sort()
                            si_exp_list.sort() ### This object contains both gene expression and SI values
                            max_rpkm = rpkms[-1]
                            half_max_rpkm = max_rpkm/2 ### Only look at genes in which there is less than a 2 fold differnce
                            s = bisect.bisect_right(rpkms,half_max_rpkm)
                            si_highGeneExp = map(lambda (rpkm,si): si, si_exp_list[s:])
                            #print si_exp_list[s:]
                            #cv = statistics.stdev(si_highGeneExp)/statistics.avg(si_highGeneExp)
                            si_highGeneExp.sort()
                            try:
                                biggest_diff = si_highGeneExp[-2]-si_highGeneExp[1]
                                #print biggest_diff
                                #print cv
                                if gene == 'ENSG00000009413':
                                    print altExon, biggest_diff
                                    print si_highGeneExp
                                if biggest_diff>2 and len(si_highGeneExp)>20:
                                    try:
                                        feature_db = gene_comparison_db[gene]
                                        feature_db[altExon] = log_folds
                                    except Exception:
                                        feature_db={}
                                        feature_db[altExon] = log_folds
                                        gene_comparison_db[gene] = feature_db
                                    
                                    #added[altExon]=[]
                                    #log_folds = map(str, log_folds)
                                    #values = string.join([altExon]+log_folds,'\t')+'\n' ### [t[3]]+ before log_folds?
                                    #if gene in gene_to_symbol: symbol = gene_to_symbol[gene][0]+"  "
                                    #else: symbol = ''
                                    #export_data.write(symbol+values)
                                    #exported_IDs+=1
                            except Exception: pass ### Occurs with less than 4 samples in the si_highGeneExp set
    
    print exported_IDs, 'exported ID values for clustering out of',ids_exmained
    export_data.close()
    return export_dir, exported_IDs

def AllGroupsNIComparison(root_dir, species, array_type):
    if array_type == 'RNASeq': avg_all_for_SS = 'yes'
    else: avg_all_for_SS = 'no'
    agglomerate_inclusion_probesets = 'no'
    #calculateNormalizedIntensities(root_dir, species, array_type, avg_all_for_SS = avg_all_for_SS, probeset_type = 'core')
    ### This analysis is designed for datasets without known variables (e.g., single cell seq)
    
    from build_scripts import JunctionArray
    exon_db, constitutive_probeset_db = AltAnalyze.importSplicingAnnotations(array_type,species,probeset_type,avg_all_for_SS,root_dir)
    alt_junction_db,critical_exon_db,exon_dbase,exon_inclusion_db,exon_db = JunctionArray.getPutativeSpliceEvents(species,array_type,exon_db,agglomerate_inclusion_probesets,root_dir)
    print 'Number of Genes with Examined Splice Events:',len(alt_junction_db)
 
    for affygene in alt_junction_db:
        for event in alt_junction_db[affygene]:
            event.InclusionProbeset()
            event.ExclusionProbeset()

def createExpressionSQLdb(species,platform,expFile):
    """ Store junction/exon RPKMs or probesets expression in a SQL database"""
    start=time.time()
    from import_scripts import SQLInterface
    DBname = 'FeatureExpression'
    schema_text ='''-- Schema for species specific AltAnalyze junction/exon expression data.

-- Genes store general information on each Ensembl gene ID
create table ExonExp (
    uid          text primary key,
    gene         text,
    expression     text
);
'''
    conn = SQLInterface.populateSQLite(species,platform,DBname,schema_text=schema_text) ### conn is the database connnection interface
    
    ### Populate the database
    print 'importing', expFile
    fn=filepath(expFile)
    for line in open(fn,'r').xreadlines():
        data = line.strip()
        t = string.split(data,'\t')
        uid = t[0]; expression = string.join(t[1:],'\t')
        try: gene = string.split(uid,':')[0]
        except Exception: print 'not RNASeq - function not supported';kill
        #print exonID,gene,sequence
        ### Store this data in the SQL database
        command = """insert into ExonExp (uid, gene, expression)
        values ('%s', '%s','%s')""" % (uid,gene,expression)
        conn.execute(command)
        
    conn.commit() ### Needed to commit changes
    conn.close()
    time_diff = str(round(time.time()-start,1))
    print 'Exon/Junction Expression Data added to SQLite database in %s seconds' % time_diff

def logratio(list):
    return list[0] - list[1]

def matchAndCorrelate(prime, secondary, output_dir, rho_cutoff):
    ### Take two files and correlate their IDs to any matching
    
    export_object = export.ExportFile(output_dir[:-4]+'-'+str(rho_cutoff)+'.txt')
    export_object.write('Feature1\tFeature2\trho\n')
    
    firstLine = True; prime_db={}
    for line in open(prime,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            firstLine = False
        else:
            prime_db[t[0]] = map(float,t[1:])
            
    firstLine = True; secondary_db={}; key_db={}
    for line in open(secondary,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            firstLine = False
        else:
            try: gene_id,probe_id = string.split(t[0],':')
            except Exception: gene_id = t[0]; probe_id = t[0]
            try: secondary_db[gene_id].append((map(float,t[1:]),probe_id))
            except Exception: secondary_db[gene_id] = [(map(float,t[1:]),probe_id)]
        
    from  scipy import stats
    top_correlated={}
    for gene in prime_db:
        prime_profile = prime_db[gene]
        if gene in secondary_db:
            for (secondary_db_profile, probe_id) in secondary_db[gene]:
                rho,p = stats.pearsonr(prime_profile,secondary_db_profile)
                if rho > rho_cutoff or rho < -1*rho_cutoff:
                    #print gene, '\t',probe_id, '\t',rho
                    export_object.write(gene+'\t'+probe_id+'\t'+str(rho)+'\n')
    export_object.close()
    
def getMax(values):
    values2=[]
    for i in values:
        try: values2.append(float(i))
        except Exception: pass
    return max(values2)

def replaceNulls(x,avg):
        if x=='':
            return '0'
        else:
            return str(float(x)-avg)
        
def nullReplace(dpsi_values,combined,min_exp): ### Don't count un-detected genes in later stats
        null_replaced=[]
        i=0
        for v in combined:
            if v<(min_exp): null_replaced.append('') #changed July3 2017(min_exp+1)-Meenakshi
            else: null_replaced.append(str(dpsi_values[i]))
            i+=1
        return null_replaced
    
def averageWithNulls(values):
        avg_vals=[]
        for i in values:
            try: avg_vals.append(float(i))
            except Exception: pass
        avg = statistics.avg(avg_vals)
        return avg
    
def expressionSortImport(filename,filter_db=None):
    firstLine = True; exp_db={}; lines=0; max_var = 3
    for line in open(filename,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            headers = t[1:]
            header_ind = map(lambda x: (x,headers.index(x)),headers) ### store header index values
            header_ind.sort()
            #print header_ind
            headers_ind = map(lambda (x,i): i,header_ind)
            firstLine = False
        else:
            try: exp_data = map(float,t[1:])
            except Exception:
                exp_data=[]
                for value in t[1:]:
                    try: value = float(value)
                    except Exception: pass
                    exp_data.append(value)
                
            exp_data = map(lambda i: exp_data[i],headers_ind)
            if filter_db != None:
                key = t[0]
                if ':' in key:
                    key = string.split(key,':')[0]
                    max_var = 0
                if key in filter_db:
                    if max(exp_data)>max_var:
                        exp_db[key] = exp_data
            else:
                exp_db[t[0]] = exp_data
        lines+=1
    print len(exp_db),'IDs imported from', export.findFilename(filename)
    return exp_db
            
def featureCorrelate(species,query_dir,feature_dir,output_dir,feature_type):
    #python ExpressionBuilder.py --species Hs --i /Volumes/SEQ-DATA/SingleCell-Churko/Filtered/Unsupervised-AllExons/AltResults/Unbiased/junctions-2-5/top_alt_junctions_clust-TTN_all_selected.txt --additional /Volumes/SEQ-DATA/SingleCell-Churko/Filtered/Unsupervised-AllExons/ExpressionInput/exp.CM-TTN-steady-state.txt --analysis featureCorrelate --var "splicing regulator"
    ### Correlate features in a file to feature-specific gene expression data (e.g., "splicing regulator")
    try: export_object = export.ExportFile(output_dir[:-4]+'-'+feature_type+'.txt')
    except Exception: export_object = export.ExportFile(output_dir[:-4]+'-None.txt')
    export_object.write('UID\tFeature\trho\n')
    
    import RNASeq; import ExpressionBuilder
    biological_categories = RNASeq.importBiologicalRelationships(species)
    gene_to_symbol_db = ExpressionBuilder.importGeneAnnotations(species)
    if feature_type != None:
        filter_genes = biological_categories[feature_type]
    else:
        filter_genes=None
    try: print len(filter_genes),feature_type,'genes imported for comparison...'
    except Exception: pass
    query_db = expressionSortImport(query_dir)
    feature_db = expressionSortImport(feature_dir,filter_genes)
        
    from  scipy import stats
    top_correlated={}
    for uid in query_db:
        query_exp_profile = query_db[uid]
        for gene in feature_db:
            feature_exp_profile = feature_db[gene]
            try: rho,p = stats.pearsonr(query_exp_profile,feature_exp_profile)
            except Exception:
                ### If missing values are present, only correlate to where the values are present
                query_exp_profile2=[]
                feature_exp_profile2=[]
                i=0
                for v in query_exp_profile:
                    if v!='':
                        query_exp_profile2.append(query_exp_profile[i])
                        feature_exp_profile2.append(feature_exp_profile[i])
                    i+=1
                if len(feature_exp_profile2)>20:
                    rho,p = stats.pearsonr(query_exp_profile2,feature_exp_profile2)
                else:
                    rho = 0
            try: symbol = gene_to_symbol_db[gene]
            except Exception: symbol = gene
            try: top_correlated[uid].append([abs(rho),symbol[0],rho])
            except Exception: top_correlated[uid]=[[abs(rho),symbol[0],rho]]

    for uid in top_correlated:
        res = top_correlated[uid]
        res.sort()
        feature = res[-1][1]
        rho = res[-1][-1]
        export_object.write(uid+'\t'+feature+'\t'+str(rho)+'\n')
    export_object.close()
    
def lncRNANeighborCorrelationAnalysis(dataset_dir):
    ### dataset_dir is the ExpressionOuput DATASET file location
    
    #Get all analyzed genes and coordinates
    print 'Importing the DATASET file'
    global gene_symbol_db
    fn=filepath(dataset_dir); gene_coordinate_db={}; all_lncRNA_db={}; coord_gene_db={}; gene_symbol_db={}
    chr_coord_list=[]; positive_strand=[]; negative_strand=[]
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        geneID = t[0]
        symbol = t[2]
        protein_class = t[8]
        chr = t[9]
        strand = t[10]
        coordinates = tuple(string.split(t[11],'-'))
        coord_list = chr,coordinates,strand
        gene_coordinate_db[geneID]=coord_list
        coord_gene_db[coord_list] = geneID
        gene_symbol_db[geneID] = symbol
        chr_coord_list.append(coord_list)
        if '+' in strand:
            positive_strand.append(coord_list)
        else:
            negative_strand.append(coord_list)
        if 'lincRNA' in protein_class or 'lncRNA' in protein_class:
            all_lncRNA_db[geneID]=[] 
    chr_coord_list.sort(); positive_strand.sort(); negative_strand.sort()

    useClusterFile = False
    #Get all significantly differentially expressed genes
    if useClusterFile:
        cluster_file = string.replace(dataset_dir,'ExpressionOutput','ExpressionOutput/Clustering/')
        cluster_file = string.replace(cluster_file,'DATASET-','SampleLogFolds-')
    else:
        cluster_file = string.replace(dataset_dir,'ExpressionOutput','ExpressionInput')
        cluster_file = string.replace(cluster_file,'DATASET-','exp.')
        cluster_file = string.replace(cluster_file,'.txt','-steady-state.txt')
    
    print 'Importing the cluster file'
    fn=filepath(cluster_file); differentially_exp_db={}; lncRNA_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        uid = string.split(data,'\t')[0]
        if ' ' in uid:
            uid = string.split(uid,' ')[0]
        differentially_exp_db[uid]=[]
        if uid in all_lncRNA_db:
            lncRNA_db[uid]=[]

    #import random
    #lncRNA_db = random.sample(differentially_exp_db,len(lncRNA_db))
    
    print 'Number of lncRNAs regulated in clusters:',len(lncRNA_db)
    #Get the MarkerFinder cluster assignments of all analyzed genes
    root_dir = string.split(dataset_dir,'ExpressionOutput')[0]
    markerfinder = root_dir+'ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt'
    
    print 'Importing the MarkerFinder file'
    fn=filepath(markerfinder); cluster_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        geneID = t[0]
        cluster = t[-1]
        if geneID in differentially_exp_db:
            cluster_db[geneID] = cluster

    cluster_regulated_lncRNAs={}
    for geneID in lncRNA_db:
        try:
            cluster_regulated_lncRNAs[cluster_db[geneID]]+=1
        except Exception:
            try: cluster_regulated_lncRNAs[cluster_db[geneID]]=1
            except Exception: pass
    for cluster in cluster_regulated_lncRNAs:
        print cluster, cluster_regulated_lncRNAs[cluster]
    
    print 'Searching for lncRNA positional correlations'
    direction_list=['both','forward','reverse']
    print '\tExamining both strands'
    for direction in direction_list:
        nmc = searchUpOrDownstreamGenes(lncRNA_db,gene_coordinate_db,cluster_db,coord_gene_db,chr_coord_list,direction)
        #print len(nmc),len(chr_coord_list)
    print '\tExamining the positive strand'
    for direction in direction_list:
        nmc = searchUpOrDownstreamGenes(lncRNA_db,gene_coordinate_db,cluster_db,coord_gene_db,positive_strand,direction)
        #print len(nmc),len(positive_strand)
    print '\tExamining the negative strand'
    for direction in direction_list:
        nmc = searchUpOrDownstreamGenes(lncRNA_db,gene_coordinate_db,cluster_db,coord_gene_db,negative_strand,direction)
        #print len(nmc),len(negative_strand)
        
def searchUpOrDownstreamGenes(lncRNA_db,gene_coordinate_db,cluster_db,coord_gene_db,coord_list,direction):
    neighbor_matching_cluster_db={}; multiple_lncRNAs=0; number_of_neighbors=0
    for geneID in lncRNA_db:
        coordinates = gene_coordinate_db[geneID]
        if coordinates in coord_list: ### strand dependent
            rank_index = coord_list.index(coordinates)
            if geneID in cluster_db:
                cluster = cluster_db[geneID]
                if direction == 'forward':
                    search_pos = [4, 3,2,1]
                    search_pos = [1]
                elif direction == 'reverse':
                    search_pos = [-4, -3,-2,-1]
                    search_pos = [-1]
                else:
                    search_pos = [4,3,2,1,-3,-2,-1, -4]
                    search_pos = [1,-1]
                for oi in search_pos:
                    i = coord_list[rank_index-oi]
                    neighbor_gene = coord_gene_db[i]
                    symbol = gene_symbol_db[neighbor_gene]
                    if neighbor_gene in cluster_db and neighbor_gene not in lncRNA_db and neighbor_gene != geneID and '.' not in symbol:
                        ncluster = cluster_db[neighbor_gene]
                        if cluster == ncluster:
                            if neighbor_gene in lncRNA_db:
                                multiple_lncRNAs+=1
                            try: neighbor_matching_cluster_db[geneID]+=1; number_of_neighbors+=1
                            except Exception: neighbor_matching_cluster_db[geneID]=1; number_of_neighbors+=1
                            print cluster,gene_symbol_db[geneID],gene_symbol_db[neighbor_gene]
        #print 'multiple_lncRNAs:', multiple_lncRNAs, number_of_neighbors
    return neighbor_matching_cluster_db

def getHighestExpressingGenes(input_file,output_dir,topReported):
    ### Sorts genes based on RPKM (ignore read counts)
    bisectValues = False
    if topReported<100:
        bisectValues = True
    firstLine = True
    sampleExpression_db={}
    for line in open(input_file,'rU').xreadlines():
        data = cleanUpLine(line)
        values = string.split(data,'\t')
        if firstLine:
            headers = values[1:]
            for i in headers:
                sampleExpression_db[i]=[]
            firstLine = False
            print len(values)
        else:
            gene = values[0]
            i=0
            for rpkm in values[1:]:
                sampleExpression_db[headers[i]].append((float(rpkm),gene))

                i+=1
    for sample in sampleExpression_db:
        Sample = string.replace(sample,'.bed','')
        Sample = string.replace(Sample,'.cel','')
        Sample = string.replace(Sample,'.CEL','')
        Sample = string.replace(Sample,':','-')
        export_object = export.ExportFile(output_dir+'/'+Sample+'-top_'+str(topReported)+'.txt')
        export_object.write('Genes\tSystemCode\tChanged\n')
        sampleExpression_db[sample].sort()
        if bisectValues:
            rpkms =  map(lambda x: x[0], sampleExpression_db[sample])
            print rpkms[-5:]
            s = bisect.bisect_right(rpkms,float(topReported))
            topExpGenes = map(lambda x: str(x[1]), sampleExpression_db[sample][-1*(len(rpkms)-s):])
            print Sample,len(topExpGenes), s
        else:
            topExpGenes = map(lambda x: str(x[1]), sampleExpression_db[sample][-1*topReported:])
        for gene in topExpGenes:
            if 'ENS' in gene or 'ENF' in gene: system = 'En'
            else: system = 'Sy'
            export_object.write(gene+'\t'+system+'\t1\n')
        export_object.close()
    print 'The top',topReported,'expressing genes have been exported to',output_file
    
def returnRowHeaderForMaxEntry(filename,top):
    ### Used for enrichment analysis matrices to find the most significant term for each comparison/group/sample
    output_file = filename[:-4]+'_top%d.txt' % top
    export_object = export.ExportFile(output_file)
    from visualization_scripts import clustering; import numpy
    matrix, column_header, row_header, dataset_name, group_db = clustering.importData(filename,reverseOrder=False)
    matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
    column_header, row_header = row_header, column_header
    
    x=0
    for row in matrix:
        comparison = row_header[x]
        copied_row_values = list(row)
        copied_row_values.sort()
        max_vals = copied_row_values[-1*top:]
        max_vals.reverse()
        term = column_header[list(row).index(max_vals[0])]
        term+= '('+str(max_vals[0])[:4]+')|'
        
        if top>1:
            term+= column_header[list(row).index(max_vals[1])]
            term+= '('+str(max_vals[1])[:4]+')|'
        if top>2:
            term+= column_header[list(row).index(max_vals[2])]
            term+= '('+str(max_vals[2])[:4]+')|'
        if top>3:
            term+= column_header[list(row).index(max_vals[3])]
            term+= '('+str(max_vals[3])[:4]+')|'
        if top>4:
            term+= column_header[list(row).index(max_vals[4])]
            term+= '('+str(max_vals[4])[:4]+')|'
        if top>5:
            term+= column_header[list(row).index(max_vals[5])]
            term+= '('+str(max_vals[5])[:4]+')|'
        if top>6:
            term+= column_header[list(row).index(max_vals[6])]
            term+= '('+str(max_vals[6])[:4]+')|'
        if top>7:
            term+= column_header[list(row).index(max_vals[7])]
            term+= '('+str(max_vals[7])[:4]+')|'
        if top>8:
            term+= column_header[list(row).index(max_vals[8])]
            term+= '('+str(max_vals[8])[:4]+')|'
        if top>9:
            term+= column_header[list(row).index(max_vals[9])]
            term+= '('+str(max_vals[9])[:4]+')|'
        if top>10:
            term+= column_header[list(row).index(max_vals[10])]
            term+= '('+str(max_vals[10])[:4]+')|'
    
        #print comparison, term
        export_object.write(comparison+'\t'+term+'\n')
        x+=1
    export_object.close()

def orderHeatmapByMarkerFinderOrder(clustered_file):
    output_file = clustered_file[:-4]+'_MarkerFinderOrdered.txt'
    export_object = export.ExportFile(output_file)
    
    firstLine = True
    geneOrder=[]
    arrayOrder={}
    for line in open(input_file,'rU').xreadlines():
        data = line[:-1]
        values = string.split(data,'\t')
        if firstLine:
            headers = values[1:]
            for i in headers:
                group,sample = string.split(i,':')
                try: arrayOrder[group].append(i)
                except Exception: arrayOrder[group] = [i]
            firstLine = False
        else:
            gene = values[0]
            
def exportSorted(filename, sort_col, excludeHeader=True):
    ### efficient method to sort a big file without storing everything in memory
    ### http://stackoverflow.com/questions/7079473/sorting-large-text-data
    ouput_file = filename[:-4]+'-sorted' ### temporary
    index = []
    f = open(filename)
    firstLine = True
    while True:
        offset = f.tell()
        line = f.readline()
        if not line: break
        length = len(line)
        col = line.split('\t')[sort_col].strip()
        if firstLine:
            header = line
            firstLine = False
            if excludeHeader == False:
                index.append((col, offset, length))
        else:
            index.append((col, offset, length))
    f.close()
    index.sort()
    
    o = open(ouput_file,'w')
    f = open(filename)
    if excludeHeader:
        o.write(header)
    for col, offset, length in index:
        #print col, offset, length
        f.seek(offset)
        o.write(f.read(length))
    o.close()
    try:
        ### Error occurs when the file can't be deleted due to system permissions
        os.remove(filename)
        os.rename(ouput_file,filename)
        return filename
    except Exception:
        return ouput_file
    
def importJunctionPositions(species,array_type):
    ### Look up the junction coordinates for the region
    if array_type == 'RNASeq':
        probesets = 'junctions'
    else:
        probesets = 'probesets'
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_'+probesets+'.txt'
    fn=filepath(filename)
    region_db = {}
    firstRow=True
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow: firstRow = False
        else:
            probeset = t[0]
            gene = t[2]
            chr = t[4]

            if '|' in t[13]:
                region1 = string.split(t[13],'|')[1]
                region2 = string.split(t[14],'|')[0]
                junction_coord = region1+'-'+region2
                region_db[probeset] = junction_coord
                region_db[junction_coord] = gene+':'+t[12],chr ### Junction region for this version of Ensembl
    return region_db

def convertArrayReciprocalJunctionToCoordinates(species,array_type,dir_path,start_version,end_version):
    """ Script for taking junction array defined ASPIRE or LinearRegression junction pairs, extracting the region coordinates
    and exporting those coordinates with end_version EnsMart block and region IDs"""
    
    UI.exportDBversion(start_version) ### Database EnsMart version
    
    region_db = importJunctionPositions(species,array_type)
    
    comparison_db={}
    dir_list = UI.read_directory(dir_path)
    for filename in dir_list:
        if '.txt' in filename:
            comparsion = string.split(filename,'.')[0]
            proceed = False
            if ('ASPIRE' in filename or 'egress' in filename) and ('GENE' not in filename and 'inclusion' in filename): proceed = True ### Need to create a special analysis just for reciprocal junctions
            if proceed: ### Don't include splicing-index results for RNA-Seq
                comparison_db[comparsion] = {}
                fn=filepath(dir_path+'/'+filename)
                x=0
                for line in open(fn,'rU').xreadlines():
                    data = cleanUpLine(line)
                    t = string.split(data,'\t')
                    if x==0:
                        p1 = t.index('probeset1')
                        p2 = t.index('probeset2')
                        reg_call = t.index('regulation_call')
                        e1 = t.index('exons1')
                        x=1
                    else:
                        if '-' in t[e1]:
                            jc1 = region_db[t[p1]]
                            jc2 = region_db[t[p2]]
                            chr = region_db[jc1][1]
                            db = comparison_db[comparsion]
                            db[jc1,jc2] = t[reg_call],chr

    UI.exportDBversion(end_version) ### Database EnsMart version
    converted_comparison_db = {}
    region_db2 = importJunctionPositions(species,'RNASeq')
    
    eo = export.ExportFile(dir_path+'/converted_junction_events.txt')
    succeed=0; fail=0
    for comparison in comparison_db:
        for (j1,j2) in comparison_db[comparison]:
            reg_call,chr = comparison_db[comparison][(j1,j2)]
            if j1 in region_db2 and j2 in region_db2:
                junction1_id,chr = region_db2[j1]
                junction2_id,chr = region_db2[j2]
                #print junction1_id, junction2_id, j1,j2, comparison;sys.exit()
            else:
                junction1_id=''
                junction2_id=''
            j1=chr+':'+j1
            j2=chr+':'+j2
            values = string.join([comparison,j1,j2,junction1_id,junction2_id,reg_call],'\t')+'\n'
            eo.write(values)
    eo.close()
        
def convertPSIJunctionIDsToPositions(psi_file,regulated_file):
    """ Links up PSI genomic positions with IDs in a significantly differentially regulated PSI results file """
    
    fn=filepath(psi_file)
    x=0
    coord_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            symbol = t.index('Symbol')
            minor = t.index('Examined-Junction')
            major = t.index('Background-Major-Junction')
            coord = t.index('Coordinates')
            x=1
        else:
            uid = t[symbol]+':'+t[minor]+'|'+t[major]
            coordinates = t[coord]
            coord_db[uid] = coordinates
    
    dir_path = export.findParentDir(regulated_file)
    comparison = export.findFilename(regulated_file)
    eo = export.ExportFile(dir_path+'/coordinate_PSI_events.txt')
    fn=filepath(regulated_file) 
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        event = t[0]
        event = string.replace(event,'@',':')
        event = string.replace(event,'&',':')
        event = string.replace(event,'__','|')
        regulation = t[1]
        if event in coord_db:
            coordinates = coord_db[event]
            values = string.join([comparison,event,coordinates,regulation],'\t')+'\n'
            eo.write(values)
    eo.close()

if __name__ == '__main__':
    #predictSplicingEventTypes('ENSG00000123352:E15.4-E16.1','ENSG00000123352:E15.3-E16.1');sys.exit()
    test=False
    if test:
        directory = '/Volumes/SEQ-DATA/AML_junction/AltResults/AlternativeOutput/'
        dir_list = read_directory(directory)
        for file in dir_list:
          if 'PSI-clust' in file: 
            filename = meanCenterPSI(directory+'/'+file)
            #filterJunctionExpression(filename,minPercentPresent=0.75)
        #exportHeatmap('/Volumes/My Passport/AML-LAML/LAML1/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI-clust-filt.txt',color_gradient='yellow_black_blue',columnMethod='hopach')
        #sys.exit()
    #convertPSIJunctionIDsToPositions('/Volumes/SEQ-DATA/Grimeslab/TopHat/AltResults/AlternativeOutput/Mm_RNASeq_top_alt_junctions-PSI.txt','/Users/saljh8/Documents/1-dataAnalysis/SplicingFactors/Grimes-MarkerFinder-v2.txt.txt')
    #convertArrayReciprocalJunctionToCoordinates('Hs','junction','/Volumes/Time Machine Backups/dataAnalysis/SplicingFactor/Hs/hglue/Marto/AltResults/AlternativeOutput','EnsMart65','EnsMart72')
    #sys.exit()
    
    fold = 2
    pval = 0.05
    ptype = 'rawp'
    species = 'Hs'
    analysis = 'goelite'
    array_type = "3'array"
    norm = 'RPKM'
    graphic_links=[]
    additional = None
    use_downregulated_labels=True
    excludeGenes = None
    expression_data_format = 'non-log'
    expression_data_format = 'log'
    var = None
    
    ################  Comand-line arguments ################
    import getopt
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a tab-delimited input expression file in the command-line"
        print "Example: python ExpressionBuilder.py --i '/Users/me/GSEXXX/ExpressionOutput' --p 0.05 --f 1.5 --ptype rawp --analysis summary --direction up --platform RNASeq"
        sys.exit()
        
        #Building GO-Elite inputs and running GO-Elite in batch
        #python ExpressionBuilder.py --i /Users/saljh8/Desktop/C4-hESC/ExpressionOutput  --p 0.05 --f 2 --ptype adjp --analysis goelite --direction down --platform gene --species Hs --additional goelite

        #Generating signatures
        #python ExpressionBuilder.py --i /Users/saljh8/Desktop/C4-hESC/GO-Elite/upregulated/ --analysis signature --inputSource Ensembl --outputSource EntrezGene
        
        #Filtering expression datasets
        #python ExpressionBuilder.py --i /Users/saljh8/Desktop/C4-hESC/ExpressionOutput --analysis filter
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','o=','f=','p=','a=','platform=',
                                                    'ptype=','analysis=','species=','direction=',
                                                    'inputSource=','outputSource=', 'additional=',
                                                    'excludeGenes=','var='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i': directory=arg
            elif opt == '--o': output_file=arg
            elif opt == '--f': fold=float(arg)
            elif opt == '--p': pval=float(arg)
            elif opt == '--ptype': ptype=arg
            elif opt == '--analysis' or opt == '--a': analysis=arg
            elif opt == '--species': species=arg
            elif opt == '--platform': array_type=arg
            elif opt == '--inputSource': input_source=arg
            elif opt == '--outputSource': output_source=arg
            elif opt == '--additional': additional=arg
            elif opt == '--excludeGenes': excludeGenes=arg ### File location for text file with genes to exclude
            elif opt == '--var': var=arg
            elif opt == '--direction':
                if 'own' in arg:
                    use_downregulated_labels = True
                else:
                    use_downregulated_labels = False
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

    ### Allow for import of genes to exclude (e.g., sex-associated or pseudogenes)
    try: genesToExclude = excludeGenesImport(excludeGenes)
    except Exception: genesToExclude = {}
    print analysis
    if array_type == 'RNASeq':
        gene_exp_threshold = 50
        gene_rpkm_threshold = 3
    if analysis == 'matchAndCorrelate':
        matchAndCorrelate(directory, var, output_source, additional)
    if analysis == 'returnRowHeaderForMaxEntry':
        ### Used primarily for combining LineageProfiler z-scores to report the top categories across compendiums
        try: returnRowHeaderForMaxEntry(directory,int(var))
        except Exception: pass
    if analysis == 'featureCorrelate':
        try: output_file = output_file
        except Exception: output_file=directory
        featureCorrelate(species,directory,additional,output_file,var)
    if analysis == 'MarkerFinderOrder':
        ### Used for combining the all gene MarkerFinder ordered results with already clustered results (e.g., significant gene)
        ### to return significantly differentially expressed genes (expressed sufficiently) and cluster samples within classes
        ### but order by MarkerFinder correlations and groups
        orderHeatmapByMarkerFinderOrder(directory)
    if analysis == 'unbiased':
        #python ExpressionBuilder.py --species Hs --platform RNASeq --i "/Volumes/My Passport/salomonis2/SRP042161_GBM-single-cell/bams/" --a unbiased --additional "/Volumes/My Passport/salomonis2/SRP042161_GBM-single-cell/bams/ExpressionInput/counts.GBM_scRNA-Seq.txt"
        import RNASeq
        #export_dir = '/Volumes/SEQ-DATA/Grimes/14018_gmp-pro/Lattice/Full/AltResults/Unbiased/DataPlots/Clustering-myeloblast-hierarchical_euclidean_euclidean.txt'
        #export_dir = '/Volumes/SEQ-DATA/SingleCell-Churko/AltResults/Unbiased/DataPlots/Clustering-CM-hierarchical_euclidean_euclidean.txt'
        #calculateNormalizedIntensities(directory, species, array_type, analysis_type = 'raw', expFile = additional)
        var = unbiasedComparisonSpliceProfiles(directory,species,array_type,expFile=additional,min_events=1,med_events=1)
        #export_dir, exported_IDs = var
        #print export_dir
        #RNASeq.correlateClusteredGenes(export_dir)
    if analysis == 'highest-expressing':
        getHighestExpressingGenes(directory,output_file,float(var))
    if analysis == 'lncRNA':
        lncRNANeighborCorrelationAnalysis(directory)
    if analysis == 'NI':
        calculateNormalizedIntensities(directory,species,array_type)
    if analysis == 'AltExonConfirmed':
        ### Grab the alternative exons in the AltExonConfirmed GO-Elite folder, combine them and filter the splicing-index raw table
        input_dir = directory+'/AltExonConfirmed/'
        cluster_file, rows_in_file = buildAltExonClusterInputs(input_dir,species,array_type,dataType='AltExonConfirmed')
        if rows_in_file < 7000:
            exportHeatmap(cluster_file,size=rows_in_file)
    if analysis == 'goelite' or analysis == 'summary':
        #python ExpressionBuilder.py --f 2 --p 0.05 --ptype adjp --analysis summary --i /inputs
        buildCriterion(fold, pval, ptype, directory+'/',analysis,UseDownRegulatedLabel=use_downregulated_labels,genesToExclude=genesToExclude)
        if additional == 'goelite':
            import multiprocessing as mlp
            runGOElite(species,directory)
    if analysis == 'filter':
        filterDatasetFile(directory+'/')
    if analysis == 'signature':
        import gene_associations
        directory+='/'; gene_conversion_db={}
        dir_list = read_directory(directory)
        for file in dir_list:
            filename = directory+'/'+file
            db,input_data_db = gene_associations.IDconverter(filename, species, input_source, output_source,analysis=analysis)
            gene_conversion_db[file] = db,input_data_db
        exportSignatures(gene_conversion_db,directory,species)
    if analysis == 'QC':
        graphic_links = visualizeQCPlots(directory)
    elif analysis == 'LineageProfiler':
        graphic_links = performLineageProfiler(directory,graphic_links)

