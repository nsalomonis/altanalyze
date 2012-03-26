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
import unique
import statistics
import math
import reorder_arrays
import ExonArray
import export
import time
import UI
import BuildAffymetrixAssociations; reload(BuildAffymetrixAssociations)
import FilterDabg; reload(FilterDabg)

use_Tkinter = 'no'
try:
    from Tkinter import *
    use_Tkinter = 'yes'
except ImportError: use_Tkinter = 'yes'; print "\nPmw or Tkinter not found... Tkinter print out not available";
debug_mode = 'no'

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

def calculate_expression_measures(expr_input_dir,expr_group_dir,experiment_name,comp_group_dir,probeset_db,annotate_db):
    print "Processing the expression file:",expr_input_dir
    
    global array_fold_headers; global summary_filtering_stats; global raw_data_comp_headers; global array_folds
    fn1=filepath(expr_input_dir)
    x = 0; y = 0; d = 0
    array_folds={}
    for line in open(fn1,'rU').xreadlines():             
      data = cleanUpLine(line)
      if data[0] != '#':
        fold_data = string.split(data,'\t'); arrayid = fold_data[0]
        #if 'counts.' in expr_input_dir: arrayid,coordinates = string.split(arrayid,'=') ### needed for exon-level analyses only
        ### differentiate data from column headers
        if x == 1:
            fold_data = fold_data[1:]; fold_data2=[]
            for fold in fold_data:
                fold = string.replace(fold,'"','')
                try:
                    if len(fold)>0: fold = float(fold); fold_data2.append(fold)
                    else: null = float('a') ###Force a ValueError since no data is present in this cell
                except ValueError: 
                    print_out = 'WARNING!!! The ID'+arrayid+ 'has an invalid expression value:'+fold+'\n. Correct and re-run'
                    try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); sys.exit()
                    except NameError: print print_out; sys.exit()
            if expression_data_format == 'non-log':
                fold_data3=[] ###Convert numeric expression to log fold (previous to version 2.05 added 1)
                for fold in fold_data2:
                    try:
                        log_fold = math.log((float(fold)),2) ### changed from - log_fold = math.log((float(fold)+1),2) - version 2.05
                        fold_data3.append(log_fold)
                    except ValueError:  ###Not an ideal situation: Value is negative - Convert to zero
                        if float(fold)<=0: math.log(1,2); fold_data3.append(log_fold)
                        else:
                            print_out = 'WARNING!!! The ID'+arrayid+ 'has an invalid expression value:'+fold+'\n. Correct and re-run'
                            try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); sys.exit()
                            except NameError: print print_out; sys.exit()
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
    
    try:
        array_folds, array_fold_headers, summary_filtering_stats,raw_data_comp_headers = reorder_arrays.reorder(array_folds,array_names,expr_group_list,comp_group_list,probeset_db,include_raw_data,array_type,normalization_method,fl)
    except Exception: 
        print_out = 'AltAnalyze encountered an error with the format of the expression file.\nIf the data was designated as log intensities and it is not, then re-run as non-log.'
        try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); root.destroy(); force_exit ### Forces the error log to pop-up
        except NameError: print print_out; sys.exit()

    ### Integrate maximum counts for each gene for the purpose of filtering (RNASeq data only)
    if array_type == 'RNASeq' and 'counts.' not in expr_input_dir: addMaxReadCounts(expr_input_dir)
    
    ### Export these results to a DATASET statistics and annotation results file
    if 'exp.' in expr_input_dir:
        if array_type == 'RNASeq' and norm == 'RPKM': filterRNASeq(count_statistics_db)
        exportAnalyzedData(comp_group_list2,expr_group_db)
    
        ### Export formatted results for input as an expression dataset into GenMAPP or PathVisio
        if data_type == 'expression':
            if include_raw_data == 'yes': headers = removeRawData(array_fold_headers)
            else: headers = array_fold_headers
            exportDataForGenMAPP(headers)
            
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
    max_count_db = RNASeq.importGeneCounts(filename)
    for gene in summary_filtering_stats:
        gs = summary_filtering_stats[gene]
        gs.setMaxCount(max_count_db[gene]) ### Shouldn't cause an error, but we want to get an exception if it does (something is wrong with the analysis)        
        
def importArrayGroups(expr_group_dir,array_linker_db):
    new_index_order = 0    
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
                    except KeyError:
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
        
    except Exception:
        exportArrayHeaders(expr_group_dir,array_linker_db)
        print_out = 'No groups or comps files found for'+expr_group_dir+'... exiting program.'
        try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); root.destroy(); sys.exit()
        except Exception: print print_out; sys.exit()
        
    return expr_group_list,expr_group_db

def exportArrayHeaders(expr_group_dir,array_linker_db):
    new_file = string.replace(expr_group_dir,'groups.','arrays.')
    new_file = string.replace(new_file,'exp.','arrays.')
    new_file = string.replace(new_file,'counts.','arrays.')
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
            
def exportDataForGenMAPP(headers):
    ###Export summary columns for GenMAPP analysis
    systems = importSystemCodes()
    GenMAPP_file = expression_dataset_output_dir + 'GenMAPP-'+experiment_name+'.txt'
    try: genmapp = export.createExportFile(GenMAPP_file,expression_dataset_output_dir[:-1])
    except RuntimeError:
        export.isFileOpen(GenMAPP_file,expression_dataset_output_dir[:-1])
        genmapp = export.createExportFile(GenMAPP_file,expression_dataset_output_dir[:-1])
    if array_type == "3'array":
        if vendor == 'Affymetrix': system_code = 'X'
        if vendor == 'Illumina': system_code = 'Il'
        if vendor == 'Agilent': system_code = 'Ag'
        if vendor == 'Codelink': system_code = 'Co'
    elif array_type != 'AltMouse': system_code = 'En'
    else:
        try: system_code = systems[vendor]
        except Exception: system_code = 'X'
    genmapp_title = ['GeneID','SystemCode'] + headers
    genmapp_title = string.join(genmapp_title,'\t')+'\t'+'ANOVA-rawp'+'\t'+'ANOVA-adjp'+'\t'+'largest fold'+'\n'
    genmapp.write(genmapp_title)

    for probeset in array_folds:
        data_val = probeset+'\t'+system_code
        for value in array_folds[probeset]: data_val += '\t'+ str(value)
        gs = summary_filtering_stats[probeset]
        data_val += '\t'+ str(gs.Pval()) +'\t'+ str(gs.AdjP()) +'\t'+ str(gs.LogFold()) +'\n'
        genmapp.write(data_val)
    genmapp.close()
    exportGOEliteInput(headers,system_code)
    print 'Exported GO-Elite input files...'

def buildCriterion(ge_fold_cutoffs, ge_pvalue_cutoffs, ge_ptype, main_output_folder, operation):
    global array_folds; global m_cutoff; global p_cutoff; global expression_dataset_output_dir
    global ptype_to_use
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
                    array_folds[probeset] = values
            if operation == 'summary':
                exportGeneRegulationSummary(headers,system_code)
            else:
                input_files_exported = exportGOEliteInput(headers,system_code)
    array_folds=[]
    
def exportGOEliteInput(headers,system_code):
    ### Filter statistics based on user-defined thresholds as input for GO-Elite analysis
    criterion_db={}; denominator_geneids={}; index = 0; ttest=[]
    for column in headers:
        if ptype_to_use in column and 'ANOVA' not in column: ttest.append(index)
        index+=1
        
    for probeset in array_folds:
        index = 0; af = array_folds[probeset]
        for value in array_folds[probeset]:
            denominator_geneids[probeset]=[]
            if index in ttest:
                criterion_name = headers[index][5:]
                try: log_fold = float(af[index-2])
                except Exception: log_fold = 0 ### Occurs when a fold change is annotated as 'Insufficient Expression'
                try: p_value = float(value)
                except Exception: p_value = 1 ### Occurs when a p-value is annotated as 'Insufficient Expression'
                if abs(log_fold)>m_cutoff and (p_value<p_cutoff or p_value==1):
                    try: criterion_db[criterion_name].append((probeset,log_fold,p_value))
                    except KeyError: criterion_db[criterion_name] = [(probeset,log_fold,p_value)]
                    if log_fold>0:
                        try: criterion_db[criterion_name+'-upregulated'].append((probeset,log_fold,p_value))
                        except KeyError: criterion_db[criterion_name+'-upregulated'] = [(probeset,log_fold,p_value)]
                    else:
                        try: criterion_db[criterion_name+'-downregulated'].append((probeset,log_fold,p_value))
                        except KeyError: criterion_db[criterion_name+'-downregulated'] = [(probeset,log_fold,p_value)]                    
            index += 1
    if len(criterion_db)>0:
        ### Export denominator gene IDs
        input_files_exported = 'yes'
        expression_dir = string.replace(expression_dataset_output_dir,'ExpressionOutput/','')
        goelite_file = expression_dir +'GO-Elite/denominator/GE.denominator.txt'
        goelite = export.createExportFile(goelite_file,expression_dir+'GO-Elite/denominator')        
        goelite_title = ['GeneID','SystemCode']
        goelite_title = string.join(goelite_title,'\t')+'\n'; goelite.write(goelite_title)
        for probeset in denominator_geneids:
            values = string.join([probeset,system_code],'\t')+'\n'; goelite.write(values)
        goelite.close()

        ### Export criterion gene IDs and minimal data     
        for criterion_name in criterion_db:
            if criterion_name[-1] == ' ': criterion_file_name = criterion_name[:-1]
            else: criterion_file_name = criterion_name
            goelite_file = expression_dir + 'GO-Elite/input/GE.'+criterion_file_name+'.txt'
            goelite = export.createExportFile(goelite_file,expression_dir+'GO-Elite/input')        
            goelite_title = ['GeneID','SystemCode',criterion_name+'-log_fold',criterion_name+'-p_value']
            goelite_title = string.join(goelite_title,'\t')+'\n'; goelite.write(goelite_title)
            for (probeset,log_fold,p_value) in criterion_db[criterion_name]:
                values = string.join([probeset,system_code,str(log_fold),str(p_value)],'\t')+'\n'
                goelite.write(values)
            goelite.close()
    else: input_files_exported = 'no'
    return input_files_exported

def exportGeneRegulationSummary(headers,system_code):
    """Perform a series of targetted queries to report the number of coding and non-coding genes expressed
    along with various regulation and annotation parameters"""
    
    criterion_db={}; detected_exp_db={}; denominator_geneids={}; index = 0; ttest=[]; avg_columns=[]; all_criterion=[]; all_groups=[]
    search_miR = 'miR-1('
    coding_types = ['protein_coding','ncRNA']
    
    for column in headers:
        if ptype_to_use in column and 'ANOVA' not in column: ttest.append(index)
        index+=1
        if 'Protein Classes' in column: pc = index-1
        if 'microRNA' in column: mi = index-1
        if 'avg-' in column: avg_columns.append(index-1)
        
    for probeset in array_folds:
        index = 0; af = array_folds[probeset]
        for value in array_folds[probeset]:
            denominator_geneids[probeset]=[]
            if index in avg_columns:
                group_name = headers[index][4:]
                try: protein_class = af[pc]
                except Exception: protein_class = 'NULL'
                if float(af[index])>0:
                    if group_name not in all_groups: all_groups.append(group_name)
                    if 'protein_coding' in protein_class:
                        try: detected_exp_db[group_name,'protein_coding']+=1
                        except KeyError: detected_exp_db[group_name,'protein_coding']=1
                    else:
                        try: detected_exp_db[group_name,'ncRNA']+=1
                        except KeyError: detected_exp_db[group_name,'ncRNA']=1
            if index in ttest:
                criterion_name = headers[index][5:]
                try: log_fold = float(af[index-2])
                except Exception: log_fold = 0 ### Occurs when a fold change is annotated as 'Insufficient Expression'
                try: p_value = float(value)
                except Exception: p_value = 1 ### Occurs when a p-value is annotated as 'Insufficient Expression'
                try: protein_class = af[pc]
                except Exception: protein_class = 'NULL'
                if abs(log_fold)>m_cutoff and (p_value<p_cutoff or p_value==1):
                    if criterion_name not in all_criterion: all_criterion.append(criterion_name)
                    try: criterion_db[criterion_name]+=1
                    except KeyError: criterion_db[criterion_name] = 1
                    
                    if 'protein_coding' in protein_class:
                        if log_fold>0:
                            try: criterion_db[criterion_name,'upregulated','protein_coding']+=1
                            except KeyError: criterion_db[criterion_name,'upregulated','protein_coding'] = 1
                            if 'miR-1(' in af[mi]:
                                try: criterion_db[criterion_name,'upregulated','protein_coding',search_miR[:-1]]+=1
                                except KeyError: criterion_db[criterion_name,'upregulated','protein_coding',search_miR[:-1]] = 1
                        else:
                            try: criterion_db[criterion_name,'downregulated','protein_coding']+=1
                            except KeyError: criterion_db[criterion_name,'downregulated','protein_coding'] = 1
                            if 'miR-1(' in af[mi]:
                                try: criterion_db[criterion_name,'downregulated','protein_coding',search_miR[:-1]]+=1
                                except KeyError: criterion_db[criterion_name,'downregulated','protein_coding',search_miR[:-1]] = 1
                    else:
                        if log_fold>0:
                            try: criterion_db[criterion_name,'upregulated','ncRNA']+=1
                            except KeyError: criterion_db[criterion_name,'upregulated','ncRNA'] = 1
                            if 'miR-1(' in af[mi]:
                                try: criterion_db[criterion_name,'upregulated','ncRNA',search_miR[:-1]]+=1
                                except KeyError: criterion_db[criterion_name,'upregulated','ncRNA',search_miR[:-1]] = 1
                        else:
                            try: criterion_db[criterion_name,'downregulated','ncRNA']+=1
                            except KeyError: criterion_db[criterion_name,'downregulated','ncRNA'] = 1
                            if 'miR-1(' in af[mi]:
                                try: criterion_db[criterion_name,'downregulated','ncRNA',search_miR[:-1]]+=1
                                except KeyError: criterion_db[criterion_name,'downregulated','ncRNA',search_miR[:-1]] = 1

            index += 1

    if len(criterion_db)>0:

        summary_path = expression_dataset_output_dir +'summary.txt'
        export_data = export.ExportFile(summary_path)
        print 'Export summary gene expression results to:',summary_path

        ### Output Number of Expressed Genes        
        title = ['Biological group']
        for group_name in all_groups: title.append(group_name)
        title = string.join(title,'\t')+'\n'; export_data.write(title)

        for coding_type in coding_types:
            if coding_type == 'protein_coding': values = ['Expressed protein-coding genes']
            else: values = ['Expressed ncRNAs']
            for group in all_groups:
                for group_name in detected_exp_db:
                    if group in group_name and coding_type in group_name:
                        values.append(str(detected_exp_db[group_name]))
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

def exportAnalyzedData(comp_group_list2,expr_group_db):
    report = 'multiple'; report = 'single'
    try: ensembl_microRNA_db = importMicrornaAssociations(species,report)
    except IOError: ensembl_microRNA_db={}
    if array_type != "AltMouse" and array_type != "3'array":
        try:
            import EnsemblImport
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
        if array_type != "AltMouse" and array_type != "3'array" :
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
            elif array_type != 'AltMouse' and array_type != "3'array":
                try:
                    try: definition = annotate_db[arrayid][0]; symbol = annotate_db[arrayid][1]; rna_processing = annotate_db[arrayid][2]
                    except TypeError: print arrayid, annotate_db[arrayid]; kill
                except KeyError: definition=''; symbol=''; rna_processing=''
                report = 'all'
                try: miRs = ensembl_microRNA_db[arrayid]
                except KeyError: miRs = ''
                trans_cluster = probeset_db[arrayid][0]
                exon_ids = probeset_db[arrayid][1]
                probesets = probeset_db[arrayid][2]
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
            elif array_type == "3'array":
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

def avg(array):
    denominator = len(array)
    total = float(sum(array))
    average = total/denominator
    return average

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

def performTissueProfiling(expr_input_dir):
    try:
        import TissueProfiler
        compendium_platform = 'gene'
        exp_output = expression_dataset_output_dir + 'DATASET-'+experiment_name+'.txt'
        TissueProfiler.runTissueProfiler(species,array_type,expr_input_dir,exp_output,'protein_coding',compendium_platform)
    except IOError: null=[] ### Analysis may not be supported for species or data is incompatible
    
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
  global arrayCode; arrayCode = 0; global probability_statistic; global fl

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
  
  print "Beginning to Process the",species,array_type,'dataset'
  
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
  elif array_type == "3'array":
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
          conventional_array_db = BuildAffymetrixAssociations.getArrayAnnotationsFromGOElite(conventional_array_db,species,vendor,use_go)
      except Exception: print "Additional annotation import failed"
      print len(conventional_array_db), "Array IDs with annotations from",vendor,"annotation files imported."
  elif array_type != "AltMouse":
      probeset_db = []; annotate_db = []; constitutive_db = []; conventional_array_db = []
      ### The below function gathers GO annotations from the GO-Elite database (not Affymetrix as the module name implies)
      conventional_array_db = BuildAffymetrixAssociations.getEnsemblAnnotationsFromGOElite(species)
                
  altanalyze_files = []; datasets_with_all_necessary_files=0
  for dataset in exp_file_location_db:
      experiment_name = string.replace(dataset,'exp.',''); experiment_name = string.replace(experiment_name,'.txt','')
      fl = exp_file_location_db[dataset]
      expr_input_dir = fl.ExpFile()
      stats_input_dir = fl.StatsFile()
      expr_group_dir = fl.GroupsFile()
      comp_group_dir = fl.CompsFile()
      try: norm = fl.FeatureNormalization()
      except Exception: norm = 'NA'
      try: probability_statistic = fl.ProbabilityStatistic()
      except Exception: probability_statistic = 'unpaired t-test'
      residuals_input_dir = string.replace(expr_input_dir,'exp.','residuals.')
      root_dir = fl.RootDir()
      datasets_with_all_necessary_files +=1
      checkArrayHeaders(expr_input_dir,expr_group_dir)
      expression_dataset_output_dir = root_dir+"ExpressionOutput/"
      if array_type != "3'array": #array_type != 'AltMouse' and 
          try: probeset_db,annotate_db,comparison_filename_list = ExonArray.getAnnotations(fl,array_type,dabg_p,expression_threshold,data_source,vendor,constitutive_source,species,avg_all_for_ss,filter_by_dabg,perform_alt_analysis,expression_data_format)
          except Exception, e:
              print_out = 'Error ecountered for the '+species+', '+array_type+' dataset. Check to ensure that:\n(1) the correct platform and species were selected and\n(2) some expression values are present in ExpressionInput/exp.YourDataset.txt'
              try: UI.WarningWindow(print_out,'Critical Error - Exiting Program!!!'); root.destroy(); sys.exit()
              except Exception: print print_out; sys.exit()
          if array_type != 'AltMouse': expr_input_dir = expr_input_dir[:-4]+'-steady-state.txt'
          else: probeset_db = original_probeset_db; annotate_db = original_annotate_db
          for file in comparison_filename_list: altanalyze_files.append(file)
          residual_file_status = ExonArray.verifyFile(residuals_input_dir)
          ### Separate residual file into comparison files for AltAnalyze (if running FIRMA)
          if residual_file_status == 'found': ExonArray.processResiduals(fl,Array_type,Species,perform_alt_analysis)
      if norm == 'RPKM':
        ### Separately analyze steady-state counts first, to replace fold changes
        counts_expr_dir = string.replace(expr_input_dir,'exp.','counts.')
        count_statistics_db, count_statistics_headers = calculate_expression_measures(counts_expr_dir,expr_group_dir,experiment_name,comp_group_dir,probeset_db,annotate_db)
      calculate_expression_measures(expr_input_dir,expr_group_dir,experiment_name,comp_group_dir,probeset_db,annotate_db)
      #performTissueProfiling(expr_input_dir) ### Correlate gene-level expression values with known cells and tissues
  annotate_db={}; probeset_db={}; constitutive_db={}; array_fold_db={}; raw_data_comps={}; conventional_array_db=[]
  clearObjectsFromMemory(conventional_array_db); conventional_array_db=[]
  try: clearObjectsFromMemory(summary_filtering_stats); summary_filtering_stats=[]
  except Exception: null=[]
  try: clearObjectsFromMemory(array_folds); array_folds=[]
  except Exception: null=[]
  try: clearObjectsFromMemory(count_statistics_db); count_statistics_db=[]
  except Exception: null=[]
  
  #print 'after deleted'; returnLargeGlobalVars()

  ### Code in progress for version 2.1
  try: buildCriterion(GE_fold_cutoffs, p_cutoff, ptype_to_use, root_dir+'/ExpressionOutput/','summary') ###Outputs a summary of the dataset and all comparisons to ExpressionOutput/summary.txt
  except Exception: null=[]
  
  if datasets_with_all_necessary_files == 0:
      ###Thus no files were found with valid inputs for all file types
      print 'WARNING....No propperly named datasets were found. ExpressionBuilder requires that there are at least 3 files with the prefixes "exp.", "groups." and "comps.", with the following dataset name being identical with all three files.'
      print "...check these file names before running again."
      inp = sys.stdin.readline(); sys.exit()
  altanalyze_files = unique.unique(altanalyze_files) ###currently not used, since declaring altanalyze_files a global is problematic (not available from ExonArray... could add though)
  if array_type != "3'array" and perform_alt_analysis != 'expression':
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

if __name__ == '__main__':
  buildCriterion(2, 0.05, 'rawp', 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/Kathy-RNASeq/RPKM-normalization/ExpressionOutput/','summary')
  """
  dabg_p = 0.75
  a = "3'array"; b = "exon"; c = "AltMouse"; e = "custom"; Array_type = c
  l = 'log'; n = 'non-log'; Expression_data_format = l
  w = 'Agilent'; x = 'Affymetrix'; y = 'Ensembl'; z = 'default'; data_source = y; constitutive_source = z; vendor = x
  hs = 'Hs'; mm = 'Mm'; dr = 'Dr'; rn = 'Rn'; Species = mm
  Include_raw_data = 'yes'  
  expression_threshold = 0 ### Based on suggestion from BMC Genomics. 2006 Dec 27;7:325. PMID: 17192196, for hu-exon 1.0 st array
  avg_all_for_ss = 'no'  ###Default is 'no' since we don't want all probes averaged for the exon arrays
  
  remoteExpressionBuilder(Species,Array_type,dabg_p,expression_threshold,avg_all_for_ss,Expression_data_format,vendor,constitutive_source,data_source,Include_raw_data)
  """
