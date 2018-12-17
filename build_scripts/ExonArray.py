###ExonArray
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

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os.path
import unique
from stats_scripts import statistics
import math
from build_scripts import EnsemblImport; reload(EnsemblImport)
from build_scripts import ExonArrayEnsemblRules; reload(ExonArrayEnsemblRules)
from build_scripts import ExonArrayAffyRules
import ExpressionBuilder
import reorder_arrays
import time
import export
import traceback

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    #add in code to prevent folder names from being included
    dir_list2 = []
    for file in dir_list:
        if '.txt' in file: dir_list2.append(file)
    return dir_list2

################# Begin Analysis

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
            if len(globals()[var])>500:
                print var, len(globals()[var])
        except Exception: null=[]

def clearObjectsFromMemory(db_to_clear):
    db_keys={}
    for key in db_to_clear: db_keys[key]=[]
    for key in db_keys: del db_to_clear[key]
    
def exportMetaProbesets(array_type,species):
    import AltAnalyze; reload(AltAnalyze)
    import export
    probeset_types = ['core','extended','full']
    if array_type == 'junction': probeset_types = ['all']
    for probeset_type in probeset_types:
        exon_db,null = AltAnalyze.importSplicingAnnotations(array_type,species,probeset_type,'yes','')
        gene_db={}; null=[]
        for probeset in exon_db:
            ### At this point, exon_db is filtered by the probeset_type (e.g., core)
            ensembl_gene_id = exon_db[probeset].GeneID()
            try: gene_db[ensembl_gene_id].append(probeset)
            except Exception: gene_db[ensembl_gene_id] = [probeset]
            
        exon_db=[]; uid=0
        output_dir = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_'+array_type+'_'+probeset_type+'.mps'
        #output_cv_dir = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Conversion_'+array_type+'_'+probeset_type+'.txt'
        #data_conversion = export.ExportFile(output_cv_dir)
        data = export.ExportFile(output_dir)
        data.write('probeset_id\ttranscript_cluster_id\tprobeset_list\tprobe_count\n')
        print "Exporting",len(gene_db),"to",output_dir
        for ensembl_gene_id in gene_db:
            probeset_strlist = string.join(gene_db[ensembl_gene_id],' '); uid+=1
            line = string.join([str(uid),str(uid),probeset_strlist,str(len(gene_db[ensembl_gene_id])*4)],'\t')+'\n'
            data.write(line)
            #conversion_line = string.join([str(uid),ensembl_gene_id],'\t')+'\n'; data_conversion.write(conversion_line)
        data.close(); #data_conversion.close()

def adjustCounts(exp_vals):
    exp_vals2=[]
    for i in exp_vals: exp_vals2.append(int(i)+1) ### Increment the rwcounts by 1
    return exp_vals

def remoteExonProbesetData(filename,import_these_probesets,import_type,platform):
    global array_type
    array_type = platform
    results = importExonProbesetData(filename,import_these_probesets,import_type)
    return results

def importExonProbesetData(filename,import_these_probesets,import_type):
    """This is a powerfull function, that allows exon-array data import and processing on a line-by-line basis,
    allowing the program to immediately write out data or summarize it without storing large amounts of data."""

    fn=filepath(filename); start_time = time.time()
    exp_dbase={}; filtered_exp_db={}; ftest_gene_db={}; filtered_gene_db={}; probeset_gene_db={}; biotypes={}
    d = 0; x = 0
    
    if 'stats.' in filename: filetype = 'dabg'
    else:
        filetype = 'expression'
        if 'FullDatasets' in filename and import_type == 'filterDataset':
            output_file = string.replace(filename,'.txt','-temp.txt')
            temp_data = export.ExportFile(output_file)       
                                         
    ###Import expression data (non-log space)
    
    if 'counts.' in filename: counts = 'yes'
    else: counts = 'no'
    try:
        for line in open(fn,'rU').xreadlines():             
          data = cleanUpLine(line);
          if len(data)==0: null=[]
          elif data[0] != '#' and x == 1:   ###Grab expression values
            tab_delimited_data = string.split(data,'\t')
            probeset = tab_delimited_data[0]
            if '=' in probeset: probeset = string.split(probeset,'=')[0]
            if '' in tab_delimited_data or '0' in tab_delimited_data and counts == 'no':
                None ### Rare GEO datasets remove values from the exon-level data (exclude the whole probeset)
            elif import_type == 'raw':
                try:
                    null = import_these_probesets[probeset]; exp_vals = tab_delimited_data[1:]; exp_dbase[probeset] = exp_vals
                except KeyError: null = [] ###Don't import any probeset data
            elif import_type == 'filterDataset':
                try:
                    null = import_these_probesets[probeset]; temp_data.write(line)
                except KeyError: null = [] ###Don't import any probeset data                
            elif import_type == 'reorderFilterAndExportAll':
                if '-' in probeset: biotypes['junction'] = []
                else: biotypes['exon'] = []
                try:
                    ###For filtering, don't remove re-organized entries but export filtered probesets to another file
                    if exp_analysis_type == 'expression': null = import_these_probesets[probeset]
                    exp_vals = tab_delimited_data[1:]
                    #if counts == 'yes': exp_vals = adjustCounts(exp_vals)
                    filtered_exp_db={}; filtered_exp_db[probeset] = exp_vals
                    reorderArraysOnly(filtered_exp_db,filetype,counts) ###order and directly write data
                except KeyError: null = [] ###Don't import any probeset data
          elif data[0] != '#' and x == 0:  ###Grab labels
              array_names = []; array_linker_db = {}; z = 0
              tab_delimited_data = string.split(data,'\t')
              for entry in tab_delimited_data:
                  if z != 0: array_names.append(entry)
                  z += 1
              for array in array_names: #use this to have an orignal index order of arrays
                  array = string.replace(array,'\r','') ###This occured once... not sure why
                  array_linker_db[array] = d; d +=1
              x += 1
              ### Process and export expression dataset headers
              if import_type == 'reorderFilterAndExportAll':
                  if filetype == 'expression':
                      headers = tab_delimited_data[1:]; probeset_header = tab_delimited_data[0]
                      filtered_exp_db[probeset_header] = headers
                      reorderArraysHeader(filtered_exp_db)
              elif import_type == 'filterDataset': temp_data.write(line)
    except IOError:
        #print traceback.format_exc()
        print filename, 'not found.'
        null=[]

    end_time = time.time(); time_diff = int(end_time-start_time)
    print "Exon data imported in %d seconds" % time_diff
    if array_type == 'RNASeq': id_name = 'junction IDs'
    else: id_name = 'array IDs'
    if import_type == 'filterDataset': temp_data.close()
    if import_type == 'arraynames': return array_linker_db,array_names
    if import_type == 'raw':
        print len(exp_dbase),id_name,"imported with expression values"
        return exp_dbase
    else:
        return biotypes

def exportGroupedComparisonProbesetData(filename,probeset_db,data_type,array_names,array_linker_db,perform_alt_analysis):
        """This function organizes the raw expression data into sorted groups, exports the organized data for all conditions and comparisons
        and calculates which probesets have groups that meet the user defined dabg and expression thresholds."""
        #comparison_filename_list=[]
        #if perform_alt_analysis != 'expression': ### User Option (removed in version 2.0 since the option prevented propper filtering)
        comparison_filename_list=[]
        probeset_dbase={}; exp_dbase={}; constitutive_gene_db={}; probeset_gene_db={} ### reset databases to conserve memory
        global expr_group_list; global comp_group_list; global expr_group_db
        if data_type == 'residuals':
            expr_group_dir = string.replace(filename,'residuals.','groups.')
            comp_group_dir = string.replace(filename,'residuals.','comps.')
        elif data_type == 'expression':
            expr_group_dir = string.replace(filename,'exp.','groups.')
            comp_group_dir = string.replace(filename,'exp.','comps.')
            if 'counts.' in filename:
                expr_group_dir = string.replace(expr_group_dir,'counts.','groups.')
                comp_group_dir = string.replace(comp_group_dir,'counts.','comps.')
                data_type = 'counts'
        elif data_type == 'dabg':
            expr_group_dir = string.replace(filename,'stats.','groups.')
            comp_group_dir = string.replace(filename,'stats.','comps.')

        comp_group_list, comp_group_list2 = ExpressionBuilder.importComparisonGroups(comp_group_dir)
        expr_group_list,expr_group_db = ExpressionBuilder.importArrayGroups(expr_group_dir,array_linker_db)

        print "Reorganizing expression data into comparison groups for export to down-stream splicing analysis software"
        ###Do this only for the header data
        group_count,raw_data_comp_headers = reorder_arrays.reorderArrayHeaders(array_names,expr_group_list,comp_group_list,array_linker_db)

        ###Export the header info and store the export write data for reorder_arrays
        global comparision_export_db; comparision_export_db={}; array_type_name = 'Exon'
        if array_type == 'junction': array_type_name = 'Junction'
        elif array_type == 'RNASeq': array_type_name = 'RNASeq'
        if data_type != 'residuals': AltAnalzye_input_dir = root_dir+"AltExpression/pre-filtered/"+data_type+'/'
        else: AltAnalzye_input_dir = root_dir+"AltExpression/FIRMA/residuals/"+array_type+'/'+species+'/' ### These files does not need to be filtered until AltAnalyze.py

        for comparison in comp_group_list2: ###loop throught the list of comparisons
            group1 = comparison[0]; group2 = comparison[1]
            group1_name = expr_group_db[group1]; group2_name = expr_group_db[group2]
            comparison_filename = species+'_'+array_type_name+'_'+ group1_name + '_vs_' + group2_name + '.txt'
                
            new_file = AltAnalzye_input_dir + comparison_filename; comparison_filename_list.append(comparison_filename)
            data = export.createExportFile(new_file,AltAnalzye_input_dir[:-1])

            try: array_names = raw_data_comp_headers[comparison]
            except KeyError: print raw_data_comp_headers;kill
            title = ['UID']+array_names; title = string.join(title,'\t')+'\n'; data.write(title)
            comparision_export_db[comparison] = data ###store the export file write data so we can write after organizing
        #print filename, normalize_feature_exp
        biotypes = importExonProbesetData(filename,probeset_db,'reorderFilterAndExportAll')
        
        if normalize_feature_exp == 'RPKM': ### Add the gene-level RPKM data (this is in addition to the counts. file)
            exp_gene_db={}
            for i in probeset_db: exp_gene_db[probeset_db[i][0]]=[]
            filename = string.replace(filename,'.txt','-steady-state.txt')
            #print filename, normalize_feature_exp, 'here'
            importExonProbesetData(filename,exp_gene_db,'reorderFilterAndExportAll')
            
        for comparison in comparision_export_db:
            data = comparision_export_db[comparison]; data.close()
        print "Pairwise comparisons for AltAnalyze exported..."
        try: fulldataset_export_object.close()
        except Exception: null=[]
        return comparison_filename_list, biotypes

def filterExpressionData(filename,pre_filtered_db,constitutive_gene_db,probeset_db,data_type,array_names,perform_alt_analysis):
    """Probeset level data and gene level data are handled differently by this program for exon and tiling based arrays.
    Probeset level data is sequentially filtered to reduce the dataset to a minimal set of expressed probesets (exp p-values) that
    that align to genes with multiple lines of evidence (e.g. ensembl) and show some evidence of regulation for any probesets in a
    gene (transcript_cluster). Gene level analyses are handled under the main module of the program, exactly like 3' arrays, while
    probe level data is reorganized, filtered and output from this module."""
          
    ###First time we import, just grab the probesets associated
    if array_names != 'null': ### Then there is only an expr. file and not a stats file
        ###Identify constitutive probesets to import (sometimes not all probesets on the array are imported)
        possible_constitutive_probeset={}; probeset_gene_db={}
        for probeset in probeset_db:
            try:
                probe_data = probeset_db[probeset]
                gene = probe_data[0]; affy_class = probe_data[-1]; external_exonid = probe_data[-2]
                if affy_class == 'core' or len(external_exonid)>2: ### These are known exon only (e.g., 'E' probesets)
                    proceed = 'yes'
                    if array_type == 'RNASeq' and 'exon' in biotypes: ### Restrict the analysis to exon RPKM or count data for constitutive calculation
                        if '-' in probeset: proceed = 'no'
                    elif array_type == 'RNASeq' and 'junction' in biotypes:
                        if '-' not in probeset: proceed = 'no' ### Use this option to override 
                    if proceed == 'yes':
                        try: probeset_gene_db[gene].append(probeset)
                        except KeyError: probeset_gene_db[gene] = [probeset]
                        possible_constitutive_probeset[probeset] = []
            except KeyError: null = []

        ### Only import probesets that can be used to calculate gene expression values OR link to gene annotations (which have at least one dabg p<0.05 for all samples normalized)
        constitutive_exp_dbase = importExonProbesetData(filename,possible_constitutive_probeset,'raw') 
        generateConstitutiveExpression(constitutive_exp_dbase,constitutive_gene_db,probeset_gene_db,pre_filtered_db,array_names,filename)
        
        if array_type != 'RNASeq':
            ### Repeat the analysis to export gene-level DABG reports for LineageProfiler (done by re-calling the current function for RNASeq and counts)
            try:
                constitutive_exp_dbase = importExonProbesetData(stats_input_dir,possible_constitutive_probeset,'raw')
                generateConstitutiveExpression(constitutive_exp_dbase,constitutive_gene_db,probeset_gene_db,pre_filtered_db,array_names,stats_input_dir)
            except Exception:
                print 'No dabg p-value specified for analysis (must be supplied in command-line or GUI mode - e.g., stats.experiment.txt)'
                None ### Occurs when no stats_input_dir specified
        
        constitutive_exp_dbase = {}; possible_constitutive_probeset={}; pre_filtered_db={}; probeset_db={}; constitutive_gene_db={} ### reset databases to conserve memory
        probeset_gene_db={}
        
        """
        print 'global vars'
        returnLargeGlobalVars()
        print 'local vars'
        all = [var for var in locals() if (var[:2], var[-2:]) != ("__", "__")]
        for var in all:
                try:
                    if len(locals()[var])>500: print var, len(locals()[var])
                except Exception: null=[]
        """
        
def reorderArraysHeader(filtered_exp_db):
    ### These are all just headers from the first line
    for probeset in filtered_exp_db:
        grouped_ordered_array_list = {}; group_list = []
        for x in expr_group_list:
            y = x[1]; group = x[2]  ### this is the new first index
            try: new_item = filtered_exp_db[probeset][y]
            except Exception:
                print y,group, probeset
                print 'Prior counts.YourExperiment-steady-state.txt exists and has different samples... delete before proceeding\n'
                bad_exit
            try: grouped_ordered_array_list[group].append(new_item)
            except KeyError: grouped_ordered_array_list[group] = [new_item]
        for group in grouped_ordered_array_list: group_list.append(group)
        group_list.sort(); combined_value_list=[]
        for group in group_list:
            group_name = expr_group_db[str(group)]
            g_data2 = []; g_data = grouped_ordered_array_list[group]
            for header in g_data: g_data2.append(group_name+':'+header)
            combined_value_list+=g_data2
        values = string.join([probeset]+combined_value_list,'\t')+'\n'
        fulldataset_export_object.write(values)
        
def reorderArraysOnly(filtered_exp_db,filetype,counts): 
    ###array_order gives the final level order sorted, followed by the original index order as a tuple                   
    ###expr_group_list gives the final level order sorted, followed by the original index order as a tuple
    for probeset in filtered_exp_db:
        grouped_ordered_array_list = {}; group_list = []
        for x in expr_group_list:
            y = x[1]; group = x[2]  ### this is the new first index
            ### for example y = 5, therefore the filtered_exp_db[probeset][5] entry is now the first
            try:
                try: new_item = filtered_exp_db[probeset][y]
                except TypeError: print y,x,expr_group_list; kill
            except IndexError: print probeset,y,x,expr_group_list,'\n',filtered_exp_db[probeset];kill
            ###Used for comparision analysis
            try: grouped_ordered_array_list[group].append(new_item)
            except KeyError: grouped_ordered_array_list[group] = [new_item]
            
        ### For the exon-level expression data, export the group pair data for all pairwise comparisons to different comp files
        ###*******Include a database with the raw values saved for permuteAltAnalyze*******
        for info in comp_group_list:
            group1 = int(info[0]); group2 = int(info[1]); comp = str(info[0]),str(info[1])
            g1_data = grouped_ordered_array_list[group1]
            g2_data = grouped_ordered_array_list[group2]
            #print probeset, group1, group2, g1_data, g2_data, info;kill
            data = comparision_export_db[comp]
            values = [probeset]+g2_data+g1_data; values = string.join(values,'\t')+'\n' ###groups are reversed since so are the labels
            #raw_data_comps[probeset,comp] = temp_raw
            data.write(values)
            
        ### Export all values grouped from the array
        for group in grouped_ordered_array_list: group_list.append(group)
        group_list.sort(); combined_value_list=[]; avg_values=[]
        for group in group_list:
            g_data = grouped_ordered_array_list[group]
            if exp_analysis_type == 'expression':
                try: avg_gdata = statistics.avg(g_data); avg_values.append(avg_gdata)
                except Exception:
                    print g_data
                    print avg_values
                    kill
            combined_value_list+=g_data
        
        if exp_data_format == 'non-log' and counts == 'no':
            try: combined_value_list = logTransform(combined_value_list)
            except Exception:
                print probeset, combined_value_list,comp_group_list,expr_group_list
                print filtered_exp_db[probeset]; kill

        if filetype == 'expression':
            ### Export the expression values for all samples grouped (if meeting the above thresholds)
            values = string.join([probeset]+combined_value_list,'\t')+'\n'
            fulldataset_export_object.write(values) ### Don't need this for dabg data

        if exp_analysis_type == 'expression':
            avg_values.sort() ### Sort to get the lowest dabg and largest average expression
            if filetype == 'dabg':
                if avg_values[0]<=dabg_p_threshold: dabg_summary[probeset]=[] ### store probeset if the minimum p<user-threshold
            else:
                #if 'ENSMUSG00000018263:' in probeset: print probeset,[avg_values[-1],expression_threshold]
                if avg_values[-1]>=expression_threshold:
                    expression_summary[probeset]=[] ### store probeset if the minimum p<user-threshold
                    
def logTransform(exp_values):
    ### This code was added in version 1.16 in conjunction with a switch from logstatus to
    ### non-log in AltAnalyze to prevent "Process AltAnalyze Filtered" associated errors
    exp_values_log2=[]
    
    for exp_val in exp_values:
        exp_values_log2.append(str(math.log(float(exp_val),2))) ### changed from - log_fold = math.log((float(exp_val)+1),2) - version 2.05
    return exp_values_log2
              
def generateConstitutiveExpression(exp_dbase,constitutive_gene_db,probeset_gene_db,pre_filtered_db,array_names,filename):
    """Generate Steady-State expression values for each gene for analysis in the main module of this package"""
    steady_state_db={}; k=0; l=0
    remove_nonexpressed_genes = 'no' ### By default set to 'no'

    ###1st Pass: Identify probesets for steady-state calculation
    for gene in probeset_gene_db:
        if avg_all_probes_for_steady_state == 'yes': average_all_probesets[gene] = probeset_gene_db[gene] ### These are all exon aligning (not intron) probesets
        else:
            if gene not in constitutive_gene_db: average_all_probesets[gene] = probeset_gene_db[gene]
            else:
                constitutive_probeset_list = constitutive_gene_db[gene]
                constitutive_filtered=[] ###Added this extra code to eliminate constitutive probesets not in exp_dbase (gene level filters are more efficient when dealing with this many probesets)
                for probeset in constitutive_probeset_list:
                    if probeset in probeset_gene_db[gene]: constitutive_filtered.append(probeset)
                if len(constitutive_filtered)>0: average_all_probesets[gene] = constitutive_filtered
                else: average_all_probesets[gene] = probeset_gene_db[gene]

    ###2nd Pass: Remove probesets that have no detected expression (keep all if none are expressed)
    if excludeLowExpressionExons:
        non_expressed_genes={} ### keep track of these for internal QC
        for gene in average_all_probesets:
            gene_probe_list=[]; x = 0
            for probeset in average_all_probesets[gene]:
                if probeset in pre_filtered_db: gene_probe_list.append(probeset); x += 1
            ###If no constitutive and there are probes with detected expression: replace entry
            if x >0: average_all_probesets[gene] = gene_probe_list
            elif remove_nonexpressed_genes == 'yes': non_expressed_genes[gene]=[]   

    if remove_nonexpressed_genes == 'yes':
        for gene in non_expressed_genes: del average_all_probesets[gene]
    ###3rd Pass: Make sure the probesets are present in the input set (this is not typical unless a user is loading a pre-filtered probeset expression dataset)
    for gene in average_all_probesets:
        v=0
        for probeset in average_all_probesets[gene]:
            try: null = exp_dbase[probeset]; v+=1
            except KeyError: null =[] ###occurs if the expression probeset list is missing some of these probesets
            if v==0: ###Therefore, no probesets were found that were previously predicted to be best constitutive
                try: average_all_probesets[gene] = probeset_gene_db[gene] ###expand the average_all_probesets to include any exon linked to the gene
                except KeyError: print gene, probeset, len(probeset_gene_db), len(average_all_probesets);kill
    
    for probeset in exp_dbase:
        array_count = len(exp_dbase[probeset]); break

    try: null = array_count
    except Exception:
        print 'WARNING...CRITICAL ERROR. Make sure the correct array type is selected and that all input expression files are indeed present (array_count ERROR).'; forceError
        
    ###Calculate avg expression for each array for each probeset (using constitutive values)
    gene_count_db={}
    for gene in average_all_probesets:
        x = 0 ###For each array, average all probeset expression values
        gene_sum=0
        probeset_list = average_all_probesets[gene]#; k+= len(average_all_probesets[gene])
        if array_type != 'RNASeq': ### Just retain the list of probesets for RNA-seq
            while x < array_count:
                exp_list=[] ### average all exp values for constituitive probesets for each array
                for probeset in probeset_list:
                    try:
                        exp_val = exp_dbase[probeset][x]
                        exp_list.append(exp_val)
                    except KeyError: null =[] ###occurs if the expression probeset list is missing some of these probesets
                try:
                    if len(exp_list)==0:                
                        for probeset in probeset_list:
                            try:
                                exp_val = exp_dbase[probeset][x]
                                exp_list.append(exp_val)
                            except KeyError: null =[] ###occurs if the expression probeset list is missing some of these probesets
                    avg_const_exp=statistics.avg(exp_list)
                    ### Add only one avg-expression value for each array, this loop
                    try: steady_state_db[gene].append(avg_const_exp)
                    except KeyError: steady_state_db[gene] = [avg_const_exp]
                except ZeroDivisionError: null=[] ### Occurs when processing a truncated dataset (for testing usually) - no values for the gene should be included
                x += 1

    l = len(probeset_gene_db) - len(steady_state_db)
    steady_state_export = filename[0:-4]+'-steady-state.txt'
    steady_state_export = string.replace(steady_state_export,'counts.','exp.')
    fn=filepath(steady_state_export); data = open(fn,'w'); title = 'Gene_ID'
    
    if array_type == 'RNASeq':
        import RNASeq
        steady_state_db, pre_filtered_db = RNASeq.calculateGeneLevelStatistics(steady_state_export,species,average_all_probesets,normalize_feature_exp,array_names,UserOptions,excludeLowExp=excludeLowExpressionExons)
        ### This "pre_filtered_db" replaces the above since the RNASeq module performs the exon and junction-level filtering, not ExonArray (RPKM and count based)
        ### Use pre_filtered_db to exclude non-expressed features for multi-group alternative exon analysis
        removeNonExpressedProbesets(pre_filtered_db,full_dataset_export_dir)
        reload(RNASeq)
    
    for array in array_names: title = title +'\t'+ array
    data.write(title+'\n')
    for gene in steady_state_db:
        ss_vals = gene
        for exp_val in steady_state_db[gene]:
            ss_vals = ss_vals +'\t'+ str(exp_val)
        data.write(ss_vals+'\n')
    data.close()
    exp_dbase={}; steady_state_db={}; pre_filtered_db ={}
    #print k, "probesets were not found in the expression file, that could be used for the constitutive expression calculation"
    #print l, "genes were also not included that did not have such expression data"
    print "Steady-state data exported to",steady_state_export
    
def permformFtests(filtered_exp_db,group_count,probeset_db):
    ###Perform an f-test analysis to filter out low significance probesets
    ftest_gene_db={}; filtered_gene_db={}; filtered_gene_db2={}
    for probeset in filtered_exp_db:
        len_p = 0; ftest_list=[]
        try: gene_id = probeset_db[probeset][0]
        except KeyError: continue
        for len_s in group_count:
            index = len_s + len_p
            exp_group = filtered_exp_db[probeset][len_p:index]
            ftest_list.append(exp_group); len_p = index
        fstat,df1,df2 = statistics.Ftest(ftest_list)
        if fstat > f_cutoff:
            ftest_gene_db[gene_id] = 1
        try: filtered_gene_db[gene_id].append(probeset)
        except KeyError: filtered_gene_db[gene_id] = [probeset]
    ###Remove genes with no significant f-test probesets
    #print "len(ftest_gene_db)",len(filtered_gene_db)
    for gene_id in filtered_gene_db:
        if gene_id in ftest_gene_db:
            filtered_gene_db2[gene_id] = filtered_gene_db[gene_id]
    #print "len(filtered_gene_db)",len(filtered_gene_db2)
    return filtered_gene_db2

################# Resolve data annotations and filters

def eliminateRedundant(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def filtereLists(list1,db1):
    ###Used for large lists where string searches are computationally intensive
    list2 = []; db2=[]
    for key in db1:
        for entry in db1[key]:
            id = entry[1][-1]; list2.append(id)
    for key in list1: db2.append(key)
    for entry in list2: db2.append(entry)
    temp={}; combined = []
    for entry in db2:
        try:temp[entry] += 1
        except KeyError:temp[entry] = 1
    for entry in temp:
        if temp[entry] > 1:combined.append(entry)
    return combined

def makeGeneLevelAnnotations(probeset_db):
    transcluster_db = {}; exon_db = {}; probeset_gene_db={}
    #probeset_db[probeset] = gene,transcluster,exon_id,ens_exon_ids,exon_annotations,constitutitive
    ### Previously built the list of probesets for calculating steady-stae
    for gene in average_all_probesets:
        for probeset in average_all_probesets[gene]:
            gene_id,transcluster,exon_id,ens_exon_ids,affy_class = probeset_db[probeset]
            try: transcluster_db[gene].append(transcluster)
            except KeyError: transcluster_db[gene] = [transcluster]
            ens_exon_list = string.split(ens_exon_ids,'|')
            for exon in ens_exon_list:
                try: exon_db[gene].append(ens_exon_ids)
                except KeyError: exon_db[gene] = [ens_exon_ids]
    transcluster_db = eliminateRedundant(transcluster_db)
    exon_db = eliminateRedundant(exon_db)

    for gene in average_all_probesets:
        transcript_cluster_ids = string.join(transcluster_db[gene],'|')
        probeset_ids = string.join(average_all_probesets[gene],'|')
        exon_ids = string.join(exon_db[gene],'|')
        probeset_gene_db[gene] = transcript_cluster_ids,exon_ids,probeset_ids
    return probeset_gene_db

def getAnnotations(fl,Array_type,p_threshold,e_threshold,data_source,manufacturer,constitutive_source,Species,avg_all_for_ss,filter_by_DABG,perform_alt_analysis,expression_data_format):
    global species; species = Species; global average_all_probesets; average_all_probesets={}
    global avg_all_probes_for_steady_state; avg_all_probes_for_steady_state = avg_all_for_ss; global filter_by_dabg; filter_by_dabg = filter_by_DABG
    global dabg_p_threshold; dabg_p_threshold = float(p_threshold); global root_dir; global biotypes; global normalize_feature_exp
    global expression_threshold; global exp_data_format; exp_data_format = expression_data_format; global UserOptions; UserOptions = fl
    global full_dataset_export_dir; global excludeLowExpressionExons

    """
    try: exon_exp_threshold = fl.ExonExpThreshold()
    except Exception: exon_exp_threshold = 0
    try: exon_rpkm_threshold = fl.ExonRPKMThreshold()
    except Exception: exon_rpkm_threshold = 0
    try: gene_rpkm_threshold = fl.RPKMThreshold()
    except Exception: gene_rpkm_threshold = 0
    try: gene_exp_threshold = fl.GeneExpThreshold()
    except Exception: gene_exp_threshold = 0
    """
    
    ### The input expression data can be log or non-log. If non-log, transform to log in FilterDABG prior to the alternative exon analysis - v.1.16    
    if expression_data_format == 'log':
        try: expression_threshold = math.log(float(e_threshold),2)
        except Exception: expression_threshold = 0 ### Applies to RNASeq datasets
    else:
        expression_threshold = float(e_threshold)
    
    process_from_scratch = 'no' ###internal variables used while testing
    global dabg_summary; global expression_summary; dabg_summary={};expression_summary={}
    global fulldataset_export_object; global array_type; array_type = Array_type
    global exp_analysis_type; exp_analysis_type = 'expression'
    global stats_input_dir
    expr_input_dir = fl.ExpFile(); stats_input_dir = fl.StatsFile(); root_dir = fl.RootDir()
    try: normalize_feature_exp = fl.FeatureNormalization()
    except Exception: normalize_feature_exp = 'NA'
    try: excludeLowExpressionExons = fl.excludeLowExpressionExons()
    except Exception: excludeLowExpressionExons = True
    try:
        useJunctionsForGeneExpression = fl.useJunctionsForGeneExpression()
        if useJunctionsForGeneExpression:
            print 'Using known junction only to estimate gene expression!!!'
    except Exception: useJunctionsForGeneExpression = False
    
    source_biotype = 'mRNA'
    if array_type == 'gene': source_biotype = 'gene'
    elif array_type == 'junction': source_biotype = 'junction'
    ###Get annotations using Affymetrix as a trusted source or via links to Ensembl

    if array_type == 'AltMouse':
        probeset_db,constitutive_gene_db = ExpressionBuilder.importAltMerge('full'); annotate_db={}
        source_biotype = 'AltMouse'
    elif manufacturer == 'Affymetrix' or array_type == 'RNASeq':
        if array_type == 'RNASeq':
            source_biotype = array_type, root_dir

        probeset_db,annotate_db,constitutive_gene_db,splicing_analysis_db = ExonArrayEnsemblRules.getAnnotations(process_from_scratch,constitutive_source,source_biotype,species)

    ### Get all file locations and get array headers
    #print len(splicing_analysis_db),"genes included in the splicing annotation database (constitutive only containing)"
    stats_file_status = verifyFile(stats_input_dir)
    array_linker_db,array_names = importExonProbesetData(expr_input_dir,{},'arraynames')
    input_dir_split = string.split(expr_input_dir,'/')
    full_dataset_export_dir = root_dir+'AltExpression/FullDatasets/ExonArray/'+species+'/'+string.replace(input_dir_split[-1],'exp.','')
    if array_type == 'gene': full_dataset_export_dir = string.replace(full_dataset_export_dir,'ExonArray','GeneArray')
    if array_type == 'junction': full_dataset_export_dir = string.replace(full_dataset_export_dir,'ExonArray','JunctionArray')
    if array_type == 'AltMouse': full_dataset_export_dir = string.replace(full_dataset_export_dir,'ExonArray','AltMouse')
    if array_type == 'RNASeq': full_dataset_export_dir = string.replace(full_dataset_export_dir,'ExonArray','RNASeq')
    try: fulldataset_export_object = export.ExportFile(full_dataset_export_dir)
    except Exception:
        print 'AltAnalyze is having trouble creating the directory:\n',full_dataset_export_dir
        print 'Report this issue to the AltAnalyze help desk or create this directory manually (Error Code X1).'; force_exception
    ### Organize arrays according to groups and export all probeset data and any pairwise comparisons
    data_type = 'expression'
    if array_type == 'RNASeq':
        expr_input_dir = string.replace(expr_input_dir,'exp.','counts.') ### Filter based on the counts file and then replace values with the normalized as the last step
    comparison_filename_list,biotypes = exportGroupedComparisonProbesetData(expr_input_dir,probeset_db,data_type,array_names,array_linker_db,perform_alt_analysis)
    if useJunctionsForGeneExpression:
        if 'junction' in biotypes:
            if 'exon' in biotypes: del biotypes['exon']
    if filter_by_dabg == 'yes' and stats_file_status == 'found':
        data_type = 'dabg'
        exportGroupedComparisonProbesetData(stats_input_dir,probeset_db,data_type,array_names,array_linker_db,perform_alt_analysis)
    ###Filter expression data based on DABG and annotation filtered probesets (will work without DABG filtering as well) - won't work for RNA-Seq (execute function later)
    filtered_exon_db = removeNonExpressedProbesets(probeset_db,full_dataset_export_dir)
    filterExpressionData(expr_input_dir,filtered_exon_db,constitutive_gene_db,probeset_db,'expression',array_names,perform_alt_analysis)
    constitutive_gene_db={}; probeset_gene_db = makeGeneLevelAnnotations(probeset_db)

    if array_type == 'RNASeq':
        fulldataset_export_object = export.ExportFile(full_dataset_export_dir)
        data_type = 'expression' ### Repeat with counts and then with exp. to add gene-level estimates to both
        exportGroupedComparisonProbesetData(expr_input_dir,probeset_db,data_type,array_names,array_linker_db,perform_alt_analysis)
        fulldataset_export_object = export.ExportFile(full_dataset_export_dir)
        expr_input_dir = string.replace(expr_input_dir,'counts.','exp.')
        exportGroupedComparisonProbesetData(expr_input_dir,probeset_db,data_type,array_names,array_linker_db,perform_alt_analysis)
        
    try: clearObjectsFromMemory(average_all_probesets); clearObjectsFromMemory(expression_summary); clearObjectsFromMemory(splicing_analysis_db)
    except Exception: null=[]
    filtered_exon_db=[]; probeset_db={}; average_all_probesets={}; expression_summary={}; splicing_analysis_db={}
    #filtered_exp_db,group_count,ranked_array_headers = filterExpressionData(expr_input_dir,filtered_exon_db,constitutive_gene_db,probeset_db)
    #filtered_gene_db = permformFtests(filtered_exp_db,group_count,probeset_db)

    """    
    pre_filtered_db=[]
    print 'global vars'
    returnLargeGlobalVars()    
    print 'local vars'
    all = [var for var in locals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all:
            try:
                if len(locals()[var])>500: print var, len(locals()[var])
            except Exception: null=[]
    """
    return probeset_gene_db, annotate_db, comparison_filename_list

def processResiduals(fl,Array_type,Species,perform_alt_analysis):
    global species; species = Species; global root_dir; global fulldataset_export_object
    global array_type; array_type = Array_type; global exp_analysis_type; exp_analysis_type = 'residual'
    
    ### Get all file locations and get array headers
    expr_input_dir = fl.ExpFile(); root_dir = fl.RootDir()
    array_linker_db,array_names = importExonProbesetData(expr_input_dir,{},'arraynames')
    input_dir_split = string.split(expr_input_dir,'/')
    full_dataset_export_dir = root_dir+'AltExpression/FIRMA/FullDatasets/'+array_type+'/'+species+'/'+string.replace(input_dir_split[-1],'exp.','')
    expr_input_dir = string.replace(expr_input_dir,'exp.','residuals.') ### Wait to change this untile the above processes are finished
    
    fulldataset_export_object = export.ExportFile(full_dataset_export_dir)
    #print full_dataset_export_dir
    ### Organize arrays according to groups and export all probeset data and any pairwise comparisons
    comparison_filename_list = exportGroupedComparisonProbesetData(expr_input_dir,{},'residuals',array_names,array_linker_db,perform_alt_analysis)

def removeNonExpressedProbesets(probeset_db,full_dataset_export_dir):
    combined_db={}
    if array_type == 'RNASeq':
        id_name = 'junction IDs'
        combined_db = probeset_db
        print len(combined_db), id_name,'after detection RPKM and read count filtering.'
    else:
        id_name = 'array IDs'
        print len(expression_summary), 'expression and',len(dabg_summary),'detection p-value filtered '+id_name+' out of', len(probeset_db)
        for probeset in probeset_db:
            if len(dabg_summary)>0:
                try:
                    n = expression_summary[probeset]
                    s = dabg_summary[probeset]
                    combined_db[probeset]=[]
                except Exception: null=[]
            else:
                try:
                    n = expression_summary[probeset]
                    combined_db[probeset]=[]
                except Exception: null=[]
        print len(combined_db), id_name,'after detection p-value and expression filtering.'
    rewriteOrganizedExpressionFile(combined_db,full_dataset_export_dir)
    return combined_db

def rewriteOrganizedExpressionFile(combined_db,full_dataset_export_dir):
    importExonProbesetData(full_dataset_export_dir,combined_db,'filterDataset')
    temp_dir = string.replace(full_dataset_export_dir,'.txt','-temp.txt')
    import shutil
    try:
        shutil.copyfile(temp_dir, full_dataset_export_dir) ### replace unfiltered file
        os.remove(temp_dir)
    except Exception: null=[]
    
def inputResultFiles(filename,file_type):
    fn=filepath(filename)
    gene_db={}; x = 0
    ###Import expression data (non-log space)
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            if file_type == 'gene': geneid = t[0]; gene_db[geneid]=[]
            else: probeset = t[7]; gene_db[probeset]=[]
    return gene_db
                
def grabExonIntronPromoterSequences(species,array_type,data_type,output_types):
    ### output_types could be adjacent intron sequences, adjacent exon sequences, targets exon sequence or promoter
    sequence_input_dir_list=[]
    if data_type == 'probeset': sequence_input_dir = '/AltResults/AlternativeOutput/'+array_type+'/sequence_input'
    if data_type == 'gene': sequence_input_dir = '/ExpressionOutput/'+array_type+'/sequence_input'
    
    dir_list = read_directory(sequence_input_dir)
    for input_file in dir_list:
        filedir = sequence_input_dir[1:]+'/'+input_file
        filter_db = inputResultFiles(filedir,data_type)
        export_exon_filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_probesets.txt'        
        ensembl_probeset_db = ExonArrayEnsemblRules.reimportEnsemblProbesetsForSeqExtraction(export_exon_filename,data_type,filter_db)
        """for gene in ensembl_probeset_db:
            if gene == 'ENSG00000139737':
                for x in ensembl_probeset_db[gene]:
                    exon_id,((probe_start,probe_stop,probeset_id,exon_class,transcript_clust),ed) = x
                    print gene, ed.ExonID()
        kill"""
        analysis_type = 'get_sequence'
        dir = 'AltDatabase/ensembl/'+species+'/'; gene_seq_filename = dir+species+'_gene-seq-2000_flank'
        ensembl_probeset_db = EnsemblImport.import_sequence_data(gene_seq_filename,ensembl_probeset_db,species,analysis_type)

        """
        critical_exon_file = 'AltDatabase/'+species+'/'+ array_type + '/' + array_type+'_critical-exon-seq.txt'
        if output_types == 'all' and data_type == 'probeset':
            output_types = ['alt-promoter','promoter','exon','adjacent-exons','adjacent-introns']
        else: output_types = [output_types]
        
        for output_type in output_types:
            sequence_input_dir = string.replace(sequence_input_dir,'_input','_output')
            filename = sequence_input_dir[1:]+'/ExportedSequence-'+data_type+'-'+output_type+'.txt'
            exportExonIntronPromoterSequences(filename, ensembl_probeset_db,data_type,output_type)
        """
        if output_types == 'all' and data_type == 'probeset':
            output_types = ['alt-promoter','promoter','exon','adjacent-exons','adjacent-introns']
        else: output_types = [output_types]
        
        for output_type in output_types:
            sequence_input_dir2 = string.replace(sequence_input_dir,'_input','_output')
            filename = sequence_input_dir2[1:]+'/'+input_file[:-4]+'-'+data_type+'-'+output_type+'.txt'
            exportExonIntronPromoterSequences(filename, ensembl_probeset_db,data_type,output_type)

def exportExonIntronPromoterSequences(filename,ensembl_probeset_db,data_type,output_type):
    exon_seq_db_filename = filename
    fn=filepath(exon_seq_db_filename); data = open(fn,'w'); gene_data_exported={}; probe_data_exported={}

    if data_type == 'gene' or output_type == 'promoter': seq_title = 'PromoterSeq'
    elif output_type == 'alt-promoter': seq_title = 'PromoterSeq'
    elif output_type == 'exon': seq_title = 'ExonSeq'
    elif output_type == 'adjacent-exons': seq_title = 'PrevExonSeq\tNextExonSeq'
    elif output_type == 'adjacent-introns': seq_title = 'PrevIntronSeq\tNextIntronSeq'

    title = ['Ensembl_GeneID','ExonID','ExternalExonIDs','ProbesetID',seq_title]
    title = string.join(title,'\t')+'\n'; #data.write(title)
    for ens_gene in ensembl_probeset_db:
        for probe_data in ensembl_probeset_db[ens_gene]:
            exon_id,((probe_start,probe_stop,probeset_id,exon_class,transcript_clust),ed) = probe_data
            ens_exon_list = ed.ExonID(); ens_exons = string.join(ens_exon_list,' ')
            proceed = 'no'
            try:
                if data_type == 'gene' or output_type == 'promoter':
                    try: seq = [ed.PromoterSeq()]; proceed = 'yes'
                    except AttributeError: proceed = 'no'
                elif output_type == 'alt-promoter':
                    if 'alt-N-term' in ed.AssociatedSplicingEvent() or 'altPromoter' in ed.AssociatedSplicingEvent():
                        try: seq = [ed.PrevIntronSeq()]; proceed = 'yes'
                        except AttributeError: proceed = 'no'
                elif output_type == 'exon':
                    try: seq = [ed.ExonSeq()]; proceed = 'yes'
                    except AttributeError: proceed = 'no'
                elif output_type == 'adjacent-exons':
                    if len(ed.PrevExonSeq())>1 and len(ed.NextExonSeq())>1:
                        try: seq = ['(prior-exon-seq)'+ed.PrevExonSeq(),'(next-exon-seq)'+ed.NextExonSeq()]; proceed = 'yes'
                        except AttributeError: proceed = 'no'
                elif output_type == 'adjacent-introns':
                    if len(ed.PrevIntronSeq())>1 and len(ed.NextIntronSeq())>1:
                        try: seq = ['(prior-intron-seq)'+ed.PrevIntronSeq(), '(next-intron-seq)'+ed.NextIntronSeq()]; proceed = 'yes'
                        except AttributeError: proceed = 'no'
            except AttributeError: proceed = 'no'
            if proceed == 'yes':
                gene_data_exported[ens_gene]=[]
                probe_data_exported[probeset_id]=[]
                if data_type == 'gene':
                    values = ['>'+ens_gene]
                else: values = ['>'+ens_gene,exon_id,ens_exons,probeset_id]
                values = string.join(values,'|')+'\n';data.write(values)
                for seq_data in seq:
                    i = 0; e = 100
                    while i < len(seq_data):
                        seq_line = seq_data[i:e]; data.write(seq_line+'\n')
                        i+=100; e+=100

    #print len(gene_data_exported), 'gene entries exported'
    #print len(probe_data_exported), 'probeset entries exported'
    data.close()
    #print exon_seq_db_filename, 'exported....'

def verifyFile(filename):
    status = 'not found'
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = 'found';break
    except Exception: status = 'not found'
    return status

if __name__ == '__main__':
    m = 'Mm'
    h = 'Hs'
    r = 'Rn'
    Species = h
    Array_type = 'junction'
    exportMetaProbesets(Array_type,Species);sys.exit()
    Data_type = 'probeset'
    Output_types = 'promoter'
    Output_types = 'all'
    
    grabExonIntronPromoterSequences(Species,Array_type,Data_type,Output_types)
    sys.exit()
    #"""
    avg_all_for_ss = 'yes'
    import_dir = '/AltDatabase/'+Species+ '/exon'
    expr_file_dir = 'ExpressionInput\exp.HEK-confluency.plier.txt'
    dagb_p = 0.001
    f_cutoff = 2.297
    exons_to_grab = "core"
    x = 'Affymetrix'
    y = 'Ensembl'
    z = 'default'
    data_source = y
    constitutive_source = z
    filename = expr_file_dir; p = dagb_p
    getAnnotations(expr_file_dir,dagb_p,exons_to_grab,data_source,constitutive_source,Species)
    global species; species = Species
    process_from_scratch = 'no'
    ###Get annotations using Affymetrix as a trusted source or via links to Ensembl
    if data_source == 'Affymetrix':
        annotation_dbases = ExonArrayAffyRules.getAnnotations(exons_to_grab,constitutive_source,process_from_scratch)
        probe_association_db,constitutive_gene_db,exon_location_db, trans_annotation_db, trans_annot_extended = annotation_dbases
    else:
        probeset_db,annotate_db,constitutive_gene_db,splicing_analysis_db = ExonArrayEnsemblRules.getAnnotations(process_from_scratch,constitutive_source,species,avg_all_for_ss)

    filterExpressionData(filename,filtered_exon_db,constitutive_gene_db,probeset_db,data_type)
    #filtered_gene_db = permformFtests(filtered_exp_db,group_count,probeset_db)
    
    
