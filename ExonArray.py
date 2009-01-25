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

import sys, string
import os.path
import unique
import statistics
import math
import EnsemblImport; reload(EnsemblImport)
import ExonArrayEnsemblRules; reload(ExonArrayEnsemblRules)
import ExonArrayAffyRules
import ExpressionBuilder
import reorder_arrays
import time
import export

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

def getFilteredExons(filename,p,probeset_db):
    start_time = time.time()
    ###Get statistics file based on expression file's name
    #Re-write to take into account groups
    #filename= string.replace(filename,'plier','dabg')
    #filename = string.replace(filename,'exp.','stats.')
    fn=filepath(filename)
    filtered_exon_list={}; total = 0; x=0
    try:
        for line in open(fn,'rU').xreadlines():
          data = cleanUpLine(line)
          if data[0] != '#' and x == 1:
            tab_delimited_data = string.split(data,'\t')
            probeset = tab_delimited_data[0]
            try:
                null = probeset_db[probeset] ### Verify the probeset linked to gene annotations
                stats = tab_delimited_data[1:]
                number_significant = 0
                for statistic in stats:
                    try:
                        if float(statistic) < p: number_significant += 1
                    except ValueError: print line,'\n',data,'\n',tab_delimited_data;kill
                if number_significant > 0: ###Make this as simple as possible... if at least one value has a p<0.05, accept it
                    filtered_exon_list[probeset]=[] #float(number_significant)/len(stats)
                total += 1
            except KeyError: null = []
          if data[0] != '#' and x == 0: x=1
        end_time = time.time(); time_diff = int(end_time-start_time)
        print len(filtered_exon_list),"probesets remaining after dabg filtering, out of",total,'(data analyzed in %d seconds)' % time_diff  
    except IOError:
        filtered_exon_list = filtered_exon_list
    return filtered_exon_list

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def importExonProbesetData(filename,import_these_probesets,import_type):
    fn=filepath(filename); start_time = time.time()
    exp_dbase={}; filtered_exp_db={}; ftest_gene_db={}; filtered_gene_db={}; probeset_gene_db={}; d = 0; x = 0
    ###Import expression data (non-log space)
    try:
        for line in open(fn,'rU').xreadlines():             
          data = cleanUpLine(line)
          if len(data)==0: null = ''
          elif data[0] != '#' and x == 1:   ###Grab expression values
            tab_delimited_data = string.split(data,'\t')
            probeset = tab_delimited_data[0]
            if import_type == 'probesets': ###Occurs on the first round
                exp_vals = []; exp_dbase[probeset] = exp_vals
            elif import_type == 'raw':
                try:
                    null = import_these_probesets[probeset]; exp_vals = tab_delimited_data[1:]; exp_dbase[probeset] = exp_vals
                except KeyError: null = [] ###Don't import any probeset data
            elif import_type == 'comparisons':
                try:
                    null = import_these_probesets[probeset]; exp_vals = tab_delimited_data[1:]
                    filtered_exp_db={}; filtered_exp_db[probeset] = exp_vals
                    reorderArraysOnly(filtered_exp_db) ###order and directly write data
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
    except IOError: exp_dbase = exp_dbase

    end_time = time.time(); time_diff = int(end_time-start_time)
    print "Exon probeset data imported in %d seconds" % time_diff
    
    if import_type == 'raw': print len(exp_dbase),"probesets imported with expression values"
    try:
        if import_type != 'comparisons': return exp_dbase,array_linker_db,array_names
    except UnboundLocalError: ### Occurs if no p-value file exists
        return 'null', 'null','null'
    
def filterExpressionData(filename,pre_filtered_db,constitutive_gene_db,probeset_db,data_type,perform_alt_analysis):
    """Probeset level data and gene level data are handled differently by this program for exon and tiling based arrays.
    Probeset level data is sequentially filtered to reduce the dataset to a minimal set of expressed probesets (exp p-values) that
    that align to genes with multiple lines of evidence (e.g. ensembl) and show some evidence of regulation for any probesets in a
    gene (transcript_cluster). Gene level analyses are handled under the main module of the program, exactly like 3' arrays, while
    probe level data is reorganized, filtered and output from this module."""

    comparison_filename_list=[]
    ###First time we import, just grab the probesets associated
    exp_dbase,array_linker_db,array_names = importExonProbesetData(filename,{},'probesets')
    if array_names != 'null': ### Then there is only an expr. file and not a stats file
        constitutive_probeset_db = {}
        for gene in constitutive_gene_db:
            for probeset in constitutive_gene_db[gene]:
                constitutive_probeset_db[probeset] = gene

        ###Identify constitutive probesets to import (sometimes not all probesets on the array are imported)
        possible_genes={}; possible_constitutive_probeset={}; genes_with_constitutive_probesets={}; probeset_gene_db={}
        for probeset in exp_dbase:
            try:
                gene = probeset_db[probeset][0]; exon_id = probeset_db[probeset][2]
                try: probeset_gene_db[gene].append(probeset)
                except KeyError: probeset_gene_db[gene] = [probeset]
                if gene not in constitutive_gene_db: ###Thus there no are constitutive probesets for the gene
                    if 'E' in exon_id:
                        possible_constitutive_probeset[probeset] = []; genes_with_constitutive_probesets[gene] = []
                elif probeset in constitutive_probeset_db: possible_constitutive_probeset[probeset] = []; genes_with_constitutive_probesets[gene] = []
                possible_genes[gene] = []
            except KeyError: null = []
        if len(possible_genes) != genes_with_constitutive_probesets: ###Thus all constitutive probesets for a gene are missing from the input file.  Need to use all other probesets for the gene
            for gene in genes_with_constitutive_probesets: del possible_genes[gene] ###only leave unique unaccounted for genes (probably more efficient that a On2 lookup
            for probeset in exp_dbase:
                try:
                    gene = probeset_db[probeset][0]
                    if gene in possible_genes: possible_constitutive_probeset[probeset] = []
                except KeyError: null = [] ###gene not in annotation dbase

        ### Only import probesets that can be used to calculate gene expression values OR link to gene annotations (which have at least one dabg p<0.05 for all samples normalized)
        import_these_probesets={}
        for probeset in possible_constitutive_probeset: import_these_probesets[probeset] = []
        #for probeset in pre_filtered_db: import_these_probesets[probeset] = []
        
        print "Processing Expression Dataset..."
        constitutive_exp_dbase,array_linker_db,array_names = importExonProbesetData(filename,import_these_probesets,'raw') 
        """constitutive_exp_dbase={} ### Create a new dictionary containing just constitutive probesets
        for probeset in possible_constitutive_probeset:
            constitutive_exp_dbase[probeset] = filtered_exp_db[probeset]"""

        generateConstitutiveExpression(constitutive_exp_dbase,constitutive_gene_db,probeset_gene_db,pre_filtered_db,array_names,filename)
        probeset_db = {}; constitutive_exp_dbase = {}; possible_constitutive_probeset={} ### reset databases to conserve memory

        if perform_alt_analysis != 'expression': ### User Option
            exp_dbase={}; constitutive_gene_db={}; probeset_gene_db={} ### reset databases to conserve memory
            global expr_group_list; global comp_group_list
            if data_type == 'expression':
                expr_group_dir = string.replace(filename,'exp.','groups.')
                comp_group_dir = string.replace(filename,'exp.','comps.')
            else:
                expr_group_dir = string.replace(filename,'stats.','groups.')
                comp_group_dir = string.replace(filename,'stats.','comps.')
            comp_group_list, comp_group_list2 = ExpressionBuilder.importComparisonGroups(comp_group_dir)
            expr_group_list,expr_group_db = ExpressionBuilder.importArrayGroups(expr_group_dir,array_linker_db)

            print "Reorganizing array data into comparison groups for export to down-stream splicing analysis software"
            ###Do this only for the header data
            group_count,raw_data_comp_headers = reorder_arrays.reorderArrayHeaders(array_names,expr_group_list,comp_group_list,array_linker_db)

            ###Export the header info and store the export write data for reorder_arrays
            global comparision_export_db; comparision_export_db={}
            AltAnalzye_input_dir = "AltExpression/pre-filtered/"+data_type+'/'; array_type_name = 'Exon'
            for comparison in comp_group_list2: ###loop throught the list of comparisons
                group1 = comparison[0]; group2 = comparison[1]
                group1_name = expr_group_db[group1]; group2_name = expr_group_db[group2]
                comparison_filename = species+'_'+array_type_name+'_'+ group1_name + '_vs_' + group2_name + '.txt'
                
                new_file = AltAnalzye_input_dir + comparison_filename; comparison_filename_list.append(comparison_filename)
                data = export.createExportFile(new_file,AltAnalzye_input_dir[:-1])

                try: array_names = raw_data_comp_headers[comparison]
                except KeyError: print raw_data_comp_headers;kill
                title = ['Probesets']+array_names; title = string.join(title,'\t')+'\n'; data.write(title)
                comparision_export_db[comparison] = data ###store the export file write data so we can write after organizing

            #reorder_arrays.reorderArraysOnly(filtered_exp_db,expr_group_list,comp_group_list)
            print len(pre_filtered_db), 'pre_filtered_db probesets being selected from expression file'
            importExonProbesetData(filename,pre_filtered_db,'comparisons')
            pre_filtered_db={}
            
            for comparison in comparision_export_db:
                data = comparision_export_db[comparison]
                data.close()
            print "Pairwise comparisons for AltAnalyze exported..."
            
            #export_summary_stats='no'; array_type = "exon-array"
            #ExpressionBuilder.exportSplicingInput(species,array_type,expr_group_db,raw_data_comp_headers,comp_group_list2,raw_data_comps,export_summary_stats,data_type)
            """
            null,filename = string.split(filename,'\\')
            filtered_exp_export = 'AltExpression/pre-filtered/'+species+'_Exon-'+filename[4:]
            fn=filepath(filtered_exp_export); data = open(fn,'w'); title = 'Probeset'
            for array in ranked_array_headers: title = title +'\t'+ array
            data.write(title+'\n')
            for probeset in filtered_exp_db:
                exp_vals = probeset
                for exp_val in filtered_exp_db[probeset]:
                    exp_vals = exp_vals +'\t'+ str(exp_val)
                data.write(exp_vals+'\n')
            data.close()
            #return filtered_exp_db,group_count,ranked_array_headers"""
    return comparison_filename_list

def reorderArraysOnly(filtered_exp_db): 
    ###array_order gives the final level order sorted, followed by the original index order as a tuple                   
    ###expr_group_list gives the final level order sorted, followed by the original index order as a tuple                   
    for probeset in filtered_exp_db:
        grouped_ordered_array_list = {}
        for x in expr_group_list:
            y = x[1]  ### this is the new first index
            group = x[2]   
            ### for example y = 5, therefore the filtered_exp_db[probeset][5] entry is now the first
            try:
                try: new_item = filtered_exp_db[probeset][y]
                except TypeError: print y,x,expr_group_list; kill
            except IndexError: print probeset,y,x,expr_group_list,'\n',filtered_exp_db[probeset];kill
            ###Used for comparision analysis
            try: grouped_ordered_array_list[group].append(new_item)
            except KeyError: grouped_ordered_array_list[group] = [new_item]
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

def generateConstitutiveExpression(exp_dbase,constitutive_gene_db,probeset_gene_db,pre_filtered_db,array_names,filename):
    """Generate Steady-State expression values for each gene for analysis in the main module of this package"""
    steady_state_db={}; k=0; l=0
    ###1st Pass: Identify probesets for steady-state calculation
    for gene in probeset_gene_db:
        if avg_all_probes_for_steady_state == 'yes': average_all_probesets[gene] = probeset_gene_db[gene]
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
    for gene in average_all_probesets:
        gene_probe_list=[]; x = 0
        if filter_by_dabg == 'yes':
            for probeset in average_all_probesets[gene]:
                if probeset in pre_filtered_db: gene_probe_list.append(probeset); x += 1
        else:
            for probeset in average_all_probesets[gene]: gene_probe_list.append(probeset); x += 1
        ###If no constitutive and there are probes with detected expression: replace entry
        if x >0: average_all_probesets[gene] = gene_probe_list

    ###3rd Pass: Make sure the probesets are present in the input set
    for gene in average_all_probesets:
        v=0
        for probeset in average_all_probesets[gene]:
            try: null = exp_dbase[probeset]; v+=1
            except KeyError: null =[] ###occurs if the expression probeset list is missing some of these probesets
            if v==0: ###Therefore, no probesets were found that were previously predicted to be best constitutive
                try: average_all_probesets[gene] = probeset_gene_db[gene] ###expand the average_all_probesets to include any exon linked to the gene
                except KeyError: print gene, probeset, len(probeset_gene_db), len(average_all_probesets);kill
                
    for probeset in exp_dbase: array_count = len(exp_dbase[probeset]); break
    ###Calculate avg expression for each array for each probeset (using constitutive values)
    for gene in average_all_probesets:
        x = 0 ###For each array, average all probeset expression values
        probeset_list = average_all_probesets[gene]#; k+= len(average_all_probesets[gene])
        while x < array_count:
            exp_list=[] ### average all exp values for constituitive probesets for each array
            for probeset in probeset_list:
                try: exp_val = exp_dbase[probeset][x]; exp_list.append(exp_val)
                except KeyError: null =[] ###occurs if the expression probeset list is missing some of these probesets
            avg_const_exp=statistics.avg(exp_list)
            ### Add only one avg-expression value for each array, this loop
            try: steady_state_db[gene].append(avg_const_exp)
            except KeyError: steady_state_db[gene] = [avg_const_exp]
            x += 1
            
    l = len(probeset_gene_db) - len(steady_state_db)
    steady_state_export = filename[0:-4]+'-steady-state.txt'
    fn=filepath(steady_state_export); data = open(fn,'w'); title = 'Gene_ID'
    for array in array_names: title = title +'\t'+ array
    data.write(title+'\n')
    for gene in steady_state_db:
        ss_vals = gene
        for exp_val in steady_state_db[gene]:
            ss_vals = ss_vals +'\t'+ str(exp_val)
        data.write(ss_vals+'\n')
    data.close()
    exp_dbase={}; steady_state_db={}
    print k, "probesets were not found in the expression file, that could be used for the constitutive expression calculation"
    print l, "genes were also not included that did not have such expression data"
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
    print "len(ftest_gene_db)",len(filtered_gene_db)
    for gene_id in filtered_gene_db:
        if gene_id in ftest_gene_db:
            filtered_gene_db2[gene_id] = filtered_gene_db[gene_id]
    print "len(filtered_gene_db)",len(filtered_gene_db2)
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

def filterFilteredData(db1,db2):
    """Filter out probesets that were filtered based on expression by whether those probesets are annotated as
    corresponding to exons/introns/utr regions of known gene"""
    combined_db={}
    if len(db2)<1:
       for probeset in db1: combined_db[probeset]=[]
    else:
        for probeset in db1:
                try:
                    null = db2[probeset] ###faster than find (or if statement)
                    combined_db[probeset]=[]
                except KeyError: null = []
    print 'Number of filtered probesets (expression p-value and choosen annotation method):', len(combined_db)
    return combined_db

def makeGeneLevelAnnotations(probeset_db):
    transcluster_db = {}; exon_db = {}; probeset_gene_db={}
    #probeset_db[probeset] = gene,transcluster,exon_id,ens_exon_ids,exon_annotations,constitutitive
    ### Previously built the list of probesets for calculating steady-stae
    for gene in average_all_probesets:
        for probeset in average_all_probesets[gene]:
            gene_id,transcluster,exon_id,ens_exon_ids = probeset_db[probeset]
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

def getAnnotations(fl,stats_input_dir,p,exons_to_grab,data_source,manufacturer,constitutive_source,Species,avg_all_for_ss,filter_by_DABG,perform_alt_analysis):
    global species; species = Species; global average_all_probesets; average_all_probesets={}; filtered_exon_list=[]
    global avg_all_probes_for_steady_state; avg_all_probes_for_steady_state = avg_all_for_ss; global filter_by_dabg; filter_by_dabg = filter_by_DABG
    process_from_scratch = 'no'; constitutive_analyzed = 'no'  ###internal variables used while testing
    source_biotype = 'mRNA'
    ###Get annotations using Affymetrix as a trusted source or via links to Ensembl
    if data_source == 'Affymetrix':
        annotation_dbases = ExonArrayAffyRules.getAnnotations(exons_to_grab,constitutive_source,process_from_scratch)
        probe_association_db,constitutive_gene_db,exon_location_db, trans_annotation_db, trans_annot_extended = annotation_dbases
    else:
        if manufacturer == 'Affymetrix': probeset_db,annotate_db,constitutive_gene_db,splicing_analysis_db = ExonArrayEnsemblRules.getAnnotations(process_from_scratch,constitutive_source,source_biotype,species)
        if constitutive_analyzed == 'yes':
            return probeset_gene_db,annotate_db
    if filter_by_dabg == 'yes':
        print len(splicing_analysis_db),"genes included in the splicing annotation database (constitutive only containing)"
        ###Step 1) generate a list of probesets with a DABG p<threshold
        expr_input_dir = fl.ExpFile(); stats_input_dir = fl.StatsFile()
        filtered_exon_list = getFilteredExons(stats_input_dir,p,probeset_db)
        #filtered_exon_list=[]
        ###Step 2) find the intersection of probesets from step 1 and probesets with valid gene annotations (for export to splicing-analyses)
        filtered_exon_list = filterFilteredData(probeset_db,filtered_exon_list)

    ###Filter expression data based on DABG and annotation filtered probesets (will work without DABG filtering as well)
    data_type = 'expression'
    comparison_filename_list = filterExpressionData(expr_input_dir,filtered_exon_list,constitutive_gene_db,probeset_db,data_type,perform_alt_analysis)
    data_type = 'dabg' ###repeat, but for dabg statistics rather than expression data
    comparison_filename_list2 = filterExpressionData(stats_input_dir,filtered_exon_list,constitutive_gene_db,probeset_db,data_type,perform_alt_analysis)
    constitutive_gene_db={}
    probeset_gene_db = makeGeneLevelAnnotations(probeset_db)
    
    filtered_exon_list=[]; probeset_db={}
    #filtered_exp_db,group_count,ranked_array_headers = filterExpressionData(expr_input_dir,filtered_exon_list,constitutive_gene_db,probeset_db)
    #filtered_gene_db = permformFtests(filtered_exp_db,group_count,probeset_db)
    return probeset_gene_db,annotate_db, comparison_filename_list

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

    print len(gene_data_exported), 'gene entries exported'
    print len(probe_data_exported), 'probeset entries exported'
    data.close()
    print exon_seq_db_filename, 'exported....'
        
if __name__ == '__main__':
    m = 'Mm'
    h = 'Hs'
    Species = h
    Array_type = 'exon'
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
    filtered_exon_list = {}
    filtered_exon_list = getFilteredExons(filename,p)
    filtered_exon_list = filterFilteredData(splicing_analysis_db,filtered_exon_list)
    #filtered_exon_list = filterFilteredData(splicing_analysis_db,filtered_exon_list)
    #"""
    filterExpressionData(filename,filtered_exon_list,constitutive_gene_db,probeset_db,data_type)
    #filtered_gene_db = permformFtests(filtered_exp_db,group_count,probeset_db)
    
    