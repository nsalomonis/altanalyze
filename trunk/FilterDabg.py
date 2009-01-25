###FilterDabg
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

import math
import statistics
import sys, string
import export
import os.path
import unique
import time

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

########### Begin Analyses ###########

def import_altmerge(filename,array_type):
    global exon_db
    fn=filepath(filename)
    altmerge_constituitive ={}; constituitive_gene={}; constitutive_original={}
    exon_db={}   #use this as a general annotation database
    count = 0; x = 0
    original_probesets_add = 0
    if array_type != 'exon':
        for line in open(fn,'rU').xreadlines():             
            probeset_data = cleanUpLine(line)  #remove endline
            probeset,affygene,exons,transcript_num,transcripts,probe_type_call,ensembl,block_exon_ids,block_structure,comparison_info = string.split(probeset_data,'\t')
            ###note: currently exclude comparison_info since not applicable for existing analyses
            if x == 0: x = 1
            else:
                probe_data = affygene,exons,ensembl,block_exon_ids,block_structure
                if exons[-1] == '|': exons = exons[0:-1]
                if affygene[-1] == '|': affygene = affygene[0:-1]
                exon_db[probeset] = affygene
                if probe_type_call == 'gene': #looked through the probe annotations and the gene seems to be the most consistent constituitive feature
                    altmerge_constituitive[probeset] = affygene
    else:
        for line in open(fn,'rU').xreadlines():             
            probeset_data = cleanUpLine(line)  #remove endline
            if x == 0: x = 1
            else:
                probeset_id, exon_id, ensembl_gene_id, transcript_cluster_id, chromosome, strand, probeset_start, probeset_stop, affy_class, constitutitive_probeset, ens_exon_ids, ens_const_exons, exon_region, exon_region_start, exon_region_stop, splicing_event, splice_junctions = string.split(probeset_data,'\t')
                #probeset_id, exon_id, ensembl_gene_id, transcript_cluster_id, chromosome, strand, probeset_start, probeset_stop, affy_class, constitutitive_probeset, ens_exon_ids, exon_annotations = string.split(line,'\t')
                probe_data = ensembl_gene_id
                exon_db[probeset_id] = probe_data
                original_constitutitive_probeset_call = constitutitive_probeset
                if len(splicing_event)>0: constitutitive_probeset = 'no'
                if constitutitive_probeset == 'yes':
                    altmerge_constituitive[probeset_id] = ensembl_gene_id
                    constituitive_gene[ensembl_gene_id]=[]
                if original_constitutitive_probeset_call == 'yes':
                    try: constitutive_original[ensembl_gene_id].append(probeset_id)
                    except KeyError: constitutive_original[ensembl_gene_id] = [probeset_id]

        ###If no constitutive probesets for a gene as a result of additional filtering (removing all probesets associated with a splice event), add these back
        for gene in constitutive_original:
            if gene not in constituitive_gene:
                original_probesets_add +=1
                for probeset in constitutive_original[gene]: altmerge_constituitive[probeset] = gene
    print original_probesets_add, 'genes not viewed as constitutive as a result of filtering probesets based on splicing evidence, added back'
    return altmerge_constituitive,exon_db

def parse_input_data(filename,data_type):
    fn=filepath(filename); first_line = 1; array_group_db = {}; array_group_name_db = {}; z=0
    print "Reading",filename
    for line in open(fn,'rU').xreadlines():
      data = cleanUpLine(line); data2 = string.split(data,'\t'); probeset = data2[0]; z+=1
      if first_line == 1:
          first_line = 0 #makes this value null for the next loop of actual array data
          ###Below ocucrs if the data is raw opposed to precomputed
          if data_type == 'export':
              if array_type == 'exon': folder = 'ExonArray'+'/'+species + '/'
              else: folder = array_type + '/'
              output_file = root_dir+'AltExpression/'+folder + altanalzye_input[0:-4] + '.p' + str(int(100*p)) +'_'+ filter_method+'.txt'
              print "...Exporting",output_file
              export_data = export.createExportFile(output_file,root_dir+'AltExpression/'+folder)
              fn=filepath(output_file); export_data = open(fn,'w');
              export_data.write(line)
          if ':' in data2[1]:
              array_group_list = []; x=0 ###gives us an original index value for each entry in the group
              for entry in data2[1:]:
                  array_group,array_name = string.split(entry,':')
                  try:
                      array_group_db[array_group].append(x)
                      array_group_name_db[array_group].append(array_name)
                  except KeyError:
                      array_group_db[array_group] = [x]
                      array_group_name_db[array_group] = [array_name]
                      ### below only occurs with a new group addition
                      array_group_list.append(array_group) #use this to generate comparisons in the below linked function
                  x += 1
      elif len(probeset)>0 and data_type != 'export':
          ###Use the index values from above to assign each expression value to a new database
          temp_group_array={}; array_index_list = []  ###Use this list for permutation analysis
          for group in array_group_db:
              #array_index_list.append(array_group_db[group])
              group_values = []
              for array_index in array_group_db[group]:
                  try: exp_val = float(data2[array_index+1])
                  except IndexError: print data2, z,'\n',array_index,'\n',group, probeset;kill
                  ###If non-log array data
                  if exp_data_format == 'non-log': exp_val = math.log(exp_val+1,2) #add 1 since probeset is the first column
                  group_values.append(exp_val)
              avg_stat = statistics.avg(group_values)

              if data_type == 'expression':
                  if avg_stat>log_expression_threshold: k=1
                  else: k=0
                  try: expression_status_db[probeset].append(k)
                  except KeyError: expression_status_db[probeset] = [k]
              if data_type == 'p-value':
                  if avg_stat<p: k=1
                  else: k=0
                  if 'G7216513_a_at' in probeset: print k, avg_stat
                  try: pvalue_status_db[probeset].append(k)
                  except KeyError: pvalue_status_db[probeset] = [k]
      elif data_type == 'export':
          try: null = export_db[probeset]; export_data.write(line)
          except KeyError: null = [] ### occurs if not a probeset to include in the filtered results export file
    if data_type == 'export': export_data.close()
    
def expr_analysis(filename,filename2,altmerge_constituitive,exon_db,analyze_dabg):
    """import list of expression values for arrayids and calculates statistics"""
    constitutive_keep={}; keep_probesets={}; keep_genes={}

    global expression_status_db; global pvalue_status_db; global export_db
    expression_status_db={}; pvalue_status_db={}; export_db={}

    parse_input_data(filename,'expression') ### Parse expression file
    if analyze_dabg == 'yes':
        parse_input_data(filename2,'p-value') ### Parse DABG p-value file

    for probeset in expression_status_db:
        proceed = 'no'
        if probeset in pvalue_status_db: proceed = 'yes' ### Indicates there are both expression and dabg files with the same probe sets
        elif analyze_dabg == 'no': proceed = 'yes' ### Indicates there is only an expression file and no dabg
        if proceed == 'yes':
            exp_stat1,exp_stat2 = expression_status_db[probeset]
            if analyze_dabg == 'yes': p_stat1,p_stat2 = pvalue_status_db[probeset]
            else: p_stat1=1; p_stat2=1 ### Automatically assigned an "expressed" call.
            affygene = exon_db[probeset]
            if exp_stat1 == 1 and p_stat1 == 1: k = 1 ### Thus it is "expressed"
            else: k = 0
            if exp_stat2 == 1 and p_stat2 == 1: b = 1 ### Thus it is "expressed"
            else: b = 0
            #if 'G7216513_a_at' in probeset: print b,k,affygene,(probeset,exp_stat1,p_stat1, exp_stat2, p_stat2),pvalue_status_db[probeset];kill
            if probeset in altmerge_constituitive:
                if b == 1 or k == 1:
                    keep_probesets[probeset] = [] ### If either have an "expressed" call, keep... but evaluate constitutive below
                    try: keep_genes[affygene].append(probeset)
                    except KeyError: keep_genes[affygene] = [probeset]
                                
                    if b==1 and k==1: constitutive_keep[affygene] = [] ### This will only keep a gene in the analysis if at least one probeset is 'expressed' in both conditions
            else:
                if b == 0 and k == 0: null =  ''
                else:    
                    keep_probesets[probeset] = []
                    try: keep_genes[affygene].append(probeset)
                    except KeyError: keep_genes[affygene] = [probeset]
    for gene in constitutive_keep:
        probeset_list = keep_genes[gene]
        for probeset in probeset_list: export_db[probeset]=[]
    parse_input_data(filename,'export') ### Parse expression file
    expression_status_db={}; pvalue_status_db={}; export_db={}

def combine_profiles(profile_list):
    profile_group_sizes={}
    for db in profile_list:
        for key in db: profile_group_sizes[key] = len(db[key])
        break

    new_profile_db={}
    for key in profile_group_sizes:
        x = profile_group_sizes[key] ###number of elements in list for key
        new_val_list=[]; i = 0
        while i<x:
            temp_val_list=[]
            for db in profile_list:
                if key in db: val = db[key][i]; temp_val_list.append(val)
            i+=1; val_avg = statistics.avg(temp_val_list); new_val_list.append(val_avg)
        new_profile_db[key] = new_val_list
    return new_profile_db

########### Misc. Functions ###########
def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def remoteRun(Species,Array_type,expression_threshold,filter_method_type,p_val,express_data_format,altanalyze_file_list,Root_dir):
  start_time = time.time()
  global p; global filter_method; global exp_data_format; global array_type; global species; global root_dir
  aspire_output_list=[]; aspire_output_gene_list=[]
  filter_method = filter_method_type
  altanalyze_files = altanalyze_file_list
  p = p_val; species = Species; array_type = Array_type
  exp_data_format = express_data_format
  root_dir = Root_dir

  if 'exon' in array_type: array_type = 'exon' ###In AnalayzeExpressionDataset module, this is named 'exon-array'
  
  global log_expression_threshold
  try: log_expression_threshold = math.log(expression_threshold,2)
  except OverflowError: log_expression_threshold = 0 ###Occurs if expression_threshold == 0
  
  import_dir = '/AltExpression/pre-filtered/expression/'; import_dir_dabg = '/AltExpression/pre-filtered/dabg/'
  
  dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
  try: dir_list2 = read_directory(import_dir_dabg)
  except Exception: dir_list2=[]

  if len(altanalyze_files) == 0: altanalyze_files = dir_list  ###if no filenames input

  if array_type == 'exon': altmerge_db = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_probesets.txt"
  else: altmerge_db = "AltDatabase/"+species+"/"+array_type+"/MASTER-probeset-transcript.txt"
  ###Import probe-level associations
  altmerge_constituitive,exon_db = import_altmerge(altmerge_db,array_type)
        
  global altanalzye_input
  for altanalzye_input in dir_list:    #loop through each file in the directory to output results
    if altanalzye_input in altanalyze_files:
        if altanalzye_input in dir_list2: analyze_dabg = 'yes'
        else: analyze_dabg = 'no'
        ind_start_time = time.time()
        array_db = import_dir + "/"+ altanalzye_input
        dabg_db = import_dir_dabg + "/"+ altanalzye_input
        array_db = array_db[1:] #not sure why, but the '\' needs to be there while reading initally but not while accessing the file late
        dabg_db = dabg_db[1:]
        dataset_name = altanalzye_input[0:-4] + '-'
        print "Begining to process",dataset_name[0:-1]
        print "Array type is:",array_type
        print "Species is:", species
        print "Expression format is:",exp_data_format
        print "DABG p-value cut-off is:",p
        print "Filter method is:",filter_method
        print "Log2 expression cut-off is:",log_expression_threshold

        ###Import expression data and stats
        try: expr_analysis(array_db,dabg_db,altmerge_constituitive,exon_db,analyze_dabg)    #filter the expression data based on fold and p-value OR expression threshold
        except KeyError: print "Impropper array type (",dataset_name[0:-1],") for",array_type,species,'. Skipping array.'
        ind_end_time = time.time(); time_diff = int(ind_end_time-ind_start_time)
        print dataset_name,"filtering finished in %d seconds" % time_diff
  end_time = time.time(); time_diff = int(end_time-start_time)
  print "Filtering complete for all files in %d seconds" % time_diff
  exon_db={}; altmerge_constituitive={}
  
if __name__ == '__main__':
  m = 'Mm'; h = 'Hs'
  species = m
  x = 'AltMouse'; y = 'exon'
  array_type = y
  l = 'log'; n = 'non-log'
  exp_data_format = l  
  p = 0.75
  filter_method = 'average'
  expression_threshold = 0
  
  print "Filter Data For:"
  print "1) Human 1.0 ST exon data\n2) AltMouse"
  inp = sys.stdin.readline(); inp = inp.strip()
  if inp == "1": species = h; p = 0.05; array_type = y; expression_threshold = 70
  elif inp == "2": species = m; p = 0.75; expression_threshold = 0

  altanalyze_files = []
  remoteRun(species,array_type,expression_threshold,filter_method,p,exp_data_format,altanalyze_files)