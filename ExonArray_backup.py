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

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

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
    

    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line);

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
    
    
