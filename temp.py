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

def importAnnotations(filename):
    firstLine = True
    fn = filepath(filename)
    rows = 0
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line);
        tab_delimited_data = string.split(data,'\t')
        if rows > 10: sys.exit()
        print tab_delimited_data#;sys.exit()
        rows+=1

def importMethylationData(filename,betaLow=0.4,betaHigh=0.6,counts=1):
    
    firstLine = True
    rows=0
    for line in open(filename,'rU').xreadlines():             
        data = cleanUpLine(line);
        t = string.split(data,'\t')
        rows+=1
        if rows<10:
            print [line]
        else:
            sys.exit()
        rows+=1

def conFloat(x,betaValues):
    try: x = float(x)
    except Exception:  x=None
    if x== None or x == 0:
        floats=[]
        for i in betaValues:
            if i=='': pass
            elif float(i)==0: pass
            else: floats.append(float(i))
        try: return min(floats)
        except Exception: print betaValues;sys.exit()
    else:
        return x

def betaHighCount(x,betaHigh):
    if x>betaHigh:
        return 1
    else: return 0
    
def betaLowCount(x,betaLow):
    if x<betaLow:
        return 1
    else: return 0


if __name__ == '__main__':
    filename = 'AltDatabase/ucsc/Hs/wgEncodeHaibMethyl450CpgIslandDetails.txt'
    input_file = '/Volumes/SEQ-DATA/Kamath/BedFiles/ExpressionInput/URSA/ERX011182.pcl'
    Output_types = 'promoter'
    Output_types = 'all'
    Species = 'Hs'
    importMethylationData(input_file); sys.exit()
    importAnnotations(filename);sys.exit()
    
    
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
    
    
