###ProcessAgilentArrays
#Copyright 2012 J. David Gladstone Institutes, San Francisco California
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
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
from stats_scripts import statistics
import os.path
import unique
import UI
import export
import time
import traceback
import RNASeq
    
def agilentSummarize(exp_file_location_db):
    print 'Agilent array import started'
    
    global red_channel_db
    global green_channel_db
    red_channel_db={}
    green_channel_db={}
    
    for dataset in exp_file_location_db: ### Instance of the Class ExpressionFileLocationData
        fl = exp_file_location_db[dataset]
        output_dir = fl.OutputDir()
        array_dir=fl.CELFileDir()
        group_dir = fl.GroupsFile() ### provides the list of array_files
        channel_to_extract = fl.ChannelToExtract()
        expression_file = fl.ExpFile()
        array_group_list = UI.importArrayGroupsSimple(group_dir,[])[0]
        normalization_method = fl.NormMatrix()

    arrays = map(lambda agd: agd.Array(), array_group_list) ### Pull the array names out of this list of objects

    dir_list = unique.read_directory(array_dir)
    count=0
    for array in dir_list:
        if array in arrays: ### Important since other text files may exist in that directory
            count+=1
            filename = array_dir+'/'+array
            importAgilentExpressionValues(filename,array,channel_to_extract)
            if count == 50:
                print '' ### For progress printing
                count = 0
            
    if len(green_channel_db)>0:
        filename = output_dir+ '/'+ 'gProcessed/gProcessed-'+dataset+'-raw.txt'
        exportExpressionData(filename,green_channel_db)
        if 'quantile' in normalization_method:
            print '\nPerforming quantile normalization on the green channel...'
            green_channel_db = RNASeq.quantileNormalizationSimple(green_channel_db)
            filename = output_dir+ '/'+ 'gProcessed/gProcessed-'+dataset+'-quantile.txt'
            exportExpressionData(filename,green_channel_db)
        final_exp_db = green_channel_db
        
    if len(red_channel_db)>0:
        filename = output_dir+ '/'+ 'rProcessed/rProcessed-'+dataset+'-raw.txt'
        exportExpressionData(filename,red_channel_db)
        if 'quantile' in normalization_method:
            print '\nPerforming quantile normalization on the red channel...'
            red_channel_db = RNASeq.quantileNormalizationSimple(red_channel_db)
            filename = output_dir+ '/'+ 'rProcessed/rProcessed-'+dataset+'-quantile.txt'
            exportExpressionData(filename,red_channel_db)
        final_exp_db = red_channel_db

    if len(red_channel_db)>0 and len(green_channel_db)>0:
        if channel_to_extract == 'green/red ratio':
            final_exp_db = calculateRatios(green_channel_db,red_channel_db)
        elif channel_to_extract == 'red/green ratio':
            final_exp_db = calculateRatios(red_channel_db,green_channel_db)
            
    exportExpressionData(expression_file,final_exp_db)
    print 'Exported expression input file to:',expression_file
            
def calculateRatios(db1,db2):
    ratio_db={}
    for array in db1:
        exp_ratios={}
        exp_db = db1[array]
        for probe_name in exp_db:
            exp_ratios[probe_name] = str(float(exp_db[probe_name])-float(db2[array][probe_name])) ### log2 ratio
        ratio_db[array]=exp_ratios
    return ratio_db

def importAgilentExpressionValues(filename,array,channel_to_extract):
    """ Imports Agilent Feature Extraction files for one or more channels """
    print '.',
    red_expr_db={}
    green_expr_db={}
    parse=False
    fn=unique.filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = UI.cleanUpLine(line)
        if parse==False:
            if 'ProbeName' in data:
                headers = string.split(data,'\t')
                pn = headers.index('ProbeName')
                try: gc = headers.index('gProcessedSignal')
                except Exception: pass
                try: rc = headers.index('rProcessedSignal')
                except Exception: pass
                parse = True
        else:
            t = string.split(data,'\t')
            probe_name = t[pn]
            try: green_channel = math.log(float(t[gc])+1,2) #min is 0
            except Exception: pass
            try: red_channel = math.log(float(t[rc])+1,2) #min is 0
            except Exception: pass
            if 'red' in channel_to_extract:
                red_expr_db[probe_name] = red_channel
            if 'green' in channel_to_extract:
                green_expr_db[probe_name] = green_channel

    if 'red' in channel_to_extract:
        red_channel_db[array] = red_expr_db
    if 'green' in channel_to_extract:
        green_channel_db[array] = green_expr_db
                
def exportExpressionData(filename,sample_db):
    export_text = export.ExportFile(filename)
    all_genes_db = {}
    sample_list=[]
    for sample in sample_db:
        sample_list.append(sample)
        gene_db = sample_db[sample]
        for geneid in gene_db:
            all_genes_db[geneid]=[]
    sample_list.sort() ### Organize these alphabetically rather than randomly
    column_header = string.join(['ProbeName']+sample_list,'\t')+'\n' ### format column-names for export
    export_text.write(column_header)

    for geneid in all_genes_db:
        values=[]
        for sample in sample_list:
            try: values.append(sample_db[sample][geneid]) ### protein_expression
            except Exception: values.append(0)
        export_text.write(string.join([geneid]+map(str, values),'\t')+'\n')
    export_text.close()
    
if __name__ == '__main__':
    None