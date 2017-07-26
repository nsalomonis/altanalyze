###NormalizeDataset
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
import ExpressionBuilder
    
def normalizeDataset(filename,output = None, normalization='quantile',platform="3'array"):
    """ Perform Quantile Normalization on an input expression dataset """
    
    if output==None:
        output = filename
        moved_exp_dir = export.findParentDir(filename)+'Non-Normalized/'+export.findFilename(filename)
        try:
            export.copyFile(filename, moved_exp_dir)
            print 'Moved original expression file to:'
            print '\t'+moved_exp_dir
        except Exception: None
        
    if normalization == 'Quantile' or normalization == 'quantile':
        print "Importing data..."
        sample_expression_db = importExpressionValues(filename)
        print "Performing quantile normalization..."    
        sample_expression_db = RNASeq.quantileNormalizationSimple(sample_expression_db)
        exportExpressionData(output,sample_expression_db)
    elif normalization == 'group':
        performGroupNormalization(moved_exp_dir,filename,platform)    
    print 'Exported expression input file to:',output
            
def importGroups(fn):
    try: group_db=collections.OrderedDict()
    except Exception:
        try:
            import ordereddict
            group_db=ordereddict.OrderedDict()
        except Exception: group_db={}
    for line in open(fn,'rU').xreadlines():
        data = ExpressionBuilder.cleanUpLine(line)
        sample_filename,group_number,group_name = string.split(data,'\t')
        try: group_db[group_name].append(sample_filename)
        except Exception: group_db[group_name] = [sample_filename]
    return group_db

def performGroupNormalization(filename,export_dir,platform):
    expressionDataFormat,increment,convertNonLogToLog = ExpressionBuilder.checkExpressionFileFormat(filename)
    groups_dir = string.replace(export_dir,'exp.','batch.')
    fn=unique.filepath(filename); row_number=0; exp_db={}; relative_headers_exported = False
    group_db = importGroups(groups_dir)
    export_data = export.ExportFile(export_dir)
    for line in open(fn,'rU').xreadlines():
        data = ExpressionBuilder.cleanUpLine(line)
        t = string.split(data,'\t')
        if data[0]=='#' and row_number==0: row_number = 0
        elif row_number==0:
            sample_list = t[1:]
            new_sample_list = []
            for group in group_db:
                group_samples = group_db[group]
                try:
                    sample_index_list = map(lambda x: sample_list.index(x), group_samples)
                    group_db[group] = sample_index_list
                    new_sample_list+=group_samples
                except Exception:
                    missing=[]
                    for x in sample_list:
                        if x not in t[1:]: missing.append(x)
                    print 'missing:',missing
                    print t
                    print sample_list
                    print filename, groups_dir
                    print 'Unknown Error!!! Skipping cluster input file build (check column and row formats for conflicts)'; forceExit
            title = string.join([t[0]]+new_sample_list,'\t')+'\n' ### output the new sample order (group file order)
            export_data.write(title)
            row_number=1
        else:
            gene = t[0]
            if expressionDataFormat == 'non-log' and (convertNonLogToLog or platform == 'RNASeq'):
                ### Convert to log2 RPKM values - or counts
    
                try: all_values = map(lambda x: math.log(float(x)+increment,2), t[1:])
                except Exception:
                    all_values = ExpressionBuilder.logTransformWithNAs(t[1:],increment)
            else:
                try: all_values = map(float,t[1:])
                except Exception:
                    all_values = ExpressionBuilder.logTransformWithNAs(t[1:],increment)
            row_number+=1 ### Keep track of the first gene as to write out column headers for the relative outputs
            gene_log_folds = []

            for group in group_db:
                sample_index_list = group_db[group]
                ### Calculate log-fold values relative to the mean of all sample expression values
                try: values = map(lambda x: all_values[x], sample_index_list) ### simple and fast way to reorganize the samples
                except Exception:
                    print len(values), sample_index_list;kill
                try: avg = statistics.avg(values)
                except Exception:
                    values2=[]
                    for v in values:
                        try: values2.append(float(v))
                        except Exception: pass
                    values = values2
                    try: avg = statistics.avg(values)
                    except Exception:
                        if len(values)>0: avg = values[0]
                        else: avg = 0
                try: log_folds = map(lambda x: (x-avg), values)
                except Exception: 
                    log_folds=[]
                    for x in values:
                        try: log_folds.append(x-avg)
                        except Exception: log_folds.append('')
                gene_log_folds+=log_folds                            
            gene_log_folds = map(lambda x: str(x),gene_log_folds)
            export_data.write(string.join([gene]+gene_log_folds,'\t')+'\n')
    export_data.close()

def calculateRatios(db1,db2):
    ratio_db={}
    for array in db1:
        exp_ratios={}
        exp_db = db1[array]
        for probe_name in exp_db:
            exp_ratios[probe_name] = str(float(exp_db[probe_name])-float(db2[array][probe_name])) ### log2 ratio
        ratio_db[array]=exp_ratios
    return ratio_db

def importExpressionValues(filename):
    """ Imports tab-delimited expression values"""

    header = True
    sample_expression_db={}
    fn=unique.filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = UI.cleanUpLine(line)
        if header:
            sample_names = string.split(data,'\t')
            header = False
        else:
            exp_values = string.split(data,'\t')
            gene = exp_values[0]
            index=1
            for value in exp_values[1:]:
                sample_name = sample_names[index]
                if sample_name in sample_expression_db:
                    gene_expression_db = sample_expression_db[sample_name]
                    gene_expression_db[gene] = value
                else:
                    gene_expression_db={}
                    gene_expression_db[gene] = value
                    sample_expression_db[sample_name] = gene_expression_db
                index+=1
    return sample_expression_db
                
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
    filename = "/Volumes/salomonis1/projects/Beena2/fastq/ExpressionInput/exp.AF.txt'"
    normalizeDataset(filename,normalization='Group')
    #normalizeDataset(filename)