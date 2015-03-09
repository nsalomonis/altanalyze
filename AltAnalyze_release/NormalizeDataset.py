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
import statistics
import sys, string
import os.path
import unique
import UI
import export
import time
import traceback
import RNASeq
    
def normalizeDataset(filename,output = None):
    """ Perform Quantile Normalization on an input expression dataset """
    
    print "Importing data..."
    sample_expression_db = importExpressionValues(filename)
    print "Performing quantile normalization..."
    sample_expression_db = RNASeq.quantileNormalizationSimple(sample_expression_db)
    
    if output==None:
        output = filename
        moved_exp_dir = export.findParentDir(filename)+'Non-Quantile/'+export.findFilename(filename)
        try:
            export.copyFile(filename, moved_exp_dir)
            print 'Moved original expression file to:'
            print '\t'+moved_exp_dir
        except Exception: None
    
    exportExpressionData(output,sample_expression_db)
    print 'Exported expression input file to:',output
            
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
    filename = "/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/NODAT raw exp/Normalize/diabetes_type2.txt"
    filename = '/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/qPCR/deltaCT/LabMeeting/ExpressionInput/exp.Adult_ABI_dCT_viia7-Ensembl-STA.txt'
    normalizeDataset(filename)