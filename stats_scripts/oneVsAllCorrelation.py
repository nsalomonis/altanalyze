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
import statistics
import export
import copy
from visualization_scripts import clustering
from scipy import stats

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def getSimpleCorrelations(filename,ref_gene):    
    X, column_header, row_header, dataset_name, group_db = clustering.importData(filename)
    i=0
    
    ri = row_header.index(ref_gene)
    ref_values = X[ri]
    
    dir = export.findParentDir(filename)
    eos = export.ExportFile(dir+'/results/'+ref_gene+'.txt')
    
    for ls in X:
        try: rho,p = stats.pearsonr(ls,ref_values)
        except: rho = -1
        if rho>0.2:
            gene = row_header[i]
            eos.write(string.join([gene,ref_gene,str(rho)],'\t')+'\n')
        i+=1
    eos.close()
    
if __name__ == '__main__':
    
    ################  Comand-line arguments ################
    import getopt

    #python stats_scripts/ADT.py --i /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/CPTT/exp.Biolegend-ADT-1.txt --a /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/feature_reference.txt --g /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/CPTT/AltAnalyze/cellHarmony/cellHarmony/QueryGroups.cellHarmony_captures-filtered.txt --s /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Isolation-Strategy/HSCP.txt
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print 'WARNING!!!! Too commands supplied.'
        
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','g='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i':
                exp_file = arg
            if opt == '--g':
                gene = arg
                
    getSimpleCorrelations(exp_file,gene)