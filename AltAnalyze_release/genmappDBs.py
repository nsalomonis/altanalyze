###genmapp
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
import math
import os.path
import unique
import copy
import time
import export

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def makeUnique(item):
    db1={}; list1=[]; k=0
    for i in item:
        try: db1[i]=[]
        except TypeError: db1[tuple(i)]=[]; k=1
    for i in db1:
        if k==0: list1.append(i)
        else: list1.append(list(i))
    list1.sort()
    return list1

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def getAnnotations(species):
    exportChromosomeStrandCoordinates(species)

def exportChromosomeStrandCoordinates(species):
    import EnsemblImport
    gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,'RNASeq','key_by_array')

    import ExpressionBuilder
    gene_biotype_db = ExpressionBuilder.importTranscriptBiotypeAnnotations(species)
    export_path = 'GenMAPPDBs/'+species+'/chr_gene_locations.txt'
    export_data = export.ExportFile(export_path)

    import ExonAnalyze_module
    gene_annotation_file = "AltDatabase/ensembl/"+species+"/"+species+"_Ensembl-annotations.txt"
    annotate_db = ExonAnalyze_module.import_annotations(gene_annotation_file,'RNASeq')
      
    print 'Annotations for',len(gene_location_db),'genes imported'
    
    sorted_list=[]; protein_coding=0 
    for gene in gene_location_db:
        chr,strand,start,end = gene_location_db[gene]
        if gene in gene_biotype_db:
            biotype = gene_biotype_db[gene][-1]
            if biotype == 'protein_coding': protein_coding+=1
                
        else: biotype = 'NA'
        if len(chr)<7:
            sorted_list.append([chr,strand,int(start),int(end),gene,biotype])
        #else: print chr;sys.exit()
    print len(sorted_list),'genes for typical chromosomes present'
    print protein_coding, 'protein coding genes present'
    sorted_list.sort()        
    for values in sorted_list:
        chr,strand,start,end,gene,biotype=values
        try: symbol = annotate_db[gene].Symbol()
        except Exception: symbol = ''
        values = [gene,symbol,chr,strand,str(start),str(end),biotype]
        export_data.write(string.join(values,'\t')+'\n')
    export_data.close()
    print species, 'chromosome locations exported to:\n',export_path
    
if __name__ == '__main__':
    getAnnotations('Sc')
    
