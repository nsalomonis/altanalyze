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
import collections
import sys,string,os
import statistics
import export
import unique
import copy
from visualization_scripts import clustering
import UI
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

def reformatOrganizedDifferentials(fn,datasetName,organizedDiffGene_db):
    firstRow = True
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if firstRow:
            firstRow=False
        else:
            t = string.split(data,'\t')
            uid = t[0]
            geneSymbol = string.split(uid,":")[1]
            if geneSymbol not in organizedDiffGene_db:
                organizedDiffGene_db[geneSymbol]=uid
    return organizedDiffGene_db

def importFolds(fn,organizedDiffGene_db,fold_db,datasetName,header_db):
    firstRow = True
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            header_temp = t[1:]
            header=[]
            for h in header_temp:
                #cell_type = string.join(string.split(h,'_')[:-1],'_')
                cell_type = string.split(h,'_')[0]
                header.append(cell_type)
                try: header_db[cell_type]+=1
                except: header_db[cell_type]=1
            firstRow=False
        else:
            geneSymbol = t[0]
            folds = t[1:]
            if geneSymbol in organizedDiffGene_db:
                i=0
                fold_objects=collections.OrderedDict()
                for fold in folds:
                    celltype = header[i]
                    fold_objects[celltype]=fold
                    i+=1
                if geneSymbol in fold_db:
                    dataset_folds_db = fold_db[geneSymbol]
                    dataset_folds_db[datasetName] = fold_objects
                else:
                    dataset_folds_db = collections.OrderedDict()
                    dataset_folds_db[datasetName] = fold_objects
                    fold_db[geneSymbol] = dataset_folds_db
    return fold_db,header_db
        
def getDatasetName(other_files_dir):
    files = unique.read_directory(other_files_dir+'/')
    for file in files:
        if 'exp.' in file and 'AllCells-folds.txt' in file:
            folds_file = other_files_dir+'/'+file
            datasetName = file[4:-19]
    return folds_file, datasetName
        
def combine(organized_differentials,species,output_dir):
    datasets=[]
    fold_files=[]
    organizedDiffGene_db = collections.OrderedDict()
    
    for OD in organized_differentials:
        cellHarmony_dir = export.findParentDir(OD)
        folds_file, datasetName = getDatasetName(cellHarmony_dir+'/OtherFiles/')
        organizedDiffGene_db = reformatOrganizedDifferentials(OD,datasetName,organizedDiffGene_db)
        fold_files.append(folds_file)
        datasets.append(datasetName)
    
    i=0
    fold_db = collections.OrderedDict()
    header_db = collections.OrderedDict()
    for FF in fold_files:
        datasetName = datasets[i]
        fold_db,header_db = importFolds(FF,organizedDiffGene_db,fold_db,datasetName,header_db)
        i+=1
        
    final_celltypes = []
    for celltype in header_db:
        print celltype, header_db[celltype]
        if header_db[celltype]==i:
            final_celltypes.append(celltype) ### Hence the celltype is present in all dataset comparisons
    print final_celltypes
    header=['UID']
    groups=[]
    export_header=True
    export_file = output_dir+'/exp.combinedCellHarmony.txt'
    groups_file = output_dir+'/groups1.combinedCellHarmony.txt'
    eo = export.ExportFile(export_file)
    eog = export.ExportFile(groups_file)

    for gene in organizedDiffGene_db:
        folds=[]
        fullGeneName = organizedDiffGene_db[gene]
        dataset_folds_db = fold_db[gene]
        for celltype in final_celltypes:
            for datasetName in datasets:
                fold = dataset_folds_db[datasetName][celltype]
                folds.append(fold)
                header.append(celltype+':'+celltype+'-'+datasetName)
                groups.append(celltype+'-'+datasetName+'\t'+datasetName+'\t'+datasetName+'\n')
        if export_header:
            export_header=False
            eo.write(string.join(header,'\t')+'\n')
            for g in groups: ### export a groups file to denote which comparisons each fold derives from
                eog.write(g)
        eo.write(string.join([fullGeneName]+folds,'\t')+'\n')
    eo.close()
    eog.close()
    
    gsp = UI.GeneSelectionParameters(species,'RNASeq','RNASeq')
    gsp.setPathwaySelect('None Selected')
    gsp.setGeneSelection('')
    gsp.setOntologyID('')
    gsp.setGeneSet('None Selected')
    gsp.setJustShowTheseIDs('')                
    gsp.setTranspose(False)
    gsp.setNormalize('NA')
    gsp.setGeneSelection('')
    #gsp.setClusterGOElite('GeneOntology')
    gsp.setClusterGOElite('PathwayCommons')
    row_method = None; row_metric = 'correlation'; column_method = None; column_metric = 'cosine'; color_gradient = 'yellow_black_blue'
    transpose = False; Normalize='NA'

    print 'Producing a heatmap'
    graphics = clustering.runHCexplicit(export_file, [], row_method, row_metric,
                column_method, column_metric, color_gradient, gsp, Normalize=Normalize,
                contrast=7, display=False)

if __name__ == '__main__':
    
    ################  Comand-line arguments ################
    import getopt

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