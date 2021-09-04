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

def iterativeMarkerFinder(root_dir,dataType='PSI',geneRPKM=0):
    """ Iteratively perform MarkerFinder and combine the results """
    import ExpressionBuilder, shutil
    species = 'Hs'
    
    if dataType == 'PSI':
        platform = 'exon'
    else:
        platform = dataType
    import markerFinder
    import collections
    consolidated_MarkerFinder = collections.OrderedDict()
    ordered_groups = []
    files = UI.read_directory(root_dir+'/ExpressionInput')
    files.sort()
    for file in files:
        if 'exp.' in file:
            print file
            marker_dir = root_dir+'/ExpressionInput/'+file
            group_dir = marker_dir.replace('exp.','groups.')
            comps_dir = marker_dir.replace('exp.','comps.')
            sample_group_db = ExpressionBuilder.simplerGroupImport(group_dir)
            for sample in sample_group_db:
                group = sample_group_db[sample]
                if group not in ordered_groups: ordered_groups.append(group)
                
            fl = UI.ExpressionFileLocationData(marker_dir,'','',''); fl.setOutputDir(root_dir)
            fl.setSpecies(species); fl.setVendor(dataType); fl.setVendor(dataType)
            fl.setRPKMThreshold(geneRPKM)
            fl.setCorrelationDirection('up')
            logTransform = False
            markerFinder.analyzeData(marker_dir,species,platform,'protein_coding',geneToReport=50,correlateAll=True,AdditionalParameters=fl,logTransform=logTransform)
            ICGS_State_ranked = importMarkerFinderHits(root_dir+'/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt',dataType)
            shutil.copy(root_dir+'/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt',root_dir+'/ExpressionOutput/MarkerFinder/'+file[:-4]+'-up.txt')
            consolidated_MarkerFinder[file[:-4],'up']=ICGS_State_ranked
            fl.setCorrelationDirection('down')
            markerFinder.analyzeData(marker_dir,species,platform,'protein_coding',geneToReport=50,correlateAll=True,AdditionalParameters=fl,logTransform=logTransform)
            ICGS_State_ranked = importMarkerFinderHits(root_dir+'/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt',dataType)
            shutil.copy(root_dir+'/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt',root_dir+'/ExpressionOutput/MarkerFinder/'+file[:-4]+'-down.txt')
            consolidated_MarkerFinder[file[:-4],'down']=ICGS_State_ranked
            #graphics_mf = markerFinder.generateMarkerHeatMaps(fl,dataType,convertNonLogToLog=logTransform,Species=species)
    
    ### Reorganize groups
    try:
        if 'del' in ordered_groups[3] or 'Del' in ordered_groups[3]:
            deleted = ordered_groups[3]
            del ordered_groups[3]
            ordered_groups.append(deleted)
    except: pass
    organizeConsolidatedMarkerFinder(consolidated_MarkerFinder,ordered_groups,root_dir,marker_dir)
    
def organizeConsolidatedMarkerFinder(consolidated_MarkerFinder,ordered_groups,root_dir,marker_dir):
    organized_patterns={}
    for (file,direction) in consolidated_MarkerFinder:
        ICGS_State_ranked = consolidated_MarkerFinder[(file,direction)]
        for cell_state in ICGS_State_ranked:
            for (rho,gene,symbol) in ICGS_State_ranked[cell_state]:
                try: organized_patterns[gene].append([rho,cell_state,direction])
                except: organized_patterns[gene] = [[rho,cell_state,direction]]
    ranked_cell_state_genes={}
    count=0

    for gene in organized_patterns:
        gene_pattern = unique.unique(organized_patterns[gene])
        gene_pattern.sort()
        rho, cell_state, direction = gene_pattern[-1]
        try:
            rho2, cell_state2, direction2 = gene_pattern[-2]
            if rho == rho2:
                ### Occurs with oppositive patterns solving the same objective
                if direction2 == 'up':
                    rho, cell_state, direction = rho2, cell_state2, direction2
                    if direction2 != direction:
                        print gene, rho, cell_state, direction, rho2, cell_state2, direction2
        except: pass ### No other instances with a rho>cutoff
        try: ranked_cell_state_genes[direction,cell_state].append([rho,gene])
        except Exception: ranked_cell_state_genes[direction,cell_state] = [[rho,gene]]
        count+=1
    
    ordered_patterns=[]
    for cell_state in ordered_groups:
        try: ranked_cell_state_genes['up',cell_state].sort()
        except:
            print 'up',cell_state, '---failed'
            continue
        ranked_cell_state_genes['up',cell_state].reverse()
        for (rho,gene) in ranked_cell_state_genes['up',cell_state]:
            ordered_patterns.append(['up-'+cell_state,gene])
    for cell_state in ordered_groups:
        try: ranked_cell_state_genes['down',cell_state].sort()
        except:
            print 'down',cell_state, '---failed'
            continue
        ranked_cell_state_genes['down',cell_state].reverse()
        for (rho,gene) in ranked_cell_state_genes['down',cell_state]:
            ordered_patterns.append(['down-'+cell_state,gene])
    #print 'Number of events matching the MarkerFinder and stastical cutoffs:',len(ordered_patterns)
    
    export_file = root_dir+'/DataPlots/Consolidated-MarkerFinder.txt'
    eo = export.ExportFile(export_file)
    matrix, column_header, row_header, dataset_name, group_db = clustering.importData(marker_dir)
    revised_column_headers = ['UID']
    for i in column_header:
        revised_column_headers.append(group_db[i][0] + ':' + i)
    eo.write(string.join(revised_column_headers,'\t')+'\n')
    for (pattern,gene) in ordered_patterns:
        i = row_header.index(gene)
        gene = string.replace(gene,':','__')
        eo.write(string.join([pattern+':'+gene]+map(str,matrix[i]),'\t')+'\n')
    eo.close()

    row_method = None; row_metric = 'correlation'; column_method = None; column_metric = 'cosine'; color_gradient = 'yellow_black_blue'
    transpose = False; Normalize='median'; #gsp.setClusterGOElite('PathwayCommons')

    graphics = clustering.runHCexplicit(export_file, [], row_method, row_metric,
                column_method, column_metric, color_gradient, transpose, Normalize=Normalize,
                contrast=10, display=False)
    
def importMarkerFinderHits(fn,dataType):
    if dataType == 'PSI': cutoff = 0.5
    elif dataType == 'ADT': cutoff = 0.2
    elif dataType == 'RNASeq': cutoff = 0.7
    else: cutoff = 0.6
    print "Using a MarkerFinder Pearson rho >",cutoff
    genes={}
    ICGS_State_ranked={}
    skip=True
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if skip: skip=False
        else:
            try:
                gene,symbol,rho,ICGS_State = string.split(data,'\t')
            except Exception:
                gene,symbol,rho,rho_p,ICGS_State = string.split(data,'\t')
            if float(rho)>cutoff:
                try: ICGS_State_ranked[ICGS_State].append([float(rho),gene,symbol])
                except Exception: ICGS_State_ranked[ICGS_State] = [[float(rho),gene,symbol]]
    return ICGS_State_ranked

if __name__ == '__main__':
    
    ################  Comand-line arguments ################
    import getopt

    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print 'WARNING!!!! Too commands supplied.'
        
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','p='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i':
                input_dir = arg
            if opt == '--p':
                dataType = arg
                
    iterativeMarkerFinder(input_dir,dataType='PSI',geneRPKM=0);sys.exit()