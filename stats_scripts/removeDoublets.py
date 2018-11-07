### removeDoublets
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

import sys,string,os,copy, path
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

command_args = string.join(sys.argv,' ')
if len(sys.argv[1:])>0 and '--' in command_args: commandLine=True
else: commandLine=False

display_label_names = True

import traceback

def removeMarkerFinderDoublets(heatmap_file,diff=1):
    matrix, column_header, row_header, dataset_name, group_db, priorColumnClusters, priorRowClusters = remoteImportData(heatmap_file)
    
    priorRowClusters.reverse()
    if len(priorColumnClusters)==0:
        for c in column_header:
            cluster = string.split(c,':')[0]
            priorColumnClusters.append(cluster)
        for r in row_header:
            cluster = string.split(r,':')[0]
            priorRowClusters.append(cluster)

    import collections
    cluster_db = collections.OrderedDict()
    i=0
    for cluster in priorRowClusters:
        try: cluster_db[cluster].append(matrix[i])
        except: cluster_db[cluster] = [matrix[i]]
        i+=1
    
    transposed_data_matrix=[]
    clusters=[]
    for cluster in cluster_db:
        cluster_cell_means = numpy.mean(cluster_db[cluster],axis=0)
        cluster_db[cluster] = cluster_cell_means
        transposed_data_matrix.append(cluster_cell_means)
        if cluster not in clusters:
            clusters.append(cluster)
    transposed_data_matrix = zip(*transposed_data_matrix)
    
    i=0
    cell_max_scores=[]
    cell_max_score_db = collections.OrderedDict()

    for cell_scores in transposed_data_matrix:
        cluster = priorColumnClusters[i]
        cell = column_header[i]
        ci = clusters.index(cluster)
        #print ci, cell, cluster, cell_scores;sys.exit()
        cell_state_score = cell_scores[ci] ### This is the score for that cell for it's assigned MarkerFinder cluster
        alternate_state_scores=[]
        for score in cell_scores:
            if score != cell_state_score:
                alternate_state_scores.append(score)
        alt_max_score = max(alternate_state_scores)
        alt_sum_score = sum(alternate_state_scores)
        cell_max_scores.append([cell_state_score,alt_max_score,alt_sum_score]) ### max and secondary max score - max for the cell-state should be greater than secondary max
        try: cell_max_score_db[cluster].append(([cell_state_score,alt_max_score,alt_sum_score]))
        except: cell_max_score_db[cluster] = [[cell_state_score,alt_max_score,alt_sum_score]]
        i+=1
    
    
    for cluster in cell_max_score_db:
        cluster_cell_means = numpy.median(cell_max_score_db[cluster],axis=0)
        cell_max_score_db[cluster] = cluster_cell_means ### This is the cell-state mean score for all cells in that cluster and the alternative max mean score (difference gives you the threshold for detecting double)
    i=0
    print len(cell_max_scores)
    keep=['row_clusters-flat']
    keep_alt=['row_clusters-flat']
    for (cell_score,alt_score,alt_sum) in cell_max_scores:
        cluster = priorColumnClusters[i]
        cell = column_header[i]
        ref_max, ref_alt, ref_sum = cell_max_score_db[cluster]
        ci = clusters.index(cluster)
        ref_diff= math.pow(2,(ref_max-ref_alt))*diff #1.1
        ref_alt = math.pow(2,(ref_alt))
        cell_diff = math.pow(2,(cell_score-alt_score))
        cell_score = math.pow(2,cell_score)
        if cell_diff>ref_diff and cell_diff>diff: #cell_score cutoff removes some, but cell_diff is more crucial
                #if alt_sum<cell_score:
                assignment=0 #1.2
                keep.append(cell)
                keep_alt.append(string.split(cell,':')[1]) ### if prefix added
        else:
            assignment=1

        #print assignment
        i+=1
    print len(keep)
    from import_scripts import sampleIndexSelection
    input_file=heatmap_file
    output_file = heatmap_file[:-4]+'-ReOrdered.txt'
    try: sampleIndexSelection.filterFile(input_file,output_file,keep)
    except: sampleIndexSelection.filterFile(input_file,output_file,keep_alt)

if __name__ == '__main__':
    import getopt
    filter_rows=False
    filter_file=None
    genome = 'hg19'
    dataset_name = '10X_filtered'
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Insufficient options provided";sys.exit()
        #Filtering samples in a datasets
        #python 10XProcessing.py --i /Users/test/10X/outs/filtered_gene_bc_matrices/ --g hg19 --n My10XExperiment
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','g=','n='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i': matrices_dir=arg
            elif opt == '--g': genome=arg
            elif opt == '--n': dataset_name=arg
    removeMarkerFinderDoublets('/Users/saljh8/Desktop/dataAnalysis/Collaborative/Ratner-10X/2mo-NF/Mm_Run2050/Final.8.11.2018/DataPlots/MarkerFinder/Clustering-MarkerGenes_correlations-ReplicateBased-20180812-091909_euclidean_cosine.txt',diff=diff);sys.exit()
    