import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
from scipy import sparse, io
import numpy
import LineageProfilerIterate
import cluster_corr

""" cellHarmony with Louvain Clustering Funcitons """

def manage_louvain_alignment(species,platform,query_exp_file,exp_output,
            customMarkers=False,useMulti=False,fl=None,customLabels=None):

    if 'ICGS' in customMarkers or 'MarkerGene' in customMarkers:
        """ When performing cellHarmony, build an ICGS expression reference with log2 TPM values rather than fold """
        print 'Converting ICGS folds to ICGS expression values as a reference first...'
        try: customMarkers = LineageProfilerIterate.convertICGSClustersToExpression(customMarkers,query_exp_file,returnCentroids=False,species=species)
        except:
            print "Using the supplied reference file only (not importing raw expression)...Proceeding without differential expression analsyes..."
            pass

    gene_list = None
    if species != None:
        gene_list = cluster_corr.read_gene_list(customMarkers)
        
    export_directory = os.path.abspath(os.path.join(query_exp_file, os.pardir))
    dataset_name = string.replace(string.split(query_exp_file,'/')[-1][:-4],'exp.','')
    try: os.mkdir(export_directory+'/CellClassification/')
    except: pass
    output_classification_file = export_directory+'/CellClassification/'+dataset_name+'-CellClassification.txt'
    reference = customMarkers
    query = query_exp_file
    louvain_results = cluster_corr.find_nearest_cells(reference,
                    query_exp_file,
                    gene_list=gene_list,
                    num_neighbors=10,
                    num_trees=100,
                    louvain_level=0,
                    min_cluster_correlation=-1,
                    genome=species)
    cluster_corr.write_results_to_file(louvain_results, output_classification_file, labels=customLabels)
    
    LineageProfilerIterate.harmonizeClassifiedSamples(species, reference, query_exp_file, output_classification_file,fl=fl)
    
    return True

if __name__ == '__main__':

    from argparse import ArgumentParser
    parser = ArgumentParser(description="find the cells in reference_h5 that are most similar to the cells in query_h5")
    parser.add_argument("reference_h5", help="a CellRanger h5 file")
    parser.add_argument("query_h5", help="a CellRanger h5 file")
    parser.add_argument("output", help="the result file to write")
    parser.add_argument("-g", "--genes", default=None, help="an ICGS file with the genes to use")
    parser.add_argument("-s", "--genome", default=None, help="genome aligned to")
    parser.add_argument("-k", "--num_neighbors", type=int, default=10,
                        help="number of nearest neighbors to use in clustering, default: %(default)s")
    parser.add_argument("-t", "--num_trees", type=int, default=100,
                        help="number of trees to use in random forest for approximating nearest neighbors, default: %(default)s")
    parser.add_argument("-l", "--louvain", type=int, default=0,
                        help="what level to cut the clustering dendrogram.  0 is the most granular, -1 the least.  Default: %(default)s")
    parser.add_argument("-m", "--min_correlation", type=float, default=-1,
                        help="the lowest correlation permissible between clusters.  Any clusters in query that don't correlate to ref at least this well will be skipped.  Default: %(default)s")
    parser.add_argument("-d", "--diff_expression", type=float, default=True,
                        help="perform differential expression analyses.  Default: %(default)s")
    parser.add_argument("-b", "--labels", type=str, default=None, help = "a tab-delimited text file with two columns (reference cell barcode and cluster name)")

    args = parser.parse_args()

    genome = None
    if args.genes != None:
        genome = args.genome

    if args.labels != None:
        labels = cluster_corr.read_labels_dictionary(args.labels)

    platform = 'RNASeq'
    manage_louvain_alignment(genome,platform,args.query_h5,None,
            customMarkers=args.reference_h5,useMulti=False,fl=None,customLabels=labels)
    