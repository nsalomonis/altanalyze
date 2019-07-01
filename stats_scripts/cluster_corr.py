from __future__ import print_function

from cell_collection import CellCollection
from annoy import AnnoyIndex
import community # python-louvain package, imported as community
import networkx
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) # import parent dir dependencies
import numpy as np
import time

def read_gene_list(filename):
    """
    Reads the gene list from a file 
    """
    gene_list = []
    with open(filename, 'rU') as f:
        for line in f:
            gene = line.split('\t', 1)[0]
            if ' ' in gene:
                gene = string.split(gene.rstrip(),' ')[0]
            if ':' in gene:
                gene_list.append((gene.rstrip().split(':'))[1])
            else:
                gene_list.append(gene.rstrip())
    return gene_list

def read_labels_dictionary(filename):
    """
    Reads the labels assigned by the user to each barcode (first and last columns)
    """
    label_dictionary = {}
    with open(filename, 'rU') as f:
        header_lines = 2
        for line in f:
            barcode = line.split('\t', 1)[0]
            label = line.split('\t', 1)[-1]
            label_dictionary[barcode]=label.rstrip()
    return label_dictionary

def data_check(ref_h5_filename, query_h5_filename):
    """ If h5 files, are both h5? """
    if 'h5' in ref_h5_filename and 'h5' not in query_h5_filename:
        return False
    else:
        return True
    
def find_nearest_cells(ref_h5_filename, query_h5_filename, gene_list=None, genome=None,
                        num_neighbors=10, num_trees=100, louvain_level=-1, min_cluster_correlation=-1):
    """
    For every cell in query_h5_filename, identifies the most similar cell in ref_h5_filename.

    For parameter definitions, see partition_h5_file() and find_closest_cluster(), this is a convenience
    function that calls them (among others)
    """

    ### Do the two input file formats match?
    matching = data_check(ref_h5_filename, query_h5_filename)
    
    startT = time.time()
    if '.mtx' not in ref_h5_filename and '.h5' not in ref_h5_filename:
        """ Assumes partial overlapping gene lists present """
        gene_list = find_shared_genes(ref_h5_filename,genome=genome,gene_list=gene_list)
        gene_list = find_shared_genes(query_h5_filename,genome=genome,gene_list=gene_list)

    ref_partition = partition_h5_file(ref_h5_filename, gene_list=gene_list, num_neighbors=num_neighbors, 
                            num_trees=num_trees, louvain_level=louvain_level,genome=genome)
    query_partition = partition_h5_file(query_h5_filename, gene_list=gene_list, num_neighbors=num_neighbors, 
                            num_trees=num_trees, louvain_level=louvain_level,genome=genome)
    best_match = find_closest_cluster(query_partition, ref_partition, min_correlation=min_cluster_correlation)
    result = {}
    for query_part_id, ref_part_id in best_match:
        ref = ref_partition[ref_part_id]
        query = query_partition[query_part_id]
        for idx in range(query.num_cells()):
            q_barcode = query.get_barcode(idx)
            best_bc, best_cor = ref.find_best_correlated(query.get_cell_expression_vector(idx))
            result[q_barcode] = {'barcode': best_bc, 
                                 'correlation': best_cor, 
                                 'query_partition': query_part_id,
                                 'ref_partition': ref_part_id}
            
    print('cellHarmony-community alignment complete in %s seconds' % str(time.time()-startT))
    return result

def write_results_to_file(results, filename, labels=None):
    def add_labels(barcode):
        if barcode in labels:
            return labels[barcode]
        if ':' in barcode:
            return labels[barcode.split(':')[1]]
        else:
            return 'NA'
            
    if labels == None:
        with open(filename, 'w') as f:
            print("\t".join( ("Query Barcode", "Ref Barcode", "Correlation", "Query Partition", "Ref Partition") ), file=f)
            for q in results.keys():
                print("\t".join( (q, 
                    results[q]['barcode'], 
                    str(results[q]['correlation']), 
                    str(results[q]['query_partition']),
                    str(results[q]['ref_partition'])) ), file=f)
    else:
        with open(filename, 'w') as f:
            print("\t".join( ("Query Barcode", "Ref Barcode", "Correlation", "Query Partition", "Ref Partition", "Label") ), file=f)
            for q in results.keys():
                print("\t".join( (q, 
                    results[q]['barcode'], 
                    str(results[q]['correlation']), 
                    str(results[q]['query_partition']),
                    str(results[q]['ref_partition']),
                    add_labels(results[q]['barcode'])) ), file=f)

def nearest_neighbors(collection, num_neighbors=10, n_trees=100):
    """
    Finds the num_neighbors nearest neighbors to each cell in the sparse matrix

    Return result is a dictionary of lists, where the key is an index into the cells, 
    and the value is the neighbors of that cell
    """
    nn_idx = AnnoyIndex(collection.num_genes())
    # Add the elements in reverse order because Annoy allocates the memory based on
    # the value of the element added - so adding in increasing order will trigger
    # lots of allocations
    for i in range(collection.num_cells()-1, -1, -1):
        nn_idx.add_item(i, collection.get_cell_expression_vector(i))
    nn_idx.build(n_trees)
    return { i: nn_idx.get_nns_by_item(i, num_neighbors) for i in range(collection.num_cells()) }

def identify_clusters(graph, louvain_level=-1):
    """
    Identifies clusters in the given NetworkX Graph by Louvain partitioning.
    
    The parameter louvain_level controls the degree of partitioning.  0 is the most granular
    partition, and granularity decreases as louvain_level increases.  Since the number of
    levels can't be known a priori, negative values "count down" from the max - ie, -1
    means to use the maximum possible value and thus get the largest clusters
    """
    dendrogram = community.generate_dendrogram(graph)
    if louvain_level < 0:
        louvain_level = max(0, len(dendrogram) + louvain_level)
    if louvain_level >= len(dendrogram):
        #print("Warning [identify_clusters]: louvain_level set to {}, max allowable is {}.  Resetting".format(louvain_level, len(dendrogram)-1), file=sys.stderr)
        louvain_level = len(dendrogram) - 1
    #print("Cutting the Louvain dendrogram at level {}".format(louvain_level), file=sys.stderr)
    return community.partition_at_level(dendrogram, louvain_level)

def find_shared_genes(h5_filename,genome=None,gene_list=None):
    """
    Selects genes shared by the reference, query and gene_list
    for filtering genes
    """    
    
    if gene_list !=None:
        if 'h5' in h5_filename:
            genes = CellCollection.from_cellranger_h5(h5_filename,returnGenes=True)
        elif 'txt' in h5_filename:
            try:
                genes = CellCollection.from_tsvfile_alt(h5_filename,genome,returnGenes=True,gene_list=gene_list)
            except:
                genes = CellCollection.from_tsvfile(h5_filename,genome,returnGenes=True,gene_list=gene_list)
        else:
            genes = CellCollection.from_cellranger_mtx(h5_filename,genome,returnGenes=True)
        gene_list = list(set(genes) & set(gene_list))
        
    return gene_list

def partition_h5_file(h5_filename, gene_list=None, num_neighbors=10, num_trees=100,
                    louvain_level=-1,genome=None):
    """
    Reads a CellRanger h5 file and partitions it by clustering on the k-nearest neighbor graph

    Keyword arguments:
    gene_list - restricts the analysis to the specified list of gene symbols.  Default is to not restrict
    num_neighbors - the number of nearest neighbors to compute for each cell
    num_trees - the number of trees used in the random forest that approximates the nearest neighbor calculation
    louvain_level - the level of the Louvain clustering dendrogram to cut at.  Level 0 is the lowest (most granular) 
                    level, and higher levels get less granular.  The highest level is considered the "best" set
                    of clusters, but the number of levels is not known a priori.  Hence, negative values will
                    count down from the highest level, so -1 will always be the "best" clustering, regardless of
                    the actual number of levels in the dendrogram

    Return Result: A dictionary, where the keys are partition ids and the values are the CellCollection for that partition
    """
    if 'h5' in h5_filename:
        collection = CellCollection.from_cellranger_h5(h5_filename)
        data_type = 'h5'
    elif 'txt' in h5_filename:
        try:
            collection = CellCollection.from_tsvfile_alt(h5_filename,genome,gene_list=gene_list)
        except:
            collection = CellCollection.from_tsvfile(h5_filename,genome)
        data_type = 'txt'
    else:
        collection = CellCollection.from_cellranger_mtx(h5_filename,genome)
        data_type = 'mtx'
    if gene_list != None:
        collection.filter_genes_by_symbol(gene_list,data_type)
    neighbor_dict = nearest_neighbors(collection, num_neighbors=num_neighbors, n_trees=num_trees)
    cluster_definition = identify_clusters(networkx.from_dict_of_lists(neighbor_dict), louvain_level=louvain_level)
    return collection.partition(cluster_definition)

def compute_centroids(collection_dict):
    """Returns (centroid matrix, ID list) for the given dictionary of CellCollections"""
    centroids = np.concatenate( [ p.centroid() for p in collection_dict.values() ], axis=1 )
    return centroids, collection_dict.keys()

def find_closest_cluster(query, ref, min_correlation=-1):
    """
    For each collection in query, identifies the collection in ref that is most similar

    query and ref are both dictionaries of CellCollections, keyed by a "partition id"

    Returns a list containing the best matches for each collection in query that meet the 
    min_correlation threshold.  Each member of the list is itself a list containing the 
    id of the query collection and the id of its best match in ref
    """
    query_centroids, query_ids = compute_centroids(query)
    ref_centroids, ref_ids = compute_centroids(ref)
    print('number of reference partions %d, number of query partions %d' % (len(ref_ids),len(query_ids)))
    all_correlations = np.corrcoef(np.concatenate((ref_centroids, query_centroids), axis=1), rowvar=False)

    # At this point, we have the correlations of everything vs everything.  We only care about query vs ref
    # Extract the top-right corner of the matrix
    nref = len(ref)
    corr = np.hsplit(np.vsplit(all_correlations, (nref, ))[0], (nref,))[1]
    best_match = zip(range(corr.shape[1]), np.argmax(corr, 0))
    # At this point, best_match is: 1) using indices into the array rather than ids, 
    # and 2) not restricted by the threshold.  Fix before returning
    return ( (query_ids[q], ref_ids[r]) for q, r in best_match if corr[r,q] >= min_correlation )

if __name__ == "__main__":
    # genes = read_gene_list('data/FinalMarkerHeatmap_all.txt')
    # results = find_nearest_cells('data/reference-filtered_gene_bc_matrices_h5.h5', 
    #                             'data/query-filtered_gene_bc_matrices_h5.h5',
    #                             gene_list=genes, louvain_level=-1)
    # write_results_to_file(results, 'temp.txt')

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
    parser.add_argument("-b", "--labels", type=str, default=None, help = "a tab-delimited text file with two columns (reference cell barcode and cluster name)")

    args = parser.parse_args()

    gene_list = None
    genome = None
    labels = None
    if args.genes != None:
        gene_list = read_gene_list(args.genes)
    if args.labels != None:
        labels = read_labels_dictionary(args.labels)
    if args.genome != None:
        genome = args.genome

    results = find_nearest_cells(args.reference_h5,
                                 args.query_h5,
                                 gene_list=gene_list,
                                 num_neighbors=args.num_neighbors,
                                 num_trees=args.num_trees,
                                 louvain_level=args.louvain,
                                 min_cluster_correlation=args.min_correlation,
                                 genome=genome)
    write_results_to_file(results, args.output,labels=labels)