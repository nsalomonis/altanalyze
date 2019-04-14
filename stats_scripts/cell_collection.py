from __future__ import print_function

import h5py
import scipy.io
from scipy import sparse, stats
import numpy as np
import sys, string, os, csv, math
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

def index_items(universe, itemset):
    """
    Returns a list of indices to the items in universe that match items in itemset
    """
    return [ idx for idx, item in enumerate(universe) if item in itemset ]

class CellCollection:
    """
    Encapsulates a cohort of cells, ie from a CellRanger run

    Expression values are stored in a sparse matrix, and barcodes/gene identifiers are
    maintained in parallel arrays.  Construct by calling CellCollection.from_cellranger_h5()
    """

    @staticmethod
    def from_cellranger_h5(h5_filename, genome=None, returnGenes=False):
        """
        Creates a CellCollection from the contents of an H5 file created by CellRanger.

        If genome is None, the only genome in the H5 will be loaded.  If genome is None and multiple genomes
        are in the H5, an Exception is raised.
        """
        coll = CellCollection()
        f = h5py.File(h5_filename, 'r')
        if genome == None:
            possible_genomes = f.keys()
            if len(possible_genomes) != 1:
                raise Exception("{} contains multiple genomes ({}).  Explicitly select one".format(h5_filename, ", ".join(possible_genomes)))
            genome = possible_genomes[0]
            #print("Auto-selecting genome {}".format(genome), file=sys.stderr)

        coll._gene_names = f[genome]['gene_names']
        coll._gene_ids = f[genome]['genes']
        coll._barcodes = f[genome]['barcodes']
        coll._barcodes = map(lambda x: string.replace(x,'-1',''), coll._barcodes)
        
        if returnGenes:
            """ Do not import the matrix at this point """
            return list(coll._gene_names)
        
        coll._matrix = sparse.csc_matrix((f[genome]['data'], f[genome]['indices'], f[genome]['indptr']))

        print('sparse matrix data imported from h5 file...')
        return coll

    @staticmethod
    def from_cellranger_mtx(matrix_directory, genome=None, returnGenes=False):
        """
        Creates a CellCollection from the contents of an mtx directory created by CellRanger.
        """
        
        coll = CellCollection()

        genes_path = os.path.join(matrix_directory, "genes.tsv")
        if os.path.isfile(genes_path)==False:
            genes_path = os.path.join(matrix_directory, "features.tsv")
        gene_ids = [row[0] for row in csv.reader(open(genes_path), delimiter="\t")]

        try:
            coll._gene_names = np.array([row[1] for row in csv.reader(open(genes_path), delimiter="\t")])
        except:
            coll._gene_names = np.array([row[0] for row in csv.reader(open(genes_path), delimiter="\t")])
        if returnGenes:
            """ Do not import the matrix at this point """
            return list(coll._gene_names)
        
        coll._gene_ids = coll._gene_names
        barcodes_path = os.path.join(matrix_directory, "barcodes.tsv")
        barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]
        coll._barcodes = np.array(map(lambda x: string.replace(x,'-1',''), barcodes))        
        coll._matrix = scipy.io.mmread(os.path.join(matrix_directory, "matrix.mtx"))
        
        print('sparse matrix data imported from mtx file...')
        return coll
    
    @staticmethod
    def from_tsvfile(tsv_file, genome=None, returnGenes=False, gene_list=None):
        """
        Creates a CellCollection from the contents of a tab-separated text file.
        """
        
        coll = CellCollection()
        
        UseDense=False
        header=True
        skip=False
        for line in open(tsv_file,'rU').xreadlines():
            if header:
                delimiter = ',' # CSV file
                start = 1
                if 'row_clusters' in line:
                    start=2 # An extra column and row are present from the ICGS file
                    skip=True
                if '\t' in line:
                    delimiter = '\t' # TSV file

                coll._barcodes=string.split(line.rstrip(),delimiter)[start:]
                coll._gene_names=[]
                data_array=[]
                header=False
            elif skip:
                skip=False # Igore the second row in the file that has cluster info
            else:
                values = line.rstrip().split(delimiter)
                gene = values[0]
                if gene_list!=None:
                    if gene not in gene_list:
                        continue
                if ' ' in gene:
                    gene = string.split(gene,' ')[0]
                if ':' in gene:
                    coll._gene_names.append((gene.rstrip().split(':'))[1])
                else:
                    coll._gene_names.append(gene)

                """ If the data (always log2) is a float, increment by 0.5 to round up """
                if returnGenes==False:
                    if UseDense:
                        data_array.append(map(float,values[start:]))
                    else:
                        #data_array.append(map(lambda x: round(math.pow(2,float(x))),values[start:]))
                        data_array.append(map(float,values[start:]))
                    
        if returnGenes:
            """ Do not import the matrix at this point """
            return list(coll._gene_names)
        
        if UseDense:
            coll._matrix = np.array(data_array)
        else:
            """ Convert to a sparse matrix """
            coll._matrix = sparse.csc_matrix(np.array(data_array))

        coll._barcodes = np.array(coll._barcodes)
        coll._gene_names = np.array(coll._gene_names)
        coll._gene_ids = coll._gene_names
        print('sparse matrix data imported from TSV file...')
        return coll
    
    def __init__(self):
        self._matrix = sparse.csc_matrix((0,0), dtype=np.int8)
        self._barcodes = ()
        self._gene_names = ()
        self._gene_ids = ()
    
    def __getattr__(self, name):
        """
        Methods/attributes not explicitly defined in the CellCollection are passed down
        to the matrix
        """
        return getattr(self._matrix, name)

    def num_genes(self):
        return len(self._gene_ids)

    def num_cells(self):
        return len(self._barcodes)

    def get_barcode(self, cell_index):
        return self._barcodes[cell_index]
    
    def get_cell_expression_vector(self, cell_index):
        """
        Returns a (standard, non-sparse) sequence of expression values for a given cell
        """
        try:
            return self._matrix.getcol(cell_index).todense()
        except:
            return self._matrix[:,cell_index] # ith column for existing dense matrix

    def centroid(self):
        """
        Returns the centroid of this collection as a (standard, non-sparse) sequence.

        The centroid is defined as the mean expression of each gene
        """
        return self._matrix.mean(axis=1)
    
    def partition(self, partition_dict):
        """
        Returns a dictionary of CellCollections, each a distinct subset (by cell) of self.
        partition_dict is a dictionary of cell index => set id, as generated by 
        the python-louvain methods
        """
        partitions = {}
        for k, v in partition_dict.items():
            if v not in partitions: partitions[v] = []
            partitions[v].append(k)
        result = {}
        for part_id in partitions.keys():
            result[part_id] = self.subset_by_cell_index(partitions[part_id])
        return result

    def find_best_correlated(self, query):
        """
        Identifies the cell in this collection that has the highest Pearson's correlation
        with query (a sequence of expression values in the same order as in this collection)

        Returns the pair of (barcode, r^2 value) for the best match in ref
        """
        best_cor = -2
        best_bc = "<None>"
        for idx in range(self.num_cells()):
            r = self.get_cell_expression_vector(idx)
            cor = stats.pearsonr(query, r)[0][0] # pearsonr returns the pair (r^2, p-val), and for some reason the r^2 is a list
            if cor > best_cor:
                best_cor = cor
                best_bc = self.get_barcode(idx)
        return best_bc, best_cor

    def filter_by_cell_index(self, cell_index):
        self._matrix = self._matrix[:, cell_index]
        self._barcodes = self._barcodes[cell_index]

    def subset_by_cell_index(self, cell_index):
        """
        Returns a new CellCollection containing only chosen cells from self
        """
        cc = CellCollection()
        cc._gene_ids = self._gene_ids
        cc._gene_names = self._gene_names
        cc._matrix = self._matrix[:, cell_index]
        cc._barcodes = self._barcodes[cell_index]
        return cc

    def filter_barcodes(self, barcode_list):
        """
        Reduces the CellCollection in-place to only contain the barcodes requested
        """
        barcode_subset = set(barcode_list)
        #print("Selecting {} barcodes".format(len(barcode_subset)), file=sys.stderr)
        barcode_index = index_items(self._barcodes, barcode_subset)
        self.filter_by_cell_index(barcode_index)
    
    def subset_barcodes(self, barcode_list):
        barcode_subset = set(barcode_list)
        barcode_index = index_items(self._barcodes, barcode_subset)
        return self.subset_by_cell_index(barcode_index)

    def _filter_genes_by_index(self, gene_index):
        #print(gene_index);sys.exit()
        self._matrix = self._matrix[gene_index, :]
        self._gene_ids = self._gene_ids[gene_index]
        self._gene_names = self._gene_names[gene_index]
        #mat_array_original = self._matrix.toarray()
        #print(len(mat_array_original))

    def filter_genes_by_symbol(self, symbol_list):
        """
        Reduces the CellCollection in-place to only contain the genes requested.

        Note that gene symbols could be non-unique, and thus more genes may remain in the
        filtered collection than were requested. The order of the genes in the h5 may also
        differ and the same genes may not be present in the different sets
        """
        gene_subset = set(symbol_list)
        #print("Selecting {} genes".format(len(gene_subset)), file=sys.stderr)
        #gene_index = index_items(self._gene_names, gene_subset) # will output genes in the full dataset order
        gene_index=[]
        gene_names = list(self._gene_names)
        for gene in gene_subset:
            if gene in gene_names:
                gene_index.append(gene_names.index(gene))
        self._filter_genes_by_index(gene_index)

    def filter_genes_by_id(self, id_list):
        """
        Reduces the CellCollection in-place to only contain the genes requested.
        """
        gene_subset = set(id_list)
        #print("Selecting {} genes".format(len(gene_subset)), file=sys.stderr)
        gene_index = index_items(self._gene_ids, gene_subset)
        self._filter_genes_by_index(gene_index)
    
    
