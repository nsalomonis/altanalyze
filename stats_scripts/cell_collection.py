from __future__ import print_function

try:
    import h5py
    from h5py import defs, utils, h5ac, _proxy # for py2app
except:
    print ('Missing the h5py library (hdf5 support)...')

import gzip
import scipy.io
from scipy import sparse, stats, io
import numpy as np
import sys, string, os, csv, math
import time
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
    maintained in parallel arrays.  Construct by calling CellCollection.from_file(), or one
    of the other specialized static constructors
    """

    @staticmethod
    def from_cellranger_h5(h5_filename, genome=None, returnGenes=False):
        """
        Creates a CellCollection from the contents of an H5 file created by CellRanger.

        The meaning of the genome parameter differs depending on the version of CellRanger that created the h5.
        
        For CellRanger version 2, the genome parameters specifies the matrix to load.  If genome is None, the
        single matrix present will be loaded (using genome==None when multiple genomes are present in the file
        is an error and will cause an exception).

        For CellRanger version 3, genome is now specified as an attribute of the features (typically genes).
        In this version, specifying a genome will filter the matrix to only include features from that genome.
        Whether a genome is specified or not, non-gene features will be removed
        """
        start = time.time()
        coll = CellCollection()
        f = h5py.File(h5_filename, 'r')
        if 'matrix' in f:
            # CellRanger v3
            coll._barcodes = f['matrix']['barcodes']
            coll._gene_ids = f['matrix']['features']['id']
            coll._gene_names = f['matrix']['features']['name']
            
            if returnGenes:
                """ Do not import the matrix at this point """
                return list(coll._gene_names)
            
            coll._matrix = sparse.csc_matrix((f['matrix']['data'], f['matrix']['indices'], f['matrix']['indptr']), shape=f['matrix']['shape'])
            indices = np.flatnonzero(np.array(f['matrix']['features']['genome']) != '') if \
                genome == None else \
                np.flatnonzero(np.array(f['matrix']['features']['genome']) == genome)
            coll._filter_genes_by_index(indices.tolist())
        else:
            # CellRanger v2
            if genome == None:
                possible_genomes = f.keys()
                if len(possible_genomes) != 1:
                    raise Exception("{} contains multiple genomes ({}).  Explicitly select one".format(h5_filename, ", ".join(possible_genomes)))
                genome = possible_genomes[0]
                #print("Auto-selecting genome {}".format(genome), file=sys.stderr)

            coll._gene_names = f[genome]['gene_names']
            if returnGenes:
                """ Do not import the matrix at this point """
                return list(coll._gene_names)

            coll._matrix = sparse.csc_matrix((f[genome]['data'], f[genome]['indices'], f[genome]['indptr']))   
            coll._barcodes = f[genome]['barcodes']
            coll._gene_ids = f[genome]['genes']

        print('sparse matrix data imported from h5 file in %s seconds' % str(time.time()-start))
        return coll

    @staticmethod
    def from_cellranger_mtx(mtx_directory, genome=None, returnGenes=False):

        """
        Creates a CellCollection from a sparse matrix (.mtx and associated files) exported by CellRanger

        Recognize directories from CellRanger version 2 (files: matrix.mtx, genes.tsv, barcodes.tsv) and
        CellRanger v3 (files: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)
        """

        start = time.time()
        coll = CellCollection()
        cellranger_version = 2
        if '.mtx' in mtx_directory:
            mtx_file = mtx_directory ### Hence an mtx file was directly supplied
            mtx_directory = os.path.abspath(os.path.join(mtx_file, os.pardir))
        else:
            mtx_file = os.path.join(mtx_directory, "matrix.mtx")
            
        if not os.path.exists(mtx_file):
            cellranger_version = 3
            mtx_file = mtx_file + ".gz"
            if not os.path.exists(mtx_file):
                raise Exception("Directory {} does not contain a recognizable matrix file".format(mtx_directory))
        if '.gz' in mtx_file:
            cellranger_version = 3
        sparse_matrix = io.mmread(mtx_file)
        coll._matrix = sparse_matrix.tocsc()
        coll._gene_ids = np.empty((coll._matrix.shape[0], ), np.object)
        coll._gene_names = np.empty((coll._matrix.shape[0], ), np.object)
        
        if cellranger_version == 2:
            with open(os.path.join(mtx_directory, "genes.tsv"), "rU") as f:
                idx = 0
                for line in f:
                    i, n = line.rstrip().split("\t")
                    coll._gene_ids[idx] = i
                    coll._gene_names[idx] = n
                    idx += 1
            with open(os.path.join(mtx_directory, "barcodes.tsv"), "rU") as f:
                coll._barcodes = np.array( [ line.rstrip() for line in f ] )
        else:
            with gzip.open(os.path.join(mtx_directory, "features.tsv.gz"), "rt") as f:
                idx = 0
                indices = []
                for line in f:
                    i, n, t = line.rstrip().split("\t")
                    coll._gene_ids[idx] = i
                    coll._gene_names[idx] = n
                    if t == 'Gene Expression':
                        indices.append(idx)
                    idx += 1
                coll._filter_genes_by_index(indices)
            with gzip.open(os.path.join(mtx_directory, "barcodes.tsv.gz"), "rt") as f:
                coll._barcodes = np.array( [ line.rstrip() for line in f ] )

        if returnGenes:
            """ Do not import the matrix at this point """
            return list(coll._gene_names)
            
        print('sparse matrix data imported from mtx file in %s seconds' % str(time.time()-start))
        return coll

    @staticmethod
    def from_tsvfile_alt(tsv_file, genome=None, returnGenes=False, gene_list=None):
        """
        Creates a CellCollection from the contents of a tab-separated text file.
        """
        startT = time.time()
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
                barcodes = string.split(line.rstrip(),delimiter)[start:]
                if ':' in line:
                    barcodes = map(lambda x:x.split(':')[1],barcodes)

                coll._barcodes=barcodes
                coll._gene_names=[]
                data_array=[]
                header=False
            elif skip:
                skip=False # Igore the second row in the file that has cluster info
            else:
                values = line.rstrip().split(delimiter)
                gene = values[0]
                if ' ' in gene:
                    gene = string.split(gene,' ')[0]
                if ':' in gene:
                    gene = (gene.rstrip().split(':'))[1]
                if gene_list!=None:
                    if gene not in gene_list:
                        continue
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

        print('sparse matrix data imported from TSV file in %s seconds' % str(time.time()-startT))
        #print (len(coll._gene_ids),len(coll._barcodes))
        return coll
    
    @staticmethod
    def from_tsvfile(tsv_filename, genome=None, returnGenes=False, gene_list=None):
        """
        Generates a CellCollection from a (dense) tab-separated file, where cells are in
        columns and
        """
        
        start = time.time()
        coll = CellCollection()
        with open(tsv_filename, "rU") as f:
            try:
                line = next(f)
            except StopIteration:
                raise Exception("TSV file {} is empty".format(tsv_filename))
            
            ### Check formatting
            skip=False
            if '\t' in line:
                delimiter = '\t' # TSV file
            else:
                delimiter = ','
            col_start = 1
            if 'row_clusters' in line:
                col_start=2 # An extra column and row are present from the ICGS file
                skip=True
            ### Check formatting end
                    
            coll._barcodes = np.array(line.rstrip().split(delimiter)[col_start:])
            sparse_matrix = sparse.lil_matrix((50000, len(coll._barcodes)), dtype=np.float_)
            coll._gene_names = np.empty((sparse_matrix.shape[0], ), np.object)
            row = 0
            for line in f:
                if row==0 and skip:
                    skip = False
                    continue
                vals = line.rstrip().split(delimiter)
                coll._gene_names[row] = vals[0]
                if returnGenes==False:
                    for i in range(col_start, len(vals)):
                        if vals[i] != "0":
                            sparse_matrix[row, i-col_start] = float(vals[i])
                    if row == sparse_matrix.shape[0]-1:
                        sparse_matrix.resize(sparse_matrix.shape + (10000, 0))
                        coll._gene_names.resize(coll._gene_names.shape + (10000, 0))
                row += 1

        coll._gene_names.resize((row, ))
        if returnGenes:
            """ Do not import the matrix at this point """
            return list(coll._gene_names)
        
        sparse_matrix.resize((row, len(coll._barcodes)))
        coll._matrix = sparse_matrix.tocsc()
        coll._gene_ids = coll._gene_names
    
        #print('matrix shape: {}'.format(coll._matrix.shape))
        print('sparse matrix data imported from TSV file in %s seconds' % str(time.time()-start))
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
        #try:
        return self._matrix.getcol(cell_index).todense()
        #except:
        #    return self._matrix[:,cell_index] # ith column for existing dense matrix

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

    def filter_genes_by_symbol(self, symbol_list, data_type):
        """
        Reduces the CellCollection in-place to only contain the genes requested.

        Note that gene symbols could be non-unique, and thus more genes may remain in the
        filtered collection than were requested. The order of the genes in the h5 may also
        differ and the same genes may not be present in the different sets
        """
        gene_subset = set(symbol_list)
        #print("Selecting {} genes".format(len(gene_subset)), file=sys.stderr)
        gene_index=[]
        gene_names = list(self._gene_names)
        if data_type == 'txt':
            ### below code is problematic for h5 and probably sparse matrix files
            for gene in gene_subset:
                if gene in gene_names:
                    gene_index.append(gene_names.index(gene))
        else:
            gene_index = index_items(self._gene_names, gene_subset) # will output genes in the full dataset order
            
        self._filter_genes_by_index(gene_index)

    def filter_genes_by_id(self, id_list):
        """
        Reduces the CellCollection in-place to only contain the genes requested.
        """
        gene_subset = set(id_list)
        #print("Selecting {} genes".format(len(gene_subset)), file=sys.stderr)
        gene_index = index_items(self._gene_ids, gene_subset)
        self._filter_genes_by_index(gene_index)
    
    
