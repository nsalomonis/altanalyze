import os, sys, string
from scipy import sparse, io
import numpy
import math
import time

""" Converts a tab-delimited text file to a sparse matrix """

class SparseMatrix:
    def __init__(self,barcodes,features,data_matrix):
        self.barcodes = barcodes
        self.features = features
        self.data_matrix = data_matrix
    def Barcodes(self): return self.barcodes
    def Features(self): return self.features
    def Matrix(self): return self.data_matrix
    
def import_filter_genes(fn):
    gene_filter = []
    for line in open(fn,'rU').xreadlines():
        if '\t' in line:
            gene = string.split(line.rstrip(),'\t')[0]
        else:
            gene = string.split(line.rstrip(),',')[0]
        gene_filter.append(gene)
    return gene_filter

def covert_table_to_matrix(fn,delimiter='\t',gene_filter=None,Export=True):
    header=True
    skip=False
    start_time = time.time()
    for line in open(fn,'rU').xreadlines():
        if header:
            delimiter = ',' # CSV file
            start = 1
            if 'row_clusters' in line:
                start=2 # An extra column and row are present from the ICGS file
                skip=True
            if '\t' in line:
                delimiter = '\t' # TSV file
                
            t = string.split(line.rstrip(),delimiter)[start:]
            
            """ Optionally write out the matrix """
            if Export:
                export_directory = os.path.abspath(os.path.join(fn, os.pardir))+'/sparse/'
                try: os.mkdir(export_directory)
                except: pass
                barcodes_export = open(export_directory+'barcodes.tsv', 'w')
                features_export = open(export_directory+'features.tsv', 'w')
                matrix_export = open(export_directory+'matrix.mtx', 'w')
                barcodes = t
                barcodes_export.write(string.join(t,'\n'))
                barcodes_export.close()
            genes=[]
            data_array=[]
            header=False
        elif skip:
            skip=False # Igore the second row in the file that has cluster info
        else:
            values = string.split(line.rstrip(),'\t')
            gene = values[0]
            if gene_filter!=None:
                """ Exclude the gene from the large input matrix if not in the filter list """
                if gene not in gene_filter:
                    continue
            if ' ' in gene:
                gene = string.split(gene,' ')[0]
            if ':' in gene:
                genes.append((gene.rstrip().split(':'))[1])
            else:
                genes.append(gene)
            """ If the data is a float, increment by 0.5 to round up """
            values = map(float,values[start:])
            def convert(x):
                if x==0:
                    return 0
                else:
                    return int(math.pow(2,x)-1)
            values = map(lambda x: convert(x),values) 
            data_array.append(values)

    data_array = sparse.csr_matrix(numpy.array(data_array))

    end_time = time.time()
    print 'Sparse matrix conversion in',end_time-start_time,'seconds.'
    
    if Export:
        features_export.write(string.join(genes,'\n'))
        features_export.close()
        io.mmwrite(matrix_export, data_array)
    else:
        sm = SparseMatrix(barcodes,genes,data_array)
        return sm

if __name__ == '__main__':
    import getopt
    
    gene_filter = None
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','f='])
        for opt, arg in options:
            if opt == '--i': fn=arg
            elif opt == '--f': gene_filter=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    if gene_filter != None:
         gene_filter = import_filter_genes(gene_filter)
    covert_table_to_matrix(fn,gene_filter=gene_filter)