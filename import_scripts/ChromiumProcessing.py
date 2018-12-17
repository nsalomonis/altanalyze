import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import csv
import scipy.io
import numpy
import time
import math

def import10XSparseMatrix(matrices_dir,genome,dataset_name, expFile=None, log=True):
    start_time = time.time()
    human_matrix_dir = os.path.join(matrices_dir, genome)
    mat = scipy.io.mmread(os.path.join(human_matrix_dir, "matrix.mtx"))
    genes_path = os.path.join(human_matrix_dir, "genes.tsv")
    gene_ids = [row[0] for row in csv.reader(open(genes_path), delimiter="\t")]
    gene_names = [row[1] for row in csv.reader(open(genes_path), delimiter="\t")]
    barcodes_path = os.path.join(human_matrix_dir, "barcodes.tsv")
    barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]
    barcodes = map(lambda x: string.replace(x,'-1',''), barcodes)

    ### Write out raw data matrix
    counts_path = matrices_dir+'/'+dataset_name+'_matrix.txt'
    if expFile!=None:
        if 'exp.' in expFile:
            counts_path = string.replace(expFile,'exp.','counts.')
            
    ### Efficiently write the data to an external file (fastest way)
    mat_array_original = mat.toarray() ### convert sparse matrix to numpy array
    mat_array = numpy.ndarray.tolist(mat_array_original) ### convert to non-numpy list
    i=0
    updated_mat=[['UID']+barcodes]
    for ls in mat_array:
        updated_mat.append([gene_names[i]]+ls); i+=1
    mat_array = numpy.array(mat_array)
    updated_mat = numpy.array(updated_mat)
    
    numpy.savetxt(counts_path,updated_mat,fmt='%s',delimiter='\t')
    del updated_mat
    print 'Raw-counts written to file:'
    print counts_path 
    #print time.time()-start_time
    
    ### Write out CPM normalized data matrix
    if expFile==None:
        norm_path = matrices_dir+'/'+dataset_name+'_matrix_CPTT.txt'
    else:
        norm_path = expFile ### AltAnalyze designated ExpressionInput file (not counts)

    print 'Normalizing gene counts to counts per ten thousand (CPTT)'
    
    barcode_sum = numpy.sum(mat_array,axis=0) ### Get the sum counts for all barcodes, 0 is the y-axis in the matrix
    
    start_time = time.time()
    mat_array = mat_array.transpose()
    def calculateCPTT(val,barcode_sum):
        if val==0:
            return '0'
        else:
            if log:
                return math.log((10000.00*val/barcode_sum)+1.0,2) ### convert to log2 expression
            else:
                return 10000.00*val/barcode_sum

    vfunc = numpy.vectorize(calculateCPTT)
    norm_mat_array=[]
    i=0
    for vector in mat_array:
        norm_mat_array.append(vfunc(vector,barcode_sum[i]))
        i+=1
    mat_array = numpy.array(norm_mat_array)
    mat_array = mat_array.transpose()
    
    mat_array = numpy.ndarray.tolist(mat_array) ### convert to non-numpy list
    i=0
    updated_mat=[['UID']+barcodes]
    for ls in mat_array:
        updated_mat.append([gene_names[i]]+ls); i+=1
    updated_mat = numpy.array(updated_mat)
    del mat_array
    numpy.savetxt(norm_path,updated_mat,fmt='%s',delimiter='\t')
    
    #print time.time()-start_time
    print 'CPTT written to file:',
    print norm_path

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
    import10XSparseMatrix(matrices_dir,genome,dataset_name)