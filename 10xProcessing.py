import sys,string
import csv
import os
import scipy.io

def import10XSparseMatrix(matrices_dir,genome,dataset_name):
    human_matrix_dir = os.path.join(matrices_dir, genome)
    mat = scipy.io.mmread(os.path.join(human_matrix_dir, "matrix.mtx"))
    genes_path = os.path.join(human_matrix_dir, "genes.tsv")
    gene_ids = [row[0] for row in csv.reader(open(genes_path), delimiter="\t")]
    gene_names = [row[1] for row in csv.reader(open(genes_path), delimiter="\t")]
    barcodes_path = os.path.join(human_matrix_dir, "barcodes.tsv")
    barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]
        
    ### Write out raw data matrix
    counts_path = matrices_dir+'/'+dataset_name+'_matrix.txt'
    outfile = open(counts_path, 'w')
    outfile.write(string.join(['UID']+barcodes,'\t')+'\n')
    mat_array = mat.toarray()
    
    i=0
    for k in mat_array:
        gene = gene_names[i]
        values = map(str,mat_array[i])
        outfile.write(string.join([gene]+values,'\t')+'\n')
        i+=1
    outfile.close()
    print 'Raw-counts written to file:'
    print counts_path
    ### Write out CPM normalized data matrix
    norm_path = matrices_dir+'/'+dataset_name+'_matrix_CPTT.txt'
    outfile = open(norm_path, 'w')
    outfile.write(string.join(['UID']+barcodes,'\t')+'\n')
    
    mat_array_t = mat_array.transpose()

    print 'Normalizing gene counts to counts per ten thousand (CPTT)'
    barcode_sum=[]
    for k in mat_array_t:
        barcode_sum.append(sum(k))
    
    i=0
    for k in mat_array:
        gene = gene_names[i]
        l=0; cpms=[]
        for x in k:
            if x!= 0:
                value = (float(x)/barcode_sum[l])*10000
                cpms.append(str(value))
                l+=1
            else:
                cpms.append('0')
        outfile.write(string.join([gene]+cpms,'\t')+'\n')
        i+=1
    
    outfile.close()
    print 'CPTT written to file:'
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