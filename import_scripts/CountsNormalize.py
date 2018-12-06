import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import csv
import scipy.io
import numpy
import time
import math

def normalizeDropSeqCounts(expFile,log=True):
    start_time = time.time()
    
    print 'log2 conversion equals',log
    firstLine = True
    mat_array=[]
    gene_names = []
    for line in open(expFile,'rU').xreadlines():
        line = string.replace(line,'"','')
        data = line.rstrip('\n')
        t = string.split(data,'\t')
        if firstLine:
            barcodes = t[1:]
            firstLine = False
        else:
            gene_names.append(t[0])
            mat_array.append(map(float,t[1:]))
    mat_array = numpy.array(mat_array)
    
    ### Write out CPM normalized data matrix
    norm_path = expFile[:-4]+'_matrix_CPTT.txt' ### AltAnalyze designated ExpressionInput file (not counts)

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
    log=True
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Insufficient options provided";sys.exit()
        #Filtering samples in a datasets
        #python DropSeqProcessing.py --i dropseq.txt
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','log='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i': matrices_dir=arg
            if opt == '--log':
                if string.lower(arg) == 'true' or string.lower(arg) == 'yes':
                    pass
                else:
                    log = False
    normalizeDropSeqCounts(matrices_dir,log=log)