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
    
    print '... scaling completed in ',time.time()-start_time, 'seconds'
    print 'CPTT written to file:',
    print norm_path

def calculateCPTT(val,barcode_sum,log=True):
        if val==0:
            return '0'
        else:
            if log:
                return math.log((10000.00*val/barcode_sum)+1.0,2) ### convert to log2 expression
            else:
                return 10000.00*val/barcode_sum

def normalizeDropSeqCountsMemoryEfficient(expFile,log=True):
    """ A more memory efficient function than the above for scaling scRNA-Seq, line-by-line """
    
    start_time = time.time()
    
    print 'log2 conversion equals',log
    firstLine = True
    mat_array=[]
    gene_names = []
    for line in open(expFile,'rU').xreadlines():
        line = string.replace(line,'"','')
        data = line.rstrip('\n')
        if '\t' in data:
            t = string.split(data,'\t')
        else:
            t = string.split(data,',')
        if firstLine:
            barcodes = t[1:]
            firstLine = False
            count_sum_array=[0]*len(barcodes)
        else:
            gene_names.append(t[0])
            values = map(float,t[1:])
            count_sum_array = [sum(value) for value in zip(*[count_sum_array,values])]

    ### Import the expression dataset again and now scale
    output_file = expFile[:-4]+'_CPTT-log2.txt'
    export_object = open(output_file,'w')
    firstLine=True
    for line in open(expFile,'rU').xreadlines():
        line = string.replace(line,'"','')
        data = line.rstrip('\n')
        if '\t' in data:
            t = string.split(data,'\t')
        else:
            t = string.split(data,',')
        if firstLine:
            export_object.write(string.join(t,'\t')+'\n')
            firstLine = False
        else:
            gene = t[0]
            if 'ENS' in gene and '.' in gene:
                gene = string.split(gene,'.')[0]
            values = map(float,t[1:])
            index=0
            cptt_values = []
            for barcode in barcodes:
                barcode_sum = count_sum_array[index]
                val = values[index]
                cptt_val = calculateCPTT(val,barcode_sum)
                cptt_values.append(cptt_val)
                index+=1
            values = string.join([gene]+map(lambda x: str(x)[:7], cptt_values),'\t')
            export_object.write(values+'\n')
    export_object.close()
    print '... scaling completed in ',time.time()-start_time, 'seconds'
    return output_file
    
def CSVformat(matrices_dir,cells_dir,emptydrops=True,target_organ=None,expressionCutoff=500):
    """ Process an HCA DCP csv dense matrix format """
    
    ### Import and process the cell metadata
    export_object = open(cells_dir[:-4]+'_'+target_organ+'-filtered.csv','w')
    cell_ids_to_retain={}
    firstLine=True
    for line in open(cells_dir,'rU').xreadlines():
        line = string.replace(line,'"','')
        data = line.rstrip('\n')
        t = string.split(data,',')
        if firstLine:
            index=0
            for i in t:
                if i == 'emptydrops_is_cell': edi = index
                elif i == 'genes_detected': gdi = index
                elif i == 'barcode': bi = index
                elif i == 'derived_organ_label': oi = index
                index+=1
            export_object.write(line)
            firstLine = False
        else:
            cell_id = t[0]
            empty_drops = t[edi]
            genes_detected = int(t[gdi])
            barcode = t[bi]
            organ = t[oi]
            proceed=True
            if emptydrops:
                if empty_drops == 'f':
                    proceed = False
            if target_organ != None:
                if target_organ != organ:
                    proceed = False
            if genes_detected<expressionCutoff:
                proceed = False
            if proceed:
                cell_ids_to_retain[cell_id]=barcode
                export_object.write(line)
    print len(cell_ids_to_retain), 'IDs matching the user filters'
    export_object.close()

    if target_organ != None:
        export_object = open(matrices_dir[:-4]+'_'+target_organ+'-filtered.txt','w')
    else:
        export_object = open(matrices_dir[:-4]+'-filtered.txt','w')
        
    ### Import and filter the flat expression data
    firstLine=True
    count=0
    print 'Increments of 10,000 cells exported:',
    for line in open(matrices_dir,'rU').xreadlines():
        line = string.replace(line,'"','')
        data = line.rstrip('\n')
        t = string.split(data,',')
        if firstLine:
            export_object.write(string.join(t,'\t')+'\n')
            firstLine = False
        else:
            cell_id = t[0]
            if cell_id in cell_ids_to_retain:
                export_object.write(string.join(t,'\t')+'\n')
                if count == 10000:
                    count = 0
                    print '*',
                count+=1
    export_object.close()
    
if __name__ == '__main__':
    import getopt
    log=True
    expressionCutoff = 500
    organ = None
    emptydrops = True
    cells_dir = None
    
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Insufficient options provided";sys.exit()
        #Filtering samples in a datasets
        #python DropSeqProcessing.py --i dropseq.txt
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','log=','csv=','organ=','expressionCutoff=',
                                                             'emptydrops=','cells='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i': matrices_dir=arg
            elif opt == '--csv': matrices_dir=arg
            elif opt == '--cells': cells_dir=arg
            elif opt == '--organ': target_organ=arg
            elif opt == '--expressionCutoff':
                expressionCutoff=float(arg)
            elif opt == '--emptydrops':
                if 'f' in arg or 'F' in arg:
                    emptydrops=False
            elif opt == '--log':
                if string.lower(arg) == 'true' or string.lower(arg) == 'yes':
                    pass
                else:
                    log = False
    if cells_dir != None:
        CSVformat(matrices_dir,cells_dir,emptydrops=emptydrops,target_organ=target_organ,expressionCutoff=expressionCutoff)
    else:
        normalizeDropSeqCountsMemoryEfficient(matrices_dir)
        #normalizeDropSeqCounts(matrices_dir,log=log)