import sys,string,os,math
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

def getFiles(sub_dir):
    dir_list = os.listdir(sub_dir); dir_list2 = []
    for entry in dir_list:
        if '.csv' in entry: dir_list2.append(entry)
    return dir_list2

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def combineCSVFiles(root_dir):
    files = getFiles(root_dir)
    import export
    folder = export.getParentDir(root_dir+'/blah.txt')
    first_file = True
    genes=[]
    cells=[]
    data_matrices=[]
    for file in files:
        cells.append(file[:-4])
        matrix=[]
        for line in open(root_dir+'/'+file,'rU').xreadlines():
            data = cleanUpLine(line)
            gene,count = string.split(data,',')
            if first_file:
                genes.append(gene)
            matrix.append(float(count))
        data_matrices.append(matrix)
        first_file=False
        
    data_matrices = zip(*data_matrices)
    export_object = open(root_dir+'/'+folder+'-counts.txt','w')
    headers = string.join(['UID']+cells,'\t')+'\n'
    export_object.write(headers)
    index=0
    for geneID in genes:
        values = string.join([geneID]+map(str,list(data_matrices[index])),'\t')+'\n'
        export_object.write(values)
        index+=1
    export_object.close()

    data_matrices = zip(*data_matrices)
    
    index=0
    for cell in cells:
        matrix = data_matrices[index]
        barcode_sum = sum(matrix)
        data_matrices[index] = map(lambda val: math.log((10000.00*val/barcode_sum)+1.0,2), matrix)
        index+=1
        
    data_matrices = zip(*data_matrices)

    export_object = open(root_dir+'/'+folder+'-CPTT.txt','w')
    headers = string.join(['UID']+cells,'\t')+'\n'
    export_object.write(headers)
    index=0
    for geneID in genes:
        values = string.join([geneID]+map(str,list(data_matrices[index])),'\t')+'\n'
        export_object.write(values)
        index+=1
    export_object.close()
    
if __name__ == '__main__':
    ################  Comand-line arguments ################
    import getopt
    options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','GeneType=', 'IDType='])
    for opt, arg in options:
        if opt == '--i': input_dir=arg
            
    combineCSVFiles(input_dir)

