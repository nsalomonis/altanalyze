import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

def importFPKMFile(input_file):
    added_key={}
    firstLine = True
    for line in open(input_file,'rU').xreadlines():
        data = line.rstrip('\n')
        t = string.split(data,'\t')
        if firstLine:
            firstLine = False
        else:
            try: geneID= t[0]; symbol=t[1]; position=t[2]; fpkm = t[5]
            except Exception: print t;sys.exit()
            try: fpkm_db[geneID].append(fpkm)
            except Exception: fpkm_db[geneID] = [fpkm]
            added_key[geneID]=[]
    
    for i in fpkm_db:
        if i not in added_key:
            fpkm_db[i].append('0.00')
    
def importCufflinksFPKMFileEXCEL(filename):
    from xlrd import open_workbook
    print filename
    wb = open_workbook(filename)
    for w in wb.sheets():
        print w.name;sys.exit()
    sys.exit()
    rows=[]
    for row in range(worksheet.nrows):
        print row
        values = []
        for col in range(worksheet.ncols):
            try: values.append(str(s.cell(row,col).value))
            except Exception: pass
        rows.append(values)

def getFiles(sub_dir,directories=True):
    dir_list = os.listdir(sub_dir); dir_list2 = []
    for entry in dir_list:
        if directories:
            if '.' not in entry: dir_list2.append(entry)
        else:
            if '.' in entry: dir_list2.append(entry)
    return dir_list2

def combineCufflinks(root_dir):
    export_object = open(root_dir+'/RSEM.txt','w')
    global fpkm_db
    global headers
    fpkm_db={}; headers=['GeneID']
    files = getFiles(root_dir,False)
    for file in files:
        filename = root_dir+'/'+file
        if ('.genes.results' in file and '.gz' not in file) or ('.txt' in file and '.gz' not in file and 'RSEM.txt' not in file):
            importFPKMFile(filename)
            headers.append(file)
        if '.xls' in file:
            importCufflinksFPKMFileEXCEL(filename)
            headers.append(file)

                
                
    for i in fpkm_db:
        fpkm_db[i]=[]

    for file in files:
        filename = root_dir+'/'+file
        if ('.genes.results' in file and '.gz' not in file) or ('.txt' in file and '.gz' not in file and 'RSEM.txt' not in file):
            importFPKMFile(filename)
    
    headers = string.join(headers,'\t')+'\n'
    export_object.write(headers)
    for geneID in fpkm_db:
        values = map(str,fpkm_db[geneID])
        values = string.join([geneID]+values,'\t')+'\n'
        export_object.write(values)
    export_object.close()
    
                
def gunzipfiles(root_dir):
    import gzip
    import shutil;
    folders = getFiles(root_dir,True)
    for folder in folders:
        l1 = root_dir+'/'+folder
        files = getFiles(l1,False)
        for file in files:
            filename = l1+'/'+file
            if 'genes.fpkm_tracking' in filename:
                content = gzip.GzipFile(filename, 'rb')
                decompressed_filepath = string.replace(filename,'.gz','')
                data = open(decompressed_filepath,'wb')
                shutil.copyfileobj(content,data) 

if __name__ == '__main__':
    #gunzipfiles('/Users/saljh8/Downloads/6b_CUFFLINKS_output/');sys.exit()
    
    import getopt
    options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','f='])
    #print sys.argv[1:]
    for opt, arg in options:
        if opt == '--i': input_file=arg
            
    if '/' in input_file: delim = '/'
    else: delim = '\\'
    combineCufflinks(input_file);sys.exit()

