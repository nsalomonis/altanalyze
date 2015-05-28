import sys,string
import os
import unique

def importFPKMFile(input_file,folder):
    firstLine = True
    print folder
    for line in open(input_file,'rU').xreadlines():
        data = line.rstrip('\n')
        t = string.split(data,'\t')
        if firstLine:
            firstLine = False
        else:
            try: geneID= t[3]; fpkm = float(t[9])
            except Exception: print t;sys.exit()
            try: fpkm_db[geneID].append(fpkm)
            except Exception: fpkm_db[geneID] = [fpkm]

def getFiles(sub_dir,directories=True):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    for entry in dir_list:
        if directories:
            if '.' not in entry: dir_list2.append(entry)
        else:
            if '.' in entry: dir_list2.append(entry)
    return dir_list2

def combineCufflinks(root_dir):
    export_object = open(root_dir+'/Cufflinks.txt','w')
    global fpkm_db
    global headers
    fpkm_db={}; headers=['UID']
    folders = getFiles(root_dir,True)
    for folder in folders:
        l1 = root_dir+'/'+folder
        files = getFiles(l1,False)
        for file in files:
            filename = l1+'/'+file
            if 'genes.fpkm_tracking' in file and '.gz' not in file:
                importFPKMFile(filename,folder)
                headers.append(folder)
                
    headers = string.join(headers,'\t')+'\n'
    export_object.write(headers)
    for gene in fpkm_db:
        values = map(str,fpkm_db[gene])
        values = string.join([gene]+values,'\t')+'\n'
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
    combineCufflinks('/Users/saljh8/Downloads/6b_CUFFLINKS_output/');sys.exit()
    
    ################  Comand-line arguments ################
    import getopt
    options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','f='])
    #print sys.argv[1:]
    for opt, arg in options:
        if opt == '--i': input_file=arg
            
    if '/' in input_file: delim = '/'
    else: delim = '\\'
    root_dir = string.join(string.split(input_file,delim)[:-1], delim)
    processFile(input_file,root_dir)

