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
            try: geneID= t[0]; symbol=t[4]; position=t[6]; fpkm = t[9]
            except Exception: geneID= t[0]; symbol=t[1]; position=t[2]; fpkm = t[4]
            if 'chr' not in position: position = 'chr'+position
            try:
                null=position_symbol_db[position]
                if len(symbol)>1: position_symbol_db[position].append(symbol)
            except Exception: position_symbol_db[position] = [symbol]
            try:
                null=position_gene_db[position]
                if '.' not in geneID: position_gene_db[position].append(geneID)
            except Exception: position_gene_db[position] = [geneID]

            if IDType == 'symbol':
                if symbol!= '-': geneID = symbol
            if IDType == 'position':
                geneID = position
            if getData:
                try: fpkm_db[geneID].append(fpkm)
                except Exception: fpkm_db[geneID] = [fpkm]
                added_key[geneID]=[]
            else:
                fpkm_db[geneID]=[]
                added_key[geneID]=[]
    
    for i in fpkm_db:
        if i not in added_key:
            fpkm_db[i].append('0.00')
    
def getFiles(sub_dir,directories=True):
    dir_list = os.listdir(sub_dir); dir_list2 = []
    for entry in dir_list:
        if directories:
            if '.' not in entry: dir_list2.append(entry)
        else:
            if '.' in entry: dir_list2.append(entry)
    return dir_list2

def combineCufflinks(root_dir,gene_type,id_type):
    global fpkm_db
    global headers
    global IDType; IDType = id_type
    global position_symbol_db; position_symbol_db={}
    global position_gene_db; position_gene_db={}
    global getData; getData=False
    fpkm_db={}; headers=['GeneID']
    folders = getFiles(root_dir,True)
    for folder in folders:
        l1 = root_dir+'/'+folder
        """  
        files = getFiles(l1,False)
        for file in files:
            
            filename = l1+'/'+file
            if gene_type in file and '.gz' not in file:
                print filename
                importFPKMFile(filename)
                headers.append(folder)
                break
        break
        """
        folders2 = getFiles(l1,True)
        for folder in folders2:
            l2 = l1+'/'+folder
            files = getFiles(l2,False)
            for file in files:
                filename = l2+'/'+file
                if gene_type in file and '.gz' not in file:
                    print filename
                    importFPKMFile(filename)
                    headers.append(folder)    
    for i in fpkm_db:
        fpkm_db[i]=[]

    getData=True
    headers=['UID']
    if IDType == 'position':
        headers=['Position','GeneID','Symbol']
    for folder in folders:
        l1 = root_dir+'/'+folder
        """
        files = getFiles(l1,False)
        for file in files:
            
            filename = l1+'/'+file
            if gene_type in file and '.gz' not in file:
                print filename
                importFPKMFile(filename)
                headers.append(folder)
                break
        break
        """
        folders2 = getFiles(l1,True)
        for folder in folders2:
            l2 = l1+'/'+folder
            files = getFiles(l2,False)
            for file in files:
                filename = l2+'/'+file
                if gene_type in file and '.gz' not in file:
                    print filename
                    importFPKMFile(filename)
                    headers.append(folder)
                    
    if IDType == 'position':
        for position in position_gene_db:
            gene = '-'
            #print position,position_gene_db[position]
            for i in position_gene_db[position]:
                if '.' not in i: gene = i
            position_gene_db[position] = gene
            #print position,position_gene_db[position],'\n'
            #print position,position_symbol_db[position]
            symbol = '-'
            for i in position_symbol_db[position]:
                if i != '-': symbol = i
            position_symbol_db[position] = symbol
            #print position,position_symbol_db[position]
    #print position_symbol_db[position]

    export_object = open(root_dir+'/'+gene_type+'-Cufflinks.txt','w')
    headers = string.join(headers,'\t')+'\n'
    export_object.write(headers)
    for geneID in fpkm_db:
        values = map(str,fpkm_db[geneID])
        if IDType != 'position':
            values = string.join([geneID]+values,'\t')+'\n'
        else:
            values = string.join([geneID,position_gene_db[geneID],position_symbol_db[geneID]]+values,'\t')+'\n'
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
    #combineCufflinks('/Users/saljh8/Downloads/6b_CUFFLINKS_output/');sys.exit()
    type = 'genes'
    id_type = 'position'
    ################  Comand-line arguments ################
    import getopt
    options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','GeneType=', 'IDType='])
    for opt, arg in options:
        if opt == '--i': input_dir=arg
        if opt == '--GeneType': type=arg
        if opt == '--IDType': id_type=arg
            
    combineCufflinks(input_dir,type,id_type)

