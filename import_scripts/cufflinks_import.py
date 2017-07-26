#!/usr/local/bin/python2.6

###cufflinks_import
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import gzip
import getopt

def read_directory(sub_dir):
    dir_list = os.listdir(sub_dir)
    dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def getFiles(sub_dir):
    dir_list = os.listdir(sub_dir)
    dir_list2 = []
    ###Only get folder names
    for entry in dir_list:
        dir_list2.append(entry)
    return dir_list2

def zipDirectory(dir):
    #http://www.testingreflections.com/node/view/8173
    import zipfile
    zip_file = dir+'.zip'
    p = string.split(dir,'/'); top=p[-1]
    zip = zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for f in files:
            fullpath = os.path.join(root, f)
            archive_name = os.path.join(top+archive_root, f)
            zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()
    return zip_file

def unzipFiles(filename,dir):
    import zipfile
    output_filepath = dir+filename
    try:
        zfile = zipfile.ZipFile(output_filepath)
        for name in zfile.namelist():
            if name.endswith('/'):null=[] ### Don't need to export
            else:
                try: outfile = open(dir+name,'w')
                except Exception:
                    outfile = open(dir+name[1:],'w')
                outfile.write(zfile.read(name)); outfile.close()
        #print 'Zip extracted to:',output_filepath
        status = 'completed'
    except Exception, e:
        try:
            ### Use the operating system's unzip if all else fails
            extracted_path = string.replace(output_filepath,'.zip','')
            try: os.remove(extracted_path) ### This is necessary, otherwise the empty file created above will require user authorization to delete
            except Exception: null=[]
            subprocessUnzip(dir,output_filepath)
            status = 'completed'
        except IOError:
            print e
            print 'WARNING!!!! The zip file',output_filepath,'does not appear to be a valid zip archive file or is currupt.'
            status = 'failed'
    return status

def gunzip():
    import gzip; content = gzip.GzipFile(gz_filepath, 'rb')
    data = open(decompressed_filepath,'wb')
    import shutil; shutil.copyfileobj(content,data)

def importCufflinksDir(directory):
    root_dir_files = getFiles(directory)
    global sample_FPKM_db
    sample_FPKM_db={}
    for sample in root_dir_files: ###
        if 'fpkm_tracking' in sample:
            x_db,transcript_db = readFPKMs(directory+'/'+sample)
            sample_FPKM_db.update(x_db)
        ### Below occurs if it is a directory of FPKM results with the folder name as the sample
        try:
            files = getFiles(directory+'/'+sample)
            for file in files:
                if 'fpkm_tracking' in file:
                    x_db,transcript_db = readFPKMs(directory+'/'+sample+'/'+file)
                    sample_FPKM_db.update(x_db)
        except Exception: pass
    
    ### Get the samples
    samples = sample_FPKM_db.keys()
    samples.sort()
    ### Get the transcripts
    gene_fpkm_db={}
    for sample in samples:
        fpkm_db = sample_FPKM_db[sample]
        for gene in fpkm_db:
            fpkm = fpkm_db[gene]
            try: gene_fpkm_db[gene].append(fpkm)
            except Exception: gene_fpkm_db[gene] = [fpkm]
        
    export_object = open(directory+'/Cufflinks.txt','w')
    headers = string.join(['UID']+samples,'\t')+'\n'
    export_object.write(headers)
    for geneID in gene_fpkm_db:
        values = map(str,gene_fpkm_db[geneID])
        values = string.join([geneID]+values,'\t')+'\n'
        export_object.write(values)
    export_object.close()
    
def findFilename(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'//','/') ### If /// present
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1 ### get just the parent directory
    return filename[x:]

def readFPKMs(path):
    if '.gz' in path:
        f=gzip.open(path,'rb')
    else:
        f=open(path,"rU")
    file_content=f.read()
    fpkm_data = string.split(file_content,'\n')
    sample = findFilename(path)
    if 'fpkm_tracking' in sample:
        sample = string.split(sample,'.fpkm_tracking')[0]
        sample = string.replace(sample,'.sorted.genes','')
    fpkm_db={}
    transcript_db={}
    firstLine=True
    row_count=0
    for line in fpkm_data:
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            try:
                track_i = t.index('tracking_id')
                gene_i = t.index('gene_id')
                fpkm_i = t.index('FPKM')
            except Exception:
                fpkm_i = 9
                gene_i = 3
                row_count = 1
            firstLine = False
        if firstLine == False and row_count>0:
            if len(t)>1:
                geneID = t[gene_i]
                transcriptID = t[gene_i]
                fpkm = t[fpkm_i]
                fpkm_db[transcriptID] = float(fpkm)
                transcript_db[transcriptID] = geneID
        row_count+=1
    sample_FPKM_db[sample] = fpkm_db
    return sample_FPKM_db,transcript_db

if __name__ == "__main__":
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a directory of fpkm_tracking file as input in the command-line"
        print "Example: python cufflinks_import.py --i /Users/cufflinks/"
        sys.exit()
    else:
        Species = None
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','species='])
        for opt, arg in options:
            if opt == '--i': dir=arg ### full path of a BAM file
            elif opt == '--species': Species=arg ### full path of a BAM file
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
            
    try: importCufflinksDir(dir)
    except ZeroDivisionError:
        print [sys.argv[1:]],'error'; error

