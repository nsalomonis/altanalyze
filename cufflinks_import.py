#!/usr/local/bin/python2.6

###update
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

"""This module contains generic methods for downloading and decompressing files designated
online files and coordinating specific database build operations for all AltAnalyze supported
gene ID systems with other program modules. """

import os
import sys
import unique
import string
import export
import gzip

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
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
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Only get folder names
    for entry in dir_list:
        dir_list2.append(entry)
    return dir_list2

def zipDirectory(dir):
    #http://www.testingreflections.com/node/view/8173
    import zipfile
    dir = filepath(dir); zip_file = dir+'.zip'
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
    output_filepath = filepath(dir+filename)
    try:
        zfile = zipfile.ZipFile(output_filepath)
        for name in zfile.namelist():
            if name.endswith('/'):null=[] ### Don't need to export
            else:
                try: outfile = export.ExportFile(filepath(dir+name))
                except Exception: outfile = export.ExportFile(filepath(dir+name[1:]))
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
    sample_FPKM_db={}
    for sample in root_dir_files: ### e.g.,
        try:
            files = getFiles(directory+'/'+sample)
            for file in files:
                if file == 'isoforms.fpkm_tracking.gz':
                    x_db,transcript_db = readFPKMs(directory+'/'+sample+'/'+file)
                    sample_FPKM_db.uppdate(x_db)
        except Exception: pass
    
    ### Get the samples
    samples = sample_FPKM_db.keys()
    samples.sort()
    ### Get the transcripts
    fpkm_db = sample_FPKM_db[samples[0]]
    transcripts = fpkm_db.keys()                
    for transcript in transcripts:
        for sample in samples:
            fpkm = sample_FPKM_db[sample][transcript]
            fpkm_list.append(fpkm)

def readFPKMs(path):
    f=gzip.open(path,'rb')
    file_content=f.read()
    fpkm_data = string.split(file_content,'\n')
    sample = export.findFilename(path)
    fpkm_db={}
    transcript_db={}
    firstLine=True
    for line in fpkm_data:
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            track_i = t.index('tracking_id')
            gene_i = t.index('gene_id')
            fpkm_i = t.index('FPKM')
            firstLine = False
        else:
            geneID = t[gene_i]
            transcriptID = t[gene_i]
            fpkm = t[fpkm_i]
            fpkm_db[transcriptID] = float(fpkm)
            transcript_db[transcriptID] = geneID
    sample_FPKM_db[sample] = fpkm_db
    return sample_FPKM_db,transcript_db
    
if __name__ == '__main__':
    dir = '/Volumes/SEQ-DATA/Cufflinks'
    importCufflinksDir(dir)

