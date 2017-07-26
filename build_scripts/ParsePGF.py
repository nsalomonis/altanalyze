###AltAnalyze
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
import os.path
import platform

def filepath(filename):
    #fn = unique.filepath(filename)
    return filename

def read_directory(sub_dir):
    dir_list = os.listdir(sub_dir)
    return dir_list

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        ls = list(set(database[key]))
        ls.sort()
        db1[key] = ls
    return db1

class GrabFiles:
    def setdirectory(self,value):
        self.data = value
    def display(self):
        print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        file_dir,file = getDirectoryFiles(self.data,str(search_term))
        if len(file)<1: print search_term,'not found'
        return file_dir,file
    
def getDirectoryFiles(import_dir, search_term):
    exact_file = ''; exact_file_dir=''
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    import_dir = filepath(import_dir)
    for data in dir_list:    #loop through each file in the directory to output results
        if (':' in import_dir) or ('/Users/' == import_dir[:7]) or ('Linux' in platform.system()): affy_data_dir = import_dir+'/'+data
        else: affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir: exact_file_dir = affy_data_dir; exact_file = data
    return exact_file_dir,exact_file

def makeUnique(item):
    db1={}; list1=[]; k=0
    for i in item:
        try: db1[i]=[]
        except TypeError: db1[tuple(i)]=[]; k=1
    for i in db1:
        if k==0: list1.append(i)
        else: list1.append(list(i))
    list1.sort()
    return list1

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def importPGF(dir,species,filename):
    fn=filepath(filename); probe_db = {}; x=0
    psr_file = dir+'/'+species+'/'+array_type+'/'+species+'_probeset-psr.txt'
    psr_file = string.replace(psr_file,'affymetrix/LibraryFiles/','')
    try: eo = open(filepath(psr_file),'w')
    except Exception: eo = open(filepath(psr_file[1:]),'w')
    for line in open(fn,'rU').xreadlines():
        if line[0] != '#':
            data = cleanUpLine(line); x+=1
            t = string.split(data,'\t')
            if len(t)==2 or len(t)==3:
                if len(t[0])>0:
                    probeset = t[0]; type = t[1]
                    eo.write(probeset+'\t'+t[-1]+'\n') ### Used for HTA array where we need to have PSR to probeset IDs
            else:
                try:
                    probe = t[2]
                    #if probeset == '10701621': print probe
                    try: probe_db[probeset].append(probe)
                    except KeyError: probe_db[probeset] = [probe]
                except Exception: null=[]
    eo.close() 
    new_file = dir+'/'+species+'/'+array_type+'/'+species+'_probeset-probes.txt'
    new_file = string.replace(new_file,'affymetrix/LibraryFiles/','')
    headers = 'probeset\t' + 'probe\n'; n=0
    try: data = open(filepath(new_file),'w')
    except Exception: data = open(filepath(new_file[1:]),'w')
    data.write(headers)
    for probeset in probe_db:
        for probe in probe_db[probeset]:
            data.write(probeset+'\t'+probe+'\n'); n+=1
    data.close()
    print n, 'Entries exported for', new_file
    
if __name__ == '__main__':
    skip_intro = 'yes'
    array_type = 'gene'
    #array_type = 'exon'
    #array_type = 'junction'
    array_type = 'gene'
    species = 'Mm'
    parent_dir = 'AltDatabase/'+species+'/'+array_type+'/library'
    parent_dir = '/AltDatabase/affymetrix/LibraryFiles'
    e = GrabFiles(); e.setdirectory(parent_dir)
    pgf_dir,pgf_file = e.searchdirectory('MoGene-2_0-st.pgf')
    importPGF(parent_dir,species,pgf_dir)
