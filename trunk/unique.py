###unique
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

"""This module contains instructions for recalling operating system file and directory paths,
eliminating redundant list entries, removing unecessary file paths from py2app or py2exe
and reading the propper Ensembl database version to allow for version specific access."""

import sys, string
import os.path, platform
import unique ### Import itself as a reference to it's location
dirfile = unique

py2app_adj = '/GO_Elite.app/Contents/Resources/Python/site-packages.zip'
py2app_adj1 = '/GO_Elite.app/Contents/Resources/lib/python2.4/site-packages.zip'
py2app_adj2 = '/GO_Elite.app/Contents/Resources/lib/python2.5/site-packages.zip'
py2app_adj3 = '/GO_Elite.app/Contents/Resources/lib/python2.6/site-packages.zip'
py2app_adj4 = '/GO_Elite.app/Contents/Resources/lib/python2.7/site-packages.zip'
py2exe_adj = '\\library.zip' ###py2exe
py2app_ge_dirs = [py2app_adj,py2exe_adj,py2app_adj1,py2app_adj2,py2app_adj3,py2app_adj4]

py2app_adj = '/AltAnalyze.app/Contents/Resources/Python/site-packages.zip'
py2app_adj1 = '/AltAnalyze.app/Contents/Resources/lib/python2.4/site-packages.zip'
py2app_adj2 = '/AltAnalyze.app/Contents/Resources/lib/python2.5/site-packages.zip'
py2app_adj3 = '/AltAnalyze.app/Contents/Resources/lib/python2.6/site-packages.zip'
py2app_adj4 = '/AltAnalyze.app/Contents/Resources/lib/python2.7/site-packages.zip'
py2exe_adj = '\\library.zip' ###py2exe
py2app_aa_dirs = [py2app_adj,py2app_adj1,py2exe_adj,py2app_adj2,py2app_adj3,py2app_adj4]
py2app_dirs = py2app_ge_dirs + py2app_aa_dirs

def filepath(filename):
    dir=os.path.dirname(dirfile.__file__)       #directory file is input as a variable under the main
    filename = string.replace(filename,'//','/'); filename = string.replace(filename,'\\','/')
    folders = string.split(filename,'/'); top_dir = folders[0]
    try: top_dirs = returnDirectories(dir) ### useful for Linux
    except Exception: top_dirs = [] ### Error arrises with PC compiled version
    if len(filename)==0: fn = dir+'/'
    elif top_dir in top_dirs: fn=os.path.join(dir,filename)
    elif (':' in filename) or ('/Users/' == filename[:7]) or ('/Volumes/' in filename) or ('Linux' in platform.system()): fn = filename #':' is for Windows dirs and '/Users/' for Mac
    else: fn=os.path.join(dir,filename)
    if '/Volumes/' in filename: filenames = string.split(filename,'/Volumes/'); fn = '/Volumes/'+filenames[-1]
    for py2app_dir in py2app_dirs: fn = string.replace(fn,py2app_dir,'')
    if 'Databases' in fn or 'AltDatabase' in fn:
        getCurrentGeneDatabaseVersion()
        fn = correctGeneDatabaseDir(fn)
    return fn

def read_directory(sub_dir):
    dir=os.path.dirname(dirfile.__file__)
    for py2app_dir in py2app_dirs: dir = string.replace(dir,py2app_dir,'')
    if 'Databases' in sub_dir or 'AltDatabase' in sub_dir:
        getCurrentGeneDatabaseVersion()
        sub_dir = correctGeneDatabaseDir(sub_dir)
    try:
        if (':' in sub_dir) or ('/Users/' == sub_dir[:7]) or ('/Volumes/' in sub_dir) or ('Linux' in platform.system()): dir_list = os.listdir(dir+sub_dir) ### Thus, the whole path is provided already
        else:
            try: dir_list = os.listdir(dir + sub_dir)
            except Exception: dir_list = os.listdir(sub_dir) ### For linux
    except Exception:
        if (':' in sub_dir) or ('/Users/' == sub_dir[:7]) or ('/Volumes/' in sub_dir) or ('Linux' in platform.system()): dir_list = os.listdir(sub_dir) ### Thus, the whole path is provided already
        else:
            try: dir_list = os.listdir(dir + sub_dir)
            except Exception: dir_list = os.listdir(sub_dir) ### For linux
    return dir_list
    
def returnDirectories(sub_dir):
    dir=os.path.dirname(dirfile.__file__)
    if 'Databases' in sub_dir or 'AltDatabase' in sub_dir:
        getCurrentGeneDatabaseVersion()
        sub_dir = correctGeneDatabaseDir(sub_dir)
    for py2app_dir in py2app_dirs:
        dir = string.replace(dir,py2app_dir,'')
        try: dir_list = os.listdir(dir + sub_dir)
        except Exception: dir_list = os.listdir(sub_dir) ### For linux
    return dir_list

def returnDirectoriesNoReplace(sub_dir):
    dir=os.path.dirname(dirfile.__file__)
    for py2app_dir in py2app_dirs:
        dir = string.replace(dir,py2app_dir,'')
    try: dir_list = os.listdir(dir + sub_dir)
    except Exception: dir_list = os.listdir(sub_dir) ### For linux
    return dir_list

def refDir():
    reference_dir=os.path.dirname(dirfile.__file__)       #directory file is input as a variable under the main            
    for py2app_dir in py2app_dirs: 
        reference_dir = string.replace(reference_dir,py2app_adj,'')
    return reference_dir

def whatProgramIsThis():
    reference_dir = refDir()
    if 'AltAnalyze' in reference_dir: type = 'AltAnalyze'; database_dir = 'AltDatabase/goelite'
    elif 'GO-Elite' in reference_dir: type = 'GO-Elite'; database_dir = 'Databases'
    if 'Linux' in platform.system(): database_dir += '/'
    return type,database_dir

def correctGeneDatabaseDir(fn):
    try:
        if gene_database_dir not in fn and 'EnsMart' not in fn:
            fn = string.replace(fn,'Databases','Databases/'+gene_database_dir)
            if 'AltDatabase/affymetrix' not in fn and 'NoVersion' not in fn:
                fn = string.replace(fn,'AltDatabase','AltDatabase/'+gene_database_dir)
        fn = string.replace(fn,'NoVersion','') ### When the text 'NoVersion' is in a filepath, is tells the program to ignore it for adding the database version
    except Exception: null = ''
    return fn

def getCurrentGeneDatabaseVersion():
    global gene_database_dir
    try:
        filename = 'Config/version.txt'; fn=filepath(filename)
        for line in open(fn,'r').readlines():
            gene_database_dir, previous_date = string.split(line,'\t')
    except Exception: gene_database_dir=''
    return gene_database_dir
    
def unique(s):
    #we need to remove duplicates from a list, unsuccessfully tried many different methods
    #so I found the below function at: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52560
        n = len(s)
        if n == 0: return []    
        u = {}
        try:
            for x in s: u[x] = 1
        except TypeError: del u  # move on to the next method
        else: return u.keys()
        try: t = list(s); t.sort()
        except TypeError: del t  # move on to the next method
        else:
            assert n > 0
            last = t[0]; lasti = i = 1
            while i < n:
                if t[i] != last: t[lasti] = last = t[i]; lasti += 1
                i += 1
            return t[:lasti]
        u = []
        for x in s:
            if x not in u: u.append(x)
        return u
    
def dictionary(s):
    d={}
    for i in s:
        try: d[i]=[]
        except TypeError: d[tuple(i)]=[]
    return d

def unique_db(s):
    d={}; t=[]
    for i in s:
        try: d[i]=[]
        except TypeError: d[tuple(i)]=[]
    for i in d: t.append(i)
    return t

def list(d):
    t=[]
    for i in d: t.append(i)
    return t

if __name__ == '__main__':
    fn = filepath('Databases/Mm/uid-gene/null.txt')
    print fn;kill
    unique_db([1,2,3,4,4,4,5])
