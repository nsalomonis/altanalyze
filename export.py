import os
import sys
import string
import unique
dirfile = unique
py2app_adj = '/AltAnalyze.app/Contents/Resources/Python/site-packages.zip'

def filepath(filename):
    dir=os.path.dirname(dirfile.__file__)       #directory file is input as a variable under the main            
    fn=os.path.join(dir,filename)
    fn = string.replace(fn,py2app_adj,'')
    fn = string.replace(fn,'\\library.zip','') ###py2exe on some systems, searches for all files in the library file, eroneously
    return fn

def createExportFile(new_file,dir):
    try:
        fn=filepath(new_file); file_var = open(fn,'w')
    except IOError:
        print "IOError", fn
        fn = filepath(dir)
        try:
            os.mkdir(fn) ###Re-Create directory if deleted
            #print fn, 'written'
        except OSError: createExportDir(new_file,dir) ###Occurs if the parent directory is also missing
        fn=filepath(new_file); file_var = open(fn,'w')
    return file_var

"""def createExportDir(fn,dir):
    dir_ls = string.split(dir,'//')
    if len(dir_ls) == 1: ### Thus, '//' not found
        dir_ls = string.split(dir,'/')
    if len(dir_ls) != 1: parent_dir = string.join(dir_ls[:-1],'/')
    else: #print "Parent directory not found locally for", dir_ls; sys.exit()
    pfn = filepath(parent_dir)
    print "Trying to create the directory:",parent_dir
    try:
        os.mkdir(pfn) ###Re-Create directory if deleted
        try: os.mkdir(fn) ###Re-Create directory if deleted
        except OSError: createExportFile(fn,dir) ###Occurs if the parent directory is also missing        
    except OSError: createExportDir(fn,pfn) ###Occurs if the parent directory is also missing"""

def createExportDir(new_file,dir):
    dir_ls = string.split(dir,'//')
    if len(dir_ls) == 1: ### Thus, '//' not found
        dir_ls = string.split(dir,'/')
    if len(dir_ls) != 1:
        index = 1
        while index < (len(dir_ls)+1):
            parent_dir = string.join(dir_ls[:index],'/')
            index+=1
            print "Trying to create the directory:",parent_dir
            try:
                pfn = filepath(parent_dir)
                os.mkdir(pfn)
                print parent_dir, 'written' #, index-1, len(dir_ls)
            except OSError:
                #print "Can not write this dir"
                continue
        createExportFile(new_file,dir)
    else: print "Parent directory not found locally for", dir_ls; sys.exit()