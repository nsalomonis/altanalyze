###export
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

import os
import sys
import string
import unique
import UI

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def findParentDir(filename):
    ### :: reverses string
    x = string.find(filename[::-1],'\\')*-1 ### get just the parent directory
    if x == 1:
        x = string.find(filename[::-1],'//')*-1 ### get just the parent directory
    if x == 1:
        x = string.find(filename[::-1],'/')*-1 ### get just the parent directory
    return filename[:x]

def findFilename(filename):
    x = string.find(filename[::-1],'\\')*-1 ### get just the parent directory
    if x == 1:
        x = string.find(filename[::-1],'//')*-1 ### get just the parent directory
    if x == 1:
        x = string.find(filename[::-1],'/')*-1 ### get just the parent directory
    return filename[x:]

def ExportFile(filename):
    filename = string.replace(filename,'//','/')
    dir = findParentDir(filename)
    file_var = createExportFile(filename,dir)
    return file_var

def isFileOpen(new_file,dir):
    try: 
        file_open = 'yes'
        dir_list = read_directory(dir)
        if len(dir_list)>0:
            while file_open == 'yes':
                file_open = 'no'
                for file in dir_list:
                    if file in new_file: ###Thus the file is open
                        try:
                            fn=filepath(new_file); file_var = open(fn,'w')
                            file_open = 'no'
                        except IOError:
                            print_out = 'Results file: '+new_file+ '\nis open...can not re-write.\nPlease close file and select "OK".'
                            try: UI.WarningWindow(print_out,' OK ');
                            except Exception:
                                print print_out; print 'Please correct (hit return to continue)'
                                inp = sys.stdin.readline()
                            file_open = 'yes'
    except OSError: null = []
            
def createExportFile(new_file,dir):
    try:
        isFileOpen(new_file,dir)
        fn=filepath(new_file); file_var = open(fn,'w')
    except IOError:
        fn = filepath(dir)
        try:
            os.mkdir(fn) ###Re-Create directory if deleted
            print fn, 'written'
        except OSError:
            createExportDir(new_file,dir) ###Occurs if the parent directory is also missing
        fn=filepath(new_file); file_var = open(fn,'w')
    return file_var

def createExportDir(new_file,dir):
    dir = string.replace(dir,'//','/')
    dir = string.replace(dir,'\\','/')
    dir = string.replace(dir,'\\','/')
    dir_ls = string.split(dir,'/')
    if len(dir_ls) != 1:
        index = 1
        while index < (len(dir_ls)+1):
            parent_dir = string.join(dir_ls[:index],'/')
            index+=1
            print "Trying to create the directory:",parent_dir
            try:
                pfn = filepath(parent_dir)
                print pfn
                os.mkdir(pfn)
                print parent_dir, 'written' #, index-1, len(dir_ls)
            except OSError:
                print "Can not write this dir"
                #break
                continue
        createExportFile(new_file,dir)
    else: print "Parent directory not found locally for", dir_ls; sys.exit()

def createExportFolder(dir):
    dir = string.replace(dir,'//','/')
    dir = string.replace(dir,'\\','/')
    dir = string.replace(dir,'\\','/')
    dir_ls = string.split(dir,'/')
    if len(dir_ls) != 1:
        index = 1
        while index < (len(dir_ls)+1):
            parent_dir = string.join(dir_ls[:index],'/')
            index+=1
            print "Trying to create the directory:",parent_dir
            try:
                pfn = filepath(parent_dir)
                print pfn
                os.mkdir(pfn)
                print parent_dir, 'written' #, index-1, len(dir_ls)
            except OSError:
                print "Can not write this dir"
                #break
                continue
    else: print "Parent directory not found locally for", dir_ls; sys.exit()    

def deleteFolder(dir):
    try:
        dir = filepath(dir); dir_list = read_directory(dir) ### Get all files in directory
        for file in dir_list:
            fn = filepath(dir+'/'+file)
            os.remove(fn)
        os.removedirs(dir)
        return 'success'
    except OSError:
        return 'failed'

if __name__ == '__main__':
    fn = '/Users/nsalomonis/Desktop/GSE13297_RAW//AltExpression/ExonArray/Hs/Hs_Exon_CS_vs_hESC.p5_average.txt'
    createExportFile(fn,'/Users/nsalomonis/Desktop/GSE13297_RAW//AltExpression/ExonArray/Hs')
    ExportFile(fn)