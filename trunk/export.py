###export
#Copyright 2005-2008 J. Davide Gladstone Institutes, San Francisco California
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

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def createExportFile(new_file,dir):
    try:
        fn=filepath(new_file); file_var = open(fn,'w')
    except IOError:
        #print "IOError", fn
        fn = filepath(dir)
        try:
            os.mkdir(fn) ###Re-Create directory if deleted
            #print fn, 'written'
        except OSError: createExportDir(new_file,dir) ###Occurs if the parent directory is also missing
        fn=filepath(new_file); file_var = open(fn,'w')
    return file_var

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