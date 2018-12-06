import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os.path
import unique
import itertools
dirfile = unique

############ File Import Functions #############

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def returnDirectories(sub_dir):
    dir_list = unique.returnDirectories(sub_dir)
    return dir_list

class GrabFiles:
    def setdirectory(self,value): self.data = value
    def display(self): print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        files = getDirectoryFiles(self.data,search_term)
        if len(files)<1: print 'files not found'
        return files
    def returndirectory(self):
        dir_list = getAllDirectoryFiles(self.data)
        return dir_list

def getAllDirectoryFiles(import_dir):
    all_files = []
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        data_dir = import_dir[1:]+'/'+data
        if '.' in data:
            all_files.append(data_dir)
    return all_files

def getDirectoryFiles(import_dir,search_term):
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    matches=[]
    for data in dir_list:    #loop through each file in the directory to output results
        data_dir = import_dir[1:]+'/'+data
        if search_term not in data_dir and '.' in data: matches.append(data_dir)
    return matches

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data
    
def combineAllLists(files_to_merge,original_filename,includeColumns=False):
    headers =[]; all_keys={}; dataset_data={}; files=[]; unique_filenames=[]
    count=0
    for filename in files_to_merge:
        duplicates=[]
        count+=1
        fn=filepath(filename); x=0; combined_data ={}
        if '/' in filename:
            file = string.split(filename,'/')[-1][:-4]
        else:
            file = string.split(filename,'\\')[-1][:-4]
        ### If two files with the same name being merged
        if file in unique_filenames:
            file += str(count)
        unique_filenames.append(file)
        print file
        files.append(filename)
        for line in open(fn,'rU').xreadlines():         
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0:
                if data[0]!='#':
                    x=1
                    try: t = t[1:]
                    except Exception: t = ['null']
                    if includeColumns==False:
                        for i in t:    
                            headers.append(i+'.'+file)
                            #headers.append(i)
                    else:
                        headers.append(t[includeColumns]+'.'+file)
            else: #elif 'FOXP1' in data or 'SLK' in data or 'MBD2' in data:
                key = t[0]
                if includeColumns==False:
                    try: values = t[1:]
                    except Exception: values = ['null']
                    try:
                        if original_filename in filename and len(original_filename)>0: key = t[1]; values = t[2:]
                    except IndexError: print original_filename,filename,t;kill
                else:
                    values = [t[includeColumns]]
                #key = string.replace(key,' ','')
                if len(key)>0 and key != ' ' and key not in combined_data: ### When the same key is present in the same dataset more than once
                    try: all_keys[key] += 1
                    except KeyError: all_keys[key] = 1
                if permform_all_pairwise == 'yes':
                    try: combined_data[key].append(values); duplicates.append(key)
                    except Exception: combined_data[key] = [values]
                else:
                    combined_data[key] = values
        #print duplicates
        dataset_data[filename] = combined_data
    for i in dataset_data:
        print len(dataset_data[i]), i
    ###Add null values for key's in all_keys not in the list for each individual dataset
    combined_file_data = {}
    for filename in files:
        combined_data = dataset_data[filename]
        ###Determine the number of unique values for each key for each dataset
        null_values = []; i=0
        for key in combined_data: number_of_values = len(combined_data[key][0]); break
        while i<number_of_values: null_values.append('0'); i+=1
        for key in all_keys:
            include = 'yes'
            if combine_type == 'intersection':
                if all_keys[key]>(len(files_to_merge)-1): include = 'yes'
                else: include = 'no'
            if include == 'yes':
                try: values = combined_data[key]
                except KeyError:
                    values = null_values
                    if permform_all_pairwise == 'yes':
                        values = [null_values]
                if permform_all_pairwise == 'yes':
                    try:
                        val_list = combined_file_data[key]
                        val_list.append(values)
                        combined_file_data[key] = val_list
                    except KeyError: combined_file_data[key] = [values]
                else:
                    try:
                        previous_val = combined_file_data[key]
                        new_val = previous_val + values
                        combined_file_data[key] = new_val
                    except KeyError: combined_file_data[key] = values

    original_filename = string.replace(original_filename,'1.',  '1.AS-')
    export_file = output_dir+'/MergedFiles.txt'
    fn=filepath(export_file);data = open(fn,'w')
    title = string.join(['uid']+headers,'\t')+'\n'; data.write(title)
    for key in combined_file_data:
        #if key == 'ENSG00000121067': print key,combined_file_data[key];kill
        new_key_data = string.split(key,'-'); new_key = new_key_data[0]
        if permform_all_pairwise == 'yes':
            results = getAllListCombinations(combined_file_data[key])
            for result in results:
                merged=[]
                for i in result: merged+=i
                values = string.join([key]+merged,'\t')+'\n'; data.write(values) ###removed [new_key]+ from string.join
        else:
            try:
                values = string.join([key]+combined_file_data[key],'\t')+'\n'; data.write(values) ###removed [new_key]+ from string.join
            except Exception: print combined_file_data[key];sys.exit()
    data.close()
    print "exported",len(dataset_data),"to",export_file

def customLSDeepCopy(ls):
    ls2=[]
    for i in ls: ls2.append(i)
    return ls2

def getAllListCombinationsLong(a):
    ls1 = ['a1','a2','a3']
    ls2 = ['b1','b2','b3']
    ls3 = ['c1','c2','c3']
    ls = ls1,ls2,ls3
    
    list_len_db={}
    for ls in a:
        list_len_db[len(x)]=[]
    print len(list_len_db), list_len_db;sys.exit()
    if len(list_len_db)==1 and 1 in list_len_db:
        ### Just simply combine non-complex data
        r=[]
        for i in a:
            r+=i
    else:
        #http://code.activestate.com/recipes/496807-list-of-all-combination-from-multiple-lists/
        r=[[]]
        for x in a:
            t = []
            for y in x:
                for i in r:
                    t.append(i+[y])
            r = t
    return r

def combineUniqueAllLists(files_to_merge,original_filename):
    headers =[]; all_keys={}; dataset_data={}; files=[]
    for filename in files_to_merge:
        print filename
        fn=filepath(filename); x=0; combined_data ={}; files.append(filename)
        if '/' in filename:
            file = string.split(filename,'/')[-1][:-4]
        else:
            file = string.split(filename,'\\')[-1][:-4]
        for line in open(fn,'rU').xreadlines():         
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0:
                if data[0]!='#':
                    x=1
                    try: t = t[1:]
                    except Exception: t = ['null']
                    for i in t:    
                        headers.append(i+'.'+file)
            if x==0:
                if data[0]!='#':
                    x=1;
                    headers+=t[1:] ###Occurs for the header line
                    headers+=['null']
            else: #elif 'FOXP1' in data or 'SLK' in data or 'MBD2' in data:
                key = t[0]
                try: values = t[1:]
                except Exception: values = ['null']
                try:
                    if original_filename in filename and len(original_filename)>0: key = t[1]; values = t[2:]
                except IndexError: print original_filename,filename,t;kill
                #key = string.replace(key,' ','')
                combined_data[key] = values
                if len(key)>0 and key != ' ':
                    try: all_keys[key] += 1
                    except KeyError: all_keys[key] = 1
        dataset_data[filename] = combined_data

    ###Add null values for key's in all_keys not in the list for each individual dataset
    combined_file_data = {}
    for filename in files:
        combined_data = dataset_data[filename]
        ###Determine the number of unique values for each key for each dataset
        null_values = []; i=0
        for key in combined_data: number_of_values = len(combined_data[key]); break
        while i<number_of_values: null_values.append('0'); i+=1
        for key in all_keys:
            include = 'yes'
            if combine_type == 'intersection':
                if all_keys[key]>(len(files_to_merge)-1): include = 'yes'
                else: include = 'no'
            if include == 'yes':
                try: values = combined_data[key]
                except KeyError: values = null_values
                try:
                    previous_val = combined_file_data[key]
                    new_val = previous_val + values
                    combined_file_data[key] = new_val
                except KeyError: combined_file_data[key] = values

    original_filename = string.replace(original_filename,'1.',  '1.AS-')          
    export_file = output_dir+'/MergedFiles.txt'
    fn=filepath(export_file);data = open(fn,'w')
    title = string.join(['uid']+headers,'\t')+'\n'; data.write(title)
    for key in combined_file_data:
        #if key == 'ENSG00000121067': print key,combined_file_data[key];kill
        new_key_data = string.split(key,'-'); new_key = new_key_data[0]
        values = string.join([key]+combined_file_data[key],'\t')+'\n'; data.write(values) ###removed [new_key]+ from string.join
    data.close()
    print "exported",len(dataset_data),"to",export_file
    
def getAllListCombinations(a):
    #http://www.saltycrane.com/blog/2011/11/find-all-combinations-set-lists-itertoolsproduct/
    """ Nice code to get all combinations of lists like in the above example, where each element from each list is represented only once """
    list_len_db={}
    for x in a:
        list_len_db[len(x)]=[]
    if len(list_len_db)==1 and 1 in list_len_db:
        ### Just simply combine non-complex data
        r=[]
        for i in a:
            r.append(i[0])
        return [r]
    else:
        return list(itertools.product(*a))
    
def joinFiles(files_to_merge,CombineType,unique_join,outputDir):
    """ Join multiple files into a single output file """
    global combine_type
    global permform_all_pairwise
    global output_dir
    output_dir = outputDir
    combine_type = string.lower(CombineType)
    permform_all_pairwise = 'yes'

    print 'combine type:',combine_type
    print 'join type:', unique_join
    #g = GrabFiles(); g.setdirectory(import_dir)
    #files_to_merge = g.searchdirectory('xyz') ###made this a term to excluded
    
    if unique_join:
        combineUniqueAllLists(files_to_merge,'')
    else:
        combineAllLists(files_to_merge,'')
        
    return output_dir+'/MergedFiles.txt'

if __name__ == '__main__':
    dirfile = unique
    includeColumns=-2
    includeColumns = False
    output_dir = filepath('output')
    combine_type = 'union'
    permform_all_pairwise = 'yes'
    print "Analysis Mode:"
    print "1) Batch Analysis"
    print "2) Single Output"
    inp = sys.stdin.readline(); inp = inp.strip()
    if inp == "1": batch_mode = 'yes'
    elif inp == "2": batch_mode = 'no'
    
    print "Combine Lists Using:"
    print "1) Grab Union"
    print "2) Grab Intersection"
    inp = sys.stdin.readline(); inp = inp.strip()
    if inp == "1": combine_type = 'union'
    elif inp == "2": combine_type = 'intersection'

    if batch_mode == 'yes': import_dir = '/batch/general_input'
    else: import_dir = '/input'
    g = GrabFiles(); g.setdirectory(import_dir)
    files_to_merge = g.searchdirectory('xyz') ###made this a term to excluded

    if batch_mode == 'yes':
        second_import_dir = '/batch/primary_input'
        g = GrabFiles(); g.setdirectory(second_import_dir)
        files_to_merge2 = g.searchdirectory('xyz') ###made this a term to excluded
        for file in files_to_merge2:
            temp_files_to_merge = customLSDeepCopy(files_to_merge)
            original_filename = string.split(file,'/'); original_filename = original_filename[-1]
            temp_files_to_merge.append(file)
            if '.' in file:
                combineAllLists(temp_files_to_merge,original_filename)
    else:
        combineAllLists(files_to_merge,'',includeColumns=includeColumns)
    print "Finished combining lists. Select return/enter to exit"; inp = sys.stdin.readline()
