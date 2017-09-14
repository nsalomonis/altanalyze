###GO-Elite
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

"""This module contains instructions importing existing over-representation analysis (ORA)
results, pruning redundant Gene Ontology paths, reading in command-line options and coordinating
all GO-Elite analysis and update functions found in other modules."""

try:
    try: import traceback
    except Exception: None
    import sys, string
    import os.path, platform
    import unique; reload(unique)
    import math
    import copy
    import time
    import shutil
    import webbrowser
    import gene_associations; reload(gene_associations)
    try: from import_scripts import OBO_import
    except Exception: pass
    import export
    import UI
    import mappfinder; reload(mappfinder)
    import datetime
    from visualization_scripts import WikiPathways_webservice
except Exception:
    print_out = "\nWarning!!! Critical Python incompatiblity encoutered.\nThis can occur if the users calls the GO-Elite "
    print_out += "python source code WITHIN A COMPILED version directory which results in critical conflicts between "
    print_out += "the compiled distributed binaries and the operating systems installed version of python. "
    print_out += "If this applies, either:\n(A) double-click on the executable GO-Elite file (GO-Elite or GO-Elite.exe) "
    print_out += "or\n(B) download the Source-Code version of GO-Elite and run this version (after installing any needed "
    print_out += "dependencies (see our Wiki or Documentation). Otherwise, contact us at: altanalyze@gmail.com\n\n"
    print_out += "Installation Wiki: http://code.google.com/p/go-elite/wiki/Installation\n"
    print print_out
    
    try:
        ### Create a log report of this
        import unique
        try: log_file = unique.filepath('GO-Elite_error-report.log')
        except Exception: log_file = unique.filepath('/GO-Elite_error-report.log')
        log_report = open(log_file,'a')
        log_report.write(print_out)
        try: log_report.write(traceback.format_exc())
        except Exception: None
        log_report.close()
        print 'See log file for additional technical details:',log_file
        ### Open this file
        if os.name == 'nt':
            try: os.startfile('"'+log_file+'"')
            except Exception:  os.system('open "'+log_file+'"')
        elif 'darwin' in sys.platform: os.system('open "'+log_file+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+log_file+'/"')   
    except Exception: None
    sys.exit()

debug_mode = 'no'
Tkinter_failure = False
start_time = time.time()

try: command_args = string.join(sys.argv,' ')
except Exception: command_args = ''
if len(sys.argv[1:])>1 and '-' in command_args:
    use_Tkinter = 'no'
else:
    try:
        import Tkinter
        from Tkinter import *
        from visualization_scripts import PmwFreeze
        use_Tkinter = 'yes'
    except ImportError:
        use_Tkinter = 'yes'
        Tkinter_failure = True
        
###### File Import Functions ######
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

def read_directory_extended(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv" or entry[-4:] == ".tab": dir_list2.append(entry)
    return dir_list2

def refDir():
    reference_dir=unique.refDir()
    return reference_dir

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database: list = unique.unique(database[key]); list.sort(); db1[key] = list
    return db1

###### Classes ######
class GrabFiles:
    def setdirectory(self,value):
        self.data = value
    def display(self):
        print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        files,file_dir,file = gene_associations.getDirectoryFiles(self.data,str(search_term))
        if len(file)<1: print search_term,'not found'
        return file_dir
    def searchdirectory_start(self,search_term):
        #self is an instance while self.data is the value of the instance
        files,file_dir,file = gene_associations.getDirectoryFiles(self.data,str(search_term))
        match = ''
        for file in files:
            split_strs = string.split(file,'.')
            split_strs = string.split(split_strs[0],'/')
            if search_term  == split_strs[-1]: match = file
        return match
    
def getDirectoryFiles(import_dir, search_term):
    exact_file = ''; exact_file_dir=''
    dir_list = read_directory_extended(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        if (':' in import_dir) or ('/Users/' == import_dir[:7]) or ('Linux' in platform.system()): affy_data_dir = import_dir+'/'+data
        else: affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir: exact_file = affy_data_dir
    return exact_file

def getAllDirectoryFiles(import_dir, search_term):
    exact_files = []
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        if ':' in import_dir or ('/Users/' == import_dir[:7]) or ('Linux' in platform.system()): affy_data_dir = import_dir+'/'+data
        else: affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir: exact_files.append(affy_data_dir)
    return exact_files

###### GO-Elite Functions ######

def import_path_index():
    global path_id_to_goid;  global full_path_db; global all_nested_paths
    try:  built_go_paths, path_goid_db, path_dictionary = OBO_import.remoteImportOntologyTree(ontology_type)
    except IOError:
        ### This exception was added in version 1.2 and replaces the code in OBO_import.buildNestedOntologyTree which resets the OBO version to 0/0/00 and re-runs (see unlisted_variable = kill)
        print_out = "Unknown error encountered during Gene Ontology Tree import.\nPlease report to altanalyze@gmail.com if this error re-occurs."
        try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy(); sys.exit()
        except Exception: print print_out; sys.exit()
    print 'Imported',ontology_type,'tree-structure for pruning'
    full_path_db = built_go_paths
    path_id_to_goid = path_goid_db
    #all_nested_paths = path_dictionary 
    #print len(full_path_db), 'GO paths imported'
    
class MAPPFinderResults:
    def __init__(self,goid,go_name,zscore,num_changed,onto_type,permute_p,result_line):
        self._goid = goid; self._zscore = zscore; self._num_changed = num_changed
        self._onto_type = onto_type; self._permute_p = permute_p; self._go_name= go_name
        self._result_line = result_line
    def GOID(self): return self._goid
    def GOName(self): return self._go_name
    def ZScore(self): return self._zscore
    def NumChanged(self): return self._num_changed
    def OntoType(self): return self._onto_type
    def PermuteP(self): return self._permute_p
    def ResultLine(self): return self._result_line
    def Report(self):
        output = self.GOID()+'|'+self.GOName()
        return output
    def __repr__(self): return self.Report()

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def moveMAPPFinderFiles(input_dir):
    ###Create new archive directory
    import datetime
    dir_list = read_directory_extended(input_dir)
    input_dir2 = input_dir
    if len(dir_list)>0: ###if there are recently run files in the directory
        today = str(datetime.date.today()); today = string.split(today,'-'); today = today[0]+''+today[1]+''+today[2]
        time_stamp = string.replace(time.ctime(),':','')
        time_stamp = string.replace(time_stamp,'  ',' ')
        time_stamp = string.split(time_stamp,' ') ###Use a time-stamp as the output dir (minus the day)
        time_stamp = today+'-'+time_stamp[3]
        if 'archived-' not in input_dir2: ### We don't want users to select an archived directory and have the results moved to a child archive directory.
            archive_dir = input_dir2 +'/archived-'+time_stamp
            fn = filepath(archive_dir)
            try: os.mkdir(fn) ###create new directory
            except Exception:
                try: archive_dir = archive_dir[1:]; fn = filepath(archive_dir); os.mkdir(fn)
                except Exception: null = [] ### directory already exists
            ###Read all files in old directory and copy to new directory
            m = GrabFiles(); m.setdirectory(input_dir)
            #print 'Moving MAPPFinder results to the archived directory',archive_dir
            for mappfinder_input in dir_list: 
                mappfinder_input_dir = m.searchdirectory(mappfinder_input)
                fn = filepath(mappfinder_input_dir)
                fn2 = string.split(fn,'/')
                destination_dir = string.join(fn2[:-1]+['/archived-'+time_stamp+'/']+[fn2[-1]],'/')
                #destination_dir = string.replace(fn,'/'+mappfinder_dir+'/','/'+mappfinder_dir+'/archived-'+time_stamp+'/');
                try: shutil.copyfile(fn, destination_dir)
                except Exception:
                    print_out = "WARNING!!! Unable to move ORA results to an archived directory."
                    try: UI.WarningWindow(print_out,' Exit ')
                    except Exception: print print_out
                    sys.exit()
                proceed = 'no'
                while proceed == 'no':
                    try: os.remove(fn); proceed = 'yes'
                    except Exception:
                        print 'Tried to move the file',mappfinder_input,'to an archived folder, but it is currently open.'
                        print 'Please close this file and hit return or quit GO-Elite'
                        inp = sys.stdin.readline()

def checkPathwayType(filename):
    type='GeneSet'
    ontology_type = string.split(filename,'-')[-1][:-4]
    if ontology_type == 'local': ontology_type = 'WikiPathways'
    if ontology_type == 'GO': ontology_type = 'GeneOntology'
    fn=filepath(filename); x=20; y=0
    for line in open(fn,'rU').readlines():
        y+=1
        if y<x:
            if 'Ontology-ID' in line: type = 'Ontology'
    return ontology_type,type

def importMAPPFinderResults(filename):
    zscore_changed_path_db = {}; go_full = {}; zscore_goid = {}
    global dummy_db; dummy_db = {}; global file_headers; file_headers = []
    global full_go_name_db; full_go_name_db = {}; filtered_mapps=[]; filtered_mapp_list = []
    run_mod = ''; run_source=''; gene_identifier_file = ''
    fn=filepath(filename); x=0; y=0
    #print "Importing MAPPFinder data for:",filename
    for line in open(fn,'rU').readlines():
        data = cleanUpLine(line)
        if x == 0 and y == 0:
            if 'Ontology-ID' in data:
                x = 1; go_titles = data + '\n'; input_file_type = 'Ontology'
            elif 'Gene-Set Name' in data:
                y = 1; mapp_titles = data + '\n'; input_file_type = 'GeneSet'
            elif data in species_list: ###Check species from results file
                if data in species_list: species = data; species_code = species_codes[species]
            elif 'source identifiers supplied' in data:
                space_delimited = string.split(data,' ')
                gene_input_file_data = string.split(data,'supplied in the input file:')
                putative_source_type = space_delimited[1]  ###Extract out uid from GO/MAPP results file
                if putative_source_type in source_types: run_source = putative_source_type
                elif putative_source_type in mod_types: run_source = putative_source_type
                if len(gene_input_file_data)>1: gene_identifier_file = gene_input_file_data[-1]
            elif 'source identifiers linked to' in data:
                space_delimited = string.split(data,' ')
                putative_mod_type = space_delimited[-2] ###Extract out MOD from GO/MAPP results file
                if putative_mod_type in mod_types: run_mod = putative_mod_type
            elif 'Probes linked to a ' in data:
                space_delimited = string.split(data,' ')
                putative_mod_type = space_delimited[-2] ###Extract out MOD from GO/MAPP results file (derived from GenMAPP MAPPFinder results)
                if putative_mod_type in mod_types: run_mod = putative_mod_type
                else:
                    print_out = "WARNING!!! The MOD "+putative_mod_type+"\nused by GenMAPP MAPPFinder, is not\navailable in GO-Elite. Errors may\nresult in deriving propper gene\nassociations since a different MOD has to be\nused."
                    try: UI.WarningWindow(print_out,' Exit ')
                    except Exception: print print_out
        elif x == 1:
            try: goid, go_name, go_type, number_changed, number_measured, number_go, percent_changed, percent_present, z_score, permutep, adjustedp = string.split(data,'\t')
            except ValueError:
                t = string.split(data,'\t')
                print_out = "WARNING...GeneOntology results file has too few or too many columns.\nExpected 16 columns, but read"+ str(len(t))
                try: UI.WarningWindow(print_out,' Exit ')
                except Exception:
                    print print_out; print t; print 'Please correct and re-run'
                sys.exit()                
            if goid != 'GO': ###Exclude the top level of the GO tree (not an ID number MAPPFinder output)
                #go_id_full = 'GO:'+goid
                ###Check to make sure that the GOID is the right length (usually proceeded by zeros which can be removed)
                if len(goid)<7:
                    extended_goid=goid
                    while len(extended_goid)< 7: extended_goid = '0'+extended_goid
                    if extended_goid in full_path_db: goid = extended_goid
                if goid in full_path_db:
                    for path in full_path_db[goid]:
                        #goid = int(goid); path = path_data.PathTuple() ###Path is a tuple of integers
                        z_score = float(z_score); permutep = float(permutep); number_changed = int(number_changed)
                        #path_list = string.split(path,'.')#; path_list=[]; full_tree_index[path] = path_list ###used to created nested gene associations to GO-elite terms
                        full_go_name_db[goid] = go_name
                        #enrichment_type
                        if number_changed > change_threshold and permutep < p_val_threshold and int(number_go) <= max_member_count:
                            if (z_score > z_threshold and enrichment_type=='ORA') or (z_score < (z_threshold*-1) and enrichment_type=='URA'):
                                #tree_index[path] = path_list; zscore_changed_path_db[path] = z_score, number_changed; zscore_goid[goid] = z_score, go_type; go_terms[path] = goid
                                z = MAPPFinderResults(goid,go_name,z_score,number_changed,input_file_type,permutep,data)
                                zscore_changed_path_db[path] = z; zscore_goid[goid] = z
                                go_full[goid]= data + '\n'
            ### for all entries, save a dummy entry.  This is needed if we want to re-examine in MAPPFinder but don't want to visualize non-GO-elite selected terms
            dummy_data = str(goid) +"\t"+ go_name +"\t"+ go_type
            dummy_data = dummy_data+ "\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"1"+"\t"+"1" +"\n"
            dummy_db[goid] = dummy_data    
        elif y == 1:
                ###CODE USED TO PROCESS LOCAL RESULT FILES FOR MAPPS
                try:
                    try:
                        t = string.split(data,'\t')
                        mapp_name = t[0]; num_changed = t[1]; zscore = t[6]; permutep = t[7]; num_on_mapp = t[3]
                        num_changed = int(num_changed); zscore = float(zscore); permutep = float(permutep)
                        #print [mapp_name, num_changed, zscore,permutep]
                        if num_changed > change_threshold and permutep < p_val_threshold and int(num_on_mapp) <= max_member_count:
                            if (zscore > z_threshold and enrichment_type=='ORA') or (zscore < (z_threshold*-1) and enrichment_type=='URA'):
                                filtered_mapps.append((zscore, line, mapp_name))
                                #zscore_goid[mapp_name] = zscore, 'MAPP'
                                z = MAPPFinderResults(mapp_name,mapp_name,zscore,num_changed,input_file_type,permutep,data)
                                zscore_goid[mapp_name] = z
                                filtered_mapp_list.append(mapp_name) ###Use this for downstream analyses
                    except ValueError: continue
                except IndexError: continue
    if x == 1:
        try: species_code = species_code ###Test to see if there is a proper species name found
        except UnboundLocalError:
            #print 'No proper species name found in GO results. Please change and re-run'
            #inp = sys.stdin.readline(); sys.exit()
            ###Pick arbitrary species code
            for species in species_codes: species_code = species_codes[species]
        return run_mod,run_source,zscore_changed_path_db, go_full, go_titles, zscore_goid, input_file_type, gene_identifier_file, species_code
    if y == 1:
        filtered_mapps.sort(); filtered_mapps.reverse()
        return run_mod,run_source,filtered_mapp_list, filtered_mapps, mapp_titles, zscore_goid, input_file_type, gene_identifier_file, species_code

def buildFullTreePath(go_index_list):
    """Find common parent nodes among children and parents and collapse into a single branch node"""
    global path_dictionary; path_dictionary = {}; path_dictionary_filtered={}
    for tuple_index in go_index_list: path_dictionary[tuple_index] = [] ###store the original top level path-indeces
    for tuple_index in go_index_list:
        path_len = len(tuple_index); i=-1
        while path_len+i > 0:
            parent_path = tuple_index[:i]
            if parent_path in path_dictionary:
                path_dictionary[parent_path].append(tuple_index) ###if a parent path of a
            i-=1
    for parent_path in path_dictionary:
        if len(path_dictionary[parent_path])>0:
            path_dictionary[parent_path].sort()
            path_dictionary_filtered[parent_path] = path_dictionary[parent_path]
        ### All GOIDs NOT in path_dictionary_filtered (no significant parents) are added back later 
    return path_dictionary_filtered

def buildOrderedTree(go_index_list):
    """Find common parent nodes among children and parents and collapse into a single branch node"""
    path_dictionary = {}; path_dictionary_filtered={}; organized_tree_db={}
    original_path_dictionary={}
    for tuple_index in go_index_list:
        path_dictionary[tuple_index] = [] ###store the original top level path-indeces
        original_path_dictionary[tuple_index] = []
    for tuple_index in go_index_list:
        path_len = len(tuple_index); i=-1
        while path_len+i > 0:
            parent_path = tuple_index[:i]
            if parent_path in path_dictionary:
                path_dictionary[tuple_index].append(parent_path)
                original_path_dictionary[parent_path].append(tuple_index)
            i-=1
    for parent_path in original_path_dictionary:
        if len(original_path_dictionary[parent_path])>0:
            original_path_dictionary[parent_path].sort()
            path_dictionary_filtered[parent_path] = original_path_dictionary[parent_path]
        ### All GOIDs NOT in path_dictionary_filtered (no significant parents) are added back later

    keys_to_delete={}   
    for child in path_dictionary:
        if len(path_dictionary[child])>0:
            path_dictionary[child].sort()
            new_path = path_dictionary[child][1:]
            for path in new_path:
                keys_to_delete[path] = []
                
    for path in keys_to_delete:
        try: del path_dictionary[path]
        except Exception: null=[]

    parent_tree_db={}
    for child in path_dictionary:
        if len(path_dictionary[child])>0:
            path_dictionary[child].sort()
            new_path = path_dictionary[child][1:]+[child]
            #print new_path, path_dictionary[child][0], path_dictionary[child];kill
            try: parent_tree_db[path_dictionary[child][0]].append(new_path)
            except Exception: parent_tree_db[path_dictionary[child][0]] = [new_path]
        else:
            organized_tree_db['node',child]=['null']
            
    for path in keys_to_delete:
        try: del parent_tree_db[path]
        except Exception: null=[]

    for path in parent_tree_db:
        del organized_tree_db['node',path]

    finished = 'no'
    while finished == 'no':
        parent_tree_db2={}; count=0
        for parent_path in parent_tree_db:
            if len(parent_tree_db[parent_path])>1:
                for childset in parent_tree_db[parent_path]:
                    top_child = childset[0]
                    if len(childset)>1:
                        if 'node' not in parent_path:
                            new_parent = 'node',parent_path,top_child
                        else: new_parent = tuple(list(parent_path)+[top_child])
                        try: parent_tree_db2[new_parent].append(childset[1:])
                        except Exception: parent_tree_db2[new_parent] = [childset[1:]]
                        count+=1
                    elif len(childset)==1:
                        if 'node' not in parent_path:
                            new_parent = 'node',parent_path,top_child
                        else: new_parent = tuple(list(parent_path)+[top_child])
                        organized_tree_db[new_parent] = ['null']
            else:
                childset = parent_tree_db[parent_path][0]
                if 'node' not in parent_path:
                    new_parent = 'node',parent_path
                else: new_parent = parent_path
                #if len(childset)>1: print new_parent, childset;kill
                organized_tree_db[new_parent] = childset
        if count == 0: finished = 'yes'
        parent_tree_db = parent_tree_db2

    possible_sibling_paths={}; siblings_reverse_lookup={}
    for path in organized_tree_db:
        if len(path)>2 and organized_tree_db[path]== ['null']:
            try: possible_sibling_paths[path[-1][:-1]].append(path)
            except Exception: possible_sibling_paths[path[-1][:-1]] = [path]
            ### Since some deepest children will have siblings, along with cousins on different branches, we must exclude these
            for node in path[1:-1]:
                try: siblings_reverse_lookup[node].append(path)
                except KeyError: siblings_reverse_lookup[node] = [path]

    for node in siblings_reverse_lookup:
        try:
            if len(siblings_reverse_lookup[node])!=len(possible_sibling_paths[node]):
                try: del possible_sibling_paths[node]
                except Exception: null=[]
        except Exception: null=[]

    for common_path in possible_sibling_paths:
        if len(possible_sibling_paths[common_path])>1:
            sibling_paths=[]
            for key in possible_sibling_paths[common_path]:
                new_key = tuple(['sibling1']+list(key[1:]))
                organized_tree_db[new_key] =  ['siblings_type1']
                sibling_paths.append(key[-1])
                del organized_tree_db[key]
            sibling_paths = unique.unique(sibling_paths)
            sibling_paths.sort()
            new_key = tuple(['sibling2']+list(key[1:-1]))
            organized_tree_db[new_key]=sibling_paths
    """
    for i in organized_tree_db:
        print i, organized_tree_db[i]
    kill"""
    return organized_tree_db,path_dictionary_filtered

def link_score_to_all_paths(all_paths,zscore_changed_path_db):
    temp_list = []
    parent_highest_score = {}
    for entry in all_paths:
        if entry in zscore_changed_path_db:
            if filter_method == 'z-score': new_entry = zscore_changed_path_db[entry].ZScore(), entry
            elif filter_method == 'gene number':new_entry = zscore_changed_path_db[entry].NumChanged(), entry
            elif filter_method == 'combination':
                z_score_val = zscore_changed_path_db[entry].ZScore()
                gene_num = math.log(zscore_changed_path_db[entry].NumChanged(),2)
                new_entry = gene_num*z_score_val, entry
            else:
                print 'error, no filter method selected!!'
                break
        for item in all_paths[entry]:
            if filter_method == 'z-score':
                new_item = zscore_changed_path_db[item].ZScore(), item
                temp_list.append(new_item)
            elif filter_method == 'gene number':
                new_item = zscore_changed_path_db[item].NumChanged(), item
                temp_list.append(new_item)
            elif filter_method == 'combination':
                z_score_val = zscore_changed_path_db[item].ZScore()
                gene_num = math.log(zscore_changed_path_db[item].NumChanged(),2)
                new_item = gene_num*z_score_val, item
                temp_list.append(new_item)
            else:
                print 'error, no filter method selected!!'
                break
        max_child_score = max(temp_list)[0]
        max_child_path = max(temp_list)[1]
        parent_score = new_entry[0]
        parent_path = new_entry[1]
        if max_child_score <= parent_score:
            parent_highest_score[parent_path] = all_paths[entry]   #include the elite parent and it's children
            #print "parent_path",parent_path,parent_score,"max_child_path", max_child_path, max_child_score
        temp_list=[]
    """for entry in parent_highest_score:
        print entry,':',parent_highest_score[entry]"""
    #print "Number of parents > child in tree:",len(parent_highest_score)
    return parent_highest_score

def calculate_score_for_children(tree, zscore_changed_path_db):
    temp_list = []
    child_highest_score = {}
    for entry in tree:  ###The tree is a key with one or more nodes representing a branch, chewed back, with or without children
        for item in tree[entry]: #get the children of the tree
            if item in zscore_changed_path_db:
                if filter_method == 'z-score':
                    new_item = zscore_changed_path_db[item].ZScore(), item
                    temp_list.append(new_item)
                elif filter_method == 'gene number':
                    new_item = zscore_changed_path_db[item].NumChanged(), item
                    temp_list.append(new_item)
                elif filter_method == 'combination':
                    z_score_val = zscore_changed_path_db[item].ZScore()
                    gene_num = math.log(zscore_changed_path_db[item].NumChanged(),2)
                    new_item = gene_num*z_score_val, item
                    temp_list.append(new_item)
                    #print new_item,z_score_val
                else:
                    print 'error, no filter method selected!!'
                    break
            """elif item == 'siblings_type1': #this will only occur if an parent had multiple children, none of them with nested children
                 if filter_method == 'z-score':
                    parent_item = entry[-1]
                    new_item = zscore_changed_path_db[parent_item][0], parent_item
                    temp_list.append(new_item)
                elif filter_method == 'gene number':
                    new_item = zscore_changed_path_db[parent_item][1], parent_item
                    temp_list.append(new_item)
                elif filter_method == 'combination':
                    z_score_val = zscore_changed_path_db[parent_item][0]
                    gene_num = math.log(zscore_changed_path_db[parent_item][1],2)
                    new_item = gene_num*z_score_val, parent_item
                    temp_list.append(new_item)"""
        #print new_item,z_score_val               
        if len(temp_list)> 0:  #if there is at least one nested go_path
            max_child_score = max(temp_list)[0]
            max_child_path = max(temp_list)[1]
            child_highest_score[entry]=max_child_path
            temp_list = []
    """for entry in child_highest_score:
        print entry,':',child_highest_score[entry]"""
    return child_highest_score

def collapse_tree(parent_highest_score,child_highest_score,tree):
    collapsed = {}
    collapsed_parents = {}
    collapsed_children = {}
    for entry in tree:
        count=0 #see next comment
        for item in entry:
            if item in parent_highest_score and count==0: #ensures that only the top parent is added
                collapsed_parents[entry] = item, parent_highest_score[item] #include the children of the elite parent to exclude later
                ### aka elite_parents[parent_paths] = elite_parent_path, non_elite_children
                count +=1  #see last comment
                break
    for entry in tree:
        if entry not in collapsed_parents:  ###if the highest score does not occur for one of the parent terms in the key of the tree-branch
            if entry in child_highest_score: ###Entry is the parent path stored and the value is the child path with the max score (greater than parent)
                max_child_path = child_highest_score[entry]
                collapsed_children[entry] = max_child_path, tree[entry] #include the children to exclude from other entries later (includes the max_child_path though)
    """for entry in collapsed_children:
        print entry,':',collapsed_children[entry]"""
    temp_list = []
    for entry in collapsed_parents:
        node = collapsed_parents[entry][0]
        children = collapsed_parents[entry][1]
        collapsed['parent',node] = children
        """if entry in tree:
            for item in entry:
                temp_list.append(item)
            for item in tree[entry]:
                temp_list.append(item)
            temp_list.sort()
            temp_list2 = unique.unique(temp_list)
            collapsed['parent',node] = temp_list2
            temp_list=[]"""
    for entry in collapsed_children:
        #previously we included the blocked out code which also
        #would later exclude entries that contain key entry items in the elite entries
        node = collapsed_children[entry][0]
        siblings = collapsed_children[entry][1]
        collapsed['child',node] = siblings
        """if entry in tree:
            for item in entry: temp_list.append(item)
            for item in tree[entry]: temp_list.append(item)
            temp_list.sort()
            temp_list2 = unique.unique(temp_list)
            collapsed['child',node] = temp_list2
            temp_list=[]"""
    for entry in tree:
        nested_path = tree[entry]
        if nested_path == ['null']:
            collapsed['unique',entry[-1]] = 'unique',entry[-1]
    #print "The length of the collapsed list is:",len(collapsed)
    ###This code is used to check whether or not a tree with multiple siblings listed as children was added to collapsed because the parent term(s) had a higher score
    ###if not, we must individually add the child sibling terms
    for entry in tree:
        if tree[entry] == ['siblings_type1']:
            parent_root_node = 'parent',entry[-2]
            if parent_root_node not in collapsed:
                #use 'child' since we have to over-write the child entry created if the 'parent' was'nt more significant (thus syblings would eroneously be considered children)
                 collapsed['child',entry[-1]] = 'node',entry[-1] 
    #for entry in collapsed:
        #print entry,collapsed[entry]
    return collapsed

def link_goid(tree,zscore_changed_path_db,all_paths):
    """Starting with the collapsed tree, convert to goids"""
    temp1=[]; temp2=[]; new_tree = {}; count = 0

    for entry in tree:
        for top_node in entry:
            try: goid = zscore_changed_path_db[top_node].GOID(); temp1.append(goid)
            except KeyError: not_valid_goid = ''
        for child_node in tree[entry]:
            try: goid= zscore_changed_path_db[child_node].GOID(); temp2.append(goid)
            except KeyError: not_valid_goid = ''
        temp2.sort()
        new_top = tuple(temp1) #one elite go-term in goid form
        new_child = temp2   #list of children (value items) in dictionary, in goid form
        count +=1; temp_list=[]
        if new_top in new_tree:  #if the new dictionary, new_tree already has this key added
            top_node = new_top[0]
            for node in new_tree[new_top]:
                if node != top_node: temp_list.append(node)
            for node in new_child:
                if node != top_node: temp_list.append(node)
            temp_list.sort(); new_list = unique.unique(temp_list)
            new_tree[new_top] = new_list; temp_list=[]
        else:
            #print new_child[0], new_top[0]
            if new_child[0] == new_top[0]:  #if the first dictionary (which should be the only) value equals the key
                new_tree[new_top] = []
            else:  ###note: if two parents, both more significant than all children, share a child, both will still be reported. This is fine, since we don't know the parent's relationship
                temp_list2 = []
                for node in new_child:  
                    if node != new_top[0]: #get rid of values that equal the key!!!
                        temp_list2.append(node)
                temp_list2.sort()
                new_child = temp_list2
                new_tree[new_top] = new_child #this should be the unique term
                #print new_top, new_child
        temp1 = []; temp2 = []
    ###remove entries that are parents or children of the most 'elite' entry, in dictionary
    ###note: parents of the entry shouldn't be selected, since collapsing chooses
    ###selects the highest scoring node.
    if exclude_related == 'yes':
        exclude = {}
        for entry in new_tree:
            for item in new_tree[entry]:
                exclude[item] = ''
                #print item, entry, new_tree[entry]
        new_tree2={}               
        for entry in new_tree:
            if entry[0] not in exclude:
                new_tree2[entry[0]] = ''
        new_tree = new_tree2
    
    #for entry in tree: print entry,':',tree[entry]
    #for entry in new_tree: print entry,':',new_tree[entry] 
    
    #print "Length of input tree/output tree",len(tree),'/',len(new_tree),'count:',count
    return new_tree
                            
def exportGOResults(go_full,go_titles,collapsed_go_tree,zscore_goid,go_gene_annotation_db,go_values_db,value_headers,goids_with_redundant_genes):
    reference_dir = refDir()
    if len(go_elite_output_folder) == 0: output_folder = "output"; ref_dir = reference_dir +'/'+output_folder
    else: output_folder = go_elite_output_folder; ref_dir = go_elite_output_folder
    output_results = output_folder+'/'+mappfinder_input[0:-4]+ '_' +filter_method +'_elite.txt'
    if ontology_type == 'WikiPathways':
        output_results = string.replace(output_results, 'WikiPathways','local') ### Makes the output filename compatible with GenMAPP-CS plugin filenames
    if ontology_type == 'GO':
        output_results = string.replace(output_results, 'GeneOntology','GO') ### Makes the output filename compatible with GenMAPP-CS plugin filenames

    entries_added=[]  #keep track of added GOIDs - later, add back un-added for MAPPFinder examination
    print_out = 'Results file: '+output_results+ '\nis open...can not re-write.\nPlease close file and select "OK".'
    try: raw = export.ExportFile(output_results)
    except Exception:
        try: UI.WarningWindow(print_out,' OK ')
        except Exception:
            print print_out; print 'Please correct (hit return to continue)'; inp = sys.stdin.readline()
    if export_user_summary_values == 'yes' and len(value_headers)>0:
        value_headers2=[]; stdev_headers=[]
        for i in value_headers: value_headers2.append('AVG-'+i); stdev_headers.append('STDEV-'+i)
        value_headers = '\t'+string.join(value_headers2+stdev_headers,'\t')
    else: value_headers = ''
    if len(go_gene_annotation_db)>0: go_titles = go_titles[:-1]+'\t'+'redundant with terms'+'\t'+'inverse redundant'+'\t'+'gene symbols'+value_headers+'\n' ###If re-outputing results after building gene_associations
    raw.write(go_titles)
    combined_results[output_results] = [((ontology_type,'Ontology'),mappfinder_input,go_titles)]
    
    #append data to a list, to sort it by go_cateogry and z-score
    collapsed_go_list = []; collapsed_go_list2 = []; goids_redundant_with_others={}
    for entry in collapsed_go_tree:
        goid = entry
        entries_added.append(goid)
        if goid in zscore_goid:
            z_score_val = zscore_goid[goid].ZScore()
            go_type = zscore_goid[goid].OntoType()
            if sort_only_by_zscore == 'yes': info = z_score_val,z_score_val,goid
            else: info = go_type,z_score_val,goid
            collapsed_go_list.append(info); collapsed_go_list2.append(goid)
    collapsed_go_list.sort()
    collapsed_go_list.reverse()
    for (go_type,z_score_val,goid) in collapsed_go_list:
        if goid in go_full:
            data = go_full[goid]
            if len(go_gene_annotation_db)>0:
                symbol_ls = []; goid_vals = ''
                if goid in go_gene_annotation_db:
                    for s in go_gene_annotation_db[goid]:
                        if len(s.Symbol())>0: symbol_ls.append(s.Symbol()) #s.GeneID()
                        else: symbol_ls.append(s.GeneID())
                if goid in go_values_db:
                    goid_vals = string.join(go_values_db[goid][0]+go_values_db[goid][1],'\t') ###mean of values from input file, for each GO term
                symbol_ls = unique.unique(symbol_ls); symbol_ls.sort()
                symbols = string.join(symbol_ls,'|')
                try:
                    rr = goids_with_redundant_genes[goid]
                    redundant_with_terms = rr.RedundantNames(); inverse_redundant = rr.InverseNames()
                    if len(redundant_with_terms)>1: goids_redundant_with_others[goid]=[]
                except KeyError: redundant_with_terms = ' '; inverse_redundant =' '
                if export_user_summary_values == 'yes' and len(goid_vals)>0: goid_vals = '\t'+goid_vals
                else: goid_vals = ''
                data = data[:-1]+'\t'+redundant_with_terms+'\t'+inverse_redundant+'\t'+symbols+goid_vals+'\n'  ###If re-outputing results after building gene_associations
            raw.write(data)
            combined_results[output_results].append(((ontology_type,'Ontology'),mappfinder_input,data))
    raw.close()

    combined_results[output_results].append((('',''),'',''))
    summary_data_db['redundant_go_term_count'] = len(goids_redundant_with_others)
    #print "New MAPPFinder Elite file", output_results, "written...."
    return collapsed_go_list2

def exportLocalResults(filtered_mapps,mapp_titles,mapp_gene_annotation_db,mapp_values_db,mapp_value_headers,mapps_with_redundant_genes):
    if len(go_elite_output_folder) == 0: output_folder = "output"
    else: output_folder = go_elite_output_folder

    output_results = output_folder+'/'+mappfinder_input[0:-4]+ '_' +filter_method + '_elite.txt'
    try: raw = export.ExportFile(output_results)
    except Exception:
        try: UI.WarningWindow(print_out,' OK ')
        except Exception:
            print print_out; print 'Please correct (hit return to continue)'; inp = sys.stdin.readline()

    if export_user_summary_values == 'yes' and len(mapp_value_headers)>0:
        mapp_value_headers2=[]; stdev_headers=[]
        for i in mapp_value_headers: mapp_value_headers2.append('AVG-'+i); stdev_headers.append('STDEV-'+i)
        mapp_value_headers = '\t'+string.join(mapp_value_headers2+stdev_headers,'\t')
    else: mapp_value_headers = ''
    if len(mapp_gene_annotation_db)>0: mapp_titles = mapp_titles[:-1]+'\t'+'redundant with terms'+'\t'+'inverse redundant'+'\t'+'gene symbols'+mapp_value_headers+'\n'
    raw.write(mapp_titles)
    combined_results[output_results] = [((ontology_type,'GeneSet'),mappfinder_input,mapp_titles)]

    filtered_mapp_list = []; mapps_redundant_with_others={}
    for (zscore,line,mapp_name) in filtered_mapps:
        if len(mapp_gene_annotation_db)>0:
            symbol_ls = []; mapp_vals=''
            if mapp_name in mapp_gene_annotation_db:
                for s in mapp_gene_annotation_db[mapp_name]:
                    if len(s.Symbol())>0: symbol_ls.append(s.Symbol()) #s.GeneID()
                    else: symbol_ls.append(s.GeneID())
                symbol_ls = unique.unique(symbol_ls); symbol_ls.sort()
                symbols = string.join(symbol_ls,'|')
                if mapp_name in mapp_values_db:
                    mapp_vals = string.join(mapp_values_db[mapp_name][0]+mapp_values_db[mapp_name][1],'\t') ###mean of values from input file, for each MAPP
                try:
                    rr = mapps_with_redundant_genes[mapp_name]
                    redundant_with_terms = rr.RedundantNames(); inverse_redundant = rr.InverseNames()
                    if len(redundant_with_terms)>1: mapps_redundant_with_others[mapp_name]=[]
                except KeyError: redundant_with_terms = ' '; inverse_redundant = ' '

                if export_user_summary_values == 'yes' and len(mapp_vals)>0: mapp_vals = '\t'+mapp_vals
                else: mapp_vals = ''                
                line = line[:-1]+'\t'+redundant_with_terms+'\t'+inverse_redundant+'\t'+symbols+mapp_vals+'\n'  ###If re-outputing results after building gene_associations
        raw.write(line)
        combined_results[output_results].append(((ontology_type,'GeneSet'),mappfinder_input,line))
    raw.close()
    combined_results[output_results].append((('',''),'',''))
    #print "Local Filtered MAPPFinder file", output_results, "written...."
    summary_data_db['redundant_mapp_term_count'] = len(mapps_redundant_with_others)
    
def importORASimple(ora_dir,elite_pathway_db,file_type):
    """ This function imports pathway data for any elite pathway from the unfiltered results """
    summary_results_db={}
    ontology_name = file_type
    if file_type == 'GeneOntology': file_type = 'GO.'
    if file_type == 'WikiPathways': file_type = 'local'
    dir_list = read_directory(ora_dir)
    for file in dir_list:
        if '.txt' in file and file_type in file:
            fn=filepath(ora_dir+'/'+file)
            for line in open(fn,'rU').readlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                try:
                    ### Applies to Ontology ORA files
                    pathway_id = t[0]; pathway_name = t[1]; num_changed = t[3]; num_measured = t[4]; zscore = t[8]; permutep = t[9]; pathway_type = t[2]
                except IndexError:
                    ### Applies to Local pathway files
                    try:
                        pathway_name = t[0]; num_changed = t[1]; num_measured = t[2]; zscore = t[6]; permutep = t[7]; pathway_id = pathway_name; pathway_type = 'GeneSet'
                    except Exception: pathway_name='null'
                    
                ### If the pathway/ontology term was considered a pruned term (elite)

                if (pathway_name,ontology_name) in elite_pathway_db:
                    try:
                        int(num_changed)
                        ### Encode this data immediately for combination with the other input files and for export
                        key = (pathway_id,pathway_name,ontology_name,pathway_type)
                        values = [num_changed,num_measured,zscore,permutep]
                        try: summary_results_db[pathway_name,key].append([file,values])
                        except Exception: summary_results_db[pathway_name,key] = [[file,values]]
                    except Exception: null=[]
    return summary_results_db

def outputOverlappingResults(combined_results,ora_dir):
    """For any elite term/pathway from a user analysis, output all associated original ORA data"""
    
    output_folder = go_elite_output_folder
    output_results = output_folder+'/'+'overlapping-results_' +filter_method + '_elite.txt'
    proceed = 'no'
    while proceed == 'no':
        try: raw = export.ExportFile(output_results); proceed = 'yes'
        except Exception:
            print_out = output_results, '\nis open. Please close and select "Continue".'
            try: UI.WarningWindow(print_out,' OK ')
            except Exception:
                print print_out; print 'Please correct (hit return to continue)'; inp = sys.stdin.readline()
    filename_list=[]; criterion={}
    for filename in combined_results:
        filename_list.append(filename)
        criterion[string.join(string.split(filename,'-')[:-2],'-')]=None
    filename_list.sort() ### rank the results alphebetically so local and GO results are sequential
    ontology_to_filename={}; pathway_to_filename={}
    
    ### determine which Elite pathways are in which files and store

    if len(criterion)>1:
        ontology_files={}; local_files={}
        for filename in filename_list:
            file = export.findFilename(filename)[:-4]
            for (pathway_type,mapp_name,line) in combined_results[filename]:
                specific_pathway_type,pathway_type = pathway_type
                t = string.split(line,'\t')
                if pathway_type == 'Ontology':
                    go_name = t[1]
                    try: ontology_to_filename[go_name,specific_pathway_type].append(file)
                    except Exception: ontology_to_filename[go_name,specific_pathway_type]=[file]
                    if specific_pathway_type not in ontology_files:
                        ontology_files[specific_pathway_type] = [file]
                    elif file not in ontology_files[specific_pathway_type]:
                        ontology_files[specific_pathway_type].append(file)
                elif pathway_type == 'GeneSet':
                    pathway_name = t[0]
                    try: pathway_to_filename[pathway_name,specific_pathway_type].append(file)
                    except Exception: pathway_to_filename[pathway_name,specific_pathway_type]=[file]
                    if specific_pathway_type not in local_files:
                        local_files[specific_pathway_type] = [file]
                    elif file not in local_files[specific_pathway_type]:
                        local_files[specific_pathway_type].append(file)
                        
        headers = ['number_changed','number_measured',' z_score','permuteP']
        header_ontology = ['Ontology-ID','Ontology-term','Ontology-name','Ontology-type','elite_term_in']
        header_local = ['GeneSet-ID','GeneSet-term','GeneSet-name','GeneSet-type','elite_term_in']
        for ontology_type in ontologies_analyzed:
            input_file_type = ontologies_analyzed[ontology_type]
            if input_file_type == 'Ontology':
                ontology_elite_db = importORASimple(ora_dir,ontology_to_filename,ontology_type)
                writeOverlapLine(raw,ontology_files,headers,header_ontology,ontology_elite_db,ontology_to_filename,ontology_type)
            else:
                local_elite_db = importORASimple(ora_dir,pathway_to_filename,ontology_type)
                writeOverlapLine(raw,local_files,headers,header_local,local_elite_db,pathway_to_filename,ontology_type)
            raw.write('\n')
        raw.close()
        
        ### Move to the root directory
        fn = filepath(output_results)
        fn2 = string.replace(fn,'CompleteResults/ORA_pruned','')
        try:
            export.customFileMove(fn,fn2)
            from visualization_scripts import clustering
            clustering.clusterPathwayZscores(fn2) ### outputs the overlapping results as a heatmap
        except Exception,e:
            #print e
            #print fn; print fn2; print "OverlapResults failed to be copied from CompleteResults (please see CompleteResults instead)"
            pass
    
def writeOverlapLine(raw,ontology_files,headers,header_ontology,ontology_elite_db,ontology_to_filename,ontology_type):
    file_headers=[]
    filenames = ontology_files[ontology_type]
    filenames.sort() ### Sort by filename
    for filename in filenames:
        file_headers += updateHeaderName(headers,filename)
    title = string.join(header_ontology+file_headers,'\t')+'\n'
    raw.write(title)
    for (pathway_name,key) in ontology_elite_db:
        elite_files = string.join(ontology_to_filename[pathway_name,ontology_type],'|')
        row = list(key)+[elite_files]
        scores_data = ontology_elite_db[(pathway_name,key)]
        scores_data.sort() ### Sort by filename to match above
        for (file,values) in scores_data:
            row += values
        raw.write(string.join(row,'\t')+'\n')
        
def updateHeaderName(headers,filename):
    headers_copy = copy.deepcopy(headers)
    i=0
    for header in headers_copy:
        headers_copy[i]=header+'.'+filename
        i+=1
    return headers_copy

def output_combined_results(combined_results):
    if len(go_elite_output_folder) == 0: output_folder = "output"
    else: output_folder = go_elite_output_folder
    output_results = output_folder+'/'+'pruned-results_' +filter_method + '_elite.txt'
    proceed = 'no'
    while proceed == 'no':
        try: raw = export.ExportFile(output_results); proceed = 'yes'
        except Exception:
            print_out = output_results, '\nis open. Please close and select "Continue".'
            try: UI.WarningWindow(print_out,' OK ')
            except Exception:
                print print_out; print 'Please correct (hit return to continue)'; inp = sys.stdin.readline()
    filename_list=[]
    for filename in combined_results: filename_list.append(filename)
    filename_list.sort() ### rank the results alphebetically so local and GO results are seqeuntial
    for filename in filename_list:
        for (pathway_type,mapp_name,line) in combined_results[filename]:
            specific_pathway_type,pathway_type = pathway_type
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if pathway_type == 'Ontology':
                goid = t[0]
                del t[0]; #del t[2:7]
                go_type = t[1]; del t[1]; specific_pathway_type = go_type
                t[0] = t[0]+'('+goid+')'
            t.reverse(); t.append(specific_pathway_type); t.append(mapp_name); t.reverse()
            vals = string.join(t,'\t'); vals = vals + '\n'
            raw.write(vals)
        #print filename, len(combined_results[filename])
    raw.close()

def identifyGeneFiles(import_dir, mappfinder_input):
    e = GrabFiles(); e.setdirectory(import_dir)
    if len(gene_identifier_file)>0: ###global variable set upon MAPPFinder result import
        gene_file_dir = e.searchdirectory(gene_identifier_file)
        return gene_file_dir
    else:
        mappfinder_input = string.replace(mappfinder_input,'-GO.txt','') ###the prefix in most cases will be the same for the MAPPFinder results and gene input filename
        mappfinder_input = string.replace(mappfinder_input,'-local.txt','')
        split_name = string.split(mappfinder_input,'.')
        
        gene_file_dir = e.searchdirectory(mappfinder_input)
        if len(gene_file_dir)>0: return gene_file_dir
        else:
            try:
                index = int(split_name[0])
                index_str = str(index)
                gene_file_dir = e.searchdirectory_start(index_str)
            except ValueError: gene_file_dir =''
            return gene_file_dir

def grabAllNestedGOIDs(collapsed_go_tree,all_paths):
    ###Updated all_paths contains all possible paths listed in GO
    nested_collapsed_go_tree={}
    for goid_tuple in collapsed_go_tree:
        for goid in goid_tuple:
            child_goids=[]
            path_id_list = full_path_db[goid]
            for path_id in path_id_list: ###examine all possible paths for that goid (and different possible children)
                #path_id = path_id_data.PathTuple()
                if path_id in all_paths:
                    child_path_id_list = all_paths[path_id]
                    for child_path_id in child_path_id_list:
                        child_goid = path_id_to_goid[child_path_id]
                        child_goids.append(child_goid) ###append all possible children terms to a new list
        nested_collapsed_go_tree[goid] = child_goids
    nested_collapsed_go_tree = eliminate_redundant_dict_values(nested_collapsed_go_tree)
    return nested_collapsed_go_tree

def checkGOEliteSpecies(species):
    ### For an external program, see if the species data is supported
    speciesData()
    if species in species_names: return 'yes'
    else: return 'no'

def speciesData():
    importSpeciesData()
    if len(species_names)==0:
        UI.remoteSpeciesInfo('no') ### will get the missing data from a backup file
        importSpeciesData()
        
def importSpeciesData():
    program_type,database_dir = unique.whatProgramIsThis()
    if program_type == 'GO-Elite': filename = 'Config/species.txt'
    else: filename = 'Config/goelite_species.txt'
    x=0
    fn=filepath(filename);global species_list; species_list=[]; global species_codes; species_codes={}
    global species_names; species_names={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        abrev=t[0]; species=t[1]
        if x==0: x=1
        else:
            species_list.append(species)
            species_codes[species] = abrev
            species_names[abrev] = species

def testFileLength(fn):
    x=0
    for line in open(fn,'rU').readlines(): x+=1
    return x

def sourceDataCheck():
    filename = 'Config/source_data.txt'; fn=filepath(filename)    
    file_length = testFileLength(fn)
    if file_length <2:
        fn2 = string.replace(fn,'.txt','_archive.txt')
        import shutil; shutil.copyfile(fn2,fn) ### Bad file was downloaded (with warning)
    
def getSourceData():
    sourceDataCheck(); sourceData()
    return system_codes,source_types,mod_types

def getSourceDataNames():
    sourceDataCheck(); sourceData()
    system_code_names={}
    for sc in system_codes:
        name = system_codes[sc]
        system_code_names[name]=sc
    return system_code_names,source_types,mod_types

def sourceData():
    filename = 'Config/source_data.txt'; x=0
    fn=filepath(filename)
    global source_types; source_types=[]
    global system_codes; system_codes={}
    global mod_types; mod_types=[]
    for line in open(fn,'rU').readlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t'); source=t[0]
        try: system_code=t[1]
        except IndexError: system_code = 'NuLL'
        if x==0: x=1
        else:
            if len(t)>2: ### Therefore, this ID system is a potential MOD
                if t[2] == 'MOD': mod_types.append(source)
            if source not in mod_types: source_types.append(source) 
            system_codes[system_code] = source ###Used when users include system code data in their input file

def buildFullPathForCollapsedTree(collapsed_tree):
    filtered_tree_index=[]
    for (type,path) in collapsed_tree: filtered_tree_index.append(path)
    all_nested_paths = buildFullTreePath(filtered_tree_index)
    return all_nested_paths

def getNonEliteTerms(collapsed_goid_list,full_goid_db):
    non_elite_db=[]
    for goid in full_goid_db:
        if goid not in collapsed_goid_list: non_elite_db.append(goid)
        #non_elite_db.append(goid)
    print "Out of",len(collapsed_goid_list),"GO-Elite terms for",len(full_goid_db),"terms examined,",len(non_elite_db),"terms are not-elite"
    return non_elite_db

def countGOFullGenes(go_full,go_to_mod_genes):
    ### Count unique genes associated with each filtered pathway
    unique_go_full_genes={}
    for goid in go_full:
        if goid in go_to_mod_genes:
            for gene in go_to_mod_genes[goid]: unique_go_full_genes[gene] = []
    return len(unique_go_full_genes)
                    
def reorganizeResourceList(pathway_list):
    
    ### Cluster by resource type to minimize import of annotation and hierarchies
    geneset_types = {}
    for pathway in pathway_list:
        geneset_type = string.split(pathway,'-')[-1][:-4]
        try: geneset_types[geneset_type].append(pathway)
        except Exception: geneset_types[geneset_type] = [pathway]
        
    ### Save as a list
    pathway_list=[]
    for geneset_type in geneset_types:
        pathway_list+=geneset_types[geneset_type]
    
    ### Make sure that WikiPathways and GO are analyzed last, so that gene results are also reported last to GO_Elite.py
    add_pathway=[]
    pathway_list_reorganized=[]
    for pathway in pathway_list:
        if '-local.txt' in pathway: add_pathway.append(pathway)
        elif '-GO.txt' in pathway: add_pathway.append(pathway)
        else: pathway_list_reorganized.append(pathway)
    pathway_list_reorganized+=add_pathway
    
    return pathway_list_reorganized

def getAvaialbleResources(species_code):
    program_type,database_dir = unique.whatProgramIsThis()
    import_dir1 = '/'+database_dir+species_code+'/gene-mapp'
    import_dir2 = '/'+database_dir+species_code+'/gene-go'

    default_resources=[]
    try:
        gene_mapp_list = read_directory(import_dir1)
        gene_mapp_list.sort()
        for file in gene_mapp_list:
            resource = string.split(file,'-')[-1][:-4]
            if resource != 'MAPP' and resource not in default_resources and '.txt' in file:
                default_resources.append(resource)
            if resource == 'MAPP' and 'Pathways' not in default_resources and '.txt' in file:
                default_resources.append('Pathways')
    except Exception: pass
    try:
        gene_go_list = read_directory(import_dir2)
        gene_go_list.sort()
        for file in gene_go_list:
            resource = string.split(file,'-')[-1][:-4]
            if resource not in default_resources and 'version' not in resource and '.txt' in file:
                default_resources.append(resource)
    except Exception: pass

    return default_resources

def multiMappfinder(species,species_code,source_data,mod,system_codes,permutations,resources_to_analyze,file_dirs,root):
    """ Run in multiprocessing mode to run all jobs in parallel """
    
    multiProcessing = True
    multiprocessing_pipe = True
    if multiprocessing_pipe:
        queue = mlp.Queue()
    
    if type(resources_to_analyze) is list: resources = resources_to_analyze
    elif resources_to_analyze=='all': resources = getAvaialbleResources(species_code)
    elif resources_to_analyze == 'both': resources = ['GeneOntology','Pathways']
    elif isinstance(resources_to_analyze, list): resources = resources_to_analyze
    else: resources = [resources_to_analyze]
        
    variable_list=[]; incl=True
    for resource in resources:
        criterion_input_folder, criterion_denom_folder, output_dir, custom_sets_folder = file_dirs
        input_dir_files = readDirText(criterion_input_folder)
        for input_file in input_dir_files:
            if incl: custom_sets_folder = custom_sets_folder; x=1 ### This file should only be processed once
            else: custom_sets_folder = ''
            new_file_dirs = input_file, criterion_denom_folder, output_dir, custom_sets_folder 
            variables = species,species_code,source_data,mod,system_codes,permutations,resource,new_file_dirs,''
            variable_list.append(variables)

    processes=mlp.cpu_count()
    ### This is approach less efficeint than pool but we print out the progress
    s = processes; b=0
    list_of_vl=[]
    while s<len(variable_list):
        list_of_vl.append(variable_list[b:s])
        b+=processes; s+=processes
    list_of_vl.append(variable_list[b:s])
    
    start_time = time.time()
    print 'Performing Over Representation Analyses (ORA) on input datasets... (be patient)'
    go_ora_genes = {}
    local_ora_genes = {}
    ontology_ora_genes = {}
    geneset_ora_genes = {}
    
    if multiProcessing:
        if multiprocessing_pipe:
            for variable_list in list_of_vl:
                procs=list()
                for varset in variable_list:
                    species,species_code,source_data,mod,system_codes,permutations,resource,new_file_dirs,null = varset
                    proc = mlp.Process(target=mappfinder.generateMAPPFinderScores,args=(species,species_code,source_data,mod,system_codes,permutations,resource,new_file_dirs,None,True,queue))
                    procs.append(proc)
                    proc.start()
                    
                for _ in procs:
                    vals= queue.get()
                    if len(vals)==1:
                        print_out = vals[0]
                        if root==None:
                            print '\nWarning!!!',print_out,'\n'; sys.exit()
                        else:
                            try: UI.WarningWindow(print_out,'Critical Error!'); sys.exit()
                            except Exception: print '\nWarning!!!',print_out,'\n'; sys.exit()
                    else:
                        go_to_mod_genes, mapp_to_mod_genes, timediff, mappfinder_input, resource = vals
                        if 'GeneOntology' in resource:
                            go_ora_genes[mappfinder_input] = go_to_mod_genes
                        elif 'Pathways' in resource or 'MAPP' in resource:
                            local_ora_genes[mappfinder_input] = mapp_to_mod_genes
                        elif len(go_to_mod_genes)>0:
                            ontology_ora_genes[mappfinder_input] = go_to_mod_genes
                        elif len(mapp_to_mod_genes)>0:
                            geneset_ora_genes[mappfinder_input] = mapp_to_mod_genes
                        print '\tORA of',mappfinder_input[:-4],resource,'complete in %s seconds' % timediff
                    
                for proc in procs:
                    proc.join()
            
            ### Select a GO and local set (or other) to report for the last mappfinder_input considered
            input_file = export.findFilename(input_file)
            try:
                if len(go_ora_genes)>0:
                    go_to_mod_genes = go_ora_genes[input_file]
                elif len(ontology_ora_genes)>0:
                    go_to_mod_genes = ontology_ora_genes[input_file]
                if len(go_ora_genes)>0:
                    mapp_to_mod_genes = local_ora_genes[input_file]
                elif len(ontology_ora_genes)>0:
                    mapp_to_mod_genes = geneset_ora_genes[input_file]
            except Exception: pass
                
        else:
            pool_size = mlp.cpu_count()
            pool = mlp.Pool(processes=pool_size)
            results = pool.map(runMappfinderProcesses,variable_list)
            ### go_to_mod_genes and mapp_to_mod_genes are only needed for the last criterion (local and GO) for screen reporting
            if len(results)==1:
                print results[0]; forceError
            go_to_mod_genes, mapp_to_mod_genes, timediff, mappfinder_input, resource = results[-1]
    else:
        for variables in variable_list:
            species,species_code,source_data,mod,system_codes,permutations,resource,file_dirs,r = variables
            go_to_mod_genes, mapp_to_mod_genes, timediff, mappfinder_input, resource = mappfinder.generateMAPPFinderScores(species,species_code,source_data,mod,system_codes,permutations,resource,file_dirs,r,Multi=mlp)
    
    try: pool.close(); pool.join(); pool = None
    except Exception: pass
    end_time = time.time()
    time_diff = mappfinder.formatTime(start_time,end_time)
    print 'ORA analyses finished in %s seconds' % time_diff,'\n'
    return go_to_mod_genes, mapp_to_mod_genes

def runMappfinderProcesses(variables):
    try:
        species,species_code,source_data,mod,system_codes,permutations,resources_to_analyze,file_dirs,root = variables
        go_to_mod_genes, mapp_to_mod_genes, time_diff, mappfinder_input, resource = mappfinder.generateMAPPFinderScores(species,species_code,source_data,mod,system_codes,permutations,resources_to_analyze,file_dirs,root,poolVar=True)
    except Exception:
        return [traceback.format_exc()]
    return go_to_mod_genes, mapp_to_mod_genes, input_file, time_diff, mappfinder_input, resource

def runMappfinderPipe(variables):
    try:
        species,species_code,source_data,mod,system_codes,permutations,resources_to_analyze,file_dirs,queue = variables
        mappfinder.generateMAPPFinderScores(species,species_code,source_data,mod,system_codes,permutations,resources_to_analyze,file_dirs,root,poolVar=True)
    except Exception:
        root.put([traceback.format_exc()])

def readDirText(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt": dir_list2.append(sub_dir+'/'+entry)
    return dir_list2

def runGOElite(mod):

  print 'Running GO-Elite version 1.2.6\n'
  #UI.WarningWindow('test','test')    
  source_data = mod; speciesData(); sourceDataCheck(); sourceData()
  try: species_code = species_codes[species]
  except KeyError: species_code=''

  global exclude_related; global mappfinder_input; global nested_collapsed_go_tree; global uid_to_go; global collapsed_go_tree
  global path_dictionary_include; global organized_tree; global parent_highest_score; global gene_identifier_file; global main_output_folder
  global ontology_type; global ontologies_analyzed; ontologies_analyzed = {}
  exclude_related = 'yes'
  
  global program_type; global log_file
  program_type,database_dir = unique.whatProgramIsThis()
  global combined_results; combined_results={}; global go_to_mod_genes; global mapp_to_mod_genes; global zscore_goid
  global combined_associations; combined_associations={}; global combined_gene_ranking; combined_gene_ranking={}

  ### Add empty values for when there is no data on these values
  summary_keys = ['filtered_go_term_count','elite_go_term_count','redundant_go_term_count','filtered_local_count']
  summary_keys +=['redundant_mapp_term_count','go_gene_full_count','go_gene_elite_count','mapp_gene_elite_count']
  for key in summary_keys: summary_data_db[key]='0'
  
  print_items=[]
  print_items.append("Primary GO-Elite Parameters")
  print_items.append('\t'+'commandLineMode'+': '+commandLineMode)
  print_items.append('\t'+'species'+': '+species_code)
  print_items.append('\t'+'results-directory'+': '+main_output_folder)
  print_items.append('\t'+'filter_method'+': '+filter_method)
  print_items.append('\t'+'z_threshold'+': '+str(z_threshold))
  print_items.append('\t'+'p_val_threshold'+': '+str(p_val_threshold))
  print_items.append('\t'+'change_threshold'+': '+str(change_threshold))
  print_items.append('\t'+'sort_only_by_zscore'+': '+sort_only_by_zscore)
  print_items.append('\t'+'analysis_method'+': '+str(analysis_method))
  print_items.append('\t'+'criterion_input_folder'+': '+criterion_input_folder)
  print_items.append('\t'+'custom_sets_folder'+': '+custom_sets_folder)
  print_items.append('\t'+'enrichment_type'+': '+enrichment_type)

  universalPrintFunction(print_items)

  if run_mappfinder=='yes':
      global path_id_to_goid;  global full_path_db; global all_nested_paths; global collapsed_tree; export_databases = 'no'
      file_dirs = criterion_input_folder, criterion_denom_folder, main_output_folder, custom_sets_folder
      denom_search_dir = criterion_denom_folder
      
      multiProcessing = False ### For GO-Elite as a stand-alone app, be more conservative and allow for more printouts (less files typically run - slightly slower)
      if program_type == 'AltAnalyze': multiProcessing = False ### Faster in speed tests thus far (likely memory issues for many input sets run simultaneously)
      
      resources = resources_to_analyze ### Have to rename due to global conflict reassignment 
      if type(resources) is list:
        if 'all' in resources: resources = 'all'
        elif len(resources)==1: resources = resources[0]
        else: 
            multiProcessing = True ### Allows for distinct resources to be run in parallel
    
      multiProcessing = False
      permute = permutations ### Have to rename due to global conflict reassignment 
      if denom_search_dir==None:
          permute = 'FisherExactTest'
          
      print_items = []
      print_items.append("ORA Parameters")
      print_items.append('\t'+'mod'+': '+mod)
      print_items.append('\t'+'permutations'+': '+str(permute))
      print_items.append('\t'+'resources_to_analyze'+': '+str(resources))
      universalPrintFunction(print_items)
      
      if denom_search_dir==None:
        print '\nUPDATED PARAMETERS - Forcing ORA algorithm to be FisherExactTest (required with no supplied denominator)\n'
        denom_search_dir='' ### Need this to be a string for downstream functions
      ### go_to_mod_genes and mapp_to_mod_genes are only needed for the last criterion (local and GO) for screen reporting
      if multiProcessing == True:
           go_to_mod_genes, mapp_to_mod_genes = multiMappfinder(species,species_code,source_data,mod,system_codes,permute,resources,file_dirs,root)
      else:
           try:go_to_mod_genes, mapp_to_mod_genes, timediff, mappfinder_input, resource = mappfinder.generateMAPPFinderScores(species,species_code,source_data,mod,system_codes,permute,resources,file_dirs,root,Multi=mlp)
           except Exception: print traceback.format_exc()

      reload(mappfinder) ### Clear memory of any retained objects

  else: denom_search_dir = ''
  global export_user_summary_values; global go_elite_output_folder
  get_non_elite_terms_only = 'no' ### Used for internal benchmarking and QC
  export_user_summary_values = 'yes'
  user_defined_source = source_data
  user_defined_mod = mod

  if len(main_output_folder) == 0: import_dir = '/input/MAPPFinder'; go_elite_output_folder = ''
  else:
      if run_mappfinder == 'no':
          import_dir = main_output_folder
          output_dir = main_output_folder
          if 'GO-Elite_results/CompleteResults/ORA_pruned' in output_dir: output_dir = string.replace(output_dir,'/GO-Elite_results/CompleteResults/ORA_pruned','')
          main_output_folders = string.split(output_dir,'/')
          go_elite_output_folder = string.join(main_output_folders[:-1],'/') + '/GO-Elite_results/CompleteResults/ORA_pruned' ###Get only the parent directory
      else:
          if 'GO-Elite_results/CompleteResults/ORA' not in main_output_folder:
              import_dir = main_output_folder + '/GO-Elite_results/CompleteResults/ORA'
              go_elite_output_folder = main_output_folder+ '/GO-Elite_results/CompleteResults/ORA_pruned'
          else:
              import_dir = main_output_folder
              root_import_dir,child_dirs = string.split(main_output_folder,'/GO-Elite_results/CompleteResults/ORA')
              go_elite_output_folder = root_import_dir+ '/GO-Elite_results/CompleteResults/ORA_pruned'
          
  if len(criterion_input_folder)==0: gene_input_dir = '/input/GenesToQuery/'+species_code
  else: gene_input_dir = criterion_input_folder
  m = GrabFiles(); m.setdirectory(import_dir)
  try: dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
  except Exception:
      try: dir_list = read_directory(filepath('Databases/'+species_code+'/nested'))
      except Exception: dir_list = []
      log_report = open(log_file,'a')
      d = "File in the folder nested are: "+str(dir_list); log_report.write(d+'\n')
      log_report.write('('+filepath('Databases/'+species_code+'/nested')+')\n')
      mappfinder_db_input_dir = '/'+species_code+'/nested/'
      #try: print os.name, platform.win32_ver()[0], platform.architecture(), platform.mac_ver(), platform.libc_ver(), platform.platform()
      #except Exception: print os.name
      #exportLog(log_report)
      ### This exception was added in version 1.2 and replaces the code in OBO_import.buildNestedOntologyTree which resets the OBO version to 0/0/00 and re-runs (see unlisted_variable = kill)
      print_out = "Unknown error encountered during data processing.\nPlease see logfile in:\n\n"+log_file+"\nand report to altanalyze@gmail.com."
      if root != None:
        program,program_dir = unique.whatProgramIsThis()
        if program!= 'AltAnalyze':
            try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy()
            except Exception: print print_out
      """
      if os.name == 'nt':
            try: os.startfile('"'+log_file+'"')
            except Exception:  os.system('open "'+log_file+'"')
      elif 'darwin' in sys.platform: os.system('open "'+log_file+'"')
      elif 'linux' in sys.platform: os.system('xdg-open "'+log_file+'/"')   
      sys.exit()
      """
      
  global dataset_name; go_ora_files=0; local_ora_files=0; oraDirTogeneDir={}
  dir_list = reorganizeResourceList(dir_list)
  prior_ontology_type = None

  for mappfinder_input in dir_list:    #loop through each file in the directory to output results
    try:
        dataset_name = string.join(string.split(mappfinder_input,'-')[:-1],'-')
        print 'Pruning',mappfinder_input
        source_data = user_defined_source ###Resets the variables back to the user defined when missing from the MAPPFinder files
        mod = user_defined_mod
        mappfinder_input_dir = m.searchdirectory(mappfinder_input)
        ontology_type,input_file_type = checkPathwayType(mappfinder_input_dir)
        ontologies_analyzed[ontology_type] = input_file_type ### Use when examining overlapping terms between criterion

        if input_file_type == 'Ontology' and prior_ontology_type != ontology_type: ### Ensures we don't reimport an already imported tree-structure
            import_path_index()  ### need_batch_efficient
        try:
            run_mod,run_source,zscore_changed_path_db,go_full,go_titles,zscore_goid,input_file_type,gene_identifier_file,species_code = importMAPPFinderResults(mappfinder_input_dir)
        except Exception:
            print_out = 'Impropper ORA file format! Make sure to\n select the correct pre-processed input files.'
            if root != None:
                try: print "Analysis Failed"; UI.WarningWindow(print_out,'Critical Error!'); root.destroy()
                except IOError: print "Analysis Failed\n"
            else:
                print "Analysis Failed\n"
            sys.exit()
        if len(run_mod)>1: mod = run_mod
        if len(run_source)>1: source_data = run_source
        if input_file_type == 'Ontology':
            go_ora_files+=1
            if run_mappfinder == 'yes' and 'GO.txt' in mappfinder_input_dir:
                ### Used for gene-level summary reporting
                go_gene_full_count = countGOFullGenes(go_full,go_to_mod_genes); summary_data_db['go_gene_full_count'] = go_gene_full_count
            organized_tree,all_paths = buildOrderedTree(zscore_changed_path_db) ### Improved GO-tree reconstruction method - implemented in version 1.22
            parent_highest_score = link_score_to_all_paths(all_paths,zscore_changed_path_db)
            child_highest_score = calculate_score_for_children(organized_tree,zscore_changed_path_db)
            collapsed_tree = collapse_tree(parent_highest_score,child_highest_score,organized_tree)
            collapsed_go_tree = link_goid(collapsed_tree,zscore_changed_path_db,all_paths)
            collapsed_go_list = exportGOResults(go_full,go_titles,collapsed_go_tree,zscore_goid,{},{},{},{})
            summary_data_db['elite_go_term_count'] =  len(collapsed_go_list) ### All Elite GO-terms
            summary_data_db['filtered_go_term_count'] = len(go_full)
            ###For testing purposes, we can export all non-Elite terms
            if get_non_elite_terms_only == 'yes': collapsed_go_list = getNonEliteTerms(collapsed_go_list,go_full)
    
            ###Identify gene lists for the corresponding GO-elite results and generate nested associations
            try: gene_file_dir = identifyGeneFiles(gene_input_dir, mappfinder_input)
            except Exception: gene_file_dir=''
            if len(gene_file_dir) > 0:
                #print "Begining to re-annotate pruned results...",
                nested_paths_stored = species_code
                if len(denom_search_dir)==0: search_dir = export.findParentDir(gene_input_dir)
                else: search_dir = denom_search_dir ### Occurs when analyzing an input list and existing pruned results
                if prior_ontology_type != ontology_type:
                    try:
                        uid_to_go, uid_system, gene_annotations = gene_associations.grabNestedGeneToOntologyAssociations(species_code,
                                                                mod,source_data,system_codes,search_dir,ontology_type)
                    except Exception,e:
                        if root != None:
                            print_out = e
                            try: UI.WarningWindow(print_out,' Continue '); sys.exit()
                            except Exception: print print_out; sys.exit()
                #print 'Annotations imported'
                try:
                    vals = gene_associations.matchInputIDsToGOEliteTerms(gene_file_dir,
                                            go_elite_output_folder,system_codes,mappfinder_input_dir,
                                            collapsed_go_list,uid_to_go,gene_annotations,full_go_name_db,
                                            uid_system,combined_associations,combined_gene_ranking)
                    combined_associations,combined_gene_ranking,go_gene_annotation_db,go_values_db,value_headers,goids_with_redundant_genes,go_gene_elite_count = vals
                    summary_data_db['go_gene_elite_count'] = go_gene_elite_count
                    unique_genes={}
                    for goid in go_gene_annotation_db:
                        for s in go_gene_annotation_db[goid]: unique_genes[s.GeneID()]=[]
                    #print len(unique_genes), "unique genes associated with GO-Elite terms"          
                    ###Re-output results, now with gene annotation data
                    collapsed_go_list = exportGOResults(go_full,go_titles,collapsed_go_list,zscore_goid,go_gene_annotation_db,go_values_db,value_headers,goids_with_redundant_genes)    
                    exportFilteredSIF(mod,species_code,collapsed_go_list,mappfinder_input_dir,None)
                except Exception:
                    continue
        else:
            local_ora_files+=1
            if run_mappfinder == 'yes': mapp_gene_full_count = countGOFullGenes(zscore_goid,mapp_to_mod_genes); summary_data_db['mapp_gene_full_count'] = mapp_gene_full_count
            filtered_mapp_list = zscore_changed_path_db
            exportLocalResults(go_full,go_titles,{},{},{},{})
            summary_data_db['filtered_local_count'] =  len(go_full) ### All Elite GO-terms
            ###Identify gene lists for the corresponding GO-elite results and generate nested associations
            try: gene_file_dir = identifyGeneFiles(gene_input_dir, mappfinder_input)
            except Exception: gene_file_dir = ''
            oraDirTogeneDir[mappfinder_input] = gene_file_dir ### Store the corresponding gene file for each ORA file
            if len(gene_file_dir) > 0:
                nested_paths_stored = species_code
                if prior_ontology_type != ontology_type:
                    uid_to_mapp, uid_system, gene_annotations = gene_associations.grabNestedGeneToPathwayAssociations(species_code,
                                                        mod,source_data,system_codes,custom_sets_folder,denom_search_dir,ontology_type)
                    #print 'Annotations imported'
                if len(uid_to_mapp)>0: ### alternative occurs if analyzing a custom_gene_set result without referencing it again (only should occur during testing)
                    try:
                        vals = gene_associations.matchInputIDsToMAPPEliteTerms(gene_file_dir,
                                                go_elite_output_folder,system_codes,mappfinder_input_dir,
                                                uid_to_mapp,filtered_mapp_list,gene_annotations,uid_system,
                                                combined_associations,combined_gene_ranking)
                        combined_associations,combined_gene_ranking,mapp_gene_annotation_db,mapp_values_db,mapp_value_headers,mapps_with_redundant_genes,mapp_gene_elite_count = vals
                        summary_data_db['mapp_gene_elite_count'] = mapp_gene_elite_count
                    except Exception:
                        continue
    
                exportLocalResults(go_full,go_titles,mapp_gene_annotation_db,mapp_values_db,mapp_value_headers,mapps_with_redundant_genes)
                exportFilteredSIF(mod,species_code,mapp_gene_annotation_db,mappfinder_input_dir,oraDirTogeneDir)
            if program_type != 'GO-Elite' and mappfinder_input[:3] == 'AS.':
                local_filename = go_elite_output_folder+'/'+mappfinder_input[0:-4]+ '_'+filter_method+'_elite.txt'
                print 'Copying GO-Elite results to DomainGraph folder...'
                fn = filepath(local_filename)
                fn2 = string.replace(fn,'GO-Elite_results','DomainGraph')
                fn2 = string.replace(fn2,'GO-Elite','AltResults')
                fn2 = string.replace(fn2,'AS.','')
                fn2 = string.split(fn2,'-local'); fn2=fn2[0]+'-pathways-DomainGraph.txt'
                #shutil.copyfile(fn,fn2)
        prior_ontology_type = ontology_type
    except Exception,e:
        print traceback.format_exc()
        print 'Error encountered in GO-Elite results pruning for this gene-set type'
  print 'gene associations assigned'

  if '/app' in filepath(import_dir): webservice = 'yes'
  else: webservice = 'no'

  if len(combined_results)>0:
      ora_files = [local_ora_files,go_ora_files]
      output_combined_results(combined_results)

      if max(ora_files)>1 and 'Heatmap' not in analysis_method: ### When GO-Elite is run from the heatmap function in clustering, calling clustering below will cause extra matplotlib show occurances
        outputOverlappingResults(combined_results,import_dir) ### Compare all input file results
        
      gene_associations.exportCombinedAssociations(combined_associations,go_elite_output_folder,'gene-associations')
      gene_associations.exportCombinedAssociations(combined_gene_ranking,go_elite_output_folder,'gene-ranking')
      
      ### Copy over ss from CompleteResults
      combined_ora_summary = go_elite_output_folder+'/'+'pruned-results_' +filter_method + '_elite.txt'
      fn = filepath(combined_ora_summary)
      fn2 = string.replace(fn,'CompleteResults/ORA_pruned','')
      try: export.customFileMove(fn,fn2)
      except Exception: print fn; print fn2; print "SummaryResults failed to be copied from CompleteResults (please see CompleteResults instead)"
      combined_ora_gene = go_elite_output_folder+'/gene_associations/pruned-gene-associations.txt'
      fn = filepath(combined_ora_gene)
      fn2 = string.replace(fn,'CompleteResults/ORA_pruned/gene_associations','')
      try: export.customFileMove(fn,fn2)
      except IOError: print fn; print fn2; print "SummaryResults failed to be copied from CompleteResults (please see CompleteResults instead)"
      
      visualizePathways(species_code,oraDirTogeneDir,combined_results)
      
      if webservice == 'no':
          try: moveMAPPFinderFiles(import_dir)
          except Exception,e:  print "Could not move ORA results... this will not impact this analysis:",e
      end_time = time.time(); time_diff = int(end_time-start_time)
      print_out = 'Analysis completed. GO-Elite results\nexported to the specified GO-Elite\noutput and ORA output directories.'
      try:
          if (program_type == 'GO-Elite' or analysis_method == 'UI') and use_Tkinter == 'yes' and root != None:
              print 'Analysis completed in %d seconds... Exiting GO-Elite' % time_diff
              UI.InfoWindow(print_out,'Analysis Completed!')
              tl = Toplevel(); SummaryResultsWindow(tl,main_output_folder)
          else:
              ### Record results
              result_list = recordGOEliteStats(); log_report.append('\n')
      except Exception:
          try: print 'Analysis completed in %d seconds... Exiting GO-Elite' % time_diff
          except Exception: null=[]
  else:
      end_time = time.time(); time_diff = int(end_time-start_time)
      print_out = 'No input files to summarize!'
      if program_type == 'GO-Elite' and use_Tkinter == 'yes' and root != None:
          print 'Analysis completed in %d seconds... Exiting GO-Elite' % time_diff
          UI.WarningWindow(print_out,'Analysis Completed!')

  #exportLog(log_report)
  run_parameter = 'skip'
  end_time = time.time(); time_diff = int(end_time-start_time)
  if program_type == 'GO-Elite' and use_Tkinter == 'yes': importGOEliteParameters(run_parameter)

def exportFilteredSIF(mod,species_code,collapsed_term_list,mappfinder_input_dir,oraDirTogeneDir):
    gene_association_file = string.replace(mappfinder_input_dir,'.txt','-associations.tab')
    sif_output = string.replace(mappfinder_input_dir,'ORA','ORA_pruned')
    sif_output = string.replace(sif_output,'.txt','.sif')
    sif = export.ExportFile(sif_output)
    
    ### Import the gene associations for all GO/pathways analyzed and output Elite filtered results as a SIF
    fn=filepath(gene_association_file); x = 0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x=1
        else:
            ID,symbol,term = string.split(data,'\t')
            if term in collapsed_term_list:
                if term in full_go_name_db: term = full_go_name_db[term]
                sif.write(string.join([term,'pr',ID+':'+symbol],'\t')+'\n')
    sif.close()
    
    try:
        from visualization_scripts import clustering
        try:
            criterion_name = export.findFilename(mappfinder_input_dir)
            ora_input_dir = oraDirTogeneDir[criterion_name] ### This is ONLY needed for transcription factor graph visualization
        except Exception: ora_input_dir = None
        clustering.buildGraphFromSIF(mod,species_code,sif_output,ora_input_dir)
    except Exception:
        #print traceback.format_exc()
        pass #Export from PyGraphViz not supported
        
class SummaryResultsWindow:
    def __init__(self,tl,output_dir):
        def showLink(event):
            idx= int(event.widget.tag_names(CURRENT)[1])
            webbrowser.open(LINKS[idx])
        url = 'http://www.genmapp.org/go_elite/help_main.htm'
        LINKS=(url,'')
        self.LINKS = LINKS
        tl.title('GO-Elite version 1.2.6'); self.tl = tl
        
        #"""
        filename = 'Config/icon.gif'
        fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(tl); can.pack(side='top'); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        #"""
        use_scroll = 'no'
        
        label_text_str = 'GO-Elite Result Summary'; height = 300; width = 510
        if os.name != 'nt': width = 550
        elif platform.win32_ver()[0] == 'Vista': width = 550
        self.sf = PmwFreeze.ScrolledFrame(tl,
            labelpos = 'n', label_text = label_text_str,
            usehullsize = 1, hull_width = width, hull_height = height)
        self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
        self.frame = self.sf.interior()
        
        txt=Text(self.frame,bg='gray')                    
        txt.pack(expand=True, fill="both")
        txt.insert(END, 'Primary Analysis Finished....\n')
        txt.insert(END, '\nResults saved to:\n'+output_dir+'\n')
        txt.insert(END, '\nFor more information see the ')
        txt.insert(END, "GO-Elite Help", ('link', str(0)))
        txt.insert(END, '\n\n')

        try: result_list = recordGOEliteStats()
        except Exception: result_list=[]
        for d in result_list: txt.insert(END, d+'\n');
            
        txt.tag_config('link', foreground="blue", underline = 1)
        txt.tag_bind('link', '<Button-1>', showLink)

        open_results_folder = Button(self.tl, text = 'Results Folder', command = self.openDirectory)
        open_results_folder.pack(side = 'left', padx = 5, pady = 5);

        self.dg_url = 'http://www.wikipathways.org'
        text_button = Button(self.tl, text='WikiPathways website', command=self.WPlinkout)
        text_button.pack(side = 'right', padx = 5, pady = 5)
        self.output_dir = output_dir
        self.whatNext_url = 'http://code.google.com/p/go-elite/wiki/Tutorial_GUI_version#Downstream_Analyses'
        whatNext_pdf = 'Documentation/what_next_goelite.pdf'; whatNext_pdf = filepath(whatNext_pdf); self.whatNext_pdf = whatNext_pdf
        
        what_next = Button(self.tl, text='What Next?', command=self.whatNextlinkout)
        what_next.pack(side = 'right', padx = 5, pady = 5)
        quit_buttonTL = Button(self.tl,text='Close View', command=self.close)
        quit_buttonTL.pack(side = 'right', padx = 5, pady = 5)
        
        continue_to_next_win = Button(text = 'Continue', command = root.destroy)
        continue_to_next_win.pack(side = 'right', padx = 10, pady = 10)
        quit_button = Button(root,text='Quit', command=self.quit)
        quit_button.pack(side = 'right', padx = 5, pady = 5)

        button_text = 'Help'; url = 'http://www.genmapp.org/go_elite/help_main.htm'; self.help_url = url      
        pdf_help_file = 'Documentation/GO-Elite_Manual.pdf'; pdf_help_file = filepath(pdf_help_file); self.pdf_help_file = pdf_help_file
        help_button = Button(root, text=button_text, command=self.Helplinkout)
        help_button.pack(side = 'left', padx = 5, pady = 5); root.mainloop()
        
        tl.mainloop() ###Needed to show graphic
    def openDirectory(self):
        if os.name == 'nt':
            try: os.startfile('"'+self.output_dir+'"')
            except Exception:  os.system('open "'+self.output_dir+'"')
        elif 'darwin' in sys.platform: os.system('open "'+self.output_dir+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+self.output_dir+'"')  

    def Helplinkout(self): self.GetHelpTopLevel(self.help_url,self.pdf_help_file)
    def whatNextlinkout(self):
        #self.GetHelpTopLevel(self.whatNext_url,self.whatNext_pdf) ### Provides the option to open a URL or PDF
        webbrowser.open(self.whatNext_url) ### Just defaults to the URL
    def GetHelpTopLevel(self,url,pdf_file):
        self.pdf_file = pdf_file; self.url = url
        message = ''; self.message = message; self.online_help = 'Online Documentation'; self.pdf_help = 'Local PDF File'
        tl = Toplevel(); self._tl = tl; nulls = '\t\t\t\t'; tl.title('Please select one of the options')
        self.sf = PmwFreeze.ScrolledFrame(self._tl,
                labelpos = 'n', label_text = '', usehullsize = 1, hull_width = 220, hull_height = 150)
        self.sf.pack(padx = 10, pady = 10, fill = 'both', expand = 1)
        self.frame = self.sf.interior()
        group = PmwFreeze.Group(self.sf.interior(),tag_text = 'Options')
        group.pack(fill = 'both', expand = 1, padx = 20, pady = 10)
        l1 = Label(group.interior(), text=nulls);  l1.pack(side = 'bottom')
        text_button2 = Button(group.interior(), text=self.online_help, command=self.openOnlineHelp); text_button2.pack(side = 'top', padx = 5, pady = 5) 
        try: text_button = Button(group.interior(), text=self.pdf_help, command=self.openPDFHelp); text_button.pack(side = 'top', padx = 5, pady = 5)
        except Exception: text_button = Button(group.interior(), text=self.pdf_help, command=self.openPDFHelp); text_button.pack(side = 'top', padx = 5, pady = 5)
        tl.mainloop()
    def openPDFHelp(self):
        if os.name == 'nt':
            try: os.startfile('"'+self.pdf_file+'"')
            except Exception:  os.system('open "'+self.pdf_file+'"')
        elif 'darwin' in sys.platform: os.system('open "'+self.pdf_file+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+self.pdf_file+'"')   
        self._tl.destroy()
    def openOnlineHelp(self):
        try: webbrowser.open(self.url)
        except Exception: null=[]
        self._tl.destroy()
    def WPlinkout(self): webbrowser.open(self.dg_url)
    def quit(self):
        #exportLog(log_report)
        try: self.tl.quit(); self.tl.destroy()
        except Exception: null=[]
        root.quit(); root.destroy(); sys.exit()
    def close(self):
        self.tl.quit()
        self.tl.destroy()

def recordGOEliteStats():
    if commandLineMode == 'no': log_report = open(log_file,'a') 
    result_list=[]
    for key in summary_data_db: summary_data_db[key] = str(summary_data_db[key])
    d = 'Dataset name: '+ dataset_name; result_list.append(d+'\n')
    d = summary_data_db['filtered_go_term_count']+':\tReported Filtered GO-terms'; result_list.append(d)
    d = summary_data_db['elite_go_term_count']+':\tReported Elite GO-terms'; result_list.append(d)
    d = summary_data_db['redundant_go_term_count']+':\tElite GO-terms redundant with other Elite GO-terms'; result_list.append(d)
    d = summary_data_db['filtered_local_count']+':\tReported Filtered Pathways'; result_list.append(d)
    d = summary_data_db['redundant_mapp_term_count']+':\tFiltered Pathways redundant with other Filtered Pathways'; result_list.append(d)

    if run_mappfinder == 'yes':
        try: 
            if mod == 'HMDB': name = 'metabolites'
            else: name = 'genes'
        except Exception: name = 'genes'
    else: name = 'genes'
    if run_mappfinder == 'yes':
        d = summary_data_db['go_gene_full_count']+':\tNumber of '+name+' associated with Filtered GO-terms'; result_list.append(d)
    d = summary_data_db['go_gene_elite_count']+':\tNumber of '+name+' associated with Elite GO-terms'; result_list.append(d)
    d = summary_data_db['mapp_gene_elite_count']+':\tNumber of '+name+' associated with Filtered Pathways'; result_list.append(d)
        
    if commandLineMode == 'no':
        log_report.write('\n')
    for d in result_list:
        if commandLineMode == 'no': log_report.write(d+'\n')
        else: print d
    if commandLineMode == 'no': log_report.close()
    return result_list
    
class StatusWindow:
    def __init__(self,root,mod):
            root.title('GO-Elite 1.2.6 - Status Window')
            self._parent = root
            
            statusVar = StringVar() ### Class method for Tkinter. Description: "Value holder for strings variables."

            if os.name == 'nt': height = 450; width = 500
            else: height = 490; width = 640
            self.sf = PmwFreeze.ScrolledFrame(self._parent,
                    labelpos = 'n', label_text = 'Results Status Window',
                    usehullsize = 1, hull_width = width, hull_height = height)
            self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
            self.frame = self.sf.interior()

            group = PmwFreeze.Group(self.sf.interior(),tag_text = 'Output')
            group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)
                
            Label(group.interior(),width=180,height=452,justify=LEFT, bg='black', fg = 'white',anchor=NW,padx = 5,pady = 5, textvariable=statusVar).pack(fill=X,expand=Y)

            status = StringVarFile(statusVar,root) ### Likely captures the stdout
            sys.stdout = status; root.after(100, runGOElite(mod))

            try: self._parent.mainloop()
            except Exception: null=[]
            try: self._parent.destroy()
            except Exception: null=[]

    def deleteWindow(self):
        #tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)
        try: self._parent.destroy(); sys.exit() ### just quit instead
        except Exception: sys.exit()
    def quit(self):
        try: self._parent.quit(); self._parent.destroy(); sys.exit()
        except Exception: sys.exit()
    
class StringVarFile:
    def __init__(self,stringVar,window):
        self.__newline = 0; self.__stringvar = stringVar; self.__window = window
    def write(self,s):
        log_report = open(log_file,'a')
        log_report.write(s); log_report.close() ### Variable to record each print statement
        new = self.__stringvar.get()
        for c in s:
            #if c == '\n': self.__newline = 1
            if c == '\k': self.__newline = 1### This should not be found and thus results in a continous feed rather than replacing a single line
            else:
                if self.__newline: new = ""; self.__newline = 0
                new = new+c
        self.set(new)
    def set(self,s):
        try: self.__stringvar.set(s); self.__window.update()
        except Exception: sys.exit() ### When the application is closed or exited by force
    def get(self): return self.__stringvar.get()
    def flush(self): pass
                
def timestamp():
    import datetime
    today = str(datetime.date.today()); today = string.split(today,'-'); today = today[0]+''+today[1]+''+today[2]
    time_stamp = string.replace(time.ctime(),':','')
    time_stamp = string.replace(time_stamp,'  ',' ')
    time_stamp = string.split(time_stamp,' ') ###Use a time-stamp as the output dir (minus the day)
    time_stamp = today+'-'+time_stamp[3]
    return time_stamp

def TimeStamp():
    time_stamp = time.localtime()
    year = str(time_stamp[0]); month = str(time_stamp[1]); day = str(time_stamp[2])
    if len(month)<2: month = '0'+month
    if len(day)<2: day = '0'+day
    return year+month+day

def universalPrintFunction(print_items):
    if commandLineMode == 'no':
        log_report = open(log_file,'a')
    for item in print_items:
        if commandLineMode == 'no': ### Command-line has it's own log file write method (Logger)
            log_report.write(item+'\n')
        else: print item
    if commandLineMode == 'no':
        log_report.close()

class Logger(object):
    def __init__(self,null): 
        self.terminal = sys.stdout
        self.log = open(log_file, "w")

    def write(self, message):
        self.log = open(log_file, "a")
        self.terminal.write(message)
        self.log.write(message)
        self.log.close()
        
    def flush(self): pass
    
def importGOEliteParameters(run_parameter):
    global max_member_count; global filter_method; global z_threshold; global p_val_threshold; global change_threshold
    global sort_only_by_zscore; global permutations; global species; global analysis_method; global resources_to_analyze
    global run_mappfinder; global criterion_input_folder; global criterion_denom_folder; global main_output_folder
    global log_file; global summary_data_db; summary_data_db = {}; global custom_sets_folder; global mod; global returnPathways
    global imageType; imageType = True; global commandLineMode; commandLineMode = 'no'; global enrichment_type

    parameters = UI.getUserParameters(run_parameter)
    species, run_mappfinder, mod, permutations, filter_method, z_threshold, p_val_threshold, change_threshold, resources_to_analyze, max_member_count, returnPathways, file_dirs, enrichment_type = parameters
    change_threshold = int(change_threshold); p_val_threshold = float(p_val_threshold); z_threshold = float(z_threshold)
    criterion_input_folder, criterion_denom_folder, main_output_folder, custom_sets_folder = file_dirs
    time_stamp = timestamp() 
    log_file = filepath(main_output_folder+'/GO-Elite_report-'+time_stamp+'.log')
    log_report = open(log_file,'w'); log_report.close()

    sort_only_by_zscore = 'yes'; analysis_method = 'UI'
    
    global root; root=''
    if use_Tkinter == 'yes' and debug_mode == 'no':
        root = Tk()
        StatusWindow(root,mod)
        try: root.destroy()
        except Exception: null=[]
    else:
        sys.stdout = Logger('')
        runGOElite(mod)

def remoteAnalysis(variables,run_type,Multi=None):
    global max_member_count; global filter_method; global z_threshold; global p_val_threshold; global change_threshold
    global sort_only_by_zscore; global permutations; global species; global root; global custom_sets_folder
    global run_mappfinder; global criterion_input_folder; global criterion_denom_folder; global main_output_folder
    global summary_data_db; summary_data_db = {}; global analysis_method; global log_file; global mod
    global resources_to_analyze; global returnPathways; global imageType; imageType = 'both'
    global commandLineMode; commandLineMode = 'no'; global enrichment_type; enrichment_type = 'ORA'
    global mlp; mlp = Multi
    
    try: species_code,mod,permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,returnPathways,file_dirs,enrichment_type,parent = variables
    except Exception: species_code,mod,permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,returnPathways,file_dirs,parent = variables
    #print variables
    try: permutations = int(permutations)
    except Exception: permutations = permutations ### For Fisher Exact
    change_threshold = int(change_threshold); p_val_threshold = float(p_val_threshold); z_threshold = float(z_threshold)
    if resources_to_analyze == 'WikiPathways' or resources_to_analyze == 'local':
        resources_to_analyze = 'Pathways'
    speciesData()
    species = species_names[species_code]
    max_member_count = 10000; sort_only_by_zscore = 'yes'; run_mappfinder = 'yes'
    criterion_input_folder, criterion_denom_folder, main_output_folder = file_dirs; custom_sets_folder = '' ### Variable not currently used for AltAnalyze
    
    time_stamp = timestamp() 
    log_file = filepath(main_output_folder+'/GO-Elite_report-'+time_stamp+'.log')
    
    if 'non-UI' in run_type:
        analysis_method = run_type
        commandLineMode = 'yes'
        sys.stdout = Logger('')
        root = parent; runGOElite(mod)
    elif run_type == 'UI':
        log_report = open(log_file,'a'); log_report.close()
        analysis_method = run_type
        root = Tk()
        StatusWindow(root,mod)
        try: root.destroy()
        except Exception: null=[]

def visualizePathways(species_code,oraDirTogeneDir,combined_results):
    """ Sends all over-represented pathways to the WikiPathways API for visualization """
    wp_start_time = time.time()
    try:
        failed=[]
        if returnPathways != None and returnPathways != 'None' and returnPathways != 'no':
            ### If only the top X pathways should be returned, get this number
            returnNumber = None
            if 'top' in returnPathways:
                returnNumber = int(string.replace(returnPathways,'top',''))
                
            count=0
            filename_list=[]
            for filename in combined_results: filename_list.append(filename)
            filename_list.sort() ### rank the results alphebetically so local and GO results are sequential
            wp_to_visualize = []
            for filename in filename_list:
                if 'local' in filename:
                    criterion_name = string.split(export.findFilename(filename)[:-4],'-local')[0]+'-local.txt'
                for (pathway_type,mapp_name,line) in combined_results[filename]:
                    if ':WP' in line:
                        gene_file_dir = oraDirTogeneDir[criterion_name]
                        wpid = 'WP'+string.split(line,':WP')[-1]
                        wpid = string.split(wpid,'\t')[0]
                        wp_to_visualize.append([gene_file_dir,wpid])
            if len(wp_to_visualize)>0:
                print 'Exporting pathway images for %s Wikipathways (expect %s minute runtime)' % (str(len(wp_to_visualize)),str(len(wp_to_visualize)/15)+'-'+str(len(wp_to_visualize)/2))
                
                if imageType == 'png' or imageType == 'pdf':
                    print '...restricting to',imageType

                try: wp_to_visualize = getProducedPathways(wp_to_visualize); #print len(wp_to_visualize) ### Retain pathways from older runs
                except Exception: pass
                try: poolWPVisualization(wp_to_visualize,species_code,mod,imageType)
                except Exception,e:
                    #print e
                    try: wp_to_visualize = getProducedPathways(wp_to_visualize); #print len(wp_to_visualize)
                    except Exception: pass
                    print 'Error encountered with Multiple Processor Mode or server timeout... trying single processor mode'
                    for (gene_file_dir,wpid) in wp_to_visualize:
                        ### graphic_link is a dictionary of PNG locations
                        try:
                            graphic_link = WikiPathways_webservice.visualizePathwayAssociations(gene_file_dir,species_code,mod,wpid,imageExport=imageType)
                            count+=1
                            print '.',
                            if returnNumber != None:
                                if returnNumber == count:
                                    break
                        except Exception:
                            error_report = traceback.format_exc()
                            failed.append(wpid)
                    print 'exported to the folder "WikiPathways"'
                if len(failed)>0:
                    print len(failed),'Wikipathways failed to be exported (e.g.,',failed[0],')'
                    print error_report
    except Exception:
        pass
    wp_end_time = time.time(); time_diff = int(wp_end_time-wp_start_time)
    print "Wikipathways output in %d seconds" % time_diff
    
def makeVariableList(wp_to_visualize,species_code,mod,imageType):
    variable_list=[]
    for (gene_file_dir,wpid) in wp_to_visualize:
        variables = gene_file_dir,species_code,mod,wpid,imageType
        variable_list.append(variables)
    return variable_list

def poolWPVisualization(wp_to_visualize,species_code,mod,imageType):

    wp_n = len(wp_to_visualize)
    variable_list = makeVariableList(wp_to_visualize,species_code,mod,imageType)
        
    pool_size = int(mlp.cpu_count() * 1)
    pool = mlp.Pool(processes=pool_size)
    try: results = pool.map(runWPProcesses,variable_list) ### try once
    except Exception:
        if 'forceTimeError' in traceback.format_exc():
            try: pool.close(); pool.join(); pool = None
            except Exception: pass
            forceTimeError
        try:
            wp_to_visualize = getProducedPathways(wp_to_visualize)
            variable_list = makeVariableList(wp_to_visualize,species_code,mod,imageType)
        except Exception: pass #Occurs when the WikiPathways directory is not created yet
        #if wp_n == len(wp_to_visualize): force_pool_error
        #else:
        try: results = pool.map(runWPProcesses,variable_list) ### try again with only failed pathways
        except Exception:
            try:
                wp_to_visualize = getProducedPathways(wp_to_visualize)
                variable_list = makeVariableList(wp_to_visualize,species_code,mod,imageType)
            except Exception: pass #Occurs when the WikiPathways directory is not created yet
            results = pool.map(runWPProcesses,variable_list)
            
    try: pool.close(); pool.join(); pool = None
    except Exception: pass
    return results

def getProducedPathways(wp_to_visualize):
    """ If a pooled run fails, try again for those pathways not output """
    
    for (filename,wpid) in wp_to_visualize: break
    root_dir = export.findParentDir(filename)
    if 'GO-Elite/input' in root_dir:
        root_dir = string.replace(root_dir,'GO-Elite/input','WikiPathways')
    else:
        root_dir+='WikiPathways/'
    
    possible_pathways={} ### Create a list of all possible created pathways
    for (gene_file_dir,wpid) in wp_to_visualize:
        dataset_name = export.findFilename(gene_file_dir)[:-4]
        try: possible_pathways[wpid].append((dataset_name,gene_file_dir,wpid))
        except Exception: possible_pathways[wpid] = [(dataset_name,gene_file_dir,wpid)]

    pathway_list = unique.read_directory(root_dir)
    pathway_db={}
    for pathway in pathway_list:
        wpid = string.split(pathway,'_')[0]
        try: pathway_db[wpid].append(pathway)
        except Exception: pathway_db[wpid] = [pathway] ### pathway name contains the wpid, pathway name, dataset (pathway name from webservice)
    
    wp_to_visualize=[]
    for wpid in possible_pathways:
        if wpid in pathway_db:
            pathway_images = pathway_db[wpid]
            for (dataset_name,gene_file_dir,wpid) in possible_pathways[wpid]:
                proceed = False
                for pathway in pathway_images:
                    if dataset_name in pathway: proceed = True#; print dataset_name, wpid
                if proceed==False:
                    wp_to_visualize.append((gene_file_dir,wpid))    
        else:
            for (dataset_name,gene_file_dir,wpid) in possible_pathways[wpid]:
                wp_to_visualize.append((gene_file_dir,wpid))    
    return wp_to_visualize
    
def runWPProcesses(variables):
    gene_file_dir,species_code,mod,wpid,imageType = variables
    st = time.time()
    try: graphic_link = WikiPathways_webservice.visualizePathwayAssociations(gene_file_dir,species_code,mod,wpid,imageExport=imageType)
    except Exception:
        try: graphic_link = WikiPathways_webservice.visualizePathwayAssociations(gene_file_dir,species_code,mod,wpid,imageExport=imageType)
        except Exception: pass
    if (time.time() - st)> 100: forceTimeError ### Too slow this way
    return graphic_link

def displayHelp():
    fn=filepath('Documentation/commandline.txt')
    print '\n################################################\nGO-Elite Command-Line Help'
    for line in open(fn,'rU').readlines():
        print cleanUpLine(line)
    print '\n################################################'
    sys.exit()
    
class SpeciesData:
    def __init__(self, abrev, species, systems, taxid):
        self._abrev = abrev; self._species = species; self._systems = systems; self._taxid = taxid
    def SpeciesCode(self): return self._abrev
    def SpeciesName(self): return self._species
    def Systems(self): return self._systems
    def TaxID(self): return self._taxid
    def __repr__(self): return self.SpeciesCode()+'|'+self.SpeciesName()
        
def importSpeciesInfo():
    filename = 'Config/species_all_archive.txt'
    fn=filepath(filename); global species_list; species_list=[]; global species_codes; species_codes={}; x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        try:
            abrev,species,taxid,compatible_mods = string.split(data,'\t')
        except Exception:
            if '!DOCTYPE': print_out = "A internet connection could not be established.\nPlease fix the problem before proceeding."
            else: print_out = "Unknown file error encountered."
            raw = export.ExportFile(fn); raw.close(); GO_Elite.importGOEliteParameters('skip'); sys.exit()
        if x==0: x=1
        else:
            compatible_mods = string.split(compatible_mods,'|')
            species_list.append(species)
            sd = SpeciesData(abrev,species,compatible_mods,taxid)
            species_codes[species] = sd
    return species_codes

def returnDirectoriesNoReplace(dir):
    dir_list = unique.returnDirectoriesNoReplace(dir); dir_list2 = []
    for entry in dir_list:
        if '.' not in entry: dir_list2.append(entry)
    return dir_list2

###### Command Line Functions (AKA Headless Mode) ######
def commandLineRun():
    import getopt
    #python GO_Elite.py --species Mm --mod Ensembl --permutations 2000 --method "combination" --zscore 1.96 --pval 0.05 --num 2 --input "C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample/input_list_small" --denom "C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample/denominator" --output "C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample"
    #python GO_Elite.py --species Mm --mod Ensembl --permutations 200 --method "combination" --zscore 1.96 --pval 0.05 --num 2 --input "C:/input" --denom "C:/denominator" --output "C:/output"
    #python GO_Elite.py --species Mm --input "C:/input" --denom "C:/denominator" --output "C:/output" --mod Ensembl --permutations 200
    global max_member_count; global filter_method; global z_threshold; global p_val_threshold; global change_threshold
    global sort_only_by_zscore; global permutations; global species; global root; global resources_to_analyze
    global run_mappfinder; global criterion_input_folder; global criterion_denom_folder; global main_output_folder; global custom_sets_folder
    global log_file; global summary_data_db; summary_data_db = {}; global analysis_method; global returnPathways; global mod; global imageType
    global commandLineMode; commandLineMode = 'yes'; global enrichment_type
    
    ###optional
    permutations = 2000; filter_method = 'z-score'; mod = 'Ensembl'; z_threshold = 1.96; enrichment_type = 'ORA'
    p_val_threshold = 0.05; change_threshold = 3; main_output_folder = ''; resources_to_analyze = 'all'

    ###required
    species_code = None
    species_full = None
    criterion_input_folder = None
    criterion_denom_folder = None
    custom_sets_folder = ''
    main_output_folder = None
    update_dbs = None
    update_method = []
    process_affygo='no'
    speciesfull=None
    update_ensrel = []
    replaceDB = 'no'
    species_taxid = ''
    incorporate_previous = 'yes'
    new_species = 'no'
    force = 'yes'
    ensembl_version = 'current'
    goelite_db_version = ''
    remove_download_files = 'no'
    export_versions_info = 'no'
    archive = 'no'
    additional_resources = [None]
    buildNested = 'no'
    continue_build = 'no'
    OBOurl = ''
    buildLocalWPFiles = 'yes'
    GOtype = 'GeneOntology'
    input_var = sys.argv[1:]
    nestTerms = 'no'
    wpid = None
    image_export = None
    returnPathways = None
    imageType = True
    resources=[]
    
    #input_var = ['--update', 'EntrezGene', '--species' ,'Mm'] ###used for testing

    if '--help' in input_var or '--h' in input_var:
        displayHelp()
    if '--version' in input_var or '--v' in input_var:
        print 'GO-Elite version 1.2.6 (http://genmapp.org/go_elite)'
        
    try:
        options, remainder = getopt.getopt(input_var,'', ['species=', 'mod=','permutations=',
                                    'method=','zscore=','pval=','num=','input=','denom=','output=',
                                    'update=','uaffygo=','system=', 'replaceDB=', 'speciesfull=',
                                    'taxid=', 'addspecies=', 'force=', 'version=', 'delfiles=',
                                    'archive=','exportinfo=','dataToAnalyze=','buildNested=',
                                    'customSet=','OBOurl=','GOtype=','nestTerms=','additional=',
                                    'image=','wpid=','returnPathways=','imageType=','enrichment='])
    except Exception, e: print e;sys.exit()
    
    for opt, arg in options:
        if opt == '--species': species_code=arg
        elif opt == '--mod': mod=arg
        elif opt == '--image': image_export=arg
        elif opt == '--wpid': wpid=arg
        elif opt == '--imageType':
            imageType=arg
            if imageType=='all':
                imageType = True
        elif opt == '--permutations': permutations=arg
        elif opt == '--method':
            filter_method=arg
            if filter_method == 'gene_number': filter_method = 'gene number'
            if filter_method == 'gene': filter_method = 'gene number'
            if 'z' in filter_method or 'Z' in filter_method: filter_method = 'z-score'
        elif opt == '--zscore': z_threshold=arg
        elif opt == '--pval': p_val_threshold=arg
        elif opt == '--num': change_threshold=arg
        elif opt == '--dataToAnalyze':
            if arg == 'local' or arg == 'WikiPathways':
                arg = 'Pathways'
            resources.append(arg)
        elif opt == '--input': criterion_input_folder=arg
        elif opt == '--denom': criterion_denom_folder=arg
        elif opt == '--output': main_output_folder=arg
        elif opt == '--enrichment':
            if string.lower(arg)=='ura' or 'under' in string.lower(arg):
                enrichment_type='URA'
        elif opt == '--update': update_dbs='yes'; update_method.append(arg)
        elif opt == '--additional':
            if additional_resources[0] == None:
                additional_resources=[]
                additional_resources.append(arg)
            else:
                additional_resources.append(arg)
        elif opt == '--uaffygo': process_affygo=arg
        elif opt == '--system': update_ensrel.append(arg) ### This is the related system. Multiple args to this flag are valid
        elif opt == '--version': ensembl_version = arg
        elif opt == '--replaceDB': replaceDB=arg
        elif opt == '--speciesfull': species_full=arg
        elif opt == '--taxid': species_taxid=arg
        elif opt == '--addspecies': update_dbs='yes'; new_species = 'yes'
        elif opt == '--force': force=arg
        elif opt == '--delfiles': remove_download_files = arg
        elif opt == '--archive': archive = arg
        elif opt == '--exportinfo': export_versions_info = arg
        elif opt == '--buildNested': buildNested = arg
        elif opt == '--customSet': custom_sets_folder = arg
        elif opt == '--OBOurl': OBOurl = arg
        elif opt == '--GOtype': GOtype = arg
        elif opt == '--nestTerms': nestTerms = arg
        elif opt == '--returnPathways': returnPathways = arg

    """ Build Database Outline:
    python GO_Elite.py --update Ensembl --system all --version 72 --species Hs --update WikiPathways --system EntrezGene
    GO_Elite.py --update EntrezGene --version 72 --species Hs
    python GO_Elite.py --update metabolites --version 72 --species Hs --force no
    ython GO_Elite.py --update Affymetrix --update WikiPathways --species Hs --replaceDB no --force no
    """

    if len(resources)>1: resources_to_analyze = resources
    elif len(resources)>0: resources_to_analyze = resources[0]
    
    species_full_original = species_full; species_code_original = species_code
    if image_export != None:
        if image_export == 'WikiPathways':
            #python GO_Elite.py --input /users/test/input/criterion1.txt --image WikiPathways --mod Ensembl --system arrays --species Hs --wpid WP536
            if wpid==None:
                print 'Please provide a valid WikiPathways ID (e.g., WP1234)';sys.exit()
            if species_code==None:
                print 'Please provide a valid species ID for an installed database (to install: --update Official --species Hs --version EnsMart62Plus)';sys.exit()
            if criterion_input_folder==None:
                print 'Please provide a valid file location for your input IDs (also needs to inlcude system code and value column)';sys.exit()
            from visualization_scripts import WikiPathways_webservice
            try:
                print 'Attempting to output a WikiPathways colored image from user data'
                print 'mod:',mod
                print 'species_code:',species_code
                print 'wpid:',wpid
                print 'imageType',imageType
                print 'input GO-Elite ID file:',criterion_input_folder
                graphic_link = WikiPathways_webservice.visualizePathwayAssociations(criterion_input_folder,species_code,mod,wpid,imageExport=imageType)
            except Exception,e:
                print traceback.format_exc()
                if 'force_no_matching_error' in traceback.format_exc():
                    print 'None of the input IDs mapped to this pathway' 
                elif 'IndexError' in traceback.format_exc():
                    print 'Input ID file does not have at least 3 columns, with the second column being system code'
                elif 'ValueError' in traceback.format_exc():
                    print 'Input ID file error. Please check that you do not have extra rows with no data' 
                elif 'source_data' in traceback.format_exc():
                    print 'Input ID file does not contain a valid system code' 
                else:
                    print 'Error generating the pathway "%s"' % wpid
        try: print 'Finished exporting visualized pathway to:',graphic_link['WP']
        except Exception: None ### Occurs if nothing output
        sys.exit()
        
    if 'EnsMart' in ensembl_version:
        import UI; UI.exportDBversion(ensembl_version)
    elif 'Plant' in ensembl_version:
        import UI; UI.exportDBversion(string.replace(ensembl_version,'Plant','EnsMart'))
    elif 'Bacteria' in ensembl_version:
        import UI; UI.exportDBversion(string.replace(ensembl_version,'Bacteria','EnsMart'))
    elif 'Fung' in ensembl_version:
        import UI; UI.exportDBversion(string.replace(ensembl_version,'Fungi','EnsMart'))

    program_type,database_dir = unique.whatProgramIsThis()
    if program_type == 'AltAnalyze': database_dir = '/AltDatabase'; goelite = '/goelite'
    else: database_dir = '/Databases'; goelite = ''

    if archive == 'yes':
        import update; import UI
        
        db_versions = UI.returnDirectoriesNoReplace(database_dir)

        for version in db_versions:
            print database_dir[1:]+'/'+version+goelite
            species_dirs = returnDirectoriesNoReplace(database_dir+'/'+version+goelite)
            print species_dirs
            for i in species_dirs:
                update.zipDirectory(database_dir[1:]+'/'+version+goelite+'/'+i); print 'Zipping',i

    if export_versions_info == 'yes':
        import UI; species_archive_db={}
        speciesData()
        db_versions = UI.returnDirectoriesNoReplace(database_dir)
        ### Export species names for each Official database version based on zip files in each folder
        #print db_versions
        for version in db_versions:
            #print version
            species_file_dirs = UI.returnFilesNoReplace(database_dir+'/'+version+goelite)
            #print species_dirs
            for file in species_file_dirs:
                if '.zip' in file:
                    file = string.replace(file,'.zip','')
                    if file in species_names:
                        species_name = species_names[file]
                        try: species_archive_db[species_name].append(version)
                        except Exception: species_archive_db[species_name] = [version]
        print 'len(species_archive_db)',len(species_archive_db)
        if len(species_archive_db)>0: UI.exportSpeciesVersionInfo(species_archive_db)

        ### Export array systems for each species for each Official database version
        species_array_db={}
        for version in db_versions:
            #print version
            species_dirs = UI.returnDirectoriesNoReplace(database_dir+'/'+version+goelite)
            for species_dir in species_dirs:
                supported_arrays=[]
                if species_dir in species_names:
                    species_name = species_names[species_dir]
                    species_file_dirs = UI.returnFilesNoReplace(database_dir+'/'+version+goelite+'/'+species_dir+'/uid-gene')
                    for file in species_file_dirs:
                        if 'Affymetrix' in file: supported_arrays.append('Affymetrix')
                        if 'MiscArray' in file: supported_arrays.append('MiscArray')
                        if 'Codelink' in file: supported_arrays.append('Codelink')
                        if 'Illumina' in file: supported_arrays.append('Illumina')
                        if 'Agilent' in file: supported_arrays.append('Agilent')
                    if len(supported_arrays)>0:
                        species_array_db[species_name,version] = supported_arrays
        print 'len(species_array_db)',len(species_array_db)
        if len(species_array_db)>0: UI.exportArrayVersionInfo(species_array_db)
    
    if replaceDB == 'yes': incorporate_previous = 'no'
    if update_dbs == 'yes' and ((species_code == None and species_full == None) and (update_method != ['Ontology'])) and update_method != ['metabolites']:
        print '\nInsufficient flags entered (requires --species or --speciesfull)'; sys.exit()
    elif (update_dbs == 'yes' or buildNested == 'yes') and (species_code != None or species_full != None or update_method == ['Ontology']):
    
        from import_scripts import BuildAffymetrixAssociations; import update; from build_scripts import EnsemblSQL; import UI
        file_location_defaults = UI.importDefaultFileLocations()
        speciesData()
        species_codes = UI.importSpeciesInfo()
        species_code_list=[]
        if len(species_codes) == 0:
            UI.remoteSpeciesInfo('no')
            species_codes = importSpeciesInfo() ### Gets the information from the backup version
    
        if ensembl_version != 'current' and 'release-' not in ensembl_version and 'EnsMart' not in ensembl_version:
            if 'Plant' not in ensembl_version and 'Fungi' not in ensembl_version:
                try: version_int = int(ensembl_version); ensembl_version = 'release-'+ensembl_version
                except ValueError: print 'The Ensembl version number is not formatted corrected. Please indicate the desired version number to use (e.g., "55").'; sys.exit()                
        
        if update_method == ['Ontology']: species_code_list=[]
        elif species_code == 'all':
            ### Add all species from the current database
            for species_code in species_names: species_code_list.append(species_code)
        elif species_full != None and species_full != 'all':
            species_full = [species_full] ###If the Ensembl species is not '' but is defined
        elif species_full == 'all':
            species_full = []
            child_dirs, ensembl_species, ensembl_versions = EnsemblSQL.getCurrentEnsemblSpecies(ensembl_version)
            for ens_species in ensembl_species:
                ens_species = string.replace(ens_species,'_',' ')
                species_full.append(ens_species)
        else: species_code_list = [species_code]
        if 'Official' in update_method:
            existing_species_codes = importSpeciesInfo()
            if len(existing_species_codes) == 0:
                UI.remoteSpeciesInfo('no')
                existing_species_codes = importSpeciesInfo() ### Gets the information from the backup version
                
            UI.getOnlineEliteConfig(file_location_defaults,'')
            if species_code != None:
                ### Integrate the speciescode
                speciesData()
                species_names_temp = UI.remoteSpeciesInfo('yes')
                if species_full == 'all': species_code_ls = species_names_temp
                else: species_code_ls = [species_code]
                for species_code in species_names_temp:
                    sd = species_names_temp[species_code]
                    existing_species_codes[sd.SpeciesName()] = sd
                UI.exportSpeciesInfo(existing_species_codes)
                speciesData()
        elif species_full != None:
            try: ensembl_species = ensembl_species
            except Exception: child_dirs, ensembl_species, ensembl_versions = EnsemblSQL.getCurrentEnsemblSpecies(ensembl_version)
            for species_ens in species_full:
                if species_ens in ensembl_species:
                    genus,species_code = string.split(species_ens,' ')
                    species_code = genus[0]+species_code[0]
                    taxid = ''; compatible_mods = ['En']
                    if species_code in species_names:
                        species = species_names[species_code]
                        sd = species_codes[species]
                        compatible_mods = sd.Systems()
                        taxid = sd.TaxID()
                        if species != species_ens: species_code = genus[:2]; species_names[species_code] = species_ens
                        elif 'En' not in compatible_mods: compatible_mods.append('En')
                    sd = UI.SpeciesData(species_code,species_ens,compatible_mods,taxid)
                    species_codes[species_ens] = sd
                    species_code_list.append(species_code) ### Add all Ensembl species (either all or single specified)
                else: print "The Ensembl species",species_ens,"was not found."; sys.exit()
            UI.exportSpeciesInfo(species_codes)
        species_code_list = unique.unique(species_code_list)

        if continue_build == 'yes' and (species_code_original == 'all' or species_full_original == 'all'):
            ### Only analyze species NOT in the directory
            current_species_ls = unique.read_directory('/Databases'); species_code_list2=[]
            for sc in species_code_list:
                if sc not in current_species_ls: species_code_list2.append(sc)
            species_code_list = species_code_list2
        
        species_iteration=0
        speciesData(); species_codes = importSpeciesInfo() ### Re-import the species data updated above
        if len(species_codes) == 0:
            UI.remoteSpeciesInfo('no'); species_codes = importSpeciesInfo() ### Gets the information from the backup version
        #print species_code_list, update_ensrel
        species_code_list.sort(); species_code_list.reverse()
        
        if 'WikiPathways' in update_method and buildLocalWPFiles == 'yes':
            import gene_associations; all_species = 'no'
            try:
                gene_associations.convertAllGPML(species_code_list,all_species) ### Downloads GPMLs and builds flat files
                null=[]
            except Exception:
                print 'Unable to connect to http://www.wikipathways.org'; sys.exit()

        for species_code in species_code_list:
            species_iteration +=1
            system_codes,system_list,mod_list = UI.remoteSystemInfo()
            try: species = species_names[species_code]
            except Exception: print 'Species code %s not found. Please add the species to the database.' % species_code; sys.exit()
            print "Starting to update databases for",species, string.join(update_method,',')
            
            ### Update EntrezGene-GO Databases
            if 'EntrezGene' in update_method:
                ncbi_go_file = file_location_defaults['EntrezGO'].Location(); status = 'null'
                if species_iteration == 1:
                    if force == 'yes': fln,status = update.download(ncbi_go_file,'BuildDBs/Entrez/Gene2GO/','txt')
                    else:
                        file_found = UI.verifyFile('BuildDBs/Entrez/Gene2GO/gene2go.txt')
                        if file_found == 'no':
                            fln,status = update.download(ncbi_go_file,'BuildDBs/Entrez/Gene2GO/','txt')
                if 'Internet' in status: print status
                else:
                    try: sd = species_codes[species]; species_taxid = sd.TaxID()
                    except KeyError: species_taxid = species_taxid
                    try: run_status = BuildAffymetrixAssociations.parseGene2GO(species_taxid,species_code,'over-write previous',incorporate_previous)
                    except Exception: run_status = 'no'
                    if run_status == 'run': print 'Finished building EntrezGene-GeneOntology associations files.'
                    else: print 'Gene2GO file not found. Select download to obtain database prior to extraction'
                    if remove_download_files == 'yes' and len(species_code_list)==1: export.deleteFolder('BuildDBs/Entrez/Gene2GO')
                    
                import gene_associations
                try: gene_associations.swapAndExportSystems(species_code,'Ensembl','EntrezGene') ### Allows for analysis of Ensembl IDs with EntrezGene based GO annotations (which can vary from Ensembl)
                except Exception: null=[] ### Occurs if EntrezGene not supported
                try: gene_associations.augmentEnsemblGO(species_code)
                except Exception: null=[] ### Occurs if EntrezGene not supported

            if new_species == 'yes':
                try:
                    species = species_full[0]
                    species_code = species_code
                    species_taxid = species_taxid
                except Exception:
                    print_out = 'Additional values are needed to add a new species'
                    print_out +='(e.g., --speciesfull "Homo sapiens" --species Hs --taxid 9606)'
                    print print_out; sys.exit()
                compatible_mods = ['L']
                sd = UI.SpeciesData(species_code,species,compatible_mods,species_taxid)
                species_codes[new_species_name] = sd
                UI.exportSpeciesInfo(species_codes)
                                        
            ### Download Official Databases
            if 'Official' in update_method:
                import UI
                try: os.remove(filepath('Databases/'+species_code+'/nested/version.txt'))
                except Exception: null=[] ### Remove old nested file
                
                buildNested = 'yes'
                UI.getOnlineEliteConfig(file_location_defaults,'')
                db_versions = UI.importOnlineDatabaseVersions(); db_version_list=[]
                for version in db_versions: db_version_list.append(version)
                db_version_list.sort(); db_version_list.reverse(); select_version = db_version_list[0]
                db_versions[select_version].sort()
                if ensembl_version != 'current':
                    if ensembl_version not in db_versions:
                        print ensembl_version, 'is not a valid version of Ensembl, while',select_version, 'is.'; sys.exit()
                    else: select_version = ensembl_version
                if species not in db_versions[select_version]:
                    print species, ': This species is not available for this version %s of the Official database.' % select_version
                else:
                    base_url = file_location_defaults['url'].Location()
                    #print base_url+'Databases/'+select_version+'/'+species_code+'.zip'
                    fln,status = update.download(base_url+'Databases/'+select_version+'/'+species_code+'.zip','Databases/','')
                    
                    ### Creates gene-Symbol.txt, EntrezGene-Ensembl.txt and augements gene-GO tables for between system analyses
                    UI.buildInferrenceTables(species_code)
                    
                    ### Attempt to download additional Ontologies and GeneSets
                    update_method.append('AdditionalResources')
                    
            ### Attempt to download additional Ontologies and GeneSets
            if 'AdditionalResources' in update_method:
                try:
                    from build_scripts import GeneSetDownloader
                    print 'Adding supplemental GeneSet and Ontology Collections'
                    if 'all' in additional_resources:
                        additionalResources = UI.importResourceList() ### Get's all additional possible resources
                    else: additionalResources = additional_resources
                    GeneSetDownloader.buildAccessoryPathwayDatabases([species_code],additionalResources,force)
                    print 'Finished adding additional analysis resources.'
                except Exception: print 'Download error encountered for additional Ontologies andGeneSets...\nplease try again later.'
                if 'Official' in update_method:
                    if 'Internet' not in status:
                        print 'Finished downloading the latest species database files.'
                    
            ### Download Ensembl Database
            if 'Ensembl' in update_method:
                externalDBName_list=[]
                try: child_dirs, ensembl_species, ensembl_versions = EnsemblSQL.getCurrentEnsemblSpecies(ensembl_version)
                except Exception: print "\nPlease try a different version. This one does not appear to be valid."; sys.exit()
                try: ensembl_sql_dir,ensembl_sql_description_dir = child_dirs[species]
                except Exception:
                    print traceback.format_exc()
                    print species,'species not supported in Ensembl'; continue
                ### Download the latest version of Ensembl
                EnsemblSQL.updateFiles(ensembl_sql_dir,'Config/','external_db.txt','yes')

                raw = export.ExportFile('Config/array.txt'); raw.close() ### Delete existing
                try: EnsemblSQL.updateFiles(string.replace(ensembl_sql_dir,'core','funcgen'),'Config/','array.txt','yes')
                except Exception: raw = export.ExportFile('Config/array.txt'); raw.close()
                external_dbs, external_system, array_db, external_ids = UI.importExternalDBs(species)
                if len(update_ensrel) == 0:
                    print "\nPlease indicate the system to update (e.g., --system all --system arrays --system Entrez)."; sys.exit()
                for i in update_ensrel:
                    i = string.replace(i,'\x93',''); i = string.replace(i,'\x94','') ### Occurs when there is a forward slash in the system name
                    if i in external_system: externalDBName_list.append(i)
                    elif i != 'all' and i != 'arrays':
                        print '\nEnsembl related system',[i], 'not found!!! Check the file Config/external_db.txt for available valid system names before proceeding.'; sys.exit()
                if 'all' in update_ensrel:
                    for i in external_system:  ### Build for all possible systems
                        if '\\N_' not in i: externalDBName_list.append(i)
                        #print [i]
                elif 'arrays' in update_ensrel:
                    externalDBName_list=[]
                    for array in array_db:
                        if '\\N_' not in array:
                            if 'ProbeLevel' in update_method:
                                if 'AFFY' in array: externalDBName_list.append(array) ### Ensures probe-level genomic coordinate assingment is restricted to Affy
                            else: externalDBName_list.append(array)

                import datetime
                today = str(datetime.date.today()); today = string.split(today,'-'); today = today[1]+'/'+today[2]+'/'+today[0]
                dirversion = string.replace(ensembl_version,'release-','EnsMart')
                OBO_import.exportVersionData(dirversion,today,'Config/')
                
                overwrite_previous = 'over-write previous'
                configType = 'Basic'; iteration=0
                from build_scripts import EnsemblSQL; reload(EnsemblSQL)
                if 'arrays' not in update_ensrel:
                    try: all_external_ids = EnsemblSQL.buildGOEliteDBs(species_code,ensembl_sql_dir,ensembl_sql_description_dir,'GO',configType,'GeneAndExternal',overwrite_previous,replaceDB,external_system,force); iteration+=1
                    except Exception, e:
                        print traceback.format_exc()
                        print 'Critical Error!!!! Exiting GO-Elite...'; sys.exit()
                    externalDBName_list_updated = UI.filterExternalDBs(all_external_ids,externalDBName_list,external_ids,array_db)
                    
                ###Add additional systems not in our annotated Config file if the user specified parsing of all systems or all array systems
                if 'arrays' in update_ensrel or 'all' in update_ensrel:
                    externalDBName_list = externalDBName_list_updated
                    
                for externalDBName in externalDBName_list:
                    if externalDBName != ' ':
                        if force == 'yes' and iteration == 1: force = 'no'
                        from build_scripts import EnsemblSQL; reload(EnsemblSQL)
                        
                        output_dir = 'BuildDBs/EnsemblSQL/'+species_code+'/'
                        if force == 'yes': ### Delete any existing data in the destination directory that can muck up tables from a new Ensembl build
                            export.deleteFolder(output_dir)
                            
                        if externalDBName in array_db:
                            #python GO_Elite.py --update Ensembl --update ProbeLevel --system arrays --version 65 --species Dr --force no
                            analysisType = 'FuncGen'
                            print [externalDBName], analysisType
                            if 'ProbeLevel' in update_method:
                                analysisType = 'ProbeLevel'
                                EnsemblSQL.buildGOEliteDBs(species_code,ensembl_sql_dir,ensembl_sql_description_dir,externalDBName,configType,analysisType,overwrite_previous,replaceDB,external_system,force); iteration+=1
                                #except Exception,e: print e;sys.exit()
                            else:
                                EnsemblSQL.buildGOEliteDBs(species_code,ensembl_sql_dir,ensembl_sql_description_dir,externalDBName,configType,analysisType,overwrite_previous,replaceDB,external_system,force); iteration+=1
                                #except Exception,e: print e;sys.exit()
                        else:
                            analysisType = 'GeneAndExternal'
                            print [externalDBName], analysisType
                            EnsemblSQL.buildGOEliteDBs(species_code,ensembl_sql_dir,ensembl_sql_description_dir,externalDBName,configType,analysisType,overwrite_previous,replaceDB,external_system,force); iteration+=1
                            #except Exception,e: print e;sys.exit()
                            
                try: swapAndExportSystems(species_code,'Ensembl','EntrezGene') ### Allows for analysis of Ensembl IDs with EntrezGene based GO annotations (which can vary from Ensembl)
                except Exception: null=[] ### Occurs if EntrezGene not supported
                
                if remove_download_files == 'yes': export.deleteFolder('BuildDBs/EnsemblSQL/'+species_code)
                
            if 'Affymetrix' in update_method or 'WikiPathways' in update_method:
                continue_analysis = 'no'
                if 'WikiPathways' in update_method:
                    relationship_types = ['native','mapped']
                    for relationship_type in relationship_types:
                        print 'Processing',relationship_type,'relationships'
                        index=0
                        if buildLocalWPFiles == 'yes':
                            status = '' ### Used when building a flat file from GPML zip file
                        else:
                            ### This is used when downloading a pre-built flat file from WikiPathways
                            UI.deleteWPFiles() ### Remove prior WP files
                            date = UI.TimeStamp(); file_type = ('wikipathways_'+date+'_'+species+'.tab','.txt')
                            if relationship_type == 'mapped':
                                url = 'http://wikipathways.org/wpi/cache/wikipathways_data_'+string.replace(species,' ','%20')+'.tab'
                                print url
                            else:
                                url = 'http://wikipathways.org/wpi/cache/wikipathways_native_data_'+string.replace(species,' ','%20')+'.tab'
                                print url
                            if force == 'yes':
                                fln,status = update.download(url,'BuildDBs/wikipathways/',file_type)
                                if 'Internet' in status: ### This file now loads for a minute or two, so this is common
                                    print 'Connection timed out... trying again in 15 seconds.'
                                    start_time = time.time(); end_time = time.time()
                                    while (end_time-start_time)<15: end_time = time.time()
                                    fln,status = update.download(url,'BuildDBs/wikipathways/',file_type)
                            else: status = 'null'
                        if 'Internet' not in status:
                            if 'Affymetrix' not in update_method: integrate_affy_associations = 'no'
                            else: integrate_affy_associations = 'yes'
                            counts = BuildAffymetrixAssociations.importWikipathways(system_codes,incorporate_previous,process_affygo,species,species_code,integrate_affy_associations,relationship_type,'over-write previous')
                            index+=1
                            if counts == 0: print 'No Affymetrix annotation files found and thus NO new results.'
                            else: print 'Finished parsing the latest Wikipathways and Affymetrix annotations.'
                        else: print status
                else:
                    try:
                        dir = '/BuildDBs/Affymetrix'; dir_list = UI.getFolders(dir)
                        if species_code in dir_list: continue_analysis = 'yes'
                    except IOError: continue_analysis = 'yes'
                    if continue_analysis == 'yes':
                        BuildAffymetrixAssociations.buildAffymetrixCSVAnnotations(species_code,incorporate_previous,process_affygo,'no','yes','over-write previous')
            if buildNested == 'yes':
                #try: os.remove(filepath('OBO/version.txt')) ### Maybe not a good idea?
                #except KeyError: null=[] ### Remove old nested file
                
                species_names_temp = UI.remoteSpeciesInfo('yes')
                if species_full == 'all' or species_code == 'all': species_code_ls = species_names_temp
                elif species_code != None: species_code_ls = [species_code]
                elif species_full != None:
                    species1,species2 = string.split(species_full,' ')
                    species_code_ls = [species1[0]+species2[0]]
                for species_code in species_code_ls:
                    current_species_dirs = unique.returnDirectories('/Databases')
                    if species_code in current_species_dirs:
                        try: export.deleteFolder(filepath('Databases/'+species_code+'/nested')) ### Delete the existing nested folder (will force rebuilding it)
                        except Exception: null=[]
                        ### Creates a nested GO (and stores objects in memory, but not needed
                        export_databases = 'no'; genmapp_mod = 'Ensembl'; sourceData()
                        print species_code,'Building Nested for mod types:',mod_types
                        avaialble_ontologies = OBO_import.findAvailableOntologies(species_code,mod_types)
                        for ontology_type in avaialble_ontologies:
                            try: full_path_db,path_id_to_goid,null = OBO_import.buildNestedOntologyAssociations(species_code,mod_types,ontology_type)
                            except Exception: None
                    try: UI.buildInferrenceTables(species_code)
                    except Exception: pass
                    
            if 'GORelationships' in update_method:            
                import gene_associations
                import datetime
                today = string.replace(str(datetime.date.today()),'-','')
                go_annotations = OBO_import.importPreviousOntologyAnnotations(GOtype)
                try:
                    #if GOtype != 'GOslim':
                    if nestTerms == 'yes': nested = 'nested'
                    else: nested = 'null'
                    ens_goslim = gene_associations.importGeneToOntologyData(species_code,'Ensembl',nested,ontology_type) ### The second to last last argument is gotype which can be null or nested
                    #else: ens_goslim = gene_associations.importUidGeneSimple(species_code,'Ensembl-goslim_goa')
                    print len(ens_goslim),'Ensembl-'+GOtype+ ' relationships imported...'
                    output_results = 'OBO/builds/'+species_code+'_'+GOtype+'_'+today+'.txt'
                    export_data = export.ExportFile(output_results)
                    title = string.join(['Ensembl','annotation.GO BIOLOGICAL_PROCESS','annotation.GO CELLULAR_COMPONENT','annotation.GO MOLECULAR_FUNCTION'],'\t')+'\n'
                    export_data.write(title)
                except Exception: ens_goslim=[] ### GO slim relationships not supported
                for gene in ens_goslim:
                    mf=[]; bp=[]; cc=[]
                    counts=0
                    for goid in ens_goslim[gene]:
                        try:
                            go_type = go_annotations[goid].GOType()
                            if go_type == 'molecular_function': mf.append(goid); counts+=1
                            elif go_type == 'biological_process': bp.append(goid); counts+=1
                            elif go_type == 'cellular_component': cc.append(goid); counts+=1
                        except Exception: null=[]
                    if counts>0:
                        mf = string.join(mf,','); bp = string.join(bp,','); cc = string.join(cc,',')
                        values = string.join([gene,bp,cc,mf],'\t')+'\n'
                        export_data.write(values)
                try: export_data.close()
                except Exception: null=[] ### If relationships not supported
              
        if 'Ontology' in update_method:
            UI.updateOBOfiles(file_location_defaults,'yes',OBOurl,'')

    elif criterion_input_folder != None and main_output_folder != None and species_code != None: #and criterion_denom_folder!= None
        file_dirs = criterion_input_folder, criterion_denom_folder, main_output_folder, custom_sets_folder
        parent = None
    
        time_stamp = timestamp() 
        log_file = filepath(main_output_folder+'/GO-Elite_report-'+time_stamp+'.log')
        try: log_report = open(log_file,'a'); log_report.close(); sys.stdout = Logger('')
        except Exception:
            print "Warning! The output directory must exist before running GO-Elite. Since it is easy to accidently specify an invalid output directory, it is best that GO-Elite does not create this for you.\n"
            sys.exit()
    
        if ensembl_version != 'current':
            ### Select database version otherwise use default
            goelite_db_version = ensembl_version
            import UI
            species_dirs = UI.returnDirectoriesNoReplace('/Databases')
            if goelite_db_version in species_dirs:
                import UI; UI.exportDBversion(goelite_db_version)
                print 'Using database version',goelite_db_version
                
        try: permutations = int(permutations);change_threshold = int(change_threshold)
        except Exception: permutations = permutations
        p_val_threshold = float(p_val_threshold); z_threshold = float(z_threshold)
        
        try: speciesData(); species = species_names[species_code]
        except Exception: print 'Species code not found. Please add the species to the database.'; sys.exit()
        
        max_member_count = 10000; sort_only_by_zscore = 'yes'; run_mappfinder = 'yes'
        #criterion_input_folder, criterion_denom_folder, main_output_folder, custom_sets_folder = file_dirs
        analysis_method = 'non-UI'
        change_threshold = change_threshold-1
        
        try:
            import UI
            species_dirs = UI.returnDirectoriesNoReplace('/Databases')
        except Exception:
            print '\nPlease install a species database (to install: python GO_Elite.py --update Official --species Hs --version EnsMart62Plus)';sys.exit()

        print ''
        root = parent; runGOElite(mod)
    else:
        print '\nInsufficient flags entered (requires --species, --input and --output)'; sys.exit()
    if 'metabolites' in update_method:
        import MetabolomicsParser
        MetabolomicsParser.buildMetabolomicsDatabase(force) ### will update any installed species
    sys.exit()
    
if __name__ == '__main__':
    try:
        import multiprocessing as mlp
        mlp.freeze_support()
    except Exception:
        print 'Note: Multiprocessing not supported for this verison python.'
        mlp = None
        
    run_parameter = 'intro'
    if Tkinter_failure == True:
        print "\nPmw or Tkinter not found... Tkinter print out not available"
        
    try: 
        ###### Determine Command Line versus GUI Control ######
        program_type,database_dir = unique.whatProgramIsThis()
        command_args = string.join(sys.argv,' ')
        if len(sys.argv[1:])>0 and '-' in command_args:
            ### Run in headless mode (Tkinter not required)
            commandLineRun()

        if use_Tkinter == 'yes':
            ### Run in GUI mode (Tkinter required)
            importGOEliteParameters(run_parameter)
    except Exception, exception:
            try:
                print_out = "Operating System Info: "+os.name +', '+str(platform.win32_ver()[0])+', '+str(platform.architecture())+', '+str(platform.mac_ver())+', '+str(platform.libc_ver())+', '+str(platform.platform())
            except Exception:
                print_out = "Operating System Info: "+os.name
            
            trace =  traceback.format_exc()
            print '\nWARNING!!! Critical error encountered (see below details)\n'
            print trace
            print "\n...exiting GO-Elite due to unexpected error (contact altanalyze@gmail.com for assistance)."
            time_stamp = timestamp()
            if commandLineMode == 'no':
                try: log_file = log_file
                except Exception:
                    try: log_file = filepath(main_output_folder+'/GO-Elite_report-error_'+time_stamp+'.log')
                    except Exception:
                        try: log_file = filepath('GO-Elite_report-error_'+time_stamp+'.log')
                        except Exception: log_file = filepath('/GO-Elite_report-error_'+time_stamp+'.log')
                            
                try: log_report = open(log_file,'a') ### append
                except Exception: log_report = open(log_file,'a') ### write
                log_report.write(print_out+'\n')
                log_report.write(trace)
                
                print_out = "Unknown error encountered during data processing.\nPlease see logfile in:\n\n"+log_file+"\nand report to altanalyze@gmail.com."
                if use_Tkinter == 'yes':
                    program,program_dir = unique.whatProgramIsThis()
                    if program!= 'AltAnalyze':
                        try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy()
                        except Exception: print print_out
                else: print print_out
                log_report.close()
                if len(log_file)>0:
                    if os.name == 'nt':
                          try: os.startfile('"'+log_file+'"')
                          except Exception:  os.system('open "'+log_file+'"')
                    elif 'darwin' in sys.platform: os.system('open "'+log_file+'"')
                    elif 'linux' in sys.platform: os.system('xdg-open "'+log_file+'/"')   
            sys.exit()
