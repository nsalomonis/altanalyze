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

import sys, string
import os.path, platform
import unique; reload(unique)
import math
import copy
import time
import shutil
import webbrowser
import gene_associations; reload(gene_associations)
import OBO_import; reload(OBO_import)
import export
import UI
import mappfinder; reload(mappfinder)
import traceback
use_Tkinter = 'no'
debug_mode = 'no'
start_time = time.time()

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
        file = getDirectoryFiles(self.data,str(search_term))
        if len(file)<1: print search_term,'not found'
        return file
    def searchdirectory_start(self,search_term):
        #self is an instance while self.data is the value of the instance
        files = getAllDirectoryFiles(self.data,str(search_term))
        match = ''
        for file in files:
            split_strs = string.split(file,'.')
            split_strs = string.split(split_strs[0],'/')
            if search_term  == split_strs[-1]: match = file
        return match
    
def getDirectoryFiles(import_dir, search_term):
    exact_file = ''; exact_file_dir=''
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
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
    try: built_go_paths, path_goid_db, path_dictionary = OBO_import.buildNestedGOTree('no','no')
    except Exception:
        ### This exception was added in version 1.2 and replaces the code in OBO_import.buildNestedGOTree which resets the OBO version to 0/0/00 and re-runs (see unlisted_variable = kill)
        print_out = "Unknown error encountered during Gene Ontology Tree import.\nPlease report to genmapp@gladstone.ucsf.edu if this error re-occurs."
        try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy(); sys.exit()
        except Exception: print print_out; sys.exit()
    
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
    dir_list = read_directory(input_dir)
    if (':' not in input_dir) and ('/Users/' not in input_dir[:7]) and ('Linux' in platform.system()): input_dir2 = input_dir
    else: input_dir2 = input_dir
    if len(dir_list)>0: ###if there are recently run files in the directory
        import datetime
        today = str(datetime.date.today()); today = string.split(today,'-'); today = today[0]+''+today[1]+''+today[2]
        
        time_stamp = string.replace(time.ctime(),':','')
        time_stamp = string.split(time_stamp,' ') ###Use a time-stamp as the output dir (minus the day)
        time_stamp = today+'-'+time_stamp[3]
        if 'archived-' not in input_dir2: ### We don't want users to select an archived directory and have the results moved to a child archive directory.
            archive_dir = input_dir2 +'/archived-'+time_stamp
            try: fn = filepath(archive_dir); os.mkdir(fn) ###create new directory
            except Exception: archive_dir = archive_dir[1:]; fn = filepath(archive_dir); os.mkdir(fn)
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
            if 'GOID' in data: x = 1; go_titles = data + '\n'; input_file_type = 'GO'
            elif 'MAPP Name' in data: y = 1; mapp_titles = data + '\n'; input_file_type = 'local'
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
            try: goid, go_name, go_type, local_changed, local_measured, local_go, percent_changed_local, percent_present_local, number_changed, number_measured, number_go, percent_changed, percent_present, z_score, permutep, adjustedp = string.split(data,'\t')
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
                        if z_score > z_threshold and number_changed > change_threshold and permutep < p_val_threshold:
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
                        mapp_name = t[0]; num_changed = t[1]; zscore = t[6]; permutep = t[7]
                        num_changed = int(num_changed); zscore = float(zscore); permutep = float(permutep)
                        #print [mapp_name, num_changed, zscore,permutep]
                        if zscore > z_threshold and num_changed > change_threshold and permutep < p_val_threshold:
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
    return path_dictionary_filtered

def comparison(zscore_changed_path_db):
    """This function performs the following operations:
    1) builds a database of all possible parent and child paths (e.g. nodes) from the list of filtered GO terms.
    2) identifies the most top-level nodes from 1 that will contain all nodes (where there's at one parent and one child).
    3) identifies all node without a parent or child in zscore_changed_path_db (paths filtered upon import)"""
    ###For each node in the hierarchy that has a child, list the parent as the key and any and all children (not specific branches) as the values.
    ###Each parent could be a child if it has a parent in another dictionary.  Nodes not in path_dictionary have no parents or children regulated.
    path_dictionary = buildFullTreePath(zscore_changed_path_db)  
    print '*',
    """Find the top parent nodes generated in the above path_dict query that are not in other dictionary
    items.  These represent unique top level parent nodes... the rest are redundant, since they will be children of the top parent"""
    path_dictionary_exclude = {}
    for redundant_path in path_dictionary:
        for unique_path in path_dictionary:
            if redundant_path!= unique_path:  ###If not looking at the same path
                if redundant_path in path_dictionary[unique_path]:
                    ####Therefore unique_path is the actual parent of redundant_path, example: redundant_path_index = (0, 3, 24, 30, 8) and unique_path_index = (0, 3)
                    redundant_children_paths = path_dictionary[redundant_path]
                    ###Store these so we can remove both the parent and child paths
                    path_dictionary_exclude[redundant_path] = redundant_children_paths  ###Probably don't need to store the redundant_children_paths (see next paragraph)
    print '* * * *',
    path_dictionary_include = {}; path_dictionary_include2 = {}
    ###Find the top level parents by removeing the rendundant entries from the path_dictionary (excluding these entries from a new dictionary)
    for parent_path in path_dictionary:
        if parent_path not in path_dictionary_exclude:
            unique_children_paths = path_dictionary[parent_path]
            path_dictionary_include['node',parent_path] = unique_children_paths ### The node annotations is used later on to destinguish this as a known top level parent node
            path_dictionary_include2[parent_path] = unique_children_paths  

    ###Now, identify those top-level nodes not found in the above paragraph, since they are regulated but have no other significant parents or
    ###child terms in the zscore_changed_path_db dictionary.  This will be likely listed in the output file.
    exclude = {}; zscore_changed_path_db_unique = {}
    for filtered_path in zscore_changed_path_db:
        for parent_path in path_dictionary_include2:
            if filtered_path in path_dictionary_include2[parent_path]: exclude[filtered_path]=filtered_path   
        if filtered_path in path_dictionary_include2: exclude[filtered_path]=filtered_path         
    for filtered_path in zscore_changed_path_db:
        if filtered_path not in exclude: zscore_changed_path_db_unique['node',filtered_path] = ['null']
    print '* * *',
    s = "\nRound 1: length of original GO index: "+str(len(zscore_changed_path_db))+'\n'; log_report.append(s)
    s = "Non-nested GO terms (unique): "+str(len(zscore_changed_path_db_unique))+'\n'; log_report.append(s)
    s = "Length of initial GO tree build: "+str(len(path_dictionary))+'\n'; log_report.append(s)
    s = "Length of GO tree node dictionary excluded: "+str(len(path_dictionary_exclude))+'\n'; log_report.append(s)
    s = "Length of GO tree node dictionary included: "+str(len(path_dictionary_include))+'\n'; log_report.append(s)

    iteration_count = 0  ###The next function iteratures many of the queries in this function, until optimized tree-hierarchies of nodes are created.  This states that no interations have occured yet

    return path_dictionary_include,zscore_changed_path_db_unique,iteration_count, path_dictionary
    
def comparison_iterate(dbase, finished_node_dbase,iteration_count):
    """This function is almost identical to "comparison" but handles two input databases from comparison and
    allows for multiple items in the key field of the input and output tables. NOTE: comparison1 could be
    adapted to serve this purpose"""
    #This code finds common parent nodes among children and parents and collapses into a single branch
    #node (i.e. dictionary).    
    path_dictionary = {}; path_dictionary_temp = {}
    for item in dbase:
        for entry in dbase[item]:   #unlike the 1st comparison function, query the dictionary definitions
            path_char = len(entry)  #get the number of characters to compare
            if item in dbase:        #compare back to the same entry
                for entry2 in dbase[item]:
                    if entry != entry2:  #looking for children of entry
                        base_path = entry2[0:path_char]
                        #print entry, "base_path:",base_path
                        if base_path == entry:
                            new_node_set=[]
                            for object in item:      #grab all individual key entries 
                                new_node_set.append(object)  #append the next child node in the new key
                            new_node_set.append(entry)
                            new_node_set = tuple(new_node_set)   #need to convert the list to a tuple to be a dictionary key
                            path_dictionary_temp[item[-1]] = ['temp']
                            try: path_dictionary[new_node_set].append(entry2) #append child of entry
                            except KeyError: path_dictionary[new_node_set] = [entry2]
    #need to only look at the last entry of dbase to find last unique entry
    dbase_temp = {}
    for entry in dbase: dbase_temp[entry[-1]] = ['temp']

    #find those nodes that had no more children this run and append to unique entries from last run
    dbase_unique = {}
    for entry in dbase:
        entry2 = entry[-1]   #we added the word "node" into the key at [0]
        #this adds entries that had no more nesting(one entry) to the finished dictionary as well
        #as parent entries added to path_dictionary
        if entry2 not in path_dictionary_temp:
            if len(dbase[entry]) == 1: #this should have no sibling entries below the parent
                finished_node_dbase[entry] = dbase[entry]
            else:  #therefore, there are siblings that need to be separated into distinct entries
                for node in dbase[entry]:
                    entry3 = entry[1:]
                    temp_list=[]
                    temp_list.append('sibling1')
                    for item in entry3:
                        temp_list.append(item)
                    temp_list.append(node)
                    new_key = tuple(temp_list)
                    finished_node_dbase[new_key] = ['siblings_type1']   #type1 specifies that the parent is unique, therefore no nested siblings of the parent term, just loners=> can analyze in calculate_score_for_children function
                    entry3 = entry[1:]
                    temp_list=[]
                    temp_list.append('sibling2')
                    for item in entry3:
                        temp_list.append(item)
                    new_key = tuple(temp_list)
                    finished_node_dbase[new_key] = dbase[entry] #also add an entry which contains both children under the true parent: note: instead of node, sibling is listed now to distinguish

    #find parent nodes generated in the above path_dict query that are not in other dictionary
    #items of path_dict.  These represent unique top level parent nodes... the rest are redundant
    path_dictionary_exclude = {}
    for item in path_dictionary:
        #item_parent = item[-2]
        item_child = item[-1]
        for item2 in path_dictionary:
            #item2_parent = item2[-2]
            #item2_child = item2[-1]
            if item != item2:
                if item_child in path_dictionary[item2]: #if the parent node is in the list of another entry, then it was prematurely added and should be excluded
                    parent_path = path_dictionary[item] #not really necessary to add the parent_path?
                    path_dictionary_exclude[item] = parent_path #exclude the parent item, prematurely added
    path_dictionary_include = {}
    for item in path_dictionary:
        if item not in path_dictionary_exclude:
            parent_path = path_dictionary[item]
            path_dictionary_include[item] = parent_path
    
    #print "\nRound 2: length of original GO tree:",len(dbase)
    #print "length of initial GO tree node dictionary:", len(path_dictionary)
    #print "number of non-nested GO terms:",len(finished_node_dbase)
    #print "length of GO tree node dictionary excluded:",len(path_dictionary_exclude)
    #print "length of GO tree node dictionary included:",len(path_dictionary_include)
            
    ### added this code to account for nodes that have a parent, no children, but other siblings with and without children
    ### previously these would have been eroneously eliminated
    nodes_examined = []
    for entry in path_dictionary_include:
        for node in path_dictionary_include[entry]:
             nodes_examined.append(node)
        for node in entry:
             if node != 'node':
                nodes_examined.append(node)
    for entry in finished_node_dbase:
        for node in finished_node_dbase[entry]:
             nodes_examined.append(node)
        for node in entry:
             if node != 'node':
                nodes_examined.append(node)
    nodes_examined = unique.unique(nodes_examined)
    for nodes in dbase: ###This is where it switched
        for entry in dbase[nodes]:
            if entry not in nodes_examined:
                temp_list = []
                for item in nodes:
                    temp_list.append(item)
                temp_list.append(entry)
                new_key = tuple(temp_list)
                #print "new_key",new_key
                finished_node_dbase[new_key] = ['null']
    if iteration_count < 2: print '* * * *',
    elif iteration_count < 3: print '* * *',
    if len(path_dictionary_include) > 0:
        iteration_count +=1
        """print "path_dictionary_include"
        for entry in path_dictionary_include:
            print entry, ":", path_dictionary_include[entry]
        print "finished_node_dbase"
        for entry in finished_node_dbase:
            print entry, ":", finished_node_dbase[entry] """
        return comparison_iterate(path_dictionary_include, finished_node_dbase, iteration_count)
    else:
        """for i in finished_node_dbase:
            if (0, 3, 24) in i: print i, finished_node_dbase[i]
        kill"""
        #print "Total number of iterations:",iteration_count
        print '* GO-tree collapsing finished'
        return chew_back(finished_node_dbase)

def chew_back(dbase):
    """This function takes a database of go nodes (e.g. [0.1.2,0.1.2.1,0.1.2.1.2]:[0.1.2.1.2.23]) and checks
    to see if the second to last parent node is found in any other dictionary keys.  If it is, then there
    is no need to transfer parent nodes to the children list, if no, then we must chew back."""
    iterate_dbase ={}; finished_dbase = {}
    #for entry in dbase: print entry,":", dbase[entry]
    
    original_dbase = copy.deepcopy(dbase)
    #print "Lenth of input chewback database:",len(dbase)
    #print "Looking for additional redundant terms across multiple branches...(please be patient)"
    loop = 0; count = 0; iter=0
    while len(dbase) > 0:
        loop +=1; original_increment = int(len(dbase)/20); increment = original_increment
        for entry in dbase:
            count = 0; iter+=1; parent_node = entry[-2]; unique_branch_node = entry[-1]
            if loop == 1:
                if iter == increment: increment+=original_increment; print '*',
            if len(entry)<3: #if there is only one parent node
                finished_dbase[entry] = dbase[entry]; count = 1
            else:
                for entry2 in original_dbase:
                    if entry != entry2 and unique_branch_node not in entry2:
                        if parent_node in entry2: finished_dbase[entry] = dbase[entry]; count = 1
            if count == 0:   #if the parent entry was not found in any of the dictionary keys, it's unique and thus, we need to chew back
                new_entry = entry[0:-1]
                children = dbase[entry]; last_entry = entry[-1]
                updated_children = [last_entry]+children
                iterate_dbase[new_entry] = updated_children
        dbase = iterate_dbase
        iterate_dbase = {}
        #print "number of chew-back iterations:",loop
        
        
    #print "Length of output chewback database:",len(finished_dbase)
    print "Database chewback complete"
    """for entry in finished_dbase: print entry,':',finished_dbase[entry] """
    return finished_dbase             

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
    if include_headers_in_output == 'yes': output_results = output_folder+'/'+mappfinder_input[0:-4]+ '_' +filter_method +'_elite-headers.txt'
    else: output_results = output_folder+'/'+mappfinder_input[0:-4]+ '_' +filter_method +'_elite.txt'
    entries_added=[]  #keep track of added GOIDs - later, add back un-added for MAPPFinder examination
    print_out = 'Results file: '+output_results+ '\nis open...can not re-write.\nPlease close file and select "OK".'
    try: raw = export.ExportFile(output_results)
    except Exception:
        try: UI.WarningWindow(print_out,' OK ')
        except Exception:
            print print_out; print 'Please correct (hit return to continue)'; inp = sys.stdin.readline()
    if include_headers_in_output == 'yes':
        for row in file_headers:
            row = row + '\n'
            raw.write(row) 
    if export_user_summary_values == 'yes' and len(value_headers)>0:
        value_headers2=[]; stdev_headers=[]
        for i in value_headers: value_headers2.append('AVG-'+i); stdev_headers.append('STDEV-'+i)
        value_headers = '\t'+string.join(value_headers2+stdev_headers,'\t')
    else: value_headers = ''
    if len(go_gene_annotation_db)>0: go_titles = go_titles[:-1]+'\t'+'redundant with terms'+'\t'+'inverse redundant'+'\t'+'gene symbols'+value_headers+'\n' ###If re-outputing results after building gene_associations
    raw.write(go_titles)
    combined_results[output_results] = [('go',mappfinder_input,go_titles)]
    
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
            combined_results[output_results].append(('go',mappfinder_input,data))
    if include_headers_in_output == 'yes':
        for goid in dummy_db:
            if goid not in entries_added:
                dummy_data = dummy_db[goid]
                raw.write(dummy_data)
    raw.close()

    combined_results[output_results].append(('','',''))
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
    combined_results[output_results] = [('mapp',mappfinder_input,mapp_titles)]

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
        combined_results[output_results].append(('mapp',mappfinder_input,line))
    raw.close()
    combined_results[output_results].append(('','',''))
    #print "Local Filtered MAPPFinder file", output_results, "written...."
    summary_data_db['redundant_mapp_term_count'] = len(mapps_redundant_with_others)
    
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
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if pathway_type == 'go':
                del t[0]; del t[2:7]
                go_type = t[1]; del t[1]; specific_pathway_type = go_type
            elif pathway_type == 'mapp': specific_pathway_type = 'local'
            else: specific_pathway_type = ''
            t.reverse(); t.append(specific_pathway_type); t.append(mapp_name); t.reverse()
            vals = string.join(t,'\t'); vals = vals + '\n'
            raw.write(vals)
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
    program_type,database_dir = unique.whatProgramIsThis()
    if program_type == 'GO-Elite': filename = 'Config/species.txt'
    else: filename = 'Config/goelite_species.txt'
    x=0
    fn=filepath(filename);global species_list; species_list=[]; global species_codes; species_codes={}
    global species_names; species_names={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t'); abrev=t[0]; species=t[1]
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
                    
def runGOElite(mod):
  #UI.WarningWindow('test','test')    
  source_data = mod; speciesData(); sourceDataCheck(); sourceData()
  try: species_code = species_codes[species]
  except KeyError: species_code=''

  global exclude_related; global mappfinder_input; global nested_collapsed_go_tree; global uid_to_go; global collapsed_go_tree
  global path_dictionary_include; global organized_tree; global parent_highest_score; global gene_identifier_file; global main_output_folder
  exclude_related = 'yes'
  nested_paths_stored = '' ###Variable used to determine if gene databases for GO analysis are imported for that species
  
  global program_type; global log_file
  program_type,database_dir = unique.whatProgramIsThis()
  global combined_results; combined_results={}; global go_to_mod_genes; global mapp_to_mod_genes; global zscore_goid
  global combined_associations; combined_associations={}; global combined_gene_ranking; combined_gene_ranking={}
  log_file = main_output_folder+'/GO-Elite_report.log'

  ### Add empty values for when there is no data on these values
  summary_keys = ['filtered_go_term_count','elite_go_term_count','redundant_go_term_count','filtered_local_count']
  summary_keys +=['redundant_mapp_term_count','go_gene_full_count','go_gene_elite_count','mapp_gene_elite_count']
  for key in summary_keys: summary_data_db[key]='0'

  if run_mappfinder=='yes':
          global path_id_to_goid;  global full_path_db; global all_nested_paths; global collapsed_tree; export_databases = 'no'
          genmapp_mod = 'Ensembl' ### default
          try:
              if resources_to_analyze != 'Pathways':
                  full_path_db,path_id_to_goid,null = OBO_import.buildNestedGOAssociations(species_code,export_databases,mod_types,genmapp_mod)
          except Exception:
              ### This exception was added in version 1.2 and replaces the code in OBO_import.buildNestedGOTree which resets the OBO version to 0/0/00 and re-runs (see unlisted_variable = kill)
              print_out = "Unknown error encountered during Gene Ontology Tree import.\nPlease report to genmapp@gladstone.ucsf.edu if this error re-occurs."
              try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy(); sys.exit()
              except Exception: print print_out; sys.exit()
          file_dirs = criterion_input_folder, criterion_denom_folder, main_output_folder
          go_to_mod_genes, mapp_to_mod_genes = mappfinder.generateMAPPFinderScores(species,species_code,source_data,mod,system_codes,permutations,resources_to_analyze,file_dirs,root)      
  else: import_path_index()
  global export_user_summary_values; global go_elite_output_folder
  get_non_elite_terms_only = 'no'
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
      d = "File in the folder nested are: "+str(dir_list); log_report.append(d+'\n')
      log_report.append('('+filepath('Databases/'+species_code+'/nested')+')\n')
      mappfinder_db_input_dir = '/'+species_code+'/nested/'
      nested_version = OBO_import.remoteImportVersionData(mappfinder_db_input_dir)
      #except Exception: nested_version = 'null'
      try: print os.name, platform.win32_ver()[0], platform.architecture(), platform.mac_ver(), platform.libc_ver(), platform.platform(), sys.getwindowsversion()[-1]
      except Exception: print os.name
      d = "Version of the nested files is: "+nested_version; log_report.append(d+'\n')
      nested_version = OBO_import.remoteImportVersionData('OBO/')
      d = "Version of the OBO files are: "+nested_version; log_report.append(d+'\n')
      exportLog(log_report); log_file = main_output_folder+'/GO-Elite_report.log'
      ### This exception was added in version 1.2 and replaces the code in OBO_import.buildNestedGOTree which resets the OBO version to 0/0/00 and re-runs (see unlisted_variable = kill)
      print_out = "Unknown error encountered during data processing.\nPlease see logfile in:\n\n"+log_file+"\nand report to genmapp@gladstone.ucsf.edu."
      try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy()
      except Exception: print print_out
      try: os.startfile('"'+log_file+'"')
      except Exception: null=[]
      sys.exit()

  global dataset_name; dataset_report=[]
  for mappfinder_input in dir_list:    #loop through each file in the directory to output results
    dataset_name = string.replace(mappfinder_input,'-GO.txt',''); dataset_name = string.replace(dataset_name,'-local.txt','')
    if dataset_name not in dataset_report:
        print 'Begining GO-Elite pruning for:',dataset_name; dataset_report.append(dataset_name)
    source_data = user_defined_source ###Resets the variables back to the user defined when missing from the MAPPFinder files
    mod = user_defined_mod
    mappfinder_input_dir = m.searchdirectory(mappfinder_input)
    try: run_mod,run_source,zscore_changed_path_db,go_full,go_titles,zscore_goid,input_file_type,gene_identifier_file,species_code = importMAPPFinderResults(mappfinder_input_dir)
    except Exception:
        print_out = 'Impropper ORA file format! Make sure to\n select the correct pre-processed input files.'
        try: print "Analysis Complete"; UI.WarningWindow(print_out,'Critical Error!')
        except Exception: print "Analysis Complete\n"
    if len(run_mod)>1: mod = run_mod
    if len(run_source)>1: source_data = run_source
    if input_file_type == 'GO':
        if run_mappfinder == 'yes':
            go_gene_full_count = countGOFullGenes(go_full,go_to_mod_genes); summary_data_db['go_gene_full_count'] = go_gene_full_count
        path_dictionary_include,zscore_changed_path_db_unique,iteration_count,all_paths = comparison(zscore_changed_path_db)
        organized_tree = comparison_iterate(path_dictionary_include,zscore_changed_path_db_unique,iteration_count)
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
            print "Begining to re-annotate pruned results..."
            if species_code != nested_paths_stored:  ###Otherwise these relationships are built already
                nested_paths_stored = species_code
                uid_to_go, uid_to_mapp, uid_system, gene_annotations = gene_associations.grabNestedGeneAssociations(species_code,mod,source_data)
            #else: print 'Apply Existing Gene Association Databases to GO Files'  ### Relevant species databases already exist (uid_to_go,gene_annotations,all_nested_paths)
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

            except KeyError:
                print_out = "\nWARNING: Could not delete parent GOID\n"+parent_goid+ "from go_count_db."
                print_out += '\nTypically occurs if the MOD or Source Type\nspecified in the MAPPFinder file is wrong (check this)\nor does not match the user selected MOD.'
                print_out += '\kReport bug to the GO-Elite help desk.'
                try: UI.WarningWindow(print_out,' Continue ')
                except Exception: print print_out
    else:
        if run_mappfinder == 'yes': mapp_gene_full_count = countGOFullGenes(zscore_goid,mapp_to_mod_genes); summary_data_db['mapp_gene_full_count'] = mapp_gene_full_count
        filtered_mapp_list = zscore_changed_path_db
        exportLocalResults(go_full,go_titles,{},{},{},{})
        summary_data_db['filtered_local_count'] =  len(go_full) ### All Elite GO-terms
        ###Identify gene lists for the corresponding GO-elite results and generate nested associations
        try: gene_file_dir = identifyGeneFiles(gene_input_dir, mappfinder_input)
        except Exception: gene_file_dir = ''
        if len(gene_file_dir) > 0:
            if species_code != nested_paths_stored:  ###Otherwise these relationships are built already
                nested_paths_stored = species_code
                uid_to_go, uid_to_mapp, uid_system, gene_annotations = gene_associations.grabNestedGeneAssociations(species_code,mod,source_data)
            #else: print 'Apply Existing Gene Association Databases to Local Files' ### Relevant species databases already exist (uid_to_go,gene_annotations,all_nested_paths)
            try:
                vals = gene_associations.matchInputIDsToMAPPEliteTerms(gene_file_dir,
                                        go_elite_output_folder,system_codes,mappfinder_input_dir,
                                        uid_to_mapp,filtered_mapp_list,gene_annotations,uid_system,
                                        combined_associations,combined_gene_ranking)
                combined_associations,combined_gene_ranking,mapp_gene_annotation_db,mapp_values_db,mapp_value_headers,mapps_with_redundant_genes,mapp_gene_elite_count = vals
                summary_data_db['mapp_gene_elite_count'] = mapp_gene_elite_count
                exportLocalResults(go_full,go_titles,mapp_gene_annotation_db,mapp_values_db,mapp_value_headers,mapps_with_redundant_genes)
            except KeyError:
                print_out = "\nWARNING: Could not delete parent GOID\n"+parent_goid+ "from go_count_db."
                print_out += '\nTypically occurs if the MOD or Source Type\nspecified in the MAPPFinder file is wrong (check this)\nor does not match the user selected MOD.'
                print_out += '\kReport bug to the GO-Elite help desk.'
                try: UI.WarningWindow(print_out,' Continue ')
                except Exception: print print_out
        
        if program_type != 'GO-Elite' and mappfinder_input[:3] == 'AS.':
            local_filename = go_elite_output_folder+'/'+mappfinder_input[0:-4]+ '_'+filter_method+'_elite.txt'
            #print 'Copying GO-Elite results to DomainGraph folder...'
            fn = filepath(local_filename)
            fn2 = string.replace(fn,'GO-Elite_results','DomainGraph')
            fn2 = string.replace(fn2,'GO-Elite','AltResults')
            fn2 = string.replace(fn2,'AS.','')
            fn2 = string.split(fn2,'-local'); fn2=fn2[0]+'-pathways-DomainGraph.txt'
            #shutil.copyfile(fn,fn2)
                
    ### Record results
    try:
        result_list = recordGOEliteStats()
        for d in result_list: log_report.append(d+'\n')
        log_report.append('\n')
    except Exception: null=[]
  print 'gene associations assigned'
  if len(combined_results)>0:
      output_combined_results(combined_results)
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

      try: moveMAPPFinderFiles(import_dir)
      except Exception,e:  print "Could not move ORA results... this will not impact this analysis:",e
      end_time = time.time(); time_diff = int(end_time-start_time)
      print_out = 'Analysis complete. GO-Elite results\nexported to the specified GO-Elite\noutput and MAPPFinder output directories.'
      try:
          if (program_type == 'GO-Elite' or analysis_method == 'UI') and use_Tkinter == 'yes':
              print 'Analysis Complete in %d seconds... Exiting GO-Elite' % time_diff
              UI.InfoWindow(print_out,'Analysis Completed!')
              tl = Toplevel(); SummaryResultsWindow(tl,main_output_folder)           
      except Exception:
          try: print 'Analysis Complete in %d seconds... Exiting GO-Elite' % time_diff
          except Exception: null=[]
  else:
      end_time = time.time(); time_diff = int(end_time-start_time)
      print_out = 'No input files to summarize!'
      if program_type == 'GO-Elite' and use_Tkinter == 'yes':
          print 'Analysis Complete in %d seconds... Exiting GO-Elite' % time_diff
          UI.WarningWindow(print_out,'Analysis Completed!')
          
  exportLog(log_report)
  run_parameter = 'skip'
  end_time = time.time(); time_diff = int(end_time-start_time)
  if program_type == 'GO-Elite' and use_Tkinter == 'yes': importGOEliteParameters(run_parameter)
  else:
      #try: import AltAnalyze; AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
      try: print 'Analysis Complete in %d seconds... Exiting GO-Elite' % time_diff
      except Exception: null=[]

class SummaryResultsWindow:
    def __init__(self,tl,output_dir):
        def showLink(event):
            idx= int(event.widget.tag_names(CURRENT)[1])
            webbrowser.open(LINKS[idx])
        url = 'http://www.genmapp.org/go_elite/help_main.htm'; url = filepath(url)
        LINKS=(url,'')
        self.LINKS = LINKS
        tl.title('GO-Elite 1.2 beta'); self.tl = tl
        
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

        result_list = recordGOEliteStats()
        for d in result_list: txt.insert(END, d+'\n');
            
        txt.tag_config('link', foreground="blue", underline = 1)
        txt.tag_bind('link', '<Button-1>', showLink)

        open_results_folder = Button(self.tl, text = 'Results Folder', command = self.openDirectory)
        open_results_folder.pack(side = 'left', padx = 5, pady = 5);

        self.dg_url = 'http://www.wikipathways.org'
        text_button = Button(self.tl, text='View WikiPathways', command=self.WPlinkout)
        text_button.pack(side = 'right', padx = 5, pady = 5)
        self.output_dir = output_dir
        self.whatNext_url = 'http://groups.google.com/group/go-elite/web/after-running-go-elite'
        whatNext_pdf = 'Documentation/what_next_goelite.pdf'; whatNext_pdf = filepath(whatNext_pdf); self.whatNext_pdf = whatNext_pdf
        
        what_next = Button(self.tl, text='What Next?', command=self.whatNextlinkout)
        what_next.pack(side = 'right', padx = 5, pady = 5)
        quit_buttonTL = Button(self.tl,text='Close View', command=self.close)
        quit_buttonTL.pack(side = 'right', padx = 5, pady = 5)
        
        continue_to_next_win = Button(text = 'Continue', command = root.destroy)
        continue_to_next_win.pack(side = 'right', padx = 10, pady = 10)
        quit_button = Button(root,text='Quit', command=self.quit)
        quit_button.pack(side = 'right', padx = 5, pady = 5)

        button_text = 'Help'; url = 'http://www.genmapp.org/go_elite/help_main.htm'; self.help_url = filepath(url)        
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
    def whatNextlinkout(self): self.GetHelpTopLevel(self.whatNext_url,self.whatNext_pdf)
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
        exportLog(log_report)
        try: self.tl.quit(); self.tl.destroy()
        except Exception: null=[]
        root.quit(); root.destroy(); sys.exit()
    def close(self):
        self.tl.quit()
        self.tl.destroy()

def recordGOEliteStats():
        result_list=[]
        for key in summary_data_db: summary_data_db[key] = str(summary_data_db[key])
        d = 'Dataset name: '+ dataset_name; result_list.append(d+'\n')
        d = summary_data_db['filtered_go_term_count']+':\tReported Filtered GO-terms'; result_list.append(d)
        d = summary_data_db['elite_go_term_count']+':\tReported Elite GO-terms'; result_list.append(d)
        d = summary_data_db['redundant_go_term_count']+':\tElite GO-terms redundant with other Elite GO-terms'; result_list.append(d)
        d = summary_data_db['filtered_local_count']+':\tReported Filtered Pathways'; result_list.append(d)
        d = summary_data_db['redundant_mapp_term_count']+':\tFiltered Pathways redundant with other Filtered Pathways'; result_list.append(d)

        if run_mappfinder == 'yes': 
            d = summary_data_db['go_gene_full_count']+':\tNumber of genes associated with Filtered GO-terms'; result_list.append(d)
            d = summary_data_db['go_gene_elite_count']+':\tNumber of genes associated with Elite GO-terms'; result_list.append(d)
            d = summary_data_db['mapp_gene_elite_count']+':\tNumber of genes associated with Filtered Pathways'; result_list.append(d)
        return result_list
    
def exportLog(log_report):
    log_file = main_output_folder+'/GO-Elite_report.log'
    fn=filepath(log_file); data = open(fn,'w')
    log_report = string.join(log_report,'')
    data.write(log_report); data.close()
    
class StatusWindow:
    def __init__(self,root,mod):

            self._parent = root
            root.title('GO-Elite 1.2 Beta - Status Window')
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
                
            Label(group.interior(),width=180,height=152,justify=LEFT, bg='black', fg = 'white',anchor=NW,padx = 5,pady = 5, textvariable=statusVar).pack(fill=X,expand=Y)

            status = StringVarFile(statusVar,root) ### Likely captures the stdout
            sys.stdout = status; root.after(100, runGOElite(mod))

            try: self._parent.mainloop()
            except Exception: null=[]
            try: self._parent.destroy()
            except Exception: null=[]

    def deleteWindow(self): tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)
    def quit(self): self._parent.quit(); self._parent.destroy(); sys.exit()
    
class StringVarFile:
    def __init__(self,stringVar,window):
        self.__newline = 0; self.__stringvar = stringVar; self.__window = window
    def write(self,s):
        log_report.append(s) ### Variable to record each print statement
        new = self.__stringvar.get()
        for c in s:
            #if c == '\n': self.__newline = 1
            if c == '\k': self.__newline = 1### This should not be found and thus results in a continous feed rather than replacing a single line
            else:
                if self.__newline: new = ""; self.__newline = 0
                new = new+c
        self.set(new)
    def set(self,s): self.__stringvar.set(s); self.__window.update()   
    def get(self): return self.__stringvar.get()
                
def importGOEliteParameters(run_parameter):
    global include_headers_in_output; global filter_method; global z_threshold; global p_val_threshold; global change_threshold
    global sort_only_by_zscore; global permutations; global species; global analysis_method; global resources_to_analyze
    global run_mappfinder; global criterion_input_folder; global criterion_denom_folder; global main_output_folder
    global log_report; log_report=[]; global summary_data_db; summary_data_db = {}
    parameters = UI.getUserParameters(run_parameter)
    species, run_mappfinder, mod, permutations, filter_method, z_threshold, p_val_threshold, change_threshold, resources_to_analyze, include_headers_in_output, file_dirs = parameters
    change_threshold = int(change_threshold); p_val_threshold = float(p_val_threshold); z_threshold = float(z_threshold)
    criterion_input_folder, criterion_denom_folder, main_output_folder = file_dirs
    try: permutations = int(permutations)
    except ValueError: permutations = 0
    sort_only_by_zscore = 'yes'; analysis_method = 'UI'
    
    global root; root=''
    if use_Tkinter == 'yes' and debug_mode == 'no':
        root = Tk()
        StatusWindow(root,mod)
        try: root.destroy()
        except Exception: null=[]
    else: runGOElite(mod)

def remoteAnalysis(variables,run_type):
    global include_headers_in_output; global filter_method; global z_threshold; global p_val_threshold; global change_threshold
    global sort_only_by_zscore; global permutations; global species; global root; global resources_to_analyze
    global run_mappfinder; global criterion_input_folder; global criterion_denom_folder; global main_output_folder
    global log_report; log_report=[]; global summary_data_db; summary_data_db = {}; global analysis_method

    species_code,mod,permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,file_dirs,parent = variables
    #print variables
    permutations = int(permutations); change_threshold = int(change_threshold); p_val_threshold = float(p_val_threshold)
    z_threshold = float(z_threshold); change_threshold = change_threshold-1
    speciesData()
    species = species_names[species_code]
    include_headers_in_output = 'no'; sort_only_by_zscore = 'yes'; run_mappfinder = 'yes'
    criterion_input_folder, criterion_denom_folder, main_output_folder = file_dirs
    if run_type == 'non-UI':
        analysis_method = run_type
        root = parent; runGOElite(mod)
    elif run_type == 'UI':
        analysis_method = run_type
        root = Tk()
        StatusWindow(root,mod)
        try: root.destroy()
        except Exception: null=[]

###### Command Line Functions (AKA Headless Mode) ######
def commandLineRun():
    import getopt
    #python GO_Elite.py --species Mm --mod Ensembl --permutations 2000 --method "combination" --zscore 1.96 --pval 0.05 --num 2 --input "C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample/input_list_small" --denom "C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample/denominator" --output "C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample"
    #python GO_Elite.py --species Mm --mod Ensembl --permutations 200 --method "combination" --zscore 1.96 --pval 0.05 --num 2 --input "C:/input" --denom "C:/denominator" --output "C:/output"
    #python GO_Elite.py --species Mm --input "C:/input" --denom "C:/denominator" --output "C:/output" --mod Ensembl --permutations 200
    global include_headers_in_output; global filter_method; global z_threshold; global p_val_threshold; global change_threshold
    global sort_only_by_zscore; global permutations; global species; global root; global resources_to_analyze
    global run_mappfinder; global criterion_input_folder; global criterion_denom_folder; global main_output_folder
    global log_report; log_report=[]; global summary_data_db; summary_data_db = {}; global analysis_method

    ###optional
    permutations = 2000; filter_method = 'z-score'; mod = 'EntrezGene'; z_threshold = 1.96
    p_val_threshold = 0.05; change_threshold = 3; main_output_folder = ''; resources_to_analyze = 'both'

    ###required
    species_code = None
    species_full = None
    criterion_input_folder = None
    criterion_denom_folder = None
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
    
    input_var = sys.argv[1:]
    #input_var = ['--update', 'EntrezGene', '--species' ,'Mm'] ###used for testing

    try:
        options, remainder = getopt.getopt(input_var,'', ['species=', 'mod=','permutations=',
                                    'method=','zscore=','pval=','num=','input=','denom=','output=',
                                    'update=','uaffygo=','system=', 'replaceDB=', 'speciesfull=',
                                    'taxid=', 'addspecies=', 'force=', 'version=', 'delfiles=',
                                    'archive=','exportinfo=','resources_to_analyze='])
    except Exception, e: print e;sys.exit()
    
    for opt, arg in options:
        if opt == '--species': species_code=arg
        elif opt == '--mod': mod=arg
        elif opt == '--permutations': permutations=arg
        elif opt == '--method': filter_method=arg
        elif opt == '--zscore': z_threshold=arg
        elif opt == '--pval': p_val_threshold=arg
        elif opt == '--num': change_threshold=arg
        elif opt == '--dataToAnalyze': resources_to_analyze=arg
        elif opt == '--input': criterion_input_folder=arg
        elif opt == '--denom': criterion_denom_folder=arg
        elif opt == '--output': main_output_folder=arg

        elif opt == '--update': update_dbs='yes'; update_method.append(arg)        
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
        
    print ''

    if 'EnsMart' in ensembl_version:
        import UI; UI.exportDBversion(ensembl_version)
        
    if archive == 'yes':
        import update; import UI
        db_versions = UI.returnDirectoriesNoReplace('/Databases')
        for version in db_versions:
            species_dirs = UI.returnDirectoriesNoReplace('/Databases/'+version)
            for i in species_dirs: update.zipDirectory('Databases/'+version+'/'+i); print 'Zipping',i

    if export_versions_info == 'yes':
        import UI; species_archive_db={}
        speciesData()
        db_versions = UI.returnDirectoriesNoReplace('/Databases')
        ### Export species names for each Official database version based on zip files in each folder
        print db_versions
        for version in db_versions:
            print version
            species_file_dirs = UI.returnFilesNoReplace('/Databases/'+version)
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
            print version
            species_dirs = UI.returnDirectoriesNoReplace('/Databases/'+version)
            for species_dir in species_dirs:
                supported_arrays=[]
                if species_dir in species_names:
                    species_name = species_names[species_dir]
                    species_file_dirs = UI.returnFilesNoReplace('/Databases/'+version+'/'+species_dir+'/uid-gene')
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
    if update_dbs == 'yes' and ((species_code == None and species_full) or update_method == ['GO']) == None:
        print '\nInsufficient flags entered (requires --species or --speciesfull)'; sys.exit()
    elif update_dbs == 'yes' and (species_code != None or species_full != None or update_method == ['GO']):
        import BuildAffymetrixAssociations; import update; import EnsemblSQL; import UI
        file_location_defaults = UI.importDefaultFileLocations()
        speciesData(); species_codes = UI.importSpeciesInfo(); species_code_list=[]

        if ensembl_version != 'current' and 'release-' not in ensembl_version and 'EnsMart' not in ensembl_version:
            try: version_int = int(ensembl_version); ensembl_version = 'release-'+ensembl_version
            except ValueError: print 'The Ensembl version number is not formatted corrected. Please indicate the desired version number to use (e.g., "55").'; sys.exit()                
        
        if update_method == ['GO']: species_code_list=[]
        elif species_code == 'all':
            ### Add all species from the current database
            for species_code in species_names: species_code_list.append(species_code)
        elif species_full != None and species_full != 'all':
            species_full = [species_full] ###If the Ensembl species is not 'all' but is defined
        elif species_full == 'all':
            species_full = []
            child_dirs, ensembl_species, ensembl_versions = EnsemblSQL.getCurrentEnsemblSpecies(ensembl_version)
            for ens_species in ensembl_species:
                ens_species = string.replace(ens_species,'_',' ')
                species_full.append(ens_species)
        else: species_code_list = [species_code]
        if 'Official' in update_method:
            existing_species_codes = UI.importSpeciesInfo()
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

        """
        ### Only analyze species NOT in the directory
        current_species_ls = unique.read_directory('/Databases'); species_code_list2=[]
        for sc in species_code_list:
            if sc not in current_species_ls: species_code_list2.append(sc)
        species_code_list = species_code_list2
        """
        species_iteration=0
        speciesData(); species_codes = UI.importSpeciesInfo() ### Re-import the species data updated above
        for species_code in species_code_list:
            species_iteration +=1
            system_codes,system_list,mod_list = UI.remoteSystemInfo()
            try: species = species_names[species_code]
            except Exception: print 'Species code %s not found. Please add the species to the database.' % species_code; sys.exit()
            print "Starting to update databases for",species
            
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
                
            ### Update Ensembl Databases
            if 'Official' in update_method:
                import UI
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
                    if 'Internet' not in status:
                        print 'Finished downloading the latest species database files.'
            
            if 'Ensembl' in update_method:
                externalDBName_list=['GO']
                try: child_dirs, ensembl_species, ensembl_versions = EnsemblSQL.getCurrentEnsemblSpecies(ensembl_version)
                except Exception: print "\nPlease try a different version. This one does not appear to be valid."; sys.exit()
                ensembl_sql_dir,ensembl_sql_description_dir = child_dirs[species]
                ### Download the latest version of Ensembl
                EnsemblSQL.updateFiles(ensembl_sql_dir,'Config/','external_db.txt','yes')
                try: EnsemblSQL.updateFiles(string.replace(ensembl_sql_dir,'core','funcgen'),'Config/','array.txt','yes')
                except Exception: raw = export.ExportFile('Config/array.txt'); raw.close()
                external_dbs, external_system, array_db = UI.importExternalDBs(species)
                for i in update_ensrel:
                    i = string.replace(i,'\x93',''); i = string.replace(i,'\x94','') ### Occurs when there is a forward slash in the system name
                    if i in external_system: externalDBName_list.append(i)
                    elif i != 'all' and i != 'arrays': print 'Ensembl related system',i, 'not found.'
                if 'all' in update_ensrel:
                    for i in external_system:  ### Build for all possible systems
                        if '\\N_' not in i: externalDBName_list.append(i)
                elif 'arrays' in update_ensrel:
                    externalDBName_list=[]
                    for array in array_db:
                        if '\\N_' not in array: externalDBName_list.append(array)
                
                import datetime
                today = str(datetime.date.today()); today = string.split(today,'-'); today = today[1]+'/'+today[2]+'/'+today[0]
                dirversion = string.replace(ensembl_version,'release-','EnsMart')
                OBO_import.exportVersionData(dirversion,today,'Config/')
                
                overwrite_previous = 'over-write previous'
                configType = 'Basic'; iteration=0
                for externalDBName in externalDBName_list:
                    if externalDBName != ' ':
                        if force == 'yes' and iteration == 1: force = 'no'
                        import EnsemblSQL; reload(EnsemblSQL)
                        if externalDBName in array_db:
                            analysisType = 'FuncGen'                        
                            EnsemblSQL.buildGOEliteDBs(species_code,ensembl_sql_dir,ensembl_sql_description_dir,externalDBName,configType,analysisType,overwrite_previous,replaceDB,external_system,force); iteration+=1
                            #except Exception,e: print e;sys.exit()
                        else:
                            analysisType = 'GeneAndExternal'                       
                            EnsemblSQL.buildGOEliteDBs(species_code,ensembl_sql_dir,ensembl_sql_description_dir,externalDBName,configType,analysisType,overwrite_previous,replaceDB,external_system,force); iteration+=1
                            #except Exception,e: print e;sys.exit()

                if remove_download_files == 'yes': export.deleteFolder('BuildDBs/EnsemblSQL/'+species_code)
                
            if 'Affymetrix' in update_method or 'WikiPathways' in update_method:
                continue_analysis = 'no'
                if 'WikiPathways' in update_method:
                    date = UI.TimeStamp(); file_type = ('wikipathways_'+date+'.tab','.txt')
                    url = 'http://www.wikipathways.org/wpi/pathway_content_flatfile.php?output=tab'
                    #url = 'http://137.120.14.24/wikipathways-test/wpi/pathway_content_flatfile.php?species=Homo%20sapiens&output=tab'
                    if species_iteration == 1 and force == 'yes':
                        fln,status = update.download(url,'BuildDBs/wikipathways/',file_type)
                        if 'Internet' in status: ### This file now loads for a minute or two, so this is common
                            fln,status = update.download(url,'BuildDBs/wikipathways/',file_type)
                    else: status = 'null'
                    if 'Internet' not in status:
                        if 'Affymetrix' not in update_method: integrate_affy_associations = 'no'
                        else: integrate_affy_associations = 'yes'
                        counts = BuildAffymetrixAssociations.importWikipathways(system_codes,incorporate_previous,process_affygo,species,species_code,integrate_affy_associations,'over-write previous')
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
        if 'GO' in update_method:
            UI.updateOBOfiles(file_location_defaults,'yes','')
    elif criterion_input_folder != None and criterion_denom_folder!= None and main_output_folder != None and species_code != None:
        file_dirs = criterion_input_folder, criterion_denom_folder, main_output_folder
        parent = None
    
        if ensembl_version != 'current':
            ### Select database version otherwise use default
            goelite_db_version = ensembl_version
            import UI
            species_dirs = UI.returnDirectoriesNoReplace('/Databases')
            if goelite_db_version in species_dirs:
                import UI; UI.exportDBversion(goelite_db_version)
                print 'Using database version',goelite_db_version
                
        permutations = int(permutations);change_threshold = int(change_threshold)
        p_val_threshold = float(p_val_threshold); z_threshold = float(z_threshold)
        
        try: speciesData(); species = species_names[species_code]
        except Exception: print 'Species code not found. Please add the species to the database.'; sys.exit()
        
        include_headers_in_output = 'no'; sort_only_by_zscore = 'yes'; run_mappfinder = 'yes'
        criterion_input_folder, criterion_denom_folder, main_output_folder = file_dirs
        analysis_method = 'non-UI'
        change_threshold = change_threshold-1
        print ''
        root = parent; runGOElite(mod)

    else:
        print '\nInsufficient flags entered (requires --species, --input, --denom and --output)'; sys.exit()
    sys.exit()
    
###### Determine Command Line versus GUI Control ######
command_args = string.join(sys.argv,' ')
if len(sys.argv[1:])>1 and '-' in command_args and 'AltAnalyze.py' not in command_args: commandLineRun()
else:
    try:
        import Tkinter
        from Tkinter import *
        import PmwFreeze
        use_Tkinter = 'yes'
    except ImportError: use_Tkinter = 'yes'; print "\nPmw or Tkinter not found... Tkinter print out not available";
    
if __name__ == '__main__':
    run_parameter = 'intro'
    if use_Tkinter == 'yes':
        try: importGOEliteParameters(run_parameter)
        except Exception, exception:

            try: print "Operating System Info:",os.name, platform.win32_ver()[0], platform.architecture(), platform.mac_ver(), platform.libc_ver(), platform.platform(), sys.getwindowsversion()[-1]
            except Exception: print os.name
            
            print traceback.format_exc()
            exportLog(log_report)
            print "\n...exiting GO-Elite due to unexpected error"
            try: log_file = log_file
            except Exception: log_file = ''            
            print_out = "Unknown error encountered during data processing.\nPlease see logfile in:\n\n"+log_file+"\nand report to genmapp@gladstone.ucsf.edu."
            try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy()
            except Exception: print print_out
            if len(log_file)>0:
                try: os.startfile('"'+log_file+'"')
                except Exception: null=[]
            sys.exit()
