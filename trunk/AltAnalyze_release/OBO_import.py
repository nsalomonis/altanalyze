###OBO_import
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

"""This module contains methods for reading in OBO format Gene Ontology files and building
numeric nested hierarchy paths (e.g., reconstructing the directed acyclic graph), importing
prebuilt hiearchy paths, creating nested GO associations from existing gene-GO files and
building GenMAPP format database tables (only by running this module directly)."""

import sys, string
import export
import os.path, platform
import unique
import math
import shutil
import time
import gene_associations
import copy

################# Parse directory files
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv" or ".ontology" in entry or '.obo' in entry: dir_list2.append(entry)
    return dir_list2

###### Classes ######
class GrabFiles:
    def setdirectory(self,value):
        self.data = value
    def display(self):
        print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        files = getDirectoryFiles(self.data,str(search_term))
        return files
    
def getDirectoryFiles(import_dir, search_term):
    matching_files = []
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir: matching_files.append(affy_data_dir)
    return matching_files

################# Import and Annotate Data
def eliminate_redundant_dict_values(database):
    db1={}
    for key in database: list = unique.unique(database[key]); list.sort(); db1[key] = list
    return db1

class GOPath:
    def __init__(self,goid,go_name,current_level,rank,path,GO_type):
        self._goid = goid; self._go_name = go_name; self._current_level = current_level
        self._rank = rank; self._path = path; self._GO_type = GO_type
    def GOID(self): return self._goid
    def GOIDStr(self): return self._goid[3:]
    def GOName(self): return self._go_name
    def GOLevel(self): return self._current_level
    def GOType(self): return self._GO_type
    def Rank(self): return self._rank
    def PathStr(self):
        path_index = pathToString(self.PathList())
        return path_index
    def PathList(self): return self._path
    def PathTuple(self): return tuple(self._path)
    def Report(self):
        output = self.GOID()+'|'+self.GOName()
        return output
    def __repr__(self): return self.Report()
    
def pathToString(path_list):
    path_str=[]
    for path_int in path_list: path_str.append(str(path_int))
    path_index = string.join(path_str,'.')
    return path_index

class GOTree:
    def __init__(self,goid,go_name,GO_type):
        self._goid = goid; self._go_name = go_name; self._GO_type = GO_type
    def GOID(self): return self._goid
    def GOName(self): return self._go_name
    def GOType(self): return self._GO_type
    def setGOType(self,gotype): self._GO_type=gotype
    def Report(self):
        output = self.GOID()+'|'+self.GOName()
        return output
    def __repr__(self): return self.Report()

class GOTreeDetailed(GOTree):
    ###Class not currently used
    def __init__(self,goid,go_name,GO_type,parent_goid,relation):
        self._goid = goid; self._go_name = go_name; self._GO_type = GO_type
        self._parent_goid = parent_goid; self._relation = relation
    def ParentGOID(self): return self._parent_goid
    def Relation(self): return self._relation
    def Report(self):
        output = self.GOID()+'|'+self.GOName()
        return output
    def __repr__(self): return self.Report()

###################################### UPDATED OBO CODE - BEGIN

def nestTree(parent_node,path,export_data,count_nodes):
    ### export_data,count_nodes are used for QC only
    children = edges[parent_node]
    path.append(0)
    for child in children.keys():
        tuple_path = tuple(path)
        #count_nodes+=1
        #try: temp = string.join(edges[child].keys(),'|')
        #except Exception: temp = ''
        #export_data.write(str(tuple_path)+'\t'+child+'\t'+temp+'\n')
        p = list(path)  ### Otherwise, the same path somehow gets used (alternative to copy.deepcopy())
        if child in edges:
            count_nodes = nestTree(child,p,export_data,count_nodes)
        #if count_nodes==1000: kill

        path_goid_db[tuple_path] = child
        if child not in built_go_paths:
            built_go_paths[child] = [tuple_path]
        elif tuple_path not in built_go_paths[child]:
            built_go_paths[child].append(tuple_path)
        path[-1]+=1
    return count_nodes
                                 
def importOBONew(filedir,path,GO_type,rank):
    if GO_type == '': discover_root = 'yes'
    else: discover_root = 'no'
    global edges
    #print [discover_root,GO_type,path]
    fn=filepath(filedir); x=0; stored={}; edges={}; category = 'null'; all_children={}; go_annotations_extra={}
    goid=''; go_name=''; edge_count=0; root_node = None
    for line in open(fn,'r').xreadlines():             
        data = cleanUpLine(line)
        s = string.split(data,' '); d = string.split(data,':')
        if x == 0:
            ###First few lines of files
            if 'date:' in data:
                build_date,version = extractDate(s[1])
                #print 'Current version of GO being parsed',build_date
                if version == previous_version:  ###If the current version is different (allows for any version differences to trigger a rebuild)
                    return 'old-version','old-version','old-version','old-version','old-version','old-version'
                    break
                else: exportVersionData(version,build_date,'OBO/')
                x=1
        if x > 0:
            #if s[0]=='def:': definition = d[1]
            if s[0]=='id:':
                try:
                    goid = s[1]
                    #goid=string.split(goid,':')[1]
                    category = 'null'
                except Exception: null=[]; goid = ''; go_name = ''
            if s[0]=='namespace:': category = s[1]
            if s[0]=='name:':
                go_name = d[1][1:]
                if go_name == GO_type:
                    root_node = goid
            if category == GO_type or discover_root=='yes':
                if s[0]=='is_a:': ### Note: sometimes there are multiple parents indicated for a single child
                    parent = s[1] ### immediate parent node
                    #parent=string.split(parent,':')[1]
                    if parent in edges: ### each child has one parent, one parent can have many children
                        children = edges[parent]
                        children[goid]=[]
                    else: children = {goid:[]}
                    edges[parent]=children
                    edge_count+=1
                    if discover_root=='yes': all_children[goid] = []
                    if goid not in go_annotations:
                        gt = GOTree(goid,go_name,GO_type)  
                        go_annotations[goid] = gt
                elif root_node == goid: ### For example, biological process
                    gt = GOTree(goid,go_name,GO_type)  
                    go_annotations[goid] = gt
                elif goid != '' and go_name != '':
                    gt = GOTree(goid,go_name,GO_type)
                    go_annotations_extra[goid] = gt
                    
    if discover_root=='yes':
        ### The root node should not exist as a child node
        for parent in edges:
            if parent not in all_children: root_node = parent
        GO_type = go_annotations_extra[root_node].GOName()
        print 'Parent node assigned to:',GO_type
        ### Assing the root_node name as the GO-Type
        for goid in go_annotations:
            go_annotations[goid].setGOType(GO_type)
    if root_node == None:
        print 'NO ROOT NODE IDENTIFIED... SHOULD BE:', GO_type; kill
    
    if len(path)==0: path.append(0); path_goid_db[tuple(path)] = root_node; return_path = list(path); #print [tuple(path)]
    else: path = [path[0]+1]; path_goid_db[tuple(path)] = root_node; return_path = list(path); #print [tuple(path)]
    print edge_count,'edges imported'
    print len(go_annotations), 'GO annotations'
    
    #export_data = export.ExportFile('OBO/test.txt')
    export_data=''
    nestTree(root_node,path,export_data,0)
    #export_data.close()
    #print 'Tree built'
    if run_mappfinder == 'yes':
        for path in path_goid_db:
            if mappfinder_version != current_version:
                path_dictionary[path]=[path]
                ###Build nested Path-index
                path_len = len(path); i=-1
                while path_len+i > 0:
                    parent_path = path[:i]
                    if parent_path in path_dictionary: path_dictionary[parent_path].append(path)
                    i-=1
                    
    print GO_type, 'parsed....'
    #print [[[return_path]]]
    return path_goid_db,built_go_paths,go_annotations,path_dictionary,return_path,rank
    
###################################### UPDATED OBO CODE - END

###### Begin Gene-Relationship Parsing ######
def importGOTree(filedir,path,rank):
    """Use the code designed to import OBO gene ontology DAG tree data for the GO-Elite mappfinder analysis to build
    the GenMAPP GO relationship tables.  At this stage, there is no species filtering"""
    fn=filepath(filedir); x=0; stored={}
    GO_type_name = string.replace(string.split(filedir,'/')[-1],'.ontology','')
    GO_type = string.upper(GO_type_name[0])
    if len(built_go_paths)<1: last_level=-1
    else: last_level=2
    for line in open(fn,'r').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x == 0:
            ###First few lines of files are annotations
            if data[0] != '!': x = 1
            ###Determine the database version (if newer than previous build, proceed)
            if '!date' in data:
                build_date,version = extractDate(data)
                #print 'Current version of GO being parsed',build_date
                if version == previous_version:  ###If the current version is different (allows for any version differences to trigger a rebuild)
                    return 'old-version','old-version','old-version','old-version','old-version','old-version'
                    break
                else: exportVersionData(version,build_date,'OBO/')
        if x > 0:
            goid,go_name,current_level,relation = parseLine(data)
            if current_level != -1:  ###Ignore very top level term (Gene_Ontology): Start number with next level down
                if current_level>last_level:
                    path.append(0) ###If the current level is lower than the previous, add a deeper (child) path
                    #print path, '1'
                elif current_level==last_level:
                    path[-1]+=1 ###If the current level is the same as the previous, increment number (distinct sibling to previous)
                    #print path, '2'
                elif current_level<last_level:
                    ###If the current level is higher than the previous, increment number of parent and remove child information
                    try: path = path[:current_level+1]; path[-1]+=1 ###remove child data and increment
                    except IndexError: print path, current_level, last_level, goid, go_name;kill
                    #print path, '3'
                #if go_name == 'mitochondrion': print goid,go_name;kill
                ###s = GOPath(goid,go_name,current_level,rank,path,GO_type)
                path = tuple(path)
                goid = goid[3:]
                if goid not in go_annotations:
                    s = GOTree(goid,go_name,GO_type)  
                    go_annotations[goid] = s
                try: built_go_paths[goid].append(path)
                except KeyError: built_go_paths[goid] = [path]
                path_goid_db[path] = goid
                
                if export_databases == 'yes':
                    ###identify the current goid's parent goid (needed for GenMAPP database export file)
                    try:parent_path = path[:-1]; parent_goid = path_goid_db[parent_path]
                    except KeyError: parent_goid = '' ###Occurs for all three top GO levels
                    if (goid,parent_goid) not in stored: ###don't add the same relationship twice
                        go_summary = [goid,go_name,GO_type,parent_goid,relation,'',build_date,'']
                        go_summary = string.join(go_summary,'\t')+'\n'; gofo.write(go_summary)
                        stored[(goid,parent_goid)] = []
                    tree_summary =  [str(rank+1),str(current_level+1), goid, go_name]
                    tree_summary = string.join(tree_summary,'\t')+'\n'; treefo.write(tree_summary)
                    
                if run_mappfinder == 'yes':
                    if mappfinder_version != current_version:
                        path_dictionary[path]=[path]
                        ###Build nested Path-index
                        path_len = len(path); i=-1
                        while path_len+i > 0:
                            parent_path = path[:i]
                            if parent_path in path_dictionary: path_dictionary[parent_path].append(path)
                            i-=1
                path = list(path)
                rank+=1; last_level = current_level
                
    print GO_type_name, 'GO terms parsed'
    return built_go_paths,path_goid_db,go_annotations,path_dictionary,path,rank

def exportGOCounts(built_go_paths):
    ###Make the GenMAPP GO-Counts table
    for goid in built_go_paths:
        goid_count = str(len(built_go_paths[goid]))
        count_summary =  [goid,goid_count]
        count_summary = string.join(count_summary,'\t')+'\n'; countfo.write(count_summary)
    countfo.close()
        
def parseLine(c):
    y=[]
    x = string.find(c,'%')
    if x>-1: y.append(x)
    x = string.find(c,'$')
    if x>-1: y.append(x)
    x = string.find(c,'<')
    if x>-1: y.append(x)
    y.sort(); level = y[0] ###indicates the number of spaces and thus the level of the GO-term
    c_list = string.split(c[level+1:],' ; '); goid = c_list[1]; go_name = c_list[0]
    relation = c[level]
    goid_list = string.split(goid,',') ###sometimes, there are more than one GOIDs per GO-term (synonym IDs). Only take the first
    goid = goid_list[0]
    goid_list = string.split(goid,' ') ### goid can often contain additional nested GO-terms downstream. Remove these
    goid = goid_list[0]
    return goid,go_name,level-1,relation ###-1 Since we eliminate Gene_Ontology as a level (for simplicity when joining GO Categories)

def extractDate(string_val):
    if ' ' in string_val:
        month_db = {}
        month_db['Jan'] = '1'; month_db['Feb'] = '2'; month_db['Mar'] = '3'; month_db['Apr'] = '4'; month_db['May'] = '5'
        month_db['Jun'] = '6'; month_db['Jul'] = '7'; month_db['Aug'] = '8'; month_db['Sep'] = '9'; month_db['Oct'] = '10'
        month_db['Nov'] = '11'; month_db['Dec'] = '12'
        date_info = string.split(string_val,' ')
        month = date_info[-5]; month = month_db[month]
        year = date_info[-1]
        day = date_info[-4]
        avg_days_in_month = 365/12.0
        version = int(((int(month)*avg_days_in_month)-avg_days_in_month)+int(day)+(int(year)*365.0)) ###The previous "version" was not precise. Instead use the date as # of days starting from Jan 1, 0000.
        date = month+'/'+day+'/'+year
    else:
        #19:05:2011
        day,month,year=string.split(string_val,':')
        month=str(int(month))
        date = month+'/'+day+'/'+year
        avg_days_in_month = 365/12.0
        version = int(((int(month)*avg_days_in_month)-avg_days_in_month)+int(day)+(int(year)*365.0))
    return date,version

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def remoteImportVersionData(dir):
    importVersionData(dir)
    return str(previous_version)

def importVersionData(dir):
    program_type,database_dir = unique.whatProgramIsThis()
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    elif 'OBO' in dir: parent_dir = ''
    else: parent_dir = database_dir
    original_dir = dir
    dir = parent_dir+dir
    global previous_version; global previous_date
    filename = dir+'version.txt'; fn=filepath(filename)
    try:
        for line in open(fn,'r').readlines():
            data = cleanUpLine(line)
            previous_version, previous_date = string.split(data,'\t')
            #print 'Last version of GO parsed',previous_date
    except Exception:
        #print original_dir,dir
        exportVersionData(0,'0/0/0',original_dir) ###Occurs if the version file is deleted... created a new one that requires a new build
        importVersionData(original_dir)
    previous_version = int(previous_version)

def exportVersionData(version,version_date,dir):
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    elif 'OBO' in dir or 'Config' in dir: parent_dir = ''
    else: parent_dir = database_dir
    dir = parent_dir+dir    
    global current_version; current_version = version
    global current_version_date; current_version_date = version_date
    new_file = dir+'version.txt'
    data = export.ExportFile(new_file)
    data.write(str(version)+'\t'+str(version_date)+'\n'); data.close()
    
def swapKeyValues(db):
    swapped={}
    for key in db:
        values = list(db[key]) ###If the value is not a list, make a list
        for value in values:
            try: swapped[value].append(key)
            except KeyError: swapped[value] = [key]
    swapped = eliminate_redundant_dict_values(swapped)
    return swapped

def exportCurrentGOBuild(path_goid_db,go_annotations):
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    new_file = parent_dir+'OBO/builds/built_go_paths.txt'
    try: fn=filepath(new_file); data = open(fn,'w')
    except Exception:
        new_dir = parent_dir+'OBO/builds'; fn = filepath(new_dir)
        os.mkdir(fn) ###Re-Create directory if deleted
        fn=filepath(new_file); data = open(fn,'w')
    data.write('Path'+'\t'+'GOID'+'\n')
    for path in path_goid_db:
        goid = path_goid_db[path]; path = pathToString(path)
        data.write(path +'\t'+ goid +'\n')
    data.close()

    new_file = parent_dir+'OBO/builds/go_annotations.txt'
    fn=filepath(new_file); data = open(fn,'w')
    data.write('GOID'+'\t'+'GO Name'+'\t'+'GO Type'+'\n')    
    for goid in go_annotations:
        s = go_annotations[goid]
        data.write(goid +'\t'+ s.GOName() +'\t'+ s.GOType() +'\n')
    data.close()

def convertStrListToIntList(path):
    path_int=[]
    for str in path: path_int.append(int(str))
    return path_int

def importPreviousGOAnnotations():
    go_annotations={}
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    filename = parent_dir+'OBO/builds/go_annotations.txt'; fn=filepath(filename); x=0
    for line in open(fn,'r').xreadlines():
        if x==0: x=1 ###Skip the title line
        else:
            data = cleanUpLine(line)
            goid,go_name,go_type = string.split(data,'\t')
            if ':' not in goid: goid = 'GO:'+goid
            if go_name[0]== ' ': go_name = go_name[1:]
            s = GOTree(goid,go_name,go_type)
            go_annotations[goid] = s
    return go_annotations

def importPreviousGOBuild():
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    filename = parent_dir+'OBO/builds/built_go_paths.txt'; fn=filepath(filename); x=0; count=0
    importVersionData('OBO/'); global current_version; current_version = previous_version
    global current_version_date; current_version_date = previous_date ### Needed for buildNestedGOAssociations
    for line in open(fn,'r').xreadlines(): count+=1
    original_increment = int(count/20); increment = original_increment
    
    try: ### This reduces run-time for the typical analysis where the databases are in sync and up-to-date
        if run_mappfinder == 'yes':
            if (mappfinder_version != current_version and current_version != 0) or verified_nested == 'no':
                build_nestedDB='yes'
            else: build_nestedDB = 'no'
        else: build_nestedDB = 'no'
    except Exception: build_nestedDB = 'no'
            
    for line in open(fn,'r').xreadlines():
        if x==0: x+=1 ###Skip the title line
        else:
            x+=1
            if x == increment: increment+=original_increment; print '*',    
            data = cleanUpLine(line)
            path,goid = string.split(data,'\t')
            path = tuple(map(int,string.split(path,'.')))
            #path = string.split(path_str,'.'); path = convertStrListToIntList(path); path = tuple(path)
            #s = GOPath(goid,'','','',path,''); s = GOPathAbr(goid,path)
            if ':' not in goid: goid = 'GO:'+goid
            path_goid_db[path] = goid
            try: built_go_paths[goid].append(path)
            except KeyError: built_go_paths[goid] = [path]
            if build_nestedDB == 'yes':
                if mappfinder_version != current_version:
                    path_dictionary[path]=[path]
    ###All of the paths need to be added before  
    if build_nestedDB == 'yes':
        if mappfinder_version != current_version:
            for path in path_dictionary:
                ###Build nested Path-index
                path_len = len(path); i=-1
                while path_len+i > 0:
                    parent_path = path[:i]
                    try: path_dictionary[parent_path].append(path)
                    except Exception: null=[]
                    i-=1    
    
#### Import gene data and associate with Nested GO
def grabNestedGOIDs():
    nested_go_tree={}
    for path in path_dictionary:
        parent_goid = path_goid_db[path]
        child_goid_list=[]
        for child_path in path_dictionary[path]:
            child_goid = path_goid_db[child_path]; child_goid_list.append(child_goid)
        child_goid_list = unique.unique(child_goid_list)
        nested_go_tree[parent_goid] = child_goid_list
    return nested_go_tree
    
def linkGenesToNestedGO(go_to_gene):
    nested_go_genes={}; made_unique={}; x=0
    original_increment = int(len(nested_go_tree)/20); increment = original_increment    

    for parent_goid in nested_go_tree:
        x+=1
        if x == increment: increment+=original_increment; print '*',
        for child_goid in nested_go_tree[parent_goid]: ### This list of goids includes the parent, since it is the first entry in path_dictionary
            if child_goid in go_to_gene:
                ensembls=go_to_gene[child_goid]
                for ensembl in ensembls:
                    try:
                        ens_db = nested_go_genes[parent_goid]
                        ens_db[ensembl] = ''
                    except KeyError:
                        ens_db = {}; ens_db[ensembl] = ''; e = ens_db
                        nested_go_genes[parent_goid] = e 
    return nested_go_genes

def exportGORelationships(nested_go_gene,gene_to_source_id,mod,source_type):
    program_type,database_dir = unique.whatProgramIsThis()
    new_file = database_dir+'/'+species_code+'/nested/'+mod+'_to_Nested-GO.txt'
    data = export.ExportFile(new_file)
    title = [mod,'GOID']; title_str = string.join(title,'\t')
    data.write(title_str+'\n')
    for goid in nested_go_gene:
        for gene in nested_go_gene[goid]:
            output_list = [gene,goid]
            output_str = string.join(output_list,'\t')
            data.write(output_str+'\n')
    data.close()
    print new_file, 'saved to disk'

def setGenMAPPexportGlobals():
    global gofo; global treefo; global countfo
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    of = parent_dir+'OBO/builds/genmapp_go.txt';of=filepath(of); gofo = open(of,'w') ###create output for GenMAPP GeneOntology table database
    oft = parent_dir+'OBO/builds/genmapp_go_tree.txt'; oft=filepath(oft); treefo = open(oft,'w') ###create output for GenMAPP GeneOntology table database
    ofc = parent_dir+'OBO/builds/genmapp_go_count.txt'; ofc=filepath(ofc);countfo = open(ofc,'w') ###create output for GenMAPP GeneOntology table database
    gofo_title = ['ID', 'Name', 'Type', 'Parent', 'Relation', 'Species', 'Date', 'Remarks']
    treefo_title = ['OrderNo', 'Level', 'ID', 'Name']
    countfo_title = ['ID', 'Count']
    gofo_title = string.join(gofo_title,'\t')+'\n'; gofo.write(gofo_title)
    treefo_title = string.join(treefo_title,'\t')+'\n'; treefo.write(treefo_title)
    countfo_title = string.join(countfo_title,'\t')+'\n'; countfo.write(countfo_title)
    return of,oft,ofc

def reOrganizeGenMAPPFiles(of,oft,ofc):
    ###copy files to different filename
    new_of = string.replace(of,'go','go_table'); shutil.copyfile(of, new_of)
    new_oft = string.replace(oft,'go_tree','go_tree_table'); shutil.copyfile(oft, new_oft)
    new_ofc = string.replace(ofc,'go_count','go_count_table'); shutil.copyfile(ofc, new_ofc)

#### Main functions that grab data from above functions
def buildNestedGOTree(mappfinder,export_dbases):
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    
    global run_mappfinder; run_mappfinder = mappfinder
    global export_databases; export_databases = export_dbases ###declaring the global here makes it applicable everywhere else
    ###Import all the OBO GO tree information from http:/www.geneontology.org/
    import_dir = '/'+parent_dir+'OBO'; global path_dictionary; path_dictionary={}; global built_go_paths; built_go_paths={}
    global go_annotations; go_annotations={} ###global built_go_paths; built_go_paths={};
    global path_goid_db; path_goid_db={}; global GO_version; path=[]; rank=0
    c = GrabFiles(); c.setdirectory(import_dir); file_dirs = c.searchdirectory('.ontology'); file_dirs.reverse()
    x = file_dirs[1:]+file_dirs[0:1] ###Reorganize to mimic Steve Lawlors program
    if len(file_dirs)==0:
        file_dirs = c.searchdirectory('.obo'); file_dirs.reverse(); x = file_dirs[1:]+file_dirs[0:1]
    importVersionData('OBO/')
    if export_databases == 'yes': of,oft,ofc = setGenMAPPexportGlobals()
    start_time = time.time()
    for file_dir in file_dirs:
        ###Import the 3 main GO files and index them so that the first path corresponds to the Ontology type - Software checks the date before parsing
        if '.obo' in file_dir or '.ontology' in file_dir:
            if 'gene_ontology' in file_dir or 'goslim' in file_dir:
                ### This is the latest supported format
                path_goid_db,built_go_paths,go_annotations,path_dictionary,path,rank = importOBONew(file_dir,path,'biological_process',rank)
                #print len(path_goid_db),len(built_go_paths),len(go_annotations)
                try: path_goid_db,built_go_paths,go_annotations,path_dictionary,path,rank = importOBONew(file_dir,path,'molecular_function',rank)
                except Exception: null=[]
                #print len(path_goid_db),len(built_go_paths),len(go_annotations)
                path_goid_db,built_go_paths,go_annotations,path_dictionary,path,rank = importOBONew(file_dir,path,'cellular_component',rank)
                #print len(path_goid_db),len(built_go_paths),len(go_annotations)
            else:
                path_goid_db,built_go_paths,go_annotations,path_dictionary,path,rank = importOBONew(file_dir,path,'',rank)
        else:
            ### This is the old flat file format which is no longer supported
            built_go_paths,path_goid_db,go_annotations,path_dictionary,path,rank = importGOTree(file_dir,path,rank)
        if built_go_paths == 'old-version': break
    if built_go_paths == 'old-version' or program_type == 'AltAnalyze' or len(file_dirs)==0:
        GO_version = 'old'
        print "Importing previously parsed version of GO"; go_annotations={}; built_go_paths={}; path_dictionary={}; path_goid_db={}
        try:
            importPreviousGOBuild()
            if export_databases == 'yes':
                gofo.close(); treefo.close(); countfo.close()
                os.remove(of);os.remove(oft);os.remove(ofc) ###If export file is built with a new version, don't need this empty file when parsing an old version
        except IOError: ### Occurs when the OBO/builds directory is deleted but the version file still exists
            """
            print "build directory not found. Re-importing GO annotations"
            exportVersionData(0,'0/0/0','OBO/'); current_version = -1
            buildNestedGOTree(mappfinder,export_dbases) ### Re-initialize current function
            """
            unlisted_variable = kill
    else:
        GO_version = 'new'
        exportCurrentGOBuild(path_goid_db,go_annotations)
        if export_databases == 'yes':
            gofo.close(); treefo.close()
            exportGOCounts(built_go_paths)
            reOrganizeGenMAPPFiles(of,oft,ofc)
            os.remove(of);os.remove(oft);os.remove(ofc) ###If export file is built with a new version, don't need this empty file when parsing an old version
    end_time = time.time(); time_diff = int(end_time-start_time)

    print "GO categories imported and nested in %d seconds" % time_diff
    return built_go_paths, path_goid_db, path_dictionary

def verifyFileLength(filename):
    count = 0
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            count+=1
            if count>3: break
    except Exception: null=[]
    return count

def verifyNestedFileCreation(species,mod_types):
    ### Determine which mods are present for GO
    mods_present = []; nested_present=[]; verified = 'no'
    for mod in mod_types:
        go_file = 'Databases/'+species+'/gene-go/'+mod+'-GeneOntology.txt'
        count = verifyFileLength(go_file) ### See if there are lines present in the file (if present)
        if count>1: mods_present.append(mod)
    if len(mods_present)>0:
        for mod in mods_present:
            go_file = 'Databases/'+species+'/nested/'+mod+'_to_Nested-GO.txt'
            count = verifyFileLength(go_file) ### See if there are lines present in the file (if present)
            if count>1: nested_present.append(mod)
        if len(nested_present) == len(mods_present): verified = 'yes'
    return verified
        
def buildNestedGOAssociations(species,export_dbases,mod_types,genmapp_mod):
    global species_code; species_code = species; global verified_nested
    ###Check to see if a current version of the nested gene databases exists
    if ('Linux' in platform.system()): mappfinder_db_input_dir = species_code+'/nested/'
    else: mappfinder_db_input_dir = '/'+species_code+'/nested/'    
    verified_nested = verifyNestedFileCreation(species,mod_types) ### This is in addition to looking at the version files (in case a problem occurs with these)
    if verified_nested == 'no':
        exportVersionData(0,'0/0/0',mappfinder_db_input_dir)
    importVersionData(mappfinder_db_input_dir)
    
    global mappfinder_version; mappfinder_version = previous_version
    built_go_paths, path_goid_db, path_dictionary = buildNestedGOTree('yes',export_dbases)
    
    if (mappfinder_version != current_version and current_version != 0) or verified_nested == 'no':  ### modified this code such that any version change warrants a rebuild and if reset by BuildEntrezAffymetrixAssociations or other, that it triggers a rebuild
        print 'Building GO nested gene association files for',species_code
        ###Build Gene to GO associations for all MODs and export these for re-import by
        ###the MAPPFinder module
        global nested_go_tree; nested_go_tree = grabNestedGOIDs()
        for mod in mod_types:
            try:
                start_time = time.time()
                mod_to_go = gene_associations.importGeneGOData(species_code,mod,'null')
                go_to_mod = swapKeyValues(mod_to_go); total_gene_count = len(mod_to_go); mod_to_go=[]
                ###Obtain a database of GOIDs with all nested gene associations
                nested_go_mod = linkGenesToNestedGO(go_to_mod)
                exportGORelationships(nested_go_mod,{},mod,'')
                end_time = time.time(); time_diff = int(end_time-start_time)
                exportVersionData(current_version,current_version_date,mappfinder_db_input_dir)
                print "GO Nested Lists Process/Created in %d seconds" % time_diff
                if export_databases == 'yes' and mod == genmapp_mod:
                    ###Read-In GenMAPP GO databases (built in ImportGOTree) and filter to make species specific using nested GO-association
                    ###Do this only using Ensembl associations, since these the MOD for GenMAPP
                    print 'exporting GenMAPP-GeneOntology tables for the gene ID system',mod,'for',species
                    exportAndFilterGenMAPPGOTables(species,go_to_mod,nested_go_mod,total_gene_count)
            except Exception:
                print mod, 'associated files not present!'
                null = 'null' ###No MOD relational files present
    return built_go_paths, path_goid_db, path_dictionary

def exportAndFilterGenMAPPGOTables(species,go_to_ensembl,nested_go_gene,total_gene_count):
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'

    import_dir = '/'+parent_dir+'OBO/builds'
    c = GrabFiles(); c.setdirectory(import_dir)
    go_tables = c.searchdirectory('go_table'); go_trees = c.searchdirectory('go_tree_table'); go_counts = c.searchdirectory('go_count_table')
    go_tables=filepath(go_tables[0]);go_trees=filepath(go_trees[0]);go_counts=filepath(go_counts[0]); x=0
    filterGenMAPPGOTables(go_tables,nested_go_gene,species)
    filterGenMAPPGOTables(go_trees,nested_go_gene,species)
    filterGenMAPPGOTables(go_counts,nested_go_gene,species)
    ###The fourth file to be output is the ensembl-GO count table built using nested and non-nested GO-gene counts
    program_type,database_dir = unique.whatProgramIsThis()
    export_data_fn = database_dir+'/'+species+'/nested/'+'genmapp_Ens-GOCounts.txt';
    export_data_fn=filepath(export_data_fn); export_data = open(export_data_fn,'w')
    title = ['GO','Count','Total']; title = string.join(title,'\t')+'\n'; export_data.write(title)
    values = string.join(['GO','0',str(total_gene_count)],'\t')+'\n'; export_data.write(values)
    for goid in nested_go_gene:
        ###Some terms in go_to_ensembl are not in nested_go_gene (less than 0.1% of terms). These are retired IDs that are in the Ensembl-GO table
        try: gene_count = str(len(go_to_ensembl[goid]))
        except KeyError: gene_count = '0'
        nested_gene_count = str(len(nested_go_gene[goid]))
        values = [goid,gene_count,nested_gene_count]
        values = string.join(values,'\t')+'\n'; export_data.write(values)
    export_data.close()

def filterGenMAPPGOTables(fn,nested_go_gene,species):
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
        
    export_data_fn = string.replace(fn,parent_dir+'OBO/builds',database_dir+'/'+species+'/nested/')
    print 'Writing species specific GenMAPP GO databases to',export_data_fn
    export_data = open(export_data_fn,'w'); x=0
    for line in open(fn,'r').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x == 0:
            y = 0; x=1; export_data.write(line)
            if 'OrderNo' in line:
                reassign_order = 'yes' ###The order is reassigned without breaks
            else: reassign_order = 'no'
            for i in t: ###Since the ID field varies in location, find it in the header
                if i == 'ID': id_index = y
                y+=1
        else:
            goid = t[id_index] ###If the GOID in the nested database, that term should be included (since nested account for all child terms linked to GO)
            if goid in nested_go_gene:
                if reassign_order == 'yes':
                    values = [str(x)]+t[1:]; line = string.join(values,'\t')+'\n'; x+=1
                export_data.write(line)
    export_data.close()

def speciesData():
    program_type,database_dir = unique.whatProgramIsThis()
    filename = 'Config/species.txt'
    fn=filepath(filename); global species_list; species_list=[]; global species_codes; species_codes={}
    for line in open(fn,'r').readlines():             
        data = cleanUpLine(line)
        abrev,species = string.split(data,'\t')
        species_list.append(species)
        species_codes[species] = abrev

def sourceData():
    program_type,database_dir = unique.whatProgramIsThis()
    filename = 'Config/source_data.txt'
    fn=filepath(filename)
    global source_types; source_types=[]
    global system_codes; system_codes={}
    global mod_types; mod_types=[]
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t'); source=t[0]
        try: system_code=t[1]
        except IndexError: system_code = 'NuLL'
        if len(t)>2: ### Therefore, this ID system is a potential MOD
            if t[2] == 'MOD': mod_types.append(source)
        if source not in mod_types: source_types.append(source) 
        system_codes[system_code] = source ###Used when users include system code data in their input file
                
if __name__ == '__main__':
    """This module imports GO hierarchy data, nests it, outputs it to GO-elite and associates
    gene level data with nested GO terms for MAPPFinder"""
    buildNestedGOTree('no','no');sys.exit()
    genmapp_mod = 'Ensembl'; species_code = 'Mm'; export_databases = 'no'; mod_types = ['Ensembl','EntrezGene']
    buildNestedGOAssociations(species_code,export_databases,mod_types,genmapp_mod)
    kill
    sourceData(); speciesData() ###finds species for which available information is stored locally   
    x = 1; print "Select species for GO table extraction"
    for species_name in species_list: print str(x)+')',species_name; x+=1
    inp = sys.stdin.readline(); inp = int(inp.strip())
    species = species_list[inp-1]; species_code = species_codes[species]
    export_databases = 'yes'

    x = 1; print "Specify MOD for GenMAPP-GeneOntology database export"
    for mod_name in mod_types: print str(x)+')',mod_name; x+=1
    print str(x)+') NA'
    inp = sys.stdin.readline(); inp = int(inp.strip())
    genmapp_mod = mod_types[inp-1]
    
    print 'Building GenMAPP tables for',species,species_code   
    #buildNestedGOTree('no') ###always no unless running via the buildNestedGOAssociations function
    mod_types = ['Ensembl','EntrezGene']
    buildNestedGOAssociations(species_code,export_databases,mod_types,genmapp_mod)
    print "Build Complete (hit any key to exit)"; sys.stdin.readline()
    
#!/usr/bin/python
###########################
#Program:	GO-elite.py
#Author:	Nathan Salomonis
#Date:		12/12/06
#Website:	http://www.genmapp.org
#Email:	nsalomonis@gmail.com
###########################


        
    
