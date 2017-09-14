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
prebuilt hiearchy paths, creating nested Ontology associations from existing gene-Ontology files."""

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
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
        try: files = getDirectoryFiles(self.data,str(search_term))
        except Exception:
            files = [] ### directory doesn't exist
            #print self.data, "doesn't exist"
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

class OntologyPath:
    def __init__(self,ontology_id,ontology_term,current_level,rank,path,specific_type):
        self._ontology_id = ontology_id; self._ontology_term = ontology_term; self._current_level = current_level
        self._rank = rank; self._path = path; self._specific_type = specific_type
    def OntologyID(self): return self._ontology_id
    def OntologyIDStr(self): return self._ontology_id[3:]
    def OntologyTerm(self): return self._ontology_term
    def OntologyLevel(self): return self._current_level
    def OntologyType(self): return self._specific_type
    def Rank(self): return self._rank
    def PathStr(self):
        path_index = pathToString(self.PathList())
        return path_index
    def PathList(self): return self._path
    def PathTuple(self): return tuple(self._path)
    def Report(self):
        output = self.OntologyID()+'|'+self.OntologyTerm()
        return output
    def __repr__(self): return self.Report()
    
def pathToString(path_list):
    path_str=[]
    for path_int in path_list: path_str.append(str(path_int))
    path_index = string.join(path_str,'.')
    return path_index

class OntologyTree:
    def __init__(self,ontology_id,ontology_term,ontology_type):
        self._ontology_id = ontology_id; self._ontology_term = ontology_term; self._ontology_type = ontology_type
    def OntologyID(self): return self._ontology_id
    def OntologyTerm(self): return self._ontology_term
    def OntologyType(self): return self._ontology_type
    def setOntologyType(self,ontology_type): self._ontology_type=ontology_type
    def Report(self):
        output = self.OntologyID()+'|'+self.OntologyTerm()
        return output
    def __repr__(self): return self.Report()

class OntologyTreeDetailed(OntologyTree):
    ###Class not currently used
    def __init__(self,ontology_id,ontology_term,ontology_type,parent_ontology_id,relation):
        self._ontology_id = ontology_id; self._ontology_term = ontology_term; self._ontology_type = ontology_type
        self._parent_ontology_id = parent_ontology_id; self._relation = relation
    def ParentOntologyID(self): return self._parent_ontology_id
    def Relation(self): return self._relation
    def Report(self):
        output = self.OntologyID()+'|'+self.OntologyTerm()
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

        path_ontology_db[tuple_path] = child
        if child not in built_ontology_paths:
            built_ontology_paths[child] = [tuple_path]
        elif tuple_path not in built_ontology_paths[child]:
            built_ontology_paths[child].append(tuple_path)
        path[-1]+=1
    return count_nodes
                                 
def importOBONew(filedir,path,specific_type,rank):
    if specific_type == '': discover_root = 'yes'
    else: discover_root = 'no'
    global edges
    #print [discover_root,specific_type,path]
    fn=filepath(filedir); x=0; stored={}; edges={}; category = 'null'; all_children={}; ontology_annotations_extra={}
    ontology_id=''; ontology_term=''; edge_count=0; root_node = None
    for line in open(fn,'r').xreadlines():             
        data = cleanUpLine(line)
        s = string.split(data,' '); d = string.split(data,':')
        if x == 0:
                x=1
        if x > 0:
            #if s[0]=='def:': definition = d[1]
            if s[0]=='id:':
                try:
                    ontology_id = s[1]
                    #ontology_id=string.split(ontology_id,':')[1]
                    category = 'null'
                except Exception: null=[]; ontology_id = ''; ontology_term = ''
            if s[0]=='namespace:': category = s[1]
            if s[0]=='name:':
                ontology_term = d[1][1:]
                if ontology_term == specific_type:
                    root_node = ontology_id
            if category == specific_type or discover_root=='yes':
                if s[0]=='is_a:': ### Note: sometimes there are multiple parents indicated for a single child
                    parent = s[1] ### immediate parent node
                    #parent=string.split(parent,':')[1]
                    if parent in edges: ### each child has one parent, one parent can have many children
                        children = edges[parent]
                        children[ontology_id]=[]
                    else: children = {ontology_id:[]}
                    edges[parent]=children
                    edge_count+=1
                    if discover_root=='yes': all_children[ontology_id] = []
                    if ontology_id not in ontology_annotations:
                        gt = OntologyTree(ontology_id,ontology_term,specific_type)  
                        ontology_annotations[ontology_id] = gt
                elif root_node == ontology_id: ### For example, biological process
                    gt = OntologyTree(ontology_id,ontology_term,specific_type)  
                    ontology_annotations[ontology_id] = gt
                elif ontology_id != '' and ontology_term != '':
                    gt = OntologyTree(ontology_id,ontology_term,specific_type)
                    ontology_annotations_extra[ontology_id] = gt
                    
    if discover_root=='yes':
        ### The root node should not exist as a child node
        for parent in edges:
            if parent not in all_children: root_node = parent
        specific_type = ontology_annotations_extra[root_node].OntologyTerm()
        #print 'Parent node assigned to:',specific_type
        ### Assing the root_node name as the Ontology-Type
        for ontology_id in ontology_annotations:
            ontology_annotations[ontology_id].setOntologyType(specific_type)
    if root_node == None:
        print 'NO ROOT NODE IDENTIFIED... SHOULD BE:', specific_type
        print filedir; kill
    
    if len(path)==0: path.append(0); path_ontology_db[tuple(path)] = root_node; return_path = list(path); #print [tuple(path)]
    else: path = [path[0]+1]; path_ontology_db[tuple(path)] = root_node; return_path = list(path); #print [tuple(path)]
    
    #export_data = export.ExportFile('OBO/test.txt')
    export_data=''
    nestTree(root_node,path,export_data,0)
    #export_data.close()
    #print 'Tree built'
    for path in path_ontology_db:
        path_dictionary[path]=[path]
        ###Build nested Path-index
        path_len = len(path); i=-1
        while path_len+i > 0:
            parent_path = path[:i]
            if parent_path in path_dictionary: path_dictionary[parent_path].append(path)
            i-=1
                    
    print edge_count,'edges and',len(ontology_annotations), 'Ontology annotations imported for',specific_type
    #print [[[return_path]]]
    return path_ontology_db,built_ontology_paths,ontology_annotations,path_dictionary,return_path,rank
###################################### UPDATED OBO CODE - END
        
def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def swapKeyValues(db):
    swapped={}
    for key in db:
        values = list(db[key]) ###If the value is not a list, make a list
        for value in values:
            try: swapped[value].append(key)
            except KeyError: swapped[value] = [key]
    swapped = eliminate_redundant_dict_values(swapped)
    return swapped

def exportCurrentOntologyBuild(path_ontology_db,ontology_annotations,ontology_type, display=False):
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    new_file = parent_dir+'OBO/builds/built_'+ontology_type+'_paths.txt'
    try: fn=filepath(new_file); data = open(fn,'w')
    except Exception:
        new_dir = parent_dir+'OBO/builds'; fn = filepath(new_dir)
        os.mkdir(fn) ###Re-Create directory if deleted
        fn=filepath(new_file); data = open(fn,'w')
    data.write('Path'+'\t'+'ontology_id'+'\n')
    for path in path_ontology_db:
        ontology_id = path_ontology_db[path]; path = pathToString(path)
        data.write(path +'\t'+ ontology_id +'\n')
    data.close()

    new_file = parent_dir+'OBO/builds/'+ontology_type+'_annotations.txt'
    fn=filepath(new_file); data = open(fn,'w')
    data.write('ontology_id'+'\t'+'Ontology Name'+'\t'+'Ontology Type'+'\n')    
    for ontology_id in ontology_annotations:
        s = ontology_annotations[ontology_id]
        data.write(ontology_id +'\t'+ s.OntologyTerm() +'\t'+ s.OntologyType() +'\n')
    data.close()

def convertStrListToIntList(path):
    path_int=[]
    for str in path: path_int.append(int(str))
    return path_int

def importPreviousOntologyAnnotations(target_ontology_type):
    ontology_annotations={}
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    if target_ontology_type == 'GeneOntology': target_ontology_type = 'go'
    filename = parent_dir+'OBO/builds/'+target_ontology_type+'_annotations.txt'; fn=filepath(filename); x=0
    for line in open(fn,'r').xreadlines():
        if x==0: x=1 ###Skip the title line
        else:
            data = cleanUpLine(line)
            ontology_id,ontology_name,ontology_type = string.split(data,'\t')
            if ':' not in ontology_id: ontology_id = 'GO:'+ontology_id
            if ontology_name[0]== ' ': ontology_name = ontology_name[1:]
            s = OntologyTree(ontology_id,ontology_name,ontology_type)
            ontology_annotations[ontology_id] = s
    return ontology_annotations

def importPreviousOntologyBuild(ontology_type,display=True):
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    if ontology_type == 'GeneOntology': ontology_type = 'go'
    filename = parent_dir+'OBO/builds/built_'+ontology_type+'_paths.txt'; fn=filepath(filename); x=0; count=0

    for line in open(fn,'r').xreadlines(): count+=1
    original_increment = int(count/10); increment = original_increment
    
    try: ### This reduces run-time for the typical analysis where the databases are in sync and up-to-date
        if run_mappfinder == 'yes':
            if verified_nested == 'no':
                build_nestedDB='yes'
            else: build_nestedDB = 'no'
        else: build_nestedDB = 'no'
    except Exception: build_nestedDB = 'yes'
            
    for line in open(fn,'r').xreadlines():
        if x==0: x+=1 ###Skip the title line
        else:
            x+=1
            if x == increment and display: increment+=original_increment; print '*',    
            data = cleanUpLine(line)
            path,ontology_id = string.split(data,'\t')
            path = tuple(map(int,string.split(path,'.')))
            #path = string.split(path_str,'.'); path = convertStrListToIntList(path); path = tuple(path)
            #s = OntologyPath(ontology_id,'','','',path,''); s = OntologyPathAbr(ontology_id,path)
            if ':' not in ontology_id: ontology_id = 'GO:'+ontology_id
            path_ontology_db[path] = ontology_id
            try: built_ontology_paths[ontology_id].append(path)
            except KeyError: built_ontology_paths[ontology_id] = [path]
            if build_nestedDB == 'yes':
                path_dictionary[path]=[path]
    ###All of the paths need to be added before  
    if build_nestedDB == 'yes':
        if build_nestedDB == 'yes':
            for path in path_dictionary:
                ###Build nested Path-index
                path_len = len(path); i=-1
                while path_len+i > 0:
                    parent_path = path[:i]
                    try: path_dictionary[parent_path].append(path)
                    except Exception: null=[]
                    i-=1    

#### Import gene data and associate with Nested Ontology
def grabNestedOntologyIDs():
    nested_ontology_tree={}
    for path in path_dictionary:
        parent_ontology_id = path_ontology_db[path]
        child_ontology_list=[]
        for child_path in path_dictionary[path]:
            child_ontology_id = path_ontology_db[child_path]; child_ontology_list.append(child_ontology_id)
        child_ontology_list = unique.unique(child_ontology_list)
        nested_ontology_tree[parent_ontology_id] = child_ontology_list
    return nested_ontology_tree
    
def linkGenesToNestedOntology(ontology_to_gene):
    nested_ontology_genes={}; made_unique={}; x=0
    original_increment = int(len(nested_ontology_tree)/10); increment = original_increment    

    for parent_ontology_id in nested_ontology_tree:
        x+=1
        if x == increment: increment+=original_increment; print '*',
        for child_ontology_id in nested_ontology_tree[parent_ontology_id]: ### This list of ontology_ids includes the parent, since it is the first entry in path_dictionary
            if child_ontology_id in ontology_to_gene:
                ensembls=ontology_to_gene[child_ontology_id]
                for ensembl in ensembls:
                    try:
                        ens_db = nested_ontology_genes[parent_ontology_id]
                        ens_db[ensembl] = ''
                    except KeyError:
                        ens_db = {}; ens_db[ensembl] = ''; e = ens_db
                        nested_ontology_genes[parent_ontology_id] = e 
    return nested_ontology_genes

def exportVersionData(version,version_date,dir):
    ### Used by the module UI
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
    
def exportOntologyRelationships(nested_ontology_gene,gene_to_source_id,mod,source_type,ontology_type):
    program_type,database_dir = unique.whatProgramIsThis()
    if ontology_type == 'GeneOntology': ontology_type = 'GO'
    new_file = database_dir+'/'+species_code+'/nested/'+mod+'_to_Nested-'+ontology_type+'.txt'
    data = export.ExportFile(new_file)
    title = [mod,'ontology_id']; title_str = string.join(title,'\t')
    data.write(title_str+'\n')
    for ontology_id in nested_ontology_gene:
        for gene in nested_ontology_gene[ontology_id]:
            output_list = [gene,ontology_id]
            output_str = string.join(output_list,'\t')
            data.write(output_str+'\n')
    data.close()
    print new_file, 'saved to disk'

#### Main functions that grab data from above functions
def remoteImportOntologyTree(ontology_type):
    global built_ontology_paths; global path_ontology_db; global path_dictionary
    built_ontology_paths={}; path_ontology_db={}; path_dictionary={}
    importPreviousOntologyBuild(ontology_type)
    return built_ontology_paths, path_ontology_db, path_dictionary
    
def buildNestedOntologyTree(mappfinder,display=True):
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    
    global run_mappfinder; run_mappfinder = mappfinder
    ###Import all the OBO Ontology tree information from http:/www.geneontology.org/
    import_dir = '/'+parent_dir+'OBO'; global Ontology_version; path=[]; rank=0
    c = GrabFiles(); c.setdirectory(import_dir)
    file_dirs = c.searchdirectory('.ontology')
    file_dirs += c.searchdirectory('.obo')
    file_dirs.reverse()
    x = file_dirs[1:]+file_dirs[0:1] ###Reorganize to mimic GenMAPP order
    start_time = time.time()
    ontology_type = ''
    #print file_dirs
    for file_dir in file_dirs:
        try:
            if '.obo' in file_dir or '.ontology' in file_dir:
                if 'gene_ontology' in file_dir or 'goslim' in file_dir:
                    ontology_type = 'GeneOntology'
                    if 'goslim' in file_dir: ontology_type = 'GOSlim'
                    ###Import the 3 main Ontology files and index them so that the first path corresponds to the Ontology type - Software checks the date before parsing
                    path_ontology_db,built_ontology_paths,ontology_annotations,path_dictionary,path,rank = importOBONew(file_dir,path,'biological_process',rank)
                    try: path_ontology_db,built_ontology_paths,ontology_annotations,path_dictionary,path,rank = importOBONew(file_dir,path,'molecular_function',rank)
                    except Exception: null=[] ### Sometimes missing from GO-Slim
                    path_ontology_db,built_ontology_paths,ontology_annotations,path_dictionary,path,rank = importOBONew(file_dir,path,'cellular_component',rank)
                else:
                    ontology_type = getOntologyType(file_dir)
                    path_ontology_db,built_ontology_paths,ontology_annotations,path_dictionary,path,rank = importOBONew(file_dir,path,'',rank)
                deleteNestedOntologyFiles(ontology_type) ### Necessary to trigger an update for all species
            else:
                if display: print 'The ontology format present in',file_dir,'is no longer supported.'
            exportCurrentOntologyBuild(path_ontology_db,ontology_annotations,ontology_type,display=display)
        except Exception:
            pass ### If an Ontology file fails download, it still may create an empty file that will screw up the processing of other obo files - just skip it
    end_time = time.time(); time_diff = int(end_time-start_time)
    
    if display: print "Ontology categories imported and nested in %d seconds" % time_diff
        
def getOntologyType(file_dir):
    ontology_type = string.split(file_dir,'/')[-1]
    if '_' in ontology_type:
        ontology_type = string.split(ontology_type,'_')[0]+'Ontology'
    else:
        ontology_type = string.split(ontology_type,'.')[0]+'Ontology'
    return ontology_type

def deleteNestedOntologyFiles(ontology_type):
    program_type,database_dir = unique.whatProgramIsThis()
    current_species_dirs = unique.read_directory('/'+database_dir)
    for species_code in current_species_dirs:
        c = GrabFiles(); c.setdirectory('/'+database_dir+'/'+species_code+'/nested')
        if ontology_type == 'GeneOntology': ontology_type = 'GO'
        file_dirs = c.searchdirectory('-'+ontology_type) ### list all nested files referencing the Ontology type
        for file in file_dirs:
            try: os.remove(filepath(database_dir+'/'+species_code+'/nested/'+file))
            except Exception: null=[]
    
def verifyFileLength(filename):
    count = 0
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            count+=1
            if count>3: break
    except Exception: null=[]
    return count

def verifyNestedFileCreation(species,mod_types,ontology_type):
    ### Determine which mods are present for Ontology
    program_type,database_dir = unique.whatProgramIsThis()
    mods_present = []; nested_present=[]; verified = 'no'
    for mod in mod_types:
        ontology_file = database_dir+'/'+species+'/gene-go/'+mod+'-'+ontology_type+'.txt'
        count = verifyFileLength(ontology_file) ### See if there are lines present in the file (if present)
        if count>1: mods_present.append(mod)
    if len(mods_present)>0:
        for mod in mods_present:
            if ontology_type == 'GeneOntology': ontology_type = 'GO'
            ontology_file = database_dir+'/'+species+'/nested/'+mod+'_to_Nested-'+ontology_type+'.txt'
            count = verifyFileLength(ontology_file) ### See if there are lines present in the file (if present)
            if count>1: nested_present.append(mod)
        if len(nested_present) == len(mods_present): verified = 'yes'
    return verified
        
def findAvailableOntologies(species,mod_types):
    program_type,database_dir = unique.whatProgramIsThis()
    c = GrabFiles(); c.setdirectory('/'+database_dir+'/'+species+'/gene-go'); file_dirs=[]
    for mod in mod_types:
        file_dirs+= c.searchdirectory(mod+'-')
    avaialble_ontologies=[]
    for filedir in file_dirs:
        ontology_type = string.split(filedir,'-')[-1][:-4] ### remove the .txt
        avaialble_ontologies.append(ontology_type)
    avaialble_ontologies = unique.unique(avaialble_ontologies)
    return avaialble_ontologies

def moveOntologyToArchiveDir(display=True):
    ### Move any existing OBO files to an archived directory as to not combine new with old annotations
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    c = GrabFiles()
    c.setdirectory('/'+parent_dir+'OBO')
    file_dirs = c.searchdirectory('.ontology')+c.searchdirectory('.obo')
    
    for file_dir in file_dirs:
        new_file_dir = string.replace(file_dir,parent_dir+'OBO/',parent_dir+'OBO/archive/')
        if display: print 'Moving:',file_dir,'to:',new_file_dir
        export.customFileMove(file_dir,new_file_dir)

    if len(file_dirs)==0:
        c.setdirectory('/'+'OBO')
        file_dirs = c.searchdirectory('.ontology')+c.searchdirectory('.obo')
        for file_dir in file_dirs:
            new_file_dir = string.replace(file_dir,'OBO/',parent_dir+'OBO/')
            if display: print 'Moving:',file_dir,'to:',new_file_dir
            export.customFileMove(file_dir,new_file_dir)
                
def buildNestedOntologyAssociations(species,mod_types,target_ontology_type,display=True):
    global species_code; species_code = species; global verified_nested
    global path_dictionary; path_dictionary={}
    global built_ontology_paths; built_ontology_paths={}
    global ontology_annotations; ontology_annotations={}
    global path_ontology_db; path_ontology_db={}
    
    if ('Linux' in platform.system()): mappfinder_db_input_dir = species_code+'/nested/'
    else: mappfinder_db_input_dir = '/'+species_code+'/nested/'
            
    buildNestedOntologyTree('yes') ### Checks the OBO directory to process new ontology files (if there)
    moveOntologyToArchiveDir(display=display) ### Move any new read ontology files to

    avaialble_ontologies = findAvailableOntologies(species,mod_types)

    verified_nested_db={}
    for ontology_type in avaialble_ontologies:
        ### This module verifies that the nested files are present (no longer considers database versions)
        verified_nested = verifyNestedFileCreation(species,mod_types,ontology_type) 
        verified_nested_db[ontology_type] = verified_nested
    verified_nested = verified_nested_db[target_ontology_type]
    importPreviousOntologyBuild(target_ontology_type,display=display) ### populates the global variables we return below
    if verified_nested == 'no':  ### modified this code such that any version change warrants a rebuild and if reset by BuildEntrezAffymetrixAssociations or other, that it triggers a rebuild
        if display: print 'Building %s Ontology nested gene association files for %s' % (target_ontology_type,species_code)
        ###Build Gene to Ontology associations for all MODs and export these for re-import by the MAPPFinder module
        global nested_ontology_tree
        nested_ontology_tree = grabNestedOntologyIDs()
        for mod in mod_types:
            try:
                start_time = time.time()
                mod_to_ontology = gene_associations.importGeneToOntologyData(species_code,mod,'null',target_ontology_type)
                ontology_to_mod = swapKeyValues(mod_to_ontology); total_gene_count = len(mod_to_ontology); mod_to_ontology=[]
                ###Obtain a database of ontology_ids with all nested gene associations
                nested_ontology_mod = linkGenesToNestedOntology(ontology_to_mod)
                exportOntologyRelationships(nested_ontology_mod,{},mod,'',target_ontology_type)
                end_time = time.time(); time_diff = int(end_time-start_time)
                if display: print "Ontology Nested Lists Process/Created in %d seconds" % time_diff
            except Exception:
                if mod != 'HMDB':
                    None ### optionally indicate if a MOD doesn't have local files supporting the creation of a nested set
                    #print mod, 'associated files not present!'
    return built_ontology_paths, path_ontology_db, path_dictionary

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
    """This module imports Ontology hierarchy data, nests it, outputs it to GO-Elite and associates
    gene level data with nested Ontology terms for MAPPFinder"""
    species_code = 'Hs'; mod_types = ['Ensembl','EntrezGene']; ontology_type = 'MPhenoOntology'
    buildNestedOntologyAssociations(species_code,mod_types,ontology_type); sys.exit()
    
#!/usr/bin/python
###########################
#Program:	GO-elite.py
#Author:	Nathan Salomonis
#Date:		12/12/06
#Website:	http://www.genmapp.org
#Email:	nsalomonis@gmail.com
###########################


        
    
