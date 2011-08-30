###gene_associations
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

"""This module contains fairly generic methods for reading different GO-Elite gene relationship
files, user input and denominator files, and summarizing quantitative and gene annotation data."""

import sys, string
import os.path
import unique
import math
import statistics
import gene_associations
import OBO_import
import export

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

###### Classes ######
class GrabFiles:
    def setdirectory(self,value):
        self.data = value
    def display(self):
        print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        file_dir,file = getDirectoryFiles(self.data,str(search_term))
        if len(file)<1: print search_term,'not found',self.data
        return file_dir,file
    
def getDirectoryFiles(import_dir, search_term):
    exact_file = ''; exact_file_dir=''
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        if (':' in import_dir) or ('/Users/' == import_dir[:7]): affy_data_dir = import_dir+'/'+data
        else: affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir:
            if 'version.txt' not in affy_data_dir: exact_file_dir = affy_data_dir; exact_file = data
    return exact_file_dir,exact_file

class GeneRelationships:
    def __init__(self,uid,gene,terms,uid_system,gene_system):
        self._uid = uid; self._gene = gene; self._terms = terms
        self._uid_system = uid_system; self._gene_system = gene_system
    def UniqueID(self): return self._uid
    def GeneID(self): return self._gene
    def GOIDs(self): return self._terms
    def GOIDInts(self):
        goid_list = self._terms; goid_list_int = []
        for goid in goid_list: goid_list_int.append(int(goid))
        return goid_list_int
    def GOIDStrings(self):
        goid_list = self._terms; goid_list_str = []
        for goid in goid_list: goid_list_str.append('GO:'+goid)
        return goid_list_str
    def MAPPs(self):
        mapp_list = self._terms
        return mapp_list
    def UIDSystem(self): return self.uid_system
    def GeneIDSystem(self): return self._gene_system
    def setUIDValues(self,uid_values):
        self._uid_values = uid_values
    def UIDValues(self): return self._uid_values
    def Report(self):
        output = self.UniqueID()
        return output
    def __repr__(self): return self.Report()

class GeneAnnotations:
    def __init__(self,gene,symbol,name,gene_system):
        self._gene = gene; self._symbol = symbol; self._name = name; 
        self._gene_system = gene_system
    def GeneID(self): return self._gene
    def Symbol(self): return self._symbol
    def Description(self): return self._name
    def GeneIDSystem(self): return self._gene_system
    def Report(self):
        output = self.GeneID()+'|'+self.Symbol()
        return output
    def __repr__(self): return self.Report()

class RedundantRelationships:
    def __init__(self,redundant_names,inverse_names):
        self._redundant_names = redundant_names; self._inverse_names = inverse_names
    def RedundantNames(self):
        redundant_with_terms = string.join(self._redundant_names,'|')
        return redundant_with_terms        
    def InverseNames(self):
        redundant_with_terms = string.join(self._inverse_names,'|')
        return redundant_with_terms
    def Report(self):
        output = self.RedundantNames()+'|'+self.InverseNames()
        return output
    def __repr__(self): return self.Report()
 
def cleanUpLine(line):
    data = string.replace(line,'\n','')
    data = string.replace(data,'\c','')
    data = string.replace(data,'\r','')
    data = string.replace(data,'"','')
    return data
    
###### Begin Gene-Relationship Parsing ######
def importGeneData(species_code,mod):
    ###Pass in species_code, mod
    program_type,database_dir = unique.whatProgramIsThis()
    gene_import_dir = '/'+database_dir+'/'+species_code+'/gene'
    g = GrabFiles(); g.setdirectory(gene_import_dir)
    filedir,filename = g.searchdirectory(mod) ###Identify gene files corresponding to a particular MOD
    
    system_type = string.replace(filename,'.txt','')
    fn=filepath(filedir); gene_annotations={}; x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x == 0: x = 1
        else:
            gene = t[0]; symbol = t[1]; name = t[2]
            s = GeneAnnotations(gene,symbol,name,system_type)
            gene_annotations[gene] = s
    return gene_annotations

def buildNestedAssociations(species_code):
    print "Nested gene association, not yet build for",species_code
    print "Proceeding to build nested associations for this species\n"
    export_databases = 'no'; genmapp_mod = 'Ensembl' ### default

    import GO_Elite; system_codes,source_types,mod_types = GO_Elite.getSourceData()
    OBO_import.buildNestedGOAssociations(species_code,export_databases,mod_types,genmapp_mod)
        
def importGeneGOData(species_code,mod,gotype):
    program_type,database_dir = unique.whatProgramIsThis()
    if gotype == 'nested': geneGO_import_dir = '/'+database_dir+'/'+species_code+'/nested'
    elif gotype == 'null': geneGO_import_dir = '/'+database_dir+'/'+species_code+'/gene-go'
    gg = GrabFiles(); gg.setdirectory(geneGO_import_dir)
    try: filedir,file = gg.searchdirectory(mod) ###Identify gene files corresponding to a particular MOD
    except Exception:
        if gotype == 'nested':
            buildNestedAssociations(species_code)
            filedir,file = gg.searchdirectory(mod)
    global gene_to_go; x=0
    
    fn=filepath(filedir); gene_to_go={}
    for line in open(fn,'rU').readlines():
        if gotype == 'nested' and x==0: x = 1
        else:
            if 'GeneOntology' not in line and 'GO ID' not in line: ### Header included in input from AP or BioMart
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if len(t)>1:
                    try: gene = t[0]; goid = t[1]
                    except IndexError: print [line],t, 'GO-gene relationship not valid...please clean up',filedir;sys.exit()
                    goid = string.replace(goid,'GO:','') ### Included in input from AP
                    ### If supplied without proceeding zeros
                    if len(goid)<7 and len(goid)>0:
                        extended_goid=goid
                        while len(extended_goid)< 7: extended_goid = '0'+extended_goid
                        goid = extended_goid
                    if len(gene)>0 and len(goid)>0:
                        try: gene_to_go[gene].append(goid)
                        except KeyError: gene_to_go[gene]= [goid]
    return gene_to_go

def importGeneMAPPData(species_code,mod):
    program_type,database_dir = unique.whatProgramIsThis()
    geneMAPP_import_dir = '/'+database_dir+'/'+species_code+'/gene-mapp'
    gm = GrabFiles(); gm.setdirectory(geneMAPP_import_dir)
    filedir,file = gm.searchdirectory(mod) ### Identify gene files corresponding to a particular MOD
    global gene_to_mapp; x = 0
    
    fn=filepath(filedir); gene_to_mapp={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            gene = t[0]; mapp = t[2]
            try: gene_to_mapp[gene].append(mapp)
            except KeyError: gene_to_mapp[gene]= [mapp]
    return gene_to_mapp

def grabFileRelationships(filename):
    filename = string.replace(filename,'.txt','')
    system1,system2 = string.split(filename,'-')
    return system1,system2

def importUidGeneSimple(species_code,mod_source):
    program_type,database_dir = unique.whatProgramIsThis()
    geneUID_import_dir = '/'+database_dir+'/'+species_code+'/uid-gene'
    ug = GrabFiles(); ug.setdirectory(geneUID_import_dir)
    filedir,file = ug.searchdirectory(mod_source) ### Identify gene files corresponding to a particular MOD
    gene_to_uid={}; x = 0
    
    fn=filepath(filedir); gene_to_mapp={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            gene = t[0]; uid = t[1]
            try: gene_to_uid[gene].append(uid)
            except KeyError: gene_to_uid[gene]= [uid]
    return gene_to_uid

def importUidGeneData(species_code,mod_source,gene_to_go,gene_to_mapp):
    program_type,database_dir = unique.whatProgramIsThis()
    geneUID_import_dir = '/'+database_dir+'/'+species_code+'/uid-gene'
    ug = GrabFiles(); ug.setdirectory(geneUID_import_dir)
    #print "begining to parse",mod_source, 'from',geneGO_import_dir
    filedir,filename = ug.searchdirectory(mod_source) ###Identify gene files corresponding to a particular MOD
    fn=filepath(filedir); uid_to_go={}; uid_to_mapp={}; count=0; x=0
    uid_system,gene_system = grabFileRelationships(filename)
    for line in open(fn,'r').xreadlines(): count+=1
    original_increment = int(count/20); increment = original_increment
    for line in open(fn,'rU').readlines():           
        data = cleanUpLine(line); x+=1
        if program_type == 'GO-Elite':
            if x == increment: increment+=original_increment; print '*',
        t = string.split(data,'\t')
        uid = t[1]; gene = t[0] 
        if gene in gene_to_go:
            goid = gene_to_go[gene]
            y = GeneRelationships(uid,gene,goid,uid_system,gene_system)
            try: uid_to_go[uid].append(y)
            except KeyError: uid_to_go[uid] = [y]
        if gene in gene_to_mapp:
            mapp_name = gene_to_mapp[gene]
            y = GeneRelationships(uid,gene,mapp_name,uid_system,gene_system)
            try: uid_to_mapp[uid].append(y)
            except KeyError: uid_to_mapp[uid] = [y]            
    return uid_to_go,uid_to_mapp,uid_system

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database: list = unique.unique(database[key]); list.sort(); db1[key] = list
    return db1

def getGeneToUid(species_code,mod_source):
    program_type,database_dir = unique.whatProgramIsThis()
    import_dir = '/'+database_dir+'/'+species_code+'/uid-gene'
    ug = GrabFiles(); ug.setdirectory(import_dir)
    filedir,filename = ug.searchdirectory(mod_source) ###Identify gene files corresponding to a particular MOD
    
    fn=filepath(filedir); gene_to_uid={}; count = 0; x=0
    uid_system,gene_system = grabFileRelationships(filename)

    for line in open(fn,'r').xreadlines(): count+=1
    original_increment = int(count/20); increment = original_increment
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line); x+=1
        if x == increment: increment+=original_increment; print '*',
        t = string.split(data,'\t')
        uid = t[1]; gene = t[0]
        try: gene_to_uid[gene].append(uid)
        except KeyError: gene_to_uid[gene] = [uid]
    gene_to_uid = eliminate_redundant_dict_values(gene_to_uid)
    return gene_to_uid

def getGeneToUidNoExon(species_code,mod_source):
    gene_to_uid = getGeneToUid(species_code,mod_source)
    try:
        probeset_db = simpleExonImporter(species_code); filtered_db={}
        for gene in gene_to_uid:
            for uid in gene_to_uid[gene]:
                try: probeset_db[uid]
                except KeyError: ### Only inlcude probesets not in the exon database
                    try: filtered_db[gene].append(uid)
                    except KeyError: filtered_db[gene] = [uid]
        return filtered_db
    except Exception: return gene_to_uid

def simpleExonImporter(species_code):
    filename = 'AltDatabase/'+species_code+'/exon/'+species_code+'_Ensembl_probesets.txt'
    fn=filepath(filename); probeset_db={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        data = string.split(data,'\t')
        probeset_db[data[0]]=[]
    return probeset_db

def predictIDSource(id,system_codes):
    au = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    nm = ['1','2','3','4','5','6','7','8','9']
    affy_suffix = '_at'; ensembl_prefix = 'ENS'; source_data = ''; id_type = ''
    if len(id)>3:
        if affy_suffix == id[-3:]: id_type = 'Affymetrix'
        elif ensembl_prefix == id[:3]: id_type = 'Ensembl'
        elif id[2] == '.': id_type = 'Unigene'
        elif (id[0] in au and id[1] in nm) or '_' in id: id_type = 'UniProt'
        else:
            try:
                int_val = int(id)
                if len(id) == 7: id_type = 'Affymetrix' ###All newer Affymetrix transcript_cluster_ids and probesets (can be mistaken for EntrezGene IDs)
                else: id_type = 'EntrezGene'
            except ValueError: null = []
    ###If the user changes the names of the above id_types, we need to verify that the id_type is in the database
    if len(id_type)>0:
        for code in system_codes:
            if system_codes[code] == id_type: source_data = id_type
    return source_data

def addNewCustomSystem(filedir,system,save_option,species_code):
    gene_annotations={}
    if 'update' in save_option:
        print "Importing existing",species_code,system,"relationships."
        gene_annotations = importGeneData(species_code,mod)
    fn=filepath(filedir); length3=0; lengthnot3=0
    for line in open(fn,'rU').readlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        data_length = len(t)
        ### Check the length of the input file
        if data_length==3:
            length3+=1
            if len(t[0])>0:
                gene = t[0]; symbol = t[1]; name = t[2]
                s = GeneAnnotations(gene,symbol,name,'')
                gene_annotations[gene] = s
        else: lengthnot3+=1
        
    if length3>lengthnot3:
        ###Ensure that the input file is almost completely 2 columns
        print 'Writing new annotation file:',species_code,system
        export_dir = 'Databases/'+species_code+'/gene/'+system+'.txt'
        data = export.ExportFile(export_dir)
        data.write('ID\tSymbol\tDescription\n')
        for gene in gene_annotations:
            s = gene_annotations[gene]
            data.write(string.join([gene,s.Symbol(),s.Description()],'\t')+'\n')
        data.close()
        return 'exported'
    else: return 'not-exported'
    
def addNewCustomRelationships(filedir,relationship_file,save_option,species_code):
    mod,data_type = string.split(relationship_file,'-'); mod_id_to_related={}
    if 'update' in save_option:
        print "Importing existing",species_code,relationship_file,"relationships."
        if data_type == 'MAPP':
            mod_id_to_related=importGeneMAPPData(species_code,mod)
        elif data_type == 'GeneOntology':
            mod_id_to_related=importGeneGOData(species_code,mod,'null')
            ###We have to reset the version to trigger a rebuild if we change any values in this table
            OBO_import.exportVersionData(0,'0/0/0','/'+species_code+'/nested/')
        else: mod_id_to_related=importUidGeneSimple(species_code,mod+'-'+data_type)
    fn=filepath(filedir); length2=0; lengthnot2=0
    for line in open(fn,'rU').readlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        data_length = len(t)
        ### Check the length of the input file
        if data_length==2:
            length2+=1
            if len(t[0])>0 and len(t[1])>0:
                try: mod_id_to_related[t[0]].append(t[1])
                except KeyError: mod_id_to_related[t[0]] = [t[1]]
        else: lengthnot2+=1
        
    if length2>lengthnot2:
        ###Ensure that the input file is almost completely 2 columns
        print 'Writing new relationship file:',species_code,relationship_file
        if data_type == 'MAPP': export_dir = 'Databases/'+species_code+'/gene-mapp/'+relationship_file+'.txt'
        elif data_type == 'GeneOntology': export_dir = 'Databases/'+species_code+'/gene-go/'+relationship_file+'.txt'
        else: export_dir = 'Databases/'+species_code+'/uid-gene/'+relationship_file+'.txt'
        data = export.ExportFile(export_dir)
        data.write('ID\tRelated ID\n')
        for mod_id in mod_id_to_related:
            for related_id in mod_id_to_related[mod_id]:
                if data_type == 'MAPP': ###Pathway files have 3 rather than 2 columns
                    data.write(string.join([mod_id,'',related_id],'\t')+'\n')
                else: data.write(string.join([mod_id,related_id],'\t')+'\n')
        data.close()
        if data_type == 'GeneOntology':
            OBO_import.exportVersionData(0,'0/0/0','/'+species_code+'/nested/')  ### Erase the existing file so that the database is re-remade
        return 'exported'
    else: return 'not-exported'
                
def importUIDsForMAPPFinderQuery(filedir,system_codes,return_uid_values):
    fn=filepath(filedir); x=0; uid_list = {}; uid_value_db = {}; source_data = ''; value_len_list = []; system_code_found = 'no'; first_uid=''; first_line = 'no lines in file'
    source_data_db={}; error=''
    for line in open(fn,'rU').readlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1; value_headers = t[1:]
        else:
            first_uid = t[0]; first_line = t
            x+=1; source_data_read = ''
            uid = t[0]
            uid = string.replace(uid,'---','') ###occurs for Affymetrix gene ID input
            uids = string.split(uid,' /// ') ###occurs for Affymetrix gene IDs where there are multiple associations
            for uid in uids:  
                if len(uid)>0: uid_list[uid]=[]
            if len(t)>1: ###Therefore there is either gene values or system data associated with imported primary IDs
                if t[1] in system_codes:
                    #if len(source_data)<1: source_data = system_codes[t[1]]; system_code_found = 'yes'
                    source_data = system_codes[t[1]]; system_code_found = 'yes'; source_data_read = source_data
                    if len(t)==2: values = [] ###Therefore, there is system code data but no associated values
                    else: values = t[2:]
                else: values = t[1:]
                value_len_list.append(len(values))
                if len(values)>0: uid_value_db[uid] = values
            
            if len(source_data_read)<1:
                for uid in uids:
                    source_data_read = predictIDSource(uid,system_codes)
                if len(source_data_read)>0: source_data = source_data_read
            if len(source_data_read)>0: source_data_db[source_data_read]=[]

    first_uid = string.replace(first_uid,'Worksheet','!!!!')
    filenames = string.split(filedir,'/'); filename = filenames[-1]
    if '!!!!' in first_uid:
        error = 'WARNING!!! There appears to be a formatting file issue with the file:\n"'+filename+'"\nPlease correct and re-run (should be tab-delimited text).'
    if len(source_data_db)>1:
        error = 'WARNING!!! There is more than one gene system (e.g., Ensembl and EntrezGene) in the file:\n"'+filename+'"\nPlease correct and re-run.'
    if source_data == '':
        error = 'WARNING!!! No gene SystemCode identified in:\n"'+filename+'"\nPlease provide in input text file and re-run.\n If SystemCode column present, the file format may be incorrect.'

    if return_uid_values == 'yes' and len(value_len_list)>0:
        value_len_list.sort(); longest_value_list = value_len_list[-1]
        for uid in uid_value_db:
            values = uid_value_db[uid]
            if len(values) != longest_value_list:
                while len(values) != longest_value_list: values.append('')
        if system_code_found == 'yes': value_headers = value_headers[1:]###Therfore, column #2 has system code data
        return uid_list,uid_value_db,value_headers,error
    elif return_uid_values == 'yes': return uid_list,uid_value_db,value_headers,error
    else: return uid_list,source_data,error

def grabNestedGeneAssociations(species_code,mod,source_data):
    program_type,database_dir = unique.whatProgramIsThis()
    global printout; printout = 'yes'
    gene_annotations = importGeneData(species_code,mod)
    try: gene_to_go = importGeneGOData(species_code,mod,'nested')
    except Exception: gene_to_go={}; print "NOTE: Skipping Gene-Ontology annotation inclusion for "+mod+"\n("+database_dir+'/'+species_code+') due to missing annotations.'
    try: gene_to_mapp = importGeneMAPPData(species_code,mod)
    except Exception: gene_to_mapp = {}
    if source_data != mod:
        mod_source = mod+'-'+source_data
        uid_to_go,uid_to_mapp,uid_system = importUidGeneData(species_code,mod_source,gene_to_go,gene_to_mapp)
        return uid_to_go, uid_to_mapp, uid_system, gene_annotations
    else:
        gene_to_go = convertGeneDBtoGeneObjects(gene_to_go,mod)
        gene_to_mapp = convertGeneDBtoGeneObjects(gene_to_mapp,mod)
    return gene_to_go, gene_to_mapp, mod, gene_annotations

def convertGeneDBtoGeneObjects(gene_to_pathway,mod):
    gene_to_pathway_objects={}
    for gene in gene_to_pathway:
        pathways = gene_to_pathway[gene]
        y = GeneRelationships(gene,gene,pathways,mod,mod) ###Note: in the analagous function for unique_ids, pathways is called pathway but is still a list (correct this at some point)
        try: gene_to_pathway_objects[gene].append(y)
        except KeyError: gene_to_pathway_objects[gene] = [y]
    return gene_to_pathway_objects
            
def matchInputIDsToGOEliteTerms(gene_file_dir,go_elite_output_dir,system_codes,mappfinder_file,nested_collapsed_go_tree,uid_to_go,gene_annotations,full_go_name_db,uid_system,prev_combined_associations,prev_combined_gene_ranking):
    global combined_associations ###This variable exists as a global under GO-Elite main as well, but needs to be re-initialized here as well
    global combined_gene_ranking; global value_headers
    combined_associations = prev_combined_associations
    combined_gene_ranking = prev_combined_gene_ranking
    try: uid_list,uid_value_db,value_headers,error = importUIDsForMAPPFinderQuery(gene_file_dir,system_codes,'yes')
    except Exception: uid_list,uid_value_db,value_headers,error = importUIDsForMAPPFinderQuery(gene_file_dir[1:],system_codes,'yes')

    filename_list = string.split(mappfinder_file,"/")
    filename = filename_list[-1]
    related_goid_db={}
    for uid in uid_list:
        if uid in uid_to_go:  ###Find probesets in input list with GO annotataions
            for y in uid_to_go[uid]: ### Different instances for each gene associated with the uid (e.g. multiple probeset to gene associations)
                related_goids = y.GOIDs()
                ###If there are numerical values associated with the primary unique ID, store these with the GeneRelationships instance data
                if uid in uid_value_db:
                    uid_values = uid_value_db[uid]
                    y.setUIDValues(uid_values)
                else: y.setUIDValues([])
                for goid in related_goids:
                    if goid in nested_collapsed_go_tree:
                        try: related_goid_db[goid].append(y)
                        except KeyError: related_goid_db[goid] = [y] ###store probesets linked to each goid
    ###Reduce related_goid_db to unique GO to gene relationships (rather than unique ids)
    for goid in related_goid_db:
        gene_db={}
        for y in related_goid_db[goid]:
            try: gene_db[y.GeneID()].append(y)
            except KeyError: gene_db[y.GeneID()]=[y]
        related_goid_db[goid] = gene_db
        
    go_elite_gene_associations={}; gene_uid_data={}
    for parent_goid in nested_collapsed_go_tree:
        go_associated = [] ###Add all genes belonging to the parent
        if parent_goid in related_goid_db:
            gene_db = related_goid_db[parent_goid]
            for gene in gene_db:  ###Too complicated to store this information in the next list: Create a database a call later
                gene_uid_data[gene] = gene_db[gene]
                go_associated.append(gene)
        go_associated = unique.unique(go_associated)
        go_elite_gene_associations[parent_goid] = go_associated
    return exportGeneGOAssociations(filename,go_elite_output_dir,go_elite_gene_associations,gene_uid_data,gene_annotations,full_go_name_db,uid_system)

def matchInputIDsToMAPPEliteTerms(gene_file_dir,go_elite_output_dir,system_codes,mappfinder_file,uid_to_mapp,zscore_mapp,gene_annotations,uid_system,prev_combined_associations,prev_combined_gene_ranking):
    global combined_associations ###This variable exists as a global under GO-Elite main as well, but needs to be re-initialized here as well
    global combined_gene_ranking; global value_headers
    combined_associations = prev_combined_associations
    combined_gene_ranking = prev_combined_gene_ranking
    try: uid_list,uid_value_db,value_headers,error = importUIDsForMAPPFinderQuery(gene_file_dir,system_codes,'yes')
    except Exception: uid_list,uid_value_db,value_headers,error = importUIDsForMAPPFinderQuery(gene_file_dir[1:],system_codes,'yes')

    filename_list = string.split(mappfinder_file,"/")
    filename = filename_list[-1]
    related_mapp_db={}
    for uid in uid_list:
        if uid in uid_to_mapp:  ###Find probesets in input list with MAPP annotataions
            for y in uid_to_mapp[uid]: 
                related_mapps = y.MAPPs()
                ###If there are numerical values associated with the primary unique ID, store these with the GeneRelationships instance data
                if uid in uid_value_db:
                    uid_values = uid_value_db[uid]
                    y.setUIDValues(uid_values)
                else: y.setUIDValues([])
                for mapp in related_mapps:
                    if mapp in zscore_mapp: ###Thus this MAPP is a GO-elite Local term
                        try: related_mapp_db[mapp].append(y)
                        except KeyError: related_mapp_db[mapp] = [y] ###store probesets linked to each MAPP
    ###Reduce related_mapp_db to unique GO to gene relationships (rather than unique ids)
    for mapp in related_mapp_db:
        gene_db={}
        for y in related_mapp_db[mapp]:
            try: gene_db[y.GeneID()].append(y)
            except KeyError: gene_db[y.GeneID()]=[y]
        related_mapp_db[mapp] = gene_db
        
    go_elite_gene_associations={}; gene_uid_data={}
    
    for mapp in zscore_mapp: ###Redundant with above, but preserving the structure of the analagous GO processing function
        mapp_associated = [] 
        if mapp in related_mapp_db:
            gene_db = related_mapp_db[mapp]
            for gene in gene_db:  ###Too complicated to store this information in the next list: Create a database a call later
                gene_uid_data[gene] = gene_db[gene]
                mapp_associated.append(gene)
        mapp_associated = unique.unique(mapp_associated)
        go_elite_gene_associations[mapp] = mapp_associated

    full_mapp_name_db={}
    return exportGeneGOAssociations(filename,go_elite_output_dir,go_elite_gene_associations,gene_uid_data,gene_annotations,full_mapp_name_db,uid_system)

def identifyRedundantPathways(go_elite_gene_associations):
    for goid_or_mapp1 in go_elite_gene_associations:
        for goid_or_mapp2 in go_elite_gene_associations:
            if goid_or_mapp1 != goid_or_mapp2:
                genes1 = go_elite_gene_associations[goid_or_mapp1]; total_genes1 = len(genes1)
                genes2 = go_elite_gene_associations[goid_or_mapp2]; total_genes1 = len(genes2)
                common_genes1 = []; common_genes2 = []
                for gene in genes1:
                    if gene in genes2: common_genes1.append(gene)
                for gene in genes2:
                    if gene in gene1: common_genes2.append(gene)
                ###Can get recipricol hits (unique_common_genes1 = unique_common_genes2), or hits in just one
    
def exportGeneGOAssociations(filename,go_elite_output_dir,go_elite_gene_associations,gene_uid_data,gene_annotations,full_go_name_db,uid_system):
    gene_annotation_filename = filename[0:-4]+ '_' +uid_system+"-gene-associations.txt"
    gene_ranking_filename = filename[0:-4]+ '_' +uid_system+"-Gene-Rankings.txt"
    
    original_gene_annotation_filename = gene_annotation_filename
    original_gene_ranking_filename = gene_ranking_filename
    if len(go_elite_output_dir) == 0:
        gene_annotation_filename = "output/gene_associations/" + gene_annotation_filename
        gene_ranking_filename = "output/gene_associations/" + gene_ranking_filename
        new_dir = 'output/gene_associations'
    else:
        gene_annotation_filename = go_elite_output_dir+'/gene_associations/'+gene_annotation_filename
        gene_ranking_filename = go_elite_output_dir+'/gene_associations/'+gene_ranking_filename
        new_dir = go_elite_output_dir
        
    go_gene_annotation_db={}; go_values_db={}
    data = export.ExportFile(gene_annotation_filename)
                 
    title = ['GeneID','Symbol','Description','UID-str','parent GOID','GO Name']+value_headers
    title = string.join(title,'\t')
    combined_associations[original_gene_annotation_filename]=['']  ###create a new key (filename) and append all entries to it
    combined_associations[original_gene_annotation_filename].append(title)
    data.write(title+'\n'); gene_centric_results={}; gene_to_goid={}
    for parent_goid in go_elite_gene_associations:
        try: go_name = full_go_name_db[parent_goid]
        except KeyError: go_name = parent_goid
        for gene in go_elite_gene_associations[parent_goid]:
            uid_data = gene_uid_data[gene]
            uid_list=[]; uid_value_db={}; indexes=[]; uid_value_list = []; uid_mean_value_list = []
            for y in uid_data:
                uid_list.append(y.UniqueID())
                index = 0
                for value in y.UIDValues():
                    try: uid_value_db[index].append(value)
                    except KeyError: uid_value_db[index] = [value]
                    indexes.append(index); index+=1
            indexes = unique.unique(indexes); indexes.sort()
            for index in indexes:
                value_str = string.join(uid_value_db[index],'|')
                try: value_mean = statistics.avg(uid_value_db[index])
                except ValueError: value_mean = 0
                uid_value_list.append(value_str); uid_mean_value_list.append(value_mean)
            uid_str = string.join(uid_list,'|')
            try:
                s = gene_annotations[gene]
                info = [s.GeneID(),s.Symbol(),s.Description(),uid_str,str(parent_goid),go_name]+uid_value_list
            except KeyError:
                s = GeneAnnotations(gene,'','','')
                info = [gene,'','',uid_str,str(parent_goid),go_name]
            ###Record these associations to include in the main GO-elite output results file
            try: go_gene_annotation_db[parent_goid].append(s)
            except KeyError: go_gene_annotation_db[parent_goid] = [s]
            try: go_values_db[parent_goid].append(uid_mean_value_list)
            except KeyError: go_values_db[parent_goid] = [uid_mean_value_list]
                
            try: gene_centric_results[gene,uid_str].append(go_name)
            except KeyError: gene_centric_results[gene,uid_str] = [go_name]
            try: gene_to_goid[gene].append(parent_goid)
            except KeyError: gene_to_goid[gene] = [parent_goid]
            try: info = string.join(info,'\t')
            except TypeError: print info;kill
            combined_associations[original_gene_annotation_filename].append(info)
            data.write(info+'\n')
    data.close()
    #print 'Nested gene associations (re-derived) written for GO-Elite results to:\n',gene_annotation_filename
    data = export.ExportFile(gene_ranking_filename); sorted_gene_ranking=[] 
    title = ['GeneID','Symbol','Description','#GO-Terms Associated','Percent GO-Terms Associated','UID-str','GO Names']
    title = string.join(title,'\t')
    combined_gene_ranking[original_gene_ranking_filename]=['']
    combined_gene_ranking[original_gene_ranking_filename].append(title)
    data.write(title+'\n')
    for (gene,uid_str) in gene_centric_results:
        go_terms = gene_centric_results[(gene,uid_str)]
        num_go_terms = len(go_terms); total_go_elite_terms = len(go_elite_gene_associations); percent = float(num_go_terms)/total_go_elite_terms
        go_name_str = string.join(go_terms,' // ')
        try:
            s = gene_annotations[gene]
            info = [s.GeneID(),s.Symbol(),s.Description(),str(num_go_terms),str(percent),uid_str,go_name_str]
        except KeyError:
            info = [gene,'','',str(num_go_terms),str(percent),uid_str,go_name_str]
        info = string.join(info,'\t')
        sorted_gene_ranking.append((num_go_terms,info))
    
    sorted_gene_ranking.sort(); sorted_gene_ranking.reverse()
    for (rank,info) in sorted_gene_ranking:
        combined_gene_ranking[original_gene_ranking_filename].append(info)
        data.write(info+'\n')
    data.close()
    ###Additionaly identify redundant GO-terms based solely on gene content
    goids_with_redundant_genes={}
    for parent_goid in go_elite_gene_associations:
        go_count_db={}
        for gene in go_elite_gene_associations[parent_goid]:
            associated_goids = gene_to_goid[gene]
            for goid in associated_goids:
                ###Count the number of genes associated with each goid associated with all genes (linked to a parent_goid)
                try: go_count_db[goid] += 1
                except KeyError: go_count_db[goid] = 1
                
        ##### MAJOR SOURCE OF KEYERRORS DUE TO MISMATCHING DATABASES OR ARRAY FILES - RAISE EXCEPTION IN GO-ELITE.PY
        del go_count_db[parent_goid] ###this goterm would otherwise be redundant

        for goid in go_count_db:
            overlap_gene_count = go_count_db[goid] ###number of overlaping genes between the parent_goid and current goid
            total_gene_count = len(go_elite_gene_associations[goid])
            if overlap_gene_count == total_gene_count:
                try: parent_go_name = full_go_name_db[parent_goid]
                except KeyError: parent_go_name = parent_goid
                try: go_name = full_go_name_db[goid]
                except KeyError: go_name = goid
                try: goids_with_redundant_genes[goid].append((parent_go_name,parent_goid))
                except KeyError: goids_with_redundant_genes[goid] = [(parent_go_name,parent_goid)]
    inverse_goids_with_redundant_genes={}
    ###Now, select the inverse of redundant GOIDs (report parent redundancies for each redundant term
    all_redundant={}
    for child_goid in goids_with_redundant_genes:
        try: child_go_name = full_go_name_db[child_goid]
        except KeyError: child_go_name = child_goid
        parent_go_names = []
        for (parent_go_name,parent_goid) in goids_with_redundant_genes[child_goid]:
            parent_go_names.append(parent_go_name) ###we just want to report the name, but we needed the goid for this function
            try: inverse_goids_with_redundant_genes[parent_goid].append(child_go_name)
            except KeyError: inverse_goids_with_redundant_genes[parent_goid] = [child_go_name]
            all_redundant[parent_goid] = []
        goids_with_redundant_genes[child_goid] = parent_go_names ###replace the existing entry with a list of names only
        all_redundant[child_goid] = []
    goids_with_redundant_genes = eliminate_redundant_dict_values(goids_with_redundant_genes)
    inverse_goids_with_redundant_genes = eliminate_redundant_dict_values(inverse_goids_with_redundant_genes)
    goids_with_redundant_genes2={} ###Combine the redundant and inverse data and store as instances of RedundantRelationships
    for goid in all_redundant:
        redundant_go_names=[' ']; inverse_go_names=[' ']
        if goid in goids_with_redundant_genes: redundant_go_names = goids_with_redundant_genes[goid]
        if goid in inverse_goids_with_redundant_genes: inverse_go_names = inverse_goids_with_redundant_genes[goid]
        rr = RedundantRelationships(redundant_go_names,inverse_go_names)
        goids_with_redundant_genes2[goid] = rr
    goids_with_redundant_genes = goids_with_redundant_genes2
    #for go_name in goids_with_redundant_genes: print go_name,goids_with_redundant_genes[go_name]
    ###For each column of values, summarize the numerical values (take the mean) for each GO/MAPP term
    for goid in go_values_db:
        index_db={}; index_ls=[]
        for vals in go_values_db[goid]:
            index = 0
            for val in vals:
                if val != 0: ###only occurs when the value was non-numeric, otherwise it's a float
                    try: index_db[index].append(val)
                    except KeyError: index_db[index] = [val]
                    index_ls.append(index) ###for sorting through
                index+=1
        index_ls = unique.unique(index_ls); index_ls.sort()
        summary_values = []; summary_stdev_values = []
        for index in index_ls:
            try:
                try: avg_val = statistics.avg(index_db[index]); summary_values.append(str(avg_val))
                except KeyError: summary_values.append('')
                try: stdev_val = statistics.stdev(index_db[index]); summary_stdev_values.append(str(stdev_val))
                except KeyError: summary_stdev_values.append('')
            except ValueError: summary_values.append(''); summary_stdev_values.append('')
        go_values_db[goid] = summary_values, summary_stdev_values
    #print 'Gene counts (re-derived) for GO-Elite results writen to:\n',gene_ranking_filename
    return combined_associations,combined_gene_ranking,go_gene_annotation_db,go_values_db,value_headers,goids_with_redundant_genes,len(gene_to_goid)

def exportCombinedAssociations(combined_associations,go_elite_output_folder,file_suffix):
    new_filename = 'pruned-'+file_suffix+'.txt'
    if len(go_elite_output_folder)==0: output_results = 'output/gene_associations/' + new_filename
    else: output_results = go_elite_output_folder +'/gene_associations/'+ new_filename
    
    data = export.ExportFile(output_results)
    for file in combined_associations:
        for line in combined_associations[file]:
            info = file +'\t'+ line +'\n'
            data.write(info)
    data.close()
    #print 'Combined gene associations (re-derived) written for GO-Elite results to:\n',output_results
                              
if __name__ == '__main__':
    species_code = 'Bt'; mod = 'Ensembl'; gotype='nested'
    filedir = 'C:/Documents and Settings/Nathan Salomonis/My Documents/GO-Elite_120beta/Databases/EnsMart56Plus/Ce/gene/EntrezGene.txt'
    system = 'Macaroni'
    addNewCustomSystem(filedir,system,'yes','Ms'); kill
    addNewCustomRelationships('test.txt','Ensembl-MAPP','update','Mm');kill
    importGeneGOData(species_code,mod,gotype);kill
    
    species_name = 'Homo Sapiens'; species_code = 'Hs'; source_data = 'EntrezGene'; mod = 'EntrezGene'
    system_codes={}; system_codes['X'] = 'Affymetrix'
    import_dir = '/input/GenesToQuery/'+species_code
    m = GrabFiles(); m.setdirectory(import_dir)
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for mappfinder_input in dir_list:    #loop through each file in the directory
        permuted_z_scores={}; original_go_z_score_data={}; original_mapp_z_score_data={}
        gene_file_dir, gene_file = m.searchdirectory(mappfinder_input)
        ###Import Input gene/source-id lists
        input_gene_list,source_data_input,error = importUIDsForMAPPFinderQuery('input/GenesToQuery/'+species_code+'/'+gene_file,system_codes,'no'); input_count = len(input_gene_list)
        
    uid_to_go, uid_system, gene_annotations = grabNestedGeneAssociations(species_code,mod)

#!/usr/bin/python
###########################
#Program:	GO-elite.py
#Author:	Nathan Salomonis
#Date:		12/12/06
#Website:	http://www.genmapp.org
#Email:	nsalomonis@gmail.com
###########################


        
    
