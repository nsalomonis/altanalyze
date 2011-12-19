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
import re

###### File Import Functions ######
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv" or ".owl" in entry or ".gpml" in entry or ".xml" in entry or ".gmt" in entry: dir_list2.append(entry)
    return dir_list2

###### Classes ######
class GrabFiles:
    def setdirectory(self,value):
        self.data = value
    def display(self):
        print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        all_matching,file_dir,file = getDirectoryFiles(self.data,str(search_term))
        #if len(file)<1: print self.data,search_term,'not found', filepath(self.data)
        return file_dir,file
    def getAllFiles(self,search_term):
        #self is an instance while self.data is the value of the instance
        all_matching,file_dir,file = getDirectoryFiles(self.data,str(search_term))
        #if len(file)<1: print search_term,'not found'
        return all_matching
    
def getDirectoryFiles(import_dir, search_term):
    exact_file = ''; exact_file_dir=''; all_matching=[]
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        if (':' in import_dir) or ('/Users/' == import_dir[:7]): affy_data_dir = import_dir+'/'+data
        else: affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir:
            if 'version.txt' not in affy_data_dir: exact_file_dir = affy_data_dir; exact_file = data; all_matching.append(exact_file_dir)
    return all_matching, exact_file_dir,exact_file

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
        for goid in goid_list:
            if ':' not in goid: goid_list_str.append('GO:'+goid)
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
    def SymbolLower(self): return string.lower(self._symbol)
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
    
def importGeneric(filename):
    fn=filepath(filename); key_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        key_db[t[0]] = t[1:]
    return key_db

###### Begin Gene-Relationship Parsing ######
def importGeneData(species_code,mod):
    if 'export' in mod: export_status,mod = mod
    else: export_status = 'no' 
    ###Pass in species_code, mod
    program_type,database_dir = unique.whatProgramIsThis()
    gene_import_dir = '/'+database_dir+'/'+species_code+'/gene'
    g = GrabFiles(); g.setdirectory(gene_import_dir)
    filedir,filename = g.searchdirectory(mod) ###Identify gene files corresponding to a particular MOD
    
    system_type = string.replace(filename,'.txt','')
    fn=filepath(filedir); gene_annotations={}; x=0
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x == 0: x = 1
        else:
            gene = t[0]; symbol = t[1]; name = t[2]
            s = GeneAnnotations(gene,symbol,name,system_type)
            gene_annotations[gene] = s
    if export_status == 'export':
        export_dir = 'Databases/'+species_code+'/uid-gene/'+mod+'-Symbol'+'.txt'
        data = export.ExportFile(export_dir)        
        for gene in gene_annotations:
            s = gene_annotations[gene]
            if len(s.Symbol())>0:
                data.write(gene+'\t'+s.Symbol()+'\n')
        data.close()
    else:
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
    for line in open(fn,'rU').xreadlines():
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
                        if ':' not in goid:
                            goid = 'GO:'+goid
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
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            gene = t[0]; mapp = t[2]
            try: gene_to_mapp[gene].append(mapp)
            except KeyError: gene_to_mapp[gene]= [mapp]
    return gene_to_mapp

def exportCustomPathwayMappings(gene_to_custom,mod,system_codes,custom_sets_folder):
    export_dir = custom_sets_folder+'/CustomGeneSets/custom_gene_set.txt'
    print 'Exporting:',export_dir
    try: data = export.ExportFile(export_dir)
    except Exception: data = export.ExportFile(export_dir[1:])
    for system_code in system_codes:
        if system_codes[system_code] == mod: mod_code = system_code
    
    data.write('GeneID\t'+mod+'\tPathway\n'); relationships=0
    for gene in gene_to_custom:
        for pathway in gene_to_custom[gene]:
            values = string.join([gene,mod_code,pathway],'\t')+'\n'
            data.write(values); relationships+=1
    data.close()
    print relationships,'Custom pathway-to-ID relationships exported...'

def exportNodeInteractions(pathway_db,mod,system_codes,custom_sets_folder):
    export_dir = custom_sets_folder+'/Interactomes/interactions.txt'
    print 'Exporting:',export_dir
    try: data = export.ExportFile(export_dir)
    except Exception: data = export.ExportFile(export_dir[1:])
    for system_code in system_codes:
        if system_codes[system_code] == mod: mod_code = system_code
    
    data.write('Symbol1\tInteractionType\tSymbol2\t'+mod+'\tGeneID1\tGeneID2\tPathway\n'); relationships=0
    for pathway_id in pathway_db:
        wpd = pathway_db[pathway_id]
        
        for itd in wpd.Interactions():
            gi1 = itd.GeneObject1()
            gi2 = itd.GeneObject2()
            try:
                values = string.join([gi1.GeneName(),itd.InteractionType(),gi2.GeneName(),gi1.MODID(),gi2.MODID(),wpd.Pathway()],'\t')+'\n'
                data.write(values); relationships+=1
            except AttributeError: null=[] ### Occurs if a MODID is not present for one of the nodes
            
    data.close()
    print relationships,'Interactions exported...'   

def importGeneCustomData(species_code,system_codes,custom_sets_folder,mod):
    print 'Importing custom pathway relationships...'
    #print 'Trying to import text data'
    gene_to_custom = importTextCustomData(species_code,system_codes,custom_sets_folder,mod)
    #print 'Trying to import gmt data'
    gmt_data = parseGMT(custom_sets_folder)
    #print 'Trying to import gpml data'
    gpml_data,pathway_db = parseGPML(custom_sets_folder)
    #print 'Trying to import biopax data'
    biopax_data = parseBioPax(custom_sets_folder)
    #print 'Unifying gene systems for biopax'
    gene_to_BioPax = unifyGeneSystems(biopax_data,species_code,mod); #print len(gene_to_BioPax)
    #print 'Unifying gene systems for WP'
    gene_to_WP = unifyGeneSystems(gpml_data,species_code,mod); #print len(gene_to_WP)
    #print 'Unifying gene systems for gmt'
    gene_to_GMT = unifyGeneSystems(gmt_data,species_code,mod); #print len(gene_to_WP)
    #print 'Combine WP-biopax'
    gene_to_xml = combineDBs(gene_to_WP,gene_to_BioPax)
    #print 'Combine xml-text'
    gene_to_custom = combineDBs(gene_to_xml,gene_to_custom)
    #print 'Combine gmt-other'
    gene_to_custom = combineDBs(gene_to_GMT,gene_to_custom)
    exportCustomPathwayMappings(gene_to_custom,mod,system_codes,custom_sets_folder)
    if len(gene_to_WP)>0:
        try: exportNodeInteractions(pathway_db,mod,system_codes,custom_sets_folder)
        except Exception: null=[]
    
    ### Combine WikiPathway associations with the custom
    try: gene_to_mapp = importGeneMAPPData(species_code,mod)
    except Exception: gene_to_mapp = {}
    for gene in gene_to_mapp:
        for mapp in gene_to_mapp[gene]:
            try: gene_to_custom[gene].append(mapp)
            except KeyError: gene_to_custom[gene]= [mapp]

    return gene_to_custom

def importTextCustomData(species_code,system_codes,custom_sets_folder,mod):
    program_type,database_dir = unique.whatProgramIsThis()
    gm = GrabFiles(); gm.setdirectory(custom_sets_folder); system = None
    filedirs = gm.getAllFiles('.txt')
    global gene_to_custom; gene_to_custom={}
    for filedir in filedirs:
        file = string.split(filedir,'/')[-1]
        print "Reading custom gene set",filedir
        fn=filepath(filedir); gene_to_custom={}; x = 0
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            if x==0: x=1
            else:
                t = string.split(data,'\t')
                try: gene = t[0]; mapp = t[2]; system_code = t[1]
                except Exception:
                    print file, 'is not propperly formatted (skipping import of relationships)'; break
                if system_code in system_codes: system = system_codes[system_code]
                else: print system_code, "is not a recognized system code. Skipping import of",file; break
                try: gene_to_custom[gene].append(mapp)
                except KeyError: gene_to_custom[gene]= [mapp]
                
        ### If the system code is not the MOD - Convert to the MOD
        if (system != mod) and (system != None):
            gene_to_custom2={}
            mod_source = mod+'-'+system
            try: gene_to_source_id = getGeneToUid(species_code,mod_source)
            except Exception: print mod_source,'relationships not found. Skipping import of',file; break
            source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
            for source_id in gene_to_custom:
                if source_id in source_to_gene:
                    for gene in source_to_gene[source_id]:
                        gene_to_custom2[gene] = gene_to_custom[source_id]
            gene_to_custom = gene_to_custom2

    return gene_to_custom

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
    
    fn=filepath(filedir)
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            gene = t[0]; uid = t[1]
            try: gene_to_uid[gene].append(uid)
            except KeyError: gene_to_uid[gene]= [uid]
    return gene_to_uid

def augmentEnsemblGO(species_code):
    entrez_ens = importUidGeneSimple(species_code,'EntrezGene-Ensembl')
    try: ens_go=importGeneGOData(species_code,'Ensembl','null')
    except Exception: ens_go = {}
    try: entrez_go=importGeneGOData(species_code,'Entrez','null')
    except Exception: entrez_go = {}
    ens_go_translated = {}
    for entrez in entrez_go:
        if entrez in entrez_ens:
            if len(entrez_ens[entrez])<3: ### Limit bad associations
                #print entrez,entrez_ens[entrez];kill
                for ens in entrez_ens[entrez]: ens_go_translated[ens] = entrez_go[entrez]
    
    ### Add these relationships to the original
    for ens in ens_go_translated:
        try: ens_go[ens] = unique.unique(ens_go[ens] + ens_go_translated[ens])
        except Exception: ens_go[ens] = ens_go_translated[ens]

    program_type,database_dir = unique.whatProgramIsThis()
    export_dir = database_dir+'/'+species_code+'/gene-go/Ensembl-GeneOntology.txt'
    
    print 'Exporting augmented:',export_dir
    try: data = export.ExportFile(export_dir)
    except Exception: data = export.ExportFile(export_dir[1:])
    
    data.write('GeneID'+'\tGO ID\n')
    for gene in ens_go:
        for goid in ens_go[gene]:
            values = string.join([gene,goid],'\t')+'\n'
            data.write(values)
    data.close()

    ### Build Nested
    export_databases = 'no'; genmapp_mod = 'Ensembl'
    full_path_db,path_id_to_goid,null = OBO_import.buildNestedGOAssociations(species_code,export_databases,['Ensembl'],genmapp_mod)

    
def importUidGeneData(species_code,mod_source,gene_to_go,gene_to_mapp,denominator_source_ids):
    program_type,database_dir = unique.whatProgramIsThis()
    geneUID_import_dir = '/'+database_dir+'/'+species_code+'/uid-gene'
    ug = GrabFiles(); ug.setdirectory(geneUID_import_dir)
    #print "begining to parse",mod_source, 'from',geneGO_import_dir
    filedir,filename = ug.searchdirectory(mod_source) ###Identify gene files corresponding to a particular MOD
    fn=filepath(filedir); uid_to_go={}; uid_to_mapp={}; count=0; x=0
    uid_system,gene_system = grabFileRelationships(filename)
    for line in open(fn,'rU').xreadlines(): count+=1
    original_increment = int(count/20); increment = original_increment
    for line in open(fn,'rU').xreadlines():           
        data = cleanUpLine(line); x+=1
        if program_type == 'GO-Elite':
            if x == increment: increment+=original_increment; print '*',
        t = string.split(data,'\t')
        uid = t[1]; gene = t[0] 
        try:
            if len(denominator_source_ids)>0:
                null=denominator_source_ids[uid] ### requires that the source ID be in the list of analyzed denominators
            goid = gene_to_go[gene]
            y = GeneRelationships(uid,gene,goid,uid_system,gene_system)
            try: uid_to_go[uid].append(y)
            except KeyError: uid_to_go[uid] = [y]
        except Exception: null=[]
        try:
            if len(denominator_source_ids)>0:
                null=denominator_source_ids[uid] ### requires that the source ID be in the list of analyzed denominators
            mapp_name = gene_to_mapp[gene]
            y = GeneRelationships(uid,gene,mapp_name,uid_system,gene_system)
            try: uid_to_mapp[uid].append(y)
            except KeyError: uid_to_mapp[uid] = [y]
        except Exception: null=[]
    return uid_to_go,uid_to_mapp,uid_system

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database: list = unique.unique(database[key]); list.sort(); db1[key] = list
    return db1

def getGeneToUid(species_code,mod_source):
    if 'hide' in mod_source: show_progress, mod_source = mod_source
    else: show_progress = 'yes'
    program_type,database_dir = unique.whatProgramIsThis()
    import_dir = '/'+database_dir+'/'+species_code+'/uid-gene'
    ug = GrabFiles(); ug.setdirectory(import_dir)
    filedir,filename = ug.searchdirectory(mod_source) ###Identify gene files corresponding to a particular MOD
    
    fn=filepath(filedir); gene_to_uid={}; count = 0; x=0
    uid_system,gene_system = grabFileRelationships(filename)

    for line in open(fn,'r').xreadlines(): count+=1
    original_increment = int(count/20); increment = original_increment
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line); x+=1
        if x == increment and show_progress == 'yes': increment+=original_increment; print '*',
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
    for line in open(fn,'rU').xreadlines():             
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
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        data_length = len(t)
        ### Check the length of the input file
        if data_length==3:
            length3+=1
            if len(t[0])>0:
                try: gene = t[0]; symbol = t[1]; name = t[2]
                except Exception:
                    print 'Unexpected error in the following line'; print t, filedir;kill
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
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        data_length = len(t)
        ### Check the length of the input file
        if data_length==2:
            length2+=1
            if len(t[0])>0 and len(t[1])>0:
                if ',' in t[0]: keys = string.split(t[0],',') ### These can be present in the MGI phenotype ontology gene association files
                else: keys = [t[0]]
                for key in keys:
                    try: mod_id_to_related[key].append(t[1])
                    except KeyError: mod_id_to_related[key] = [t[1]]
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
    source_data_db={}; error=''; bad_alphanumerics=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        tnull = string.join(t,''); ### Indicates whether the line is blank - even if tabs are present
        alphanumeric = string.join(re.findall(r"\w",data))
        if x==0:
            x=1; value_headers = t[1:]
        else:
            if len(tnull)>0:
                if len(alphanumeric)==0: bad_alphanumerics+=1
                first_uid = t[0]; first_line = t
                x+=1; source_data_read = ''
                uid = t[0]
                uid = string.replace(uid,'---','') ###occurs for Affymetrix gene ID input
                uids = string.split(uid,' /// ') ###occurs for Affymetrix gene IDs where there are multiple associations
                for uid in uids:  
                    if len(uid) > 0 and uid != ' ':
                        if uid[0]==' ': uid = uid[1:]
                        if uid[-1]==' ': uid = uid[:-1]
                        all_alphanumeric = re.findall(r"\w",uid)
                        if len(all_alphanumeric)>0: uid_list[uid]=[]
                if len(t)>1: ###Therefore there is either gene values or system data associated with imported primary IDs
                    try:
                        system = t[1]
                        if system[0]==' ': system = system[1:]
                        if system[-1]==' ': system = system[:-1]
                        if system in system_codes:
                            #if len(source_data)<1: source_data = system_codes[t[1]]; system_code_found = 'yes'
                            source_data = system_codes[system]; system_code_found = 'yes'; source_data_read = source_data
                            if len(t)==2: values = [] ###Therefore, there is system code data but no associated values
                            else: values = t[2:]
                        else: values = t[1:]
                        value_len_list.append(len(values))
                        if len(values)>0:
                            if len(uid) > 0 and uid != ' ':
                                all_alphanumeric = re.findall(r"\w",uid)
                                if len(all_alphanumeric)>0: uid_value_db[uid] = values
                    except Exception:
                        source_data_read =''
                if len(source_data_read)<1:
                    for uid in uids:
                        source_data_read = predictIDSource(uid,system_codes)
                    if len(source_data_read)>0: source_data = source_data_read
                if len(source_data_read)>0: source_data_db[source_data_read]=[]
    first_uid = string.replace(first_uid,'Worksheet','!!!!')
    filenames = string.split(filedir,'/'); filename = filenames[-1]
    if '!!!!' in first_uid:
        error = 'WARNING!!! There appears to be a formatting file issue with the file:\n"'+filename+'"\nPlease correct and re-run (should be tab-delimited text).'
    elif len(source_data_db)>1:
        error = 'WARNING!!! There is more than one gene system (e.g., Ensembl and EntrezGene) in the file:\n"'+filename+'"\nPlease correct and re-run.'
    elif source_data == '':
        error = 'WARNING!!! No gene SystemCode identified in:\n"'+filename+'"\nPlease provide in input text file and re-run.\n If SystemCode column present, the file format may be incorrect.'
    elif x>0 and bad_alphanumerics>1:
        error = 'WARNING!!! Invalid text file encoding found (invalid alphanumeric values). Please resave text file in a standard tab-delimited format'
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
    
def getAllDenominators(denom_search_dir,system_codes):
    import mappfinder; denom_all_db={}
    m = GrabFiles(); m.setdirectory(denom_search_dir)
    if len(denom_search_dir)>0: denom_dir_list = mappfinder.readDirText(denom_search_dir)
    else: denom_dir_list = []
    for denominator_file in denom_dir_list:
        gene_file_dir, gene_file = m.searchdirectory(denominator_file)
        denom_gene_list,source_data_denom,error_message = importUIDsForMAPPFinderQuery(denom_search_dir+'/'+gene_file,system_codes,'no')
        for i in denom_gene_list: denom_all_db[i]=[]
    return denom_all_db
    
def grabNestedGeneAssociations(species_code,mod,source_data,system_codes,custom_sets_folder,denom_search_dir):
    program_type,database_dir = unique.whatProgramIsThis()
    global printout; printout = 'yes'
    gene_annotations = importGeneData(species_code,mod)
    ### Filter source IDs for those in all user denominator files
    denominator_source_ids = getAllDenominators(denom_search_dir,system_codes)
    
    try: gene_to_go = importGeneGOData(species_code,mod,'nested')
    except Exception: print "Warning...the MOD you have selected:",mod,"is missing the appropriate relationship files, necessary to run GO-Elite.  Either replace the missing files ("+database_dir+'/'+species_code+') or select a different MOD at runtime. Exiting program.'; sys.stdin.readline(); sys.exit()
    if len(custom_sets_folder)>0:
        gene_to_mapp = importGeneCustomData(species_code,system_codes,custom_sets_folder,mod)
    else:
        try: gene_to_mapp = importGeneMAPPData(species_code,mod)
        except Exception: gene_to_mapp = {}
    if source_data != mod:
        mod_source = mod+'-'+source_data
        uid_to_go,uid_to_mapp,uid_system = importUidGeneData(species_code,mod_source,gene_to_go,gene_to_mapp,denominator_source_ids)
        gene_to_go = convertGeneDBtoGeneObjects(gene_to_go,mod)
        gene_to_mapp = convertGeneDBtoGeneObjects(gene_to_mapp,mod)
        ### Combine the gene and UID data in order to later filter out gene's linked to an arrayid that are not linked to the pathway (common for some pathways)
        for gene in gene_to_go: uid_to_go[gene] = gene_to_go[gene]
        for gene in gene_to_mapp: uid_to_mapp[gene] = gene_to_mapp[gene]
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
        try:
            for y in uid_to_go[uid]: ### Different instances for each gene associated with the uid (e.g. multiple probeset to gene associations)
                try: related_goids = y.GOIDs()
                except Exception: print uid, y;kill
                ###If there are numerical values associated with the primary unique ID, store these with the GeneRelationships instance data
                if uid in uid_value_db:
                    uid_values = uid_value_db[uid]
                    y.setUIDValues(uid_values)
                else: y.setUIDValues([])
                for goid in related_goids:
                    if goid in nested_collapsed_go_tree:
                        try: related_goid_db[goid].append(y)
                        except KeyError: related_goid_db[goid] = [y] ###store probesets linked to each goid
        except Exception: null=[]
    ###Reduce related_goid_db to unique GO to gene relationships (rather than unique ids)
    for goid in related_goid_db:
        gene_db={}
        for y in related_goid_db[goid]:
            if y.GeneID() in uid_to_go: ### Ensures that genes not in the GO term, but linked to an array id are excluded
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
            if y.GeneID() in uid_to_mapp: ### Ensures that genes not on the pathway, but linked to an array id are excluded
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
                       
def swapAndExportSystems(species_code,system1,system2):
    program_type,database_dir = unique.whatProgramIsThis()
    geneUID_import_dir = '/'+database_dir+'/'+species_code+'/uid-gene'
    ug = GrabFiles(); ug.setdirectory(geneUID_import_dir)
    filedir,file = ug.searchdirectory(system1+'-'+system2) ### Identify gene files corresponding to a particular MOD
    gene_to_uid={}; x = 0
    
    fn=filepath(filedir)
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            id1 = t[0]; id2 = t[1]
            try: gene_to_uid[id2].append(id1)
            except KeyError: gene_to_uid[id2]= [id1]
            
    export_dir = geneUID_import_dir[1:]+'/'+system2+'-'+system1+'.txt'
    try: data = export.ExportFile(export_dir)
    except Exception: data = export.ExportFile(export_dir[1:])
    
    data.write(system2+'\t'+system1+'\n')
    for id2 in gene_to_uid:
        for id1 in gene_to_uid[id2]:
            values = string.join([id2,id1],'\t')+'\n'
            data.write(values)
    data.close()

class GeneIDInfo:
    def __init__(self,system_name,geneID,pathway):
        self.system_name = system_name; self.geneID = geneID; self.pathway = pathway
        self.geneID = string.replace(geneID,'...','')
        #print self.geneID, system_name
        try:
            if ' ' == self.geneID[-1]: self.geneID = self.geneID[:-1]
        except Exception: null=[]
        
    def GeneID(self): return str(self.geneID)
    def System(self): return str(self.system_name)
    def Pathway(self): return str(self.pathway)
    def setGraphID(self,graphid): self.graphid = graphid
    def setGroupID(self,groupid): self.groupid = groupid
    def setGeneName(self,gene_name): self.gene_name = gene_name
    def setMODID(self,mod_id): self.mod_id = mod_id
    def GeneName(self): return self.gene_name
    def GraphID(self): return self.graphid
    def GroupID(self): return self.groupid
    def MODID(self): return self.mod_id
    def Report(self):
        output = self.GeneID()+'|'+self.System()
        return output
    def __repr__(self): return self.Report()
    
def parseGMT(custom_sets_folder):
    ### File format from the Broad Institutes MSigDB ()
    gm = GrabFiles(); gm.setdirectory(custom_sets_folder); system = None
    filedirs = gm.getAllFiles('.gmt') ### Identify gene files corresponding to a particular MOD
    gene_data=[]
    for gmt in filedirs:
        gmt=filepath(gmt)
        for line in open(gmt,'rU').xreadlines():             
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            pathway_name = t[0]; url = t[1]; geneids = t[2:]
            try: null=int(geneids[0]); system_name = 'Entrez Gene'
            except Exception: system_name = 'Symbol'
            for id in geneids:
                gi = GeneIDInfo(system_name,id,pathway_name)
                gene_data.append(gi)
    return gene_data

def exportWikiPathwayData(species_name,pathway_db,type):
    export_dir = 'BuildDBs/wikipathways/wikipathways_'+type+'_data_'+species_name+'.tab'
    export_data = export.ExportFile(export_dir)
        
    title = ['Pathway Name', 'Organism', 'Gene Ontology', 'Url to WikiPathways', 'Last Changed', 'Last Revision', 'Author', 'Count', 'Entrez Gene', 'Ensembl', 'Uniprot/TrEMBL', 'UniGene', 'RefSeq', 'MOD', 'PubChem', 'CAS', 'ChEBI']
    title = string.join(title,'\t')+'\n'
    export_data.write(title)
    
    for wpid in pathway_db:
        wpd = pathway_db[wpid]
        values = [wpd.Pathway(), wpd.Organism(), '', wpd.URL(), '', wpd.Revision(), '', wpd.OriginalCount(), wpd.EntrezGene(), wpd.Ensembl()]
        values +=[wpd.Uniprot(), wpd.Unigene(), wpd.Refseq(), wpd.MOD(), wpd.Pubchem(), wpd.CAS(), wpd.Chebi()]
        values = string.join(values,'\t')
        values = string.replace(values,'\n','') ###Included by mistake
        export_data.write(values+'\n')
    export_data.close()
    print 'WikiPathways data exported to:',export_dir
    
def clusterPrimaryPathwayGeneSystems(species_code,pathway_db):
    st = systemTranslation()
    for wpid in pathway_db:
        wpd = pathway_db[wpid]
        xml_data = wpd.PathwayGeneData()
        ensembl=[]; entrez=[]; refseq=[]; unigene=[]; uniprot=[]; cas=[]; chebi=[]; pubchem=[]; mod=[]
        ### If the MOD gene IDs are in the pathway then add these
        for gi in xml_data:
            if len(gi.GeneID())>0:
                source_data = gi.System() 
                if 'Ensembl' in source_data: source_data = 'Ensembl'
                if source_data in st: source_data = st[source_data] ### convert the name to the GO-Elite compatible name
                if source_data == '' and gi.GeneID()[:2]=='WP': source_data = 'WikiPathways'
                if source_data == 'Ensembl': ensembl.append(gi.GeneID())
                elif source_data == 'EntrezGene': entrez.append(gi.GeneID())
                elif string.lower(source_data) == 'refseq': refseq.append(gi.GeneID())
                elif string.lower(source_data) == 'unigene': unigene.append(gi.GeneID())
                elif 'uniprot' in string.lower(source_data): uniprot.append(gi.GeneID())
                elif string.lower(source_data) == 'swissprot': uniprot.append(gi.GeneID())
                elif 'trembl' in string.lower(source_data): uniprot.append(gi.GeneID())
                elif string.lower(source_data) == 'cas': cas.append(gi.GeneID())
                elif string.lower(source_data) == 'chebi': chebi.append(gi.GeneID())
                elif string.lower(source_data) == 'pubchem': pubchem.append(gi.GeneID())
                else: mod.append(gi.GeneID()+'('+source_data+')')
        wpd.setGeneData(ensembl,uniprot,refseq,unigene,entrez,mod,pubchem,cas,chebi)
        wpd.setOriginalCount(wpd.Count()) ### Since the mapped count will include all infered, only reference the original count

def convertAllGPML(specific_species,all_species):
    import update; import UI
    species_names = UI.remoteSpeciesInfo('yes')
    for species_code in species_names:
        if all_species == 'yes' or species_code in specific_species:
            species_name = species_names[species_code].SpeciesName()
            species = string.replace(species_name,' ','_')
            
            ### Clear and create the output dir
            try: export.deleteFolder('GPML')
            except Exception: null=[]
            os.mkdir(filepath('GPML'))
            
            ### Download all species GPML from .zip
            url = 'http://wikipathways.org//wpi/cache/wikipathways_'+species+'_Curation-AnalysisCollection__gpml.zip'
            fln,status = update.download(url,'GPML/','')
            
            if 'Internet' not in status:
                gene_data,pathway_db = parseGPML('/GPML') ### Get pathway associations (gene and annotation)
                clusterPrimaryPathwayGeneSystems(species_code,pathway_db)
                exportWikiPathwayData(species_name,pathway_db,'native')
    
                ### Create a "MAPPED" version (mapped will contain BAD MAPPINGS provide by the source database)!!!
                ensembl_to_WP = unifyGeneSystems(gene_data,species_code,'Ensembl') ### convert primary systems to MOD IDs
                hmdb_to_WP = unifyGeneSystems(gene_data,species_code,'HMDB') ### convert primary systems to MOD IDs
                convertBetweenSystems(species_code,pathway_db,ensembl_to_WP,hmdb_to_WP)
                exportWikiPathwayData(species_name,pathway_db,'mapped')
                #sys.exit()

def convertBetweenSystems(species_code,pathway_db,ensembl_to_WP,hmdb_to_WP):        
    WP_to_ensembl = OBO_import.swapKeyValues(ensembl_to_WP)
    WP_to_hmdb = OBO_import.swapKeyValues(hmdb_to_WP)

    ens_uniprot = getRelated(species_code,'Ensembl-'+'Uniprot')
    ens_refseq = getRelated(species_code,'Ensembl-'+'RefSeq')
    ens_unigene = getRelated(species_code,'Ensembl-'+'UniGene')
    ens_entrez = getRelated(species_code,'Ensembl-'+'EntrezGene')
        
    hmdb_pubchem = getRelated(species_code,'HMDB-'+'PubChem')
    hmdb_cas = getRelated(species_code,'HMDB-'+'CAS')
    hmdb_chebi = getRelated(species_code,'HMDB-'+'ChEBI')
    
    for wpid in pathway_db:
        wpd = pathway_db[wpid]
        
        try: ens_ids = WP_to_ensembl[wpd.Pathway()]
        except Exception: ens_ids=[]; #print wpid,len(WP_to_ensembl);sys.exit()
        try: hmdb_ids = WP_to_hmdb[wpd.Pathway()]
        except Exception: hmdb_ids=[]
        #print wpid,wpid.Pathway(),hmdb_ids;sys.exit()
        ensembl,uniprot,refseq,unigene,entrez,mod = wpd.GeneDataSystems()
        pubchem,cas,chebi = wpd.ChemSystems()
        uniprot = getConverted(ens_uniprot,ens_ids) +uniprot
        refseq = getConverted(ens_refseq,ens_ids) +refseq
        unigene = getConverted(ens_unigene,ens_ids) +unigene
        entrez = getConverted(ens_entrez,ens_ids) +entrez
        ensembl = ens_ids + ensembl
        
        pubchem = getConverted(hmdb_pubchem,hmdb_ids) +pubchem
        cas = getConverted(hmdb_cas,hmdb_ids) +cas
        chebi = getConverted(hmdb_chebi,hmdb_ids) +chebi
        
        ### Reset these
        wpd.setGeneData(ensembl,uniprot,refseq,unigene,entrez,mod,pubchem,cas,chebi)
        
def getConverted(gene_to_source_id,mod_ids):
    source_ids=[]
    for id in mod_ids:
        try:
            for i in gene_to_source_id[id]: source_ids.append(i)
        except Exception: null=[] 
    return source_ids

def getRelated(species_code,mod_source):
    try:
        gene_to_source_id = getGeneToUid(species_code,('hide',mod_source)); #print mod_source, 'relationships imported.'
    except Exception: gene_to_source_id={}
    return gene_to_source_id

            
class WikiPathwaysData:
    def __init__(self,pathway,wpid,revision,organism,gi):
        self.pathway = pathway; self.wpid = wpid; self.revision = revision; self.gi = gi
        self.organism = organism
    def Pathway(self): return self.pathway
    def WPID(self): return self.wpid
    def URL(self): return 'http://www.wikipathways.org/index.php/Pathway:'+self.wpid
    def Organism(self): return self.organism
    def Revision(self): return self.revision
    def PathwayGeneData(self): return self.gi
    def Report(self):
        output = self.GeneID()+'|'+self.System()
        return output
    def setGeneData(self,ensembl,uniprot,refseq,unigene,entrez,mod,pubchem,cas,chebi):
        self.ensembl=ensembl;self.uniprot=uniprot;self.refseq=refseq;self.unigene=unigene
        self.entrez=entrez;self.mod=mod;self.pubchem=pubchem;self.cas=cas;self.chebi=chebi
        combined = ensembl+uniprot+refseq+unigene+entrez+mod+pubchem+cas+chebi
        self.count = len(unique.unique(combined))
    def setOriginalCount(self,original_count): self.original_count = original_count
    def OriginalCount(self): return self.original_count
    def setInteractions(self,interactions): self.interactions = interactions
    def Interactions(self): return self.interactions
    def GeneDataSystems(self): return self.ensembl,self.uniprot,self.refseq,self.unigene,self.entrez,self.mod
    def ChemSystems(self): return self.pubchem,self.cas,self.chebi
    def Ensembl(self): return self.Join(self.ensembl)
    def Uniprot(self): return self.Join(self.uniprot)
    def Refseq(self): return self.Join(self.refseq)
    def Unigene(self): return self.Join(self.unigene)
    def EntrezGene(self): return self.Join(self.entrez)
    def MOD(self): return self.Join(self.mod)
    def Pubchem(self): return self.Join(self.pubchem)
    def CAS(self): return self.Join(self.cas)
    def Chebi(self): return self.Join(self.chebi)
    def Join(self,ls): return string.join(unique.unique(ls),',')
    def Count(self): return str(self.count)
    def __repr__(self): return self.Report()

class InteractionData:
    def __init__(self, gene1,gene2,int_type):
        self.gene1 = gene1; self.gene2 = gene2; self.int_type = int_type
    def GeneObject1(self): return self.gene1
    def GeneObject2(self): return self.gene2
    def InteractionType(self): return self.int_type

class EdgeData:
    def __init__(self, graphid1,graphid2,int_type):
        self.graphid1 = graphid1; self.graphid2 = graphid2; self.int_type = int_type
    def GraphID1(self): return self.graphid1
    def GraphID2(self): return self.graphid2
    def InteractionType(self): return self.int_type
    
def parseGPML(custom_sets_folder):
    import xml.dom.minidom
    from xml.dom.minidom import Node
    from xml.dom.minidom import parse, parseString

    gm = GrabFiles(); gm.setdirectory(custom_sets_folder); system = None
    filedirs = gm.getAllFiles('.gpml') ### Identify gene files corresponding to a particular MOD
    gene_data=[]; pathway_db={}
    for xml in filedirs:
        #if 'Wnt' in xml and 'Pluri' in xml:
        pathway_gene_data=[]; complexes_data={}; edge_data=[]; graph_node_data=[]
        xml=filepath(xml)
        filename = string.split(xml,'/')[-1]
        wpid = string.split(filename,'_')[-2]
        revision = string.split(filename,'_')[-1][:-5]
        dom = parse(xml)
        tags = dom.getElementsByTagName('Xref') ### gene IDs
        datanodes = dom.getElementsByTagName('DataNode') ### graph IDs and data types
        groups = dom.getElementsByTagName('Group') ### complexes
        edges = dom.getElementsByTagName('Point') ### interacting nodes/complexes
        pathway_tag = dom.getElementsByTagName('Pathway') ### pathway info
        
        for pn in pathway_tag:
            pathway_name = pn.getAttribute("Name")
            organism = pn.getAttribute("Organism")
        for ed in edges:
            ### Store internal graph data for pathway edges to build gene interaction networks later
            graphid = ed.getAttribute("GraphRef")
            edge_type = ed.getAttribute("ArrowHead")
            if edge_type == '': edge_pair = [graphid] ### either just a graphical line or the begining of a node-node edge
            else:
                edge_pair.append(graphid)
                edd = EdgeData(str(edge_pair[0]),str(edge_pair[1]),str(edge_type))
                edge_data.append(edd)
        for gd in groups:
            ### Specific to groups and complexes
            groupid = gd.getAttribute("GroupId") ### Group node ID (same as "GroupRef")
            graphid = gd.getAttribute("GraphId") ### GPML node ID
            complexes_data[str(groupid)] = str(graphid)
        for dn in datanodes:
            ### Specific to individual nodes
            graphid = dn.getAttribute("GraphId") ### GPML node ID
            groupid = dn.getAttribute('GroupRef')### Group node ID
            gene_name = dn.getAttribute('TextLabel')### Group node ID
            gene_name = dn.getAttribute('TextLabel')### Group node ID
            type = dn.getAttribute('Type') #E.g.', GeneProduct, Metabolite
            graph_node_data.append([graphid,groupid,gene_name,type])
        count=0 ### Keeps track of the same gene entry
        for i in tags:
            system_name = i.getAttribute("Database") ### gene ID type (e.g., Entrez Gene)
            id = i.getAttribute("ID") ### primary gene ID (e.g., ENS123)
            if len(id)>0:
                gi = GeneIDInfo(str(system_name),str(id),pathway_name)
                graphid,groupid,gene_name,type = graph_node_data[count]
                gi.setGraphID(str(graphid)); gi.setGroupID(str(groupid)) ### Include internal graph IDs for determining edges
                try: gi.setGeneName(str(gene_name))
                except Exception: gi.setGeneName(str(''))
                gene_data.append(gi)
                pathway_gene_data.append(gi)
            count+=1
        wpd=WikiPathwaysData(pathway_name,wpid,revision,organism,pathway_gene_data)
        pathway_db[wpid]=wpd
        interaction_data = getInteractions(complexes_data,edge_data,wpd)
        wpd.setInteractions(interaction_data)
    return gene_data,pathway_db

def getInteractions(complexes_data,edge_data,wpd):
    ### Annotate the interactions between nodes and between groups of nodes from WikiPathways
    interaction_data=[]
    wpid = wpd.WPID()
    graphID_db={}

    for gi in wpd.PathwayGeneData():
        graphID_db[gi.GraphID()] = [gi] ### graphID for the node
        if gi.GroupID() in complexes_data:
            graph_id = complexes_data[gi.GroupID()]
            try: graphID_db[graph_id].append(gi) ### graphID for the group of nodes with putative edges
            except Exception: graphID_db[graph_id] = [gi]

    for eed in edge_data:
        if eed.GraphID1() != '' and eed.GraphID2() != '':
            try:
                gi_list1 = graphID_db[eed.GraphID1()]
                gi_list2 = graphID_db[eed.GraphID2()]
                
                for gi1 in gi_list1:
                    for gi2 in gi_list2:
                        intd = InteractionData(gi1,gi2,eed.InteractionType())
                        interaction_data.append(intd)
            except KeyError: null=[] ### Typically occurs for interactions with Labels and similar objects
    return interaction_data

def parseGPML2():
    import urllib
    from xml.dom import minidom
    gpml = 'http://wikipathways.org//wpi/wpi.php?action=downloadFile&type=gpml&pwTitle=Pathway:WP201'
    #pathway = string.split(gpml,'/')[-1]; pathway = string.split(pathway,'.')[0]
    dom = minidom.parse(urllib.urlopen(gpml))
    tags = dom.getElementsByTagName('Xref')
    pathway_tag = dom.getElementsByTagName('Pathway')
    for pn in pathway_tag:
        pathway_name = pn.getAttribute("Name")
        organism = pn.getAttribute("Organism")
        print pathway_name,organism;kill
    gene_data=[]
    for i in tags:
        system_name = i.getAttribute("Database")
        id = i.getAttribute("ID")
        gi = GeneIDInfo(system_name,id,pathway_name)
        gene_data.append(gi)
    return gene_data

def parseBioPax3(custom_sets_folder):
    import xml.dom.minidom
    from xml.dom.minidom import Node
    from xml.dom.minidom import parse, parseString

    gm = GrabFiles(); gm.setdirectory(custom_sets_folder); system = None
    filedirs = gm.getAllFiles('.owl') ### Identify gene files corresponding to a particular MOD
    gene_data=[]
    for xml in filedirs:
        xml=filepath(xml)
        for line in open(xml,'rU').xreadlines():             
            data = cleanUpLine(line)
            if 'pathway rdf:ID=' in data: ##<pathway rdf:ID="EGFR1_pathway_Human">
                pathway_name = string.split(data,'pathway rdf:ID')[-1][:-1]
                
def processXMLline(data_type,line):
    data_types = ['NAME','pathway']
    for data_type in data_types:
        full_str='<bp:'+data_type+' rdf:datatype="xsd:string">' 
        if full_str in line:
            value = string.replace(line,full_str,'')
    
def parseBioPax(custom_sets_folder):
    #import urllib
    from xml.dom import minidom
    from xml.dom.minidom import Node
    
    import xml.dom.minidom
    from xml.dom.minidom import Node
    from xml.dom.minidom import parse, parseString

    gm = GrabFiles(); gm.setdirectory(custom_sets_folder); system = None
    filedirs = gm.getAllFiles('.owl') ### Identify gene files corresponding to a particular MOD
    gene_data=[]
 
    for xml in filedirs:
        #print xml
        #xml = 'http://www.netpath.org/data/biopax/NetPath_4.owl'
        #pathway = string.split(xml,'/')[-1]; pathway = string.split(pathway,'.')[0]
        #dom = minidom.parse(urllib.urlopen(xml))
        xml=filepath(xml)
        dom = parse(xml); pathway_name=[]
        tags = dom.getElementsByTagName('unificationXref')
        pathway_tag = dom.getElementsByTagName('pathway')
        for pathway_name in pathway_tag: pathway_name = pathway_name.getAttribute("rdf:ID")
        if len(pathway_name)==0:
            pathway_tag = dom.getElementsByTagName('bp:NAME')
            for pathway_name in pathway_tag: pathway_name = pathway_name.getAttribute("rdf:datatype='xsd:string'")
            #print pathway_name
        for node in tags:
            id_object =  node.getElementsByTagName("ID")
            db_object =  node.getElementsByTagName("DB")
            for node2 in id_object:
                id = ''
                for node3 in node2.childNodes:
                    if node3.nodeType == Node.TEXT_NODE: id+= node3.data
            for node2 in db_object:
                system_name = ''
                for node3 in node2.childNodes:
                    if node3.nodeType == Node.TEXT_NODE: system_name+= node3.data
            gi = GeneIDInfo(system_name,id,pathway_name)
            gene_data.append(gi)
    #print gene_data
    return gene_data

def parseBioPax2(filename):
    import urllib
    from xml.dom import minidom
    from xml.dom.minidom import Node
    
    xml = 'http://www.netpath.org/data/biopax/NetPath_4.owl'
    #pathway = string.split(xml,'/')[-1]; pathway = string.split(pathway,'.')[0]
    dom = minidom.parse(urllib.urlopen(xml))          
    tags = dom.getElementsByTagName('unificationXref')
    pathway_tag = dom.getElementsByTagName('pathway')
    for pathway_name in pathway_tag: pathway_name = pathway_name.getAttribute("rdf:ID")
    gene_data=[]
    for node in tags:
        id_object =  node.getElementsByTagName("ID")
        db_object =  node.getElementsByTagName("DB")
        for node2 in id_object:
            id = ''
            for node3 in node2.childNodes:
                if node3.nodeType == Node.TEXT_NODE: id+= node3.data
        for node2 in db_object:
            system_name = ''
            for node3 in node2.childNodes:
                if node3.nodeType == Node.TEXT_NODE: system_name+= node3.data
        gi = GeneIDInfo(system_name,id,pathway_name)
        gene_data.append(gi)
    return gene_data

def systemTranslation():
    st={}
    st['Entrez Gene']='EntrezGene'
    st['refseq']='RefSeq'
    st['uniprot']='Uniprot'
    st['TubercuList']='Ensembl'
    st['Kegg Compound']='KeggCompound'
    st['Uniprot/TrEMBL']='Uniprot'
    st['SwissProt']='Uniprot'
    return st
    
def unifyGeneSystems(xml_data,species_code,mod):
    systems_db={}
    for gi in xml_data:
        try:
            try: systems_db[gi.System()]+=1
            except KeyError: systems_db[gi.System()]=1
        except Exception: print gi;kill
        
    #for i in systems_db: print i
    
    ### Import and combine all secondary systems
    system_ids={}
    st = systemTranslation();
    for source_data in systems_db:
        original_source = source_data
        if 'Ensembl' in source_data: source_data = 'Ensembl'
        if source_data in st: source_data = st[source_data] ### convert the name to the GO-Elite compatible name
        if source_data != mod:
            mod_source = mod+'-'+source_data
            try: gene_to_source_id = getGeneToUid(species_code,('hide',mod_source)); #print mod_source, 'relationships imported.'
            except Exception: gene_to_source_id={}
            #print len(gene_to_source_id),mod_source,source_data
            source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
            system_ids[original_source]=source_to_gene
            
    for system in system_ids:
        if system == 'Symbol':
            source_to_gene = system_ids[system]
            source_to_gene = lowerAllIDs(source_to_gene)
            system_ids[system] = source_to_gene
            
    mod_pathway={}
    ### Convert source IDs to MOD IDs
    for gi in xml_data:
        if gi.System() in system_ids:
            source_to_gene = system_ids[gi.System()]
            geneID = gi.GeneID()
            if gi.System() == 'Symbol':
                geneID = string.lower(geneID)
            if geneID in source_to_gene:
                for mod_id in source_to_gene[geneID]:
                    try: mod_pathway[mod_id].append(gi.Pathway())
                    except Exception: mod_pathway[mod_id] = [gi.Pathway()]
                    gi.setMODID(mod_id)
                    #if 'Trp53' == gi.GeneName(): print gi.GeneID(),mod_id,gi.System(),mod,'a'
            else:
                null=[]
                #print [gi.GeneID()]
    
    ### If the MOD gene IDs are in the pathway then add these
    for gi in xml_data:
        source_data = gi.System() 
        if 'Ensembl' in source_data: source_data = 'Ensembl'
        if source_data in st: source_data = st[source_data] ### convert the name to the GO-Elite compatible name
        if source_data == mod:
            try: mod_pathway[gi.GeneID()].append(gi.Pathway())
            except Exception: mod_pathway[gi.GeneID()] = [gi.Pathway()]
            gi.setMODID(gi.GeneID())
            #if 'Trp53' == gi.GeneName(): print gi.GeneID(),mod_id,gi.System(),mod,'b'
            
    #print len(system_ids),len(mod_pathway),len(mod_pathway)
    return mod_pathway

def lowerAllIDs(source_to_gene):
    source_to_gene2={}
    for source in source_to_gene:
        source_to_gene2[string.lower(source)] = source_to_gene[source]
    return source_to_gene2

def combineDBs(db1,db2):
    for id in db1:
        if id in db2: db2[id]+=db1[id]
        else: db2[id]=db1[id]
    return db2

if __name__ == '__main__':
    species_code = 'Mx'; mod = 'HMDB'; gotype='nested'
    filedir = 'C:/Documents and Settings/Nathan Salomonis/My Documents/GO-Elite_120beta/Databases/EnsMart56Plus/Ce/gene/EntrezGene.txt'
    system = 'Macaroni'
    
    #biopax_data = parseBioPax('/test'); sys.exit()
    
    import GO_Elite
    system_codes,source_types,mod_types = GO_Elite.getSourceData()
    custom_sets_folder = '/test'
    #importGeneCustomData(species_code,system_codes,custom_sets_folder,mod); sys.exit()

    ### Test interactions export    
    custom_sets_folder = 'C:/Users/Nathan Salomonis/Desktop/GO-Elite/GO-Elite_v.1.2.2-Win64/GPML'
    species_code = 'Hs'; mod = 'Ensembl'
    gpml_data,pathway_db = parseGPML(custom_sets_folder)
    gene_to_WP = unifyGeneSystems(gpml_data,species_code,mod)
    exportNodeInteractions(pathway_db,mod,system_codes,custom_sets_folder)
    sys.exit()
    
    mod = 'Ensembl'
    gene_to_BioPax = unifyGeneSystems(biopax_data,species_code,mod)
    #gene_to_WP = unifyGeneSystems(gpml_data,species_code,mod)
    #gene_to_BioPax = combineDBs(gene_to_WP,gene_to_BioPax)
    for i in gene_to_BioPax:
        print i, gene_to_BioPax[i]
    print len(gene_to_BioPax); sys.exit()

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


        
    
