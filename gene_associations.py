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
from stats_scripts import statistics
try: from import_scripts import OBO_import
except Exception: pass
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
        if entry[-4:] == ".txt" or entry[-4:] == ".csv" or ".owl" in entry or ".gpml" in entry or ".xml" in entry or ".gmt" in entry or ".tab" in entry: dir_list2.append(entry)
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
    def getMatchingFolders(self,search_term):
        dir_list = unique.read_directory(self.data); matching = ''
        root_dir = filepath('')
        for folder in dir_list:
            if search_term in folder and '.' not in folder:
                matching = filepath(self.data[1:]+'/'+folder)
                if root_dir not in matching:
                    matching = filepath(root_dir +'/'+ matching) ### python bug????
        return matching
        
def getDirectoryFiles(import_dir, search_term):
    exact_file = ''; exact_file_dir=''; all_matching=[]
    if '.txt' in import_dir:
        import_dir = export.findParentDir(import_dir)
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        present = verifyFile(import_dir+'/'+data) ### Check to see if the file is present in formatting the full file_dir
        if present == True:
            affy_data_dir = import_dir+'/'+data
        else: affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir and '._' not in affy_data_dir:
            if 'version.txt' not in affy_data_dir: exact_file_dir = affy_data_dir; exact_file = data; all_matching.append(exact_file_dir)
    return all_matching, exact_file_dir,exact_file

def verifyFile(filename):
    present = False
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines(): present = True; break
    except Exception: present = False
    return present

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

def importGenericDB(filename):
    try:
        key_db=collections.OrderedDict() ### Retain order if possible
    except Exception:
        key_db={}
    fn=filepath(filename); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if data[0]=='#': x=x
        elif x==0:
            x=1
            headers = t
        else:
            try: key_db[t[0]].append(t[1:])
            except Exception: key_db[t[0]] = [t[1:]]
    return key_db, headers

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
        export_dir = database_dir+'/'+species_code+'/uid-gene/'+mod+'-Symbol'+'.txt'
        data = export.ExportFile(export_dir)
        data.write('GeneID\tSymbol\n')
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
    OBO_import.buildNestedOntologyAssociations(species_code,export_databases,mod_types,genmapp_mod)
        
def importGeneToOntologyData(species_code,mod,gotype,ontology_type):
    program_type,database_dir = unique.whatProgramIsThis()
    if gotype == 'nested':
        geneGO_import_dir = '/'+database_dir+'/'+species_code+'/nested'
        if ontology_type == 'GeneOntology': ontology_type = 'GO.'
        mod += '_to_Nested'
    elif gotype == 'null':
        geneGO_import_dir = '/'+database_dir+'/'+species_code+'/gene-go'
    gg = GrabFiles(); gg.setdirectory(geneGO_import_dir)
    try: filedir,file = gg.searchdirectory(mod+'-'+ontology_type) ###Identify gene files corresponding to a particular MOD
    except Exception:
        if gotype == 'nested':
            buildNestedAssociations(species_code)
            filedir,file = gg.searchdirectory(mod+'-'+ontology_type)
    global gene_to_go; x=0
    
    fn=filepath(filedir); gene_to_go={}
    for line in open(fn,'rU').xreadlines():
        if gotype == 'nested' and x==0: x = 1
        else:
            if 'Ontology' not in line and 'ID' not in line: ### Header included in input from AP or BioMart
                #data = cleanUpLine(line)
                data = line.strip()
                t = string.split(data,'\t')
                if len(t)>1:
                    try: gene = t[0]; goid = t[1]
                    except IndexError: print [line],t, 'Ontology-gene relationship not valid...please clean up',filedir;sys.exit()
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
    global gene_to_mapp; x = True
    
    fn=filepath(filedir); gene_to_mapp={}
    for line in open(fn,'rU').xreadlines():             
        #data = cleanUpLine(line)
        data = line.strip()
        if x: x=False
        else:
            t = string.split(data,'\t')
            gene = t[0]; mapp = t[2]
            try: gene_to_mapp[gene].append(mapp)
            except KeyError: gene_to_mapp[gene]= [mapp]
    return gene_to_mapp

def exportCustomPathwayMappings(gene_to_custom,mod,system_codes,custom_sets_folder):
    if '.txt' in custom_sets_folder:
        export_dir = custom_sets_folder
    else:
        export_dir = custom_sets_folder+'/CustomGeneSets/custom_gene_set.txt'
    #print 'Exporting:',export_dir
    try: data = export.ExportFile(export_dir)
    except Exception: data = export.ExportFile(export_dir[1:])
    for system_code in system_codes:
        if system_codes[system_code] == mod: mod_code = system_code
    
    data.write('GeneID\t'+mod+'\tPathway\n'); relationships=0
    
    gene_to_custom = eliminate_redundant_dict_values(gene_to_custom)
    for gene in gene_to_custom:
        for pathway in gene_to_custom[gene]:
            values = string.join([gene,mod_code,pathway],'\t')+'\n'
            data.write(values); relationships+=1
    data.close()
    #print relationships,'Custom pathway-to-ID relationships exported...'

def exportNodeInteractions(pathway_db,mod,custom_sets_folder):
    
    import GO_Elite
    system_codes,source_types,mod_types = GO_Elite.getSourceData()
    
    export_dir = custom_sets_folder+'/Interactomes/interactions.txt'
    print 'Exporting:',export_dir
    try: data = export.ExportFile(export_dir)
    except Exception: data = export.ExportFile(export_dir[1:])
    for system_code in system_codes:
        if system_codes[system_code] == mod: mod_code = system_code

    try: gene_to_symbol_db = getGeneToUid(species_code,('hide',mod+'-Symbol.txt')); #print mod_source, 'relationships imported.'
    except Exception: gene_to_symbol_db={}
            
    data.write('Symbol1\tInteractionType\tSymbol2\tGeneID1\tGeneID2\tPathway\n'); relationships=0

    for pathway_id in pathway_db:
        wpd = pathway_db[pathway_id]
        for itd in wpd.Interactions():
            gi1 = itd.GeneObject1()
            gi2 = itd.GeneObject2()
            if len(gene_to_symbol_db)>0:
                try:
                    if gi1.ModID()[0] in gene_to_symbol_db:
                        symbol = gene_to_symbol_db[gi1.ModID()[0]][0]
                        if len(symbol)>0:
                            gi1.setLabel(symbol) ### Replace the WikiPathways user annnotated symbol with a MOD symbol
                except Exception:
                    None
                try:
                    if gi2.ModID()[0] in gene_to_symbol_db:
                        symbol = gene_to_symbol_db[gi2.ModID()[0]][0]
                        if len(symbol)>0:
                            gi2.setLabel(symbol) ### Replace the WikiPathways user annnotated symbol with a MOD symbol
                except Exception:
                    None
            try:
                values = string.join([gi1.Label(),itd.InteractionType(),gi2.Label(),gi1.ModID()[0],gi2.ModID()[0],wpd.Pathway()],'\t')
                relationships+=1
                try:
                    values = cleanUpLine(values)+'\n' ### get rid of any end-of lines introduced by the labels
                    data.write(values)
                except Exception:
                    try: ### Occurs when 'ascii' codec can't encode characters due to UnicodeEncodeError
                        values = string.join(['',itd.InteractionType(),gi2.Label(),gi1.ModID()[0],gi2.ModID()[0],wpd.Pathway()],'\t')
                        values = cleanUpLine(values)+'\n' ### get rid of any end-of lines introduced by the labels
                        data.write(values)
                    except Exception:
                        values = string.join([gi1.Label(),itd.InteractionType(),'',gi1.ModID()[0],gi2.ModID()[0],wpd.Pathway()],'\t')
                        values = cleanUpLine(values)+'\n' ### get rid of any end-of lines introduced by the labels
                        try: data.write(values)
                        except Exception: None ### Occurs due to illegal characters
            except AttributeError,e:
                #print e
                null=[] ### Occurs if a MODID is not present for one of the nodes
            
    data.close()
    print relationships,'Interactions exported...'
    
def importGeneCustomData(species_code,system_codes,custom_sets_folder,mod):
    print 'Importing custom pathway relationships...'
    #print 'Trying to import text data'
    gene_to_custom = importTextCustomData(species_code,system_codes,custom_sets_folder,mod)
    #print len(gene_to_custom)
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
    if len(gene_to_WP)>0: ### Export all pathway interactions
        try: exportNodeInteractions(pathway_db,mod,custom_sets_folder)
        except Exception: null=[]
    """
    ### Combine WikiPathway associations with the custom
    try: gene_to_mapp = importGeneMAPPData(species_code,mod)
    except Exception: gene_to_mapp = {}
    for gene in gene_to_mapp:
        for mapp in gene_to_mapp[gene]:
            try: gene_to_custom[gene].append(mapp)
            except KeyError: gene_to_custom[gene]= [mapp]"""

    return gene_to_custom

def importTextCustomData(species_code,system_codes,custom_sets_folder,mod):
    program_type,database_dir = unique.whatProgramIsThis()
    gm = GrabFiles(); gm.setdirectory(custom_sets_folder); system = None
    filedirs = gm.getAllFiles('.txt')
    global gene_to_custom
    gene_to_custom={}
    file_gene_to_custom={}
    for filedir in filedirs:
        try:
            file = string.split(filedir,'/')[-1]
            print "Reading custom gene set",filedir
            fn=filepath(filedir); x = 1
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                if x==0: x=1
                else:
                    x+=1
                    t = string.split(data,'\t')
                    try:
                        gene = t[0]; mapp = t[2]; system_code = t[1]
                        if system_code in system_codes: system = system_codes[system_code]
                        else: 
                            if x == 3: print system_code, "is not a recognized system code. Skipping import of",file; break
                    except Exception:
                        if len(t)>0:
                            gene = t[0]
                            if len(t)==1: ### Hence, no system code is provided and only one gene-set is indicated per file
                                source_data = predictIDSource(t[0],system_codes)
                                if len(source_data)>0: system = source_data; mapp = file
                                else:
                                    if x == 3: print file, 'is not propperly formatted (skipping import of relationships)'; break
                            elif len(t)==2:
                                if t[1] in system_codes: system = system_codes[t[1]]; mapp = file ### Hence, system code is provided by only one gene-set is indicated per file
                                else:
                                    source_data = predictIDSource(t[0],system_codes)
                                    if len(source_data)>0: system = source_data; mapp = t[1]
                                    else:
                                        if x == 3: print file, 'is not propperly formatted (skipping import of relationships)'; break
                        else: continue ### Skip line
                    try: file_gene_to_custom[gene].append(mapp)
                    except KeyError: file_gene_to_custom[gene]= [mapp]
                    
            #print [system, mod, len(file_gene_to_custom)]
            ### If the system code is not the MOD - Convert to the MOD
            if (system != mod) and (system != None):
                mod_source = 'hide',mod+'-'+system+'.txt'
                try: gene_to_source_id = getGeneToUid(species_code,mod_source)
                except Exception: print mod_source,'relationships not found. Skipping import of',file; break
                source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
                if system == 'Symbol': source_to_gene = lowerAllIDs(source_to_gene)
                for source_id in file_gene_to_custom:
                    original_source_id = source_id ### necessary when Symbol
                    if system == 'Symbol': source_id = string.lower(source_id)
                    if source_id in source_to_gene:
                        for gene in source_to_gene[source_id]:
                            try: gene_to_custom[gene] += file_gene_to_custom[original_source_id]
                            except Exception: gene_to_custom[gene] = file_gene_to_custom[original_source_id]
            else:
                for gene in file_gene_to_custom:
                    try: gene_to_custom[gene] += file_gene_to_custom[gene]
                    except Exception: gene_to_custom[gene] = file_gene_to_custom[gene]
        except Exception:
            print file, 'not formatted propperly!'
        file_gene_to_custom={} ### Clear this file specific object
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
        #data = cleanUpLine(line)
        data = line.strip()
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            gene = t[0]; uid = t[1]
            try: gene_to_uid[gene].append(uid)
            except KeyError: gene_to_uid[gene]= [uid]
    return gene_to_uid

def augmentEnsemblGO(species_code):
    ontology_type = 'GeneOntology'
    entrez_ens = importUidGeneSimple(species_code,'EntrezGene-Ensembl')
    try: ens_go=importGeneToOntologyData(species_code,'Ensembl','null',ontology_type)
    except Exception: ens_go = {}
    try: entrez_go=importGeneToOntologyData(species_code,'Entrez','null',ontology_type)
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
    full_path_db,path_id_to_goid,null = OBO_import.buildNestedOntologyAssociations(species_code,export_databases,['Ensembl'],genmapp_mod)

def importOntologyUIDGeneData(species_code,mod_source,gene_to_go,denominator_source_ids):
    program_type,database_dir = unique.whatProgramIsThis()
    geneUID_import_dir = '/'+database_dir+'/'+species_code+'/uid-gene'
    ug = GrabFiles(); ug.setdirectory(geneUID_import_dir)
    #print "begining to parse",mod_source, 'from',geneGO_import_dir
    filedir,filename = ug.searchdirectory(mod_source) ###Identify gene files corresponding to a particular MOD
    fn=filepath(filedir); uid_to_go={}; count=0; x=0
    uid_system,gene_system = grabFileRelationships(filename)
    for line in open(fn,'rU').xreadlines(): count+=1
    original_increment = int(count/10); increment = original_increment
    for line in open(fn,'rU').xreadlines():           
        #data = cleanUpLine(line); x+=1
        data = line.strip(); x+=1
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
        except Exception: null=[]
    return uid_to_go,uid_system

def importGeneSetUIDGeneData(species_code,mod_source,gene_to_mapp,denominator_source_ids):
    program_type,database_dir = unique.whatProgramIsThis()
    geneUID_import_dir = '/'+database_dir+'/'+species_code+'/uid-gene'
    ug = GrabFiles(); ug.setdirectory(geneUID_import_dir)
    #print "begining to parse",mod_source, 'from',geneGO_import_dir
    filedir,filename = ug.searchdirectory(mod_source) ###Identify gene files corresponding to a particular MOD
    fn=filepath(filedir); uid_to_mapp={}; count=0; x=0
    uid_system,gene_system = grabFileRelationships(filename)
    for line in open(fn,'rU').xreadlines(): count+=1
    original_increment = int(count/10); increment = original_increment
    for line in open(fn,'rU').xreadlines():           
        #data = cleanUpLine(line); x+=1
        data = line.strip(); x+=1
        if program_type == 'GO-Elite':
            if x == increment: increment+=original_increment; print '*',
        t = string.split(data,'\t')
        uid = t[1]; gene = t[0] 
        try:
            if len(denominator_source_ids)>0:
                null=denominator_source_ids[uid] ### requires that the source ID be in the list of analyzed denominators
            mapp_name = gene_to_mapp[gene]
            y = GeneRelationships(uid,gene,mapp_name,uid_system,gene_system)
            try: uid_to_mapp[uid].append(y)
            except KeyError: uid_to_mapp[uid] = [y]
        except Exception: null=[]
    return uid_to_mapp,uid_system

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database: list = unique.unique(database[key]); list.sort(); db1[key] = list
    return db1

def getGeneToUid(species_code,mod_source,display=True):
    if 'hide' in mod_source: show_progress, mod_source = mod_source
    elif display==False: show_progress = 'no'
    else: show_progress = 'yes'

    program_type,database_dir = unique.whatProgramIsThis()
    import_dir = '/'+database_dir+'/'+species_code+'/uid-gene'
    ug = GrabFiles(); ug.setdirectory(import_dir)
    filedir,filename = ug.searchdirectory(mod_source) ###Identify gene files corresponding to a particular MOD
    fn=filepath(filedir); gene_to_uid={}; count = 0; x=0
    uid_system,gene_system = grabFileRelationships(filename)

    for line in open(fn,'r').xreadlines(): count+=1
    original_increment = int(count/10); increment = original_increment
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line); x+=1
        if x == increment and show_progress == 'yes':
            increment+=original_increment; print '*',
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
        #data = cleanUpLine(line)
        data = line.strip()
        data = string.split(data,'\t')
        probeset_db[data[0]]=[]
    return probeset_db

def predictIDSource(id,system_codes):
    au = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    nm = ['1','2','3','4','5','6','7','8','9']
    affy_suffix = '_at'; ensembl_prefix = 'ENS'; source_data = ''; id_type = 'Symbol'; id_types={}
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
            try: id_types[id_type]+=1
            except Exception: id_types[id_type]=1

    ###If the user changes the names of the above id_types, we need to verify that the id_type is in the database
    if len(id_type)>0 and len(id_types)>0:
        id_type_count=[]
        for i in id_types:
            id_type_count.append((id_types[i],i))
        id_type_count.sort()
        #print id_type_count
        id_type = id_type_count[-1][-1]
    for code in system_codes:
        if system_codes[code] == id_type: source_data = id_type
    return source_data

def predictIDSourceSimple(id):
    au = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    nm = ['1','2','3','4','5','6','7','8','9']
    affy_suffix = '_at'; ensembl_prefix = 'ENS'; source_data = ''; id_type = 'Sy'; id_types={}
    if len(id)>3:
        if affy_suffix == id[-3:]: id_type = 'X'
        elif ensembl_prefix == id[:3]:
            if ' ' in id:
                id_type = 'En:Sy'
            else:
                id_type = 'En'
        elif id[2] == '.': id_type = 'Ug'
        #elif (id[0] in au and id[1] in nm) or '_' in id: id_type = 'S'
        else:
            try:
                int_val = int(id)
                if len(id) == 7: id_type = 'X' ###All newer Affymetrix transcript_cluster_ids and probesets (can be mistaken for EntrezGene IDs)
                else: id_type = 'L'
            except ValueError: null = []
        if id_type != 'En' and ensembl_prefix in id and ':' in id:
            prefix = string.split(id,':')[0]
            if ensembl_prefix not in prefix and ' ' in id:
                id_type = '$En:Sy'
            else:
                id_type = 'Ae'
    return id_type

def addNewCustomSystem(filedir,system,save_option,species_code):
    print 'Adding new custom system (be patient)' ### Print statement here forces the status window to appear quicker otherwise stalls
    gene_annotations={}
    if 'update' in save_option:
        print "Importing existing",species_code,system,"relationships."
        gene_annotations = importGeneData(species_code,mod)
    fn=filepath(filedir); length3=0; lengthnot3=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if len(t[0])>0:
            try:
                ### Allow the user to only include one column here
                gene = t[0]
                try: symbol = t[1]
                except Exception: symbol = gene
                try: name = ''
                except Exception: name = ''
                s = GeneAnnotations(gene,symbol,name,'')
                gene_annotations[gene] = s
            except Exception:
                print 'Unexpected error in the following line'; print t, filedir;kill
        
    if len(gene_annotations)>0:
        print 'Writing new annotation file:',species_code,system
        export_dir = 'Databases/'+species_code+'/gene/'+system+'.txt'
        data = export.ExportFile(export_dir)
        data.write('UID\tSymbol\tDescription\n')
        for gene in gene_annotations:
            s = gene_annotations[gene]
            data.write(string.join([gene,s.Symbol(),s.Description()],'\t')+'\n')
        data.close()
        return 'exported'
    else: return 'not-exported'
    
def importGeneSetsIntoDatabase(source_file,species_code,mod):
    source_filename = export.findFilename(source_file)
    export.deleteFolder('BuildDBs/temp') ### Delete any previous data
    destination_dir = filepath('BuildDBs/temp/'+source_filename)
    export.customFileCopy(source_file,destination_dir)
    #print destination_dir
    custom_sets_folder = export.findParentDir(destination_dir)
    import GO_Elite; system_codes,source_types,mod_types = GO_Elite.getSourceData()
    gene_to_custom = importGeneCustomData(species_code,system_codes,custom_sets_folder,mod)
    return gene_to_custom

def addNewCustomRelationships(filedir,relationship_file,save_option,species_code):
    print 'Adding new custom relationships (be patient)' ### Print statement here forces the status window to appear quicker otherwise stalls
    relationship_filename = export.findFilename(relationship_file)
    mod,data_type = string.split(relationship_filename,'-')
    mod_id_to_related={}
    if 'update' in save_option or '.txt' not in filedir:
        if 'gene-mapp' in relationship_file:
            if '.txt' not in filedir: ### Then these are probably gmt, gpml or owl
                mod_id_to_related = importGeneSetsIntoDatabase(filedir,species_code,mod)
            else:
                mod_id_to_related = importGeneMAPPData(species_code,mod)
        elif 'gene-go' in relationship_file:
            mod_id_to_related=importGeneToOntologyData(species_code,mod,'null',data_type)
        else: mod_id_to_related=importUidGeneSimple(species_code,mod+'-'+data_type+'.txt')
    fn=filepath(filedir); length2=0; lengthnot2=0
    if '.txt' in fn:
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            data_length = len(t)
            ### Check the length of the input file
            if data_length==2 or data_length==3: #This can occur if users load a custom gene set file created by GO-Elite
                length2+=1
                if len(t[0])>0 and len(t[-1])>0:
                    if ',' in t[0]: keys = string.split(t[0],',') ### These can be present in the MGI phenotype ontology gene association files
                    else: keys = [t[0]]
                    for key in keys:
                        try: mod_id_to_related[key].append(t[-1])
                        except KeyError: mod_id_to_related[key] = [t[-1]]
            else: lengthnot2+=1
    else: length2+=1
        
    if length2>lengthnot2:
        ###Ensure that the input file is almost completely 2 columns
        print 'Writing new relationship file:',species_code,relationship_file
        if 'Databases' in relationship_file:
            export_dir = relationship_file ### Sometimes we just include the entire path for clarification
            if 'mapp' in export_dir: data_type = 'MAPP'
            else: data_type = 'Ontology'
        elif data_type == 'MAPP': export_dir = 'Databases/'+species_code+'/gene-mapp/'+relationship_file+'.txt'
        elif data_type == 'Ontology':
            export_dir = 'Databases/'+species_code+'/gene-go/'+relationship_file+'.txt'
        else: export_dir = 'Databases/'+species_code+'/uid-gene/'+relationship_file+'.txt'
        if data_type == 'Ontology':
            ### To trigger a rebuild of the ontology nested, must deleted the existing nested for this ontology
            nested_path = filepath(string.replace(relationship_file,'gene-go','nested'))
            nested_path = string.replace(nested_path,'GeneOntology','GO')
            obo=string.split(nested_path,'-')[-1]
            nested_path = string.replace(nested_path,'-'+obo,'_to_Nested-'+obo) ### Will break if there is another '-' in the name
            try: os.remove(nested_path)
            except Exception: null=[]
        data = export.ExportFile(export_dir)
        data.write('UID\t\tRelated ID\n')
        for mod_id in mod_id_to_related:
            for related_id in mod_id_to_related[mod_id]:
                if data_type == 'MAPP': ###Pathway files have 3 rather than 2 columns
                    data.write(string.join([mod_id,'',related_id],'\t')+'\n')
                else: data.write(string.join([mod_id,related_id],'\t')+'\n')
        data.close()
        return 'exported'
    else: return 'not-exported'
                
def importUIDsForMAPPFinderQuery(filedir,system_codes,return_uid_values):
    fn=filepath(filedir); x=0; uid_list = {}; uid_value_db = {}; source_data = ''; value_len_list = []
    system_code_found = 'no'; first_uid=''; first_line = 'no lines in file'; underscore = False
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
                        system = string.replace(system,' ','')
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
                                if '_' in uid: underscore=True
                    except Exception:
                        source_data_read =''
                if len(source_data_read)<1:
                    for uid in uids:
                        source_data_read = predictIDSource(uid,system_codes)
                    if len(source_data_read)>0: source_data = source_data_read
                try: source_data_db[source_data_read]+=1
                except Exception: source_data_db[source_data_read]=1
    first_uid = string.replace(first_uid,'Worksheet','!!!!')
    filenames = string.split(filedir,'/'); filename = filenames[-1]
    if x==1:
        error = 'No results in input file:'+filename
        print error
    elif '!!!!' in first_uid:
        error = 'WARNING!!! There appears to be a formatting file issue with the file:\n"'+filename+'"\nPlease correct and re-run (should be tab-delimited text).'
    elif len(source_data_db)>1:
        #error = 'WARNING!!! There is more than one gene system (e.g., Ensembl and EntrezGene) in the file:\n"'+filename+'"\nPlease correct and re-run.'
        sources = []
        for s in source_data_db:
            sources.append([source_data_db[s],s])
        sources.sort(); source_data = sources[-1][1]
        #print 'Using the system code:', source_data, "(multiple systems present)"
    elif source_data == '':
        error = 'WARNING!!! No System Code identified in:\n"'+filename+'"\nPlease provide in input text file and re-run.\n If System Code column present, the file format may be incorrect.'
        try: error +='Possible system code: '+system+' not recognized.'
        except Exception: None
    elif x>0 and bad_alphanumerics>1:
        error = 'WARNING!!! Invalid text file encoding found (invalid alphanumeric values). Please resave text file in a standard tab-delimited format'
    if underscore and system == 'S':
        error += '\nSwissProt IDs of the type P53_HUMAN (P04637) are not recognized as propper UniProt IDs. Consider using alternative compatible IDs (e.g., P04637).'
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

def grabNestedGeneToOntologyAssociations(species_code,mod,source_data,system_codes,denom_search_dir,ontology_type):
    program_type,database_dir = unique.whatProgramIsThis()
    global printout; printout = 'yes'
    gene_annotations = importGeneData(species_code,mod)
    ### Filter source IDs for those in all user denominator files
    denominator_source_ids = getAllDenominators(denom_search_dir,system_codes)
    try: gene_to_go = importGeneToOntologyData(species_code,mod,'nested',ontology_type)
    except Exception:
        print "Warning...the MOD you have selected:",mod,ontology_type,"is missing the appropriate relationship files",
        print "necessary to run GO-Elite.  Either replace the missing files ("+database_dir+'/'+species_code+') or select a different MOD at runtime.'
        print 'Exiting program.'; forceExit
    if source_data != mod:
        mod_source = mod+'-'+source_data+'.txt'
        uid_to_go,uid_system = importOntologyUIDGeneData(species_code,mod_source,gene_to_go,denominator_source_ids)
        gene_to_go = convertGeneDBtoGeneObjects(gene_to_go,mod)
        for gene in gene_to_go: uid_to_go[gene] = gene_to_go[gene]
        return uid_to_go, uid_system, gene_annotations
    else:
        gene_to_go = convertGeneDBtoGeneObjects(gene_to_go,mod)
        return gene_to_go, mod, gene_annotations

def grabNestedGeneToPathwayAssociations(species_code,mod,source_data,system_codes,custom_sets_folder,denom_search_dir,ontology_type):
    program_type,database_dir = unique.whatProgramIsThis()
    global printout; printout = 'yes'
    gene_annotations = importGeneData(species_code,mod)
    ### Filter source IDs for those in all user denominator files
    denominator_source_ids = getAllDenominators(denom_search_dir,system_codes)
    try:
        if ontology_type == 'WikiPathways': ontology_type = 'MAPP'
        gene_to_mapp = importGeneMAPPData(species_code,mod+'-'+ontology_type)
    except Exception:
        ### If a custom collection, shouldn't be found in the gene-mapp directory
        try: gene_to_mapp = importGeneCustomData(species_code,system_codes,custom_sets_folder,mod)
        except Exception: gene_to_mapp = {}
    if source_data != mod:
        mod_source = mod+'-'+source_data+'.txt'
        uid_to_mapp,uid_system = importGeneSetUIDGeneData(species_code,mod_source,gene_to_mapp,denominator_source_ids)
        gene_to_mapp = convertGeneDBtoGeneObjects(gene_to_mapp,mod)
        ### Combine the gene and UID data in order to later filter out gene's linked to an arrayid that are not linked to the pathway (common for some pathways)
        for gene in gene_to_mapp: uid_to_mapp[gene] = gene_to_mapp[gene]
        return uid_to_mapp, uid_system, gene_annotations
    else:
        gene_to_mapp = convertGeneDBtoGeneObjects(gene_to_mapp,mod)
        return gene_to_mapp, mod, gene_annotations

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

def simpleGenePathwayImport(species,geneset_type,pathway,OntologyID,directory):
    ###Import a gene-set and only store geneIDs for that pathway
    associated_IDs={}
    if geneset_type == 'WikiPathways': geneset_type = 'MAPP'
    filename = 'AltDatabase/goelite/'+species+'/'+directory+'/'+'Ensembl-'+geneset_type+'.txt'
    if directory == 'nested':
        if geneset_type == 'GeneOntology': geneset_type = 'GO'
        if len(OntologyID)>1: pathway = OntologyID
        else: pathway = lookupOntologyID(geneset_type,pathway) ### Translates from name to ID
        filename = 'AltDatabase/goelite/'+species+'/'+directory+'/'+'Ensembl_to_Nested-'+geneset_type+'.txt'
    
    fn=filepath(filename)
    ### Imports a geneset category and stores pathway-level names
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        #print [t[-1]],[pathway,OntologyID]
        if t[-1] == pathway:
            associated_IDs[t[0]] = None
    return associated_IDs

def filterGeneToUID(species_code,mod,source,filter_db):
    mod_source = mod+'-'+source
    if 'hide' in mod_source: show_progress, mod_source = mod_source
    else: show_progress = 'yes'
    program_type,database_dir = unique.whatProgramIsThis()
    import_dir = '/'+database_dir+'/'+species_code+'/uid-gene'
    ug = GrabFiles(); ug.setdirectory(import_dir)
    filedir,filename = ug.searchdirectory(mod_source) ###Identify gene files corresponding to a particular MOD
    fn=filepath(filedir); gene_to_uid={}; count = 0; x=0
    uid_system,gene_system = grabFileRelationships(filename)

    for line in open(fn,'rU').xreadlines():             
        #data = cleanUpLine(line)
        data = line.strip()
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            if t[0] in filter_db or t[1] in filter_db:
                uid = t[1]; gene = t[0]
                try: gene_to_uid[uid].append(gene)
                except KeyError: gene_to_uid[uid] = [gene]
            elif len(filter_db)==0:
                uid = t[1]; gene = t[0]
                try: gene_to_uid[uid].append(gene)
                except KeyError: gene_to_uid[uid] = [gene]
    gene_to_uid = eliminate_redundant_dict_values(gene_to_uid)
    return gene_to_uid

def lookupOntologyID(geneset_type,ontology_name,type='name'):
    if geneset_type == 'GeneOntology': geneset_type = 'go'
    filename = 'AltDatabase/goelite/OBO/builds/'+geneset_type+'_annotations.txt'

    fn=filepath(filename)
    ### Imports a geneset category and stores pathway-level names
    i=0
    for line in open(fn,'rU').xreadlines():
        if i==0: i=1 ### Skip the header
        else:
            #data = cleanUpLine(line)
            data = line.strip()
            t = string.split(data,'\t')
            geneset_category = t[1]
            ontologyID = t[0]
            if type=='ID':
                if ontologyID == ontology_name:
                    ontology_id = geneset_category; break
            elif geneset_category == ontology_name:
                ontology_id = t[0]; break
    return ontology_id

def swapKeyValues(db):
    swapped={}
    for key in db:
        values = list(db[key]) ###If the value is not a list, make a list
        for value in values:
            try: swapped[value].append(key)
            except KeyError: swapped[value] = [key]
    swapped = eliminate_redundant_dict_values(swapped)
    return swapped

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
        #data = cleanUpLine(line)
        data = line.strip()
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
    def Pathway(self):
        try: return str(self.pathway)
        except Exception: return 'invalid'
    def setGraphID(self,graphID): self.graphID = graphID
    def setGroupID(self,groupid): self.groupid = groupid
    def GraphID(self): return self.graphID
    def GroupID(self): return self.groupid
    def setLabel(self,label): self.label = label
    def Label(self): return self.label
    def setModID(self,mod_list): self.mod_list = mod_list
    def ModID(self): return self.mod_list
    def Report(self):
        try: output = self.GeneID()+'|'+self.System()
        except Exception: print self.Label()
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
        try: export_data.write(values+'\n')
        except Exception: pass
    export_data.close()
    #print 'WikiPathways data exported to:',export_dir
    
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
    global species_code
    try: species_names = UI.remoteSpeciesInfo('yes') ### GO-Elite version
    except Exception: species_names = UI.remoteSpeciesAlt() ### AltAnalyze version
    for species_code in species_names:
        if all_species == 'yes' or species_code in specific_species:
            species_name = species_names[species_code].SpeciesName()
            species = string.replace(species_name,' ','_')
            
            ### Clear and create the output dir
            try: export.deleteFolder('GPML')
            except Exception: null=[]
            try: os.mkdir(filepath('GPML'))
            except Exception: null=[]
            
            ### Download all species GPML from .zip
            url = 'http://wikipathways.org//wpi/cache/wikipathways_'+species+'_Curation-AnalysisCollection__gpml.zip'
            fln,status = update.download(url,'GPML/','')
            
            if 'Internet' not in status:
                if len(specific_species) == 1:
                    print 'Including the latest WikiPathways associations'
                gene_data,pathway_db = parseGPML('/GPML') ### Get pathway associations (gene and annotation)
                clusterPrimaryPathwayGeneSystems(species_code,pathway_db)
                exportWikiPathwayData(species_name,pathway_db,'native')
  
                ### Create a "MAPPED" version (mapped will contain BAD MAPPINGS provide by the source database)!!!
                ensembl_to_WP = unifyGeneSystems(gene_data,species_code,'Ensembl') ### convert primary systems to MOD IDs                
                hmdb_to_WP = unifyGeneSystems(gene_data,species_code,'HMDB') ### convert primary systems to MOD IDs
                convertBetweenSystems(species_code,pathway_db,ensembl_to_WP,hmdb_to_WP) ### convert MOD IDs to all related primary systems (opposite as last method)
                exportWikiPathwayData(species_name,pathway_db,'mapped')
                if len(pathway_db)>0: ### Export all pathway interactions
                    try: exportNodeInteractions(pathway_db,'Ensembl','GPML')
                    except Exception: null=[]

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
        output = self.WPID()
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
        pathway_gene_data=[]
        complexes_data={}
        edge_data=[]
        #graph_node_data=[]
        pathway_type = 'GPML'
        xml=filepath(xml)
        filename = string.split(xml,'/')[-1]
        try: wpid = string.split(filename,'_')[-2]
        except Exception: wpid = filename[:-5]
        revision = string.split(filename,'_')[-1][:-5]
        try: dom = parse(xml)
        except Except: continue
        tags = dom.getElementsByTagName('Xref')
        data_node_tags = dom.getElementsByTagName('DataNode')
        groups = dom.getElementsByTagName('Group') ### complexes
        edges = dom.getElementsByTagName('Point') ### interacting nodes/complexes
        pathway_tag = dom.getElementsByTagName('Pathway')
        comment_tag = dom.getElementsByTagName('Comment')
        #print comment_tag.nodeValue; sys.exit()
        for pn in pathway_tag:
            pathway_name = pn.getAttribute("Name")
            organism = pn.getAttribute("Organism")
            data_source = pn.getAttribute("Data-Source")
        for ed in edges:
            ### Store internal graph data for pathway edges to build gene interaction networks later
            graphid = ed.getAttribute("GraphRef")
            edge_type = ed.getAttribute("ArrowHead")
            if edge_type == '': edge_pair = [graphid] ### either just a graphical line or the begining of a node-node edge
            else:
                try:
                    edge_pair.append(graphid)
                    edd = EdgeData(str(edge_pair[0]),str(edge_pair[1]),str(edge_type))
                    edge_data.append(edd)
                except Exception:
                    None ### Can happen with some pathways
        for gd in groups:
            ### Specific to groups and complexes
            groupID = gd.getAttribute("GroupId") ### Group node ID (same as "GroupRef")
            graphID = gd.getAttribute("GraphId") ### GPML node ID
            complexes_data[str(groupID)] = str(graphID)
        for cm in comment_tag:
            pathway_type = cm.getAttribute("Source")
            ### This is an issue with KEGG pathways who's names are truncated
            try:
                extended_comment_text = cm.childNodes[0].nodeValue
                if 'truncated' in extended_comment_text:
                    pathway_name = string.split(extended_comment_text,': ')[-1]
            except IndexError: null=[]
        if 'Kegg' in pathway_type:
            if '?' in data_source:
                pathway_id = string.split(data_source,'?')[-1]
                pathway_name += ':KEGG-'+pathway_id
                if 'WP' not in wpid:
                    wpid = pathway_id
        for i in data_node_tags:
            #print i.toxml()
            id = ''
            system_name = ''
            for x in i.childNodes: ### DataNode is the parent attribute with children GraphId, TextLabel, Xref
                if x.nodeName == 'Xref': ### Since the attributes we want are children of these nodes, must find the parents first
                    system_name = x.getAttribute("Database") ### System Code
                    id = x.getAttribute("ID") ### Gene or metabolite ID
            label = i.getAttribute("TextLabel") ### Gene or metabolite label
            type = i.getAttribute('Type') #E.g.', GeneProduct, Metabolite
            graphID = i.getAttribute("GraphId") ### WikiPathways graph ID
            groupID = i.getAttribute('GroupRef')### Group node ID
            #graph_node_data.append([graphID,groupID,label,type])
            gi = GeneIDInfo(str(system_name),str(id),pathway_name)
            gi.setGroupID(str(groupID)) ### Include internal graph IDs for determining edges
            gi.setGraphID(graphID)
            gi.setLabel(label)
            if len(id)>0 or 'Tissue' in pathway_name: ### Applies to the Lineage Profiler pathway which doesn't have IDs
                gene_data.append(gi)
                pathway_gene_data.append(gi)
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

def getGPMLGraphData(custom_sets_folder,species_code,mod):
    """ Calls methods to import GPML data and retrieve MOD IDs (Ensembl/Entrez) for GPML graphIDs """
    gpml_data,pathway_db = parseGPML(custom_sets_folder)
    gene_to_WP = unifyGeneSystems(gpml_data,species_code,mod)
    return pathway_db ### This object contains all pathways and all ID associations
    
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
                #test   
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
    if 'Symbol' not in systems_db: systems_db['Symbol']=1
    st = systemTranslation();
    for source_data in systems_db:
        original_source = source_data
        if 'Ensembl' in source_data: source_data = 'Ensembl'
        if source_data in st: source_data = st[source_data] ### convert the name to the GO-Elite compatible name
        if source_data != mod:
            mod_source = mod+'-'+source_data+'.txt'
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
    mapped=0
    for gi in xml_data:
        if gi.System() in system_ids:
            source_to_gene = system_ids[gi.System()]
            geneID = gi.GeneID()
            if gi.System() == 'Symbol':
                geneID = string.lower(geneID)
            if geneID in source_to_gene:
                for mod_id in source_to_gene[geneID]:
                    mapped+=1
                    try: mod_pathway[mod_id].append(gi.Pathway())
                    except Exception: mod_pathway[mod_id] = [gi.Pathway()]
                gi.setModID(source_to_gene[geneID]) ### Update this object to include associated MOD IDs (e.g., Ensembl or Entrez)
            else:
                source_to_gene = system_ids['Symbol'] ### Assume the missing ID is a symbol
                geneID = string.lower(geneID)
                if geneID in source_to_gene:
                    for mod_id in source_to_gene[geneID]:
                        mapped+=1
                        try: mod_pathway[mod_id].append(gi.Pathway())
                        except Exception: mod_pathway[mod_id] = [gi.Pathway()]
                    gi.setModID(source_to_gene[geneID]) ### Update this object to include associated MOD IDs (e.g., Ensembl or Entrez)
    #print mapped;sys.exit()
         
    ### If the MOD gene IDs are in the pathway then add these
    for gi in xml_data:
        source_data = gi.System() 
        if 'Ensembl' in source_data: source_data = 'Ensembl'
        if source_data in st: source_data = st[source_data] ### convert the name to the GO-Elite compatible name
        if source_data == mod:
            try: mod_pathway[gi.GeneID()].append(gi.Pathway())
            except Exception: mod_pathway[gi.GeneID()] = [gi.Pathway()]
            gi.setModID([gi.GeneID()])
        
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

def IDconverter(filename,species_code,input_system_name, output_system_name,analysis=None):
    """ This is a function built to convert the IDs in an input file from one system to another while preserving the original members """
    
    if 'HMDB' in input_system_name or 'HMDB' in output_system_name: mod = 'HMDB'
    elif 'ChEBI' in input_system_name or 'ChEBI' in output_system_name: mod = 'HMDB'
    elif 'KeggCompound' in input_system_name or 'KeggCompound' in output_system_name: mod = 'HMDB'
    elif 'CAS' in input_system_name or 'CAS' in output_system_name: mod = 'HMDB'
    elif 'PubChem' in input_system_name or 'PubChem' in output_system_name: mod = 'HMDB'
    else: mod = 'Ensembl'
        
    print 'Attempting to convert IDs from ', export.findFilename(filename)
    
    input_data_db, headers = importGenericDB(filename)
    if input_system_name == mod: ### This is or MOD
        source1_to_gene={}
        for id in input_data_db:
            source1_to_gene[id] = [id] ### make the primary ensembl ID = ensembl ID
    else:
        gene_to_source1 = getGeneToUid(species_code,mod+'-'+input_system_name)
        source1_to_gene = OBO_import.swapKeyValues(gene_to_source1)
    
    if output_system_name == mod: ### This is or MOD
        gene_to_source2={}
        gene_annotations = importGeneData(species_code,mod)
        for gene in gene_annotations:
            if 'LRG_' not in gene: ### Bad gene IDs from Ensembl
                gene_to_source2[gene] = [gene] ###import and use all availalable Ensembl IDs
    else:
        gene_to_source2 = getGeneToUid(species_code,mod+'-'+output_system_name)
        
    converted=0
    converted_ids={}
    for id in input_data_db:
        secondary_ids = ''
        genes = ''
        if id in source1_to_gene:
            genes = source1_to_gene[id]
            for gene in genes:
                if gene in gene_to_source2:
                    secondary_ids = gene_to_source2[gene]
                    secondary_ids = string.join(secondary_ids,'|')
            genes = string.join(genes,'|')
            converted+=1
        converted_ids[id] = secondary_ids, genes
    
    if analysis != 'signature':
        if '.txt' in filename:
            filename = string.replace(filename,'.txt','-'+output_system_name+'.txt')
        else:
            filename = filename[:-4]+'-'+output_system_name+'.txt'
        export_data = export.ExportFile(filename)
        headers = string.join([output_system_name,mod+' IDs']+headers,'\t')+'\n'
        export_data.write(headers)
        for id in input_data_db:
            secondary_ids, genes = converted_ids[id]
            for t in input_data_db[id]:
                export_values = string.join([secondary_ids,genes]+[id]+t,'\t')+'\n'
                export_data.write(export_values)
        export_data.close()
        
        print ''
        print converted, 'input',input_system_name,'IDs converted to',output_system_name,'out of',len(input_data_db)
        filename = export.findFilename(filename)
        return filename
    else:
        #print len(input_data_db), len(converted_ids)
        return converted_ids, input_data_db

if __name__ == '__main__':
    species_code = 'Hs'; mod = 'Ensembl'; gotype='nested'
    filedir = 'C:/Documents and Settings/Nathan Salomonis/My Documents/GO-Elite_120beta/Databases/EnsMart56Plus/Ce/gene/EntrezGene.txt'
    system = 'Macaroni'
    gene_annotations = importGeneData('Hs','EntrezGene')
    for i in gene_annotations:
        print i, gene_annotations[i].Symbol(); break
    print len(gene_annotations)
    sys.exit()
    import GO_Elite
    system_codes,source_types,mod_types = GO_Elite.getSourceData()
    #custom_sets_folder = '/test'
    #importGeneCustomData(species_code,system_codes,custom_sets_folder,mod); sys.exit()
    
    ### Test interactions export    
    custom_sets_folder = 'GPML'
    species_code = 'Hs'; mod = 'Ensembl'

    gene_to_symbol_db = getGeneToUid(species_code,('hide',mod+'-Symbol.txt')); #print mod_source, 'relationships imported.'
                            
    gpml_data,pathway_db = parseGPML(custom_sets_folder)
    gene_to_WP = unifyGeneSystems(gpml_data,species_code,mod)
    exportNodeInteractions(pathway_db,mod,custom_sets_folder)
    sys.exit()
    biopax_data = parseBioPax('/test'); sys.exit()
    

    
    gene_data,pathway_db = parseGPML('/test')
    print len(gpml_data)
    mod = 'Ensembl'
    gene_to_BioPax = unifyGeneSystems(biopax_data,species_code,mod)
    #gene_to_WP = unifyGeneSystems(gpml_data,species_code,mod)
    #gene_to_BioPax = combineDBs(gene_to_WP,gene_to_BioPax)
    for i in gene_to_BioPax:
        print i, gene_to_BioPax[i]
    print len(gene_to_BioPax); sys.exit()

    addNewCustomSystem(filedir,system,'yes','Ms'); kill
    addNewCustomRelationships('test.txt','Ensembl-MAPP','update','Mm');kill
    importGeneToOntologyData(species_code,mod,gotype,ontology_type);kill
    
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
        

#!/usr/bin/python
###########################
#Program:	GO-elite.py
#Author:	Nathan Salomonis
#Date:		12/12/06
#Website:	http://www.genmapp.org
#Email:	nsalomonis@gmail.com
###########################


        
    
