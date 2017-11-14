###mappfinder
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

"""This module contains methods for performing over-representation analysis (ORA) on input gene
lists provided by the user relative to denominator gene lists for nested Gene Ontology relationships
and WikiPathway biological pathways. These methods include a permutation based analysis and multiple
hypthesis correction."""

import sys, string
import os.path, platform
import unique
import math
import time
import gene_associations; reload(gene_associations)
try: from import_scripts import OBO_import
except Exception: pass
import GO_Elite
from stats_scripts import statistics
import random
import UI
import export; reload(export)
import re
from stats_scripts import fishers_exact_test
import traceback
import warnings

try:
    from scipy import stats
except Exception:
    pass ### scipy is not required but is used as a faster implementation of Fisher Exact Test when present

################# Parse directory files
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    try:
        dir_list = unique.read_directory(sub_dir)
    except Exception:
        dir_list=[] ### Directory does not exist
    dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

def readDirText(sub_dir):
    if sub_dir == None or sub_dir == '':
        dir_list2 = [] ### For automatically assigned denominators
    elif '.txt' in sub_dir:
        dir_list2 = [export.findFilename(sub_dir)] ### for pooled analyses - analyze only the specific file direction provided
    else:
        dir_list = unique.read_directory(sub_dir); dir_list2 = []
        ###Code to prevent folder names from being included
        for entry in dir_list:
            if entry[-4:] == ".txt" and '._' not in entry: dir_list2.append(entry)
    return dir_list2

###### Classes ######
class GrabFiles:
    def setdirectory(self,value):
        if '.txt' in value:
            value = export.findParentDir(value) ### For pooled analyses where only one file is submitted
        self.data = value
    def display(self):
        print self.data
    def directory(self):
        return self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        all_matching,file_dir,file = gene_associations.getDirectoryFiles(self.data,str(search_term))
        #if len(file)<1: print search_term,'not found'
        return file_dir,file
    def getAllFiles(self,search_term):
        #self is an instance while self.data is the value of the instance
        all_matching,file_dir,file = gene_associations.getDirectoryFiles(self.data,str(search_term))
        #if len(file)<1: print search_term,'not found'
        return all_matching
    
################# Import and Annotate Data
def eliminate_redundant_dict_values(database):
    db1={}
    for key in database: list = unique.unique(database[key]); list.sort(); db1[key] = list
    return db1

def swapKeyValues(db):
    swapped={}
    for key in db:
        values = db[key]
        for value in values:
            try: swapped[value].append(key)
            except KeyError: swapped[value] = [key]
    return swapped

def identifyGeneFiles(import_dir,gene_file):
    split_name = string.split(gene_file,'.')
    e = GrabFiles(); e.setdirectory(import_dir)
    dir_files = read_directory(import_dir)
    if len(split_name)>2:
        prefix_id = split_name[0]+'.'
        denominator_file_dir,denominator_file = e.searchdirectory(prefix_id)
    else: denominator_file_dir =''
    if len(dir_files)==1 or denominator_file_dir=='':
        try: denominator_file_dir,denominator_file = e.searchdirectory(dir_files[0])
        except IndexError:
            print_out = "WARNING: No denominator file included in\nthe GeneQuery/DenominatorGenes directory.\nTo proceed, place all denominator\nIDs in a file in that directory."
            if PoolVar: q.put([print_out]); return None
            ForceCriticalError(print_out)
            
    return denominator_file_dir

def associateInputSourceWithGene(source_to_gene,source_id_list):
    gene_db={}; count_null = 0
    for source_id in source_id_list:
        try:
            gene_ids = source_to_gene[source_id]
            for gene_id in gene_ids:
                try: gene_db[gene_id].append(source_id)
                except KeyError: gene_db[gene_id] = [source_id]
        except KeyError: count_null+=1
    #print count_null, 'source IDs not imported'
    return gene_db

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def importVersionData(dir):
    global OBO_date
    program_type,database_dir = unique.whatProgramIsThis(); parent_dir = ''
    if program_type == 'AltAnalyze': parent_dir = 'AltDatabase/goelite/'
    dir = parent_dir+dir
    
    filename = dir+'version.txt'; fn=filepath(filename)
    for line in open(fn,'r').readlines():
        data = cleanUpLine(line)
        OBO_version, OBO_date = string.split(data,'\t')
    return OBO_date

def checkDenominatorMatchesInput(input_gene_list,denominator_gene_list,gene_file):
    for id in input_gene_list:
        try: null = denominator_gene_list[id] ###this object was changed from a list to a dictionary for efficiency
        except KeyError: ###Only occurs if an input ID is NOT found in the denominator
            all_alphanumeric = string.join(re.findall(r"\w",id))
            print_out = 'Identifier: '+'"'+id+'"'+' not found in Denominator set '+str(len(id))+' '+all_alphanumeric+' '+str(len(input_gene_list))+' '+str(len(denominator_gene_list))+'\n'
            print_out = 'WARNING!!! Job stopped... Denominator gene list\ndoes not match the input gene list for\n%s' % gene_file
            if PoolVar: q.put([print_out]); return None
            ForceCriticalError(print_out)

def formatTime(start_time,end_time):
    intgr,decim = string.split(str(end_time-start_time),'.')
    ### Alternatively, use - round(end_time-start_time,1)
    return intgr+'.'+decim[0]
  
def importGeneData():
    ### Simple dictionary for verifying ID sets (important when wrong species selected)
    program_type,database_dir = unique.whatProgramIsThis()
    gene_import_dir = database_dir+'/'+species_code+'/gene/'+mod+'.txt'
    fn=filepath(gene_import_dir); gene_annotations={}; x=0
    for line in open(fn,'rU').xreadlines():             
        data = line.strip()
        t = string.split(data,'\t')
        if x == 0: x = 1
        else:
            gene_annotations[t[0]] = t[0]
    return gene_annotations
            
def checkCorrectSystemAndSpecies(input_gene_list):
    ### Verify that the input IDs are the correct system and species
    gene_annotations = importGeneData()

    missing=True
    while missing:
        for i in input_gene_list:
            try:
                null=gene_annotations[i]
                missing=False
            except Exception: pass
        break
    return missing

def generateMAPPFinderScores(species_title,species_id,source,mod_db,system_Codes,permute,resources_to_analyze,file_dirs,parent_root,poolVar=False,Q=None,Multi=None):
    global mappfinder_output_dir; global custom_sets_folder; global root; root = parent_root
    global mapp_to_mod_genes; global ontology_to_mod_genes; global system_codes; system_codes = system_Codes
    global q; q=Q; global mlp; mlp = Multi; global display; poolVar = True
    
    criterion_input_folder, criterion_denom_folder, output_dir, custom_sets_folder = file_dirs
    previous_denominator_file_dir = ''
    ontology_to_mod_genes={}; mapp_to_mod_genes={}; global test; test = 'no'; global PoolVar; PoolVar=poolVar
    program_type,database_dir = unique.whatProgramIsThis()
    if resources_to_analyze == 'Gene Ontology': resources_to_analyze = 'GeneOntology'

    if len(output_dir) == 0: mappfinder_output_dir = 'input/MAPPFinder'
    else: mappfinder_output_dir = output_dir + '/GO-Elite_results/CompleteResults/ORA'
    denom_folder = criterion_denom_folder
    
    global source_data; source_data = source; global mod; mod = mod_db
    global species_code; species_code = species_id
    global species_name; species_name = species_title; global gene_to_mapp
    global permutations; permutations = permute
    global eliminate_redundant_genes; eliminate_redundant_genes = 'yes'
    global permuted_z_scores; global ontology_annotations
    global original_ontology_z_score_data; global original_mapp_z_score_data
    global input_gene_list; global denominator_gene_list
    global gene_file; global denom_file_status
    global input_count; global denom_count; global gene_annotations
    global source_to_gene; global use_FET
    if permutations == "FisherExactTest":
        use_FET = 'yes' ### Use Fisher's Exact test instead of permutation-based p-values
        permutations = 0
    else:
        use_FET = 'no'
    if poolVar: display = False
    else: display = True
    start_time = time.time()   
    
    gene_annotations = gene_associations.importGeneData(species_code,mod)

    OBO_date = importVersionData('OBO/')

    if len(criterion_input_folder) == 0: import_dir = '/input/GenesToQuery/'+species_code
    else: import_dir = criterion_input_folder
    m = GrabFiles(); m.setdirectory(import_dir)
    import_dir_alt = m.directory()
    try: dir_list = readDirText(import_dir)  #send a sub_directory to a function to identify all files in a directory
    except Exception:
        print_out = 'Warning! Input directory location is not a valid folder. Exiting GO-Elite.'
        try:
            if PoolVar: q.put([print_out]); return None
        except Exception: pass
        ForceCriticalError(print_out)
    try: denom_dir_list = readDirText(criterion_denom_folder)
    except Exception:
        print_out = 'Warning! Denominator directory location is not a valid folder. Exiting GO-Elite.'
        try:
            if PoolVar: q.put([print_out]); return None
        except Exception: pass
        ForceCriticalError(print_out)
    if len(dir_list)==0:
        error_message = 'No files with the extension ".txt" found in the input directory.'
        try:
            if PoolVar: q.put([print_out]); return None
        except Exception: pass
        ForceCriticalError(error_message)
    if len(denom_dir_list)==0:
        skipDenomImport=True #### No denominator supplied by the user
        #error_message = 'No files with the extension ".txt" found in the denominator directory.'
        #if PoolVar: q.put([print_out]); return None
        #ForceCriticalError(error_message)
    else:
        skipDenomImport=False
    inputs_analyzed=0
    for mappfinder_input in dir_list:    #loop through each file in the directory
        permuted_z_scores={}; original_ontology_z_score_data={}; original_mapp_z_score_data={}
        if PoolVar==False:
            print 'Performing over-representation analysis (ORA) on',mappfinder_input
        gene_file_dir, gene_file = m.searchdirectory(mappfinder_input)
        ###Import Input gene/source-id lists
        input_gene_list,source_data_input,error_message = gene_associations.importUIDsForMAPPFinderQuery(import_dir_alt+'/'+gene_file,system_codes,'no'); input_count = len(input_gene_list)
        if 'No results' in error_message:
            continue
        if 'WARNING!!!' in error_message: ### Warn the user about SwissProt issues when importing the denominator
            try:
                if PoolVar: q.put([print_out]); return None
            except Exception: pass
            ForceCriticalError(error_message)
            continue
        if skipDenomImport==False: 
            denominator_file_dir = identifyGeneFiles(denom_folder,gene_file) ###input is in input\Genes, denominator in
            try:
                denominator_file_dir = identifyGeneFiles(denom_folder,gene_file) ###input is in input\Genes, denominator in
                denominator_file = string.split(denominator_file_dir,'/')[-1]
                if PoolVar==False:
                    print 'Using:', denominator_file,'for the denominator.'
            except Exception:
                print_out = "WARNING: No denominator file included in\nthe Denominator directory.\nTo proceed, place all denominator\nIDs in a file in that directory."
                try:
                    if PoolVar: q.put([print_out]); return None
                except Exception: pass
                ForceCriticalError(print_out)
        else:
            denominator_file_dir = None
        if denominator_file_dir == previous_denominator_file_dir: denom_file_status = 'old'
        else: denom_file_status = 'new'
        if denom_file_status == 'new':
            if skipDenomImport==False:
                previous_denominator_file_dir = denominator_file_dir
                denominator_gene_list,source_data_denom,error_message = gene_associations.importUIDsForMAPPFinderQuery(denominator_file_dir,system_codes,'no')
                denom_count = len(denominator_gene_list)
                if 'SwissProt' in error_message and 'WARNING!!!' not in error_message:
                    if len(input_gene_list)==0:
                        error_message+='\nNo valid input IDs found. Exiting GO-Elite.'
                        if PoolVar: q.put([print_out]); return None
                        ForceCriticalError(error_message)
                    else:
                        if PoolVar: q.put([print_out]); return None
                        ForceCriticalError(error_message)
                elif len(error_message)>0:
                    if PoolVar: q.put([print_out]); return None
                    ForceCriticalError(error_message)
                if len(denominator_gene_list) == len(input_gene_list):
                    print_out = 'Input and Denominator lists have identical counts.\nPlease load a propper denominator set (containing\nthe input list with all assayed gene IDs) before proceeding.'
                    if PoolVar: q.put([print_out]); return None
                    ForceCriticalError(print_out)
                original_denominator_gene_list=[]
                for id in denominator_gene_list: original_denominator_gene_list.append(id) ###need this to be a valid list not dictionary for permutation analysis
            else:
                ### Pull in all MOD IDs as a surrogate denominator
                denominator_gene_list = importGeneData()
                original_denominator_gene_list=[]
                for id in denominator_gene_list: original_denominator_gene_list.append(id)
                denom_count = len(denominator_gene_list)
        if len(source_data_input)>0: source_data = source_data_input ###over-ride source_data if a source was identified from the input file
        if source_data != mod:
            if denom_file_status == 'new':
                mod_source = mod+'-'+source_data+'.txt'
                try:
                    gene_to_source_id = gene_associations.getGeneToUid(species_code,mod_source,display=display)
                    if PoolVar==False:
                        print mod_source, 'imported'
                except Exception:
                    try:
                        if mod=='EntrezGene': mod = 'Ensembl'
                        else: mod = 'EntrezGene'
                        if PoolVar==False:
                            print 'The primary system (MOD) has been switched from',mod_db,'to',mod,'\n('+mod_db,'not supported for the %s ID system).' % source_data
                        mod_source = mod+'-'+source_data+'.txt'
                        gene_to_source_id = gene_associations.getGeneToUid(species_code,mod_source,display=display)
                    except Exception:
                        print_out = "WARNING: The primary gene ID system '"+mod+"'\ndoes not support relationships with '"+ source_data +"'.\nRe-run using a supported primary ID system."
                        try:
                            if PoolVar: q.put([print_out]); return None
                        except Exception: pass
                        ForceCriticalError(print_out)
                source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
                if skipDenomImport==False:
                    denominator_gene_list = associateInputSourceWithGene(source_to_gene,denominator_gene_list)
                    ### Introduced the below method in version 1.2.1 to improve permutation speed (no longer need to search all source IDs)
                    ### Only includes source ID to gene relationships represented in the denominator file (needed for Affymetrix)
                    source_to_gene = OBO_import.swapKeyValues(denominator_gene_list)
            ###Replace input lists with corresponding MOD IDs
            input_gene_list = associateInputSourceWithGene(source_to_gene,input_gene_list)
        else:
            if len(input_gene_list)>3:
                missing = checkCorrectSystemAndSpecies(input_gene_list)
                if missing:
                    print_out = "None of the input IDs match the selected mod (%s) or species (%s)." % (mod,species_name)
                    print_out += '\nVerify the correct ID system is indicated in the input file and the correct species seleced.'
                    try:
                        if PoolVar: q.put([print_out]); return None
                    except Exception: pass
                    ForceCriticalError(print_out)
        if skipDenomImport==False:
            checkDenominatorMatchesInput(input_gene_list,denominator_gene_list,gene_file) ###This is for only the associated MOD IDs

        gd = GrabFiles(); gd.setdirectory('/'+database_dir+'/'+species_code+'/gene-mapp')
        available_genesets = reorganizeResourceList(gd.getAllFiles(mod))
        od = GrabFiles(); od.setdirectory('/'+database_dir+'/'+species_code+'/gene-go')
        available_ontologies = reorganizeResourceList(od.getAllFiles(mod))
        
        input_gene_count = len(input_gene_list) ###Count number of genes associated with source input IDs

        if len(input_gene_list)==0 or len(denominator_gene_list)==0:
            if len(input_gene_list)==0:
                try:
                    print_out = 'WARNING!!!! None of the input IDs provided map to genes for '+mappfinder_input+'. Check to make sure the selected species is correct.'
                    print_out += '\nSelected species: '+species_name
                    print_out += '\nInput ID system: '+str(source_data_input)
                    print_out += '\nPrimary ID system (MOD): '+str(mod)
                    if PoolVar: q.put([print_out]); return None 
                except Exception:
                    pass
                ForceCriticalError(print_out)
            if len(denominator_gene_list)==0:
                try:
                    print_out = 'WARNING!!!! None of the denominator IDs provided map to genes for '+denominator_file_dir+'. Check to make sure the selected species is correct.'
                    print_out += '\nSelected species: '+species_name
                    print_out += '\nDenominator ID system: '+str(source)
                    print_out += '\nPrimary ID system (MOD):'+str(mod)
                    if PoolVar: q.put([print_out]); return None
                    ForceCriticalError(print_out)   
                except Exception:
                    pass

        elif len(available_ontologies) == 0 and len(available_genesets) == 0:
            print_out = 'WARNING!!!! No Ontology or GeneSets appear to be available for this species. Please supply and re-analyze.'
            try:
                if PoolVar: q.put([print_out]); return None
            except Exception: pass
            ForceCriticalError(print_out)
        else:
            """ Perform permutation analysis and ORA on available GeneSets or Ontologies"""
            inputs_analyzed+=1
            
            global permute_inputs; permute_inputs=[]
            if permutations != 0 or use_FET == 'no':
                buildPermutationDatabase(original_denominator_gene_list,input_count)
            
            run_status = 0
            ### Analyzed ontologies
            if len(available_ontologies)>0:
                if PoolVar==False:
                    print '    Analyzing input ID list with available ontologies'

            for ontology_dir in available_ontologies:
                ontology_type = getResourceType(ontology_dir)
                permuted_z_scores={}; original_ontology_z_score_data={}
                #print ontology_type, resources_to_analyze
                if (resources_to_analyze == ontology_type) or (resources_to_analyze == 'all') or (resources_to_analyze == 'both' and ontology_type == 'GeneOntology'):
                    ontology_annotations = importOntologyAnnotations(species_code,ontology_type)
                    if ontology_annotations!=None: ### Occurs when the files are named or formatted correctly
                        status, ontology_to_mod_genes = performOntologyORA(ontology_dir)
                        run_status += status
                    
            ### Analyzed gene-sets
            if len(available_genesets)>0:
                if PoolVar==False:
                    print '    Analyzing input ID list with available gene-sets'
            for geneset_dir in available_genesets:
                geneset_type = getResourceType(geneset_dir)
                permuted_z_scores={}; original_mapp_z_score_data={}
                #print geneset_type, resources_to_analyze;sys.exit()
                if geneset_type in resources_to_analyze or resources_to_analyze == 'all' or (resources_to_analyze == 'both' and geneset_type == 'Pathways'):
                    status, mapp_to_mod_genes = performGeneSetORA(geneset_dir)
                    run_status += status
            if len(custom_sets_folder)>0:
                ### Hence - Analyze User Supplied GeneSets
                permuted_z_scores={}; original_mapp_z_score_data={}
                run_status += performGeneSetORA('UserSuppliedAssociations')[0]
            
            permute_inputs=[]; permute_mapp_inputs=[]
            ontology_input_gene_count=[]; mapp_input_gene_count=[]

            if run_status == 0:
                ### Returns the number of successfully analyzed gene-set databases
                program_type,database_dir = unique.whatProgramIsThis()
                print_out = "Warning!!! Either the MOD you have selected: "+mod+"\nis missing the appropriate relationship files necessary to run GO-Elite\nor you have selected an invalid resource to analyze.  Either replace\nthe missing MOD files in "+database_dir+'/'+species_code+' sub-directories or\nselect a different MOD at run-time.'
                if PoolVar: q.put([print_out]); return None
                ForceCriticalError(print_out)
                
    end_time = time.time()
    time_diff = formatTime(start_time,end_time)
    if PoolVar==False:
        print 'ORA analyses finished in %s seconds' % time_diff
    try:
        q.put([ontology_to_mod_genes, mapp_to_mod_genes,time_diff,mappfinder_input,resources_to_analyze])
    except Exception:
        q = [ontology_to_mod_genes, mapp_to_mod_genes,time_diff,mappfinder_input,resources_to_analyze]
        return q ###Return the MOD genes associated with each GO term and MAPP

def importOntologyAnnotations(species_code,ontology_type):
    try:
        system_codes,source_types,mod_types = GO_Elite.getSourceData()
        verified_nested = OBO_import.verifyNestedFileCreation(species_code,mod_types,ontology_type)
        if verified_nested == 'no': force_error
        ontology_annotations = OBO_import.importPreviousOntologyAnnotations(ontology_type)
    except Exception:
        try:
            ### Occurs when the annotation file isn't built yet - if so try to build
            OBO_import.buildNestedOntologyAssociations(species_code,mod_types,ontology_type,display=display)
            ontology_annotations = OBO_import.importPreviousOntologyAnnotations(ontology_type)
        except Exception:
            ### Occurs for non-ontologies
            #print 'error 2',traceback.format_exc()
            ontology_annotations=None
    return ontology_annotations

def getResourceType(pathway_dir):
    pathway_type = string.split(pathway_dir,'-')[-1][:-4]
    if pathway_type == 'MAPP':
        pathway_type = 'Pathways'
    return pathway_type

def reorganizeResourceList(pathway_list):
    ### Make sure that WikiPathways and GO are analyzed last, so that gene results are also reported last to GO_Elite.py
    add_pathway=[]
    pathway_list_reorganized=[]
    for pathway in pathway_list:
        if '-MAPP.txt' in pathway: add_pathway.append(pathway)
        elif '-GeneOntology.txt' in pathway: add_pathway.append(pathway)
        else: pathway_list_reorganized.append(pathway)
    pathway_list_reorganized+=add_pathway
    return pathway_list_reorganized

def ForceCriticalError(print_out):
    print '\n'+print_out+'\n'
    """
    if len(sys.argv[1:])<2: ### Don't create a Tkinter window if command-line options supplied
        try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); #sys.exit()
        except Exception: None
    else:
        #sys.exit()
        pass
    """
    
def buildPermutationDatabase(original_denominator_gene_list,input_count):
    if PoolVar==False:
        print "Building %d permuted ID sets" % permutations,
    
    global k; k=0; x=0
    try: original_increment = int(permutations/10); increment = original_increment
    except Exception: null=None
    if PoolVar==False:
        if permutations!=0: print '*',
    start_time = time.time() ### Build Permutation Identifier Database
            
    while x<permutations:
        if x == increment and PoolVar==False:
            increment+=original_increment; print '*',
        try: permute_input_list = random.sample(original_denominator_gene_list,input_count); x+=1
        except ValueError:
            print_out = 'Input count>Denominator '+str(len(original_denominator_gene_list))+' '+str(input_count)+'\n terminating'
            if PoolVar: q.put([print_out]); return None
            ForceCriticalError(print_out)
        #permute_input_list = random.sample(denominator_gene_list,len(input_gene_list)); x+=1
        #permute_input_list = random.shuffle(original_denominator_gene_list); x+=1; permute_input_list = permute_input_list[:input_count]
        if source_data!=mod:  ###Store the randomly choosen input lists for GenMAPP MAPP Permutation analysis
            permute_input_list = associateInputSourceWithGene(source_to_gene,permute_input_list)
            if len(permute_input_list)>len(input_gene_list): k+=1
        permute_inputs.append(permute_input_list)

    end_time = time.time()
    time_diff = formatTime(start_time,end_time)
    if PoolVar==False:
        print 'completed in %s seconds' % time_diff
            
def swapKeyValues(db):
    rev_db = {}
    for i, k_list in db.iteritems():
        for k in k_list:
            try: rev_db[k].append(i)
            except Exception: rev_db[k] = [i]
    return rev_db

def performGeneSetORA(geneset_dir):
    """ Perform over-representation analysis (ORA) on any provided Gene Set """
    
    start_time = time.time()
    geneset_type = getResourceType(geneset_dir)
    #permuted_z_scores={}; original_mapp_z_score_data={}
    
    if geneset_type == 'Pathways': geneset_type = 'WikiPathways'
    ### Since MAPP tables can be provided by the user, allow the file to be missing
    if geneset_dir == 'UserSuppliedAssociations':
        gene_to_mapp = gene_associations.importGeneCustomData(species_code,system_codes,custom_sets_folder,mod)
        geneset_type = geneset_dir
    else:
        try: gene_to_mapp = gene_associations.importGeneMAPPData(species_code,geneset_dir)
        except Exception: gene_to_mapp = {}

    mapp_to_gene = swapKeyValues(gene_to_mapp)
    
    if len(gene_to_mapp)==0:
        return 0, None
    else:
        
        ###Calculate primary z-scores for GeneSets
        #print len(input_gene_list), len(gene_to_mapp)
        mapp_to_mod_genes = getGenesInPathway(input_gene_list,gene_to_mapp) ### For summary reporting
        mapp_input_gene_count,Rm,input_linked_mapp = countGenesInPathway(input_gene_list,gene_to_mapp,'yes')

        mapp_denominator_gene_count,Nm,denom_linked_mapp = countGenesInPathway(denominator_gene_list,gene_to_mapp,'yes')
        #print Nm,"unique genes, linked to GeneSets and in dataset and", Rm, "unique GeneSets\n linked genes matching criterion."
        
        #zstart_time = time.time()
        
        try:
            #exhaustiveMultiZScores(mapp_input_gene_count,mapp_denominator_gene_count,Nm,Rm,mapp_to_gene,'MAPP')
            if PoolVar==False:
                try: multiZScores(mapp_input_gene_count,mapp_denominator_gene_count,Nm,Rm,mapp_to_gene,'MAPP')
                except Exception: calculateZScores(mapp_input_gene_count,mapp_denominator_gene_count,Nm,Rm,mapp_to_gene,'MAPP')
            else:
                calculateZScores(mapp_input_gene_count,mapp_denominator_gene_count,Nm,Rm,mapp_to_gene,'MAPP')
        except Exception:
            calculateZScores(mapp_input_gene_count,mapp_denominator_gene_count,Nm,Rm,mapp_to_gene,'MAPP')
        
        #print len(permuted_z_scores), len(original_mapp_z_score_data)
        #print time.time()-zstart_time, 'seconds...';sys.exit()
        if use_FET == 'no':
            permute_mapp_inputs=[]
            ###Begin GeneSets Permutation Analysis
            try: original_increment = int(permutations/10); increment = original_increment
            except Exception: null=None
            x=0
            if PoolVar==False:
                if permutations!=0: print '*',
            for permute_input_list in permute_inputs:
                if PoolVar==False:
                    if x == increment: increment+=original_increment; print '*',
                x+=1
                permute_mapp_input_gene_count,null,null = countGenesInPathway(permute_input_list,gene_to_mapp,'no')
                permute_mapp_inputs.append(permute_mapp_input_gene_count)
            calculatePermuteZScores(permute_mapp_inputs,mapp_denominator_gene_count,Nm,Rm)
            calculatePermuteStats(original_mapp_z_score_data)
        adjustPermuteStats(original_mapp_z_score_data)
    
        mapp_headers = formatHeaders(gene_file,input_count,input_linked_mapp,denom_count,denom_linked_mapp,Rm,Nm,'MAPP',OBO_date)
        exportPathwayData(original_mapp_z_score_data,gene_file,mapp_headers,geneset_type,'local')
                                
        ### Export all gene associations (added in version 1.21)
        exportPathwayToGeneAssociations(mapp_to_mod_genes,mod,gene_file,gene_annotations,geneset_type,'local')
    
        end_time = time.time()
        time_diff = formatTime(start_time,end_time)
        if PoolVar==False:
            print "Initial results for %s calculated in %s seconds" % (geneset_type,time_diff)
        permute_mapp_inputs=[]
        
        return 1, mapp_to_mod_genes
    
def performOntologyORA(ontology_dir):
    """ Perform over-representation analysis (ORA) on any provided Ontology """
    
    start_time = time.time()
    ontology_type = getResourceType(ontology_dir)
    
    ######### Import Gene-to-Nested-Ontology #########
    gene_to_ontology = gene_associations.importGeneToOntologyData(species_code,mod,'nested',ontology_type)
    ontology_to_gene = OBO_import.swapKeyValues(gene_to_ontology)   
    if len(gene_to_ontology)==0:
        return 0, None
    else:           
        ######### Calculate primary z-scores for GO terms
        #a = time.time()
        ontology_to_mod_genes = getGenesInPathway(input_gene_list,gene_to_ontology) ### For summary gene reporting
        #b = time.time(); print 'a',b-a
        ontology_input_gene_count,Rg,input_linked_ontology = countGenesInPathway(input_gene_list,gene_to_ontology,'yes')
        #c = time.time(); print 'b',c-b
        ontology_denominator_gene_count,Ng,denom_linked_ontology = countGenesInPathway(denominator_gene_list,gene_to_ontology,'yes')
        #d = time.time(); print 'c',d-c
        #print Ng,"unique genes, linked to GO and in dataset and", Rg, "unique GO linked genes matching criterion."
        try:
            if PoolVar==False:
                multiZScores(ontology_input_gene_count,ontology_denominator_gene_count,Ng,Rg,ontology_to_gene,'Ontology')
            else:
                calculateZScores(ontology_input_gene_count,ontology_denominator_gene_count,Ng,Rg,ontology_to_gene,'Ontology')
        except Exception: calculateZScores(ontology_input_gene_count,ontology_denominator_gene_count,Ng,Rg,ontology_to_gene,'Ontology')
        #e = time.time(); print 'd',e-d; sys.exit()
        if use_FET == 'no':
            ###Begining Ontology Permutation Analysis
            try: original_increment = int(permutations/10); increment = original_increment
            except Exception: null=None
            x=0
            permute_ontology_inputs=[]
            if PoolVar==False:
                if permutations!=0: print '*',
            for permute_input_list in permute_inputs:
                ### http://docs.python.org/library/multiprocessing.html
                if PoolVar==False:
                    if x == increment: increment+=original_increment; print '*',
                x+=1
                permute_ontology_input_gene_count,null,null = countGenesInPathway(permute_input_list,gene_to_ontology,'no'); permute_input_list=[]
                permute_ontology_inputs.append(permute_ontology_input_gene_count)
            if PoolVar==False:
                if permutations !=0: print 'Ontology finished'
            calculatePermuteZScores(permute_ontology_inputs,ontology_denominator_gene_count,Ng,Rg)
            calculatePermuteStats(original_ontology_z_score_data)
        adjustPermuteStats(original_ontology_z_score_data)
        go_headers = formatHeaders(gene_file,input_count,input_linked_ontology,denom_count,denom_linked_ontology,Rg,Ng,'Ontology',OBO_date)
        exportPathwayData(original_ontology_z_score_data,gene_file,go_headers,ontology_type,'Ontology')

        ### Export all gene associations (added in version 1.21)          
        exportPathwayToGeneAssociations(ontology_to_mod_genes,mod,gene_file,gene_annotations,ontology_type,'Ontology')
        end_time = time.time()
        time_diff = formatTime(start_time,end_time)
        if PoolVar==False:
            print "Initial results for %s calculated in %s seconds" % (ontology_type,time_diff)
        permute_ontology_inputs=[]
        return 1, ontology_to_mod_genes
            
def exportPathwayToGeneAssociations(pathway_to_mod_genes,mod,gene_file,gene_annotations,resource_name,pathway_type):
    headers = string.join([mod,'symbol',resource_name],'\t')+'\n'
    if resource_name == 'GeneOntology': resource_name = 'GO' ### Makes the output filename compatible with GenMAPP-CS plugin filenames
    if resource_name == 'WikiPathways': resource_name = 'local' ### Makes the output filename compatible with GenMAPP-CS plugin filenames
    new_file = mappfinder_output_dir+'/'+gene_file[:-4]+'-'+resource_name+'-associations.tab'
    data = export.ExportFile(new_file); data.write(headers)
    
    for pathway in pathway_to_mod_genes:
        for gene in pathway_to_mod_genes[pathway]:
            try: symbol = gene_annotations[gene].Symbol()
            except Exception: symbol = ''
            if pathway_type == 'Ontology' and ':' not in pathway: pathway = 'GO:'+ pathway
            values = string.join([gene,symbol,pathway],'\t')+'\n'
            data.write(values)
    data.close()
    
def formatHeaders(gene_file,input_count,input_linked,denom_count,denom_linked,R,N,pathway_type,OBO_date):
    headers = []
    headers.append('GO-Elite ORA Results')
    headers.append('File:')
    headers.append('Table:')
    if pathway_type == 'Ontology':
        headers.append('Database: Based on OBO-Database version: '+OBO_date)
    headers.append('colors:')
    t = time.localtime(); dt = str(t[1])+'/'+str(t[2])+'/'+str(t[0])
    headers.append(dt)
    headers.append(species_name)
    headers.append('Pvalues = true')
    headers.append('Calculation Summary:')
    headers.append(str(input_count)+' '+source_data+' source identifiers supplied in the input file:'+gene_file)
    headers.append(str(input_linked)+' source identifiers meeting the filter linked to a '+mod+' ID.')
    headers.append(str(R)+' genes meeting the criterion linked to a term.')
    headers.append(str(denom_count)+' source identifiers in this dataset.')
    headers.append(str(denom_linked)+' source identifiers linked to a '+mod+' ID.')
    headers.append(str(N)+' Genes linked to a term.')
    headers.append('The z score is based on an N of '+str(N)+' and a R of '+str(R)+' distinct genes in all terms.\n')
    
    if use_FET == 'yes': prob = "FisherExactP"
    else: prob = "PermuteP"
    if pathway_type == 'Ontology':
        title = ['Ontology-ID','Ontology Name','Ontology Type','Number Changed','Number Measured','Number in Ontology','Percent Changed','Percent Present','Z Score',prob,'AdjustedP']
        title = string.join(title,'\t'); headers.append(title)
    else:
        title = ['Gene-Set Name','Number Changed','Number Measured','Number in Gene-Set','Percent Changed','Percent Present','Z Score',prob,'AdjustedP']
        title = string.join(title,'\t'); headers.append(title)
    header_str = string.join(headers,'\n')
    return header_str+'\n'
                                                               
def exportPathwayData(original_pathway_z_score_data,gene_file,headers,resource_name,pathway_type,altOutputDir=None):
    if resource_name == 'GeneOntology': resource_name = 'GO' ### Makes the output filename compatible with GenMAPP-CS plugin filenames
    if resource_name == 'WikiPathways': resource_name = 'local' ### Makes the output filename compatible with GenMAPP-CS plugin filenames
    new_file = mappfinder_output_dir+'/'+gene_file[:-4]+'-'+resource_name+'.txt'
    
    global sort_results
    data = export.ExportFile(new_file); data.write(headers); sort_results=[]
    #print "Results for",len(original_pathway_z_score_data),"pathways exported to",new_file
    for pathway in original_pathway_z_score_data:
        zsd=original_pathway_z_score_data[pathway]
        try: results = [zsd.Changed(), zsd.Measured(), zsd.InPathway(), zsd.PercentChanged(), zsd.PercentPresent(), zsd.ZScore(), zsd.PermuteP(), zsd.AdjP()]
        except AttributeError:
            return traceback.format_exc()
            #print pathway,len(permuted_z_scores[pathway]);kill
        try: ###This is unnecessary, unless using the non-nested GO associations (which can have out of sync GOIDs)
            if pathway_type == 'Ontology':
                s = ontology_annotations[pathway]
                annotations = [s.OntologyID(),s.OntologyTerm(),s.OntologyType()]; results = annotations + results 
            else:
                results = [pathway] + results
            results = string.join(results,'\t') + '\n'
            sort_results.append([float(zsd.ZScore()),-1/float(zsd.Measured()),results])
        except KeyError: null = []
    sort_results.sort(); sort_results.reverse()
    for values in sort_results:
        results = values[2]
        data.write(results)
    data.close()
    
def swapKeyValuesTuple(db):
    swapped={}
    for key in db:
        values = tuple(db[key]) ###If the value is not a list, make a list
        swapped[values] = [key]
    swapped = eliminate_redundant_dict_values(swapped)
    return swapped

class ZScoreData:
    def __init__(self,pathway,changed,measured,zscore,null_z,in_pathway):
        self._pathway = pathway; self._changed = changed; self._measured = measured
        self._zscore = zscore; self._null_z = null_z; self._in_pathway = in_pathway
    def PathwayID(self): return self._pathway
    def Changed(self): return str(int(self._changed))
    def Measured(self): return str(int(self._measured))
    def InPathway(self): return str(self._in_pathway)
    def ZScore(self): return str(self._zscore)
    def SetP(self,p): self._permute_p = p
    def PermuteP(self): return str(self._permute_p)
    def SetAdjP(self,adjp): self._adj_p = adjp
    def AdjP(self): return str(self._adj_p)
    def setAssociatedIDs(self,ids): self.ids = ids
    def AssociatedIDs(self): return self.ids
    def PercentChanged(self):
        try: pc = float(self.Changed())/float(self.Measured())*100
        except Exception: pc = 0
        return str(pc)
    def PercentPresent(self):
        try: pp = float(self.Measured())/float(self.InPathway())*100
        except Exception: pp = 0
        return str(pp)
    def NullZ(self): return self._null_z
    def Report(self):
        output = self.PathwayID()
        return output
    def __repr__(self): return self.Report()
    
def workerExhaustive(queue,pathway,genes_in_pathway,pathway_input_gene_count,pathway_denominator_gene_count,N,R,pathway_type):
    permuted_z_scores_instance={}; original_mapp_z_score_data_instance={}; original_ontology_z_score_data_instance={}
    """ Exhaustive multiprocessing execution """
    try:
        n = pathway_denominator_gene_count[pathway]
        try: r = pathway_input_gene_count[pathway]
        except Exception: r = 0.0000
    except Exception: n = 0.0000; r = 0.0000
    if n != 0:
        try: z = Zscore(r,n,N,R)
        except ZeroDivisionError: z = 0.0000
        try: null_z = Zscore(0,n,N,R)
        except ZeroDivisionError: null_z = 0.000
        zsd = ZScoreData(pathway,r,n,z,null_z,genes_in_pathway)
        if pathway_type == 'Ontology': original_ontology_z_score_data_instance[pathway] = zsd
        else: original_mapp_z_score_data_instance[pathway] = zsd
        permuted_z_scores_instance[pathway] = [z]
        #if '06878' in pathway: print pathway, z, null_z, r,n, N, R;kill
        if use_FET == 'yes':
            ### Alternatively calculate p using the Fisher's Exact Test
            p = FishersExactTest(r,n,R,N)
            zsd.SetP(p)
    queue.put((permuted_z_scores_instance,original_mapp_z_score_data_instance,original_ontology_z_score_data_instance)) 

def exhaustiveMultiZScores(pathway_input_gene_count,pathway_denominator_gene_count,N,R,pathway_db,pathway_type):
    """ Exhaustive multiprocessing - create a process for each entry in the dictionary for Z-score calcualtion """
    procs=list()
    
    queue = mlp.Queue()
    for pathway in pathway_db:
        genes_in_pathway = len(pathway_db[pathway])
        p = mlp.Process(target=workerExhaustive, args=(queue,pathway,genes_in_pathway,pathway_input_gene_count,pathway_denominator_gene_count,N,R,pathway_type))
        procs.append(p)
        p.start()
    
    permuted_z_scores_list=[]; original_mapp_z_score_data_list=[]; original_ontology_z_score_data_list=[]
    for _ in procs:
        val = queue.get()
        permuted_z_scores_list.append(val[0]); original_mapp_z_score_data_list.append(val[1])
        original_ontology_z_score_data_list.append(val[2])
        
    for p in procs:
        p.join()
     
    for i in permuted_z_scores_list: permuted_z_scores.update(i)
    for i in original_mapp_z_score_data_list: original_mapp_z_score_data.update(i)
    for i in original_ontology_z_score_data_list: original_ontology_z_score_data.update(i)
    
def multiZScores(pathway_input_gene_count,pathway_denominator_gene_count,N,R,pathway_db,pathway_type):
    """ Create a finate pool of processes (4) for Z-score calculation """
    if mlp.cpu_count() < 3:
        processors = mlp.cpu_count()
    else: processors = 4
    pool = mlp.Pool(processes=processors)
    si = (len(pathway_db)/processors)
    s = si; b=0
    db_ls=[]
    if len(pathway_db)<10: forceError ### will si to be zero and an infanite loop
    while s<len(pathway_db):
        db_ls.append(dict(pathway_db.items()[b:s]))
        b+=si; s+=si
    db_ls.append(dict(pathway_db.items()[b:s]))

    ### Create an instance of MultiZscoreWorker (store the variables to save memory)
    workerMulti = MultiZscoreWorker(pathway_input_gene_count,pathway_denominator_gene_count,N,R,pathway_type,use_FET)
    
    results = pool.map(workerMulti,db_ls)
    pool.close(); pool.join(); pool = None
    permuted_z_scores_list=[]; original_mapp_z_score_data_list=[]; original_ontology_z_score_data_list=[]
    for (a,b,c) in results:
        permuted_z_scores_list.append(a); original_mapp_z_score_data_list.append(b)
        original_ontology_z_score_data_list.append(c)
        
    for i in permuted_z_scores_list: permuted_z_scores.update(i)
    for i in original_mapp_z_score_data_list: original_mapp_z_score_data.update(i)
    for i in original_ontology_z_score_data_list: original_ontology_z_score_data.update(i)

class MultiZscoreWorker:
    def __init__(self,pathway_input_gene_count,pathway_denominator_gene_count,N,R,pathway_type,use_FET):
        self.pathway_input_gene_count = pathway_input_gene_count
        self.pathway_denominator_gene_count = pathway_denominator_gene_count
        self.N = N
        self.R = R
        self.pathway_type = pathway_type
        self.use_FET = use_FET
        
    def __call__(self,pathway_db):
        N = self.N; R = self.R; use_FET = self.use_FET
        zt=0; nzt=0; zsdt=0; rt=0; ft=0

        permuted_z_scores_instance={}; original_mapp_z_score_data_instance={}; original_ontology_z_score_data_instance={}
        for pathway in pathway_db:
            genes_in_pathway = len(pathway_db[pathway])
            try:
                n = self.pathway_denominator_gene_count[pathway]
                #r1 = time.time()
                try: r = self.pathway_input_gene_count[pathway]
                except Exception: r = 0.0000
                #rt += time.time() - r1
            except Exception: n = 0.0000; r = 0.0000
            if n != 0:
                z1 = time.time()
                try: z = Zscore(r,n,N,R)
                except ZeroDivisionError: z = 0.0000
                #zt += time.time() - z1
                #nz1 = time.time()
                try: null_z = Zscore(0,n,N,R)
                except ZeroDivisionError: null_z = 0.000
                #nzt+= time.time() - nz1
                #zsd1 = time.time()
                zsd = ZScoreData(pathway,r,n,z,null_z,genes_in_pathway)
                #zsdt+= time.time() - zsd1
                if self.pathway_type == 'Ontology': original_ontology_z_score_data_instance[pathway] = zsd
                else: original_mapp_z_score_data_instance[pathway] = zsd
                permuted_z_scores_instance[pathway] = [z]
                #if '06878' in pathway: print pathway, z, null_z, r,n, N, R;kill
                if use_FET == 'yes':
                    ### Alternatively calculate p using the Fisher's Exact Test
                    #ft1 = time.time()
                    p = FishersExactTest(r,n,R,N)
                    #ft+= time.time() - ft1
                    zsd.SetP(p)
        #print zt,nzt,zsdt,rt,ft ### Used for efficiency evaluation
        return permuted_z_scores_instance,original_mapp_z_score_data_instance, original_ontology_z_score_data_instance

def calculateZScores(pathway_input_gene_count,pathway_denominator_gene_count,N,R,pathway_db,pathway_type):
    """where N is the total number of genes measured: 
    R is the total number of genes meeting the criterion:
    n is the total number of genes in this specific MAPP: 
    r is the number of genes meeting the criterion in this MAPP: """
    for pathway in pathway_db:
        try:
            n = pathway_denominator_gene_count[pathway]
            try: r = pathway_input_gene_count[pathway]
            except Exception: r = 0.0000
        except Exception: n = 0.0000; r = 0.0000
        if n != 0:
            try: z = Zscore(r,n,N,R)
            except ZeroDivisionError: z = 0.0000
            try: null_z = Zscore(0,n,N,R)
            except ZeroDivisionError: null_z = 0.000
            genes_in_pathway = len(pathway_db[pathway])
            zsd = ZScoreData(pathway,r,n,z,null_z,genes_in_pathway)
            if pathway_type == 'Ontology': original_ontology_z_score_data[pathway] = zsd
            else: original_mapp_z_score_data[pathway] = zsd
            permuted_z_scores[pathway] = [z]
            #if '06878' in pathway: print pathway, z, null_z, r,n, N, R;kill
            if use_FET == 'yes':
                ### Alternatively calculate p using the Fisher's Exact Test
                p = FishersExactTest(r,n,R,N)
                zsd.SetP(p)

def Zscore(r,n,N,R):
    N=float(N) ### This bring all other values into float space
    z = (r - n*(R/N))/math.sqrt(n*(R/N)*(1-(R/N))*(1-((n-1)/(N-1))))
    return z

def calculatePermuteZScores(permute_pathway_inputs,pathway_denominator_gene_count,N,R):
    for pathway_input_gene_count in permute_pathway_inputs:  
        for pathway in pathway_input_gene_count:
            r = pathway_input_gene_count[pathway]
            n = pathway_denominator_gene_count[pathway]
            try: z = statistics.zscore(r,n,N,R)
            except ZeroDivisionError: z = 0
            permuted_z_scores[pathway].append(abs(z))
            #if pathway == '0005488':
            #a.append(r)
                
def calculatePermuteStats(original_pathway_z_score_data):
    for pathway in original_pathway_z_score_data:
        zsd = original_pathway_z_score_data[pathway]
        z = abs(permuted_z_scores[pathway][0])
        permute_scores = permuted_z_scores[pathway][1:] ###Exclude the true value
        nullz = zsd.NullZ()
        if abs(nullz) == z: ###Only add the nullz values if they can count towards the p-value (if equal to the original z)
            null_z_to_add = permutations - len(permute_scores)
            permute_scores+=[abs(nullz)]*null_z_to_add ###Add null_z's in proportion to the amount of times there were not genes found for that pathway
        if len(permute_scores)>0: p = permute_p(permute_scores,z)  
        else: p = 0
        #if p>1: p=1
        zsd.SetP(p)

def adjustPermuteStats(original_pathway_z_score_data):
 
    #1. Sort ascending the original input p value vector.  Call this spval.  Keep the original indecies so you can sort back.
    #2. Define a new vector called tmp.  tmp= spval.  tmp will contain the BH p values.
    #3. m is the length of tmp (also spval)
    #4. i=m-1
    #5  tmp[ i ]=min(tmp[i+1], min((m/i)*spval[ i ],1)) - second to last, last, last/second to last
    #6. i=m-2
    #7  tmp[ i ]=min(tmp[i+1], min((m/i)*spval[ i ],1))
    #8  repeat step 7 for m-3, m-4,... until i=1
    #9. sort tmp back to the original order of the input p values.
    
    global spval; spval=[]; adj_p_list=[]
    for pathway in original_pathway_z_score_data:
        zsd = original_pathway_z_score_data[pathway]
        p = float(zsd.PermuteP())
        spval.append([p,pathway])
        
    spval.sort(); tmp = spval; m = len(spval); i=m-2; x=0 ###Step 1-4
    l=0
    
    while i > -1:
        adjp = min(tmp[i+1][0], min((float(m)/(i+1))*spval[i][0],1))
        tmp[i]=adjp,tmp[i][1]; i -= 1
        if adjp !=0: adj_p_list.append(adjp) ### get the minimum adjp
        
    for (adjp,pathway) in tmp:
        try:
            if adjp == 0: adjp = min(adj_p_list)
        except Exception: null=[]
        zsd = original_pathway_z_score_data[pathway]
        zsd.SetAdjP(adjp)
   
def adjustPermuteStatsTemp(pval_db):
    
    global spval; spval=[]
    for element in pval_db:
        zsd = pval_db[element]
        try:
            try: p = float(zsd.PermuteP())
            except AttributeError: p = float(zsd[0]) ### When values are indeces rather than objects
        except Exception: p = 1
        spval.append([p,element])

    spval.sort(); tmp = spval; m = len(spval); i=m-2; x=0 ###Step 1-4       
    #spval.sort(); tmp = spval; m = len(spval)-1; i=m-1; x=0 ###Step 1-4
    
    while i > -1:
        tmp[i]=min(tmp[i+1][0], min((float(m)/(i+1))*spval[i][0],1)),tmp[i][1]; i -= 1
        
    for (adjp,element) in tmp:
        zsd = pval_db[element]
        try: zsd.SetAdjP(adjp)
        except AttributeError: zsd[1] = adjp ### When values are indeces rather than objects
     
def permute_p(null_list,true_value):
    y = 0; z = 0; x = permutations
    for value in null_list:
        if value >= true_value: y += 1
    #if true_value > 8: global a; a = null_list; print true_value,y,x;kill
    return (float(y)/float(x))  ###Multiply probabilty x2?

def FishersExactTest(r,n,R,N):
    """
    N is the total number of genes measured (Ensembl linked from denom) (total number of ) (number of exons evaluated)
    R is the total number of genes meeting the criterion (Ensembl linked from input) (number of exonic/intronic regions overlaping with any CLIP peeks)
    n is the total number of genes in this specific MAPP (Ensembl denom in MAPP) (number of exonic/intronic regions associated with the SF)
    r is the number of genes meeting the criterion in this MAPP (Ensembl input in MAPP) (number of exonic/intronic regions with peeks overlapping with the SF)
    
    With these values, we must create a 2x2 contingency table for a Fisher's Exact Test
    that reports:
     
    +---+---+    a is the # of IDs in the term regulated
    | a | b |    b is the # of IDs in the term not-regulated 
    +---+---+    c is the # of IDs not-in-term and regulated
    | c | d |    d is the # of IDs not-in-term and not-regulated
    +---+---+
    
    If we know r=20, R=80, n=437 and N=14480
    +----+-----+    
    | 20 | 417 |  437   
    +----+-----+    
    | 65 |13978| 14043
    +----+-----+
      85  14395  14480
    """

    a = r; b = n-r; c=R-r; d=N-R-b
    table = [[int(a),int(b)], [int(c),int(d)]]

    """
    print a,b; print c,d
    from stats_scripts import fishers_exact_test; table = [[a,b], [c,d]]
    ft = fishers_exact_test.FishersExactTest(table)
    print ft.probability_of_table(table); print ft.two_tail_p()
    print ft.right_tail_p(); print ft.left_tail_p()
    """
    
    try: ### Scipy version - cuts down rutime by ~1/3rd the time
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides import warnings
            oddsratio, pvalue = stats.fisher_exact(table)
        return pvalue
    except Exception:
        ft = fishers_exact_test.FishersExactTest(table)
        return ft.two_tail_p() 

def getGenesInPathway(gene_list,gene_to_pathway):
    ###This function is similar to countGenesInPathway, but is used to return the genes associated with a pathway
    ### Can be used to improve downstream annotation speed when this file is present rather than re-derive
    pathway_to_gene={}
    for gene in gene_list:
        if gene in gene_to_pathway:
            pathways = gene_to_pathway[gene]
            for pathway in pathways:
                try: pathway_to_gene[pathway].append(gene)
                except KeyError: pathway_to_gene[pathway] = [gene]
    return pathway_to_gene
            
def countGenesInPathway(gene_list,gene_to_pathway,count_linked_source):
    pathway_count={}; associated_genes={}; linked_source={}
    ### Add genes to a dictionary of pathways to get unique counts (could count directly, but biased by redundant source-id associations with MOD)
    for gene in gene_list:
        if source_data != mod and eliminate_redundant_genes == 'yes':
            gene_id = tuple(gene_list[gene]) ### switches gene with list of source_ids (if made unique, decreased redundant association)
            if count_linked_source == 'yes':
                for id in gene_id: linked_source[id] = []
        else: gene_id = gene; linked_source[gene_id] = []
        try:
            pathways = gene_to_pathway[gene]
            associated_genes[gene_id] = []
            for pathway in pathways:
                try: pathway_count[pathway][gene_id]=None ### more efficent storage for length determination below
                except: pathway_count[pathway] = {gene_id:None}
        except Exception: pass
        
    ### Count unique gene or source set associations per pathway
    unique_associated_gene_count = len(associated_genes)
    linked_count = len(linked_source)
    for pathway in pathway_count:
        pathway_count[pathway] = len(pathway_count[pathway])
    return pathway_count, unique_associated_gene_count, linked_count

if __name__ == '__main__':
    #r=20; R=85; n=437; N=14480
    species_name = 'Galus galus'; species_code = 'Gg'; source_data = 'EnsTranscript'; mod = 'Ensembl'
    species_name = 'Mus musculus'; species_code = 'Mm'; source_data = 'EntrezGene'; mod = 'EntrezGene'
    species_name = 'Homo sapiens'; species_code = 'Hs'; source_data = 'Ensembl'; mod = 'Ensembl'

    system_codes={}; system_codes['L'] = 'EntrezGene'; system_codes['En'] = 'Ensembl'; system_codes['X'] = 'Affymetrix'
    file_dirs = 'C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample/input_list_small','C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample/denominator','C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample'
    file_dirs = '/Users/nsalomonis/Desktop/GOElite-test/input','/Users/nsalomonis/Desktop/GOElite-test/denom','/Users/nsalomonis/Desktop/GOElite-test','/Users/nsalomonis/Desktop/GOElite-test/miR'
    permute = 20000
    permute = 'FisherExactTest'
    generateMAPPFinderScores(species_name,species_code,source_data,mod,system_codes,permute,'all',file_dirs,'')

#!/usr/bin/python
###########################
#Program:	GO-elite.py
#Author:	Nathan Salomonis
#Date:		12/12/06
#Website:	http://www.genmapp.org
#Email:	nsalomonis@gmail.com
###########################


        
    
