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
import OBO_import
import GO_Elite
import statistics
import random
import UI
import export; reload(export)
import re

################# Parse directory files
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

def readDirText(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt": dir_list2.append(entry)
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
        if len(file)<1: print search_term,'not found'
        return file_dir,file
    
def getDirectoryFiles(import_dir, search_term):
    exact_file = ''; exact_file_dir=''
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        if (':' in import_dir) or ('/Users/' == import_dir[:7]) or ('Linux' in platform.system()): affy_data_dir = import_dir+'/'+data
        else: affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir: exact_file_dir = affy_data_dir; exact_file = data
    return exact_file_dir,exact_file

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
            try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
            except Exception: print print_out; sys.exit()
            
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
            all_alphanumeric = re.findall(r"\w",id)
            print 'Identifier:','"'+id+'"', 'not found in Denominator set',len(id),all_alphanumeric,len(input_gene_list),len(denominator_gene_list)
            print_out = 'WARNING!!! Job stopped... Denominator gene list\ndoes not match the input gene list for\n%s' % gene_file
            try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
            except Exception:
                print '\nWARNING!!!! Job stopped... Denominator gene list does not match the input gene list for',gene_file
                print 'Please use an appropriate denominator list (e.g. all probe sets from the microarray experiment or all known species gene IDs).'; sys.exit()

def generateMAPPFinderScores(species_title,species_id,source,mod_db,system_codes,permute,resources_to_analyze,file_dirs,parent_root):
    global mappfinder_output_dir; global root; root = parent_root; resource_type = ''
    criterion_input_folder, criterion_denom_folder, output_dir, custom_sets_folder = file_dirs
    go_to_mod_genes={}; mapp_to_mod_genes={}
    if len(output_dir) == 0: mappfinder_output_dir = 'input/MAPPFinder'
    else: mappfinder_output_dir = output_dir + '/GO-Elite_results/CompleteResults/ORA'

    global source_data; source_data = source; global mod; mod = mod_db
    global species_code; species_code = species_id
    global species_name; species_name = species_title; global gene_to_mapp
    global permutations; permutations = permute; previous_denominator_file_dir = ''
    global eliminate_redundant_genes; eliminate_redundant_genes = 'yes'
    global permuted_z_scores; global go_annotations
    global original_go_z_score_data; global original_mapp_z_score_data 
    
    start_time = time.time()   
    go_annotations = OBO_import.importPreviousGOAnnotations()
    gene_annotations = gene_associations.importGeneData(species_code,mod)
    
    OBO_date = importVersionData('OBO/')
    if len(criterion_input_folder) == 0: import_dir = '/input/GenesToQuery/'+species_code; import_dir_alt = import_dir[1:]
    else: import_dir = criterion_input_folder; import_dir_alt = criterion_input_folder
    m = GrabFiles(); m.setdirectory(import_dir)
    dir_list = readDirText(import_dir)  #send a sub_directory to a function to identify all files in a directory
    denom_dir_list = readDirText(criterion_denom_folder)
    if len(dir_list)==0:
        error_message = 'No files with the extension ".txt" found in the input directory.'
        try: UI.WarningWindow(error_message,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
        except Exception: print error_message; sys.exit()
    if len(denom_dir_list)==0:
        error_message = 'No files with the extension ".txt" found in the denominator directory.'
        try: UI.WarningWindow(error_message,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
        except Exception: print error_message; sys.exit()
        
    inputs_analyzed=0
    for mappfinder_input in dir_list:    #loop through each file in the directory
        permuted_z_scores={}; original_go_z_score_data={}; original_mapp_z_score_data={}
        print 'Performing mappfinder analysis on',mappfinder_input
        gene_file_dir, gene_file = m.searchdirectory(mappfinder_input)
        ###Import Input gene/source-id lists
        input_gene_list,source_data_input,error_message = gene_associations.importUIDsForMAPPFinderQuery(import_dir_alt+'/'+gene_file,system_codes,'no'); input_count = len(input_gene_list)
        if len(error_message)>0:
            try: UI.WarningWindow(error_message,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
            except Exception: print error_message; sys.exit()
        if len(criterion_denom_folder)==0: denom_folder = '/input/GenesToQuery/'+species_code+'/DenominatorGenes'
        else: denom_folder = criterion_denom_folder
        error_warning = "\nThe directory\n"+'['+denom_folder+']'+"\nis not found. Please create the directory\nand place an appropriate denominator file\nor files in it."
        denominator_file_dir = identifyGeneFiles(denom_folder,gene_file) ###input is in input\Genes, denominator in
        try: denominator_file_dir = identifyGeneFiles(denom_folder,gene_file) ###input is in input\Genes, denominator in
        except Exception:
            print_out = "WARNING: No denominator file included in\nthe Denominator directory.\nTo proceed, place all denominator\nIDs in a file in that directory."
            try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
            except Exception: print error_warning; sys.exit()
        if denominator_file_dir == previous_denominator_file_dir: denom_file_status = 'old'
        else: denom_file_status = 'new'
        if denom_file_status == 'new':
            denominator_gene_list,source_data_denom,error_message = gene_associations.importUIDsForMAPPFinderQuery(denominator_file_dir,system_codes,'no'); denom_count = len(denominator_gene_list)
            if len(error_message)>0:
                try: UI.WarningWindow(error_message,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
                except Exception: print error_message; sys.exit()
            if len(denominator_gene_list) == len(input_gene_list):
                try: UI.WarningWindow('Input and Denominator lists have identical counts.\nPlease load a propper denominator set (containing\nthe input list with all assayed gene IDs) before proceeding.','Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
                except Exception: print print_out; sys.exit()    
            original_denominator_gene_list=[]
            for id in denominator_gene_list: original_denominator_gene_list.append(id) ###need this to be a valid list not dictionary for permutation analysis
        if len(source_data_input)>0: source_data = source_data_input ###over-ride source_data if a source was identified from the input file
        if source_data != mod:
            if denom_file_status == 'new':
                mod_source = mod+'-'+source_data
                #checkDenominatorMatchesInput(input_gene_list,denominator_gene_list,gene_file) ###This is checked for the source IDs not associated MOD IDs
                try: gene_to_source_id = gene_associations.getGeneToUid(species_code,mod_source)
                except Exception:
                    try:
                        if mod=='EntrezGene': mod = 'Ensembl'
                        else: mod = 'EntrezGene'
                        print 'The primary system (MOD) has been switched from',mod_db,'to',mod,'('+mod_db,'not supported for the input ID system).'
                        mod_source = mod+'-'+source_data
                        gene_to_source_id = gene_associations.getGeneToUid(species_code,mod_source)
                    except Exception:
                        print_out = "WARNING: The primary gene ID system '"+mod+"'\ndoes not support relationships with '"+ source_data +"'.\nRe-run using a supported primary ID system."
                        try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
                        except Exception: print print_out; sys.exit()
                source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
                denominator_gene_list = associateInputSourceWithGene(source_to_gene,denominator_gene_list)
                ### Introduced the below method in version 1.21 to improve permutation speed (no longer need to search all source IDs)
                ### Only includes source ID to gene relationships represented in the denominator file (needed for Affymetrix)
                source_to_gene = OBO_import.swapKeyValues(denominator_gene_list)
            ###Replace input lists with corresponding MOD IDs
            input_gene_list = associateInputSourceWithGene(source_to_gene,input_gene_list)
        checkDenominatorMatchesInput(input_gene_list,denominator_gene_list,gene_file) ###This is for only the associated MOD IDs

        if inputs_analyzed == 0:
            ### Import Gene-to-Nested-GO and Gene-to-MAPP Associations
            """Rather than perform these import steps at the begining of the function, we do this here
            in case the mod needs to be changed. This occurs when the user designated mod does not link
            to the source IDs determined by GO-Elite"""
            try:
                if resources_to_analyze != 'Pathways':
                    gene_to_go = gene_associations.importGeneGOData(species_code,mod,'nested')
                    go_to_gene = OBO_import.swapKeyValues(gene_to_go)
                    resource_type += 'Gene Ontology'
            except Exception:
                program_type,database_dir = unique.whatProgramIsThis()
                print_out = "Warning!!! The MOD you have selected: "+mod+"\nis missing the appropriate relationship\nfiles necessary to run GO-Elite.  Either\nreplace the missing MOD files in\n"+database_dir+'/'+species_code+' sub-directories or\nselect a different MOD at run-time.'
                try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
                except Exception: print print_out
                sys.exit()        
            ### Since MAPP tables can be provided by the user, allow the file to be missing
            if resources_to_analyze != 'GeneOntology':
                if len(custom_sets_folder)>0:
                    gene_to_mapp = gene_associations.importGeneCustomData(species_code,system_codes,custom_sets_folder,mod)
                else:
                    try: gene_to_mapp = gene_associations.importGeneMAPPData(species_code,mod)
                    except Exception: gene_to_mapp = {}
                mapp_to_gene = OBO_import.swapKeyValues(gene_to_mapp)
                if len(resource_type)==0: resource_type += 'Pathways'
                else: resource_type += '/Pathways'
                
        inputs_analyzed+=1
        input_gene_count = len(input_gene_list) ###Count number of genes associated with source input IDs
        if len(input_gene_list)>0 and len(denominator_gene_list)>0:
            ###Calculate primary z-scores for GO terms
            
            ### The values below are for gene summary reporting (only the last file will be summarized)
            if resources_to_analyze != 'Pathways': go_to_mod_genes = getGenesInPathway(input_gene_list,gene_to_go)
            if resources_to_analyze != 'GeneOntology': mapp_to_mod_genes = getGenesInPathway(input_gene_list,gene_to_mapp)
            
            if resources_to_analyze != 'Pathways':
                go_input_gene_count,Rg,input_linked_go = countGenesInPathway(input_gene_list,gene_to_go,'yes')
                if denom_file_status == 'new':
                    go_denominator_gene_count,Ng,denom_linked_go = countGenesInPathway(denominator_gene_list,gene_to_go,'yes')
                #print Ng,"unique genes, linked to GO and in dataset and", Rg, "unique GO linked genes matching criterion."
                calculateZScores(go_input_gene_count,go_denominator_gene_count,Ng,Rg,go_to_gene,'GO')
            if resources_to_analyze != 'GeneOntology':
                ###Calculate primary z-scores for GenMAPP MAPPs
                mapp_input_gene_count,Rm,input_linked_mapp = countGenesInPathway(input_gene_list,gene_to_mapp,'yes')
                mapp_denominator_gene_count,Nm,denom_linked_mapp = countGenesInPathway(denominator_gene_list,gene_to_mapp,'yes')
                #print Nm,"unique genes, linked to GenMAPP MAPPs and in dataset and", Rm, "unique GenMAPP\nMAPPs linked genes matching criterion."
                calculateZScores(mapp_input_gene_count,mapp_denominator_gene_count,Nm,Rm,mapp_to_gene,'MAPP')
            end_time = time.time(); time_diff = int(end_time-start_time)
            print "Primary "+resource_type+" Z-scores calculated in %d seconds." % time_diff
            print "Begining permutation analysis..."

            ###Begin GO Permute Analysis
            global k; k=0
            start_time = time.time(); #print 'Permuting the Gene Ontology analysis %d times' % permutations
            x=0; permute_inputs=[]; permute_go_inputs=[] #; denominator_gene_list = unique.unique_db(denominator_gene_list)
            try: original_increment = int(permutations/20); increment = original_increment
            except Exception: null=None
            if permutations!=0: print '*',
            while x<permutations:
                if x == increment: increment+=original_increment; print '*',
                try: permute_input_list = random.sample(original_denominator_gene_list,input_count); x+=1
                except ValueError: print 'Input count>Denominator',len(original_denominator_gene_list), input_count,'\n','terminating'; sys.stdin.readline(); sys.exit()
                #permute_input_list = random.sample(denominator_gene_list,len(input_gene_list)); x+=1
                #permute_input_list = random.shuffle(original_denominator_gene_list); x+=1; permute_input_list = permute_input_list[:input_count]
                if source_data!=mod:  ###Store the randomly choosen input lists for GenMAPP MAPP Permutation analysis
                    permute_input_list = associateInputSourceWithGene(source_to_gene,permute_input_list)
                    if len(permute_input_list)>len(input_gene_list): k+=1
                permute_inputs.append(permute_input_list) 
                if resources_to_analyze != 'Pathways':
                    permute_go_input_gene_count,null,null = countGenesInPathway(permute_input_list,gene_to_go,'no'); permute_input_list=[]
                    permute_go_inputs.append(permute_go_input_gene_count)
            if resources_to_analyze != 'Pathways':
                if permutations !=0: print 'Gene Ontology finished'
                calculatePermuteZScores(permute_go_inputs,go_denominator_gene_count,Ng,Rg)
                calculatePermuteStats(original_go_z_score_data)
                adjustPermuteStats(original_go_z_score_data)
                end_time = time.time(); time_diff = int(end_time-start_time)
                print "Permuted p-values for GO calculated in %d seconds" % time_diff
                go_headers = formatHeaders(gene_file,input_count,input_linked_go,denom_count,denom_linked_go,Rg,Ng,'GO',OBO_date)
                exportPathwayData(original_go_z_score_data,gene_file,go_headers,'GO')

            if resources_to_analyze != 'GeneOntology':
                start_time = time.time(); #print 'Permuting the GenMAPP MAPP analysis %d times' % permutations
                permute_mapp_inputs=[]
                ###Begin GenMAPP MAPP Permute Analysis
                try: original_increment = int(permutations/20); increment = original_increment
                except Exception: null=None
                x=0
                if permutations!=0: print '*',
                for permute_input_list in permute_inputs:
                    if x == increment: increment+=original_increment; print '*',
                    x+=1
                    permute_mapp_input_gene_count,null,null = countGenesInPathway(permute_input_list,gene_to_mapp,'no')
                    permute_mapp_inputs.append(permute_mapp_input_gene_count)
                calculatePermuteZScores(permute_mapp_inputs,mapp_denominator_gene_count,Nm,Rm)
                calculatePermuteStats(original_mapp_z_score_data)
                adjustPermuteStats(original_mapp_z_score_data)
                if permutations !=0: print 'Pathways finished'
                end_time = time.time(); time_diff = int(end_time-start_time)
                print "Permuted p-values for Pathways calculated in %d seconds" % time_diff
                mapp_headers = formatHeaders(gene_file,input_count,input_linked_mapp,denom_count,denom_linked_mapp,Rm,Nm,'MAPP',OBO_date)
                exportPathwayData(original_mapp_z_score_data,gene_file,mapp_headers,'local')
                
                previous_denominator_file_dir = denominator_file_dir
                permute_inputs=[]; permute_go_inputs=[]; permute_mapp_inputs=[]
                go_input_gene_count=[]; mapp_input_gene_count=[]
        else:
            if len(input_gene_list)==0:
                print_out = 'WARNING!!!! The length of the input file\n '+ mappfinder_input +' is zero.\nLength of input ID list is '+str(input_count)+' for system: '+str(source_data_input)+' and mod: '+str(mod)
                print_out += '. This warning is likely due to the selection of an incorrect species.\nCurrent species is '+ species_name
                try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
                except Exception: print print_out; sys.exit()   
            if len(denominator_gene_list)==0:
                print_out = 'WARNING!!!! The length of the denominator file\n '+denominator_file_dir+' is zero.\nLength of denominator ID list is '+str(denom_count)+' for system: ' + str(source_data_denom)+' '+str(source)+' '+str(mod)
                try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('yes'); sys.exit()
                except Exception: print print_out; sys.exit()
      
        ### Export all gene associations (added in version 1.21)          
        exportPathwayToGeneAssociations(go_to_mod_genes,mod,gene_file,gene_annotations,'GO')
        exportPathwayToGeneAssociations(mapp_to_mod_genes,mod,gene_file,gene_annotations,'local')
    
    return go_to_mod_genes, mapp_to_mod_genes ###Return the MOD genes associated with each GO term and MAPP

def exportPathwayToGeneAssociations(pathway_to_mod_genes,mod,gene_file,gene_annotations,pathway_type):
    new_file = mappfinder_output_dir+'/'+gene_file[:-4]+'-'+pathway_type+'-associations.tab'
    headers = string.join([mod,'symbol',pathway_type],'\t')+'\n'
    data = export.ExportFile(new_file); data.write(headers)
    
    for pathway in pathway_to_mod_genes:
        for gene in pathway_to_mod_genes[pathway]:
            try: symbol = gene_annotations[gene].Symbol()
            except Exception: symbol = ''
            if pathway_type == 'GO' and ':' not in pathway: pathway = 'GO:'+ pathway
            values = string.join([gene,symbol,pathway],'\t')+'\n'
            data.write(values)
    data.close()
    
def formatHeaders(gene_file,input_count,input_linked,denom_count,denom_linked,R,N,pathway_type,OBO_date):
    headers = []
    headers.append('GO-Elite MAPPFinder Results')
    headers.append('File:')
    headers.append('Table:')
    headers.append('Database: Based on OBO-Database version: '+OBO_date)
    headers.append('colors:')
    t = time.localtime(); dt = str(t[1])+'/'+str(t[2])+'/'+str(t[0])
    headers.append(dt)
    headers.append(species_name)
    headers.append('Pvalues = true')
    headers.append('Calculation Summary:')
    headers.append(str(input_count)+' '+source_data+' source identifiers supplied in the input file:'+gene_file)
    headers.append(str(input_linked)+' source identifiers meeting the filter linked to a '+mod+' ID.')
    headers.append(str(R)+' genes meeting the criterion linked to a GO term.')
    headers.append(str(denom_count)+' source identifiers in this dataset.')
    headers.append(str(denom_linked)+' source identifiers linked to a '+mod+' ID.')
    headers.append(str(N)+' Genes linked to a GO term.')
    headers.append('The z score is based on an N of '+str(N)+' and a R of '+str(R)+' distinct genes in the GO.\n')
    if pathway_type == 'GO':
        title = ['GOID','GO Name','GO Type','Number Changed Local','Number Measured Local','Number in GO Local','Percent Changed Local','Percent Present Local','Number Changed','Number Measured','Number in GO','Percent Changed','Percent Present','Z Score','PermuteP','AdjustedP']
        title = string.join(title,'\t'); headers.append(title)
    else:
        title = ['MAPP Name','Number Changed','Number Measured','Number on MAPP','Percent Changed','Percent Present','Z Score','PermuteP','AdjustedP']
        title = string.join(title,'\t'); headers.append(title)
    header_str = string.join(headers,'\n')
    return header_str+'\n'
                                                               
def exportPathwayData(original_pathway_z_score_data,gene_file,headers,pathway_type):
    new_file = mappfinder_output_dir+'/'+gene_file[:-4]+'-'+pathway_type+'.txt'
    
    global sort_results
    data = export.ExportFile(new_file); data.write(headers); sort_results=[]
    #print "Results for",len(original_pathway_z_score_data),"pathways exported to",new_file
    for pathway in original_pathway_z_score_data:
        zsd=original_pathway_z_score_data[pathway]
        try: results = [zsd.Changed(), zsd.Measured(), zsd.InPathway(), zsd.PercentChanged(), zsd.PercentPresent(), zsd.ZScore(), zsd.PermuteP(), zsd.AdjP()]
        except AttributeError: print pathway,len(permuted_z_scores[pathway]);kill
        try: ###This is unnecessary, unless using the non-nested GO associations (which can have out of sync GOIDs)
            if pathway_type == 'GO':
                s = go_annotations[pathway]
                annotations = [s.GOID(),s.GOName(),s.GOType()]; annotations += ['']*5; results = annotations + results 
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
    def Changed(self): return str(self._changed)
    def Measured(self): return str(self._measured)
    def InPathway(self): return str(self._in_pathway)
    def ZScore(self): return str(self._zscore)
    def SetP(self,p): self._permute_p = p
    def PermuteP(self): return str(self._permute_p)
    def SetAdjP(self,adjp): self._adj_p = adjp
    def AdjP(self): return str(self._adj_p)
    def PercentChanged(self):
        try: pc = float(self.Changed())/float(self.Measured())*100
        except ZeroDivisionError: pc = 0
        return str(pc)
    def PercentPresent(self):
        try: pp = float(self.Measured())/float(self.InPathway())*100
        except ZeroDivisionError: pp = 0
        return str(pp)
    def NullZ(self): return self._null_z
    def Report(self):
        output = self.PathwayID()
        return output
    def __repr__(self): return self.Report()
    
def calculateZScores(pathway_input_gene_count,pathway_denominator_gene_count,N,R,pathway_db,pathway_type):
    """where N is the total number of genes measured: 
    R is the total number of genes meeting the criterion:
    n is the total number of genes in this specific MAPP: 
    r is the number of genes meeting the criterion in this MAPP: """
    for pathway in pathway_db:
        if pathway in pathway_denominator_gene_count:
            n = pathway_denominator_gene_count[pathway]
            if pathway in pathway_input_gene_count: r = pathway_input_gene_count[pathway]
            else: r = 0
        else: n = 0; r = 0
        if n != 0:
            try: z = statistics.zscore(r,n,N,R)
            except ZeroDivisionError: z = 0
            try: null_z = statistics.zscore(0,n,N,R)
            except ZeroDivisionError: null_z = 0
            genes_in_pathway = len(pathway_db[pathway])
            zsd = ZScoreData(pathway,r,n,z,null_z,genes_in_pathway)
            if pathway_type == 'GO': original_go_z_score_data[pathway] = zsd
            else: original_mapp_z_score_data[pathway] = zsd
            permuted_z_scores[pathway] = [z]
            #if '06878' in pathway: print pathway, z, null_z, r,n, N, R;kill

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
        
def permute_p(null_list,true_value):
    y = 0; z = 0; x = permutations
    for value in null_list:
        if value >= true_value: y += 1
    #if true_value > 8: global a; a = null_list; print true_value,y,x;kill
    return (float(y)/float(x))  ###Multiply probabilty x2?

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
        if gene in gene_to_pathway:
            associated_genes[gene_id] = []
            pathways = gene_to_pathway[gene]
            for pathway in pathways:
                try: pathway_count[pathway].append(gene_id)
                except: pathway_count[pathway] = [gene_id]
    ### Count unique gene or source set associations per pathway
    #global pathway_redundancies; pathway_redundancies={}
    linked_count=0
    for pathway in pathway_count:
        unique_genes = unique.dictionary(pathway_count[pathway])  ###more efficient code for unique
        """###Code used to assess redundancy in probeset assignments
        if source_data != mod:
            redundant_source_ids={}
            for source_id_set in unique_genes:
                for source in source_id_set:
                    try: redundant_source_ids[source].append(source_id_set)
                    except KeyError: redundant_source_ids[source] = [source_id_set]
            for id in redundant_source_ids:
                redundant_source_ids_filtered={}
                if len(redundant_source_ids[id])>1:
                    redundant_source_ids_filtered[id] = redundant_source_ids[id]
                    try: pathway_redundancies[pathway].append([redundant_source_ids_filtered])
                    except KeyError: pathway_redundancies[pathway] = [redundant_source_ids_filtered]
        """
        count = len(unique_genes)
        pathway_count[pathway] = count
    unique_associated_gene_count = len(associated_genes)
    linked_count = len(linked_source)
    return pathway_count, unique_associated_gene_count, linked_count

if __name__ == '__main__':
    species_name = 'Galus galus'; species_code = 'Gg'; source_data = 'EnsTranscript'; mod = 'Ensembl'
    species_name = 'Homo sapiens'; species_code = 'Hs'; source_data = 'Ensembl'; mod = 'Ensembl'
    species_name = 'Mus musculus'; species_code = 'Mm'; source_data = 'EntrezGene'; mod = 'EntrezGene'
    system_codes={}; system_codes['L'] = 'EntrezGene'; system_codes['En'] = 'Ensembl'; system_codes['X'] = 'Affymetrix'
    file_dirs = 'C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample/input_list_small','C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample/denominator','C:/Documents and Settings/Nathan/Desktop/GenMAPP/Mm_sample'
    generateMAPPFinderScores(species_name,species_code,source_data,mod,system_codes,0,file_dirs,'')
    
#!/usr/bin/python
###########################
#Program:	GO-elite.py
#Author:	Nathan Salomonis
#Date:		12/12/06
#Website:	http://www.genmapp.org
#Email:	nsalomonis@gmail.com
###########################


        
    
