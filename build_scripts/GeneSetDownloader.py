###MetabolomicsParser
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

"""This module contains methods for reading the HMDB and storing relationships"""

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os.path
import unique
import export
import time
import update; reload(update)
from import_scripts import OBO_import
import gene_associations
import traceback

############# Common file handling routines ############# 
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def lowerSymbolDB(source_to_gene):
    source_to_gene2={}
    for symbol in source_to_gene:
        source_to_gene2[string.lower(symbol)]=source_to_gene[symbol]
    return source_to_gene2
    
def verifyFile(filename):
    fn=filepath(filename); file_found = 'yes'
    try:
        for line in open(fn,'rU').xreadlines():break
    except Exception: file_found = 'no'
    return file_found

def importSpeciesData():
    if program_type == 'GO-Elite': filename = 'Config/species_all.txt' ### species.txt can be cleared during updating
    else: filename = 'Config/goelite_species.txt'
    x=0
    fn=filepath(filename);global species_list; species_list=[]; global species_codes; species_codes={}
    global species_names; species_names={}
    global species_taxids; species_taxids={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t'); abrev=t[0]; species=t[1]
        try: taxid = t[2]
        except Exception: taxid = None
        if x==0: x=1
        else:
            species_list.append(species)
            species_codes[species] = abrev
            species_names[abrev] = species
            species_taxids[abrev] = taxid
           
def getSourceData():
    filename = 'Config/source_data.txt'; x=0
    fn=filepath(filename)
    global source_types; source_types={}
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
            source_types[source]=system_code
            system_codes[system_code] = source ###Used when users include system code data in their input file
            
############# File download/extraction #############           
def downloadPAZARAssocations():
    url = 'http://www.pazar.info/tftargets/tftargets.zip'
    print 'Downloading Transcription Factor to Target associations'
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/tftargets/','')
    return 'raw'

def downloadPAZARAssocationsOld():
    """ This works fine, but is redundant with the new zip file that contains all files"""
    
    base_url = 'http://www.pazar.info/tftargets/'
    filenames = getPAZARFileNames()
    print 'Downloading Transcription Factor to Target associations'
    source = 'raw'
    r = 4; k = -1
    
    for resource in filenames:
        filename = filenames[resource]
        url = base_url+filename
        start_time = time.time()
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/tftargets/','')
        end_time = time.time()
        if (end_time-start_time)>3: ### Hence the internet connection is very slow (will take forever to get everything)
            downloadPreCompiledPAZAR() ### Just get the compiled symbol data instead
            print '...access to source PAZAR files too slow, getting pre-compiled from genmapp.org'
            source = 'precompiled'
            break
        k+=1
        if r==k:
            k=0
            print '*',
    print ''
    return source

def downloadPreCompiledPAZAR():
    """ Downloads the already merged symbol to TF file from PAZAR files """
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/tf-target.txt'
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/tftargets/symbol/','')
     
def downloadAmadeusPredictions():
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/symbol-Metazoan-Amadeus.txt'
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Amadeus/','')
    
def downloadBioMarkers(output_dir):
    ensembl_version = unique.getCurrentGeneDatabaseVersion()
    #url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/Hs_exon_tissue-specific_protein_coding.zip'
    url = 'http://altanalyze.org/archiveDBs/AltDatabase/updated/'+ensembl_version+'/Hs_LineageProfiler.zip'
    print 'Downloading BioMarker associations'
    fln,status = update.downloadSuppressPrintOuts(url,output_dir,'')
    
    #url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/Mm_gene_tissue-specific_protein_coding.zip'
    url = 'http://altanalyze.org/archiveDBs/AltDatabase/updated/'+ensembl_version+'/Mm_LineageProfiler.zip'
    fln,status = update.downloadSuppressPrintOuts(url,output_dir,'')
    
def downloadKEGGPathways(species):
    print "Integrating KEGG associations for "+species
    url = 'http://www.genmapp.org/go_elite/Databases/KEGG/'+species+'-KEGG_20110518.zip'
    ### This is a fixed date resource since KEGG licensed their material after this date
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/KEGG/','')

def downloadDomainAssociations(selected_species):
    paths=[]
    if selected_species != None: ### Restrict to selected species only
        current_species_dirs=selected_species
    else:
        current_species_dirs = unique.read_directory('/'+database_dir)
    for species in current_species_dirs:
        url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/Domains/'+species+'_Ensembl-Domain.gz'
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Domains/','txt')
        if 'Internet' not in status:
            paths.append((species,fln))
    return paths

def downloadPhenotypeOntologyOBO():
    print 'Downloading Phenotype Ontology structure and associations'
    url = 'ftp://ftp.informatics.jax.org/pub/reports/MPheno_OBO.ontology'
    
    fln,status = update.downloadSuppressPrintOuts(url,program_dir+'OBO/','')

def downloadPhenotypeOntologyGeneAssociations():
    url = 'ftp://ftp.informatics.jax.org/pub/reports/HMD_HumanPhenotype.rpt'
    #url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/HMD_HumanPhenotype.rpt'
    ### Mouse and human gene symbols and gene IDs (use the gene symbols)
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Pheno/','')

def downloadBioGRIDAssociations():
    print 'Downloading BioGRID associations'
    url = 'http://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-ALL-LATEST.tab2.zip'
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/BioGRID/','')

def downloadDrugBankAssociations():
    print 'Downloading DrugBank associations'
    url = 'http://www.drugbank.ca/system/downloads/current/drugbank.txt.zip'
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/DrugBank/','')
    
def downloadPathwayCommons():
    print 'Downloading PathwayCommons associations'
    url = 'http://www.pathwaycommons.org/pc-snapshot/current-release/gsea/by_species/homo-sapiens-9606-gene-symbol.gmt.zip'
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/PathwayCommons/','')
    
def downloadDiseaseOntologyOBO():
    print 'Downloading Disease Ontology structure and associations'
    
    """ Unfortunately, we have to download versions that are not as frequently updated, since RGDs server
        reliability is poor """
    #url = 'ftp://rgd.mcw.edu/pub/data_release/ontology_obo_files/disease/CTD.obo'
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/CTD.obo'
    
    ### Includes congenital and environmental diseases - http://ctdbase.org/detail.go?type=disease&acc=MESH%3aD002318
    fln,status = update.downloadSuppressPrintOuts(url,program_dir+'OBO/','')
    
def downloadDiseaseOntologyGeneAssociations(selected_species):
    if selected_species == None: sc = []
    else: sc = selected_species
    """ Unfortunately, we have to download versions that are not as frequently updated, since RGDs server
        reliability is poor """
        
    if 'Hs' in sc or len(sc)==0:
        #url = 'ftp://rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/homo_genes_do'
        url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/homo_genes_do'
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Disease/','')

    if 'Mm' in sc or len(sc)==0:
        #url = 'ftp://rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/mus_genes_do'
        url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/mus_genes_do'
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Disease/','')

    if 'Rn' in sc or len(sc)==0:
        #url = 'ftp://rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/rattus_genes_do'
        url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/rattus_genes_do'
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Disease/','')

def downloadMiRDatabases(species):
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/'+species+'_microRNA-Ensembl-GOElite_strict.txt'
    selected = ['Hs','Mm','Rn'] ### these are simply zipped where the others are not
    ### These files should be updated on a regular basis
    if species in selected:
        url = string.replace(url,'.txt','.zip')
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/microRNATargets/','')
    else:
        ### Where strict is too strict
        url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/'+species+'_microRNA-Ensembl-GOElite_lax.txt'
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/microRNATargets/','')
    fln = string.replace(fln,'.zip','.txt')
    return fln

def downloadRvistaDatabases(species):
    ### Source files from  http://hazelton.lbl.gov/pub/poliakov/wgrvista_paper/
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/'+species+'_RVista_factors.zip'
    ### These files should be updated on a regular basis
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/RVista/','')
    fln = string.replace(fln,'.zip','.txt')
    return fln

def remoteDownloadEnsemblTranscriptAssocations(species):
    global program_dir
    program_type,database_dir = unique.whatProgramIsThis()
    if program_type == 'AltAnalyze':
        program_dir = database_dir
    downloadEnsemblTranscriptAssociations(species)
    
def downloadEnsemblTranscriptAssociations(species):
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/Transcripts/'+species+'/Ensembl-EnsTranscript.txt'
    ### These files should be updated on a regular basis
    fln,status = update.downloadSuppressPrintOuts(url,program_dir+species+'/uid-gene/','')

def downloadGOSlimOBO():
    url = 'http://geneontology.org/ontology/subsets/goslim_pir.obo'
    #url = 'http://www.geneontology.org/GO_slims/goslim_generic.obo'  ### Missing 
    fln,status = update.downloadSuppressPrintOuts(url,database_dir+'OBO/','')
    
def importUniProtAnnotations(species_db,force):
    base_url = 'http://www.altanalyze.org/archiveDBs/'
    uniprot_ensembl_db={}
    for species in species_db:
        url = base_url+species+'/custom_annotations.txt'
        if force=='yes':
            fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/UniProt/'+species+'/','')
        else:
            fln = 'BuildDBs/UniProt/'+species+'/custom_annotations.txt'
        for line in open(fln,'rU').xreadlines():
            data = cleanUpLine(line)
            try:
                ens_gene,compartment,function,symbols,full_name,uniprot_name,uniprot_ids,unigene = string.split(data,'\t')
                symbols = string.split(string.replace(symbols,'; Synonyms=',', '),', ')
                uniprot_ensembl_db[species,uniprot_name] = ens_gene
                species_extension = string.split(uniprot_name,'_')[-1]
                full_name = string.split(full_name,';')[0]
                if 'Transcription factor' in full_name:
                    symbols.append(string.split(full_name,'Transcription factor ')[-1]) ### Add this additional synonym to symbols
                ### Extend this database out to account for weird names in PAZAR
                for symbol in symbols:
                    new_name = string.upper(symbol)+'_'+species_extension
                    if new_name not in uniprot_ensembl_db:
                        uniprot_ensembl_db[species,symbol+'_'+species_extension] = ens_gene
                    uniprot_ensembl_db[species,string.upper(symbol)] = ens_gene
            except Exception:
                None
    return uniprot_ensembl_db

############# Import/processing/export #############    
def getPAZARFileNames():
    """ Filenames are manually and periodically downloaded from: http://www.pazar.info/cgi-bin/downloads_csv.pl"""
    fn = filepath('Config/PAZAR_list.txt')
    x=0
    filenames = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            resource, filename = string.split(data,'\t')
            filenames[resource]=filename
    return filenames

class TFTargetInfo:
    def __init__(self,tf_name,ens_gene,project,pmid,analysis_method):
        self.tf_name=tf_name
        self.ens_gene=ens_gene
        self.project=project
        self.pmid=pmid
        self.analysis_method=analysis_method
    def TFName(self): return self.tf_name
    def Ensembl(self): return self.ens_gene
    def Project(self):
        if self.project[-1]=='_':
            return self.project[:-1]
        else:
            return self.project
    def PMID(self): return self.pmid
    def AnalysisMethod(self): return self.analysis_method
    def __repr__(self): return self.TFName()
    
def importPAZARAssociations(force):
    pazar_files = unique.read_directory('/BuildDBs/tftargets')
    species_db={}
    tf_to_target={}
    tf_name_to_ID_db={}
    for file in pazar_files:
        if '.csv' in file:
            name = string.join(string.split(file,'_')[1:-1],'_')
            fn = filepath('BuildDBs/tftargets/'+file)
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                try:
                    ### Each line contains the following 11 tab-delim fields:
                    ### Fields are: <PAZAR TF ID>  <TF Name>  <PAZAR Gene ID>  <ensembl gene accession>  <chromosome>  <gene start coordinate>  <gene end coordinate>  <species>  <project name>  <PMID>  <analysis method> 
                    pazar_tf_id, ens_tf_transcript, tf_name, pazar_geneid, ens_gene, chr, gene_start,gene_end,species,project,pmid,analysis_method = string.split(data,'\t')
                    if ens_tf_transcript == 'ENSMUST00000105345' and 'Pluripotency' in project:
                        ### This is a specific error (TCF3 corresponds to TCF7L1 but poor nomenclature resulted in a mis-annotation here)
                        ens_tf_transcript = 'ENSMUST00000069536'
                    species,genus = string.split(species,' ')
                    species = species[0]+genus[0]
                    tft=TFTargetInfo(tf_name,ens_gene,project,pmid,analysis_method)
                    try: tf_to_target[species,tf_name].append(tft)
                    except Exception: tf_to_target[species,tf_name] = [tft]
                    species_db[species]=[]
                    tf_name_to_ID_db[tf_name] = ens_tf_transcript ### This is an Ensembl transcript ID -> convert to gene ID
                except Exception:
                    None ### Occurs due to file formatting issues (during an update?)
    
    if determine_tf_geneids == 'yes':
        """ The below code is probably most useful for creation of complex regulatory inference networks in Cytoscape """
        uniprot_ensembl_db = importUniProtAnnotations(species_db,force) ### Get UniProt IDs that often match Pazar TF names
        transcript_to_gene_db={}
        gene_to_symbol_db={}
        #species_db={}
        #species_db['Mm']=[]
        for species in species_db:
            try:
                try: gene_to_transcript_db = gene_associations.getGeneToUid(species,('hide','Ensembl-EnsTranscript')); #print mod_source, 'relationships imported.'
                except Exception:
                    try:
                        downloadEnsemblTranscriptAssociations(species)
                        gene_to_transcript_db = gene_associations.getGeneToUid(species,('hide','Ensembl-EnsTranscript'))
                    except Exception: gene_to_transcript_db={}
                try: gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
                except Exception:
                    ### Download this
                    base_url = 'http://www.altanalyze.org/archiveDBs/'
                    url = base_url+species+'/Ensembl-Symbol.txt'
                    update.downloadSuppressPrintOuts(url,'AltDatabase/goelite/'+species+'/uid-gene/','')
                    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
            except Exception: gene_to_transcript_db={}; gene_to_symbol={}
            #print len(gene_to_transcript_db), species
            for gene in gene_to_transcript_db:
                for transcript in gene_to_transcript_db[gene]:
                    transcript_to_gene_db[transcript]=gene
            for gene in gene_to_symbol:
                gene_to_symbol_db[gene] = gene_to_symbol[gene]

        missing=[]
        for (species,tf_name) in tf_to_target:
            original_tf_name = tf_name
            try:
                ens_tf_transcript = tf_name_to_ID_db[tf_name]
                ens_gene = transcript_to_gene_db[ens_tf_transcript]
                #if 'ENSMUST00000025271' == ens_tf_transcript: print ens_gene;kill
                #print gene_to_symbol_db[ens_gene];sys.exit()
                symbol = string.lower(gene_to_symbol_db[ens_gene][0]) ### covert gene ID to lower case symbol ID
                ### Store the original TF name and Ensembl symbol (different species TFs can have different symbols - store all)
                try: tf_to_target_symbol[original_tf_name].append(symbol)
                except Exception: tf_to_target_symbol[original_tf_name] = [symbol]
            except Exception:
                try:
                    #ens_tf_transcript = tf_name_to_ID_db[tf_name]
                    #print species, tf_name, ens_tf_transcript;sys.exit()
                    ens_gene = uniprot_ensembl_db[species,tf_name]
                    symbol = string.lower(gene_to_symbol_db[ens_gene][0]) ### covert gene ID to lower case symbol ID
                    try: tf_to_target_symbol[original_tf_name].append(symbol)
                    except Exception: tf_to_target_symbol[original_tf_name] = [symbol]
                except Exception:
                    try:
                        tf_name = string.split(tf_name,'_')[0]
                        ens_gene = uniprot_ensembl_db[species,tf_name]
                        symbol = string.lower(gene_to_symbol_db[ens_gene][0]) ### covert gene ID to lower case symbol ID
                        try: tf_to_target_symbol[original_tf_name].append(symbol)
                        except Exception: tf_to_target_symbol[original_tf_name] = [symbol]
                    except Exception:
                        try:
                            tf_names=[]
                            if '/' in tf_name:
                                tf_names = string.split(tf_name,'/')
                            elif ' ' in tf_name:
                                tf_names = string.split(tf_name,' ')
                            for tf_name in tf_names:
                                ens_gene = uniprot_ensembl_db[species,tf_name]
                                symbol = string.lower(gene_to_symbol_db[ens_gene][0]) ### covert gene ID to lower case symbol ID
                                try: tf_to_target_symbol[original_tf_name].append(symbol)
                                except Exception: tf_to_target_symbol[original_tf_name] = [symbol]
                        except Exception: missing.append((tf_name,species))
        print 'Ensembl IDs found for transcript or UniProt Transcription factor names:',len(tf_to_target_symbol),'and missing:', len(missing)
        #print missing[:20]
        
    ### Translate all species data to gene symbol to export for all species
    species_tf_targets={}
    for (species,tf_name) in tf_to_target:
        try:
            tf_db = species_tf_targets[species]
            tf_db[tf_name] = tf_to_target[species,tf_name]
        except Exception:
            tf_db = {}
            tf_db[tf_name] = tf_to_target[species,tf_name]
            species_tf_targets[species] = tf_db
        
    tf_dir = 'BuildDBs/tftargets/symbol/tf-target.txt'
    tf_data = export.ExportFile(tf_dir)
    tf_to_symbol={}
    #print 'Exporting:',tf_dir
    #print len(species_tf_targets)
    for species in species_tf_targets:
        try: gene_to_source_id = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
        except Exception: gene_to_source_id={}
        tf_db = species_tf_targets[species]
        for tf_name in tf_db:
            for tft in tf_db[tf_name]:
                try:
                    for symbol in gene_to_source_id[tft.Ensembl()]:
                        symbol = string.lower(symbol)
                        tf_id = tf_name+'(Source:'+tft.Project()+'-PAZAR'+')'
                        tf_data.write(tf_id+'\t'+symbol+'\n')
                        try: tf_to_symbol[tf_id].append(symbol)
                        except Exception: tf_to_symbol[tf_id] = [symbol]
                except Exception: null=[]; 
    tf_data.close()
    tf_to_symbol = gene_associations.eliminate_redundant_dict_values(tf_to_symbol)
    return tf_to_symbol

def importPAZARcompiled():
    """ Skips over the above function when these tf-target file is downlaoded directly """
    tf_dir = 'BuildDBs/tftargets/symbol/tf-target.txt'
    tf_to_symbol={}
    fn = filepath(tf_dir)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        tf_id,symbol = string.split(data,'\t')
        try: tf_to_symbol[tf_id].append(symbol)
        except Exception: tf_to_symbol[tf_id] = [symbol]
    tf_to_symbol = gene_associations.eliminate_redundant_dict_values(tf_to_symbol)
    return tf_to_symbol

def importPhenotypeOntologyGeneAssociations():
    x=0
    pheno_symbol={}; phen=[]
    fn = filepath('BuildDBs/Pheno/HMD_HumanPhenotype.rpt')
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        hs_symbol=t[0]; hs_entrez=t[1]; mm_symbol=t[2]; mgi=t[3]; pheno_ids=t[4]
        if 'MP' not in pheno_ids:
            try: pheno_ids = t[5]
            except Exception: pass
        hs_symbol = string.lower(hs_symbol)
        mm_symbol = string.lower(mm_symbol)
        symbols = [mm_symbol,hs_symbol]
        pheno_ids = string.split(pheno_ids,' '); phen+=pheno_ids
        for pheno_id in pheno_ids:
            if len(pheno_id)>0:
                for symbol in symbols:
                    try: pheno_symbol[pheno_id].append(symbol)
                    except Exception: pheno_symbol[pheno_id]=[symbol]
    phen = unique.unique(phen)
    pheno_symbol = gene_associations.eliminate_redundant_dict_values(pheno_symbol)
    return pheno_symbol

def importAmandeusPredictions(force):
    if force == 'yes':
        downloadAmadeusPredictions()
        
    x=0
    tf_symbol_db={}
    fn = filepath('BuildDBs/Amadeus/symbol-Metazoan-Amadeus.txt')
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            symbol,system,tf_name = string.split(data,'\t')
            symbol = string.lower(symbol)
            if tf_name == 'Oct4': tf_name = 'Pou5f1' ### Known annotation issue
            try: tf_symbol_db[tf_name].append(symbol)
            except Exception: tf_symbol_db[tf_name]=[symbol]
    tf_symbol_db = gene_associations.eliminate_redundant_dict_values(tf_symbol_db)
    
    if determine_tf_geneids == 'yes':
        """ ### Since this data is species independent (not indicated in the current file) can't apply this yet
        uniprot_ensembl_db = importUniProtAnnotations(species_db,force)
        try: gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
        except Exception: gene_to_symbol={}
        symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
        symbol_to_gene = lowerSymbolDB(symbol_to_gene)
        uniprot_ensembl_db = lowerSymbolDB(uniprot_ensembl_db)
        """
        ### Build a TFname to Ensembl gene symbol database for Amadeus
        for tf_name in tf_symbol_db:
            symbol = string.lower(string.split(tf_name,'(')[0])
            if 'miR' not in tf_name and 'let' not in tf_name: ### Exlude miRNAs
                try: tf_to_target_symbol[tf_name].append(symbol)
                except Exception: tf_to_target_symbol[tf_name] = [symbol]
    return tf_symbol_db
    
def importDiseaseOntologyGeneAssocations():
    disease_ontology_files = unique.read_directory('/BuildDBs/Disease')
    symbol_to_DO={}
    for file in disease_ontology_files:
        if '_do' in file:
            fn = filepath('BuildDBs/Disease/'+file)
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if len(t)>1:
                    symbol=string.lower(t[2]); doid = t[4]
                    try: symbol_to_DO[doid].append(symbol)
                    except Exception: symbol_to_DO[doid]=[symbol]
    return symbol_to_DO
    
def exportSymbolRelationships(pathway_to_symbol,selected_species,pathway_type,type):    
    if selected_species != None: ### Restrict to selected species only
        current_species_dirs=selected_species
    else:
        current_species_dirs = unique.read_directory('/'+database_dir)
    
    for species in current_species_dirs:
        if '.' not in species:
            ens_dir = database_dir+'/'+species+'/gene-'+type+'/Ensembl-'+pathway_type+'.txt'
            ens_data = export.ExportFile(ens_dir)
            try:
                if determine_tf_geneids == 'yes':
                    ### When TF data is exported and TF gene-gene interactions are defined, export them
                    tf_network_dir = 'BuildDBs/TF_Interactions/'+species+'/interactions.txt'
                    tf_network_data = export.ExportFile(tf_network_dir)
                    tf_network_data.write('Symbol1\tInteractionType\tSymbol2\tGeneID1\tGeneID2\tSource\n')
                    interaction = 'transcriptional_target'
                    print 'Exporting TF-Target gene-gene relationships to',tf_network_dir
                    tf_count=0; unique_interactions={}
                    try:
                        gene_chr_db = gene_associations.getGeneToUid(species,('hide','Ensembl-chr'))
                    except Exception: 
                        try:
                            ensembl_version = unique.getCurrentGeneDatabaseVersion()
                            from build_scripts import EnsemblSQL
                            EnsemblSQL.getChrGeneOnly(species,'Basic',ensembl_version,'yes')
                            gene_chr_db = gene_associations.getGeneToUid(species,('hide','Ensembl-chr'))
                        except Exception: gene_chr_db={}
                    
            except Exception: None
            if 'mapp' in type: ens_data.write('GeneID\tSystem\tGeneSet\n')
            else: ens_data.write('GeneID\tGeneSet\n')
            try: ens_to_entrez = gene_associations.getGeneToUid(species,('hide','Ensembl-EntrezGene'))
            except Exception: ens_to_entrez ={}
            if len(ens_to_entrez)>0:
                entrez_dir = database_dir+'/'+species+'/gene-'+type+'/EntrezGene-'+pathway_type+'.txt'
                entrez_data = export.ExportFile(entrez_dir)
                if 'mapp' in type: entrez_data.write('GeneID\tSystem\tGeneSet\n')
                else: entrez_data.write('GeneID\tGeneSet\n')
            #print 'Exporting '+pathway_type+' databases for:',species
            try: gene_to_source_id = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
            except Exception: gene_to_source_id={}
            source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
            source_to_gene = lowerSymbolDB(source_to_gene)
            for pathway in pathway_to_symbol:
                for symbol in pathway_to_symbol[pathway]:
                    try:
                        genes = source_to_gene[symbol]
                        for gene in genes:
                            if len(genes)<5: ### don't propagate redundant associations
                                if 'mapp' in type: ens_data.write(gene+'\tEn\t'+pathway+'\n')
                                else: ens_data.write(gene+'\t'+pathway+'\n')
                            if gene in ens_to_entrez:
                                for entrez in ens_to_entrez[gene]:
                                    if 'mapp' in type: entrez_data.write(entrez+'\tL\t'+pathway+'\n')
                                    else: entrez_data.write(entrez+'\t'+pathway+'\n')
                            try:
                                if determine_tf_geneids == 'yes':
                                    if '(' in pathway:
                                        source_name = string.split(pathway,'(')[0]
                                    else:
                                        source_name =  pathway
                                    proceed = True
                                    if gene in gene_chr_db:
                                        if len(gene_chr_db[gene][0])>2: ### not a valid chromosome (e.g., HSCHR6_MHC_COX)
                                            proceed = False
                                    if proceed ==True or proceed == False:
                                        try: symbols = tf_to_target_symbol[source_name]
                                        except Exception: symbols = tf_to_target_symbol[pathway]
                                        for symbol in symbols:
                                            tf_gene = source_to_gene[symbol][0] ### meta-species converted symbol -> species Ensembl gene
                                            tf_symbol = gene_to_source_id[tf_gene][0] ### species formatted symbol for TF
                                            symbol = gene_to_source_id[gene][0] ### species formatted symbol for target
                                            if (tf_symbol,symbol,pathway) not in unique_interactions:
                                                tf_network_data.write(string.join([tf_symbol,interaction,symbol,tf_gene,gene,pathway],'\t')+'\n')
                                                unique_interactions[tf_symbol,symbol,pathway]=[]
                                                try: merged_tf_interactions[tf_symbol].append(string.lower(symbol))
                                                except Exception: merged_tf_interactions[tf_symbol] = [string.lower(symbol)]
                                                tf_count+=1
                            except Exception: None
                    except Exception: null=[]
            ens_data.close()
            try: entrez_data.close()
            except Exception: null=[]
            try:
                if determine_tf_geneids == 'yes':
                    tf_network_data.close()
                    print tf_count,'TF to target interactions exported..'
            except Exception: None

def translateBioMarkersBetweenSpecies(input_dir,species):
    ### Convert the species Ensembl primary key IDs from the source to 
    try:
        biomarker_files = unique.read_directory(input_dir)
    except Exception:
        biomarker_files = unique.read_directory(input_dir)
        
    try: gene_to_source_id = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    except Exception: gene_to_source_id={}
    source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
    source_to_gene = lowerSymbolDB(source_to_gene)
    print len(source_to_gene)

    marker_symbol_db={}
    for file in biomarker_files:
        if 'tissue-specific' in file and species not in file: ### Don't re-translate an already translated file
            export_dir = 'AltDatabase/ensembl/'+species+'/'+species+file[2:]
            export_data = export.ExportFile(export_dir)
            fn = filepath(input_dir+'/'+file)
            x=0
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if x==0:
                    x = 1; y=0
                    for i in t:
                        if 'Symbol' in i: sy = y
                        y+=1
                    export_data.write(line)
                else:
                    ensembl = t[0]; symbol = string.lower(t[sy])
                    #print symbol, len(source_to_gene);sys.exit()
                    if symbol in source_to_gene:
                        species_gene = source_to_gene[symbol][0] ### should only be one gene per symbol (hopefully)
                        values = string.join([species_gene]+t[1:],'\t')+'\n'
                        export_data.write(values)
            export_data.close()
            print export_dir,'written...'
         
def extractKEGGAssociations(species,mod,system_codes):
    import_dir = filepath('/BuildDBs/KEGG')
    g = gene_associations.GrabFiles(); g.setdirectory(import_dir)
    filedir = g.getMatchingFolders(species)
    gpml_data,pathway_db = gene_associations.parseGPML(filepath(filedir))
    gene_to_WP = gene_associations.unifyGeneSystems(gpml_data,species,mod)
    gene_associations.exportCustomPathwayMappings(gene_to_WP,mod,system_codes,filepath(database_dir+'/'+species+'/gene-mapp/'+mod+'-KEGG.txt'))
    if len(gene_to_WP)>0 and mod == 'Ensembl': ### Export all pathway interactions
        try: gene_associations.exportNodeInteractions(pathway_db,mod,filepath(filedir))
        except Exception: null=[]
    elif len(gene_to_WP)>0 and mod == 'HMDB': ### Export all pathway interactions
        try: gene_associations.exportNodeInteractions(pathway_db,mod,filepath(filedir),appendToFile=True)
        except Exception: null=[]
    return filedir

def extractGMTAssociations(species,mod,system_codes,data_type):
    if mod != 'HMDB':
        import_dir = filepath('/BuildDBs/'+data_type)
        gmt_data = gene_associations.parseGMT(import_dir)
        gene_to_custom = gene_associations.unifyGeneSystems(gmt_data,species,mod)
        gene_associations.exportCustomPathwayMappings(gene_to_custom,mod,system_codes,filepath(database_dir+'/'+species+'/gene-mapp/'+mod+'-'+data_type+'.txt'))
        
def transferGOSlimGeneAssociations(selected_species):
    if selected_species != None: ### Restrict to selected species only
        current_species_dirs=selected_species
    else:
        current_species_dirs = unique.read_directory('/'+database_dir)
    for species_code in current_species_dirs:
        try:
            ens_go_file_dir = filepath(database_dir+'/'+species_code+'/gene-go/Ensembl-GOSlim.txt')
            goslim_ens_file = filepath(database_dir+'/'+species_code+'/uid-gene/Ensembl-goslim_goa.txt')
            export.copyFile(goslim_ens_file,ens_go_file_dir)
            translateToEntrezGene(species_code,ens_go_file_dir)
        except Exception: null=[]

def translateToEntrezGene(species,filename):
    x=0; type = 'pathway'
    try: ens_to_entrez = gene_associations.getGeneToUid(species,('hide','Ensembl-EntrezGene'))
    except Exception: ens_to_entrez ={}

    if len(ens_to_entrez)>0:
        export_file = string.replace(filename,'Ensembl','EntrezGene')
        export_data = export.ExportFile(export_file)
        export_data.write('EntrezGene\tOntologyID\n')
        fn = filepath(filename)
        for line in open(fn,'rU').xreadlines():
            if x==0: x=1
            else:
                data = cleanUpLine(line)
                try:
                    ensembl,pathway = string.split(data,'\t')
                    type = 'ontology'
                except Exception:
                    ensembl,null,pathway = string.split(data,'\t')
                try:
                    entrezs = ens_to_entrez[ensembl]
                    for entrez in entrezs:
                        if type == 'ontology':
                            export_data.write(entrez+'\t'+pathway+'\n')
                        else:
                            export_data.write(entrez+'\tEn\t'+pathway+'\n')
                except Exception:
                    null=[]
        export_data.close()

def importRVistaGeneAssociations(species_code,source_path):
    x=0; tf_symbol_db={}
    fn = filepath(source_path)
    
    TF_symbol_db={}
    if species_code == 'Dm' or species_code == 'Mm': ### Dm IDs are Ensembl
        gene_to_symbol = gene_associations.getGeneToUid(species_code,('hide','Ensembl-Symbol'))
    increment = 10000
    print 'Importing symbol-R Vista TF assocations (be patient)'
    x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        x+=1
        #if x==increment: x=0; print '*',
        try:
            symbol,TF = string.split(data,'\t')
            if species_code == 'Dm' or species_code == 'Mm':
                ensembls = [symbol]
                try: symbol = gene_to_symbol[symbol][0]
                except Exception: forceError
            try: TF_symbol_db[TF].append(string.lower(symbol))
            except Exception: TF_symbol_db[TF]=[string.lower(symbol)]
        except Exception: None
    TF_symbol_db = gene_associations.eliminate_redundant_dict_values(TF_symbol_db)
    return TF_symbol_db
    
def importMiRGeneAssociations(species_code,source_path):
    try:
        destination_path = filepath(database_dir+'/'+species_code+'/gene-mapp/Ensembl-microRNATargets.txt')
        export.copyFile(source_path,destination_path)
        translateToEntrezGene(species_code,destination_path)
    except Exception: null=[]            
        
def importBioMarkerGeneAssociations(input_dir):
    try:
        biomarker_folder = unique.read_directory(input_dir)
    except Exception:
        biomarker_folder = unique.read_directory(input_dir)
    marker_symbol_db={}
    for folder in biomarker_folder:
        if '.' in folder: continue ### This is a file not a folder
        else: biomarker_files = unique.read_directory(input_dir+'/'+folder)
        for file in biomarker_files:
            x=0
            if '.txt' in file and '-correlation' not in file and 'AltExon' not in file:
                fn = filepath('BuildDBs/BioMarkers/'+folder+'/'+file)
                for line in open(fn,'rU').xreadlines():
                    data = cleanUpLine(line)
                    t = string.split(data,'\t')
                    if x==0:
                        x = 1; y=0
                        for i in t:
                            if 'marker-in' in i: mi = y
                            if 'Symbol' in i: sy = y
                            y+=1
                    ensembl = t[0]; symbol = string.lower(t[sy]); marker = t[mi]
                    markers = string.split(marker,'|')
                    for marker in markers:
                        try: marker_symbol_db[marker].append(symbol)
                        except Exception: marker_symbol_db[marker]=[symbol] 
    marker_symbol_db = gene_associations.eliminate_redundant_dict_values(marker_symbol_db)
    return marker_symbol_db

def importDomainGeneAssociations(species_code,source_path):
    try:
        destination_path = filepath(database_dir+'/'+species_code+'/gene-mapp/Ensembl-Domains.txt')
        export.copyFile(source_path,destination_path)
        translateToEntrezGene(species_code,destination_path)
    except Exception: null=[]            

def importBioGRIDGeneAssociations(taxid,species):
    model_mammal_tax = {}
    model_mammal_tax['9606'] = 'Hs'
    model_mammal_tax['10090'] = 'Mm'
    model_mammal_tax['10116'] = 'Rn'
    
    filtered=considerOnlyMammalian([species]) ### See if the species is a mammal
    if len(filtered)==0: model_mammal_tax={} ### Don't inlcude the gold standard mammals if not a mammal
    
    model_mammal_tax[taxid]=species
    
    biogrid_files = unique.read_directory('/BuildDBs/BioGRID')
    latest_file = biogrid_files[-1]
    fn = filepath('BuildDBs/BioGRID/'+latest_file)
    
    ens_dir = database_dir+'/'+species+'/gene-interactions/Ensembl-BioGRID.txt'
    ens_data = export.ExportFile(ens_dir)
    ens_data.write('Symbol1\tInteractionType\tSymbol2\tGeneID1\tGeneID2\tSource\n')
    
    try: gene_to_source_id = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    except Exception: gene_to_source_id={}
    source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
    source_to_gene = lowerSymbolDB(source_to_gene)
            
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        species_tax = t[15]
        if species_tax in model_mammal_tax:
            symbol1 = t[7]; symbol2 = t[8]
            source_exp = t[11]; interaction_type = t[12]
            try:
                ens1 = source_to_gene[string.lower(symbol1)][0]
                ens2 = source_to_gene[string.lower(symbol2)][0]
                values = string.join([symbol1,interaction_type,symbol2,ens1,ens2,source_exp],'\t')+'\n'
                ens_data.write(values)
            except Exception:
                None

def importDrugBankAssociations(species):
    fn = filepath('BuildDBs/DrugBank/drugbank.txt')
    
    ens_dir = database_dir+'/'+species+'/gene-interactions/Ensembl-DrugBank.txt'
    ens_data = export.ExportFile(ens_dir)
    ens_data.write('DrugName\tMechanism\tGeneSymbol\tDrugBankDB-ID\tGeneID\tSource\n')
    
    try: gene_to_source_id = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    except Exception: gene_to_source_id={}
    source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
    source_to_gene = lowerSymbolDB(source_to_gene)
    
    getCAS=False
    getGenericName=False
    getMechanim = False
    getGeneName = False
    geneNames=[]
    mechanism=''
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        
        if data == '# Primary_Accession_No:': getCAS=True ### switched this from CAS to Drug bank ID for consistency and namining simplicity
        elif getCAS: casID = data; getCAS=False
        
        if data == '# Generic_Name:': getGenericName=True
        elif getGenericName: genericName = data; getGenericName=False
        
        if data == '# Mechanism_Of_Action:': getMechanim=True
        elif getMechanim:
            if len(data)>0: mechanism += data + ' '
            else: getMechanim=False
            
        if '# Drug_Target_' in data and '_Gene_Name' in data: getGeneName=True
        elif getGeneName: geneNames.append(data); getGeneName=False
        
        if '#END_DRUGCARD' in data:
            for symbol in geneNames:
                try:
                    ens = source_to_gene[string.lower(symbol)][0]
                    values = string.join([genericName,mechanism,symbol,casID,ens,'DrugBank'],'\t')+'\n'
                    ens_data.write(values)
                except Exception:
                    None
            casID=''; genericName=''; mechanism=''; geneNames=[]
    
############# Central buid functions #############

def importWikiPathways(selected_species,force):
    if selected_species == None:
        selected_species = unique.read_directory('/'+database_dir)
    importSpeciesData()
    getSourceData()
    all_species = 'no'
    if force == 'yes':
        try:
            gene_associations.convertAllGPML(selected_species,all_species) ### Downloads GPMLs and builds flat files
            for species_code in selected_species:
                interaction_file = 'GPML/Interactomes/interactions.txt'
                moveInteractionsToInteractionsDir(interaction_file,species_code,'WikiPathways')
            status = 'built'
        except IOError:
            print 'Unable to connect to http://www.wikipathways.org'
            status = 'failed'
    status = 'built'
    if status == 'built':
        from import_scripts import BuildAffymetrixAssociations

        for species_code in selected_species:
            species_name = species_names[species_code]
            if status == 'built':          
                relationship_types = ['native','mapped']
                for relationship_type in relationship_types:
                    #print 'Processing',relationship_type,'relationships'
                    index=0
                    integrate_affy_associations = 'no'
                    incorporate_previous = 'yes'
                    process_affygo = 'no'
                    counts = BuildAffymetrixAssociations.importWikipathways(source_types,incorporate_previous,process_affygo,species_name,species_code,integrate_affy_associations,relationship_type,'over-write previous')
                    index+=1
    print 'Finished integrating updated WikiPathways'

def moveInteractionsToInteractionsDir(source_file,species,name):
    destination_file = filepath('AltDatabase/goelite/'+species+'/gene-interactions/Ensembl-'+name+'.txt')
    source_file = filepath(source_file)
    try: export.copyFile(source_file,destination_file)
    except Exception: None ### No file to move
    
def importKEGGAssociations(selected_species,force):
    supported_databases = ['Ag','At','Ce','Dm','Dr','Hs','Mm','Os','Rn']
    getSourceData()
    
    if selected_species != None: ### Restrict by selected species
        supported_databases2=[]
        for species in selected_species:
            if species in supported_databases:
                supported_databases2.append(species)
        supported_databases = supported_databases2

    mod_types_list=[]
    for i in mod_types:  mod_types_list.append(i)
    mod_types_list.sort()
    
    for species in supported_databases:
        if force == 'yes':
            downloadKEGGPathways(species)
        for mod in mod_types_list:
            buildDB_dir = extractKEGGAssociations(species,mod,system_codes)
            
        interaction_file = buildDB_dir+'/Interactomes/interactions.txt'
        moveInteractionsToInteractionsDir(interaction_file,species,'KEGG')
                
def importPathwayCommons(selected_species,force):
    original_species = selected_species
    selected_species = considerOnlyMammalian(selected_species)
    if len(selected_species) == 0:
        print 'PLEASE NOTE: %s does not support PathwayCommons update.' % string.join(original_species,',')
    else:
        if force == 'yes':
            downloadPathwayCommons()
    
        getSourceData()
    
        for species in selected_species:
            for mod in mod_types:
                extractGMTAssociations(species,mod,system_codes,'PathwayCommons')
    
def importTranscriptionTargetAssociations(selected_species,force):
    original_species = selected_species
    selected_species = considerOnlyMammalian(selected_species)
    x=[]
    if len(selected_species) == 0:
        print 'PLEASE NOTE: %s does not support Transcription Factor association update.' % string.join(original_species,',')
    else:
        global determine_tf_geneids
        global tf_to_target_symbol ### Used for TF-target interaction networks
        global merged_tf_interactions
        tf_to_target_symbol={}
        source = 'raw'#'precompiled'
        merged_tf_interactions={} ### Stores the final merged PAZAR-Amadeus merged data
        determine_tf_geneids = 'yes'
        
        ### No need to specify a species since the database will be added only to currently installed species
        if force == 'yes':
            source = downloadPAZARAssocations()
        if source == 'raw':
            x = importPAZARAssociations(force) ### builds the PAZAR TF-symbol associations from resource.csv files
        if source == 'precompiled' or len(x)==0:
            x = importPAZARcompiled() ### imports from pre-compiled/downloaded TF-symbol associations
            
        y = importAmandeusPredictions(force)
        z = dict(x.items() + y.items())
        geneset = 'TFTargets'
        if determine_tf_geneids == 'yes':
            tf_to_target_symbol = gene_associations.eliminate_redundant_dict_values(tf_to_target_symbol)
            exportSymbolRelationships(z,selected_species,geneset,'mapp')
            determine_tf_geneids = 'no'
            geneset = 'MergedTFTargets'
            z = merged_tf_interactions
        exportSymbolRelationships(merged_tf_interactions,selected_species,geneset,'mapp')
        
        for species in selected_species:
            interaction_file = 'BuildDBs/TF_Interactions/'+species+'/interactions.txt'
            moveInteractionsToInteractionsDir(interaction_file,species,'TFTargets')

def importPhenotypeOntologyData(selected_species,force):
    original_species = selected_species
    selected_species = considerOnlyMammalian(selected_species)
    if len(selected_species) == 0:
        print 'PLEASE NOTE: %s does not support Phenotype Ontology update.' % string.join(original_species,',')
    else:
        ### No need to specify a species since the database will be added only to currently installed species    
        if force == 'yes':
            downloadPhenotypeOntologyOBO()
            downloadPhenotypeOntologyGeneAssociations()
        x = importPhenotypeOntologyGeneAssociations()
        exportSymbolRelationships(x,selected_species,'MPhenoOntology','go')
    
def importDiseaseOntologyAssociations(selected_species,force):
    original_species = selected_species
    selected_species = considerOnlyMammalian(selected_species)
    if len(selected_species) == 0:
        print 'PLEASE NOTE: %s does not support Disease Ontology update.' % string.join(original_species,',')
    else:
        if force == 'yes':
            downloadDiseaseOntologyOBO()
            downloadDiseaseOntologyGeneAssociations(selected_species)
        x = importDiseaseOntologyGeneAssocations()
        exportSymbolRelationships(x,selected_species,'CTDOntology','go')

def importGOSlimAssociations(selected_species,force):
    if force == 'yes':
        downloadGOSlimOBO()
    transferGOSlimGeneAssociations(selected_species)
    
def importRVistaAssocations(selected_species,force):
    supported_databases = ['Hs','Mm','Dm']
    selected_supported_databases=[]
    for species in selected_species:
        if species in supported_databases:
            selected_supported_databases.append(species)
            
    missing_Rvista_associations=[]
    found_Rvista_associations=[]
    for species in selected_supported_databases:
        if force == 'yes':
            try:
                fn = downloadRvistaDatabases(species)
                found_Rvista_associations.append((species,fn))
            except Exception:
                missing_Rvista_associations.append(species)
        else:
            fn = filepath('BuildDBs/RVista/'+species+'_RVista_factors.txt')
            found_Rvista_associations.append((species,fn))
   
    for (species,fn) in found_Rvista_associations:
        TF_symbol_db = importRVistaGeneAssociations(species,fn)
        exportSymbolRelationships(TF_symbol_db,[species],'RVista_TFsites','mapp')
        
def importMiRAssociations(selected_species,force):
    supported_databases = unique.read_directory('/'+database_dir)
    if selected_species != None: ### Restrict by selected species
        supported_databases=selected_species

    missing_miR_associations=[]
    found_miR_associations=[]
    for species in supported_databases:
        if force == 'yes':
            try:
                fn = downloadMiRDatabases(species)
                found_miR_associations.append((species,fn))
            except Exception:
                missing_miR_associations.append(species)
                
    for (species,fn) in found_miR_associations:
        importMiRGeneAssociations(species,fn)
        
        interaction_file = 'AltDatabase/goelite/'+species+'/gene-mapp/Ensembl-microRNATargets.txt'
        moveInteractionsToInteractionsDir(interaction_file,species,'microRNATargets')
            
def importBioMarkerAssociations(selected_species,force):
    original_species = selected_species
    selected_species = considerOnlyMammalian(selected_species)
    if len(selected_species) == 0:
        print 'PLEASE NOTE: %s does not support BioMarker association update.' % string.join(original_species,',')
    else:
        if force == 'yes':
            downloadBioMarkers('BuildDBs/BioMarkers/')
        x = importBioMarkerGeneAssociations('BuildDBs/BioMarkers')
        exportSymbolRelationships(x,selected_species,'BioMarkers','mapp')

def importDrugBank(selected_species,force):
    if force == 'yes':
        downloadDrugBankAssociations()
    for species in selected_species:
        importDrugBankAssociations(species)
                
def importBioGRID(selected_species,force):
    if force == 'yes':
        downloadBioGRIDAssociations()
    for species in selected_species:
        importSpeciesData() ### Creates the global species_taxids
        if species in species_taxids:
            taxid = species_taxids[species]
            importBioGRIDGeneAssociations(taxid,species)
      
def importDomainAssociations(selected_species,force):
    if force == 'yes':
        paths = downloadDomainAssociations(selected_species)
        for (species,path) in paths:
            path = string.replace(path,'.gz','.txt')
            importDomainGeneAssociations(species, path)
    
def considerOnlyMammalian(selected_species):
    supported_mammals = ['Am','Bt', 'Cf', 'Ch', 'Cj', 'Cp', 'Do', 'Ec', 'Ee', 'Et', 'Fc', 'Gg', 'Go', 'Hs',
                         'La', 'Ma', 'Md', 'Me', 'Mi', 'Ml', 'Mm', 'Oa', 'Oc','Og', 'Op', 'Pc', 'Pp',
                         'Pt', 'Pv', 'Rn', 'Sa', 'Ss', 'St', 'Tb', 'Tn', 'Tr', 'Ts', 'Tt', 'Vp']
    filtered_species=[]
    if selected_species == None:
        selected_species = unique.read_directory('/'+database_dir)
        
    for i in selected_species:
        if i in supported_mammals:
            filtered_species.append(i)
    return filtered_species

def buildInferrenceTables(selected_species):
    for species_code in selected_species:
        file_found = verifyFile(database_dir+'/'+species_code+'/uid-gene/Ensembl-Symbol'+'.txt') ### If file is present, the below is not needed
        if file_found == 'no':
            try: gene_associations.swapAndExportSystems(species_code,'Ensembl','EntrezGene') ### Allows for analysis of Ensembl IDs with EntrezGene based GO annotations (which can vary from Ensembl)
            except Exception: null=[] ### Occurs if EntrezGene not supported
    
            ### Build out these symbol association files
            try: gene_associations.importGeneData(species_code,('export','Ensembl'))
            except Exception: null=[] ### Occurs if EntrezGene not supported
            try: gene_associations.importGeneData(species_code,('export','EntrezGene'))
            except Exception: null=[] ### Occurs if EntrezGene not supported

def exportBioTypes(selected_species):
    for species in selected_species:
        eo = export.ExportFile('AltDatabase/goelite/'+species+'/gene-mapp/Ensembl-Biotypes.txt')
        eo.write('Ensembl\tSystemCode\tClass\n')
        filename = 'AltDatabase/uniprot/'+species+'/custom_annotations.txt'
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            ens_gene,compartment,custom_class = t[:3]
            if 'GPCR' in custom_class:
                custom_class = ['GPCR']
            else:
                custom_class = string.split(custom_class,'|')
            custom_class = string.split(compartment,'|')+custom_class
            for c in custom_class:
                if len(c)>0:
                    eo.write(ens_gene+'\t\t'+c+'\n')
    eo.close() 

def buildAccessoryPathwayDatabases(selected_species,additional_resources,force):
    global database_dir
    global program_dir
    global program_type
    program_type,database_dir = unique.whatProgramIsThis()
    if program_type == 'AltAnalyze':
        program_dir = database_dir
    else:
        program_dir=''

    buildInferrenceTables(selected_species) ### Make sure these tables are present first!!!
        
    print 'Attempting to update:', string.join(additional_resources,',')
    if 'Latest WikiPathways' in additional_resources:
        try: importWikiPathways(selected_species,force)
        except Exception: print 'WikiPathways import failed (cause unknown)'
    if 'KEGG' in additional_resources:
        try: importKEGGAssociations(selected_species,force)
        except Exception: print 'KEGG import failed (cause unknown)'
    if 'Transcription Factor Targets' in additional_resources:
        try: importTranscriptionTargetAssociations(selected_species,force)
        except Exception:
            #print traceback.format_exc()
            print 'Transcription Factor Targets import failed (cause unknown)'
    if 'Phenotype Ontology' in additional_resources:
        try: importPhenotypeOntologyData(selected_species,force)
        except Exception: print 'Phenotype Ontology import failed (cause unknown)'
    if 'Disease Ontology' in additional_resources:
        try: importDiseaseOntologyAssociations(selected_species,force)
        except Exception: print 'Disease Ontology import failed (cause unknown)'
    if 'GOSlim' in additional_resources:
        try: importGOSlimAssociations(selected_species,force)
        except Exception: print 'GOSlim import failed (cause unknown)'
    if 'miRNA Targets' in additional_resources:
        try: importMiRAssociations(selected_species,force)
        except Exception: print 'miRNA Targets import failed (cause unknown)'
    if 'BioMarkers' in additional_resources:
        try: importBioMarkerAssociations(selected_species,force)
        except Exception:
            print 'BioMarkers import failed (cause unknown)'#,traceback.format_exc()
    if 'Domains2' in additional_resources: ### Currently disabled since it's utility is likely low but analysis time is long
        try: importDomainAssociations(selected_species,force)
        except Exception: print 'Domains import failed (cause unknown)'
    if 'PathwayCommons' in additional_resources:
        try: importPathwayCommons(selected_species,force)
        except Exception: print 'PathwayCommons import failed (cause unknown)'
    if 'RVista Transcription Factor Sites' in additional_resources:
        try: importRVistaAssocations(selected_species,force)
        except Exception: print 'R Vista Transcription Factor Site import failed (cause unknown)'
    if 'BioGRID' in additional_resources:
        try: importBioGRID(selected_species,force)
        except Exception: print 'BioGRID import failed (cause unknown)'
    if 'DrugBank' in additional_resources:
        try: importDrugBank(selected_species,force)
        except Exception: print 'Drug Bank import failed (cause unknown)'
    try: exportBioTypes(selected_species)
    except Exception: pass
    
if __name__ == '__main__':
    selected_species = ['Mm']
    force = 'yes'
    download_species = 'Hs'
    species = 'Hs'
    species = 'Mm'
    #exportBioTypes(selected_species);sys.exit()
    additional_resources=['Latest WikiPathways']#'KEGG','BioGRID','DrugBank','miRNA Targets','Transcription Factor Targets']
    #translateBioMarkersBetweenSpecies('AltDatabase/ensembl/'+download_species,species);sys.exit()
    additional_resources=['Latest WikiPathways','PathwayCommons','Transcription Factor Targets','Domains','BioMarkers']
    additional_resources+=['miRNA Targets','GOSlim','Disease Ontology','Phenotype Ontology','KEGG','RVista Transcription Factor Sites']
    additional_resources=['Latest WikiPathways']
    buildAccessoryPathwayDatabases(selected_species,additional_resources,force)
    