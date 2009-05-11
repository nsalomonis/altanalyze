###update
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

import os
import sys
import unique
import string
import export

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

############## Update Databases ##############
def buildAltMouseExonAnnotations(species,array_type,force,genomic_build):
    """Code required to:
    1) Extract out Affymetrix provided exon sequence (probeset sequence extracted from "probeset_sequence_reversed.txt", derived
       directly from the Affymetrix AltMouse probe annotation file), from the "SEQUENCE-transcript-dbase.txt" (built using
       dump-chip1 .gff sequence and AltMerge-Peptide Informatics script "sequence_analysis_AltMouse_refseq.py").
    2) Once exported, grab full length exon sequences using exon/intron coordinates matches to full-length gene sequences with 2kb
       flanking sequence to efficiently predict microRNA binding site exclusion (reAnnotateCriticalExonSequences) and later for
       coordinate mapping to get exons aligning with UCSC annotated splicing annotations and exons. This sequence data replaced
       the previous file (don't need to re-run this - see rederive_exonseq == 'yes' below for reference).
    3) Match the updated exon sequences to the most recent genomic coordinates and build the exact equivalent of the exon array
       Mm_Ensembl_probeset.txt database (same structure and ExonArrayEnsemblRules.py code). This involves running EnsemblImport.
    This code should be run before the exon array location build code since the "Mm_Ensembl_probeset.txt" is created and then re-
    written as "Mm_AltMouse_Ensembl_probeset.txt".
    """
    
    import JunctionArray
    import JunctionArrayEnsemblRules    
    rederive_exonseq = 'no'
    ### Only needs to be run once, to export exon sequence for AltMouse array the original (1 and 2 above)
    if rederive_exonseq == 'yes':
        import AltAnalyze
        import ExonAnnotate_module
        import ExonAnalyze_module
        agglomerate_inclusion_probesets = 'no'; onlyAnalyzeJunctions='no'
        probeset_annotations_file = "AltDatabase/"+species+"/"+array_type+"/"+"MASTER-probeset-transcript.txt"
        exon_db={}; filtered_arrayids={};filter_status='no'
        constituitive_probeset_db,exon_db,genes_being_analyzed = AltAnalyze.importSplicingAnnotationDatabase(probeset_annotations_file,array_type,filtered_arrayids,filter_status)
        alt_junction_db,critical_exon_db,exon_dbase,exon_inclusion_db,exon_db = ExonAnnotate_module.identify_putative_splice_events(exon_db,constituitive_probeset_db,{},agglomerate_inclusion_probesets,onlyAnalyzeJunctions)
        ExonAnnotate_module.exportJunctionComparisons(alt_junction_db,critical_exon_db,exon_dbase)
        print "Finished exporting junctions used in AltMouse array comparisons."

        ExonAnalyze_module.exportAltMouseExonSequence()
        JunctionArray.reAnnotateCriticalExonSequences(species,array_type)

    ### Get UCSC associations (download databases if necessary)
    import UCSCImport
    mRNA_Type = 'mrna'; run_from_scratch = 'yes'
    export_all_associations = 'no' ### YES only for protein prediction analysis
    UCSCImport.runUCSCEnsemblAssociations(species,mRNA_Type,export_all_associations,run_from_scratch,force)
    
    if genomic_build == 'new':
        ### Need to run with every new genomic build (match up new coordinates
        JunctionArray.identifyCriticalExonLocations(species,array_type)
    JunctionArrayEnsemblRules.getAnnotations(species,array_type,force)

def buildExonArrayExonAnnotations(species,force):

    ### Get UCSC associations (download databases if necessary)
    import UCSCImport
    mRNA_Type = 'mrna'; run_from_scratch = 'yes'
    export_all_associations = 'no' ### YES only for protein prediction analysis
    UCSCImport.runUCSCEnsemblAssociations(species,mRNA_Type,export_all_associations,run_from_scratch,force)
    
    import ExonArrayEnsemblRules; reload(ExonArrayEnsemblRules)
    process_from_scratch='yes'
    constitutive_source='default'
    ### Build the databases and return the variables (not used here)
    source_biotype = 'mRNA'
    probeset_db,annotate_db,constitutive_gene_db,splicing_analysis_db = ExonArrayEnsemblRules.getAnnotations(process_from_scratch,constitutive_source,source_biotype,species)

def buildUniProtFunctAnnotations(species,force):
    import UI
    file_location_defaults = UI.importDefaultFileLocations()
    """Identify the appropriate download location for the UniProt database for the selected species"""
    uis = file_location_defaults['UniProt']
    tis = file_location_defaults['TrEMBL']
    
    for ui in uis:
        if species in ui.Species(): uniprot_filename_url = ui.Location()
    for ti in tis:
        if species in ti.Species(): trembl_filename_url = ti.Location()
    species_codes = importSpeciesInfo(); species_full = species_codes[species].SpeciesName()
    
    import ExtractUniProtFunctAnnot; reload(ExtractUniProtFunctAnnot)
    ExtractUniProtFunctAnnot.runExtractUniProt(species,species_full,uniprot_filename_url,trembl_filename_url,force)

class SpeciesData:
    def __init__(self, abrev, species):
        self._abrev = abrev; self._species = species
    def SpeciesCode(self): return self._abrev
    def SpeciesName(self): return self._species
    def __repr__(self): return self.Report()
    
def importSpeciesInfo():
    filename = 'Config/species.txt'; x=0
    fn=filepath(filename); species_codes={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        abrev,species,algorithms = string.split(data,'\t')
        if x==0: x=1
        else:
            sd = SpeciesData(abrev,species)
            species_codes[abrev] = sd
    return species_codes
    
def download(url,dir,file_type):
    try: gz_filepath, status = download_protocol(url,dir,file_type)
    except Exception: gz_filepath='failed'; status = "Internet connection not established. Re-establsih and try again."
    if status == 'remove':
        print "\nRemoving zip file:",gz_filepath
        try: os.remove(gz_filepath); status = 'removed'
        except OSError: status = 'not-removed' ### Not sure why this error occurs since the file is not open
    return gz_filepath, status

def download_protocol(url,dir,file_type):
    """Copy the contents of a file from a given URL to a local file."""
    print "Downloading the following file:",url
    import urllib2
    webFile = urllib2.urlopen(url)
    filename = url.split('/')[-1]
    output_filepath_object = export.createExportFile(dir+filename,dir[:-1])
    output_filepath = filepath(dir+filename)
    localFile = open(output_filepath, 'wb')
    localFile.write(webFile.read())
    webFile.close()
    localFile.close()
    print "File downloaded to:",output_filepath
    if '.zip' in filename:
        import zipfile
        zfile = zipfile.ZipFile(output_filepath)
        for name in zfile.namelist():
            outfile = open(filepath(dir+name), 'wb')
            outfile.write(zfile.read(name)); outfile.close()
        gz_filepath = filepath(output_filepath); status = 'remove'
        print "zip file extracted..."
    if '.gz' in filename:
        gz_filepath = output_filepath
        if len(file_type)==0: extension = '.gz'
        else: extension = 'gz'
        decompressed_filepath = string.replace(gz_filepath,extension,file_type)
        try:
            ### Below code can be too memory intensive
            file_size = os.path.getsize(output_filepath)
            megabtyes = file_size/1000000.00
            if megabtyes>40: force_error ### force_error is an undefined variable which causes an exception
            import gzip; zfile = gzip.GzipFile(gz_filepath); content = zfile.read(); zfile.close()
            data = open(decompressed_filepath,'w')
            print "\nExtracting downloaded file:",gz_filepath; data.write(content); data.close()
            status = 'remove'
        except Exception:
            ### If file(s) will be over-written, user will be prompted
            try: os.remove(decompressed_filepath) ### So OS does not prompt you to see if you want to over-write existing
            except OSError: null=[]
            gz_filepath = string.replace('"'+output_filepath+'"','\\','/')
            os.system("gunzip " + gz_filepath)
            gz_filepath = filepath(output_filepath)
            status = 'remove'
            print "zip file extracted..."
    else: gz_filepath = ''; status = 'NA'
    return gz_filepath, status

def downloadCurrentVersionUI(filename,secondary_dir,file_type,root):
    continue_analysis = downloadCurrentVersion(filename,secondary_dir,file_type)
    if continue_analysis == 'no':
        import UI
        root.destroy(); UI.getUserParameters('no'); sys.exit()
    root.destroy()
        
def downloadCurrentVersion(filename,secondary_dir,file_type):
    import UI
    file_location_defaults = UI.importDefaultFileLocations()

    uds = file_location_defaults['url'] ### Get the location of the download site from Config/default-files.csv
    for ud in uds: url_dir = ud.Location() ### Only one entry
    
    dir = export.findParentDir(filename)  
    filename = export.findFilename(filename)
    url = url_dir+secondary_dir+'/'+filename
    
    file,status = download(url,dir,file_type); continue_analysis = 'yes'
    if 'Internet' in status:
        import UI
        print_out = "File could not be found on server\nor internet connection is unavailable."
        UI.WarningWindow(print_out,'WARNING!!!')
        continue_analysis = 'no'
    elif status == 'remove':
        try: os.remove(file) ### Not sure why this works now and not before
        except OSError: status = status
    return continue_analysis
            
def updateDBs(species,array_type):
    update_uniprot='no'; update_ensembl='no'; update_probeset_to_ensembl='no'; update_domain='no'; update_all='no'; update_miRs='no'
    proceed = 'no'

    while proceed == 'no':
        print "\n*****Select Species*****"
        print "1) Human"
        print "2) Mouse"
        print "3) Rat"
        print "4) Quit"
        inp = sys.stdin.readline(); inp = inp.strip()
        if inp  == '1': species = 'Hs'; proceed = 'yes'
        elif inp == '2': species = 'Mm'; proceed = 'yes'
        elif inp == '3': species = 'Rn'; proceed = 'yes'
        elif inp == '4': sys.exit()
        else: print "Sorry... that command is not an option\n"  

    proceed = 'no'
    while proceed == 'no':
        print "\n*****Select Array Type*****"
        print "1) Exon"
        print "2) AltMouse"
        print "3) Quit"
        inp = sys.stdin.readline(); inp = inp.strip()
        if inp  == '1': array_type = 'exon'; proceed = 'yes'
        elif inp == '2': array_type = 'AltMouse'; proceed = 'yes'
        elif inp == '3': sys.exit()
        else: print "Sorry... that command is not an option\n"

    proceed = 'no'        
    while proceed == 'no':
        print "\n*****Update AltAnalyze Database Options***** \n"
        print "Program Options: Updates for species",species,"and array type",array_type
        print "NOTE: some files will needed be manually updated before continuing (e.g., Affymetrix and Ensembl files). See ReadMe for more info."
        print "1) Update UniProt"
        print "2) Update Ensembl"
        print "3) Update Probeset-Ensembl associations"
        print "4) Update Protein/Domain associations"
        print "5) Update miRNA aligments"
        print "6) Update All"
        print "7) Quit"
        inp = sys.stdin.readline(); inp = inp.strip()
        if inp  == '1': update_uniprot = 'yes'; proceed = 'yes'
        elif inp == '2': update_ensembl = 'yes'; proceed = 'yes'
        elif inp == '3': update_probeset_to_ensembl = 'yes'; proceed = 'yes'
        elif inp == '4': update_domain = 'yes'; proceed = 'yes'
        elif inp == '5': update_miRs = 'yes'; proceed = 'yes'
        elif inp == '6': update_all = 'yes'; proceed = 'yes'
        elif inp == '7': sys.exit()
        else: print "Sorry... that command is not an option\n"

    proceed = 'no'
    while proceed == 'no':
        print "\n*****Update AltAnalyze Database Options*****"
        print "1) Force download of all files (whether present locally or not) "
        print "2) Use local or download if necessary"
        print "3) Quit"
        inp = sys.stdin.readline(); inp = inp.strip()
        if inp  == '1': force = 'yes'; proceed = 'yes'
        elif inp == '2': force = 'no'; proceed = 'yes'
        elif inp == '3': sys.exit()
        else: print "Sorry... that command is not an option\n"        

    if array_type == 'AltMouse':
        proceed = 'no'
        while proceed == 'no':
            print "\n*****Update AltAnalyze Database Options*****"
            print "1) New genomic build."
            print "2) Genomic build is identical to prior database build."
            print "3) Quit"
            inp = sys.stdin.readline(); inp = inp.strip()
            if inp  == '1': genomic_build = 'new'; proceed = 'yes'
            elif inp == '2': genomic_build = 'old'; proceed = 'yes'
            elif inp == '3': sys.exit()
            else: print "Sorry... that command is not an option\n"
        
    if update_all == 'yes':
        update_uniprot='yes'; update_ensembl='yes'; update_probeset_to_ensembl='yes'; update_domain='yes'; update_miRs = 'yes'
        
    if update_ensembl == 'yes':
        import EnsemblSQL; reload(EnsemblSQL)

        """ Used to grab all essential Ensembl annotations previously obtained via BioMart"""        
        configType = 'Advanced'; analysisType = 'AltAnalyzeDBs'; externalDBName = ''
        EnsemblSQL.buildEnsemblRelationalTablesFromSQL(species,configType,analysisType,externalDBName,force)
        
        """ Used to grab Ensembl-to-External gene associations"""
        configType = 'Basic'; analysisType = 'ExternalOnly'; externalDBName = 'Uniprot/SWISSPROT'
        EnsemblSQL.buildEnsemblRelationalTablesFromSQL(species,configType,analysisType,externalDBName,force)

    if update_uniprot == 'yes':            
        ###Might need to delete the existing versions of downloaded databases or force download
        buildUniProtFunctAnnotations(species,force)
                
    if update_probeset_to_ensembl == 'yes':
        if species == 'Mm' and array_type == 'AltMouse':
            buildAltMouseExonAnnotations(species,array_type,force,genomic_build)
        else: buildExonArrayExonAnnotations(species,force)

    if update_domain == 'yes':

        ### Get UCSC associations for all Ensembl linked genes (download databases if necessary)
        import UCSCImport
        mRNA_Type = 'mrna'; run_from_scratch = 'yes'
        export_all_associations = 'yes' ### YES only for protein prediction analysis
        UCSCImport.runUCSCEnsemblAssociations(species,mRNA_Type,export_all_associations,run_from_scratch,force)        

        if species == 'Mm' and array_type == 'AltMouse':
            """Imports and re-exports array-Ensembl annotations"""
            import JunctionArray
            null = JunctionArray.importArrayAnnotations(species,array_type); null={}
            """Performs probeset sequence aligment to Ensembl and UCSC transcripts. To do: Need to setup download if files missing"""
            import mRNASeqAlign
            mRNASeqAlign.alignProbesetsToTranscripts(species,array_type,force)
       
        import IdentifyAltIsoforms; run_seqcomp = 'no'
        IdentifyAltIsoforms.runProgram(species,array_type,force,run_seqcomp)

        import FeatureAlignment
        FeatureAlignment.findDomainsByGenomeCoordinates(species,array_type)
            
    if update_miRs == 'yes':
        import MatchMiRTargetPredictions; only_add_sequence_to_previous_results = 'no'
        MatchMiRTargetPredictions.runProgram(species,force,only_add_sequence_to_previous_results)

        if array_type == 'exon':        
            import ExonSeqModule
            stringency = 'strict'; process_microRNA_predictions = 'yes'; mir_source = 'multiple'
            ExonSeqModule.runProgram(species,array_type,process_microRNA_predictions,mir_source,stringency)
            stringency = 'lax'
            ExonSeqModule.runProgram(species,array_type,process_microRNA_predictions,mir_source,stringency)
        else:
            import JunctionSeqModule
            stringency = 'strict'; mir_source = 'multiple'
            JunctionSeqModule.runProgram(species,array_type,mir_source,stringency,force)
            stringency = 'lax'
            JunctionSeqModule.runProgram(species,array_type,mir_source,stringency,force)
            
if __name__ == '__main__':
    #buildUniProtFunctAnnotations('Hs',force='no')

    url='http://altanalyze.org/archiveDBs/LibraryFiles/MoGene-1_0-st-v1.r3.zip'
    dir='AltDatabase/affymetrix/LibraryFiles/'
    file_type = ''
    download(url,dir,file_type);kill
    
    species='Hs'; array_type = 'exon'
    updateDBs(species,array_type); sys.exit()
 
    