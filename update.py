###update
#Copyright 2005-2008 J. Davide Gladstone Institutes, San Francisco California
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
def buildAltMouseExonAnnotations():
    """Code required to:
    1) Extract out Affymetrix provided exon sequence (probeset sequence extracted from "probeset_sequence_reversed.txt", derived
       directly from the Affymetrix AltMouse probe annotation file), from the "SEQUENCE-transcript-dbase.txt" (built using
       dump-chip1 .gff sequence and AltMerge-Peptide Informatics script "sequence_analysis_AltMouse_refseq.py").
    2) Once exported, grab full length exon sequences using exon/intron coordinates matches to full-length gene sequences with 2kb
       flanking sequence to efficiently predict microRNA binding site exclusion (reAnnotateCriticalExonSequences) and later for
       coordinate mapping to get exons aligning with UCSC annotated splicing annotations and exons. This sequence data replaced
       the previous file (don't need to re-run this - see commented out code below for reference).
    3) Match the updated exon sequences to the most recent genomic coordinates and build the exact equivalent of the exon array
       Mm_Ensembl_probeset.txt database (same structure and ExonArrayEnsemblRules.py code). This involves running EnsemblImport.
    This code should be run before the exon array location build code since the "Mm_Ensembl_probeset.txt" is created and then re-
    written as "Mm_AltMouse_Ensembl_probeset.txt".
    """

    import JunctionArray
    import JunctionArrayEnsemblRules    
    species = 'Mm'; array_type = 'AltMouse'
    rederive_exonseq == 'no'
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
    
    ### Need to run with every new genomic build (match up new coordinates
    JunctionArray.identifyCriticalExonLocations(species,array_type)
    JunctionArrayEnsemblRules.getAnnotations(species,array_type)

def buildExonArrayExonAnnotations(species):
    import ExonArrayEnsemblRules
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
    if '.gz' in filename:
        gz_filepath = output_filepath
        decompressed_filepath = string.replace(gz_filepath,'gz',file_type)
        import gzip
        zfile = gzip.GzipFile(gz_filepath)
        content = zfile.read()
        zfile.close()
        data = open(decompressed_filepath,'w')
        print "\nExtracting downloaded file:",gz_filepath
        data.write(content); data.close()
        print "\nRemoving gz file:",gz_filepath
        try: os.remove(gz_filepath); status = 'removed'
        except OSError: status = 'not-removed' ### Not sure why this error occurs since the file is not open
    else: gz_filepath = ''; status = 'removed'
    return gz_filepath, status

def downloadCurrentVersion(filename,dir,secondary_dir,file_type):
    import UI
    file_location_defaults = UI.importDefaultFileLocations()    
    url_dir = file_location_defaults['url'].Location() ### Get the location of the download site from Config/default-files.csv
    url = url_dir+secondary_dir+'/'+filename
    file,status = download(url,dir,file_type)

def updateDBs(species,array_type):
    proceed = 'no'
    while proceed == 'no':
        print "\n*****Update AltAnalyze Database Options***** \n"
        print "Program Options: Updates for species",species,"and array type",array_type
        print "NOTE: some files will needed be manually updated before continuing (e.g., Affymetrix and Ensembl files). See ReadMe for more info."
        print "1) Download all external databases from source websites (UCSC, NCBI, and UniProt). "
        print "2) Perform step 1 and immediately build databases (must have other files downloaded),"
        print "3) Build databases will local files present."
        print "4) Quit"
        inp = sys.stdin.readline(); inp = inp.strip()
        if inp  == '1': download_new_dbs = 'yes'; rebuild_altanalyze_dbs = 'no'; proceed = 'yes'
        elif inp == '2': download_new_dbs = 'yes'; rebuild_altanalyze_dbs = 'yes'; proceed = 'yes'
        elif inp == '3': download_new_dbs = 'no'; rebuild_altanalyze_dbs = 'yes'; proceed = 'yes'
        elif inp == '4': sys.exit()
        else: print "Sorry... that command is not an option\n"

    force = 'no'
    if download_new_dbs == 'yes':
        force = 'no'
        
    if rebuild_altanalyze_dbs == 'yes':            
        ###Might need to delete the existing versions of downloaded databases or force download
        #buildUniProtFunctAnnotations(species,force)

        if species == 'Mm' and array_type == 'AltMouse':
            buildAltMouseExonAnnotations()
        else: buildExonArrayExonAnnotations(species)
            
if __name__ == '__main__':
    #buildUniProtFunctAnnotations('Hs',force='no')
    species='Hs'; array_type = 'exon'
    updateDBs(species,array_type)
    