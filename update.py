#!/usr/local/bin/python2.6

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

"""This module contains generic methods for downloading and decompressing files designated
online files and coordinating specific database build operations for all AltAnalyze supported
gene ID systems with other program modules. """

import os
import sys
import unique
import string
import export
import traceback
try: from build_scripts import ExonArray
except Exception: null=[]
import traceback

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

def zipDirectory(dir):
    #http://www.testingreflections.com/node/view/8173
    import zipfile
    dir = filepath(dir); zip_file = dir+'.zip'
    p = string.split(dir,'/'); top=p[-1]
    zip = zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for f in files:
            fullpath = os.path.join(root, f)
            archive_name = os.path.join(top+archive_root, f)
            zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()
    return zip_file

def unzipFiles(filename,dir):
    import zipfile
    output_filepath = filepath(dir+'/'+filename)
    try:
        zfile = zipfile.ZipFile(output_filepath)
        for name in zfile.namelist():
            if name.endswith('/'):null=[] ### Don't need to export
            else:
                if 'EnsMart' in name and 'EnsMart' in dir:
                    dir = export.findParentDir(dir[:-1]) ### Remove EnsMart suffix directory
                try: outfile = export.ExportFile(filepath(dir+name))
                except Exception: outfile = export.ExportFile(filepath(dir+name[1:]))
                outfile.write(zfile.read(name)); outfile.close()
        #print 'Zip extracted to:',output_filepath
        status = 'completed'
    except Exception, e:
        try:
            ### Use the operating system's unzip if all else fails
            extracted_path = string.replace(output_filepath,'.zip','')
            try: os.remove(extracted_path) ### This is necessary, otherwise the empty file created above will require user authorization to delete
            except Exception: null=[]
            subprocessUnzip(dir,output_filepath)
            status = 'completed'
        except IOError:
            print e
            print 'WARNING!!!! The zip file',output_filepath,'does not appear to be a valid zip archive file or is currupt.'
            status = 'failed'
    return status

def subprocessUnzip(dir,output_filepath):
    import subprocess
    dir = filepath(dir)
    subprocess.Popen(["unzip", "-d", dir, output_filepath]).wait()
    
############## Update Databases ##############
def buildJunctionExonAnnotations(species,array_type,specific_array_type,force,genomic_build):
    ### Get UCSC associations (download databases if necessary)
    mRNA_Type = 'mrna'; run_from_scratch = 'yes'; force='no'
    export_all_associations = 'no' ### YES only for protein prediction analysis
    #buildUCSCAnnoationFiles(species,mRNA_Type,export_all_associations,run_from_scratch,force)

    ### Get genomic locations and initial annotations for exon sequences (exon pobesets and junctions)    
    from build_scripts import JunctionArray
    from build_scripts import JunctionArrayEnsemblRules
    """ The following functions:
    1) Extract transcript cluster-to-gene annotations
    2) Extract exon sequences for junctions and exon probesets from the Affymetrix annotation file (version 2.0),
    3) Map these sequences to Ensembl gene sequences (build specific) plus and minus 2KB, upstream and downstream
    4) Obtain AltAnalyze exon region annotations and obtain full-length exon sequences for each exon probeset
    5) Consoladate these into an Ensembl_probeset.txt file (rather than Ensembl_junction_probeset.txt) with junctions
       having a single probeset identifier.
    6) Determine which junctions and junction-exons represent recipricol junctions using:
       a) AltAnalyze identified recipricol junctions from Ensembl and UCSC and
       b) Affymetrix suggested recipricol junctions based on common exon cluster annotations, creating
          Mm_junction_comps_updated.txt.
       c) De novo comparison of all exon-junction region IDs for all junctions using the EnsemblImport method compareJunctions().
    """
    ### Steps 1-3
    JunctionArray.getJunctionExonLocations(species,array_type,specific_array_type)
    ### Step 4
    JunctionArrayEnsemblRules.getAnnotations(species,array_type,'yes',force)
    ### Step 5-6
    JunctionArray.identifyJunctionComps(species,array_type,specific_array_type)
    
    
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
    
    from build_scripts import JunctionArray
    from build_scripts import JunctionArrayEnsemblRules    
    rederive_exonseq = 'no'
    ### Only needs to be run once, to export exon sequence for AltMouse array the original (1 and 2 above)
    if rederive_exonseq == 'yes':
        import AltAnalyze
        from import_scripts import ExonAnnotate_module
        from build_scripts import ExonAnalyze_module
        agglomerate_inclusion_probesets = 'no'; onlyAnalyzeJunctions='no'
        probeset_annotations_file = "AltDatabase/"+species+"/"+array_type+"/"+"MASTER-probeset-transcript.txt"
        verifyFile(probeset_annotations_file,array_type) ### Will force download if missing
        exon_db={}; filtered_arrayids={};filter_status='no'
        constituitive_probeset_db,exon_db,genes_being_analyzed = AltAnalyze.importSplicingAnnotationDatabase(probeset_annotations_file,array_type,filtered_arrayids,filter_status)
        alt_junction_db,critical_exon_db,exon_dbase,exon_inclusion_db,exon_db = ExonAnnotate_module.identifyPutativeSpliceEvents(exon_db,constituitive_probeset_db,{},agglomerate_inclusion_probesets,onlyAnalyzeJunctions)
        ExonAnnotate_module.exportJunctionComparisons(alt_junction_db,critical_exon_db,exon_dbase)
        print "Finished exporting junctions used in AltMouse array comparisons."

        ExonAnalyze_module.exportAltMouseExonSequence()
        JunctionArray.reAnnotateCriticalExonSequences(species,array_type)

    ### Get UCSC associations (download databases if necessary)
    mRNA_Type = 'mrna'; run_from_scratch = 'yes'
    export_all_associations = 'no' ### YES only for protein prediction analysis
    buildUCSCAnnoationFiles(species,mRNA_Type,export_all_associations,run_from_scratch,force)

    reannotate_exon_seq = 'yes'
    print 'genomic_build', genomic_build
    if genomic_build == 'new':
        ### Need to run with every new genomic build (match up new coordinates
        print "Begining to derive exon sequence from new genomic build"
        JunctionArray.identifyCriticalExonLocations(species,array_type)
        reannotate_exon_seq = 'yes'
    JunctionArrayEnsemblRules.getAnnotations(species,array_type,reannotate_exon_seq,force)
    
    ### Download files required during AltAnalyze analysis but not during the database build process
    filename = "AltDatabase/"+species+"/"+array_type+"/"+"MASTER-probeset-transcript.txt"
    verifyFile(filename,array_type) ### Will force download if missing
    filename = "AltDatabase/"+species+'/'+ array_type+'/'+array_type+"_annotations.txt"
    verifyFile(filename,array_type) ### Will force download if missing
    
def buildExonArrayExonAnnotations(species, array_type, force):

    ### Get UCSC associations (download databases if necessary)
    mRNA_Type = 'mrna'; run_from_scratch = 'yes'
    export_all_associations = 'no' ### YES only for protein prediction analysis
    buildUCSCAnnoationFiles(species,mRNA_Type,export_all_associations,run_from_scratch,force)
    from build_scripts import ExonArrayEnsemblRules; reload(ExonArrayEnsemblRules)
    process_from_scratch='yes'
    constitutive_source='default'
    ### Build the databases and return the variables (not used here)
    source_biotype = 'mRNA'
    if array_type == 'gene': source_biotype = 'gene'
    probeset_db,annotate_db,constitutive_gene_db,splicing_analysis_db = ExonArrayEnsemblRules.getAnnotations(process_from_scratch,constitutive_source,source_biotype,species)

def getFTPData(ftp_server,subdir,filename_search):
    ### This is a more generic function for downloading FTP files based on input directories and a search term
    
    from ftplib import FTP
    print 'Connecting to',ftp_server
    ftp = FTP(ftp_server); ftp.login()
    ftp.cwd(subdir)
    
    ftpfilenames = []; ftp.dir(ftpfilenames.append); ftp.quit()
    matching=[]

    for line in ftpfilenames:
        line = string.split(line,' '); file_dir = line[-1]
        dir = 'ftp://'+ftp_server+subdir+'/'+file_dir
        #print dir
        if filename_search in dir and '.md5' not in dir: matching.append(dir)
    if len(matching)==1:
        return matching[0]
    elif len(matching)==0:
        print filename_search, 'not found at',ftp_server+subdir
        return string.replace(filename_search,'.gz','')
    else:
        return matching
    
def verifyFile(filename,server_folder):
    fn=filepath(filename); counts=0
    try:
        for line in open(fn,'rU').xreadlines():
            if len(line)>2: counts+=1  ### Needed for connection error files
            if counts>10: break
    except Exception:
        counts=0
    if server_folder == 'counts': ### Used if the file cannot be downloaded from http://www.altanalyze.org
        return counts
    elif counts == 0:
        if server_folder == None: server_folder = 'AltMouse'
        continue_analysis = downloadCurrentVersion(filename,server_folder,'')
        if continue_analysis == 'no' and 'nnot' not in filename:
            print 'The file:\n',filename, '\nis missing and cannot be found online. Please save to the designated directory or contact AltAnalyze support.'
    else:
        return counts

def getFileLocations(species_code,search_term):
    ### The three supplied variables determine the type of data to obtain from default-files.csv
    import UI
    file_location_defaults = UI.importDefaultFileLocations()
    sfl = file_location_defaults[search_term]
    
    for sf in sfl:
        if species_code in sf.Species(): filename = sf.Location()
    return filename
    
def buildUniProtFunctAnnotations(species,force):
    import UI
    file_location_defaults = UI.importDefaultFileLocations()
    """Identify the appropriate download location for the UniProt database for the selected species"""
    uis = file_location_defaults['UniProt']
    trembl_filename_url=''
    
    for ui in uis:
        if species in ui.Species(): uniprot_filename_url = ui.Location()
    species_codes = importSpeciesInfo(); species_full = species_codes[species].SpeciesName()
    
    from build_scripts import ExtractUniProtFunctAnnot; reload(ExtractUniProtFunctAnnot)
    ExtractUniProtFunctAnnot.runExtractUniProt(species,species_full,uniprot_filename_url,trembl_filename_url,force)

class SpeciesData:
    def __init__(self, abrev, species):
        self._abrev = abrev; self._species = species
    def SpeciesCode(self): return self._abrev
    def SpeciesName(self): return self._species
    def __repr__(self): return self.Report()
    
def importSpeciesInfo():
    filename = 'Config/species_all.txt'; x=0
    fn=filepath(filename); species_codes={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        abrev = t[0]; species = t[1]
        if x==0: x=1
        else:
            sd = SpeciesData(abrev,species)
            species_codes[abrev] = sd
    return species_codes
    
def downloadSuppressPrintOuts(url,dir,file_type):
    global Suppress_Printouts
    Suppress_Printouts = 'yes'
    return download(url,dir,file_type)
    
def download(url,dir,file_type):
    global suppress_printouts
    try: suppress_printouts = Suppress_Printouts
    except Exception: suppress_printouts = 'no'
    try: dp = download_protocol(url,dir,file_type); output_filepath, status  = dp.getStatus(); fp = output_filepath
    except Exception:
        #print traceback.format_exc()
        try:
            dir = unique.filepath(dir) ### Can result in the wrong filepath exported for AltDatabase RNA-Seq zip files (don't include by default)
            dp = download_protocol(url,dir,file_type); output_filepath, status  = dp.getStatus(); fp = output_filepath
        except Exception:
            output_filepath='failed'; status = "Internet connection not established. Re-establish and try again."
            fp = filepath(dir+url.split('/')[-1]) ### Remove this empty object if saved
    if 'Internet' not in status:
        if '.zip' in fp or '.gz' in fp or '.tar' in fp:
            #print "\nRemoving zip file:",fp
            try: os.remove(fp); status = 'removed'
            except Exception: null=[] ### Not sure why this error occurs since the file is not open
            #print "\nRemoving zip file:",string.replace(fp,'.gz','')
            if '.tar' in fp:
                try: os.remove(string.replace(fp,'.gz',''))
                except Exception: null=[]
    return output_filepath, status

class download_protocol:
    def __init__(self,url,dir,file_type):
        try: self.suppress = suppress_printouts
        except Exception: self.suppress = 'no'
        """Copy the contents of a file from a given URL to a local file."""
        filename = url.split('/')[-1]; self.status = ''
        #print [url, dir, file_type]
        #dir = unique.filepath(dir) ### Can screw up directory structures
        if file_type == None: file_type =''
        if len(file_type) == 2: filename, file_type = file_type ### Added this feature for when a file has an invalid filename
        output_filepath = unique.filepath(dir+filename, force='application-path')
        dir = export.findParentDir(output_filepath)
        output_filepath_object = export.createExportFile(output_filepath,dir[:-1])

        self.output_filepath = output_filepath
        
        if self.suppress == 'no':
            print "Downloading the following file:",filename,' ',
        self.original_increment = 5
        self.increment = 0
        import urllib
        from urllib import urlretrieve
        #if 'gene.txt.gz' in url: print [self.reporthookFunction];sys.exit()
        try:
            try: webfile, msg = urlretrieve(url,output_filepath,reporthook=self.reporthookFunction)
            except IOError:
                if 'Binary' in traceback.format_exc(): #IOError: [Errno ftp error] 200 Switching to Binary mode.
                    ### https://bugs.python.org/issue1067702 - some machines the socket doesn't close and causes an error - reload to close the socket
                    reload(urllib)
                    webfile, msg = urlretrieve(url,output_filepath,reporthook=self.reporthookFunction)
                    reload(urllib)
        except:
            print 'Unknown URL error encountered...'; forceURLError

        if self.suppress == 'no': print ''
        self.testFile()
        if self.suppress == 'no': print self.status

        if 'Internet' not in self.status:
            if '.zip' in filename:
                if self.suppress == 'no': print "Extracting zip file...",
                try: decompressZipStackOverflow(filename,dir); status = 'completed'
                except Exception:
                    #print 'Native unzip not present...trying python unzip methods...'
                    status = unzipFiles(filename,dir)
                    if status == 'failed': print 'zip extraction failed!'
                self.gz_filepath = filepath(output_filepath); self.status = 'remove'
                if self.suppress == 'no': print "zip file extracted"
            elif '.gz' in filename:
                self.gz_filepath = output_filepath
                if len(file_type)==0: extension = '.gz'
                else: extension = 'gz'
                decompressed_filepath = string.replace(self.gz_filepath,extension,file_type)
                ### Below code can be too memory intensive
                #file_size = os.path.getsize(output_filepath)
                #megabtyes = file_size/1000000.00
                #if megabtyes>5000: force_error ### force_error is an undefined variable which causes an exception
                import gzip; content = gzip.GzipFile(self.gz_filepath, 'rb')
                data = open(decompressed_filepath,'wb')
                #print "\nExtracting downloaded file:",self.gz_filepath
                import shutil; shutil.copyfileobj(content,data)
                # http://pythonicprose.blogspot.com/2009/10/python-extract-or-unzip-tar-file.html
                os.chdir(filepath(dir))
                if '.tar' in decompressed_filepath:
                    import tarfile
                    tfile = tarfile.open(decompressed_filepath)
                    tfile.extractall()
                    tfile.close()
                    tar_dir = string.replace(decompressed_filepath,'.tar','')
                self.status = 'remove'
            else: self.gz_filepath = ''; self.status = 'remove'
    def testFile(self):
        fn=filepath(self.output_filepath)
        try:
            for line in open(fn,'rU').xreadlines():
                if '!DOCTYPE' in line: self.status = "Internet connection not established. Re-establish and try again."
                break
        except Exception: null=[]
        
    def getStatus(self): return self.output_filepath, self.status
    def reporthookFunction(self, blocks_read, block_size, total_size):
        if not blocks_read:
            if self.suppress == 'no':
                print 'Connection opened. Downloading (be patient)'
        if total_size < 0:
            # Unknown size
            if self.suppress == 'no':
                print 'Read %d blocks' % blocks_read
        else:
            amount_read = blocks_read * block_size
            percent_read = ((amount_read)*1.00/total_size)*100
            if percent_read>self.increment:
                #print '%d%% downloaded' % self.increment
                if self.suppress == 'no':
                    print '*',
                self.increment += self.original_increment
            #print 'Read %d blocks, or %d/%d' % (blocks_read, amount_read, (amount_read/total_size)*100.000)

def createExportFile(new_file,dir):
    try:
        fn=filepath(new_file); file_var = open(fn,'w')
    except IOError:
        #print "IOError", fn
        fn = filepath(dir)
        try:
            os.mkdir(fn) ###Re-Create directory if deleted
            #print fn, 'written'
        except OSError: createExportDir(new_file,dir) ###Occurs if the parent directory is also missing
        fn=filepath(new_file); file_var = open(fn,'w')
    return file_var

def downloadCurrentVersionUI(filename,secondary_dir,file_type,root):
    continue_analysis = downloadCurrentVersion(filename,secondary_dir,file_type)
    if continue_analysis == 'no' and 'nnot' not in filename:
        import UI
        try: root.destroy(); UI.getUserParameters('no'); sys.exit()
        except Exception: sys.exit()
    try: root.destroy()
    except Exception(): null=[] ### Occurs when running from command-line
        
def downloadCurrentVersion(filename,secondary_dir,file_type):
    import UI
    file_location_defaults = UI.importDefaultFileLocations()

    ud = file_location_defaults['url'] ### Get the location of the download site from Config/default-files.csv
    url_dir = ud.Location() ### Only one entry
    
    dir = export.findParentDir(filename)
    dir = string.replace(dir,'hGlue','')  ### Used since the hGlue data is in a sub-directory
    filename = export.findFilename(filename)
    url = url_dir+secondary_dir+'/'+filename
    file,status = download(url,dir,file_type); continue_analysis = 'yes'
    if 'Internet' in status and 'nnot' not in filename: ### Exclude for Affymetrix annotation files
        print_out = "File:\n"+url+"\ncould not be found on the server or an internet connection is unavailable."
        if len(sys.argv)<2:
            try:
                UI.WarningWindow(print_out,'WARNING!!!')
                continue_analysis = 'no'
            except Exception:
                print 'cannot be downloaded';force_error
        else: print 'cannot be downloaded';force_error
    elif status == 'remove' and ('.zip' in file or '.tar' in file or '.gz' in file):
        try: os.remove(file) ### Not sure why this works now and not before
        except Exception: status = status
    return continue_analysis

def decompressZipStackOverflow(original_zip_file,dir):
    zip_file = filepath(dir+original_zip_file)
    #print 'Using OS native unzip software to extract'
    ###http://stackoverflow.com/questions/339053/how-do-you-unzip-very-large-files-in-python
    import zipfile
    import zlib
    src = open(zip_file,"rb")
    zf = zipfile.ZipFile(src)
    for m in zf.infolist():
        # Examine the header
        #print m.filename, m.header_offset, m.compress_size, repr(m.extra), repr(m.comment)
        src.seek(m.header_offset)
        src.read(30) # Good to use struct to unpack this.
        nm= src.read(len(m.filename))
        if len(m.extra) > 0: ex= src.read(len(m.extra))
        if len(m.comment) > 0: cm=src.read(len(m.comment)) 

        # Build a decompression object
        decomp= zlib.decompressobj(-15)

        # This can be done with a loop reading blocks
        out=open(filepath(dir+m.filename), "wb")
        result=decomp.decompress(src.read(m.compress_size))
        out.write(result); result=decomp.flush()
        out.write(result); out.close()
    zf.close()
    src.close()
    
def buildUCSCAnnoationFiles(species,mRNA_Type,export_all_associations,run_from_scratch,force):
    ### Test whether files already exist and if not downloads/builds them
    if export_all_associations == 'no': ### Applies to EnsemblImport.py analyses
        filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_mrna.txt'
    else: ### Applies to the file used for Domain-level analyses
        filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_COMPLETE-mrna.txt'
    counts = verifyFile(filename,'counts')
    if counts<9:
        from build_scripts import UCSCImport
        try: UCSCImport.runUCSCEnsemblAssociations(species,mRNA_Type,export_all_associations,run_from_scratch,force)
        except Exception: UCSCImport.exportNullDatabases(species) ### used for species not supported by UCSC
        
def executeParameters(species,array_type,force,genomic_build,update_uniprot,update_ensembl,update_probeset_to_ensembl,update_domain,update_miRs,update_all,update_miR_seq,ensembl_version):    
    if '|' in array_type: array_type, specific_array_type = string.split(array_type,'|') ### To destinguish between array sub-types, like the HJAY and hGlue
    else: specific_array_type = array_type
    
    if update_all == 'yes':
        update_uniprot='yes'; update_ensembl='yes'; update_probeset_to_ensembl='yes'; update_domain='yes'; update_miRs = 'yes'
        
    if update_ensembl == 'yes':
        from build_scripts import EnsemblSQL; reload(EnsemblSQL)

        """ Used to grab all essential Ensembl annotations previously obtained via BioMart"""        
        configType = 'Advanced'; analysisType = 'AltAnalyzeDBs'; externalDBName = ''
        EnsemblSQL.buildEnsemblRelationalTablesFromSQL(species,configType,analysisType,externalDBName,ensembl_version,force)
        
        """ Used to grab Ensembl-to-External gene associations"""
        configType = 'Basic'; analysisType = 'ExternalOnly'; externalDBName = 'Uniprot/SWISSPROT'
        EnsemblSQL.buildEnsemblRelationalTablesFromSQL(species,configType,analysisType,externalDBName,ensembl_version,force)
        
        """ Used to grab Ensembl full gene sequence plus promoter and 3'UTRs """
        if array_type == 'AltMouse' or array_type == 'junction' or array_type == 'RNASeq':
            EnsemblSQL.getFullGeneSequences(ensembl_version,species)
            
    if update_uniprot == 'yes':            
        ###Might need to delete the existing versions of downloaded databases or force download
        buildUniProtFunctAnnotations(species,force)
                
    if update_probeset_to_ensembl == 'yes':
        if species == 'Mm' and array_type == 'AltMouse':
            buildAltMouseExonAnnotations(species,array_type,force,genomic_build)
        elif array_type == 'junction':
            buildJunctionExonAnnotations(species,array_type,specific_array_type,force,genomic_build)
        elif array_type == 'RNASeq':
            import RNASeq; test_status = 'no'; data_type = 'mRNA'
            RNASeq.getEnsemblAssociations(species,data_type,test_status,force)
        else: buildExonArrayExonAnnotations(species,array_type,force)

    if update_domain == 'yes':
        if array_type == 'RNASeq':
            only_rely_on_coordinate_mapping = True ### This will provide more accurate results as many junctions have missing sequences
        else:
            only_rely_on_coordinate_mapping = False

        from build_scripts import FeatureAlignment
        from build_scripts import JunctionArray
        from build_scripts import mRNASeqAlign
        from build_scripts import IdentifyAltIsoforms
        
        ### Get UCSC associations for all Ensembl linked genes (download databases if necessary)        if species == 'Mm' and array_type == 'AltMouse':
        mRNA_Type = 'mrna'; run_from_scratch = 'yes'
        export_all_associations = 'yes' ### YES only for protein prediction analysis
        buildUCSCAnnoationFiles(species,mRNA_Type,export_all_associations,run_from_scratch,force)

        if (species == 'Mm' and array_type == 'AltMouse'):
            """Imports and re-exports array-Ensembl annotations"""
            null = JunctionArray.importArrayAnnotations(species,array_type); null={}
        if (species == 'Mm' and array_type == 'AltMouse') or array_type == 'junction' or array_type == 'RNASeq':
            if only_rely_on_coordinate_mapping == False:
                """Performs probeset sequence aligment to Ensembl and UCSC transcripts. To do: Need to setup download if files missing"""
                analysis_type = 'reciprocal'
                mRNASeqAlign.alignProbesetsToTranscripts(species,array_type,analysis_type,force)
    
        run_seqcomp = 'no'
        if only_rely_on_coordinate_mapping == False:
            IdentifyAltIsoforms.runProgram(species,array_type,'null',force,run_seqcomp)
            FeatureAlignment.findDomainsByGenomeCoordinates(species,array_type,'null')
        
        if array_type == 'junction' or array_type == 'RNASeq':
            if only_rely_on_coordinate_mapping == False:
                ### For junction probeset sequences from mRNASeqAlign(), find and assess alternative proteins - export to the folder 'junction'
                mRNASeqAlign.alignProbesetsToTranscripts(species,array_type,'single',force)
                IdentifyAltIsoforms.runProgram(species,array_type,'junction',force,run_seqcomp)
                FeatureAlignment.findDomainsByGenomeCoordinates(species,array_type,'junction')
                ### For exon probesets (and junction exons) align and assess alternative proteins - export to the folder 'exon'
                IdentifyAltIsoforms.runProgram(species,array_type,'exon',force,run_seqcomp)
                FeatureAlignment.findDomainsByGenomeCoordinates(species,array_type,'exon') # not needed
            
            """ Repeat above with CoordinateBasedMatching = True """ 
            ### Peform coordinate based junction mapping to transcripts (requires certain sequence files built in IdentifyAltIosofmrs)
            analysis_type = 'reciprocal'
            mRNASeqAlign.alignProbesetsToTranscripts(species,array_type,analysis_type,force,CoordinateBasedMatching = True)
            IdentifyAltIsoforms.runProgram(species,array_type,'null',force,run_seqcomp)
            FeatureAlignment.findDomainsByGenomeCoordinates(species,array_type,'null')
            mRNASeqAlign.alignProbesetsToTranscripts(species,array_type,'single',force,CoordinateBasedMatching = True)
            IdentifyAltIsoforms.runProgram(species,array_type,'junction',force,run_seqcomp)
            FeatureAlignment.findDomainsByGenomeCoordinates(species,array_type,'junction')
            IdentifyAltIsoforms.runProgram(species,array_type,'exon',force,run_seqcomp)
            if array_type == 'RNASeq':
                JunctionArray.combineExonJunctionAnnotations(species,array_type)
                
    if update_miRs == 'yes':
        if update_miR_seq == 'yes':
            from build_scripts import MatchMiRTargetPredictions; only_add_sequence_to_previous_results = 'no'
            MatchMiRTargetPredictions.runProgram(species,force,only_add_sequence_to_previous_results)

        if array_type == 'exon' or array_type == 'gene':        
            from build_scripts import ExonSeqModule
            stringency = 'strict'; process_microRNA_predictions = 'yes'; mir_source = 'multiple'
            ExonSeqModule.runProgram(species,array_type,process_microRNA_predictions,mir_source,stringency)
            stringency = 'lax'
            ExonSeqModule.runProgram(species,array_type,process_microRNA_predictions,mir_source,stringency)
            ExonArray.exportMetaProbesets(array_type,species) ### Export metaprobesets for this build
        else:
            from build_scripts  import JunctionSeqModule
            stringency = 'strict'; mir_source = 'multiple'
            JunctionSeqModule.runProgram(species,array_type,mir_source,stringency,force)
            stringency = 'lax'
            JunctionSeqModule.runProgram(species,array_type,mir_source,stringency,force)

    if array_type == 'junction':
        try:
            from build_scripts import JunctionArray; from build_scripts import JunctionArrayEnsemblRules
            JunctionArray.filterForCriticalExons(species,array_type)
            JunctionArray.overRideExonEntriesWithJunctions(species,array_type)
            JunctionArrayEnsemblRules.annotateJunctionIDsAsExon(species,array_type)
            ExonArray.exportMetaProbesets(array_type,species) ### Export metaprobesets for this build
        except IOError: print 'No built junction files to analyze';sys.exit()
    if array_type == 'RNASeq' and (species == 'Hs' or species == 'Mm' or species == 'Rn'):
        from build_scripts import JunctionArray; from build_scripts import JunctionArrayEnsemblRules
        try: JunctionArrayEnsemblRules.annotateJunctionIDsAsExon(species,array_type)
        except IOError: print 'No Ensembl_exons.txt file to analyze';sys.exit()
    
    try:
        filename = 'AltDatabase/'+species+'/SequenceData/miRBS-combined_gene-targets.txt'; ef=filepath(filename)
        er = string.replace(ef,species+'/SequenceData/miRBS-combined_gene-targets.txt','ensembl/'+species+'/'+species+'_microRNA-Ensembl.txt')
        import shutil; shutil.copyfile(ef,er)
    except Exception: null=[]
    if array_type != 'RNASeq':
        ### Get the probeset-probe relationships from online - needed for FIRMA analysis
        filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_probeset-probes.txt'
        if array_type == 'junction' and 'lue' in specific_array_type:
            server_folder = 'junction/hGlue'
            verifyFile(filename,server_folder) ### Will force download if missing
            verifyFile('AltDatabase/'+species+'/'+array_type+'/platform.txt',server_folder) ### Will force download if missing
        elif array_type != 'AltMouse': verifyFile(filename,array_type) ### Will force download if missing
        if (array_type == 'exon' or array_type == 'AltMouse') and species != 'Rn':
            try:
                ### Available for select exon-arrays and AltMouse
                probeset_to_remove_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_probes_to_remove.txt'
                verifyFile(probeset_to_remove_file,array_type)
            except Exception: null=[]
    
def temp(array_type,species):
    specific_array_type = 'hGlue'
    ExonArray.exportMetaProbesets(array_type,species) ### Export metaprobesets for this build
    if array_type != 'RNASeq':
        ### Get the probeset-probe relationships from online - needed for FIRMA analysis
        filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_probeset-probes.txt'
        if array_type == 'junction' and 'lue' in specific_array_type:
            server_folder = 'junction/hGlue'
            verifyFile(filename,server_folder) ### Will force download if missing
            verifyFile('AltDatabase/'+species+'/'+array_type+'/platform.txt',server_folder) ### Will force download if missing
        elif array_type != 'AltMouse': verifyFile(filename,array_type) ### Will force download if missing
        if (array_type == 'exon' or array_type == 'AltMouse') and species != 'Rn':
            try:
                ### Available for select exon-arrays and AltMouse
                probeset_to_remove_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_probes_to_remove.txt'
                verifyFile(probeset_to_remove_file,array_type)
            except Exception: null=[]
            
if __name__ == '__main__':
    
    dir = '/home/nsalomonis/software/AltAnalyze_v.2.0.5-Py/AltDatabase/EnsMart65/miRBS/'; filename = 'UTR_Sequences.txt.zip'
    array_type = 'gene'; species = 'Mm'
    filename = "AltDatabase/"+species+'/'+ array_type+'/'+array_type+"_annotations.txt"
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_probeset-probes.txt'
    #verifyFile(filename,array_type);sys.exit()
    #temp(array_type,species); sys.exit()
    
    #unzipFiles(filename,dir); kill
    #zipDirectory('Rn/EnsMart49');kill
    #unzipFiles('Rn.zip', 'AltDatabaseNoVersion/');kill    
    #filename = 'http://altanalyze.org/archiveDBs/LibraryFiles/Mouse430_2.zip'
    #filename = 'AltDatabase/affymetrix/LibraryFiles/Mouse430_2.zip'
    filename = 'AltDatabase/Mm_RNASeq.zip'; dir = 'AltDatabase/updated/EnsMart72'
    #downloadCurrentVersionUI(filename,dir,'','')
    import update
    dp = update.download_protocol('ftp://ftp.ensembl.org/pub/release-72/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.72.cdna.all.fa.gz','AltDatabase/Hs/SequenceData/','');sys.exit()
    #dp = update.download_protocol('ftp://ftp.ensembl.org/pub/release-72/mysql/macaca_mulatta_core_72_10/gene.txt.gz','AltDatabase/ensembl/Ma/EnsemblSQL/','');sys.exit()

    #kill
    #target_folder = 'Databases/Ci'; zipDirectory(target_folder)
    filename = 'Mm_RNASeq.zip'; dir = 'AltDatabase/EnsMart72/'
    unzipFiles(filename, dir);sys.exit()
    #buildUniProtFunctAnnotations('Hs',force='no')

    species = 'Hs'; array_type = 'junction'; force = 'yes'; run_seqcomp = 'no'
    
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_probeset-probes.txt'
    filename = string.replace(filename,'junction','junction/hGlue') ### Due this also for the hGlue array
    #verifyFile(filename,'junction/hGlue')
    #print filename; sys.exit()
    #from build_scripts import IdentifyAltIsoforms
    #IdentifyAltIsoforms.runProgram(species,array_type,'exon',force,run_seqcomp); sys.exit()
            
    import UI
    #date = UI.TimeStamp(); file_type = ('wikipathways_'+date+'.tab','.txt')
    url  ='http://www.wikipathways.org/wpi/pathway_content_flatfile.php?output=tab'
    output = 'BuildDBs/wikipathways/'
    file_type = ''
    url = 'http://www.genmapp.org/go_elite/Databases/EnsMart56/Cs.zip'
    output = 'Databases/'
    url = 'http://www.altanalyze.org/archiveDBs/Cytoscape/cytoscape.tar.gz'
    output = ''
    dp = download_protocol(url,output,file_type);sys.exit()
    fln,status = download(url,output,file_type)
    print status;sys.exit()
        
    url='http://altanalyze.org/archiveDBs/LibraryFiles/MoEx-1_0-st-v1.r2.antigenomic.bgp.gz'
    dir='AltDatabase/affymetrix/LibraryFiles/'
    file_type = ''
    download(url,dir,file_type)
    
    species='Hs'; array_type = 'junction'
    
    #"""
    array_type = ['RNASeq']; species = ['Mm']; force = 'yes'; ensembl_version = '60'
    update_uniprot='no'; update_ensembl='no'; update_probeset_to_ensembl='yes'; update_domain='yes'; update_all='no'; update_miRs='yes'
    proceed = 'yes'; genomic_build = 'old'; update_miR_seq = 'yes'

    for specific_species in species:
        for platform_name in array_type:
            if platform_name == 'AltMouse' and specific_species == 'Mm': proceed = 'yes'
            elif platform_name == 'exon' or platform_name == 'gene' or platform_name == 'junction': proceed = 'yes'
            else: proceed = 'no'
            if proceed == 'yes':
                print "Analyzing", specific_species, platform_name
                executeParameters(specific_species,platform_name,force,genomic_build,update_uniprot,update_ensembl,update_probeset_to_ensembl,update_domain,update_miRs,update_all,update_miR_seq,ensembl_version)

    #"""    
