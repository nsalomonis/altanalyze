###EnsemblSQL
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

"""This module contains methods for reading the Ensembl SQL FTP directory structure to
identify appropriate species and systems for download for various versions of Ensembl,
download the necessary database files to reconstruct any BioMart annotation files and
determine genomic coordinates for the start and end positions of protein domains."""

try: clearall()
except NameError: null = [] ### Occurs when re-running the script to clear all global variables

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os.path
import unique
import export
import traceback
import update; reload(update)

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

def getChrGeneOnly(Species,configType,ensembl_version,force):
    global species; species = Species
    global rewrite_existing; rewrite_existing = 'yes'
    print 'Downloading Ensembl flat files for parsing from Ensembl SQL FTP server...'
    
    global ensembl_build
    ensembl_sql_dir, ensembl_sql_description_dir = getEnsemblSQLDir(ensembl_version)
    print ensembl_sql_dir
    ensembl_build = string.split(ensembl_sql_dir,'core')[-1][:-1]
    sql_file_db,sql_group_db = importEnsemblSQLInfo(configType) ###Import the Config file with the files and fields to parse from the downloaded SQL files

    filtered_sql_group_db={} ### Get only these tables
    filtered_sql_group_db['Primary'] = ['gene.txt','gene_stable_id.txt','seq_region.txt']
        
    filtered_sql_group_db['Description'] = [ensembl_sql_description_dir]
    sq = EnsemblSQLInfo(ensembl_sql_description_dir, 'EnsemblSQLDescriptions', 'Description', '', '', '')
    sql_file_db['Primary',ensembl_sql_description_dir] = sq        
    output_dir = 'AltDatabase/ensembl/'+species+'/EnsemblSQL/'
    importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,filtered_sql_group_db,sql_file_db,output_dir,'Primary',force) ###Download and import the Ensembl SQL files
    
    program_type,database_dir = unique.whatProgramIsThis()
    if program_type == 'AltAnalyze':
        parent_dir = 'AltDatabase/goelite/'
    else:
        parent_dir = 'Databases/'
    output_dir = parent_dir+species+'/uid-gene/Ensembl-chr.txt'
            
    headers = ['Ensembl Gene ID', 'Chromosome']
    values_list=[]
    for geneid in gene_db:
        gi = gene_db[geneid]
        try: ens_gene = gene_db[geneid].StableId()
        except Exception: ens_gene = gene_stable_id_db[geneid].StableId() 
        seq_region_id = gi.SeqRegionId()
        chr = seq_region_db[seq_region_id].Name()
        values = [ens_gene,chr]
        values_list.append(values)
    exportEnsemblTable(values_list,headers,output_dir)
    
def getGeneTranscriptOnly(Species,configType,ensembl_version,force):
    global species; species = Species
    global rewrite_existing; rewrite_existing = 'yes'
    print 'Downloading Ensembl flat files for parsing from Ensembl SQL FTP server...'
    
    global ensembl_build
    ensembl_sql_dir, ensembl_sql_description_dir = getEnsemblSQLDir(ensembl_version)
    print ensembl_sql_dir
    ensembl_build = string.split(ensembl_sql_dir,'core')[-1][:-1]
    sql_file_db,sql_group_db = importEnsemblSQLInfo(configType) ###Import the Config file with the files and fields to parse from the downloaded SQL files

    filtered_sql_group_db={} ### Get only these tables
    filtered_sql_group_db['Primary'] = ['gene.txt','gene_stable_id.txt','transcript.txt','transcript_stable_id.txt']
        
    filtered_sql_group_db['Description'] = [ensembl_sql_description_dir]
    sq = EnsemblSQLInfo(ensembl_sql_description_dir, 'EnsemblSQLDescriptions', 'Description', '', '', '')
    sql_file_db['Primary',ensembl_sql_description_dir] = sq        
    output_dir = 'AltDatabase/ensembl/'+species+'/EnsemblSQL/'
    additionalFilter = ['transcript','gene']
    importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,filtered_sql_group_db,sql_file_db,output_dir,'Primary',force) ###Download and import the Ensembl SQL files
    
    program_type,database_dir = unique.whatProgramIsThis()
    if program_type == 'AltAnalyze':
        parent_dir = 'AltDatabase/goelite/'
    else:
        parent_dir = 'Databases/'
    output_dir = parent_dir+species+'/uid-gene/Ensembl-EnsTranscript.txt'
            
    headers = ['Ensembl Gene ID', 'Ensembl Transcript ID']
    values_list=[]
    for transcript_id in transcript_db:
        ti = transcript_db[transcript_id]
        try: ens_gene = gene_db[ti.GeneId()].StableId()
        except Exception: ens_gene = gene_stable_id_db[ti.GeneId()].StableId()
        try: ens_transcript = transcript_db[transcript_id].StableId()
        except Exception: ens_transcript = transcript_stable_id_db[transcript_id].StableId()
        values = [ens_gene,ens_transcript]
        values_list.append(values)
    exportEnsemblTable(values_list,headers,output_dir)
    
def getEnsemblSQLDir(ensembl_version):
    if 'Plus' in ensembl_version:
        ensembl_version = string.replace(ensembl_version,'Plus','')
    if 'EnsMart' in ensembl_version:
        ensembl_version = string.replace(ensembl_version, 'EnsMart','')
    original_version = ensembl_version
    if 'Plant' in ensembl_version:
        ensembl_version = string.replace(ensembl_version, 'Plant','')
    if 'Bacteria' in ensembl_version:
        ensembl_version = string.replace(ensembl_version, 'Bacteria','')
    if 'Fungi' in ensembl_version:
        ensembl_version = string.replace(ensembl_version, 'Fungi','')
    try:
        check = int(ensembl_version)    
        import UI; UI.exportDBversion('EnsMart'+ensembl_version)
        ensembl_version = 'release-'+ensembl_version
    except Exception: print 'EnsMart version name is incorrect (e.g., should be "60")'; sys.exit()
    
    import UI; species_names = UI.getSpeciesInfo()
    try:
        species_full = species_names[species]
        ens_species = string.replace(string.lower(species_full),' ','_')
    except Exception:
        ens_species = species
        species_full = species

    try: child_dirs, ensembl_species, ensembl_versions = getCurrentEnsemblSpecies(original_version)
    except Exception:
        print "\nPlease try a different version. This one does not appear to be valid."
        print traceback.format_exc();sys.exit()
    ensembl_sql_dir,ensembl_sql_description_dir = child_dirs[species_full]
    return ensembl_sql_dir, ensembl_sql_description_dir

def buildEnsemblRelationalTablesFromSQL(Species,configType,analysisType,externalDBName,ensembl_version,force,buildCommand=None):
    import UI; import update; global external_xref_key_db; global species; species = Species
    global system_synonym_db; system_synonym_db={} ### Currently only used by GO-Elite
    global rewrite_existing; rewrite_existing = 'yes'; global all_external_ids; all_external_ids={}; global added_systems; added_systems={}
    original_version = ensembl_version
    print 'Downloading Ensembl flat files for parsing from Ensembl SQL FTP server...'
    ### Get version directories for Ensembl

    global ensembl_build
    ensembl_sql_dir, ensembl_sql_description_dir = getEnsemblSQLDir(ensembl_version)
    ensembl_build = string.split(ensembl_sql_dir,'core')[-1][:-1]

    sql_file_db,sql_group_db = importEnsemblSQLInfo(configType) ###Import the Config file with the files and fields to parse from the downloaded SQL files
    sql_group_db['Description'] = [ensembl_sql_description_dir]
    sq = EnsemblSQLInfo(ensembl_sql_description_dir, 'EnsemblSQLDescriptions', 'Description', '', '', '')
    sql_file_db['Primary',ensembl_sql_description_dir] = sq
    output_dir = 'AltDatabase/ensembl/'+species+'/EnsemblSQL/'
    importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'Primary',force) ###Download and import the Ensembl SQL files

    if analysisType != 'ExternalOnly':
        ### Build and export the basic Ensembl gene annotation table
        xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'Xref',force)
        buildEnsemblGeneAnnotationTable(species,xref_db)

    if analysisType == 'ExternalOnly':
        ###Export data for Ensembl-External gene system
        buildFilterDBForExternalDB(externalDBName)
        xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'Xref',force)
        external_xref_key_db = xref_db; #resetExternalFilterDB()
        object_xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'Object-Xref',force)
        export_status='yes';output_id_type='Gene'
        buildEnsemblExternalDBRelationshipTable(externalDBName,xref_db,object_xref_db,output_id_type,species)

    if analysisType == 'AltAnalyzeDBs':
        ###Export data for Ensembl-External gene system
        buildTranscriptStructureTable()
        if buildCommand == 'exon':
            return None
        ###Export transcript biotype annotations (e.g., NMD)
        exportTranscriptBioType()
        
        ###Get Interpro AC display name
        buildFilterDBForExternalDB('Interpro')
        xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'Xref',force)
        getDomainGenomicCoordinates(species,xref_db)

        if force == 'yes':
            ### Download Transcript seqeunces
            getEnsemblTranscriptSequences(original_version,species)
     
def getEnsemblTranscriptSequences(ensembl_version,species,restrictTo=None):
    ### Download Seqeunce data
    import UI; species_names = UI.getSpeciesInfo()
    species_full = species_names[species]
    ens_species = string.replace(string.lower(species_full),' ','_')
    
    if 'release-' not in ensembl_version:
        ensembl_version = 'release-'+ensembl_version
    
    if restrictTo == None or restrictTo == 'protein':
        dirtype = 'fasta/'+ens_species+'/pep'
        if 'Plant' in ensembl_version or 'Bacteria' in ensembl_version or 'Fungi' in ensembl_version:
            ensembl_protseq_dir = getCurrentEnsemblGenomesSequences(ensembl_version,dirtype,ens_species)
        else:
            ensembl_protseq_dir = getCurrentEnsemblSequences(ensembl_version,dirtype,ens_species)

        output_dir = 'AltDatabase/ensembl/'+species + '/'
        gz_filepath, status = update.download(ensembl_protseq_dir,output_dir,'')
        if status == 'not-removed':
            try: os.remove(gz_filepath) ### Not sure why this works now and not before
            except OSError: status = status

    if restrictTo == None or restrictTo == 'cDNA':
        dirtype = 'fasta/'+ens_species+'/cdna'
        if 'Plant' in ensembl_version or 'Bacteria' in ensembl_version or 'Fungi' in ensembl_version:
            ensembl_cdnaseq_dir = getCurrentEnsemblGenomesSequences(ensembl_version,dirtype,ens_species)
        else:
            ensembl_cdnaseq_dir = getCurrentEnsemblSequences(ensembl_version,dirtype,ens_species)
            
        output_dir = 'AltDatabase/'+species + '/SequenceData/'
        #print [ensembl_cdnaseq_dir]
        gz_filepath, status = update.download(ensembl_cdnaseq_dir,output_dir,'')
        if status == 'not-removed':
            try: os.remove(gz_filepath) ### Not sure why this works now and not before
            except OSError: status = status
                
def getFullGeneSequences(ensembl_version,species):
    import UI; species_names = UI.getSpeciesInfo()
    try:
        species_full = species_names[species]
        ens_species = string.replace(string.lower(species_full),' ','_')
    except Exception:
        ens_species = species
        species_full = species

    if 'EnsMart' in ensembl_version:
        ensembl_version = string.replace(ensembl_version, 'EnsMart','')
    if 'release-' not in ensembl_version:
        ensembl_version = 'release-'+ensembl_version
        
    dirtype = 'fasta/'+ens_species+'/dna'
    if 'Plant' in ensembl_version or 'Bacteria' in ensembl_version or 'Fungi' in ensembl_version:
        ensembl_dnaseq_dirs = getCurrentEnsemblGenomesSequences(ensembl_version,dirtype,ens_species)
    else:
        ensembl_dnaseq_dirs = getCurrentEnsemblSequences(ensembl_version,dirtype,ens_species)

    output_dir = 'AltDatabase/'+species + '/SequenceData/chr/'
    global dna
    dna = export.ExportFile(output_dir+species+'_gene-seq-2000_flank.fa')
    print 'Writing gene sequence file plus 2000 flanking base-pairs'
        
    from build_scripts import EnsemblImport
    chr_gene_db,location_gene_db = EnsemblImport.getEnsemblGeneLocations(species,'','key_by_chromosome')
    
    for (chr_file,chr_url) in ensembl_dnaseq_dirs:
        #if 'MT' in chr_file:
        if 'toplevel' in chr_file:
            gz_filepath, status = update.download(chr_url,output_dir,'')
            dna_db = divideUpTopLevelChrFASTA(output_dir+chr_file[:-3],output_dir,'store') ### Would work nice, but too file intensive
            print len(dna_db), 'scaffolds read into memory'
            for chr in chr_gene_db:
                if chr in dna_db:
                    #print chr,'found'
                    parseChrFASTA(dna_db[chr],chr,chr_gene_db[chr],location_gene_db)
                else: print chr,'not found'
            try: os.remove(gz_filepath) ### Not sure why this works now and not before
            except Exception: null=[]
            try: os.remove(gz_filepath[:-3]) ### Delete the chromosome sequence file
            except Exception: null=[]
        else:   
            chr=string.split(chr_file,'.')[-3]
            if chr in chr_gene_db:
                gz_filepath, status = update.download(chr_url,output_dir,'')
                parseChrFASTA(output_dir+chr_file[:-3],chr,chr_gene_db[chr],location_gene_db)
                try: os.remove(gz_filepath) ### Not sure why this works now and not before
                except Exception: null=[]
                try: os.remove(gz_filepath[:-3]) ### Delete the chromosome sequence file
                except Exception: null=[]
    dna.close()

def divideUpTopLevelChrFASTA(filename,output_dir,action):
    """ When individual chromosome files are not present, parse the TopLevel file which contains segments divided by scaffold rather
    than individual chromosomes. Since the Scaffolds are reported rather than the chromosomes for genes in the other Ensembl annotation
    files, the scaffolds can serve as a surrogate for individual chromosome files when parsed out."""
    fn=filepath(filename); fasta=[]; sdna=None
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if '>' == data[0]:
            ### Save values from previous scaffold
            if sdna != None:
                if action == 'write':
                    for ln in fasta: sdna.write(ln)
                    sdna.close()
                else: sdna[scaffold_id] = [scaffold_line]+fasta
                #if scaffold_id == 'scaffold_753': print scaffold_id, scaffold_line, fasta;kill
                fasta=[]
            else: sdna = {} ### performed once
            ### Create new export object
            scaffold_id=string.split(data[1:],' ')[0]; scaffold_line = line
            if action == 'write':
                export_dir = string.replace(filename,'toplevel','chromosome.'+scaffold_id)
                sdna = export.ExportFile(export_dir)
                sdna.write(line)
        else: fasta.append(line)

    if action == 'write':
        for ln in fasta: sdna.write(ln)
        sdna.close()
    else: sdna[scaffold_id] = [scaffold_line]+fasta
    return sdna

def parseChrFASTA(filename,chr,chr_gene_locations,location_gene_db):
    """ This function parses FASTA formatted chromosomal sequence, storing only sequence needed to get the next gene plus any buffer seqence, in memory.
    Note: not straight forward to follow initially, as some tricks are needed to get the gene plus buffer sequence at the ends of the chromosome"""
    genes_to_be_examined={}
    for (ch,cs,ce) in location_gene_db:
        gene,strand = location_gene_db[ch,cs,ce]
        if ch==chr:
            genes_to_be_examined[gene]=[ch,cs,ce]
        
    from build_scripts import EnsemblImport

    if ('/' in filename) or ('\\' in filename): ### Occurs when parsing a toplevel chromosome sequence file
        print "Parsing chromosome %s sequence..." % chr
        fn=filepath(filename)
        sequence_data = open(fn,'rU').xreadlines()
        readtype = 'file'
    else:
        """The scaffold architecture of some species sequence files creates a challenge to integrate non-chromosomal sequence,
        however, we still want to analyze it in the same way as chromosomal sequence. One alternative is to write out each scaffold
        it's own fasta file. The problem with this is that it creates dozens of thousands of files which is time consuming and messy
        - see divideUpTopLevelChrFASTA(). To resolve this, we instead store the sequence in a dictionary as sequence lines, but read
        in the sequence lines just like they are directly read from the fasta sequence, allowing both pipelines to work the same
        (rather than maintain two complex functions that do the same thing)."""
        scaffold_line_db = filename 
        sequence_data = scaffold_line_db
        #print len(scaffold_line_db)
        readtype = 'dictionary'
        
    chr_gene_locations.sort()
    cs=chr_gene_locations[0][0]; ce=chr_gene_locations[0][1]; buffer=2000; failed=0; gene_count=0; terminate_chr_analysis='no'
    max_ce = chr_gene_locations[-1][1]; genes_analyzed={}
    gene,strand = location_gene_db[chr,cs,ce]
    x=0; sequence=''; running_seq_length=0; tr=0 ### tr is the total number of nucleotides removed (only store sequence following the last gene)
    for line in sequence_data:
        data = cleanUpLine(line)
        if x==0:
            #>MT dna:chromosome chromosome:GRCh37:MT:1:16569:1
            max_chr_length = int(string.split(data,':')[-2])-70
            x=1
        else:
            iterate = 'yes' ### Since the buffer region can itself contain multiple genes, we may need to iterature through the last stored sequence set for multiple genes
            sequence += data; running_seq_length+=len(data)
            if (ce+buffer)>=max_chr_length: target_seq_length = max_chr_length
            else: target_seq_length = ce+buffer
            if running_seq_length>=target_seq_length:
                while iterate == 'yes':
                    internal_buffer = buffer
                    adj_start = cs-tr-buffer-1; adj_end = ce-tr+buffer; original_adj_start=adj_start
                    if adj_start<0: original_adj_start = abs(adj_start); internal_buffer=buffer-(adj_start*-1); adj_start=0
                    try:
                        gene_seq = sequence[adj_start:adj_end]
                        ### Add "N" to the sequence, before and after, if the gene starts to close the chromosome end for the buffer
                        if adj_start==0: gene_seq=original_adj_start*'N'+gene_seq
                        elif len(gene_seq) != (adj_end-adj_start): gene_seq+=((adj_end-adj_start)-len(gene_seq))*'N'
                        start=str(cs-buffer); end=str(ce+buffer)
                        ### write this out
                        if strand=='-': gene_seq=EnsemblImport.reverse_orientation(gene_seq)
                        header = string.join(['>'+gene,chr,str(cs),str(ce)],'|') #['>'+gene,chr,start,end]
                        #print header, cs, ce, strand,adj_start,adj_end,len(sequence),len(gene_seq)
                        #print x,[gene_seq]
                        if gene not in genes_analyzed: ### Indicates something is iterating where it shouldn't be, but doesn't seem to create a signficant data handling issue
                            dna.write(header+'\n'); dna.write(gene_seq+'\n'); x+=1
                            genes_analyzed[gene]=[]
                    except KeyError: failed+=1; kill ### gene is too close to chromosome end to add buffer
                    seq_start = adj_start; tr = cs-internal_buffer-1; sequence=sequence[seq_start:]; gene_count+=1
                    try:
                        cs=chr_gene_locations[gene_count][0]; ce=chr_gene_locations[gene_count][1]
                        gene,strand = location_gene_db[chr,cs,ce]
                    except IndexError: terminate_chr_analysis = 'yes'; iterate='no'; break ### no more genes to examine on this chromosome
                    if len(data)<60: iterate='yes' ### Forces re-iteration through the last stored sequence set
                    else: iterate='no'
                if terminate_chr_analysis == 'yes': break
    if ('/' in filename) or ('\\' in filename):
        print "Sequence for",len(genes_analyzed),"out of",len(genes_to_be_examined),"genes exported..." 
        """
        for gene in genes_to_be_examined:
            if gene not in genes_analyzed:
                print gene, genes_to_be_examined[gene]
                sys.exit()"""
    elif len(genes_analyzed) != len(genes_to_be_examined):
        print len(genes_to_be_examined)-len(genes_analyzed), 'not found for',chr
    
def buildTranscriptStructureTable():
    ### Function for building Ensembl transcription annotation file
    temp_values_list = []; output_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    headers = ['Ensembl Gene ID', 'Chromosome', 'Strand', 'Exon Start (bp)', 'Exon End (bp)', 'Ensembl Exon ID', 'Constitutive Exon', 'Ensembl Transcript ID']
    for eti in exon_transcript_db:
        exonid = eti.ExonId(); transcript_id = eti.TranscriptId(); rank = eti.Rank()
        ei = exon_db[exonid]; seq_region_start = ei.SeqRegionStart(); seq_region_end = ei.SeqRegionEnd()
        try: constitutive_call = str(ei.IsConstitutive())
        except Exception: constitutive_call = '0'
        try: ens_exon = exon_db[exonid].StableId()
        except Exception: ens_exon = exon_stable_id_db[exonid].StableId()
        ti = transcript_db[transcript_id]; geneid = ti.GeneId()
        try: ens_transcript = transcript_db[transcript_id].StableId()
        except Exception: ens_transcript = transcript_stable_id_db[transcript_id].StableId()
        gi = gene_db[geneid]; seq_region_id = gi.SeqRegionId(); strand = gi.SeqRegionStrand()
        chr = seq_region_db[seq_region_id].Name()
        try: ens_gene = gene_db[geneid].StableId()
        except Exception: ens_gene = gene_stable_id_db[geneid].StableId()   
        values = [geneid,transcript_id,rank,ens_gene,chr,str(strand),str(seq_region_start),str(seq_region_end),ens_exon,constitutive_call,ens_transcript]
        temp_values_list.append(values)
        
    values_list=[]; temp_values_list.sort() ###Make sure the gene, transcripts and exons are grouped and ranked appropriately
    for values in temp_values_list: values_list.append(values[3:])
    temp_values_list=[]
    exportEnsemblTable(values_list,headers,output_dir)

def exportTranscriptBioType():
    ### Function for extracting biotype annotations for transcripts
    transcript_protein_id={}
    for protein_id in translation_db:
        ti = translation_db[protein_id]; transcript_id = ti.TranscriptId()
        transcript_protein_id[transcript_id] = protein_id
        
    values_list = []; output_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-biotypes.txt'
    headers = ['Ensembl Gene ID','Ensembl Translation ID', 'Biotype']
    for transcript_id in transcript_db:
        ti = transcript_db[transcript_id]
        
        try: ens_transcript = transcript_db[transcript_id].StableId()
        except Exception: ens_transcript = transcript_stable_id_db[transcript_id].StableId()
        try:
            try: protein_id = transcript_protein_id[transcript_id]; ens_protein = translation_db[protein_id].StableId()
            except Exception: protein_id = transcript_protein_id[transcript_id]; ens_protein = translation_stable_id_db[protein_id].StableId()
        except Exception: ens_protein = ens_transcript+'-PEP'
        geneid = ti.GeneId(); geneid = ti.GeneId()
        try: ens_gene = gene_db[geneid].StableId()
        except Exception: ens_gene = gene_stable_id_db[geneid].StableId()
        values = [ens_gene,ens_protein,ti.Biotype()]
        values_list.append(values)
    exportEnsemblTable(values_list,headers,output_dir)
    
def aaselect(nt_length):
    ### Convert to protein length
    if (float(nt_length)/3) == (int(nt_length)/3):
        return nt_length/3
    else:
        return int(string.split(str(float(nt_length)/3),'.')[0])+1
    
def getDomainGenomicCoordinates(species,xref_db):
    ### Get all transcripts relative to genes to export gene, transcript and protein ID relationships
    transcript_protein_id={}            
    for protein_id in translation_db:
        ti = translation_db[protein_id]; transcript_id = ti.TranscriptId()
        transcript_protein_id[transcript_id] = protein_id
        
    output_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_Protein_'+ensembl_build+'.tab'
    headers = ['Gene', 'Trans', 'Protein']; values_list=[]
    for transcript_id in transcript_db:
        geneid = transcript_db[transcript_id].GeneId()
        try: ens_gene = gene_db[geneid].StableId()
        except Exception: ens_gene = gene_stable_id_db[geneid].StableId()
        try:
            try: protein_id = transcript_protein_id[transcript_id]; ens_protein = translation_db[protein_id].StableId()
            except Exception: protein_id = transcript_protein_id[transcript_id]; ens_protein = translation_stable_id_db[protein_id].StableId()
        except KeyError: ens_protein = ''
        try: ens_transcript = transcript_db[transcript_id].StableId()
        except Exception: ens_transcript = transcript_stable_id_db[transcript_id].StableId()
        values_list.append([ens_gene,ens_transcript,ens_protein])
    exportEnsemblTable(values_list,headers,output_dir)
    
    ### Function for building Ensembl transcription annotation file
    exon_transcript_db2={} ### exon_transcript_db has data stored as a list, so store as a ranked db
    for eti in exon_transcript_db:
        transcript_id = eti.TranscriptId(); rank = eti.Rank()
        try: exon_transcript_db2[transcript_id].append((rank,eti))
        except KeyError: exon_transcript_db2[transcript_id] = [(rank,eti)]

    interpro_annotation_db={}
    for xref_id in xref_db:
        interprot_id = xref_db[xref_id].DbprimaryAcc(); display_label = xref_db[xref_id].DisplayLabel()
        interpro_annotation_db[interprot_id] = display_label

    ### Get the protein coding positions for each exon, to later determine the genomic position of non-InterPro features (way downstream)
    ### This code was adapted from the following paragraph to get InterPro locations (this is simpiler)
    output_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_ProteinCoordinates_build'+ensembl_build+'.tab'
    headers = ['protienID', 'exonID', 'AA_NT_Start', 'AA_NT_Stop', 'Genomic_Start', 'Genomic_Stop']; values_list=[]
    cds_output_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_EnsemblTranscriptCDSPositions.tab'
    cds_headers = ['transcriptID', 'CDS_Start', 'CDS_Stop']; cds_values_list=[]
    protein_coding_exon_db={}; protein_length_db={}
    ### Get the bp position (relative to the exon not genomic) for each transcript, where protein translation begins and ends. Store with the exon data.
    for protein_id in translation_db:
        ti = translation_db[protein_id]; transcript_id = ti.TranscriptId()
        seq_start = ti.SeqStart(); start_exon_id = ti.StartExonId(); seq_end = ti.SeqEnd(); end_exon_id = ti.EndExonId()
        eti_list = exon_transcript_db2[transcript_id]; eti_list.sort()
        ### Get info for exporting exon protein coordinate data
        try: ens_protein = translation_db[protein_id].StableId()
        except Exception: ens_protein = translation_stable_id_db[protein_id].StableId()
        try: ens_transcript = transcript_db[transcript_id].StableId()
        except Exception: ens_transcript = transcript_stable_id_db[transcript_id].StableId()
        #if ens_protein == 'ENSDARP00000087122':
        cumulative_exon_length = 0
        tis = transcript_db[transcript_id]; geneid = tis.GeneId(); gi = gene_db[geneid]; strand = gi.SeqRegionStrand() #Get strand
        if strand == '-1': start_correction = -3; end_correction = 2
        else: start_correction = 0; end_correction = 0
        translation_pos1=0; translation_pos2=0 ### Indicate amino acid positions in nt space, begining and ending in the exon
        coding_exons = []; ce=0
        for (rank,eti) in eti_list:
            exonid = eti.ExonId()
            try: ens_exon = exon_db[exonid].StableId()
            except Exception: ens_exon = exon_stable_id_db[exonid].StableId()
            ei = exon_db[exonid]; genomic_exon_start = ei.SeqRegionStart(); genomic_exon_end = ei.SeqRegionEnd()
            exon_length = abs(genomic_exon_end-genomic_exon_start)+1
            #print exonid,start_exon_id,end_exon_id
            if exonid == start_exon_id: ### Thus we are in the start exon for this transcript
                coding_bp_in_exon = exon_length - seq_start; eti.setCodingBpInExon(coding_bp_in_exon) ### -1 since coding can't start at 0, but rather 1 at a minimum
                #print 'start',genomic_exon_start,genomic_exon_end,exon_length,seq_start;kill
                coding_exons.append(eti); ce+=1; translation_pos2=coding_bp_in_exon+1
                if strand == -1:
                    genomic_translation_start = ei.SeqRegionStart()+coding_bp_in_exon+start_correction
                    genomic_exon_end = ei.SeqRegionStart()+end_correction
                else:               
                    genomic_translation_start = ei.SeqRegionEnd()-coding_bp_in_exon
                    genomic_exon_end = ei.SeqRegionEnd()
                #print 'trans1:',float(translation_pos1)/3, float(translation_pos2)/3
                values_list.append([ens_protein,ens_exon,1,aaselect(translation_pos2),genomic_translation_start,genomic_exon_end])
                cds_start = seq_start+cumulative_exon_length ### start position in this exon plus last exon cumulative length
            elif exonid == end_exon_id: ### Thus we are in the end exon for this transcript        
                coding_bp_in_exon = seq_end; eti.setCodingBpInExon(coding_bp_in_exon)
                coding_exons.append(eti); ce = 0
                translation_pos1=translation_pos2+1
                translation_pos2=translation_pos1+coding_bp_in_exon-4
                if strand == -1:
                    genomic_translation_start = ei.SeqRegionEnd()+end_correction+start_correction
                    genomic_exon_end = ei.SeqRegionEnd()-coding_bp_in_exon+end_correction
                else:               
                    genomic_translation_start = ei.SeqRegionStart()
                    genomic_exon_end = ei.SeqRegionStart()+coding_bp_in_exon
                #print 'trans1:',float(translation_pos1)/3, float(translation_pos2)/3
                values_list.append([ens_protein,ens_exon,aaselect(translation_pos1),aaselect(translation_pos2),genomic_translation_start,genomic_exon_end])
                cds_end = seq_end+cumulative_exon_length ### start position in this exon plus last exon cumulative length
                cds_values_list.append([ens_transcript,cds_start,cds_end])
                #ens_exon = exon_stable_id_db[exonid].StableId()  
            elif ce != 0: ###If we are in coding exons
                eti.setCodingBpInExon(exon_length)
                coding_exons.append(eti)
                translation_pos1=translation_pos2+1 ### 1 nucleotide difference from the last exon position
                translation_pos2=translation_pos1+exon_length-1
                if strand == -1:
                    genomic_translation_start = ei.SeqRegionEnd()+start_correction
                    genomic_exon_end = ei.SeqRegionStart()+end_correction
                else:               
                    genomic_translation_start = ei.SeqRegionStart()
                    genomic_exon_end = ei.SeqRegionEnd()
                #print 'trans1:',float(translation_pos1)/3,float(translation_pos2)/3
                values_list.append([ens_protein,ens_exon,aaselect(translation_pos1),aaselect(translation_pos2),genomic_translation_start,genomic_exon_end])
            cumulative_exon_length += exon_length
            #ti = transcript_db[transcript_id]; geneid = ti.GeneId(); ens_gene = gene_stable_id_db[geneid].StableId()
            #ens_exon = exon_stable_id_db[exonid].StableId()
            #if ens_gene == 'ENSG00000106785':
            #if ens_exon == 'ENSE00001381077':
                #print exon_length, seq_start, coding_bp_in_exon;kill
                #print 'start',ens_gene,rank,len(eti_list),exonid,ens_exon,start_exon_id,end_exon_id,genomic_exon_start,genomic_exon_end,exon_length,seq_start,coding_bp_in_exon,seq_end
        protein_coding_exon_db[protein_id] = coding_exons
    print 'Exporting protein-to-exon genomic position translations',len(values_list)
    exportEnsemblTable(values_list,headers,output_dir)
    exportEnsemblTable(cds_values_list,cds_headers,cds_output_dir)
    #sys.exit()
    
    ### Using the exon coding positions, determine InterPro coding genomic locations
    output_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_ProteinFeatures_build'+ensembl_build+'.tab'
    headers = ['ID', 'AA_Start', 'AA_Stop', 'Start', 'Stop', 'Name', 'Interpro', 'Description']
    interprot_match=0; values_list=[]
    #print len(protein_feature_db),len(translation_db),len(interpro_db),len(interpro_annotation_db);sys.exit()
    for protein_feature_id in protein_feature_db:
        pfi = protein_feature_db[protein_feature_id]; protein_id = pfi.TranslationId(); evalue = pfi.Evalue()
        try: ens_protein = translation_db[protein_id].StableId()
        except Exception: ens_protein = translation_stable_id_db[protein_id].StableId()
        seq_start = pfi.SeqStart(); seq_end = pfi.SeqEnd(); hit_id = pfi.HitId() ### hit_id is the domain accession which maps to 'id' in interpro
        seq_start = seq_start*3-2; seq_end = seq_end*3 ###convert to transcript basepair positions
        coding_exons = protein_coding_exon_db[protein_id]
        cumulative_coding_length = 0; last_exon_cumulative_coding_length = 1
        ti = translation_db[protein_id]; transcript_id = ti.TranscriptId() ### Get transcript ID, to look up strand
        ti = transcript_db[transcript_id]; geneid = ti.GeneId(); gi = gene_db[geneid]; strand = gi.SeqRegionStrand() #Get strand
        if strand == '-1': start_correction = -3; end_correction = 2
        else: start_correction = 0; end_correction = 0
        if hit_id in interpro_db and evalue<1: ### Only analyze domain-level features that correspond to known protein feature IDs
            interpro_ac = interpro_db[hit_id].InterproAc(); interprot_match+=1
            if interpro_ac in interpro_annotation_db:
                interpro_name = interpro_annotation_db[interpro_ac] ###Built this annotation database using the xref database
                #print interpro_name,ens_protein,seq_start,seq_end
                genomic_domain_start = 0; genomic_domain_end = 0; domain_in_first_exon = 'no'; non_coding_seq_len = 0
                for eti in coding_exons:
                    domain_in_first_exon = 'no'
                    coding_bp_in_exon = eti.CodingBpInExon(); cumulative_coding_length += coding_bp_in_exon
                    if seq_start <= cumulative_coding_length and seq_start >= last_exon_cumulative_coding_length: ### Thus, domain starts in this exon
                        var = 'start'
                        exonid = eti.ExonId(); ei = exon_db[exonid]
                        if abs(ei.SeqRegionEnd()-ei.SeqRegionStart()+1) != coding_bp_in_exon:
                            ### Can occur in the first exon
                            if last_exon_cumulative_coding_length<2:
                                domain_in_first_exon = 'yes'
                                non_coding_seq_len = abs(ei.SeqRegionEnd()-ei.SeqRegionStart()+1) - coding_bp_in_exon
                                genomic_bp_exon_offset = seq_start - last_exon_cumulative_coding_length + non_coding_seq_len
                            else:
                                genomic_bp_exon_offset = seq_start - last_exon_cumulative_coding_length
                        else: genomic_bp_exon_offset = seq_start - last_exon_cumulative_coding_length
                        if strand == -1:
                            genomic_exon_start = ei.SeqRegionEnd() ### This needs to be reversed if reverse strand
                            genomic_domain_start = genomic_exon_start-genomic_bp_exon_offset+start_correction ### This is what we want! (minus 3 so that we start at the first bp of that codon, not the first of the next (don't count the starting coding as 3bp)
                        else:
                            genomic_exon_start = ei.SeqRegionStart()                        
                            genomic_domain_start = genomic_bp_exon_offset+genomic_exon_start+start_correction ### This is what we want! (minus 3 so that we start at the first bp of that codon, not the first of the next (don't count the starting coding as 3bp)
                        #print genomic_exon_start,last_exon_cumulative_coding_length,genomic_domain_start,genomic_bp_exon_offset;kill
                        #pfi.setGenomicStart(genomic_domain_start)
                    if seq_end <= cumulative_coding_length and seq_end >= last_exon_cumulative_coding_length: ### Thus, domain ends in this exon
                        var = 'end'
                        exonid = eti.ExonId(); ei = exon_db[exonid]
                        genomic_bp_exon_offset = seq_end - last_exon_cumulative_coding_length
                        if (abs(ei.SeqRegionEnd()-ei.SeqRegionStart()+1) != coding_bp_in_exon) and domain_in_first_exon == 'yes': genomic_bp_exon_offset += non_coding_seq_len ### If the domain starts/ends in the first exon
                        if strand == -1:
                            genomic_exon_start = ei.SeqRegionEnd() ### This needs to be reversed if reverse strand      
                            genomic_domain_end = genomic_exon_start-genomic_bp_exon_offset+end_correction ### This is what we want!
                        else:
                            genomic_exon_start = ei.SeqRegionStart()
                            genomic_domain_end = genomic_bp_exon_offset+genomic_exon_start+end_correction ### This is what we want!
                        #pfi.setGenomicEnd(genomic_domain_end)
                    #"""
                    #ens_exon = exon_stable_id_db[eti.ExonId()].StableId()  
                    #if cumulative_coding_length == seq_end and strand == -1 and seq_start == 1 and domain_in_first_exon == 'yes':
                        #print interpro_name,protein_id,eti.ExonId(),ens_protein,ens_exon,seq_end,genomic_domain_start,genomic_domain_end;kill
                    
                    if ens_protein == 'ENSMUSP00000097740':
                        print interpro_name, genomic_domain_start, genomic_domain_end, last_exon_cumulative_coding_length, seq_end, seq_start, non_coding_seq_len
                        #print 'coding_bp_in_exon, cumulative_coding_length, genomic_bp_exon_offset',exon_db[exonid].StableId(), coding_bp_in_exon, cumulative_coding_length, genomic_bp_exon_offset
                        if var == 'start':
                            print interpro_name, var,genomic_exon_start,genomic_bp_exon_offset,start_correction, ei.SeqRegionStart(), ei.SeqRegionEnd()
                        else:
                            print interpro_name, var,genomic_exon_start,genomic_bp_exon_offset,end_correction, ei.SeqRegionStart(), ei.SeqRegionEnd()
                        #print 'exon',ens_exon,eti.ExonId(),ei.SeqRegionStart(), ei.SeqRegionEnd()#"""
                        #print non_coding_seq_len, domain_in_first_exon, coding_bp_in_exon, genomic_exon_start, genomic_domain_start, genomic_bp_exon_offset, start_correction
                        #print seq_start, interpro_name,seq_end, cumulative_coding_length,last_exon_cumulative_coding_length, genomic_domain_start, genomic_domain_end, ei.SeqRegionStart(), ei.SeqRegionEnd()
                    last_exon_cumulative_coding_length = cumulative_coding_length + 1
                if genomic_domain_start !=0 and genomic_domain_end !=0:
                    values = [ens_protein,(seq_start/3)+1,seq_end/3,genomic_domain_start,genomic_domain_end,hit_id,interpro_ac,interpro_name]
                    values_list.append(values)

    print 'interprot_matches:',interprot_match
    exportEnsemblTable(values_list,headers,output_dir)

def buildFilterDBForExternalDB(externalDBName):
    ### Generic function for pulling out any specific type of Xref without storing the whole database in memory
    global external_filter_db
    external_filter_db={}; key_filter_db={} ### Reset key_filter_db, otherwise this would cause importPrimaryEnsemblSQLTables to import the wrong entries
    for external_db_id in external_db_db:
        db_name = external_db_db[external_db_id].DbName()
        if db_name == externalDBName: external_filter_db[external_db_id]=[]

def buildFilterDBForArrayDB(externalDBName):
    ### Generic function for pulling out any specific type of Xref without storing the whole database in memory
    global external_filter_db
    external_filter_db={}; key_filter_db={} ### Reset key_filter_db, otherwise this would cause importPrimaryEnsemblSQLTables to import the wrong entries
    for array_chip_id in array_chip_db:
        array_id = array_chip_db[array_chip_id].ArrayID()
        try:
            name = array_db[array_id].Name()
            vendor = array_db[array_id].Vendor()
            format = array_db[array_id].Format()
            if name == externalDBName and (format != 'TILED' and '\\N' not in vendor):
                external_filter_db[array_chip_id]=[vendor]
        except KeyError: null=[]

def resetExternalFilterDB():
    global external_filter_db; external_filter_db={}
    external_filter_db=external_filter_db
    
def buildEnsemblGeneAnnotationTable(species,xref_db):
    ### Get Definitions and Symbol for each Ensembl gene ID
    values_list = []; output_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations_simple.txt'
    headers = ['Ensembl Gene ID','Description','Gene name']
    for geneid in gene_db:
        gi = gene_db[geneid]; display_xref_id = gi.DisplayXrefId(); description = gi.Description()
        if description == '\\N': description = ''
        try: symbol = xref_db[display_xref_id].DisplayLabel()
        except KeyError: symbol = ''
        try: ens_gene = gene_db[geneid].StableId()
        except Exception: ens_gene = gene_stable_id_db[geneid].StableId()
        values = [ens_gene,description,symbol]
        values_list.append(values)
        ###Also, create a filter_db that allows for the creation of Ensembl-external gene db tables
    exportEnsemblTable(values_list,headers,output_dir)
    
def buildEnsemblExternalDBRelationshipTable(external_system,xref_db,object_xref_db,output_id_type,species,writeToGOEliteFolder=False):
    ### Get xref annotations (e.g., external geneID) for a set of Ensembl IDs (e.g. gene IDs)
    external_system = string.replace(external_system,'/','-') ###Some external relationship systems have a slash in them which will create unwanted directories
    
    program_type,database_dir = unique.whatProgramIsThis()
    if program_type == 'AltAnalyze' and writeToGOEliteFolder == False:
        output_dir = 'AltDatabase/'+external_system+'/'+species+'/'+species+'_Ensembl-'+external_system+'.txt'
    else:
        if program_type == 'AltAnalyze':
            parent_dir = 'AltDatabase/goelite' ### Build directly with the AltAnalyze database
        elif 'over-write previous' in overwrite_previous: parent_dir = 'Databases'
        else: parent_dir = 'NewDatabases'
        if external_system == 'GO':
            output_dir = parent_dir+'/'+species+'/gene-go/Ensembl-GeneOntology.txt'
        elif 'meta' in external_system:
            output_dir = parent_dir+'/'+species+'/uid-gene/Ensembl_EntrezGene-meta.txt'
        elif external_system in system_synonym_db:
            system_name = system_synonym_db[external_system]
            output_dir = parent_dir+'/'+species+'/uid-gene/Ensembl-'+system_name+'.txt'
        elif 'Uniprot' in external_system: ### Needed for AltAnalyze
            output_dir = parent_dir+'/'+species+'/uid-gene/Ensembl-Uniprot.txt'
        else:
            output_dir = parent_dir+'/'+species+'/uid-gene/Ensembl-'+external_system+'.txt'
        version_info = species+' Ensembl relationships downloaded from EnsemblSQL server, build '+ensembl_build

    try: exportVersionInfo(output_dir,version_info)
    except Exception: null=[]
    
    headers = ['Ensembl ID',external_system + ' ID']; id_type_db={}; index=0
    id_type_db={}; gene_relationship_db={}; transcript_relationship_db={}; translation_relationship_db={}
    for xref_id in object_xref_db:
        for ox in object_xref_db[xref_id]:
            ens_numeric_id = ox.EnsemblId()
            try:
                index+=1
                try: dbprimary_acc = xref_db[xref_id].DbprimaryAcc()
                except Exception: print external_system, xref_id, xref_db[xref_id], len(xref_db);sys.exit() 
                external_db_id = xref_db[xref_id].ExternalDbId()
                ens_object_type = ox.EnsemblObjectType()
                #print len(xref_db),len(object_xref_db),len(external_filter_db),xref_id,dbprimary_acc,external_db_id,ens_object_type;kill
                #if '13065181' == external_db_id or '13065181' == xref_id or '13065181' == dbprimary_acc: print 'blah';kill
                if external_db_id in external_filter_db: ###Make sure other gene systems are not included
                    ### For testing determine the most likely system linked to by the xref (e.g. gene, transcript or translation ID).
                    try: id_type_db[ox.EnsemblObjectType()]+=1
                    except KeyError: id_type_db[ox.EnsemblObjectType()]=1
                        
                    if ens_object_type == 'Gene': gene_relationship_db[ens_numeric_id,dbprimary_acc]=[]
                    if ens_object_type == 'Transcript': transcript_relationship_db[ens_numeric_id,dbprimary_acc]=[]
                    if ens_object_type == 'Translation': translation_relationship_db[ens_numeric_id,dbprimary_acc]=[]
            except KeyError: null=[]       
    ids=['ID types linked to '+external_system+' are: ']
    for id_type in id_type_db: ### Tells us which ID types are most connected to a given external reference ID (don't often know)
        ids+=[str(id_type),': ',str(id_type_db[id_type]),'\t']
    ids = string.join(ids,''); print ids
    
    values_list = convertBetweenEnsemblIDTypes(output_id_type,transcript_relationship_db,translation_relationship_db,gene_relationship_db)
    if 'meta' in external_system:
        values_list2=[]
        for values in values_list: values_list2.append([values[1],values[0]]); values_list2.append(values)
        values_list = values_list2
    if len(values_list)>0:
        exportEnsemblTable(values_list,headers,output_dir)
        added_systems[external_system]=[]
    return len(values_list)

def exportVersionInfo(dir,version_info):
    dirs = string.split(dir,'/')
    dir = string.join(dirs[:-1],'/') ### Remove the filename
    data = export.ExportFile(dir+'/Ensembl_version.txt')
    data.write(version_info+'\n')
    
def convertBetweenEnsemblIDTypes(output_id_type,transcript_relationship_db,translation_relationship_db,gene_relationship_db):
    ### Starting with gene, transcript or protein relationships to an xref, convert to a single bio-type
    values_list=[]
    
    ### Get all proteins relative to transcripts
    transcript_to_protein_db = {}
    for protein_id in translation_db:
        transcript_id = translation_db[protein_id].TranscriptId()
        transcript_to_protein_db[transcript_id] = protein_id

    ### Get all transcripts relative to genes
    gene_to_transcript_db={}
    for transcript_id in transcript_db:
        geneid = transcript_db[transcript_id].GeneId()
        try: gene_to_transcript_db[geneid].append(transcript_id)
        except KeyError: gene_to_transcript_db[geneid] = [transcript_id]

    for (ens_numeric_id,dbprimary_acc) in transcript_relationship_db:
        if output_id_type == 'Gene':
            geneid = transcript_db[ens_numeric_id].GeneId();
            try:
                try: ens_id = gene_db[geneid].StableId()
                except Exception: ens_id = gene_stable_id_db[geneid].StableId()
            except KeyError: null = [] ### Again, this occurs in version 47
        elif output_id_type == 'Transcription':
            try: ens_id = transcript_db[ens_numeric_id].StableId()
            except Exception: ens_id = transcript_stable_id_db[ens_numeric_id].StableId()
        elif output_id_type == 'Translation':
            try:
                try: protein_id = transcript_to_protein_db[ens_numeric_id]; ens_id = translation_db[protein_id].StableId()
                except Exception: protein_id = transcript_to_protein_db[ens_numeric_id]; ens_id = translation_stable_id_db[protein_id].StableId()
            except KeyError: null = []
        try: values = [ens_id,dbprimary_acc]; values_list.append(values)
        except NameError: null = []
        
    for (ens_numeric_id,dbprimary_acc) in translation_relationship_db:
        if output_id_type == 'Gene':
            transcript_id = translation_db[ens_numeric_id].TranscriptId()
            geneid = transcript_db[transcript_id].GeneId()
            try:
                try: ens_id = gene_db[geneid].StableId()
                except Exception: ens_id = gene_stable_id_db[geneid].StableId()
            except KeyError: null = [] ### Again, this occurs in version 47
        elif output_id_type == 'Transcription':
            transcript_id = translation_db[ens_numeric_id].TranscriptId()
            try: ens_id = transcript_db[transcript_id].StableId()
            except Exception: ens_id = transcript_stable_id_db[transcript_id].StableId()
        elif output_id_type == 'Translation': ens_id = translation_stable_id_db[ens_numeric_id].StableId()
        try: values = [ens_id,dbprimary_acc]; values_list.append(values)
        except NameError: null = []
        
    for (ens_numeric_id,dbprimary_acc) in gene_relationship_db:
        if output_id_type == 'Gene':
            try: ens_id = gene_db[ens_numeric_id].StableId()
            except Exception: ens_id = gene_stable_id_db[ens_numeric_id].StableId()
            values = [ens_id,dbprimary_acc]; values_list.append(values)
        elif output_id_type == 'Transcription':
            transcript_ids = gene_to_transcript_db[ens_numeric_id]
            for transcript_id in transcript_ids:
                try: ens_id = transcript_db[transcript_id].StableId()
                except Exception: ens_id = transcript_stable_id_db[transcript_id].StableId()
                values = [ens_id,dbprimary_acc]; values_list.append(values)
        elif output_id_type == 'Translation':
            transcript_ids = gene_to_transcript_db[ens_numeric_id]
            for transcript_id in transcript_ids:
                try: ### Translate between transcripts to protein IDs
                    protein_id = transcript_to_protein_db[ens_numeric_id]
                    try: ens_id = translation_db[protein_id].StableId()
                    except Exception: ens_id = translation_stable_id_db[protein_id].StableId()
                    values = [ens_id,dbprimary_acc]; values_list.append(values)
                except KeyError: null = []
    values_list = unique.unique(values_list)    
    return values_list

def exportListstoFiles(values_list,headers,output_dir,rewrite):
    global rewrite_existing; rewrite_existing = rewrite
    exportEnsemblTable(values_list,headers,output_dir)
    
def exportEnsemblTable(values_list,headers,output_dir):
    if rewrite_existing == 'no':
        print 'Appending new data to',output_dir
        try:values_list = combineEnsemblTables(values_list,output_dir) ###Combine with previous
        except Exception: null=[]
    
    data = export.ExportFile(output_dir)
    if len(headers)>0:
        headers = string.join(headers,'\t')+'\n'
        data.write(headers)
    for values in values_list:
        try: values = string.join(values,'\t')+'\n'
        except TypeError:
            values_temp = values; values = []
            for value in values_temp: values.append(str(value))
            values = string.join(values,'\t')+'\n'
        data.write(values)
    data.close()
    print 'File:',output_dir,'exported.'

def combineEnsemblTables(values_list,filename):
    fn=filepath(filename); x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            values_list.append(t)
    values_list = unique.unique(values_list)
    return values_list

class EnsemblSQLInfo:
    def __init__(self, filename, url_type, importGroup, key, values, comments):
        self._filename = filename; self._url_type = url_type; self._importGroup = importGroup
        self._key = key; self._values = values; self._comments = comments
    def Filename(self): return self._filename
    def URLType(self): return self._url_type
    def ImportGroup(self): return self._importGroup
    def Key(self): return self._key
    def Values(self):
        value_list = string.split(self._values,'|')
        return value_list
    def FieldNames(self): ### These are fields in the description file to parse from downloaded files
        field_names = self.Values()
        field_names.append(self.Key())
        return field_names
    def setIndexDB(self,index_db): self._index_db = index_db
    def IndexDB(self): return self._index_db
    def Comments(self): return self._comments
    def Report(self): return self.Key()+'|'+self.Values()
    def __repr__(self): return self.Report()

class EnsemblSQLEntryData:
    def setSQLValue(self,header,value):
        if header == "seq_region_start": self.seq_region_start = value
        elif header == "transcript_id": self.transcript_id = value
        elif header == "stable_id": self.stable_id = value
        elif header == "gene_id": self.gene_id = value
        elif header == "biotype": self.biotype = value
        elif header == "translation_id": self.translation_id = value
        elif header == "id": self.id = value
        elif header == "synonym": self.synonym = value
        elif header == "external_db_id": self.external_db_id = value
        elif header == "object_xref_id": self.object_xref_id = value
        elif header == "db_name": self.db_name = value
        elif header == "seq_end": self.seq_end = value
        elif header == "end_exon_id": self.end_exon_id = value
        elif header == "description": self.description = value
        elif header == "hit_id": self.hit_id = value
        elif header == "hit_name": self.hit_id = value
        elif header == "ensembl_object_type": self.ensembl_object_type = value
        elif header == "start_exon_id": self.start_exon_id = value
        elif header == "seq_region_end": self.seq_region_end = value
        elif header == "dbprimary_acc": self.dbprimary_acc = value
        elif header == "seq_start": self.seq_start = value
        elif header == "display_xref_id": self.display_xref_id = value
        elif header == "display_label": self.display_label = value
        elif header == "ensembl_id": self.ensembl_id = value
        elif header == "seq_region_strand": self.seq_region_strand = value
        elif header == "rank": self.rank = value
        elif header == "seq_region_id": self.seq_region_id = value
        elif header == "name": self.name = value
        elif header == "exon_id": self.exon_id = value
        elif header == "is_constitutive": self.is_constitutive = value
        elif header == "interpro_ac": self.interpro_ac = value
        elif header == "xref_id": self.xref_id = value
        elif header == "evalue":
            try: self.evalue = float(value)
            except Exception: self.evalue = 0 ### For yeast, can be NA (likely a problem with Ensembl)
        elif header == "vendor": self.vendor = value
        elif header == "array_id": self.array_id = value
        elif header == "array_chip_id": self.array_chip_id = value
        elif header == "probe_set_id": self.probe_set_id = value
        elif header == "format": self.format = value
        
        else: ###Shouldn't occur, unless we didn't account for an object type
            print 'Warning!!! An object type has been imported which does not exist in this class'
            print 'Object type =',header;sys.exit()
    ### Create objects specified in the SQL Description and Configuration file
    def SeqRegionStart(self): return self.seq_region_start
    def TranscriptId(self): return self.transcript_id
    def StableId(self): return self.stable_id
    def GeneId(self): return self.gene_id
    def Biotype(self): return self.biotype
    def TranslationId(self): return self.translation_id
    def Id(self): return self.id
    def Synonym(self): return self.synonym
    def ExternalDbId(self): return self.external_db_id
    def ObjectXrefId(self): return self.object_xref_id
    def DbName(self): return self.db_name
    def ExonId(self): return self.exon_id
    def IsConstitutive(self): return self.is_constitutive
    def SeqEnd(self): return self.seq_end
    def EndExonId(self): return self.end_exon_id
    def Description(self): return self.description
    def HitId(self): return self.hit_id
    def EnsemblObjectType(self): return self.ensembl_object_type
    def StartExonId(self): return self.start_exon_id
    def SeqRegionEnd(self): return self.seq_region_end
    def DbprimaryAcc(self): return self.dbprimary_acc
    def SeqStart(self): return self.seq_start
    def DisplayXrefId(self): return self.display_xref_id
    def DisplayLabel(self): return self.display_label
    def EnsemblId(self): return self.ensembl_id
    def SeqRegionStrand(self): return self.seq_region_strand
    def Rank(self): return self.rank
    def Name(self): return self.name
    def SeqRegionId(self): return self.seq_region_id
    ### Create new objects designated by downstream custom code
    def setCodingBpInExon(self,coding_bp_in_exon): self.coding_bp_in_exon = coding_bp_in_exon
    def CodingBpInExon(self): return self.coding_bp_in_exon
    def setGenomicStart(self,genomic_start): self.genomic_start = genomic_start
    def GenomicStart(self): return self.genomic_start
    def setGenomicEnd(self,genomic_end): self.genomic_end = genomic_end
    def GenomicEnd(self): return self.genomic_end
    def InterproAc(self): return self.interpro_ac
    def XrefId(self): return self.xref_id
    def Evalue(self): return self.evalue
    def Vendor(self): return self.vendor
    def ArrayID(self): return self.array_id
    def ArrayChipID(self): return self.array_chip_id
    def ProbeSetID(self): return self.probe_set_id
    def Format(self): return self.format
    def setProbeSetID(self,probe_set_id): self.probe_set_id = probe_set_id
        
def importEnsemblSQLInfo(configType):
    filename = 'Config/EnsemblSQL.txt'
    fn=filepath(filename); sql_file_db={}; sql_group_db={}; x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        filename, url_type, importGroup, configGroup, key, values, comments = string.split(data,'\t')
        if x==0: x=1
        else:
            sq = EnsemblSQLInfo(filename, url_type, importGroup, key, values, comments)
            sql_file_db[importGroup,filename] = sq
            ### To conserve memory, only import necessary files for each run (results in 1/5th the memory usage)
            if configType == 'Basic' and configGroup == 'Basic': proceed = 'yes'
            elif configType == 'Basic': proceed = 'no'
            else: proceed = 'yes'
            if proceed == 'yes':
                try: sql_group_db[importGroup].append(filename)
                except KeyError: sql_group_db[importGroup] = [filename]
    return sql_file_db,sql_group_db

def importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_description_dir,sql_group_db,sql_file_db,output_dir,import_group,force):
    if 'Primary' in import_group:
        global exon_db; global transcript_db; global translation_stable_id_db
        global exon_transcript_db; global transcript_stable_id_db; global gene_db
        global exon_stable_id_db; global translation_db; global gene_stable_id_db
        global protein_feature_db; global interpro_db; global external_synonym_db
        global external_db_db; global key_filter_db; global external_filter_db
        global seq_region_db; global array_db; global array_chip_db
        global probe_set_db; global probe_db; global probe_feature_db
        key_filter_db={}; external_filter_db={}; probe_set_db={}

    ###Generic function for importing all EnsemblSQL tables based on the SQLDescription tabl
    for filename in sql_group_db['Description']: ### Gets the SQL description table
        try: sql_filepaths = updateFiles(ensembl_sql_description_dir,output_dir,filename,force)
        except Exception: ### Try again... can be a server problem
            sql_filepaths = updateFiles(ensembl_sql_description_dir,output_dir,filename,force)
            #print ensembl_sql_description_dir,output_dir,filename; print e; sys.exit()
        for sql_filepath in sql_filepaths:
            sql_file_db = importSQLDescriptions(import_group,sql_filepath,sql_file_db)

    for filename in sql_group_db[import_group]:
        sql_filepaths = updateFiles(ensembl_sql_dir,output_dir,filename,force)
        if sql_filepaths == 'stable-combined-version':
            ### Hence, the file was not downloaded
            print filename, 'not present in this version of Ensembl' ### Stable files are merged with the indexed files in Ensembl 65 and later
        else:
            #if 'object' in filename and 'Func' in output_dir: sql_filepaths = ['BuildDBs/EnsemblSQL/Dr/FuncGen/object_xref.001.txt','BuildDBs/EnsemblSQL/Dr/FuncGen/object_xref.002.txt']
            for sql_filepath in sql_filepaths: ### If multiple files for a single table, run all files (retain orginal filename)
                print 'Processing:',filename
                try:
                    try: key_value_db = importPrimaryEnsemblSQLTables(sql_filepath,filename,sql_file_db[import_group,filename])
                    except Exception:
                        if key_value_db == 'stable-combined-version':
                            sql_filepaths == 'stable-combined-version'  ### This is a new issue in version 2.1.1 - temporary fix
                            continue
                    """except IOError:
                        sql_filepaths = updateFiles(ensembl_sql_dir,output_dir,filename,'yes')
                        key_value_db = importPrimaryEnsemblSQLTables(sql_filepath,filename,sql_file_db[import_group,filename])"""
                    if filename == "exon.txt": exon_db = key_value_db
                    elif filename == "exon_transcript.txt": exon_transcript_db = key_value_db
                    elif filename == "exon_stable_id.txt": exon_stable_id_db = key_value_db
                    elif filename == "transcript.txt": transcript_db = key_value_db
                    elif filename == "transcript_stable_id.txt": transcript_stable_id_db = key_value_db
                    elif filename == "translation.txt": translation_db = key_value_db
                    elif filename == "translation_stable_id.txt": translation_stable_id_db = key_value_db
                    elif filename == "gene.txt": gene_db = key_value_db
                    elif filename == "gene_stable_id.txt": gene_stable_id_db = key_value_db
                    elif filename == "protein_feature.txt": protein_feature_db = key_value_db
                    elif filename == "interpro.txt": interpro_db = key_value_db
                    elif filename == "external_synonym.txt": external_synonym_db = key_value_db
                    elif filename == "external_db.txt": external_db_db = key_value_db
                    elif filename == "seq_region.txt": seq_region_db = key_value_db
                    elif filename == "array.txt": array_db = key_value_db
                    elif filename == "array_chip.txt":
                        ### Add chip type and vendor to external_filter_db to filter probe.txt
                        array_chip_db = key_value_db; buildFilterDBForArrayDB(externalDBName)
                    elif filename == "xref.txt":
                        try: xref_db = combineDBs(xref_db,key_value_db,'string')
                        except Exception: xref_db = key_value_db
                        if '.0' in sql_filepath: print 'Entries in xref_db', len(xref_db)
                    elif filename == "object_xref.txt":
                        try: object_xref_db = combineDBs(object_xref_db,key_value_db,'list')
                        except Exception: object_xref_db = key_value_db
                        if '.0' in sql_filepath: print 'Entries in object_xref_db', len(object_xref_db)
                    elif filename == "probe.txt":
                        try: probe_db = combineDBs(probe_db,key_value_db,'list')
                        except Exception: probe_db = key_value_db
                        if '.0' in sql_filepath: print 'Entries in probe_db', len(probe_db)
                        if 'AFFY' in manufacturer and 'ProbeLevel' in analysisType:
                            for probe_id in probe_db:
                                for pd in probe_db[probe_id]:
                                    probe_set_db[pd.ProbeSetID()]=[] ### this dictionary is introduced prior to prob_set.txt import when annotating probe genome location
                    elif filename == "probe_set.txt":
                        if 'AFFY' in manufacturer and 'ProbeLevel' in analysisType:
                            for probe_id in probe_db:
                                for pi in probe_db[probe_id]:
                                    pi.setProbeSetID(key_value_db[pi.ProbeSetID()].Name()) ### Add the probeset annotation name
                            del probe_set_db ### this object is only temporarily created during probe_db addition - used to restrict to probesets in probe_db
                        elif 'AFFY' in manufacturer:
                            probe_set_db = key_value_db
                            del probe_db
                        else:
                            probe_set_db = probe_db
                            del probe_db
                    elif filename == "probe_feature.txt":
                        try: probe_feature_db = key_value_db
                        except Exception: key_value_db = key_value_db
                    else: ###Shouldn't occur, unless we didn't account for a file type
                        print 'Warning!!! A file has been imported which does not exist in the program space'
                        print 'Filename =',filename;sys.exit()
                except IOError,e:
                    print e
                    print '...Likely due to present SQL tables from a prior Ensembl version run that are no longer supported for this version. Ignoring and proceeding..'
    if sql_filepaths != 'stable-combined-version':
        if import_group == 'Primary':
            key_filter_db={}
            for geneid in gene_db:
                gi = gene_db[geneid]; display_xref_id = gi.DisplayXrefId()
                key_filter_db[display_xref_id]=[]
        elif 'Object-Xref' in import_group:
            return object_xref_db
        elif 'Xref' in import_group:
            return xref_db
        elif 'PrimaryFunc' in import_group:
            return xref_db

def combineDBs(db1,db2,type):
    if type == 'string':
        for key in db2: db1[key] = db2[key] ### No common keys ever with string
    if type == 'list':
        for key in db2:
            try: db1[key] = db1[key]+db2[key] ### Occurs when same key exists, but different lists
            except KeyError: db1[key] = db2[key]
    return db1

def updateFiles(ensembl_sql_dir,output_dir,filename,force):
    if force == 'no':
        file_found = verifyFile(output_dir+filename)
        #print file_found,output_dir+filename
        if file_found == 'no':
            index=1; sql_filepaths = []
            while index<10:
                filename_new = string.replace(filename,'.txt','.00'+str(index)+'.txt')
                file_found = verifyFile(output_dir+filename_new); index+=1
                if file_found == 'yes':
                    sql_filepaths.append(output_dir + filename_new)
            if len(sql_filepaths)<1: force = 'yes'
        else: sql_filepaths = [output_dir + filename]
    if force == 'yes': ### Download the files, rather than re-use existing
        ftp_url = ensembl_sql_dir+filename + '.gz'
        try: gz_filepath, status = update.download(ftp_url,output_dir,'')
        except Exception:
            if 'stable' in filename:
                sql_filepaths=[]
        if 'Internet' in status:
            index=1; status = ''; sql_filepaths=[]
            while index<10:
                if 'Internet' not in status: ### sometimes, instead of object_xref.txt the file will be name of object_xref.001.txt
                    ftp_url_new = string.replace(ftp_url,'.txt','.00'+str(index)+'.txt')
                    gz_filepath, status = update.download(ftp_url_new,output_dir,'')
                    if status == 'not-removed':
                        try: os.remove(gz_filepath) ### Not sure why this works now and not before
                        except OSError: status = status
                    sql_filepaths.append(gz_filepath[:-3])
                    index+=1
                    #print [[[sql_filepaths]]]
                else:
                    sql_filepaths=sql_filepaths[:-1]; print ''
                    #print [[sql_filepaths]]
                    break
        else:
            if status == 'not-removed':
                try: os.remove(gz_filepath) ### Not sure why this works now and not before
                except OSError: status = status
            sql_filepaths = [gz_filepath[:-3]]
    #print [sql_filepaths]
    if len(sql_filepaths)==0:
        if 'stable' in filename: sql_filepaths = 'stable-combined-version' ### For Ensembl 65 and later (changed table organization from previous versions)
        else:
            print '\nThe file:',filename, 'is missing from Ensembl FTP directory for this version of Ensembl. Contact altanalyze@gmail.com to inform our developers of this change to the Ensembl database table structure or download the latest version of this software.'; force_exit
    return sql_filepaths

def importPrimaryEnsemblSQLTables(sql_filepath,filename,sfd):
    fn=filepath(sql_filepath)
    index=0; key_value_db={}; data_types={}
    if len(key_filter_db)>0: key_filter = 'yes'
    else: key_filter = 'no'
    
    if len(external_filter_db)>0: external_filter = 'yes'
    else: external_filter = 'no'

    try:
        if len(external_xref_key_db)>0: external_xref_key_filter = 'yes'
        else: external_xref_key_filter = 'no'
    except NameError: external_xref_key_filter = 'no'
    
    #print filename, key_filter, external_filter, external_xref_key_filter
    try: index_db = sfd.IndexDB(); key_name = sfd.Key(); entries=0
    except Exception:
        return 'stable-combined-version'
    if len(key_name)<1: key_value_db=[] ### No key, so store data in a list
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line); data = string.split(data,'\t')
        ese = EnsemblSQLEntryData(); skip = 'no'; entries+=1
        if len(line)>3: ### In the most recent version (plant-16) encountered lines with just a single /
            for index in index_db:
                header_name,value_type = index_db[index]
                try: value = data[index]
                except Exception:
                    if 'array.' in filename: skip = 'yes'; value = ''
                    """
                    ### This will often occur due to prior lines having a premature end of line (bad formatting in source data) creating a bad line (just skip it)
                    else:
                        ### Likely occurs when Ensembl adds new line types to their SQL file (bastards)
                        print index_db,data, index, len(data),header_name,value_type
                        kill"""
                        
                try:
                    if value_type == 'integer': value = int(value) # Integers will take up less space in memory
                except ValueError:
                    if value == '\\N': value = value
                    elif 'array.' in filename: skip = 'yes' ### There can be formatting issues with this file when read in python
                    else: skip = 'yes'; value = '' #print filename,[header_name], index,value_type,[value];kill
                ###Although we are setting each index to value, the distinct headers will instruct
                ###EnsemblSQLEntryData to assign the value to a distinct object
                if header_name != key_name:
                    ese.setSQLValue(header_name,value)
                else:
                    key = value
            ### Filtering primarily used for all Xref, since this database is very large
            #if filename == 'xref.txt':print len(key_name), key_filter,external_filter, external_xref_key_filter;kill
            if skip == 'yes': null = []
            elif 'probe.' in filename:
                ### For all array types (Affy & Other) - used to select specific array types - also stores non-Affy probe-data
                ### Each array has a specific ID - determine if this ID is the one we set in external_filter_db
                if ese.ArrayChipID() in external_filter_db:
                    vendor = external_filter_db[ese.ArrayChipID()][0]
                    if vendor == 'AFFY' and analysisType != 'ProbeLevel':
                        key_value_db[ese.ProbeSetID()] = [] ### key is the probe_id - used for Affy probe ID location analysis
                    else:
                        try: key_value_db[key].append(ese) ### probe_set.txt only appears to contain AFFY IDs, the remainder are saved in probe.txt under probe_id
                        except KeyError: key_value_db[key] = [ese]
            elif 'probe_set.' in filename:
                if analysisType == 'ProbeLevel':
                    if key in probe_set_db: key_value_db[key] = ese
                else:
                    if key in probe_db: key_value_db[key] = [ese]
            elif 'probe_feature.' in filename:
                if key in probe_db: key_value_db[key] = ese
            elif len(key_name)<1: key_value_db.append(ese)
            elif key_filter == 'no' and external_filter == 'no': key_value_db[key] = ese
            elif external_xref_key_filter == 'yes':
                if key in external_xref_key_db:
                    try: key_value_db[key].append(ese)
                    except KeyError: key_value_db[key] = [ese]
            elif 'object_xref' in filename:
              try:
                if key in probe_set_db:
                    if (manufacturer == 'AFFY' and ese.EnsemblObjectType() == 'ProbeSet') or (manufacturer != 'AFFY' and ese.EnsemblObjectType() == 'Probe'):
                        try: key_value_db[key].append(ese)
                        except KeyError: key_value_db[key] = [ese]
                        data_types[ese.EnsemblObjectType()]=[]
              except Exception: print external_xref_key_filter, key_filter, filename;kill
            elif key in key_filter_db: key_value_db[key] = ese
            elif 'seq_region.' in filename: key_value_db[key] = ese ### Specifically applies to ProbeLevel analyses
            elif external_filter == 'yes': ### For example, if UniGene's dbase ID is in the external_filter_db (when parsing xref)
                #if 'xref' in filename: print filename,[header_name], index,value_type,[value],ese.ExternalDbId(),external_filter_db;kill
                try:
                    #print key, [ese.ExternalDbId()], [ese.DbprimaryAcc()], [ese.DisplayLabel()];kill
                    #if key == 1214287: print [ese.ExternalDbId()], [ese.DbprimaryAcc()], [ese.DisplayLabel()]
                    if ese.ExternalDbId() in external_filter_db: key_value_db[key] = ese
                    all_external_ids[ese.ExternalDbId()]=[]
                except AttributeError:
                    print len(external_filter_db),len(key_filter_db)
                    print 'index_db',index_db,'\n'
                    print 'key_name',key_name
                    print len(external_filter_db)
                    print len(key_value_db);kill

    #key_filter_db={}; external_filter_db={}; external_xref_key_db={}
                
    if 'object_xref' in filename:
        print external_xref_key_filter, key_filter, filename
        print "Extracted",len(key_value_db),"out of",entries,"entries for",filename
        print 'Data types linked to probes:',
        for data_type in data_types: print data_type,
        print ''
    return key_value_db
        
def importSQLDescriptions(import_group,filename,sql_file_db):
    fn=filepath(filename)
    index_db={}; index=0
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        if 'CREATE TABLE' in data:
            ###New file descriptions
            file_broken = string.split(data,'`'); filename = file_broken[1]+".txt"
            if (import_group,filename) in sql_file_db:
                sql_data = sql_file_db[import_group,filename]
                target_field_names = sql_data.FieldNames()
            else: target_field_names=[]
        elif data[:3] == '  `':
            field_broken = string.split(data,'`'); field_name = field_broken[1]
            if 'int(' in data: type = 'integer'
            else: type = 'string'
            if field_name in target_field_names: index_db[index]=field_name,type
            index+=1
        elif len(data)<2:
            ###Write out previous data here and clear entries
            if len(index_db)>0: ### Thus fields in the Config file are found in the description file
                sql_data.setIndexDB(index_db)
            index_db = {}; index=0 ### re-set
        elif '/*' in data or 'DROP TABLE IF EXISTS' in data:
            None ### Not sure what this line is that has recently been added
        else: index+=1
    if len(index_db)>0: asql_data.setIndexDB(index_db)
    return sql_file_db

def storeFTPDirs(ftp_server,subdir,dirtype):
    from ftplib import FTP
    ftp = FTP(ftp_server); ftp.login()

    try: ftp.cwd(subdir)
    except Exception:
        subdir = string.replace(subdir,'/mysql','') ### Older version don't have this subdir 
        ftp.cwd(subdir)
    data = []; child_dirs={};species_list=[]
    ftp.dir(data.append); ftp.quit()
    for line in data:
        line = string.split(line,' '); file_dir = line[-1]
        if dirtype in file_dir:
            species_name_data = string.split(file_dir,dirtype)
            species_name = species_name_data[0]
            if 'bacteria' not in species_name: ### Occurs for Bacteria Genomes
                species_name = string.replace(string.upper(species_name[0])+species_name[1:],'_',' ')
            ensembl_sql_dir = 'ftp://'+ftp_server+subdir+'/'+file_dir+'/'
            ensembl_sql_description_dir = file_dir+'.sql'
            child_dirs[species_name] = ensembl_sql_dir,ensembl_sql_description_dir
            species_list.append(species_name)
    species_list.sort()
    return child_dirs,species_list

def getEnsemblVersions(ftp_server,subdir):
    from ftplib import FTP
    ftp = FTP(ftp_server); ftp.login()
    ftp.cwd(subdir)
    data = []; ensembl_versions=[]
    ftp.dir(data.append); ftp.quit()
    for line in data:
        line = string.split(line,' '); file_dir = line[-1]
        if 'release' in file_dir and '/' not in file_dir:
            version_number = int(string.replace(file_dir,'release-',''))
            if version_number>46: ###Before this version, the SQL FTP folder structure differed substantially
                ensembl_versions.append(file_dir)
    return ensembl_versions

def clearall():
    all = [var for var in globals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all: del globals()[var]

def clearvar(varname):
    all = [var for var in globals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all:
        if var == varname: del globals()[var]
    
def getCurrentEnsemblSpecies(version):
    original_version = version
    if 'EnsMart' in version:
        version = string.replace(version,'EnsMart','') ### User may enter EnsMart65, but we want 65 (release- is added before the number)
        version = string.replace(version,'Plus','')
    if 'Plant' in version:
        version = string.replace(version,'Plant','')
    if 'Fungi' in version:
        version = string.replace(version,'Fungi','')
    if 'Bacteria' in version:
        version = string.replace(version,'Bacteria','')
    if 'release' not in version and 'current' not in version:
        version = 'release-'+version
    if 'Plant' in original_version or 'Bacteria' in original_version or 'Fungi' in original_version:
        ftp_server = 'ftp.ensemblgenomes.org'
    else:
        ftp_server = 'ftp.ensembl.org'
    if version == 'current':
        subdir = '/pub/current_mysql'
    elif 'Plant' in original_version:
        subdir = '/pub/plants/'+version+'/mysql'
    elif 'Bacteria' in original_version:
        subdir = '/pub/'+version+'/bacteria/mysql'
    elif 'Fungi' in original_version:
        subdir = '/pub/'+version+'/fungi/mysql'
    else:
        subdir = '/pub/'+version+'/mysql'
    dirtype = '_core_'
    ensembl_versions = getEnsemblVersions(ftp_server,'/pub')
    child_dirs, species_list = storeFTPDirs(ftp_server,subdir,dirtype)
    return child_dirs, species_list, ensembl_versions

def getCurrentEnsemblSequences(version,dirtype,species):
    ftp_server = 'ftp.ensembl.org'
    if version == 'current': subdir = '/pub/current_'+dirtype
    else: subdir = '/pub/'+version+'/'+dirtype
    seq_dir = storeSeqFTPDirs(ftp_server,species,subdir,dirtype)
    return seq_dir

def getCurrentEnsemblGenomesSequences(version,dirtype,species):
    original_version = version    
    if 'Fungi' in version:
        version = string.replace(version,'Fungi','')
    if 'Plant' in version:
        version = string.replace(version,'Plant','')
    if 'Bacteria' in version:
        version = string.replace(version,'Bacteria','')
    ftp_server = 'ftp.ensemblgenomes.org'
    if version == 'current': subdir = '/pub/current_'+dirtype
    elif 'Bacteria' in original_version:
        subdir = '/pub/'+version+'/bacteria/'+dirtype
    elif 'Fungi' in original_version:
        subdir = '/pub/'+version+'/fungi/'+dirtype
    else: subdir = '/pub/plants/'+version+'/'+dirtype
    seq_dir = storeSeqFTPDirs(ftp_server,species,subdir,dirtype)
    return seq_dir

def storeSeqFTPDirs(ftp_server,species,subdir,dirtype):
    from ftplib import FTP
    ftp = FTP(ftp_server); ftp.login()
    #print subdir;sys.exit()
    try: ftp.cwd(subdir)
    except Exception:
        subdir = string.replace(subdir,'/'+dirtype,'') ### Older version don't have this subdir 
        ftp.cwd(subdir)
    data = []; seq_dir=[]; ftp.dir(data.append); ftp.quit()
    for line in data:
        line = string.split(line,' '); file_dir = line[-1]
        if species[1:] in file_dir:
            if '.fa' in file_dir and '.all' in file_dir: seq_dir = 'ftp://'+ftp_server+subdir+'/'+file_dir
            elif '.fa' in file_dir and 'dna.chromosome' in file_dir:
                try: seq_dir.append((file_dir,'ftp://'+ftp_server+subdir+'/'+file_dir))
                except Exception: seq_dir=[]; seq_dir.append((file_dir,'ftp://'+ftp_server+subdir+'/'+file_dir))
            elif '.fa' in file_dir and 'dna.' in file_dir and 'nonchromosomal' not in file_dir:
                ### This is needed when there are numbered chromosomes (e.g., ftp://ftp.ensembl.org/pub/release-60/fasta/anolis_carolinensis/dna/)
                try: seq_dir2.append((file_dir,'ftp://'+ftp_server+subdir+'/'+file_dir))
                except Exception: seq_dir2=[]; seq_dir2.append((file_dir,'ftp://'+ftp_server+subdir+'/'+file_dir))
    if len(seq_dir)==0: seq_dir = seq_dir2
    return seq_dir
    
"""
def getExternalDBs(Species,ensembl_sql_dir,ensembl_sql_description_dir):
    global species; species = Species;
    configType = 'Basic'
    sql_file_db,sql_group_db = importEnsemblSQLInfo(configType) ###Import the Config file with the files and fields to parse from the donwloaded SQL files

    sql_group_db['Description'] = [ensembl_sql_description_dir]
    info = string.split(ensembl_sql_description_dir,'_'); ensembl_build = info[-2]
    sq = EnsemblSQLInfo(ensembl_sql_description_dir, 'EnsemblSQLDescriptions', 'Description', '', '', '')
    sql_file_db['Primary',ensembl_sql_description_dir] = sq
    
    output_dir = 'BuildDBs/EnsemblSQL/'+species+'/'
    importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'Primary',force) ###Download and import the Ensembl SQL files
"""

class SystemData:
    def __init__(self, syscode, sysname, mod):
        self._syscode = syscode; self._sysname = sysname; self._mod = mod
    def SystemCode(self): return self._syscode
    def SystemName(self): return self._sysname
    def MOD(self): return self._mod
    def __repr__(self): return self.SystemCode()+'|'+self.SystemName()+'|'+self.MOD()
    
def importSystemInfo():
    filename = 'Config/source_data.txt'; x=0
    system_list=[]; system_codes={}
    fn=filepath(filename); mod_list=[]
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if '!DOCTYPE' in data:
            fn2 = string.replace(fn,'.txt','_archive.txt')
            import shutil; shutil.copyfile(fn2,fn) ### Bad file was downloaded (with warning)
            importSystemInfo(); break
        else:
            try: sysname=t[0];syscode=t[1]
            except Exception: sysname=''
        try: mod = t[2]
        except Exception: mod = ''
        if x==0: x=1
        elif sysname != '':
            system_list.append(sysname)
            ad = SystemData(syscode,sysname,mod)
            if len(mod)>1: mod_list.append(sysname)
            system_codes[sysname] = ad
    return system_codes,system_list,mod_list


def buildGOEliteDBs(Species,ensembl_sql_dir,ensembl_sql_description_dir,ExternalDBName,configType,analysis_type,Overwrite_previous,Rewrite_existing,external_system_db,force):
    global external_xref_key_db; global species; species = Species; global overwrite_previous; overwrite_previous = Overwrite_previous
    global rewrite_existing; rewrite_existing = Rewrite_existing; global ensembl_build; global externalDBName; externalDBName = ExternalDBName
    global ensembl_build; ensembl_build = string.split(ensembl_sql_dir,'core')[-1][:-1]; global analysisType; analysisType = analysis_type
    global manufacturer; global system_synonym_db; global added_systems; added_systems={}; global all_external_ids; all_external_ids={}

    ### Get System Code info and create DBs for automatically formatting system output names
    ### This is necessary to ensure similiar systems with different names are saved and referenced
    ### by the same system name and system code

    import UI;
    try: system_codes,source_types,mod_types = UI.remoteSystemInfo()
    except Exception: system_codes,source_types,mod_types = importSystemInfo()
    system_synonym_db = {}; system_code_db={}; new_system_codes={}
    for system_name in system_codes: system_code_db[system_codes[system_name].SystemCode()] = system_name
    if externalDBName != 'GO':
        system_code = external_system_db[externalDBName]
        try: system_name = system_code_db[system_code]
        except Exception:
            system_name = externalDBName
            ad = UI.SystemData(system_code,system_name,'')
            new_system_codes[externalDBName] = ad ### Add these to the source_data file if relationships added (see below)
            system_code_db[system_code] = system_name ### Make sure that only one new system name for the code is added
        system_synonym_db[externalDBName] = system_name

    ### Get Ensembl Data            

    sql_file_db,sql_group_db = importEnsemblSQLInfo(configType) ###Import the Config file with the files and fields to parse from the donwloaded SQL files

    sql_group_db['Description'] = [ensembl_sql_description_dir]
    info = string.split(ensembl_sql_description_dir,'_'); ensembl_build = info[-2]
    sq = EnsemblSQLInfo(ensembl_sql_description_dir, 'EnsemblSQLDescriptions', 'Description', '', '', '')
    sql_file_db['Primary',ensembl_sql_description_dir] = sq
    
    output_dir = 'BuildDBs/EnsemblSQL/'+species+'/'
    importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'Primary',force) ###Download and import the Ensembl SQL files
    
    export_status='yes';output_id_type='Gene'
    if analysisType == 'GeneAndExternal':
        if 'over-write previous' in overwrite_previous: output_ens_dir = 'Databases/'+species+'/gene/Ensembl.txt'
        else: output_ens_dir = 'NewDatabases/'+species+'/gene/Ensembl.txt'
        file_found = verifyFile(output_ens_dir)
        revert_force = 'no'
        if file_found == 'no':
            if force == 'yes': force = 'no'; revert_force = 'yes'
            xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'Xref',force)
            buildEnsemblGeneGOEliteTable(species,xref_db,overwrite_previous)        
        ###Export data for Ensembl-External gene system
        buildFilterDBForExternalDB(externalDBName)
        xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'Xref',force)

        external_xref_key_db = xref_db
        if revert_force == 'yes': force = 'yes'
        object_xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'Object-Xref',force)
        num_relationships_exported = buildEnsemblExternalDBRelationshipTable(externalDBName,xref_db,
                                                object_xref_db,output_id_type,species,writeToGOEliteFolder=True)
        if 'EntrezGene' in externalDBName: ### Make a meta DB table for translating WikiPathway primary IDs
            buildEnsemblExternalDBRelationshipTable('meta',xref_db,object_xref_db,
                                                output_id_type,species,writeToGOEliteFolder=True)

        ### Add New Systems to the source_data relationships file (only when valid relationships found)
    
        for externalDBName in new_system_codes:
            if externalDBName in added_systems:
                ad = new_system_codes[externalDBName]
                system_codes[ad.SystemName()] = ad
        UI.exportSystemInfoRemote(system_codes)

    if analysisType == 'FuncGen' or analysisType == 'ProbeLevel':
        sl = string.split(externalDBName,'_'); manufacturer=sl[0];externalDBName=string.replace(externalDBName,manufacturer+'_','')
        ensembl_sql_dir = string.replace(ensembl_sql_dir,'_core_', '_funcgen_')
        ensembl_sql_description_dir = string.replace(ensembl_sql_description_dir,'_core_', '_funcgen_')
        sql_file_db,sql_group_db = importEnsemblSQLInfo('FuncGen') ###Import the Config file with the files and fields to parse from the donwloaded SQL files
        sql_group_db['Description'] = [ensembl_sql_description_dir]
        info = string.split(ensembl_sql_description_dir,'_'); ensembl_build = info[-2]
        sq = EnsemblSQLInfo(ensembl_sql_description_dir, 'EnsemblSQLFuncDescriptions', 'Description', '', '', '')
        sql_file_db['PrimaryFunc',ensembl_sql_description_dir] = sq
        output_dir = 'BuildDBs/EnsemblSQL/'+species+'/FuncGen/'   
        xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'PrimaryFunc',force) ###Download and import the Ensembl SQL files
            
        if analysisType == 'FuncGen':
            #print len(probe_set_db), len(transcript_stable_id_db), len(xref_db)
            object_xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'Object-XrefFunc',force)
            buildEnsemblArrayDBRelationshipTable(xref_db,object_xref_db,output_id_type,externalDBName,species)
        if analysisType == 'ProbeLevel':
            ### There is alot of unnecessary data imported for this specific analysis but this strategy is likely the most straight forward code addition
            sql_file_db['ProbeFeature',ensembl_sql_description_dir] = sq
            importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_dir,sql_group_db,sql_file_db,output_dir,'ProbeFeature',force)
            outputProbeGenomicLocations(externalDBName,species)
    return all_external_ids

def buildEnsemblArrayDBRelationshipTable(xref_db,object_xref_db,output_id_type,externalDBName,species):
    ### Get xref annotations (e.g., external geneID) for a set of Ensembl IDs (e.g. gene IDs)
    external_system = manufacturer
    external_system = external_system[0]+string.lower(external_system[1:])
    
    print 'Exporting',externalDBName
    program_type,database_dir = unique.whatProgramIsThis()
    if program_type == 'AltAnalyze':
        parent_dir = 'AltDatabase/goelite'
    elif 'over-write previous' in overwrite_previous: parent_dir = 'Databases'
    else: parent_dir = 'NewDatabases'
    if 'Affy' in external_system:
        output_dir = parent_dir+'/'+species+'/uid-gene/Ensembl-Affymetrix.txt'
    elif 'Agilent' in external_system:
        output_dir = parent_dir+'/'+species+'/uid-gene/Ensembl-Agilent.txt'
    elif 'Illumina' in external_system:
        output_dir = parent_dir+'/'+species+'/uid-gene/Ensembl-Illumina.txt'
    else: output_dir = parent_dir+'/'+species+'/uid-gene/Ensembl-MiscArray.txt'
        
    version_info = species+' Ensembl relationships downloaded from EnsemblSQL server, build '+ensembl_build
    exportVersionInfo(output_dir,version_info)
    
    headers = ['Ensembl ID',external_system + ' ID']; id_type_db={}; index=0
    id_type_db={}; gene_relationship_db={}; transcript_relationship_db={}; translation_relationship_db={}

    ### Get Ensembl gene-transcript relationships from core files
    ens_transcript_db={}
    for transcript_id in transcript_db:
        ti = transcript_db[transcript_id]
        try: ens_transcript = transcript_db[transcript_id].StableId()
        except Exception: ens_transcript = transcript_stable_id_db[transcript_id].StableId()
        ens_transcript_db[ens_transcript] = ti

    values_list=[]; nulls=0; type_errors=[]
    ### Get array ID-gene relationships from exref and probe dbs
    if len(probe_set_db)>0:
        for probe_set_id in object_xref_db:
            try:
                for ox in object_xref_db[probe_set_id]:
                    try:
                        transcript_id = ox.XrefId()
                        ens_transcript = xref_db[transcript_id].DbprimaryAcc()
                        ti = ens_transcript_db[ens_transcript]; geneid = ti.GeneId()
                        try: ens_gene = gene_db[geneid].StableId()
                        except Exception: ens_gene = gene_stable_id_db[geneid].StableId()
                        for ps in probe_set_db[probe_set_id]:
                            probeset = ps.Name() ### The same probe ID shouldn't exist more than once, but for some arrays it does (associaed with multiple probe names)
                            values_list.append([ens_gene,probeset])
                    except Exception: nulls+=1
            except TypeError,e:
                null=[]
                print probe_set_id, len(probe_set_db), len(probe_set_db[probe_set_id])
                for ps in probe_set_db[probe_set_id]:
                    probeset = ps.Name()
                    print probeset
                print 'ERROR ENCOUNTERED.',e; sys.exit()

        values_list = unique.unique(values_list)
        print 'ID types linked to '+external_system+' are: Gene:',len(values_list),'... not linked:',nulls
        if len(values_list)>0:
            exportEnsemblTable(values_list,headers,output_dir)
            
    else: print "****** Unknown problem encountered with the parsing of:", externalDBName          

def outputProbeGenomicLocations(externalDBName,species):
    external_system = manufacturer
    external_system = external_system[0]+string.lower(external_system[1:])
    
    print 'Exporting',externalDBName
    output_dir = 'Affymetrix/'+species+'/'+externalDBName+'.txt'
    
    if 'Affy' in external_system:
        headers = ['ProbeXY','ProbesetName','Chr','Start','End','Strand']; values_list=[]
        print len(probe_db), ':Probes in Ensembl', len(probe_feature_db),':Probes with genomic locations'
        for probe_id in probe_db:
            for pi in probe_db[probe_id]: ### probe ID
                probe_set_name = pi.ProbeSetID() ### probe_set_id is only a superfluous annotation - not used downstream
                probe_name = pi.Name()
                try:
                    pl = probe_feature_db[probe_id] ### probe genomic location data
                    chr = seq_region_db[pl.SeqRegionId()].Name()
                    seq_region_start = pl.SeqRegionStart(); seq_region_end = pl.SeqRegionEnd(); strand = pl.SeqRegionStrand()
                    values_list.append([probe_name,probe_set_name,chr,seq_region_start,seq_region_end,strand])
                except Exception: null=[]
        values_list = unique.unique(values_list)
        print 'ID types linked to '+external_system+' are: Probe:',len(values_list)
        if len(values_list)>0:
            
            exportEnsemblTable(values_list,headers,output_dir)
    
def buildEnsemblGeneGOEliteTable(species,xref_db,overwrite_previous):
    ### Get Definitions and Symbol for each Ensembl gene ID
    values_list = []
    program_type,database_dir = unique.whatProgramIsThis()
    if program_type == 'AltAnalyze':
        output_dir = 'AltDatabase/goelite/'+species+'/gene/Ensembl.txt'
    elif 'over-write previous' in overwrite_previous:
        output_dir = 'Databases/'+species+'/gene/Ensembl.txt'
    else:
        output_dir = 'NewDatabases/'+species+'/gene/Ensembl.txt'
    headers = ['ID','Symbol','Description']
    for geneid in gene_db:
        gi = gene_db[geneid]; display_xref_id = gi.DisplayXrefId(); description = gi.Description()
        if description == '\\N': description = ''
        try: symbol = xref_db[display_xref_id].DisplayLabel()
        except KeyError: symbol = ''
        try:
            try: ens_gene = gene_db[geneid].StableId()
            except Exception: ens_gene = gene_stable_id_db[geneid].StableId()
            values = [ens_gene,symbol,description]
            values_list.append(values)
        except KeyError: null=[] ### In version 47, discovered that some geneids are not in the gene_stable - inclear why this is
        ###Also, create a filter_db that allows for the creation of Ensembl-external gene db tables
    exportEnsemblTable(values_list,headers,output_dir)

def verifyFile(filename):
    fn=filepath(filename); file_found = 'yes'
    try:
        for line in open(fn,'rU').xreadlines():break
    except Exception: file_found = 'no'
    return file_found
            
if __name__ == '__main__':
    import update
    dp = update.download_protocol('ftp://ftp.ensembl.org/pub/release-72/mysql/macaca_mulatta_core_72_10/gene.txt.gz','AltDatabase/ensembl/Ma/EnsemblSQL/','')
    reload(update)
    dp = update.download_protocol('ftp://ftp.ensembl.org/pub/release-72/mysql/macaca_mulatta_core_72_10/gene.txt.gz','AltDatabase/ensembl/Ma/EnsemblSQL/','');sys.exit()

    getGeneTranscriptOnly('Gg','Basic','EnsMart65','yes');sys.exit()
    #getChrGeneOnly('Hs','Basic','EnsMart65','yes');sys.exit()
    analysisType = 'GeneAndExternal'; externalDBName_list = ['Ens_Gg_transcript']
    force = 'yes'; configType = 'Basic'; overwrite_previous = 'no'; iteration=0; version = 'current'
    print 'proceeding'
    analysisType = 'ExternalOnly'
    
    ensembl_version = '65'
    species = 'Gg'
    
    #ensembl_version = 'Fungi27'
    #species = 'Nc'
    #print string.replace(unique.getCurrentGeneDatabaseVersion(),'EnsMart','');sys.exit()
    #getEnsemblTranscriptSequences(ensembl_version,species,restrictTo='cDNA');sys.exit()
    
    #getFullGeneSequences('Bacteria18','bacteria_1_collection'); sys.exit()
    #for i in child_dirs: print child_dirs[i]
    #"""
    ### WON'T WORK FOR MORE THAN ONE EXTERNAL DATABASE -- WHEN RUN WITHIN THIS MOD
    species_full = 'Neurospora crassa'
    species_full = 'Gallus gallus'
    species_full = 'Mus musculus'; ensembl_version = '72'; force = 'no'; species = 'Mm'; analysisType = 'AltAnalyzeDBs'; configType = 'Advanced'
    #child_dirs, ensembl_species, ensembl_versions = getCurrentEnsemblSpecies(ensembl_version)
    #genus,species = string.split(species_full,' '); species = genus[0]+species[0]
    #ensembl_sql_dir,ensembl_sql_description_dir = child_dirs[species_full]
    rewrite_existing = 'no'
    external_system_db = {'Ens_Gg_transcript':'Et'}
    for externalDBName in externalDBName_list:
        if force == 'yes' and iteration == 1: force = 'no'
        #buildGOEliteDBs(species,ensembl_sql_dir,ensembl_sql_description_dir,externalDBName,configType,analysisType,overwrite_previous,rewrite_existing,external_system_db,force); iteration+=1
        buildEnsemblRelationalTablesFromSQL(species,configType,analysisType,externalDBName,ensembl_version,force)

