import matplotlib
import matplotlib.pyplot as pylab
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy

import string
import time
import random
import math
import sys, os
import sqlite3
import export

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filepath(filename):
    try:
        import unique ### local to AltAnalyze
        fn = unique.filepath(filename)
    except Exception:
        ### Should work fine when run as a script with this (AltAnalyze code is specific for packaging with AltAnalyze)
        dir=os.path.dirname(dirfile.__file__)
        try: dir_list = os.listdir(filename); fn = filename ### test to see if the path can be found (then it is the full path)
        except Exception: fn=os.path.join(dir,filename)
    return fn

##### SQLite Database Access ######

def createSchemaTextFile():
    schema_filename = 'relational_databases/AltAnalyze_schema.sql'
    export_data = export.ExportFile(schema_filename)
    schema_text ='''-- Schema for species specific AltAnalyze transcript data.

-- Genes store general information on each Ensembl gene ID
create table genes (
    id          text primary key,
    name        text,
    description text,
    chr         text,
    strand      text
);

-- Exon_regions store coordinates and annotations for AltAnalyze defined unique exon regions
create table exons (
    id           text,
    start        integer,
    end          integer,
    gene         text not null references genes(id)
);

-- Junctions store coordinates (5' donor and 3' acceptor splice site) and annotations
create table junctions (
    id           text,
    start        integer,
    end          integer,
    gene         text not null references genes(id)
);

-- Transcripts store basic transcript data
create table transcripts (
    id           text primary key,
    gene         text not null references genes(id)
);

-- CDS store basic transcript data
create table cds (
    id           text primary key references transcripts(id),
    start        integer,
    end          integer
);

-- Transcript_exons store coordinates for Ensembl and UCSC defined transcript exons
create table transcript_exons (
    start        integer,
    end          integer,
    transcript   text not null references transcripts(id)
);

-- Proteins store basic protein info
create table proteins (
    id           text primary key,
    genome_start integer,
    genome_end   integer,
    transcript   text not null references transcripts(id)
);

-- Protein_features store InterPro and UniProt protein feature info
create table protein_feature (
    id           text,
    name         text,
    genome_start integer,
    genome_end   integer,
    protein      text not null references proteins(id)
);
'''
    ### We will need to augment the database with protein feature annotations for 
    export_data.write(schema_text)
    export_data.close()
   
def populateSQLite(species,platform):
    global conn
    """ Since we wish to work with only one gene at a time which can be associated with a lot of data
    it would be more memory efficient transfer this data to a propper relational database for each query """
    
    db_filename = filepath('AltDatabase/'+species+'/'+platform+'/AltAnalyze.db') ### store in user directory
    schema_filename = filepath('AltDatabase/'+species+'/'+platform+'/AltAnalyze_schema.sql')
    
    ### Check to see if the database exists already and if not creat it
    db_is_new = not os.path.exists(db_filename)

    with sqlite3.connect(db_filename) as conn:
        if db_is_new:
            createSchemaTextFile()
            print 'Creating schema'
            with open(schema_filename, 'rt') as f:
                schema = f.read()
            conn.executescript(schema)
    
            print 'Inserting initial data'
            species = 'Hs'
            importAllTranscriptData(species)
        else:
            print 'Database exists, assume schema does too.'
            retreiveDatabaseFields()
            sys.exit()
    #return conn
    
def retreiveDatabaseFields():
    """ Retreive data from specific fields from the database """
    cursor = conn.cursor()
    
    id = 'ENSG00000114127'
    query = "select id, name, description, chr, strand from genes where id = ?"
    cursor.execute(query,(id,)) ### In this way, don't have to use %s and specify type
    
    #if only one entry applicable
    #id, name, description, chr, strand = cursor.fetchone()
    #print '%s %s %s %s %s' % (id, name, description, chr, strand);sys.exit()
    
    for row in cursor.fetchall():
        id, name, description, chr, strand = row
        print '%s %s %s %s %s' % (id, name, description, chr, strand)
  
def bulkLoading():
    import csv
    import sqlite3
    import sys
    
    db_filename = 'todo.db'
    data_filename = sys.argv[1]
    
    SQL = """insert into task (details, priority, status, deadline, project)
             values (:details, :priority, 'active', :deadline, :project)
          """
    
    with open(data_filename, 'rt') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        
        with sqlite3.connect(db_filename) as conn:
            cursor = conn.cursor()
            cursor.executemany(SQL, csv_reader)
        
def verifyFile(filename):
    fn=filepath(filename)
    try:
        for line in open(fn,'rU').xreadlines(): found = True; break
    except Exception: found = False
    return  found

def isoformViewer():
    """
    Make a "broken" horizontal bar plot, ie one with gaps
    """    
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.broken_barh([ (110, 30), (150, 10) ] , (10, 5), facecolors=('gray','blue')) # (position, length) - top row
    ax.broken_barh([ (10, 50), (100, 20),  (130, 10)] , (20, 5),
                    facecolors=('red', 'yellow', 'green')) # (position, length) - next row down

    ### Straight line
    pylab.plot((140,150),(12.5,12.5),lw=2,color = 'red') ### x coordinates of the line, y-coordinates of the line, line-thickness - iterate through a list of coordinates to do this
    
    ### Curved line
    verts = [
        (140, 15),  # P0
        (145, 20), # (x coordinate, half distance to this y coordinate)
        (150, 15), # P2
        ]
    
    codes = [Path.MOVETO,
             Path.CURVE4,
             Path.CURVE4,
             ]
    
    path = Path(verts, codes)
    
    patch = patches.PathPatch(path, facecolor='none', lw=2, edgecolor = 'green')
    ax.add_patch(patch)
    
    #midpt = cubic_bezier(pts, .5)
    ax.text(142, 17.7, '25 reads')
    
    ### End-curved line
    
    ax.set_ylim(5,35)
    ax.set_xlim(0,200)
    ax.set_xlabel('Transcript Exons')
    ax.set_yticks([15,25])
    ax.set_yticklabels(['isoform A', 'isoform B'])
    ax.grid(True)
    """
    ax.annotate('alternative splice site', (61, 25),
                xytext=(0.8, 0.9), textcoords='axes fraction',
                arrowprops=dict(facecolor='black', shrink=0.05),
                fontsize=16,
                horizontalalignment='right', verticalalignment='top')
    """
    pylab.show()
  
def countReadsGenomicInterval():
    ### Could use an interval tree implementation - http://informatics.malariagen.net/2011/07/07/using-interval-trees-to-query-genome-annotations-by-position/
    """ other options -
    http://www.biostars.org/post/show/99/fast-interval-intersection-methodologies/
    PMID: 17061921
    http://pypi.python.org/pypi/fastinterval
    """
    a = numpy.array([0,1,2,3,4,5,6,7,8,9]) ### could be genomic start of reads from a chromosome
    count = ((25 < a) & (a < 100)).sum() ### gives the total count

class FeatureData:
    def __init__(self,start,end,annotation):
        self.start = start; self.end = end; self.annotation = annotation
    def Start(self): return self.start
    def End(self): return self.end
    def Annotation(self): return self.annotation

class GeneData:
    def __init__(self,chr,strand):
        self.chr = chr; self.strand = strand
    def setAnnotations(self,symbol,description):
        self.symbol = symbol; self.description = description
    def Chr(self): return self.chr
    def Strand(self): return self.strand
    def Symbol(self): return self.symbol
    def Description(self): return self.description
    
class TranscriptData:
    def __init__(self,strand,start,stop):
        self.start = start; self.stop = stop; self.strand = strand
    def Strand(self): return self.strand
    def Start(self): return self.start
    def Stop(self): return self.stop
            
def importExonAnnotations(species,type):
    start_time = time.time()
    if 'exons' in type:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
    else:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_junction.txt'

    fn=filepath(filename); x=0; exon_annotation_db={}; gene_annotation_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            gene, exonid, chr, strand, start, end, constitutive_call, ens_exon_ids, splice_events, splice_junctions = t; proceed = 'yes'
            if proceed == 'yes': #if gene == 'ENSG00000111671':
                if type == 'junctions': 
                    exon1_start,exon1_end = string.split(start,'|')
                    exon2_start,exon2_end = string.split(end,'|')
                    if strand == '-':
                        exon1_end,exon1_start = exon1_start,exon1_end
                        exon2_end,exon2_start = exon2_start,exon2_end
                        start = int(exon1_end); end = int(exon2_start)
                    else:
                        start = int(exon1_end); end = int(exon2_start)
                else:
                    start = int(start); end = int(end)
                    if gene not in gene_annotation_db:
                        gd = GeneData(chr,strand)
                        gene_annotation_db[gene]=gd
                        
                ### Store this data in the SQL database
                command = """insert into %s (id, start, end, gene)
                values ('%s', %d, %d, '%s')""" % (type,exonid,start,end,gene)
                conn.execute(command)
                        
                #gi = FeatureData(start,end,exonid)
                #try: exon_annotation_db[gene].append(gi)
                #except KeyError: exon_annotation_db[gene]=[gi]

    time_diff = str(round(time.time()-start_time,1))
    print 'Dataset import in %s seconds' % time_diff
    
    if type == 'exons':
        return gene_annotation_db

def importGeneAnnotations(species,gene_annotation_db):
    start_time = time.time()
    gene_annotation_file = "AltDatabase/ensembl/"+species+"/"+species+"_Ensembl-annotations_simple.txt"
    fn=filepath(gene_annotation_file)
    count = 0            
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if count == 0: count = 1
        else:
            gene, description, symbol = string.split(data,'\t')
            description = string.replace(description,"'","") ### single ' will cause problems
            #gene_annotation_db[gene].setAnnotations(symbol, description) ### don't need to store this
            chr = gene_annotation_db[gene].Chr()
            strand = gene_annotation_db[gene].Strand()
            
            ### Store this data in the SQL database
            command = """insert into genes (id, name, description, chr, strand)
            values ('%s', '%s', '%s', '%s', '%s')""" % (gene,symbol,description,chr,strand)
            try: conn.execute(command)
            except Exception:
                print [command];sys.exit()
                
    del gene_annotation_db
    time_diff = str(round(time.time()-start_time,1))
    #print 'Dataset import in %s seconds' % time_diff

def importProcessedSpliceData(filename):
    start_time = time.time()
    fn=filepath(filename)
    count = 0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if count == 0: count = 1
        else:
            gene, symbol, description = string.split(data,'\t')
            gene_annotation_db[gene].setAnnotations(symbol, description)
            
    time_diff = str(round(time.time()-start_time,1))
    print 'Dataset import in %s seconds' % time_diff
    
def importEnsExonStructureData(filename,option):
    start_time = time.time()
    fn=filepath(filename); count=0; last_transcript = ''
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if count==0: count=1
        else:
            gene, chr, strand, start, end, ens_exonid, constitutive_exon, transcript = t
            end = int(end); start = int(start)

            if option == 'SQL':
                ### Store this data in the SQL database
                command = """insert into transcripts (id, gene)
                values ('%s', '%s')""" % (transcript, gene)
                try: conn.execute(command)
                except Exception: None ### Occurs when transcript in database already
                
                command = """insert into transcript_exons (start, end, transcript)
                values ('%d', '%d', '%s')""" % (start, end, transcript)
                conn.execute(command)
            elif option == 'transcript':
                ### Need to separate the two types for in silico translation
                if 'Ensembl' in filename:
                    type = 'Ensembl'
                else:
                    type = 'GenBank'
                try: gene_transcript_db[gene].append((transcript,type))
                except Exception: gene_transcript_db[gene] = [(transcript,type)]
            elif option == 'exon':
                if transcript in cds_location_db and transcript not in cds_genomic_db: # and strand == '-1'
                    cds_start, cds_stop = cds_location_db[transcript]
                    cds_start, cds_stop = int(cds_start), int(cds_stop)
                    if transcript != last_transcript:
                        cumulative = 0
                        last_transcript = transcript
                        
                    if strand == '-1': start_correction = -3; end_correction = 2
                    else: start_correction = 0; end_correction = 0
                    diff1 = cumulative-cds_start
                    diff2 = ((end-start+cumulative) - cds_start)+2
                    diff3 = cumulative-cds_stop 
                    diff4 = ((end-start+cumulative) - cds_stop)+2

                    if diff1 <= 0 and diff2 > 0: ### CDS start is in the first exon
                        exon_length = abs(end-start)+1
                        coding_bp_in_exon = exon_length - cds_start
                        if strand == '-1':
                            cds_genomic_start = start + coding_bp_in_exon + start_correction
                        else:
                            cds_genomic_start = end + coding_bp_in_exon
                    if diff3 < 0 and diff4 >= 0: ### CDS start is in the first exon
                        coding_bp_in_exon = cds_stop
                        if strand == '-1':
                            cds_genomic_stop = end - coding_bp_in_exon + start_correction
                        else:
                            cds_genomic_stop = start + coding_bp_in_exon
                            
                        try:
                            cds_genomic_db[transcript] = cds_genomic_start,cds_genomic_stop
                        except Exception:
                            print cds_start, cds_stop, transcript, cds_genomic_stop; sys.exit()
                        if transcript == 'ENST00000326513':
                            print chr+':'+str(cds_genomic_stop)+'-'+str(cds_genomic_start)
                        del cds_genomic_stop
                        del cds_genomic_start
                    if transcript == 'ENST00000326513':
                        print diff1,diff2,diff3,diff4
                        print end,start,cumulative,cds_start,cds_stop
                    cumulative += (end-start)
                    """
                    cumulative = 1
                        last_transcript = transcript
                        
                    diff1 = cumulative-cds_start
                    diff2 = ((end-start+cumulative) - cds_start)+1
                    diff3 = cumulative-cds_stop
                    diff4 = ((end-start+cumulative) - cds_stop)+1

                    if diff1 <= 0 and diff2 > 0: ### CDS start is in the first exon
                        if strand == '1':
                            cds_genomic_start = end - diff2 + 1
                        else:
                            cds_genomic_start = start + diff2
                    if diff3 < 0 and diff4 >= 0: ### CDS start is in the first exon
                        if strand == '1':
                            cds_genomic_stop = end - diff4
                        else:
                            cds_genomic_stop = start + diff4 - 1
                        try:
                            cds_genomic_db[transcript] = cds_genomic_start,cds_genomic_stop
                        except Exception:
                            print cds_start, cds_stop, transcript, cds_genomic_stop; sys.exit()
                        if transcript == 'ENST00000436739':
                            print chr+':'+str(cds_genomic_stop)+'-'+str(cds_genomic_start);sys.exit()
                        del cds_genomic_stop
                        del cds_genomic_start
                    if transcript == 'ENST00000436739':
                        print diff1,diff2,diff3,diff4
                        print end,start,cumulative,cds_stop
                    cumulative += (end-start)
                    """
                            
    #sys.exit()
    time_diff = str(round(time.time()-start_time,1))
    print 'Dataset import in %s seconds' % time_diff

def importProteinFeatureData(species,cds_location_db):
    filename = 'AltDatabase/ensembl/'+species+'/ProteinFeatureIsoform_complete.txt'
    db={}
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        transcript, protein, residue_start, residue_stop, feature_annotation, feature_id = t
        if transcript in cds_location_db:
            cds_start,cds_end = cds_location_db[transcript]
            cds_feature_start = cds_start+(int(residue_start)*3)
            cds_feature_end = cds_start+(int(residue_stop)*3)
            #print residue_start, residue_stop
            #print cds_start,cds_end
            #print transcript, protein, feature_annotation, feature_id, cds_feature_start, cds_feature_end;sys.exit()
            if 'blah' in transcript:
                command = """insert into transcripts (id, gene)
                values ('%s', '%s')""" % (transcript, gene)
                try: conn.execute(command)
                except Exception: None ### Occurs when transcript in database already
            ### There are often many features that overlap within a transcript, so consistently pick just one
            if transcript in transcript_feature_db:
                db = transcript_feature_db[transcript]
                db[cds_feature_start, cds_feature_end].append([feature_annotation, feature_id])
            else:
                db={}
                db[cds_feature_start, cds_feature_end]=[[feature_annotation, feature_id]]
                transcript_feature_db[transcript] = db
                
    #for transcript in transcript_feature_db:
        

    print len(db)
    
def importAllTranscriptData(species):
    """
    gene_annotation_db = importExonAnnotations(species,'exons') ### stores gene annotations and adds exon data to SQL database
    importExonAnnotations(species,'junctions') ### adds junction data to SQL database
    importGeneAnnotations(species,gene_annotation_db) ### adds gene annotations to SQL database
    
    option = 'SQL'
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    importEnsExonStructureData(filename,option)
    filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_mrna.txt'
    try: importEnsExonStructureData(filename,option)
    except Exception: None ### Not available for all species - needs to be built prior to transcript model creation
    """
    cds_location_db = importCDSsimple(species) ### Import all CDS coordinates
    importProteinFeatureData(species,cds_location_db)

    processed_splice_file = '/Users/nsalomonis/Desktop/AltExonViewer/ProcessedSpliceData/Hs_RNASeq_H9_ES_vs_BJ1_Fibroblast.ExpCutoff-2.0_average-splicing-index-ProcessedSpliceData.txt'

def importCDSsimple(species):
    """ Reads in the combined mRNA CDS coordinates compiled in the IdentifyAltIsoforms.importCDScoordinates() """
    
    cds_location_db={}
    filename = 'AltDatabase/ensembl/'+species +'/AllTranscriptCDSPositions.txt'
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():             
        line_data = cleanUpLine(line)  #remove endline
        mRNA_AC,start,end = string.split(line_data,'\t')
        start,end = int(start),int(end)
        
        #command = """insert into cds (id, start, end)
        #values ('%s', '%s', '%s')""" % (transcript, start, end)
        #try: conn.execute(command)
        #except Exception: None ### Occurs when transcript in database already

        cds_location_db[mRNA_AC] = start, end
    return cds_location_db
                                                       
def alignAllDomainsToTranscripts(species,platform):
    """ This function is only run during the database build process to create files available for subsequent download.
    This recapitulates several functions executed during the database build process but does so explicitely for each
    isoform with the goal of obtained genomic coordinates of each protein feature post de novo sequence alignment.
    This includes all Ensembl proteins, UCSC mRNAs and in silico translated RNAs """
    
    ### Import all transcript to gene associations for Ensembl and UCSC transcripts
    global gene_transcript_db
    gene_transcript_db={}
    option = 'transcript'
    print 'Importing transcript data into memory'
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    importEnsExonStructureData(filename,option)
    filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_mrna.txt'
    try: importEnsExonStructureData(filename,option)
    except Exception: None ### Not available for all species - needs to be built prior to transcript model creation
    
    from build_scripts import FeatureAlignment
    ucsc_transcripts={}
    gene_db = {}
    gene_transcript_db = FeatureAlignment.eliminateRedundant(gene_transcript_db)
    for gene in gene_transcript_db:
        for (ac,type) in gene_transcript_db[gene]:
            if type != 'Ensembl':
                ucsc_transcripts[ac]=[] ### Store all the untranslated UCSC mRNAs
        gene_db[gene] = [gene] ### mimics the necessary structure for FeatureAlignment
    ### Identify untranslated Ensembl transcripts
    
    print 'Importing Ensembl transcript to protein'
    ens_transcript_protein_db = importEnsemblTranscriptAssociations(species)
    
    ### Import protein ID and protein sequence into a dictionary
    #global protein_sequence_db
    #protein_sequence_db = FeatureAlignment.remoteEnsemblProtSeqImport(species) ### All Ensembl protein sequences
    
    """This code imports all protein sequences (NCBI, Ensembl, in silico translated) associated with optimal isoform pairs,
    however, not all isoforms analyzed in the database are here, hence, this should be considered a subset of in silico
    translated Ensembl mRNAs, UCSC ,RNAs, and known analyzed UCSC proteins"""
    #ucsc_transcripts={}
    #ucsc_transcripts['BC065499']=[]
    #ucsc_transcripts['AK309510']=[] ### in silico translated
    #ens_transcript_protein_db={}
    ### Download or translate ANY AND ALL mRNAs considered by AltAnalyze via in silico translation
    from build_scripts import IdentifyAltIsoforms
    analysis_type = 'fetch_new' # analysis_type = 'fetch' ???

    #IdentifyAltIsoforms.remoteTranslateRNAs(species,ucsc_transcripts,ens_transcript_protein_db,analysis_type)
    ### Derive all protein ID, domain and genomic coordinate data from Ensembl and UniProt
    """ This data is available for Ensembl and UniProt isoforms but we re-derive the associations based on sequence for completeness """

    ### Get the domain sequences and genomic coordinates
    """
    # for testing
    gt = {}; y=0
    for gene in gene_db:
        if y < 20:
            gt[gene] = gene_db[gene]
        else: break
        y+=1
    """
    protein_ft_db,domain_gene_counts = FeatureAlignment.grab_exon_level_feature_calls(species,platform,gene_db)
    from build_scripts import ExonAnalyze_module
    seq_files, mRNA_protein_seq_db = IdentifyAltIsoforms.importProteinSequences(species,'getSequence') ### Import all available protein sequences (downloaded or in silico)
    coordinate_type = 'genomic'; #coordinate_type = 'protein'
    ExonAnalyze_module.getFeatureIsoformGenomePositions(species,protein_ft_db,mRNA_protein_seq_db,gene_transcript_db,coordinate_type)

    ### We may need to augment the above domain coordinate to isoform information with the Ensembl and UniProt files (if seq alignment failed for some reason - see grab_exon_level_feature_calls)!
    
def importEnsemblTranscriptAssociations(species):
    """ Import all protein ID-gene relationships used in AltAnalyze. Requires accessing multiple flat files """
    ens_transcript_protein_db={}
    
    ### Import Ensembl protein IDs
    import gene_associations
    gene_import_dir = '/AltDatabase/ensembl/'+species
    g = gene_associations.GrabFiles()
    g.setdirectory(gene_import_dir)
    filedir,filename = g.searchdirectory('Ensembl_Protein__')

    fn=filepath(filedir)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        gene,transcript,protein = string.split(data,'\t')
        if len(protein)==0:
            ens_transcript_protein_db[transcript] = transcript ### Infer protein sequence (in silico translation)
        else:
            ens_transcript_protein_db[transcript] = protein
    return ens_transcript_protein_db

def getCodingGenomicCoordinates(species):
    global cds_location_db
    global cds_genomic_db
    from build_scripts import IdentifyAltIsoforms
    cds_location_db = IdentifyAltIsoforms.importCDScoordinates(species)
    #print cds_location_db['ENST00000436739'];sys.exit()
    cds_genomic_db={}
    
    option = 'exon'
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    importEnsExonStructureData(filename,option)
    filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_mrna.txt'
    try: importEnsExonStructureData(filename,option)
    except Exception: None ### Not available for all species - needs to be built prior to transcript model creation
    
def buildAltExonDatabases(species,platform):
    alignAllDomainsToTranscripts(species,platform)
    
if __name__ == '__main__':
    #isoformViewer();sys.exit()
    species = 'Mm'; type = 'exon'; platform = 'RNASeq'
    #importAllTranscriptData(species); sys.exit()
    #getCodingGenomicCoordinates(species)
    #importProteinFeatureData(species)
    buildAltExonDatabases(species,platform)
    sys.exit()
    populateSQLite()
    #sys.exit()
    #importAllTranscriptData(species);sys.exit()
    #test()
    isoformViewer()
    