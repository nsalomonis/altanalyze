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

try: clearall()
except NameError: null = [] ### Occurs when re-running the script to clear all global variables

import sys, string
import os.path
import unique
import export
import update

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

def buildEnsemblRelationalTablesFromSQL(species,configType,analysisType,externalDBName,force):
    import UI; import update; global external_xref_key_db
    file_location_defaults = UI.importDefaultFileLocations()
    """Identify the appropriate download location for the Ensembl SQL databases for the selected species"""
    fld = file_location_defaults['EnsemblSQL']
    fldd = file_location_defaults['EnsemblSQLDescriptions']
    
    for fl in fld:
        if species in fl.Species(): ensembl_sql_dir = fl.Location()
    fl = fldd[0]; ensembl_sql_description_dir = fl.Location() ### Since the name for the description file varies for species, automatically use human

    global ensembl_build; ensembl_build = string.split(ensembl_sql_dir,'core')[-1][:-1]
        
    sql_file_db,sql_group_db = importEnsemblSQLInfo(configType) ###Import the Config file with the files and fields to parse from the donwloaded SQL files
    output_dir = 'AltDatabase/ensembl/'+species+'/EnsemblSQL/'
    importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_description_dir,sql_group_db,sql_file_db,output_dir,'Primary',force) ###Download and import the Ensembl SQL files

    if analysisType != 'ExternalOnly':
        ### Build and export the basic Ensembl gene annotation table
        xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_description_dir,sql_group_db,sql_file_db,output_dir,'Xref',force)
        buildEnsemblGeneAnnotationTable(species,xref_db)

    if analysisType == 'ExternalOnly':
        ###Export data for Ensembl-External gene system
        buildFilterDBForExternalDB(externalDBName)
        xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_description_dir,sql_group_db,sql_file_db,output_dir,'Xref',force)
        external_xref_key_db = xref_db; #resetExternalFilterDB()
        object_xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_description_dir,sql_group_db,sql_file_db,output_dir,'Object-Xref',force)
        export_status='yes';output_id_type='Gene'
        buildEnsemblExternalDBRelationshipTable(externalDBName,xref_db,object_xref_db,output_id_type,species)

    if analysisType == 'AltAnalyzeDBs':
        ###Export data for Ensembl-External gene system
        buildTranscriptStructureTable()
        
        ###Get Interpro AC display name
        buildFilterDBForExternalDB('Interpro')
        xref_db = importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_description_dir,sql_group_db,sql_file_db,output_dir,'Xref',force)
        getDomainGenomicCoordinates(species,xref_db)
                
def buildTranscriptStructureTable():
    ### Function for building Ensembl transcription annotation file
    temp_values_list = []; output_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    headers = ['Ensembl Gene ID', 'Chromosome', 'Strand', 'Exon Start (bp)', 'Exon End (bp)', 'Ensembl Exon ID', 'Constitutive Exon', 'Ensembl Transcript ID']
    for eti in exon_transcript_db:
        exonid = eti.ExonId(); transcript_id = eti.TranscriptId(); rank = eti.Rank()
        ei = exon_db[exonid]; seq_region_start = ei.SeqRegionStart(); seq_region_end = ei.SeqRegionEnd()
        ens_exon = exon_stable_id_db[exonid].StableId()  
        ti = transcript_db[transcript_id]; geneid = ti.GeneId(); ens_transcript = transcript_stable_id_db[transcript_id].StableId()  
        gi = gene_db[geneid]; seq_region_id = gi.SeqRegionId(); strand = gi.SeqRegionStrand()
        chr = seq_region_db[seq_region_id].Name()
        ens_gene = gene_stable_id_db[geneid].StableId()        
        values = [geneid,transcript_id,rank,ens_gene,chr,str(strand),str(seq_region_start),str(seq_region_end),ens_exon,'0',ens_transcript]
        temp_values_list.append(values)
        
    values_list=[]; temp_values_list.sort() ###Make sure the gene, transcripts and exons are grouped and ranked appropriately
    for values in temp_values_list: values_list.append(values[3:])
    temp_values_list=[]
    exportEnsemblTable(values_list,headers,output_dir)

def getDomainGenomicCoordinates(species,xref_db):
    ### Function for building Ensembl transcription annotation file
    output_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_ProteinFeatures_build'+ensembl_build+'.tab'
    headers = ['ID', 'AA_Start', 'AA_Stop', 'Start', 'Stop', 'Name', 'Interpro', 'Description']
    exon_transcript_db2={} ### exon_transcript_db has data stored as a list, so store as a ranked db
    values_list=[]
    for eti in exon_transcript_db:
        transcript_id = eti.TranscriptId(); rank = eti.Rank()
        try: exon_transcript_db2[transcript_id].append((rank,eti))
        except KeyError: exon_transcript_db2[transcript_id] = [(rank,eti)]

    interpro_annotation_db={}
    for xref_id in xref_db:
        interprot_id = xref_db[xref_id].DbprimaryAcc(); display_label = xref_db[xref_id].DisplayLabel()
        interpro_annotation_db[interprot_id] = display_label
        
    protein_coding_exon_db={}            
    for protein_id in translation_db:
        ti = translation_db[protein_id]; transcript_id = ti.TranscriptId(); 
        seq_start = ti.SeqStart(); start_exon_id = ti.StartExonId(); seq_end = ti.SeqEnd(); end_exon_id = ti.EndExonId()
        eti_list = exon_transcript_db2[transcript_id]; eti_list.sort()
        ti = transcript_db[transcript_id]; geneid = ti.GeneId(); gi = gene_db[geneid]; strand = gi.SeqRegionStrand() #Get strand 
        coding_exons = []; ce=0
        for (rank,eti) in eti_list:
            exonid = eti.ExonId()
            #print exonid,start_exon_id,end_exon_id
            if exonid == start_exon_id: ### Thus we are in the start exon for this transcript
                ei = exon_db[exonid]; genomic_exon_start = ei.SeqRegionStart(); genomic_exon_end = ei.SeqRegionEnd()
                exon_length = abs(genomic_exon_end-genomic_exon_start)
                coding_bp_in_exon = exon_length - seq_start; eti.setCodingBpInExon(coding_bp_in_exon)
                #print 'start',genomic_exon_start,genomic_exon_end,exon_length,seq_start;kill
                coding_exons.append(eti); ce+=1
            elif ce != 0: ###If we are in coding exons
                ei = exon_db[exonid]; genomic_exon_start = ei.SeqRegionStart(); genomic_exon_end = ei.SeqRegionEnd()
                exon_length = abs(genomic_exon_end-genomic_exon_start)
                eti.setCodingBpInExon(exon_length)
                coding_exons.append(eti)
            elif exonid == end_exon_id: ### Thus we are in the end exon for this transcript        
                ei = exon_db[exonid]; genomic_exon_start = ei.SeqRegionStart(); genomic_exon_end = ei.SeqRegionEnd()
                exon_length = abs(genomic_exon_end-genomic_exon_start)
                coding_bp_in_exon = exon_length - seq_end; eti.setCodingBpInExon(coding_bp_in_exon)
                coding_exons.append(eti); ce = 0
        protein_coding_exon_db[protein_id] = coding_exons

    interprot_match=0
    for protein_feature_id in protein_feature_db:
        pfi = protein_feature_db[protein_feature_id]; protein_id = pfi.TranslationId(); evalue = pfi.Evalue()
        ens_protein = translation_stable_id_db[protein_id].StableId()
        seq_start = pfi.SeqStart(); seq_end = pfi.SeqEnd(); hit_id = pfi.HitId() ### hit_id is the domain accession which maps to 'id' in interpro
        seq_start = seq_start*3; seq_end = seq_end*3 ###convert to transcript basepair positions
        coding_exons = protein_coding_exon_db[protein_id]
        cummulative_coding_length = 0; last_exon_cummulative_coding_length = 1

        if hit_id in interpro_db and evalue<1: ### Only analyze domain-level features that correspond to known protein feature IDs
            interpro_ac = interpro_db[hit_id].InterproAc(); interprot_match+=1
            if interpro_ac in interpro_annotation_db:
                interpro_name = interpro_annotation_db[interpro_ac] ###Built this annotation database using the xref database
                #print interpro_name,ens_protein,seq_start,seq_end
                for eti in coding_exons:
                    coding_bp_in_exon = eti.CodingBpInExon(); cummulative_coding_length += coding_bp_in_exon
                    if seq_start <= cummulative_coding_length and seq_start >= last_exon_cummulative_coding_length: ### Thus, domain starts in this exon
                        exonid = eti.ExonId(); ei = exon_db[exonid]; genomic_exon_start = ei.SeqRegionStart()
                        genomic_bp_exon_offset = seq_start - last_exon_cummulative_coding_length
                        genomic_domain_start = genomic_bp_exon_offset+genomic_exon_start-3 ### This is what we want! (minus 3 so that we start at the first bp of that codon, not the first of the next (don't count the starting coding as 3bp)
                        #print genomic_exon_start,last_exon_cummulative_coding_length,genomic_domain_start,genomic_bp_exon_offset;kill
                        #pfi.setGenomicStart(genomic_domain_start)
                    if seq_end <= cummulative_coding_length and seq_end >= last_exon_cummulative_coding_length: ### Thus, domain starts in this exon
                        exonid = eti.ExonId(); ei = exon_db[exonid]; genomic_exon_start = ei.SeqRegionStart()
                        genomic_bp_exon_offset = seq_end - last_exon_cummulative_coding_length
                        genomic_domain_end = genomic_bp_exon_offset+genomic_exon_start-1 ### This is what we want!
                        #pfi.setGenomicEnd(genomic_domain_end)
                    last_exon_cummulative_coding_length = cummulative_coding_length + 1
                values = [ens_protein,seq_start/3,seq_end/3,genomic_domain_start,genomic_domain_end,hit_id,interpro_ac,interpro_name]
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
        ens_gene = gene_stable_id_db[geneid].StableId()
        values = [ens_gene,description,symbol]
        values_list.append(values)
        ###Also, create a filter_db that allows for the creation of Ensembl-external gene db tables
    exportEnsemblTable(values_list,headers,output_dir)
    
def buildEnsemblExternalDBRelationshipTable(external_system,xref_db,object_xref_db,output_id_type,species):
    ### Get xref annotations (e.g., external geneID) for a set of Ensembl IDs (e.g. gene IDs)
    external_system = string.replace(external_system,'/','-') ###Some external relationship systems have a slash in them which will create unwanted directories
    output_dir = 'AltDatabase/'+external_system+'/'+species+'/'+species+'_Ensembl-'+external_system+'.txt'
    headers = ['Ensembl ID',external_system]; id_type_db={}
    id_type_db={}; gene_relationship_db={}; transcript_relationship_db={}; translation_relationship_db={}
    for xref_id in object_xref_db:
        for ox in object_xref_db[xref_id]:
            ens_numeric_id = ox.EnsemblId()
            if xref_id in xref_db:
                dbprimary_acc = xref_db[xref_id].DbprimaryAcc()
                external_db_id = xref_db[xref_id].ExternalDbId()
                ens_object_type = ox.EnsemblObjectType()
                if external_db_id in external_filter_db: ###Make sure other gene systems are not included
                    ### For testing determine the most likely system linked to by the xref (e.g. gene, transcript or translation ID).
                    try: id_type_db[ox.EnsemblObjectType()]+=1
                    except KeyError: id_type_db[ox.EnsemblObjectType()]=1
                        
                    if ens_object_type == 'Gene': gene_relationship_db[ens_numeric_id,dbprimary_acc]=[]
                    if ens_object_type == 'Transcript': transcript_relationship_db[ens_numeric_id,dbprimary_acc]=[]
                    if ens_object_type == 'Translation': translation_relationship_db[ens_numeric_id,dbprimary_acc]=[]
                            
    ids=['ID types linked to '+external_system+' are:']
    for id_type in id_type_db: ### Tells us which ID types are most connected to a given external reference ID (don't often know)
        ids+=[str(id_type),': ',str(id_type_db[id_type]),'\t']
    ids = string.join(ids,''); print ids

    values_list = convertBetweenEnsemblIDTypes(output_id_type,transcript_relationship_db,translation_relationship_db,gene_relationship_db)
    exportEnsemblTable(values_list,headers,output_dir)
            
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
            geneid = transcript_db[ens_numeric_id].GeneId(); ens_id = gene_stable_id_db[geneid].StableId()
        elif output_id_type == 'Transcription': ens_id = transcript_stable_id_db[ens_numeric_id].StableId()
        elif output_id_type == 'Translation':
            try: protein_id = transcript_to_protein_db[ens_numeric_id]; ens_id = translation_stable_id_db[protein_id].StableId()
            except KeyError: null = []
        try: values = [ens_id,dbprimary_acc]; values_list.append(values)
        except NameError: null = []
        
    for (ens_numeric_id,dbprimary_acc) in translation_relationship_db:
        if output_id_type == 'Gene':
            transcript_id = translation_db[ens_numeric_id].TranscriptId()
            geneid = transcript_db[transcript_id].GeneId(); ens_id = gene_stable_id_db[geneid].StableId()
        elif output_id_type == 'Transcription':
            transcript_id = translation_db[ens_numeric_id].TranscriptId()
            ens_id = transcript_stable_id_db[transcript_id].StableId()
        elif output_id_type == 'Translation': ens_id = translation_stable_id_db[ens_numeric_id].StableId()
        values = [ens_id,dbprimary_acc]; values_list.append(values)
        
    for (ens_numeric_id,dbprimary_acc) in gene_relationship_db:
        if output_id_type == 'Gene':
            ens_id = gene_stable_id_db[ens_numeric_id].StableId()
            values = [ens_id,dbprimary_acc]; values_list.append(values)
        elif output_id_type == 'Transcription':
            transcript_ids = gene_to_transcript_db[ens_numeric_id]
            for transcript_id in transcript_ids:
                ens_id = transcript_stable_id_db[transcript_id].StableId()
                values = [ens_id,dbprimary_acc]; values_list.append(values)
        elif output_id_type == 'Translation':
            transcript_ids = gene_to_transcript_db[ens_numeric_id]
            for transcript_id in transcript_ids:
                try: ### Translate between transcripts to protein IDs
                    protein_id = transcript_to_protein_db[ens_numeric_id]
                    ens_id = translation_stable_id_db[protein_id].StableId()
                    values = [ens_id,dbprimary_acc]; values_list.append(values)
                except KeyError: null = []
    values_list = unique.unique(values_list)    
    return values_list
                
def exportEnsemblTable(values_list,headers,output_dir):
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
        elif header == "interpro_ac": self.interpro_ac = value
        elif header == "xref_id": self.xref_id = value
        elif header == "evalue": self.evalue = float(value)
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
        
def importEnsemblSQLInfo(configType):
    filename = 'Config/EnsemblSQL.txt'
    fn=filepath(filename); sql_file_db={}; sql_group_db={}; x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        filename, url_type, importGroup, configGroup, key, values, comments = string.split(data,'\t')
        if x==0: x=1
        else:
            sq = EnsemblSQLInfo(filename, url_type, importGroup, key, values, comments)
            sql_file_db[filename] = sq
            ### To conserve memory, only import necessary files for each run (results in 1/5th the memory usage)
            if configType == 'Basic' and configGroup == 'Basic': proceed = 'yes'
            elif configType == 'Basic': proceed = 'no'
            else: proceed = 'yes'
            if proceed == 'yes':
                try: sql_group_db[importGroup].append(filename)
                except KeyError: sql_group_db[importGroup] = [filename]
    return sql_file_db,sql_group_db

def importEnsemblSQLFiles(ensembl_sql_dir,ensembl_sql_description_dir,sql_group_db,sql_file_db,output_dir,import_group,force):
    if import_group == 'Primary':
        global exon_db; global transcript_db; global translation_stable_id_db
        global exon_transcript_db; global transcript_stable_id_db; global gene_db
        global exon_stable_id_db; global translation_db; global gene_stable_id_db
        global protein_feature_db; global interpro_db; global external_synonym_db
        global external_db_db; global key_filter_db; global external_filter_db
        global seq_region_db
        key_filter_db={}; external_filter_db={}        

    ###Generic function for importing all EnsemblSQL tables based on the SQLDescription table        
    for filename in sql_group_db['Description']: ### Gets the SQL description table
        sql_filepath = updateFiles(ensembl_sql_description_dir,output_dir,filename,force)
        sql_file_db = importSQLDescriptions(sql_filepath,sql_file_db)
    for filename in sql_group_db[import_group]:
        sql_filepath = updateFiles(ensembl_sql_dir,output_dir,filename,force)
        try: key_value_db = importPrimaryEnsemblSQLTables(sql_filepath,filename,sql_file_db[filename])
        except IOError:
            sql_filepath = updateFiles(ensembl_sql_dir,output_dir,filename,'yes')
            key_value_db = importPrimaryEnsemblSQLTables(sql_filepath,filename,sql_file_db[filename])
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
        elif filename == "xref.txt": xref_db = key_value_db
        elif filename == "object_xref": object_xref = key_value_db
        elif filename == "external_db.txt": external_db_db = key_value_db
        elif filename == "object_xref.txt": object_xref_db = key_value_db
        elif filename == "seq_region.txt": seq_region_db = key_value_db
        else: ###Shouldn't occur, unless we didn't account for a file type
            print 'Warning!!! A file has been imported which does not exist in the program space'
            print 'Filename =',filename;sys.exit()
    if import_group == 'Primary':
        key_filter_db={}
        for geneid in gene_db:
            gi = gene_db[geneid]; display_xref_id = gi.DisplayXrefId()
            key_filter_db[display_xref_id]=[]
    elif import_group == 'Xref':
        return xref_db
    elif import_group == 'Object-Xref':
        return object_xref_db
        
def updateFiles(ensembl_sql_dir,output_dir,filename,force):
    if force == 'yes': ### Download the files, rather than re-use existing
        ftp_url = ensembl_sql_dir+filename + '.gz'
        gz_filepath, status = update.download(ftp_url,output_dir,'')
        if status == 'not-removed':
            try: os.remove(gz_filepath) ### Not sure why this works now and not before
            except OSError: status = status
        sql_filepath = gz_filepath[:-3]
    else: sql_filepath = output_dir + filename
    return sql_filepath

def importPrimaryEnsemblSQLTables(sql_filepath,filename,sfd):
    fn=filepath(sql_filepath)
    index=0; key_value_db={}
    if len(key_filter_db)>0: key_filter = 'yes'
    else: key_filter = 'no'
    
    if len(external_filter_db)>0: external_filter = 'yes'
    else: external_filter = 'no'

    try:
        if len(external_xref_key_db)>0: external_xref_key_filter = 'yes'
        else: external_xref_key_filter = 'no'
    except NameError: external_xref_key_filter = 'no'

    index_db = sfd.IndexDB(); key_name = sfd.Key()
    if len(key_name)<1: key_value_db=[] ### No key, so store data in a list
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line); data = string.split(data,'\t')
        ese = EnsemblSQLEntryData()
        for index in index_db:
            header_name,value_type = index_db[index]
            value = data[index]
            try:
                if value_type == 'integer': value = int(value) # Integers will take up less space in memory
            except ValueError:
                if value == '\\N': value = value
                else: print filename,[header_name],[value];kill
            ###Although we are setting each index to value, the distinct headers will instruct
            ###EnsemblSQLEntryData to assign the value to a distinct object
            if header_name != key_name:
                ese.setSQLValue(header_name,value)
            else:
                key = value
        ### Filtering primarily used for all Xref, since this database is very large
        #if filename == 'xref.txt':print len(key_name), key_filter,external_filter, external_xref_key_filter;kill
        if len(key_name)<1: key_value_db.append(ese)
        elif key_filter == 'no' and external_filter == 'no': key_value_db[key] = ese
        elif external_xref_key_filter == 'yes':
            if key in external_xref_key_db:
                try: key_value_db[key].append(ese)
                except KeyError: key_value_db[key] = [ese]
        elif key in key_filter_db: key_value_db[key] = ese
        elif external_filter == 'yes': ### For example, if UniGene's dbase ID is in the external_filter_db (when parsing xref)
            try:
                #print key, [ese.ExternalDbId()], [ese.DbprimaryAcc()], [ese.DisplayLabel()];kill
                #if key == 1214287: print [ese.ExternalDbId()], [ese.DbprimaryAcc()], [ese.DisplayLabel()]
                if ese.ExternalDbId() in external_filter_db: key_value_db[key] = ese
            except AttributeError:
                print len(external_filter_db),len(key_filter_db)
                print 'index_db',index_db,'\n'
                print 'key_name',key_name
                print len(external_filter_db)
                print len(key_value_db);kill

    #key_filter_db={}; external_filter_db={}; external_xref_key_db={}
    print "Extracted",len(key_value_db),"entries for",filename
    return key_value_db
        
def importSQLDescriptions(filename,sql_file_db):
    fn=filepath(filename)
    index_db={}; index=0
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        if 'CREATE TABLE' in data:
            ###New file descriptions
            file_broken = string.split(data,'`'); filename = file_broken[1]+".txt"
            if filename in sql_file_db:
                sql_data = sql_file_db[filename]
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
        else: index+=1
    if len(index_db)>0:
        sql_data.setIndexDB(index_db)
    return sql_file_db

def clearall():
    all = [var for var in globals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all: del globals()[var]
    
if __name__ == '__main__':
    analysisType = 'AltAnalyzeDBs'; analysisType = 'GeneAndExternal'
    
    #species='Hs'; force = 'no'; configType = 'Basic'; analysisType = 'ExternalOnly'; externalDBName = 'GO'
    #buildEnsemblRelationalTablesFromSQL(species,configType,analysisType,externalDBName,force)

    species='Mm'; force = 'yes'; configType = 'Advanced'; analysisType = 'AltAnalyzeDBs'; externalDBName = ''
    buildEnsemblRelationalTablesFromSQL(species,configType,analysisType,externalDBName,force)

    