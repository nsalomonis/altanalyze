###FeatureAlignment
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

import sys, string
import os.path
import unique
import ExonAnalyze_module
import copy

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt"or entry[-4:] == ".tab" or entry[-4:] == ".csv" or '.fa' in entry: dir_list2.append(entry)
    return dir_list2

class GrabFiles:
    def setdirectory(self,value):
        self.data = value
    def display(self):
        print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        file_dir,file = getDirectoryFiles(self.data,str(search_term))
        if len(file)<1: print search_term,'not found',self.data
        return file_dir
    
def getDirectoryFiles(import_dir, search_term):
    exact_file = ''; exact_file_dir=''
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir: exact_file_dir = affy_data_dir; exact_file = data
    return exact_file_dir,exact_file

########## End generic file import ##########

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

class ProteinFunctionalSeqData:
    def __init__(self, protein_accession,primary_annotation, secondary_annotation, ft_start_pos, ft_end_pos, ft_seq):
        self._protein_accession = protein_accession; self._primary_annotation = primary_annotation
        self._secondary_annotation = secondary_annotation; self._ft_start_pos = ft_start_pos
        self._ft_end_pos = ft_end_pos; self._ft_seq = ft_seq
    def ProteinID(self): return self._protein_accession
    def PrimaryAnnot(self): return self._primary_annotation
    def SecondaryAnnot(self): return self._secondary_annotation
    def CombinedAnnot(self): return self.PrimaryAnnot()+'-'+self.SecondaryAnnot()
    def DomainStart(self): return int(self._ft_start_pos)
    def DomainEnd(self): return int(self._ft_end_pos)
    def DomainSeq(self): return self._ft_seq
    def DomainLen(self):
        domain_len = self.DomainEnd()-self.DomainStart()
        return domain_len
    def SummaryValues(self):
        output = self.PrimaryAnnot()+'|'+self.SecondaryAnnot()
        return output
    def __repr__(self): return self.SummaryValues()

class FullProteinSeqData:
    def __init__(self, primary_id, secondary_ids, sequence, type):
        self._primary_id = primary_id; self._secondary_ids = secondary_ids
        self._sequence = sequence; self._type = type
    def PrimaryID(self): return self._primary_id
    def SecondaryID(self): return self._secondary_ids
    def Sequence(self): return self._sequence
    def SequenceLength(self): return len(self._sequence)
    def AccessionType(self): return self._type
    def SummaryValues(self):
        output = self.PrimaryID()+'|'+self.SecondaryID()+'|'+self.AccessionType()
        return output
    def __repr__(self): return self.SummaryValues()    

######Import ArrayID to Protein/Gene Relationships
def import_arrayid_ensembl(filename):
    fn=filepath(filename); x = 0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        gene_id,ensembl_gene_id = string.split(data,'\t')
        try: ensembl_arrayid_db[ensembl_gene_id].append(gene_id)
        except KeyError: ensembl_arrayid_db[ensembl_gene_id] = [gene_id]

def findDomainsByGenomeCoordinates(species,array_type):
    ### Grab Ensembl relationships from a custom Ensembl Perl script or BioMart
    protein_relationship_file,protein_feature_file,protein_seq_fasta = getEnsemblRelationshipDirs(species)
    ens_transcript_protein_db = importEnsemblRelationships(protein_relationship_file,'transcript') ### From Perl script to Ensembl API
    ens_protein_gene_db = importEnsemblRelationships(protein_relationship_file,'gene') ### From Perl script to Ensembl API
    exon_protein_db = importEnsExonStructureDataSimple(species,ens_transcript_protein_db) ### From BioMart

    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'    
    first_last_exon_coord_db = importEnsExonStructureDataCustom(filename,species,{})
    ### Add UCSC transcript data to ens_transcript_exon_db and ens_gene_transcript_db
    filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_mrna.txt' ### Use the non-filtered database to propperly analyze exon composition 
    first_last_exon_coord_db = importEnsExonStructureDataCustom(filename,species,first_last_exon_coord_db)
 
    if array_type == 'exon': ens_probeset_file = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_probesets.txt"    
    else: ens_probeset_file = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_"+array_type+"_probesets.txt"
    protein_probeset_db,gene_probeset_db = importSplicingAnnotationDatabase(ens_probeset_file,exon_protein_db) ### Derived from ExonArrayEnsemblRules
    matchEnsemblDomainCoordinates(protein_feature_file,species,array_type,protein_probeset_db,ens_protein_gene_db,gene_probeset_db,first_last_exon_coord_db)
    
def matchEnsemblDomainCoordinates(filename,species,array_type,protein_probeset_db,ens_protein_gene_db,gene_probeset_db,first_last_exon_coord_db):
    fn=filepath(filename); x=0; probeset_domain_match={}; probeset_domain_indirect_match={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x=1 ###Don't extract the headers
        else:
            coor_list=[]
            ensembl_prot, aa_start, aa_stop, genomic_start, genomic_stop, name, interpro_id, description = string.split(data,'\t')
            coor_list.append(int(genomic_start)); coor_list.append(int(genomic_stop)); coor_list.sort()
            genomic_start, genomic_stop = coor_list
            if ensembl_prot in protein_probeset_db and len(description)>1:
                probeset_data = protein_probeset_db[ensembl_prot]
                for sad in probeset_data:
                    proceed = 'no'
                    if ((genomic_start <= sad.ProbeStart()) and (genomic_stop >= sad.ProbeStart())) and ((genomic_stop >= sad.ProbeStop()) and (genomic_start <= sad.ProbeStop())):
                        overlap_area = int(abs(sad.ProbeStop()-sad.ProbeStart())/3)
                        proceed = 'yes'
                    elif ((genomic_start <= sad.ProbeStart()) and (genomic_stop >= sad.ProbeStart())):
                        overlap_area = int(abs(genomic_stop - sad.ProbeStart())/3)
                        proceed = 'yes'
                    elif ((genomic_stop >= sad.ProbeStop()) and (genomic_start <= sad.ProbeStop())):
                        overlap_area = int(abs(sad.ProbeStop() - genomic_start)/3)
                        proceed = 'yes'
                    if proceed == 'yes':
                        #if sad.Probeset() == '3217131': print ensembl_prot, sad.ProbeStart(),sad.ProbeStop(),genomic_start,genomic_stop,interpro_id,description;kill

                        ipd = description+'-'+interpro_id
                        try: probeset_domain_match[sad.Probeset()].append(ipd)
                        except KeyError: probeset_domain_match[sad.Probeset()] = [ipd]
            ### Grab all gene associated probesets that are not in the first or last exon of any transcript (make a db of the first and last exon coordiantes of all analyzed transcripts also no UTR exons OR just remove those that have an alt-N, Alt-C or UTR annotation)
            if ensembl_prot in ens_protein_gene_db:
                ens_gene = ens_protein_gene_db[ensembl_prot]
                if ens_gene in gene_probeset_db:
                    probeset_data = gene_probeset_db[ens_gene]
                    for sad in probeset_data:
                        if ((genomic_start <= sad.ProbeStart()) and (genomic_stop >= sad.ProbeStart())) or ((genomic_stop >= sad.ProbeStop()) and (genomic_start <= sad.ProbeStop())):  
                            #if sad.Probeset() == '3217131': print ensembl_prot, sad.ProbeStart(),sad.ProbeStop(),genomic_start,genomic_stop,interpro_id,description;kill

                            ipd = description+'-'+interpro_id
                            try: probeset_domain_indirect_match[sad.Probeset()].append(ipd)
                            except KeyError: probeset_domain_indirect_match[sad.Probeset()] = [ipd]                  

    probeset_domain_indirect_match = eliminateRedundant(probeset_domain_indirect_match)
    probeset_domain_match = eliminateRedundant(probeset_domain_match) ###Remove redundant matches
    print len(probeset_domain_match),'probesets with associated protein domains'
    
    probeset_domain_indirect_match2={}; probeset_domain_indirect_match3={}
    for probeset in probeset_domain_indirect_match:
        if probeset not in probeset_domain_match: ### Only have probesets that don't directly map to domains
            probeset_domain_indirect_match2[probeset] = probeset_domain_indirect_match[probeset]

    probesets_to_exclude={}
    for gene in gene_probeset_db:
        ### Remove probesets from database that overlap with the first or last exon of any associated transcript (not high confidence for domain alignment)
        probeset_data = gene_probeset_db[gene]
        for sad in probeset_data:
            if sad.Probeset() in probeset_domain_indirect_match2:
                exon_coordinates = first_last_exon_coord_db[gene]
                for exon_coor in exon_coordinates:
                    exon_coor.sort()
                    genomic_start,genomic_stop = exon_coor
                    if ((genomic_start <= sad.ProbeStart()) and (genomic_stop >= sad.ProbeStart())) and ((genomic_stop >= sad.ProbeStop()) and (genomic_start <= sad.ProbeStop())):  
                        probesets_to_exclude[sad.Probeset()] = []
                    #if sad.Probeset() == '3217131': print gene, sad.ProbeStart(),sad.ProbeStop(),genomic_start,genomic_stop;kill

    for probeset in probeset_domain_indirect_match2:
        if probeset not in probesets_to_exclude:
            probeset_domain_indirect_match3[probeset] = probeset_domain_indirect_match2[probeset]

    print len(probeset_domain_indirect_match2), len(probesets_to_exclude), len(probeset_domain_indirect_match3)   
    
    exportProbesetDomainMappings(species,array_type,'',probeset_domain_match)
    exportProbesetDomainMappings(species,array_type,'indirect_',probeset_domain_indirect_match3)
    
def importEnsExonStructureDataCustom(filename,species,first_last_exon_coord_db):
    fn=filepath(filename); x=0; ens_transcript_exon_db={}; ens_gene_transcript_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            gene, chr, strand, exon_start, exon_end, ens_exonid, constitutive_exon, ens_transcriptid = t
            exon_start = int(exon_start); exon_end = int(exon_end)
            try: ens_transcript_exon_db[ens_transcriptid].append([exon_start,exon_end])
            except KeyError: ens_transcript_exon_db[ens_transcriptid] = [[exon_start,exon_end]]
            ens_gene_transcript_db[ens_transcriptid] = gene

    for transcript in ens_transcript_exon_db:
        gene = ens_gene_transcript_db[transcript]
        try: first_last_exon_coord_db[gene].append(ens_transcript_exon_db[transcript][0])
        except KeyError: first_last_exon_coord_db[gene] = [ens_transcript_exon_db[transcript][0]]
        try: first_last_exon_coord_db[gene].append(ens_transcript_exon_db[transcript][-1])
        except KeyError: first_last_exon_coord_db[gene] = [ens_transcript_exon_db[transcript][-1]]
    first_last_exon_coord_db = eliminateRedundant(first_last_exon_coord_db)
    
    return first_last_exon_coord_db

def exportProbesetDomainMappings(species,array_type,indirect_mapping,probeset_domain_match):            
    import export
    export_file = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_"+indirect_mapping+"domain_aligning_probesets.txt"                       
    data = export.createExportFile(export_file,"AltDatabase/"+species+"/"+array_type)
    data.write('Probeset\tInterPro-Description\n')
    for probeset in probeset_domain_match:
        domain_info_list = probeset_domain_match[probeset]
        for ipd in domain_info_list: data.write(probeset+'\t'+ipd+'\n')
    data.close()
    print "Direct probeset to domain associations exported to:",export_file
    
def importSplicingAnnotationDatabase(filename,exon_protein_db):
    fn=filepath(filename); x=0; protein_probeset_db={}; gene_probeset_db={}
    for line in open(fn,'rU').xreadlines():             
        probeset_data = cleanUpLine(line)  #remove endline
        if x == 0: x = 1
        else:
            t=string.split(probeset_data,'\t'); probeset=t[0]; exon_id = t[1]; ens_gene=t[2]; probeset_start=t[6]; probeset_stop=t[7]; external_exonid=t[10]; splicing_event=t[-2]; 
            sad = SplicingAnnotationData(probeset,probeset_start,probeset_stop)
            if 'U' not in exon_id:
                #if 'alt-C-term' not in splicing_event and 'alt-N-term' not in splicing_event and 'altPromoter' not in splicing_event:
                try: gene_probeset_db[ens_gene].append(sad)
                except KeyError: gene_probeset_db[ens_gene] = [sad]
            if 'ENS' in external_exonid: ### We only need to examine probesets linked to Ensembl transcripts
                external_exonid_list = string.split(external_exonid,'|'); ens_exonid_list=[]
                for exon in external_exonid_list:
                    if 'ENS' in exon: ens_exonid_list.append(exon)
                for ens_exon in ens_exonid_list:
                    if ens_exon in exon_protein_db:
                        ens_proteins = exon_protein_db[ens_exon]
                        for ens_protein_id in ens_proteins:
                            try: protein_probeset_db[ens_protein_id].append(sad)
                            except KeyError: protein_probeset_db[ens_protein_id] = [sad]
                                            
    print len(protein_probeset_db),'Ensembl proteins with associated probesets'
    return protein_probeset_db,gene_probeset_db

class SplicingAnnotationData:
    def __init__(self,probeset,start,stop):
        self._probeset = probeset; self._start = start; self._stop = stop
    def Probeset(self): return self._probeset
    def ProbeStart(self): return int(self._start)
    def ProbeStop(self): return int(self._stop)
    
def importEnsExonStructureDataSimple(species,ens_transcript_protein_db):
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    fn=filepath(filename); x=0; exon_protein_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            gene, chr, strand, exon_start, exon_end, ens_exonid, constitutive_exon, ens_transcriptid = t
            if ens_transcriptid in ens_transcript_protein_db:
                ens_protein_id = ens_transcript_protein_db[ens_transcriptid]
                try: exon_protein_db[ens_exonid].append(ens_protein_id)
                except KeyError: exon_protein_db[ens_exonid]=[ens_protein_id]
    print len(exon_protein_db),'Ensembl exons with associated protein IDs'
    return exon_protein_db

def import_ensembl_ft_data(species,filename,ensembl_arrayid_db,array_type):
    """Over the lifetime of this program, the inputs for protein sequences and relationships have changed.
    To support multiple versions, this function now routes the data to two different possible functions,
    grabbing the same type of data (InterPro relationships and protein sequence) from different sets of files"""
    try: ensembl_protein_seq_db,ensembl_ft_db,domain_gene_counts = importCombinedEnsemblFTdata(filename,ensembl_arrayid_db,array_type)
    except IOError:
        protein_relationship_file,protein_feature_file,protein_seq_fasta = getEnsemblRelationshipDirs(species)
        ensembl_protein_seq_db = importEnsemblProtSeqFasta(protein_seq_fasta)
        ensembl_protein_gene_db = importEnsemblRelationships(protein_relationship_file,'gene')
        ensembl_ft_db, domain_gene_counts = importEnsemblFTdata(protein_feature_file,ensembl_arrayid_db,array_type,ensembl_protein_seq_db,ensembl_protein_gene_db)
        
    return ensembl_protein_seq_db,ensembl_ft_db,domain_gene_counts

def getEnsemblRelationshipDirs(species):
    import_dir = '/AltDatabase/ensembl/'+species
    m = GrabFiles(); m.setdirectory(import_dir)
    protein_relationship_file = m.searchdirectory(species+'_Ensembl_Protein_')
    protein_seq_fasta = m.searchdirectory('pep.all')
    protein_feature_file = m.searchdirectory(species+'_ProteinFeatures_')
    return protein_relationship_file,protein_feature_file,protein_seq_fasta

def importEnsemblProtSeqFasta(filename):
    print "Begining generic fasta import of",filename
    fn=filepath(filename); ensembl_protein_seq_db={}; sequence = ''
    for line in open(fn,'r').xreadlines():
        try: data, newline= string.split(line,'\n')
        except ValueError: continue
        try:
            if data[0] == '>':
                if len(sequence) > 0:
                    seq_data = FullProteinSeqData(ensembl_prot,[ensembl_prot],sequence,'EnsProt')
                    ensembl_protein_seq_db[ensembl_prot] = seq_data
                ### Parse new line                      
                t= string.split(data[1:],' '); sequence=''
                ensembl_prot = t[0]
        except IndexError: continue
        try:
            if data[0] != '>': sequence = sequence + data
        except IndexError:  continue
    seq_data = FullProteinSeqData(ensembl_prot,[ensembl_prot],sequence,'EnsProt')
    ensembl_protein_seq_db[ensembl_prot] = seq_data
    return ensembl_protein_seq_db

def importEnsemblRelationships(filename,type):
    fn=filepath(filename); ensembl_protein_gene_db={}; x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x=1 ###Don't extract the headers
        else:
            ensembl_gene, ensembl_transcript, ensembl_protein = string.split(data,'\t')
            if type == 'gene': ensembl_protein_gene_db[ensembl_protein] =  ensembl_gene
            if type == 'transcript': ensembl_protein_gene_db[ensembl_transcript] =  ensembl_protein
    return ensembl_protein_gene_db

def importEnsemblFTdata(filename,ensembl_arrayid_db,array_type,ensembl_protein_seq_db,ensembl_protein_gene_db):
    global arrayid_ensembl_protein_db; arrayid_ensembl_protein_db={}; x=0
    fn=filepath(filename); ensembl_ft_db = {}; ensembl_ft_summary_db = {}# Use the last database for summary statistics
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x=1 ###Don't extract the headers
        else: 
            ensembl_prot, aa_start, aa_stop, start, stop, name, interpro_id, description = string.split(data,'\t')
            sd = ensembl_protein_seq_db[ensembl_prot]; protein_sequence = sd.Sequence()
            ensembl_gene = ensembl_protein_gene_db[ensembl_prot]
            ft_start_pos = aa_start; ft_end_pos = aa_stop
            """Below code is ALMOST identical to importCombinedEnsemblFTdata (original code commented out)"""
            ###If array_type is exon, ensembl is used as the primary gene ID and thus no internal arrayid is needed. Get only Ensembl's that are currently being analyzed (for over-representation analysis)
            if ensembl_gene in ensembl_arrayid_db:
                id_list = ensembl_arrayid_db[ensembl_gene]
                for gene_id in id_list: 
                    try: arrayid_ensembl_protein_db[gene_id].append(ensembl_prot)
                    except KeyError: arrayid_ensembl_protein_db[gene_id] = [ensembl_prot]
            #for entry in ft_info_list:
            """try: peptide_start_end, gene_start_end, feature_source, interpro_id, description = string.split(entry,' ')
            except ValueError: continue
            ###142-180 3015022-3015156 Pfam IPR002050 Env_polyprotein
            ft_start_pos, ft_end_pos = string.split(peptide_start_end,'-')"""
            pos1 = int(ft_start_pos); pos2 = int(ft_end_pos)
            ft_length = pos2-pos1
            if ft_length > 6: pos_1 = pos1; pos_2 = pos2
            else:
                if ft_length < 3: pos_1 = pos1 - 3; pos_2 = pos2 + 3
                else: pos_1 = pos1 - 1; pos_2 = pos2 + 1
            sequence_fragment = protein_sequence[pos_1:pos_2]  ###We will search for this sequence, so have this expanded if too small (see above code)
            if len(description)>1 or len(interpro_id)>1:
                #ft_info = ProteinFunctionalSeqData(ensembl_prot,description,interpro_id,pos1,pos2,sequence_fragment)
                ft_info = ensembl_prot,description,interpro_id,pos1,pos2,sequence_fragment ###don't store as an instance yet... wait till we eliminate duplicates
                if ensembl_gene in ensembl_arrayid_db:  ###If the ensembl gene is connected to microarray identifiers
                    arrayids = ensembl_arrayid_db[ensembl_gene]
                    for arrayid in arrayids: ###This file differs in structure to the UniProt data 
                        try: ensembl_ft_db[arrayid].append(ft_info)
                        except KeyError: ensembl_ft_db[arrayid] = [ft_info]
                        
    ensembl_ft_db2 = {}                        
    ensembl_ft_db = eliminateRedundant(ensembl_ft_db) ###duplicate interprot information is typically present
    for arrayid in ensembl_ft_db:
        for ft_info in ensembl_ft_db[arrayid]:
            ensembl_prot,description,interpro_id,pos1,pos2,sequence_fragment = ft_info
            ft_info2 = ProteinFunctionalSeqData(ensembl_prot,description,interpro_id,pos1,pos2,sequence_fragment)
            try: ensembl_ft_db2[arrayid].append(ft_info2)
            except KeyError: ensembl_ft_db2[arrayid] = [ft_info2]
            
    domain_gene_counts = summarizeEnsDomainData(ensembl_ft_db2)
    return ensembl_ft_db2,domain_gene_counts
    
def importCombinedEnsemblFTdata(filename,ensembl_arrayid_db,array_type):
    global arrayid_ensembl_protein_db; arrayid_ensembl_protein_db={}
    fn=filepath(filename); x = 0; ensembl_ft_db = {}; ensembl_ft_summary_db = {}# Use the last database for summary statistics
    ensembl_protein_seq_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 1:
            try: ensembl_gene, chr, mgi, uniprot, ensembl_prot, protein_sequence, position_info = string.split(data,'\t')
            except ValueError: continue
            ft_info_list = string.split(position_info,' | ')
            seq_data = FullProteinSeqData(ensembl_prot,[ensembl_prot],protein_sequence,'EnsProt')
            ensembl_protein_seq_db[ensembl_prot] =  seq_data   ###use this db as input for the UniProt exon based ft search below
            ###If array_type is exon, ensembl is used as the primary gene ID and thus no internal arrayid is needed. Get only Ensembl's that are currently being analyzed (for over-representation analysis)
            if ensembl_gene in ensembl_arrayid_db:
                id_list = ensembl_arrayid_db[ensembl_gene]
                for gene_id in id_list: 
                    try: arrayid_ensembl_protein_db[gene_id].append(ensembl_prot)
                    except KeyError: arrayid_ensembl_protein_db[gene_id] = [ensembl_prot]
            for entry in ft_info_list:
                try: peptide_start_end, gene_start_end, feature_source, interpro_id, description = string.split(entry,' ')
                except ValueError: continue
                ###142-180 3015022-3015156 Pfam IPR002050 Env_polyprotein
                ft_start_pos, ft_end_pos = string.split(peptide_start_end,'-')
                pos1 = int(ft_start_pos); pos2 = int(ft_end_pos)
                ft_length = pos2-pos1
                if ft_length > 6: pos_1 = pos1; pos_2 = pos2
                else:
                    if ft_length < 3: pos_1 = pos1 - 3; pos_2 = pos2 + 3
                    else: pos_1 = pos1 - 1; pos_2 = pos2 + 1
                sequence_fragment = protein_sequence[pos_1:pos_2]  ###We will search for this sequence, so have this expanded if too small (see above code)
                if len(description)>1 or len(interpro_id)>1:
                    ft_info = ProteinFunctionalSeqData(ensembl_prot,description,interpro_id,pos1,pos2,sequence_fragment)
                    if ensembl_gene in ensembl_arrayid_db:  ###If the ensembl gene is connected to microarray identifiers
                        arrayids = ensembl_arrayid_db[ensembl_gene]
                        for arrayid in arrayids: ###This file differs in structure to the UniProt data 
                            try: ensembl_ft_db[arrayid].append(ft_info)
                            except KeyError: ensembl_ft_db[arrayid] = [ft_info]
        else:
            if data[0:6] == 'GeneID': x = 1
            
    domain_gene_counts = summarizeEnsDomainData(ensembl_ft_db)
    return ensembl_protein_seq_db,ensembl_ft_db,domain_gene_counts

def summarizeEnsDomainData(ensembl_ft_db):
    """This is a function because the data can be extracted from different functions, using different file formats"""
    ensembl_ft_db = eliminateRedundant(ensembl_ft_db)
    domain_gene_counts = {}; domain_gene_counts2 = {}
    ###Count the number of domains present in all genes (count a domain only once per gene)
    for gene in ensembl_ft_db:
        for ft_info in ensembl_ft_db[gene]:
            try: domain_gene_counts[ft_info.PrimaryAnnot(),ft_info.SecondaryAnnot()].append(gene)
            except KeyError: domain_gene_counts[ft_info.PrimaryAnnot(),ft_info.SecondaryAnnot()] = [gene]
    domain_gene_counts = eliminateRedundant(domain_gene_counts)
    
    for (primary,secondary) in domain_gene_counts:
        if len(secondary)>0: key = primary+'-'+secondary
        else: key = primary
        domain_gene_counts2[key] = len(domain_gene_counts[(primary,secondary)])
    domain_gene_counts = domain_gene_counts2
    print "Number of Ensembl genes, linked to array genes with domain annotations:",len(ensembl_ft_db)
    print "Number of Ensembl domains:",len(domain_gene_counts)
    return domain_gene_counts

def import_arrayid_uniprot(filename):
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        arrayid, uniprot = string.split(data,'\t')
        try: arrayid_uniprot_db[arrayid].append(uniprot)
        except KeyError: arrayid_uniprot_db[arrayid] = [uniprot]
        try: uniprot_arrayid_db[uniprot].append(arrayid)
        except KeyError: uniprot_arrayid_db[uniprot] = [arrayid]
                  
def importUniProtSeqeunces(species,ensembl_arrayid_db,array_type):
    filename = 'AltDatabase/uniprot/'+species+'/'+'uniprot_sequence.txt'           
    fn=filepath(filename); uniprot_protein_seq_db = {}; external_transcript_to_uniprot_protein_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        id=t[0];ac=t[1];ensembls=t[4];seq=t[2];type=t[6];unigenes=t[7];embls=t[9]
        ac=string.split(ac,','); ensembls=string.split(ensembls,','); embls=string.split(embls,','); unigenes=string.split(unigenes,',')
        y = FullProteinSeqData(id,ac,seq,type)
        if type=='swissprot': uniprot_protein_seq_db[id] = y
        if array_type == 'exon':
            for ensembl in ensembls:
                if len(ensembl)>0 and ensembl in ensembl_arrayid_db:  ###remove genes not being analyzed now
                    ###This database replaces the empty arrayid_uniprot_db
                    try: arrayid_uniprot_db[ensembl].append(id)
                    except KeyError: arrayid_uniprot_db[ensembl] = [id]
                    try: uniprot_arrayid_db[id].append(ensembl)
                    except KeyError: uniprot_arrayid_db[id] = [ensembl]
    return uniprot_protein_seq_db
    return uniprot_protein_seq_db

######## End - Derive protein predictions for Exon array probesets
    
def import_uniprot_ft_data(species,filename,domain_gene_counts,ensembl_arrayid_db,array_type):
    uniprot_protein_seq_db = importUniProtSeqeunces(species,ensembl_arrayid_db,array_type)

    filename = 'AltDatabase/uniprot/'+species+'/'+'uniprot_feature_file.txt'
    fn=filepath(filename); uniprot_ft_db = {}        
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        try:  primary_uniprot_id,ac,ft,pos1,pos2,annotation = t
        except ValueError:
            try: primary_uniprot_id,ac,ft,pos1,pos2 = t
            except ValueError:
                ###Not sure why, but an extra \t is in at least one description.
                primary_uniprot_id = t[0]; ac=t[1]; ft=t[2];pos1=t[3];pos2=t[4];annotation=string.join(t[-2:])
            try: pos2,annotation = string.split(pos2,'/')
            except ValueError: annotation = ''
            annotation = '/' + annotation
        try:
            annotation = annotation[0:-1]
            if '(By similarity)' in annotation: annotation,null = string.split(annotation,'(By similarity)')
            if '(Potential)' in annotation: annotation,null = string.split(annotation,'(Potential)')
            if 'By similarity' in annotation: annotation,null = string.split(annotation,'By similarity')
            if 'Potential' in annotation: annotation,null = string.split(annotation,'Potential')
            try:
                if ' ' == annotation[-1]:  annotation = annotation[0:-1]
            except IndexError: annotation = annotation
            if '.' in annotation: annotation,null = string.split(annotation,'.')
            pos1 = int(pos1); pos2 = int(pos2)
            pos1 = pos1-1
            ft_length = pos2-pos1
            if ft_length > 6: pos_1 = pos1; pos_2 = pos2
            else:
                if ft_length < 3: pos_1 = pos1 - 3; pos_2 = pos2 + 3
                else: pos_1 = pos1 - 1; pos_2 = pos2 + 1
                
            if primary_uniprot_id in uniprot_protein_seq_db:
                full_prot_seq = uniprot_protein_seq_db[primary_uniprot_id].Sequence()
                sequence_fragment = full_prot_seq[pos_1:pos_2] ###We will search for this sequence, so have this expanded if too small (see above code)
                if ft != 'CHAIN' and ft != 'CONFLICT' and ft != 'VARIANT' and ft != 'VARSPLIC' and ft != 'VAR_SEQ' and '>' not in annotation: ###exlcludes variant, splice variant SNP and conflict info
                    ft_info = ProteinFunctionalSeqData(primary_uniprot_id,ft,annotation,pos1,pos2,sequence_fragment)
                    ###Store the primary ID as the arrayid (gene accession number)
                    if primary_uniprot_id in uniprot_arrayid_db:
                        arrayids = uniprot_arrayid_db[primary_uniprot_id]
                        for arrayid in arrayids:
                            try: uniprot_ft_db[arrayid].append(ft_info)
                            except KeyError: uniprot_ft_db[arrayid] = [ft_info]
            else:
                ###Occurs for non-SwissProt ft_data (e.g. TrEMBL)
                continue
        except ValueError: continue
        
    domain_gene_count_temp={}
    for geneid in uniprot_ft_db:
        for ft_info in uniprot_ft_db[geneid]:
            try: domain_gene_count_temp[ft_info.PrimaryAnnot(),ft_info.SecondaryAnnot()].append(geneid)
            except KeyError: domain_gene_count_temp[ft_info.PrimaryAnnot(),ft_info.SecondaryAnnot()] = [geneid]
            
    domain_gene_count_temp = eliminateRedundant(domain_gene_count_temp)
    
    for (primary,secondary) in domain_gene_count_temp:
        if len(secondary)>0: key = primary+'-'+secondary
        else: key = primary
        domain_gene_counts[key] = len(domain_gene_count_temp[(primary,secondary)])
    
    print "Number of species uniprot entries imported", len(uniprot_protein_seq_db)
    print "Number of species feature containing entries imported", len(uniprot_ft_db)
    return uniprot_protein_seq_db,uniprot_ft_db,domain_gene_counts

def customDeepCopy(db):
    db2={}
    for i in db:
        try: ###occurs when the contents of the dictionary are an item versus a list
            for e in db[i]:
                try: db2[i].append(e)
                except KeyError: db2[i]=[e]
        except TypeError:
            db2[i] = db[i]
    return db2

def eliminateRedundant(database):
    for key in database:
        try:
            list = makeUnique(database[key])
            list.sort()
        except Exception: list = unique.unique(database[key])
        database[key] = list
    return database

def makeUnique(item):
    db1={}; list1=[]
    for i in item: db1[i]=[]
    for i in db1: list1.append(i)
    list1.sort()
    return list1

def grab_exon_level_feature_calls(species,array_type,genes_analyzed):
    arrayid_uniprot_file = 'AltDatabase/uniprot/'+species+'/'+'arrayid-uniprot.txt'    
    arrayid_ensembl_file = 'AltDatabase/ensembl/'+species+'/'+array_type+'/'+array_type+'-Ensembl.txt'
    ensembl_ft_file = 'AltDatabase/ensembl/'+species+'/'+'DomainFile_All.txt'
    uniprot_feature_file = 'AltDatabase/uniprot/'+species+'/'+'uniprot_feature_file.txt'

    global uniprot_arrayid_db; uniprot_arrayid_db = {}; global arrayid_uniprot_db; arrayid_uniprot_db = {}
    global ensembl_arrayid_db; ensembl_arrayid_db={}
    if array_type != 'exon':
        import_arrayid_uniprot(arrayid_uniprot_file)
        import_arrayid_ensembl(arrayid_ensembl_file)
        ###Otherwise, these databases can be built on-the-fly in downstream methods, since Ensembl will be used as the array gene id
    else: ensembl_arrayid_db = genes_analyzed ###ensembl to ensembl for those being analyzed in the program
    ensembl_protein_seq_db,ensembl_ft_db,domain_gene_counts = import_ensembl_ft_data(species,ensembl_ft_file,ensembl_arrayid_db,array_type) ###Import function domain annotations for Ensembl proteins
    uniprot_protein_seq_db,uniprot_ft_db,domain_gene_counts = import_uniprot_ft_data(species,uniprot_feature_file,domain_gene_counts,ensembl_arrayid_db,array_type)  ###" " " " UniProt "

    arrayid_ft_db = combineDatabases(uniprot_ft_db,ensembl_ft_db)  ###arrayid relating to classes of functional domain attributes and associated proteins (ensembl and uniprot)
    return arrayid_ft_db,domain_gene_counts

def combineDatabases(x,y):
    db1 = customDeepCopy(x); db2 = customDeepCopy(y); db3={}
    for entry in db1: db3[entry] = db1[entry]
    for entry in db2:
        if entry in db3: db3[entry]+=db2[entry]
        else: db3[entry]=db2[entry]
    return db3

def clearall():
    all = [var for var in globals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all: del globals()[var]

if __name__ == '__main__':
    species = 'Mm'; array_type = 'AltMouse'
    findDomainsByGenomeCoordinates(species,array_type); sys.exit()
    
    matchEnsemblDomainCoordinates(protein_feature_file,species,array_type,protein_probeset_db,ens_protein_gene_db,gene_probeset_db)
    kill
    protein_relationship_file,protein_feature_file,protein_seq_fasta = getEnsemblRelationshipDirs(species)
    ens_transcript_protein_db = importEnsemblRelationships(protein_relationship_file,'transcript') ### From Perl script to Ensembl API
    ens_protein_gene_db = importEnsemblRelationships(protein_relationship_file,'gene') ### From Perl script to Ensembl API
    exon_protein_db = importEnsExonStructureDataSimple(species,ens_transcript_protein_db) ### From BioMart
    #kill
    if array_type == 'exon': ens_probeset_file = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_probesets.txt"    
    else: ens_probeset_file = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_"+array_type+"_probesets.txt"
    protein_probeset_db,gene_probeset_db = importSplicingAnnotationDatabase(ens_probeset_file,exon_protein_db) ### Derived from ExonArrayEnsemblRules
    #matchEnsemblDomainCoordinates(protein_feature_file,species,array_type,protein_probeset_db,ens_protein_gene_db,gene_probeset_db)
    kill
