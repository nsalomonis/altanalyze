###ExonAnalyze_module
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

import math
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os.path
import unique
from stats_scripts import statistics
import update

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list
            
def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def grabJunctionData(species,array_type,key_type,root_dir):
    if array_type == 'AltMouse': filename = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'_junction-comparisons.txt'
    else:
        filename = 'AltDatabase/' + species + '/'+array_type+'/'+ species + '_junction_comps_updated.txt'
        if array_type == 'RNASeq': filename = root_dir+filename ### This is a dataset specific file
    fn=filepath(filename); critical_exon_junction_db = {}; x=0
    for line in open(fn,'rU').xreadlines():
        if x==0: x=1
        else:
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if array_type == 'AltMouse':
                geneid,probeset1,probeset2,critical_exons = t
                probesets = [probeset1, probeset2]; critical_exons = string.split(critical_exons,'|')
            else:
                geneid,critical_exon,excl_junction,incl_junction,excl_junction_probeset,incl_junction_probeset,source = t
                probeset1 = incl_junction_probeset; probeset2 = excl_junction_probeset
                probesets = [incl_junction_probeset, excl_junction_probeset]; critical_exons = [critical_exon]
            for exon in critical_exons:
                if key_type == 'gene-exon':
                    key = geneid+':'+exon
                    ###Record for each probeset what critical junctions it is associated with
                    try: critical_exon_junction_db[key].append((probeset1, probeset2))
                    except KeyError: critical_exon_junction_db[key] = [(probeset1, probeset2)]
                    #print key,(probeset1, probeset2)
                else: critical_exon_junction_db[(probeset1, probeset2)] = critical_exons
    return critical_exon_junction_db

class GeneAnnotationData:
    def __init__(self, geneid, description, symbol, external_geneid, rna_processing_annot):
        self._geneid = geneid; self._description = description; self._symbol = symbol
        self._rna_processing_annot = rna_processing_annot; self._external_geneid = external_geneid
    def GeneID(self): return self._geneid
    def Description(self): return self._description
    def Symbol(self): return self._symbol
    def ExternalGeneID(self): return self._external_geneid
    def TranscriptClusterIDs(self):
        if self.GeneID() in gene_transcript_cluster_db:
            transcript_cluster_list = gene_transcript_cluster_db[self.GeneID()]
            return transcript_cluster_list
        else:
            try: return transcript_cluster_list
            except UnboundLocalError: return [''] ### Occurs for SplicingIndex method with AltMouse data
    def RNAProcessing(self): return self._rna_processing_annot
    def Report(self):
        output = self.GeneID() +'|'+ self.Symbol() +'|'+ self.Description()
        return output
    def __repr__(self): return self.Report()
    
def import_annotations(filename,array_type,keyBySymbol=False):
    fn=filepath(filename); annotate_db = {}; x = 0
    if array_type == 'AltMouse':
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            if x == 0: x = 1
            else:
                try: affygene, description, ll_id, symbol, rna_processing_annot = string.split(data,'\t')
                except ValueError: affygene, description, ll_id, symbol = string.split(data,'\t'); rna_processing_annot = ''
                if '"' in description: null,description,null = string.split(description,'"')
                y = GeneAnnotationData(affygene, description, symbol, ll_id, rna_processing_annot)
                if keyBySymbol:
                    annotate_db[symbol] = y
                else:
                    annotate_db[affygene] = y
    else:
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            rna_processing_annot=''
            try: ensembl, symbol, description, rna_processing_annot = string.split(data,'\t')
            except ValueError: ensembl, description, symbol = string.split(data,'\t')
            y = GeneAnnotationData(ensembl, description, symbol, ensembl, rna_processing_annot)
            annotate_db[ensembl] = y
            if keyBySymbol:
                annotate_db[symbol] = y
            else:
                annotate_db[ensembl] = y
    return annotate_db

def importmicroRNADataExon(species,array_type,exon_db,microRNA_prediction_method,explicit_data_type,root_dir):
    filename = "AltDatabase/"+species+"/"+array_type+"/"+species+"_probeset_microRNAs_"+microRNA_prediction_method+".txt"
    if array_type == 'junction': filename = string.replace(filename,'.txt','-filtered.txt')

    fn=filepath(filename); microRNA_full_exon_db={}; microRNA_count_db={}; gene_microRNA_denom={}
    if array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null'):
        critical_exon_junction_db = grabJunctionData(species,array_type,'gene-exon',root_dir)
        if array_type == 'AltMouse':
            gene_annotation_file = "AltDatabase/"+species+"/"+array_type+"/"+array_type+"_gene_annotations.txt"
            annotate_db = import_annotations(gene_annotation_file,array_type)

    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        probeset=t[0];microRNA=t[1]
        try: mir_seq=t[2];mir_sources=t[3]
        except IndexError: mir_seq='';mir_sources=''
        try:
            if array_type == 'exon' or array_type == 'gene' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type != 'null'):
                ed = exon_db[probeset]; probeset_list = [probeset]; exonid = ed.ExonID(); symbol = ed.Symbol(); geneid = ed.GeneID()
            else:
                probeset_list = []; uid = probeset ###This is the gene with the critical exon it's mapped to
                for junctions in critical_exon_junction_db[uid]:
                    if junctions in exon_db:
                        geneid = exon_db[junctions].GeneID()
                        if array_type == 'AltMouse':
                            if geneid in annotate_db: symbol = annotate_db[geneid].Symbol()
                            else: symbol = geneid
                        else:
                            try: symbol = exon_db[junctions].Symbol()
                            except Exception: symbol=''
                    else: geneid = ''
                    #if junctions in exon_db:
                    probeset_list.append(junctions)
            if len(geneid)>0:
                microRNA_info = microRNA, symbol, mir_seq, mir_sources
                for probeset in probeset_list:
                    try: microRNA_full_exon_db[geneid,probeset].append(microRNA_info)
                    except KeyError: microRNA_full_exon_db[geneid,probeset] = [microRNA_info]
                    
                    try: microRNA_count_db[microRNA].append(geneid)
                    except KeyError: microRNA_count_db[microRNA] = [geneid]
                    #if 'ENS' in microRNA: print [data,t];kill
                    gene_microRNA_denom[geneid] = []
        except KeyError: null=[]
    microRNA_count_db = eliminate_redundant_dict_values(microRNA_count_db)

    ###Replace the actual genes with the unique gene count per microRNA
    for microRNA in microRNA_count_db: microRNA_count_db[microRNA] = len(microRNA_count_db[microRNA])
    if array_type == 'RNASeq': id_name = 'junction IDs'
    else: id_name = 'array IDs'
    print len(gene_microRNA_denom),"genes with a predicted microRNA binding site aligning to a",id_name
    return microRNA_full_exon_db,microRNA_count_db,gene_microRNA_denom
        
def filterMicroRNAProbesetAssociations(microRNA_full_exon_db,exon_hits):
    filtered_microRNA_exon_db = {}
    microRNA_gene_count_db = {}
    for key in microRNA_full_exon_db:
        array_geneid = key[0]
        if key in exon_hits:
            filtered_microRNA_exon_db[key] = microRNA_full_exon_db[key]
            for (microRNA,gene_symbol,seq,source) in microRNA_full_exon_db[key]:
                gene_info = array_geneid, gene_symbol
                try: microRNA_gene_count_db[microRNA].append(gene_info)
                except KeyError: microRNA_gene_count_db[microRNA] = [gene_info]
    return filtered_microRNA_exon_db
                        
class ExonSequenceData:
    def __init__(self, gene_id, exon_id, exon_seq, complete_mRNA_length,strand):
        self._gene_id = gene_id; self._exon_id = exon_id; self._strand = strand
        self._exon_seq = exon_seq; self._complete_mRNA_length = complete_mRNA_length
    def GeneID(self): return self._gene_id
    def ExonID(self): return self._exon_id
    def ExonSeq(self): return self._exon_seq
    def RNALen(self): return int(self._complete_mRNA_length)
    def Strand(self): return self._strand
    def SummaryValues(self):
        output = self.GeneID()+'|'+self.SecondaryID()+'|'+self.ExonID()
        return output
    def __repr__(self): return self.SummaryValues()

class ExonProteinAlignmentData:
    def __init__(self,geneid,probeset,exonid,protein_hit_id,protein_null_id):
        self._gene_id = geneid; self._exon_id = exonid; self._probeset = probeset
        self._protein_hit_id = protein_hit_id; self._protein_null_id = protein_null_id
    def GeneID(self): return self._gene_id
    def ExonID(self): return self._exon_id
    def Probeset(self): return self._probeset
    def HitProteinID(self): return self._protein_hit_id  
    def NullProteinID(self): return self._protein_null_id
    def RecordHitProteinExonIDs(self,hit_exon_ids): self._hit_exon_ids = hit_exon_ids
    def RecordNullProteinExonIDs(self,null_exon_ids): self._null_exon_ids = null_exon_ids
    def HitProteinExonIDs(self):
        exon_list_str = string.join(self._hit_exon_ids,',')
        return exon_list_str
    def NullProteinExonIDs(self):
        exon_list_str = string.join(self._null_exon_ids,',')
        return exon_list_str
    def SummaryValues(self):
        output = self.GeneID()+'|'+self.ExonID()+'|'+self.Probeset()
        return output
    def __repr__(self): return self.SummaryValues()

def importExonSequenceBuild(filename,exon_db):
    ###Parse protein sequence
    fn=filepath(filename); protein_sequence_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        protein_id,protein_seq = string.split(data,'\t')
        protein_sequence_db[protein_id] = protein_seq
    
    ###Parse protein/probeset data (built in FeatureAlignment)        
    filename = string.replace(filename,'SEQUENCE','probeset')
    fn=filepath(filename); probeset_protein_db = {}; positive_protein_associations = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line); t = string.split(data,'\t')
        try: probeset,protein_hit_id,protein_null_id = t
        except ValueError: gene,probeset,exonid,protein_hit_id,protein_null_id = t
        try:
            ed = exon_db[probeset] ###update this information every time you run (in case a new exon database is used)
            ep = ExonProteinAlignmentData(ed.GeneID(),probeset,ed.ExonID(),protein_hit_id,protein_null_id)
            probeset_protein_db[probeset] = ep
            try: positive_protein_associations[protein_hit_id].append((probeset,ed.ExonID())) ### Record which probesets actually match to which proteins (e.g. informs you which align to which nulls too)
            except KeyError: positive_protein_associations[protein_hit_id] = [(probeset,ed.ExonID())]
        except KeyError: null=[]
    ###Determine for each null and hit proteins, what exons are they typically associated with (probably not too informative but record anyways)
    for probeset in probeset_protein_db:
        ep = probeset_protein_db[probeset]
        try:
            null_tuple_exonid_list = positive_protein_associations[ep.NullProteinID()]
            null_exonid_list = formatExonLists(null_tuple_exonid_list)
        except KeyError: null_exonid_list = []
        hit_tuple_exonid_list = positive_protein_associations[ep.HitProteinID()]
        hit_exonid_list = formatExonLists(hit_tuple_exonid_list)
        ep.RecordHitProteinExonIDs(hit_exonid_list)
        ep.RecordNullProteinExonIDs(null_exonid_list)
        #a = ep.HitProteinExonIDs()
        #b = ep.NullProteinExonIDs()
        #print ep.HitProteinExonIDs()
        #print ep.NullProteinExonIDs();kill
    return probeset_protein_db,protein_sequence_db

def formatExonLists(exonid_tuple_lists):
    exonid_list=[]; exonid_tuple_lists.sort()
    for (probeset,exonid) in exonid_tuple_lists: exonid_list.append(exonid)
    return exonid_list

def import_existing_sequence_build(filename):
    fn=filepath(filename); transcript_cdna_sequence_dbase = {}; exon_sequence_database = {}; transcript_associations = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        try: array_geneid, altmerge, strand, dna_exon_seq, exon_structure, coding_seq, peptide_length, reference_seq_name,reference_seq, ref_AC, full_length, start_exon, stop_exon, null,null,null,null,refseq_mRNA,null,null,null,null,null,null,null,nmd_call = string.split(data,'\t')       
        except ValueError: t = string.split(data,'\t');print len(t);kill
        if array_geneid == 'array_geneid': title = data
        else:
            dna_exon_seq = dna_exon_seq[0:-3]  #there is a "***" at the end of each transcript sequence
            dna_exon_seq = string.split(dna_exon_seq,'('); dna_exon_seq = dna_exon_seq[1:] #get rid of the first blank generated by parsing
            ### separate the exon number and exon sequence into a tuple and store in order within 'exon_seq_tuple_list'
            mRNA_length = 0; exon_seq_tuple_list = []
            for exon_seq in dna_exon_seq:   #exon_seq : E20)tccccagctttgggtggtgg
                exon_num, dna_seq = string.split(exon_seq,')'); mRNA_length += len(dna_seq)
                if mRNA_length > 18:
                    esd = ExonSequenceData(array_geneid,exon_num,dna_seq,mRNA_length,strand)
                    exon_sequence_database[array_geneid,exon_num] = esd
            transcript_data = [int(peptide_length), altmerge,[start_exon, stop_exon], exon_structure, nmd_call, full_length, coding_seq, strand, reference_seq,ref_AC,float(mRNA_length),refseq_mRNA]
            try: transcript_cdna_sequence_dbase[array_geneid].append(transcript_data)
            except KeyError: transcript_cdna_sequence_dbase[array_geneid] = [transcript_data]  
            transcript_associations[array_geneid,exon_structure] = altmerge
    print "\nNumber of exon sequences imported",len(exon_sequence_database)
    print "Number of transcript sequences imported: ",len(transcript_cdna_sequence_dbase)
    return transcript_cdna_sequence_dbase, transcript_associations, exon_sequence_database

def exportAltMouseExonSequence():
    probeset_exon_db={}; x=0
    species = 'Mm'; array_type = 'AltMouse'

    critical_exon_import_file = 'AltDatabase/Mm/AltMouse/AltMouse_junction-comparisons.txt'
    update.verifyFile(critical_exon_import_file,array_type)
    critical_exon_db={}; critical_probesets={}
    fn=filepath(critical_exon_import_file)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        gene,probeset1,probeset2,critical_exons=string.split(data,'\t')
        critical_exons= string.split(critical_exons,'|')
        for exon in critical_exons:
            try: critical_exon_db[gene,exon].append(probeset1+'-'+probeset2)
            except KeyError: critical_exon_db[gene,exon] = [probeset1+'-'+probeset2]
            critical_probesets[probeset1]=[]; critical_probesets[probeset2]=[]
            
    probeset_annotations_file = "AltDatabase/Mm/AltMouse/MASTER-probeset-transcript.txt"
    update.verifyFile(probeset_annotations_file,array_type)
    fn=filepath(probeset_annotations_file)
    for line in open(fn,'rU').xreadlines():             
        probeset_data = cleanUpLine(line)  #remove endline
        if x==0: x=1
        else:
            probeset,affygene,exons,transcript_num,transcripts,probe_type_call,ensembl,block_exon_ids,block_structure,comparison_info = string.split(probeset_data,'\t')
            if probeset in critical_probesets:
                exons = exons[:-1]; exons = string.split(exons,'-')
                affygene = affygene[:-1]
                if '|' in exons: print exons;kill
                probeset_exon_db[probeset,affygene]=exons

    exon_protein_sequence_file = "AltDatabase/Mm/AltMouse/SEQUENCE-transcript-dbase.txt"
    update.verifyFile(exon_protein_sequence_file,array_type)
    transcript_cdna_sequence_dbase,transcript_associations,exon_sequence_database = import_existing_sequence_build(exon_protein_sequence_file)
    
    critical_exon_seq_export = 'AltDatabase/Mm/AltMouse/AltMouse_critical-exon-seq.txt'
    update.verifyFile(critical_exon_seq_export,array_type)
    fn=filepath(critical_exon_seq_export)
    data = open(fn,'w')
    title = ['Affygene:exon','critical_exon-num','critical-probeset-comps']; title = string.join(title,'\t')+'\n'; data.write(title)    
    for (gene,exon_num) in critical_exon_db:
        probeset_comp_list = critical_exon_db[(gene,exon_num)]; probeset_comp_list = string.join(probeset_comp_list,'|')
        try: ###Restrict export to previously exported critical exons (ExonAnnotate_module)
            exon_sequence_database[(gene,exon_num)]; esd = exon_sequence_database[(gene,exon_num)]
            exon_seq = esd.ExonSeq()
            exon_data = string.join([gene+':'+exon_num,probeset_comp_list,exon_seq],'\t')+'\n'
            data.write(exon_data)
        except KeyError: null=[]
    data.close()

    probeset_seq_file = 'AltDatabase/Mm/AltMouse/probeset_sequence_reversed.txt'
    update.verifyFile(probeset_seq_file,array_type)
    probeset_seq_db={}; x=0
    fn=filepath(probeset_seq_file)
    for line in open(fn,'rU').xreadlines():
        if x == 0: x=1
        else:
            data = cleanUpLine(line); t = string.split(data,'\t')
            probeset = t[0]
            probeset_seq_list = t[1:]
            probeset_seq_db[probeset] = probeset_seq_list
            
    critical_junction_seq_export = 'AltDatabase/Mm/AltMouse/AltMouse_critical-junction-seq.txt'
    update.verifyFile(critical_junction_seq_export,array_type)
    fn=filepath(critical_junction_seq_export)
    data = open(fn,'w'); x=0; k=0;l=0
    title = ['probeset','probeset-seq','junction-seq']; title = string.join(title,'\t')+'\n'; data.write(title)
    for (probeset,gene) in probeset_exon_db:
        junction_seq = []; y=0; positions=[]
        try:
            probeset_seq_list = probeset_seq_db[probeset]
            for exon_num in probeset_exon_db[(probeset,gene)]: 
                try: ###Restrict export to previously exported critical exons (ExonAnnotate_module)
                    exon_sequence_database[(gene,exon_num)]; esd = exon_sequence_database[(gene,exon_num)]
                    exon_seq = esd.ExonSeq(); strand = esd.Strand()
                    junction_seq.append(exon_seq); y+=1
                    #exon_data = string.join([gene+':'+exon_num,probeset_comp_list,exon_seq],'\t')+'\n'
                    #data.write(exon_data)
                except KeyError: null=[]
            #if 'E5' in probeset_exon_db[(probeset,gene)]:
            if y>0:
                if strand == '-': junction_seq.reverse()
                junction_seq_str = string.join(junction_seq,'')
                junction_seq_str = string.upper(junction_seq_str)
                not_found = 0
                for probeset_seq in probeset_seq_list:
                    #probeset_seq = reverse_string(probeset_seq)
                    probeset_seq_rev = reverse_orientation(probeset_seq)
                    if probeset_seq in junction_seq_str:
                        f = string.find(junction_seq_str,probeset_seq)
                        positions.append((f,len(probeset_seq)))
                        k+=1
                    else:
                        not_found = 1
                        x+=1
                if not_found == 1:
                    new_probeset_seq = probeset_seq_list[0] ###pick the first probe sequence found
                if len(positions)>0:
                    positions.sort()
                    new_probeset_seq = junction_seq_str[positions[0][0]:positions[-1][0]+positions[-1][1]]
                    #print new_probeset_seq,positions, probeset,probeset_exon_db[(probeset,gene)],probeset_seq_list,junction_seq;kill
                junction_seq = string.join(junction_seq,'|') ###indicate where the junction is
                probe_seq_data = string.join([probeset,new_probeset_seq,junction_seq],'\t')+'\n'
                data.write(probe_seq_data)
        except KeyError: null=[]
    data.close()
    print k,x

def reverse_orientation(sequence):
    """reverse the orientation of a sequence (opposite strand)"""
    exchange = []
    for nucleotide in sequence:
        if nucleotide == 'A': nucleotide = 'T'
        elif nucleotide == 'T': nucleotide = 'A'
        elif nucleotide == 'G': nucleotide = 'C'
        elif nucleotide == 'C': nucleotide = 'G'
        exchange.append(nucleotide)
    complementary_sequence = reverse_string(exchange)
    return complementary_sequence

def reverse_string(astring):
    "http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/65225"
    revchars = list(astring)        # string -> list of chars
    revchars.reverse()              # inplace reverse the list
    revchars = ''.join(revchars)    # list of strings -> string
    return revchars

def compareProteinFeatures(protein_ft,neg_coding_seq,pos_coding_seq):
    ###Parse out ft-information. Generate ft-fragment sequences for querying
    ###This is a modification of the original script from FeatureAlignment but simplified for exon analysis
    protein_ft_unique=[]; new_ft_list = []
    for ft_data in protein_ft:
        ft_name = ft_data.PrimaryAnnot(); domain_seq = ft_data.DomainSeq(); annotation = ft_data.SecondaryAnnot()
        protein_ft_unique.append((ft_name,annotation,domain_seq))
    ###Redundant entries that are class objects can't be eliminated, so save to a new list and eliminate redundant entries
    protein_ft_unique = unique.unique(protein_ft_unique)
    for (ft_name,annotation,domain_seq) in protein_ft_unique:
        ft_length = len(domain_seq)
        new_ft_data = 'null',domain_seq,ft_name,annotation
        new_ft_list.append(new_ft_data)
    new_ft_list = unique.unique(new_ft_list)
    pos_ft = []; neg_ft = []; all_fts = []
    for (pos,seq,ft_name,annot) in new_ft_list:
        if seq in pos_coding_seq:
            pos_ft.append([pos,seq,ft_name,annot]); all_fts.append([pos,seq,ft_name,annot])
        if seq in neg_coding_seq:
            neg_ft.append([pos,seq,ft_name,annot]); all_fts.append([pos,seq,ft_name,annot])
    all_fts = unique.unique(all_fts)
    pos_ft_missing=[]; neg_ft_missing=[]
    for entry in all_fts:
        if entry not in pos_ft: pos_ft_missing.append(entry)
        if entry not in neg_ft: neg_ft_missing.append(entry)
    pos_ft_missing2=[]; neg_ft_missing2=[]
    for entry in pos_ft_missing: entry[1] = ''; pos_ft_missing2.append(entry)
    for entry in neg_ft_missing: entry[1] = ''; neg_ft_missing2.append(entry)
        
    pos_ft_missing2 = unique.unique(pos_ft_missing2)
    neg_ft_missing2 = unique.unique(neg_ft_missing2)  
    return neg_ft_missing2,pos_ft_missing2

def getFeatureIsoformGenomePositions(species,protein_ft_db,mRNA_protein_seq_db,gene_transcript_db,coordinate_type):
    """ Adapted from compareProteinFeatures but for one isoform and returns genomic coordinates for each feature
    This function is designed to export all unique isoforms rather than just comparison isoforms """
    
    import export
    export_file = 'AltDatabase/ensembl/'+species+'/ProteinFeatureIsoform_complete.txt'                
    export_data = export.ExportFile(export_file)

    failed = 0
    worked = 0
    failed_ac=[]
    for gene in protein_ft_db:
        transcript_feature_db={}
        for ft in protein_ft_db[gene]:
            try:
                ft_name = ft.PrimaryAnnot(); annotation = ft.SecondaryAnnot()
                for (mRNA,type) in gene_transcript_db[gene]:
                    try:
                        protein,protein_seq = mRNA_protein_seq_db[mRNA]
                        error = False
                    except Exception:
                        failed_ac.append(mRNA)
                        error = True
                    if error == False:
                        if ft.DomainSeq() in protein_seq:
                            #if coordinate_type == 'genomic':
                            pos1_genomic = ft.GenomicStart(); pos2_genomic = ft.GenomicStop()
                            #else:
                            pos1 = str(ft.DomainStart()); pos2 = str(ft.DomainEnd())
    
                            ### There are often many features that overlap within a transcript, so consistently pick just one
                            if mRNA in transcript_feature_db:
                                db = transcript_feature_db[mRNA]
                                if (pos1,pos2) in db:
                                    db[pos1, pos2].append([pos1_genomic, pos2_genomic, protein,ft_name,annotation])
                                else:
                                    db[pos1, pos2]=[[pos1_genomic, pos2_genomic, protein,ft_name,annotation]]
                            else:
                                db={}
                                db[pos1, pos2]=[[pos1_genomic, pos2_genomic, protein,ft_name,annotation]]
                                transcript_feature_db[mRNA] = db
                                
                            #values = [mRNA, protein, pos1, pos2,ft_name,annotation]; unique_entries.append(values)
                            worked+=1
            except IOError:
                failed+=1

        for transcript in transcript_feature_db:
            db = transcript_feature_db[transcript]
            for (pos1,pos2) in db:
                db[pos1,pos2].sort() ### Pick the alphabetically listed first feature
                pos1_genomic, pos2_genomic, protein,ft_name,annotation = db[pos1,pos2][0]
                values = [transcript, protein, pos1, pos2,pos1_genomic, pos2_genomic, ft_name,annotation]
                export_data.write(string.join(values,'\t')+'\n')
                
    export_data.close()
    print failed,'features failed to have corresponding aligned genomic locations out of', worked+failed
    failed_ac = unique.unique(failed_ac)
    print len(failed_ac),'mRNAs without identified/in silico derived proteins'  ### Appear to be ncRNAs without ATGs
    print failed_ac[:20]
    
def identifyAltIsoformsProteinComp(probeset_gene_db,species,array_type,protein_domain_db,compare_all_features,data_type):
    """ This function is used by the module IdentifyAltIsoforms to run 'characterizeProteinLevelExonChanges'"""
    global protein_ft_db; protein_ft_db = protein_domain_db; protein_domain_db=[]
    exon_db={} ### Create a simplified version of the exon_db dictionary with probesets that map to a match and null protein
    for probeset in probeset_gene_db:
        gene, exon_id = probeset_gene_db[probeset]
        ep = ExonProteinAlignmentData(gene,probeset,exon_id,'',''); exon_db[probeset] = ep
    global protein_sequence_db
    if compare_all_features == 'yes': type = 'seqcomp'
    else: type = 'exoncomp'
    if (array_type == 'junction' or array_type == 'RNASeq') and data_type != 'null':
        exon_protein_sequence_file = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/'+'SEQUENCE-protein-dbase_'+type+'.txt'
    else:
        exon_protein_sequence_file = 'AltDatabase/'+species+'/'+array_type+'/'+'SEQUENCE-protein-dbase_'+type+'.txt'
    probeset_protein_db,protein_sequence_db = importExonSequenceBuild(exon_protein_sequence_file,exon_db)
    
    exon_hits={}
    for probeset in probeset_protein_db:
        gene = probeset_protein_db[probeset].GeneID()
        exon_hits[gene,probeset]=[]

    include_sequences = 'no' ### Sequences for comparisons are unnecessary to store. List array-type as exon since AltMouse data has been re-organized, later get rid of AltMouse specific functionality in this function
    functional_attribute_db,protein_features = characterizeProteinLevelExonChanges(species,exon_hits,probeset_protein_db,'exon',include_sequences)

    if (array_type == 'junction' or array_type == 'RNASeq') and data_type != 'null':
        export_file = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/probeset-domain-annotations-'+type+'.txt' 
    else:
        export_file = 'AltDatabase/'+species+'/'+array_type+'/probeset-domain-annotations-'+type+'.txt' 
    formatAttributeForExport(protein_features,export_file)

    if (array_type == 'junction' or array_type == 'RNASeq') and data_type != 'null':
        export_file = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/probeset-protein-annotations-'+type+'.txt' 
    else:
        export_file = 'AltDatabase/'+species+'/'+array_type+'/probeset-protein-annotations-'+type+'.txt' 
    formatAttributeForExport(functional_attribute_db,export_file)

def formatAttributeForExport(attribute_db,filename):
    from build_scripts import IdentifyAltIsoforms
    export_db={}
    for (gene,probeset) in attribute_db:
        attribute_list = attribute_db[(gene,probeset)]; attribute_list2=[]
        for (attribute,direction) in attribute_list:
            attribute = string.replace(attribute,'|',' ')
            attribute_list2.append(attribute+'|'+direction)
        export_db[probeset]=attribute_list2
    print 'Exporting:',filename
    IdentifyAltIsoforms.exportSimple(export_db,filename,'')
        
def importTranscriptBiotypeAnnotations(species):
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-biotypes.txt'
    accepted_biotypes = ['nonsense_mediated_decay','non_coding','retained_intron','retrotransposed']
    fn=filepath(filename); biotype_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        gene,transcript,biotype = string.split(data,'\t')
        if biotype in accepted_biotypes:
            biotype_db[transcript] = biotype
    return biotype_db

def characterizeProteinLevelExonChanges(species,exon_hits,probeset_protein_db,array_type,include_sequences):
    """Examines the the two reciprocal isoforms for overal changes in protein sequence and general sequence composition"""
    
    try: biotype_db = importTranscriptBiotypeAnnotations(species)
    except Exception: biotype_db={}
    
    functional_attribute_db={}; protein_features={}
    for (array_geneid,uid) in exon_hits: ###uid is probeset or (probeset1,probeset2) value, depending on array_type
            if array_type != 'exon' and array_type != 'gene': probeset,probeset2 = uid
            else: probeset = uid
            hv=0
            try:
                pp = probeset_protein_db[probeset]; hv=1
                pos_ref_AC = pp.HitProteinID()
                neg_ref_AC = pp.NullProteinID()
                if array_type != 'exon' and array_type != 'gene': ### Instead of using the null_hit, use the second junction probeset
                    try: np = probeset_protein_db[probeset2]; neg_ref_AC = np.HitProteinID()
                    except KeyError: null =[] ###just use the existing null                
            except KeyError: ###occurs if there is not protein associated with the first probeset (NA for exon arrays)
                if array_type != 'exon' and array_type != 'gene': ###Reverse of above
                    try:
                        np = probeset_protein_db[probeset2]; hv=2
                        pos_ref_AC = np.NullProteinID()
                        neg_ref_AC = np.HitProteinID()
                    except KeyError: hv=0
            if hv!=0:
                neg_coding_seq = protein_sequence_db[neg_ref_AC]
                pos_coding_seq = protein_sequence_db[pos_ref_AC]
                
                pos_biotype=None; neg_biotype=None
                if pos_ref_AC in biotype_db: pos_biotype = biotype_db[pos_ref_AC]
                if neg_ref_AC in biotype_db: neg_biotype = biotype_db[neg_ref_AC]
                
                neg_length = len(neg_coding_seq)
                pos_length = len(pos_coding_seq)
                pos_length = float(pos_length);  neg_length = float(neg_length)
                if array_geneid in protein_ft_db:
                        protein_ft = protein_ft_db[array_geneid]
                        neg_ft_missing,pos_ft_missing = compareProteinFeatures(protein_ft,neg_coding_seq,pos_coding_seq)
                        for (pos,blank,ft_name,annotation) in pos_ft_missing:
                            call_x = '-'
                            if len(annotation)>0:
                                data_tuple = ft_name,call_x 
                                data_tuple = ft_name +'-'+annotation,call_x
                            else: data_tuple = ft_name,call_x
                            try: protein_features[array_geneid,uid].append(data_tuple)
                            except KeyError: protein_features[array_geneid,uid] = [data_tuple]
                        for (pos,blank,ft_name,annotation) in neg_ft_missing:
                            ###If missing from the negative list, it is present in the positive state
                            call_x = '+'
                            if len(annotation)>0:
                                data_tuple = ft_name,call_x 
                                data_tuple = ft_name +'-'+annotation,call_x
                            else: data_tuple = ft_name,call_x
                            try: protein_features[array_geneid,uid].append(data_tuple)
                            except KeyError: protein_features[array_geneid,uid] = [data_tuple]
                call=''
                if pos_biotype != None or neg_biotype !=None:
                    ### For example, if one or both transcript(s) are annotated as nonsense_mediated_decay
                    if (neg_length - pos_length)>0: call = '-'
                    elif (pos_length - neg_length)>0: call = '+'
                    else: call = '~'
                    if pos_biotype != None:
                        data_tuple = pos_biotype,call
                        try: functional_attribute_db[array_geneid,uid].append(data_tuple)
                        except KeyError: functional_attribute_db[array_geneid,uid]= [data_tuple]
                    if neg_biotype != None:
                        data_tuple = neg_biotype,call
                        try: functional_attribute_db[array_geneid,uid].append(data_tuple)
                        except KeyError: functional_attribute_db[array_geneid,uid]= [data_tuple]
                if pos_coding_seq[:5] != neg_coding_seq[:5]:
                    function_var = 'alt-N-terminus'
                    if (neg_length - pos_length)>0: call = '-'
                    elif (pos_length - neg_length)>0: call = '+'
                    else: call = '~'
                    data_tuple = function_var,call
                    try: functional_attribute_db[array_geneid,uid].append(data_tuple)
                    except KeyError: functional_attribute_db[array_geneid,uid]= [data_tuple]
                if pos_coding_seq[-5:] != neg_coding_seq[-5:]:
                    #print probeset,(1-((neg_length - pos_length)/neg_length)), ((neg_length - pos_length)/neg_length),[neg_length],[pos_length]
                    #print probeset,(1-((pos_length - neg_length)/pos_length)),((pos_length - neg_length)/pos_length)
                    if (1-((neg_length - pos_length)/neg_length)) < 0.5 and ((neg_length - pos_length)/neg_length) > 0 and (pos_coding_seq[:5] == neg_coding_seq[:5]):
                        function_var = 'truncated'
                        call = '+'; data_tuple = function_var,call
                        try: functional_attribute_db[array_geneid,uid].append(data_tuple)
                        except KeyError: functional_attribute_db[array_geneid,uid]=[data_tuple]                        
                    elif (1-((pos_length - neg_length)/pos_length)) < 0.5 and ((pos_length - neg_length)/pos_length) > 0 and (pos_coding_seq[:5] == neg_coding_seq[:5]): 
                        function_var = 'truncated'
                        call = '-'; data_tuple = function_var,call
                        try: functional_attribute_db[array_geneid,uid].append(data_tuple)
                        except KeyError: functional_attribute_db[array_geneid,uid]=[data_tuple]
                    else:
                        function_var = 'alt-C-terminus'
                        if (neg_length - pos_length)>0: call = '-'
                        elif (pos_length - neg_length)>0: call = '+'
                        else: call = '~'
                        data_tuple = function_var,call
                        try: functional_attribute_db[array_geneid,uid].append(data_tuple)
                        except KeyError: functional_attribute_db[array_geneid,uid]=[data_tuple]
                if call == '':
                    if pos_coding_seq != neg_coding_seq:
                        function_var = 'alt-coding'
                        if (neg_length - pos_length)>0: call = '-'
                        elif (pos_length - neg_length)>0: call = '+'
                        else: call = '~'
                        data_tuple = function_var,call
                        try: functional_attribute_db[array_geneid,uid].append(data_tuple)
                        except KeyError: functional_attribute_db[array_geneid,uid]= [data_tuple]
                ### Record change in peptide size
                if neg_length > pos_length: fcall = '-'
                elif neg_length < pos_length: fcall = '+'
                elif neg_length == pos_length: fcall = '~'
                if len(pos_ref_AC)<1: pos_ref_AC = 'NULL'
                if len(neg_ref_AC)<1: neg_ref_AC = 'NULL'
                if fcall == '-':  
                    function_var1 = 'AA:' + str(int(pos_length))+'('+pos_ref_AC+')' +'->'+ str(int(neg_length))+'('+neg_ref_AC+')'
                    if include_sequences == 'yes': function_var2 = 'sequence: ' +'('+pos_ref_AC+')'+pos_coding_seq +' -> '+ '('+neg_ref_AC+')'+neg_coding_seq
                else:
                    function_var1 = 'AA:' +str(int(neg_length))+ '('+neg_ref_AC+')' +'->'+ str(int(pos_length))+'('+pos_ref_AC+')'
                    if include_sequences == 'yes': function_var2 = 'sequence: ' +'('+neg_ref_AC+')'+neg_coding_seq +' -> '+'('+pos_ref_AC+')'+pos_coding_seq
                data_tuple1 = function_var1,fcall
                try: functional_attribute_db[array_geneid,uid].append(data_tuple1)
                except KeyError: functional_attribute_db[array_geneid,uid]= [data_tuple1]
                ### Record sequence change
                if include_sequences == 'yes':
                    data_tuple2 = function_var2,fcall
                    try: functional_attribute_db[array_geneid,uid].append(data_tuple2)
                    except KeyError: functional_attribute_db[array_geneid,uid]= [data_tuple2]
    print len(functional_attribute_db),'Genes with affected functional attributes'
    return functional_attribute_db,protein_features
                                              
def combine_databases(db1,db2):
    for key in db2:
        if key not in db1:
            db1[key] = db2[key]
    return db1

if __name__ == '__main__':
    species = 'Mm'; array_type='AltMouse'
 
    probeset_protein_db = exon_sequence_database
    protein_sequence_db = uniprot_seq_string
    exon_hits={}
    exon_hits[('ENSG00000130939', 'E9-3|')] =  [['E9-3', '2319586']]
    functional_attribute_db,protein_features = characterizeProteinLevelExonChanges(exon_hits,probeset_protein_db)
    sys.exit()

    exon_protein_sequence_file = "AltDatabase/"+species+"/"+array_type+"/"+"SEQUENCE-protein-dbase.txt"
    transcript_cdna_sequence_dbase,uniprot_seq_string = importExonSequenceBuild(exon_protein_sequence_file)
    sys.exit()
