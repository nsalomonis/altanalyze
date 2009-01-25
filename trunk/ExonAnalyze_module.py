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
import sys, string
import os.path
import unique
import statistics

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

def grabJunctionData(species,array_type,key_type):
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'_junction-comparisons.txt'
    fn=filepath(filename); critical_exon_junction_db = {}; x=0
    for line in open(fn,'rU').xreadlines():
        if x==0: x=1
        else:
            data = cleanUpLine(line)
            geneid,probeset1,probeset2,critical_exons = string.split(data,'\t')
            probesets = [probeset1, probeset2]; critical_exons = string.split(critical_exons,'|')
            for exon in critical_exons:
                if key_type == 'gene-exon':
                    key = geneid+':'+exon
                    ###Record for each probeset what critical junctions it is associated with
                    try: critical_exon_junction_db[key].append((probeset1, probeset2))
                    except KeyError: critical_exon_junction_db[key] = [(probeset1, probeset2)]
                    #print key,(probeset1, probeset2)
                else: critical_exon_junction_db[(probeset1, probeset2)] = critical_exons
    return critical_exon_junction_db

def importmicroRNADataExon(species,array_type,exon_db,microRNA_prediction_method):
    filename = "AltDatabase/"+species+"/"+array_type+"/"+species+"_probeset_microRNAs_"+microRNA_prediction_method+".txt"
    fn=filepath(filename); microRNA_full_exon_db={}; microRNA_count_db={}; gene_microRNA_denom={}
    if array_type != 'exon': critical_exon_junction_db = grabJunctionData(species,array_type,'gene-exon')
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        probeset=t[0];microRNA=t[1]
        try: mir_seq=t[2];mir_sources=t[3]
        except IndexError: mir_seq='';mir_sources=''
        try:
            if array_type == 'exon': ed = exon_db[probeset]; probeset_list = [probeset]; exonid = ed.ExonID()
            else:
                probeset_list = []; uid = probeset ###This is the gene with the critical exon it's mapped to
                for junctions in critical_exon_junction_db[uid]: probeset_list.append(junctions); ed = exon_db[junctions[0]]
            microRNA_info = microRNA, ed.Symbol(), mir_seq, mir_sources
            for probeset in probeset_list:
                try: microRNA_full_exon_db[ed.GeneID(),probeset].append(microRNA_info)
                except KeyError: microRNA_full_exon_db[ed.GeneID(),probeset] = [microRNA_info]
                
                try: microRNA_count_db[microRNA].append(ed.GeneID())
                except KeyError: microRNA_count_db[microRNA] = [ed.GeneID()]
                gene_microRNA_denom[ed.GeneID()] = []
        except KeyError: null=[]
    microRNA_count_db = eliminate_redundant_dict_values(microRNA_count_db)

    ###Replace the actual genes with the unique gene count per microRNA
    for microRNA in microRNA_count_db: microRNA_count_db[microRNA] = len(microRNA_count_db[microRNA])
    print len(gene_microRNA_denom),"genes with a predicted microRNA binding site alinging to a probeset"
    return microRNA_full_exon_db,microRNA_count_db,gene_microRNA_denom
        
def link_microRNA_exon_to_decon_db(microRNA_full_exon_db,exon_hits):
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
    transcript_cdna_sequence_dbase,transcript_associations,exon_sequence_database = import_existing_sequence_build(exon_protein_sequence_file)
    
    critical_exon_seq_export = 'AltDatabase/Mm/AltMouse/AltMouse_critical-exon-seq.txt'
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

def identifyAltIsoformsProteinComp(probeset_gene_db,species,array_type,protein_domain_db,compare_all_features):
    """ This function is used by the module IdentifyAltIsoforms to run 'characterizeProteinLevelExonChanges'"""
    global protein_ft_db; protein_ft_db = protein_domain_db; protein_domain_db=[]
    exon_db={} ### Create a simplified version of the exon_db dictionary with probesets that map to a match and null protein
    for probeset in probeset_gene_db:
        gene, exon_id = probeset_gene_db[probeset]
        ep = ExonProteinAlignmentData(gene,probeset,exon_id,'',''); exon_db[probeset] = ep
    global protein_sequence_db
    if compare_all_features == 'yes': type = 'seqcomp'
    else: type = 'exoncomp'
    exon_protein_sequence_file = 'AltDatabase/'+species+'/'+array_type+'/'+'SEQUENCE-protein-dbase_'+type+'.txt'
    probeset_protein_db,protein_sequence_db = importExonSequenceBuild(exon_protein_sequence_file,exon_db)
    
    exon_hits={}
    for probeset in probeset_protein_db:
        gene = probeset_protein_db[probeset].GeneID()
        exon_hits[gene,probeset]=[]

    include_sequences = 'no' ### Sequences for comparisons are unnecessary to store. List array-type as exon since AltMouse data has been re-organized, later get rid of AltMouse specific functionality in this function
    functional_attribute_db,protein_features = characterizeProteinLevelExonChanges(exon_hits,probeset_protein_db,'exon',include_sequences)

    export_file = 'AltDatabase/'+species+'/'+array_type+'/probeset-domain-annotations-'+type+'.txt' 
    formatAttributeForExport(functional_attribute_db,export_file)
    
    export_file = 'AltDatabase/'+species+'/'+array_type+'/probeset-protein-annotations-'+type+'.txt' 
    formatAttributeForExport(protein_features,export_file)

def formatAttributeForExport(attribute_db,filename):
    import IdentifyAltIsoforms
    export_db={}
    for (gene,probeset) in attribute_db:
        attribute_list = attribute_db[(gene,probeset)]; attribute_list2=[]
        for (attribute,direction) in attribute_list: attribute_list2.append(attribute+'|'+direction)
        export_db[probeset]=attribute_list2
    IdentifyAltIsoforms.exportSimple(export_db,filename,'')
        
def characterizeProteinLevelExonChanges(exon_hits,probeset_protein_db,array_type,include_sequences):    
    functional_attribute_db={}; protein_features={}
    for (array_geneid,uid) in exon_hits: ###uid is probeset or (probeset1,probeset2) value, depending on array_type
            if array_type != 'exon': probeset,probeset2 = uid
            else: probeset = uid
            hv=0
            try:
                pp = probeset_protein_db[probeset]; hv=1
                pos_ref_AC = pp.HitProteinID()
                neg_ref_AC = pp.NullProteinID()
                if array_type != 'exon': ### Instead of using the null_hit, use the second juncition probeset
                    try: np = probeset_protein_db[probeset2]; neg_ref_AC = np.HitProteinID()
                    except KeyError: null =[] ###just use the existing null                
            except KeyError: ###occurs if there is not protein associated with the first probeset (NA for exon arrays)
                if array_type != 'exon': ###Reverse of above
                    try:
                        np = probeset_protein_db[probeset2]; hv=2
                        pos_ref_AC = np.NullProteinID()
                        neg_ref_AC = np.HitProteinID()
                    except KeyError: hv=0
            if hv!=0:
                neg_coding_seq = protein_sequence_db[neg_ref_AC]
                pos_coding_seq = protein_sequence_db[pos_ref_AC]
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
    return functional_attribute_db,protein_features
                                              
def combine_databases(db1,db2):
    for key in db2:
        if key not in db1:
            db1[key] = db2[key]
    return db1

def function_exon_analysis_main(array_type,exon_hits,dataset_name_original,microRNA_full_exon_db,protein_sequence_dbase,probeset_protein_db,protein_ft_dbase):
    ###This function was created so other modules could import ExonAnalyze code and run specific functions
    global dataset_name; global new_results_dir; global refseq_start_stop_seq_db
    global uniprot_ft_db; global protein_ft_db
    dataset_name = dataset_name_original
    new_results_dir = "Exon-Analyze/Additional_results"
    
    protein_ft_db = protein_ft_dbase
    
    global protein_sequence_db
    protein_sequence_db = protein_sequence_dbase
    include_sequences = 'yes'
    functional_attribute_db,protein_features = characterizeProteinLevelExonChanges(exon_hits,probeset_protein_db,array_type,include_sequences)
    exon_attribute_db={}; gene_top_attribute_db={}
    filtered_microRNA_exon_db = link_microRNA_exon_to_decon_db(microRNA_full_exon_db,exon_hits)
        
    return filtered_microRNA_exon_db, exon_attribute_db, gene_top_attribute_db, functional_attribute_db, protein_features

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