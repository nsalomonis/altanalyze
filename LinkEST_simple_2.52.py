import sys, string
import os.path
import unique
import statistics
import copy
import time

################# Parse directory files
dirfile = unique

def read_directory(sub_dir):
    dir=os.path.dirname(dirfile.__file__)
    dir_list = os.listdir(dir + sub_dir)
    #add in code to prevent folder names from being included
    dir_list2 = [] 
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".all" or entry[-5:] == ".data" or entry[-3:] == ".fa":
            dir_list2.append(entry)
    return dir_list2

def filepath(filename):
    dir=os.path.dirname(dirfile.__file__)       #directory file is input as a variable under the main            
    fn=os.path.join(dir,filename)
    return fn

def returnDirectories(sub_dir):
    dir=os.path.dirname(dirfile.__file__)
    dir_list = os.listdir(dir + sub_dir)
    ###Below code used to prevent FILE names from being included
    dir_list2 = []
    for entry in dir_list:
        if "." not in entry: dir_list2.append(entry)
    return dir_list2

def find(s, p):
    # find first occurrence of p in s
    n = len(s)
    m = len(p)
    skip = delta1(p)[p[m-1]]
    i = 0
    while i <= n-m:
        if s[i+m-1] == p[m-1]: # (boyer-moore)
            # potential match
            if s[i:i+m-1] == p[:m-1]:
                return i
            if s[i+m] not in p:
                i = i + m + 1 # (sunday)
            else:
                i = i + skip # (horspool)
        else:
            # skip
            if s[i+m] not in p:
                i = i + m + 1 # (sunday)
            else:
                i = i + 1
    return -1 # not found

class GrabFiles:
    def setdirectory(self,value): self.data = value
    def display(self): print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        files = getDirectoryFiles(self.data,str(search_term))
        if len(files)<1: print search_term,'not found'
        return files
    def returndirectory(self):
        dir_list = getAllDirectoryFiles(self.data)
        return dir_list

def getAllDirectoryFiles(import_dir):
    all_files = []
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        data_dir = import_dir[1:]+'/'+data
        all_files.append(data_dir)
    return all_files

def getDirectoryFiles(import_dir, search_term):
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    matches=[]
    for data in dir_list:    #loop through each file in the directory to output results
        data_dir = import_dir[1:]+'/'+data
        if search_term in data_dir: matches.append(data_dir)
    return matches

def eliminateRedundant(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

################# Classes #################

class UnigeneAssociations:
    def __init__(self,entrez,unigene_id,symbol,accession_numbers):
        self._entrez = entrez; self._unigene = unigene_id; self._symbol = symbol; self._accession = accession_numbers
    def Entrez(self): return self._entrez
    def Unigene(self): return self.unigene
    def Symbol(self): return self._symbol
    def Accession(self): return self._accession
    def OutputValues(self):
        output = self.Entrez()+'|'+self.Symbol()
        return output
    def __repr__(self): return self.OutputValues()

class SequenceRetrieval:
    def __init__(self,accession,unigene_id,seq_type,mRNA_length,call):
        self._accession = accession; self._mRNA_length = mRNA_length; self._unigene_id = unigene_id
        self._seq_type = seq_type; self._call = call
    def UnigeneAC(self): return self._unigene_id
    def Accession(self): return self._accession
    def SetSequenceData(self,seqeunce): self._sequence = seqeunce
    def SequenceData(self): return self._sequence
    def SeqType(self): return self._seq_type
    def mRNALength(self): return self._mRNA_length
    def AlignmentCall(self): return self._call
    def OutputValues(self):
        output = self.Accession()
        return output
    def __repr__(self): return self.OutputValues()

################# Probeset Annoation Parsers #################

def import_annotations(filename,array_type):
    fn=filepath(filename)
    annotation_db = {}; x = 0
    if array_type == 'AltMouse':
        for line in open(fn,'r').readlines():
            data,null = string.split(line,'\n')
            if x == 0: x = 1
            else:
                try: affygene, description, ll_id, symbol, rna_processing_annot = string.split(data,'\t')
                except ValueError:
                    affygene, description, ll_id, symbol = string.split(data,'\t')
                    splicing_annotation = ''
                if '"' in description: null,description,null = string.split(description,'"')
                annotation_db[affygene] = description, symbol,ll_id,rna_processing_annot
    else:
        for line in open(fn,'r').readlines():
            data,null = string.split(line,'\n')
            if x == 0: x = 1
            else:
                try: ensembl, description, symbol, rna_processing_annot = string.split(data,'\t')
                except ValueError: ensembl, description, symbol = string.split(data,'\t')
                annotation_db[ensembl] = description, symbol, ensembl, rna_processing_annot
    return annotation_db

################# EST Sequence Parsers #################

def importESTannotations(filename):
    start_time = time.time()
    output_file = 'AltDatabase/'+species+'/SequenceData/output/'+species+'/'+species+'_cDNA-alignments.txt'
    fn1=filepath(output_file)
    dataw = open(fn1,'w'); previous_ug = ''; probeset_match_db_temp={}; unigene_id=''

    output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/'+species+'_EST_sequences.txt'
    fn2=filepath(output_file)
    datar = open(fn2,'w')    
    """
    output_file = 'AltDatabase/'+species+'/SequenceData/output/'+species+'_probeset-mRNA_relationships.txt'
    fn3=filepath(output_file)
    datamr = open(fn3,'w')
    """
    print "Begining generic fasta import of",filename
    #'>gnl|UG|Mm#S10859962 Mus musculus 12 days embryo spinal ganglion cDNA, RIKEN full-length enriched library, clone:D130006G06 product:unclassifiable, full insert sequence /gb=AK051143 /gi=26094349 /ug=Mm.1 /len=2289']
    #'ATCGTGGTGTGCCCAGCTCTTCCAAGGACTGCTGCGCTTCGGGGCCCAGGTGAGTCCCGC'
    fn=filepath(filename); sequence = '|'; k=0
    for line in open(fn,'r').xreadlines():
        try: data, newline= string.split(line,'\n')
        except ValueError: continue
        if len(data)>0:
            if data[0] != '#':
                try:
                    if data[0] == '>':
                        if test == 'yes':
                            if test_unigene == unigene_id: proceed = 'yes'
                            else: proceed = 'no'
                        else: proceed = 'yes'
                        if proceed == 'yes' and len(sequence) > 1:
                            if unigene_id in unigene_ensembl:
                                    ensembl_list = unigene_ensembl[unigene_id]
                                    for ens_gene in ensembl_list:
                                      if ens_gene in splice_event_db:
                                        sequence = string.upper(sequence)
                                        probeset_seq_data = splice_event_db[ens_gene]
                                        #y = SequenceRetrieval(accession,unigene_id,sequence[1:],seq_type)
                                        #fasta[accession] = y
                                        est_seq = sequence[1:]
                                        mRNA_length = len(est_seq)
                                        results = simpleSeqMatchProtocol(probeset_seq_data,est_seq)
                                        for (call,probeset) in results:
                                            #y = SequenceRetrieval(accession,unigene_id,seq_type,mRNA_length,call)
                                            #try: probeset_match_db_temp[probeset].append((call,mRNA_length,y))
                                            #except KeyError: probeset_match_db_temp[probeset] = [(call,mRNA_length,y)]
                                            dataw.write(string.join([probeset,str(call),str(mRNA_length),accession,seq_type],'\t')+'\n')
                            try:
                                ###Save all sequences to the disk rather than store these in memory. Just select the optimal sequences later.
                                values = [accession,est_seq]
                                values = string.join(values,'\t')+'\n'
                                datar.write(values); k+=1
                            except UnboundLocalError: null=[]
                            if 'NM_' in data: seq_type = 'full-length'
                            elif 'mRNA' in data: seq_type = 'full-length'
                            elif 'full-length' in data: seq_type = 'full-length'
                            elif 'complete cds' in data: seq_type = 'full-length'
                            elif 'cDNA' in data: seq_type = 'full-length'
                            elif 'clone' not in data: seq_type = 'full-length'
                            else: seq_type = 'clone'
                            previous_ug = unigene_id
                            z = string.find(data,'ug='); ug_info = string.split(data[z+3:],' /'); unigene_id = ug_info[0]
                            x = string.find(data,'gb='); gb_info = string.split(data[x+3:],' /'); accession = gb_info[0]
                            sequence = '|'; current_ug = unigene_id
                            """
                            ###Analyze data for one Unigene cluster if it differs from the previous
                            if current_ug != previous_ug:
                                for probeset in probeset_match_db_temp:
                                    try: probeset_match_db_temp[probeset] += probeset_match_db[probeset]
                                    except KeyError: null = ''
                                identifyCorrespondingESTs(probeset_match_db_temp,datamr)
                                probeset_match_db_temp={} ###re-set this database
                                """
                        else:   ###only used for the first entry, then the above if is always used
                            if 'NM_' in data: seq_type = 'full-length'
                            elif 'clone' not in data: seq_type = 'full-length'
                            else: seq_type = 'clone'; previous_ug = unigene_id
                            z = string.find(data,'ug='); ug_info = string.split(data[z+3:],' /'); unigene_id = ug_info[0]; current_ug = unigene_id
                            x = string.find(data,'gb='); gb_info = string.split(data[x+3:],' /'); accession = gb_info[0];  sequence = '|'
                            continue
                except IndexError:
                                if unigene_id in unigene_ensembl:
                                    ensembl_list = unigene_ensembl[unigene_id]
                                    for ens_gene in ensembl_list:
                                      if ens_gene in splice_event_db:
                                        probeset_seq_data = splice_event_db[ens_gene]
                                        #y = SequenceRetrieval(accession,unigene_id,sequence[1:],seq_type)
                                        #fasta[accession] = y
                                        est_seq = sequence[1:]
                                        mRNA_length = len(est_seq)
                                        results = simpleSeqMatchProtocol(probeset_seq_data,est_seq)
                                        for (call,probeset) in results:
                                            #y = SequenceRetrieval(accession,unigene_id,seq_type,mRNA_length,call)
                                            #try: probeset_match_db[probeset].append((call,mRNA_length,y))
                                            #except KeyError: probeset_match_db[probeset] = [(call,mRNA_length,y)]
                                            dataw.write(string.join([probeset,str(call),str(mRNA_length),accession,seq_type],'\t')+'\n')
                                try:
                                    ###Save all sequences to the disk rather than store these in memory. Just select the optimal sequences later.
                                    values = [accession,est_seq]
                                    values = string.join(values,'\t')+'\n'
                                    datar.write(values)
                                except UnboundLocalError: null=[]
                                """
                                ###Analyze data for one Unigene cluster if it differs from the previous
                                if current_ug != previous_ug:
                                    for probeset in probeset_match_db_temp:
                                        try: probeset_match_db_temp[probeset] += probeset_match_db[probeset]
                                        except KeyError: null = ''
                                    identifyCorrespondingESTs(probeset_match_db_temp,datamr)
                                    probeset_match_db_temp={} ###re-set this database
                                continue"""
                try:
                    if data[0] != '>': sequence = sequence + data
                except IndexError:
                    print dog; continue
    datar.close(); dataw.close() #; datamr.close()
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "Unigene mRNA sequences analyzed in %d seconds" % time_diff    
    print "Number of imported sequences:", k

def combineDBs(db1,db2):
    for i in db2:
        try: db1[i]+=db2[i]
        except KeyError: db1[i]=db2[i]
    return db1

def importEnsemblTranscriptSequence(filename):
    start_time = time.time()
    output_file = 'AltDatabase/'+species+'/SequenceData/output/'+species+'_Ens-alignments.txt'
    fn1=filepath(output_file)
    dataw = open(fn1,'w')
    
    output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/'+species+'_Ens_sequences.txt'
    fn2=filepath(output_file)
    datar = open(fn2,'w')
    
    print "Begining generic fasta import of",filename
    fn=filepath(filename)
    sequence = ''; x = 0; count = 0; global gene_not_found; gene_not_found=[]; genes_found={}
    for line in open(fn,'r').xreadlines():
        #for line in open(fn)(): ###try this out, redundant
        exon_start=1; exon_stop=1
        try: data, newline= string.split(line,'\n')
        except ValueError: continue
        try:
            if data[0] == '>':
                    if len(sequence) > 0:
                        gene_found = 'no'
                        if test == 'yes':
                            if ensembl_id == test_ensembl: proceed = 'yes'
                            else: proceed = 'no'
                        else: proceed = 'yes'
                        if proceed == 'yes':
                            if ensembl_id in splice_event_db:
                                gene_found = 'yes'; genes_found[ensembl_id]=[]
                                seq_type = 'full-length'
                                probeset_seq_data = splice_event_db[ensembl_id]
                                cDNA_seq = sequence[1:]
                                mRNA_length = len(cDNA_seq)
                                results = simpleSeqMatchProtocol(probeset_seq_data,cDNA_seq)
                                for (call,probeset) in results:
                                        #y = SequenceRetrieval(transid,ensembl_id,seq_type,mRNA_length,call)
                                        #try: probeset_match_db[probeset].append((call,mRNA_length,y))
                                        #except KeyError: probeset_match_db[probeset] = [(call,mRNA_length,y)]
                                        dataw.write(string.join([probeset,str(call),str(mRNA_length),transid,seq_type],'\t')+'\n')
                        if gene_found == 'yes':
                            ###Save all sequences to the disk rather than store these in memory. Just select the optimal sequences later. 
                            try: values = [transid,cDNA_seq]
                            except UnboundLocalError:
                                print ensembl_id
                                if ensembl_id in splice_event_db: print 'yes'
                                else: print 'no'
                                print transid;kill
                            values = string.join(values,'\t')+'\n'
                            datar.write(values); x+=1
                        else: gene_not_found.append(ensembl_id)
                        sequence = ''; t= string.split(data,'|')
                        try: ensembl_id,chr,strand,transid,prot_id = t
                        except ValueError: ensembl_id,chr,strand,transid = t
                        ensembl_id = ensembl_id[1:]
                    else: ###Occurs for the first entry
                        t= string.split(data,'|')
                        try: ensembl_id,chr,strand,transid,prot_id = t
                        except ValueError: ensembl_id,chr,strand,transid = t
                        ensembl_id = ensembl_id[1:]
        except IndexError: continue
        try:
            if data[0] != '>': sequence = sequence + data
        except IndexError:  continue
    datar.close(); dataw.close()
    end_time = time.time(); time_diff = int(end_time-start_time)
    gene_not_found = unique.unique(gene_not_found)
    print len(genes_found), 'genes associated with Ensembl transcripts'
    print len(gene_not_found), "gene's not found in the Ensembl-probeset database (should be there unless conflict present)"
    print gene_not_found[0:10],'not found examples'
    print "Ensembl transcript sequences analyzed in %d seconds" % time_diff

def importEnsemblProteinSeq(export_seq):
    filename = 'AltDatabase/'+species+'/SequenceData/'+species+'_EnsProt_sequence.txt'
    fn=filepath(filename); ensembl_prot_seq_db = {}; ensembl_prot_seq_len_db = {}
    sequence = ''
    for line in open(fn,'r').xreadlines():
        try: data, newline= string.split(line,'\n')
        except ValueError: continue
        try:
            if data[0] == '>':
                    if len(sequence) > 0:
                        if len(ens_prot_id)>0: ###occurs once in the file
                            ensembl_prot_seq_db[ens_prot_id] = sequence[:-1]
                            ensembl_prot_seq_len_db[ens_prot_id] = len(sequence[:-1])
                        ens_prot_id = data[1:]
                        sequence = ''
                    else:
                        sequence = ''
                        ens_prot_id = data[1:] ###Occurs for the first entry  
        except IndexError: continue
        try:
            if data[0] != '>': sequence = sequence + data
        except IndexError:  continue
    ensembl_prot_seq_db[ens_prot_id] = sequence[:-1]  ###Need this for the last entry
    ensembl_prot_seq_len_db[ens_prot_id] = len(sequence[:-1])
    if export_seq == 0: return ensembl_prot_seq_len_db
    else: return ensembl_prot_seq_db

def importNCBIProteinSeq(protein_linked_ACs):
    import_dir = '/AltDatabase/SequenceData' ### Multi-species fiel
    g = GrabFiles(); g.setdirectory(import_dir)
    seq_files = g.searchdirectory('fsa_aa')

    try: filename = seq_files[0]
    except IndexError: filename = ''
    fn=filepath(filename); ncbi_prot_seq_db = {}
    sequence = ''
    for line in open(fn,'r').xreadlines():
        try: data, newline= string.split(line,'\n')
        except ValueError: continue
        try:
            if data[0] == '>':
                if len(sequence) > 0:
                    if len(prot_id)>0 and prot_id in protein_linked_ACs: ###occurs once in the file
                        ncbi_prot_seq_db[prot_id] = sequence
                    """else:
                        print prot_id
                        for i in protein_linked_ACs:
                            print i;kill"""
                    t= string.split(data,'|')
                    prot_id = t[3][:-2]
                    sequence = ''
                else:
                    sequence = ''
                    t= string.split(data,'|')
                    prot_id = t[3][:-2] ###Occurs for the first entry  
        except IndexError: continue
        try:
            if data[0] != '>': sequence = sequence + data
        except IndexError:  continue
    if len(prot_id)>0 and prot_id in protein_linked_ACs:
        ncbi_prot_seq_db[prot_id] = sequence  ###Need this for the last entry
    print len(ncbi_prot_seq_db), "genbank protein ACs associated with sequence, out of", len(protein_linked_ACs)
    return ncbi_prot_seq_db
    
def import_ensembl_unigene(filename):
    fn=filepath(filename)
    for line in open(fn,'r').xreadlines():        
        data, newline= string.split(line,'\n')
        ensembl,unigene = string.split(data,'\t')
        try: unigene_ensembl[unigene].append(ensembl)
        except KeyError: unigene_ensembl[unigene] = [ensembl]
    
def simpleSeqMatchProtocol(probeset_seq_data,est_seq):
    results=[]
    for y in probeset_seq_data:
        try: exon_seq = y.ExonSeq(); probeset = y.Probeset()
        except AttributeError: probeset = y.Probeset(); print [probeset],'seq missing';kill
        """if probeset == 'G7203664@J946608_RC@j_at':
            print exon_seq
            print est_seq;kill"""
        if exon_seq in est_seq: call=1
        elif array_type == 'exon':
            if exon_seq[:25] in est_seq: call=1
            elif exon_seq[-25:] in est_seq: call=1
            else: call = 0
        else: call = 0
        #else: junction_seq = y.JunctionSeq()
        results.append((call,probeset))
    return results

def importUCSCSequences(determine_transcript_length,mRNA_probesets_to_analyze,match_type,var):
    start_time = time.time()
    filename = 'AltDatabase/'+species+'/SequenceData/mrna.fa'
    output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/'+species+'_UCSC_'+str(var)+'_sequences.txt'
    fn2=filepath(output_file)
    datar = open(fn2,'w')    
    
    print "Begining generic fasta import of",filename
    #'>gnl|ENS|Mm#S10859962 Mus musculus 12 days embryo spinal ganglion cDNA, RIKEN full-length enriched library, clone:D130006G06 product:unclassifiable, full insert sequence /gb=AK051143 /gi=26094349 /ens=Mm.1 /len=2289']
    #'ATCGTGGTGTGCCCAGCTCTTCCAAGGACTGCTGCGCTTCGGGGCCCAGGTGAGTCCCGC'
    fn=filepath(filename); sequence = '|'; ucsc_mRNA_hit_len={}; ucsc_probeset_null_hits={}; k=0
    for line in open(fn,'r').xreadlines():
        try: data, newline= string.split(line,'\n')
        except ValueError: continue
        if len(data)>0:
            if data[0] != '#':
                try:
                    if data[0] == '>':
                        """if test == 'yes':
                            if test_ensembl == ensembl_id: proceed = 'yes'
                            else: proceed = 'no'"""
                        proceed = 'yes'
                        if proceed == 'yes' and len(sequence) > 1:
                            if accession in mRNA_probesets_to_analyze:
                                sequence = string.upper(sequence)
                                mRNA_seq = sequence[1:]; mRNA_length = len(mRNA_seq)
                                if accession in mRNA_probesets_to_analyze:
                                    k+=1
                                    probeset_seq_data = mRNA_probesets_to_analyze[accession]
                                    results = simpleSeqMatchProtocol(probeset_seq_data,mRNA_seq)
                                    for (call,probeset) in results:
                                        #y = SequenceRetrieval(transid,ensembl_id,seq_type,mRNA_length,call)
                                        #try: probeset_match_db[probeset].append((call,mRNA_length,y))
                                        #except KeyError: probeset_match_db[probeset] = [(call,mRNA_length,y)]
                                        dataw.write(string.join([probeset,str(call),str(mRNA_length),transid,seq_type],'\t')+'\n')
                            if accession in determine_transcript_length:
                                sequence = string.upper(sequence)
                                mRNA_seq = sequence[1:]; mRNA_length = len(mRNA_seq)
                                ucsc_mRNA_hit_len[accession] = mRNA_length ###record the length of the mRNA, specifically for this entry
                                values = [accession,mRNA_seq]; values = string.join(values,'\t')+'\n'
                                datar.write(values)
                        values = string.split(data,' '); accession = values[0][1:]
                        sequence = '|'; continue
                except IndexError: null = []       
                try:
                    if data[0] != '>': sequence = sequence + data
                except IndexError: print kill; continue
    datar.close()
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "UCSC mRNA sequences analyzed in %d seconds" % time_diff
    print "The mRNA sequence length for",len(ucsc_mRNA_hit_len),"accession numbers obtained (only for multiple hits inferred from location match)"
    print len(ucsc_probeset_null_hits),"probesets found that have new null hits identified via sequence match from UCSC mRNAs"
    return ucsc_mRNA_hit_len,probeset_match_db

def importEnsExonTranscriptAssociations():
    ###This function is used to extract out EnsExon to EnsTranscript relationships to find out directly
    ###which probesets associate with which transcripts and then which proteins
    filename = 'AltDatabase/'+species+'/SequenceData/'+species+'_Ensembl_transcript-annotations.txt'
    fn=filepath(filename); ens_exon_to_transcript={}; ens_gene_to_transcript={}
    for line in open(fn,'r').readlines():
        data = string.replace(line,'\n',''); data = string.replace(data,'\r','')
        t = string.split(data,'\t')
        ens_geneid,chr,strand,exon_start,exon_end,ens_exonid,constitutive_exon,ens_transcript_id = t
        try: ens_exon_to_transcript[ens_exonid].append(ens_transcript_id)
        except KeyError:ens_exon_to_transcript[ens_exonid] = [ens_transcript_id]
        try: ens_gene_to_transcript[ens_geneid].append(ens_transcript_id)
        except KeyError:ens_gene_to_transcript[ens_geneid] = [ens_transcript_id]
    ens_exon_to_transcript = eliminateRedundant(ens_exon_to_transcript)
    ens_gene_to_transcript = eliminateRedundant(ens_gene_to_transcript)  
    return ens_exon_to_transcript,ens_gene_to_transcript

def importUCSCTranscriptAssociations():
    ###This function is used to extract out EnsExon to EnsTranscript relationships to find out directly
    ###which probesets associate with which transcripts and then which proteins
    filename = 'AltDatabase/'+species+'/SequenceData/'+species+'_UCSC-accession-to-gene_mrna.txt'
    fn=filepath(filename); ucsc_mrna_to_gene={}
    for line in open(fn,'r').readlines():
        data = string.replace(line,'\n',''); data = string.replace(data,'\r','')
        t = string.split(data,'\t')
        mRNA_ac, ens_genes, num_ens_exons = t
        #if mRNA_ac not in accession_index: ###only add mRNAs not examined in UniGene
        ens_genes = string.split(ens_genes,'|')
        ucsc_mrna_to_gene[mRNA_ac] = ens_genes
    return ucsc_mrna_to_gene

def characterizeProteinLevelExonChanges(uid,array_geneid,hit_data,null_data,array_type):    
            domains=[]; functional_attribute=[]
            if array_type != 'exon': probeset,probeset2 = uid
            else: probeset = uid
            hv=0
            try:
                hv=1
                pos_ref_AC = hit_data[0]
                neg_ref_AC = null_data[0]
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
                neg_coding_seq = null_data[1]
                pos_coding_seq = hit_data[1]
                neg_length = len(neg_coding_seq)
                pos_length = len(pos_coding_seq)
                pos_length = float(pos_length);  neg_length = float(neg_length)
                if array_geneid in protein_ft_db:
                        protein_ft = protein_ft_db[array_geneid]
                        neg_ft_missing,pos_ft_missing = compareProteinFeatures(protein_ft,neg_coding_seq,pos_coding_seq)
                        for (pos,blank,ft_name,annotation) in pos_ft_missing:
                            domains.append(ft_name)
                        for (pos,blank,ft_name,annotation) in neg_ft_missing:
                            ###If missing from the negative list, it is present in the positive state
                            domains.append(ft_name)
                if pos_coding_seq[0:11] != neg_coding_seq[0:11]:
                    function_var = 'alt-N-terminus'
                    functional_attribute.append(function_var)
                if pos_coding_seq[-11:] != neg_coding_seq[-11:]:
                    function_var = 'alt-C-terminus'
                    functional_attribute.append(function_var)
                ### Record change in peptide size
                protein_length_diff = abs(pos_length-neg_length)
            results = [len(functional_attribute),len(domains),protein_length_diff,pos_ref_AC,neg_ref_AC] ###number of domains, number N or C-terminal differences, protein length difference
            return results

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

def reAnalyzeRNAProbesetMatches(align_files,probeset_seq_file):
    ###Note: what is missing are transcripts that overlap with Ensembl genes that are not in Unigene.
    ###wether to include these or not is unclear, but could do it be searching specifically for alignments to UCSC
    ###transcripts from the external-exon field of Hs_probeset-Ensembl.txt.
    probeset_index = {}
    ###determine which probesets had either no hit or no null and further analyze these
    accession_index = {}; match_db = {}; null_match_db = {}
    print "Begining to parse ~200,000,000 entries"
    align_files.sort(); align_files = [align_files[0]]
    for filename in align_files:
        print 'Reading',filename
        if 'Ens' in filename: filetype = 'Ensembl'
        else: filetype = 'EST'
        start_time = time.time(); ens_probeset_match_db={}; ens_probeset_null_match_db={}
        fn=filepath(filename); x = 0; z = 400000
        for line in open(fn,'r').xreadlines():
            values = string.replace(line,'\n',''); x+=1 
            probeset,call,mRNA_length,accession,seq_type = string.split(values,'\t')
            mRNA_length = int(mRNA_length)
            if x == z:
                print z, 'entries parsed'
                z = z+400000#; break
            if call == '1': #and probeset == '2523723'
                try: ens_probeset_match_db[probeset].append(accession)
                except KeyError: ens_probeset_match_db[probeset] = [accession]
                probeset_index[probeset]=[]
            elif call == '0':
                try: ens_probeset_null_match_db[probeset].append(accession)
                except KeyError: ens_probeset_null_match_db[probeset] = [accession]
                probeset_index[probeset]=[]

        end_time = time.time(); time_diff = int(end_time-start_time)
        print "databases built in %d seconds" % time_diff

    ### Get probesets only that are either only a null or only a hit
    match_in_both_db = {}; match_in_one_only_db={}
    for probeset in probeset_index:
        if probeset in ens_probeset_null_match_db and probeset in ens_probeset_match_db: match_in_both_db[probeset] = []
        else: match_in_one_only_db[probeset]=[]

    ###To mirror the previous probeset-protein database structure, we need additional exon information data
    if array_type == 'exon':
        probeset_annotations_file = 'AltDatabase/'+species+'/SequenceData/'+species+'_Ensembl_probesets.txt'
        import ExonSeqModule
        exon_db = ExonSeqModule.importSplicingAnnotationDatabase(probeset_annotations_file)
        #exon_db = ExonSeqModule.getExonAnnotationsAndSequence(probeset_seq_file,species)
    else:
        import JunctionSeqModule; exon_db = JunctionSeqModule.importSplicingAnnotationDatabaseAndSequence(species,array_type,'probeset')
       
    ### Do the same for any probeset not matched at all by Ensembl
    for probeset in exon_db:
        if probeset not in probeset_index: match_in_one_only_db[probeset] = []

    ### Get Ensembl protein sequence and transcript-protein   
    ensembl_prot_seq_db = importEnsemblProteinSeq(1)
    ens_transcript_to_protein_db,null = importEnsProtRelationships(species); null=[]

    minimal_effect_db={}; accession_seq_db={};start_time = time.time()
    print len(match_in_both_db), "probesets versus", len(match_in_one_only_db), "probesets, that have a match and null in Ensembl"
    for probeset in match_in_both_db:
        geneid = exon_db[probeset].GeneID()
        match_list = ens_probeset_match_db[probeset]; null_list = ens_probeset_null_match_db[probeset]; prot_match_list=[]; prot_null_list=[]
        for match_ac in match_list:
            if match_ac in ens_transcript_to_protein_db:
                protein_ac = ens_transcript_to_protein_db[match_ac]
                protein_seq = ensembl_prot_seq_db[protein_ac]
                prot_match_list.append([protein_ac,protein_seq])
        for null_ac in null_list:
            if null_ac in ens_transcript_to_protein_db:
                protein_ac = ens_transcript_to_protein_db[null_ac]
                protein_seq = ensembl_prot_seq_db[protein_ac]
                prot_null_list.append([protein_ac,protein_seq])
        if len(prot_match_list)>0 and len(prot_null_list)>0:
            results_list=[]
            for mi in prot_match_list:
                for ni in prot_null_list:
                    results = characterizeProteinLevelExonChanges(probeset,geneid,mi,ni,array_type)
                    results_list.append(results)
            results_list.sort()
            hit_ac = results_list[0][-2]
            null_ac = results_list[0][-1]
            minimal_effect_db[probeset] = [hit_ac,null_ac]
            accession_seq_db[hit_ac] = [ensembl_prot_seq_db[hit_ac]]
            accession_seq_db[null_ac] = [ensembl_prot_seq_db[null_ac]]
            #print results_list[0][-2],results_list[0][-1]
            #print probeset,len(prot_match_list),len(prot_null_list),results_list;kill
        else: match_in_one_only_db[probeset] = []

    print len(minimal_effect_db),"optimal pairs of probesets linked to Ensembl found."
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "databases built in %d seconds" % time_diff

    exportSimple(minimal_effect_db,species+'_Ensembl_protein_pairs.txt')
    exportSimple(accession_seq_db,'sequences/'+species+'_match_null_Ensembl_seq.txt')
    kill
    k=0; l=0
    for probeset in exon_db:
        ###See if there are Ensembl transcript hits by location and not by sequence
        proceed = []
        if (probeset not in probeset_match_ENS): proceed.append(1)
        if (probeset in probeset_match_null_ENS): proceed.append(2)
        if (probeset in probeset_match_ENS): proceed.append(3)
        if len(proceed)>0:
            ens_transcript_matches = []; ens_transcript_matches2 = []
            if array_type == 'exon':
                pd = exon_db[probeset]
                exons = pd.ExternalExonIDList()
                for exon in exons:
                    if 'ENS' in exon:
                        ens_transcripts = ens_exon_to_transcript[exon]
                        if probeset == '2594513': print '111-AAAAA',exons, exon, ens_transcripts, proceed
                        ens_transcript_matches+=ens_transcripts
            else:
                try: ens_transcript_matches = ens_probeset_match_db[probeset]  ###use this as an alternative for junction arrays, where probesets don't match to a single exon (could use for exon analysis as well)
                except KeyError: ens_transcript_matches=[]
            ens_transcript_matches = unique.unique(ens_transcript_matches)
            if len(ens_transcript_matches)>0: ###Thus, there are Ensembl matches
                if (1 in proceed) or (3 in proceed):                    
                    for ens_transcript in ens_transcript_matches:
                        if ens_transcript in ens_transcript_len:
                            mRNA_len = ens_transcript_len[ens_transcript] ###length recorded when importing all probeset to transcript associations
                            ens_transcript_matches2.append((mRNA_len,ens_transcript)) ###only add it if it's in the database, if not, there can be problems down-the-line
                    ens_transcript_matches2.sort()
                    if len(ens_transcript_matches2)>0:
                        ens_protein_matches = [] ###The longest transcript does not always mean the longest protein (e.g. probeset: 2594513, ENST00000337057 is the appropriate match not ENST00000368127 at the protein level)
                        if probeset == '2594513': print 'AAAA',ens_transcript_matches2, ens_transcript_matches, ens_transcript_matches2
                        for (ml,ens_transcript) in ens_transcript_matches2:
                            if ens_transcript in ens_transcript_to_protein_db:
                                ens_protein_id = ens_transcript_to_protein_db[ens_transcript]
                                if probeset == '2594513': print '!!!!!AAAA',[ens_protein_id]
                                if ens_protein_id in ensembl_prot_seq_len_db:
                                    pl = ensembl_prot_seq_len_db[ens_protein_id]
                                    ens_protein_matches.append((pl,ens_transcript))
                                    if probeset == '2594513': print 'BBBB',ens_protein_matches
                        ens_protein_matches.sort()
                        if probeset == '2594513': print 'CCCCC',ens_protein_matches
                        if len(ens_protein_matches)>0:
                            ens_transcript_matches2 = ens_protein_matches
                            if probeset == '2594513': print 'DDDD',ens_transcript_matches2
                        if probeset == '2594513': print 'EEEE',ens_transcript_matches2
                        ###Create a new entry with that probeset information
                        probeset_match_ENS[probeset]=ens_transcript_matches2[-1]; k+=1
                        match_db[probeset] = []
                if 2 in proceed:
                    ens_transcript_null_matches = []
                    ml, ens_transcript = probeset_match_null_ENS[probeset]
                    if array_type == 'exon': geneid = pd.GeneID()
                    else: geneid = ens_transcript_to_gene[ens_transcript_matches[0]] ###Any of the transcript matches is a correct gene match as well (used for non-exon arrays where a one-to-one exon association may not be present)
                    for ens_transcript in ens_gene_to_transcript[geneid]:
                        if ens_transcript not in ens_transcript_matches: ####Find transcripts NOT in the matches
                            if ens_transcript in ens_transcript_len:
                                mRNA_len = ens_transcript_len[ens_transcript] ###length recorded when importing all probeset to transcript associations
                                ens_transcript_null_matches.append((mRNA_len,ens_transcript))  ###This is really NULL matches
                    ens_transcript_null_matches.sort()
                    ens_protein_matches = [] ###The longest transcript does not always mean the longest protein (e.g. probeset: 2594513, ENST00000337057 is the appropriate match not ENST00000368127 at the protein level)
                    for (ml,ens_transcript) in ens_transcript_null_matches:
                        if ens_transcript in ens_transcript_to_protein_db:
                            ens_protein_id = ens_transcript_to_protein_db[ens_transcript]
                            if ens_protein_id in ensembl_prot_seq_len_db:
                                pl = ensembl_prot_seq_len_db[ens_protein_id]
                                ens_protein_matches.append((pl,ens_transcript))  ###NULL matching transcripts
                    ens_protein_matches.sort()
                    
                    if ens_transcript not in ens_transcript_matches and len(ens_protein_matches)>0:
                        ###Check to see if there are valid proteins that are longer than the null selected (can occur if the mRNA is longer)
                        if probeset == '2594513': print 'KKKK1',ens_protein_matches
                        probeset_match_null_ENS[probeset]=ens_protein_matches[-1]
                        
                    ###check to see if the transcript WHICH SHOULD NOT BE MATCH THE PROBESET, matches the probeset
                    if ens_transcript in ens_transcript_matches:
                        ###Thus the entry is not valid (probeset matches to transcript - happended for (probset: 2594513-ENSG00000215573)
                        if probeset == '2594513': print 'FFFF1',ens_transcript_matches
                        if probeset == '2594513': print 'FFFF2',ens_transcript_null_matches
                        if len(ens_transcript_null_matches)>0:
                            ###Then one of these is the NEW NULL Ens Transcript
                            if len(ens_protein_matches)>0: ens_transcript_null_matches = ens_protein_matches
                            if probeset == '2594513': print 'GGGG',ens_transcript_null_matches
                            ###REPLACE EXISTING ENTRY with transcripts that DON'T MATCH
                            probeset_match_null_ENS[probeset]=ens_transcript_null_matches[-1]; l+=1
                        else:
                            del probeset_match_null_ENS[probeset]
                            if (probeset not in probeset_match_null_FL) and (probeset not in probeset_match_null_EST):
                                del null_match_db[probeset]  ###allows the to be further analysed in the context of the UCSC data

    print k, "probesets associated with Ensembl transcripts by location and not sequence"
    print l, "probesets replaced, where the transcript did match WHEN PREDICTED NULL"
    ###Find probesets with no null or no hits
    for probeset in null_match_db:
        try: h = match_db[probeset]
        except KeyError: not_aligned_hit_probesets[probeset] = 'null'
    for probeset in match_db:
        try: h = null_match_db[probeset]
        except KeyError: not_aligned_null_probesets[probeset] = 'hit'

    probeset = '2594513'
    if probeset in not_aligned_hit_probesets:  print '%$$$$$$$$$$$$-',probeset,'not_aligned_hit_probesets',not_aligned_hit_probesets[probeset]
    if probeset in not_aligned_null_probesets: print '%$$$$$$$$$$$$-',probeset,'not_aligned_null_probesets',not_aligned_null_probesets[probeset]

    if probeset == '2594513':
        if probeset in probeset_match_FL: print 'FL',probeset_match_FL[probeset]
        if probeset in probeset_match_EST: print 'FL',probeset_match_EST[probeset]    
    ###for probesets with no matching transcript, check to see if we previously associated the probeset with an exon
    determine_transcript_length={}; determine_transcript_length_probesets={}; ensembl_probeset_no_seq={}; k=0
    for probeset in exon_db: ###look at all rather than just in 'not_aligned_hit_probesets'
        pd = exon_db[probeset]
        exons = pd.ExternalExonIDList()
        if probeset == '2594513': print len(exons), exons[0], len(exons[0])
        """if len(exons)==1:
            exon = exons[0]
            if len(exon)>1: ###otherwise there is no exon associations
                if '-' in exon:
                    ucsc_transcript,null = string.split(exon,'-')
                    probeset_match_FL[probeset] = 10000,ucsc_transcript; k+=1  #length is arbitrary
                    
                else:
                    if probeset in not_aligned_hit_probesets:
                        ###need to search for probeset alignments, even if aligned to Ensembl, when sequence is missing from our database.
                        ensembl_probeset_no_seq[probeset]=[]
                        if probeset == '2594513': print len(exons), exons[0], len(exons[0]),ensembl_probeset_no_seq[probeset]"""
        if len(exons)>0:  ###This was for multiple exon alignments, but if we replaced the previous FL with a new one, without sequence this caused problems - therefore require all have sequence before using
            ucsc_transcript_matches = []; no_sequence_ensembls=[]
            for exon in exons:
                if len(exon)>1:
                    if '-' in exon: ucsc_transcript_matches.append(exon[:-2])
                    else:
                        if probeset in not_aligned_hit_probesets: ensembl_probeset_no_seq[probeset]=[]
            if len(ucsc_transcript_matches)>0:
                ###Look up the length when importing all sequences and record it
                determine_transcript_length_probesets[probeset] = ucsc_transcript_matches
                for transcript in ucsc_transcript_matches: determine_transcript_length[transcript] = []
        ###This effectively takes care of probesets with no hit (although we will have to determine the length of some mRNAs and add these back after sequence import)
        ###However, Ensembl exons with no transcript links CURRENTLY IN OUR DATABASE need to be linked to UCSC accession numbers NOT IN ENSEMBL (see sequence match below)
    
    if '2594513' in ensembl_probeset_no_seq: print '***2594513 in ensembl_probeset_no_seq',ensembl_probeset_no_seq['2594513']

    print k, "missing probeset mRNA hits, added to database based on UCSC location based associations"
    print len(determine_transcript_length), 'for which length needs to be determined'

    ###Need to determine which UCSC transcripts DO NOT contain this probeset.  The best way to do this is by sequence matching.
    ###First, determine which transcripts we need to look at (since we've looked at alot already, we can eliminate these). This is done by only looking
    ###at transcripts NOT in "accession_index", that are associated with Ensembl gene ids (exported from the ExpressionBuilder module UCSC).
    ucsc_mrna_to_gene = importUCSCTranscriptAssociations()     
    ###Next, determine which probesets should be compared to which UCSC mRNA transcripts
    ucsc_gene_to_mrna={}; mRNA_probesets_to_analyze={}; ucsc_gene_to_mrna_non_optimal={}; k=0
    for mRNA in ucsc_mrna_to_gene:
        gene_ids = ucsc_mrna_to_gene[mRNA]
        if len(gene_ids)==1: ###thus the mRNA doesn't overlap with multiple Ensembl gene IDs (some cases where this might cause problems)
            geneid = gene_ids[0]
            k+=1
            try: ucsc_gene_to_mrna[geneid].append(mRNA)
            except KeyError: ucsc_gene_to_mrna[geneid] = [mRNA]
        else: 
            for geneid in gene_ids: ###Use these for a second round of sequence query if the first round fails
                try: ucsc_gene_to_mrna_non_optimal[geneid].append(mRNA)
                except KeyError: ucsc_gene_to_mrna_non_optimal[geneid] = [mRNA]
    #print 'XXXXXX ucsc_gene_to_mrna[ENSG00000080298]',ucsc_gene_to_mrna['ENSG00000080298']

    for probeset in not_aligned_null_probesets:
        pd = exon_db[probeset]
        if array_type == 'exon': geneids = [pd.GeneID()]
        else: geneids = pd.GeneID()
        for geneid in geneids:
            if geneid in ucsc_gene_to_mrna:
                mRNAs = ucsc_gene_to_mrna[geneid]
                for mRNA in mRNAs:
                    if probeset in not_sequence_characterized:
                        try: mRNA_probesets_to_analyze[mRNA].append(pd)
                        except KeyError: mRNA_probesets_to_analyze[mRNA] = [pd]
                    else:  ###only look at mRNA's not examined by ESTMatch if looked at already
                        if mRNA not in accession_index: 
                            try: mRNA_probesets_to_analyze[mRNA].append(pd)
                            except KeyError: mRNA_probesets_to_analyze[mRNA] = [pd]

    #print 'XXXXXX mRNA_probesets_to_analyze[AK093692]',mRNA_probesets_to_analyze['AK093692']    
    for mRNA in mRNA_probesets_to_analyze:
        for pd in mRNA_probesets_to_analyze[mRNA]:
            if pd.Probeset() == '2594513': print '###### 2594513 mRNA examined is',mRNA

    print k, "mRNAs from UCSC with only one gene associated"
    print len(not_aligned_hit_probesets), 'probesets and',len(mRNA_probesets_to_analyze),'mRNAs being compared to identify null probeset mRNA matches from UCSC'
    ###get the length of the first and second sets and hit type for the second (only care if it doesn't match).
    match_type = 0 ###Find accession that LACK the probeset sequence
    importUCSCSequences(determine_transcript_length,mRNA_probesets_to_analyze,match_type)

    ###if the above sequence match queries didn't work, send back using sub-optimal gene associations
    mRNA_probesets_to_analyze={}
    for probeset in not_aligned_null_probesets:
        if probeset not in probeset_match_null_FL:
            ###probesets that only had sequence matches associated after UCSC sequence query
            pd = exon_db[probeset]
            if array_type == 'exon': geneids = [pd.GeneID()]
            else: geneids = pd.GeneID()
            for geneid in geneids:
                if geneid in ucsc_gene_to_mrna_non_optimal:
                    mRNAs = ucsc_gene_to_mrna_non_optimal[geneid]
                    for mRNA in mRNAs:
                        if probeset in not_sequence_characterized:
                            try: mRNA_probesets_to_analyze[mRNA].append(pd)
                            except KeyError: mRNA_probesets_to_analyze[mRNA] = [pd]
                        else:  ###only look at mRNA's not examined by ESTMatch if looked at already
                            if mRNA not in accession_index: 
                                try: mRNA_probesets_to_analyze[mRNA].append(pd)
                                except KeyError: mRNA_probesets_to_analyze[mRNA] = [pd]
    k=0
    for probeset in determine_transcript_length_probesets:
        longest_mRNA_ls=[]
        for mRNA in determine_transcript_length_probesets[probeset]:
            try: mRNA_len = ucsc_mRNA_hit_len[mRNA]
            except KeyError: mRNA_len = 1000
            longest_mRNA_ls.append((mRNA_len,mRNA))
        longest_mRNA_ls.sort(); longest_mRNA = longest_mRNA_ls[-1][1]
        probeset_match_FL[probeset] = 10000,longest_mRNA; k+=1
    print k,'probesets added to analysis with more than one mRNA hit based on UCSC location analysis'

    output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/'+species+'_probeset-mRNA_relationships.txt'
    fn=filepath(output_file);data = open(fn,'w')
    filtered_probeset_mRNA_db={}
    for probeset in exon_db:
        longest_mRNA = ''; longest_null_mRNA = ''
        if probeset in probeset_match_ENS: l1,longest_mRNA = probeset_match_ENS[probeset]
        elif probeset in probeset_match_FL: l1,longest_mRNA = probeset_match_FL[probeset]
        elif probeset in probeset_match_EST: l1,longest_mRNA = probeset_match_EST[probeset]
        #if probeset == '2594513': print '$$$$$$$$$$$$',longest_mRNA
        if len(longest_mRNA)>1:
            if probeset in probeset_match_null_ENS: l2,longest_null_mRNA = probeset_match_null_ENS[probeset]
            elif probeset in probeset_match_null_FL:
                l2,longest_null_mRNA = probeset_match_null_FL[probeset]
                #if probeset == '2594513': print '$$$$$$$$$$$$',longest_null_mRNA
                if probeset in probeset_match_null_EST: ###check to see if there is a more optimal EST (more than 20% longer than the FL)
                    l2b,longest_null_mRNAb = probeset_match_null_EST[probeset]
                    if (float(l2b)/float(l2))>1.2:
                        #if longest_null_mRNA == 'X06369' and probeset == '2328854': print '2328854 from',longest_null_mRNA,'to',longest_null_mRNAb
                        l2,longest_null_mRNA = l2b,longest_null_mRNAb
            elif probeset in probeset_match_null_EST: l2,longest_null_mRNA = probeset_match_null_EST[probeset]
            if len(longest_null_mRNA)>0:
                if probeset == '2594513': print '%$$$$$$$$$$$$',longest_mRNA
                if probeset == '2594513': print '%$$$$$$$$$$$$',longest_null_mRNA
                data.write(string.join([probeset,longest_mRNA,longest_null_mRNA],'\t')+'\n')
    data.close()
    
def identifyCorrespondingESTs(probeset_match_db_temp,data):
    accession = ''; unigene_id='nullX'; seq_type = 'null'
    l = SequenceRetrieval(accession,unigene_id,seq_type,0,0)

    splice_event_results = []
    global missing; missing=[]; global found; found=[]; splice_event_results={}
    
    for probeset in probeset_match_db_temp:
            matches = probeset_match_db_temp[probeset]
            matches.sort(); matches.reverse()
            call_db={}
            for (call,mRNA_length,y) in matches:
                try: call_db[call].append(y)
                except KeyError: call_db[call] = [y]
            try:
                longest_null = call_db[0][0]; ens=0
                for y in call_db[0]:
                    if 'ENS' in y.Accession(): longest_null = y; ens=1; break
                if ens==0:
                    for y in call_db[0]:
                        if y.SeqType() == 'full-length': longest_null = y; break
            except KeyError: longest_null = l
            try:
                longest_pos = call_db[1][0]; ens=0
                for y in call_db[1]:
                    if 'ENS' in y.Accession(): longest_pos = y; ens=1; break
                if ens==0:
                    for y in call_db[1]:
                        if y.SeqType() == 'full-length': longest_pos = y; break
            except KeyError: longest_pos = l
            lp = longest_pos; ln = longest_null
            values = string.join([probeset,lp.Accession(),ln.Accession()],'\t')+'\n'
            data.write(values)

def exportSimple(db,filename):
    output_file = 'AltDatabase/'+species+'/SequenceData/output/'+filename
    fn=filepath(output_file)
    data = open(fn,'w')
    for key in db:
        values = string.join([key] + db[key],'\t')+'\n'; data.write(values)
    data.close()


def exportAssociationsSimple(splice_event_results):
    output_file = 'AltDatabase/'+species+'/SequenceData/output/'+species+'_probeset-mRNA_accession.txt'
    fn=filepath(output_file)
    data = open(fn,'w')
    accession_seq_db={}
    for key in splice_event_results:
        y,gene_id,x1,x2,hits1,hits = splice_event_results[key]
        
        values = [y.Probeset(),x1.Accession(),x2.Accession()]
        accession_seq_db[x1.Accession()] = x1.SequenceData()
        accession_seq_db[x2.Accession()] = x2.SequenceData()
        values = string.join(values,'/t')+'\n'; data.write(values)

    data.close()

    if len(accession_seq_db)>0:
        output_file = 'AltDatabase/'+species+'/SequenceData/output/'+species+'_probeset-mRNA_sequence.txt'
        fn=filepath(output_file); data = open(fn,'w')
        for ac in accession_seq_db:
            values = [ac,seq]
            values = string.join(values,'/t')+'\n'; data.write(values)
        data.close()
                
def compareSequences(seq1,seq2,seq1_short):
    if seq1 in seq2: call = 1
    else: call = 0
    if call == 0:
        seq1 = reverse_orientation(seq1)
        if seq1 in seq2: call = 1
    if call == 0:
        if seq1_short in seq2:
            call = 1
    return call

def compareSequencesSimple(seq1,seq2):
    if seq1 in seq2: call = 1
    else: call = 0
    
    if call == 0:
        seq1 = reverse_orientation(seq1)
        if seq1 in seq2: call = 1
    return call


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

def importBestAccessionAlignments(seq_files):
    ensembl_prot_seq_db = importEnsemblProteinSeq(1)
    refseq_mRNA_db,refseq_protein_seq_db = importRefProtRelationships(species)
    ens_transcript_to_protein_db,external_transcript_to_ens_prot_db = importEnsProtRelationships(species)
    uniprot_protein_seq_db,external_transcript_to_uniprot_protein_db,n_terminal_seq,c_terminal_seq = importUniProtSeqeunces(species,'sequence')
    importNCBISequenceAssociations(species)
    
    print len(external_transcript_to_uniprot_protein_db), 'EMBL IDs linked to UniProt'
    relationship_file = '/AltDatabase/'+species+'/SequenceData/output'
    g = GrabFiles(); g.setdirectory(relationship_file)
    relationship_files = g.searchdirectory('probeset-mRNA_relationships')

    relationship_db = {}    
    for filename in relationship_files:
        print 'Reading',filename
        fn=filepath(filename)
        for line in open(fn,'r').xreadlines():
            values, newline= string.split(line,'\n')
            probeset,hit,null = string.split(values,'\t')
            if probeset == '2594513': print '$$$$$$$$$$$-2594513',[hit],[null]
            try: relationship_db[probeset].append((hit,null))
            except KeyError: relationship_db[probeset]= [(hit,null)]
            
    relationship_db = eliminateRedundant(relationship_db)
    x=0; y=0; z=0; relationship_unique_db={}
    for probeset in relationship_db:
        if len(relationship_db[probeset])>1:
            ###Eliminate entries with the least optimal overall matches (no blanks versus one blank)
            keepers=[]; a=[]; b=[]; c=[]
            for pair in relationship_db[probeset]:
                if pair[0]!='' and pair[1]!='': a.append(pair)
                elif pair[0]!='' and pair[1]=='': b.append(pair)
                else: c.append(pair)
            if len(a)>0: keepers = a
            elif len(b)>0: keepers = b
            elif len(c)>0: keepers = c
            if len(keepers)>1:
                y+=1
                ###Compare null AC(s)
                index1=0; ens_index1=0; rs_index1=0; rx_index1=0
                ens1=0; rs1=0; rx1=0
                for pair in keepers:
                    if 'ENS' in pair[1]: ens_index1 = index1; ens1=1
                    if 'NM_' in pair[1]: rs_index1 = index1; rs1=1
                    if 'XM_' in pair[1] or 'XR_' in pair[1] or 'NR_' in pair[1] : rx_index1 = index1; rx1=1
                    index1+=1
                ###Pick the best ones
                index2=0; ens_index2=0; rs_index2=0; rx_index2=0
                ens2=0; rs2=0; rx2=0
                for pair in keepers:
                    if 'ENS' in pair[0]: ens_index2 = index2; ens2=1
                    if 'NM_' in pair[0]: rs_index2 = index2; rs2=1
                    if 'XM_' in pair[0] or 'XR_' in pair[1] or 'NR_' in pair[1] : rx_index2 = index2; rx2=1
                    index2+=1
                if ens2==1 and ens1==1: keep = keepers[ens_index2]
                elif ens2==1: keep = keepers[ens_index2]
                elif rs2==1 and (rs1==1 or ens1==1 or rx1==1): keep = keepers[rs_index2]
                elif rs2==1: keep = keepers[rs_index2]
                elif rx2==1 and (rs1==1 or ens1==1 or rx1==1): keep = keepers[rs_index2]
                elif rx2==1: keep = keepers[rx_index2]
                elif len(a)>0 or len(b)>0: ###thus hit is not blank
                    keep = keepers[0] ###just pick the first entry, if no other optimization present
                    ###COULD IDENTIFY THE mRNA WITH A CORRESPONDING PROTEIN HERE!!!!!!
                else:
                    if ens1==1: keep = keepers[ens_index1]
                    elif rs1==1: keep = keepers[rs_index1]
                    elif rx1==1: keep = keepers[rx_index1]
                    else:
                        keep = keepers[0] ###just pick the first entry, if no other optimization present
                        ###COULD IDENTIFY THE mRNA WITH A CORRESPONDING PROTEIN HERE!!!!!!
                relationship_unique_db[probeset]=keep
            else:
                relationship_unique_db[probeset]=keepers[0]
                x+=1
        else:
            relationship_unique_db[probeset]=relationship_db[probeset][0]
            x+=1
    print x, 'single mRNA matches'
    print y, 'multiple mRNA matches'

    ###Write out the relationship_unique_db
    output_file = 'AltDatabase/'+species+'/SequenceData/output/'+species+'_unique-probeset-mRNAs.txt'
    fn1=filepath(output_file); data = open(fn1,'w')
    for probeset in relationship_unique_db:
        hit_ac,null_ac = relationship_unique_db[probeset]
        data.write(probeset+'\t'+hit_ac+'\t'+null_ac+'\n')
    data.close()           
        
    mRNA_ACs={}; protein_linked_ACs={}; full_length_ACs={}; non_protein_linked_ACs={}; protein_id_not_found={}; x=0; y=0
    miss_assigned = {}
    ###Link each mRNA accession with proteins or store unlinked mRNAs for later translation
    ###Also check to see if mRNAs link to Ensembl transcrpts, even when they didn't align before
    for probeset in relationship_unique_db:
        hit,null = relationship_unique_db[probeset]
        if len(hit)>2 and len(null)>2:
            y+=1
            for ac in relationship_unique_db[probeset]:
                if 'ENS' not in ac and '_' not in ac:
                    mRNA_ACs[ac]=[]
                    if ac in external_transcript_to_ens_prot_db: ###shouldn't occur since these should already be Ensembl transcript ids
                        ###Decided that we should trust the alignments over Ensembl's associations, since unclear associations are made (aligning sequence occurs outside of that list by Ensembl)
                        #full_length_ACs[ac] = external_transcript_to_ens_prot_db[ac][0]  ###Don't think there should be more than 1
                        try: miss_assigned[ac,tuple(external_transcript_to_ens_prot_db[ac])].append(probeset)
                        except KeyError: miss_assigned[ac,tuple(external_transcript_to_ens_prot_db[ac])]= [probeset] ###shouldn't occur since these should already be Ensembl transcript ids
                    if ac in external_transcript_to_uniprot_protein_db: 
                        protein_linked_ACs[ac] = external_transcript_to_uniprot_protein_db[ac]
                    elif ac in mRNA_protein_db: ###two databases from NCBI and UCSC
                        protein_linked_ACs[ac] = mRNA_protein_db[ac]
                        #elif x<40:
                        #print ac; x+=1
                    else: non_protein_linked_ACs[ac] = '' ###For in silico translation
                else:
                    if 'ENS' in ac:
                        if ac in ens_transcript_to_protein_db:
                            ensprot_id = ens_transcript_to_protein_db[ac]
                            full_length_ACs[ac] = ensprot_id   ###relate ensembl transcripts to proteins
                        else:
                            protein_id_not_found[ac] = []
                            non_protein_linked_ACs[ac] = ''
                    else:
                        if ac in refseq_mRNA_db:
                            refprot_id = refseq_mRNA_db[ac]
                            full_length_ACs[ac] = refprot_id ###relate refseq transcripts to proteins
                        elif ac in external_transcript_to_uniprot_protein_db:
                            protein_linked_ACs[ac] = external_transcript_to_uniprot_protein_db[ac]
                        elif ac in mRNA_protein_db: ###two databases from NCBI and UCSC
                            protein_linked_ACs[ac] = mRNA_protein_db[ac]
                        else:
                            protein_id_not_found[ac] = []
                            non_protein_linked_ACs[ac] = '' ###For in silico translation
    probeset1 = '3935913'                
    if probeset1 in protein_linked_ACs: print '$$$$$$$$$$$-',probeset1,'protein_linked_ACs',[protein_linked_ACs[probeset1]]
    if probeset1 in non_protein_linked_ACs: print '$$$$$$$$$$$-',probeset1,'non_protein_linked_ACs',[non_protein_linked_ACs[probeset1]]
    
    print y, 'probesets linked to complete mRNA data'                   
    print len(full_length_ACs), "unique Ensembl/RefSeq accession numbers linked to probesets"
    print len(mRNA_ACs), "unique non-Ensembl accession numbers linked to probesets"
    print len(protein_linked_ACs), "protein_linked_ACs - unique non-Ensembl accession numbers linked to UniProt ACs"
    print len(non_protein_linked_ACs), "non_protein_linked_ACs - mRNA ACs without protein links"
    print len(protein_id_not_found), "protein_id_not_found (should be 0)"
    print len(miss_assigned), "linked indirectly through Ensembl to EnsProt, but not by direct seq match: shouldn't occur"
    #('AB040879', ('ENSP00000347301',)) probesets - ['3579614', '3579612', '3579611'], ENST00000353598
    for i in miss_assigned:
        print i,miss_assigned[i]; break
    protein_linked_ACs = eliminateRedundant(protein_linked_ACs)

    ###Now we have a several different pieces of information about protein ac association to mRNA
    ###1) mRNA ID's linked to only one protein ID (Ensembl or RefSeq)
    ###2) mRNA ID's linked to one or more proteins (see next segment of code)
    ###3) protein IDs from 1 and 2 with sequence in parsed files and without (these need to be added to the in silico list)
    
    x=0; z=[]; linked_protein_ids = {}
    ###In cases where there is more than one linked protein ID, pick the best
    for ac in protein_linked_ACs:
        if len(protein_linked_ACs[ac])>1:
            sp = []; trembl = []
            for i in protein_linked_ACs[ac]:
                if '_' in i: sp.append(i) ###If matching a SwissProt ID
                elif i[0] == 'Q' or i[0] == 'P': trembl.append(i)  ###If matching a TrEMBL ID
            if len(sp)>0: protein_linked_ACs[ac] = sp[0]; linked_protein_ids[sp[0]]=[]; a = sp
            elif len(trembl)>0: protein_linked_ACs[ac] = trembl[0]; linked_protein_ids[trembl[0]]=[]; b= trembl
            else: protein_linked_ACs[ac] = protein_linked_ACs[ac][0]; linked_protein_ids[protein_linked_ACs[ac]]=[]; c = protein_linked_ACs[ac]
        else: protein_linked_ACs[ac] = protein_linked_ACs[ac][0]; linked_protein_ids[protein_linked_ACs[ac]]=[]; d = protein_linked_ACs[ac]

    #print a,b,c,d
    print len(linked_protein_ids), len(protein_linked_ACs)

    ### Grab protein sequences from NCBI that correspond to those IDs linking to mRNAs 
    try: ncbi_protein_seq_db = importNCBIProteinSeq(linked_protein_ids)
    except IOError: ncbi_protein_seq_db={}
    
    ###Output initial protein associations and filtered mRNA to probeset relationships
        
    ###Where there aren't any protein linked mRNA ac's, perform in silico translation
    ###Also, build the mRNA-AC to protein sequence db
    mRNA_protein_seq_db={}
    ###first look at the non-full length sequences then the full-length, to over-ride duplicates with the higher quality data
    for mrna_ac in protein_linked_ACs:
        seq_found = 0
        protein_id = protein_linked_ACs[mrna_ac]
        if protein_id in uniprot_protein_seq_db: sequence = uniprot_protein_seq_db[protein_id]; seq_found = 1
        elif protein_id in ncbi_protein_seq_db: sequence = ncbi_protein_seq_db[protein_id]; seq_found = 1
        else: non_protein_linked_ACs[mrna_ac] = protein_id ###For in silico translation, but record protein ID
        if seq_found == 1: mRNA_protein_seq_db[mrna_ac] = protein_id,sequence    
    for mrna_ac in full_length_ACs: #protein_linked_ACs
        seq_found = 0
        protein_id = full_length_ACs[mrna_ac]
        if protein_id in ensembl_prot_seq_db: sequence = ensembl_prot_seq_db[protein_id]; seq_found = 1
        elif protein_id in ensembl_prot_seq_db: sequence = refseq_protein_seq_db[protein_id]; seq_found = 1
        else: non_protein_linked_ACs[mrna_ac] = protein_id ###For in silico translation, but record protein ID
        if seq_found == 1: mRNA_protein_seq_db[mrna_ac] = protein_id,sequence

    print "Revised number of non_protein_linked_ACs is", len(non_protein_linked_ACs)

    if 'AK025306' in protein_linked_ACs: print '$##$$$$$$$$$-AK025306-protein_linked_ACs',[protein_linked_ACs['AK025306']]
    if 'AK025306' in full_length_ACs: print '$##$$$$$$$$$-AK025306-full_length_ACs',[full_length_ACs['AK025306']]
    if 'AK025306' in non_protein_linked_ACs: print '$##$$$$$$$$$-AK025306-non_protein_linked_ACs',[non_protein_linked_ACs['AK025306']]
    
    ###clear memory of the stored full protein sequence databases
    ensembl_prot_seq_db=[]; refseq_protein_seq_db=[]; ncbi_protein_seq_db=[]; uniprot_protein_seq_db=[]
    ens_transcript_to_protein_db=[]; relationship_db = []
    external_transcript_to_ens_prot_db=[]; refseq_mRNA_db=[]; external_transcript_to_uniprot_protein_db=[]

    ###First obtain the mRNA sequences for all these unliked ac's (collected from the previous section and before), from non_protein_linked_ACs
    output_file = 'AltDatabase/'+species+'/SequenceData/output/'+species+'_filtered-no-protein_mRNA-sequences.txt'
    input_file = output_file
    start_time = time.time()
    fn1=filepath(output_file)
    data = open(fn1,'w')

    if 'AK025306' in non_protein_linked_ACs: print '****$##$$$$$$$$$-AK025306-non_protein_linked_ACs',[non_protein_linked_ACs['AK025306']]
    
    for filename in seq_files:
        print 'Reading',filename
        fn=filepath(filename)
        for line in open(fn,'r').xreadlines():
            values, newline= string.split(line,'\n')
            ac,seq = string.split(values,'\t')
            if 'AK025306' in line: print '&&&&&&&&&&&&&&&&&&&&&',[ac]
            if ac in non_protein_linked_ACs:
                protein_id = non_protein_linked_ACs[ac]
                if 'AK025306' in line: print '&&&&&&&&&&&&&&&&&&&&&',[ac], [seq[:40]],[protein_id]
                try: data.write(ac+'\t'+protein_id+'\t'+seq+'\n')
                except TypeError: print ac, protein_id, seq;kill
    data.close()           
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "Non-protein linked mRNA sequences filtered in %d seconds" % time_diff 
    mRNA_seq_db = importFilteredmRNAs(input_file) ### mRNA sequences from non_protein_linked_ACs

    translation_db = BuildInSilicoTranslations(mRNA_seq_db,n_terminal_seq,c_terminal_seq)

    ###Merge the predicted and downloaded translations
    for mRNA_ac in translation_db:
        mRNA_protein_seq_db[mRNA_ac] = translation_db[mRNA_ac]
        
    print len(mRNA_protein_seq_db), "protein sequences parsed, to be written to file"

    ###To mirror the previous probeset-protein database structure, we need additional exon information data    
    probeset_annotations_file = 'AltDatabase/'+species+'/SequenceData/'+species+'_Ensembl_probesets.txt'
    if array_type == 'exon':
        import ExonSeqModule; exon_db = ExonSeqModule.importSplicingAnnotationDatabase(probeset_annotations_file)
    else:
        import JunctionSeqModule; exon_db = JunctionSeqModule.importSplicingAnnotationDatabaseAndSequence(species,array_type,'probeset')

    """for i in mRNA_protein_seq_db:
        print '----',[i]"""
    if 'AK025306' in mRNA_protein_seq_db:
        print 'AK025306',mRNA_protein_seq_db['AK025306'][0],mRNA_protein_seq_db['AK025306'][1][:10]
    if 'AK025306' in mRNA_protein_seq_db:
        print 'AK025306',mRNA_protein_seq_db['AK025306'][0],mRNA_protein_seq_db['AK025306'][1][:10]
    if '2594513' in relationship_unique_db:
        print '2594513', relationship_unique_db['2594513']
        
    k=0    
    output_file = 'AltDatabase/'+species+'/SequenceData/output/'+species+'/probeset-protein-dbase.txt'
    fn=filepath(output_file); data = open(fn,'w'); filtered_protein_seq_data={}
    title = ['GeneID','Probeset','ExonID','Aligned protein_id','Non-aligned protein_id']
    data.write(string.join(title,'\t')+'\n')
    for probeset in relationship_unique_db:
        hit_ac,null_ac = relationship_unique_db[probeset]
        if probeset == '2594513': print '####$$$$$$$$#########-2594513', [hit_ac], [null_ac]
        hit_prot_ac = ''; null_prot_ac = ''
        if hit_ac in mRNA_protein_seq_db: hit_prot_ac, seq = mRNA_protein_seq_db[hit_ac]#; print [hit_prot_ac], len(seq)
        if null_ac in mRNA_protein_seq_db: null_prot_ac, seq = mRNA_protein_seq_db[null_ac]#; print [null_prot_ac], len(seq)
        if len(hit_prot_ac)>1 and len(null_prot_ac)>1:
            if probeset == '2594513': print '%%%%%%%%%%-2594513', [hit_ac], [null_ac]
            k+=1
            pd = exon_db[probeset]
            if array_type == 'exon': geneid = pd.GeneID()
            else: geneid = pd.ExternalGeneID()
            data.write(geneid+'\t'+probeset+'\t'+pd.ExonID()+'\t'+hit_prot_ac+'\t'+null_prot_ac+'\n')
            filtered_protein_seq_data[hit_ac] = mRNA_protein_seq_db[hit_ac]
            filtered_protein_seq_data[null_ac] = mRNA_protein_seq_db[null_ac]
    data.close()           

    output_file = 'AltDatabase/'+species+'/SequenceData/output/'+species+'/SEQUENCE-protein-dbase.txt'
    fn=filepath(output_file); data = open(fn,'w')
    for mRNA_ac in filtered_protein_seq_data:
        protein_id, sequence = filtered_protein_seq_data[mRNA_ac]
        data.write(protein_id+'\t'+sequence+'\n')
    data.close()
    print k,"probesets with both a positive and negatively aligning optimal protein (or theoretical translation) exported"
    
def importFilteredmRNAs(filename):
    print "Importing filtered RNA sequences"
    fn=filepath(filename); mRNA_seq_db = {}
    for line in open(fn,'r').xreadlines():
        values, newline= string.split(line,'\n')
        ac,protein_id,seq = string.split(values,'\t')
        mRNA_seq_db[ac]=protein_id,seq
    return mRNA_seq_db
           
def BuildInSilicoTranslations(mRNA_db,n_terminal_seq,c_terminal_seq):
    """
    BI517798
    Seq('ATGTGGCCAGGAGACGCCACTGGAGAACATGCTGTTCGCCTCCTTCTACCTTCTGGATTT ...', IUPACUnambiguousDNA())
    213 MWPGDATGEHAVRLLLPSGFYPGFSWQYPGSVAFHPRPQVRDPGQRVPDASGRGRLVVRAGPAHPPGLPLLWEPLAIWGNRMPSHRLPLLPQHVRQHLLPHLHQRRPFPGHCAPGQVPQAPQAPLRTPGLCLPVGGGGCGHGPAAGEPTDRADKHTVGLPAAVPGEGSNMPGVPWQWPSLPVHHQVTCTVIIRSCGRPRVEKALRTRQGHESP

    211 MNGLEVAPPGLITNFSLATAEQCGQETPLENMLFASFYLLDFILALVGNTLALWLFIRDHKSGTPANVFLMHLAVADLSCVLVLPTRLVYHFSGNHWPFGEIACRLTGFLFYLNMYASIYFLTCISADRFLAIVHPVKSLKLRRPLYAHLACAFLWVVVAVAMAPLLVSPQTVQTNTRWVCLQLYREKAPTCLVSLGSGLHFPFITRSRVL
    """    
    translation_db={}
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    from Bio import Translate
    first_time = 0
    for mRNA_AC in mRNA_db:
        if mRNA_AC == 'AK025306': print '@@@@@@@@@@^^^^AK025306...attempting in silico translation'
        temp_protein_list=[]; y=0
        protein_id,sequence = mRNA_db[mRNA_AC]
        if protein_id == '': protein_id = mRNA_AC+'-PEP'
        sequence = string.upper(sequence)
        while (string.find(sequence,'ATG')) != -1:  #while there is still a methionine in the DNA sequence, reload this DNA sequence for translation: find the longest ORF
            x = string.find(sequence,'ATG') #before doing this, need to find the start codon ourselves
            y += x  #maintain a running count of the sequence position
            sequence_met = sequence[x:]  #x gives us the position where the first Met* is.
            seq_type = IUPAC.unambiguous_dna
            #seq_type = IUPAC.ambiguous_dna       
            dna_seq = Seq(sequence_met,seq_type)
            standard_translator = Translate.unambiguous_dna_by_id[1]
            prot_seq = standard_translator.translate_to_stop(dna_seq) #convert the dna to protein sequence
            #prot_seq2 = standard_translator.translate(dna_seq)  #include stop codons
            prot_seq_string = prot_seq.tostring()
            prot_seq_tuple = len(prot_seq_string),y,prot_seq_string,dna_seq  #added DNA sequence to determine which exon we're in later  
            temp_protein_list.append(prot_seq_tuple) #create a list of protein sequences to select the longest one
            sequence = sequence_met[3:]  # sequence_met is the sequence after the first or proceeduring methionine, reset the sequence for the next loop
        if len(temp_protein_list) == 0:
            continue
        else:
            #temp_protein_list = pick_optimal_peptide_seq(temp_protein_list) ###Used a more complex method in the original code to determine the best selection
            temp_protein_list.sort(); temp_protein_list.reverse()
            peptide_len1 = temp_protein_list[0][0]
            prot_seq_string = temp_protein_list[0][2] #extract out the protein sequence string portion of the tuple
            coding_dna_seq_string = temp_protein_list[0][3]
            pos1 = temp_protein_list[0][1] ###position in DNA sequence where the translation starts
            n_term1 = prot_seq_string[:5]; c_term1 = prot_seq_string[-5:]
            ###Check the top protein sequences and see if there are frame differences
            choose = 0
            for protein_data in temp_protein_list[1:]: ###exlcude the first entry
                peptide_len2 = protein_data[0]; pos2= protein_data[1]
                percent_of_top = (float(peptide_len1)/peptide_len2)*100
                if (percent_of_top>70) and (peptide_len2>20):
                    prot_seq_string2 = protein_data[2]
                    n_term2 = prot_seq_string2[:5]; c_term2 = prot_seq_string2[-5:]
                    frame_shift = check4FrameShifts(pos1,pos2)
                    if frame_shift == 'yes':
                        ###determine which prediction is better to use
                        if n_term1 in n_terminal_seq: choose = 1
                        elif n_term2 in n_terminal_seq: choose = 2
                        elif c_term1 in c_terminal_seq: choose = 1
                        elif c_term2 in c_terminal_seq: choose = 2
                        if choose == 2:
                            prot_seq_string = protein_data[2]
                            coding_dna_seq_string = protein_data[3]
                            alt_prot_seq_string = temp_protein_list[0][2]
                            alt_coding_dna_seq_string = temp_protein_list[0][3]
                            if first_time == 0:
                                print mRNA_AC
                                print coding_dna_seq_string
                                print len(prot_seq_string),prot_seq_string
                                print alt_coding_dna_seq_string
                                print len(alt_prot_seq_string),alt_prot_seq_string                                
                                first_time = 1
                            ###write this data out in the future
                    else: break
                else: break ###do not need to keep looking
            dl = (len(prot_seq_string))*3  #figure out what the DNA coding sequence length is
            dna_seq_string_coding_to_end = coding_dna_seq_string.tostring()
            coding_dna_seq_string = dna_seq_string_coding_to_end[0:dl]
            ### Determine if a stop codon is in the sequence or there's a premature end
            coding_diff = len(dna_seq_string_coding_to_end) - len(coding_dna_seq_string)
            if coding_diff > 4: stop_found = 'stop-codon-present'
            else: stop_found = 'stop-codon-absent'
            #print [mRNA_AC],[protein_id],prot_seq_string[0:10]
            if mRNA_AC == 'AK025306': print '*********AK025306',[protein_id],prot_seq_string[0:10]
            translation_db[mRNA_AC] = protein_id,prot_seq_string
    return translation_db

def check4FrameShifts(pos1,pos2):
    pos_diff = abs(pos2 - pos1)
    str_codon_diff = str(float(pos_diff)/3)
    value1,value2 = string.split(str_codon_diff,'.')
    if value2 == '0': frame_shift = 'no'
    else: frame_shift = 'yes'
    return frame_shift

def determineNMD():
        nmd_call = 'NA'
        nmd_length = 0
        x = 0
        y = 0
        n = 0
        if affygene in affygene_end_seq_db2:
            nmd_call = 'none-detected'
            three_coding_exon_reference = affygene_end_seq_db2[affygene]
            if three_coding_exon_reference in exon_annotation_string: #exon_annotation_string is a list of exons
                #print three_coding_exon_reference, last_coding_exon, dog
                if three_coding_exon_reference != last_coding_exon:
                    for exon_info in exon_sequence_list:
                        exon_id = exon_info[0]
                        exon_sequence_data = exon_info[1]
                        iterim_transcript_length = exon_info[2]
                        if three_coding_exon_reference == exon_id:
                            y = 1
                        if x == 1 and y == 0:
                            nmd_length += len(exon_sequence_data)
                            n += 1
                        if last_coding_exon == exon_id:
                            x = 1
                            nmd_length += (int(iterim_transcript_length) - int(three_coding_location))
                            n += 1
        if nmd_length > 50 and n>1:
            nmd_call = 'NMD-suggested'

def importRefProtRelationships(species):
    filename = 'AltDatabase/'+species+'/SequenceData/'+species+'_refseq_mRNA-to-protein.txt'
    fn=filepath(filename); refseq_mRNA_db={};refseq_protein_seq_db={} 
    for line in open(fn,'r').readlines():
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t')
        refseq_protein, refseq_mRNA, prot_seq = t
        refseq_mRNA_db[refseq_mRNA[:-2]] = refseq_protein
        refseq_protein_seq_db[refseq_protein] = prot_seq
    return refseq_mRNA_db,refseq_protein_seq_db

def importEnsProtRelationships(species):
    ###This function is primarily used to extract out EnsTranscript to EnsProt and EMBL to EnsProt relationships for
    ###making protein predictions for exon array probesets
    filename = 'AltDatabase/'+species+'/SequenceData/'+species+'_EnsProt-Annotations.txt'
    fn=filepath(filename); ens_transcript_to_protein_db={}; external_transcript_to_ens_prot_db={}
    for line in open(fn,'r').readlines():
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t')
        try: ens_geneid, ens_transcript_id, ens_peptide_id, refseq_dna_id, refseq_predicted_dna_id, embl = t
        except ValueError: print t;kill
        if len(ens_transcript_id)>0: ens_transcript_to_protein_db[ens_transcript_id] = ens_peptide_id
        dna_ids = [refseq_dna_id,refseq_predicted_dna_id,embl]
        for dna_id in dna_ids:
            if len(dna_id)>0 and len(ens_peptide_id)>0:
                try: external_transcript_to_ens_prot_db[dna_id].append(ens_peptide_id)
                except KeyError: external_transcript_to_ens_prot_db[dna_id] = [ens_peptide_id]
    external_transcript_to_ens_prot_db = eliminateRedundant(external_transcript_to_ens_prot_db)
    return ens_transcript_to_protein_db,external_transcript_to_ens_prot_db

def importUniProtSeqeunces(species,parse_type):
    n_terminal_seq={}; c_terminal_seq={}
    filename = 'AltAnalyze/'+species+'/SequenceData/'+'uniprot_trembl_sequence.txt'           
    fn=filepath(filename); uniprot_protein_seq_db = {}; external_transcript_to_uniprot_protein_db={}
    unigene_ensembl_up = {}
    for line in open(fn,'r').readlines():
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t')
        id=t[0];ac=t[1];ensembls=t[4];seq=t[2];type=t[6];unigenes=t[7];embls=t[9]
        ac=string.split(ac,','); ensembls=string.split(ensembls,','); embls=string.split(embls,','); unigenes=string.split(unigenes,',')
        if parse_type != 'unigene':
            for embl in embls:
                if (len(embl)>0) and type == 'fragment':  ###NOTE: Fragment annotated seem to be the only protein IDs that contain direct references to a specific mRNA rather than promiscous (as opposed to Swissprot and Variant)
                    try: external_transcript_to_uniprot_protein_db[embl].append(id)
                    except KeyError: external_transcript_to_uniprot_protein_db[embl] = [id]
            uniprot_protein_seq_db[id] = seq
            n_terminal_seq[seq[:5]] = []
            c_terminal_seq[seq[-5:]] = []
        else:
            for unigene in unigenes:
                for ensembl in ensembls:
                    if len(unigene)>1 and len(ensembl)>1:
                        try: unigene_ensembl_up[unigene].append(ensembl)
                        except KeyError: unigene_ensembl_up[unigene] = [ensembl]
    if parse_type != 'unigene':
        external_transcript_to_uniprot_protein_db = eliminateRedundant(external_transcript_to_uniprot_protein_db)
        return uniprot_protein_seq_db,external_transcript_to_uniprot_protein_db,n_terminal_seq,c_terminal_seq
    else: return unigene_ensembl_up

def importNCBISequenceAssociations(species):
    filename = 'AltDatabase/SequenceData/gene2accession.txt'   ###Multi-species file        
    fn=filepath(filename)
    for line in open(fn,'r').readlines():
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t')
        try: mRNA_AC=t[3][:-2];prot_AC=t[5][:-2]
        except IndexError: mRNA_AC='';prot_AC=''
        if len(mRNA_AC)>2 and len(prot_AC)>2:
            try: mRNA_protein_db[mRNA_AC].append(prot_AC)
            except KeyError: mRNA_protein_db[mRNA_AC] = [prot_AC]
    print len(mRNA_protein_db), "mRNA to protein accession relationships parsed from EntrezGene data"
    importUCSCSequenceAssociations(species)

def importUCSCSequenceAssociations(species):
    filename = 'AltDatabase/SequenceData/'+species+'/'+'kgXref.txt'           
    fn=filepath(filename)
    for line in open(fn,'r').readlines():
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t')
        try: mRNA_AC=t[1];prot_AC=t[3]  ###Information actually seems largely redundant with what we get from UniProt direclty
        except IndexError: mRNA_AC='';prot_AC=''
        if len(mRNA_AC)>2 and len(prot_AC)>2:
            try: mRNA_protein_db[mRNA_AC].append(prot_AC)
            except KeyError: mRNA_protein_db[mRNA_AC] = [prot_AC]
    print len(mRNA_protein_db), "mRNA to protein accession relationships parsed from EntrezGene and UCSC data"

if __name__ == '__main__':
    a = 'AltMouse'; b = 'exon'
    array_type = b
    h = 'Hs'; m = 'Mm'
    species = h

    match_primers_against_all = 'no' ###Variable used to see if we even want to look at primer sequences (used for AltMouse only)
    parse_from_scratch = 'yes'  ###If using a new Ensembl-probeset.txt file (new genome build)
    kill_after_Ensembl = 'no'

    test = 'yes'
    test_ensembl = 'ENSMUSG00000000902'
    test_unigene = 'Mm.279751'
    re_analyze_mRNA_probeset_matches = 'yes'
    
    print "******Analysis Parameters*******\nChoose Array Type"
    print "1) Affymetrix Exon 1.0 ST Array\n2) Affymetrix AltMouseA Junction Array"; inp = sys.stdin.readline(); inp = inp.strip()
    if inp == '1': array_type = b
    if inp == '2': array_type = a; species = m

    if array_type == b:
        print "Choose Species For Analysis"
        print "1) Human\n2) Mouse"; inp = sys.stdin.readline(); inp = inp.strip()
        if inp == '1': species = h
        if inp == '2': species = m        

    print "******Analysis Parameters*******\nChoose Run Options"
    print "1) Derive probeset-mRNA associations from scratch\n2) Analyze existing"; inp = sys.stdin.readline(); inp = inp.strip()
    if inp == '1': parse_from_scratch = 'yes'
    if inp == '2': parse_from_scratch = 'no'

    print "******Analysis Parameters*******\nChoose Run Options"
    print "1) Run program\n2) Test only"; inp = sys.stdin.readline(); inp = inp.strip()
    if inp == '1': test = 'no'
    if inp == '2': test = 'yes'
    
    print "Analyzing array type:",array_type,species
    print "parse_from_scratch:",parse_from_scratch
      
    global probeset_match_db; probeset_match_db={}

    import_dir = '/AltDatabase'+'/'+species+'/SequenceData'
    filedir = import_dir[1:]+'/'
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for input_file in dir_list:    #loop through each file in the directory to output results
        #if '.seq.all' in input_file: est_seq_annot_file = filedir+input_file
        #if '.data' in input_file: unigene_annot_file = filedir+input_file
        if 'cDNA.fasta' in input_file: trans_file = filedir+input_file
        #if 'Ensembl-Unigene' in input_file: ensembl_unigene_file = filedir+input_file
        if '.probeset.fa' in input_file: probeset_seq_file = filedir+input_file
            
    if parse_from_scratch == 'yes':
        re_analyze_mRNA_probeset_matches = 'yes' ###this is needed when ever you have just built new results from scratch (the previous sort methods don't work well)
        #sys.exit()
        """
        global unigene_ensembl; unigene_ensembl={}
        print "Extracting Unigene relationships..."
        import_ensembl_unigene(ensembl_unigene_file)
        unigene_ensembl_up = importUniProtSeqeunces(species,'unigene')
        
        for unigene in unigene_ensembl_up:
            ensembls = unigene_ensembl_up[unigene]
            for ensembl in ensembls:
                try: unigene_ensembl[unigene].append(ensembl)
                except KeyError: unigene_ensembl[unigene] = [ensembl]
        unigene_ensembl = eliminateRedundant(unigene_ensembl)
        """
        
        print "Extracting microarray probesets sequences..."        
        if array_type == 'AltMouse':
            import JunctionSeqModule; data_type = 'junctions'; probeset_seq_file=''
            splice_event_db = JunctionSeqModule.getParametersAndExecute(probeset_seq_file,array_type,species,data_type)
        elif array_type == 'exon':
            import ExonSeqModule
            splice_event_db = ExonSeqModule.getParametersAndExecute(probeset_seq_file,array_type,species)
        
        importEnsemblTranscriptSequence(trans_file)
        #importESTannotations(est_seq_annot_file)
    global mRNA_protein_db; mRNA_protein_db={}

    import_dir = '/AltDatabase/'+species+'/SequenceData/output/sequences'
    g = GrabFiles(); g.setdirectory(import_dir)
    seq_files = g.searchdirectory('sequences')

    probeset_annotations_file = 'AltDatabase/'+species+'/SequenceData/'+species+'_Ensembl_probesets.txt'
    if array_type == 'exon':
        import ExonSeqModule; exon_db = ExonSeqModule.importSplicingAnnotationDatabase(probeset_annotations_file)
    else:
        import JunctionSeqModule; exon_db = JunctionSeqModule.importSplicingAnnotationDatabaseAndSequence(species,array_type,'probeset')
        
    genes_being_analyzed={}
    for probeset in exon_db:
        gene = exon_db[probeset].GeneID(); genes_being_analyzed[gene] = [gene]

    import FeatureAlignment; global protein_ft_db        
    protein_ft_db,domain_gene_counts = FeatureAlignment.grab_exon_level_feature_calls(species,array_type,genes_being_analyzed)
    if re_analyze_mRNA_probeset_matches == 'yes':
        import_dir = '/AltDatabase/'+species+'/SequenceData/output'
        c = GrabFiles(); c.setdirectory(import_dir)
        align_files = c.searchdirectory('alignments.txt')
        reAnalyzeRNAProbesetMatches(align_files,probeset_seq_file)

    importBestAccessionAlignments(seq_files)
    #resolve issue with probesets matches to EnsProt sequence idirectly.
    #after finding same size, different frames, check N-terminus and C-terminus and always mark as multiple frame predictions in file with NMD
    #generate microRNA alignments
    #export multiple files (mRNA to probeset, protein to probeset, etc.)... may need to correct mRNAs that now get linked to EnsProt indirectly.