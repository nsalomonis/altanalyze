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
import export
import time
from Bio import Entrez

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt"or entry[-4:] == ".tab" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

class GrabFiles:
    def setdirectory(self,value):
        self.data = value
    def display(self):
        print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        file_dirs = getDirectoryFiles(self.data,str(search_term))
        if len(file_dirs)<1: print search_term,'not found',self.data
        return file_dirs
    
def getDirectoryFiles(import_dir, search_term):
    exact_file = ''; exact_file_dirs=[]
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir: exact_file_dirs.append(affy_data_dir)
    return exact_file_dirs

########## End generic file import ##########

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def importSplicingAnnotationDatabase(filename):
    fn=filepath(filename); x=0; probeset_gene_db={}
    for line in open(fn,'rU').xreadlines():             
        probeset_data = cleanUpLine(line)  #remove endline
        if x == 0: x = 1
        else:
            try: t=string.split(probeset_data,'\t'); probeset=t[0]; ens_gene=t[2]
            except IndexError: print t;kill
            probeset_gene_db[probeset]=ens_gene
    print 'Probeset to Ensembl gene data imported'
    return probeset_gene_db

def importProbesetExonCoordinate(species,array_type):
    if array_type == 'exon': filename = "AltDatabase/"+species+"/"+array_type+"/"+species+"_probeset-exon-align-coord.txt"    
    else: filename = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_"+array_type+"_probeset-exon-align-coord.txt"
    fn=filepath(filename); x=0; probeset_exon_coor_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            probeset, exon_start, exon_stop = t
            exon_start = int(exon_start); exon_stop = int(exon_stop)
            #if probeset == '2991483':
            try: probeset_exon_coor_db[probeset].append([exon_start, exon_stop])
            except KeyError: probeset_exon_coor_db[probeset] = [[exon_start, exon_stop]]
    print 'Probeset exon coordiante data imported'
    return probeset_exon_coor_db

def importEnsExonStructureDataSimple(filename,species,ens_transcript_exon_db,ens_gene_transcript_db):
    fn=filepath(filename); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            gene, chr, strand, exon_start, exon_end, ens_exonid, constitutive_exon, ens_transcriptid = t
            exon_start = int(exon_start); exon_end = int(exon_end)
            try: ens_transcript_exon_db[ens_transcriptid].append([exon_start,exon_end])
            except KeyError: ens_transcript_exon_db[ens_transcriptid] = [[exon_start,exon_end]]
            try: ens_gene_transcript_db[gene].append(ens_transcriptid)
            except KeyError: ens_gene_transcript_db[gene] = [ens_transcriptid]
            
    ens_gene_transcript_db = eliminateRedundant(ens_gene_transcript_db)
    print 'Exon/transcript data imported'
    return ens_transcript_exon_db,ens_gene_transcript_db

def eliminateRedundant(database):
    for key in database:
        list = makeUnique(database[key])
        list.sort()
        database[key] = list
    return database

def makeUnique(item):
    db1={}; list1=[]
    for i in item: db1[i]=[]
    for i in db1: list1.append(i)
    list1.sort()
    return list1

def compareExonComposition(species,array_type):
    if array_type == 'exon': ens_probeset_file = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_probesets.txt"    
    else: ens_probeset_file = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_"+array_type+"_probesets.txt"
    probeset_gene_db = importSplicingAnnotationDatabase(ens_probeset_file)

    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'    
    ens_transcript_exon_db,ens_gene_transcript_db = importEnsExonStructureDataSimple(filename,species,{},{})
    ### Add UCSC transcript data to ens_transcript_exon_db and ens_gene_transcript_db
    filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_mrna.txt' ### Use the non-filtered database to propperly analyze exon composition 
    ens_transcript_exon_db,ens_gene_transcript_db = importEnsExonStructureDataSimple(filename,species,ens_transcript_exon_db,ens_gene_transcript_db)
    ### Import exon coordinates that align to probesets
    probeset_exon_coor_db = importProbesetExonCoordinate(species,array_type)

    export_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_all-transcript-matches.txt'                
    data = export.ExportFile(export_file)
    
    ### Identifying isforms containing and not containing the probeset
    probeset_transcript_db={}; match_pairs_missing=0; valid_transcript_pairs={}; ok_transcript_pairs={}; probesets_not_found=0
    for probeset in probeset_exon_coor_db:
        if probeset in probeset_gene_db:
            geneid = probeset_gene_db[probeset]
            transcripts = ens_gene_transcript_db[geneid]
            #print probeset_exon_coor_db[probeset];kill
            exon_match_data_probeset=[]
            for coordinates in probeset_exon_coor_db[probeset]: ### Multiple exons may align
                matching={}; not_matching={}
                for transcript in transcripts:
                    if coordinates in ens_transcript_exon_db[transcript]:
                        #print 'match', transcript
                        other_coordinate_list=[]
                        for other_coordinates in ens_transcript_exon_db[transcript]:
                            if coordinates != other_coordinates: ### Add exon coordinates for all exons that DO NOT aligning to the probeset
                                other_coordinate_list.append(other_coordinates)
                        if len(other_coordinate_list)>0: other_coordinate_list.sort();other_coordinate_list[0][0] = 0; other_coordinate_list[-1][-1] = 0
                        other_coordinate_list = convertListsToTuple(other_coordinate_list)
                        matching[tuple(other_coordinate_list)]=transcript
                    else:
                        #print 'not matched', transcript
                        other_coordinate_list=[]
                        for other_coordinates in ens_transcript_exon_db[transcript]:
                            if coordinates != other_coordinates: ### Add exon coordinates for all exons that DO NOT aligning to the probeset
                                other_coordinate_list.append(other_coordinates)
                        if len(other_coordinate_list)>0: other_coordinate_list.sort();other_coordinate_list[0][0] = 0; other_coordinate_list[-1][-1] = 0
                        other_coordinate_list = convertListsToTuple(other_coordinate_list)
                        not_matching[tuple(other_coordinate_list)]=transcript
                #if probeset == '2429323': print '\n',matching, not_matching
                if len(matching)>0 and len(not_matching)>0:
                    perfect_match_found='no'; exon_match_data=[] ### List is only used if more than a single cassette exon difference between transcripts
                    matching_transcripts=[]; not_matching_transcripts=[]
                    for exon_list in matching:
                        matching_transcript = matching[exon_list]; matching_transcripts.append(matching_transcript)
                        if exon_list in not_matching:
                            not_matching_transcript = not_matching[exon_list]; not_matching_transcripts.append(not_matching_transcript)
                            valid_transcript_pairs[probeset] = matching_transcript, not_matching_transcript
                            perfect_match_found = 'yes'
                            #print probeset,matching_transcript, not_matching_transcript
                            #print ens_transcript_exon_db[matching_transcript],'\n'
                            #print ens_transcript_exon_db[not_matching_transcript]; kill
                        else:
                            unique_exon_count_db={} ### Determine how many exons are common and which are transcript distinct for a pair-wise comp
                            for exon_coor in exon_list:
                                try: unique_exon_count_db[exon_coor] +=1
                                except KeyError: unique_exon_count_db[exon_coor] =1
                            for exon_list in not_matching:
                                not_matching_transcript = not_matching[exon_list]; not_matching_transcripts.append(not_matching_transcript)
                                for exon_coor in exon_list:
                                    try: unique_exon_count_db[exon_coor] +=1
                                    except KeyError: unique_exon_count_db[exon_coor] =1
                                exon_count_db={}
                                for exon_coor in unique_exon_count_db:
                                    num_trans_present = unique_exon_count_db[exon_coor]
                                    try: exon_count_db[num_trans_present]+=1
                                    except KeyError: exon_count_db[num_trans_present]=1
                                try: exon_count_results = [exon_count_db[1],-1*(exon_count_db[2]),matching_transcript,not_matching_transcript]
                                except KeyError: null =[] ###Occurs if no exons are found in common (2)
                                exon_match_data.append(exon_count_results)
                    if perfect_match_found == 'no':
                        exon_match_data.sort()
                        exon_match_data_probeset.append(exon_match_data[0]) ### Since multiple exons may match, ultimately compare the stats of all top hits
                    ###Export transcript comparison sets to an external file for different analyses
                    matching_transcripts = unique.unique(matching_transcripts)
                    not_matching_transcripts = unique.unique(not_matching_transcripts)
                    matching_transcripts=string.join(matching_transcripts,'|')
                    not_matching_transcripts=string.join(not_matching_transcripts,'|')
                    values = string.join([probeset,matching_transcripts,not_matching_transcripts],'\t')+'\n'
                    data.write(values)
                else: match_pairs_missing+=1
            if len(exon_match_data_probeset)>0:
                if probeset not in valid_transcript_pairs:
                    ### Thus, no perfect match found and one or more pairs of top transcript pairs found
                    exon_match_data_probeset.sort()
                    matching_transcript = exon_match_data_probeset[0][-2]
                    not_matching_transcript = exon_match_data_probeset[0][-1]
                    ok_transcript_pairs[probeset] = matching_transcript, not_matching_transcript
                    #if exon_match_data_probeset[0][0] <3:
                    #print probeset, exon_match_data_probeset;kill
        else: probesets_not_found+=1
    print probesets_not_found,'probesets not found'
    print match_pairs_missing,'probesets missing either an alinging or non-alinging transcript'
    print len(valid_transcript_pairs),'probesets with a single exon difference aligning to two isoforms'
    print len(ok_transcript_pairs),'probesets with more than one exon difference aligning to two isoforms'
    data.close()

    export_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_top-transcript-matches.txt'                
    data = export.ExportFile(export_file)    

    for probeset in valid_transcript_pairs:
        matching_transcript,not_matching_transcript = valid_transcript_pairs[probeset]
        values = string.join([probeset,matching_transcript,not_matching_transcript],'\t')+'\n'
        data.write(values)

    for probeset in ok_transcript_pairs:
        matching_transcript,not_matching_transcript = ok_transcript_pairs[probeset]
        values = string.join([probeset,matching_transcript,not_matching_transcript],'\t')+'\n'
        data.write(values)
    data.close()
    
def compareProteinComposition(species,array_type,translate):
    probeset_transcript_db,unique_ens_transcripts,unique_transcripts = importProbesetTranscriptMatches(species,array_type)
    if translate == 'yes':
        translateRNAs(unique_transcripts,unique_ens_transcripts,'fetch')
        
def importProteinSequences(species,just_get_ids):
    transcript_protein_seq_db={}
    import_dir = '/AltDatabase/'+species+'/SequenceData/output/sequences' ### Multi-species fiel
    g = GrabFiles(); g.setdirectory(import_dir)
    seq_files = g.searchdirectory('Transcript-')
    for seq_file in seq_files:
        fn=filepath(seq_file)
        for line in open(fn,'rU').xreadlines():             
            probeset_data = cleanUpLine(line)  #remove endline
            mRNA_AC,protein_AC,protein_seq = string.split(probeset_data,'\t')
            if just_get_ids == 'yes': transcript_protein_seq_db[mRNA_AC]=[]
            else: transcript_protein_seq_db[mRNA_AC] = protein_AC,protein_seq
    return transcript_protein_seq_db

def importProbesetTranscriptMatches(species,array_type):
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_all-transcript-matches.txt'
    fn=filepath(filename); probeset_transcript_db={}; unique_transcripts={}; unique_ens_transcripts={}
    for line in open(fn,'rU').xreadlines():             
        probeset_data = cleanUpLine(line)  #remove endline
        try: probeset,match_transcripts,not_match_transcripts = string.split(probeset_data,'\t')
        except IndexError: print t;kill
        match_transcripts = string.split(match_transcripts,'|')
        not_match_transcripts = string.split(not_match_transcripts,'|')
        ### Multiple transcript comparison sets can exist, depending on how many unique exon coordiantes the probeset aligns to
        try: probeset_transcript_db[probeset].append([match_transcripts,not_match_transcripts])
        except KeyError : probeset_transcript_db[probeset] = [[match_transcripts,not_match_transcripts]]

        ###Store unique transcripts        
        for transcript in match_transcripts:
            if 'ENS' in transcript: unique_ens_transcripts[transcript]=[]
            else: unique_transcripts[transcript]=[]
        for transcript in not_match_transcripts:
            if 'ENS' in transcript: unique_ens_transcripts[transcript]=[]
            else: unique_transcripts[transcript]=[]
    print 'len(unique_ens_transcripts)',len(unique_ens_transcripts)
    print 'len(unique_transcripts)',len(unique_transcripts)
    return probeset_transcript_db,unique_ens_transcripts,unique_transcripts

def translateRNAs(unique_transcripts,unique_ens_transcripts,analysis_type):
    if analysis_type == 'local':
        ### Get protein ACs for UCSC transcripts if provided by UCSC
        mRNA_protein_db,missing_protein_ACs_UCSC = importUCSCSequenceAssociations(species,unique_transcripts)
        ### For missing_protein_ACs, check to see if they are in UniProt. If so, export the protein sequence
        missing_protein_ACs_UniProt = importUniProtSeqeunces(species,mRNA_protein_db,missing_protein_ACs_UCSC)
        ### For transcripts with protein ACs, see if we can find sequence from NCBI
        #missing_protein_ACs_NCBI = importNCBIProteinSeq(mRNA_protein_db)
        ### Combine missing_protein_ACs_NCBI and missing_protein_ACs_UniProt
        #for mRNA_AC in missing_protein_ACs_NCBI: missing_protein_ACs_UniProt[mRNA_AC] = missing_protein_ACs_NCBI[mRNA_AC]
        ### Import mRNA sequences for mRNA ACs with no associated protein sequence and submit for in silico translation
        missing_protein_ACs_inSilico = importUCSCSequences(missing_protein_ACs_NCBI)
        
    if analysis_type != 'local': missing_protein_ACs_UniProt = importUniProtSeqeunces(species,{},{})
        
    ### Export Ensembl protein sequences for matching isoforms and identify transcripts without protein seqeunce
    #ensembl_protein_seq_db, missing_protein_ACs_Ensembl = importEnsemblProteinSeqData(species,unique_ens_transcripts)
    ### Import Ensembl mRNA sequences for mRNA ACs with no associated protein sequence and submit for in silico translation
    #missing_Ens_protein_ACs_inSilico = importEnsemblTranscriptSequence(missing_protein_ACs_Ensembl)
    
    if analysis_type != 'local': 
        ac_list = []
        for ac in unique_transcripts: ac_list.append(ac)
        postEntrez(ac_list[:50],'nucleotide',1)

    ###Import protein sequences
    just_get_ids = 'yes'
    transcript_protein_db = importProteinSequences(species,just_get_ids)
    print len(unique_ens_transcripts)+len(unique_transcripts), 'examined transcripts.'
    print len(transcript_protein_db),'transcripts with associated protein sequence.'
    
    missing_protein_data=[]
    for ac in unique_transcripts:
        if ac not in transcript_protein_db: missing_protein_data.append(ac)

    ###Search NCBI using the esearch not efetch function to get valid GIs (some ACs make efetch crash)    
    try: missing_gi_list = searchEntrez(missing_protein_data,'nucleotide')
    except IOError: missing_gi_list = searchEntrez(missing_protein_data,'nucleotide')

    fetchSeq(missing_gi_list,'nucleotide',2)
    transcript_protein_seq_db = importProteinSequences(species,'no')
    print len(unique_ens_transcripts)+len(unique_transcripts), 'examined transcripts.'
    print len(transcript_protein_seq_db),'transcripts with associated protein sequence.'
    
def searchEntrez(accession_list,bio_type):
    start_time = time.time()
    Entrez.email = "nsalomonis@gmail.com" # Always tell NCBI who you are
    index=0; gi_list=[]
    while index<len(accession_list)+20:
        try: new_accession_list = accession_list[index:index+20]
        except IndexError: new_accession_list = accession_list[index:]
        if len(new_accession_list)<1: break
        search_handle = Entrez.esearch(db=bio_type,term=string.join(new_accession_list,','))
        search_results = Entrez.read(search_handle)
        gi_list += search_results["IdList"]
        index+=20
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "finished in %s seconds" % time_diff
    return gi_list

def postEntrez(accession_list,bio_type,version):
    output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-Protein_sequences_'+str(version)+'.txt'
    fn2=filepath(output_file)
    datar = open(fn2,'w')

    print len(accession_list), "mRNA Accession numbers submitted to eUTILs."   
    start_time = time.time()
    Entrez.email = "nsalomonis@gmail.com" # Always tell NCBI who you are
    index=0

    search_handle = Entrez.epost(db=bio_type,term=string.join(accession_list,','))
    search_results = Entrez.read(search_handle)
    search_handle.close()
    gi_list = search_results["IdList"]
    print len(accession_list);kill

        
def fetchSeq(accession_list,bio_type,version):
    output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-Protein_sequences_'+str(version)+'.txt'
    fn2=filepath(output_file)
    datar = open(fn2,'w')

    print len(accession_list), "mRNA Accession numbers submitted to eUTILs."   
    start_time = time.time()
    Entrez.email = "nsalomonis@gmail.com" # Always tell NCBI who you are
    index=0
    while index<len(accession_list)+200:
        try: new_accession_list = accession_list[index:index+200]
        except IndexError: new_accession_list = accession_list[index:]
        if len(new_accession_list)<1: break
        ### epost doesn't concatenate everything into a URL  where esearch 
        
        try:
            search_handle = Entrez.efetch(db=bio_type,id=string.join(new_accession_list,','),retmode="xml") #,usehistory="y"
            search_results = Entrez.read(search_handle)
        except Exception:
            try: ### Make sure it's not due to an internet connection issue
                search_handle = Entrez.efetch(db=bio_type,id=string.join(new_accession_list,','),retmode="xml") #,usehistory="y"
                search_results = Entrez.read(search_handle)
            except Exception: index+=200; continue
        for a in search_results: ### list of dictionaries
            accession = a['GBSeq_primary-accession']
            if bio_type == 'nucleotide':
                mRNA_seq=''; proceed = 'no'
                ### Get translated protein sequence if available
                try:
                    for entry in a["GBSeq_feature-table"]:
                        for key in entry: ### Key's in dictionary
                            if key == 'GBFeature_quals': ### composed of a list of dictionaries
                                for i in entry[key]:
                                    if i['GBQualifier_name'] == 'protein_id': protein_id = i['GBQualifier_value']
                                    if i['GBQualifier_name'] == 'translation': protein_sequence = i['GBQualifier_value']; proceed = 'yes'
                except KeyError: null = []
                alt_seq='no'
                if proceed == 'yes':
                    if protein_sequence[0] != 'M': proceed = 'no'; alt_seq = 'yes'
                    else:
                        values = [accession,protein_id,protein_sequence]; values = string.join(values,'\t')+'\n'; datar.write(values)                    
                if proceed == 'no': ### Translate mRNA seq to protein
                    mRNA_seq = a['GBSeq_sequence']
                    mRNA_db={}; mRNA_db[accession] = '',mRNA_seq
                    translation_db = BuildInSilicoTranslations(mRNA_db)
                    for mRNA_AC in translation_db: ### Export in silico protein predictions
                        protein_id, protein_sequence = translation_db[mRNA_AC]
                        values = [mRNA_AC,protein_id,protein_sequence]; values = string.join(values,'\t')+'\n'; datar.write(values)
                    if len(translation_db)==0 and alt_seq == 'yes': ### If no protein sequence starting with an "M" found, write the listed seq
                        values = [accession,protein_id,protein_sequence]; values = string.join(values,'\t')+'\n'; datar.write(values)  
            else: protein_sequence = a['GBSeq_sequence']
        index+=200
        """print 'accession',accession    
        print 'protein_id',protein_id
        print 'protein_sequence',protein_sequence
        print 'mRNA_seq',mRNA_seq,'\n' """

    end_time = time.time(); time_diff = int(end_time-start_time)
    print "finished in %s seconds" % time_diff
    datar.close()
    
def importUCSCSequences(missing_protein_ACs):
    start_time = time.time()
    filename = 'AltDatabase/'+species+'/SequenceData/mrna.fa'
    output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-InSilicoProt_sequences.txt'
    fn2=filepath(output_file)
    datar = open(fn2,'w')    
    
    print "Begining generic fasta import of",filename
    #'>gnl|ENS|Mm#S10859962 Mus musculus 12 days embryo spinal ganglion cDNA, RIKEN full-length enriched library, clone:D130006G06 product:unclassifiable, full insert sequence /gb=AK051143 /gi=26094349 /ens=Mm.1 /len=2289']
    #'ATCGTGGTGTGCCCAGCTCTTCCAAGGACTGCTGCGCTTCGGGGCCCAGGTGAGTCCCGC'
    fn=filepath(filename); sequence = '|'; translated_mRNAs={}
    for line in open(fn,'r').xreadlines():
        try: data, newline= string.split(line,'\n')
        except ValueError: continue
        if len(data)>0:
            if data[0] != '#':
                try:
                    if data[0] == '>':
                        if len(sequence) > 1:
                            if accession in missing_protein_ACs:
                                ### Perform in silico translation
                                mRNA_db = {}; mRNA_db[accession] = '',sequence[1:]
                                translation_db = BuildInSilicoTranslations(mRNA_db)
                                for mRNA_AC in translation_db: ### Export in silico protein predictions
                                    protein_id, protein_seq = translation_db[mRNA_AC]
                                    values = [mRNA_AC,protein_id,protein_seq]; values = string.join(values,'\t')+'\n'; datar.write(values)
                                    translated_mRNAs[mRNA_AC]=[]
                        ### Parse new line
                        values = string.split(data,' '); accession = values[0][1:]; sequence = '|'; continue
                except IndexError: null = []       
                try:
                    if data[0] != '>': sequence = sequence + data
                except IndexError: print kill; continue
    datar.close()
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "UCSC mRNA sequences analyzed in %d seconds" % time_diff

    missing_protein_ACs_inSilico={}
    for mRNA_AC in missing_protein_ACs:
        if mRNA_AC not in translated_mRNAs: missing_protein_ACs_inSilico[mRNA_AC]=[]
    print len(missing_protein_ACs_inSilico), 'mRNAs without mRNA sequence that were NOT translated into protein'
    return missing_protein_ACs_inSilico

def importEnsemblTranscriptSequence(missing_protein_ACs):
    filename = 'AltDatabase/'+species+'/SequenceData/'+species+'_ensembl_cDNA.fasta.txt'    
    start_time = time.time()
    output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-EnsInSilicoProt_sequences.txt'
    fn2=filepath(output_file)
    datar = open(fn2,'w')
    
    print "Begining generic fasta import of",filename
    fn=filepath(filename); translated_mRNAs={}; sequence = ''
    for line in open(fn,'r').xreadlines():
        try: data, newline= string.split(line,'\n')
        except ValueError: continue
        try:
            if data[0] == '>':
                    if len(sequence) > 0:
                        if transid in missing_protein_ACs:
                            ### Perform in silico translation
                            mRNA_db = {}; mRNA_db[transid] = '',sequence[1:]
                            translation_db = BuildInSilicoTranslations(mRNA_db)
                            for mRNA_AC in translation_db: ### Export in silico protein predictions
                                protein_id, protein_seq = translation_db[mRNA_AC]
                                values = [mRNA_AC,protein_id,protein_seq]; values = string.join(values,'\t')+'\n'; datar.write(values)
                                translated_mRNAs[mRNA_AC]=[]
                    ### Parse new line                      
                    t= string.split(data,'|'); sequence=''
                    try: ensembl_id,chr,strand,transid,prot_id = t
                    except ValueError: ensembl_id,chr,strand,transid = t
        except IndexError: continue
        try:
            if data[0] != '>': sequence = sequence + data
        except IndexError:  continue
    datar.close()
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "Ensembl transcript sequences analyzed in %d seconds" % time_diff

    missing_protein_ACs_inSilico={}
    for mRNA_AC in missing_protein_ACs:
        if mRNA_AC not in translated_mRNAs: missing_protein_ACs_inSilico[mRNA_AC]=[]
    print len(missing_protein_ACs_inSilico), 'Ensembl mRNAs without mRNA sequence that were NOT translated into protein'    
    
def importNCBIProteinSeq(mRNA_protein_db):
    export_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-NCBIProt_sequences.txt'
    export_data = export.ExportFile(export_file)
    
    ### Reverse the values and key of the dictionary
    protein_linked_ACs = {}
    for mRNA_AC in mRNA_protein_db:
        protein_AC = mRNA_protein_db[mRNA_AC]
        try: protein_linked_ACs[protein_AC].append(mRNA_AC)
        except KeyError: protein_linked_ACs[protein_AC] = [mRNA_AC]
        
    import_dir = '/AltDatabase/SequenceData' ### Multi-species fiel
    g = GrabFiles(); g.setdirectory(import_dir)
    seq_files = g.searchdirectory('fsa_aa')
    for filename in seq_files:
        fn=filepath(filename); ncbi_prot_seq_db = {}
        sequence = ''
        for line in open(fn,'rU').xreadlines():
            try: data, newline= string.split(line,'\n')
            except ValueError: continue
            try:
                if data[0] == '>':
                    if len(sequence) > 0:
                        if len(prot_id)>0 and prot_id in protein_linked_ACs:
                            mRNA_ACs = protein_linked_ACs[protein_AC]
                            for mRNA_AC in mRNA_ACs:
                                values = string.join([mRNA_AC,prot_id,sequence],'\t')+'\n'
                                export_data.write(values)
                                ncbi_prot_seq_db[mRNA_AC] = [] ###occurs once in the file
                        t= string.split(data,'|'); prot_id = t[3][:-2]; sequence = ''
                    else: sequence = ''; t= string.split(data,'|'); prot_id = t[3][:-2] ###Occurs for the first entry  
            except IndexError: continue
            try:
                if data[0] != '>': sequence = sequence + data
            except IndexError:  continue
            
        export_data.close()
        if len(prot_id)>0 and prot_id in protein_linked_ACs:
            mRNA_ACs = protein_linked_ACs[protein_AC]
            for mRNA_AC in mRNA_ACs:
                values = string.join([mRNA_AC,prot_id,sequence],'\t')+'\n'
                export_data.write(values)            
                ncbi_prot_seq_db[mRNA_AC] = []  ###Need this for the last entry

    missing_protein_ACs={}
    for mRNA_AC in mRNA_protein_db:
        if mRNA_AC not in ncbi_prot_seq_db: missing_protein_ACs[mRNA_AC]=[]
            
    print len(ncbi_prot_seq_db), "genbank protein ACs associated with sequence, out of", len(mRNA_protein_db)
    print len(missing_protein_ACs), 'missing protein ACs for associated NCBI mRNAs and', len(ncbi_prot_seq_db), 'found.'  
    return missing_protein_ACs

def importEnsemblProteinSeqData(species,unique_ens_transcripts):
    import FeatureAlignment
    protein_relationship_file,protein_seq_file,protein_feature_file = FeatureAlignment.getEnsemblRelationshipDirs(species)
    ens_transcript_protein_db = FeatureAlignment.importEnsemblRelationships(protein_relationship_file,'transcript')

    unique_ens_proteins = {}
    for transcript in ens_transcript_protein_db:
        if transcript in unique_ens_transcripts:
            protein_id = ens_transcript_protein_db[transcript]
            unique_ens_proteins[protein_id] = transcript
    ensembl_protein_seq_db = importEnsemblProtSeq(protein_seq_file,unique_ens_proteins)

    transcript_with_prot_seq = {}
    for protein_id in ensembl_protein_seq_db:
        if protein_id in unique_ens_proteins:
            transcript = unique_ens_proteins[protein_id]
            transcript_with_prot_seq[transcript]=[]
            
    missing_ens_proteins={}
    for transcript in unique_ens_transcripts:
        if transcript not in transcript_with_prot_seq: missing_ens_proteins[transcript]=[]

    print len(ensembl_protein_seq_db),'Ensembl transcripts linked to protein sequence and',len(missing_ens_proteins), 'transcripts missing protein sequence.'
    return ensembl_protein_seq_db, missing_ens_proteins
            
def importEnsemblProtSeq(filename,unique_ens_proteins):
    export_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-EnsProt_sequences.txt'
    export_data = export.ExportFile(export_file)
    
    fn=filepath(filename); ensembl_protein_seq_db={}; x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x=1 ###Don't extract the headers
        else: 
            ensembl_prot, aa_start, aa_stop, protein_sequence = string.split(data,'\t')
            if ensembl_prot in unique_ens_proteins:
                mRNA_AC = unique_ens_proteins[ensembl_prot]
                values = string.join([mRNA_AC,ensembl_prot,protein_sequence],'\t')+'\n'
                export_data.write(values)
                ensembl_protein_seq_db[ensembl_prot] =  []
    export_data.close()
    return ensembl_protein_seq_db

def importUCSCSequenceAssociations(species,transcripts_to_analyze):
    filename = 'AltDatabase/'+species+'/SequenceData/kgXref.txt'           
    fn=filepath(filename); mRNA_protein_db={}
    for line in open(fn,'r').readlines():
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t')
        try: mRNA_AC=t[1];prot_AC=t[2]  ###Information actually seems largely redundant with what we get from UniProt direclty
        except IndexError: mRNA_AC='';prot_AC=''
        #if mRNA_AC == 'BC152405': print [prot_AC];kill
        if len(mRNA_AC)>2 and len(prot_AC)>2:
            if mRNA_AC in transcripts_to_analyze:
                try: mRNA_protein_db[prot_AC].append(mRNA_AC)
                except KeyError: mRNA_protein_db[prot_AC] = [mRNA_AC]
    print len(mRNA_protein_db), "mRNA to protein ACs parsed from UCSC data"

    missing_protein_ACs={} ### These will require in Silico translation
    for mRNA_AC in transcripts_to_analyze:
        if mRNA_AC not in mRNA_protein_db: missing_protein_ACs[mRNA_AC]=[]

    print len(missing_protein_ACs), 'missing protein ACs for associated UCSC mRNAs and', len(mRNA_protein_db), 'found.'  
    return mRNA_protein_db,missing_protein_ACs

def importUniProtSeqeunces(species,transcripts_with_uniprots,transcripts_to_analyze):
    global n_terminal_seq; global c_terminal_seq
    n_terminal_seq={}; c_terminal_seq={}

    export_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-UniProt_sequences.txt'
    export_data = export.ExportFile(export_file)
    filename = 'AltDatabase/'+species+'/SequenceData/'+'uniprot_trembl_sequence.txt'           
    fn=filepath(filename); transcript_to_uniprot={}
    unigene_ensembl_up = {}
    for line in open(fn,'r').readlines():
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t')
        id=t[0];ac=t[1];ensembls=t[4];seq=t[2];type=t[6];unigenes=t[7];embls=t[9]
        ac=string.split(ac,','); embls=string.split(embls,',') #; ensembls=string.split(ensembls,','); unigenes=string.split(unigenes,',')
        if type != 'swissprot':
            ### Note: These associations are based on of UCSC, which in some cases don't look correct: see AY429540	and Q75N08 from the KgXref file.
            ### Possibly exclude
            ac = ac[0]
            if ac in transcripts_with_uniprots:
                mRNA_ACs = transcripts_with_uniprots[ac]
                for mRNA_AC in mRNA_ACs:
                    transcript_to_uniprot[mRNA_AC] = []
                    values = string.join([mRNA_AC,ac,seq],'\t')+'\n'; export_data.write(values)
        for embl in embls:
            proceed = 'no'
            if (len(embl)>0) and type == 'fragment':  ###NOTE: Fragment annotated seem to be the only protein IDs that contain direct references to a specific mRNA rather than promiscous (as opposed to Swissprot and Variant)
                if embl in transcripts_to_analyze: proceed = 'yes'
                elif embl in transcripts_with_uniprots: proceed = 'yes'
            if proceed == 'yes':
                if embl not in transcript_to_uniprot:
                    transcript_to_uniprot[embl] = []
                    values = string.join([embl,id,seq],'\t')+'\n'; export_data.write(values)
        n_terminal_seq[seq[:5]] = []
        c_terminal_seq[seq[-5:]] = []
    export_data.close()

    missing_protein_ACs={}
    for mRNA_AC in transcripts_to_analyze:
        if mRNA_AC not in transcript_to_uniprot: missing_protein_ACs[mRNA_AC]=[]

    for protein_AC in transcripts_with_uniprots:
        mRNA_ACs = transcripts_with_uniprots[protein_AC]
        for mRNA_AC in mRNA_ACs:
            if mRNA_AC not in transcript_to_uniprot: missing_protein_ACs[mRNA_AC]=[]

    print len(missing_protein_ACs), 'missing protein ACs for associated UniProt mRNAs and', len(transcript_to_uniprot), 'found.'      
    return missing_protein_ACs
           
def BuildInSilicoTranslations(mRNA_db):
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
    first_time = 1
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

def convertListsToTuple(list_of_lists):
    list_of_tuples=[]
    for ls in list_of_lists:
        list_of_tuples.append(tuple(ls))
    return list_of_tuples

if __name__ == '__main__':
    species = 'Hs'; array_type = 'exon'; translate='yes'
    #compareExonComposition(species,array_type)
    compareProteinComposition(species,array_type,translate)