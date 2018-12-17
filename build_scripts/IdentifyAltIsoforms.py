#!/usr/local/bin/python2.6

###IdentifyAltIsoforms
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

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os.path
import unique
import export
import time
import copy
from build_scripts import FeatureAlignment; reload(FeatureAlignment)
from Bio import Entrez
from build_scripts import ExonAnalyze_module
from build_scripts import mRNASeqAlign
import traceback
import requests

requests.packages.urllib3.disable_warnings()
import ssl
try: _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    # Legacy Python that doesn't verify HTTPS certificates by default
    pass
else:
    # Handle target environment that doesn't support HTTPS verification
    ssl._create_default_https_context = _create_unverified_https_context
    
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if (entry[-4:] == ".txt"or entry[-4:] == ".tab" or entry[-4:] == ".csv" or '.fa' in entry) and '.gz' not in entry and '.zip' not in entry: dir_list2.append(entry)
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
    try: dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    except Exception:
        print unique.filepath(import_dir)
        export.createDirPath(unique.filepath(import_dir[1:]))
        dir_list = read_directory(import_dir)
    dir_list.sort() ### Get the latest files
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

def importSplicingAnnotationDatabase(array_type,species,import_coordinates):

    if array_type == 'exon' or array_type == 'gene': filename = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_probesets.txt"    
    elif array_type == 'junction': filename = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_"+array_type+"_probesets.txt"
    elif array_type == 'RNASeq': filename = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_exons.txt"
    
    fn=filepath(filename); x=0; probeset_gene_db={}; probeset_coordinate_db={}
    for line in open(fn,'rU').xreadlines():             
        probeset_data = cleanUpLine(line)  #remove endline
        if x == 0: x = 1
        else:
            t=string.split(probeset_data,'\t'); probeset=t[0]; exon_id=t[1]; ens_gene=t[2]; start=int(t[6]); stop=int(t[7])
            if array_type != 'RNASeq': 
                if ':' in probeset: probeset = string.split(probeset,':')[1]
            start_stop = [start,stop]; start_stop.sort(); start,stop = start_stop; proceed = 'yes'
            #if probeset[-2] != '|': ### These are the 5' and 3' exons for a splice-junction (don't need to analyze)
            if test == 'yes':
                if ens_gene not in test_genes: proceed = 'no'
            if proceed == 'yes':
                probeset_gene_db[probeset]=ens_gene,exon_id
                if import_coordinates == 'yes': probeset_coordinate_db[probeset] = start,stop
    print 'Probeset to Ensembl gene data imported'
    if import_coordinates == 'yes': return probeset_gene_db,probeset_coordinate_db
    else: return probeset_gene_db

def importAltMouseJunctionDatabase(species,array_type):
    ### Import AffyGene to Ensembl associations (e.g., AltMouse array)    
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'-Ensembl.txt'
    fn=filepath(filename); array_ens_db={}; x = 0
    for line in open(fn,'r').xreadlines():
        data, newline = string.split(line,'\n'); t = string.split(data,'\t')
        if x==0: x=1
        else: 
            array_gene,ens_gene = t
            array_ens_db[array_gene]=ens_gene
            #try: array_ens_db[array_gene].append(ens_gene)
            #except KeyError: array_ens_db[array_gene]=[ens_gene]
    
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'_junction-comparisons.txt'
    fn=filepath(filename); probeset_gene_db={}; x = 0
    for line in open(fn,'r').xreadlines():
        data, newline = string.split(line,'\n'); t = string.split(data,'\t')
        if x==0: x=1
        else: 
            array_gene,probeset1,probeset2,critical_exons = t #; critical_exons = string.split(critical_exons,'|')
            if array_gene in array_ens_db:
                ensembl_gene_ids = array_ens_db[array_gene]
            else: ensembl_gene_ids=[]
            probeset_gene_db[probeset1+'|'+probeset2] = array_gene,critical_exons
            probeset_gene_db[probeset1] = array_gene,critical_exons
            probeset_gene_db[probeset2] = array_gene,critical_exons
    return probeset_gene_db

def importJunctionDatabase(species,array_type):
    if array_type == 'junction': filename = 'AltDatabase/' + species + '/'+array_type+'/'+ species + '_junction_comps_updated.txt'
    else: filename = 'AltDatabase/' + species + '/'+array_type+'/'+ species + '_junction_comps.txt'
    fn=filepath(filename); probeset_gene_db={}; x=0
    for line in open(fn,'rU').xreadlines():
        if x==0: x=1
        else:
            data = cleanUpLine(line)
            gene,critical_exons,excl_junction,incl_junction,excl_junction_probeset,incl_junction_probeset,source = string.split(data,'\t')
            probeset_gene_db[incl_junction_probeset+'|'+excl_junction_probeset] = gene,critical_exons
            probeset_gene_db[excl_junction_probeset] = gene,critical_exons
            probeset_gene_db[incl_junction_probeset] = gene,critical_exons
    return probeset_gene_db
            
def importEnsExonStructureDataSimple(filename,species,ens_transcript_exon_db,ens_gene_transcript_db,ens_gene_exon_db,filter_transcripts):
    fn=filepath(filename); x=0
    print fn
    """
    if 'Ensembl' not in filename:
        original_ensembl_transcripts={} ### Keep track of the original ensembl transcripts
        for i in ens_transcript_exon_db: original_ensembl_transcripts[i]=[]"""
        
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            gene, chr, strand, exon_start, exon_end, ens_exonid, constitutive_exon, ens_transcriptid = t
            exon_start = int(exon_start); exon_end = int(exon_end); proceed = 'yes'
            if test == 'yes':
                if gene not in test_genes: proceed = 'no'
            if (ens_transcriptid in filter_transcripts or len(filter_transcripts)==0) and proceed == 'yes':
                ### Only import transcripts that are in the set to analyze
                try: ens_transcript_exon_db[ens_transcriptid][exon_start,exon_end]=[]
                except KeyError: ens_transcript_exon_db[ens_transcriptid] = db = {(exon_start,exon_end):[]}
                try: ens_gene_transcript_db[gene][ens_transcriptid]=[]
                except KeyError: ens_gene_transcript_db[gene] = db = {ens_transcriptid:[]}
                try: ens_gene_exon_db[gene][exon_start,exon_end]=[]
                except KeyError: ens_gene_exon_db[gene] = db = {(exon_start,exon_end):[]}

    """
    ### Some transcripts will have the same coordinates - get rid of these (not implemented yet)
    if 'Ensembl' not in filename:
        for gene in ens_gene_transcript_db:
            exon_structure_db={}
            for transcript in ens_gene_transcript_db[gene]:
                ls=[]
                for coord in ens_transcript_exon_db[transcript]: ls.append(coord)
                ls.sort()
                try: exon_structure_db[tuple(ls)].append(transcript)
                except Exception: exon_structure_db[tuple(ls)] = [transcript]"""
    #ens_gene_transcript_db = eliminateRedundant(ens_gene_transcript_db)
    #ens_gene_exon_db = eliminateRedundant(ens_gene_exon_db)
    #if 'ENSG00000240173' in ens_gene_exon_db: print 'here'
    print 'Exon/transcript data imported for %s transcripts' % len(ens_transcript_exon_db), len(ens_gene_exon_db)
    return ens_transcript_exon_db,ens_gene_transcript_db,ens_gene_exon_db

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

def getProbesetExonCoordinates(probeset_coordinate_db,probeset_gene_db,ens_gene_exon_db):
    probeset_exon_coor_db={}
    for probeset in probeset_gene_db:
        gene = probeset_gene_db[probeset][0]
        start,stop = probeset_coordinate_db[probeset]
        coor_list = ens_gene_exon_db[gene]
        proceed = 'yes'
        if test == 'yes':
            if gene not in test_genes: proceed = 'no'
        if proceed == 'yes':
            for (t_start,t_stop) in coor_list:
                if ((start >= t_start) and (start <= t_stop)) and ((stop >= t_start) and (stop <= t_stop)):
                    ### Thus the probeset is completely in this exon
                    try: probeset_exon_coor_db[probeset].append((t_start,t_stop))
                    except KeyError: probeset_exon_coor_db[probeset]=[(t_start,t_stop)]
                    #if '20.8' in probeset: print (t_start,t_stop),start,stop
    #sys.exit()
    return probeset_exon_coor_db

def compareExonComposition(species,array_type):
    probeset_gene_db,probeset_coordinate_db = importSplicingAnnotationDatabase(array_type, species,'yes')
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'    
    ens_transcript_exon_db,ens_gene_transcript_db,ens_gene_exon_db = importEnsExonStructureDataSimple(filename,species,{},{},{},{})
    ### Add UCSC transcript data to ens_transcript_exon_db and ens_gene_transcript_db
    try:
        filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_COMPLETE-mrna.txt' ### Use the non-filtered database to propperly analyze exon composition 
        ens_transcript_exon_db,ens_gene_transcript_db,ens_gene_exon_db = importEnsExonStructureDataSimple(filename,species,ens_transcript_exon_db,ens_gene_transcript_db,ens_gene_exon_db,{})
    except Exception:pass
    
    ### Derive probeset to exon associations De Novo
    probeset_exon_coor_db = getProbesetExonCoordinates(probeset_coordinate_db,probeset_gene_db,ens_gene_exon_db)

    if (array_type == 'junction' or array_type == 'RNASeq') and data_type == 'exon':
        export_file = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/'+species+'_all-transcript-matches.txt'  
    else:
        export_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_all-transcript-matches.txt'                
    data = export.ExportFile(export_file)

    ### Identifying isforms containing and not containing the probeset
    probeset_transcript_db={}; match_pairs_missing=0; valid_transcript_pairs={}; ok_transcript_pairs={}; probesets_not_found=0

    print len(probeset_exon_coor_db)
    print len(ens_transcript_exon_db)
    print len(ens_gene_transcript_db)
    print len(ens_gene_exon_db)
    for probeset in probeset_exon_coor_db:
        geneid = probeset_gene_db[probeset][0]
        transcripts = ens_gene_transcript_db[geneid]
        matching_transcripts=[]; not_matching_transcripts=[]; matching={}; not_matching={}
        for coordinates in probeset_exon_coor_db[probeset]: ### Multiple exons may align
            for transcript in transcripts:
                ### Build a cursory list of matching and non-matching transcripts

                if coordinates in ens_transcript_exon_db[transcript]:
                    matching_transcripts.append(transcript)
                else:
                    not_matching_transcripts.append(transcript)
                
                    
        ### Filter large non-matching transcript sets to facilate processing
        not_matching_transcripts = mRNASeqAlign.filterNullMatch(not_matching_transcripts,matching_transcripts)

        ### Re-analyze the filtered list for exon content
        transcripts = unique.unique(not_matching_transcripts+matching_transcripts)
        matching_transcripts=[]; not_matching_transcripts=[]
        for coordinates in probeset_exon_coor_db[probeset]: ### Multiple exons may align
            for transcript in transcripts:
                if coordinates in ens_transcript_exon_db[transcript]:
                    matching_transcripts.append(transcript)
                    other_coordinate_list=[]
                    for other_coordinates in ens_transcript_exon_db[transcript]:
                        if coordinates != other_coordinates: ### Add exon coordinates for all exons that DO NOT aligning to the probeset
                            other_coordinate_list.append(other_coordinates)
                    if len(other_coordinate_list)>1:
                        other_coordinate_list.sort()
                        ### Instead of replacing the values in place, we need to do this (otherwise the original object will be modified)
                        other_coordinate_list = [[0,other_coordinate_list[0][1]]]+other_coordinate_list[1:-1]+[[other_coordinate_list[-1][0],0]]
                    #other_coordinate_list[0][0] = 0; other_coordinate_list[-1][-1] = 0
                    other_coordinate_list = convertListsToTuple(other_coordinate_list)
                    matching[tuple(other_coordinate_list)]=transcript
                else:
                    not_matching_transcripts.append(transcript)
                    other_coordinate_list=[]
                    for other_coordinates in ens_transcript_exon_db[transcript]:
                        if coordinates != other_coordinates: ### Add exon coordinates for all exons that DO NOT aligning to the probeset
                            other_coordinate_list.append(other_coordinates)
                    if len(other_coordinate_list)>1:
                        other_coordinate_list.sort()
                        other_coordinate_list = [[0,other_coordinate_list[0][1]]]+other_coordinate_list[1:-1]+[[other_coordinate_list[-1][0],0]]
                    other_coordinate_list = convertListsToTuple(other_coordinate_list)
                    not_matching[tuple(other_coordinate_list)]=transcript
                        
        #print '\n',len(matching_transcripts), len(not_matching_transcripts);kill
        ### Can't have transcripts in matching and not matching
        not_matching_transcripts2=[]; not_matching2={} 
        for transcript in not_matching_transcripts:
            if transcript not in matching_transcripts: not_matching_transcripts2.append(transcript)
        for coord in not_matching:
            transcript = not_matching[coord]
            if transcript  in not_matching_transcripts2: not_matching2[coord] = not_matching[coord]
        not_matching = not_matching2; not_matching_transcripts = not_matching_transcripts2
            
        #if probeset == '3431530': print '3431530a', matching_transcripts,not_matching_transcripts
        if len(matching)>0 and len(not_matching)>0:
            perfect_match_found='no'; exon_match_data=[] ### List is only used if more than a single cassette exon difference between transcripts
            for exon_list in matching:
                matching_transcript = matching[exon_list]
                if exon_list in not_matching:
                    not_matching_transcript = not_matching[exon_list]
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
                        not_matching_transcript = not_matching[exon_list]
                        for exon_coor in exon_list:
                            try: unique_exon_count_db[exon_coor] +=1
                            except KeyError: unique_exon_count_db[exon_coor] =1
                        exon_count_db={}
                        for exon_coor in unique_exon_count_db:
                            num_trans_present = unique_exon_count_db[exon_coor]
                            try: exon_count_db[num_trans_present]+=1
                            except KeyError: exon_count_db[num_trans_present]=1
                        try:
                            exon_count_results = [exon_count_db[1],-1*(exon_count_db[2]),matching_transcript,not_matching_transcript]
                            exon_match_data.append(exon_count_results)
                        except KeyError:
                            null =[] ###Occurs if no exons are found in common (2)
                        #if probeset == '3431530': print '3431530b', exon_count_results
            if perfect_match_found == 'no' and len(exon_match_data)>0:
                exon_match_data.sort()
                matching_transcript = exon_match_data[0][-2]
                not_matching_transcript = exon_match_data[0][-1]
                ok_transcript_pairs[probeset] = matching_transcript, not_matching_transcript
                #if probeset == '3431530': print '3431530', matching_transcript, not_matching_transcript
            else: match_pairs_missing+=1
                
            ###Export transcript comparison sets to an external file for different analyses
            matching_transcripts = unique.unique(matching_transcripts)
            not_matching_transcripts = unique.unique(not_matching_transcripts)
            matching_transcripts=string.join(matching_transcripts,'|')
            not_matching_transcripts=string.join(not_matching_transcripts,'|')
            values = string.join([probeset,matching_transcripts,not_matching_transcripts],'\t')+'\n'
            data.write(values)
    print match_pairs_missing,'probesets missing either an alinging or non-alinging transcript'
    print len(valid_transcript_pairs),'probesets with a single exon difference aligning to two isoforms'
    print len(ok_transcript_pairs),'probesets with more than one exon difference aligning to two isoforms'

    data.close()

    if (array_type == 'junction' or array_type == 'RNASeq') and data_type != 'null': 
        export_file = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/'+species+'_top-transcript-matches.txt'
    else:
        export_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_top-transcript-matches.txt'                
    export_data = export.ExportFile(export_file)    

    for probeset in valid_transcript_pairs:
        matching_transcript,not_matching_transcript = valid_transcript_pairs[probeset]
        values = string.join([probeset,matching_transcript,not_matching_transcript],'\t')+'\n'
        export_data.write(values)

    for probeset in ok_transcript_pairs:
        matching_transcript,not_matching_transcript = ok_transcript_pairs[probeset]
        values = string.join([probeset,matching_transcript,not_matching_transcript],'\t')+'\n'
        export_data.write(values)
        #if probeset == '3431530': print '3431530d', matching_transcript,not_matching_transcript
    export_data.close()

def compareExonCompositionJunctionArray(species,array_type):
    print 'Finding optimal isoform matches for splicing events.'
    ###Import sequence aligned transcript associations for individual or probeset-pairs. Pairs provided for match-match
    probeset_transcript_db,unique_ens_transcripts,unique_transcripts,all_transcripts = importProbesetTranscriptMatches(species,array_type,'yes')
    
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'    
    ens_transcript_exon_db,ens_gene_transcript_db,ens_gene_exon_db = importEnsExonStructureDataSimple(filename,species,{},{},{},all_transcripts)
    ### Add UCSC transcript data to ens_transcript_exon_db and ens_gene_transcript_db
    try:
        filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_COMPLETE-mrna.txt' ### Use the non-filtered database to propperly analyze exon composition 
        ens_transcript_exon_db,ens_gene_transcript_db,ens_gene_exon_db = importEnsExonStructureDataSimple(filename,species,ens_transcript_exon_db,ens_gene_transcript_db,ens_gene_exon_db,all_transcripts)
    except Exception: pass
    
    """ Used to prioritize isoform pairs for further analysis. This file is re-written as AltAnalyze
    builds it's database, hence, protein coding potential information accumulates as it runs."""
    seq_files, protein_coding_isoforms = importProteinSequences(species,just_get_ids=False,just_get_length=True)

    def isoformPrioritization(exon_match_data):
        """ Prioritize based on the composition and whether the isoform is an Ensembl protein coding isoform.
        Introduced 11/26/2017(2.1.1)."""

        exon_match_data.sort(); exon_match_data.reverse()
        coding_match=[]
        ens_match=[]
        alt_match=[]
        for (unique_count,common_count,match,non) in exon_match_data:
            ### First pair is the best match
            alt_match.append([match,non])
            if 'ENS' in match and 'ENS' in non:
                try: ens_match_count = protein_coding_isoforms[match]
                except: ens_match_count = 10000
                try: ens_non_count = protein_coding_isoforms[non]
                except: ens_non_count = -10000
                coding_diff = abs(ens_match_count-ens_non_count)
                count_ratio = (coding_diff*1.00)/max(ens_match_count,ens_non_count)
                ens_match.append([count_ratio,coding_diff,ens_match_count,match,non])
            if match in protein_coding_isoforms and non in protein_coding_isoforms:
                ### Rank by AA differing
                match_count = protein_coding_isoforms[match]
                non_count = protein_coding_isoforms[non]
                coding_diff = abs(match_count-non_count)
                count_ratio = (coding_diff*1.00)/max(match_count,non_count)
                coding_match.append([count_ratio,coding_diff,match_count,match,non])
        coding_match.sort(); ens_match.sort()
        if len(coding_match)>0: ### Prioritize minimal coding differences   
            count_ratio,diff,count,final_match,final_non = coding_match[0]
        elif len(ens_match)>0:
            count_ratio,diff,count,final_match,final_non = ens_match[0]
        else:
            final_match,final_non = alt_match[0]
        """
        if len(coding_match)>0 and len(ens_match)>0:
            if len(coding_match)!= len(ens_match):
                print coding_match
                print ens_match
                print alt_match
                if 'ENS' not in final_match:
                    print final_match,final_non
                    sys.exit()"""
        return final_match,final_non

    print len(probeset_transcript_db), "probesets with multiple pairs of matching-matching or matching-null transcripts."
    ### Identifying isoforms containing and not containing the probeset
    global transcripts_not_found; marker = 5000; increment = 5000
    match_pairs_missing=0; ok_transcript_pairs={}; transcripts_not_found={}
    for probesets in probeset_transcript_db: 
        exon_match_data=[]; start_time = time.time()
        ### Examine all valid matches and non-matches
        match_transcripts,not_match_transcripts = probeset_transcript_db[probesets]
        match_transcripts = unique.unique(match_transcripts)
        not_match_transcripts = unique.unique(not_match_transcripts)
        for matching_transcript in match_transcripts:
            if matching_transcript in ens_transcript_exon_db: ### If in the Ensembl or UCSC transcript database
                transcript_match_coord = ens_transcript_exon_db[matching_transcript] ### Get the coordinates
                for not_matching_transcript in not_match_transcripts:
                    if not_matching_transcript in ens_transcript_exon_db:
                        transcript_null_coord = ens_transcript_exon_db[not_matching_transcript]
                        unique_exon_count_db={} ### Determine how many exons are common and which are transcript distinct for a pair-wise comp
                        for exon_coor in transcript_match_coord:
                            try: unique_exon_count_db[tuple(exon_coor)] +=1
                            except KeyError: unique_exon_count_db[tuple(exon_coor)] =1
                        for exon_coor in transcript_null_coord:
                            try: unique_exon_count_db[tuple(exon_coor)] +=1
                            except KeyError: unique_exon_count_db[tuple(exon_coor)] =1
                        exon_count_db={}
                        for exon_coor in unique_exon_count_db:
                            num_trans_present = unique_exon_count_db[exon_coor]
                            try: exon_count_db[num_trans_present]+=1
                            except KeyError: exon_count_db[num_trans_present]=1
                            
                        """ 11/26/2017(2.1.1): The below required that for two isoforms (matching/non-matching to an event),
                        that at least one exon was unique to a transcript and that at least one exon was in common to both.
                        This is not always the case, for example with complete intron retention, the coordinates in the
                        intron retained transcript will be distinct from the non-retained. The requirement was eleminated."""
                        try: iso_unique_exon_count = exon_count_db[1]
                        except: iso_unique_exon_count = 0
                        try: iso_common_exon_count = exon_count_db[2]
                        except: iso_common_exon_count = 0
                        exon_count_results = [iso_unique_exon_count,-1*iso_common_exon_count,matching_transcript,not_matching_transcript]
                        exon_match_data.append(exon_count_results)
                        #if probeset == '3431530': print '3431530b', exon_count_results
            else:
                ### Should rarely occur since accession would be missing from the databases
                transcripts_not_found[matching_transcript]=[]
                
        """ Compare each isoform matching-nonmatching pair in terms of composition """
        try:
            matching_transcript,not_matching_transcript = isoformPrioritization(exon_match_data)
            ok_transcript_pairs[probesets] = matching_transcript, not_matching_transcript
        except Exception: pass
        end_time = time.time()
        if len(ok_transcript_pairs) == marker:
            marker+= increment; print '*',
                
    print match_pairs_missing,'junctions missing either an alinging or non-alinging transcript'
    print len(transcripts_not_found),'transcripts not found'
    #print len(ok_transcript_pairs),'probesets with more than one exon difference aligning to two isoforms'

    if (array_type == 'junction' or array_type == 'RNASeq') and data_type != 'null': 
        export_file = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/'+species+'_top-transcript-matches.txt'  
    else: export_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_top-transcript-matches.txt'                
    export_data = export.ExportFile(export_file)    

    for probeset in ok_transcript_pairs:
        matching_transcript,not_matching_transcript = ok_transcript_pairs[probeset]
        values = string.join([probeset,matching_transcript,not_matching_transcript],'\t')+'\n'
        export_data.write(values)
        #if probeset == '3431530': print '3431530d', matching_transcript,not_matching_transcript
    export_data.close()
    
def compareProteinComposition(species,array_type,translate,compare_all_features):
    """ Objectives 11/26/2017(2.1.1):
    1) Import ALL isoform matches for each junction or reciprocal junction pair
    2) Obtain protein translations for ALL imported isoforms
    3) With the protein isoform information, compare exon composition and protein coding length to pick two representative isoforms
    4) Find protein compositional differences between these pairs
    
    This is a modification of the previous workflow which returned fewer hits that were likely less carefully selected."""
    
    probeset_transcript_db,unique_ens_transcripts,unique_transcripts,all_transcripts = importProbesetTranscriptMatches(species,array_type,'yes')
    all_transcripts=[]

    """if translate == 'yes': ### Used if we want to re-derive all transcript-protein sequences - set to yes1 when we want to disable this option
        transcript_protein_seq_db = translateRNAs(unique_transcripts,unique_ens_transcripts,'fetch')
    else: """
    
    ### Translate all mRNA seqeunces to proteins where possible
    transcript_protein_seq_db = translateRNAs(unique_transcripts,unique_ens_transcripts,'fetch_new')
    
    ### Find the best isoform pairs for each splicing event
    compareExonCompositionJunctionArray(species,array_type)
    
    ### Re-import isoform associations for the top events exported from the above function (top versus all isoforms)
    probeset_transcript_db,unique_ens_transcripts,unique_transcripts,all_transcripts = importProbesetTranscriptMatches(species,array_type,'no')
    
    probeset_protein_db,protein_seq_db = convertTranscriptToProteinAssociations(probeset_transcript_db,transcript_protein_seq_db,'no')
    transcript_protein_seq_db=[]

    global protein_ft_db
    if array_type == 'exon' or array_type == 'gene' or data_type == 'exon':
        probeset_gene_db = importSplicingAnnotationDatabase(array_type,species,'no')
    elif (array_type == 'junction' or array_type == 'RNASeq'):
        probeset_gene_db = importJunctionDatabase(species,array_type)
    elif array_type == 'AltMouse':
        probeset_gene_db = importAltMouseJunctionDatabase(species,array_type)

    genes_being_analyzed={} ### Need to tell 'grab_exon_level_feature_calls' what genes to analyze
    for probeset in probeset_gene_db:
        if probeset in probeset_protein_db: ### Filter for genes that have probesets with matching and non-matching proteins
            gene = probeset_gene_db[probeset][0]; genes_being_analyzed[gene]=[gene]
    protein_ft_db,domain_gene_counts = FeatureAlignment.grab_exon_level_feature_calls(species,array_type,genes_being_analyzed)
    
    if compare_all_features == 'yes': ### Used when comparing all possible PROTEIN pairs to find minimal domain changes      
        compareProteinFeaturesForPairwiseComps(probeset_protein_db,protein_seq_db,probeset_gene_db,species,array_type)

    ExonAnalyze_module.identifyAltIsoformsProteinComp(probeset_gene_db,species,array_type,protein_ft_db,compare_all_features,data_type)

def remoteTranslateRNAs(Species,unique_transcripts,unique_ens_transcripts,analysis_type):
    global species
    species = Species
    translateRNAs(unique_transcripts,unique_ens_transcripts,analysis_type)
    
def translateRNAs(unique_transcripts,unique_ens_transcripts,analysis_type):
    if analysis_type == 'local':
        ### Get protein ACs for UCSC transcripts if provided by UCSC (NOT CURRENTLY USED BY THE PROGRAM!!!!)
        mRNA_protein_db,missing_protein_ACs_UCSC = importUCSCSequenceAssociations(species,unique_transcripts)
        ### For missing_protein_ACs, check to see if they are in UniProt. If so, export the protein sequence
        try: missing_protein_ACs_UniProt = importUniProtSeqeunces(species,mRNA_protein_db,missing_protein_ACs_UCSC)
        except Exception: null=[]
        ### For transcripts with protein ACs, see if we can find sequence from NCBI
        #missing_protein_ACs_NCBI = importNCBIProteinSeq(mRNA_protein_db)
        ### Combine missing_protein_ACs_NCBI and missing_protein_ACs_UniProt
        #for mRNA_AC in missing_protein_ACs_NCBI: missing_protein_ACs_UniProt[mRNA_AC] = missing_protein_ACs_NCBI[mRNA_AC]
        ### Import mRNA sequences for mRNA ACs with no associated protein sequence and submit for in silico translation
        missing_protein_ACs_inSilico = importUCSCSequences(missing_protein_ACs_NCBI)
        
    else:
        try: missing_protein_ACs_UniProt = importUniProtSeqeunces(species,{},{})
        except Exception, e:
            print e
            null=[]

    ### Export Ensembl protein sequences for matching isoforms and identify transcripts without protein seqeunce
    ensembl_protein_seq_db, missing_protein_ACs_Ensembl = importEnsemblProteinSeqData(species,unique_ens_transcripts)
    ### Import Ensembl mRNA sequences for mRNA ACs with no associated protein sequence and submit for in silico translation
    missing_Ens_protein_ACs_inSilico = importEnsemblTranscriptSequence(missing_protein_ACs_Ensembl)
    
    if analysis_type == 'fetch': 
        ac_list = []
        for ac in unique_transcripts: ac_list.append(ac)

        try: ### Get rid of the second file which will not be immediately over-written and reading before regenerating
            output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-Protein_sequences_'+str(2)+'.txt'
            fn=filepath(output_file); os.remove(fn)
        except Exception: null=[]
        try: fetchSeq(ac_list,'nucleotide',1)
        except Exception:
            print traceback.format_exc(),'\n'
            print 'WARNING!!!!! NCBI webservice connectivity failed due to the above error!!!!!\n'
        
    ###Import protein sequences
    seq_files, transcript_protein_db = importProteinSequences(species,just_get_ids=True)
    print len(unique_ens_transcripts)+len(unique_transcripts), 'examined transcripts.'
    print len(transcript_protein_db),'transcripts with associated protein sequence.'
    
    missing_protein_data=[]
    #transcript_protein_db={} ### Use to override existing mRNA-protein annotation files
    for ac in unique_transcripts:
        if ac not in transcript_protein_db: missing_protein_data.append(ac)

    ###Search NCBI using the esearch not efetch function to get valid GIs (some ACs make efetch crash)    
    try: missing_gi_list = searchEntrez(missing_protein_data,'nucleotide')
    except Exception,e:
        print 'Exception encountered:',e
        try: missing_gi_list = searchEntrez(missing_protein_data,'nucleotide')
        except Exception:
            print 'Exception encountered:',e
            try: missing_gi_list = searchEntrez(missing_protein_data,'nucleotide')
            except Exception:
                print traceback.format_exc(),'\n'
                print 'WARNING!!!!! NCBI webservice connectivity failed due to the above error!!!!!\n'

    try: fetchSeq(missing_gi_list,'nucleotide',len(seq_files)-2)
    except Exception:
        print traceback.format_exc(),'\n'
        print 'WARNING!!!!! NCBI webservice connectivity failed due to the above error!!!!!\n'
    seq_files, transcript_protein_seq_db = importProteinSequences(species)
    print len(unique_ens_transcripts)+len(unique_transcripts), 'examined transcripts.'
    print len(transcript_protein_seq_db),'transcripts with associated protein sequence.'
    return transcript_protein_seq_db

def convertTranscriptToProteinAssociations(probeset_transcript_db,transcript_protein_seq_db,compare_all_features):
    ### Convert probeset-transcript match db to probeset-protein match db
    probeset_protein_db={}; compared_protein_ids={}
    for probeset in probeset_transcript_db:
        match_transcripts,not_match_transcripts = probeset_transcript_db[probeset]
        match_proteins=[]; not_match_proteins=[]
        for transcript in match_transcripts:
            if transcript in transcript_protein_seq_db:
                protein_id = transcript_protein_seq_db[transcript][0]; match_proteins.append(protein_id); compared_protein_ids[protein_id]=[]
        for transcript in not_match_transcripts:
            if transcript in transcript_protein_seq_db:
                protein_id = transcript_protein_seq_db[transcript][0]; not_match_proteins.append(protein_id); compared_protein_ids[protein_id]=[]
        if len(match_proteins)>0 and len(not_match_proteins)>0:
            probeset_protein_db[probeset]=[match_proteins,not_match_proteins]

    protein_seq_db={}
    for transcript in transcript_protein_seq_db:
        protein_id, protein_seq = transcript_protein_seq_db[transcript]
        if protein_id in compared_protein_ids: protein_seq_db[protein_id] = [protein_seq]
    transcript_protein_seq_db=[]

    if compare_all_features == 'no': ### If yes, all pairwise comps need to still be examined
        title_row = 'Probeset\tAligned protein_id\tNon-aligned protein_id'
        if (array_type == 'junction' or array_type == 'RNASeq') and data_type != 'null':
            export_file = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/probeset-protein-dbase_exoncomp.txt'
        else: export_file = 'AltDatabase/'+species+'/'+array_type+'/probeset-protein-dbase_exoncomp.txt'
        exportSimple(probeset_protein_db,export_file,title_row)
        if (array_type == 'junction' or array_type == 'RNASeq') and data_type != 'null':
            export_file = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/SEQUENCE-protein-dbase_exoncomp.txt' 
        else: export_file = 'AltDatabase/'+species+'/'+array_type+'/SEQUENCE-protein-dbase_exoncomp.txt'     
        exportSimple(protein_seq_db,export_file,'')
    
    return probeset_protein_db,protein_seq_db
            
def importProteinSequences(species,just_get_ids=False,just_get_length=False):
    transcript_protein_seq_db={}
    import_dir = '/AltDatabase/'+species+'/SequenceData/output/sequences' ### Multi-species fiel
    g = GrabFiles(); g.setdirectory(import_dir)
    seq_files = g.searchdirectory('Transcript-')
    for seq_file in seq_files:
        fn=filepath(seq_file)
        for line in open(fn,'rU').xreadlines():             
            probeset_data = cleanUpLine(line)  #remove endline
            mRNA_AC,protein_AC,protein_seq = string.split(probeset_data,'\t')
            if just_get_ids or just_get_length:
                if just_get_length:
                    transcript_protein_seq_db[mRNA_AC]=len(protein_seq)
                else:
                    transcript_protein_seq_db[mRNA_AC]=protein_AC
            else: transcript_protein_seq_db[mRNA_AC] = protein_AC,protein_seq
    return seq_files, transcript_protein_seq_db
        
def importCDScoordinates(species):
    """ Read in the CDS locations for Ensembl proteins, NCBI proteins and in silico derived proteins and then export to a single file """
    cds_location_db={}
    errors = {}
    import_dir = '/AltDatabase/'+species+'/SequenceData/output/details' 
    import_dir2 = '/AltDatabase/ensembl/'+species 
    g = GrabFiles(); g.setdirectory(import_dir)
    g2 = GrabFiles(); g2.setdirectory(import_dir2)
    seq_files = g.searchdirectory('Transcript-')
    seq_files += g2.searchdirectory('EnsemblTranscriptCDSPositions')

    output_file = 'AltDatabase/ensembl/'+species +'/AllTranscriptCDSPositions.txt'
    dataw = export.ExportFile(output_file)
    
    for seq_file in seq_files:
        fn=filepath(seq_file)
        for line in open(fn,'rU').xreadlines():             
            line_data = cleanUpLine(line)  #remove endline
            line_data = string.replace(line_data,'>','') ### occurs for some weird entries
            line_data = string.replace(line_data,'<','') ### occurs for some weird entries
            if 'join' in line_data:
                ### Occures when the cds entry looks like this: AL137661 - join(203..1267,1267..2187) or 'AL137661\tjoin(203\t1267,1267\t2187)'
                t = string.split(line_data,'\t')
                mRNA_AC = t[0]; start = t[1][5:]; stop = t[-1][:-1]
            else:
                line_data = string.replace(line_data,'complement','') ### occurs for some weird entries
                line_data = string.replace(line_data,')','') ### occurs for some weird entries
                line_data = string.replace(line_data,'(','') ### occurs for some weird entries
                try:
                    mRNA_AC,start,stop = string.split(line_data,'\t')
                    try:
                        cds_location_db[mRNA_AC] = int(start),int(stop)
                        dataw.write(string.join([mRNA_AC,start,stop],'\t')+'\n')
                    except Exception:
                        errors[line_data]=[]
                        #print line_data;sys.exit()
                except Exception:
                    errors[line_data]=[]
                    #print line_data;sys.exit()
    print len(errors),'errors...out of',len(cds_location_db)
    dataw.close()
    return cds_location_db

def importProbesetTranscriptMatches(species,array_type,compare_all_features):
    if compare_all_features == 'yes': ### Used when comparing all possible PROTEIN pairs    
        if (array_type == 'junction' or array_type == 'RNASeq') and data_type != 'null':
            filename = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/'+species+'_all-transcript-matches.txt'
        else: filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_all-transcript-matches.txt'
    else: ### Used after comparing all possible TRANSCRIPT STRUCTURE pairs
        if (array_type == 'junction' or array_type == 'RNASeq') and data_type != 'null': 
            filename = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/'+species+'_top-transcript-matches.txt'
        else: filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_top-transcript-matches.txt'

    print 'Imported:',filename
    fn=filepath(filename); probeset_transcript_db={}; unique_transcripts={}; unique_ens_transcripts={}; all_transcripts={}
    for line in open(fn,'rU').xreadlines():             
        probeset_data = cleanUpLine(line)  #remove endline
        try: probeset,match_transcripts,not_match_transcripts = string.split(probeset_data,'\t')
        except IndexError: print t;kill
        #if probeset == '2991483':
        match_transcripts = string.split(match_transcripts,'|')
        not_match_transcripts = string.split(not_match_transcripts,'|')
        ### Multiple transcript comparison sets can exist, depending on how many unique exon coordiantes the probeset aligns to
        probeset_transcript_db[probeset] = match_transcripts,not_match_transcripts

        consider_all_ensembl = 'no'
        if array_type == 'RNASeq': ### Needed for Ensembl gene IDs that don't have 'ENS' in them.
            if 'ENS' not in probeset: consider_all_ensembl = 'yes'
        ###Store unique transcripts        
        for transcript in match_transcripts:
            if 'ENS' in transcript or consider_all_ensembl == 'yes': unique_ens_transcripts[transcript]=[]
            else: unique_transcripts[transcript]=[]
            if consider_all_ensembl == 'yes': unique_transcripts[transcript]=[] ### This redundant, but is the most unbiased way to examine all non-Ensembls as well
            all_transcripts[transcript]=[]
        for transcript in not_match_transcripts:
            if 'ENS' in transcript or consider_all_ensembl == 'yes': unique_ens_transcripts[transcript]=[]
            else: unique_transcripts[transcript]=[]
            if consider_all_ensembl == 'yes': unique_transcripts[transcript]=[] ### This redundant, but is the most unbiased way to examine all non-Ensembls as well
            all_transcripts[transcript]=[]
    print 'len(unique_ens_transcripts)',len(unique_ens_transcripts)
    print 'len(unique_transcripts)',len(unique_transcripts)
    return probeset_transcript_db,unique_ens_transcripts,unique_transcripts,all_transcripts

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
    
def fetchSeq(accession_list,bio_type,version):
    output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-Protein_sequences_'+str(version)+'.txt'
    datar = export.ExportFile(output_file)
    
    output_file = 'AltDatabase/'+species+'/SequenceData/output/details/Transcript-Protein_sequences_'+str(version)+'.txt'
    datad = export.ExportFile(output_file)

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
                mRNA_seq=''; proceed = 'no'; get_location = 'no'
                ### Get translated protein sequence if available
                try:
                    for entry in a["GBSeq_feature-table"]:
                        for key in entry: ### Key's in dictionary
                            if key == 'GBFeature_quals': ### composed of a list of dictionaries
                                for i in entry[key]:
                                    if i['GBQualifier_name'] == 'protein_id': protein_id = i['GBQualifier_value']
                                    if i['GBQualifier_name'] == 'translation': protein_sequence = i['GBQualifier_value']; proceed = 'yes'
                            if key == 'GBFeature_key':
                                if entry[key] == 'CDS':
                                    get_location = 'yes'
                                else: get_location = 'no' ### Get CDS location (occurs only following the GBFeature_key CDS, but not after)
                            if key == 'GBFeature_location' and get_location == 'yes':
                                cds_location = entry[key]
                except KeyError: null = []
                alt_seq='no'
                if proceed == 'yes':
                    if protein_sequence[0] != 'M': proceed = 'no'; alt_seq = 'yes'
                    else:
                        values = [accession,protein_id,protein_sequence]; values = string.join(values,'\t')+'\n'; datar.write(values)
                        datad.write(accession+'\t'+string.replace(cds_location,'..','\t')+'\n')
                else: ### Translate mRNA seq to protein
                    mRNA_seq = a['GBSeq_sequence']
                    mRNA_db={}; mRNA_db[accession] = '',mRNA_seq
                    translation_db = BuildInSilicoTranslations(mRNA_db)
                    for mRNA_AC in translation_db: ### Export in silico protein predictions
                        protein_id, protein_sequence,cds_location = translation_db[mRNA_AC]
                        values = [mRNA_AC,protein_id,protein_sequence]; values = string.join(values,'\t')+'\n'; datar.write(values)
                        datad.write(mRNA_AC+'\t'+string.replace(cds_location,'..','\t')+'\n')
                    if len(translation_db)==0 and alt_seq == 'yes': ### If no protein sequence starting with an "M" found, write the listed seq
                        values = [accession,protein_id,protein_sequence]; values = string.join(values,'\t')+'\n'; datar.write(values)
                        datad.write(accession+'\t'+string.replace(cds_location,'..','\t')+'\n')
            else: protein_sequence = a['GBSeq_sequence']
        index+=200
        """print 'accession',accession    
        print 'protein_id',protein_id
        print 'protein_sequence',protein_sequence
        print 'mRNA_seq',mRNA_seq,'\n' """

    end_time = time.time(); time_diff = int(end_time-start_time)
    print "finished in %s seconds" % time_diff
    datar.close()
    datad.close()

def importEnsemblTranscriptSequence(missing_protein_ACs):

    import_dir = '/AltDatabase/'+species+'/SequenceData' ### Multi-species fiel
    g = GrabFiles(); g.setdirectory(import_dir)
    seq_files = g.searchdirectory('.cdna.all.fa'); seq_files.sort(); filename = seq_files[-1]
    #filename = 'AltDatabase/'+species+'/SequenceData/'+species+'_ensembl_cDNA.fasta.txt'    
    start_time = time.time()
    output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-EnsInSilicoProt_sequences.txt'
    datar = export.ExportFile(output_file)

    output_file = 'AltDatabase/'+species+'/SequenceData/output/details/Transcript-EnsInSilicoProt_sequences.txt'
    datad = export.ExportFile(output_file)
    
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
                                protein_id, protein_seq, cds_location = translation_db[mRNA_AC]
                                values = [mRNA_AC,protein_id,protein_seq]; values = string.join(values,'\t')+'\n'; datar.write(values)
                                datad.write(mRNA_AC+'\t'+string.replace(cds_location,'..','\t')+'\n')
                                translated_mRNAs[mRNA_AC]=[]
                    ### Parse new line
                    #>ENST00000400685 cdna:known supercontig::NT_113903:9607:12778:1 gene:ENSG00000215618
                    t= string.split(data[1:],':'); sequence=''
                    transid_data = string.split(t[0],' '); transid = transid_data[0]
                    if '.' in transid:
                        transid = string.split(transid,'.')[0]
                    #try: ensembl_id,chr,strand,transid,prot_id = t
                    #except ValueError: ensembl_id,chr,strand,transid = t
        except IndexError: continue
        try:
            if data[0] != '>': sequence = sequence + data
        except IndexError:  continue

    #### Add the last entry
    if len(sequence) > 0:
        if transid in missing_protein_ACs:
            ### Perform in silico translation
            mRNA_db = {}; mRNA_db[transid] = '',sequence[1:]
            translation_db = BuildInSilicoTranslations(mRNA_db)
            for mRNA_AC in translation_db: ### Export in silico protein predictions
                protein_id, protein_seq, cds_location = translation_db[mRNA_AC]
                values = [mRNA_AC,protein_id,protein_seq]; values = string.join(values,'\t')+'\n'; datar.write(values)
                datad.write(mRNA_AC+'\t'+string.replace(cds_location,'..','\t')+'\n')
                translated_mRNAs[mRNA_AC]=[]
    datar.close()
    datad.close()
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "Ensembl transcript sequences analyzed in %d seconds" % time_diff

    missing_protein_ACs_inSilico=[]
    for mRNA_AC in missing_protein_ACs:
        if mRNA_AC not in translated_mRNAs:
            missing_protein_ACs_inSilico.append(mRNA_AC)
    print len(missing_protein_ACs_inSilico), 'Ensembl mRNAs without mRNA sequence NOT in silico translated (e.g., lncRNAs)', missing_protein_ACs_inSilico[:10]

def importEnsemblProteinSeqData(species,unique_ens_transcripts):
    from build_scripts import FeatureAlignment
    protein_relationship_file,protein_feature_file,protein_seq_fasta,null = FeatureAlignment.getEnsemblRelationshipDirs(species)
    ens_transcript_protein_db = FeatureAlignment.importEnsemblRelationships(protein_relationship_file,'transcript')

    unique_ens_proteins = {}
    for transcript in ens_transcript_protein_db:
        if transcript in unique_ens_transcripts:
            protein_id = ens_transcript_protein_db[transcript]
            if len(protein_id)>1:
                unique_ens_proteins[protein_id] = transcript
    ensembl_protein_seq_db = importEnsemblProtSeq(protein_seq_fasta,unique_ens_proteins)

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
    
    fn=filepath(filename); ensembl_protein_seq_db={}; sequence = ''; y=0
    for line in open(fn,'r').xreadlines():
        try: data, newline= string.split(line,'\n'); y+=1
        except ValueError: continue
        try:
            if data[0] == '>':
                if len(sequence) > 0:
                    try: ensembl_prot = ensembl_prot
                    except Exception: print data,y,t;kill
                    if ensembl_prot in unique_ens_proteins:
                        mRNA_AC = unique_ens_proteins[ensembl_prot]
                        values = string.join([mRNA_AC,ensembl_prot,sequence],'\t')+'\n'
                        export_data.write(values); ensembl_protein_seq_db[ensembl_prot] =  []
                ### Parse new line
                t= string.split(data[1:],' '); sequence=''
                ensembl_prot = t[0]
                if '.' in ensembl_prot:
                    ensembl_prot = string.split(ensembl_prot,'.')[0] ### Added to Ensembl after version 77
        except IndexError: continue
        try:
            if data[0] != '>': sequence = sequence + data
        except IndexError:  continue
    
    if ensembl_prot in unique_ens_proteins:
        mRNA_AC = unique_ens_proteins[ensembl_prot]
        values = string.join([mRNA_AC,ensembl_prot,sequence],'\t')+'\n'
        export_data.write(values); ensembl_protein_seq_db[ensembl_prot] =  []
    export_data.close()
    return ensembl_protein_seq_db

def importUniProtSeqeunces(species,transcripts_with_uniprots,transcripts_to_analyze):
    global n_terminal_seq; global c_terminal_seq
    n_terminal_seq={}; c_terminal_seq={}

    export_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-UniProt_sequences.txt'
    export_data = export.ExportFile(export_file)
    #filename = 'AltDatabase/'+species+'/SequenceData/'+'uniprot_trembl_sequence.txt'
    filename = 'AltDatabase/uniprot/'+species+'/uniprot_sequence.txt'    
    fn=filepath(filename); transcript_to_uniprot={}
    unigene_ensembl_up = {}
    for line in open(fn,'r').readlines():
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t')
        id=t[0];ac=t[1];ensembls=t[4];seq=t[2];type=t[6];unigenes=t[7];embls=t[9]
        ac=string.split(ac,','); embls=string.split(embls,',') #; ensembls=string.split(ensembls,','); unigenes=string.split(unigenes,',')
        if type != 'swissprot1': ### unclear why this condition was excluding swissprot so added 1 - version 2.1.1
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
            
    if len(transcripts_to_analyze)>0: ### Have to submitt ACs to report them
        print len(missing_protein_ACs), 'missing protein ACs for associated UniProt mRNAs and', len(transcript_to_uniprot), 'found.'
    print len(n_terminal_seq),len(c_terminal_seq),'N and C terminal, respectively...'
    return missing_protein_ACs
           
def BuildInSilicoTranslations(mRNA_db):
    """
    BI517798
    Seq('ATGTGGCCAGGAGACGCCACTGGAGAACATGCTGTTCGCCTCCTTCTACCTTCTGGATTT ...', IUPACUnambiguousDNA())
    213 MWPGDATGEHAVRLLLPSGFYPGFSWQYPGSVAFHPRPQVRDPGQRVPDASGRGRLVVRAGPAHPPGLPLLWEPLAIWGNRMPSHRLPLLPQHVRQHLLPHLHQRRPFPGHCAPGQVPQAPQAPLRTPGLCLPVGGGGCGHGPAAGEPTDRADKHTVGLPAAVPGEGSNMPGVPWQWPSLPVHHQVTCTVIIRSCGRPRVEKALRTRQGHESP

    211 MNGLEVAPPGLITNFSLATAEQCGQETPLENMLFASFYLLDFILALVGNTLALWLFIRDHKSGTPANVFLMHLAVADLSCVLVLPTRLVYHFSGNHWPFGEIACRLTGFLFYLNMYASIYFLTCISADRFLAIVHPVKSLKLRRPLYAHLACAFLWVVVAVAMAPLLVSPQTVQTNTRWVCLQLYREKAPTCLVSLGSGLHFPFITRSRVL
    """    
    translation_db={}
    
    from Bio.Seq import Seq

    ### New Biopython methods - http://biopython.org/wiki/Seq
    from Bio.Alphabet import generic_dna
    
    ### Deprecated
    #from Bio.Alphabet import IUPAC
    #from Bio import Translate ### deprecated

    #print 'Begining in silco translation for',len(mRNA_db),'sequences.'
    
    def cleanSeq(input_seq):
        """Wrapper for Biopython translate function.  Bio.Seq.translate will complain if input sequence is 
        not a mulitple of 3.  This wrapper function passes an acceptable input to Bio.Seq.translate in order to
        avoid this warning."""
        #https://github.com/broadinstitute/oncotator/pull/265/commits/94b20aabff48741a92b3f9e608e159957af6af30
        trailing_bases = len(input_seq) % 3
        if trailing_bases:
            input_seq = ''.join([input_seq, 'NN']) if trailing_bases == 1 else ''.join([input_seq, 'N'])
        return input_seq
 
    first_time = 1
    for mRNA_AC in mRNA_db:
        if mRNA_AC == 'AK025306': print '@@@@@@@@@@^^^^AK025306...attempting in silico translation'
        temp_protein_list=[]; y=0
        protein_id,sequence = mRNA_db[mRNA_AC]
        if protein_id == '': protein_id = mRNA_AC+'-PEP'
        original_seqeunce = sequence
        sequence = string.upper(sequence)
        loop=0
        while (string.find(sequence,'ATG')) != -1:  #while there is still a methionine in the DNA sequence, reload this DNA sequence for translation: find the longest ORF
            x = string.find(sequence,'ATG') #before doing this, need to find the start codon ourselves
            y += x  #maintain a running count of the sequence position
            if loop!=0: y+=3 ### This accounts for the loss in sequence_met
            #if y<300: print original_seqeunce[:y+2], x
            sequence_met = sequence[x:]  #x gives us the position where the first Met* is.
            
            ### New Biopython methods - http://biopython.org/wiki/Seq
            dna_clean = cleanSeq(sequence_met)
            dna_seq = Seq(dna_clean, generic_dna)
            prot_seq = dna_seq.translate(to_stop=True)
            
            ### Deprecated code
            #seq_type = IUPAC.unambiguous_dna
            #dna_seq = Seq(sequence_met,seq_type)
            #standard_translator = Translate.unambiguous_dna_by_id[1]
            #prot_seq = standard_translator.translate_to_stop(dna_seq) #convert the dna to protein sequence
            
            #prot_seq_string = prot_seq.tostring()
            prot_seq_string = str(prot_seq)
            prot_seq_tuple = len(prot_seq_string),y,prot_seq_string,dna_seq  #added DNA sequence to determine which exon we're in later  
            temp_protein_list.append(prot_seq_tuple) #create a list of protein sequences to select the longest one
            sequence = sequence_met[3:]  # sequence_met is the sequence after the first or proceeduring methionine, reset the sequence for the next loop
            loop+=1
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
                            pos1 = protein_data[1]
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
            #dna_seq_string_coding_to_end = coding_dna_seq_string.tostring()
            dna_seq_string_coding_to_end = str(coding_dna_seq_string)
            coding_dna_seq_string = dna_seq_string_coding_to_end[0:dl]
            cds_location = str(pos1+1)+'..'+str(pos1+len(prot_seq_string)*3+3)
            ### Determine if a stop codon is in the sequence or there's a premature end
            coding_diff = len(dna_seq_string_coding_to_end) - len(coding_dna_seq_string)
            if coding_diff > 4: stop_found = 'stop-codon-present'
            else: stop_found = 'stop-codon-absent'
            #print [mRNA_AC],[protein_id],prot_seq_string[0:10]
            if mRNA_AC == 'AK025306': print '*********AK025306',[protein_id],prot_seq_string[0:10]
            translation_db[mRNA_AC] = protein_id,prot_seq_string,cds_location
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

def compareProteinFeaturesForPairwiseComps(probeset_protein_db,protein_seq_db,probeset_gene_db,species,array_type):

    if (array_type == 'junction' or array_type == 'RNASeq') and data_type != 'null': 
        export_file = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/probeset-protein-dbase_seqcomp.txt'
    else: export_file = 'AltDatabase/'+species+'/'+array_type+'/probeset-protein-dbase_seqcomp.txt'  
    fn=filepath(export_file); export_data1 = open(fn,'w')
    title_row = 'Probeset\tAligned protein_id\tNon-aligned protein_id\n'; export_data1.write(title_row)
    
    minimal_effect_db={}; accession_seq_db={};start_time = time.time()
    print "Comparing protein features for all pair-wise comparisons"
    for probeset in probeset_protein_db:
        geneid = probeset_gene_db[probeset][0] ###If one probeset
        match_list,null_list = probeset_protein_db[probeset]
        prot_match_list=[]; prot_null_list=[]
        for protein_ac in match_list:
            protein_seq = protein_seq_db[protein_ac][0]
            prot_match_list.append([protein_ac,protein_seq])
        for protein_ac in null_list:
            protein_seq = protein_seq_db[protein_ac][0]
            prot_null_list.append([protein_ac,protein_seq])
            
        ### Compare all possible protein combinations to find those with the minimal change in (1) sequence and (2) domain compositions
        if len(prot_match_list)>0 and len(prot_null_list)>0:
            results_list=[]
            for mi in prot_match_list:
                for ni in prot_null_list:
                    results = characterizeProteinLevelExonChanges(probeset,geneid,mi,ni,array_type)
                    results_list.append(results)
            results_list.sort()
            hit_ac = results_list[0][-2]; null_ac = results_list[0][-1]
            values = string.join([probeset,hit_ac,null_ac],'\t')+'\n'; export_data1.write(values)
            
            #minimal_effect_db[probeset] = [hit_ac,null_ac]
            accession_seq_db[hit_ac] = [protein_seq_db[hit_ac][0]]
            accession_seq_db[null_ac] = [protein_seq_db[null_ac][0]]
            #print results_list[0][-2],results_list[0][-1]
            #print probeset,len(prot_match_list),len(prot_null_list),results_list;kill

    print len(minimal_effect_db),"optimal pairs of probesets linked to Ensembl found."
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "databases built in %d seconds" % time_diff
    export_data1.close()
    
    """
    title_row = 'Probeset\tAligned protein_id\tNon-aligned protein_id'
    export_file = 'AltDatabase/'+species+'/'+array_type+'/probeset-protein-dbase_seqcomp.txt'     
    exportSimple(minimal_effect_db,export_file,title_row)"""

    if (array_type == 'junction' or array_type == 'RNASeq') and data_type != 'null':
        export_file = 'AltDatabase/'+species+'/'+array_type+'/'+data_type+'/SEQUENCE-protein-dbase_seqcomp.txt' 
    else: export_file = 'AltDatabase/'+species+'/'+array_type+'/SEQUENCE-protein-dbase_seqcomp.txt'     
    exportSimple(accession_seq_db,export_file,'')

def exportSimple(db,output_file,title_row):
    data = export.ExportFile(output_file)
    if len(title_row)>0: data.write(title_row+'\n')
    for key in db:
        try: values = string.join([key] + db[key],'\t')+'\n'; data.write(values)
        except TypeError:
            list2=[]
            for list in db[key]: list2+=list
            values = string.join([key] + list2,'\t')+'\n'; data.write(values)
    data.close()
    print output_file, 'exported...'
    
def characterizeProteinLevelExonChanges(uid,array_geneid,hit_data,null_data,array_type):    
            domains=[]; functional_attribute=[]; probeset_num = 1
            if '|' in uid:
                probeset,probeset2 = string.split(uid,'|')
                probeset_num = 2
            else: probeset = uid
            hv=1; pos_ref_AC = hit_data[0]; neg_ref_AC = null_data[0]             
            if hv!=0:
                neg_coding_seq = null_data[1]; pos_coding_seq = hit_data[1]
                neg_length = len(neg_coding_seq); pos_length = len(pos_coding_seq)
                pos_length = float(pos_length); neg_length = float(neg_length)
                if array_geneid in protein_ft_db:
                        protein_ft = protein_ft_db[array_geneid]
                        neg_ft_missing,pos_ft_missing = ExonAnalyze_module.compareProteinFeatures(protein_ft,neg_coding_seq,pos_coding_seq)
                        for (pos,blank,ft_name,annotation) in pos_ft_missing: domains.append(ft_name)
                        for (pos,blank,ft_name,annotation) in neg_ft_missing:
                            ###If missing from the negative list, it is present in the positive state
                            domains.append(ft_name)
                if pos_coding_seq[:5] != neg_coding_seq[:5]: function_var = 'alt-N-terminus';functional_attribute.append(function_var)
                if pos_coding_seq[-5:] != neg_coding_seq[-5:]: function_var = 'alt-C-terminus';functional_attribute.append(function_var)
                ### Record change in peptide size
                protein_length_diff = abs(pos_length-neg_length)
            results = [len(functional_attribute),len(domains),protein_length_diff,pos_ref_AC,neg_ref_AC] ###number of domains, number N or C-terminal differences, protein length difference
            return results

############# Code currently not used (LOCAL PROTEIN SEQUENCE ANALYSIS) ##############
def importUCSCSequenceAssociations(species,transcripts_to_analyze):
    """ NOTE: This method is currently not used, by the fact that the kgXref is not downloaded from the
    UCSC ftp database. This file relates mRNAs primarily to UniProt rather than GenBank and thus is less
    ideal than the direct NCBI API based method used downstream """
    
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
    
def importUCSCSequences(missing_protein_ACs):
    start_time = time.time()
    filename = 'AltDatabase/'+species+'/SequenceData/mrna.fa'
    output_file = 'AltDatabase/'+species+'/SequenceData/output/sequences/Transcript-InSilicoProt_sequences.txt'
    datar = export.ExportFile(output_file)  

    output_file = 'AltDatabase/'+species+'/SequenceData/output/details/Transcript-InSilicoProt_sequences.txt'
    datad = export.ExportFile(output_file)
    
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
                                    protein_id, protein_seq, cds_location = translation_db[mRNA_AC]
                                    values = [mRNA_AC,protein_id,protein_seq]; values = string.join(values,'\t')+'\n'; datar.write(values)
                                    datad.write(mRNA_AC+'\t'+string.replace(cds_location,'..','\t')+'\n')
                                    translated_mRNAs[mRNA_AC]=[]
                        ### Parse new line
                        values = string.split(data,' '); accession = values[0][1:]; sequence = '|'; continue
                except IndexError: null = []       
                try:
                    if data[0] != '>': sequence = sequence + data
                except IndexError: print kill; continue
    datar.close()
    datad.close()
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "UCSC mRNA sequences analyzed in %d seconds" % time_diff

    missing_protein_ACs_inSilico=[]
    for mRNA_AC in missing_protein_ACs:
        if mRNA_AC not in translated_mRNAs:
            missing_protein_ACs_inSilico.append(mRNA_AC)
    print len(missing_protein_ACs_inSilico), 'mRNAs without mRNA sequence NOT in silico translated (e.g., lncRNAs)',missing_protein_ACs_inSilico[:10]
    return missing_protein_ACs_inSilico

############# END Code currently not used (LOCAL PROTEIN SEQUENCE ANALYSIS) ##############

def runProgram(Species,Array_type,Data_type,translate_seq,run_seqcomp):
    global species; global array_type; global translate; global data_type; global test; global test_genes
    species = Species; array_type = Array_type; translate = translate_seq; data_type = Data_type
    test = 'no'; test_genes = ['ENSMUSG00000029467']
    if array_type == 'gene' or array_type == 'exon' or data_type == 'exon':
        compare_all_features = 'no' 
        print 'Begin Exon-based compareProteinComposition'
        try: compareProteinComposition(species,array_type,translate,compare_all_features)
        except Exception:
            compareExonComposition(species,array_type)
            compareProteinComposition(species,array_type,translate,compare_all_features)
        if run_seqcomp == 'yes':
            compare_all_features = 'yes'; translate = 'no'
            compareProteinComposition(species,array_type,translate,compare_all_features)
    elif array_type == 'AltMouse' or array_type == 'junction' or array_type == 'RNASeq':
        """ Modified on 11/26/2017(2.1.1) - integrated compareExonComposition and compareProteinComposition
        with the goal of getting protein sequences for ALL UCSC and Ensembl transcripts (known + in silico)
        to look at protein length differences before picking the top two compared isoforms. """
        compare_all_features = 'no' 
        print 'Begin Junction-based compareProteinComposition'
        compareProteinComposition(species,array_type,translate,compare_all_features) 
        if run_seqcomp == 'yes':
            compare_all_features = 'yes'; translate = 'no'
            compareProteinComposition(species,array_type,translate,compare_all_features)
    
def runProgramTest(Species,Array_type,Data_type,translate_seq,run_seqcomp):
    global species; global array_type; global translate; global data_type 
    species = Species; array_type = Array_type; translate = translate_seq; data_type = Data_type
    if array_type == 'gene' or array_type == 'exon' or data_type == 'exon':
        compare_all_features = 'no'
        compareProteinComposition(species,array_type,translate,compare_all_features)
        if run_seqcomp == 'yes':
            compare_all_features = 'yes'; translate = 'no'
            compareProteinComposition(species,array_type,translate,compare_all_features)
    elif array_type == 'AltMouse' or array_type == 'junction' or array_type == 'RNASeq':
        compare_all_features = 'no'
        compareProteinComposition(species,array_type,translate,compare_all_features)
        if run_seqcomp == 'yes':
            compare_all_features = 'yes'; translate = 'no'
            compareProteinComposition(species,array_type,translate,compare_all_features)
            
if __name__ == '__main__':
    species = 'Mm'; array_type = 'AltMouse'; translate='no'; run_seqcomp = 'no'; data_type = 'exon'
    species = 'Mm'; array_type = 'RNASeq'; translate='yes'
    #species = 'Dr'; array_type = 'RNASeq'; translate='yes'; data_type = 'junction'
    test='no'
    a = ['ENSMUST00000138102', 'ENSMUST00000193415', 'ENSMUST00000124412', 'ENSMUST00000200569']
    importEnsemblTranscriptSequence(a);sys.exit()
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'    
    ens_transcript_exon_db,ens_gene_transcript_db,ens_gene_exon_db = importEnsExonStructureDataSimple(filename,species,{},{},{},{})
    #print len(ens_transcript_exon_db), len(ens_gene_transcript_db), len(ens_gene_exon_db);sys.exit()
    #runProgramTest(species,array_type,data_type,translate,run_seqcomp)
    runProgram(species,array_type,data_type,translate,run_seqcomp)
