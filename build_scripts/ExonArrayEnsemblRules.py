###ExonArrayEnsemblRules
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
import update
import export
import math
from build_scripts import EnsemblImport; reload(EnsemblImport)
from build_scripts import ExonArrayAffyRules
from build_scripts import SubGeneViewerExport

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

################# Parse and Analyze Files
def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def grabRNAIdentifiers(mrna_assignment):
    ensembl_ids=[]; mRNA_ids=[]; mRNA_entries = string.split(mrna_assignment,' /// ')
    for entry in mRNA_entries:
        mRNA_info = string.split(entry,' // '); mrna_ac = mRNA_info[0]
        if 'GENSCAN' not in mrna_ac:
            if 'ENS' in mrna_ac: ensembl_ids.append(mrna_ac)
            else:
                try: int(mrna_ac[-3:]); mRNA_ids.append(mrna_ac)
                except ValueError: continue
    ensembl_ids = unique.unique(ensembl_ids); mRNA_ids = unique.unique(mRNA_ids)
    return ensembl_ids, mRNA_ids

def getProbesetAssociations(filename,ensembl_exon_db,ens_transcript_db,source_biotype):
    #probeset_db[probeset] = affygene,exons,probe_type_call,ensembl
    fn=filepath(filename); global transcript_cluster_data; global exon_location; global trans_annotation_db
    probe_association_db={}; transcript_cluster_data = {}; trans_annotation_db = {}
    print "Begin reading",filename; entries=0
    exon_location={}
    for line in open(fn,'rU').xreadlines():             
        data,null = string.split(line,'\n')
        data = string.replace(data,'"',''); y = 0
        t = string.split(data,',')
        if data[0] != '#' and 'probeset_id' in data:
            affy_headers = t
            for header in affy_headers:
                index = 0
                while index < len(affy_headers):
                    if 'probeset_id' == affy_headers[index]: pi = index
                    if 'start' == affy_headers[index]: st = index
                    if 'stop' == affy_headers[index]: sp = index
                    if 'level' == affy_headers[index]: lv = index
                    if 'exon_id' == affy_headers[index]: ei = index
                    if 'transcript_cluster_id' == affy_headers[index]: tc = index
                    if 'seqname' == affy_headers[index]: sn = index
                    if 'strand' == affy_headers[index]: sd = index
                    if 'mrna_assignment' == affy_headers[index]: ma = index
                    #if 'fl' == affy_headers[index]: fn = index
                    #if 'mrna' == affy_headers[index]: mr = index
                    #if 'est' == affy_headers[index]: es = index
                    #if 'ensGene' == affy_headers[index]: eg = index
                    index += 1      
        elif data[0] != '#' and 'probeset_id' not in data:
            try:
              entries+=1
              try: probeset_id=int(t[pi]);transcript_cluster_id=int(t[tc]); chr=t[sn];strand=t[sd]
              except Exception: print affy_headers; print index; sys.exit()
              if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
              if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
              start=int(t[st]);stop=int(t[sp]); exon_type=t[lv]; #fl=int(t[fn]); mRNA=int(t[mr]); est=int(t[es]); ensembl=int(t[eg])
              continue_analysis = 'no'

              if transcript_cluster_id not in trans_annotation_db:
                  mrna_assignment = t[ma]; ens_list=[]
                  if len(mrna_assignment)>4:
                      ensembl_data = string.split(mrna_assignment,' /// ')
                      for entry in ensembl_data:
                          if 'ENS' == entry[:3]:
                              ens_entry = string.split(entry,' // ')
                              if ens_entry[0] in ens_transcript_db: ens_list.append(ens_transcript_db[ens_entry[0]][0])
                      if len(ens_list)>0:
                          ens_list = unique.unique(ens_list)
                          trans_annotation_db[transcript_cluster_id] = ens_list

              if source_biotype == 'ncRNA':
                  ### transcript_cluster_ids are only informative for looking at mRNA encoding genes (can combine diff. ncRNAs in diff. introns of the same gene)
                  transcript_cluster_id = probeset_id
                    
              if test=='yes': ###used to test the program for a single gene
                   if transcript_cluster_id in test_cluster: continue_analysis='yes'
              else: continue_analysis='yes'
              if continue_analysis=='yes':
                  try: exon_location[transcript_cluster_id,chr,strand].append((start,stop,exon_type,probeset_id))
                  except KeyError: exon_location[transcript_cluster_id,chr,strand] = [(start,stop,exon_type,probeset_id)]
                  ###Assign constitutive information on a per probeset basis since per gene (unlike transcript_cluster centric analyses)
                  const_info = [transcript_cluster_id] #[ensembl,fl,est,probeset_id,transcript_cluster_id]
                  transcript_cluster_data[probeset_id] = const_info
            except ValueError:
                continue ###Control probeset with no exon information
    
    for key in exon_location: exon_location[key].sort()
    if test=='yes':
        for i in exon_location:
            print 'exon_location    ',i
            for e in exon_location[i]: print e
            
    print entries,"entries in input file"
    print len(transcript_cluster_data),"probesets,", len(exon_location),"transcript clusters"

    print "Matching up ensembl and ExonArray probesets based on chromosomal location..."
    ###Re-organize ensembl and probeset databases based on gene position (note the ensembl db is strucutred similiarly to the transcript_cluster db)
    ###however, an additional call is required to match these up.
    ###ensembl_exon_db[(geneid,chr,strand)] = [[E5,exon_info]] #exon_info = (exon_start,exon_stop,exon_id,exon_annot)
    ensembl_gene_position_pos_db, ensembl_gene_position_neg_db = getGenePositions(ensembl_exon_db,'yes')
    affy_gene_position_pos_db, affy_gene_position_neg_db = getGenePositions(exon_location,'no')

    if test=='yes':
        for chr in affy_gene_position_pos_db:
            b = affy_gene_position_pos_db[chr]
            for i in b:
                for e in b[i]:
                    for pos in e: print 'pos',chr,i,pos
        for chr in affy_gene_position_neg_db:
            b = affy_gene_position_neg_db[chr]
            for i in b:
                for e in b[i]:
                    for pos in e: print 'neg',chr,i,pos
    ###Dump memory

    """print ensembl_exon_db
    print exon_location
    print ensembl_gene_position_pos_db
    print affy_gene_position_pos_db
    killer"""
    
    exon_location={}; global ensembl_probeset_db
    print "First round (positive strand)..."
    #global merged_gene_loc
    merged_gene_loc={}; no_match_list=[]
    merged_gene_loc,no_match_list = getChromosomalOveralap(affy_gene_position_pos_db,ensembl_gene_position_pos_db,merged_gene_loc,no_match_list)
    print "Second round (negative strand)..."
    merged_gene_loc,no_match_list = getChromosomalOveralap(affy_gene_position_neg_db,ensembl_gene_position_neg_db,merged_gene_loc,no_match_list)
    merged_gene_loc = resolveProbesetAssociations(merged_gene_loc,no_match_list)
    ensembl_probeset_db = reorderEnsemblLinkedProbesets(merged_gene_loc,transcript_cluster_data,trans_annotation_db)

    if test=='yes':
        for i in merged_gene_loc:
            print 'merged_gene_loc    ',i
            for e in merged_gene_loc[i]: print e
            
    if test=='yes':
        for i in ensembl_probeset_db:
            print 'ensembl_probeset_db    ',i
            for e in ensembl_probeset_db[i]: print e

    ###Dump memory
    affy_gene_position_neg_db={};ensembl_gene_position_pos_db={}; no_match_list={}

    print "Begining to assembl constitutive annotations (Based on Ensembl/FL/EST evidence)..."
    transcript_cluster_data={}
    print "Assigning exon-level Ensembl annotations to probesets (e.g. exon number) for", len(ensembl_probeset_db),'genes.'
    ensembl_probeset_db = annotateExons(ensembl_probeset_db,exon_clusters,ensembl_exon_db,exon_region_db,intron_retention_db,intron_clusters,ucsc_splicing_annot_db)
    print "Exporting Ensembl-merged probeset database to text file for", len(ensembl_probeset_db),'genes.'
    exportEnsemblLinkedProbesets(arraytype, ensembl_probeset_db,species)
    return ensembl_probeset_db

################# Identify overlap between Ensembl and transcript_clusters
def getGenePositions(exon_db,call):
    chromosome_pos_db={}; chromosome_neg_db={}
    for key in exon_db:
        t=[]
        strand = key[-1]; chr = key[1]; index1= 0; index2 = -1; indexa = 0; indexb = 1
        if strand == '-': index1= -1; index2 = 0; indexa = 1; indexb = 0
        if call == 'yes':
            ### For Ensembl exon coordinates - get Gene Coordinates
            
            ### Doesn't work for NKX2-5 in human for exon array - first exon spans the last exon (by AltAnalyze's definitions)
            #geneStart = exon_db[key][index1][1][indexa] #first sorted exon
            #geneStop = exon_db[key][index2][1][indexb]
            for exon_data in exon_db[key]:
                t.append(exon_data[1][0])
                t.append(exon_data[1][1])
        else:
            ### For transcript cluster data (slightly different format than above) - get Gene Coordinates
            chr = chr[3:]
            if '_' in chr: c = string.split(chr,'_'); chr=c[0] ###For _unknown chr from Affy's annotation file
            #geneStart = exon_db[key][index1][indexa] #first sorted exon
            #geneStop = exon_db[key][index2][indexb]
            for exon_data in exon_db[key]:
                t.append(exon_data[0])
                t.append(exon_data[1])
        #t=[]; t.append(geneStart); t.append(geneStop); t.sort()
        t.sort(); t = [t[0],t[-1]]
        if strand == '-':
            if chr in chromosome_neg_db:
                gene_position_db = chromosome_neg_db[chr]
                gene_position_db[tuple(t)] = key,exon_db[key]
            else:
                gene_position_db={}
                gene_position_db[tuple(t)] = key,exon_db[key]
                chromosome_neg_db[chr] = gene_position_db
        if strand == '+':
            if chr in chromosome_pos_db:
                gene_position_db = chromosome_pos_db[chr]
                gene_position_db[tuple(t)] = key,exon_db[key]
            else:
                gene_position_db={}
                gene_position_db[tuple(t)] = key,exon_db[key]
                chromosome_pos_db[chr] = gene_position_db
    return chromosome_pos_db,chromosome_neg_db

def getChromosomalOveralap(affy_chr_db,ensembl_chr_db,ensembl_transcript_clusters,no_match_list):
    """Find transcript_clusters that have overlapping start positions with Ensembl gene start and end (based on first and last exons)"""
    ###exon_location[transcript_cluster_id,chr,strand] = [(start,stop,exon_type,probeset_id)]
    y = 0; l =0; multiple_ensembl_associations=[]
    ###(bp1,ep1) = (47211632,47869699); (bp2,ep2)  =  (47216942, 47240877)
    for chr in affy_chr_db:
        try:
            ensembl_db = ensembl_chr_db[chr]
            affy_db = affy_chr_db[chr]
            for (bp1,ep1) in affy_db:
                x = 0
                transcript_cluster_key = affy_db[(bp1,ep1)][0]
                for (bp2,ep2) in ensembl_db:
                    y += 1; ensembl = ensembl_db[(bp2,ep2)][0][0]
                    #print ensembl, transcript_cluster_key, (bp2,ep2),(bp1,ep1);kill
                    ###if the two gene location ranges overlapping
                    ###affy_probeset_info = (start,stop,exon_type,probeset_id)
                    if ((bp1 >= bp2)  and (ep2 >= bp1)) or ((bp2 >= bp1)  and (ep1 >= bp2)):
                        x = 1; affy_probeset_info = affy_db[(bp1,ep1)][1]; ensembl_key = ensembl_db[(bp2,ep2)][0],(bp2,ep2),affy_probeset_info
                        try:ensembl_transcript_clusters[transcript_cluster_key].append(ensembl_key)
                        except KeyError: ensembl_transcript_clusters[transcript_cluster_key] = [ensembl_key]
                        l += 1
                if x == 0: no_match_list.append(transcript_cluster_key)
        except KeyError: print chr, 'data not found'
            
    print "Transcript Clusters overlapping with Ensembl:",len(ensembl_transcript_clusters)
    print "With NO overlapp",len(no_match_list)
    return ensembl_transcript_clusters,no_match_list

def resolveProbesetAssociations(ensembl_transcript_clusters,no_match_list):
    ensembl_transcript_clusters2={}; x = 0
    for transcript_cluster_key in ensembl_transcript_clusters:
        ensembl_groups = ensembl_transcript_clusters[transcript_cluster_key]
        probeset_data_list = ensembl_groups[0][2]
        if len(ensembl_groups)>1:
            ###although the affy_transcript_cluster may overlapp, each probe does not
            ###Therefore associate the appropriate probesets with the appropriate ensembl ids
            x += 1
            for probeset_data in probeset_data_list:
                (bp1,ep1,exon_type,probeset_id) = probeset_data
                y = 0
                for data in ensembl_groups:
                    ensembl_id = data[0]
                    (bp2,ep2) = data[1]
                    if ((bp1 >= bp2)  and (ep2 >= bp1)) or ((bp2 >= bp1)  and (ep1 >= bp2)):
                        y = 1
                        ###Reduce data to Ensembl gene level (forget transcript_cluster_ids)
                        try:
                            ensembl_transcript_clusters2[ensembl_id].append(probeset_data)
                        except KeyError:
                            ensembl_transcript_clusters2[ensembl_id] = [probeset_data]
        else:
            for data in ensembl_groups:
                ensembl_id = data[0]
                probeset_data = data[2]
                try:
                    ensembl_transcript_clusters2[ensembl_id].append(probeset_data)
                except KeyError:
                    ensembl_transcript_clusters2[ensembl_id] = [probeset_data]

    print "Unique Ensembl genes linked to one or more transcript clusters:",len(ensembl_transcript_clusters2)
    print "Ensembl genes linked to more than one transcript cluster:",x
    print "Transcript clusters with NO links to Ensembl", len(no_match_list)
    print "\nNOTE: if multiple clusters linked to an Ensembl, exons outside the Ensembl overlapp with be discarded\n"
    return ensembl_transcript_clusters2

################# Annotate Exons based on location
def convertListString(list,delimiter):
    list2 = []; list = unique.unique(list); list.sort()
    for i in list:
        if len(str(i))>0: list2.append(str(i))
    str_values = string.join(list2,delimiter)
    return str_values
    
def annotateExons(exon_location,exon_clusters,ensembl_exon_db,exon_region_db,intron_retention_db,intron_clusters,ucsc_splicing_annot_db):
    ###Annotate Affymetrix probesets independent of other annotations.  A problem with this method
    ###Is that it fails to recognize distance probesets that probe different regions of the same exon.
    ###(start,stop,probeset_id,exon_class) = probeset_data
    ###exon_clusters is from EnsemblImport: exon_data = exon,(start,stop),[accessions],[types]

    ###Reverse exon entries based on strand orientation
    for key in exon_location:
        exon_location[key].sort(); strand = key[-1]
        if strand == "-": exon_location[key].reverse()
                    
    #alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','x','y','z']
    probeset_aligments={}; exon_location2={}
    x = 0; p = 0; count = 0
    for key in exon_location:
        old_index = 0; index2 = 1; index3 = 0; exon_list=[]; strand = key[-1]; last_added = 'null'
        #if key[-1] == '-':
        gene = key[0]
        new_key = (key[0],key[1][3:],key[2])
        for exon in exon_location[key]:  ###probeset location database
            start = int(exon[0]);  stop = int(exon[1]); probeset = exon[2]; probeset_region_exon_db={}
            #if count < 20: print probeset
            #else: die
            ens_exonid_values=''; ens_constitutive=''; splice_event_value=''; splice_junction_values=''; region_value=''
            y = 0; u = 0; junction = 'no'; new_exon_data = ''; exon_start='';exon_stop=''
            count+=1 ###Create the lists here so we can store data for multiple exon hits, if they occur (not good if they do)
            temp_probeset_exon_id=[]; temp_probeset_annot=[]; region_ids=[]; splice_events=[]; splice_junctions=[]; reg_start=[];reg_stop=[]
            #print probeset,start,stop;kill

            if new_key in intron_retention_db:
                for exon_data in intron_retention_db[new_key]:
                    t_start = exon_data[0]; t_stop = exon_data[1]; ed = exon_data[2]
                    if ((start >= t_start) and (stop <= t_stop)): splice_events.append('intron-retention')
                    
            ###Match up probesets specifically with UCSC exon annotations (need to do it here to precisely align probesets)
            if gene in ucsc_splicing_annot_db:
                ucsc_events = ucsc_splicing_annot_db[gene]
                for (r_start,r_stop,splice_event) in ucsc_events:
                        if ((start >= r_start) and (start < r_stop)) or ((stop > r_start) and (stop <= r_stop)): splice_events.append(splice_event)
                        elif ((r_start >= start) and (r_start <= stop)) or ((r_stop >= start) and (r_stop <= stop)): splice_events.append(splice_event)
            for exon_values in exon_clusters[key]: ###Ensembl location database
                ref_start = exon_values[1][0]; ref_stop = exon_values[1][1]; ens_exon_ida = exon_values[2][0]
                """if probeset == 'G6857086:E15' and ens_exon_ida == 'ENSMUSE00000101953':
                    print start,ref_start,stop,ref_stop,exon_values[2]
                    if ((start >= ref_start) and (stop <= ref_stop)) or ((start >= ref_start) and (start <= ref_stop)) or ((stop >= ref_start) and (stop <= ref_stop)): print 'good'
                    else: print 'bad'
                    kill"""
                if ((start >= ref_start) and (stop <= ref_stop)) or ((start >= ref_start) and (start < ref_stop)) or ((stop > ref_start) and (stop <= ref_stop)):
                    new_index = exon_values[0]
                    #probeset_exon_list.append((new_index,key,start,stop,ref_start,ref_stop))
                    ###Grab individual exon annotations
                    #temp_probeset_exon_id=[]; temp_probeset_annot=[]; region_ids=[]; splice_events=[]; splice_junctions=[]
                    for exon_data in ensembl_exon_db[new_key]:
                        t_start = exon_data[1][0]; t_stop = exon_data[1][1]; ed = exon_data[1][2]
                        if ((start >= t_start) and (start < t_stop)) or ((stop > t_start) and (stop <= t_stop)):
                            if ((start >= t_start) and (start < t_stop)) and ((stop > t_start) and (stop <= t_stop)):
                                ### Thus the probeset is completely in this exon
                                try: probeset_aligments[probeset].append([t_start,t_stop])
                                except KeyError: probeset_aligments[probeset]=[[t_start,t_stop]]
                            """if probeset == 'G6857086:E15' and ens_exon_ida == 'ENSMUSE00000101953':
                                print probeset, start,t_start,stop,t_stop,exon_values[2], 'goood', ed.ExonID()"""          
                            temp_probeset_exon_id.append(ed.ExonID())
                            block_db = exon_region_db[gene]
                            ###Annotate exon regions
                            for rd in block_db[new_index]:
                                if strand == '-': r_stop = rd.ExonStart(); r_start = rd.ExonStop()
                                else: r_start = rd.ExonStart(); r_stop = rd.ExonStop()
                                if ((start >= r_start) and (start < r_stop)) or ((stop > r_start) and (stop <= r_stop)):
                                    reg_start.append(r_start); reg_stop.append(r_stop)
                                    region_ids.append(rd.ExonRegionID()) ###only one region per probeset unless overlapping with two
                                    try: splice_events.append(rd.AssociatedSplicingEvent());splice_junctions.append(rd.AssociatedSplicingJunctions())
                                    except AttributeError: splice_events = splice_events ###occurs when the associated exon is not a critical exon
                                    temp_probeset_annot.append(rd.Constitutive())
                            """region_ids = unique.unique(region_ids)
                            if len(region_ids)>1:
                                print key, region_ids
                                for rd in block_db[new_index]: print start, rd.ExonStart(), stop, rd.ExonStop(), probeset, new_index,rd.RegionNumber();die"""
                        #elif ((start >= t_start) and (start <= t_stop)): ###Thus only the start is inside an exon
                            #if strand == '+': print probeset,key, ed.ExonID(), start, stop, t_start,t_stop;kill
                    if len(temp_probeset_exon_id)>0:
                        ens_exonid_values = convertListString(temp_probeset_exon_id,'|')
                        if 'yes' in temp_probeset_annot: ens_constitutive = 'yes'
                        elif '1' in temp_probeset_annot: ens_constitutive = '1'
                        else: ens_constitutive = convertListString(temp_probeset_annot,'|')
                        region_value = convertListString(region_ids,'|'); splice_event_value = convertListString(splice_events,'|')
                        exon_start = convertListString(reg_start,'|'); exon_stop = convertListString(reg_stop,'|')
                        splice_junction_values = convertListString(splice_junctions,'|');u = 1
                        splice_event_value = string.replace(splice_event_value,'intron-retention','exon-region-exclusion')
                        ned = EnsemblImport.ProbesetAnnotation(ens_exonid_values, ens_constitutive, region_value, splice_event_value, splice_junction_values,exon_start,exon_stop)
                        new_exon_data = exon,ned
                        #print key; exon; exon_values[2],exon_values[3],dog
                    if new_index == old_index:           
                        exon_info = 'E'+str(old_index)+'-'+str(index2); index_val = 'old'
                        y = 1; index2 += 1;  last_added = 'exon' #;last_start = start; last_stop = stop
                        if exon_values == exon_clusters[key][-1]: final_exon_included = 'yes'
                        else: final_exon_included = 'no'
                        #print exon_clusters[key][-1], exon_values, old_index, exon_info, probeset,final_exon_included;kill
                    else:
                        index2 = 1
                        exon_info = 'E'+str(new_index)+'-'+str(index2); index_val = old_index
                        old_index = new_index
                        y = 1; index2 += 1; last_added = 'exon' #;last_start = start; last_stop = stop
                        if exon_values == exon_clusters[key][-1]: final_exon_included = 'yes'
                        else: final_exon_included = 'no'
                        #print 'blah',new_index,old_index
                    """if len(exon)>7: ###Therefore, this probeset associated with distinct exon clusters (bad design)
                        index2 = (index2 - 1); x += 1"""
                    if u != 1: ###Therefore, there were specific exon associations available
                        ens_exonid_values = convertListString(exon_values[2],'|');

                        if 'yes' in exon_values[3]: ens_constitutive = 'yes'
                        elif '1' in exon_values[3]: ens_constitutive = '1'
                        else: ens_constitutive = convertListString(exon_values[3],'|')
                        ned = EnsemblImport.ProbesetAnnotation(ens_exonid_values, ens_constitutive, region_value, splice_event_value, splice_junction_values,exon_start,exon_stop)
                        new_exon_data = exon,ned
                        """if len(exon)>7: ###Therefore, this probeset associated with distinct exon clusters (bad design)
                            index2 = (index2 - 1); x += 1"""
                    continue
                #if len(probeset_exon_list)>1: print probeset,probeset_exon_list;kill
            #if probeset == '2948539': print y,last_added,final_exon_included,old_index,new_index,index2;kill
            if y == 1:
                """if len(exon)>7:
                    exon = exon[0:5]
                    if index_val != 'old':
                        old_index = index_val; y=0; junction = 'yes'
                else:"""
                exon_info = exon_info,new_exon_data; exon_list.append(exon_info)
            if y == 0: ###Thus this probeset is not in any ensembl exons
                ###Check if the first probesets are actually in an intron
                #if probeset == '2948539': print  probeset,'a',last_added,old_index,new_index,index2
                if key in intron_clusters:
                    for intron_values in intron_clusters[key]: ###Ensembl location database
                        ref_start = intron_values[1][0]; ref_stop = intron_values[1][1]
                        if ((start >= ref_start) and (stop <= ref_stop)) or ((start >= ref_start) and (start <= ref_stop)) or ((stop >= ref_start) and (stop <= ref_stop)):
                                old_index = intron_values[0]
                                if last_added == 'intron': last_added = 'intron'
                                elif last_added == 'null': last_added = 'exon'  ###utr is the default assigned at the top
                                final_exon_included = 'no'
                                #if probeset == '2948539': print probeset,'b',last_added,old_index,new_index,index2
                temp_probeset_exon_id=[]; temp_probeset_annot=[]; intron_retention_found = 'no'
                if new_key in intron_retention_db:
                    for exon_data in intron_retention_db[new_key]:
                        t_start = exon_data[0]; t_stop = exon_data[1]; ed = exon_data[2]
                        if ((start >= t_start) and (stop <= t_stop)):
                            temp_probeset_exon_id.append(ed.ExonID()); temp_probeset_annot.append(ed.Constitutive())
                    if len(temp_probeset_exon_id)>0:
                        temp_probeset_exon_id = unique.unique(temp_probeset_exon_id); temp_probeset_annot = unique.unique(temp_probeset_annot)
                        ens_exonid_values = convertListString(temp_probeset_exon_id,'|'); ens_constitutive = convertListString(temp_probeset_annot,'|')
                        #ned = EnsemblImport.ProbesetAnnotation(ens_exonid_values, ens_constitutive, 0, 'intron-retention', '')
                        splice_event_value = 'intron-retention' 
                        #new_exon_data = exon,ned; intron_retention_found = 'yes'
                #if probeset == '2948539': print  probeset,'c',last_added,old_index,index2
                ned = EnsemblImport.ProbesetAnnotation(ens_exonid_values, ens_constitutive, region_value, splice_event_value, splice_junction_values,exon_start,exon_stop)
                new_exon_data = exon,ned
                if old_index == 0: ###Means exons that are 5' to the ensembl reference exons and not in an intron
                    #if probeset == '2948539': print  probeset,'dU',last_added,old_index,new_index,index2
                    exon_info = 'U'+str(old_index)+'-'+str(index2) ###index2 will be 0 the first time by default
                    exon_info = exon_info,new_exon_data; exon_list.append(exon_info)
                    last_added = 'utr'; index2 += 1
                    #if probeset == '2948539':print probeset,exon_info, 'xl'
                else:
                    if last_added == 'exon':
                        #if probeset == '2948539': print  probeset,'eIU',last_added,old_index,new_index,index2
                        if final_exon_included == 'no':
                            index2 = 1
                            exon_info = 'I'+str(old_index)+'-'+str(index2)
                            exon_info = exon_info,new_exon_data; exon_list.append(exon_info)
                            last_added = 'intron'; index2 += 1
                            #if probeset == '2948539':print probeset,exon_info
                        else:
                            index2 = 1
                            exon_info = 'U'+str(old_index)+'-'+str(index2)
                            exon_info = exon_info,new_exon_data; exon_list.append(exon_info)
                            last_added = 'utr'; index2 += 1
                            #if probeset == '2948539':print probeset,exon_info
                    elif last_added == 'intron':
                        #if probeset == '2948539': print  probeset,'fI',last_added,old_index,new_index,index2
                        exon_info = 'I'+str(old_index)+'-'+str(index2)
                        exon_info = exon_info,new_exon_data; exon_list.append(exon_info)
                        last_added = 'intron'; index2 += 1
                        #if probeset == '2948539':print probeset,exon_info
                    elif last_added == 'utr':
                        ### Since the probeset is not in an exon and the last added was a UTR
                        ### could be in the 3'UTR (if in the 5'UTR it would have an old_index = 0) or
                        ### could be an intron, if no probesets aligned to the first exon
                        if final_exon_included == 'no': # alternative: old_index == 1: ###thus it is really in the first intron: example is 2948539
                            exon_info = 'I'+str(old_index)+'-'+str(index2)
                            exon_info = exon_info,new_exon_data; exon_list.append(exon_info)
                            last_added = 'intron'; index2 += 1
                        else:
                            exon_info = 'U'+str(old_index)+'-'+str(index2)
                            exon_info = exon_info,new_exon_data; exon_list.append(exon_info)
                            last_added = 'utr'; index2 += 1
                        #if probeset == '2948539':print probeset,exon_info;kill
            if len(exon_list)>0: exon_location2[key] = exon_list
    #kill
    #print x,p

    #exportProbesetAlignments(probeset_aligments)
    
    for key in exon_location2:
        if key[0] == 'ENSG00000138231':
            print key
            for item in exon_location2[key]: print item

    for key in exon_location2:
        if key[0] == 'ENSG00000095794':
            print key
            for item in exon_location2[key]: print item             
    return exon_location2

def reorderEnsemblLinkedProbesets(ensembl_transcript_clusters,transcript_cluster_data,trans_annotation_db):
    print len(trans_annotation_db), 'entries in old trans_annotation_db'
    ensembl_probeset_db={}; probeset_gene_redundancy={}; gene_transcript_redundancy={}; y=0; x=0; n=0; l=0; k=0
    for key in ensembl_transcript_clusters:
        geneid = key[0]; chr = 'chr'+ key[1]; strand = key[2]
        for entry in ensembl_transcript_clusters[key]:
            ###The structure of the entries varies dependant on if there were multiple ensembl's associated with a single transcript_cluster_id
            try: entry[0][0]
            except TypeError: entry = [entry]
            for probeset_info in entry:
                ###Rebuild this database since structural inconsistencies exist
                start = str(probeset_info[0]); stop = str(probeset_info[1]); exon_class = probeset_info[2]; probeset_id = probeset_info[3]
                transcript_cluster_id = transcript_cluster_data[probeset_id][-1]
                probeset_data = [start,stop,probeset_id,exon_class,transcript_cluster_id]; y += 1
                try: ensembl_probeset_db[geneid,chr,strand].append(probeset_data)
                except KeyError: ensembl_probeset_db[geneid,chr,strand]=[probeset_data]
                try: probeset_gene_redundancy[probeset_id,start,stop].append((geneid,chr,strand))
                except KeyError: probeset_gene_redundancy[probeset_id,start,stop] = [(geneid,chr,strand)]
                try: gene_transcript_redundancy[geneid].append(transcript_cluster_id)
                except KeyError: gene_transcript_redundancy[geneid] = [transcript_cluster_id]
    ###Correct incorrect Ensembl assignments based on trans-splicing etc.
    probeset_gene_redundancy = eliminateRedundant(probeset_gene_redundancy)
    gene_transcript_redundancy = eliminateRedundant(gene_transcript_redundancy)
    
    print len(ensembl_probeset_db), 'entries in old ensembl_probeset_db' 
    print len(gene_transcript_redundancy), 'entries in old gene_transcript_redundancy' 
    ### Added this new code to determine which transcript clusters uniquely detect a single gene (exon-level)
    ### Note: this is a potentially lengthy step (added in version 2.0.5)
    valid_gene_to_cluster_annotations={}; gene_transcript_redundancy2={}
    for probeset_info in probeset_gene_redundancy:
        probeset = probeset_info[0];start = int(probeset_info[1]); stop = int(probeset_info[2])
        transcript_cluster_id = transcript_cluster_data[probeset][-1]
        for ensembl_group in probeset_gene_redundancy[probeset_info]:
            pos_ens,neg_ens = alignProbesetsToEnsembl([],[],start,stop,ensembl_group) ### Indicates that at least one probeset in the TC aligns certain Ensembl genes
            ens_gene = ensembl_group[0] 
            if len(pos_ens)>0:
                try: valid_gene_to_cluster_annotations[transcript_cluster_id].append(ens_gene)
                except Exception: valid_gene_to_cluster_annotations[transcript_cluster_id] = [ens_gene]
    valid_gene_to_cluster_annotations = eliminateRedundant(valid_gene_to_cluster_annotations)
    print len(valid_gene_to_cluster_annotations), 'Valid gene-to-transcript cluser assignments based on genomic position'
    
    ### Remove probeset-gene and transcript_cluster-gene associations not supported by exon evidence
    for tc in valid_gene_to_cluster_annotations:
        for gene in valid_gene_to_cluster_annotations[tc]:
            try: gene_transcript_redundancy2[gene].append(tc)
            except Exception: gene_transcript_redundancy2[gene] = [tc]

    del_probesets = {}
    for probeset_info in probeset_gene_redundancy:
        probeset = probeset_info[0]
        transcript_cluster_id = transcript_cluster_data[probeset][-1]
        if transcript_cluster_id in valid_gene_to_cluster_annotations: ### If not, don't change the existing relationships
            keep=[]
            for ensembl_group in probeset_gene_redundancy[probeset_info]:
                ens_gene = ensembl_group[0]
                if ens_gene in valid_gene_to_cluster_annotations[transcript_cluster_id]: keep.append(ensembl_group)
            probeset_gene_redundancy[probeset_info] = keep  ### Replace the existing with this filtered list    
        else: del_probesets[probeset_info] = []
    for pi in del_probesets: del probeset_gene_redundancy[pi]
    
    trans_annotation_db = valid_gene_to_cluster_annotations
    gene_transcript_redundancy = gene_transcript_redundancy2
    print len(trans_annotation_db), 'entries in new trans_annotation_db'
    print len(gene_transcript_redundancy), 'entries in new gene_transcript_redundancy'
    print len(probeset_gene_redundancy), 'entries in probeset_gene_redundancy'

    ensembl_probeset_db2 = {}
    for (geneid,chr,strand) in ensembl_probeset_db:
        for probeset_data in ensembl_probeset_db[(geneid,chr,strand)]:
            start,stop,probeset_id,exon_class,transcript_cluster_id = probeset_data
            try:
                if (geneid,chr,strand) in probeset_gene_redundancy[probeset_id,start,stop]:
                    try: ensembl_probeset_db2[(geneid,chr,strand)].append(probeset_data)
                    except Exception: ensembl_probeset_db2[(geneid,chr,strand)] = [probeset_data]
            except KeyError: null=[]
    ensembl_probeset_db = ensembl_probeset_db2

    print len(ensembl_probeset_db), 'entries in new ensembl_probeset_db' 
    ###First check for many transcript IDs associating with one Ensembl
    remove_transcripts_clusters={}
    for geneid in gene_transcript_redundancy:
        transcript_cluster_id_list = gene_transcript_redundancy[geneid]
        if len(transcript_cluster_id_list)>1:
            for transcript_cluster_id in transcript_cluster_id_list:
                if transcript_cluster_id in trans_annotation_db: ###Check the Affymetrix transcript-Ensembl associations
                    ensembl_list = trans_annotation_db[transcript_cluster_id]
                    #[] ['3890870', '3890907', '3890909', '3890911'] ENSG00000124209
                    if transcript_cluster_id in test_cluster: print ensembl_list,transcript_cluster_id_list,geneid                   
                    if geneid not in ensembl_list:
                        try: remove_transcripts_clusters[transcript_cluster_id].append(geneid)
                        except Exception: remove_transcripts_clusters[transcript_cluster_id]=[geneid]
 
    """
    ###Perform a multi-step more refined search and remove transcript_clusters identified from above, that have no annotated Ensembl gene               
    remove_probesets={}
    for probeset_info in probeset_gene_redundancy:
        probeset = probeset_info[0]
        start = int(probeset_info[1]); stop = int(probeset_info[2])
        transcript_cluster_id = transcript_cluster_data[probeset][-1]
        ###First check to see if any probeset in the transcript_cluster_id should be eliminated 
        if transcript_cluster_id in remove_transcripts_clusters:
            try: remove_probesets[ensembl_group].append(probeset)
            except KeyError: remove_probesets[ensembl_group] = [probeset]
        if len(probeset_gene_redundancy[probeset_info])>1:
            x += 1
            ###Grab the redundant Ensembl list for each probeset
            ensembl_list1 = probeset_gene_redundancy[probeset_info]
            ###Grab the Ensembl list aligning to the transcript_cluster annotations for that probeset
            try:ensembl_list2 = trans_annotation_db[transcript_cluster_id]
            except KeyError: ensembl_list2=[]
            pos_ens=[]; neg_ens=[]
            for ensembl_group in ensembl_list1:
                ensembl = ensembl_group[0]
                if ensembl in ensembl_list2: pos_ens.append(ensembl_group)
                else: neg_ens.append(ensembl_group)
            if len(pos_ens) == 1:
                n += 1
                for ensembl_group in neg_ens:
                    try: remove_probesets[ensembl_group].append(probeset)
                    except KeyError: remove_probesets[ensembl_group] = [probeset]
            else: ###These are probesets for where the code did not identify a 'best' ensembl
                ###get probeset location and check against each exon in the exon cluster list
                for ensembl_group in ensembl_list1:  ###exon_clusters is from EnsemblImport: exon_data = exon,(start,stop),[accessions],[types]
                    if ensembl_group in exon_clusters:
                        for exon_data in exon_clusters[ensembl_group]:
                            exon_start = exon_data[1][0]
                            exon_stop = exon_data[1][1]
                            if (start >= exon_start) and (stop <= exon_stop):
                                pos_ens.append(ensembl_group)
                        if len(pos_ens) == 0:
                            neg_ens.append(ensembl_group)   
                pos_ens = unique.unique(pos_ens); neg_ens = unique.unique(neg_ens)
                if len(pos_ens) == 1:
                    l += 1
                    for ensembl_group in neg_ens:
                        try: remove_probesets[ensembl_group].append(probeset)
                        except KeyError: remove_probesets[ensembl_group] = [probeset]
                else: ###If no method for differentiating, probeset fom the database
                    k += 1
                    for ensembl_group in ensembl_list1:
                        try: remove_probesets[ensembl_group].append(probeset)
                        except KeyError: remove_probesets[ensembl_group] = [probeset]
                        """
    remove_probesets={}
    for probeset_info in probeset_gene_redundancy:
        probeset = probeset_info[0]
        start = int(probeset_info[1]); stop = int(probeset_info[2])
        transcript_cluster_id = transcript_cluster_data[probeset][-1]
        ###Grab the redundant Ensembl list for each probeset
        ensembl_list1 = probeset_gene_redundancy[probeset_info]
        ###First check to see if any probeset in the transcript_cluster_id should be eliminated
        if transcript_cluster_id in remove_transcripts_clusters:
            remove_genes = remove_transcripts_clusters[transcript_cluster_id] ### Why is this here
            pos_ens=[]; neg_ens=[]
            for ensembl_group in ensembl_list1: ###Ensembl group cooresponds to the exon_cluster dictionary key
                pos_ens,neg_ens = alignProbesetsToEnsembl(pos_ens,neg_ens,start,stop,ensembl_group)
            pos_ens = makeUnique(pos_ens); neg_ens = makeUnique(neg_ens)
            if len(pos_ens)!=1:
                ###if there are no probesets aligning to exons or probesets aligning to exons in multiple genes, remove these
                for ensembl_group in pos_ens:
                    try: remove_probesets[ensembl_group].append(probeset)
                    except KeyError: remove_probesets[ensembl_group] = [probeset]
                for ensembl_group in neg_ens:
                    try: remove_probesets[ensembl_group].append(probeset)
                    except KeyError: remove_probesets[ensembl_group] = [probeset]
            else:
                ###If a probeset is in an Ensembl exon (pos_ens), remove associations for probesets to ensembl IDs where the probeset is not in an exon for that gene
                for ensembl_group in neg_ens:
                    try: remove_probesets[ensembl_group].append(probeset)
                    except KeyError: remove_probesets[ensembl_group] = [probeset]
        elif len(probeset_gene_redundancy[probeset_info])>1:
            x += 1
            ###Grab the Ensembl list aligning to the transcript_cluster annotations for that probeset
            try: ensembl_list2 = trans_annotation_db[transcript_cluster_id]
            except KeyError: ensembl_list2=[]
            pos_ens1=[]; neg_ens1=[]
            for ensembl_group in ensembl_list1:
                ensembl = ensembl_group[0]
                if ensembl in ensembl_list2: pos_ens1.append(ensembl_group)
                else: neg_ens1.append(ensembl_group)
            pos_ens=[]; neg_ens=[]
            ###get probeset location and check against each exon in the exon cluster list
            for ensembl_group in ensembl_list1:  ###exon_clusters is from EnsemblImport: exon_data = exon,(start,stop),[accessions],[types]
                exon_found = 0
                if ensembl_group in exon_clusters:
                    for exon_data in exon_clusters[ensembl_group]:
                        exon_start = exon_data[1][0]
                        exon_stop = exon_data[1][1]
                        if (start >= exon_start) and (stop <= exon_stop):
                            pos_ens.append(ensembl_group); exon_found = 1
                    if exon_found == 0:
                        neg_ens.append(ensembl_group) 
            pos_ens = unique.unique(pos_ens); neg_ens = unique.unique(neg_ens)
            #if probeset == 3161639: print 'b',pos_ens, transcript_cluster_id, neg_ens;sys.exit()
            if len(pos_ens) == 1:
                l += 1
                for ensembl_group in neg_ens:
                    try: remove_probesets[ensembl_group].append(probeset)
                    except KeyError: remove_probesets[ensembl_group] = [probeset]
            elif len(pos_ens1) == 1:
                n += 1
                for ensembl_group in neg_ens1:
                    try: remove_probesets[ensembl_group].append(probeset)
                    except KeyError: remove_probesets[ensembl_group] = [probeset]  
            else: ###If no method for differentiating, probeset fom the database
                k += 1
                for ensembl_group in ensembl_list1:
                    try: remove_probesets[ensembl_group].append(probeset)
                    except KeyError: remove_probesets[ensembl_group] = [probeset]    
                    #if probeset == '3890871': print ensembl_group,ensembl_list1, probeset;kill
                 
    ensembl_probeset_db = removeRedundantProbesets(ensembl_probeset_db,remove_probesets)
    print "Number of Probesets linked to Ensembl entries:", y
    print "Number of Probesets occuring in multiple Ensembl genes:",x
    print "     removed by transcript_cluster editing:",n
    print "     removed by exon matching:",l
    print "     automatically removed (not associating with a clear single gene):",k
    print "     remaining redundant:",(x-(n+l+k))
    return ensembl_probeset_db

def alignProbesetsToEnsembl(pos_ens,neg_ens,start,stop,ensembl_group):
    try:
        k=0
        for exon_data in exon_clusters[ensembl_group]:
            exon_start = exon_data[1][0];exon_stop = exon_data[1][1]
            if (start >= exon_start) and (stop <= exon_stop): pos_ens.append(ensembl_group); k=1
        if k==0: neg_ens.append(ensembl_group)
        return pos_ens,neg_ens
    except KeyError: return pos_ens,neg_ens
    
def removeRedundantProbesets(ensembl_probeset_db,remove_probesets):
    ###Remove duplicate entries and record which probesets associate with multiple ensembl's (remove these ensembl's)
    #check_for_promiscuous_transcripts={}
    for ensembl_group in remove_probesets:
        new_probe_list=[]
        for probeset_data in ensembl_probeset_db[ensembl_group]:
            probeset = probeset_data[2]
            #transcript_cluster_id = transcript_cluster_data[probeset][-1]
            if probeset not in remove_probesets[ensembl_group]:
                new_probe_list.append(probeset_data)
                #try: check_for_promiscuous_transcripts[transcript_cluster_id].append(ensembl_group)
                #except KeyError: check_for_promiscuous_transcripts[transcript_cluster_id] = [ensembl_group]
        ensembl_probeset_db[ensembl_group] = new_probe_list

    ###Final Sanity check: see if probesets could have been added to more than one Ensmebl gene
    probeset_ensembl_db={}
    for ensembl_group in ensembl_probeset_db:
        for probeset_data in ensembl_probeset_db[ensembl_group]:
            probeset = probeset_data[2]
            try: probeset_ensembl_db[probeset].append(ensembl_group)
            except KeyError: probeset_ensembl_db[probeset] = [ensembl_group]
    for probeset in probeset_ensembl_db:
        if len(probeset_ensembl_db[probeset])>1: print probeset, probeset_ensembl_db[probeset]; kill
            
    """
    ###Remove Ensembl gene entries any transcript_cluster_id linking to multiple Ensembl genes
    ###This can occur if one set of probests links to one ensembl and another set (belong to the same transcript cluster ID), links to another
    ens_deleted = 0
    for transcript_cluster_id in check_for_promiscuous_transcripts:
        if len(check_for_promiscuous_transcripts[transcript_cluster_id])>1: ###Therefore, more than one ensembl gene is associated with that probeset
            for ensembl_group in check_for_promiscuous_transcripts[transcript_cluster_id]:
                ###remove these entries from the above ensembl probeset database
                try: del ensembl_probeset_db[ensembl_group]; ens_deleted += 1
                except KeyError: null = ''
    print ens_deleted,"ensembl gene entries deleted, for linking to probesets with multiple gene associations after all other filtering"
    """
    return ensembl_probeset_db

################# Select files for analysis and output results
def exportProbesetAlignments(probeset_aligments):
    ###These are probeset-transcript annotations direclty from Affymetrix, not aligned
    probeset_annotation_export = 'AltDatabase/' + species + '/'+arraytype+'/'+ species + '_probeset-exon-align-coord.txt'
    probeset_aligments = eliminateRedundant(probeset_aligments)
    
    fn=filepath(probeset_annotation_export); data = open(fn,'w')
    title = 'probeset_id'+'\t'+'genomic-start'+'\t'+'genomic-stop'+'\n'
    data.write(title)
    for probeset_id in probeset_aligments:
        for coordinates in probeset_aligments[probeset_id]:
            values = str(probeset_id) +'\t'+ str(coordinates[0]) +'\t'+ str(coordinates[1]) +'\n'
            data.write(values)
    data.close()
    
def exportEnsemblLinkedProbesets(arraytype, ensembl_probeset_db,species):
    exon_annotation_export = 'AltDatabase/'+species+'/'+arraytype+'/'+species + '_Ensembl_probesets.txt'
    subgeneviewer_export = 'AltDatabase/ensembl/'+species+'/'+species+'_SubGeneViewer_feature-data.txt'
    fn=filepath(exon_annotation_export); data = open(fn,'w')
    fn2=filepath(subgeneviewer_export); data2 = open(fn2,'w')
    title = ['probeset_id','exon_id','ensembl_gene_id','transcript_cluster_id','chromosome','strand','probeset_start','probeset_stop']
    title +=['affy_class','constitutive_probeset','ens_exon_ids','ens_constitutive_status','exon_region','exon-region-start(s)','exon-region-stop(s)','splice_events','splice_junctions']
    title2 =['probeset','gene-id','feature-id','region-id']
    title = string.join(title,'\t') + '\n'; title2 = string.join(title2,'\t') + '\n'
    data.write(title); data2.write(title2); y=0
    print 'len(ensembl_probeset_db)',len(ensembl_probeset_db)
    for key in ensembl_probeset_db:
        geneid = key[0]; chr = key[1]; strand = key[2]
        for probeset_data in ensembl_probeset_db[key]:
            exon_id,((start,stop,probeset_id,exon_class,transcript_clust),ed) = probeset_data
            try: constitutive = ed.ConstitutiveCall()
            except Exception: constitutive = 'no'
            if len(constitutive) == 0: constitutive = 'no'
            start = str(start); stop = str(stop)
            y+=1; ens_exon_list = ed.ExonID(); ens_annot_list = ed.Constitutive()
            if len(ed.AssociatedSplicingEvent())>0: constitutive = 'no'; ens_annot_list = '0' ### Re-set these if a splicing-event is associated
            values = [str(probeset_id),exon_id,geneid,str(transcript_clust),chr,strand,start,stop,exon_class,constitutive,ens_exon_list,ens_annot_list]
            values+= [str(ed.RegionNumber()),ed.ExonStart(),ed.ExonStop(),ed.AssociatedSplicingEvent(),ed.AssociatedSplicingJunctions()]
            region_num = ed.RegionNumber()
            if len(region_num)<1: b,r = string.split(exon_id,'-'); region_num = b+'-1'
            exon_id = string.replace(exon_id,'-','.'); region_num = string.replace(region_num,'-','.')
            try: values = string.join(values,'\t')+'\n';
            except TypeError:
                print exon_id
                print [probeset_id,exon_id,geneid,transcript_clust,chr,strand,start,stop,exon_class,constitutive,ens_exon_list,ens_annot_list]
                print [ed.RegionNumber(),ed.AssociatedSplicingEvent(),ed.AssociatedSplicingJunctions()];kill
            data.write(values)
            exon_regions = string.split(region_num,'|')
            for region in exon_regions:
                try:
                    if filter_sgv_output == 'yes': ### Can filter out probeset information for all probesets not linked to a probable splice event
                        if len(ed.AssociatedSplicingEvent())>1 or len(ens_exon_list)>1: proceed = 'yes'
                        else: proceed = 'no'
                    else: proceed = 'yes'
                except NameError: proceed = 'yes'
                if proceed == 'yes':
                    #print filter_sgv_output, len(ed.AssociatedSplicingEvent()),len(ens_exon_list),[ed.AssociatedSplicingEvent()],[ens_exon_list],exon_id,region
                    values2 =[str(probeset_id),geneid,exon_id,region]
                    values2 = string.join(values2,'\t')+'\n'
                    data2.write(values2)
    data.close()
    print y, "Probesets linked to Ensembl entries exported to text file:",exon_annotation_export

def eliminateRedundant(database):
    for key in database:
        try: list = makeUnique(database[key])
        except TypeError: list = unique.unique(database[key])
        list.sort()
        database[key] = list
    return database

def makeUnique(item):
    db1={}; list1=[]
    for i in item: db1[i]=[]
    for i in db1: list1.append(i)
    list1.sort()
    return list1

def testAffyAnnotationDownload(Species,array_type):
    global species; species = Species; global arraytype; arraytype = array_type
    checkDirectoryFiles()

def checkDirectoryFiles():
    """ Check to see if the probeset annotation file is present and otherwise download AltAnalyze hosted version"""
    dir = '/AltDatabase/'+species+'/'+arraytype
    probeset_annotation_file = getDirectoryFiles(dir)
    if probeset_annotation_file == None:
        filename = update.getFileLocations(species,arraytype)
        filename = dir[1:]+'/'+ filename
        update.downloadCurrentVersion(filename,arraytype,'csv')
        probeset_annotation_file = getDirectoryFiles(dir)
    if probeset_annotation_file == None: print 'No Affymetrix annotation file present for:', arraytype, species; sys.exit()
    else: print "Affymetrix annotation file found for", arraytype, species
    return filepath(probeset_annotation_file)
    
def getDirectoryFiles(dir):
    try: dir_list = read_directory(dir)  #send a sub_directory to a function to identify all files in a directory
    except Exception:
        export.createDirPath(filepath(dir[1:])) ### This directory needs to be created
        dir_list = read_directory(dir)
    probeset_annotation_file = None
    for data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = dir[1:]+'/'+data
        if '.transcript.' in affy_data_dir: transcript_annotation_file = affy_data_dir
        elif '.annoS' in affy_data_dir: probeset_transcript_file = affy_data_dir
        elif '.probeset' in affy_data_dir and '.csv' in affy_data_dir:
            if '.zip' not in affy_data_dir: probeset_annotation_file = affy_data_dir ###This file lets you grab the same info as in probeset_transcript_file, but along with mRNA associations
    return probeset_annotation_file

def reimportEnsemblProbesets(filename,probe_db=None,cs_db=None):
    fn=filepath(filename); x = 0
    #print [fn]
    if probe_db != None:
        probe_association_db=probe_db; constitutive_db=cs_db; constitutive_original_db={} ### grab these from a separate file
    else:
        probe_association_db={}; constitutive_db={}; constitutive_original_db={}
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        if x == 0:
            x=1;continue
        else:
            #t = string.split(line,'\t')
            #probeset_id=t[0];ensembl_gene_id=t[1];chr=t[2];strand=t[3];start=t[4];stop=t[5];exon_class=t[6]
            try: probeset_id, exon_id, ensembl_gene_id, transcript_cluster_id, chromosome, strand, probeset_start, probeset_stop, affy_class, constitutive_probeset, ens_exon_ids, exon_annotations,regionid,r_start,r_stop,splice_event,splice_junctions = string.split(data,'\t')
            except Exception: print data;kill
            probe_data = ensembl_gene_id,transcript_cluster_id,exon_id,ens_exon_ids,affy_class#,exon_annotations,constitutive_probeset
            proceed = 'yes'
            """
            if 'RNASeq' in filename:
                ### Restrict the analysis to exon RPKM or count data for constitutive calculation
                if 'exon' in biotypes:
                    if '-' in probeset_id: proceed = 'no'
            """
            if proceed == 'yes':
                probe_association_db[probeset_id] = probe_data
                if constitutive_probeset == 'yes':
                    try: constitutive_db[ensembl_gene_id].append(probeset_id)
                    except KeyError: constitutive_db[ensembl_gene_id] = [probeset_id]
                else: ### There was a bug here that resulted in no entries added (AltAnalyze version 1.15) --- because constitutive selection options have changed, should not have been an issue
                    try: constitutive_original_db[ensembl_gene_id].append(probeset_id)
                    except KeyError: constitutive_original_db[ensembl_gene_id] = [probeset_id]
            x+=1  
     ###If no constitutive probesets for a gene as a result of additional filtering (removing all probesets associated with a splice event), add these back
    for gene in constitutive_original_db:
        if gene not in constitutive_db: constitutive_db[gene] = constitutive_original_db[gene]
    constitutive_original_db={}

    if 'RNASeq' in filename: id_name = 'junction IDs'
    else: id_name = 'array IDs'
    print len(constitutive_db), "constitutive genes and", len(probe_association_db), id_name, "imported out of", x,"lines."
    return probe_association_db, constitutive_db

def reimportEnsemblProbesetsForSeqExtraction(filename,filter_type,filter_db):
    fn=filepath(filename); x = 0; ensembl_exon_db={}
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        if x == 0: x=1;continue
        else:
            try: probeset_id, exon_id, ensembl_gene_id, transcript_cluster_id, chr, strand, probeset_start, probeset_stop, affy_class, constitutive_probeset, ens_exon_ids, exon_annotations,regionid,r_start,r_stop,splice_event,splice_junctions = string.split(data,'\t')
            except ValueError: t = string.split(data,'\t'); print t;kill
            ens_exon_ids = string.split(ens_exon_ids,'|')
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            ed = EnsemblImport.ExonStructureData(ensembl_gene_id, chr, strand, r_start, r_stop, constitutive_probeset, ens_exon_ids, [])
            ed.setAssociatedSplicingEvent(splice_event) ###Integrate splicing annotations to look for alternative promoters
            #ed.setAssociatedSplicingJunctions(splice_junctions)
            probe_data  = exon_id,((probeset_start,probeset_stop,probeset_id,affy_class,transcript_cluster_id),ed)
            if filter_type == 'null': proceed = 'yes'
            elif filter_type == 'junctions':
                ed = EnsemblImport.ProbesetAnnotation(exon_id, constitutive_probeset, regionid, splice_event, splice_junctions,r_start,r_stop)
                probe_data = probeset_id, strand, probeset_start, probeset_stop
                probe_data = probe_data,ed
                proceed = 'yes'
            elif filter_type == 'only-junctions':
                if '-' in exon_id and '.' in exon_id:
                    try: location = chr+':'+string.split(r_start,'|')[1]+'-'+string.split(r_stop,'|')[0]
                    except Exception: location = chr+':'+probeset_start+'-'+probeset_stop
                    try: pos = [int(string.split(r_start,'|')[1]),int(string.split(r_stop,'|')[0])]
                    except Exception:
                        try: pos = [int(r_start),int(r_stop)]
                        except Exception: pos = [int(probeset_start),int(probeset_stop)]
                    probe_data = probeset_id, pos,location
                    proceed = 'yes'
                else: proceed = 'no'
            elif filter_type == 'gene':
                if ensembl_gene_id in filter_db: proceed = 'yes'
                else: proceed = 'no'
            elif filter_type == 'gene-probesets': ### get all probesets for a query gene
                if ensembl_gene_id in filter_db:
                    proceed = 'yes'
                    if '-' in exon_id and '.' in exon_id: proceed = 'no' #junction
                    else:
                        exon_id = string.replace(exon_id,'-','.')
                        try: block,region = string.split(exon_id[1:],'.')
                        except Exception: print original_exon_id;sys.exit()
                        if '_' in region:
                            region = string.split(region,'_')[0]
                        ed = EnsemblImport.ProbesetAnnotation(exon_id, constitutive_probeset, regionid, splice_event, splice_junctions,r_start,r_stop)
                        probe_data = (int(block),int(region)),ed,probeset_id ### sort by exon number            
                else: proceed = 'no'
            elif filter_type == 'probeset':
                if probeset_id in filter_db: proceed = 'yes'
                else: proceed = 'no'
            if filter_type == 'probesets':
                if ':' in probeset_id: probeset_id = string.split(probeset_id,':')[1]
                if len(regionid)<1: regionid = exon_id
                ensembl_exon_db[probeset_id]=string.replace(regionid,'-','.')
            elif filter_type == 'junction-regions':
                if '.' in regionid and '-' in regionid: ### Indicates a junction probeset
                    try: ensembl_exon_db[ensembl_gene_id].append((probeset_id,regionid))
                    except KeyError: ensembl_exon_db[ensembl_gene_id]=[(probeset_id,regionid)]
            elif proceed == 'yes':
                if filter_type == 'gene':
                    if ensembl_gene_id not in ensembl_exon_db: ensembl_exon_db[ensembl_gene_id] = [probe_data]
                else:
                    try: ensembl_exon_db[ensembl_gene_id].append(probe_data)
                    except KeyError: ensembl_exon_db[ensembl_gene_id] = [probe_data]                
                
    print len(ensembl_exon_db), "critical genes parsed with annotated exon data..."
    return ensembl_exon_db

def getSplicingAnalysisProbesets(probeset_db,constitutive_db,annotate_db):
    splicing_analysis_db={}; count=0; count2=0; count3=0
    #ensembl_annotation_db[ensembl_gene_id] = symbol,description,mRNA_processing
    for probeset in probeset_db:
        ensembl_gene = probeset_db[probeset][0]
        if ensembl_gene in constitutive_db:
            try: splicing_analysis_db[ensembl_gene].append(probeset)
            except KeyError: splicing_analysis_db[ensembl_gene] = [probeset]
            count += 1
        else:
            if ensembl_gene in annotate_db: mRNA_processing = annotate_db[ensembl_gene][2]
            else: mRNA_processing = ''
            if mRNA_processing == 'RNA_processing/binding':
                try: splicing_analysis_db[ensembl_gene].append(probeset)
                except KeyError: splicing_analysis_db[ensembl_gene] = [probeset]
                count2 += 1
            else:
                count3 += 1
    #print 'Number of probesets with constitutive probes:',count
    #print 'Number of other mRNA processing probesets:',count2
    #print 'Number of probesets excluded:',count3
    return splicing_analysis_db
    
def getAnnotations(process_from_scratch,x,source_biotype,Species):
    """Annotate Affymetrix exon array data using files Ensembl data (sync'ed to genome release)."""
    ###     probeset_db[probeset] = gene,exon_number,ensembl_exon_annotation,ensembl_exon_id
    #NEW    probeset_db[probeset] = gene,transcluster,exon_id,ens_exon_ids,exon_annotations,constitutive
    ### NA  constitutive_db[gene] = [probeset]
    ###     annotate_db[gene] = definition, symbol,rna_processing
    global species; species = Species; global export_probeset_mRNA_associationsg; global biotypes
    global test; global test_cluster; global filter_sgv_output; global arraytype
    export_probeset_mRNA_associations = 'no'; biotypes = ''
    if source_biotype == 'junction': arraytype = 'junction'; source_biotype = 'mRNA'
    elif source_biotype == 'gene': arraytype = 'gene'; reimportEnsemblProbesetsForSeqExtraction = 'mRNA'
    elif 'RNASeq' in source_biotype: arraytype,database_root_dir = source_biotype; source_biotype = 'mRNA'
    else: arraytype = 'exon'
    filter_sgv_output = 'no'
    test = 'no'
    test_cluster = [3161519, 3161559, 3161561, 3161564, 3161566, 3161706, 3161710, 3161712, 2716656, 2475411]
    test_cluster = [2887449]
    partial_process = 'no'; status = 'null'
    if process_from_scratch == 'yes':
        if partial_process == 'no':
            #"""
            global ensembl_exon_db; global ensembl_exon_db; global exon_clusters
            global exon_region_db; global intron_retention_db; global intron_clusters; global ucsc_splicing_annot_db
            global constitutive_source; constitutive_source = x
            probeset_transcript_file = checkDirectoryFiles()
            ensembl_exon_db,ensembl_annot_db,exon_clusters,intron_clusters,exon_region_db,intron_retention_db,ucsc_splicing_annot_db,ens_transcript_db = EnsemblImport.getEnsemblAssociations(species,source_biotype,test)
            ensembl_probeset_db = getProbesetAssociations(probeset_transcript_file,ensembl_exon_db,ens_transcript_db,source_biotype)
            SubGeneViewerExport.reorganizeData(species) ### reads in the data from the external generated files
            status = 'ran'
    if (process_from_scratch == 'no') or (status == 'ran'):
            if source_biotype == 'ncRNA':
                probeset_db_mRNA,constitutive_db = reimportEnsemblProbesets('AltDatabase/'+species+'/'+arraytype+'/'+species+'_Ensembl_probesets.txt')
                probeset_db_ncRNA,constitutive_db = reimportEnsemblProbesets('AltDatabase/'+species+'/'+arraytype+'/'+species+'_Ensembl_probesets_ncRNA.txt')
                probeset_db = {}
                #print len(probeset_db_mRNA), len(probeset_db_ncRNA), len(probeset_db)
                for probeset in probeset_db_ncRNA:
                    if probeset not in probeset_db_mRNA: probeset_db[probeset] = probeset_db_ncRNA[probeset]
                probeset_db_mRNA={}; probeset_db_ncRNA={}
            else:
                filename = 'AltDatabase/'+species+'/'+arraytype+'/'+species+'_Ensembl_probesets.txt'
                if arraytype != 'RNASeq':
                    probeset_db,constitutive_db = reimportEnsemblProbesets(filename)
                else:
                    exon_standard_dir = string.replace(filename,'_probesets.txt','_exons.txt')
                    probeset_db,constitutive_db = reimportEnsemblProbesets(exon_standard_dir)
                    filename = string.replace(database_root_dir+filename,'_probesets.txt','_junctions.txt')
                    probeset_db,constitutive_db = reimportEnsemblProbesets(filename,probe_db=probeset_db,cs_db=constitutive_db) ### These only include exons and junctions detected from the experiment
            annotate_db = EnsemblImport.reimportEnsemblAnnotations(species)
            splicing_analysis_db = getSplicingAnalysisProbesets(probeset_db,constitutive_db,annotate_db)
            #print "Probeset database and Annotation database reimported"
            #print "STATs: probeset database:",len(probeset_db),"probesets imported"
            #print "       annotation database:",len(annotate_db),"genes imported"
    return probeset_db,annotate_db,constitutive_db,splicing_analysis_db

if __name__ == '__main__':

    y = 'Ensembl'
    z = 'custom'
    m = 'Mm'
    h = 'Hs'
    r = 'Rn'
    source_biotype = 'mRNA'
    Species = h
    process_from_scratch = 'yes'
    export_probeset_mRNA_associations = 'no'
    constitutive_source = z ###If 'Ensembl', program won't look at any evidence except for Ensembl. Thus, not ideal
    array_type='exon'
    #getJunctionComparisonsFromExport(Species,array_type); sys.exit()
    #testAffyAnnotationDownload(Species,array_type); sys.exit()
    probeset_db,annotate_db,constitutive_db,splicing_analysis_db = getAnnotations(process_from_scratch,constitutive_source,source_biotype,Species)
    sys.exit()
    
    """Annotate Affymetrix exon array data using files Ensembl data (sync'ed to genome release)."""
    ###     probeset_db[probeset] = gene,exon_number,ensembl_exon_annotation,ensembl_exon_id
    #NEW    probeset_db[probeset] = gene,transcluster,exon_id,ens_exon_ids,exon_annotations,constitutive
    ### NA  constitutive_db[gene] = [probeset]
    ###     annotate_db[gene] = definition, symbol,rna_processing
    global species; species = Species
    global test; global test_cluster
    global filter_sgv_output
    export_probeset_mRNA_associations = 'no'
    filter_sgv_output = 'yes'
    test = 'yes'
    test_cluster = ['7958644']
    meta_test_cluster = ["3061319","3061268"]#,"3455632","3258444","3965936","2611056","3864519","3018509","3707352","3404496","3251490"]
    #test_cluster = meta_test_cluster
    partial_process = 'no'; status = 'null'
    if process_from_scratch == 'yes':
        if partial_process == 'no':
            #"""
            global ensembl_exon_db; global ensembl_exon_db; global exon_clusters
            global exon_region_db; global intron_retention_db; global ucsc_splicing_annot_db
            probeset_transcript_file = checkDirectoryFiles()
            ensembl_exon_db,ensembl_annot_db,exon_clusters,intron_clusters,exon_region_db,intron_retention_db,ucsc_splicing_annot_db,ens_transcript_db = EnsemblImport.getEnsemblAssociations(species,source_biotype,test)
            #trans_annotation_db = ExonArrayAffyRules.getTranscriptAnnotation(transcript_annotation_file,species,test,test_cluster) ###used to associate transcript_cluster ensembl's
            #"""
            ensembl_probeset_db = getProbesetAssociations(probeset_transcript_file,ensembl_exon_db,ens_transcript_db,source_biotype)
            SubGeneViewerExport.reorganizeData(species) ### reads in the data from the external generated files
            status = 'ran'
    if (process_from_scratch == 'no') or (status == 'ran'):
            probeset_db,constitutive_db = reimportEnsemblProbesets('AltDatabase/'+species+'/'+arraytype+'/'+species+'_Ensembl_probesets.txt')
            annotate_db = EnsemblImport.reimportEnsemblAnnotations(species)
            splicing_analysis_db = getSplicingAnalysisProbesets(probeset_db,constitutive_db,annotate_db)
            print "Probeset database and Annotation database reimported"
            print "STATs: probeset database:",len(probeset_db),"probesets imported"
            print "       annotation database:",len(annotate_db),"genes imported"
            

            probeset_transcript_file = checkDirectoryFiles()
            ensembl_exon_db,ensembl_annot_db,exon_clusters,intron_clusters,exon_region_db,intron_retention_db,ucsc_splicing_annot_db,ens_transcript_db = EnsemblImport.getEnsemblAssociations(species,source_biotype,test)
            #trans_annotation_db = ExonArrayAffyRules.getTranscriptAnnotation(transcript_annotation_file,species,test,test_cluster) ###used to associate transcript_cluster ensembl's
            #"""
            ensembl_probeset_db = getProbesetAssociations(probeset_transcript_file,ensembl_exon_db,ens_transcript_db,source_biotype)
            SubGeneViewerExport.reorganizeData(species) ### reads in the data from the external generated files
            status = 'ran'
    if (process_from_scratch == 'no') or (status == 'ran'):
            probeset_db,constitutive_db = reimportEnsemblProbesets('AltDatabase/'+species+'/'+arraytype+'/'+species+'_Ensembl_probesets.txt')
            annotate_db = EnsemblImport.reimportEnsemblAnnotations(species)
            splicing_analysis_db = getSplicingAnalysisProbesets(probeset_db,constitutive_db,annotate_db)
            print "Probeset database and Annotation database reimported"
            print "STATs: probeset database:",len(probeset_db),"probesets imported"
            print "       annotation database:",len(annotate_db),"genes imported"
            
