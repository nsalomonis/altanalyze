###ExonArrayAffyRules
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
from stats_scripts import statistics
import math

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

################# Parse and Analyze Files

def getConstitutiveProbesets(filename):
    """This function assumes that the probesets used by Affymetrix to calculate gene expression
    represent constitutive features.  Examination of these annotations reveal discrepencies in this
    assumption (some constitutive probesets not used and others used instead).  Note: The previous statement
    is not correct since it assumes core means constitutive, which is not the case.
    """
    fn=filepath(filename)
    constutive={}; x = 0
    for line in open(fn,'rU').xreadlines():             
        data,null = string.split(line,'\n')
        data = string.replace(data,'"','')
        if data[0] != '#' and 'probeset' not in data:
            probeset_id, transcript_cluster_id, probeset_list, probe_count = string.split(data,'\t')
            try: probeset_list = string.split(probeset_list,' ')
            except ValueError: probeset_list = [probeset_list]
            constutive[transcript_cluster_id] = probeset_list
            x += len(probeset_list)
    print "Constitutive transcript clusters:",len(constutive),"and probesets:",x
    return constutive

def getTranscriptAnnotation(filename,Species,test_status,test_clusterid):
    """This function gathers transcript_cluster annotations.  Since Affymetrix transcript_clusters appear to
    often cover more than one transcipt or have inconsistent mappings with other resources be wary. """
    global species; species = Species; global test; global test_cluster; test = test_status; test_cluster=test_clusterid
    fn=filepath(filename);trans_annotation_db={};trans_annot_extended={}
    for line in open(fn,'rU').xreadlines():             
        data,null = string.split(line,'\n')
        if data[0] != '#' and 'seqname' not in data:
            tab_data = string.split(data[1:-1],'","'); tab_data2=[]; transcript_id = tab_data[0]
            try: strand = tab_data[3]
            except IndexError:
                continue ### Problems in the source file exist if this happens... mis-placed \n, skip line
            continue_analysis = 'no'
            if test=='yes': ###used to test the program for a single gene
                if transcript_id in test_cluster: continue_analysis='yes' 
            else: continue_analysis='yes'
            if  continue_analysis=='yes':
                gene_annot = tab_data[7]; gene_annot2 = tab_data[8]
                uniprot = tab_data[9]; unigene = tab_data[10]
                #start = tab_data[4]; stop = tab_data[5]; go_b = tab_data[11]; go_c = tab_data[12]; go_f = tab_data[13] 
                try: pathway = tab_data[14]; domain = tab_data[15]; family = tab_data[16]
                except IndexError: pathway = tab_data[14]; domain = tab_data[15]; family = ''
                try: gene_annot = string.split(gene_annot,' // '); symbol=gene_annot[1]; definition=gene_annot[2]
                except IndexError: ref_seq = ''; symbol = ''; definition = ''

                uniprot_list=[]; uniprot_data = string.split(uniprot,' /// ')
                for entry in uniprot_data:
                    if ' // ' in entry: other_id,uniprot_id = string.split(entry,' // '); uniprot_list.append(uniprot_id)
                unigene_list=[]; unigene_data = string.split(unigene,' /// ')
                for entry in unigene_data:
                    if ' // ' in entry: other_id,unigene_id,tissue_exp = string.split(entry,' // '); unigene_list.append(unigene_id)                
                ensembl_list=[]; ensembl_data = string.split(gene_annot2,'gene:ENSG')
                for entry in ensembl_data:
                    if entry[0:2] == '00': ###Ensembl id's will have 0000 following the ENSG
                        ensembl_data = string.split(entry,' '); ensembl_list.append('ENSG'+ensembl_data[0])
                trans_annotation_db[transcript_id] = definition, symbol, ensembl_list
                trans_annot_extended[transcript_id] = ensembl_list,unigene_list,uniprot_list #strand,go_b,go_c,go_f,pathway,domain,family

    exportTranscriptAnnotations(trans_annot_extended)
    return trans_annotation_db

def getProbesetAssociations(filename):
    #probeset_db[probeset] = affygene,exons,probe_type_call,ensembl
    fn=filepath(filename)
    probe_association_db={}; constitutive_ranking = {}; x = 0; exon_location={}; mRNA_associations = {}
    for line in open(fn,'rU').xreadlines():             
        data,null = string.split(line,'\n')
        data = string.replace(data,'"',''); y = 0
        t = string.split(data,',')
        if data[0] != '#' and 'probeset' in data:
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
                    if 'fl' == affy_headers[index]: fn = index
                    if 'mrna' == affy_headers[index]: mr = index
                    if 'est' == affy_headers[index]: es = index
                    if 'ensGene' == affy_headers[index]: eg = index
                    if 'mrna_assignment' == affy_headers[index]: ma = index
                    index += 1      
        elif data[0] != '#' and 'probeset' not in data:
            probeset_id=t[pi];exon_cluster_id=t[ei];transcript_cluster_id=t[tc]; chr=t[sn];strand=t[sd]; mrna_assignment=t[ma]
            start=int(t[st]);stop=int(t[sp]); exon_type=t[lv]; fl=int(t[fn]); mRNA=int(t[mr]); est=int(t[es]); ensembl=int(t[eg])
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            ###Extract out the mRNA identifiers from the mRNA_annotation field
            ensembl_ids, mRNA_ids = grabRNAIdentifiers(mrna_assignment)
            mRNA_associations[probeset_id] = ensembl_ids, mRNA_ids
            evidence_count = fl+mRNA
            if strand == "-": start = int(t[7]); stop = int(t[6])
            if exons_to_grab == 'core' and exon_type == 'core': y = 1; x += 1
            elif exons_to_grab == 'extended' and (exon_type == 'core' or exon_type == 'extended'): y = 1; x += 1
            elif exons_to_grab == 'full':  ### includes 'ambiguous' exon_type
                y = 1; x += 1
            if y == 1:
                probe_association_db[probeset_id]=transcript_cluster_id,probeset_id,exon_type
                if strand == "-":
                    new_start = stop; new_stop = start; start = new_start; stop = new_stop
                try: exon_location[transcript_cluster_id,chr,strand].append((start,stop,exon_type,probeset_id))
                except KeyError: exon_location[transcript_cluster_id,chr,strand] = [(start,stop,exon_type,probeset_id)]
                const_info = [ensembl,evidence_count,est,probeset_id]
                try: constitutive_ranking[transcript_cluster_id].append(const_info)
                except KeyError: constitutive_ranking[transcript_cluster_id] = [const_info]
    print "Probeset Associations Parsed"

    ###Export probeset to transcript annotations
    exportProbesetAnnotations(mRNA_associations)
    
    ###Optionally grab constitutive annotations based on Ensembl, full-length and EST evidence only    
    if constitutive_source != 'Affymetrix':
        print "Begining to assembl constitutive annotations (Based on Ensembl/FL/EST evidence)..."
        alt_constitutive_gene_db,const_probe_count = rankConstitutive(constitutive_ranking)
        print "Number of Constitutive genes:",len(alt_constitutive_gene_db),"Number of Constitutive probesets:",const_probe_count
    
    exon_location2 = annotateExons(exon_location)
            
    print "Selected",exons_to_grab,"Probesets:",x
    return probe_association_db,exon_location2,alt_constitutive_gene_db

def grabRNAIdentifiers(mrna_assignment):
    ensembl_ids=[]; mRNA_ids=[]
    mRNA_entries = string.split(mrna_assignment,' /// ')
    for entry in mRNA_entries:
        mRNA_info = string.split(entry,' // '); mrna_ac = mRNA_info[0]
        if 'ENS' in mrna_ac: ensembl_ids.append(mrna_ac)
        else:
            try: int(mrna_ac[-3:]); mRNA_ids.append(mrna_ac)
            except ValueError: continue
    ensembl_ids = unique.unique(ensembl_ids)
    mRNA_ids = unique.unique(mRNA_ids)
    return ensembl_ids, mRNA_ids

def rankConstitutive(constitutive_ranking):
    constitutive_gene_db = {}
    for key in constitutive_ranking:
        c = constitutive_ranking[key]
        ens_list=[]; fl_list=[]; est_list=[]
        for entry in c:
            ens = entry[0]; fl = entry[1]; est = entry[2]; probeset = entry[3]
            ens_list.append(ens); fl_list.append(fl); est_list.append(est)
        ens_list.sort();ens_list.reverse(); fl_list.sort(); fl_list.reverse(); est_list.sort(); est_list.reverse()
        top_ens = ens_list[0]; top_fl = fl_list[0]; top_est = est_list[0]
        ###Remove EST only gene entries and entries with no mRNA or EST alignment variation
        if (top_ens == ens_list[-1]) and (top_fl == fl_list[-1]):
            if (top_est != est_list[-1]) and ((top_fl > 0) or (top_ens > 0)):
                for entry in c:
                    if entry[2] == top_est:
                        const_probeset = entry[3]
                        try:constitutive_gene_db[key].append(const_probeset)
                        except KeyError:constitutive_gene_db[key] = [const_probeset]
            else: continue
        elif top_ens > 1 and top_fl > 1:
            for entry in c:
                if entry[0] == top_ens:
                    const_probeset = entry[3]
                    try:constitutive_gene_db[key].append(const_probeset)
                    except KeyError:constitutive_gene_db[key] = [const_probeset]
        elif top_fl > 1: 
            for entry in c:
                if entry[1] == top_fl:
                    const_probeset = entry[3]
                    try:constitutive_gene_db[key].append(const_probeset)
                    except KeyError:constitutive_gene_db[key] = [const_probeset]
        elif top_ens > 1:
            for entry in c:
                if entry[0] == top_ens:
                    const_probeset = entry[3]
                    try:constitutive_gene_db[key].append(const_probeset)
                    except KeyError:constitutive_gene_db[key] = [const_probeset]
                    
    n=0; constitutive_probe_db={}
    for key in constitutive_gene_db:
        for entry in constitutive_gene_db[key]:
            n+=1
    return constitutive_gene_db,n

def annotateExons(exon_location):
    ###Annotate Affymetrix probesets independent of other annotations.  A problem with this method
    ###Is that it fails to recognize distance probesets that probe different regions of the same exon.

    for key in exon_location:
        exon_location[key].sort()
        strand = key[2]
        if strand == "-":
            exon_location[key].reverse()
            
    alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','x','y','z']
    exon_location2={}
    for key in exon_location:
        index = 1; index2 = 0; exon_list=[]; y = 0
        for exon in exon_location[key]:
            if y == 0:
                exon_info = 'E'+str(index),exon
                exon_list.append(exon_info); y = 1; last_start = exon[0]; last_stop = exon[1]; index += 1; index2 = 0  
            elif y == 1:
                current_start = exon[0]; current_stop = exon[1]
                if (abs(last_stop - current_start) < 20) or (abs(last_start - current_start) < 20) or (abs(current_stop - last_stop) < 20): ###Therefore, the two exons are < 20bp away
                    exon_info = 'E'+str(index-1)+alphabet[index2],exon
                    exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index2 += 1
                else:
                    exon_info = 'E'+str(index),exon
                    exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index += 1; index2 = 0                
        exon_location2[key] = exon_list
    """for key in exon_location2:
        if key[0] == '3242353':
            print key, exon_location2[key]"""
    return exon_location2

################# Select files for analysis and output results
def getDirectoryFiles():
    dir = '/AltDatabase/'+species+'/exon'
    dir_list = read_directory(dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = dir[1:]+'/'+data
        if 'transcript-annot' in affy_data_dir:
            transcript_annotation_file = affy_data_dir
        elif '.annot' in affy_data_dir:
            probeset_transcript_file = affy_data_dir
        elif '.probeset' in affy_data_dir:
            probeset_annotation_file = affy_data_dir
    return probeset_transcript_file,probeset_annotation_file,transcript_annotation_file

"""
def getDirectoryFiles(dir,exons_to_grab):
    import_dir = '/AltDatabase/'+species+'/exon'
    probeset_annotation_file=''
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = dir[1:]+'/'+data
        if 'transcript-annot' in affy_data_dir:
            transcript_annotation_file = affy_data_dir
        elif 'annot.hg' in affy_data_dir:
            probeset_transcript_file = affy_data_dir
        elif 'annot.hg' in affy_data_dir:
            probeset_annotation_file = affy_data_dir
        if exons_to_grab in affy_data_dir: ###this is the choosen constitutive list
            constitutive_probe_file = affy_data_dir
    return transcript_annotation_file,probeset_transcript_file,probeset_annotation_file,constitutive_probe_file
"""
def eliminateRedundant(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def getAnnotations(exons_to_grab,constitutive_source,process_from_scratch):
    """Annotate Affymetrix exon array data using files provided at http://www.Affymetrix.com. Filter these annotations based on
    the choice of 'core', 'extended', or 'full' annotations."""
    probeset_transcript_file,probeset_annotation_file,transcript_annotation_file = getDirectoryFiles()
    probe_association_db, exon_location_db, alt_constitutive_gene_db = getProbesetAssociations(probeset_annotation_file)
    trans_annotation_db = getTranscriptAnnotation(transcript_annotation_file)
    constitutive_db = getConstitutiveProbesets(constitutive_probe_file)
    if constitutive_source != 'Affymetrix':
        return probe_association_db,alt_constitutive_gene_db,exon_location_db, trans_annotation_db, trans_annot_extended
    else:
        return probe_association_db,constitutive_db,exon_location_db, trans_annotation_db, trans_annot_extended

def exportProbesetAnnotations(mRNA_associations):
    probeset_annotation_export = 'AltDatabase/ensembl/' + species + '/'+ species + '_probeset-mRNA_annot.txt'
    fn=filepath(probeset_annotation_export); data = open(fn,'w')
    title = 'probeset_id'+'\t'+'ensembl_transcript_ids'+'\t'+'mRNA_accession_ids'+'\n'
    data.write(title); y=0
    for probeset_id in mRNA_associations:
        ensembl_ids, mRNA_ids = mRNA_associations[probeset_id]
        if len(ensembl_ids)>0 or len(mRNA_ids)>0:
            ensembl_ids = string.join(ensembl_ids,','); mRNA_ids = string.join(mRNA_ids,',')
            values = probeset_id +'\t'+ ensembl_ids +'\t'+ mRNA_ids +'\n'
            data.write(values); y+=1
    data.close()
    print y, "Probesets linked to mRNA accession numbers exported to text file:",probeset_annotation_export

def exportTranscriptAnnotations(transcript_annotation_db):
    transcript_annotation_export = 'AltDatabase/ensembl/' + species + '/'+ species + '_Affygene-external_annot.txt'
    fn=filepath(transcript_annotation_export); data = open(fn,'w')
    title = 'transcript_cluster_id'+'\t'+'ensembl_transcript_ids'+'\t'+'mRNA_accession_ids'+'\n'
    data.write(title); y=0
    for transcript_cluster_id in transcript_annotation_db:
        ensembl_list,unigene_list,uniprot_list = transcript_annotation_db[transcript_cluster_id]
        if len(ensembl_list)>0 or len(unigene_list)>0 or len(uniprot_list)>0:
            ensembl_ids = string.join(ensembl_list,','); unigene_ids = string.join(unigene_list,','); uniprot_ids = string.join(uniprot_list,',')
            values = transcript_cluster_id +'\t'+ ensembl_ids +'\t'+ unigene_ids +'\t'+ uniprot_ids +'\n'
            data.write(values); y+=1
    data.close()
    print y, "Transcript clusters linked to other gene annotations to text file:",transcript_annotation_export    

if __name__ == '__main__':
    species = 'Hs'
    exons_to_grab = "core"
    constitutive_source = 'custom'
    process_from_scratch = 'yes'
    getAnnotations(exons_to_grab,constitutive_source,process_from_scratch)