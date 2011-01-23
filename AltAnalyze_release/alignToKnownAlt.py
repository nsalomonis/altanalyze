###alignToKnownAlt
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
import GO_parsing
import copy
import time

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def makeUnique(item):
    db1={}; list1=[]; k=0
    for i in item:
        try: db1[i]=[]
        except TypeError: db1[tuple(i)]=[]; k=1
    for i in db1:
        if k==0: list1.append(i)
        else: list1.append(list(i))
    list1.sort()
    return list1

def customDeepCopy(db):
    db2={}
    for i in db:
        for e in db[i]:
            try: db2[i].append(e)
            except KeyError: db2[i]=[e]
    return db2

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

################### Import exon coordinate/transcript data from BIOMART
def importEnsExonStructureData(species,ensembl_gene_coordinates,ensembl_annotations,exon_annotation_db):
    ucsc_gene_coordinates={}
    filename = 'AltDatabase/ucsc/'+species+'/knownAlt.txt'
    start_time = time.time()
    fn=filepath(filename); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        regionid,chr,start,stop,splice_event_call,null,strand = string.split(data,'\t')
        start = int(start)+1; stop = int(stop); chr = chr[3:] ###checked it out and all UCSC starts are -1 from the correspond Ensembl start
        try: ucsc_gene_coordinates[chr,start,stop,strand].append(splice_event_call)
        except KeyError: ucsc_gene_coordinates[chr,start,stop,strand] = [splice_event_call]

    ensembl_chr_coordinate_db={}
    for gene in ensembl_gene_coordinates:
        a = ensembl_gene_coordinates[gene]; a.sort()
        gene_start = a[0]; gene_stop = a[-1]
        chr,strand = ensembl_annotations[gene]
        if chr in ensembl_chr_coordinate_db:
            ensembl_gene_coordinates2 = ensembl_chr_coordinate_db[chr]
            ensembl_gene_coordinates2[(gene_start,gene_stop)] = gene,strand
        else:
            ensembl_gene_coordinates2={}; ensembl_gene_coordinates2[(gene_start,gene_stop)] = gene,strand
            ensembl_chr_coordinate_db[chr]=ensembl_gene_coordinates2

    ucsc_chr_coordinate_db={}
    for geneid in ucsc_gene_coordinates:
        chr,start,stop,strand = geneid
        if chr in ucsc_chr_coordinate_db:
            ucsc_gene_coordinates2 = ucsc_chr_coordinate_db[chr]
            ucsc_gene_coordinates2[(start,stop)] = geneid,strand
        else:
            ucsc_gene_coordinates2={}; ucsc_gene_coordinates2[(start,stop)] = geneid,strand
            ucsc_chr_coordinate_db[chr] = ucsc_gene_coordinates2

    ensembl_transcript_clusters,no_match_list = getChromosomalOveralap(ucsc_chr_coordinate_db,ensembl_chr_coordinate_db)

    ensembl_ucsc_splicing_event_db = {}
    for clusterid in ensembl_transcript_clusters:
        ens_geneids = ensembl_transcript_clusters[clusterid]
        if len(ens_geneids)==1: ###If a cluster ID associates with multiple Ensembl IDs
            ens_geneid = ens_geneids[0]
            annotations = ucsc_gene_coordinates[clusterid]
            try: ensembl_ucsc_splicing_event_db[ens_geneid].append((clusterid,annotations))
            except KeyError: ensembl_ucsc_splicing_event_db[ens_geneid] = [(clusterid,annotations)]


    ensembl_ucsc_splicing_annotations={}
    for ensembl in ensembl_ucsc_splicing_event_db:
        chr,strand = ensembl_annotations[ensembl]
        key = ensembl,chr,strand
        ###Look through each of the annotations (with coordinate info) for those that are specifically AltPromoters
        ###Their coordinates occur overlapping but before the exon, so we want to change the coordinates 
        for (clusterid,annotations) in ensembl_ucsc_splicing_event_db[ensembl]:
            new_coordinates = []
            if 'altPromoter' in annotations:
                chr,bp1,ep1,strand = clusterid
                if key in exon_annotation_db:
                    exon_info_ls = exon_annotation_db[key]
                    for exon_info in exon_info_ls:
                        bp2 = exon_info[0]; ep2 = exon_info[0]; add = 0  ### Changed ep2 to be the second object in the list (previously it was also the first) 4-5-08
                        if ((bp1 >= bp2)  and (ep2 >= bp1)) or ((ep1 >= bp2)  and (ep2 >= ep1)): add = 1 ###if the start or stop of the UCSC region is inside the Ensembl start and stop
                        elif ((bp2 >= bp1)  and (ep1 >= bp2)) or ((ep2 >= bp1)  and (ep1 >= ep2)): add = 1 ###opposite
                        if add == 1:
                            new_coordinates += [bp1,bp2,ep1,ep2] ###record all coordinates and take the extreme values
                new_coordinates.sort()
                if len(new_coordinates)>0:
                    new_start = new_coordinates[0]; new_stop = new_coordinates[-1]
                    clusterid = chr,new_start,new_stop,strand
            annotation_str = string.join(annotations,'|')
            ###replace with new or old information
            start = clusterid[1]; stop = clusterid[2]
            try: ensembl_ucsc_splicing_annotations[ensembl].append((start,stop,annotation_str))
            except KeyError: ensembl_ucsc_splicing_annotations[ensembl] = [(start,stop,annotation_str)]
    return ensembl_ucsc_splicing_annotations

def getChromosomalOveralap(ucsc_chr_db,ensembl_chr_db):
    print len(ucsc_chr_db),len(ensembl_chr_db); start_time = time.time()
    """Find transcript_clusters that have overlapping start positions with Ensembl gene start and end (based on first and last exons)"""
    ###exon_location[transcript_cluster_id,chr,strand] = [(start,stop,exon_type,probeset_id)]
    y = 0; l =0; ensembl_transcript_clusters={}; no_match_list=[]
    ###(bp1,ep1) = (47211632,47869699); (bp2,ep2)  =  (47216942, 47240877)
    for chr in ucsc_chr_db:
        ucsc_db = ucsc_chr_db[chr]
        try:
         for (bp1,ep1) in ucsc_db:
            #print (bp1,ep1)
            x = 0
            gene_clusterid,ucsc_strand = ucsc_db[(bp1,ep1)]
            try:
                ensembl_db = ensembl_chr_db[chr]
                for (bp2,ep2) in ensembl_db:
                    y += 1; ensembl,ens_strand = ensembl_db[(bp2,ep2)]
                    #print (bp1,ep1),(bp2,ep2);kill
                    if ucsc_strand == ens_strand:
                        ###if the two gene location ranges overlapping
                        ##########FORCE UCSC mRNA TO EXIST WITHIN THE SPACE OF ENSEMBL TO PREVENT TRANSCRIPT CLUSTER EXCLUSION IN ExonArrayEnsemblRules
                        add = 0
                        if (bp1 >= bp2) and (ep2>= ep1): add = 1 ###if the annotations reside within the gene's start and stop position
                        #if ((bp1 >= bp2)  and (ep2 >= bp1)) or ((ep1 >= bp2)  and (ep2 >= ep1)): add = 1 ###if the start or stop of the UCSC region is inside the Ensembl start and stop
                        #elif ((bp2 >= bp1)  and (ep1 >= bp2)) or ((ep2 >= bp1)  and (ep1 >= ep2)): add = 1 ###opposite
                        if add == 1:
                            #if (bp1 >= bp2) and (ep2>= ep1): a = ''
                            #else: print gene_clusterid,ensembl,bp1,bp2,ep1,ep2;kill
                            x = 1
                            try: ensembl_transcript_clusters[gene_clusterid].append(ensembl)
                            except KeyError: ensembl_transcript_clusters[gene_clusterid] = [ensembl]
                            l += 1
            except KeyError: null=[]#; print chr, 'not found'
            if x == 0: no_match_list.append(gene_clusterid)
        except ValueError:
            for y in ucsc_db: print y;kill
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "UCSC genes matched up to Ensembl in %d seconds" % time_diff            
    print "UCSC Transcript Clusters (or accession numbers) overlapping with Ensembl:",len(ensembl_transcript_clusters)
    print "With NO overlapp",len(no_match_list)
    return ensembl_transcript_clusters,no_match_list

if __name__ == '__main__':
    species = 'Hs'
    #test = 'yes'
    #test_gene = ['ENSG00000140153','ENSG00000075413']
    ensembl_ucsc_splicing_annotations = importEnsExonStructureData(species,ensembl_gene_coordinates,ensembl_annotations,exon_annotation_db)