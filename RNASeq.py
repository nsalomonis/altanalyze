###RNASeq
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

import sys, string, os
from stats_scripts import statistics
import math
import os.path
import unique
import update
import copy
import time
import export
from build_scripts import EnsemblImport; reload(EnsemblImport)
try: from build_scripts import JunctionArrayEnsemblRules
except Exception: pass ### occurs with circular imports
try: from build_scripts import JunctionArray; reload(JunctionArray)
except Exception: pass ### occurs with circular imports
try: from build_scripts import ExonArrayEnsemblRules
except Exception: pass ### occurs with circular imports
import multiprocessing
import logging
import traceback
import warnings
import bisect
from visualization_scripts import clustering; reload(clustering)
try:
    import scipy
    import scipy.cluster.hierarchy as sch
    import scipy.spatial.distance as dist
except Exception: pass
try: import numpy
except Exception: pass

LegacyMode = True

try:
    from scipy import average as Average
    from  scipy import stats
except Exception:
    try: from statistics import avg as Average
    except Exception: pass ### occurs with circular imports
    
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list_clean=[]
    dir_list = unique.read_directory(sub_dir)
    for filepath in dir_list:
        if 'log.txt' not in filepath and '.log' not in filepath:
            dir_list_clean.append(filepath)
    return dir_list_clean

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

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

######### Below code deals with building the AltDatabase #########
def collapseNoveExonBoundaries(novel_exon_coordinates,dataset_dir):
    """ Merge exon predictions based on junction measurments from TopHat. The predicted exons are
    bound by the identified splice site and the consensus length of reads in that sample"""
    
    dataset_dir = string.replace(dataset_dir,'exp.','ExpressionInput/novel.')
    
    export_data,status = AppendOrWrite(dataset_dir) ### Export all novel exons
    if status == 'not found':
        export_data.write('GeneID\tStrand\tExonID\tCoordinates\n')
    novel_gene_exon_db={}
    for (chr,coord) in novel_exon_coordinates:
        key = (chr,coord)
        ji,side,coord2 = novel_exon_coordinates[(chr,coord)]
        try:
            if side == 'left': ### left corresponds to the position of coord
                intron = string.split(string.split(ji.ExonRegionID(),'-')[1][:2],'.')[0]
            else:
                intron = string.split(string.split(ji.ExonRegionID(),'-'),'.')[0]
            ls = [coord,coord2]
            ls.sort() ### The order of this is variable
            if ji.Strand() == '-':
                coord2,coord = ls
            else: coord,coord2 = ls
            if 'I' in intron and ji.Novel() == 'side':
                #if 'ENSG00000221983' == ji.GeneID():
                try: novel_gene_exon_db[ji.GeneID(),ji.Strand(),intron].append((coord,coord2,ji,key,side))
                except Exception: novel_gene_exon_db[ji.GeneID(),ji.Strand(),intron] = [(coord,coord2,ji,key,side)]
        except Exception: pass
    
    outdatedExons={} ### merging novel exons, delete one of the two original
    for key in novel_gene_exon_db:
        firstNovel=True ### First putative novel exon coordinates examined for that gene
        novel_gene_exon_db[key].sort()
        if key[1]=='-':
            novel_gene_exon_db[key].reverse()

        for (c1,c2,ji,k,s) in novel_gene_exon_db[key]:
            if firstNovel==False:
                #print [c1,l2] #abs(c1-l2);sys.exit()
                ### see if the difference between the start position of the second exon is less than 300 nt away from the end of the last
                if abs(c2-l1) < 300 and os!=s: ### 80% of human exons are less than 200nt - PMID: 15217358
                    proceed = True
                    #if key[1]=='-':
                    if c2 in k:
                        novel_exon_coordinates[k] = ji,s,l1
                        outdatedExons[ok]=None ### merged out entry
                    elif l1 in ok:
                        novel_exon_coordinates[ok] = li,os,c2
                        outdatedExons[k]=None ### merged out entry
                    else:
                        proceed = False ### Hence, the two splice-site ends are pointing to two distinct versus one common exons
                    """
                    if c2 == 18683670 or l1 == 18683670:
                        print key,abs(c2-l1), c1, c2, l1, l2, li.ExonRegionID(), ji.ExonRegionID();
                        print k,novel_exon_coordinates[k]
                        print ok,novel_exon_coordinates[ok]
                    """
                    if proceed:
                        values = string.join([ji.GeneID(),ji.Strand(),key[2],ji.Chr()+':'+str(l1)+'-'+str(c2)],'\t')+'\n'
                        export_data.write(values)
                    
            ### For negative strand genes, c1 is larger than c2 but is the 5' begining of the exon
            l1,l2,li,ok,os = c1,c2,ji,k,s ### record the last entry
            firstNovel=False
            
    for key in outdatedExons: ### Delete the non-merged entry
        del novel_exon_coordinates[key]
    export_data.close()
    return novel_exon_coordinates

def exportNovelExonToBedCoordinates(species,novel_exon_coordinates,chr_status,searchChr=None):
    ### Export the novel exon coordinates based on those in the junction BED file to examine the differential expression of the predicted novel exon
    #bamToBed -i accepted_hits.bam -split| coverageBed -a stdin -b /home/databases/hESC_differentiation_exons.bed > day20_7B__exons-novel.bed
    bed_export_path = filepath('AltDatabase/'+species+'/RNASeq/chr/'+species + '_Ensembl_exons'+searchChr+'.bed')
    bed_data = open(bed_export_path,'w') ### Appends to existing file
    for (chr,coord) in novel_exon_coordinates:
        ji,side,coord2 = novel_exon_coordinates[(chr,coord)]
        if side == 'left': start,stop = coord,coord2
        if side == 'right': start,stop = coord2,coord
        try: gene = ji.GeneID()
        except Exception: gene = 'NA'
        if gene == None: gene = 'NA'
        if gene == None: gene = 'NA'
        if gene != 'NA': ### Including these has no benefit for AltAnalyze (just slows down alignment and piles up memory)
            if ji.Strand() == '-': stop,start=start,stop
            if chr_status == False:
                chr = string.replace(chr,'chr','') ### This will thus match up to the BAM files
            a = [start,stop]; a.sort(); start,stop = a
            bed_values = [chr,str(start),str(stop),gene,'0',str(ji.Strand())]
            bed_values = cleanUpLine(string.join(bed_values,'\t'))+'\n'
            bed_data.write(bed_values)
    bed_data.close()
    return bed_export_path

def moveBAMtoBEDFile(species,dataset_name,root_dir):
    bed_export_path = filepath('AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.bed')
    dataset_name = string.replace(dataset_name,'exp.','')
    new_fn = root_dir+'/BAMtoBED/'+species + '_'+dataset_name+'_exons.bed'
    new_fn = string.replace(new_fn,'.txt','')
    print 'Writing exon-level coordinates to BED file:'
    print new_fn
    catFiles(bed_export_path,'chr') ### concatenate the files ot the main AltDatabase directory then move
    export.customFileMove(bed_export_path,new_fn)
    return new_fn

def reformatExonFile(species,type,chr_status):
    if type == 'exon':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
        export_path = 'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.txt'
        ### Used by BEDTools to get counts per specific AltAnalyze exon region (should augment with de novo regions identified from junction analyses)
        bed_export_path = 'AltDatabase/'+species+'/RNASeq/chr/'+species + '_Ensembl_exons.bed'
        bed_data = export.ExportFile(bed_export_path)
    else:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_junction.txt'
        export_path = 'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_junctions.txt'
    print 'Writing',export_path
    export_data = export.ExportFile(export_path)
    fn=filepath(filename); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            x+=1
            export_title = ['AltAnalyzeID','exon_id','ensembl_gene_id','transcript_cluster_id','chromosome','strand','probeset_start','probeset_stop']
            export_title +=['affy_class','constitutive_probeset','ens_exon_ids','ens_constitutive_status','exon_region','exon-region-start(s)','exon-region-stop(s)','splice_events','splice_junctions']
            export_title = string.join(export_title,'\t')+'\n'; export_data.write(export_title)
        else:
            try: gene, exonid, chr, strand, start, stop, constitutive_call, ens_exon_ids, splice_events, splice_junctions = t
            except Exception: print t;kill
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention,
            if constitutive_call == 'yes': ens_constitutive_status = '1'
            else: ens_constitutive_status = '0'
            export_values = [gene+':'+exonid, exonid, gene, '', chr, strand, start, stop, 'known', constitutive_call, ens_exon_ids, ens_constitutive_status]
            export_values+= [exonid, start, stop, splice_events, splice_junctions]
            export_values = string.join(export_values,'\t')+'\n'; export_data.write(export_values)
            if type == 'exon':
                if chr_status == False:
                    chr = string.replace(chr,'chr','') ### This will thus match up to the BAM files
                bed_values = [chr,start,stop,gene+':'+exonid+'_'+ens_exon_ids,'0',strand]
                bed_values = string.join(bed_values,'\t')+'\n'; bed_data.write(bed_values)
    export_data.close()
    if type == 'exon': bed_data.close()

def importExonAnnotations(species,type,search_chr):
    if 'exon' in type:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
    else:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_junction.txt'

    fn=filepath(filename); x=0; exon_annotation_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            gene, exonid, chr, strand, start, stop, constitutive_call, ens_exon_ids, splice_events, splice_junctions = t; proceed = 'yes'
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if len(search_chr)>0:
                if chr != search_chr: proceed = 'no'
            if proceed == 'yes':
                if type == 'exon': start = int(start); stop = int(stop)
                ea = EnsemblImport.ExonAnnotationsSimple(chr, strand, start, stop, gene, ens_exon_ids, constitutive_call, exonid, splice_events, splice_junctions)
                if type == 'junction_coordinates': 
                    exon1_start,exon1_stop = string.split(start,'|')
                    exon2_start,exon2_stop = string.split(stop,'|')
                    if strand == '-':
                        exon1_stop,exon1_start = exon1_start,exon1_stop
                        exon2_stop,exon2_start = exon2_start,exon2_stop
                    #if gene == 'ENSMUSG00000027340': print chr,int(exon1_stop),int(exon2_start)
                    exon_annotation_db[chr,int(exon1_stop),int(exon2_start)]=ea              
                elif type == 'distal-exon':
                    exon_annotation_db[gene] = exonid
                else:
                    try: exon_annotation_db[gene].append(ea)
                    except KeyError: exon_annotation_db[gene]=[ea]
    return exon_annotation_db

def exportKnownJunctionComparisons(species):
    gene_junction_db = JunctionArrayEnsemblRules.importEnsemblUCSCAltJunctions(species,'standard')
    gene_intronjunction_db = JunctionArrayEnsemblRules.importEnsemblUCSCAltJunctions(species,'_intronic')
    for i in gene_intronjunction_db: gene_junction_db[i]=[]

    gene_junction_db2={}
    for (gene,critical_exon,incl_junction,excl_junction) in gene_junction_db:
        critical_exons = string.split(critical_exon,'|')
        for critical_exon in critical_exons:
            try: gene_junction_db2[gene,incl_junction,excl_junction].append(critical_exon)
            except Exception: gene_junction_db2[gene,incl_junction,excl_junction] = [critical_exon]
    gene_junction_db = gene_junction_db2; gene_junction_db2=[]
    
    junction_export = 'AltDatabase/' + species + '/RNASeq/'+ species + '_junction_comps.txt'
    fn=filepath(junction_export); data = open(fn,'w')
    print "Exporting",junction_export
    title = 'gene'+'\t'+'critical_exon'+'\t'+'exclusion_junction_region'+'\t'+'inclusion_junction_region'+'\t'+'exclusion_probeset'+'\t'+'inclusion_probeset'+'\t'+'data_source'+'\n'
    data.write(title); temp_list=[]
    for (gene,incl_junction,excl_junction) in gene_junction_db:
        critical_exons = unique.unique(gene_junction_db[(gene,incl_junction,excl_junction)])
        critical_exon = string.join(critical_exons,'|')
        temp_list.append(string.join([gene,critical_exon,excl_junction,incl_junction,gene+':'+excl_junction,gene+':'+incl_junction,'AltAnalyze'],'\t')+'\n')
    temp_list = unique.unique(temp_list)
    for i in temp_list: data.write(i)
    data.close()

def getExonAndJunctionSequences(species):
    export_exon_filename = 'AltDatabase/'+species+'/RNASeq/'+species+'_Ensembl_exons.txt'    
    ensembl_exon_db = ExonArrayEnsemblRules.reimportEnsemblProbesetsForSeqExtraction(export_exon_filename,'null',{})

    ### Import just the probeset region for mRNA alignment analysis
    analysis_type = ('region_only','get_sequence'); array_type = 'RNASeq'
    dir = 'AltDatabase/'+species+'/SequenceData/chr/'+species; gene_seq_filename = dir+'_gene-seq-2000_flank.fa'
    ensembl_exon_db = EnsemblImport.import_sequence_data(gene_seq_filename,ensembl_exon_db,species,analysis_type)

    critical_exon_file = 'AltDatabase/'+species+'/'+ array_type + '/' + array_type+'_critical-exon-seq.txt'
    getCriticalJunctionSequences(critical_exon_file,species,ensembl_exon_db)

    """    
    ### Import the full Ensembl exon sequence (not just the probeset region) for miRNA binding site analysis   
    analysis_type = 'get_sequence'; array_type = 'RNASeq'
    dir = 'AltDatabase/'+species+'/SequenceData/chr/'+species; gene_seq_filename = dir+'_gene-seq-2000_flank.fa'
    ensembl_exon_db = EnsemblImport.import_sequence_data(gene_seq_filename,ensembl_exon_db,species,analysis_type)
    """
    
    critical_exon_file = 'AltDatabase/'+species+'/'+ array_type + '/' + array_type+'_critical-exon-seq.txt'
    updateCriticalExonSequences(critical_exon_file, ensembl_exon_db)

def updateCriticalExonSequences(filename,ensembl_exon_db):
    exon_seq_db_filename = filename[:-4]+'_updated.txt'
    exonseq_data = export.ExportFile(exon_seq_db_filename)
    
    critical_exon_seq_db={}; null_count={}
    for gene in ensembl_exon_db:
        gene_exon_data={}
        for probe_data in ensembl_exon_db[gene]:
            exon_id,((probe_start,probe_stop,probeset_id,exon_class,transcript_clust),ed) = probe_data
            try: gene_exon_data[probeset_id] = ed.ExonSeq()
            except Exception: null_count[gene]=[] ### Occurs for non-chromosomal DNA (could also download this sequence though)
        if len(gene_exon_data)>0: critical_exon_seq_db[gene] = gene_exon_data
    print len(null_count),'genes not assigned sequenced (e.g.,non-chromosomal)'
    ensembl_exon_db=[]
    
    ### Export exon sequences 
    for gene in critical_exon_seq_db:
        gene_exon_data = critical_exon_seq_db[gene]
        for probeset in gene_exon_data:
            critical_exon_seq = gene_exon_data[probeset]
            values = [probeset,'',critical_exon_seq]
            values = string.join(values,'\t')+'\n'
            exonseq_data.write(values)  
    exonseq_data.close()
    print exon_seq_db_filename, 'exported....'

def getCriticalJunctionSequences(filename,species,ensembl_exon_db):
    ### Assemble and export junction sequences
    junction_seq_db_filename = string.replace(filename,'exon-seq','junction-seq')
    junctionseq_data = export.ExportFile(junction_seq_db_filename)

    critical_exon_seq_db={}; null_count={}
    for gene in ensembl_exon_db:
        gene_exon_data={}
        for probe_data in ensembl_exon_db[gene]:
            exon_id,((probe_start,probe_stop,probeset_id,exon_class,transcript_clust),ed) = probe_data
            try: gene_exon_data[probeset_id] = ed.ExonSeq()
            except Exception: null_count[gene]=[] ### Occurs for non-chromosomal DNA (could also download this sequence though)
        if len(gene_exon_data)>0: critical_exon_seq_db[gene] = gene_exon_data
    print len(null_count),'genes not assigned sequenced (e.g.,non-chromosomal)'
    ensembl_exon_db=[]
    
    junction_annotation_db = importExonAnnotations(species,'junction',[])
    for gene in junction_annotation_db:
        if gene in critical_exon_seq_db:
            gene_exon_data = critical_exon_seq_db[gene]
            for jd in junction_annotation_db[gene]:
                exon1,exon2=string.split(jd.ExonRegionIDs(),'-')
                p1=gene+':'+exon1
                p2=gene+':'+exon2
                p1_seq=gene_exon_data[p1][-15:]
                p2_seq=gene_exon_data[p2][:15]
                junction_seq = p1_seq+'|'+p2_seq
                junctionseq_data.write(gene+':'+jd.ExonRegionIDs()+'\t'+junction_seq+'\t\n')
    junctionseq_data.close()
    print junction_seq_db_filename, 'exported....'

def getEnsemblAssociations(species,data_type,test_status,force):
    ### Get UCSC associations (download databases if necessary)
    from build_scripts import UCSCImport
    mRNA_Type = 'mrna'; run_from_scratch = 'yes'
    export_all_associations = 'no' ### YES only for protein prediction analysis
    update.buildUCSCAnnoationFiles(species,mRNA_Type,export_all_associations,run_from_scratch,force)
    
    null = EnsemblImport.getEnsemblAssociations(species,data_type,test_status); null=[]
    reformatExonFile(species,'exon',True); reformatExonFile(species,'junction',True)
    exportKnownJunctionComparisons(species)
    getExonAndJunctionSequences(species)

######### Below code deals with user read alignment as opposed to building the AltDatabase #########

class ExonInfo:
    def __init__(self,start,unique_id,annotation):
        self.start = start; self.unique_id = unique_id; self.annotation = annotation
    def ReadStart(self): return self.start
    def UniqueID(self): return self.unique_id
    def Annotation(self): return self.annotation
    def setExonRegionData(self,rd): self.rd = rd
    def ExonRegionData(self): return self.rd
    def setExonRegionID(self,region_id): self.region_id = region_id
    def ExonRegionID(self): return self.region_id
    def setAlignmentRegion(self,region_type): self.region_type = region_type
    def AlignmentRegion(self): return self.region_type
    def __repr__(self): return "ExonData values"
    
class JunctionData:
    def __init__(self,chr,strand,exon1_stop,exon2_start,junction_id,biotype):
        self.chr = chr; self.strand = strand; self._chr = chr
        self.exon1_stop = exon1_stop; self.exon2_start = exon2_start
        self.junction_id = junction_id; self.biotype = biotype
        #self.reads = reads; self.condition = condition
        self.left_exon = None; self.right_exon = None; self.jd = None; self.gene_id = None
        self.trans_splicing = None
        self.splice_events=''
        self.splice_junctions=''
        self.seq_length=''
        self.uid = None
    def Chr(self): return self.chr
    def Strand(self): return self.strand
    def Exon1Stop(self): return self.exon1_stop
    def Exon2Start(self): return self.exon2_start
    def setExon1Stop(self,exon1_stop): self.exon1_stop = exon1_stop
    def setExon2Start(self,exon2_start): self.exon2_start = exon2_start
    def setSeqLength(self,seq_length): self.seq_length = seq_length
    def SeqLength(self): return self.seq_length
    def BioType(self): return self.biotype
    def checkExonPosition(self,exon_pos):
        if exon_pos == self.Exon1Stop(): return 'left'
        else: return 'right'
    ### These are used to report novel exon boundaries
    def setExon1Start(self,exon1_start): self.exon1_start = exon1_start
    def setExon2Stop(self,exon2_stop): self.exon2_stop = exon2_stop
    def Exon1Start(self): return self.exon1_start
    def Exon2Stop(self): return self.exon2_stop
    def Reads(self): return self.reads
    def JunctionID(self): return self.junction_id
    def Condition(self): return self.condition
    def setExonAnnotations(self,jd):
        self.jd = jd
        self.splice_events = jd.AssociatedSplicingEvent()
        self.splice_junctions = jd.AssociatedSplicingJunctions()
        self.exon_region = jd.ExonRegionIDs()
        self.exonid = jd.ExonID()
        self.gene_id = jd.GeneID()
        self.uid = jd.GeneID()+':'+jd.ExonRegionIDs()
    def ExonAnnotations(self): return self.jd
    def setLeftExonAnnotations(self,ld): self.gene_id,self.left_exon = ld
    def LeftExonAnnotations(self): return self.left_exon
    def setRightExonAnnotations(self,rd): self.secondary_geneid,self.right_exon = rd
    def RightExonAnnotations(self): return self.right_exon
    def setGeneID(self,geneid): self.gene_id = geneid
    def GeneID(self): return self.gene_id
    def setSecondaryGeneID(self,secondary_geneid): self.secondary_geneid = secondary_geneid
    def SecondaryGeneID(self): return self.secondary_geneid
    def setTransSplicing(self): self.trans_splicing = 'yes'
    def TransSplicing(self): return self.trans_splicing
    def SpliceSitesFound(self):
        if self.jd != None: sites_found = 'both'
        elif self.left_exon != None and self.right_exon != None: sites_found = 'both'
        elif self.left_exon != None: sites_found = 'left'
        elif self.right_exon != None: sites_found = 'right'
        else: sites_found = None
        return sites_found
    def setConstitutive(self,constitutive): self.constitutive = constitutive
    def Constitutive(self): return self.constitutive
    def setAssociatedSplicingEvent(self,splice_events): self.splice_events = splice_events
    def AssociatedSplicingEvent(self): return self.splice_events
    def setAssociatedSplicingJunctions(self,splice_junctions): self.splice_junctions = splice_junctions
    def AssociatedSplicingJunctions(self): return self.splice_junctions
    def setExonID(self,exonid): self.exonid = exonid
    def ExonID(self): return self.exonid
    def setExonRegionID(self,exon_region): self.exon_region = exon_region
    def ExonRegionID(self): return self.exon_region
    def setUniqueID(self,uid): self.uid = uid
    def UniqueID(self): return self.uid
    def setLeftExonRegionData(self,li): self.li = li
    def LeftExonRegionData(self): return self.li
    def setRightExonRegionData(self,ri): self.ri = ri
    def RightExonRegionData(self): return self.ri
    def setNovel(self, side): self.side = side
    def Novel(self): return self.side    
    def __repr__(self): return "JunctionData values"

def checkBEDFileFormat(bed_dir,root_dir):
    """ This method checks to see if the BED files (junction or exon) have 'chr' proceeding the chr number.
    It also checks to see if some files have two underscores and one has none or if double underscores are missing from all."""
    dir_list = read_directory(bed_dir)
    x=0
    break_now = False
    chr_present = False
    condition_db={}
    for filename in dir_list:
        fn=filepath(bed_dir+filename)
        #if ('.bed' in fn or '.BED' in fn): delim = 'r'
        delim = 'rU'
        if '.tab' in string.lower(filename) or '.bed' in string.lower(filename) or '.junction_quantification.txt' in string.lower(filename):
            condition_db[filename]=[]
            for line in open(fn,delim).xreadlines(): ### changed rU to r to remove \r effectively, rather than read as end-lines
                if line[0] == '#': x=0 ### BioScope
                elif x == 0: x=1 ###skip the first line
                elif x < 10: ### Only check the first 10 lines
                    if 'chr' in line: ### Need to look at multiple input formats (chr could be in t[0] or t[1])
                        chr_present = True
                    x+=1
                else:
                    break_now = True
                    break
            if break_now == True:
                break
    
    ### Check to see if exon.bed and junction.bed file names are propper or faulty (which will result in downstream errors)
    double_underscores=[]
    no_doubles=[]
    for condition in condition_db:
        if '__' in condition:
            double_underscores.append(condition)
        else:
            no_doubles.append(condition)
    
    exon_beds=[]
    junctions_beds=[] 
    if len(double_underscores)>0 and len(no_doubles)>0:
        ### Hence, a problem is likely due to inconsistent naming
        print 'The input files appear to have inconsistent naming. If both exon and junction sample data are present, make sure they are named propperly.'
        print 'For example: cancer1__exon.bed, cancer1__junction.bed (double underscore required to match these samples up)!'
        print 'Exiting AltAnalyze'; forceError
    elif len(no_doubles)>0:
        for condition in no_doubles:
            condition = string.lower(condition)
            if 'exon' in condition:
                exon_beds.append(condition)
            if 'junction' in condition:
                junctions_beds.append(condition)
        if len(exon_beds)>0 and len(junctions_beds)>0:
            print 'The input files appear to have inconsistent naming. If both exon and junction sample data are present, make sure they are named propperly.'
            print 'For example: cancer1__exon.bed, cancer1__junction.bed (double underscore required to match these samples up)!'
            print 'Exiting AltAnalyze'; forceError
    return chr_present

def getStrandMappingData(species):
    splicesite_db={}
    refExonCoordinateFile = unique.filepath('AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt')
    firstLine=True
    for line in open(refExonCoordinateFile,'rU').xreadlines():
        if firstLine: firstLine=False
        else:
            line = line.rstrip('\n')
            t = string.split(line,'\t'); #'gene', 'exon-id', 'chromosome', 'strand', 'exon-region-start(s)', 'exon-region-stop(s)', 'constitutive_call', 'ens_exon_ids', 'splice_events', 'splice_junctions'
            geneID, exon, chr, strand, start, stop = t[:6]
            splicesite_db[chr,int(start)]=strand
            splicesite_db[chr,int(stop)]=strand
    return splicesite_db

def importBEDFile(bed_dir,root_dir,species,normalize_feature_exp,getReads=False,searchChr=None,getBiotype=None,testImport=False,filteredJunctions=None):
    dir_list = read_directory(bed_dir)
    begin_time = time.time()
    if 'chr' not in searchChr:
        searchChr = 'chr'+searchChr
    condition_count_db={}; neg_count=0; pos_count=0; junction_db={}; biotypes={}; algorithms={}; exon_len_db={}; splicesite_db={}

    if testImport == 'yes': print "Reading user RNA-seq input data files"
    for filename in dir_list:
        count_db={}; rows=0
        fn=filepath(bed_dir+filename)
        condition = export.findFilename(fn)
        if '__' in condition:
            ### Allow multiple junction files per sample to be combined (e.g. canonical and non-canonical junction alignments)
            condition=string.split(condition,'__')[0]+filename[-4:]
        if ('.bed' in fn or '.BED' in fn or '.tab' in fn or '.TAB' in fn or '.junction_quantification.txt' in fn) and '._' not in condition:
            if ('.bed' in fn or '.BED' in fn): delim = 'r'
            else: delim = 'rU'
            ### The below code removes .txt if still in the filename along with .tab or .bed
            if '.tab' in fn: condition = string.replace(condition,'.txt','.tab')
            elif '.bed' in fn: condition = string.replace(condition,'.txt','.bed')
            if '.TAB' in fn: condition = string.replace(condition,'.txt','.TAB')
            elif '.BED' in fn: condition = string.replace(condition,'.txt','.BED')
            
            if testImport == 'yes': print "Reading the bed file", [fn], condition
            ### If the BED was manually created on a Mac, will neeed 'rU' - test this
            for line in open(fn,delim).xreadlines(): break
            if len(line)>500: delim = 'rU'
            for line in open(fn,delim).xreadlines(): ### changed rU to r to remove \r effectively, rather than read as end-lines
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                rows+=1
                if rows==1 or '#' == data[0]:
                    format_description = data
                    algorithm = 'Unknown'
                    if 'TopHat' in format_description: algorithm = 'TopHat'
                    elif 'HMMSplicer' in format_description: algorithm = 'HMMSplicer'
                    elif 'SpliceMap junctions' in format_description: algorithm = 'SpliceMap'
                    elif t[0] == 'E1': algorithm = 'BioScope-junction'
                    elif '# filterOrphanedMates=' in data or 'alignmentFilteringMode=' in data or '#number_of_mapped_reads=' in data:
                        algorithm = 'BioScope-exon'
                    elif '.junction_quantification.txt' in fn:
                        algorithm = 'TCGA format'
                        if 'barcode' in t: junction_position = 1
                        else: junction_position = 0
                    elif '.tab' in fn and len(t)==9:
                        try: start = float(t[1]) ### expect this to be a numerical coordinate
                        except Exception: continue
                        algorithm = 'STAR'
                        strand = '-' ### If no strand exists
                        rows=2 ### allows this first row to be processed
                        if len(splicesite_db)==0: ### get strand to pos info
                            splicesite_db = getStrandMappingData(species)   
                    if testImport == 'yes': print condition, algorithm

                if rows>1:
                    try:
                        if ':' in t[0]:
                            chr = string.split(t[0],':')[0]
                        else: chr = t[0]
                        if 'chr' not in chr:
                            chr = 'chr'+chr
                        if searchChr == chr or ('BioScope' in algorithm and searchChr == t[1]): proceed = True
                        elif searchChr == 'chrMT' and ('BioScope' not in algorithm):
                            if 'M' in chr and len(chr)<6: proceed = True ### If you don't have the length, any random thing with an M will get included
                            else: proceed = False
                        else: proceed = False
                    except IndexError:
                        print 'The input file:\n',filename
                        print 'is not formated as expected (format='+algorithm+').'
                        print 'search chromosome:',searchChr
                        print t; force_bad_exit

                    if proceed:
                        proceed = False
                        if '.tab' in fn or '.TAB' in fn:
                            ### Applies to non-BED format Junction and Exon inputs (BioScope)
                            if 'BioScope' in algorithm:
                                if algorithm == 'BioScope-exon': ### Not BED format
                                    chr,source,data_type,start,end,reads,strand,null,gene_info=t[:9]
                                    if 'chr' not in chr: chr = 'chr'+chr
                                    if data_type == 'exon': ### Can also be CDS
                                        gene_info,test,rpkm_info,null = string.split(gene_info,';')
                                        symbol = string.split(gene_info,' ')[-1]
                                        #refseq = string.split(transcript_info,' ')[-1]
                                        rpkm = string.split(rpkm_info,' ')[-1]
                                        #if normalize_feature_exp == 'RPKM': reads = rpkm ### The RPKM should be adjusted +1 counts, so don't use this
                                        biotype = 'exon'; biotypes[biotype]=[]
                                        exon1_stop,exon2_start = int(start),int(end); junction_id=''
                                        ### Adjust exon positions - not ideal but necessary. Needed as a result of exon regions overlapping by 1nt (due to build process)
                                        exon1_stop+=1; exon2_start-=1
                                        #if float(reads)>4 or getReads:
                                        proceed = True ### Added in version 2.0.9 to remove rare novel isoforms
                                        seq_length = abs(exon1_stop-exon2_start)
                                if algorithm == 'BioScope-junction':
                                    chr = t[1]; strand = t[2]; exon1_stop = int(t[4]); exon2_start = int(t[8]); count_paired = t[17]; count_single = t[19]; score=t[21]
                                    if 'chr' not in chr: chr = 'chr'+chr
                                    try: exon1_start = int(t[3]); exon2_stop = int(t[9])
                                    except Exception: pass ### If missing, these are not assigned
                                    reads = str(int(float(count_paired))+int(float(count_single))) ### Users will either have paired or single read (this uses either)
                                    biotype = 'junction'; biotypes[biotype]=[]; junction_id=''
                                    if float(reads)>4 or getReads: proceed = True ### Added in version 2.0.9 to remove rare novel isoforms
                                    seq_length = abs(float(exon1_stop-exon2_start))
                            if 'STAR' in algorithm:
                                chr = t[0]; exon1_stop = int(t[1])-1; exon2_start = int(t[2])+1; strand=''
                                if 'chr' not in chr: chr = 'chr'+chr
                                reads = str(int(t[7])+int(t[6]))
                                biotype = 'junction'; biotypes[biotype]=[]; junction_id=''
                                if float(reads)>4 or getReads: proceed = True ### Added in version 2.0.9 to remove rare novel isoforms
                                if (chr,exon1_stop) in splicesite_db:
                                    strand = splicesite_db[chr,exon1_stop]
                                elif (chr,exon2_start) in splicesite_db:
                                    strand = splicesite_db[chr,exon2_start]
                                #else: proceed = False
                                seq_length = abs(float(exon1_stop-exon2_start))
                                if strand == '-': ### switch the orientation of the positions
                                    exon1_stop,exon2_start=exon2_start,exon1_stop
                                exon1_start = exon1_stop; exon2_stop = exon2_start
                                #if 9996685==exon1_stop and 10002682==exon2_stop:
                                #print chr, strand, reads, exon1_stop, exon2_start,proceed;sys.exit()
                        else:
                            try:
                                if algorithm == 'TCGA format':
                                    coordinates = string.split(t[junction_position],',')
                                    try: chr,pos1,strand = string.split(coordinates[0],':')
                                    except Exception: print t;sys.exit()
                                    chr,pos2,strand = string.split(coordinates[1],':')
                                    if 'chr' not in chr: chr = 'chr'+chr
                                    pos2 = str(int(pos2)-1) ### This is the bed format conversion with exons of 0 length
                                    exon1_start, exon2_stop = pos1, pos2
                                    reads = t[junction_position+1]
                                    junction_id = t[junction_position]
                                    exon1_len=0; exon2_len=0
                                else: 
                                    ### Applies to BED format Junction input
                                    chr, exon1_start, exon2_stop, junction_id, reads, strand, null, null, null, null, lengths, null = t
                                    if 'chr' not in chr: chr = 'chr'+chr
                                    exon1_len,exon2_len=string.split(lengths,',')[:2]; exon1_len = int(exon1_len); exon2_len = int(exon2_len)
                                exon1_start = int(exon1_start); exon2_stop = int(exon2_stop)
                                biotype = 'junction'; biotypes[biotype]=[]
                                if strand == '-':
                                    if (exon1_len-exon2_len)==0: ### Kallisto-Splice directly reports these coordinates
                                        exon1_stop = exon1_start
                                        exon2_start = exon2_stop
                                    else:
                                        exon1_stop = exon1_start+exon1_len; exon2_start=exon2_stop-exon2_len+1
                                    ### Exons have the opposite order
                                    a = exon1_start,exon1_stop; b = exon2_start,exon2_stop
                                    exon1_stop,exon1_start = b; exon2_stop,exon2_start = a
                                else:
                                    if (exon1_len-exon2_len)==0: ### Kallisto-Splice directly reports these coordinates
                                        exon1_stop = exon1_start
                                        exon2_start= exon2_stop
                                    else:
                                        exon1_stop = exon1_start+exon1_len; exon2_start=exon2_stop-exon2_len+1
                                if float(reads)>4 or getReads: proceed = True
                                if algorithm == 'HMMSplicer':
                                    if '|junc=' in junction_id: reads = string.split(junction_id,'|junc=')[-1]
                                    else: proceed = False
                                if algorithm == 'SpliceMap':
                                    if ')' in junction_id and len(junction_id)>1: reads = string.split(junction_id,')')[0][1:]
                                    else: proceed = False
                                seq_length = abs(float(exon1_stop-exon2_start)) ### Junction distance
                            except Exception,e:
                                #print traceback.format_exc();sys.exit()
                                ### Applies to BED format exon input (BEDTools export)
                                # bamToBed -i accepted_hits.bam -split| coverageBed -a stdin -b /home/nsalomonis/databases/Mm_Ensembl_exons.bed > day0_8B__exons.bed
                                try: chr, start, end, exon_id, null, strand, reads, bp_coverage, bp_total, percent_coverage = t
                                except Exception:
                                    print 'The file',fn,'does not appear to be propperly formatted as input.'
                                    print t; force_exception
                                if 'chr' not in chr: chr = 'chr'+chr
                                algorithm = 'TopHat-exon'; biotype = 'exon'; biotypes[biotype]=[]
                                exon1_stop,exon2_start = int(start),int(end); junction_id=exon_id; seq_length = float(bp_total)
                                if seq_length == 0:
                                    seq_length = abs(float(exon1_stop-exon2_start))
                                ### Adjust exon positions - not ideal but necessary. Needed as a result of exon regions overlapping by 1nt (due to build process)
                                exon1_stop+=1; exon2_start-=1
                                #if float(reads)>4 or getReads: ### Added in version 2.0.9 to remove rare novel isoforms
                                proceed = True
                                #else: proceed = False

                        if proceed:
                            if 'chr' not in chr:
                                chr = 'chr'+chr ### Add the chromosome prefix
                            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
                            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
                            if strand == '+': pos_count+=1
                            else: neg_count+=1
                            if getReads and seq_length>0:
                                if getBiotype == biotype:
                                    if biotype == 'junction':
                                        ### We filtered for junctions>4 reads before, now we include all reads for expressed junctions
                                        if (chr,exon1_stop,exon2_start) in filteredJunctions:
                                            count_db[chr,exon1_stop,exon2_start] = reads
                                            try: exon_len_db[chr,exon1_stop,exon2_start] = seq_length
                                            except Exception: exon_len_db[chr,exon1_stop,exon2_start] = []
                                    else:
                                        count_db[chr,exon1_stop,exon2_start] = reads
                                        try: exon_len_db[chr,exon1_stop,exon2_start] = seq_length
                                        except Exception: exon_len_db[chr,exon1_stop,exon2_start] = []
                            elif seq_length>0:
                                if (chr,exon1_stop,exon2_start) not in junction_db:
                                    ji = JunctionData(chr,strand,exon1_stop,exon2_start,junction_id,biotype)
                                    junction_db[chr,exon1_stop,exon2_start] = ji
                                    try: ji.setSeqLength(seq_length) ### If RPKM imported or calculated
                                    except Exception: null=[]
                                    try: ji.setExon1Start(exon1_start);ji.setExon2Stop(exon2_stop)
                                    except Exception: null=[]
                                    key = chr,exon1_stop,exon2_start
                algorithms[algorithm]=[]
            if getReads:
                if condition in condition_count_db:
                    ### combine the data from the different files for the same sample junction alignments
                    count_db1 = condition_count_db[condition]
                    for key in count_db:
                        if key not in count_db1: count_db1[key] = count_db[key]
                        else:
                            combined_counts = int(count_db1[key])+int(count_db[key])
                            count_db1[key] = str(combined_counts)
                    condition_count_db[condition]=count_db1
                else:
                    try: condition_count_db[condition] = count_db
                    except Exception: null=[] ### Occurs for other text files in the directory that are not used for the analysis
                
    end_time = time.time()
    if testImport == 'yes': print 'Read coordinates imported in',int(end_time-begin_time),'seconds'

    if getReads:
        #print len(exon_len_db), getBiotype, 'read counts present for',algorithm
        return condition_count_db,exon_len_db,biotypes,algorithms
    else:
        if testImport == 'yes':
            if 'exon' not in biotypes and 'BioScope' not in algorithm:
                print len(junction_db),'junctions present in',algorithm,'format BED files.' # ('+str(pos_count),str(neg_count)+' by strand).'
            elif 'exon' in biotypes and 'BioScope' not in algorithm:
                print len(junction_db),'sequence identifiers present in input files.' 
            else: print len(junction_db),'sequence identifiers present in BioScope input files.'
        return junction_db,biotypes,algorithms

def importExonCoordinates(probeCoordinateFile,search_chr,getBiotype):
    probe_coordinate_db={}
    junction_db={}
    biotypes={}
    x=0
    fn=filepath(probeCoordinateFile)
    for line in open(fn,'rU').xreadlines(): ### changed rU to r to remove \r effectively, rather than read as end-lines
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            probe_id = t[0]; probeset_id=t[1]; chr=t[2]; strand=t[3]; start=t[4]; end=t[5]
            exon1_stop,exon2_start = int(start),int(end)
            seq_length = abs(float(exon1_stop-exon2_start))
            if 'chr' not in chr:
                chr = 'chr'+chr ### Add the chromosome prefix
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if search_chr == chr or search_chr == None:
                try: biotype = t[6]
                except Exception:
                    if seq_length>25:biotype = 'junction'
                    else: biotype = 'exon'
                if strand == '-':
                    exon1_stop,exon2_start = exon2_start, exon1_stop ### this is their actual 5' -> 3' orientation
                if biotype == 'junction':
                    exon1_start,exon2_stop = exon1_stop,exon2_start
                else:
                    exon1_stop+=1; exon2_start-=1
                biotypes[biotype]=[]
                if getBiotype == biotype or getBiotype == None:
                    ji = JunctionData(chr,strand,exon1_stop,exon2_start,probe_id,biotype)
                    junction_db[chr,exon1_stop,exon2_start] = ji
                    try: ji.setSeqLength(seq_length) ### If RPKM imported or calculated
                    except Exception: null=[]
                    try: ji.setExon1Start(exon1_start);ji.setExon2Stop(exon2_stop)
                    except Exception: null=[]
                    probe_coordinate_db[probe_id] = chr,exon1_stop,exon2_start ### Import the expression data for the correct chromosomes with these IDs

    return probe_coordinate_db, junction_db, biotypes
       
def importExpressionMatrix(exp_dir,root_dir,species,fl,getReads,search_chr=None,getBiotype=None):
    """ Non-RNA-Seq expression data (typically Affymetrix microarray) import and mapping to an external probe-coordinate database """
    begin_time = time.time()
            
    condition_count_db={}; neg_count=0; pos_count=0; algorithms={}; exon_len_db={}
    
    probe_coordinate_db, junction_db, biotypes = importExonCoordinates(fl.ExonMapFile(),search_chr,getBiotype)
    
    x=0
    fn=filepath(exp_dir)[:-1]
    condition = export.findFilename(fn)
    ### If the BED was manually created on a Mac, will neeed 'rU' - test this
    for line in open(fn,'rU').xreadlines(): ### changed rU to r to remove \r effectively, rather than read as end-lines
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if '#' == data[0]: None
        elif x==0:
            if 'block' in t:
                start_index = 7
            else:
                start_index = 1
            headers = t[start_index:]
            x=1
        else:
            proceed = 'yes' ### restrict by chromosome with minimum line parsing (unless we want counts instead)
            probe_id=t[0]
            if probe_id in probe_coordinate_db:
                key = probe_coordinate_db[probe_id]
                if getReads == 'no':
                    pass
                else:
                    expression_data = t[start_index:]
                    i=0
                    for sample in headers:
                        if sample in condition_count_db:
                            count_db = condition_count_db[sample]
                            count_db[key] = expression_data[i]
                            exon_len_db[key]=[]
                        else:
                            count_db={}
                            count_db[key] = expression_data[i]
                            condition_count_db[sample] = count_db
                            exon_len_db[key]=[]
                        i+=1

    algorithms['ProbeData']=[]
    end_time = time.time()
    if testImport == 'yes': print 'Probe data imported in',int(end_time-begin_time),'seconds'
    if getReads == 'yes':
        return condition_count_db,exon_len_db,biotypes,algorithms
    else:
        return junction_db,biotypes,algorithms

def adjustCounts(condition_count_db,exon_len_db):
    for key in exon_len_db:
        try:
            null=exon_len_db[key]
            for condition in condition_count_db:
                count_db = condition_count_db[condition]
                try: read_count = float(count_db[key])+1 ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
                except KeyError: read_count = 1 ###Was zero, but needs to be one for more realistic log2 fold calculations
                count_db[key] = str(read_count) ### Replace original counts with adjusted counts
        except Exception: null=[]
    return condition_count_db

def calculateRPKM(condition_count_db,exon_len_db,biotype_to_examine):
    """Determines the total number of reads in a sample and then calculates RPMK relative to a pre-determined junction length (60).
    60 was choosen, based on Illumina single-end read lengths of 35 (5 nt allowed overhand on either side of the junction)"""
    ### Get the total number of mapped reads
    mapped_reads={}
    for condition in condition_count_db:
        mapped_reads[condition]=0
        count_db = condition_count_db[condition]
        for key in count_db:
            read_count = count_db[key]
            mapped_reads[condition]+=float(read_count)
            
    ### Use the average_total_reads when no counts reported such that 0 counts are comparable            
    average_total_reads = 0
    for i in mapped_reads:
        average_total_reads+=mapped_reads[i]
        if testImport == 'yes':
            print 'condition:',i,'total reads:',mapped_reads[i]
    average_total_reads = average_total_reads/len(condition_count_db)
    if testImport == 'yes':
        print 'average_total_reads:',average_total_reads
    k=0
    c=math.pow(10.0,9.0)
    for key in exon_len_db:
        try:
            for condition in condition_count_db:
                total_mapped_reads = mapped_reads[condition]
                try: read_count = float(condition_count_db[condition][key])+1 ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
                except KeyError: read_count = 1 ###Was zero, but needs to be one for more realistic log2 fold calculations
                if biotype_to_examine == 'junction': region_length = 60.0
                else:
                    try: region_length = exon_len_db[key]
                    except Exception: continue ### This should only occur during testing (when restricting to one or few chromosomes)
                if read_count == 1: ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
                    rpkm = c*(float(read_count)/(float(average_total_reads)*region_length)) 
                try:
                    if region_length == 0:
                        region_length = abs(int(key[2]-key[1]))
                    rpkm = c*(read_count/(float(total_mapped_reads)*region_length))
                except Exception:
                    print condition, key
                    print 'Error Encountered... Exon or Junction of zero length encoutered... RPKM failed... Exiting AltAnalyze.'
                    print 'This error may be due to inconsistent file naming. If both exon and junction sample data is present, make sure they are named propperly.'
                    print 'For example: cancer1__exon.bed, cancer1__junction.bed (double underscore required to match these samples up)!'
                    print [read_count,total_mapped_reads,region_length];k=1; forceError
                condition_count_db[condition][key] = str(rpkm) ### Replace original counts with RPMK
        except Exception:
            if k == 1: kill
            null=[]
    return condition_count_db

def calculateGeneLevelStatistics(steady_state_export,species,expressed_gene_exon_db,normalize_feature_exp,array_names,fl,excludeLowExp=True,exportRPKMs=False):
    global UserOptions; UserOptions = fl
    exp_file = string.replace(steady_state_export,'-steady-state','')
    if normalize_feature_exp == 'RPKM':
        exp_dbase, all_exp_features, array_count = importRawCountData(exp_file,expressed_gene_exon_db,excludeLowExp=excludeLowExp)
        steady_state_db = obtainGeneCounts(expressed_gene_exon_db,species,exp_dbase,array_count,normalize_feature_exp,excludeLowExp=excludeLowExp); exp_dbase=[]
        exportGeneCounts(steady_state_export,array_names,steady_state_db)
        steady_state_db = calculateGeneRPKM(steady_state_db)
        if exportRPKMs:
            exportGeneCounts(steady_state_export,array_names,steady_state_db,dataType='RPKMs')
            
    else:
        exp_dbase, all_exp_features, array_count = importNormalizedCountData(exp_file,expressed_gene_exon_db)
        steady_state_db = obtainGeneCounts(expressed_gene_exon_db,species,exp_dbase,array_count,normalize_feature_exp); exp_dbase=[]
        exportGeneCounts(steady_state_export,array_names,steady_state_db)
    return steady_state_db, all_exp_features
    

def exportGeneCounts(steady_state_export,headers,gene_count_db,dataType='counts'):
    ### In addition to RPKM gene-level data, export gene level counts and lengths (should be able to calculate gene RPKMs from this file)
    if dataType=='counts':
        export_path = string.replace(steady_state_export,'exp.','counts.')
    else:
        export_path = steady_state_export
    export_data = export.ExportFile(export_path)
    
    title = string.join(['Ensembl']+headers,'\t')+'\n'
    export_data.write(title)
        
    for gene in gene_count_db:
        sample_counts=[]
        for count_data in gene_count_db[gene]:
            try: read_count,region_length = count_data
            except Exception: read_count = count_data
            sample_counts.append(str(read_count))
        sample_counts = string.join([gene]+sample_counts,'\t')+'\n'
        export_data.write(sample_counts)
    export_data.close()

def importGeneCounts(filename,import_type):
    ### Import non-normalized original counts and return the max value            
    counts_filename = string.replace(filename,'exp.','counts.')
    status = verifyFile(counts_filename)
    if status == 'not found': ### Occurs for non-normalized counts
        counts_filename = filename
    fn=filepath(counts_filename); x=0; count_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: array_names = t[1:]; x=1
        else:
            gene = t[0]
            if import_type == 'max':
                count_db[gene] = str(max(map(float,t[1:])))
            else:
                count_db[gene] = map(float,t[1:])
    return count_db,array_names

def calculateGeneRPKM(gene_count_db):
    """Determines the total number of reads in a sample and then calculates RPMK relative to a pre-determined junction length (60).
    60 was choosen, based on Illumina single-end read lengths of 35 (5 nt allowed overhand on either side of the junction)"""
    ### Get the total number of mapped reads (relative to all gene aligned rather than genome aligned exon reads)
    mapped_reads={}
    for gene in gene_count_db:
        index=0
        for (read_count,total_len) in gene_count_db[gene]:
            try: mapped_reads[index]+=float(read_count)
            except Exception: mapped_reads[index]=float(read_count)
            index+=1
            
    ### Use the average_total_reads when no counts reported such that 0 counts are comparable            
    average_total_reads = 0
    for i in mapped_reads: average_total_reads+=mapped_reads[i]
    average_total_reads = average_total_reads/(index+1) ###
    
    c=math.pow(10.0,9.0)
    for gene in gene_count_db:
        index=0; rpkms = []
        for (read_count,region_length) in gene_count_db[gene]:
            total_mapped_reads = mapped_reads[index]
            #print [read_count],[region_length],[total_mapped_reads]
            #if gene == 'ENSMUSG00000028186': print [read_count, index, total_mapped_reads,average_total_reads,region_length]
            if read_count == 0: read_count=1; rpkm = c*(float(read_count)/(float(average_total_reads)*region_length)) ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
            else:
                try: rpkm = c*(float(read_count+1)/(float(total_mapped_reads)*region_length)) ### read count is incremented +1 (see next line)
                except Exception: read_count=1; rpkm = c*(float(read_count)/(float(average_total_reads)*region_length)) ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
            #if gene == 'ENSMUSG00000028186': print rpkm,read_count,index,total_mapped_reads,average_total_reads,region_length
            #if gene == 'ENSMUSG00000026049': print gene_count_db[gene], mapped_reads[index], rpkm
            rpkms.append(rpkm)
            index+=1
        gene_count_db[gene] = rpkms ### Replace original counts with RPMK
    return gene_count_db

def deleteOldAnnotations(species,root_dir,dataset_name):
    db_dir = root_dir+'AltDatabase/'+species
    try:
        status = export.deleteFolder(db_dir)
        if status == 'success':
            print "...Previous experiment database deleted"
    except Exception: null=[]
    
    count_dir = root_dir+'ExpressionInput/Counts'
    try: status = export.deleteFolder(count_dir)
    except Exception: pass
    
    if 'exp.' not in dataset_name: dataset_name = 'exp.'+dataset_name
    if '.txt' not in dataset_name: dataset_name+='.txt'
    export_path = root_dir+'ExpressionInput/'+dataset_name
    try: os.remove(filepath(export_path))
    except Exception: null=[]
    try: os.remove(filepath(string.replace(export_path,'exp.','counts.')))
    except Exception: null=[]
    try: os.remove(filepath(string.replace(export_path,'exp.','novel.')))
    except Exception: null=[]

from copy_reg import pickle
from types import MethodType

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

def call_it(instance, name, args=(), kwargs=None):
    "indirect caller for instance methods and multiprocessing"
    if kwargs is None:
        kwargs = {}
    return getattr(instance, name)(*args, **kwargs)

def alignExonsAndJunctionsToEnsembl(species,exp_file_location_db,dataset_name,Multi=None):
    fl = exp_file_location_db[dataset_name]
    
    try: multiThreading = fl.multiThreading()
    except Exception: multiThreading = True
    print 'multiThreading:',multiThreading
    
    normalize_feature_exp = fl.FeatureNormalization()
    
    testImport='no'
    if 'demo_data' in fl.ExpFile():
        ### If the input files are in the AltAnalyze test directory, only analyze select chromosomes
        print 'Running AltAnalyze in TEST MODE... restricting to select chromosomes only!!!!!'
        testImport='yes' 
        
    rnaseq_begin_time = time.time()
    p = AlignExonsAndJunctionsToEnsembl(species,exp_file_location_db,dataset_name,testImport)
    chromosomes = p.getChromosomes()
    
    ### The following files need to be produced from chromosome specific sets later
    countsFile = p.countsFile()
    exonFile = p.exonFile()
    junctionFile = p.junctionFile()
    junctionCompFile = p.junctionCompFile()
    novelJunctionAnnotations = p.novelJunctionAnnotations()
    
    #chromosomes = ['chrMT']
    #p('chrY'); p('chr1'); p('chr2')
    #chromosomes = ['chr8','chr17']
    multiprocessing_pipe = True
    
    if 'exp.' not in dataset_name:
        dataset_name = 'exp.'+dataset_name
        if '.txt' not in dataset_name:
            dataset_name+='.txt'
                
    try:
        mlp=Multi
        pool_size = mlp.cpu_count()
        print 'Using %d processes' % pool_size
        if multiprocessing_pipe and multiThreading:
            ### This is like pool, but less efficient (needed to get print outs)
            s = pool_size; b=0
            chr_blocks=[]
            while s<len(chromosomes):
                chr_blocks.append(chromosomes[b:s])
                b+=pool_size; s+=pool_size
            chr_blocks.append(chromosomes[b:s])
    
            queue = mlp.Queue()
            results=[]
            #parent_conn, child_conn=multiprocessing.Pipe()
            for chromosomes in chr_blocks:
                procs=list()
                #print 'Block size:',len(chromosomes)
                for search_chr in chromosomes:
                    proc = mlp.Process(target=p, args=(queue,search_chr)) ### passing sys.stdout unfortunately doesn't work to pass the Tk string 
                    procs.append(proc)
                    proc.start()

                for _ in procs:
                    val = queue.get()
                    if p.AnalysisMode() == 'GUI': print '*',
                    results.append(val)
                    
                for proc in procs:
                    proc.join()
        
        elif multiThreading:
            pool = mlp.Pool(processes=pool_size)
            chr_vars=[]
            for search_chr in chromosomes:
                chr_vars.append(([],search_chr)) ### As an alternative for the pipe version above, pass an empty list rather than queue
            
            results = pool.map(p, chr_vars) ### worker jobs initiated in tandem
            try:pool.close(); pool.join(); pool = None
            except Exception: pass
        else:
            forceThreadingError
        print 'Read exon and junction mapping complete'
        
    except Exception,e:
        #print e
        print 'Proceeding with single-processor version align...'
        try: proc.close; proc.join; proc = None
        except Exception: pass
        try: pool.close(); pool.join(); pool = None
        except Exception: pass
        results=[] ### For single-thread compatible versions of Python
        for search_chr in chromosomes:
            result = p([],search_chr)
            results.append(result)
        
    results_organized=[]
    for result_set in results:
        if len(result_set[0])>0: ### Sometimes chromsomes are missing
            biotypes = result_set[0]
        results_organized.append(list(result_set[1:]))
    
    pooled_results = [sum(value) for value in zip(*results_organized)] # combine these counts
    pooled_results = [biotypes]+pooled_results
    p.setCountsOverview(pooled_results) # store as retreivable objects

    catFiles(countsFile,'Counts')
    catFiles(junctionFile,'junctions')
    catFiles(exonFile,'exons')
    catFiles(junctionCompFile,'comps')
    catFiles(novelJunctionAnnotations,'denovo')
    
    if normalize_feature_exp == 'RPKM':
        fastRPKMCalculate(countsFile) 
    
    rnaseq_end_time = time.time()
    print '...RNA-seq import completed in',int(rnaseq_end_time-rnaseq_begin_time),'seconds\n'
    biotypes = p.outputResults()
    return biotypes

def alignCoordinatesToGeneExternal(species,coordinates_to_annotate):
    chr_strand_gene_dbs,location_gene_db,chromosomes,gene_location_db = getChromosomeStrandCoordinates(species,'no')
    read_aligned_to_gene=0
    for (chr,strand) in coordinates_to_annotate:
        if (chr,strand) in chr_strand_gene_dbs:
            chr_gene_locations = chr_strand_gene_dbs[chr,strand]
            chr_reads = coordinates_to_annotate[chr,strand]
            chr_gene_locations.sort(); chr_reads.sort()
            ### Set GeneID for each coordinate object (primary and seconardary GeneIDs)
            read_aligned_to_gene=geneAlign(chr,chr_gene_locations,location_gene_db,chr_reads,'no',read_aligned_to_gene)
            ### Gene objects will be updated
            
def catFiles(outFileDir,folder):
    """ Concatenate all the chromosomal files but retain only the first header """
    root_dir = export.findParentDir(outFileDir)+folder+'/'
    dir_list = read_directory(root_dir)
    firstFile=True
    with open(filepath(outFileDir), 'w') as outfile:
        for fname in dir_list:
            chr_file = root_dir+fname
            header=True
            with open(filepath(chr_file)) as infile:
                for line in infile:
                    if header:
                        header=False
                        if firstFile:
                            outfile.write(line)
                            firstFile=False
                    else: outfile.write(line)
    export.deleteFolder(root_dir)
                
def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)

class AlignExonsAndJunctionsToEnsembl:
    def setCountsOverview(self, overview):
        self.biotypes_store, self.known_count, self.novel_junction_count, self.trans_splicing_reads, self.junctions_without_exon_gene_alignments, self.exons_without_gene_alignment_count = overview
    
    def getChromosomes(self):
        chr_list=list()
        for c in self.chromosomes:
            ### Sort chromosome by int number
            ci=string.replace(c,'chr','')
            try: ci = int(ci)
            except Exception: pass
            chr_list.append((ci,c))
        chr_list.sort()
        chr_list2=list()
        for (i,c) in chr_list: chr_list2.append(c) ### sorted
        return chr_list2
    
    def countsFile(self):
        return string.replace(self.expfile,'exp.','counts.')
    
    def junctionFile(self):
        junction_file = self.root_dir+'AltDatabase/'+self.species+'/RNASeq/'+self.species + '_Ensembl_junctions.txt'
        return junction_file

    def exonFile(self):
        exon_file = self.root_dir+'AltDatabase/'+self.species+'/RNASeq/'+self.species + '_Ensembl_exons.txt'
        return exon_file

    def junctionCompFile(self):
        junction_comp_file = self.root_dir+'AltDatabase/'+self.species+'/RNASeq/'+self.species + '_junction_comps_updated.txt'
        return junction_comp_file
    
    def novelJunctionAnnotations(self):
        junction_annotation_file = self.root_dir+'AltDatabase/ensembl/'+self.species+'/'+self.species + '_alternative_junctions_de-novo.txt'
        return junction_annotation_file
    
    def AnalysisMode(self): return self.analysisMode
    
    def __init__(self,species,exp_file_location_db,dataset_name,testImport):
        self.species = species; self.dataset_name = dataset_name
        self.testImport = testImport
        
        fl = exp_file_location_db[dataset_name]
        bed_dir=fl.BEDFileDir()
        root_dir=fl.RootDir()
        #self.stdout = fl.STDOUT()
        try: platformType = fl.PlatformType()
        except Exception: platformType = 'RNASeq'
        try: analysisMode = fl.AnalysisMode()
        except Exception: analysisMode = 'GUI'
        
        ### This occurs when run using the BAMtoBED pipeline in the GUI
        if 'exp.' not in dataset_name:
            dataset_name = 'exp.'+dataset_name
            if '.txt' not in dataset_name:
                dataset_name+='.txt'
            self.dataset_name = dataset_name
            
        ### Import experimentally identified junction splice-sites
        normalize_feature_exp = fl.FeatureNormalization()

        if platformType == 'RNASeq':
            chr_status = checkBEDFileFormat(bed_dir,root_dir) ### If false, need to remove 'chr' from the search_chr
        else:
            chr_status = True

        #self.fl = fl # Can not pass this object in pool or it breaks
        self.platformType = platformType
        self.analysisMode = analysisMode
        self.root_dir = root_dir
        self.normalize_feature_exp = normalize_feature_exp
        self.bed_dir = bed_dir
        self.chr_status = chr_status
        self.exonBedBuildStatus = fl.ExonBedBuildStatus()
        self.expfile = root_dir+'ExpressionInput/'+dataset_name
        
        if testImport == 'yes':
            print 'Chromosome annotation detected =',chr_status
        #if self.exonBedBuildStatus == 'yes':
        reformatExonFile(species,'exon',chr_status) ### exports BED format exons for exon expression extraction
        """
        Strategies to reduce memory in RNASeq:
        1)	(done)Delete old AltDatabase-local version if it exists before starting
        2)	(done)Check to see if a file exists before writing it and if so append rather than create
        3)	(done)Get counts last and normalize last in for exons and junctions separately.
        4)	(done)Delete objects explicitly before importing any new data (define a new function that just does this).
        5)	(done)Get all chromosomes first then parse exon and junction coordinate data on a per known chromosome basis.
        6)	(done)Prior to deleting all junction/exon object info for each chromsome, save the coordinate(key)-to-annotation information for the read count export file."""
        
        ### Delete any existing annotation databases that currently exist (redundant with below)
        deleteOldAnnotations(species,root_dir,dataset_name)
        
        ###Define variables to report once reads for all chromosomes have been aligned
        #global self.known_count; global self.novel_junction_count; global self.one_found; global self.not_found; global self.both_found; global self.trans_splicing_reads
        #global self.junctions_without_exon_gene_alignments; global self.exons_without_gene_alignment_count; global self.junction_simple_db; global self.chr_strand_gene_dbs

        self.known_count=0; self.novel_junction_count=0; self.one_found=0; self.not_found=0; self.both_found=0; self.trans_splicing_reads=0
        self.junctions_without_exon_gene_alignments=0; self.exons_without_gene_alignment_count=0; self.junction_simple_db={}
        
        ###Begin Chromosome specific read to exon alignments
        self.chr_strand_gene_dbs,self.location_gene_db,chromosomes,self.gene_location_db = getChromosomeStrandCoordinates(species,testImport)
        self.chromosomes = chromosomes

        print "Processing exon/junction coordinates sequentially by chromosome"
        print "Note: this step is time intensive (can be hours) and no print statements may post for a while"
        
    def outputResults(self):
        exportDatasetLinkedGenes(self.species,self.gene_location_db,self.root_dir) ### Include an entry for gene IDs to include constitutive expression for RPKM normalized data     
        chr_gene_locations=[]; self.location_gene_db=[]; self.chr_strand_gene_dbs=[] 
        
        #print 'user coordinates imported/processed'
        #print 'Importing read counts from coordinate data...'

        biotypes = self.biotypes_store
        ### Output summary statistics
        if self.normalize_feature_exp != 'none':
            print self.normalize_feature_exp, 'normalization complete'
        if 'junction' in biotypes:
            print 'Imported Junction Statistics:'
            print '    ',self.known_count, 'junctions found in Ensembl/UCSC and',self.novel_junction_count,'are novel'
            print '    ',self.trans_splicing_reads,'trans-splicing junctions found (two aligning Ensembl genes)'
            print '    ',self.junctions_without_exon_gene_alignments, 'junctions where neither splice-site aligned to a gene'
            if (float(self.known_count)*10)<float(self.novel_junction_count):
                print '\nWARNING!!!!! Few junctions aligned to known exons. Ensure that the AltAnalyze Ensembl database\nversion matches the genome build aligned to!\n'
        if 'exon' in biotypes:
            print 'Imported Exon Statistics:'
            print '    ',self.exons_without_gene_alignment_count, 'exons where neither aligned to a gene'
        print 'User databases and read counts written to:', self.root_dir[:-1]+'ExpressionInput'
        
        
        ### END CHROMOSOME SPECIFIC ANALYSES
        if self.exonBedBuildStatus == 'yes':
            bedfile = moveBAMtoBEDFile(self.species,self.dataset_name,self.root_dir)
            print 'Exon BED file updated with novel exon predictions from junction file'
            return bedfile; sys.exit()
                
        clearObjectsFromMemory(self.junction_simple_db); self.junction_simple_db=[]
        return biotypes

    def test(self, search_chr):
        print search_chr
        
    def __call__(self, queue, search_chr):
        try:
            #sys.stdout = self.stdout
            platformType = self.platformType
            testImport = self.testImport
            species = self.species
            dataset_name = self.dataset_name
            
            platformType = self.platformType
            analysisMode = self.analysisMode
            root_dir = self.root_dir
            normalize_feature_exp = self.normalize_feature_exp
            bed_dir = self.bed_dir
            chr_status = self.chr_status
            
            junction_annotations={}
            if chr_status == False:
                searchchr = string.replace(search_chr,'chr','')
            else:
                searchchr = search_chr
            if platformType == 'RNASeq':
                junction_db,biotypes,algorithms = importBEDFile(bed_dir,root_dir,species,normalize_feature_exp,searchChr=searchchr,testImport=testImport)
            else:
                normalize_feature_exp = 'quantile'
                junction_db,biotypes,algorithms = importExpressionMatrix(bed_dir,root_dir,species,fl,'no',search_chr=searchchr)
            
            self.biotypes_store = biotypes
            if len(junction_db)>0:
                    
                ### Determine which kind of data is being imported, junctions, exons or both
                unmapped_exon_db={}
                    
                if 'junction' in biotypes:
                    ### Get all known junction splice-sites
                    ens_junction_coord_db = importExonAnnotations(species,'junction_coordinates',search_chr)
                    if testImport == 'yes':
                        print len(ens_junction_coord_db),'Ensembl/UCSC junctions imported'
                ### Identify known junctions sites found in the experimental dataset (perfect match)
                novel_junction_db={}; novel_exon_db={}
                for key in junction_db:
                    ji=junction_db[key]
                    if ji.BioType()=='junction':
                        if key in ens_junction_coord_db:
                            jd=ens_junction_coord_db[key]
                            ji.setExonAnnotations(jd)
                            self.known_count+=1
                        else:
                            novel_junction_db[key]=junction_db[key]; self.novel_junction_count+=1
                            #if 75953254 in key: print key; sys.exit()
                    else:
                        unmapped_exon_db[key]=junction_db[key]
                ens_exon_db = importExonAnnotations(species,'exon',search_chr)
                
                if 'junction' in biotypes:
                    if testImport == 'yes':
                        print self.known_count,  'junctions found in Ensembl/UCSC and',len(novel_junction_db),'are novel.'
                    
                    ### Separate each junction into a 5' and 3' splice site (exon1_coord_db and exon2_coord_db)
                    exon1_coord_db={}; exon2_coord_db={}    
                    for (chr,exon1_stop,exon2_start) in ens_junction_coord_db:
                        jd = ens_junction_coord_db[(chr,exon1_stop,exon2_start)]
                        exon1_coord_db[chr,exon1_stop] = jd.GeneID(),string.split(jd.ExonRegionIDs(),'-')[0]
                        exon2_coord_db[chr,exon2_start] = jd.GeneID(),string.split(jd.ExonRegionIDs(),'-')[1]
                    clearObjectsFromMemory(ens_junction_coord_db); ens_junction_coord_db=[] ### Clear object from memory
                    
                    ### Get and re-format individual exon info
                    exon_region_db={}
                    #if 'exon' not in biotypes:
                    for gene in ens_exon_db:
                        for rd in ens_exon_db[gene]:
                            exon_region_db[gene,rd.ExonRegionIDs()]=rd
        
                    ### Add the exon annotations from the known junctions to the exons to export dictionary
                    exons_to_export={}
                    for key in junction_db:
                        ji=junction_db[key]
                        if ji.ExonAnnotations() != None:
                            jd = ji.ExonAnnotations()
                            exon1, exon2 = string.split(jd.ExonRegionIDs(),'-')
                            key1 = jd.GeneID(),exon1; key2 = jd.GeneID(),exon2
                            exons_to_export[key1] = exon_region_db[key1]
                            exons_to_export[key2] = exon_region_db[key2]
                            
                    ### For novel experimental junctions, identify those with at least one matching known 5' or 3' site
                    exons_not_identified = {}; novel_exon_coordinates={}
                    for (chr,exon1_stop,exon2_start) in novel_junction_db:
                        ji = novel_junction_db[(chr,exon1_stop,exon2_start)]
                        coord = [exon1_stop,exon2_start]; coord.sort()
                        if (chr,exon1_stop) in exon1_coord_db and (chr,exon2_start) in exon2_coord_db:
                            ### Assign exon annotations to junctions where both splice-sites are known in Ensembl/UCSC
                            ### Store the exon objects, genes and regions (le is a tuple of gene and exon region ID)
                            ### Do this later for the below un-assigned exons
                            le=exon1_coord_db[(chr,exon1_stop)]; ji.setLeftExonAnnotations(le); ji.setLeftExonRegionData(exon_region_db[le])
                            re=exon2_coord_db[(chr,exon2_start)]; ji.setRightExonAnnotations(re);  ji.setRightExonRegionData(exon_region_db[re])                                         
                            if le[0] != re[0]: ### Indicates Trans-splicing (e.g., chr7:52,677,568-52,711,750 mouse mm9)
                                ji.setTransSplicing(); #print exon1_stop,le,exon2_start,re,ji.Chr(),ji.Strand()
                            self.both_found+=1; #print 'five',(chr,exon1_stop,exon2_start),exon1_coord_db[(chr,exon1_stop)]
                        else:
                            if (chr,exon1_stop) in exon1_coord_db: ### hence, exon1_stop is known, so report the coordinates of exon2 as novel
                                le=exon1_coord_db[(chr,exon1_stop)]; ji.setLeftExonAnnotations(le)
                                self.one_found+=1; #print 'three',(chr,exon1_stop,exon2_start),exon1_coord_db[(chr,exon1_stop)]
                                novel_exon_coordinates[ji.Chr(),exon2_start] = ji,'left',ji.Exon2Stop() ### Employ this strategy to avoid duplicate exons with differing lengths (mainly an issue if analyzing only exons results)
                                ji.setNovel('side')
                            elif (chr,exon2_start) in exon2_coord_db: ### hence, exon2_start is known, so report the coordinates of exon1 as novel
                                re=exon2_coord_db[(chr,exon2_start)]; ji.setRightExonAnnotations(re) ### In very rare cases, a gene can be assigned here, even though the splice-site is on the opposite strand (not worthwhile filtering out)
                                self.one_found+=1; #print 'three',(chr,exon1_stop,exon2_start),exon1_coord_db[(chr,exon1_stop)]
                                novel_exon_coordinates[ji.Chr(),exon1_stop] = ji,'right',ji.Exon1Start()
                                ji.setNovel('side')
                            else:
                                self.not_found+=1; #if self.not_found < 10: print (chr,exon1_stop,exon2_start)
                                novel_exon_coordinates[ji.Chr(),exon1_stop] = ji,'right',ji.Exon1Start()
                                novel_exon_coordinates[ji.Chr(),exon2_start] = ji,'left',ji.Exon2Stop()
                                ji.setNovel('both')
                            ### We examine reads where one splice-site aligns to a known but the other not, to determine if trans-splicing occurs
                            try: exons_not_identified[chr,ji.Strand()].append((coord,ji))
                            except KeyError: exons_not_identified[chr,ji.Strand()] = [(coord,ji)]
                    """
                    if fl.ExonBedBuildStatus() == 'no':
                        exportNovelJunctions(species,novel_junction_db,condition_count_db,root_dir,dataset_name,'junction') ### Includes known exons
                    """
                    
                    #print self.both_found, ' where both and', self.one_found, 'where one splice-site are known out of',self.both_found+self.one_found+self.not_found
                    #print 'Novel junctions where both splice-sites are known:',self.both_found
                    #print 'Novel junctions where one splice-site is known:',self.one_found
                    #print 'Novel junctions where the splice-sites are not known:',self.not_found
                    clearObjectsFromMemory(exon_region_db); exon_region_db=[] ### Clear memory of this object
        
                    read_aligned_to_gene=0
                    for (chr,strand) in exons_not_identified:
                        if (chr,strand) in self.chr_strand_gene_dbs:
                            chr_gene_locations = self.chr_strand_gene_dbs[chr,strand]
                            chr_reads = exons_not_identified[chr,strand]
                            chr_gene_locations.sort(); chr_reads.sort()
                            ### Set GeneID for each coordinate object (primary and seconardary GeneIDs)
                            read_aligned_to_gene=geneAlign(chr,chr_gene_locations,self.location_gene_db,chr_reads,'no',read_aligned_to_gene)
                    #print read_aligned_to_gene, 'novel junctions aligned to Ensembl genes out of',self.one_found+self.not_found
                    clearObjectsFromMemory(exons_not_identified); exons_not_identified=[] ## Clear memory of this object
                    
                    for key in novel_junction_db:
                        (chr,exon1_stop,exon2_start) = key
                        ji=novel_junction_db[key]
                        if ji.GeneID() == None:
                            try:
                                if ji.SecondaryGeneID() != None:
                                    ### Occurs if mapping is to the 5'UTR of a gene for the left splice-site (novel alternative promoter)
                                    ji.setGeneID(ji.SecondaryGeneID()); ji.setSecondaryGeneID(''); #print key, ji.GeneID(), ji.Strand(), ji.SecondaryGeneID()
                            except Exception: null=[]
                        if ji.GeneID() != None:
                            geneid = ji.GeneID()
                            proceed = 'no'                
                            if ji.SpliceSitesFound() == None: proceed = 'yes'; coordinates = [exon1_stop,exon2_start]
                            elif ji.SpliceSitesFound() == 'left': proceed = 'yes'; coordinates = [exon1_stop,exon2_start]
                            elif ji.SpliceSitesFound() == 'right': proceed = 'yes'; coordinates = [exon1_stop,exon2_start]
                            if proceed == 'yes':
                                for coordinate in coordinates:
                                    if ji.TransSplicing() == 'yes':
                                        #print ji.Chr(),ji.GeneID(), ji.SecondaryGeneID(), ji.Exon1Stop(), ji.Exon2Start()
                                        self.trans_splicing_reads+=1
                                        if ji.checkExonPosition(coordinate) == 'right': geneid = ji.SecondaryGeneID()
                                    if abs(exon2_start-exon1_stop)==1: eventType = 'novel-exon-intron' ### Indicates intron-exon boundary (intron retention)
                                    else: eventType = 'novel'
                                    exon_data = (coordinate,ji.Chr()+'-'+str(coordinate),eventType)
                                    try: novel_exon_db[geneid].append(exon_data)
                                    except KeyError: novel_exon_db[geneid] = [exon_data]
                                    
                        else:
                            ### write these out
                            self.junctions_without_exon_gene_alignments+=1
            
                    ### Remove redundant exon entries and store objects    
                    for key in novel_exon_db:
                        exon_data_objects=[]
                        exon_data_list = unique.unique(novel_exon_db[key])
                        exon_data_list.sort()
                        for e in exon_data_list:
                            ed = ExonInfo(e[0],e[1],e[2])
                            exon_data_objects.append(ed)
                        novel_exon_db[key] = exon_data_objects
                        
                    #print self.trans_splicing_reads,'trans-splicing junctions found (two aligning Ensembl genes).'
                    #print self.junctions_without_exon_gene_alignments, 'junctions where neither splice-site aligned to a gene'
                    #if 'X' in search_chr: print len(ens_exon_db),len(ens_exon_db['ENSMUSG00000044424'])
                    alignReadsToExons(novel_exon_db,ens_exon_db,testImport=testImport)
                    
                    ### Link exon annotations up with novel junctions
                    junction_region_db,exons_to_export = annotateNovelJunctions(novel_junction_db,novel_exon_db,exons_to_export)
                    
                    ### Add the exon region data from known Ensembl/UCSC matched junctions to junction_region_db for recipricol junction analysis
                    for key in junction_db:
                        ji=junction_db[key]; jd = ji.ExonAnnotations()
                        try:
                            uid = jd.GeneID()+':'+jd.ExonRegionIDs();  ji.setUniqueID(uid)
                            try: junction_region_db[jd.GeneID()].append((formatID(uid),jd.ExonRegionIDs()))
                            except KeyError: junction_region_db[jd.GeneID()] = [(formatID(uid),jd.ExonRegionIDs())]
                        except AttributeError: null=[] ### Occurs since not all entries in the dictionary are perfect junction matches
                    
                    try: novel_exon_coordinates = collapseNoveExonBoundaries(novel_exon_coordinates,root_dir+dataset_name) ### Joins inferred novel exon-IDs (5' and 3' splice sites) from adjacent and close junction predictions
                    except Exception: pass ### No errors encountered before
                    #if self.exonBedBuildStatus == 'yes':
                     ### Append to the exported BED format exon coordinate file
                    bedfile = exportNovelExonToBedCoordinates(species,novel_exon_coordinates,chr_status,searchChr=searchchr)
                    
                    ### Identify reciprocol junctions and retrieve splice-event annotations for exons and inclusion junctions
                    
                    junction_annotations,critical_exon_annotations = JunctionArray.inferJunctionComps(species,('RNASeq',junction_region_db,root_dir),searchChr=searchchr)
                    clearObjectsFromMemory(junction_region_db); junction_region_db=[]   
                    
                    ### Reformat these dictionaries to combine annotations from multiple reciprocol junctions
                    junction_annotations = combineExonAnnotations(junction_annotations)
                    critical_exon_annotations = combineExonAnnotations(critical_exon_annotations)

                if 'exon' in biotypes:
                    if testImport == 'yes':
                        print len(unmapped_exon_db),'exon genomic locations imported.'
                    ### Create a new dictionary keyed by chromosome and strand
                    exons_not_aligned={}
                    for (chr,exon1_stop,exon2_start) in unmapped_exon_db:
                        ji = unmapped_exon_db[(chr,exon1_stop,exon2_start)]
                        coord = [exon1_stop,exon2_start]; coord.sort()
                        try: exons_not_aligned[chr,ji.Strand()].append((coord,ji))
                        except KeyError: exons_not_aligned[chr,ji.Strand()] = [(coord,ji)]
                                                          
                    read_aligned_to_gene=0
                    for (chr,strand) in exons_not_aligned:
                        if (chr,strand) in self.chr_strand_gene_dbs:
                            chr_gene_locations = self.chr_strand_gene_dbs[chr,strand]
                            chr_reads = exons_not_aligned[chr,strand]
                            chr_gene_locations.sort(); chr_reads.sort()
                            read_aligned_to_gene=geneAlign(chr,chr_gene_locations,self.location_gene_db,chr_reads,'no',read_aligned_to_gene)
                            
                    #print read_aligned_to_gene, 'exons aligned to Ensembl genes out of',self.one_found+self.not_found
                    
                    align_exon_db={}; exons_without_gene_alignments={}; multigene_exon=0
                    for key in unmapped_exon_db:
                        (chr,exon1_stop,exon2_start) = key
                        ji=unmapped_exon_db[key]
                        if ji.GeneID() == None:
                            try:
                                if ji.SecondaryGeneID() != None:
                                    ### Occurs if mapping outside known exon boundaries for one side of the exon
                                    ji.setGeneID(ji.SecondaryGeneID()); ji.setSecondaryGeneID(''); #print key, ji.GeneID(), ji.Strand(), ji.SecondaryGeneID()
                            except Exception: null=[]
                        else:
                            if 'ENS' in ji.JunctionID():
                                if ji.GeneID() not in ji.JunctionID(): ### Hence, there were probably two overlapping Ensembl genes and the wrong was assigned based on the initial annotations
                                    original_geneid = string.split(ji.JunctionID(),':')[0]
                                    if original_geneid in ens_exon_db: ji.setGeneID(original_geneid) #check if in ens_exon_db (since chromosome specific)
                                    
                        if ji.GeneID() != None:
                            geneid = ji.GeneID()
                            coordinates = [exon1_stop,exon2_start]
                            for coordinate in coordinates:
                                if ji.TransSplicing() != 'yes': ### This shouldn't occur for exons
                                    exon_data = (coordinate,ji.Chr()+'-'+str(coordinate),'novel')
                                    try: align_exon_db[geneid].append(exon_data)
                                    except KeyError: align_exon_db[geneid] = [exon_data]
                                else:
                                    multigene_exon+=1 ### Shouldn't occur due to a fix in the gene-alignment method which will find the correct gene on the 2nd interation
                        else: exons_without_gene_alignments[key]=ji; self.exons_without_gene_alignment_count+=1
        
                    ### Remove redundant exon entries and store objects (this step may be unnecessary)
                    for key in align_exon_db:
                        exon_data_objects=[]
                        exon_data_list = unique.unique(align_exon_db[key])
                        exon_data_list.sort()
                        for e in exon_data_list:
                            ed = ExonInfo(e[0],e[1],e[2])
                            exon_data_objects.append(ed)
                        align_exon_db[key] = exon_data_objects
                        
                    #print self.exons_without_gene_alignment_count, 'exons where neither aligned to a gene'
                    #if self.exons_without_gene_alignment_count>3000: print 'NOTE: Poor mapping of these exons may be due to an older build of\nEnsembl than the current version. Update BAMtoBED mappings to correct.'
            
                    begin_time = time.time()
                    alignReadsToExons(align_exon_db,ens_exon_db)
                    
                    end_time = time.time()
                    if testImport == 'yes':
                        print 'Exon sequences aligned to exon regions in',int(end_time-begin_time),'seconds'
        
                    ### Combine the start and end region alignments into a single exon annotation entry
                    combineDetectedExons(unmapped_exon_db,align_exon_db,novel_exon_db)
                    clearObjectsFromMemory(unmapped_exon_db); clearObjectsFromMemory(align_exon_db); clearObjectsFromMemory(novel_exon_db)
                    unmapped_exon_db=[]; align_exon_db=[]; novel_exon_db=[]
                    """            
                    if fl.ExonBedBuildStatus() == 'no':
                        exportNovelJunctions(species,exons_without_gene_alignments,condition_count_db,root_dir,dataset_name,'exon') ### Includes known exons
                    """
                    clearObjectsFromMemory(exons_without_gene_alignments); exons_without_gene_alignments=[]
                    
                ### Export both exon and junction annotations
                
                if 'junction' in biotypes:
                    ### Export the novel user exon annotations    
                    exportDatasetLinkedExons(species,exons_to_export,critical_exon_annotations,root_dir,testImport=testImport,searchChr=searchchr)            
                            
                ### Export the novel user exon-junction annotations (original junction_db objects updated by above processing)
                exportDatasetLinkedJunctions(species,junction_db,junction_annotations,root_dir,testImport=testImport,searchChr=searchchr)
        
                ### Clear memory once results are exported (don't want to delete actively used objects)
                if 'junction' in biotypes:
                    clearObjectsFromMemory(exons_to_export); clearObjectsFromMemory(critical_exon_annotations)
                    clearObjectsFromMemory(novel_junction_db); novel_junction_db=[]
                    clearObjectsFromMemory(novel_exon_coordinates); novel_exon_coordinates=[]    
                    exons_to_export=[]; critical_exon_annotations=[]
                    clearObjectsFromMemory(exon1_coord_db); clearObjectsFromMemory(exon2_coord_db)
                    exon1_coord_db=[]; exon2_coord_db=[]
                if 'exon' in biotypes:            
                    clearObjectsFromMemory(exons_not_aligned); exons_not_aligned=[]
                    clearObjectsFromMemory(ens_exon_db); ens_exon_db=[]
                                
                ### Add chromsome specific junction_db data to a simple whole genome dictionary
                for key in junction_db:
                    ji = junction_db[key]
                    if ji.GeneID()!=None and ji.UniqueID()!=None: self.junction_simple_db[key]=ji.UniqueID()
        
                #returnLargeGlobalVars()
                clearObjectsFromMemory(junction_db); clearObjectsFromMemory(junction_annotations)
                junction_db=[]; junction_annotations=[]; chr_reads=[]
    
                for biotype in biotypes:
                    ### Import Read Counts (do this last to conserve memory)
            
                    if platformType == 'RNASeq':
    
                         condition_count_db,exon_len_db,biotypes2,algorithms = importBEDFile(bed_dir,root_dir,species,normalize_feature_exp,getReads=True,searchChr=searchchr,getBiotype=biotype,testImport=testImport,filteredJunctions=self.junction_simple_db)
                    else:
                         condition_count_db,exon_len_db,biotypes2,algorithms = importExpressionMatrix(bed_dir,root_dir,species,fl,'yes',getBiotype=biotype)
                    ###First export original counts, rather than quantile normalized or RPKM
                    self.exportJunctionCounts(species,self.junction_simple_db,exon_len_db,condition_count_db,root_dir,dataset_name,biotype,'counts',searchChr=searchchr)
                    clearObjectsFromMemory(condition_count_db); clearObjectsFromMemory(exon_len_db); condition_count_db=[]; exon_len_db=[]
            if analysisMode == 'commandline':
                    print 'finished parsing data for chromosome:',search_chr ### Unix platforms are not displaying the progress in real-time
            else:
                pass #print "*",
            try: queue.put([self.biotypes_store, self.known_count, self.novel_junction_count, self.trans_splicing_reads, self.junctions_without_exon_gene_alignments, self.exons_without_gene_alignment_count])
            except Exception:
                ### If queue is not a multiprocessing object
                queue = [self.biotypes_store, self.known_count, self.novel_junction_count, self.trans_splicing_reads, self.junctions_without_exon_gene_alignments, self.exons_without_gene_alignment_count]
            return queue
        
        except Exception:
            print traceback.format_exc()
            error(traceback.format_exc())
            multiprocessing.log_to_stderr().setLevel(logging.DEBUG)
            raise
    
    def exportJunctionCounts(self,species,junction_simple_db,exon_len_db,condition_count_db,root_dir,dataset_name,biotype,count_type,searchChr=None):
        if 'exp.' not in dataset_name: dataset_name = 'exp.'+dataset_name
        if '.txt' not in dataset_name: dataset_name+='.txt'
        export_path = root_dir+'ExpressionInput/'+dataset_name
        if count_type == 'counts':
            export_path = string.replace(export_path,'exp.','counts.') ### separately export counts
            if searchChr !=None:
                export_path = string.replace(export_path,'ExpressionInput','ExpressionInput/Counts')
                export_path = string.replace(export_path,'.txt','.'+searchChr+'.txt')

        self.countsFile = export_path
        
        if self.testImport == 'yes':
            print 'Writing',export_path
        export_data,status = AppendOrWrite(export_path)
        
        if status == 'not found':
            title = ['AltAnalyze_ID']
            for condition in condition_count_db: title.append(condition)
            export_data.write(string.join(title,'\t')+'\n')
        
        for key in self.junction_simple_db:
            chr,exon1_stop,exon2_start = key
            if biotype == 'junction':
                coordinates = chr+':'+str(exon1_stop)+'-'+str(exon2_start)
            elif biotype == 'exon':
                coordinates = chr+':'+str(exon1_stop-1)+'-'+str(exon2_start+1)
            try:
                null=exon_len_db[key]
                if count_type == 'counts': values = [self.junction_simple_db[key]+'='+coordinates]
                else: values = [self.junction_simple_db[key]]
                for condition in condition_count_db: ###Memory crash here
                    count_db = condition_count_db[condition]
                    try: read_count = count_db[key]
                    except KeyError: read_count = '0'
                    values.append(read_count)
                export_data.write(string.join(values,'\t')+'\n')
            except Exception: null=[]
        export_data.close()
    
    def countsDir(self):
        return self.countsFile
    
def calculateRPKMsFromGeneCounts(filename,species,AdjustExpression):
    """ Manual way of calculating gene RPKMs from gene counts only """
    gene_lengths = getGeneExonLengths(species)
    fastRPKMCalculate(filename,GeneLengths=gene_lengths,AdjustExpression=AdjustExpression)
    
def fastRPKMCalculate(counts_file,GeneLengths=None,AdjustExpression=True):
    export_path = string.replace(counts_file,'counts.','exp.')
    export_data = export.ExportFile(export_path) ### Write this new file
    fn=filepath(counts_file); header=True
    exon_sum_array=[]; junction_sum_array=[]
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if header:
            samples = t[1:]
            header=False
            exon_sum_array=[0]*len(samples)
            junction_sum_array=[0]*len(samples)
        else:
            try: values = map(float,t[1:])
            except Exception:
                print traceback.format_exc()
                print t
                badCountsLine
            ### get the total reads/sample
            if '-' in string.split(t[0],'=')[0]:
                junction_sum_array = [sum(value) for value in zip(*[junction_sum_array,values])]
            else:
                exon_sum_array = [sum(value) for value in zip(*[exon_sum_array,values])]
                
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides warnings associated with Scipy for n=1 sample comparisons
        jatr=Average(junction_sum_array) # Average of the total maped reads
        eatr=Average(exon_sum_array) # Average of the total maped reads
    
    if AdjustExpression:
        offset = 1
    else:
        offset = 0
    header=True
    c=math.pow(10.0,9.0)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if header:
            export_data.write(line) ### Write header
            header=False
        else:
            try:
                exon_id,coordinates = string.split(t[0],'=')
                coordinates = string.split(coordinates,':')[1]
                coordinates = string.split(coordinates,'-')
                l=abs(int(coordinates[1])-int(coordinates[0])) ### read-length  
            except Exception: ### Manual way of calculating gene RPKMs from gene counts only
                exon_id = t[0]
                try: l = GeneLengths[exon_id]
                except Exception: continue #Occurs when Ensembl genes supplied from an external analysis 
                   
            try: read_counts = map(lambda x: int(x)+offset, t[1:])
            except Exception: read_counts = map(lambda x: int(float(x))+offset, t[1:])
            if '-' in exon_id:
                count_stats = zip(read_counts,junction_sum_array)
                atr = jatr
                l=60
            else:
                count_stats = zip(read_counts,exon_sum_array)
                atr = eatr
            values=[]
            
            #rpkm = map(lambda (r,t): c*(r/(t*l)), count_stats) ### Efficent way to convert to rpkm, but doesn't work for 0 counts
            for (r,t) in count_stats:
                if r == 1: ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
                    t = atr
                try:
                    rpkm = str(c*(r/(t*l)))
                    #print c,r,t,l,exon_id,rpkm;sys.exit()
                    values.append(rpkm)
                except Exception,e:
                    print e
                    print t[0]
                    print 'Error Encountered... Exon or Junction of zero length encoutered... RPKM failed... Exiting AltAnalyze.'
                    print 'This error may be due to inconsistent file naming. If both exon and junction sample data is present, make sure they are named propperly.'
                    print 'For example: cancer1__exon.bed, cancer1__junction.bed (double underscore required to match these samples up)!'
                    print [r,t,l];k=1; forceError

            values = string.join([exon_id]+values,'\t')+'\n'
            export_data.write(values)
    export_data.close()
            
def mergeCountFiles(counts_file1,counts_file2):
    ### Used internally to merge count files that are very large and too time-consuming to recreate (regenerate them)
    export_path = string.replace(counts_file2,'counts.','temp-counts.')
    export_data = export.ExportFile(export_path) ### Write this new file
    fn=filepath(counts_file1); header=True
    count_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if header:
            samples = t[1:]
            header=False
            si = samples.index('H9.102.2.5.bed')+1
            
        else:
            try: value = t[si]
            except Exception: print t; sys.exit()
            ### get the total reads/sample
            count_db[t[0]] = value

    fn=filepath(counts_file2); header=True
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if header:
            samples = t[1:]
            header=False
            si = samples.index('H9.102.2.5.bed')+1
            export_data.write(line)
        else:
            try: t[si] = count_db[t[0]]
            except Exception: pass ### keep the current value
            export_data.write(string.join(t,'\t')+'\n')
    export_data.close()
    
def getGeneExonLengths(species):
    gene_lengths={}
    filename = 'AltDatabase/'+species+'/RNASeq/'+species+'_Ensembl_exons.txt'
    fn=filepath(filename)
    firstLine=True
    for line in open(fn,'rU').xreadlines():             
        line = line.rstrip('\n')
        if firstLine:
            firstLine=False
        else:
            t = string.split(line,'\t')
            geneID = t[2]; start = int(t[6]); end = int(t[7]); exonID = t[1]
            if 'E' in exonID:
                try: gene_lengths[geneID]+=abs(end-start)
                except Exception: gene_lengths[geneID]=abs(end-start)
    return gene_lengths

def importRawCountData(filename,expressed_gene_exon_db,excludeLowExp=True):
    """ Identifies exons or junctions to evaluate gene-level expression. This function, as it is currently written:
    1) examines the RPKM and original read counts associated with all exons
    2) removes exons/junctions that do not meet their respective RPKM AND read count cutoffs
    3) returns ONLY those exons and genes deemed expressed, whether constitutive selected or all exons
    """
    
    ### Get expression values for exon/junctions to analyze
    seq_ids_to_import={}
    for gene in expressed_gene_exon_db:
        for exonid in expressed_gene_exon_db[gene]: seq_ids_to_import[exonid]=[]
    
    ### Define thresholds
    exon_exp_threshold = UserOptions.ExonExpThreshold()
    junction_exp_threshold = UserOptions.JunctionExpThreshold()
    exon_rpkm_threshold = UserOptions.ExonRPKMThreshold()

    gene_rpkm_threshold = UserOptions.RPKMThreshold()
    gene_exp_threshold = UserOptions.GeneExpThreshold()
    
    ### Import RPKM normalized expression values               
    fn=filepath(filename); x=0; rpkm_dbase={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: array_names = t[1:]; x=1
        else:
            exon_id=t[0]
            max_count=max(map(float,t[1:]))
            if max_count>=exon_rpkm_threshold or excludeLowExp==False: rpkm_dbase[exon_id]=[] ### Only retain exons/junctions meeting the RPKM threshold
            
    ### Import non-normalized original counts                
    counts_filename = string.replace(filename,'exp.','counts.')
    fn=filepath(counts_filename); x=0; exp_dbase={}
    all_exp_features={} ### Don't filter for only gene-expression reporting
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: array_names = t[1:]; x=1
        else:
            exon_id,coordinates = string.split(t[0],'=')
            coordinates = string.split(coordinates,':')[1]
            coordinates = string.split(coordinates,'-')
            length=abs(int(coordinates[1])-int(coordinates[0]))
            max_count=max(map(float,t[1:])); proceed = 'no'
            if '-' in exon_id:
                length = 60.0
                if max_count>=junction_exp_threshold or excludeLowExp==False:
                    ### Only considered when exon data is not present in the analysis
                    proceed = 'yes'
            elif max_count>=exon_exp_threshold or excludeLowExp==False: proceed = 'yes'
            if proceed == 'yes' and exon_id in rpkm_dbase: ### Ensures that the maximum sample (not group) user defined count threshold is achieved at the exon or junction-level
                all_exp_features[exon_id]=None
                if exon_id in seq_ids_to_import:### Forces an error if not in the steady-state pre-determined set (CS or all-exons) - INCLUDE HERE TO FILTER ALL FEATURES
                    exp_dbase[exon_id] = t[1:],length ### Include sequence length for normalization

    for exon in exp_dbase: array_count = len(exp_dbase[exon][0]); break
    try:null=array_count
    except Exception:
        print 'No exons or junctions considered expressed (based user thresholds). Exiting analysis.'; force_exit
    return exp_dbase, all_exp_features, array_count

def importNormalizedCountData(filename,expressed_gene_exon_db):
    ### Get expression values for exon/junctions to analyze
    seq_ids_to_import={}
    for gene in expressed_gene_exon_db:
        for exonid in expressed_gene_exon_db[gene]: seq_ids_to_import[exonid]=[]
            
    ### Define thresholds
    exon_exp_threshold = UserOptions.ExonExpThreshold()
    junction_exp_threshold = UserOptions.JunctionExpThreshold()
    exon_rpkm_threshold = UserOptions.ExonRPKMThreshold()
    
    gene_rpkm_threshold = UserOptions.RPKMThreshold()
    gene_exp_threshold = UserOptions.GeneExpThreshold()
    
    ### Import non-normalized original counts                
    fn=filepath(filename); x=0; exp_dbase={}
    all_exp_features={} ### Don't filter for only gene-expression reporting
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: array_names = t[1:]; x=1
        else:
            exon_id=t[0]; proceed = 'no'
            max_count=max(map(float,t[1:]))
            if '-' in exon_id:
                if max_count>=junction_exp_threshold: proceed = 'yes'
            elif max_count>=exon_exp_threshold: proceed = 'yes'
            if proceed == 'yes': ### Ensures that the maximum sample (not group) user defined count threshold is achieved at the exon or junction-level
                all_exp_features[exon_id]=None
                if exon_id in seq_ids_to_import: ### If a "constitutive" or exon-level feature (filter missing prior to 2.0.8 - bug)
                    exp_dbase[exon_id] = t[1:],0  ### Add the zero just to comply with the raw count input format (indicates exon length)

    for exon in exp_dbase: array_count = len(exp_dbase[exon][0]); break        
    return exp_dbase, all_exp_features, array_count

def obtainGeneCounts(expressed_gene_exon_db,species,exp_dbase,array_count,normalize_feature_exp,excludeLowExp=True):
    ###Calculate avg expression for each sample for each exon (using constitutive or all exon values)

    if excludeLowExp == False:
        gene_lengths = getGeneExonLengths(species)
         
    steady_state_db={}
    for gene in expressed_gene_exon_db:
        x = 0; gene_sum=0
        exon_list = expressed_gene_exon_db[gene]
        while x < array_count:
            exp_list=[]; len_list=[]
            for exon in exon_list:
                try:
                    exp_val = exp_dbase[exon][0][x]
                    if normalize_feature_exp == 'RPKM':
                        ### Decided to include all exons, expressed or not to prevent including lowly expressed exons that are long, that can bias the expression call
                        #if float(exp_val) != 0: ### Here, we use the original raw count data, whereas above is the adjusted quantile or raw count data
                        exp_list.append(exp_val); len_list.append(exp_dbase[exon][1]) ### This is for RNASeq -> don't include undetected exons - made in v.204
                    else: exp_list.append(exp_val) #elif float(exp_val) != 1:
                except KeyError: null =[] ###occurs if the expression exon list is missing some of these exons
            try:
                if len(exp_list)==0:
                    for exon in exon_list:
                        try:
                            exp_list.append(exp_dbase[exon][0][x]); len_list.append(exp_dbase[exon][1])
                            #kill
                        except KeyError: null=[] ### Gene entries will cause this error, since they are in the database but not in the count file
                if normalize_feature_exp == 'RPKM':
                    sum_const_exp=sum(map(float,exp_list)); gene_sum+=sum_const_exp
                    sum_length=sum(len_list) ### can have different lengths for each sample, since only expressed exons are considered
                    if excludeLowExp == False:
                        sum_length = gene_lengths[gene] ### Uses the all annotated exon lengths
                    ### Add only one avg-expression value for each array, this loop
                    try: steady_state_db[gene].append((sum_const_exp,sum_length))
                    except KeyError: steady_state_db[gene] = [(sum_const_exp,sum_length)]
                else:
                    avg_const_exp=Average(exp_list)
                    if avg_const_exp != 1: gene_sum+=avg_const_exp
                    ### Add only one avg-expression value for each array, this loop
                    try: steady_state_db[gene].append(avg_const_exp)
                    except KeyError: steady_state_db[gene] = [avg_const_exp]
            except Exception: null=[] ### Occurs when processing a truncated dataset (for testing usually) - no values for the gene should be included
            x += 1
        if gene_sum==0:
            try:
                del steady_state_db[gene] ### Hence, no genes showed evidence of expression (most critical for RNA-Seq)
            except Exception: null=[] ### Error occurs when a gene is added to the database from self.location_gene_db, but is not expressed
            
    return steady_state_db

def returnLargeGlobalVars():
    ### Prints all large global variables retained in memory (taking up space)
    all = [var for var in globals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all:
        try:
            if len(globals()[var])>1:
                print var, len(globals()[var])
        except Exception: null=[]
        
def clearObjectsFromMemory(db_to_clear):
    db_keys={}
    for key in db_to_clear: db_keys[key]=[]
    for key in db_keys:
        try: del db_to_clear[key]
        except Exception: 
            try:
                for i in key: del i ### For lists of tuples
            except Exception: del key ### For plain lists

def verifyFile(filename):
    status = 'not found'
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = 'found';break
    except Exception: status = 'not found'
    return status

def AppendOrWrite(export_path):
    export_path = filepath(export_path)
    status = verifyFile(export_path)
    if status == 'not found':
        export_data = export.ExportFile(export_path) ### Write this new file
    else:
        export_data = open(export_path,'a') ### Appends to existing file
        
    return export_data, status

def quantileNormalizationSimple(condition_count_db):
    ### Basic quantile normalization method (average ranked expression values)

    ### Get all junction or exon entries
    key_db={}
    for condition in condition_count_db:
        count_db = condition_count_db[condition]
        for key in count_db: key_db[key]=[]
            
    condition_unnormalized_db={}
    for key in key_db:
        ### Only look at the specific biotype of interest for each normalization
        for condition in condition_count_db:
            count_db = condition_count_db[condition]
            try:
                count = float(count_db[key])+1 ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
                count_db[key] = [] ### Set equal to null as a temporary measure to save memory
            except KeyError: count = 1.00 ###Was zero, but needs to be one for more realistic log2 fold calculations                
            ### store the minimal information to recover the original count and ID data prior to quantile normalization
            try: condition_unnormalized_db[condition].append([count,key])
            except Exception: condition_unnormalized_db[condition]=[[count,key]]
                    
    quantile_normalize_db={}; key_db={}
    for condition in condition_unnormalized_db:
        condition_unnormalized_db[condition].sort() ### Sort lists by count number
        rank=0 ### thus, the ID is the rank order of counts
        for (count,key) in condition_unnormalized_db[condition]:
            try: quantile_normalize_db[rank].append(count)
            except KeyError: quantile_normalize_db[rank] = [count]
            rank+=1
            
    ### Get the average value for each index
    for rank in quantile_normalize_db:
        quantile_normalize_db[rank] = Average(quantile_normalize_db[rank])

    for condition in condition_unnormalized_db:
        rank=0
        count_db = condition_count_db[condition]
        for (count,key) in condition_unnormalized_db[condition]:
            avg_count = quantile_normalize_db[rank]
            rank+=1
            count_db[key] = str(avg_count) ### re-set this value to the normalized value
    
    try:
        clearObjectsFromMemory(condition_unnormalized_db); condition_unnormalized_db =  []
        clearObjectsFromMemory(quantile_normalize_db); quantile_normalize_db = []
    except Exception: None
    
    return condition_count_db

def combineExonAnnotations(db):
    for i in db:
        list1=[]; list2=[]
        for (junctions,splice_event) in db[i]:
            list1.append(junctions); list2.append(splice_event)
        junctions = EnsemblImport.combineAnnotations(list1)
        splice_event = EnsemblImport.combineAnnotations(list2)
        db[i] = junctions,splice_event
    return db

def formatID(id):
    ### JunctionArray methods handle IDs with ":" different than those that lack this
    return string.replace(id,':','@')

def filterChromosomes(chromosome_names):
    ### If transcriptome only aligned to Ensembl reference, many chromosomes are not real
    updated_chromosomes=[]
    chr_count=0
    for chr in chromosome_names:
        if 'chr' in chr and len(chr)<7:
            chr_count+=1
            updated_chromosomes.append(chr)
    if chr_count>1:
        return updated_chromosomes
    else:
        return chromosome_names
        
def getChromosomeStrandCoordinates(species,testImport):
    ### For novel junctions with no known-splice site, map to genes
    gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,'RNASeq','key_by_array')
    
    chr_strand_gene_db = {}; location_gene_db = {}; chromosome_names={}; all_chromosomes={}
    for gene in gene_location_db:
        chr,strand,start,end = gene_location_db[gene]
        location_gene_db[chr,int(start),int(end)] = gene,strand
        try: chr_strand_gene_db[chr,strand].append((int(start),int(end)))
        except KeyError: chr_strand_gene_db[chr,strand] = [(int(start),int(end))]
        if testImport == 'yes':
            if chr=='chr1': chromosome_names[chr]=[]
            #if chr=='chr19': chromosome_names[chr]=[] ### Gene rich chromosome
            #if chr=='chrMT': chromosome_names[chr]=[] ### Gene rich chromosome
        elif len(chr)<7: chromosome_names[chr]=[]
        all_chromosomes[chr]=[]
    #chromosome_names = filterChromosomes(chromosome_names)
    ### Some organisms aren't organized into classical chromosomes (why I don't know)
    if len(chromosome_names)<10 and len(all_chromosomes)>9 and testImport=='no': chromosome_names = all_chromosomes
    return chr_strand_gene_db,location_gene_db,chromosome_names,gene_location_db

def exportDatasetLinkedExons(species,exons_to_export,critical_exon_annotations,root_dir,testImport=None,searchChr=None):
    export_path = root_dir+'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.txt'
    if searchChr != None:
        export_path = string.replace(export_path,'RNASeq/'+species,'RNASeq/exons/'+species)
        export_path = string.replace(export_path,'.txt','.'+searchChr+'.txt')
    
    if testImport == 'yes': print 'Writing',export_path
    export_data,status = AppendOrWrite(export_path)
    
    if status == 'not found': 
        export_title = ['AltAnalyzeID','exon_id','ensembl_gene_id','transcript_cluster_id','chromosome','strand','probeset_start','probeset_stop']
        export_title +=['class','constitutive_probeset','ens_exon_ids','ens_constitutive_status','exon_region','exon-region-start(s)','exon-region-stop(s)','splice_events','splice_junctions']
        export_title = string.join(export_title,'\t')+'\n'; export_data.write(export_title)

    ### We stored these in a dictionary to make sure each exon is written only once and so we can organize by gene
    exons_to_export_list=[]
    for key in exons_to_export:
        ed = exons_to_export[key]
        exons_to_export_list.append((key,ed))
    exons_to_export_list.sort()
    
    for (key,ed) in exons_to_export_list:
        constitutive_call = 'no'; ens_constitutive_status = '0'
        try:
            red = ed.ExonRegionData()
            exon_region = ed.ExonRegionID()
            start = str(ed.ReadStart()); stop = start
            if '-' not in exon_region and '_' not in exon_region: annotation = 'known'
            else: annotation = 'novel'      
        except Exception:
            red = ed ### For annotated exons, no difference in the annotations
            exon_region = ed.ExonRegionIDs()
            start = str(red.ExonStart()); stop = str(red.ExonStop())
            constitutive_call = red.Constitutive()
            if constitutive_call == 'yes': ens_constitutive_status = '1'
            annotation = 'known'            
        uid = red.GeneID()+':'+exon_region
        splice_events = red.AssociatedSplicingEvent(); splice_junctions = red.AssociatedSplicingJunctions()
        if uid in critical_exon_annotations:
            splice_junctions,splice_events = critical_exon_annotations[uid]
        export_values = [uid, exon_region, red.GeneID(), '', red.Chr(), red.Strand(), start, stop, annotation, constitutive_call, red.ExonID(), ens_constitutive_status]
        export_values+= [exon_region, str(red.ExonStart()), str(red.ExonStop()), splice_events, splice_junctions]
        export_values = string.join(export_values,'\t')+'\n'; export_data.write(export_values)
    export_data.close()
    
def exportNovelJunctions(species,novel_junction_db,condition_count_db,root_dir,dataset_name,biotype):
    
    if 'exp.' not in dataset_name: dataset_name = 'exp.'+dataset_name
    if '.txt' not in dataset_name: dataset_name+='.txt'
    dataset_name = string.replace(dataset_name,'exp','novel')
    dataset_name = string.replace(dataset_name,'.txt','.'+biotype+'.txt')
    
    export_path = root_dir+'ExpressionInput/'+dataset_name
    export_data,status = AppendOrWrite(export_path)
    
    if status == 'not found': 
        title = ['chr','strand','start','stop','start Ensembl','end Ensembl','known start', 'known end']
        for condition in condition_count_db: title.append(condition)
        export_data.write(string.join(title,'\t')+'\n')
    
    for key in novel_junction_db:
        ji = novel_junction_db[key]
        try: gene1 = str(ji.GeneID())
        except Exception: gene1=''
        try: gene2 = str(ji.SecondaryGeneID())
        except Exception: gene2 = 'None'
        try: le = str(ji.LeftExonAnnotations())
        except Exception: le = ''
        try: re = str(ji.RightExonAnnotations())
        except Exception: re = ''
        if biotype == 'junction':
            values = [ji.Chr(), ji.Strand(), str(ji.Exon1Stop()), str(ji.Exon2Start())]
        elif biotype == 'exon':
            values = [ji.Chr(), ji.Strand(), str(ji.Exon1Stop()-1), str(ji.Exon2Start()+1)] ### correct for initial adjustment
        values += [gene1,gene2,le,re]
        for condition in condition_count_db:
            count_db = condition_count_db[condition]
            try: read_count = count_db[key]
            except KeyError: read_count = '0'
            values.append(read_count)
        export_data.write(string.join(values,'\t')+'\n')
    export_data.close()

def exportDatasetLinkedGenes(species,gene_location_db,root_dir):
    """Include an entry for gene IDs to include constitutive expression for RPKM normalized data"""

    export_path = root_dir+'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_junctions.txt'
    export_data,status = AppendOrWrite(export_path)
    for gene in gene_location_db:
        chr,strand,start,end = gene_location_db[gene]
        export_values = [gene, 'E0.1',gene, '', chr, strand, str(start), str(end), 'known', 'yes', gene, '1']
        export_values+= ['E0.1', str(start), str(end), '', '']
        export_values = string.join(export_values,'\t')+'\n'; export_data.write(export_values)
    export_data.close()
               
def exportDatasetLinkedJunctions(species,junction_db,junction_annotations,root_dir,testImport=False,searchChr=None):
    export_path = root_dir+'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_junctions.txt'
    if searchChr != None:
        
        export_path = string.replace(export_path,'RNASeq/'+species,'RNASeq/junctions/'+species)
        export_path = string.replace(export_path,'.txt','.'+searchChr+'.txt')
        
    if testImport == 'yes': print 'Writing',export_path
    export_data,status = AppendOrWrite(export_path)

    if status == 'not found': 
        export_title = ['AltAnalyzeID','exon_id','ensembl_gene_id','transcript_cluster_id','chromosome','strand','probeset_start','probeset_stop']
        export_title +=['class','constitutive_probeset','ens_exon_ids','ens_constitutive_status','exon_region','exon-region-start(s)','exon-region-stop(s)','splice_events','splice_junctions']
        export_title = string.join(export_title,'\t')+'\n'; export_data.write(export_title)

    for key in junction_db:
        (chr,exon1_stop,exon2_start) = key
        ji=junction_db[key]
        #print key, ji.UniqueID(), ji.GeneID()
        if ji.GeneID()!=None and ji.UniqueID()!=None:
            if ji.UniqueID() in junction_annotations: ### Obtained from JunctionArray.inferJunctionComps()
                junctions,splice_events = junction_annotations[ji.UniqueID()]
                if ji.TransSplicing() == 'yes':
                    if len(splice_events)>0: splice_events+= '|trans-splicing'
                    else: splice_events = 'trans-splicing'
                ji.setAssociatedSplicingEvent(splice_events); ji.setAssociatedSplicingJunctions(junctions)
            elif ji.TransSplicing() == 'yes':
                ji.setAssociatedSplicingEvent('trans-splicing')
            try:
                try: constitutive_call = ji.Constitutive()
                except Exception:
                    jd = ji.ExonAnnotations()
                    constitutive_call = jd.Constitutive()
                if constitutive_call == 'yes': ens_constitutive_status = '1'
                else: ens_constitutive_status = '0'
                annotation = 'known'
            except Exception:
                constitutive_call = 'no'; ens_constitutive_status = '0'; annotation = 'novel'
            if 'I' in ji.ExonRegionID() or 'U' in ji.ExonRegionID() or '_' in ji.ExonRegionID():
                annotation = 'novel' ### Not previously indicated well (as I remember) for exon-level reads - so do this
            export_values = [ji.UniqueID(), ji.ExonRegionID(), ji.GeneID(), '', ji.Chr(), ji.Strand(), str(ji.Exon1Stop()), str(ji.Exon2Start()), annotation, constitutive_call, ji.ExonID(), ens_constitutive_status]
            export_values+= [ji.ExonRegionID(), str(ji.Exon1Stop()), str(ji.Exon2Start()), ji.AssociatedSplicingEvent(), ji.AssociatedSplicingJunctions()]
            export_values = string.join(export_values,'\t')+'\n'; export_data.write(export_values)
    export_data.close()
       
def combineDetectedExons(unmapped_exon_db,align_exon_db,novel_exon_db):
    ### Used for exon alignments (both start position and end position aligned to exon/intron/UTR regions)
    ### Reformat align_exon_db to easily lookup exon data
    aligned_exon_lookup_db={}
    for gene in align_exon_db:
        for ed in align_exon_db[gene]:
            aligned_exon_lookup_db[gene,ed.ReadStart()]=ed
            #if gene == 'ENSMUSG00000064181': print ed.ReadStart(),ed.ExonRegionID()
            
    ### Reformat novel_exon_db to easily lookup exon data - created from junction analysis (rename above exons to match novel junctions)
    novel_exon_lookup_db={}
    for gene in novel_exon_db:
        for ed in novel_exon_db[gene]:
            try:
                ### Only store exons that are found in the novel exon file
                null = aligned_exon_lookup_db[gene,ed.ReadStart()+1] ### offset introduced on import
                novel_exon_lookup_db[gene,ed.ReadStart()+1]=ed
            except Exception: null=[]
            try:
                ### Only store exons that are found in the novel exon file
                null = aligned_exon_lookup_db[gene,ed.ReadStart()-1] ### offset introduced on import
                novel_exon_lookup_db[gene,ed.ReadStart()-1]=ed
            except Exception: null=[]
                    
    ### Lookup the propper exon region ID and gene ID to format the unique ID and export coordinates
    x = 0
    for key in unmapped_exon_db:
        (chr,exon1_stop,exon2_start) = key
        ji=unmapped_exon_db[key]
        proceed = 'no'
        if ji.GeneID() != None:
            e1 = (ji.GeneID(),exon1_stop)
            e2 = (ji.GeneID(),exon2_start)
            exon_info=[]; override_annotation = None; found=[]
            try: null = aligned_exon_lookup_db[e1]; found.append(1)
            except Exception: null=[]
            try: null = aligned_exon_lookup_db[e2]; found.append(2)
            except Exception: null=[]
            try: null = novel_exon_lookup_db[e1]; override_annotation = 1
            except Exception:
                try: null = novel_exon_lookup_db[e2]; override_annotation = 2
                except Exception: null=[]          
            if len(found)>0:
                ### Below is not the simplist way to do this, but should be the fastest
                if 1 in found: exon_info.append(aligned_exon_lookup_db[e1])
                if 2 in found: exon_info.append(aligned_exon_lookup_db[e2])
                if len(exon_info) == 2: ed1,ed2 = exon_info
                else:
                    ed1 = exon_info[0]; ed2 = ed1; x+=1 ### if only one splice site aligned to a gene region (shouldn't occur)
                    if x == 2: null=[]; #print 'SOME EXONS FOUND WITH ONLY ONE ALIGNING POSITION...',key,ji.GeneID(),ed1.ExonRegionID(),e1,e2
                try: red1 = ed1.ExonRegionData(); red2 = ed2.ExonRegionData()
                except Exception:
                    """
                    print [ji.GeneID(), ji.Chr(), key]
                    print e1, e2
                    try: print ed1.ExonRegionData()
                    except Exception: 'ed1 failed'
                    try: print ed2.ExonRegionData()
                    except Exception: 'ed2 failed'
                    """
                    continue
                region1 = ed1.ExonRegionID(); region2 = ed2.ExonRegionID()
                #print region1,region2,ji.GeneID(),ji.Chr(),ji.Strand()
             
                try: splice_junctions = EnsemblImport.combineAnnotations([red1.AssociatedSplicingJunctions(),red2.AssociatedSplicingJunctions()])
                except Exception: print red1, red2;sys.exit()
                splice_events = EnsemblImport.combineAnnotations([red1.AssociatedSplicingEvent(),red2.AssociatedSplicingEvent()])
                ji.setAssociatedSplicingJunctions(splice_junctions)
                ji.setAssociatedSplicingEvent(splice_events)
                ens_exon_ids = EnsemblImport.combineAnnotations([red1.ExonID(),red2.ExonID()])
                ji.setExonID(ens_exon_ids)
                if red1.Constitutive() == 'yes' or red2.Constitutive() == 'yes': constitutive_call = 'yes'
                else: constitutive_call = 'no'
                ji.setConstitutive(constitutive_call)
                
                report_both_regions = 'no'
                try:
                    ### If the annotations are from a BED file produced by AltAnalyze, novel alternative splice sites may be present
                    ### if the below variable is not created, then this exon may over-ride the annotated exon region (e.g., E15.1 is over-written by E15.1_1234;E15.1_1256)
                    if 'ENS' in ji.JunctionID() and ':' not in ji.JunctionID(): report_both_regions = 'yes'
                except Exception: null=[]

                try:
                    ### If the annotations are from a BED file produced by AltAnalyze, it is possible for to a known exon to share a splice-site coordinate
                    ### with a novel junction exon. This will cause both to have the same override_annotation. Prevent this with the below 2nd override
                    if 'ENS' in ji.JunctionID() and ':' in ji.JunctionID(): override_annotation = None
                except Exception: null=[]
                
                if override_annotation != None:
                    if '_' in region1: region1 = string.split(region1,'_')[0]+'_'+str(int(string.split(region1,'_')[-1])-1)
                    if '_' in region2: region2 = string.split(region2,'_')[0]+'_'+str(int(string.split(region2,'_')[-1])+1)
                    if override_annotation == 1: region_id = region1 ### This forces a TopHat exon to be named for the splice-site position
                    else: region_id = region2                        
                else:
                    if report_both_regions == 'no':
                        ### Don't include specific start and end coordinates if inside a known exon
                        if ed1.AlignmentRegion() == 'exon': region1 = string.split(region1,'_')[0]
                        if ed2.AlignmentRegion() == 'exon': region2 = string.split(region2,'_')[0]
                        if ed1.AlignmentRegion() == 'full-intron' and ed2.AlignmentRegion() == 'full-intron':
                            region1 = string.split(region1,'_')[0]; region2 = string.split(region2,'_')[0]
                    
                    ### Below adjustmements need to compenstate for adjustments made upon import
                    if '_' in region1: region1 = string.split(region1,'_')[0]+'_'+str(int(string.split(region1,'_')[-1])-1)
                    if '_' in region2: region2 = string.split(region2,'_')[0]+'_'+str(int(string.split(region2,'_')[-1])+1)
                    
                ji.setExon1Stop(ji.Exon1Stop()-1); ji.setExon2Start(ji.Exon2Start()+1)
                if override_annotation != None: null=[] ### It is already assigned above
                elif region1 == region2: region_id = region1
                elif ji.Strand() == '+': region_id = region1+';'+region2
                else: region_id = region2+';'+region1 ### start and stop or genomically assigned
                uid = ji.GeneID()+':'+region_id
                #try: exon_region_db[ji.GeneID()].append((formatID(uid),region_id))
                #except KeyError: exon_region_db[ji.GeneID()]=[(formatID(uid),region_id)]
                ji.setExonRegionID(region_id)
                ji.setUniqueID(uid) ### hgu133
                ### Export format for new exons to add to the existing critical exon database (those in exon_region_db are combined with analyzed junctions)
                #exons_to_export[ji.GeneID(),region_id] = ji
            else:
                #print key, ji.GeneID(), ji.JunctionID(); sys.exit()
                null=[] ### Occurs because two genes are overlapping
    #return exons_to_export

def annotateNovelJunctions(novel_junction_db,novel_exon_db,exons_to_export):
    ### Reformat novel_exon_db to easily lookup exon data
    novel_exon_lookup_db={}
    for gene in novel_exon_db:
        for ed in novel_exon_db[gene]:
            novel_exon_lookup_db[gene,ed.ReadStart()]=ed
            
    ### Lookup the propper exon region ID and gene ID to format the unique ID and export coordinates
    junction_region_db={}
    unknown_gene_junctions={}
    for key in novel_junction_db:
        (chr,exon1_stop,exon2_start) = key
        ji=novel_junction_db[key]
        proceed = 'no'
        if ji.GeneID() != None:
            if ji.SpliceSitesFound() != 'both':
                e1 = (ji.GeneID(),exon1_stop)
                if ji.TransSplicing() == 'yes':
                    e2 = (ji.SecondaryGeneID(),exon2_start)
                else: e2 = (ji.GeneID(),exon2_start)
                if e1 in novel_exon_lookup_db and e2 in novel_exon_lookup_db:
                    proceed = 'yes'
                    try: ed1 = novel_exon_lookup_db[e1]; red1 = ed1.ExonRegionData(); gene1 = e1[0]
                    except Exception:
                        print chr, key, e1; kill
                    ed2 = novel_exon_lookup_db[e2]; red2 = ed2.ExonRegionData(); gene2 = e2[0]
                    ### If the splice-site was a match to a known junciton splice site, use it instead of that identified by exon-region location overlapp
                    if ji.LeftExonAnnotations() != None: region1 = ji.LeftExonAnnotations()
                    else: region1 = ed1.ExonRegionID(); exons_to_export[gene1,region1] = ed1
                    if ji.RightExonAnnotations() != None: region2 = ji.RightExonAnnotations()
                    else: region2 = ed2.ExonRegionID(); exons_to_export[gene2,region2] = ed2
                    #print region1,region2,ji.GeneID(),ji.Chr(),ji.Strand(), ji.LeftExonAnnotations(), ji.RightExonAnnotations()                    
            else:
                proceed = 'yes'
                region1 = ji.LeftExonAnnotations()
                region2 = ji.RightExonAnnotations()
                red1 = ji.LeftExonRegionData()
                red2 = ji.RightExonRegionData()
                ### Store the individual exons for export
                gene1 = ji.GeneID()
                if ji.TransSplicing() == 'yes': gene2 = ji.SecondaryGeneID()
                else: gene2 = ji.GeneID()
                exons_to_export[gene1,region1] = red1
                exons_to_export[gene2,region2] = red2
                
            if proceed == 'yes':
                try: splice_junctions = EnsemblImport.combineAnnotations([red1.AssociatedSplicingJunctions(),red2.AssociatedSplicingJunctions()])
                except Exception: print red1, red2;sys.exit()
                splice_events = EnsemblImport.combineAnnotations([red1.AssociatedSplicingEvent(),red2.AssociatedSplicingEvent()])
                ji.setAssociatedSplicingJunctions(splice_junctions)
                ji.setAssociatedSplicingEvent(splice_events)
                ens_exon_ids = EnsemblImport.combineAnnotations([red1.ExonID(),red2.ExonID()])
                ji.setExonID(ens_exon_ids)
                if ji.TransSplicing() == 'yes':
                    uid = ji.GeneID()+':'+region1+'-'+ji.SecondaryGeneID()+':'+region2
                    region_id = uid
                    ### When trans-splicing occurs, add the data twice to junction_region_db for the two different genes
                    ### in JunctionArray.inferJunctionComps, establish two separate gene junctions with a unique ID for the non-gene exon
                    try: junction_region_db[ji.GeneID()].append((formatID(uid),region1+'-'+'U1000.1_'+str(ji.Exon2Start())))
                    except KeyError: junction_region_db[ji.GeneID()]=[(formatID(uid),region1+'-'+'U1000.1_'+str(ji.Exon2Start()))]
                    try: junction_region_db[ji.SecondaryGeneID()].append((formatID(uid),'U0.1_'+str(ji.Exon1Stop())+'-'+region2))
                    except KeyError: junction_region_db[ji.SecondaryGeneID()]=[(formatID(uid),'U0.1_'+str(ji.Exon1Stop())+'-'+region2)]
                else:
                    uid = ji.GeneID()+':'+region1+'-'+region2
                    region_id = region1+'-'+region2
                    try: junction_region_db[ji.GeneID()].append((formatID(uid),region_id))
                    except KeyError: junction_region_db[ji.GeneID()]=[(formatID(uid),region_id)]
                ji.setExonRegionID(region_id)
                ji.setUniqueID(uid)
            else:
                unknown_gene_junctions[key]=[]
    return junction_region_db,exons_to_export

def alignReadsToExons(novel_exon_db,ens_exon_db,testImport=False):
    ### Simple method for aligning a single coordinate to an exon/intron region of an already matched gene
    examined_exons=0; aligned_exons=0
    for gene in ens_exon_db: #novel_exon_db
        try:
            region_numbers=[]; region_starts=[]; region_stops=[]
            for ed in novel_exon_db[gene]:
                examined_exons+=1; aligned_status=0; index=-1
                for rd in ens_exon_db[gene]:
                    index+=1 ### keep track of exon/intron we are in
                    region_numbers.append(int(string.split(rd.ExonRegionIDs()[1:],'.')[0]))
                    if rd.Strand() == '-': region_starts.append(rd.ExonStop()); region_stops.append(rd.ExonStart())
                    else: region_starts.append(rd.ExonStart()); region_stops.append(rd.ExonStop())
                    #print [rd.ExonStart(),rd.ExonStop(), rd.Strand()]
                    #print [ed.ReadStart(),rd.ExonStart(),rd.ExonStop()]
                    if ed.ReadStart()>=rd.ExonStart() and ed.ReadStart()<=rd.ExonStop():
                        ed.setAlignmentRegion('exon')
                        if 'I' in rd.ExonRegionIDs(): ### In an annotated intron
                            ed.setAlignmentRegion('intron')
                            ord = rd; updated = None
                            try: ### If the splice site is a novel 3' splice site then annotate as the 3' exon (less than 50nt away)
                                nrd = ens_exon_db[gene][index+1]
                                if (abs(ed.ReadStart()-nrd.ExonStart())<3) or (abs(ed.ReadStart()-nrd.ExonStop())<3):
                                    ed.setAlignmentRegion('full-intron') ### this is the start/end of intron coordinates
                                elif (abs(ed.ReadStart()-nrd.ExonStart())<50) or (abs(ed.ReadStart()-nrd.ExonStop())<50): rd = nrd; updated = 1
                            except Exception: null=[]
                            try:
                                prd = ens_exon_db[gene][index-1]
                                if (abs(ed.ReadStart()-prd.ExonStart())<3) or (abs(ed.ReadStart()-prd.ExonStop())<3):
                                    ed.setAlignmentRegion('full-intron')### this is the start/end of intron coordinates
                                elif (abs(ed.ReadStart()-prd.ExonStart())<50) or (abs(ed.ReadStart()-prd.ExonStop())<50):
                                    if updated==1: rd = ord; ###Hence the intron is too small to descriminate between alt5' and alt3' exons
                                    else: rd = prd
                            except Exception: null=[]                               
                        ed.setExonRegionData(rd); aligned_exons+=1; aligned_status=1
                        
                        if rd.ExonStop()==ed.ReadStart():
                            ed.setExonRegionID(rd.ExonRegionIDs())
                        elif rd.ExonStart()==ed.ReadStart():
                            ed.setExonRegionID(rd.ExonRegionIDs())
                        elif 'exon-intron' in ed.Annotation(): ### intron retention
                            ed.setExonRegionID(rd.ExonRegionIDs()) ### Hence there is a 1nt difference between read
                        else:
                            ed.setExonRegionID(rd.ExonRegionIDs()+'_'+str(ed.ReadStart()))
                        break
                if aligned_status == 0: ### non-exon/intron alinging sequences
                    region_numbers.sort(); region_starts.sort(); region_stops.sort()
                    if (rd.Strand() == '+' and ed.ReadStart()>=rd.ExonStop()) or (rd.Strand() == '-' and rd.ExonStop()>=ed.ReadStart()):
                        ### Applicable to 3'UTR (or other trans-splicing) aligning
                        utr_id = 'U'+str(region_numbers[-1])+'.1_'+str(ed.ReadStart())
                        ud = EnsemblImport.ExonAnnotationsSimple(rd.Chr(),rd.Strand(),region_stops[-1],region_stops[-1],gene,'','no',utr_id,'','')
                        ed.setExonRegionID(utr_id)
                    else:
                        ### Applicable to 5'UTR (or other trans-splicing) aligning
                        utr_id = 'U0.1'+'_'+str(ed.ReadStart())
                        ud = EnsemblImport.ExonAnnotationsSimple(rd.Chr(),rd.Strand(),region_starts[0],region_starts[0],gene,'','no',utr_id,'','')
                        ed.setExonRegionID(utr_id)
                    ed.setExonRegionData(ud)
                    ed.setAlignmentRegion('UTR')
        except Exception: null=[]
    if testImport == 'yes': print aligned_exons, 'splice sites aligned to exon region out of', examined_exons

def geneAlign(chr,chr_gene_locations,location_gene_db,chr_reads,switch_coord,read_aligned_to_gene):
    """ This function aligns the start or end position for each feature (junction or exon) to a gene, in two 
    steps by calling this function twice. In the second interation, the coordinates are reversed """
    
    index = 0 ### Don't examine genes already looked at
    genes_assigned = 0; trans_splicing=[]
    for (coord,ji) in chr_reads: ### junction coordinates or exon coordinates with gene object
        if index >5: index -=5 ### It is possible for some genes to overlap, so set back the index of genomically ranked genes each time
        gene_id_obtained = 'no'
        if switch_coord == 'no': rs,re=coord ### reverse the coordinates for the second iteration
        else: re,rs=coord ### first-interation coordinates (start and end)
        while index < len(chr_gene_locations):                 
            cs,ce = chr_gene_locations[index]
            #print [re,rs,cs,ce, ji.Chromosome()];sys.exit()
            ### Determine if the first listed coordinate lies within the gene
            if cs <= rs and ce >= rs:
                ### Yes, it does
                gene,strand = location_gene_db[chr,cs,ce]
                if switch_coord == 'yes': ### Only applies to coordinates, where the end-position didn't lie in the same gene as the start-position
                    if cs <= re and ce >= re:
                        ### This occurs when the first iteration detects a partial overlap, but the gene containing both coordinates is downstream
                        ### Hence, not trans-splicing
                        ji.setGeneID(gene)
                        break
                    first_geneid = ji.GeneID() ### see what gene was assigned in the first iteration (start position only)
                    #print ['trans',coord, first_geneid, gene]  ### Note: in rare cases, an exon can overlap with two genes (bad Ensembl annotations?)
                    ji.setTransSplicing()
                    side = ji.checkExonPosition(rs)
                    if side == 'left':
                        ji.setGeneID(gene)
                        ji.setSecondaryGeneID(first_geneid)
                    else:
                        ji.setSecondaryGeneID(gene)
                        #if ji.GeneID() == None: print 'B',coord, ji.GeneID(), secondaryGeneID()
                        #print ji.GeneID(), ji.SecondaryGeneID();kill
                    genes_assigned+=1; gene_id_obtained = 'yes'
                    ### Check to see if this gene represents a multi-gene spanning region (overlaps with multiple gene loci)
                    try:
                        ### This code was used to check and see if the gene is multi-spanning. Appears that the < sign is wrong > anyways, never go to the next gene unless the next read has passed it
                        #cs2,ce2 = chr_gene_locations[index+1]
                        #if cs2 < ce: index+=1 ### Continue analysis (if above is correct, the gene will have already been assigned)
                        #else: break
                        break
                    except Exception: break
                else:
                    ### First iteration, store the identified gene ID (only looking at the start position)
                    ji.setGeneID(gene);  gene_id_obtained = 'yes'
                    #print gene, rs, re, cs, ce
                    ### Check the end position, to ensure it is also lies within the gene region
                    if cs <= re and ce >= re:
                        genes_assigned+=1
                    else:
                        ### Hence, the end lies outside the gene region
                        trans_splicing.append((coord,ji))
                    ### Check to see if this gene represents a multi-gene spanning region (overlaps with multiple gene loci)
                    try:
                        ### This code was used to check and see if the gene is multi-spanning. Appears that the < sign is wrong > anyways, never go to the next gene unless the next read has passed it
                        #cs2,ce2 = chr_gene_locations[index+1]
                        #if cs2 < ce:  index+=1 ### Continue analysis (if above is correct, the gene will have already been assigned)
                        #else: break
                        break
                    except Exception: break
            else:
                if rs < ce and re < ce: break
                elif switch_coord == 'no' and cs <= re and ce >= re:
                    ### This can occur if the left junction splice site is in an exon and the other is the UTR as opposed to another gene
                    gene,strand = location_gene_db[chr,cs,ce]
                    ji.setSecondaryGeneID(gene); gene_id_obtained = 'yes'
                    #print gene, coord, ji.Strand(), ji.GeneID()
                index+=1
        if gene_id_obtained == 'no':
            ### These often appear to be genes predicted by tBLASTn at UCSC but not by Ensembl (e.g., chr17:27,089,652-27,092,318 mouse mm9)
            null=[]
            #ji.setGeneID(None) ### This is not necessary, since if one exon does not align to a gene it is still a valid alignment
            #print chr,coord

    read_aligned_to_gene += genes_assigned           
    #print genes_assigned, chr, 'Gene IDs assigned out of', len(chr_reads)
    #print len(trans_splicing),'reads with evidence of trans-splicing'
    
    ### For any coordinate-pair where the end-position doesn't lie within the same gene as the start, re-run for those to see which gene they are in
    if switch_coord == 'no' and len(trans_splicing)>0:
        read_aligned_to_gene = geneAlign(chr,chr_gene_locations,location_gene_db,trans_splicing,'yes',read_aligned_to_gene)
    return read_aligned_to_gene


def getNovelExonCoordinates(species,root_dir):
    """ Currently, any novel exon determined during initial RNA-Seq read annotation with defined start and end coordinates, only has 
    the exon-end coordinate, not start, in it's name. However, the start and stop are indicated in the counts.Experiment.txt file.
    To get this, we parse that file and only store exons with an I or U in them and then correct for this in the matching function below """    
    exp_dir = root_dir+'/ExpressionInput/'
    dir_list = read_directory(exp_dir)
    counts_file = None
    for file in dir_list:
        if 'counts.' in file and 'steady' not in file:
            counts_file = file
    
    ### Example

    #ENSG00000137076:I17.1_35718353=chr9:35718353-35718403 (novel exon coordinates - just sorted, not necessarily in the correct order)
    #ENSG00000137076:E17.1-I17.1_35718403=chr9:35718809-35718403 (5' supporting junction)
    #ENSG00000137076:I17.1_35718353-E18.1=chr9:35718353-35717783 (3' supporting junction)
    #here, once we see that I17.1_35718353 is the exon ID, we know we need to get the function with -I17.1_35718403 (always the second value)
    
    if counts_file!=None:
        fn=filepath(exp_dir+counts_file)
        print 'Reading counts file'
        novel_exon_db = parseCountFile(fn,'exons',{}) ### Get novel exons
        print 'Reading counts file'
        novel_exon_db = parseCountFile(fn,'junctions',novel_exon_db) ### Get novel exons
    return novel_exon_db
    
def getMaxCounts(fn,cutoff,filterExport=False,filterExportDir=False):
    firstLine=True
    expressed_uids={}
    
    if filterExport != False:
        eo=export.ExportFile(filterExportDir)
        
    for line in open(fn,'rU').xreadlines():
        Line = cleanUpLine(line)
        t = string.split(Line,'\t')
        key = t[0]
        if firstLine:
            firstLine = False
            if filterExport != False:
                eo.write(line)
        else:
            if filterExport != False:
                if key in filterExport:
                    eo.write(line)
            else:
                try: uid, coordinates = string.split(key,'=')
                except Exception: uid = key
                try: maxExp = max(map(lambda x: float(x), t[1:])); #print maxExp;sys.exit()
                except Exception:
                    #print t[1:];sys.exit()
                    if 'NA' in t[1:]:
                        tn = [0 if x=='NA' else x for x in t[1:]] ### Replace NAs
                        maxExp = max(map(lambda x: float(x), tn))
                    elif '' in t[1:]:
                        tn = [0 if x=='' else x for x in t[1:]] ### Replace blanks
                        maxExp = max(map(lambda x: float(x), tn))
                    else:
                        maxExp=cutoff+1
    
                #gene = string.split(uid,':')[0]
                if maxExp > cutoff:
                    expressed_uids[uid] = []
    return expressed_uids
        
def importBiologicalRelationships(species):
    ### Combine non-coding Ensembl gene annotations with UniProt functional annotations
    import ExpressionBuilder
    custom_annotation_dbase={}
    try: coding_db = ExpressionBuilder.importTranscriptBiotypeAnnotations(species)
    except Exception: coding_db = {}
    try: gene_to_symbol_db = ExpressionBuilder.importGeneAnnotations(species)
    except Exception: gene_to_symbol_db = {}
    for gene in coding_db:
        #coding_type = string.split(coding_db[gene][-1],'|')
        coding_type = coding_db[gene][-1]
        if 'protein_coding' in coding_type:
            coding_type = 'protein_coding'
        else:
            coding_type = 'ncRNA'
        if gene in gene_to_symbol_db:
            symbol = string.lower(gene_to_symbol_db[gene][0])
            ### The below genes cause issues with many single cell datasets in terms of being highly correlated
            if 'rpl'==symbol[:3] or 'rps'==symbol[:3] or 'mt-'==symbol[:3] or '.' in symbol or 'gm'==symbol[:2]:
                coding_type = 'ncRNA'
        try: gene_db = custom_annotation_dbase[coding_type]; gene_db[gene]=[]
        except Exception: custom_annotation_dbase[coding_type] = {gene:[]}
        
    filename = 'AltDatabase/uniprot/'+species+'/custom_annotations.txt'
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        ens_gene,compartment,custom_class = t[:3]
        if 'GPCR' in custom_class:
            custom_class = ['GPCR']
        else:
            custom_class = string.split(custom_class,'|')
        custom_class = string.split(compartment,'|')+custom_class
        for cc in custom_class:
            try: gene_db = custom_annotation_dbase[cc]; gene_db[ens_gene]=[]
            except Exception: custom_annotation_dbase[cc] = {ens_gene:[]}
    #custom_annotation_dbase={}
    try:
        filename = 'AltDatabase/goelite/'+species+'/gene-mapp/Ensembl-BioMarkers.txt'
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            gene,null,celltype = t[:3]
            try: gene_db = custom_annotation_dbase['BioMarker']; gene_db[gene]=[]
            except Exception: custom_annotation_dbase['BioMarker'] = {gene:[]}
        print len(custom_annotation_dbase), 'gene classes imported'
    except Exception: pass
    return custom_annotation_dbase

def importGeneSets(geneSetType,filterType=None,geneAnnotations=None):
    gene_db={}
    if 'Ontology' in geneSetType:
        filename = 'AltDatabase/goelite/'+species+'/nested/Ensembl_to_Nested-GO.txt'
        ontology=True
    else:
        filename = 'AltDatabase/goelite/'+species+'/gene-mapp/Ensembl-'+geneSetType+'.txt'
        ontology=False
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if ontology:
            gene,category = t
        else: gene,null,category = t[:3]
        if filterType==None:
            try: gene_db[gene].append(category)
            except Exception: gene_db[gene] = [category]
        elif filterType in category:
            if gene in geneAnnotations:
                gene = geneAnnotations[gene][0]
            gene_db[gene]=[]
    return gene_db
    
def singleCellRNASeqWorkflow(Species, platform, expFile, mlp, exp_threshold=5, rpkm_threshold=5, drivers=False, parameters = None, reportOnly=False):
    global species
    global rho_cutoff
    species = Species
    removeOutliers = False
    if parameters != None:
        rpkm_threshold = parameters.ExpressionCutoff()
        exp_threshold = parameters.CountsCutoff()
        rho_cutoff = parameters.RhoCutoff()
        restrictBy = parameters.RestrictBy()
        try: removeOutliers = parameters.RemoveOutliers()
        except Exception: pass
        if platform == 'exons' or platform == 'PSI':
            rpkm_threshold=0
            exp_threshold=0
    else:
        rho_cutoff = 0.4
        restrictBy = 'protein_coding'
    onlyIncludeDrivers=True
                
    if platform != 'exons' and platform != 'PSI':
        platform = checkExpressionFileFormat(expFile,platform)

    if platform != 'RNASeq':
        if rpkm_threshold>1.9999:
            rpkm_threshold = math.log(rpkm_threshold,2) ### log2 transform
        
    if removeOutliers:
        ### Remove samples with low relative number of genes expressed
        try:
            import shutil
            print '***Removing outlier samples***'
            from import_scripts import sampleIndexSelection
            reload(sampleIndexSelection)
            output_file = expFile[:-4]+'-OutliersRemoved.txt'
            sampleIndexSelection.statisticallyFilterFile(expFile,output_file,rpkm_threshold)
            if 'exp.' in expFile:
                ### move the original groups and comps files
                groups_file = string.replace(expFile,'exp.','groups.')
                groups_file = string.replace(groups_file,'-steady-state','')
                groups_filtered_file = groups_file[:-4]+'-OutliersRemoved.txt'
                #comps_file = string.replace(groups_file,'groups.','comps.')
                #comps_filtered_file = string.replace(groups_filtered_file,'groups.','comps.')
                #counts_file = string.replace(expFile,'exp.','counts.')
                #counts_filtered_file = string.replace(output_file,'exp.','counts.')
                try: shutil.copyfile(groups_file,groups_filtered_file) ### if present copy over
                except Exception: pass
                try: shutil.copyfile(comps_file,comps_filtered_file) ### if present copy over
                except Exception: pass
                #try: shutil.copyfile(counts_file,counts_filtered_file) ### if present copy over
                #except Exception: pass
            expFile = output_file
            print ''
        except Exception:
            print '***Filtering FAILED***'
            print traceback.format_exc()
            
    expressed_uids_rpkm = getMaxCounts(expFile,rpkm_threshold)
    
    try: expressed_uids_counts = getMaxCounts(string.replace(expFile,'exp.','counts.'),exp_threshold)
    except Exception: expressed_uids_counts=expressed_uids_rpkm
    
    if len(expressed_uids_counts) > 0:
        try: expressed_uids = expressed_uids_rpkm.viewkeys() & expressed_uids_counts.viewkeys() ### common
        except Exception: expressed_uids = getOverlappingKeys(expressed_uids_rpkm,expressed_uids_counts)
    else:
        expressed_uids = expressed_uids_rpkm
    print 'Genes filtered by counts:',len(expressed_uids_counts)
    print 'Genes filtered by expression:',len(expressed_uids_rpkm),len(expressed_uids)
    #expressed_uids = filterByProteinAnnotation(species,expressed_uids)
    print len(expressed_uids), 'expressed genes by RPKM/TPM (%d) and counts (%d)' % (rpkm_threshold,exp_threshold)
    #"""

    from import_scripts import OBO_import; import ExpressionBuilder
    gene_to_symbol_db = ExpressionBuilder.importGeneAnnotations(species)
    
    try: biological_categories = importBiologicalRelationships(species)
    except Exception:
        restrictBy = None
        biological_categories={}
        print 'Missing annotation file in:','AltDatabase/uniprot/'+species+'/custom_annotations.txt !!!!!'

    if restrictBy !=None:
        print 'Attempting to restrict analysis to protein coding genes only (flag --RestrictBy protein_coding)'
        genes = biological_categories['protein_coding']
        genes_temp=dict(genes)
        for gene in genes_temp:
            if gene in gene_to_symbol_db:
                genes[gene_to_symbol_db[gene][0]]=[] ### add symbols
        genes_temp={}
    else:
        genes = {}
        for i in expressed_uids: genes[i]=[]
    """
    genes.update(biological_categories['BioMarker'])
    genes.update(biological_categories['transcription regulator'])
    genes.update(biological_categories['splicing regulator'])
    genes.update(biological_categories['kinase'])
    genes.update(biological_categories['GPCR'])
    """
    expressed_uids_db={}; guide_genes={}
    for id in expressed_uids: expressed_uids_db[id]=[]
    if platform == 'exons' or platform == 'PSI': ### For splicing-index value filtering
        expressed_uids=[]
        for uid in expressed_uids_db:
            geneID = string.split(uid,':')[0]
            geneID = string.split(geneID,' ')[-1]
            if geneID in genes: expressed_uids.append(uid)
    else:
        try: expressed_uids = genes.viewkeys() & expressed_uids_db.viewkeys() ### common
        except Exception: expressed_uids = getOverlappingKeys(genes,expressed_uids_db)

        #print len(expressed_uids)
        expressed_uids_db2={}
        for id in expressed_uids: expressed_uids_db2[id]=[]
        
        if drivers != False:
            guide_genes = getDrivers(drivers)
            if onlyIncludeDrivers:
                try: expressed_uids = guide_genes.viewkeys() & expressed_uids_db2.viewkeys() ### common
                except Exception: expressed_uids = getOverlappingKeys(guide_genes,expressed_uids_db2)

    if len(expressed_uids)<10:
        print 'NOTE: The input IDs do not sufficiently map to annotated protein coding genes...',
        print 'skipping protein coding annotation filtering.'
        expressed_uids=[]
        for uid in expressed_uids_db:
            expressed_uids.append(uid)
    print len(expressed_uids), 'expressed IDs being further analyzed'
    #sys.exit()
    #Meenakshi changes
    print_out,n = findCommonExpressionProfiles(expFile,species,platform,expressed_uids,guide_genes,mlp,parameters=parameters,reportOnly=reportOnly)
    return print_out,n

def getOverlappingKeys(db1,db2):
    db3=[]
    for key in db1:
        if key in db2:
            db3.append(key)
    return db3

def getDrivers(filename):
    fn = filepath(filename)
    firstLine=True
    drivers={}
    for line in open(fn,'rU').xreadlines():
        line = line.rstrip()
        t = string.split(line,'\t')
        if firstLine: firstLine = False
        else:
            gene = t[0]
            drivers[gene]=[]
    print 'Imported %d guide genes' % len(drivers)
    return drivers
            
def filterByProteinAnnotation(species,expressed_uids):
    import ExpressionBuilder
    custom_annotation_dbase = ExpressionBuilder.importTranscriptBiotypeAnnotations(species)
    expressed_uids_protein=[]
    for gene in expressed_uids:
        if gene in custom_annotation_dbase:
            compartment,custom_class = custom_annotation_dbase[gene]
            if 'protein_coding' in custom_class:
                expressed_uids_protein.append(gene)
    if len(expressed_uids_protein)>10:
        return expressed_uids_protein
    else:
        return expressed_uids

def CoeffVar(expFile,platform,expressed_uids,fold=2,samplesDiffering=2,guideGenes=[]):
    firstLine=True
    expressed_values={}
    expressed_values_filtered={}
    cv_list=[]
    for line in open(expFile,'rU').xreadlines():
        key = string.split(line,'\t')[0]
        t = string.split(line,'\t')
        if firstLine:
            headers = line
            firstLine = False
        else:
            try: uid, coordinates = string.split(key,'=')
            except Exception: uid = key
            values = map(lambda x: float(x), t[1:])
            #gene = string.split(uid,':')[0]
            if uid in expressed_uids:
                vs = list(values); vs.sort()
                cv = statistics.stdev(values)/statistics.avg(values)
                if samplesDiffering<1: samplesDiffering=1
                if platform == 'RNASeq':
                    if (vs[-1*samplesDiffering]/vs[samplesDiffering])>fold: ### Ensures that atleast 4 samples are significantly different in the set
                        expressed_values[uid] = values
                        cv_list.append((cv,uid))
                else:
                    if (vs[-1*samplesDiffering]-vs[samplesDiffering])>fold: ### Ensures that atleast 4 samples are significantly different in the set
                        expressed_values[uid] = values
                        cv_list.append((cv,uid))
            if uid in guideGenes:
                expressed_values[uid] = values
                cv_list.append((10000,uid)) ### Very high CV
    cv_list.sort()
    cv_list.reverse()
    x=0
    for (cv,uid) in cv_list:
        x+=1
        """
        if uid == 'ENSMUSG00000003882':
            print x, 'ilr7'
        """

    for (cv,uid) in cv_list[:5000]:
        expressed_values_filtered[uid] = expressed_values[uid]
    return expressed_values_filtered, fold, samplesDiffering, headers

def determinePattern(vs):
    max_vs = max(vs)
    min_vs = min(vs)
    lower_max = max_vs - (max_vs*0.01)
    upper_min = abs(max_vs)*0.01
    s = bisect.bisect_right(vs,upper_min) ### starting low 15% index position
    e = bisect.bisect_left(vs,lower_max) ### ending upper 85% index position
    #print vs
    #print max_vs, min_vs
    #print lower_max, upper_min
    #print s, e
    avg =  statistics.avg(vs[s:e+1])
    m = bisect.bisect_left(vs,avg)
    ratio = vs[m]/vs[((e-s)/2)+s-2] ### If the ratio is close to 1, a sigmoidal or linear pattern likely exists
    print ratio
    #sys.exit()
    return ratio

def checkExpressionFileFormat(expFile,platform):
    firstLine=True
    inputMax=0; inputMin=10000
    expressed_values={}
    for line in open(expFile,'rU').xreadlines():
        key = string.split(line,'\t')[0]
        t = string.split(line,'\t')
        if firstLine:
            headers = line
            firstLine = False
        else:
            try: uid, coordinates = string.split(key,'=')
            except Exception: uid = key
            try: values = map(lambda x: float(x), t[1:])
            except Exception:
                values=[]
                for value in t[1:]:
                    try: values.append(float(value))
                    except Exception:pass
            try:
                if max(values)>inputMax: inputMax = max(values)
            except Exception:
                pass
                
    if inputMax>100: ### Thus, not log values
        platform = 'RNASeq'
    else:
        platform = "3'array"
    return platform
    
def optimizeNumberOfGenesForDiscovery(expFile,platform,expressed_uids,fold=2,samplesDiffering=2,guideGenes=[]):
    #import warnings
    #warnings.simplefilter('error', UserWarning)
    firstLine=True
    expressed_values={}
    for line in open(expFile,'rU').xreadlines():
        key = string.split(line,'\t')[0]
        t = string.split(line,'\t')
        if firstLine:
            headers = line
            firstLine = False
        else:
            try: uid, coordinates = string.split(key,'=')
            except Exception: uid = key
            try: values = map(lambda x: float(x), t[1:])
            except Exception:
                values = t[1:]
                if 'NA' in values:
                    values = [0 if x=='NA' else x for x in values] ### Replace NAs
                    values = map(lambda x: float(x), values)
                else:  
                    values=[]
                    for value in t[1:]:
                        try: values.append(float(value))
                        except Exception: values.append(-9999)
                    values = numpy.ma.masked_values(values, -9999.)
            #gene = string.split(uid,':')[0]
            #if uid == 'ENSMUSG00000041515': print 'IRF8'
            if uid in expressed_uids:
                #slope_exp_ratio = determinePattern(vs)
                #if slope_exp_ratio<2 and slope_exp_ratio>0.5:
                if platform == 'RNASeq':
                    try: values = map(lambda x: math.log(x+1,2),values)
                    except Exception:
                        if 'NA' in values:
                            values = [0 if x=='NA' else x for x in values] ### Replace NAs
                            values = map(lambda x: math.log(x+1,2),values)
                        elif '' in values:
                            values = [0 if x=='' else x for x in values] ### Replace NAs
                            values = map(lambda x: math.log(x+1,2),values)
                    vs = list(values); vs.sort()
                    if (vs[-1*samplesDiffering]-vs[samplesDiffering-1])>math.log(fold,2): ### Ensures that atleast 4 samples are significantly different in the set
                        expressed_values[uid] = values
                else:
                    vs = list(values); vs.sort()
                    if (vs[-1*samplesDiffering]-vs[samplesDiffering-1])>math.log(fold,2): ### Ensures that atleast 4 samples are significantly different in the set
                        expressed_values[uid] = values
                if uid in guideGenes:
                    expressed_values[uid] = values
                #if uid == 'ENSMUSG00000062825': print (vs[-1*samplesDiffering]-vs[samplesDiffering]),math.log(fold,2);sys.exit()

    print len(expressed_uids),'genes examined and', len(expressed_values),'genes expressed for a fold cutoff of', fold
    if len(expressed_uids)==0 or len(expressed_values)==0:
        print options_result_in_no_genes
    elif len(expressed_uids) < 50 and len(expressed_values)>0:
        return expressed_values, fold, samplesDiffering, headers
    elif len(expressed_values)>15000:
        if platform == 'exons' or platform == 'PSI':
            fold+=0.1
        else:
            fold+=1
            samplesDiffering+=1
        expressed_values, fold, samplesDiffering, headers = optimizeNumberOfGenesForDiscovery(expFile,platform,expressed_uids,fold=fold,samplesDiffering=samplesDiffering,guideGenes=guideGenes)
    elif fold == 1.2 and samplesDiffering == 1:
        return expressed_values, fold, samplesDiffering, headers
    elif len(expressed_values)<50:
        fold-=0.2
        samplesDiffering-=1
        if samplesDiffering<1: samplesDiffering = 1
        if fold < 1.1: fold = 1.2
        expressed_values, fold, samplesDiffering, headers = optimizeNumberOfGenesForDiscovery(expFile,platform,expressed_uids,fold=fold,samplesDiffering=samplesDiffering,guideGenes=guideGenes)
    else:
        return expressed_values, fold, samplesDiffering, headers
    return expressed_values, fold, samplesDiffering, headers

def intraCorrelation(expressed_values,mlp):
    if mlp.cpu_count() < 3:
        processors = mlp.cpu_count()
    else: processors = 8
    pool = mlp.Pool(processes=processors)
    si = (len(expressed_values)/processors)
    s = si; b=0
    db_ls=[]
    if len(expressed_values)<10: forceError ### will si to be zero and an infanite loop
    while s<len(expressed_values):
        db_ls.append(dict(expressed_values.items()[b:s]))
        b+=si; s+=si
    db_ls.append(dict(expressed_values.items()[b:s]))

    ### Create an instance of MultiZscoreWorker (store the variables to save memory)
    workerMulti = MultiCorrelatePatterns(expressed_values)
    
    results = pool.map(workerMulti,db_ls)
    #for i in db_ls: workerMulti(i)
    pool.close(); pool.join(); pool = None
    correlated_genes={}
    for a in results:
        for k in a: correlated_genes[k] = a[k]
    return correlated_genes

def findCommonExpressionProfiles(expFile,species,platform,expressed_uids,guide_genes,mlp,fold=2,samplesDiffering=2,parameters=None,reportOnly=False):
    use_CV=False
    global rho_cutoff

    row_metric = 'correlation'; row_method = 'average'
    column_metric = 'cosine'; column_method = 'hopach'
    original_column_metric = column_metric
    original_column_method = column_method
    color_gradient = 'yellow_black_blue'; transpose = False; graphic_links=[]
    if parameters != None:
        try: excludeGuides = parameters.ExcludeGuides() ### Remove signatures
        except Exception: excludeGuides = None
        fold = parameters.FoldDiff()
        samplesDiffering = parameters.SamplesDiffering()
        amplifyGenes = parameters.amplifyGenes()
        if 'Guide' in parameters.GeneSelection():
            amplifyGenes = False ### This occurs when running ICGS with the BOTH option, in which Guide3 genes are retained - ignore these
            parameters.setGeneSelection('')
            parameters.setClusterGOElite('')
        excludeCellCycle = parameters.ExcludeCellCycle()
        from visualization_scripts import clustering
        row_metric = 'correlation'; row_method = 'average'
        column_metric = parameters.ColumnMetric(); column_method = parameters.ColumnMethod()
        original_column_metric = column_metric
        original_column_method = column_method
        color_gradient = 'yellow_black_blue'; graphic_links=[]
        if platform == 'exons' or platform =='PSI': color_gradient = 'yellow_black_blue'
        guide_genes = parameters.JustShowTheseIDs()
        cell_cycle_id_list = []
    else:
        amplifyGenes = False
        excludeCellCycle = False
    
    if platform != 'exons'and platform !='PSI':
        platform = checkExpressionFileFormat(expFile,platform)
    else:
        if LegacyMode: pass
        else:
            fold = math.pow(2,0.5)
        fold = 1.25
    
    #"""
    if use_CV:
        expressed_values, fold, samplesDiffering, headers = CoeffVar(expFile,platform,expressed_uids,fold=2,samplesDiffering=2,guideGenes=guide_genes)
    else:
        print 'Finding an optimal number of genes based on differing thresholds to include for clustering...'
        #fold=1; samplesDiffering=1
        expressed_values, fold, samplesDiffering, headers = optimizeNumberOfGenesForDiscovery(expFile,platform,expressed_uids,fold=fold,samplesDiffering=samplesDiffering,guideGenes=guide_genes) #fold=2,samplesDiffering=2
        print 'Evaluating',len(expressed_values),'genes, differentially expressed',fold,'fold for at least',samplesDiffering*2,'samples'
    #sys.exit()
    
    from import_scripts import OBO_import; import ExpressionBuilder
    gene_to_symbol_db = ExpressionBuilder.importGeneAnnotations(species)
    symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol_db)
    
    areYouSure=False
    if (excludeCellCycle == 'strict' or excludeCellCycle == True) and areYouSure:
        cc_param = copy.deepcopy(parameters)
        cc_param.setPathwaySelect('cell cycle')
        cc_param.setGeneSet('GeneOntology')
        cc_param.setGeneSelection('amplify')
        transpose = cc_param
        filtered_file = export.findParentDir(expFile)+'/amplify/'+export.findFilename(expFile)
        writeFilteredFile(filtered_file,platform,headers,{},expressed_values,[])
        if len(expressed_values)<1000:
            row_method = 'hopach'; row_metric = 'correlation'
        if column_method != 'hopach': row_method = 'average' ### needed due to PC errors
        if len(headers)>7000: ### For very ultra-large datasets
            column_method = 'average'
        cc_graphic_links = clustering.runHCexplicit(filtered_file, graphic_links, row_method, row_metric, column_method, column_metric, color_gradient, transpose, display=False, Normalize=True, JustShowTheseIDs=guide_genes)
        cell_cycle_id_list = genericRowIDImport(string.replace(cc_graphic_links[0][-1],'.png','.txt'))
        expressed_values2 = {}
        for id in expressed_values:
            try: symbolID = gene_to_symbol_db[id][0]
            except Exception: symbolID = id
            if id not in cell_cycle_id_list and symbolID not in cell_cycle_id_list:
                expressed_values2[id]=expressed_values[id]
        print len(expressed_values)-len(expressed_values2),'cell-cycle associated genes removed for cluster discovery'
        expressed_values = expressed_values2
    print 'amplifyGenes:',amplifyGenes
    
    ### Write out filtered list to amplify and to filtered.YourExperiment.txt
    filtered_file = export.findParentDir(expFile)+'/amplify/'+export.findFilename(expFile)
    groups_file = string.replace(expFile,'exp.','groups.')
    groups_filtered_file = string.replace(filtered_file,'exp.','groups.')
    groups_file = string.replace(groups_file,'-steady-state','')
    groups_filtered_file = string.replace(groups_filtered_file,'-steady-state','')
    try: export.customFileCopy(groups_file,groups_filtered_file) ### if present copy over
    except Exception: pass
    writeFilteredFile(filtered_file,platform,headers,{},expressed_values,[])
    filtered_file_new = string.replace(expFile,'exp.','filteredExp.')
    try: export.customFileCopy(filtered_file,filtered_file_new) ### if present copy over
    except Exception: pass
    
    if reportOnly:
        print_out = '%d genes, differentially expressed %d fold for at least %d samples' % (len(expressed_values), fold, samplesDiffering*2)
        return print_out

    if len(expressed_values)<1400 and column_method == 'hopach':
        row_method = 'hopach'; row_metric = 'correlation'
    else:
        row_method = 'weighted'; row_metric = 'cosine'
    if amplifyGenes:
        transpose = parameters
        try:
            if len(parameters.GeneSelection())>0:
                parameters.setGeneSelection(parameters.GeneSelection()+' amplify')
                print 'Finding correlated genes to the input geneset(s)...'
            else:
                print 'Finding intra-correlated genes from the input geneset(s)...'
                parameters.setGeneSelection(parameters.GeneSelection()+' IntraCorrelatedOnly amplify')
        except Exception:         
            parameters.setGeneSelection(parameters.GeneSelection()+' IntraCorrelatedOnly amplify')
            print 'Finding intra-correlated genes from the input geneset(s)...'
        if column_method != 'hopach': row_method = 'average' ### needed due to PC errors
        graphic_links = clustering.runHCexplicit(filtered_file, graphic_links, row_method, row_metric, column_method, column_metric, color_gradient, transpose, display=False, Normalize=True, JustShowTheseIDs=guide_genes)
        #return graphic_links
        from visualization_scripts import clustering
        matrix, column_header, row_header, dataset_name, group_db = clustering.importData(graphic_links[-1][-1][:-4]+'.txt')
        headers = ['UID']+column_header
        expressed_values2={}
        for i in row_header: ### Filter the expressed values for the intra-correlated queried gene set and replace
            try: expressed_values2[i]=expressed_values[i]
            except Exception:
                try:
                    e = symbol_to_gene[i][0]
                    expressed_values2[e]=expressed_values[e]
                except Exception: 
                    pass
        expressed_values = expressed_values2
        
    print 'Looking for common gene expression profiles for class assignment...',
    begin_time = time.time()

    useNumpyCorr=True
    negative_rho = rho_cutoff*-1
    
    #results_file = string.replace(expFile[:-4]+'-CORRELATED-FEATURES.txt','exp.','/SamplePrediction/')
    #eo = export.ExportFile(results_file[:-4]+'-genes.txt')
    
    if useNumpyCorr:
        row_ids=[]
        x = []
        for id in expressed_values:
            row_ids.append(id)
            x.append(expressed_values[id])
            #if id== 'Bcl2l11': print expressed_values[id];sys.exit()
        D1 = numpy.corrcoef(x)
        print 'initial correlations obtained'
        i=0
        correlated_genes={}

        if 'exons' == platform or 'PSI' == platform:
            for score_ls in D1:
                proceed = True
                correlated = []
                geneID = row_ids[i]
                refgene = string.split(geneID,':')[0]
                k=0
                if excludeGuides!=None:
                    if geneID in excludeGuides: ### skip this main event
                        proceed=False
                        continue
                for v in score_ls:            
                    if v>rho_cutoff:# or v<negative_rho:
                        if refgene not in row_ids[k]:
                            correlated.append((v,row_ids[k]))
                            if excludeGuides!=None:
                                if row_ids[k] in excludeGuides: ### skip this main event
                                    proceed=False
                                    break
                    k+=1
                correlated.sort()
                if LegacyMode == False:
                    correlated.reverse()
                if proceed:
                    correlated = map(lambda x:x[1],correlated)
                    correlated_genes[geneID] = correlated
                i+=1
        else:     
            for score_ls in D1:
                correlated = []
                geneID = row_ids[i]
                k=0; temp=[]
                for v in score_ls:
                    if v>rho_cutoff:# or v<negative_rho:
                        #scores.append((v,row_ids[k]))
                        correlated.append((v,row_ids[k]))
                        #temp.append((geneID,row_ids[k],str(v)))
                    k+=1
                correlated.sort()
                if LegacyMode == False:
                    correlated.reverse()
                correlated = map(lambda x:x[1],correlated)
                if len(correlated)>0:
                    correlated_genes[geneID] = correlated
                    #for (a,b,c) in temp: eo.write(a+'\t'+b+'\t'+c+'\n')
                    
                i+=1
    else:
        ### Find common patterns now
        performAllPairwiseComparisons = True
        if performAllPairwiseComparisons:
            correlated_genes = intraCorrelation(expressed_values,mlp)
            print len(correlated_genes), 'highly correlated genes found for downstream clustering.'
        else: correlated_genes={}
        
    atleast_10={}
    if len(correlated_genes)<70: connections = 0
    elif len(correlated_genes)<110: connections = 4
    else: connections = 5

    numb_corr=[]
    for i in correlated_genes:
        if len(correlated_genes[i])>connections:
            numb_corr.append([len(correlated_genes[i]),i])
            atleast_10[i]=correlated_genes[i] ### if atleast 10 genes apart of this pattern
            x=0
            for k in correlated_genes[i]:
                if x<30: ### cap it at 30
                    try: atleast_10[k]=correlated_genes[k] ### add all correlated keys and values
                    except Exception: pass
                x+=1

    if len(atleast_10)<30:
        print 'Initial correlated set too small, getting anything correlated'
        for i in correlated_genes:
            if len(correlated_genes[i])>0:
                numb_corr.append([len(correlated_genes[i]),i])
                try: atleast_10[i]=correlated_genes[i] ### if atleast 10 genes apart of this pattern
                except Exception: pass
                for k in correlated_genes[i]:
                    try: atleast_10[k]=correlated_genes[k] ### add all correlated keys and values
                    except Exception: pass

    if len(atleast_10) == 0:
        atleast_10 = expressed_values
        
    #eo.close()
    print len(atleast_10), 'genes correlated to multiple other members (initial filtering)'
    ### go through the list from the most linked to the least linked genes, only reported the most linked partners
    if len(atleast_10)>5000:
        print_out=""
        return print_out,atleast_10
        
    removeOutlierDrivenCorrelations=True
    exclude_corr=[]
    numb_corr.sort(); numb_corr.reverse()
    numb_corr2=[]
    #print len(numb_corr)
    if removeOutlierDrivenCorrelations and samplesDiffering != 1:
        for key in numb_corr: ### key gene
            associations,gene = key
            temp_corr_matrix_db={}; rows=[]; temp_corr_matrix=[]
            gene_exp_vals = list(expressed_values[gene]) ### copy the list
            max_index = gene_exp_vals.index(max(gene_exp_vals))
            del gene_exp_vals[max_index]
            #temp_corr_matrix.append(exp_vals); rows.append(gene)
            #if 'ENSG00000016082' in correlated_genes[gene] or 'ENSG00000016082' == gene: print gene_to_symbol_db[gene],associations
            if gene not in exclude_corr:
                #print len(correlated_genes[gene])
                for k in correlated_genes[gene]:
                    exp_vals = list(expressed_values[k]) ### copy the list
                    #print exp_vals
                    del exp_vals[max_index]
                    #temp_corr_matrix.append(exp_vals); rows.append(gene)
                    #print exp_vals,'\n'
                    temp_corr_matrix_db[k]=exp_vals
                    temp_corr_matrix.append(exp_vals); rows.append(gene)
                correlated_hits = pearsonCorrelations(gene_exp_vals,temp_corr_matrix_db)
                try: avg_corr = numpyCorrelationMatrix(temp_corr_matrix,rows,gene)
                except Exception: avg_corr = 0
                #if gene_to_symbol_db[gene][0] == 'ISL1' or gene_to_symbol_db[gene][0] == 'CD10' or gene_to_symbol_db[gene][0] == 'POU3F2':

                if len(correlated_hits)>0:
                    if LegacyMode:
                        if (float(len(correlated_hits))+1)/len(correlated_genes[gene])<0.5 or avg_corr<rho_cutoff:  ### compare to the below
                            pass
                        else:
                            numb_corr2.append([len(correlated_hits),gene])
                    else:
                        if (float(len(correlated_hits))+1)/len(correlated_genes[gene])<0.5 or avg_corr<(rho_cutoff-0.1):
                            #exclude_corr.append(key)
                            #if gene == 'XXX': print len(correlated_hits),len(correlated_genes[gene]), avg_corr, rho_cutoff-0.1
                            pass
                        else:
                            numb_corr2.append([len(correlated_hits),gene])
                #print (float(len(correlated_hits))+1)/len(correlated_genes[gene]), len(correlated_genes[gene]), key
                
        numb_corr = numb_corr2
        numb_corr.sort(); numb_corr.reverse()
    
    #print len(numb_corr)
    exclude_corr={}; new_filtered_set={}
    limit=0
    for key in numb_corr: ### key gene
        associations,gene = key
        #if 'ENSG00000016082' in correlated_genes[gene] or 'ENSG00000016082' == gene: print gene_to_symbol_db[gene],associations
        if gene not in exclude_corr:
            for k in correlated_genes[gene]: 
                exclude_corr[k]=[]
                new_filtered_set[k]=[]
            new_filtered_set[gene]=[]
            limit+=1
            #print key
            #if limit==1: break
    atleast_10 = new_filtered_set
        
    addMultipleDrivers=True
    if len(guide_genes)>0 and addMultipleDrivers: ### Artificially weight the correlated genes with known biological driverse
        for gene in guide_genes:
            y=1
            while y<2:
                if y==1:
                    try: atleast_10[gene]=expressed_values[gene]
                    except Exception: break
                else:
                    try: atleast_10[gene+'-'+str(y)]=expressed_values[gene]
                    except Exception: break
                    expressed_values[gene+'-'+str(y)]=expressed_values[gene] ### Add this new ID to the database
                    #print gene+'-'+str(y)
                y+=1
    
    #atleast_10 = expressed_values

    results_file = string.replace(expFile[:-4]+'-CORRELATED-FEATURES.txt','exp.','/SamplePrediction/')
    writeFilteredFile(results_file,platform,headers,gene_to_symbol_db,expressed_values,atleast_10)
    
    print len(atleast_10),'final correlated genes'
    end_time = time.time()
    print 'Initial clustering completed in',int(end_time-begin_time),'seconds'
    
    
    results_file = string.replace(expFile[:-4]+'-CORRELATED-FEATURES.txt','exp.','/SamplePrediction/')

    if len(atleast_10)<1200 and column_method == 'hopach':
        row_method = 'hopach'; row_metric = 'correlation'
    else:
        if LegacyMode:
            row_method = 'average'; row_metric = 'euclidean'
        else:
            row_method = 'weighted'; row_metric = 'cosine'
    #print row_method, row_metric
    correlateByArrayDirectly = False
    if correlateByArrayDirectly:
        from visualization_scripts import clustering
    
        matrix, column_header, row_header, dataset_name, group_db = clustering.importData(results_file)
        
        new_column_header = map(lambda x: int(x[5:]),column_header)
        matrix = [new_column_header]+matrix
        matrix = zip(*matrix) ### transpose
        exp_sample_db={}
        for sample_data in matrix:
            exp_sample_db[sample_data[0]] = sample_data[1:]
        
        correlated_arrays = intraCorrelation(exp_sample_db,mpl)
        print len(correlated_arrays), 'highly correlated arrays from gene subsets.'
        mimum_corr_arrays={}
        for i in correlated_arrays:
            if len(correlated_arrays[i])>1:
                linked_lists=correlated_arrays[i]+[i]
                for k in correlated_arrays[i]:
                    linked_lists+=correlated_arrays[k]
                linked_lists = unique.unique(linked_lists)
                linked_lists.sort()
               # print len(linked_lists), linked_lists
    else:
        try:
            from visualization_scripts import clustering
            if platform == 'exons': color_gradient = 'yellow_black_blue'
            transpose = False
            if column_method != 'hopach': row_method = 'average' ### needed due to PC errors (possibly outside of LegacyMode)
            graphic_links = clustering.runHCexplicit(results_file, graphic_links, row_method, row_metric, column_method, column_metric, color_gradient, transpose, display=False, Normalize=True, JustShowTheseIDs=guide_genes)
            if len(graphic_links)==0:
                graphic_links = clustering.runHCexplicit(results_file, graphic_links, row_method, row_metric, column_method, column_metric, color_gradient, transpose, display=False, Normalize=True, JustShowTheseIDs=guide_genes)
            cluster_file = string.replace(graphic_links[0][1],'.png','.txt')
        except Exception: pass
        #exportGroupsFromClusters(cluster_file,expFile,platform)
    #"""
    #filtered_file = export.findParentDir(expFile)+'/amplify/'+export.findFilename(expFile)
    #graphic_links = [(1,'/Users/saljh8/Desktop/Grimes/KashishNormalization/test/ExpressionInput/SamplePrediction/DataPlots/Clustering-CombinedSingleCell_March_15_2015-CORRELATED-FEATURES-hierarchical_cosine_euclidean.txt')]

    try: graphic_links,new_results_file = correlateClusteredGenes(platform,graphic_links[-1][-1][:-4]+'.txt',numSamplesClustered=samplesDiffering,excludeCellCycle=excludeCellCycle,graphics=graphic_links,ColumnMethod=column_method)
    except Exception: print traceback.format_exc()
        
    row_metric = 'correlation'; row_method = 'hopach'
    #column_metric = 'cosine'
    #if LegacyMode: column_method = 'hopach'

    cellCycleRemove1=[]; cellCycleRemove2=[]
    try:
        newDriverGenes1, cellCycleRemove1 = correlateClusteredGenes(platform,graphic_links[-1][-1][:-4]+'.txt',stringency='strict',numSamplesClustered=samplesDiffering,excludeCellCycle=excludeCellCycle,ColumnMethod=column_method)
        newDriverGenes1_str = 'Guide1 '+string.join(newDriverGenes1.keys(),' ')+' amplify positive'
        parameters.setGeneSelection(newDriverGenes1_str) ### force correlation to these targetGenes
        parameters.setGeneSet('None Selected') ### silence this
        parameters.setPathwaySelect('None Selected')
        if column_method != 'hopach': row_method = 'average' ### needed due to PC errors
        graphic_links = clustering.runHCexplicit(filtered_file, graphic_links, row_method, row_metric, column_method, column_metric, color_gradient, parameters, display=False, Normalize=True)
        
        newDriverGenes2, cellCycleRemove2 = correlateClusteredGenes(platform,graphic_links[-1][-1][:-4]+'.txt',stringency='strict',numSamplesClustered=samplesDiffering,excludeCellCycle=excludeCellCycle,ColumnMethod=column_method)
        newDriverGenes2_str = 'Guide2 '+string.join(newDriverGenes2.keys(),' ')+' amplify positive'
        parameters.setGeneSelection(newDriverGenes2_str) ### force correlation to these targetGenes
        parameters.setGeneSet('None Selected') ### silence this
        parameters.setPathwaySelect('None Selected')
        graphic_links = clustering.runHCexplicit(filtered_file, graphic_links, row_method, row_metric, column_method, column_metric, color_gradient, parameters, display=False, Normalize=True)
        newDriverGenes3 = unique.unique(newDriverGenes1.keys()+newDriverGenes2.keys())
        cellCycleRemove=cellCycleRemove1+cellCycleRemove2 ### It is possible for a cell cycle guide-gene to be reported in both guide1 and 2, but only as cell cycle associated in one of them
        newDriverGenes3_filtered=[]
        for i in newDriverGenes3:
            if not i in cellCycleRemove:
                newDriverGenes3_filtered.append(i)
        newDriverGenes3_str = 'Guide3 '+string.join(newDriverGenes3_filtered,' ')+' amplify positive'
        parameters.setGeneSelection(newDriverGenes3_str)
        try:
            parameters.setClusterGOElite('BioMarkers')
            """
            if species == 'Mm' or species == 'Hs' or species == 'Rn':
                parameters.setClusterGOElite('BioMarkers')
            else:
                parameters.setClusterGOElite('GeneOntology')
            """
        except Exception, e: 
            print e
        graphic_links = clustering.runHCexplicit(filtered_file, graphic_links, row_method, row_metric, column_method, column_metric, color_gradient, parameters, display=False, Normalize=True)
    except Exception:
        print traceback.format_exc()

    try: copyICGSfiles(expFile,graphic_links)
    except Exception: pass
    return graphic_links,len(atleast_10)

def copyICGSfiles(expFile,graphic_links):
    if 'ExpressionInput' in expFile:
        root_dir = string.split(expFile,'ExpressionInput')[0]
    else:
        root_dir = string.split(expFile,'AltResults')[0]
    import shutil
    destination_folder = root_dir+'/ICGS'
    try: os.mkdir(destination_folder)
    except Exception: pass
    for (order,png) in graphic_links:
        file = export.findFilename(png)
        txt = string.replace(file,'.png','.txt')
        pdf = string.replace(file,'.png','.pdf')
        dest_png = destination_folder+'/'+file
        dest_txt = destination_folder+'/'+txt
        dest_pdf = destination_folder+'/'+pdf
        shutil.copy(png, dest_png)
        shutil.copy(png[:-4]+'.txt', dest_txt)
        shutil.copy(png[:-4]+'.pdf', dest_pdf)

def pearsonCorrelations(ref_gene_exp,exp_value_db):
    correlated=[]
    for gene in exp_value_db:
        rho,p = stats.pearsonr(ref_gene_exp,exp_value_db[gene])
        if rho>rho_cutoff or rho<(rho_cutoff*-1):
            if rho!= 1:
                correlated.append(gene)
    #print len(exp_value_db),len(correlated);sys.exit()
    return correlated

def numpyCorrelationMatrix(x,rows,gene):
    D1 = numpy.corrcoef(x)
    gene_correlations={}
    i=0
    scores = []
    for score_ls in D1:
        for v in score_ls:
            scores.append(v)
    return numpy.average(scores)

def numpyCorrelationMatrixCount(x,rows,cutoff=0.4,geneTypeReport=None):
    ### Find which genes are most correlated
    D1 = numpy.corrcoef(x)
    gene_correlation_counts={}
    i=0
    for score_ls in D1:
        correlated_genes=[]
        geneID = rows[i]
        k=0; genes_to_report=[]
        for rho in score_ls:
            if rho>cutoff:
                correlated_genes.append(rows[k])
                if rows[k] in geneTypeReport:
                    genes_to_report.append(rows[k])
            k+=1
        
        gene_correlation_counts[geneID]=len(correlated_genes),genes_to_report
        i+=1
        
    return gene_correlation_counts
    
def numpyCorrelationMatrixGene(x,rows,gene):
    D1 = numpy.corrcoef(x)
    
    gene_correlations={}
    i=0
    for score_ls in D1:
        scores = []
        geneID = rows[i]
        k=0
        for v in score_ls:
            scores.append((v,rows[k]))
            k+=1
        scores.sort()
        gene_correlations[geneID] = scores
        i+=1
         
    correlated_genes={}
    rho_values = map(lambda (r,g): r,gene_correlations[gene])
    genes = map(lambda (r,g): g,gene_correlations[gene])
    s1 = bisect.bisect_right(rho_values,rho_cutoff)
    s2 = bisect.bisect_left(rho_values,-1*rho_cutoff)
    correlated = genes[:s2] ### for the right bisect, remove self correlations with -1
    correlated = genes[s1:] ### for the left bisect, remove self correlations with -1
    #print len(rows), len(correlated);sys.exit()
    return len(correlated)/len(rows)
                
def numpyCorrelationMatrixGeneAlt(x,rows,genes,gene_to_symbol,rho_cutoff):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides import warnings
        D1 = numpy.ma.corrcoef(x)
    i=0
    gene_correlations={}
    for score_ls in D1:
        scores = []
        try: symbol = gene_to_symbol[rows[i]][0]
        except Exception: symbol = '$'
        if rows[i] in genes or symbol in genes:
            k=0
            for v in score_ls:
                if str(v)!='nan':
                    if v > rho_cutoff:
                        uid = rows[k]
                        if uid in gene_to_symbol: uid = gene_to_symbol[uid][0]
                        scores.append((v,uid))
                k+=1    
            scores.sort()
            scores.reverse()
            scores = map(lambda x: x[1], scores[:140]) ### grab the top 140 correlated gene symbols only
            if len(symbol)==1: symbol = rows[i]
            gene_correlations[symbol] = scores 
        i+=1
    return gene_correlations

def genericRowIDImport(filename):
    id_list=[]
    for line in open(filename,'rU').xreadlines():
        uid = string.split(line,'\t')[0]
        if ' ' in uid:
            for id in string.split(uid,' '):
                id_list.append(id)
        else:
            id_list.append(uid)
    return id_list

def writeFilteredFile(results_file,platform,headers,gene_to_symbol_db,expressed_values,atleast_10,excludeGenes=[]):
    eo = export.ExportFile(results_file)
    try: headers = string.replace(headers,'row_clusters-flat','UID')
    except Exception:
        headers = string.join(headers,'\t')+'\n'
        headers = string.replace(headers,'row_clusters-flat','UID')
    eo.write(headers)
    keep=[]; sort_genes=False
    e=0
    if len(atleast_10)==0:
        atleast_10 = expressed_values
        sort_genes = True
    for i in atleast_10:
        if i in gene_to_symbol_db:
            symbol = gene_to_symbol_db[i][0]
        else: symbol = i
        if i not in excludeGenes and symbol not in excludeGenes:
            if i not in keep:
                keep.append((symbol,i))
    
    if sort_genes:
        keep.sort(); keep.reverse()
        
    for (symbol,i) in keep:
        """
        if platform == 'RNASeq':
            values = map(lambda x: logTransform(x), expressed_values[i])
        else:
        """
        values = map(str,expressed_values[i])
        eo.write(string.join([symbol]+values,'\t')+'\n')
        e+=1
    eo.close()

def remoteGetDriverGenes(Species,platform,results_file,numSamplesClustered=3,excludeCellCycle=False,ColumnMethod='hopach'):
    global species
    species = Species
    guideGenes, cellCycleRemove = correlateClusteredGenes(platform,results_file,stringency='strict',excludeCellCycle=excludeCellCycle,ColumnMethod=ColumnMethod)
    guideGenes = string.join(guideGenes.keys(),' ')+' amplify positive'
    return guideGenes
    
def correlateClusteredGenes(platform,results_file,stringency='medium',numSamplesClustered=3,
                excludeCellCycle=False,graphics=[],ColumnMethod='hopach',rhoCutOff=0.2, transpose=False,
                includeMoreCells=False):
    
    if numSamplesClustered<1: numSamplesClustered=1
        ### Get all highly variably but low complexity differences, typically one or two samples that are really different
    if stringency == 'medium':
        new_results_file = string.replace(results_file,'.txt','-filtered.txt')
        new_results_file = string.replace(new_results_file,'.cdt','-filtered.txt')
        eo = export.ExportFile(new_results_file)
        medVarHighComplexity=[]; medVarLowComplexity=[]; highVarHighComplexity=[]; highVarLowComplexity=[]
        if transpose==False or includeMoreCells:
            medVarLowComplexity, column_header = correlateClusteredGenesParameters(results_file,rho_cutoff=0.3,hits_cutoff=3,hits_to_report=6,transpose=transpose)
            medVarHighComplexity, column_header = correlateClusteredGenesParameters(results_file,rho_cutoff=0.1,hits_cutoff=3,hits_to_report=6,transpose=transpose) #hits_cutoff=6
            highVarLowComplexity, column_header = correlateClusteredGenesParameters(results_file,rho_cutoff=0.5,hits_cutoff=1,hits_to_report=4,transpose=transpose)
            highVarHighComplexity, column_header = correlateClusteredGenesParameters(results_file,rho_cutoff=0.2,hits_cutoff=1,hits_to_report=6,filter=True,numSamplesClustered=numSamplesClustered,transpose=transpose)
        else:
            highVarLowComplexity, column_header = correlateClusteredGenesParameters(results_file,rho_cutoff=0.5,hits_cutoff=1,hits_to_report=4,transpose=transpose)
        #combined_results = dict(medVarLowComplexity.items() + medVarLowComplexity.items() + highVarLowComplexity.items() + highVarHighComplexity.items())
        combined_results={}
        
        for i in medVarLowComplexity: combined_results[i]=[]
        for i in medVarHighComplexity: combined_results[i]=[]
        for i in highVarLowComplexity: combined_results[i]=[]
        for i in highVarHighComplexity: combined_results[i]=[]
        #combined_results = highVarHighComplexity
    if stringency == 'strict':
        medVarLowComplexity, column_header = correlateClusteredGenesParameters(results_file,rho_cutoff=0.3,hits_cutoff=4,hits_to_report=50,filter=True,numSamplesClustered=numSamplesClustered)
        medVarHighComplexity, column_header = correlateClusteredGenesParameters(results_file,rho_cutoff=0.1,hits_cutoff=4,hits_to_report=50,filter=True,numSamplesClustered=numSamplesClustered) #hits_cutoff=6
        highVarLowComplexity, column_header = correlateClusteredGenesParameters(results_file,rho_cutoff=0.5,hits_cutoff=3,hits_to_report=50,filter=True,numSamplesClustered=numSamplesClustered)
        highVarHighComplexity, column_header = correlateClusteredGenesParameters(results_file,rho_cutoff=0.3,hits_cutoff=3,hits_to_report=50,filter=True,numSamplesClustered=numSamplesClustered)
        #combined_results = dict(medVarLowComplexity.items() + medVarLowComplexity.items() + highVarLowComplexity.items() + highVarHighComplexity.items())
        combined_results={}
        for i in medVarLowComplexity: combined_results[i]=[]
        for i in medVarHighComplexity: combined_results[i]=[]
        for i in highVarLowComplexity: combined_results[i]=[]
        for i in highVarHighComplexity: combined_results[i]=[]
        guideGenes, addition_cell_cycle_associated = correlateClusteredGenesParameters(results_file,rho_cutoff=rhoCutOff,hits_cutoff=0,hits_to_report=1,geneFilter=combined_results,excludeCellCycle=excludeCellCycle)
        if guideGenes == 'TooFewBlocks':
            guideGenes, addition_cell_cycle_associated = correlateClusteredGenesParameters(results_file,rho_cutoff=rhoCutOff+0.1,hits_cutoff=0,hits_to_report=1,geneFilter=combined_results,excludeCellCycle=excludeCellCycle)
            if guideGenes == 'TooFewBlocks':
                guideGenes, addition_cell_cycle_associated = correlateClusteredGenesParameters(results_file,rho_cutoff=rhoCutOff+0.2,hits_cutoff=0,hits_to_report=1,geneFilter=combined_results,excludeCellCycle=excludeCellCycle,forceOutput=True)
        if len(guideGenes)>200:
            print 'Too many guides selected (>200)... performing more stringent filtering...'
            guideGenes, addition_cell_cycle_associated = correlateClusteredGenesParameters(results_file,rho_cutoff=0.1,hits_cutoff=0,hits_to_report=1,geneFilter=combined_results,excludeCellCycle=excludeCellCycle,restrictTFs=True)
        return guideGenes, addition_cell_cycle_associated
    #B4galt6, Prom1
    for tuple_ls in combined_results:
        data_length = len(tuple_ls);break
    if data_length == len(column_header):
        eo.write(string.join(column_header,'\t')+'\n')
    else:
        eo.write(string.join(['UID']+column_header,'\t')+'\n')
    
    #combined_results = highVarHighComplexity
    for tuple_ls in combined_results:
        eo.write(string.join(list(tuple_ls),'\t')+'\n')
    eo.close()

    cluster = True
    if cluster == True and transpose==False:
        from visualization_scripts import clustering
        if ColumnMethod == 'hopach':
            row_method = 'hopach'
            column_method = 'hopach'
        else:
            column_method = ColumnMethod
            row_method = 'average'
        row_metric = 'correlation'
        column_metric = 'cosine'
        color_gradient = 'yellow_black_blue'
        if platform == 'exons': color_gradient = 'yellow_black_blue'
        transpose = False
        try:
            len(guide_genes)
        except Exception:
            guide_genes = []
 
        graphics = clustering.runHCexplicit(new_results_file, graphics, row_method, row_metric, column_method, column_metric, color_gradient, transpose, display=False, Normalize=True, JustShowTheseIDs=guide_genes)
        cluster_file = string.replace(graphics[0][1],'.png','.txt')
        #exportGroupsFromClusters(cluster_file,expFile,platform)
    return graphics, new_results_file

def exportReDefinedClusterBlocks(results_file,block_db,rho_cutoff):
    ### Re-import the matrix to get the column cluster IDs
    matrix, column_header, row_header, dataset_name, group_db, priorColumnClusters, priorRowClusters = clustering.remoteImportData(results_file)
    
    new_block_db = {}
    centroid_blocks=[]
    centroids = []
    for block in block_db:
        if len(block_db[block])>3:
            new_block_db[block] = block_db[block]  ### Keep track of the row_header indexes associated with each blcok
            data = map(lambda x: matrix[x],block_db[block])
            ### Compute an expression centroid from the block (cluster)
            centroid = [float(sum(col))/len(col) for col in zip(*data)]
            centroids.append(centroid)
            centroid_blocks.append(block)
            
    ### Compare block centroids
    D1 = numpy.corrcoef(centroids)
    i=0
    correlated_blocks=[]
    for score_ls in D1:
        scores = []
        block = centroid_blocks[i]
        k=0
        for v in score_ls:
            if str(v)!='nan' and v>0.6:
                if block !=centroid_blocks[k]:
                    blocks = [block,centroid_blocks[k]]
                    blocks.sort()
                    if blocks not in correlated_blocks:
                        correlated_blocks.append(blocks)
            k+=1    
        i+=1
    
    newBlock=0
    existing=[]
    updated_blocks={}
    correlated_blocks.sort()
    print correlated_blocks
    ### Build a tree of related blocks (based on the code in junctionGraph)
    for (block1,block2) in correlated_blocks:
        if block1 not in existing and block2 not in existing:
            newBlock=newBlock+1
            updated_blocks[newBlock]=[block1,]
            updated_blocks[newBlock].append(block2)
            existing.append(block1)
            existing.append(block2)
                
        elif block1 in existing and block2 not in existing:
            for i in updated_blocks:
                if block1 in updated_blocks[i]:
                    updated_blocks[i].append(block2)
                    existing.append(block2)
        elif block2 in existing and block1 not in existing:
            for i in updated_blocks:
                if block2 in updated_blocks[i]:
                    updated_blocks[i].append(block1)
                    existing.append(block1)   
        elif block1 in existing and block2 in existing:
            for i in updated_blocks:
                if block1 in updated_blocks[i]:
                        b1=i
                if block2 in updated_blocks[i]:
                    b2=i
            if b1!=b2:
                for b in updated_blocks[b2]:
                    if b not in updated_blocks[b1]:
                        updated_blocks[b1].append(b)
                del updated_blocks[b2]
                
    ### Add blocks not correlated to other blocks (not in correlated_blocks)
    #print len(existing),len(centroid_blocks)
    print updated_blocks
    for block in centroid_blocks:
        if block not in existing:
            newBlock+=1
            updated_blocks[newBlock]=[block]
    
    import collections
    row_order = collections.OrderedDict()
    for newBlock in updated_blocks:
        events_in_block=0
        for block in updated_blocks[newBlock]:
            for i in new_block_db[block]:
                events_in_block+=1
        if events_in_block>5:
            for block in updated_blocks[newBlock]:
                for i in new_block_db[block]:
                    row_order[i] = newBlock ### i is a row_header index - row_header[i] is a UID
                    #if newBlock==3:
                    #if row_header[i]=='TAF2&ENSG00000064313&E9.1-I9.1_120807184__ENSG00000064313&E9.1-E10.1':
                    #print row_header[i]
    print updated_blocks

    ### Non-clustered block results - Typically not used by good to refer back to when testing
    original_block_order = collections.OrderedDict()
    for block in new_block_db:
        for i in new_block_db[block]:
            original_block_order[i]=block
    #row_order = original_block_order
    
    ### Export the results
    row_header.reverse() ### Reverse order is the default
    priorColumnClusters = map(str,priorColumnClusters)
    new_results_file = results_file[:-4]+'-BlockIDs.txt'
    eo = export.ExportFile(new_results_file)
    eo.write(string.join(['UID','row_clusters-flat']+column_header,'\t')+'\n')
    eo.write(string.join(['column_clusters-flat','']+priorColumnClusters,'\t')+'\n')
    
    for i in row_order:
        cluster_number = str(row_order[i])
        uid = row_header[i]
        values = map(str,matrix[i])
        eo.write(string.join([uid,cluster_number]+values,'\t')+'\n')
    eo.close()
    
    print 'Filtered, grouped expression clusters exported to:',new_results_file

def correlateClusteredGenesParameters(results_file,rho_cutoff=0.3,hits_cutoff=4,hits_to_report=5,
            filter=False,geneFilter=None,numSamplesClustered=3,excludeCellCycle=False,restrictTFs=False,
            forceOutput=False,ReDefinedClusterBlocks=False,transpose=False):
    from visualization_scripts import clustering
    addition_cell_cycle_associated=[]
    
    if geneFilter != None:
        geneFilter_db={}
        for i in geneFilter:
            geneFilter_db[i[0]]=[]
        geneFilter=geneFilter_db
    
    matrix, column_header, row_header, dataset_name, group_db = clustering.importData(results_file,geneFilter=geneFilter)

    if transpose: ### If performing reduce cluster heterogeneity on cells rather than on genes
        #print 'Transposing matrix'
        matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
        column_header, row_header = row_header, column_header
        
    Platform = None
    for i in row_header:
        if 'ENS' in i and '-' in i and ':' in i: Platform = 'exons'
            
    #print hits_to_report
    if hits_to_report == 1:
        ### Select the best gene using correlation counts and TFs
        try:
            from import_scripts import OBO_import; import ExpressionBuilder
            gene_to_symbol_db = ExpressionBuilder.importGeneAnnotations(species)
            symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol_db)            
            try: TFs = importGeneSets('Biotypes',filterType='transcription regulator',geneAnnotations=gene_to_symbol_db)
            except Exception: TFs = importGeneSets('BioTypes',filterType='transcription regulator',geneAnnotations=gene_to_symbol_db)
            if excludeCellCycle == True or excludeCellCycle == 'strict':
                try: cell_cycle = importGeneSets('KEGG',filterType='Cell cycle:',geneAnnotations=gene_to_symbol_db)
                except Exception:
                    cell_cycle = {}
                try: cell_cycle_go = importGeneSets('GeneOntology',filterType='GO:0022402',geneAnnotations=gene_to_symbol_db)
                except Exception: cell_cycle_go={}
                for i in cell_cycle_go:
                    cell_cycle[i]=[]
                print len(cell_cycle),'cell cycle genes being considered.'
            else:
                cell_cycle={}

        except Exception:
            print traceback.format_exc()
            symbol_to_gene={}; TFs={}; cell_cycle={}
        gene_corr_counts = numpyCorrelationMatrixCount(matrix,row_header,cutoff=0.4,geneTypeReport=TFs)
        
    #try: column_header = map(lambda x: string.split(x,':')[1],column_header[1:])
    #except Exception: column_header = column_header[1:]

    i=0
    block=0
    
    if ReDefinedClusterBlocks:
        import collections
        block_db=collections.OrderedDict() ### seems benign but could alter legacy results
    else:
        block_db={}
    for row in matrix:
        if i!=0:
            rho,p = stats.pearsonr(row,matrix[i-1]) ### correlate to the last ordered row
            #if row_header[i] == 'Pax6': print [block],row_header[i-1],rho,rho_cutoff
            """
            try:
                if row_header[i] in guide_genes: print row_header[i], rho
                if row_header[i-1] in guide_genes: print row_header[i-1], rho
                if row_header[i+1] in guide_genes: print row_header[i+1], rho
            except Exception:
                pass
            """
            #if hits_to_report == 1: print [block],row_header[i], row_header[i-1],rho,rho_cutoff 
            #print rho
            if rho>0.95:
                pass ### don't store this
            elif rho>rho_cutoff:
                try:
                    block_db[block].append(i) ### store the row index
                except Exception:
                    block_db[block] = [i] ### store the row index
            else:
                block+=1
                block_db[block] = [i] ### store the row index
        else:
            block_db[block] = [i] ### store the row index
        i+=1
        
    if ReDefinedClusterBlocks:
        ### Produces a filtered-down and centroid organized heatmap text file
        exportReDefinedClusterBlocks(results_file,block_db,rho_cutoff)
        
    if hits_to_report == 1:
        if len(block_db)<4 and forceOutput==False:
            return 'TooFewBlocks', None
        guideGenes={}
        ### Select the top TFs or non-TFs with the most gene correlations
        for b in block_db:
            corr_counts_gene = []; cell_cycle_count=[]
            #print len(block_db), b, map(lambda i: row_header[i],block_db[b])
            for (gene,i) in map(lambda i: (row_header[i],i),block_db[b]):
                corr_counts_gene.append((len(gene_corr_counts[gene][1]),gene_corr_counts[gene][0],gene))
                if gene in cell_cycle:
                    cell_cycle_count.append(gene)
            corr_counts_gene.sort(); tfs=[]
            #print b, corr_counts_gene, '***',len(cell_cycle_count)
            if (len(cell_cycle_count)>1) or (len(corr_counts_gene)<4 and (len(cell_cycle_count)>0)): pass
            else:
                tf_count=0
                for (r,t, gene) in corr_counts_gene:
                    if gene in TFs:
                        if gene not in cell_cycle:
                            if restrictTFs==True and tf_count==0: pass
                            else:
                                guideGenes[gene]=[]
                                tf_count+=1
                if len(tfs)==0:
                    gene = corr_counts_gene[-1][-1]
                    if gene in cell_cycle and LegacyMode: pass
                    else:
                        guideGenes[gene]=[]
                #block_db[b]= [corr_counts_gene[-1][-1]] ### save just the selected gene indexes
        
        ### Additional filter to remove guides that will bring in cell cycle genes (the more guides the more likely)
        if excludeCellCycle == 'strict':
            #print 'guides',len(guideGenes)
            guideCorrelated = numpyCorrelationMatrixGeneAlt(matrix,row_header,guideGenes,gene_to_symbol_db,rho_cutoff)
            guideGenes={}
            for gene in guideCorrelated:
                cell_cycle_count=[]
                for corr_gene in guideCorrelated[gene]:
                    if corr_gene in cell_cycle: cell_cycle_count.append(corr_gene)
                #print gene, len(cell_cycle_count),len(guideCorrelated[gene])
                if (float(len(cell_cycle_count))/len(guideCorrelated[gene]))>.15 or (len(guideCorrelated[gene])<4 and (len(cell_cycle_count)>0)):
                    print gene, cell_cycle_count
                    addition_cell_cycle_associated.append(gene)
                    pass
                else:
                    guideGenes[gene]=[]
            print 'additional Cell Cycle guide genes removed:',addition_cell_cycle_associated
        print len(guideGenes), 'novel guide genes discovered:', guideGenes.keys()
        return guideGenes,addition_cell_cycle_associated
        
    def greaterThan(x,results_file,numSamplesClustered):
        if 'alt_junctions' not in results_file and Platform == None:
            if x>(numSamplesClustered-1): return 1
            else: return 0
        else:
            return 1
    
    max_block_size=0
    ### Sometimes the hits_cutoff is too stringent so take the largest size instead
    for block in block_db:
        indexes = len(block_db[block])
        if indexes>max_block_size: max_block_size=indexes
        
    max_block_size-=1
    retained_ids={}; final_rows = {}
    for block in block_db:
        indexes = block_db[block]
        #print [block], len(indexes),hits_cutoff,max_block_size
        if len(indexes)>hits_cutoff or len(indexes)>max_block_size: ###Increasing this helps get rid of homogenous clusters of little significance
            #if statistics.avg(matrix[indexes[0]][1:]) < -2: print statistics.avg(matrix[indexes[0]][1:]), len(indexes)
            gene_names = map(lambda i: row_header[i], indexes)
            #if 'Pax6' in gene_names or 'WNT8A' in gene_names: print '******',hits_to_report, gene_names
            indexes = indexes[:hits_to_report]
            if filter:
                new_indexes = []
                for index in indexes:
                    vs = list(matrix[index])
                    a = map(lambda x: greaterThan(x,results_file,numSamplesClustered),vs)
                    b=[1]*numSamplesClustered
                    c = [(i, i+len(b)) for i in range(len(a)) if a[i:i+len(b)] == b]
                    if len(c)>0: #http://stackoverflow.com/questions/10459493/find-indexes-of-sequence-in-list-in-python
                        new_indexes.append(index)
                        """
                        vs.sort()
                        try:
                            if abs(vs[-5]-vs[5])>6: new_indexes.append(index)
                        except Exception:
                            if abs(vs[-1]-vs[1])>6: new_indexes.append(index)"""
                    indexes = new_indexes
            #if block == 1: print map(lambda i:row_header[i],indexes)
                #print indexes;sys.exit()
            for ls in map(lambda i: [row_header[i]]+map(str,(matrix[i])), indexes):
                final_rows[tuple(ls)]=[]
            for i in indexes:
                retained_ids[row_header[i]]=[]
                
    if len(final_rows)==0:
        for block in block_db:
            indexes = block_db[block]
            if len(indexes)>hits_cutoff or len(indexes)>max_block_size:
                indexes = indexes[:hits_to_report]
            for ls in map(lambda i: [row_header[i]]+map(str,(matrix[i])), indexes):
                final_rows[tuple(ls)]=[]
    if len(final_rows)==0:
        for block in block_db:
            indexes = block_db[block]
            for ls in map(lambda i: [row_header[i]]+map(str,(matrix[i])), indexes):
                final_rows[tuple(ls)]=[]
    #print 'block length:',len(block_db), 'genes retained:',len(retained_ids)

    return final_rows, column_header

def exportGroupsFromClusters(cluster_file,expFile,platform,suffix=None):
    lineNum=1
    for line in open(cluster_file,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if lineNum==1: names = t[2:]; lineNum+=1
        elif lineNum==2: clusters = t[2:]; lineNum+=1
        else: break

    unique_clusters=[] ### Export groups
    
    new_groups_dir = string.replace(expFile,'exp.','groups.')
    new_comps_dir = string.replace(expFile,'exp.','comps.')
    if suffix != None:
        new_groups_dir = new_groups_dir[:-4]+'-'+suffix+'.txt' ###Usually end in ICGS
        new_comps_dir = new_comps_dir[:-4]+'-'+suffix+'.txt'
    out_obj = export.ExportFile(new_groups_dir)
    cluster_number=0
    cluster_db={}
    for name in names:
        cluster = clusters[names.index(name)]
        if platform == 'RNASeq':
            if 'junction_quantification' not in name and '.bed' not in name:
                name = name+'.bed'
            elif 'junction_quantification.txt' not in name and '.txt' not in name and '.bed' not in name:
                name = name+'.txt'
        if ':' in name:
            group,name = string.split(name,':')
            if group in cluster_db:
                clust_num=cluster_db[group]
            else:
                cluster_number+=1
                cluster_db[group] = cluster_number
                clust_num = cluster_number
            if cluster=='NA': cluster = group
        else:
            clust_num = cluster
        out_obj.write(name+'\t'+str(clust_num)+'\t'+cluster+'\n')
        if cluster not in unique_clusters: unique_clusters.append(cluster)
    out_obj.close()
    comps=[] #Export comps
    out_obj = export.ExportFile(new_comps_dir)
    """ ### All possible pairwise group comparisons
    for c1 in unique_clusters:
        for c2 in unique_clusters:
            temp=[int(c2),int(c1)]; temp.sort(); temp.reverse()
            if c1 != c2 and temp not in comps:
                out_obj.write(str(temp[0])+'\t'+str(temp[1])+'\n')
                comps.append(temp)
    """
    ### Simple method comparing each subsequent ordered cluster (HOPACH orders based on relative similarity)
    last_cluster = None
    for c1 in unique_clusters:
        if last_cluster !=None:
            out_obj.write(c1+'\t'+last_cluster+'\n')
        last_cluster=c1
    out_obj.close()
    return new_groups_dir

def logTransform(value):
    try: v = math.log(value,2)
    except Exception: v = math.log(0.001,2)
    return str(v)

class MultiCorrelatePatterns():
    def __init__(self,expressed_values):
        self.expressed_values = expressed_values        
    def __call__(self,features_to_correlate):
        from  scipy import stats
        correlated_genes={}
        for uid in features_to_correlate:
            ref_values = self.expressed_values[uid]
            for uid2 in self.expressed_values:
                values = self.expressed_values[uid2]
                rho,p = stats.pearsonr(values,ref_values)
                if rho>rho_cutoff or rho<-1*rho_cutoff:
                    if uid!=uid2 and rho != 1.0:
                        try: correlated_genes[uid].append(uid2)
                        except Exception: correlated_genes[uid] = [uid]
        return correlated_genes

def parseCountFile(fn,parseFeature,search_exon_db):
    novel_exon_db={}; firstLine=True
    unique_genes={}
    for line in open(fn,'rU').xreadlines():
        key = string.split(line,'\t')[0]
        #t = string.split(line,'\t')
        if firstLine: firstLine = False
        else:
            #uid, coordinates = string.split(key,'=')
            #values = map(lambda x: float(x), t[1:])
            #gene = string.split(uid,':')[0]
            #if max(values)>5: unique_genes[gene] = []
            if '_' in key: ### Only look at novel exons
                #ENSG00000112695:I2.1_75953139=chr6:75953139-75953254
                uid, coordinates = string.split(key,'=')
                gene = string.split(uid,':')[0]
                if parseFeature == 'exons':
                    if '-' not in uid:
                        chr,coordinates = string.split(coordinates,':') ### Exclude the chromosome
                        coord1,coord2 = string.split(coordinates,'-')
                        intron = string.split(uid,'_')[0]
                        intron = string.split(intron,':')[1]
                        first = intron+'_'+coord1
                        second = intron+'_'+coord2
                        proceed = True
                        if first in uid: search_uid = second ### if the first ID is already the one looked for, store the second with the exon ID
                        elif second in uid: search_uid = first
                        else:
                            proceed = False
                            #print uid, first, second; sys.exit()
                            #example: ENSG00000160785:E2.15_156170151;E2.16_156170178=chr1:156170151-156170178
                        if proceed:
                            try: novel_exon_db[gene].append((uid,search_uid))
                            except Exception: novel_exon_db[gene] = [(uid,search_uid)]  
                elif '-' in uid and 'I' in uid: ### get junctions
                    if gene in search_exon_db:
                        for (u,search_uid) in search_exon_db[gene]:
                            #if gene == 'ENSG00000137076': print u,search_uid,uid
                            if search_uid in uid:
                                novel_exon_db[uid] = u ### Relate the currently examined novel exon ID to the junction not current associated
                                #if gene == 'ENSG00000137076': print u, uid
                                #print uid;sys.exit()
                                
    #print len(unique_genes); sys.exit()
    return novel_exon_db

def getJunctionType(species,fn):
    root_dir = string.split(fn,'ExpressionInput')[0]
    fn = filepath(root_dir+'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_junctions.txt')
    firstLine=True
    junction_type_db={}; type_db={}
    for line in open(fn,'rU').xreadlines():
        t = string.split(line,'\t')
        if firstLine: firstLine = False
        else:
            id=t[0]; junction_type = t[8]
            if '-' in id:
                if 'trans-splicing' in line:
                    junction_type = 'trans-splicing'
                junction_type_db[id] = junction_type
                try: type_db[junction_type]+=1
                except Exception: type_db[junction_type]=1
    print 'Breakdown of event types'
    for type in type_db:
        print type, type_db[type]
    return junction_type_db
            
def maxCount(ls):
    c=0
    for i in ls:
        if i>0.5: c+=1
    return c
    

def getHighExpNovelExons(species,fn):
    """ Idea - if the ranking of exons based on expression changes from one condition to another, alternative splicing is occuring """
    junction_type_db = getJunctionType(species,fn)
    
    ### Possible issue detected with novel exon reads: ['ENSG00000121577'] ['119364543'] cardiac
    exon_max_exp_db={}; uid_key_db={}; firstLine=True
    novel_intronic_junctions = {}
    novel_intronic_exons = {}
    cutoff = 0.2
    read_threshold = 0.5
    expressed_junction_types={}
    features_to_export={}
    exon_coord_db={}
    for line in open(fn,'rU').xreadlines():
        t = string.split(line,'\t')
        if firstLine: firstLine = False
        else:
            key=t[0]
            #ENSG00000112695:I2.1_75953139=chr6:75953139-75953254
            try: uid, coordinates = string.split(key,'=')
            except Exception: uid = key
            gene = string.split(uid,':')[0]
            values = map(lambda x: float(x), t[1:])
            max_read_counts = max(values)
            try: exon_max_exp_db[gene].append((max_read_counts,uid))
            except Exception: exon_max_exp_db[gene] = [(max_read_counts,uid)]
            uid_key_db[uid] = key ### retain the coordinate info
            if '-' in uid and (':E' in uid or '-E' in uid):
                junction_type = junction_type_db[uid]
                if max_read_counts>read_threshold:
                    samples_expressed = maxCount(values)
                    if samples_expressed>2:
                        try: expressed_junction_types[junction_type]+=1
                        except Exception: expressed_junction_types[junction_type]=1
                        if junction_type == 'trans-splicing' and'_' not in uid:
                            try: expressed_junction_types['known transplicing']+=1
                            except Exception: expressed_junction_types['known transplicing']=1
                        elif junction_type == 'novel' and '_' not in uid:
                            try: expressed_junction_types['novel but known sites']+=1
                            except Exception: expressed_junction_types['novel but known sites']=1
                        elif junction_type == 'novel' and 'I' not in uid:
                            try: expressed_junction_types['novel but within 50nt of a known sites']+=1
                            except Exception: expressed_junction_types['novel but within 50nt of a known sites']=1
                        elif 'I' in uid and '_' in uid and junction_type!='trans-splicing':
                            #print uid;sys.exit()
                            try: expressed_junction_types['novel intronic junctions']+=1
                            except Exception: expressed_junction_types['novel intronic junctions']=1
                            coord = string.split(uid,'_')[-1]
                            if '-' in coord:
                                coord = string.split(coord,'-')[0]
                            try: novel_intronic_junctions[gene]=[coord]
                            except Exception: novel_intronic_junctions[gene].append(coord)
            elif ('I' in uid or 'U' in uid) and '_' in uid and max_read_counts>read_threshold:
                if '-' not in uid:
                    samples_expressed = maxCount(values)
                    if samples_expressed>2:
                        try: expressed_junction_types['novel intronic exon']+=1
                        except Exception: expressed_junction_types['novel intronic exon']=1
                        coord = string.split(uid,'_')[-1]
                        #print uid, coord;sys.exit()
                        #if 'ENSG00000269897' in uid: print [gene,coord]
                        try: novel_intronic_exons[gene].append(coord)
                        except Exception: novel_intronic_exons[gene]=[coord]
                        exon_coord_db[gene,coord]=uid
                
    print 'Expressed (count>%s for at least 3 samples) junctions' % read_threshold
    for junction_type in expressed_junction_types:
        print junction_type, expressed_junction_types[junction_type]
    expressed_junction_types={}
    #print len(novel_intronic_junctions)
    #print len(novel_intronic_exons)
    
    for gene in novel_intronic_junctions:
        if gene in novel_intronic_exons:
            for coord in novel_intronic_junctions[gene]:
                if coord in novel_intronic_exons[gene]:
                    try: expressed_junction_types['confirmed novel intronic exons']+=1
                    except Exception: expressed_junction_types['confirmed novel intronic exons']=1
                    uid = exon_coord_db[gene,coord]
                    features_to_export[uid]=[]
        #else: print [gene], novel_intronic_junctions[gene]; sys.exit()
                    
    for junction_type in expressed_junction_types:
        print junction_type, expressed_junction_types[junction_type]
    
    
    out_file = string.replace(fn,'.txt','-highExp.txt')
    print 'Exporting the highest expressed exons to:', out_file
    out_obj = export.ExportFile(out_file)
    ### Compare the relative expression of junctions and exons separately for each gene (junctions are more comparable)
    for gene in exon_max_exp_db:
        junction_set=[]; exon_set=[]; junction_exp=[]; exon_exp=[]
        exon_max_exp_db[gene].sort()
        exon_max_exp_db[gene].reverse()
        for (exp,uid) in exon_max_exp_db[gene]:
            if '-' in uid: junction_set.append((exp,uid)); junction_exp.append(exp)
            else: exon_set.append((exp,uid)); exon_exp.append(exp)
        if len(junction_set)>0:
            maxJunctionExp = junction_set[0][0]
            try: lower25th,median_val,upper75th,int_qrt_range = statistics.iqr(junction_exp)
            except Exception: print junction_exp;sys.exit()
            if int_qrt_range>0:
                maxJunctionExp = int_qrt_range
            junction_percent_exp = map(lambda x: (x[1],expThreshold(x[0]/maxJunctionExp,cutoff)), junction_set)
            high_exp_junctions = []
            for (uid,p) in junction_percent_exp: ### ID and percentage of expression
                if p!='NA':
                    if uid in features_to_export: ### novel exons only right now
                        out_obj.write(uid_key_db[uid]+'\t'+p+'\n') ### write out the original ID with coordinates
            
        if len(exon_set)>0:
            maxExonExp = exon_set[0][0]
            lower25th,median_val,upper75th,int_qrt_range = statistics.iqr(exon_exp)
            if int_qrt_range>0:
                maxExonExp = int_qrt_range
            exon_percent_exp = map(lambda x: (x[1],expThreshold(x[0]/maxExonExp,cutoff)), exon_set)
            high_exp_exons = []
            for (uid,p) in exon_percent_exp: ### ID and percentage of expression
                if p!='NA':
                    if uid in features_to_export:
                        out_obj.write(uid_key_db[uid]+'\t'+p+'\n') 
    out_obj.close()

def expThreshold(ratio,cutoff):
    #print [ratio,cutoff]
    if ratio>cutoff: return str(ratio)
    else: return 'NA'

def compareExonAndJunctionResults(species,array_type,summary_results_db,root_dir):
    results_dir = root_dir +'AltResults/AlternativeOutput/'
    dir_list = read_directory(results_dir)
    filtered_dir_db={}
    #"""
    try: novel_exon_junction_db = getNovelExonCoordinates(species,root_dir)
    except Exception:
        #print traceback.format_exc()
        print 'No counts file found.'
        novel_exon_junction_db={} ### only relevant to RNA-Seq analyses

    for comparison_file in summary_results_db:
        for results_file in dir_list:
            if (comparison_file in results_file and '-exon-inclusion-results.txt' in results_file) and ('comparison' not in results_file):
                try: filtered_dir_db[comparison_file].append(results_file)
                except Exception: filtered_dir_db[comparison_file] = [results_file]
                
    try: os.remove(string.split(results_dir,'AltResults')[0]+'AltResults/Clustering/Combined-junction-exon-evidence.txt')
    except Exception: pass
    for comparison_file in filtered_dir_db:
        alt_result_files = filtered_dir_db[comparison_file]
        #print alt_result_files, comparison_file
        importAltAnalyzeExonResults(alt_result_files,novel_exon_junction_db,results_dir)
    #"""
    ### Build combined clusters of high-confidence exons
    graphics2=[]; graphics=[]
    import ExpressionBuilder
    try:
        input_dir = string.split(results_dir,'AltResults')[0]+'GO-Elite/AltExonConfirmed/'
        cluster_file, rows_in_file = ExpressionBuilder.buildAltExonClusterInputs(input_dir,species,array_type,dataType='AltExonConfirmed')
        if rows_in_file > 5000: useHOPACH = False
        else: useHOPACH = True
        if rows_in_file < 12000:
            graphics = ExpressionBuilder.exportHeatmap(cluster_file,useHOPACH=useHOPACH)
    except Exception: pass
    
    try:
        input_dir = string.split(results_dir,'AltResults')[0]+'GO-Elite/AltExon/'
        cluster_file, rows_in_file = ExpressionBuilder.buildAltExonClusterInputs(input_dir,species,array_type,dataType='AltExon')
        if rows_in_file > 5000: useHOPACH = False
        else: useHOPACH = True
        if rows_in_file < 12000:
            graphics2 = ExpressionBuilder.exportHeatmap(cluster_file,useHOPACH=useHOPACH)
    except Exception: pass
    return graphics+graphics2

class SplicingData:
    def __init__(self,score,symbol,description,exonid,probesets,direction,splicing_event,external_exon,genomic_loc,gene_exp,protein_annot,domain_inferred,domain_overlap,method,dataset):
        self.score = score; self.dataset = dataset
        self.symbol = symbol;
        self.description=description;self.exonid=exonid;self.probesets=probesets;self.direction=direction
        self.splicing_event=splicing_event;self.external_exon=external_exon;self.genomic_loc=genomic_loc;
        self.gene_exp=gene_exp;self.protein_annot=protein_annot;self.domain_inferred=domain_inferred
        self.domain_overlap=domain_overlap;self.method=method
    def Score(self): return self.score
    def setScore(self,score): self.score = score
    def GeneExpression(self): return self.gene_exp
    def Dataset(self): return self.dataset
    def Symbol(self): return self.symbol
    def Description(self): return self.description
    def ExonID(self): return self.exonid
    def appendExonID(self,exonid): self.exonid+='|'+exonid
    def Probesets(self): return self.probesets
    def ProbesetDisplay(self):
        if len(self.Probesets()[1])>0:
            return string.join(self.Probesets(),'-')
        else:
            return self.Probesets()[0]
    def ProbesetsSorted(self):
        ### Don't sort the original list
        a = [self.probesets[0],self.probesets[1]]
        a.sort()
        return a
    def Direction(self): return self.direction
    def setDirection(self,direction): self.direction = direction
    def SplicingEvent(self): return self.splicing_event
    def ProteinAnnotation(self): return self.protein_annot
    def DomainInferred(self): return self.domain_inferred
    def DomainOverlap(self): return self.domain_overlap
    def Method(self): return self.method
    def setEvidence(self,evidence): self.evidence = evidence
    def Evidence(self): return self.evidence
    def GenomicLocation(self): return self.genomic_loc
    def setExonExpStatus(self, exon_expressed): self.exon_expressed = exon_expressed
    def ExonExpStatus(self): return self.exon_expressed
    
def importAltAnalyzeExonResults(dir_list,novel_exon_junction_db,results_dir):
    
    regulated_critical_exons={}; converted_db={}
    includeExonJunctionComps=True ### Allow ASPIRE comparisons with the inclusion feature as an exon to count for additive reciprocal evidence
    print "Reading AltAnalyze results file"
    root_dir = string.split(results_dir,'AltResults')[0]
    
    for filename in dir_list:
        x=0; regulated_critical_exon_temp={}
        fn=filepath(results_dir+filename)
        new_filename = string.join(string.split(filename,'-')[:-5],'-')
        if '_vs_' in filename and '_vs_' in new_filename: export_filename = new_filename
        else: export_filename = string.join(string.split(filename,'-')[:-5],'-')

        export_path = results_dir+export_filename+'-comparison-evidence.txt'
        try: os.remove(filepath(export_path)) ### If we don't do this, the old results get added to the new
        except Exception: null=[]
    
        if 'AltMouse' in filename:
            altmouse_ensembl_db = importAltMouseEnsembl()
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0: x=1; #print t[12],t[13],t[22],t[23]
            else:
                converted = False ### Indicates both junction sides were regulated
                geneid = t[0]; exonid = t[4]; probeset1 = t[6]; probeset2 = ''; score = t[1][:4]; symbol = t[2]; description = t[3]; regions = t[-4]; direction = t[5]
                genomic_loc = t[-1]; splicing_event = t[-3]; external_exon = t[-6]; gene_exp_fold = t[-8]; protein_annot = t[14]; domain_inferred = t[15]; domain_overlap = t[17]
                expressed_exon = 'NA'
                if 'RNASeq' in filename: expressed_exon = 'no' ### Set by default
                if ':' in geneid: geneid = string.split(geneid,':')[0] ### User reported that gene:gene was appearing and not sure exactly where or why but added this to address it
                if 'FIRMA' in fn: method = 'FIRMA'
                elif 'splicing-index' in fn: method = 'splicing-index'
                if 'ASPIRE' in filename or 'linearregres' in filename:
                    f1=float(t[12]); f2=float(t[13]); probeset1 = t[8]; probeset2 = t[10]; direction = t[6]; exonid2 = t[5]; splicing_event = t[-4]
                    protein_annot = t[19]; domain_inferred = t[20]; domain_overlap = t[24]; method = 'linearregres'; regions = t[-5]
                    exon1_exp=float(t[-15]); exon2_exp=float(t[-14]); fold1=float(t[12]); fold2=float(t[13])
                    if fold1<0: fold1 = 1 ### don't factor in negative changes
                    if fold2<0: fold2 = 1 ### don't factor in negative changes
                    """
                    if 'RNASeq' not in filename:
                        exon1_exp = math.pow(2,exon1_exp)
                        exon2_exp = math.log(2,exon2_exp)
                    
                    m1 = exon1_exp*fold1
                    m2 = exon2_exp*fold2
                    max_exp = max([m1,m2])
                    min_exp = min([m1,m2])
                    percent_exon_expression = str(min_exp/max_exp)
                    """
                    if 'ASPIRE' in filename: method = 'ASPIRE'; score = t[1][:5]
                    if '-' not in exonid and includeExonJunctionComps == False:
                        exonid=None ### Occurs when the inclusion just in an exon (possibly won't indicate confirmation so exclude)
                    else: exonid = exonid+' vs. '+exonid2
                    if 'AltMouse' in filename:
                        try: geneid = altmouse_ensembl_db[geneid]
                        except Exception: geneid = geneid
                    if 'RNASeq' not in filename and 'junction' not in filename: regions = string.replace(regions,'-','.')
                else:
                    if 'RNASeq' in filename and '-' not in exonid:
                        fold = float(t[10]); exon_exp = float(t[18]); gene_exp = float(t[19])
                        if fold < 0: fold = -1.0/fold
                        GE_fold = float(gene_exp_fold)
                        if GE_fold < 0: GE_fold = -1.0/float(gene_exp_fold)
                        exon_psi1 = abs(exon_exp)/(abs(gene_exp))
                        exon_psi2 = (abs(exon_exp)*fold)/(abs(gene_exp)*GE_fold)
                        max_incl_exon_exp = max([exon_psi1,exon_psi2])
                        #if max_incl_exon_exp>0.20: expressed_exon = 'yes'
                        expressed_exon = max_incl_exon_exp
                        #if 'I2.1_75953139' in probeset1:
                        #print [exon_exp,gene_exp,exon_exp*fold,gene_exp*GE_fold]
                        #print exon_psi1, exon_psi2;sys.exit()
                probesets = [probeset1,probeset2]
                if (method == 'splicing-index' or method == 'FIRMA') and ('-' in exonid) or exonid == None:
                    pass #exclude junction IDs
                else:
                    regions = string.replace(regions,';','|')
                    regions = string.replace(regions,'-','|')
                    regions = string.split(regions,'|')
                    for region in regions:
                        if len(region) == 0:
                            try: region = t[17]+t[18] ### For junction introns where no region ID exists
                            except Exception: null=[]
                        if ':' in region: region = string.split(region,':')[-1] ### User reported that gene:gene was appearing and not sure exactly where or why but added this to address it
                        if probeset1 in novel_exon_junction_db:
                            uid = novel_exon_junction_db[probeset1] ### convert the uid (alternative exon) to the annotated ID for the novel exon
                            converted_db[uid] = probeset1
                        else:
                            uid = geneid+':'+region
                        ss = SplicingData(score,symbol,description,exonid,probesets,direction,splicing_event,external_exon,genomic_loc,gene_exp_fold,protein_annot,domain_inferred,domain_overlap,method,filename)
                        ss.setExonExpStatus(str(expressed_exon))
                        try: regulated_critical_exon_temp[uid].append(ss)
                        except Exception: regulated_critical_exon_temp[uid] = [ss]
        
        #print filename, len(regulated_critical_exon_temp)
        for uid in regulated_critical_exon_temp:
            report=None
            if len(regulated_critical_exon_temp[uid])>1:
                ### We are only reporting one here and that's OK, since we are only reporting the top scores... won't include all inclusion junctions.
                scores=[]
                for ss in regulated_critical_exon_temp[uid]: scores.append((float(ss.Score()),ss))
                scores.sort()
                if (scores[0][0]*scores[-1][0])<0:
                    ss1 = scores[0][1]; ss2 = scores[-1][1]
                    if ss1.ProbesetsSorted() == ss2.ProbesetsSorted(): ss1.setDirection('mutual') ### same exons, hence, mutually exclusive event (or similiar)
                    else: ss1.setDirection('both') ### opposite directions in the same comparison-file, hence, conflicting data
                    report=[ss1]
                else:
                    if abs(scores[0][0])>abs(scores[-1][0]): report=[scores[0][1]]
                    else: report=[scores[-1][1]]
            else:
                report=regulated_critical_exon_temp[uid]
            ### Combine data from different analysis files
            try: regulated_critical_exons[uid]+=report
            except Exception: regulated_critical_exons[uid]=report
            """if 'ENSG00000204120' in uid:
                print uid,
                for i in regulated_critical_exon_temp[uid]:
                    print i.Probesets(),
                print ''
                """
            try: report[0].setEvidence(len(regulated_critical_exon_temp[uid])) ###set the number of exons demonstrating regulation of this exons
            except Exception: null=[]

    clearObjectsFromMemory(regulated_critical_exon_temp)
    
    export_data,status = AppendOrWrite(export_path)
    
    if status == 'not found': 
        header = string.join(['uid','source-IDs','symbol','description','exonids','independent confirmation','score','regulation direction','alternative exon annotations','associated isoforms','inferred regulated domains','overlapping domains','method','supporting evidence score','novel exon: high-confidence','percent exon expression of gene','differential gene-expression','genomic location'],'\t')+'\n'
        export_data.write(header)

    combined_export_path = string.split(results_dir,'AltResults')[0]+'AltResults/Clustering/Combined-junction-exon-evidence.txt'
    combined_export_data, status= AppendOrWrite(combined_export_path)
    if status == 'not found':
        header = string.join(['uid','source-IDs','symbol','description','exonids','independent confirmation','score','regulation direction','alternative exon annotations','associated isoforms','inferred regulated domains','overlapping domains','method','supporting evidence score','novel exon: high-confidence','percent exon expression of gene','differential gene-expression','genomic location','comparison'],'\t')+'\n'
        combined_export_data.write(header)
                    
    print len(regulated_critical_exons), 'regulated exon IDs imported.\n'
    print 'writing:',export_path; n=0
   # print [len(converted_db)]
    ### Check for alternative 3' or alternative 5' exon regions that were not matched to the right reciprocal junctions (occurs because only one of the exon regions is called alternative)
    regulated_critical_exons_copy={}
    for uid in regulated_critical_exons:
        regulated_critical_exons_copy[uid]=regulated_critical_exons[uid]
    u=0
    ### This is most applicable to RNA-Seq since the junction IDs correspond to the Exon Regions not the probeset Exon IDs
    for uid in regulated_critical_exons_copy: ### Look through the copied version since we can't delete entries while iterating through
        ls = regulated_critical_exons_copy[uid]
        u+=1
        #if u<20: print uid
        for jd in ls:
            if jd.Method() != 'splicing-index' and jd.Method() != 'FIRMA':
                try: ### Applicable to RNA-Seq
                    gene,exonsEx = string.split(jd.Probesets()[1],':') ### Exclusion probeset will have the exon not annotated as the critical exon (although it should be as well)
                    gene,exonsIn = string.split(jd.Probesets()[0],':')
                except Exception:
                    gene, ce = string.split(uid,':')
                    exonsIn, exonsEx = string.split(jd.ExonID(),'vs.')
                if gene !=None:
                    critical_exon = None
                    five_prime,three_prime = string.split(exonsEx,'-')
                    try: five_primeIn,three_primeIn = string.split(exonsIn,'-')
                    except Exception: five_primeIn = exonsIn; three_primeIn = exonsIn ### Only should occur during testing when a exon rather than junction ID is considered
                    #if gene == 'ENSG00000133083': print five_prime,three_prime, five_primeIn,three_primeIn
                    if five_primeIn == five_prime: ### Hence, the exclusion 3' exon should be added
                        critical_exon = gene+':'+three_prime
                        exonid = three_prime
                    elif three_primeIn == three_prime: ### Hence, the exclusion 3' exon should be added
                        critical_exon = gene+':'+five_prime
                        exonid = five_prime
                    else:
                        if ('5' in jd.SplicingEvent()) or ('five' in jd.SplicingEvent()):
                            critical_exon = gene+':'+five_prime
                            exonid = five_prime
                        elif ('3' in jd.SplicingEvent()) or ('three' in jd.SplicingEvent()):
                            critical_exon = gene+':'+three_prime
                            exonid = three_prime
                        elif ('alt-N-term' in jd.SplicingEvent()) or ('altPromoter' in jd.SplicingEvent()):
                            critical_exon = gene+':'+five_prime
                            exonid = five_prime
                        elif ('alt-C-term' in jd.SplicingEvent()):
                            critical_exon = gene+':'+three_prime
                            exonid = three_prime
                    #print critical_exon, uid, jd.ExonID(),jd.SplicingEvent(); sys.exit() 
                    if critical_exon != None:
                        if critical_exon in regulated_critical_exons:
                            #print uid, critical_exon; sys.exit()
                            if len(regulated_critical_exons[critical_exon]) == 1:
                                if len(ls)==1 and uid in regulated_critical_exons: ### Can be deleted by this method
                                    if 'vs.' not in regulated_critical_exons[critical_exon][0].ExonID() and 'vs.' not in regulated_critical_exons[critical_exon][0].ExonID():
                                        regulated_critical_exons[uid].append(regulated_critical_exons[critical_exon][0])
                                        del regulated_critical_exons[critical_exon]
                                elif uid in regulated_critical_exons: ###If two entries already exit
                                    ed = regulated_critical_exons[uid][1]
                                    ed2 = regulated_critical_exons[critical_exon][0]
                                    if 'vs.' not in ed.ExonID() and 'vs.' not in ed2.ExonID():
                                        if ed.Direction() != ed2.Direction(): ### should be opposite directions
                                            ed.appendExonID(exonid)
                                            ed.setEvidence(ed.Evidence()+1)
                                            ed.setScore(ed.Score()+'|'+ed2.Score())
                                            del regulated_critical_exons[critical_exon]
                                            
    firstEntry=True
    for uid in regulated_critical_exons:
        if uid in converted_db:
            converted = True
        else: converted = False
        #if 'ENSG00000133083' in uid: print [uid]
        exon_level_confirmation = 'no'
        ls = regulated_critical_exons[uid]
        jd = regulated_critical_exons[uid][0] ### We are only reporting one here and that's OK, since we are only reporting the top scores... won't include all inclusion junctions.
        if len(ls)>1:
            methods = []; scores = []; direction = []; exonids = []; probesets = []; evidence = 0; genomic_location = []
            junctionids=[]
            junction_data_found = 'no'; exon_data_found = 'no'
            for jd in ls:
                if jd.Method() == 'ASPIRE' or jd.Method() == 'linearregres': 
                    junction_data_found = 'yes'
                    methods.append(jd.Method())
                    scores.append(jd.Score())
                    direction.append(jd.Direction())
                    exonids.append(jd.ExonID())
                    junctionids.append(jd.ExonID())
                    probesets.append(jd.ProbesetDisplay())
                    evidence+=jd.Evidence()
                    genomic_location.append(jd.GenomicLocation())
                    ### Prefferentially obtain isoform annotations from the reciprocal analysis which is likely more accurate
                    isoform_annotations = [jd.ProteinAnnotation(), jd.DomainInferred(), jd.DomainOverlap()]
            for ed in ls:
                if ed.Method() == 'splicing-index' or ed.Method() == 'FIRMA':
                    exon_data_found = 'yes' ### pick one of them
                    methods.append(ed.Method())
                    scores.append(ed.Score())
                    direction.append(ed.Direction())
                    exonids.append(ed.ExonID())
                    probesets.append(ed.ProbesetDisplay())
                    evidence+=ed.Evidence()
                    genomic_location.append(ed.GenomicLocation())
                    #isoform_annotations = [ed.ProteinAnnotation(), ed.DomainInferred(), ed.DomainOverlap()]
            if junction_data_found == 'yes' and exon_data_found == 'yes':
                exon_level_confirmation = 'yes'
                for junctions in junctionids:
                    if 'vs.' in junctions:
                        j1 = string.split(junctions,' vs. ')[0] ### inclusion exon or junction
                        if '-' not in j1: ### not a junction, hence, may not be sufficient to use for confirmation (see below)
                            if 'I' in j1: ### intron feature
                                if '_' in j1: ### novel predicted exon
                                    exon_level_confirmation = 'no'
                                else:
                                    exon_level_confirmation = 'yes'
                            else:
                                if '_' in j1:
                                    exon_level_confirmation = 'no'
                                else:
                                    exon_level_confirmation = 'partial'
            method = string.join(methods,'|')
            unique_direction = unique.unique(direction)
            genomic_location = unique.unique(genomic_location)
            if len(unique_direction) == 1: direction = unique_direction[0]
            else: direction = string.join(direction,'|')
            score = string.join(scores,'|')
            probesets = string.join(probesets,'|')
            exonids_unique = unique.unique(exonids)
            if len(exonids_unique) == 1: exonids = exonids_unique[0]
            else: exonids = string.join(exonids,'|')
            if len(genomic_location) == 1: genomic_location = genomic_location[0]
            else: genomic_location = string.join(genomic_location,'|')
            evidence = str(evidence)
            if 'mutual' in direction: direction = 'mutual'
        if len(ls) == 1:
            probesets = jd.ProbesetDisplay()
            direction = jd.Direction()
            score = jd.Score()
            method = jd.Method()
            exonids = jd.ExonID()
            evidence = jd.Evidence()
            genomic_location = jd.GenomicLocation()
            isoform_annotations = [jd.ProteinAnnotation(), jd.DomainInferred(), jd.DomainOverlap()]
        try:
            #if int(evidence)>4 and 'I' in uid: novel_exon = 'yes' ### high-evidence novel exon
            #else: novel_exon = 'no'
            if converted == True:
                novel_exon = 'yes'
                splicing_event = 'cassette-exon'
            else:
                novel_exon = 'no'
                splicing_event = jd.SplicingEvent()
            values = [uid, probesets, jd.Symbol(), jd.Description(), exonids, exon_level_confirmation, score, direction, splicing_event]
            values += isoform_annotations+[method, str(evidence),novel_exon,jd.ExonExpStatus(),jd.GeneExpression(),genomic_location]
            values = string.join(values,'\t')+'\n'
            #if 'yes' in exon_level_confirmation:
            export_data.write(values); n+=1
            if exon_level_confirmation != 'no' and ('|' not in direction):
                geneID = string.split(uid,':')[0]
                try: relative_exon_exp = float(jd.ExonExpStatus())
                except Exception: relative_exon_exp = 1
                if firstEntry:
                    ### Also export high-confidence predictions for GO-Elite
                    elite_export_path = string.split(results_dir,'AltResults')[0]+'GO-Elite/AltExonConfirmed/'+export_filename+'-junction-exon-evidence.txt'
                    elite_export_data = export.ExportFile(elite_export_path)
                    elite_export_data.write('GeneID\tEn\tExonID\tScores\tGenomicLocation\n')
                    firstEntry = False
                if relative_exon_exp>0.10:
                    elite_export_data.write(string.join([geneID,'En',uid,score,genomic_location],'\t')+'\n')
                #if 'DNA' in isoform_annotations[-1]:
                if '2moter' not in jd.SplicingEvent() and '2lt-N' not in jd.SplicingEvent():
                    values = [uid, probesets, jd.Symbol(), jd.Description(), exonids, exon_level_confirmation, score, direction, splicing_event]
                    values += isoform_annotations+[method, str(evidence),novel_exon,jd.ExonExpStatus(),jd.GeneExpression(),genomic_location,export_filename]
                    values = string.join(values,'\t')+'\n'
                    combined_export_data.write(values)
        except Exception, e:
            #print traceback.format_exc();sys.exit()
            pass ### Unknown error - not evaluated in 2.0.8  - isoform_annotations not referenced
    print n,'exon IDs written to file.'
    export_data.close()
    try: elite_export_data.close()
    except Exception: pass
    clearObjectsFromMemory(regulated_critical_exons)
    clearObjectsFromMemory(regulated_critical_exons_copy)
    #print '!!!!Within comparison evidence'
    #returnLargeGlobalVars()
        
def FeatureCounts(bed_ref, bam_file):
    output = bam_file[:-4]+'__FeatureCounts.bed'
    import subprocess
    #if '/bin' in kallisto_dir: kallisto_file = kallisto_dir +'/apt-probeset-summarize' ### if the user selects an APT directory
    kallisto_dir= 'AltDatabase/subreads/'
    if os.name == 'nt':
        featurecounts_file = kallisto_dir + 'PC/featureCounts.exe'; plat = 'Windows'
    elif 'darwin' in sys.platform:
        featurecounts_file = kallisto_dir + 'Mac/featureCounts'; plat = 'MacOSX'
    elif 'linux' in sys.platform:
        featurecounts_file = kallisto_dir + '/Linux/featureCounts'; plat = 'linux'
    print 'Using',featurecounts_file
    featurecounts_file = filepath(featurecounts_file)
    featurecounts_root = string.split(featurecounts_file,'bin/featureCounts')[0]
    featurecounts_file = filepath(featurecounts_file)
    print [featurecounts_file,"-a", "-F", "SAF",bed_ref, "-o", output, bam_file]
    retcode = subprocess.call([featurecounts_file,"-a",bed_ref, "-F", "SAF", "-o", output, bam_file])
    
def filterFASTAFiles(fasta_files):
    filter_fasta_files=[]
    filter_dir = export.findParentDir(fasta_files[0])+'/filtered_fasta'
    try: os.mkdir(filter_dir)
    except Exception: pass
    
    for file in fasta_files:
        import shutil
        if 'filtered.fa' in file:
            filter_fasta_files.append(file)
        else:
            filtered_fasta = file[:-3]+'-filtered.fa'
            filter_fasta_files.append(filtered_fasta)
            filename = export.findFilename(file)
            eo=export.ExportFile(filtered_fasta)
            for line in open(file,'rU').xreadlines():
                if '>'==line[0]:
                    skip=False
                    ### Exclude non-standard chromosomal transcripts
                    if 'PATCH' in line or '_1_' in line or '_1:' in line or ':HSCHR' in line or 'putative' in line or 'supercontig' in line or 'NOVEL_TEST' in line:
                        skip=True
                    else:
                        eo.write(line)
                elif skip==False:
                    eo.write(line)
            eo.close()
            shutil.move(file,filter_dir+'/'+filename)
    return filter_fasta_files

def getCoordinateFile(species):
    geneCoordFile = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    geneCoordFile = unique.filepath(geneCoordFile)
    status = verifyFile(geneCoordFile)
    if status == 'not found':
        try:
            from build_scripts import EnsemblSQL
            ensembl_version = string.replace(unique.getCurrentGeneDatabaseVersion(),'EnsMart','')
            configType = 'Advanced'; analysisType = 'AltAnalyzeDBs'; externalDBName = ''; force = 'no'
            EnsemblSQL.buildEnsemblRelationalTablesFromSQL(species,configType,analysisType,externalDBName,ensembl_version,force,buildCommand='exon')     
        except Exception:
            #print traceback.format_exc()
            print 'Failed to export a transcript-exon coordinate file (similar to a GTF)!!!!\n...Proceeding with standard Kallisto (no-splicing).'
            geneCoordFile=None
    return geneCoordFile

def runKallisto(species,dataset_name,root_dir,fastq_folder,mlp,returnSampleNames=False,customFASTA=None,log_output=True):

    #print 'Running Kallisto...please be patient'
    import subprocess

    n_threads = mlp.cpu_count()
    print 'Number of threads =',n_threads
    #n_threads = 1

    kallisto_dir_objects = os.listdir(unique.filepath('AltDatabase/kallisto'))
    ### Determine version
    version = '0.43.1'
    for subdir in kallisto_dir_objects:
        if subdir.count('.')>1: version = subdir
    kallisto_dir= 'AltDatabase/kallisto/'+version+'/'

    if os.name == 'nt':
        kallisto_file = kallisto_dir + 'PC/bin/kallisto.exe'; plat = 'Windows'
    elif 'darwin' in sys.platform:
        kallisto_file = kallisto_dir + 'Mac/bin/kallisto'; plat = 'MacOSX'
    elif 'linux' in sys.platform:
        kallisto_file = kallisto_dir + '/Linux/bin/kallisto'; plat = 'linux'
    print 'Using',kallisto_file
    kallisto_file = filepath(kallisto_file)
    kallisto_root = string.split(kallisto_file,'bin/kallisto')[0]
    fn = filepath(kallisto_file)
    
    output_dir=root_dir+'/ExpressionInput/kallisto/'
    try: os.mkdir(root_dir+'/ExpressionInput')
    except Exception: pass
    try: os.mkdir(root_dir+'/ExpressionInput/kallisto')
    except Exception: pass
    fastq_folder += '/'

    dir_list = read_directory(fastq_folder)
    fastq_paths = []
    for file in dir_list:
        file_lower = string.lower(file)
        if 'fastq' in file_lower and '._' not in file[:4]: ### Hidden files
            fastq_paths.append(fastq_folder+file)
    fastq_paths,paired = findPairs(fastq_paths)

    ### Check to see if Kallisto files already exist and use these if so (could be problematic but allows for outside quantification)
    kallisto_tsv_paths=[]
    dir_list = read_directory(output_dir)
    for folder in dir_list:
        kallisto_outdir = output_dir+folder+'/abundance.tsv'
        status = os.path.isfile(kallisto_outdir)
        if status:
            kallisto_tsv_paths.append(fastq_folder+file)
            
    if returnSampleNames:
        return fastq_paths
    
    ### Store/retreive the Kallisto index in the Ensembl specific SequenceData location
    kallisto_index_root = 'AltDatabase/'+species+'/SequenceData/'
    try: os.mkdir(filepath(kallisto_index_root))
    except Exception: pass
    indexFile =  filepath(kallisto_index_root+species)
    #indexFile = filepath(kallisto_index_root + 'Hs_intron')

    indexStatus = os.path.isfile(indexFile)

    if indexStatus == False or customFASTA!=None:
        try: fasta_files = getFASTAFile(species)
        except Exception: fasta_files = []
        index_file = filepath(kallisto_index_root+species)
        if len(fasta_files)==0 and customFASTA==None:
            ###download Ensembl fasta file to the above directory
            from build_scripts import EnsemblSQL
            ensembl_version = string.replace(unique.getCurrentGeneDatabaseVersion(),'EnsMart','')
            try:
                EnsemblSQL.getEnsemblTranscriptSequences(ensembl_version,species,restrictTo='cDNA')
                fasta_files = getFASTAFile(species)
            except Exception: pass
        elif customFASTA!=None:  ### Custom FASTA file supplied by the user
            fasta_files = [customFASTA]
            indexFile = filepath(kallisto_index_root+species+'-custom')
            try: os.remove(indexFile) ### erase any pre-existing custom index
            except Exception: pass
        if len(fasta_files)>0:
            print 'Building kallisto index file...'
            arguments = [kallisto_file, "index","-i", indexFile]
            fasta_files = filterFASTAFiles(fasta_files)

            for fasta_file in fasta_files:
                arguments.append(fasta_file)

            try:
                retcode = subprocess.call(arguments)
            except Exception:
                print traceback.format_exc()

    if customFASTA!=None:
        reimportExistingKallistoOutput = False
    elif len(kallisto_tsv_paths) == len(fastq_paths):
        reimportExistingKallistoOutput = True
    elif len(kallisto_tsv_paths) > len(fastq_paths):
        reimportExistingKallistoOutput = True ### If working with a directory of kallisto results
    else:
        reimportExistingKallistoOutput = False

    if reimportExistingKallistoOutput:
        print 'NOTE: Re-import PREVIOUSLY GENERATED kallisto output:',reimportExistingKallistoOutput
        print '...To force re-analysis of FASTQ files, delete the folder "kallisto" in "ExpressionInput"'
        ### Just get the existing Kallisto output folders
        fastq_paths = read_directory(output_dir)

    kallisto_folders=[]
    try:
        import collections
        expMatrix = collections.OrderedDict()
        countMatrix = collections.OrderedDict()
        countSampleMatrix = collections.OrderedDict()
        sample_total_counts = collections.OrderedDict()
    except Exception:
        try:
            import ordereddict
            expMatrix = ordereddict.OrderedDict()
            countMatrix = ordereddict.OrderedDict()
            countSampleMatrix = ordereddict.OrderedDict()
            sample_total_counts = ordereddict.OrderedDict()
        except Exception:
            expMatrix={}
            countMatrix={}
            countSampleMatrix={}
            sample_total_counts={}
    headers=['UID']

    ### Verify, import, create and/or ignore the transcript exon coordinate file for BAM file creation
    geneCoordFile = getCoordinateFile(species)
    for n in fastq_paths:
        output_path = output_dir+n
        kallisto_folders.append(output_path)
        if reimportExistingKallistoOutput == False:
            begin_time = time.time()
            if geneCoordFile != None: ### For BAM and BED file generation
                print 'Running kallisto on:',n,'...',
                p=fastq_paths[n]
                b=[" > "+n+'.sam']
                bedFile = root_dir+ '/' + n + '__junction.bed'
                kallisto_out = open(root_dir+ '/' + n + '.bam', 'ab')
                if log_output:
                    err_out = open(output_dir + '/log.txt', 'a')
                    err_out.seek(0, 2)  # Subprocess doesn't move the file pointer when appending!
                else:
                    err_out = None
                kallisto_out.seek(0, 2)  # Subprocess doesn't move the file pointer when appending!

            if paired == 'paired':
                s=[]
            else:
                s=["--single","-l","200","-s","20"]
            
            #geneCoordFile=None - force to run simple Kallisto
            if geneCoordFile==None:
                    try: ### Without BAM and BED file generation
                        retcode = subprocess.call([kallisto_file, "quant","-i", indexFile, "-o", output_path]+s+p)
                    except Exception:
                        print traceback.format_exc()
            else: ### Attempt to export BAM and BED files with Kallisto quantification
                    kallisto_command = [kallisto_file, "quant", "-i", indexFile, "-o", output_path,
                                        "-g", geneCoordFile, "-j", bedFile, "--threads="+str(n_threads), "--sortedbam"] + s +p
                    kallisto_process = subprocess.Popen(kallisto_command, stdout=kallisto_out, stderr=err_out)
                    kallisto_process.communicate()
                    retcode = kallisto_process.returncode
                    if os.name == 'nt':
                        try:
                            sam_process = subprocess.Popen('AltDatabase\samtools\samtools.exe index ' + root_dir+ '/' + n + '.bam')
                            sam_process.communicate()
                            retcode_sam = sam_process.returncode
                        except: pass
                    #retcode = subprocess.call([kallisto_file, "quant","-i", indexFile, "-o", output_path,"--pseudobam"]+p+b)
                    #retcode = subprocess.call([kallisto_file, "quant","-i", indexFile, "-o", output_path]+p)
                    """except Exception:
                    print traceback.format_exc()
                    kill
                    retcode = subprocess.call(['kallisto', "quant","-i", indexFile, "-o", output_path]+p)"""
            if retcode == 0: print 'completed in', int(time.time()-begin_time), 'seconds'
            else: print 'kallisto failed due to an unknown error (report to altanalyze.org help).'

            #"""
        input_path = output_path+'/abundance.txt'
        try:
            try: expMatrix,countMatrix,countSampleMatrix=importTPMs(n,input_path,expMatrix,countMatrix,countSampleMatrix)
            except Exception:
                input_path = output_path+'/abundance.tsv'
                expMatrix,countMatrix,countSampleMatrix=importTPMs(n,input_path,expMatrix,countMatrix,countSampleMatrix)
            headers.append(n)
            sample_total_counts = importTotalReadCounts(n,output_path+'/run_info.json',sample_total_counts)
        except Exception:
            print traceback.format_exc()
            sys.exit()
            print n, 'TPM expression import failed'
            if paired == 'paired':
                print '\n...Make sure the paired-end samples were correctly assigned:'
                print fastq_paths
                for i in fastq_paths:
                    print 'Common name:',i,
                    for x in fastq_paths[i]:
                        print export.findParentDir(x),
                    print '\n'
                    
    ### Summarize alignment information
    for sample in countSampleMatrix:
        try: estCounts = int(float(countSampleMatrix[sample]))
        except Exception: estCounts='NA'
        try: totalCounts = sample_total_counts[sample]
        except Exception: totalCounts = 'NA'
        try: aligned = str(100*estCounts/float(totalCounts))
        except Exception: aligned = 'NA'
        try: aligned =  string.split(aligned,'.')[0]+'.'+string.split(aligned,'.')[1][:2]
        except Exception: aligned = 'NA'
        countSampleMatrix[sample] = [str(estCounts),totalCounts,aligned]
    
    dataset_name = string.replace(dataset_name,'exp.','')
    dataset_name = string.replace(dataset_name,'.txt','')
    to = export.ExportFile(root_dir+'/ExpressionInput/transcript.'+dataset_name+'.txt')
    ico = export.ExportFile(root_dir+'/ExpressionInput/isoCounts.'+dataset_name+'.txt')
    go = export.ExportFile(root_dir+'/ExpressionInput/exp.'+dataset_name+'.txt')
    co = export.ExportFile(root_dir+'/ExpressionInput/counts.'+dataset_name+'.txt')
    so = export.ExportFile(root_dir+'/ExpressionInput/summary.'+dataset_name+'.txt')
    exportMatrix(to,headers,expMatrix) ### Export transcript expression matrix
    exportMatrix(ico,headers,countMatrix,counts=True) ### Export transcript count matrix
    
    try:
        geneMatrix = calculateGeneTPMs(species,expMatrix) ### calculate combined gene level TPMs
        countsGeneMatrix = calculateGeneTPMs(species,countMatrix) ### calculate combined gene level TPMs
        exportMatrix(go,headers,geneMatrix) ### export gene expression matrix
        exportMatrix(co,headers,countsGeneMatrix,counts=True) ### export gene expression matrix
    except Exception: 
        print 'AltAnalyze was unable to summarize gene TPMs from transcripts, proceeding with transcripts.'
        export.copyFile(root_dir+'/ExpressionInput/transcript.'+dataset_name+'.txt',root_dir+'/ExpressionInput/exp.'+dataset_name+'.txt')
    exportMatrix(so,['SampleID','Estimated Counts','Total Fragments','Percent Aligned'],countSampleMatrix) ### export gene expression matrix
    
    ### Copy results to the Kallisto_Results directory
    try: os.mkdir(root_dir+'/ExpressionInput/Kallisto_Results')
    except: pass
    import shutil
    try:
        tf = root_dir+'/ExpressionInput/transcript.'+dataset_name+'.txt'
        shutil.copyfile(tf,string.replace(tf,'ExpressionInput','ExpressionInput/Kallisto_Results'))
        tf = root_dir+'/ExpressionInput/isoCounts.'+dataset_name+'.txt'
        shutil.copyfile(tf,string.replace(tf,'ExpressionInput','ExpressionInput/Kallisto_Results'))
        tf = root_dir+'/ExpressionInput/exp.'+dataset_name+'.txt'
        shutil.copyfile(tf,string.replace(tf,'ExpressionInput','ExpressionInput/Kallisto_Results'))
        tf = root_dir+'/ExpressionInput/counts.'+dataset_name+'.txt'
        shutil.copyfile(tf,string.replace(tf,'ExpressionInput','ExpressionInput/Kallisto_Results'))
        tf = root_dir+'/ExpressionInput/summary.'+dataset_name+'.txt'
        shutil.copyfile(tf,string.replace(tf,'ExpressionInput','ExpressionInput/Kallisto_Results'))
    except:
        print traceback.format_exc()
        pass
    
def calculateGeneTPMs(species,expMatrix):
    import gene_associations
    try:
        gene_to_transcript_db = gene_associations.getGeneToUid(species,('hide','Ensembl-EnsTranscript'))
        if len(gene_to_transcript_db)<10:
            raise ValueError('Ensembl-EnsTranscript file missing, forcing download of this file')
    except Exception:
        try:
            print 'Missing transcript-to-gene associations... downloading from Ensembl.'
            from build_scripts import EnsemblSQL
            db_version = unique.getCurrentGeneDatabaseVersion()
            EnsemblSQL.getGeneTranscriptOnly(species,'Basic',db_version,'yes')
            gene_to_transcript_db = gene_associations.getGeneToUid(species,('hide','Ensembl-EnsTranscript'))
        except Exception:
            from build_scripts import GeneSetDownloader
            print 'Ensembl-EnsTranscripts required for gene conversion... downloading from the web...'
            GeneSetDownloader.remoteDownloadEnsemblTranscriptAssocations(species)
            gene_to_transcript_db = gene_associations.getGeneToUid(species,('hide','Ensembl-EnsTranscript'))
    if len(gene_to_transcript_db)<10:
        print 'NOTE: No valid Ensembl-EnsTranscripts available, proceeding with the analysis of transcripts rather than genes...'
    from import_scripts import OBO_import
    transcript_to_gene_db = OBO_import.swapKeyValues(gene_to_transcript_db)
    
    gene_matrix = {}
    present_gene_transcripts={}
    for transcript in expMatrix:
        if '.' in transcript:
            transcript_alt = string.split(transcript,'.')[0]
        else:
            transcript_alt = transcript
        if transcript_alt in transcript_to_gene_db:
            gene = transcript_to_gene_db[transcript_alt][0]
            try: present_gene_transcripts[gene].append(transcript)
            except Exception: present_gene_transcripts[gene] = [transcript]
        else: pass ### could keep track of the missing transcripts
    for gene in present_gene_transcripts:
        gene_values = []
        for transcript in present_gene_transcripts[gene]:
            gene_values.append(map(float,expMatrix[transcript]))
        gene_tpms = [sum(value) for value in zip(*gene_values)] ### sum of all transcript tmp's per sample
        gene_tpms = map(str,gene_tpms)
        gene_matrix[gene] = gene_tpms
    if len(gene_matrix)>0:
        return gene_matrix
    else:
        print "NOTE: No valid transcript-gene associations available... proceeding with Transcript IDs rather than gene."
        return expMatrix
    
def exportMatrix(eo,headers,matrix,counts=False):
    eo.write(string.join(headers,'\t')+'\n')
    for gene in matrix:
        values = matrix[gene]
        if counts:
            values = map(str,map(int,map(float,values)))
        eo.write(string.join([gene]+values,'\t')+'\n')
    eo.close()

def importTPMs(sample,input_path,expMatrix,countMatrix,countSampleMatrix):
    firstLine=True
    for line in open(input_path,'rU').xreadlines():
        data = cleanUpLine(line)
        if firstLine:
            firstLine=False
            header = string.split(data,'\t')
        else:
            target_id,length,eff_length,est_counts,tpm = string.split(data,'\t')
            try: float(est_counts); 
            except Exception: ### nan instead of float found due to lack of alignment
                est_counts = '0.0'
                tpm = '0.0'
            if '.' in target_id:
                target_id = string.split(target_id,'.')[0] ### Ensembl isoform IDs in more recent Ensembl builds
            try: expMatrix[target_id].append(tpm)
            except Exception: expMatrix[target_id]=[tpm]
            try: countSampleMatrix[sample]+=float(est_counts)
            except Exception: countSampleMatrix[sample]=float(est_counts)
            try: countMatrix[target_id].append(est_counts)
            except Exception: countMatrix[target_id]=[est_counts]
    return expMatrix,countMatrix,countSampleMatrix

def importTotalReadCounts(sample,input_path,sample_total_counts):
    ### Import from Kallisto Json file
    for line in open(input_path,'rU').xreadlines():
        data = cleanUpLine(line)
        if "n_processed: " in data:
            total = string.split(data,"n_processed: ")[1]
            total = string.split(total,',')[0]
            sample_total_counts[sample]=total
    return sample_total_counts

def findPairs(fastq_paths):
    #fastq_paths = ['/Volumes/test/run0718_lane12_read1_index701=Kopan_RBP_02_14999.fastq.gz','/Volumes/run0718_lane12_read2_index701=Kopan_RBP_02_14999.fastq.gz']
    import export
    read_notation=0
    under_suffix_notation=0
    suffix_notation=0
    equal_notation=0
    suffix_db={}
    for i in fastq_paths:
        if 'read1' in i or 'read2' in i or 'pair1' in i or 'pair2' or 'R1' in i or 'R2' in i:
            read_notation+=1
        f = export.findFilename(i)
        if 'fastq' in f:
            name = string.split(f,'fastq')[0]
        elif 'FASTQ' in f:
            name = string.split(f,'FASTQ')[0]
        elif 'fq' in f:
            name = string.split(f,'fq')[0]
        if '_1.' in name or '_2.' in name:
            under_suffix_notation+=1
        elif '1.' in name or '2.' in name:
            suffix_notation+=1
            suffix_db[name[-2:]]=[]
        if '=' in name:
            equal_notation+=1
    if read_notation==0 and suffix_notation==0 and under_suffix_notation==0:
        new_names={}
        for i in fastq_paths:
            if '/' in i or '\\' in i:
                n = export.findFilename(i)
            if '=' in n:
                n = string.split(n,'=')[1]
            new_names[n] = [i]
        ### likely single-end samples
        return new_names, 'single'
    else:
        new_names={}
        paired = 'paired'
        if equal_notation==len(fastq_paths):
            for i in fastq_paths:
                name = string.split(i,'=')[-1]
                name = string.replace(name,'.fastq.gz','')
                name = string.replace(name,'.fastq','')
                name = string.replace(name,'.FASTQ.gz','')
                name = string.replace(name,'.FASTQ','')
                name = string.replace(name,'.fq.gz','')
                name = string.replace(name,'.fq','')
                if '/' in name or '\\' in name:
                    name = export.findFilename(name)
                if '=' in name:
                    name = string.split(name,'=')[1]
                try: new_names[name].append(i)
                except Exception: new_names[name]=[i]
        else:
            for i in fastq_paths:
                if suffix_notation == len(fastq_paths) and len(suffix_db)==2: ### requires that files end in both .1 and .2
                    pairs = ['1.','2.']
                else:
                    pairs = ['-read1','-read2','-pair1','-pair2','_read1','_read2','_pair1','_pair2','read1','read2','pair1','pair2','_1.','_2.','_R1','_R2','-R1','-R2','R1','R2']
                n=str(i)
                n = string.replace(n,'fastq.gz','')
                n = string.replace(n,'fastq','')
                for p in pairs: n = string.replace(n,p,'')
                if '/' in n or '\\' in n:
                    n = export.findFilename(n)
                if '=' in n:
                    n = string.split(n,'=')[1]
                if n[-1]=='.':
                    n = n[:-1] ###remove the last decimal
                try: new_names[n].append(i)
                except Exception: new_names[n]=[i]
        for i in new_names:
            if len(new_names[i])>1:
                pass
            else:
                paired = 'single'
        new_names = checkForMultipleLanes(new_names)
        return new_names, paired
    
def checkForMultipleLanes(new_names):
    """ This function further aggregates samples run across multiple flowcells """
    read_count = 0
    lane_count = 0
    updated_names={}
    for sample in new_names:
        reads = new_names[sample]
        count=0
        for read in reads:
            read_count+=1
            if '_L00' in read and '_001':
                ### assumes no more than 9 lanes/sample
                count+=1
        if len(reads) == count: ### Multiple lanes run per sample
            lane_count+=count
    if lane_count==read_count:
        for sample in new_names:
            sample_v1 = string.replace(sample,'_001','')
            sample_v1 = string.split(sample_v1,'_L00')
            if len(sample_v1[-1])==1: ### lane number
                sample_v1 = sample_v1[0]
                if sample_v1 in updated_names:
                    updated_names[sample_v1]+=new_names[sample]
                else:
                    updated_names[sample_v1]=new_names[sample]
    if len(updated_names)==0:
        updated_names = new_names
    return updated_names
        
def getFASTAFile(species):
    fasta_folder = 'AltDatabase/'+species+'/SequenceData/'
    fasta_files=[]
    dir_list = read_directory(filepath(fasta_folder))
    for file in dir_list:
        if '.fa' in file:
            fasta_files.append(filepath(fasta_folder)+file)
    return fasta_files

if __name__ == '__main__':
    samplesDiffering = 3
    column_method = 'hopach'
    species = 'Hs'
    excludeCellCycle = False
    platform = 'RNASeq'; graphic_links=[('','/Volumes/HomeBackup/CCHMC/PBMC-10X/ExpressionInput/SamplePrediction/DataPlots/Clustering-33k_CPTT_matrix-CORRELATED-FEATURES-iterFilt-hierarchical_cosine_cosine.txt')]
    """
    graphic_links,new_results_file = correlateClusteredGenes(platform,graphic_links[-1][-1][:-4]+'.txt',
                    numSamplesClustered=samplesDiffering,excludeCellCycle=excludeCellCycle,graphics=graphic_links,
                    ColumnMethod=column_method, transpose=True, includeMoreCells=True)
    """
    import UI; import multiprocessing as mlp

    #runKallisto('Mm','BoneMarrow','/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/altanalyze/Mm-FASTQ','/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/altanalyze/Mm-FASTQ',mlp);sys.exit()
    runKallisto('Hs','BreastCancer','/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/BreastCancerDemo/FASTQs/input','/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/BreastCancerDemo/FASTQs/input',mlp);sys.exit()


    results_file = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/l/July-2017/PSI/test/Clustering-exp.round2-Guide3-hierarchical_cosine_correlation.txt'
    #correlateClusteredGenesParameters(results_file,rho_cutoff=0.3,hits_cutoff=4,hits_to_report=50,ReDefinedClusterBlocks=True,filter=True)
    #sys.exit()
    #correlateClusteredGenes('exons',results_file,stringency='strict',rhoCutOff=0.6);sys.exit()
    #sys.exit()
    species='Hs'; platform = "3'array"; vendor = "3'array"
    #FeatureCounts('/Users/saljh8/Downloads/subread-1.5.2-MaxOSX-x86_64/annotation/mm10_AltAnalyze.txt', '/Users/saljh8/Desktop/Grimes/GEC14074/Grimes_092914_Cell12.bam')
    #sys.exit()
    import UI; import multiprocessing as mlp
    gsp = UI.GeneSelectionParameters(species,platform,vendor)
    gsp.setGeneSet('None Selected')
    gsp.setPathwaySelect('')
    gsp.setGeneSelection('')
    gsp.setJustShowTheseIDs('')
    gsp.setNormalize('median')
    gsp.setSampleDiscoveryParameters(1,50,4,4,
        True,'gene','protein_coding',False,'cosine','hopach',0.4)
    #expFile = '/Users/saljh8/Desktop/Grimes/KashishNormalization/test/Original/ExpressionInput/exp.CombinedSingleCell_March_15_2015.txt'
    expFile = '/Volumes/My Passport/salomonis2/SRP042161_GBM-single-cell/bams/ExpressionInput/exp.GBM_scRNA-Seq-steady-state.txt'
    #singleCellRNASeqWorkflow('Hs', "RNASeq", expFile, mlp, parameters=gsp);sys.exit()
    
    filename = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/Trumpp-HSC-2017/counts.rawTrumpp.txt'
    filename = '/Volumes/salomonis2/Erica-data/GSE98451/counts.GSE98451_uterus_single_cell_RNA-Seq_counts-Ensembl.txt'

    #fastRPKMCalculate(filename);sys.exit()
    #calculateRPKMsFromGeneCounts(filename,'Mm',AdjustExpression=False);sys.exit()
    #copyICGSfiles('','');sys.exit()
    import multiprocessing as mlp
    import UI
    species='Mm'; platform = "3'array"; vendor = 'Ensembl'
    gsp = UI.GeneSelectionParameters(species,platform,vendor)
    gsp.setGeneSet('None Selected')
    gsp.setPathwaySelect('')
    gsp.setGeneSelection('')
    gsp.setJustShowTheseIDs('')
    gsp.setNormalize('median')
    gsp.setSampleDiscoveryParameters(0,0,1.5,3,
        False,'PSI','protein_coding',False,'cosine','hopach',0.35)
    
    #gsp.setSampleDiscoveryParameters(1,1,4,3, True,'Gene','protein_coding',False,'cosine','hopach',0.5)
    filename = '/Volumes/SEQ-DATA/AML_junction/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI-clust.txt'
    #fastRPKMCalculate(filename);sys.exit()
    results_file = '/Volumes/SEQ-DATA/Grimes/14018_gmp-pro/ExpressionInput/DataPlots/400 fold for at least 4 samples/Clustering-myeloblast-steady-state-correlated-features-hierarchical_euclidean_cosine-hopach.txt'
    guideGeneFile = '/Volumes/SEQ-DATA/Grimes/14018_gmp-pro/ExpressionInput/drivingTFs-symbol.txt'

    expFile = '/Users/saljh8/Desktop/Grimes/KashishNormalization/3-25-2015/ExpressionInput/exp.CombinedSingleCell_March_15_2015.txt'
    expFile = '/Users/saljh8/Desktop/dataAnalysis/Mm_Kiddney_tubual/ExpressionInput/exp.E15.5_Adult_IRI Data-output.txt'
    expFile = '/Users/saljh8/Desktop/PCBC_MetaData_Comparisons/temp/C4Meth450-filtered-SC-3_regulated.txt'
    expFile = '/Volumes/SEQ-DATA/Grimeslab/TopHat/AltResults/AlternativeOutput/Mm_RNASeq_top_alt_junctions-PSI-clust-filter.txt'
    expFile = '/Users/saljh8/Documents/L_TargetPSIFiles/exp.TArget_psi_noif_uncorr_03-50missing-12high.txt'
    expFile = '/Volumes/BOZEMAN2015/Hs_RNASeq_top_alt_junctions-PSI-clust-filter.txt'

    singleCellRNASeqWorkflow('Hs', "exons", expFile, mlp, exp_threshold=0, rpkm_threshold=0, parameters=gsp);sys.exit()
    
    #expFile = '/Users/saljh8/Desktop/Grimes/AltSplice/Gmp-cluster-filter.txt'
    #singleCellRNASeqWorkflow('Mm', "exons", expFile, mlp, exp_threshold=0, rpkm_threshold=0, parameters=gsp);sys.exit()
    #expFile = '/Users/saljh8/Downloads/methylation/ExpressionInput/exp.female-steady-state.txt'
    
    #singleCellRNASeqWorkflow('Hs', 'RNASeq', expFile, mlp, exp_threshold=50, rpkm_threshold=5) # drivers=guideGeneFile)
    #sys.exit()
    #correlateClusteredGenes(results_file);sys.exit()
    #reformatExonFile('Hs','exon',True);sys.exit()
    filename = '/Volumes/Time Machine Backups/dataAnalysis/PCBC_Sep2013/C4-reference/ExpressionInput/counts.C4.txt'
    #fastRPKMCalculate(filename);sys.exit()
    file1 = '/Volumes/My Passport/dataAnalysis/CardiacRNASeq/BedFiles/ExpressionInput/exp.CardiacRNASeq.txt'
    file2 = '/Volumes/Time Machine Backups/dataAnalysis/PCBC_Sep2013/C4-reference/ReferenceComps/ExpressionInput/counts.C4.txt'
    #getHighExpNovelExons('Hs',file1);sys.exit()
    #mergeCountFiles(file1,file2); sys.exit()
    import UI
    test_status = 'yes'
    data_type = 'ncRNA'
    data_type = 'mRNA'
    array_type = 'RNASeq'
    array_type = 'junction'
    species = 'Hs' ### edit this
    
    summary_results_db = {}

    root_dir = '/Volumes/Time Machine Backups/dataAnalysis/Human Blood/Exon/Multiple Sclerosis/Untreated_MS-analysis/'
    #root_dir = '/Volumes/Time Machine Backups/dataAnalysis/Human Blood/Exon/Multiple Sclerosis/2-3rds_training-untreated/'
    root_dir = '/Volumes/SEQ-DATA/Grimes/14018_gmp-pro/400-original/'
    #root_dir = '/Volumes/My Passport/dataAnalysis/PCBC_Dec2013/All/bedFiles/'
    root_dir = '/Users/saljh8/Desktop/dataAnalysis/HTA2.0 Files/'
    #summary_results_db['Hs_Junction_d14_vs_d7.p5_average-ASPIRE-exon-inclusion-results.txt'] = [] ### edit this
    #summary_results_db['Hs_Junction_d14_vs_d7.p5_average-splicing-index-exon-inclusion-results.txt'] = [] ### edit this
    
    results_dir = root_dir +'AltResults/AlternativeOutput/'
    dir_list = read_directory(results_dir)
    for i in dir_list:
        if '_average' in i:
            comparison, end = string.split(i,'_average')
            if '-exon-inclusion-results.txt' in i: summary_results_db[comparison]=[]

    compareExonAndJunctionResults(species,array_type,summary_results_db,root_dir); sys.exit()
    
    fl = UI.ExpressionFileLocationData('','','',''); fl.setCELFileDir(loc); fl.setRootDir(loc)
    exp_file_location_db={}; exp_file_location_db['test']=fl
    alignJunctionsToEnsembl(species,exp_file_location_db,'test'); sys.exit()
    getEnsemblAssociations(species,data_type,test_status,'yes'); sys.exit()
