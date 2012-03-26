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

import sys, string
import statistics
import math
import os.path
import unique
import update
import copy
import time
import export
import EnsemblImport; reload(EnsemblImport)
import JunctionArrayEnsemblRules
import JunctionArray; reload(JunctionArray)
import ExonArrayEnsemblRules

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

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

######### Below code deals with building the AltDatabase #########

def exportNovelExonToBedCoordinates(species,novel_exon_coordinates):
    ### Export the novel exon coordinates based on those in the junction BED file to examine the differential expression of the predicted novel exon
    #bamToBed -i accepted_hits.bam -split| coverageBed -a stdin -b /home/databases/hESC_differentiation_exons.bed > day20_7B__exons-novel.bed
    bed_export_path = filepath('AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.bed')
    bed_data = open(bed_export_path,'a') ### Appends to existing file
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
            bed_values = [chr,str(start),str(stop),gene,'0',str(ji.Strand())]
            bed_values = cleanUpLine(string.join(bed_values,'\t'))+'\n'
            bed_data.write(bed_values)
    bed_data.close()
    return bed_export_path

def moveBAMtoBEDFile(species,dataset_name,root_dir):
    bed_export_path = filepath('AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.bed')
    new_fn = root_dir+'/BAMtoBED/'+species + '_'+dataset_name+'_exons.bed'
    print 'Writing exon-level coordinates to BED file:'
    print new_fn
    export.customFileMove(bed_export_path,new_fn)
    return new_fn

def reformatExonFile(species,type):
    if type == 'exon':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
        export_path = 'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.txt'
        ### Used by BEDTools to get counts per specific AltAnalyze exon region (should augment with de novo regions identified from junction analyses)
        bed_export_path = 'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.bed'
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
            gene, exonid, chr, strand, start, stop, constitutive_call, ens_exon_ids, splice_events, splice_junctions = t
            if constitutive_call == 'yes': ens_constitutive_status = '1'
            else: ens_constitutive_status = '0'
            export_values = [gene+':'+exonid, exonid, gene, '', chr, strand, start, stop, 'known', constitutive_call, ens_exon_ids, ens_constitutive_status]
            export_values+= [exonid, start, stop, splice_events, splice_junctions]
            export_values = string.join(export_values,'\t')+'\n'; export_data.write(export_values)
            if type == 'exon':
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
    import UCSCImport
    mRNA_Type = 'mrna'; run_from_scratch = 'yes'
    export_all_associations = 'no' ### YES only for protein prediction analysis
    update.buildUCSCAnnoationFiles(species,mRNA_Type,export_all_associations,run_from_scratch,force)
    
    null = EnsemblImport.getEnsemblAssociations(species,data_type,test_status); null=[]
    reformatExonFile(species,'exon'); reformatExonFile(species,'junction')
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
    def __repr__(self): return "JunctionData values"
    
def importBEDFile(bed_dir,root_dir,species,normalize_feature_exp,search_chr,getReads):
    dir_list = read_directory(bed_dir)
    begin_time = time.time()
            
    condition_count_db={}; neg_count=0; pos_count=0; junction_db={}; biotypes={}; algorithms={}; exon_len_db={}; getBiotype=''
    if getReads == 'yes': getBiotype = search_chr
    
    if testImport == 'yes': print "Reading user RNA-seq input data files"
    for filename in dir_list:
        count_db={}; x=0
        fn=filepath(bed_dir+filename)
        condition = export.findFilename(fn)
        if '__' in condition:
            ### Allow multiple junction files per sample to be combined (e.g. canonical and non-canonical junction alignments)
            condition=string.split(condition,'__')[0]+filename[-4:]
        if ('.bed' in fn or '.BED' in fn or '.tab' in fn or '.TAB' in fn) and '._' not in condition:
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
                if x==0 or '#' == data[0]:
                    format_description = data; x=1
                    algorithm = 'Unknown'
                    if 'TopHat' in format_description: algorithm = 'TopHat'
                    elif 'HMMSplicer' in format_description: algorithm = 'HMMSplicer'
                    elif 'SpliceMap junctions' in format_description: algorithm = 'SpliceMap'
                    elif t[0] == 'E1': algorithm = 'BioScope-junction'
                    elif '# filterOrphanedMates=' in data or 'alignmentFilteringMode=' in data or '#number_of_mapped_reads=' in data:
                        algorithm = 'BioScope-exon'
                    if testImport == 'yes': print condition, algorithm
                else:
                    proceed = 'yes' ### restrict by chromosome with minimum line parsing (unless we want counts instead)
                    if getReads == 'no':
                        try:
                            if search_chr == t[0] or search_chr == t[1]: proceed = 'yes'
                            else: proceed = 'no'
                        except IndexError:
                            print 'The input file:\n',filename
                            print 'is not formated as expected (format='+algorithm+').'
                            print 'search chromosome:',search_chr
                            print t; force_bad_exit
                    if proceed == 'yes':
                        proceed = 'no'
                        if '.tab' in fn or '.TAB' in fn:
                            ### Applies to non-BED format Junction and Exon inputs (BioScope)
                            if 'BioScope' in algorithm:
                                if algorithm == 'BioScope-exon': ### Not BED format
                                    chr,source,data_type,start,end,reads,strand,null,gene_info=t[:9]
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
                                        if float(reads)>0: proceed = 'yes'
                                        seq_length = abs(exon1_stop-exon2_start)
                                if algorithm == 'BioScope-junction':
                                    chr = t[1]; strand = t[2]; exon1_stop = int(t[4]); exon2_start = int(t[8]); count_paired = t[17]; count_single = t[19]; score=t[21]
                                    try: exon1_start = int(t[3]); exon2_stop = int(t[9])
                                    except Exception: null=[] ### If missing, these are not assigned
                                    reads = str(int(float(count_paired))+int(float(count_single))) ### Users will either have paired or single read (this uses either)
                                    biotype = 'junction'; biotypes[biotype]=[]; junction_id=''
                                    if float(reads)>0: proceed = 'yes'
                                    seq_length = abs(float(exon1_stop-exon2_start))
                        else:
                            try:
                                ### Applies to BED format Junction input
                                chr, exon1_start, exon2_stop, junction_id, reads, strand, null, null, null, null, lengths, null = t
                                exon1_len,exon2_len=string.split(lengths,',')[:2]; exon1_len = int(exon1_len); exon2_len = int(exon2_len)
                                exon1_start = int(exon1_start); exon2_stop = int(exon2_stop)
                                biotype = 'junction'; biotypes[biotype]=[]
                                if strand == '-':
                                    exon1_stop = exon1_start+exon1_len; exon2_start=exon2_stop-exon2_len+1
                                    ### Exons have the opposite order
                                    a = exon1_start,exon1_stop; b = exon2_start,exon2_stop
                                    exon1_stop,exon1_start = b; exon2_stop,exon2_start = a
                                else:
                                    exon1_stop = exon1_start+exon1_len; exon2_start=exon2_stop-exon2_len+1
                                proceed = 'yes'
                                if algorithm == 'HMMSplicer':
                                    if '|junc=' in junction_id: reads = string.split(junction_id,'|junc=')[-1]
                                    else: proceed = 'no'
                                if algorithm == 'SpliceMap':
                                    if ')' in junction_id and len(junction_id)>1: reads = string.split(junction_id,')')[0][1:]
                                    else: proceed = 'no'
                                seq_length = abs(float(exon1_stop-exon2_start))
                            except Exception:
                                ### Applies to BED format exon input (BEDTools export)
                                # bamToBed -i accepted_hits.bam -split| coverageBed -a stdin -b /home/nsalomonis/databases/Mm_Ensembl_exons.bed > day0_8B__exons.bed
                                try: chr, start, end, exon_id, null, strand, reads, bp_coverage, bp_total, percent_coverage = t
                                except Exception:
                                    print 'The file',fn,'does not appear to be propperly formatted as input.'
                                    print t; force_exception
                                algorithm = 'TopHat-exon'; biotype = 'exon'; biotypes[biotype]=[]
                                exon1_stop,exon2_start = int(start),int(end); junction_id=exon_id; seq_length = float(bp_total)
                                if seq_length == 0:
                                    seq_length = abs(float(exon1_stop-exon2_start))
                                ### Adjust exon positions - not ideal but necessary. Needed as a result of exon regions overlapping by 1nt (due to build process)
                                exon1_stop+=1; exon2_start-=1
                                if float(reads)>0:
                                    proceed = 'yes'
                                else: proceed = 'no'
                        if proceed == 'yes':
                            if strand == '+': pos_count+=1
                            else: neg_count+=1
                            if getReads == 'yes' and seq_length>0:
                                if getBiotype == biotype:
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
            if getReads == 'yes':
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
    if getReads == 'yes':
        print len(exon_len_db), getBiotype, 'read counts present for',algorithm
        return condition_count_db,exon_len_db,biotypes,algorithms
    else:
        if testImport == 'yes':
            if 'exon' not in biotypes and 'BioScope' not in algorithm:
                print len(junction_db),'junctions present in',algorithm,'format BED files.' # ('+str(pos_count),str(neg_count)+' by strand).'
            elif 'exon' in biotypes and 'BioScope' not in algorithm:
                print len(junction_db),'sequence identifiers present in input files.' 
            else: print len(junction_db),'sequence identifiers present in BioScope input files.'
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
                    rpkm = (math.pow(10.0,8.0))*(float(read_count)/(float(average_total_reads)*region_length)) 
                try:
                    if region_length == 0:
                        region_length = abs(int(key[2]-key[1]))
                    rpkm = (math.pow(10.0,8.0))*(read_count/(float(total_mapped_reads)*region_length))
                except Exception:
                    print condition, key, 'v.2.0.6 test'
                    print 'Error Encountered... Exon or Junction of zero length encoutered... RPKM failed... Exiting AltAnalyze.'
                    print [read_count,total_mapped_reads,region_length];k=1; kill
                condition_count_db[condition][key] = str(rpkm) ### Replace original counts with RPMK
        except Exception:
            if k == 1: kill
            null=[]
    return condition_count_db

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
    average_total_reads = average_total_reads/(index+1)
            
    for gene in gene_count_db:
        index=0; rpkms = []
        for (read_count,region_length) in gene_count_db[gene]:
            total_mapped_reads = mapped_reads[index]
            #print [read_count],[region_length],[total_mapped_reads]
            #if gene == 'ENSMUSG00000028186': print [read_count, index, total_mapped_reads,average_total_reads,region_length]
            if read_count == 0: read_count=1; rpkm = (math.pow(10.0,8.0))*(float(read_count)/(float(average_total_reads)*region_length)) ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
            else:
                try: rpkm = (math.pow(10.0,8.0))*(float(read_count+1)/(float(total_mapped_reads)*region_length)) ### read count is incremented +1 (see next line)
                except Exception: read_count=1; rpkm = (math.pow(10.0,8.0))*(float(read_count)/(float(average_total_reads)*region_length)) ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
            #if gene == 'ENSMUSG00000028186': print rpkm,read_count,index,total_mapped_reads,average_total_reads,region_length
            rpkms.append(rpkm)
            index+=1
        gene_count_db[gene] = rpkms ### Replace original counts with RPMK
    return gene_count_db

def calculateGeneLevelStatistics(steady_state_export,expressed_gene_exon_db,normalize_feature_exp,array_names,fl):
    global UserOptions; UserOptions = fl
    exp_file = string.replace(steady_state_export,'-steady-state','')
    if normalize_feature_exp == 'RPKM':
        exp_dbase, array_count = importRawCountData(exp_file,expressed_gene_exon_db,normalize_feature_exp)
        steady_state_db = obtainGeneCounts(expressed_gene_exon_db,exp_dbase,array_count,normalize_feature_exp); exp_dbase=[]
        exportGeneCounts(steady_state_export,array_names,steady_state_db)
        steady_state_db = calculateGeneRPKM(steady_state_db)
    else:
        exp_dbase, array_count = importNormalizedCountData(exp_file,expressed_gene_exon_db,normalize_feature_exp)
        steady_state_db = obtainGeneCounts(expressed_gene_exon_db,exp_dbase,array_count,normalize_feature_exp); exp_dbase=[]
        exportGeneCounts(steady_state_export,array_names,steady_state_db)
    return steady_state_db
    
def exportGeneCounts(steady_state_export,headers,gene_count_db):
    ### In addition to RPKM gene-level data, export gene level counts and lengths (should be able to calculate gene RPKMs from this file)
    export_path = string.replace(steady_state_export,'exp.','counts.')
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

def importGeneCounts(filename):
    ### Import non-normalized original counts and return the max value            
    counts_filename = string.replace(filename,'exp.','counts.')
    status = verifyFile(counts_filename)
    if status == 'not found': ### Occurs for non-normalized counts
        counts_filename = filename
    fn=filepath(counts_filename); x=0; max_count_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: array_names = t[1:]; x=1
        else:
            gene = t[0]
            max_count_db[gene] = str(statistics.maxval(t[1:]))
    return max_count_db

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
    
        
    for gene in gene_count_db:
        index=0; rpkms = []
        for (read_count,region_length) in gene_count_db[gene]:
            total_mapped_reads = mapped_reads[index]
            #print [read_count],[region_length],[total_mapped_reads]
            #if gene == 'ENSMUSG00000028186': print [read_count, index, total_mapped_reads,average_total_reads,region_length]
            if read_count == 0: read_count=1; rpkm = (math.pow(10.0,8.0))*(float(read_count)/(float(average_total_reads)*region_length)) ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
            else:
                try: rpkm = (math.pow(10.0,8.0))*(float(read_count+1)/(float(total_mapped_reads)*region_length)) ### read count is incremented +1 (see next line)
                except Exception: read_count=1; rpkm = (math.pow(10.0,8.0))*(float(read_count)/(float(average_total_reads)*region_length)) ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
            #if gene == 'ENSMUSG00000028186': print rpkm,read_count,index,total_mapped_reads,average_total_reads,region_length
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
    
    if 'exp.' not in dataset_name: dataset_name = 'exp.'+dataset_name
    if '.txt' not in dataset_name: dataset_name+='.txt'
    export_path = root_dir+'ExpressionInput/'+dataset_name
    try: os.remove(filepath(export_path))
    except Exception: null=[]
    try: os.remove(filepath(string.replace(export_path,'exp.','counts.')))
    except Exception: null=[]
       
def alignExonsAndJunctionsToEnsembl(species,exp_file_location_db,dataset_name):
    global testImport
    testImport = 'no'
    
    rnaseq_begin_time = time.time()
    
    fl = exp_file_location_db[dataset_name]
    bed_dir=fl.BEDFileDir()
    root_dir=fl.RootDir()
    ### Import experimentally identified junction splice-sites
    normalize_feature_exp = fl.FeatureNormalization()
    if fl.ExonBedBuildStatus() == 'yes':
        reformatExonFile(species,'exon') ### exports BED format exons for exon expression extraction
        
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
    count=0; novel_junction_count=0; one_found=0; not_found=0; both_found=0; trans_splicing_reads=0; junctions_without_exon_gene_alignments=0
    exons_without_gene_alignment_count=0; junction_simple_db={}; ps=0
    
    ###Begin Chromosome specific read to exon alignments
    chr_strand_gene_db,location_gene_db,chromosomes,gene_location_db = getChromosomeStrandCoordinates(species)
    print "Processing exon/junction coordinates sequentially by chromosome..."
    for search_chr in chromosomes:
      if 'linux' in sys.platform or testImport == 'yes': print 'parsing data for chromosome:',search_chr ### Unix platforms are not displaying the progress in real-time
      else: print "*",
      junction_annotations={}
      junction_db,biotypes,algorithms = importBEDFile(bed_dir,root_dir,species,normalize_feature_exp,search_chr,'no')
      if len(biotypes)>0: biotypes_store = biotypes
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
                    count +=1
                else: novel_junction_db[key]=junction_db[key]; novel_junction_count+=1
            else:
                unmapped_exon_db[key]=junction_db[key]

        ens_exon_db = importExonAnnotations(species,'exon',search_chr)
        
        if 'junction' in biotypes:
            if testImport == 'yes':
                print count,  'junctions found in Ensembl/UCSC and',len(novel_junction_db),'are novel.'
            
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
                    both_found+=1; #print 'five',(chr,exon1_stop,exon2_start),exon1_coord_db[(chr,exon1_stop)]
                else:
                    if (chr,exon1_stop) in exon1_coord_db:
                        le=exon1_coord_db[(chr,exon1_stop)]; ji.setLeftExonAnnotations(le)
                        one_found+=1; #print 'three',(chr,exon1_stop,exon2_start),exon1_coord_db[(chr,exon1_stop)]
                        novel_exon_coordinates[ji.Chr(),exon1_stop] = ji,'right',ji.Exon1Start() ### Employ this strategy to avoid duplicate exons with differing lengths (mainly an issue if analyzing only exons results)
                    elif (chr,exon2_start) in exon2_coord_db:
                        re=exon2_coord_db[(chr,exon2_start)]; ji.setRightExonAnnotations(re) ### In very rare cases, a gene can be assigned here, even though the splice-site is on the opposite strand (not worthwhile filtering out)
                        one_found+=1; #print 'three',(chr,exon1_stop,exon2_start),exon1_coord_db[(chr,exon1_stop)]
                        novel_exon_coordinates[ji.Chr(),exon2_start] = ji,'left',ji.Exon2Stop()
                    else:
                        not_found+=1; #if not_found < 10: print (chr,exon1_stop,exon2_start)
                        novel_exon_coordinates[ji.Chr(),exon1_stop] = ji,'right',ji.Exon1Start()
                        novel_exon_coordinates[ji.Chr(),exon2_start] = ji,'left',ji.Exon2Stop()
                    ### We examine reads where one splice-site aligns to a known but the other not, to determine if trans-splicing occurs
                    try: exons_not_identified[chr,ji.Strand()].append((coord,ji))
                    except KeyError: exons_not_identified[chr,ji.Strand()] = [(coord,ji)]
            """
            if fl.ExonBedBuildStatus() == 'no':
                exportNovelJunctions(species,novel_junction_db,condition_count_db,root_dir,dataset_name,'junction') ### Includes known exons
            """
                
            #print both_found, ' where both and', one_found, 'where one splice-site are known out of',both_found+one_found+not_found
            #print 'Novel junctions where both splice-sites are known:',both_found
            #print 'Novel junctions where one splice-site is known:',one_found
            #print 'Novel junctions where the splice-sites are not known:',not_found
            clearObjectsFromMemory(exon_region_db); exon_region_db=[] ### Clear memory of this object

            read_aligned_to_gene=0
            for (chr,strand) in exons_not_identified:
                if (chr,strand) in chr_strand_gene_db:
                    chr_gene_locations = chr_strand_gene_db[chr,strand]
                    chr_reads = exons_not_identified[chr,strand]
                    chr_gene_locations.sort(); chr_reads.sort()
                    read_aligned_to_gene=geneAlign(chr,chr_gene_locations,location_gene_db,chr_reads,'no',read_aligned_to_gene)
            #print read_aligned_to_gene, 'novel junctions aligned to Ensembl genes out of',one_found+not_found
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
                                trans_splicing_reads+=1
                                if ji.checkExonPosition(coordinate) == 'right': geneid = ji.SecondaryGeneID()
                            exon_data = (coordinate,ji.Chr()+'-'+str(coordinate),'novel')
                            try: novel_exon_db[geneid].append(exon_data)
                            except KeyError: novel_exon_db[geneid] = [exon_data]
                            
                else: junctions_without_exon_gene_alignments+=1
    
            ### Remove redundant exon entries and store objects    
            for key in novel_exon_db:
                exon_data_objects=[]
                exon_data_list = unique.unique(novel_exon_db[key])
                exon_data_list.sort()
                for e in exon_data_list:
                    ed = ExonInfo(e[0],e[1],e[2])
                    exon_data_objects.append(ed)
                novel_exon_db[key] = exon_data_objects
                
            #print trans_splicing_reads,'trans-splicing junctions found (two aligning Ensembl genes).'
            #print junctions_without_exon_gene_alignments, 'junctions where neither splice-site aligned to a gene'
            #if 'X' in search_chr: print len(ens_exon_db),len(ens_exon_db['ENSMUSG00000044424'])
            alignReadsToExons(novel_exon_db,ens_exon_db)
            
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
            
            if fl.ExonBedBuildStatus() == 'yes':
                ### Append to the exported BED format exon coordinate file
                bedfile = exportNovelExonToBedCoordinates(species,novel_exon_coordinates)
            
            ### Identify reciprocol junctions and retrieve splice-event annotations for exons and inclusion junctions
            junction_annotations,critical_exon_annotations = JunctionArray.inferJunctionComps(species,('RNASeq',junction_region_db,root_dir))
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
                if (chr,strand) in chr_strand_gene_db:
                    chr_gene_locations = chr_strand_gene_db[chr,strand]
                    chr_reads = exons_not_aligned[chr,strand]
                    chr_gene_locations.sort(); chr_reads.sort()
                    read_aligned_to_gene=geneAlign(chr,chr_gene_locations,location_gene_db,chr_reads,'no',read_aligned_to_gene)
                    
            #read_aligned_to_gene, 'exons aligned to Ensembl genes out of',one_found+not_found
                        
            align_exon_db={}; exons_without_gene_alignments={}
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
                        if ji.GeneID() not in ji.JunctionID(): ### Hence, there were probably two overlapping Ensembl genes and the wrong was assigned based on the initial annotaitons
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
                else: exons_without_gene_alignments[key]=ji; exons_without_gene_alignment_count+=1
            
            ### Remove redundant exon entries and store objects (this step may be unnecessary)
            for key in align_exon_db:
                exon_data_objects=[]
                exon_data_list = unique.unique(align_exon_db[key])
                exon_data_list.sort()
                for e in exon_data_list:
                    ed = ExonInfo(e[0],e[1],e[2])
                    exon_data_objects.append(ed)
                align_exon_db[key] = exon_data_objects
                
            #print exons_without_gene_alignment_count, 'exons where neither aligned to a gene'
            #if exons_without_gene_alignment_count>3000: print 'NOTE: Poor mapping of these exons may be due to an older build of\nEnsembl than the current version. Update BAMtoBED mappings to correct.'
    
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
            exportDatasetLinkedExons(species,exons_to_export,critical_exon_annotations,root_dir)            
                    
        ### Export the novel user exon-junction annotations (original junction_db objects updated by above processing)
        exportDatasetLinkedJunctions(species,junction_db,junction_annotations,root_dir)

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
            if ji.GeneID()!=None and ji.UniqueID()!=None: junction_simple_db[key]=ji.UniqueID()

        #returnLargeGlobalVars()
        clearObjectsFromMemory(junction_db); clearObjectsFromMemory(junction_annotations)
        junction_db=[]; junction_annotations=[]; chr_reads=[]

    exportDatasetLinkedGenes(species,gene_location_db,root_dir) ### Include an entry for gene IDs to include constitutive expression for RPKM normalized data     
    chr_gene_locations=[]; location_gene_db=[]; chr_strand_gene_db=[] 
    
    ### END CHROMOSOME SPECIFIC ANALYSES
    if fl.ExonBedBuildStatus() == 'yes':
        bedfile = moveBAMtoBEDFile(species,dataset_name,root_dir)
        print 'Exon BED file updated with novel exon predictions from junction file'
        return bedfile; sys.exit()

    print 'user coordinates imported/processed'
    print 'Importing read counts from coordinate data...'
    ### Output, normalize and annotate read count data for exons and junctions separately
    try: biotypes = biotypes_store
    except Exception: return 'failed'
    
    for biotype in biotypes:
        ### Import Read Counts (do this last to conserve memory, but do it for all)
        condition_count_db,exon_len_db,biotypes2,algorithms = importBEDFile(bed_dir,root_dir,species,normalize_feature_exp,biotype,'yes')

        ###First export original counts, rather than quantile normalized or RPKM
        exportJunctionCounts(species,junction_simple_db,exon_len_db,condition_count_db,root_dir,dataset_name,biotype,'counts')
    
        if normalize_feature_exp != 'none' and biotype == 'junction':
            if normalize_feature_exp == 'RPKM':
                print 'Normalizing junction expression (RPKM analogue - 60nt length)...',
                condition_count_db = calculateRPKM(condition_count_db,exon_len_db,'junction')
                print 'finished'
            elif normalize_feature_exp == 'quantile':
                print 'Normalizing junction expression (quantile)...',
                condition_count_db = quantileNormalizationSimple(species,condition_count_db,'junction')
                print 'finished'
        elif biotype == 'junction':
            condition_count_db = adjustCounts(condition_count_db,exon_len_db)
            
        if normalize_feature_exp != 'none' and biotype == 'exon':
            if normalize_feature_exp == 'RPKM':
                print 'Normalizing exon expression (RPKM)...',
                condition_count_db = calculateRPKM(condition_count_db,exon_len_db,'exon')
                print 'finished'
            elif normalize_feature_exp == 'quantile':
                print 'Normalizing exon expression (quantile)...',
                condition_count_db = quantileNormalizationSimple(species,condition_count_db,'exon')
                print 'finished'
        elif biotype == 'exon':
            condition_count_db = adjustCounts(condition_count_db,exon_len_db)
            
        ### export user results with junction names and counts
        exportJunctionCounts(species,junction_simple_db,exon_len_db,condition_count_db,root_dir,dataset_name,biotype,'normalized')
        clearObjectsFromMemory(condition_count_db); clearObjectsFromMemory(exon_len_db); condition_count_db=[]; exon_len_db=[]
    
    ### Output summary statistics
    if normalize_feature_exp != 'none':
        print normalize_feature_exp, 'normalization complete'
    if 'junction' in biotypes:
        print 'Imported Junction Statistics:'
        print '    ',count, 'junctions found in Ensembl/UCSC and',novel_junction_count,'are novel'
        print '    ',trans_splicing_reads,'trans-splicing junctions found (two aligning Ensembl genes)'
        print '    ',junctions_without_exon_gene_alignments, 'junctions where neither splice-site aligned to a gene'
    if 'exon' in biotypes:
        print 'Imported Exon Statistics:'
        print '    ',exons_without_gene_alignment_count, 'exons where neither aligned to a gene'
    print 'User databases and read counts written to:', root_dir[:-1]+'ExpressionInput'

    clearObjectsFromMemory(junction_simple_db); junction_simple_db=[]
    rnaseq_end_time = time.time()
    print '...RNA-seq import completed in',int(rnaseq_end_time-rnaseq_begin_time),'seconds\n'
    """
    returnLargeGlobalVars()
    print 'local vars'
    all = [var for var in locals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all:
        try:
            if len(locals()[var])>1: print var, len(locals()[var])
        except Exception: null=[]
    """
    return biotypes

def importRawCountData(filename,expressed_gene_exon_db,normalize_feature_exp):
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
    rpkm_threshold = UserOptions.RPKMThreshold()
    junction_exp_threshold = UserOptions.ExonExpThreshold()
    
    ### Import RPKM normalized expression values               
    fn=filepath(filename); x=0; rpkm_dbase={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: array_names = t[1:]; x=1
        else:
            exon_id=t[0]
            max_count=max(map(float,t[1:]))
            if max_count>=rpkm_threshold: rpkm_dbase[exon_id]=[] ### Only retain exons/junctions meeting the RPKM threshold
            
    ### Import non-normalized original counts                
    counts_filename = string.replace(filename,'exp.','counts.')
    fn=filepath(counts_filename); x=0; exp_dbase={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: array_names = t[1:]; x=1
        else:
            exon_id,coordinates = string.split(t[0],'=')
            try:
                null = seq_ids_to_import[exon_id]
                coordinates = string.split(coordinates,':')[1]
                coordinates = string.split(coordinates,'-')
                length=abs(int(coordinates[1])-int(coordinates[0]))
                max_count=max(map(float,t[1:])); proceed = 'no'
                if '-' in exon_id:
                    length = 60.0
                    if max_count>=junction_exp_threshold: proceed = 'yes'
                elif max_count>=exon_exp_threshold: proceed = 'yes'
                if proceed == 'yes' and exon_id in rpkm_dbase: ### Ensures that the maximum sample (not group) user defined count threshold is achieved at the exon or junction-level
                    exp_dbase[exon_id] = t[1:],length ### Include sequence length for normalization
            except Exception: null=[]

    for exon in exp_dbase: array_count = len(exp_dbase[exon][0]); break
    try:null=array_count
    except Exception:
        print 'No exons or junctions considered expressed (based user thresholds). Exiting analysis.'; force_exit
    return exp_dbase, array_count

def importNormalizedCountData(filename,expressed_gene_exon_db,normalize_feature_exp):
    ### Get expression values for exon/junctions to analyze
    seq_ids_to_import={}
    for gene in expressed_gene_exon_db:
        for exonid in expressed_gene_exon_db[gene]: seq_ids_to_import[exonid]=[]
            
    ### Define thresholds
    exon_exp_threshold = UserOptions.ExonExpThreshold()
    junction_exp_threshold = UserOptions.ExonExpThreshold()
    
    ### Import non-normalized original counts                
    fn=filepath(filename); x=0; exp_dbase={}
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
                exp_dbase[exon_id] = t[1:],0  ### Add the zero just to comply with the raw count input format (indicates exon length)

    for exon in exp_dbase: array_count = len(exp_dbase[exon][0]); break        
    return exp_dbase, array_count

def obtainGeneCounts(expressed_gene_exon_db,exp_dbase,array_count,normalize_feature_exp):
    ###Calculate avg expression for each sample for each exon (using constitutive or all exon values)

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
                    ### Add only one avg-expression value for each array, this loop
                    try: steady_state_db[gene].append((sum_const_exp,sum_length))
                    except KeyError: steady_state_db[gene] = [(sum_const_exp,sum_length)]
                else:
                    avg_const_exp=statistics.avg(exp_list)
                    if avg_const_exp != 1: gene_sum+=avg_const_exp
                    ### Add only one avg-expression value for each array, this loop
                    try: steady_state_db[gene].append(avg_const_exp)
                    except KeyError: steady_state_db[gene] = [avg_const_exp]
            except ZeroDivisionError: null=[] ### Occurs when processing a truncated dataset (for testing usually) - no values for the gene should be included
            x += 1
        if gene_sum==0:
            try:
                del steady_state_db[gene] ### Hence, no genes showed evidence of expression (most critical for RNA-Seq)
            except Exception: null=[] ### Error occurs when a gene is added to the database from location_gene_db, but is not expressed
            
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
            
def getChromosomeStrandCoordinates(species):
    ### For novel junctions with no known-splice site, map to genes
    gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,'RNASeq','key_by_array')
    
    chr_strand_gene_db = {}; location_gene_db = {}; chromosomes={}; all_chromosomes={}
    for gene in gene_location_db:
        chr,strand,start,end = gene_location_db[gene]
        location_gene_db[chr,int(start),int(end)] = gene,strand
        try: chr_strand_gene_db[chr,strand].append((int(start),int(end)))
        except KeyError: chr_strand_gene_db[chr,strand] = [(int(start),int(end))]
        if testImport == 'yes':
            if chr=='chr1': chromosomes[chr]=[]
            if chr=='chr19': chromosomes[chr]=[] ### Gene rich chromosome
        elif len(chr)<7: chromosomes[chr]=[]
        all_chromosomes[chr]=[]
    ### Some organisms aren't organized into classical chromosomes (why I don't know)
    if len(chromosomes)<10 and len(all_chromosomes)>9 and testImport=='no': chromosomes = all_chromosomes
    return chr_strand_gene_db,location_gene_db,chromosomes,gene_location_db
    
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

def exportJunctionCounts(species,junction_simple_db,exon_len_db,condition_count_db,root_dir,dataset_name,biotype,count_type):
    if 'exp.' not in dataset_name: dataset_name = 'exp.'+dataset_name
    if '.txt' not in dataset_name: dataset_name+='.txt'
    if count_type == 'counts': dataset_name = string.replace(dataset_name,'exp.','counts.') ### separately export counts
    
    export_path = root_dir+'ExpressionInput/'+dataset_name
    if testImport == 'yes':
        print 'Writing',export_path
    export_data,status = AppendOrWrite(export_path)
    
    if status == 'not found':
        title = ['AltAnalyze_ID']
        for condition in condition_count_db: title.append(condition)
        export_data.write(string.join(title,'\t')+'\n')
    
    for key in junction_simple_db:
        chr,exon1_stop,exon2_start = key
        if biotype == 'junction':
            coordinates = chr+':'+str(exon1_stop)+'-'+str(exon2_start)
        elif biotype == 'exon':
            coordinates = chr+':'+str(exon1_stop-1)+'-'+str(exon2_start+1)
        try:
            null=exon_len_db[key]
            if count_type == 'counts': values = [junction_simple_db[key]+'='+coordinates]
            else: values = [junction_simple_db[key]]
            for condition in condition_count_db:
                count_db = condition_count_db[condition]
                try: read_count = count_db[key]
                except KeyError: read_count = '0'
                values.append(read_count)
            export_data.write(string.join(values,'\t')+'\n')
        except Exception: null=[]
    export_data.close()
    
def quantileNormalizationSimple(species,condition_count_db,biotype):
    ### Basic quantile normalization method (average ranked expression values)

    ### Get all junction or exon entries
    key_db={}
    for condition in condition_count_db:
        count_db = condition_count_db[condition]
        for key in count_db: key_db[key]=[]
            
    import statistics
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
        quantile_normalize_db[rank] = statistics.avg(quantile_normalize_db[rank])

    for condition in condition_unnormalized_db:
        rank=0
        count_db = condition_count_db[condition]
        for (count,key) in condition_unnormalized_db[condition]:
            avg_count = quantile_normalize_db[rank]
            rank+=1
            count_db[key] = str(avg_count) ### re-set this value to the normalized value
            
    clearObjectsFromMemory(condition_unnormalized_db); condition_unnormalized_db =  []
    clearObjectsFromMemory(quantile_normalize_db); quantile_normalize_db = []
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

def exportDatasetLinkedExons(species,exons_to_export,critical_exon_annotations,root_dir):
    export_path = root_dir+'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.txt'
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
    export_path = root_dir+'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_junctions.txt'
    export_data,status = AppendOrWrite(export_path)
    for gene in gene_location_db:
        chr,strand,start,end = gene_location_db[gene]
        export_values = [gene, 'E0.1',gene, '', chr, strand, str(start), str(end), 'known', 'yes', gene, '1']
        export_values+= ['E0.1', str(start), str(end), '', '']
        export_values = string.join(export_values,'\t')+'\n'; export_data.write(export_values)
    export_data.close()
               
def exportDatasetLinkedJunctions(species,junction_db,junction_annotations,root_dir):
    export_path = root_dir+'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_junctions.txt'
    if testImport == 'yes': print 'Writing',export_path
    export_data,status = AppendOrWrite(export_path)
    
    if status == 'not found': 
        export_title = ['AltAnalyzeID','exon_id','ensembl_gene_id','transcript_cluster_id','chromosome','strand','probeset_start','probeset_stop']
        export_title +=['class','constitutive_probeset','ens_exon_ids','ens_constitutive_status','exon_region','exon-region-start(s)','exon-region-stop(s)','splice_events','splice_junctions']
        export_title = string.join(export_title,'\t')+'\n'; export_data.write(export_title)
    for key in junction_db:
        (chr,exon1_stop,exon2_start) = key
        ji=junction_db[key]
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
                except Exception: print ji.GeneID(), key;kill
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
                ji.setUniqueID(uid)
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
                    ed1 = novel_exon_lookup_db[e1]; red1 = ed1.ExonRegionData(); gene1 = e1[0]
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
    return junction_region_db,exons_to_export

def alignReadsToExons(novel_exon_db,ens_exon_db):
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
    index = 0 ### Don't examine genes already looked at
    genes_assigned = 0; trans_splicing=[]
    for (coord,ji) in chr_reads:
        if index >5: index -=5 ### It is possible for some genes to overlap, so set back the index of genomically ranked genes each time
        gene_id_obtained = 'no'
        if switch_coord == 'no': rs,re=coord
        else: re,rs=coord
        while index < len(chr_gene_locations):                 
            cs,ce = chr_gene_locations[index]
            if cs <= rs and ce >= rs:
                gene,strand = location_gene_db[chr,cs,ce]
                if switch_coord == 'yes':
                    first_geneid = ji.GeneID()
                    ji.setTransSplicing()
                    side = ji.checkExonPosition(rs)
                    if side == 'left':
                        ji.setGeneID(gene)
                        ji.setSecondaryGeneID(first_geneid)
                    else:
                        ji.setSecondaryGeneID(gene)
                        #if ji.GeneID() == None: print 'B',coord, ji.GeneID(), econdaryGeneID()
                        #print ji.GeneID(), ji.SecondaryGeneID();kill
                    genes_assigned+=1; gene_id_obtained = 'yes'
                else:
                    ji.setGeneID(gene);  gene_id_obtained = 'yes'
                    if cs <= re and ce >= re: genes_assigned+=1
                    else: trans_splicing.append((coord,ji))
                break
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
    if switch_coord == 'no' and len(trans_splicing)>0:
        read_aligned_to_gene = geneAlign(chr,chr_gene_locations,location_gene_db,trans_splicing,'yes',read_aligned_to_gene)
    return read_aligned_to_gene

def compareExonAndJunctionResults(species,array_type,summary_results_db,root_dir):
    results_dir = root_dir +'AltResults/AlternativeOutput/'
    dir_list = read_directory(results_dir)
    filtered_dir_db={}
    for comparison_file in summary_results_db:
        for results_file in dir_list:
            if (comparison_file in results_file and '-exon-inclusion-results.txt' in results_file) and ('comparison' not in results_file):
                try: filtered_dir_db[comparison_file].append(results_file)
                except Exception: filtered_dir_db[comparison_file] = [results_file]
    for comparison_file in filtered_dir_db:
        alt_result_files = filtered_dir_db[comparison_file]
        #print alt_result_files, comparison_file
        importAltAnalyzeExonResults(alt_result_files,results_dir)
        
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
    def Dataset(self): return self.dataset
    def Symbol(self): return self.symbol
    def Description(self): return self.description
    def ExonID(self): return self.exonid
    def appendExonID(self,exonid): self.exonid+='|'+exonid
    def Probesets(self): return self.probesets
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
    
def importAltAnalyzeExonResults(dir_list,results_dir):
    regulated_critical_exons={}
    print "Reading AltAnalyze results file"
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
                geneid = t[0]; exonid = t[4]; probeset1 = t[7]; probeset2 = ''; score = t[1][:4]; symbol = t[2]; description = t[3]; regions = t[-4]; direction = t[5]
                genomic_loc = t[-1]; splicing_event = t[-3]; external_exon = t[-6]; gene_exp = t[-8]; protein_annot = t[14]; domain_inferred = t[15]; domain_overlap = t[17]
                if ':' in geneid: geneid = string.split(geneid,':')[0] ### User reported that gene:gene was appearing and not sure exactly where or why but added this to address it
                if 'FIRMA' in fn: method = 'FIRMA'
                else: method = 'splicing-index'
                if 'ASPIRE' in filename or 'linearregres' in filename:
                    f1=float(t[12]); f2=float(t[13]); probeset1 = t[8]; probeset2 = t[10]; direction = t[6]; exonid2 = t[5]; splicing_event = t[-4]
                    protein_annot = t[19]; domain_inferred = t[20]; domain_overlap = t[24]; method = 'linearregres'; regions = t[-5]
                    if 'ASPIRE' in filename: method = 'ASPIRE'
                    if '-' not in exonid: exonid=None ### Occurs when the inclusion just in an exon (possibly won't indicate confirmation so exclude)
                    else: exonid = exonid+' vs. '+exonid2
                    if 'AltMouse' in filename:
                        try: geneid = altmouse_ensembl_db[geneid]
                        except Exception: geneid = geneid
                    if 'RNASeq' not in filename and 'junction' not in filename: regions = string.replace(regions,'-','.')
                probesets = [probeset1,probeset2]
                if ((method == 'splicing-index' or method == 'FIRMA') and '-' in regions) or exonid == None:
                    null=[] ### Don't consider
                else:
                    regions = string.replace(regions,';','|')
                    regions = string.replace(regions,'-','|')
                    regions = string.split(regions,'|')
                    for region in regions:
                        if len(region) == 0:
                            try: region = t[17]+t[18] ### For junction introns where no region ID exists
                            except Exception: null=[]
                        if ':' in region: region = string.split(region,':')[-1] ### User reported that gene:gene was appearing and not sure exactly where or why but added this to address it
                        uid = geneid+':'+region
                        ss = SplicingData(score,symbol,description,exonid,probesets,direction,splicing_event,external_exon,genomic_loc,gene_exp,protein_annot,domain_inferred,domain_overlap,method,filename)
                        try: regulated_critical_exon_temp[uid].append(ss)
                        except Exception: regulated_critical_exon_temp[uid] = [ss]
        for uid in regulated_critical_exon_temp:
            report=None
            if len(regulated_critical_exon_temp[uid])>1:
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
        header = string.join(['uid','symbol','description','exonids','independent confirmation','score','regulation direction','alternative exon annotations','associated isoforms','inferred regulated domains','overlapping domains','method','supporting evidence score'],'\t')+'\n'
        export_data.write(header)

    print len(regulated_critical_exons), 'regulated exon IDs imported.\n'
    print 'writing:',export_path; n=0
    
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
    for uid in regulated_critical_exons:
        #if 'ENSG00000133083' in uid: print [uid]
        exon_level_confirmation = 'no'
        ls = regulated_critical_exons[uid]
        jd = regulated_critical_exons[uid][0]
        if len(ls)>1:
            methods = []; scores = []; direction = []; exonids = []; evidence = 0
            junction_data_found = 'no'; exon_data_found = 'no'
            for jd in ls:
                if jd.Method() == 'ASPIRE' or jd.Method() == 'linearregres': 
                    junction_data_found = 'yes'
                    methods.append(jd.Method())
                    scores.append(jd.Score())
                    direction.append(jd.Direction())
                    exonids.append(jd.ExonID())
                    evidence+=jd.Evidence()
            for ed in ls:
                if ed.Method() == 'splicing-index' or ed.Method() == 'FIRMA':
                    exon_data_found = 'yes' ### pick one of them
                    methods.append(ed.Method())
                    scores.append(ed.Score())
                    direction.append(ed.Direction())
                    exonids.append(ed.ExonID())
                    evidence+=ed.Evidence()
            if junction_data_found == 'yes' and exon_data_found == 'yes': exon_level_confirmation = 'yes';
            method = string.join(methods,'|')
            unique_direction = unique.unique(direction)
            if len(unique_direction) == 1: direction = unique_direction[0]
            else: direction = string.join(direction,'|')
            score = string.join(scores,'|')
            exonids_unique = unique.unique(exonids)
            if len(exonids_unique) == 1: exonids = exonids_unique[0]
            else: exonids = string.join(exonids,'|')
            evidence = str(evidence)
            if 'mutual' in direction: direction = 'mutual'
        if len(ls) == 1:
            direction = jd.Direction()
            score = jd.Score()
            method = jd.Method()
            exonids = jd.ExonID()
            evidence = jd.Evidence()
        values = [uid, jd.Symbol(), jd.Description(), exonids, exon_level_confirmation, score, direction, jd.SplicingEvent()]
        values += [jd.ProteinAnnotation(), jd.DomainInferred(), jd.DomainOverlap(), method, str(evidence)]
        values = string.join(values,'\t')+'\n'
        export_data.write(values); n+=1
    print n,'exon IDs written to file.'
    export_data.close()
    clearObjectsFromMemory(regulated_critical_exons)
    clearObjectsFromMemory(regulated_critical_exons_copy)
    #print '!!!!Within comparison evidence'
    #returnLargeGlobalVars()

if __name__ == '__main__':
    import UI
    test_status = 'yes'
    data_type = 'ncRNA'
    data_type = 'mRNA'
    array_type = 'junction'
    species = 'Hs' ### edit this
    
    summary_results_db = {}
    root_dir = '/Users/nsalomonis/Desktop/ExonInclusions/' ### edit this
    root_dir = '/Users/nsalomonis/Desktop/hESC-splicing-factor-KD/'#HJAY_Data/HJAY_Data/'
    root_dir = '/home/blumerjb/BEDS_J/BEDS_J3/'
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
