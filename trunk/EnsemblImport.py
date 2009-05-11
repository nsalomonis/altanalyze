###EnsemblImport
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
import alignToKnownAlt

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

def exportSubGeneViewerData(exon_regions,critical_gene_junction_db,intron_region_db,intron_retention_db):
    intron_retention_db2={}
    for (gene,chr,strand) in intron_retention_db:
        for intron_info in intron_retention_db[(gene,chr,strand)]:
            pos1,pos2,ed = intron_info; pos_list=[pos1,pos2]; pos_list.sort()
            try: intron_retention_db2[gene].append(pos_list)
            except KeyError: intron_retention_db2[gene] = [pos_list]
            
    exon_annotation_export = 'AltDatabase/ensembl/'+species+'/'+species+'_SubGeneViewer_exon-structure-data.txt'
    fn=filepath(exon_annotation_export); data = open(fn,'w')
    title = ['gene','exon-id','type','block','region','constitutive','start-exon','annotation']
    title = string.join(title,'\t')+'\n'; data.write(title)
    for gene in exon_regions:
        block_db = exon_regions[gene]
        try:intron_block_db = intron_region_db[gene]; introns = 'yes'
        except KeyError: introns = 'no'
        utr_data = [gene,'U0.1','u','0','1','n','n','']
        values = string.join(utr_data,'\t')+'\n'; data.write(values)
        index=1
        for block in block_db:
            for rd in block_db[block]:
                splice_event = rd.AssociatedSplicingEvent();exon_pos = [rd.ExonStart(), rd.ExonStop()]; exon_pos.sort()
                if gene in intron_retention_db2:
                    for retained_pos in intron_retention_db2[gene]:
                        if exon_pos == retained_pos: splice_event = 'exon-region-exclusion'
                        elif exon_pos[0] == retained_pos[0] or exon_pos[1] == retained_pos[1]: splice_event = 'exon-region-exclusion'
                        elif exon_pos[0]>retained_pos[0] and exon_pos[0]<retained_pos[1] and exon_pos[1]>retained_pos[0] and exon_pos[1]<retained_pos[1]: splice_event = 'exon-region-exclusion'
                values = [gene,rd.ExonRegionID2(),'e',str(index),str(rd.RegionNumber()),rd.ConstitutiveCall(),'n',splice_event]#,str(exon_pos[0]),str(exon_pos[1])]
                values = string.join(values,'\t')+'\n'; data.write(values)
            index+=1
            if introns == 'yes':
                try:
                    intronid = rd.IntronRegionID()
                    for rd in intron_block_db[block]:
                        intron_pos = [rd.ExonStart(), rd.ExonStop()]; intron_pos.sort()
                        splice_event = rd.AssociatedSplicingEvent()
                        if gene in intron_retention_db2:
                            for retained_pos in intron_retention_db2[gene]:
                                #if '15' in intronid: print intron_pos,retained_pos;kill
                                if intron_pos == retained_pos: splice_event = 'intron-retention'
                                elif intron_pos[0] == retained_pos[0] or intron_pos[1] == retained_pos[1]: splice_event = 'intron-retention'
                                elif intron_pos[0]>retained_pos[0] and intron_pos[0]<retained_pos[1] and intron_pos[1]>retained_pos[0] and intron_pos[1]<retained_pos[1]: splice_event = 'intron-retention'                            
                        values = [gene,intronid,'i',str(index),str(rd.RegionNumber()),rd.ConstitutiveCall(),'n',splice_event]#,str(intron_pos[0]),str(intron_pos[1])]
                        values = string.join(values,'\t')+'\n'; data.write(values)
                    index+=1
                except KeyError: null=[]
        last_exon_region,null = string.split(rd.ExonRegionID2(),'.') ### e.g. E13.1	becomes E13, 1
        
        utr_data = [gene,'U'+last_exon_region[1:]+'.1','u',str(index),'1','n','n','']
        values = string.join(utr_data,'\t')+'\n'; data.write(values)
    data.close()

    critical_gene_junction_db = eliminate_redundant_dict_values(critical_gene_junction_db)
    exon_annotation_export = 'AltDatabase/ensembl/' +species+'/'+species+ '_SubGeneViewer_junction-data.txt'
    fn=filepath(exon_annotation_export); data = open(fn,'w')
    title = ['gene',"5'exon-region","3'exon-region"]
    title = string.join(title,'\t')+'\n'; data.write(title)
    for gene in critical_gene_junction_db:
        for junction_ls in critical_gene_junction_db[gene]:
            values = [gene,junction_ls[0],junction_ls[1]]
            values = string.join(values,'\t')+'\n'; data.write(values)
    data.close()

################# Begin Analysis from parsing files
class EnsemblInformation:
    def __init__(self, chr, gene_start, gene_stop, strand, ensembl_gene_id, ensembl_exon_id, exon_start, exon_stop, constitutive_exon, new_exon_start, new_exon_stop, new_gene_start, new_gene_stop):
        self._geneid = ensembl_gene_id; self._exonid = ensembl_exon_id; self._chr = chr
        self._genestart = gene_start; self._genestop = gene_stop
        self._exonstart = exon_start; self._exonstop = exon_stop
        self._constitutive_exon = constitutive_exon; self._strand = strand
        self._newgenestart = new_gene_start;  self._newgenestop = new_gene_stop
        self._newexonstart = new_exon_start; self._newexonstop = new_exon_stop
    def GeneID(self): return self._geneid
    def ExonID(self): return self._exonid
    def reSetExonID(self,exonid): self._exonid = exonid
    def Chr(self): return self._chr
    def Strand(self): return self._strand
    def GeneStart(self): return self._genestart
    def GeneStop(self): return self._genestop
    def ExonStart(self): return self._exonstart
    def ExonStop(self): return self._exonstop
    def NewGeneStart(self): return self._newgenestart
    def NewGeneStop(self): return self._newgenestop
    def NewExonStart(self): return self._newexonstart
    def NewExonStop(self): return self._newexonstop
    def Constitutive(self): return self._constitutive_exon
    def ConstitutiveCall(self):
        if self.Constitutive() == '1': call = 'y'
        else: call = 'n'
        return call
    def setSpliceData(self,splice_event,splice_junctions):
        self._splice_event = splice_event; self._splice_junctions = splice_junctions
    def setIntronDeletionStatus(self,del_status):
        self._del_status = del_status
    def setAssociatedSplicingEvent(self,splice_event): self._splice_event = splice_event
    def AssociatedSplicingEvent(self):
        try: return self._splice_event
        except AttributeError: return ''
    def IntronDeletionStatus(self):
        try: return self._del_status
        except AttributeError: return 'no'
    def AssociatedSplicingJunctions(self):
        try: return self._splice_junctions
        except AttributeError: return ''
    def ExonRegionID(self): return self._exon_region_id
    def ExonNumber(self): return self._exon_num
    def RegionNumber(self): return self._region_num
    def ExonRegionNumbers(self): return (self._exon_num,self._region_num)
    def ExonRegionID(self): return 'E'+str(self._exon_num)+'-'+str(self._region_num)
    def ExonRegionID2(self): return 'E'+str(self._exon_num)+'.'+str(self._region_num)
    def IntronRegionID(self): return 'I'+str(self._exon_num)+'.1'
    
    def setExonSeq(self,exon_seq): self._exon_seq = exon_seq
    def ExonSeq(self): return self._exon_seq
    def setPrevExonSeq(self,prev_exon_seq): self._prev_exon_seq = prev_exon_seq
    def PrevExonSeq(self):
        try: return self._prev_exon_seq
        except Exception: return ''
    def setNextExonSeq(self,next_exon_seq): self._next_exon_seq = next_exon_seq
    def NextExonSeq(self):
        try: return self._next_exon_seq
        except Exception: return ''
    def setPrevIntronSeq(self,prev_intron_seq): self._prev_intron_seq = prev_intron_seq
    def PrevIntronSeq(self): return self._prev_intron_seq
    def setNextIntronSeq(self,next_intron_seq): self._next_intron_seq = next_intron_seq
    def NextIntronSeq(self): return self._next_intron_seq
    def setPromoterSeq(self,promoter_seq): self._promoter_seq = promoter_seq
    def PromoterSeq(self): return self._promoter_seq
    
    def AllGeneValues(self):
        output = str(self.ExonID())
        return output
    def __repr__(self): return self.AllGeneValues()

class ExonStructureData(EnsemblInformation):
    def __init__(self, ensembl_gene_id, chr, strand, exon_start, exon_stop, constitutive_exon, ensembl_exon_id, transcriptid):
        self._transcriptid = transcriptid
        self._geneid = ensembl_gene_id; self._exonid = ensembl_exon_id; self._chr = chr
        self._exonstart = exon_start; self._exonstop = exon_stop
        self._constitutive_exon = constitutive_exon; self._strand = strand        
    def TranscriptID(self): return self._transcriptid

class ExonRegionData(EnsemblInformation):
    def __init__(self, ensembl_gene_id, chr, strand, exon_start, exon_stop, ensembl_exon_id, exon_region_id, exon_num, region_num,constitutive_exon):
        self._exon_region_id = exon_region_id; self._constitutive_exon = constitutive_exon
        self._geneid = ensembl_gene_id; self._exonid = ensembl_exon_id; self._chr = chr
        self._exonstart = exon_start; self._exonstop = exon_stop
        self._exon_num = exon_num; self._region_num = region_num; self._strand = strand        

class ProbesetAnnotation(EnsemblInformation):
    def __init__(self, ensembl_exon_id, constitutive_exon, exon_region_id, splice_event, splice_junctions,exon_start,exon_stop):
        self._region_num = exon_region_id; self._constitutive_exon = constitutive_exon;self._splice_event =splice_event
        self._splice_junctions = splice_junctions;self._exonid = ensembl_exon_id
        self._exonstart = exon_start; self._exonstop = exon_stop

class CriticalExonInfo:
    def __init__(self,geneid,critical_exon,splice_type,junctions):
        self._geneid = geneid; self._junctions = junctions  
        self._critical_exon = critical_exon; self._splice_type = splice_type
    def GeneID(self): return self._geneid
    def Junctions(self): return self._junctions
    def CriticalExonRegion(self): return self._critical_exon
    def SpliceType(self): return self._splice_type

class RelativeExonLocations:
    def __init__(self,exonid,pes,pee,nes,nee):
        self._exonid = exonid; self._pes = pes; self._pee = pee
        self._nes = nes; self._nee = nee
    def ExonID(self): return self._exonid
    def PrevExonCoor(self): return (self._pes,self._pee)
    def NextExonCoor(self): return (self._nes,self._nee)
    def __repr__(self): return self.ExonID()
    
################### Import exon coordinate/transcript data from BIOMART
def importEnsExonStructureDataSimple(species,type,gene_strand_db,exon_location_db,adjacent_exon_locations):
    if type == 'ensembl': filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    elif type == 'ucsc': filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_filtered_mrna.txt'
    elif type == 'ncRNA': filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_filtered_ncRNA.txt'
    start_time = time.time()
    fn=filepath(filename); x=0; k=[]; relative_exon_locations={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            if 'Chromosome' in t[0]:  type = 'old'###when using older builds of EnsMart versus BioMart
            else: type = 'current'
            x=1
        else:
            if type == 'old': chr, strand, gene, ens_transcriptid, ens_exonid, exon_start, exon_end, constitutive_exon = t
            else: gene, chr, strand, exon_start, exon_end, ens_exonid, constitutive_exon, ens_transcriptid = t
            ###switch exon-start and stop if in the reverse orientation
            if strand == '-1' or strand == '-': strand = '-'#; exon_start2 = int(exon_end); exon_end2 = int(exon_start); exon_start=exon_start2; exon_end=exon_end2
            else: strand = '+'; exon_end = int(exon_end)#; exon_start = int(exon_start)
            exon_end = int(exon_end); exon_start = int(exon_start)
            exon_data = (exon_start,exon_end,ens_exonid)
            try: relative_exon_locations[ens_transcriptid,gene,strand].append(exon_data)
            except KeyError: relative_exon_locations[ens_transcriptid,gene,strand] = [exon_data]
            gene_strand_db[gene] = strand
            exon_location_db[ens_exonid] = exon_start,exon_end

    ###Generate a list of exon possitions for adjacent exons for all exons
    first_exon_dbase={}
    for (transcript,gene,strand) in relative_exon_locations:
        relative_exon_locations[(transcript,gene,strand)].sort()
        if strand == '-': relative_exon_locations[(transcript,gene,strand)].reverse()
        i = 0
        ex_ls = relative_exon_locations[(transcript,gene,strand)]
        for exon_data in ex_ls:
            exonid = exon_data[-1]
            if i == 0: ### first exon
                pes = -1; pee = -1 ###Thus, Index should be out of range, but since -1 is valid, it won't be
                if strand == '-': ces = ex_ls[i][1]
                else: ces = ex_ls[i][0]
                try: first_exon_dbase[gene].append([ces,exonid])
                except KeyError: first_exon_dbase[gene] = [[ces,exonid]]
            else: pes = ex_ls[i-1][0]; pee = ex_ls[i-1][1] ###pes: previous exon start, pee: previous exon end
            try:  nes = ex_ls[i+1][0]; nee = ex_ls[i+1][1]
            except IndexError: nes = -1; nee = -1
            rel = RelativeExonLocations(exonid,pes,pee,nes,nee)
            """if exonid in adjacent_exon_locations:
                rel1 = adjacent_exon_locations[exonid]
                prev_exon_start,prev_exon_stop = rel1.NextExonCoor()
                next_exon_start,next_exon_stop = rel1.PrevExonCoor()
                if prev_exon_start == -1 or next_exon_start == -1:
                    adjacent_exon_locations[exonid] = rel ###Don't over-ride the exisitng entry if no exon is proceeding or following
            else: adjacent_exon_locations[exonid] = rel"""
            adjacent_exon_locations[exonid] = rel
            i+=1

    for gene in first_exon_dbase:
        first_exon_dbase[gene].sort()
        strand = gene_strand_db[gene]
        if strand == '-': first_exon_dbase[gene].reverse()
        first_exon_dbase[gene] = first_exon_dbase[gene][0][1] ### select the most 5' of the start exons for the gene
        
    #print relative_exon_locations['ENSMUST00000025142','ENSMUSG00000024293','-']; kill
    end_time = time.time(); time_diff = int(end_time-start_time)
    print filename,"parsed in %d seconds" % time_diff
    return gene_strand_db,exon_location_db,adjacent_exon_locations,first_exon_dbase

def importEnsExonStructureData(filename,species):
    start_time = time.time()
    fn=filepath(filename); x=0; k=[]
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            if 'Chromosome' in t[0]:  type = 'old'###when using older builds of EnsMart versus BioMart
            else: type = 'current'
            x=1
        else:
            if type == 'old': chr, strand, gene, ens_transcriptid, ens_exonid, exon_start, exon_end, constitutive_exon = t
            else: gene, chr, strand, exon_start, exon_end, ens_exonid, constitutive_exon, ens_transcriptid = t
            ###switch exon-start and stop if in the reverse orientation
            if strand == '-1' or strand == '-': strand = '-'; exon_start2 = int(exon_end); exon_end2 = int(exon_start); exon_start=exon_start2; exon_end=exon_end2
            else: strand = '+'; exon_end = int(exon_end); exon_start = int(exon_start)
            ens_exonid_data = string.split(ens_exonid,'.'); ens_exonid = ens_exonid_data[0]
            if abs(exon_end-exon_start)>-1:
                if abs(exon_end-exon_start)<1:
                    try: too_short[gene].append(exon_start)
                    except KeyError: too_short[gene] = [exon_start]
                exon_coordinates = [exon_start,exon_end]
                continue_analysis = 'no'
                if test=='yes': ###used to test the program for a single gene
                    if gene in test_gene: continue_analysis='yes' 
                else: continue_analysis='yes'
                if  continue_analysis=='yes':
                    ###Create temporary databases storing just exon and just trascript and then combine in the next block of code
                    initial_exon_annotation_db[ens_exonid] = gene,chr,strand,exon_start,exon_end,constitutive_exon                    
                    try: exon_transcript_db[ens_exonid].append(ens_transcriptid)
                    except KeyError: exon_transcript_db[ens_exonid] = [ens_transcriptid]
                    ###Use this database to figure out which ensembl exons represent intron retention as a special case down-stream
                    try: transcript_exon_db[gene,chr,strand,ens_transcriptid].append(exon_coordinates)
                    except KeyError: transcript_exon_db[gene,chr,strand,ens_transcriptid] = [exon_coordinates]
                    ###Store transcript data for downstream analyses
                    transcript_gene_db[ens_transcriptid] = gene,chr,strand
                    try: gene_transcript[gene].append(ens_transcriptid)
                    except KeyError: gene_transcript[gene] = [ens_transcriptid]
    end_time = time.time(); time_diff = int(end_time-start_time)
    print len(transcript_gene_db), "number of transcripts included"
    print filename,"parsed in %d seconds" % time_diff
    
def getEnsExonStructureData(species,data_type):
    start_time = time.time()
    ###Simple function to import and organize exon/transcript data
    filename1 = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    filename2 = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_filtered_mrna.txt'
    filename3 = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_filtered_ncRNA.txt'
    global initial_exon_annotation_db; initial_exon_annotation_db={}
    global exon_transcript_db; exon_transcript_db={}
    global transcript_gene_db; transcript_gene_db={}
    global gene_transcript; gene_transcript={}
    global transcript_exon_db; transcript_exon_db={}
    global initial_junction_db; initial_junction_db = {}; global too_short; too_short={}
    global ensembl_annotations; ensembl_annotations={}; global ensembl_gene_coordinates; ensembl_gene_coordinates = {}

    if data_type == 'mRNA':
        importEnsExonStructureData(filename2,species)
        importEnsExonStructureData(filename1,species)
    elif data_type == 'ncRNA':  ###Builds database based on a mix of Ensembl, GenBank and UID ncRNA IDs
        importEnsExonStructureData(filename3,species)
        
    global exon_annotation_db; exon_annotation_db={} ###Agglomerate the data from above and store as an instance of the ExonStructureData class
    for ens_exonid in initial_exon_annotation_db:
        gene,chr,strand,exon_start,exon_end,constitutive_exon = initial_exon_annotation_db[ens_exonid]
        ens_transcript_list = exon_transcript_db[ens_exonid]

        ###Record the coordiantes for matching up UCSC exon annotations with Ensembl genes
        try: ensembl_gene_coordinates[gene].append(exon_start)
        except KeyError: ensembl_gene_coordinates[gene] = [exon_start]
        ensembl_gene_coordinates[gene].append(exon_end)
        ensembl_annotations[gene] = chr,strand
        
        y = ExonStructureData(gene, chr, strand, exon_start, exon_end, constitutive_exon, ens_exonid, ens_transcript_list)
        exon_info = [exon_start,exon_end,y]
        if gene in too_short: too_short_list = too_short[gene]
        else: too_short_list=[]
        if exon_start not in too_short_list: ###Ensembl exons that are one bp (start and stop the same)
            try: exon_annotation_db[(gene,chr,strand)].append(exon_info)
            except KeyError: exon_annotation_db[(gene,chr,strand)] = [exon_info]

    initial_exon_annotation_db={}; exon_transcript_db={}
    print 'Exon and transcript data obtained for %d genes' % len(exon_annotation_db)
    ###Grab the junction data
    for (gene,chr,strand,transcript) in transcript_exon_db:
        exon_list_data = transcript_exon_db[(gene,chr,strand,transcript)];exon_list_data.sort()
        if strand == '-': exon_list_data.reverse()
        index=0
        if gene in too_short: too_short_list = too_short[gene]
        else: too_short_list=[]
        while index+1 <len(exon_list_data):
            junction_data = exon_list_data[index],exon_list_data[index+1]
            if junction_data[0][0] not in too_short_list and junction_data[0][1] not in too_short_list and junction_data[1][0] not in too_short_list and junction_data[1][1] not in too_short_list: ###Ensembl exons that are one bp (start and stop the same)
                try: initial_junction_db[(gene,chr,strand)].append(junction_data)
                except KeyError: initial_junction_db[(gene,chr,strand)] = [junction_data]
            index+=1
    print 'Exon-junction data obtained for %d genes' % len(initial_junction_db)
    
    ###Find exons that suggest intron retention: Within a junction database, find exons that overlapp with junction boundaries
    intron_retention_db={}; retained_intron_exons={}; delete_db={}
    for key in initial_junction_db:
        for junction_info in initial_junction_db[key]:
            e5_pos = junction_info[0][1]; e3_pos = junction_info[1][0] ###grab the junction coordiantes
            for exon_info in exon_annotation_db[key]:
                exon_start,exon_stop,ed = exon_info
                loc = [e5_pos,e3_pos]; loc.sort() ###The downstream functions need these two sorted
                new_exon_info = loc[0],loc[1],ed
                retained = compareExonLocations(e5_pos,e3_pos,exon_start,exon_stop)
                exonid = ed.ExonID()
                if retained == 'yes':
                    #print e5_pos,e3_pos,exon_start,exon_stop;kill
                    intron_length = abs(e5_pos-e3_pos)
                    if intron_length>500:
                        ed.setIntronDeletionStatus('yes')
                        try: delete_db[key].append(exon_info)
                        except KeyError: delete_db[key] = [exon_info]
                    try: intron_retention_db[key].append(new_exon_info)
                    except KeyError: intron_retention_db[key]=[new_exon_info]
                    retained_intron_exons[exonid]=[]
                    #print key,ed.ExonID(),len(exon_annotation_db[key])

    k=0
    exon_annotation_db2={}
    for key in exon_annotation_db:
        for exon_info in exon_annotation_db[key]:
            if key in delete_db:
                delete_info = delete_db[key] ### coordinates and objects... looks like you can match up based on object memory locations
                if exon_info not in delete_info:
                    try: exon_annotation_db2[key].append(exon_info)
                    except KeyError: exon_annotation_db2[key]=[exon_info]
                else: k+=1
            else:
                try: exon_annotation_db2[key].append(exon_info)
                except KeyError: exon_annotation_db2[key]=[exon_info]
                
    exon_annotation_db = exon_annotation_db2

    print k, 'exon entries removed from primary exon structure, which occur in predicted retained introns'    
    initial_junction_db={}
    print len(retained_intron_exons),"ensembl exons in %d genes, show evidence of being retained introns (only sequences > 500bp are removed from the database" % len(intron_retention_db)

    end_time = time.time(); time_diff = int(end_time-start_time)
    print "Primary databases build in %d seconds" % time_diff

    try: ucsc_splicing_annot_db = alignToKnownAlt.importEnsExonStructureData(species,ensembl_gene_coordinates,ensembl_annotations,exon_annotation_db) ### Should be able to exclude
    except Exception: ucsc_splicing_annot_db={}
    return exon_annotation_db,transcript_gene_db,gene_transcript,transcript_exon_db,intron_retention_db,ucsc_splicing_annot_db

def compareExonLocations(e5_pos,e3_pos,exon_start,exon_stop):
    sort_list = [e5_pos,e3_pos,exon_start,exon_stop]; sort_list.sort()
    new_sort = sort_list[1:-1]
    if e5_pos in new_sort and e3_pos in new_sort: retained = 'yes'
    else: retained = 'no'
    return retained
    
################### Import exon sequence data from BIOMART (more flexible alternative function to above)        
def import_sequence_data(filename,filter_db,species,analysis_type):
    print "Begining generic fasta import of",filename
    fn=filepath(filename);fasta = {}; exon_db = {}; gene_db = {}; cDNA_db = {};sequence = '';count = 0
    global temp_seq; temp_seq=''; damned =0; global failed; failed={}
    addition_seq_len = 2000; var = 1000
    if 'gene' in fn:
        gene_strand_db,exon_location_db,adjacent_exon_locations,null = importEnsExonStructureDataSimple(species,'ucsc',{},{},{})
        gene_strand_db,exon_location_db,adjacent_exon_locations,first_exon_db = importEnsExonStructureDataSimple(species,'ensembl',gene_strand_db,exon_location_db,adjacent_exon_locations)
        null=[]
    for line in open(fn,'r').xreadlines():
        data = cleanUpLine(line)
        try:
            if data[0] == '>':
                    if len(sequence) > 0:
                        if 'gene' in fn:
                            start = int(start); stop = int(stop)
                        else:
                            try: exon_start = int(exon_start); exon_stop = int(exon_stop)
                            except ValueError: exon_start = exon_start
                            if strand == '-1': strand = '-'
                            else: strand = '+'
                            if strand == '-': new_exon_start = exon_stop; new_exon_stop = exon_start; exon_start = new_exon_start; exon_stop = new_exon_stop #"""
                        if 'exon' in fn:
                            exon_info = [exon_start,exon_stop,exon_id,exon_annot]
                            try: exon_db[(gene,chr,strand)].append(exon_info)
                            except KeyError: exon_db[(gene,chr,strand)] = [exon_info] #exon_info = (exon_start,exon_stop,exon_id,exon_annot)
                            fasta[geneid]=description
                        if 'cDNA' in fn:
                            cDNA_info = [transid,sequence]
                            try: cDNA_db[(gene,strand)].append(cDNA_info)
                            except KeyError: cDNA_db[(gene,strand)] = [cDNA_info]                   
                        if 'gene' in fn:
                            temp_seq = sequence
                            if gene in filter_db:
                                count += 1
                                if count == var: print var,'genes examined...'; var+=1000
                                #print gene, filter_db[gene][0].ArrayGeneID();kill
                                if (len(sequence) -(stop-start)) != ((addition_seq_len*2) +1):
                                    ###multiple issues can occur with trying to grab sequence up and downstream - for now, exlcude these
                                    ###some can be accounted for by being at the end of a chromosome, but not all.
                                    damned +=1
                                    try: failed[chr]+=1
                                    except KeyError: failed[chr]=1
                                    """if chr in failed:
                                        gene_len = (stop-start); new_start = start - addition_seq_len
                                        null = sequence[(gene_len+2000):]"""
                                else:
                                    original_start = start; original_stop = stop
                                    start = start - addition_seq_len
                                    stop = stop + addition_seq_len
                                    strand = gene_strand_db[gene]
                                    first_exonid = first_exon_db[gene]
                                    fexon_start,fexon_stop = exon_location_db[first_exonid]
                                    for ed in filter_db[gene]:
                                        if analysis_type == 'get_locations':
                                            cd = seqSearch(sequence,ed.ExonSeq())
                                            if cd != -1:
                                                if strand == '-': exon_stop = stop - cd; exon_start = exon_stop - len(ed.ExonSeq()) + 1
                                                else: exon_start = start + cd; exon_stop = exon_start + len(ed.ExonSeq())-1
                                                ed.setExonStart(exon_start); ed.setExonStop(exon_stop); ed.setGeneStart(original_start); ed.setGeneStop(original_stop)
                                                #if ed.ExonID() == 'E10' and ed.ArrayGeneID() == 'G7225860':
                                                #print exon_start, exon_stop,len(ed.ExonSeq()),ed.ExonSeq();kill
                                            else:
                                                cd = seqSearch(sequence,ed.ExonSeq()[:15])
                                                #print ed.ExonSeq()[:15],ed.ExonSeq();kill
                                                if cd == -1: cd = seqSearch(sequence,ed.ExonSeq()[-15:])
                                                if cd != -1:
                                                    if strand == '-': exon_stop = stop - cd; exon_start = exon_stop - len(ed.ExonSeq()) + 1
                                                    else: exon_start = start + cd; exon_stop = exon_start + len(ed.ExonSeq())-1
                                                    ed.setExonStart(exon_start); ed.setExonStop(exon_stop); ed.setGeneStart(original_start); ed.setGeneStop(original_stop)
                                                else: null=[]#print exon_start, exon_stop, ed.ExonSeq();kill
                                        if analysis_type == 'get_sequence':
                                            exon_id,((probe_start,probe_stop,probeset_id,exon_class,transcript_clust),ed) = ed
                                            ens_exon_list = ed.ExonID()
                                            for ens_exon in ens_exon_list:
                                                if len(ens_exon)>0:
                                                    exon_start,exon_stop = exon_location_db[ens_exon]
                                                    exon_sequence = grabSeq(sequence,strand,start,stop,exon_start,exon_stop,'exon')
                                                    """Could repeat if we build another dictionary with exon->adjacent exon positions (store in a class where
                                                    you designate last and next exon posiitons for each transcript relative to that exon), to grab downstream, upsteam exon and intron sequences"""
                                                    try:
                                                        rel = adjacent_exon_locations[ens_exon]
                                                        prev_exon_start,prev_exon_stop = rel.PrevExonCoor()
                                                        next_exon_start,next_exon_stop = rel.NextExonCoor()
                                                        prev_exon_sequence = grabSeq(sequence,strand,start,stop,prev_exon_start,prev_exon_stop,'exon')
                                                        next_exon_sequence = grabSeq(sequence,strand,start,stop,next_exon_start,next_exon_stop,'exon')
                                                        seq_type = 'intron'
                                                        if strand == '-':
                                                            if 'alt-N-term' in ed.AssociatedSplicingEvent() or 'altPromoter' in ed.AssociatedSplicingEvent(): seq_type = 'promoter' ### Thus prev_intron_seq is used to designate an alternative promoter sequence
                                                            prev_intron_sequence = grabSeq(sequence,strand,start,stop,exon_stop,prev_exon_start,seq_type)
                                                            promoter_sequence = grabSeq(sequence,strand,start,stop,fexon_stop,-1,"promoter")
                                                            next_intron_sequence = grabSeq(sequence,strand,start,stop,next_exon_stop,exon_start,'intron')
                                                        else:
                                                            prev_intron_sequence = grabSeq(sequence,strand,start,stop,prev_exon_stop,exon_start,seq_type)
                                                            promoter_sequence = grabSeq(sequence,strand,start,stop,-1,fexon_start,"promoter")
                                                            next_intron_sequence = grabSeq(sequence,strand,start,stop,exon_stop,next_exon_start,'intron')
                                                        """if 'ENS' in ens_exon:
                                                            print ens_exon, strand
                                                            print '1)',exon_sequence
                                                            print '2)',prev_intron_sequence[:20],prev_intron_sequence[-20:], len(prev_intron_sequence), strand,start,stop,prev_exon_stop,exon_start,seq_type, ed.AssociatedSplicingEvent()
                                                            print '3)',next_intron_sequence[:20],next_intron_sequence[-20:], len(next_intron_sequence)
                                                            print '4)',promoter_sequence[:20],promoter_sequence[-20:], len(promoter_sequence);kill"""
                                                        ###Intron sequences can be extreemly long so just include the first and last 1kb
                                                        if len(prev_intron_sequence)>2001 and seq_type == 'promoter': prev_intron_sequence = prev_intron_sequence[-2001:]
                                                        elif len(prev_intron_sequence)>2001: prev_intron_sequence = prev_intron_sequence[:1000]+'|'+prev_intron_sequence[-1000:]
                                                        if len(next_intron_sequence)>2001: next_intron_sequence = next_intron_sequence[:1000]+'|'+next_intron_sequence[-1000:]
                                                        if len(promoter_sequence)>2001: promoter_sequence = promoter_sequence[-2001:]
                                                    except KeyError:
                                                        prev_exon_sequence=''; next_intron_sequence=''; exon_sequence=''
                                                        if strand == '-': promoter_sequence = grabSeq(sequence,strand,start,stop,fexon_stop,-1,"promoter")
                                                        else: promoter_sequence = grabSeq(sequence,strand,start,stop,-1,fexon_start,"promoter")
                                                        if len(promoter_sequence)>2001: promoter_sequence = promoter_sequence[-2001:]                                                     
                                                    ed.setExonSeq(exon_sequence) ### Use to replace the previous probeset/critical exon sequence with sequence corresponding to the full exon
                                                    if len(prev_exon_sequence)>0:
                                                        ### Use to output sequence for ESE/ISE type motif searches
                                                        ed.setPrevExonSeq(prev_exon_sequence);
                                                    if len(next_exon_sequence)>0: ed.setNextExonSeq(next_exon_sequence)
                                                    ed.setPrevIntronSeq(prev_intron_sequence[1:-1]); ed.setNextIntronSeq(next_intron_sequence[1:-1])
                                                    ed.setPromoterSeq(promoter_sequence[1:-1])
                                                else:
                                                    if strand == '-': promoter_sequence = grabSeq(sequence,strand,start,stop,fexon_stop,-1,"promoter")
                                                    else: promoter_sequence = grabSeq(sequence,strand,start,stop,-1,fexon_start,"promoter")
                                                    if len(promoter_sequence)>2001: promoter_sequence = promoter_sequence[-2001:]
                                                    ed.setPromoterSeq(promoter_sequence[1:-1])
                            sequence = ''; data2 = data[1:]; t= string.split(data2,'|'); gene,chr,start,stop = t
                    else: data2 = data[1:]; t= string.split(data2,'|'); gene,chr,start,stop = t
        except IndexError: continue
        try:
            if data[0] != '>': sequence = sequence + data
        except IndexError: continue

    print "Number of imported sequences:", len(fasta),count
    if len(exon_db) > 0: return exon_db,fasta
    elif len(cDNA_db) > 0: return cDNA_db
    elif len(fasta) > 0: return fasta
    else: return filter_db

def grabSeq(sequence,strand,start,stop,exon_start,exon_stop,type):
    proceed = 'yes'
    if type != 'promoter' and (exon_start == -1 or exon_stop == -1): proceed = 'no' ###Thus no preceeding or succedding exons and thus no seq reported
    if proceed == 'yes':
        if strand == '-':
            if exon_stop == -1: exon_stop = stop
            exon_sequence = sequence[(stop-exon_stop):(stop-exon_stop)+(exon_stop-exon_start)+1]
        else:
            #print type, exon_start,start,exon_stop
            if exon_start == -1: exon_start = start ### For last intron 
            exon_sequence = sequence[(exon_start-start):(exon_stop-start+1)]
    else: exon_sequence = ''
    return exon_sequence

def seqSearch(sequence,exon_seq):
    cd = string.find(sequence,exon_seq)
    if cd == -1:
        rev_seq = reverse_orientation(exon_seq); cd = string.find(sequence,rev_seq)
    return cd

def reverse_string(astring):
    "http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/65225"
    revchars = list(astring)        # string -> list of chars
    revchars.reverse()              # inplace reverse the list
    revchars = ''.join(revchars)    # list of strings -> string
    return revchars

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

############# First pass for annotating exons into destict, ordered regions for further annotation
def annotate_exons(exon_location):
    print "Begining to assign intial exon block and region annotations"
    ### Sort and reverse exon orientation for transcript_cluster exons
    original_neg_strand_coord={}
    ###make negative strand coordinates look like positive strand to identify overlapping exons
    for (geneid,chr,strand) in exon_location:
        exon_location[(geneid,chr,strand)].sort()
        if strand == '-':
            exon_location[(geneid,chr,strand)].reverse()
            denominator = exon_location[(geneid,chr,strand)][0][0] ###numerical coordiantes to subtract from to normalize negative strand data
            for exon_info in exon_location[(geneid,chr,strand)]:
                start,stop,ed = exon_info; ens_exon = ed.ExonID()
                coordinates = stop,start; coordinates = copy.deepcopy(coordinates)###format these in reverse for analysis in other modules
                original_neg_strand_coord[ens_exon] = coordinates
                exon_info[0] = abs(start - denominator);exon_info[1] = abs(stop - denominator)
 
    #alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','x','y','z']
    exon_location2={}; exon_temp_list=[]
    for key in exon_location:
        index = 1; index2 = 1; exon_list=[]; y = 0 #if key[-1] == '-':
        for exon in exon_location[key]: #print exon[0],exon[1],len(exon_temp_list)
            if y == 0:
                exon_info = ['E'+str(index)+'-1',exon,(index,1)]
                exon_list.append(exon_info); y = 1; last_start = exon[0]; last_stop = exon[1]; index += 1; index2 = 2; exon_temp_list =[]; exon_temp_list.append(last_start); exon_temp_list.append(last_stop)
            elif y == 1:
                current_start = exon[0];  current_stop = exon[1]
                if ((current_start >= last_start) and (current_start <= last_stop)) or ((last_start >= current_start) and (last_start <= current_stop)):
                    exon_info = ['E'+str(index-1) +'-'+ str(index2),exon,(index-1,index2)] #+alphabet[index2]
                    exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index2 += 1; exon_temp_list.append(current_start); exon_temp_list.append(current_stop)
                elif (abs(current_start - last_stop) < 1) or (abs(last_start - current_stop) < 1):
                    exon_info = ['E'+str(index-1) +'-'+ str(index2),exon,(index-1,index2)] #+alphabet[index2]
                    exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index2 += 1; exon_temp_list.append(current_start); exon_temp_list.append(current_stop)
                elif len(exon_temp_list)>3:
                    exon_temp_list.sort()
                    if (((current_start-1) > exon_temp_list[-1]) and ((current_stop-1) > exon_temp_list[-1])) or (((current_start+1) < exon_temp_list[0]) and ((current_stop+1) < exon_temp_list[0])):
                        ###Thus an overlapp with atleast one exon DOESN'T exists
                        exon_info = ['E'+str(index)+'-1',exon,(index,1)]
                        exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index += 1; index2 = 2; exon_temp_list=[]; exon_temp_list.append(current_start); exon_temp_list.append(current_stop)
                    else:
                        exon_info = ['E'+str(index-1) +'-'+ str(index2),exon,(index-1,index2)]  #+alphabet[index2]
                        exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index2 += 1; exon_temp_list.append(current_start); exon_temp_list.append(current_stop)
                else:
                    exon_info = ['E'+str(index)+'-1',exon,(index,1)]
                    exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index += 1; index2 = 2; exon_temp_list=[]; exon_temp_list.append(current_start); exon_temp_list.append(current_stop)
        exon_location2[key] = exon_list

    for key in exon_location2:
        ###Re-assign exon coordiantes back to the actual coordiantes for negative strand genes
        strand = key[-1]
        if strand == '-':
            for exon_data in exon_location2[key]:
                ed = exon_data[1][2]; ens_exon = ed.ExonID()
                ###re-assing the coordiantes
                start,stop = original_neg_strand_coord[ens_exon]
                exon_data[1][0] = start; exon_data[1][1] = stop
    """
    for key in exon_location2:
        if key[0] == 'ENSG00000129566':
            for item in exon_location2[key]: print key, item"""
    return exon_location2

def exon_clustering(exon_location):
    """for i in exon_location[('ENSG00000129566','14','-')]:print i"""
    #exon_info = exon_start,exon_end,y
    #try: exon_annotation_db[(geneid,chr,strand)].append(exon_info)
    #except KeyError: exon_annotation_db[(geneid,chr,strand)] = [exon_info]
    
    exon_clusters={}; region_locations={}; region_gene_db={}
    for key in exon_location:
        chr = 'chr'+key[1];entries = exon_location[key];strand = key[-1]; gene = key[0]
        temp_exon_db={}; temp_exon_db2={}; exon_block_db={}
        for exon in entries:
            a = string.find(exon[0],'-')  ###outdated code: if all have a '-'
            if a == -1: exon_number = int(exon[0][1:])
            else: exon_number = int(exon[0][1:a])
            ###Group overlapping exons into exon clusters
            try: temp_exon_db[exon_number].append(exon[1])
            except KeyError: temp_exon_db[exon_number] = [exon[1]]
            ###add all start and stop values to the temp database
            try: temp_exon_db2[exon_number].append(exon[1][0])
            except KeyError: temp_exon_db2[exon_number] = [exon[1][0]]
            temp_exon_db2[exon_number].append(exon[1][1])
            try: exon_block_db[exon_number].append(exon)
            except KeyError: exon_block_db[exon_number] = [exon]
        for exon in temp_exon_db:
            #intron = 'I'+str(exon)+'-1'
            exon_info = unique.unique(temp_exon_db2[exon])
            exon_info.sort();start = exon_info[0];stop = exon_info[-1];type=[]; accession=[]
            for (exon_start,exon_stop,ed) in temp_exon_db[exon]:
                exon_type = ed.Constitutive()
                exon_id = ed.ExonID()
                type.append(exon_type); accession.append(exon_id)
                #if exon_id == 'ENSE00000888906': print exon_info,temp_exon_db[exon]
            type=unique.unique(type); accession=unique.unique(accession)
            exon_data = exon,(start,stop),accession,type
            key1 = key[0],chr,strand
            try: exon_clusters[key1].append(exon_data)
            except KeyError: exon_clusters[key1] = [exon_data]
            #if len(exon_info)-len(temp_exon_db[exon])>2: print key1,len(exon_info),len(temp_exon_db[exon]);kill
            if len(exon_info)>2:
                if strand == '-': exon_info.reverse()
                index=0; exon_data_list=[]
                while index < (len(exon_info)-1):
                    region = str(index+1)#;region_locations[key,'E'+str(exon)+'-'+region] = exon_info[index:index+2]
                    if strand == '-': new_stop,new_start = exon_info[index:index+2]
                    else: new_start,new_stop = exon_info[index:index+2]
                    ned = ExonStructureData(gene, key[1], strand, new_start, new_stop, '', '', [])
                    new_exon_info = ['E'+str(exon)+'-'+region,[new_start,new_stop,ned],(exon,index+1)]
                    exon_data_list.append(new_exon_info)
                    index+=1
                    #if gene == 'ENSG00000171735':print new_exon_info, 0
                region_locations[key,exon] = exon_data_list
            else:
                exon_data_list = [exon_block_db[exon][0]]  ###There can be multiples that occur - 2 exons 1 set of coordinates
                #if gene == 'ENSG00000171735':print exon_data_list, 1
                region_locations[key,exon] = exon_data_list
                
    ###Resort and re-populated the new exon_location database where we've re-defined the region entries    
    interim_location_db={}
    for (key,exon) in region_locations:
        exon_data_list = region_locations[(key,exon)]
        for exon_data in exon_data_list:
            try: interim_location_db[key].append((exon,exon_data))
            except KeyError: interim_location_db[key] = [(exon,exon_data)]
    for key in interim_location_db:
        interim_location_db[key].sort(); new_exon_list=[]
        for (e,i) in interim_location_db[key]: new_exon_list.append(i)
        exon_location[key] = new_exon_list

    #for i in exon_location[('ENSG00000171735', '1', '+')]:print i

    """
    for i in region_locations:
        if 'ENSG00000129566' in i[0]:print i, region_locations[i]
    
    ###Transform coordinates from the source Ensembl exon to discrete regions (regions may not be biological relevant in the 3' or 5' exon of the gene due to EST assemblies).
    for (key,exon_id) in region_locations:
        gene,chr,strand = key; id_added='no'; exon_annot=[]; ens_exonids=[]; ens_transcripts=[]
        if strand == '-': stop,start = region_locations[(key,exon_id)]
        else: start,stop = region_locations[(key,exon_id)]
        ###If the number of old regions is greater than the new, delete the old
        previous_region_number = len(exon_location[key])
        new_region_number = region_gene_db[key]
        if previous_region_number>new_region_number: 
        for exon_data in exon_location[key]:
            exon_start,exon_stop,ed = exon_data[1]; ens_exon_id = ed.ExonID(); exon_annotation = ed.Constitutive()
            #if exon_id == exon_data[0]: exon_data[1][0] = start; exon_data[1][1] = stop; id_added = 'yes'
            if exon_start == start or exon_stop == stop: exon_annot.append(exon_annotation); ens_exonids.append(ens_exon_id); ens_transcripts+=ed.TranscriptID()
            if exon_stop == start or exon_start == stop: exon_annot.append(exon_annotation); ens_exonids.append(ens_exon_id)
        exon_annot = unique.unique(exon_annot); ens_exonids = unique.unique(ens_exonids); ens_transcripts = unique.unique(ens_transcripts)
        exon_annot = string.join(exon_annot,'|'); ens_exonids = string.join(ens_exonids,'|')
        for exon_data in exon_location[key]:
            if exon_id == exon_data[0]:
                ###Replace exsting entries (we're basically replacing these with completely new data, just like below, but this is easier than deleting the existing)
                y = ExonStructureData(gene, chr, strand, start, stop, exon_annot, ens_exonids, ens_transcripts)
                exon_data[1][0] = start; exon_data[1][1] = stop; exon_data[1][2] = y; id_added = 'yes'
                #start is the lower number, with either strand
                break
        if id_added == 'no': ###This occurs when a new region must be introduced from a large now broken large exon (e.g. E1-1 pos: 1-300, E1-2 pos: 20-200, now need a 3rd E1-3 200-300)
            indeces = string.split(exon_id[1:],'-')
            index1 = int(indeces[0]); index2 = int(indeces[1])
            #new_entry = ['E'+str(index-1) +'-'+ str(index2),exon,(index-1,index2)]
            #exon_info = [exon_start,exon_stop,exon_id,exon_annot]
            ###Can include inappopriate exon IDs and annotations, but not worth specializing the code
            y = ExonStructureData(gene, chr, strand, start, stop, exon_annot, ens_exonids, ens_transcripts)
            exon_info = [start,stop,y]
            new_entry = [exon_id,exon_info,(index1,index2)]
            #if key == ('ENSG00000129566', '14', '-'): print key,new_entry;kill
            exon_location[key].append(new_entry)"""
                
    exon_regions = {}
    for (gene,chr,strand) in exon_location:
        for exon_data in exon_location[(gene,chr,strand)]:
            try: exon_region_id,exon_detailed,(exon_num,region_num) = exon_data
            except ValueError: print exon_data;kill
            start,stop,ed = exon_detailed
            if strand == '+': start,stop,ed = exon_detailed
            else: stop,start,ed = exon_detailed
            y = ExonRegionData(gene, chr, strand, start, stop, ed.ExonID(), exon_region_id, exon_num, region_num, ed.Constitutive())
            try: exon_regions[gene].append(y)
            except KeyError: exon_regions[gene] = [y]
    """
    for key in exon_location:
        if key[0] == 'ENSG00000075413':
            print key
            for item in exon_location[key]: print item"""

    ###Create a corresponding database of intron and locations for the clustered block exon database
    intron_clusters={}; intron_region_db={}
    for key in exon_clusters:
        strand = key[-1]; gene = key[0]
        exon_clusters[key].sort(); index = 0
        for exon_data in exon_clusters[key]:
            try: 
                exon_num,(start,stop),null,null = exon_data
                next_exon_data = exon_clusters[key][index+1]
                en,(st,sp),null,null = next_exon_data
                intron_num = exon_num
                if strand == '+': intron_start = stop; intron_stop = st; intron_start_stop = (intron_start,intron_stop)
                else: intron_start = start; intron_stop = sp; intron_start_stop = (intron_stop,intron_start)
                index+=1
                intron_data = intron_num,intron_start_stop,'','no'
                try: intron_clusters[key].append(intron_data)
                except KeyError: intron_clusters[key] = [intron_data]
                ###This database is used for SubGeneViewer and is analagous to region_db
                intron_region_id = 'I'+str(exon)+'-1'
                rd = ExonRegionData(gene, chr, strand, intron_start, intron_stop, ed.ExonID(), intron_region_id, intron_num, 1, 0)
                if gene in intron_region_db:
                    block_db = intron_region_db[gene]
                    block_db[intron_num] = [rd]
                else:
                    block_db={}; block_db[intron_num] = [rd]
                    intron_region_db[gene] = block_db
            except IndexError: continue
            
    return exon_clusters,intron_clusters,exon_regions,intron_region_db

def eliminate_redundant_dict_values(database):
    for key in database:
        list = makeUnique(database[key])
        list.sort()
        database[key] = list
    return database

def getEnsemblAnnotations(filename,rna_processing_ensembl):
    fn=filepath(filename)
    ensembl_annotation_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if 'Description' in data: ensembl_gene_id,description,symbol = string.split(data,'\t')
        else: ensembl_gene_id,symbol,description = string.split(data,'\t')
        ensembl_annotation_db[ensembl_gene_id] = description, symbol

    for gene in ensembl_annotation_db:
        if gene in rna_processing_ensembl: mRNA_processing = 'RNA_processing/binding'
        else: mRNA_processing = ''
        index = ensembl_annotation_db[gene] 
        ensembl_annotation_db[gene] = index[0],index[1],mRNA_processing
        
    exportEnsemblAnnotations(ensembl_annotation_db)
    return ensembl_annotation_db

def exportEnsemblAnnotations(ensembl_annotation_db):
    exon_annotation_export = 'AltDatabase/ensembl/' +species+'/'+species+ '_Ensembl-annotations.txt'
    fn=filepath(exon_annotation_export); data = open(fn,'w')
    for ensembl_gene in ensembl_annotation_db:
        a = ensembl_annotation_db[ensembl_gene]
        values = ensembl_gene +'\t'+ a[0] +'\t'+ a[1] +'\t'+ a[2] + '\n'
        data.write(values)
    data.close()

def reimportEnsemblAnnotations(species):
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations.txt'
    fn=filepath(filename)
    ensembl_annotation_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        ensembl_gene_id,symbol,description,mRNA_processing = string.split(data,'\t')
        ensembl_annotation_db[ensembl_gene_id] = symbol,description,mRNA_processing
    return ensembl_annotation_db
    
def importPreviousBuildTranslation(filename,use_class_structures):
    ###When previous genome build coordinates a provided for an array - create a translation
    ###file between ensembl exon ids
    fn=filepath(filename)
    exon_coord_translation_db = {}; gene_coor_translation_db = {}; exon_db = {}
    gene_coor_translation_list = []; x=0
    for line in open(fn,'r').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x+=1 ###ignore the first line
        else:
            chr, gene_start, gene_stop, strand, ensembl_gene_id, ensembl_exon_id, exon_start, exon_stop, null, null, null,null, constitutive_exon = string.split(data,'\t')
            gene_start = int(gene_start); gene_stop = int(gene_stop)
            exon_start = int(exon_start); exon_stop = int(exon_stop)
            if strand == '-1': strand = '-'
            if strand == '1': strand = '+'
            if strand == '-': new_gene_start = gene_stop; new_gene_stop = gene_start; new_exon_start = exon_stop; new_exon_stop = exon_start 
            else: new_gene_start = gene_start; new_gene_stop = gene_stop; new_exon_start = exon_start; new_exon_stop = exon_stop 
            new_exon_start = abs(new_gene_start - new_exon_start)
            new_exon_stop = abs(new_gene_start - new_exon_stop)
            ###constitutive_exon column contains a 0 or 1: ensembl_exon_id is versioned
            ensembl_exon_id,null = string.split(ensembl_exon_id,'.')

            if use_class_structures == 'yes':
                exon_coord_info = ensembl_gene_id,(exon_start,exon_stop),(gene_start,gene_stop),int(constitutive_exon)
                ei = EnsemblInformation(chr, gene_start, gene_stop, strand, ensembl_gene_id, ensembl_exon_id, exon_start, exon_stop, constitutive_exon, new_exon_start, new_exon_stop, new_gene_start, new_gene_stop)    
                ###Also independently determine exon clusters for previous build exon data
                exon_annot=''; exon_info = (exon_start,exon_stop,ensembl_exon_id,exon_annot)
                try: exon_db[(ensembl_gene_id,chr,strand)].append(exon_info)
                except KeyError: exon_db[(ensembl_gene_id,chr,strand)] = [exon_info]
                y = ei
            else:
                ei = [chr,(gene_start,gene_stop),ensembl_gene_id,(new_gene_start,new_gene_stop)]
                y = [chr,(exon_start,exon_stop),ensembl_exon_id,(new_exon_start,new_exon_stop)]
            try: exon_coord_translation_db[ensembl_gene_id].append(y)
            except KeyError: exon_coord_translation_db[ensembl_gene_id] = [y]
            gene_coor_translation_db[(chr,strand),gene_start,gene_stop] = ei
    for key in gene_coor_translation_db:
        ei = gene_coor_translation_db[key]
        gene_coor_translation_list.append([key,ei])
    gene_coor_translation_list.sort()
    gene_coor_translation_list2={}
    for key in gene_coor_translation_list:
        chr_strand = key[0][0]
        try: gene_coor_translation_list2[chr_strand].append(key[-1])
        except KeyError: gene_coor_translation_list2[chr_strand] = [key[-1]]
    gene_coor_translation_list = gene_coor_translation_list2
    ###Determine Exon Clusters for current Ensembl genes with poor linkage properties (can't be converted to old coordiantes) 
    exon_db2 = annotate_exons(exon_db)
    exon_clusters = exon_clustering(exon_db2); exon_db2={}
    return exon_coord_translation_db, exon_db, exon_clusters

def importExonTranscriptAnnotations(filename):
    fn=filepath(filename)
    exon_trans_association_db = {}; x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x+=1
        else:
            ensembl_gene_id,ensembl_trans_id,ensembl_exon_id,constitutive = string.split(data,'\t')
            ensembl_exon_id,null = string.split(ensembl_exon_id,'.')
            try: exon_trans_association_db[ensembl_exon_id].append([ensembl_trans_id,constitutive])
            except KeyError: exon_trans_association_db[ensembl_exon_id] = [[ensembl_trans_id,constitutive]]
    return exon_trans_association_db 

def importEnsemblDomainData(filename):
    fn=filepath(filename); x = 0; ensembl_ft_db = {}; ensembl_ft_summary_db = {} # Use the last database for summary statistics
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 1:
            try: ensembl_gene, chr, mgi, uniprot, ensembl_prot, seq_data, position_info = string.split(data,'\t')
            except ValueError: continue
            ft_info_list = string.split(position_info,' | ')
            for entry in ft_info_list:
                try: peptide_start_end, gene_start_end, feature_source, interprot_id, description = string.split(entry,' ')
                except ValueError: continue
                ###142-180 3015022-3015156 Pfam IPR002050 Env_polyprotein
                ft_start_pos, ft_end_pos = string.split(peptide_start_end,'-')
                pos1 = int(ft_start_pos); pos2 = int(ft_end_pos)
                
                #sequence_fragment = seq_data[pos1:pos2]
                if len(description)>1 or len(interprot_id)>1:
                    ft_info = [description,sequence_fragment,interprot_id]
                    ft_info2 = description,interprot_id
                    ###uniprot_ft_db[id].append([ft,pos1,pos2,annotation])
                    try: ensembl_ft_db[id].append(ft_info)
                    except KeyError: ensembl_ft_db[id] = [ft_info]
                    try: ensembl_ft_summary_db[id].append(ft_info2)
                    except KeyError: ensembl_ft_summary_db[id] = [ft_info2]         
        elif data[0:6] == 'GeneID': x = 1
        
    ensembl_ft_db = eliminate_redundant_dict_values(ensembl_ft_db)
    ensembl_ft_summary_db = eliminate_redundant_dict_values(ensembl_ft_summary_db)
    domain_gene_counts = {}
    ###Count the number of domains present in all genes (count a domain only once per gene)
    for gene in ensembl_ft_summary_db:
        for domain_info in ensembl_ft_summary_db[gene]:
            try: domain_gene_counts[domain_info] += 1
            except KeyError: domain_gene_counts[domain_info] = 1
    print "Number of Ensembl genes, linked to array genes with domain annotations:",len(ensembl_ft_db)
    print "Number of Ensembl domains:",len(domain_gene_counts)
    return ensembl_ft_db,domain_gene_counts
       
def getEnsemblAssociations(Species,data_type,test_status):
    global species; species = Species
    global test; test = test_status
    global test_gene
    meta_test = ["ENSG00000215305","ENSG00000179676","ENSG00000170484","ENSG00000138180","ENSG00000100258","ENSG00000132170","ENSG00000105767","ENSG00000105865","ENSG00000108523","ENSG00000150045","ENSG00000156026"]
    test_gene = ['ENSG00000215305']
    #test_gene = meta_test
    exon_annotation_db,transcript_gene_db,gene_transcript,transcript_exon_db,intron_retention_db,ucsc_splicing_annot_db = getEnsExonStructureData(species,data_type)
    exon_annotation_db2 = annotate_exons(exon_annotation_db); ensembl_descriptions={}
    
    exon_db = customDBDeepCopy(exon_annotation_db2) ##having problems with re-writting contents of this db when I don't want to
    exon_clusters,intron_clusters,exon_regions,intron_region_db = exon_clustering(exon_db); exon_db={}
    exon_junction_db,putative_as_junction_db,exon_junction_db = processEnsExonStructureData(exon_annotation_db,exon_regions,transcript_gene_db,gene_transcript,transcript_exon_db,intron_retention_db)
    exon_regions,critical_gene_junction_db = compareJunctions(putative_as_junction_db,exon_regions)
    exportSubGeneViewerData(exon_regions,critical_gene_junction_db,intron_region_db,intron_retention_db)
    
    ###Grab rna_processing Ensembl associations
    use_exon_data='no';get_splicing_factors = 'yes'
    rna_processing_ensembl = GO_parsing.parseAffyGO(use_exon_data,get_splicing_factors,species)
    
    ensembl_annot_file = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations_simple.txt'
    ensembl_annotation_db = getEnsemblAnnotations(ensembl_annot_file,rna_processing_ensembl)
    #exportExonClusters(exon_clusters)
    return exon_annotation_db2,ensembl_annotation_db,exon_clusters,intron_clusters,exon_regions,intron_retention_db,ucsc_splicing_annot_db,transcript_gene_db

def getExonTranscriptDomainAssociations(Species):
    global species; species = Species
    import_dir = '/AltDatabase/ensembl/'+species
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for file_name in dir_list:    #loop through each file in the directory to output results
        dir_file = 'AltDatabase/ensembl/'+species+'/'+file_name
        if 'Exon_cDNA' in dir_file: exon_trans_file = dir_file
        elif 'Domain' in dir_file: domain_file = dir_file
    exon_trans_association_db = importExonTranscriptAnnotations(exon_trans_file)
    return exon_trans_association_db
    
def exportExonClusters(exon_clusters,species):
    exon_cluster_export = 'AltDatabase/ensembl/'+species+'/'+species+'-Ensembl-exon-clusters.txt'
    fn=filepath(exon_cluster_export); data = open(fn,'w')
    for key in exon_clusters:
        ensembl_gene = key[0]
        chr = key[1]
        strand = key[2]
        for entry in exon_clusters[key]:
            exon_cluster_num = str(entry[0])
            exon_start = str(entry[1][0])
            exon_stop = str(entry[1][1])
            exon_id_list = string.join(entry[2],'|')
            annotation_list = string.join(entry[3],'|')
            values = ensembl_gene +'\t'+ chr +'\t'+ strand +'\t'+ exon_cluster_num +'\t'+ exon_start +'\t'+ exon_stop +'\t'+ exon_id_list +'\t'+ annotation_list +'\n'
            data.write(values)
    data.close()

def checkforEnsemblExons(trans_exon_data):
    proceed_status = 'no'
    for (start_pos,ste,spe) in trans_exon_data:
        #print ste.ExonRegionNumbers(),spe.ExonRegionNumbers(), [ste.ExonID(),spe.ExonID()];kill
        if ste.ExonID() == '' or spe.ExonID() == '': print ste.ExonRegionNumbers(),spe.ExonRegionNumbers(), [ste.ExonID(),spe.ExonID()];kill
        if ('ENS' in ste.ExonID()) or ('ENS' in spe.ExonID()): proceed_status = 'yes'
        #else: print ste.ExonRegionNumbers(),spe.ExonRegionNumbers(), [ste.ExonID(),spe.ExonID()];kill
    return proceed_status

def processEnsExonStructureData(exon_annotation_db,exon_regions,transcript_gene_db,gene_transcript,transcript_exon_db,intron_retention_db):
    ###Parses transcript to exon relationships and links these to previously described distinct exon regions.
    ###The exon blocks and region numbers provide a simple semantic for determining where in the transcript and what event is occuring
    ###when directly comparing junctions to each other
    transcript_exon_db={}; null=[]; k=[]
    for (gene,chr,strand) in exon_annotation_db:
        for exon_info in exon_annotation_db[(gene,chr,strand)]:
            if strand == '+': exon_start,exon_end,y = exon_info ###Link the original exon-start/stop data back to the region based data
            else: exon_end,exon_start,y = exon_info ###This shouldn't be necessary, but the coordinates get reversed in annotate_exons and copying the original via deepcopy takes way too much time
            if gene in exon_regions:
                ste=''; spe=''
                for ed in exon_regions[gene]:
                    if exon_start == ed.ExonStart(): ste = ed
                    if exon_end == ed.ExonStop(): spe = ed
                    if spe!='' and ste!='': break 
                    #print ed.ExonStart(),ed.ExonStop()
                if spe!='' and ste!='':
                    values = exon_start,ste,spe #,y
                    ste.reSetExonID(y.ExonID()); spe.reSetExonID(y.ExonID())  ###Need to represent the original source ExonID to eliminate UCSC transcripts without Ensembl exons
                    for ens_transcriptid in y.TranscriptID():
                        try: transcript_exon_db[ens_transcriptid].append(values)
                        except KeyError: transcript_exon_db[ens_transcriptid] = [values]
                else:
                    #"""
                    print len(exon_regions[gene]), gene,chr,strand,[spe],[ste],exon_start,exon_end
                    #for ed in exon_regions[gene]: print [ed.ExonID()],ed.ExonStart(),ed.ExonStop(), ed.ExonRegionNumbers()
                    for ed in exon_regions[gene]:
                        #if ed.ExonID() == ens_exonid: print exon_start,exon_end,ens_exonid, ed.ExonStart(),ed.ExonStop(), ed.ExonID(), ed.ExonRegionNumbers();kill
                        print [ed.ExonID()],[exon_start],[exon_end], [ed.ExonStart()],[ed.ExonStop()], [ed.ExonRegionNumbers()]

                    kill#"""
                    null.append(gene)
            else: k.append(gene)
    null = unique.unique(null)      
    print len(null)
    print len(k), 'genes, not matched up to region database'
    print len(transcript_gene_db), 'transcripts from Ensembl being analyzed'
    
    transcript_exon_db = eliminate_redundant_dict_values(transcript_exon_db)
    gene_transcript = eliminate_redundant_dict_values(gene_transcript)
    
    gene_transcript_multiple={}; tc=0
    for gene in gene_transcript:
        if len(gene_transcript[gene])>1: tc+=1; gene_transcript_multiple[gene]=len(gene_transcript[gene])
    print tc,"genes with multiple transcripts associated from Ensembl"

    ###If a retained intron is present we must ignore all exons downstream of where we DELETED that exon information (otherwise there is false junciton information introduced)    
    ###Here we simply delete the information for that transcript all together
    td=0
    for key in intron_retention_db:
        transcripts_to_delete = {}
        for (s1,s2,y) in intron_retention_db[key]:
            if y.IntronDeletionStatus() == 'yes': ###only do this for exon's deleted upon import
                for transcript in y.TranscriptID(): transcripts_to_delete[transcript] = []
        for transcript in transcripts_to_delete:
            ###may not be present if the exon deleted constituted the whole transcript
            if transcript in transcript_exon_db: del transcript_exon_db[transcript]; td+=1            

    print td, "transcripts deleted with intron retention. Required for down-stream junction analysis"
    exon_junction_db={}; junction_transcript_db={}; rt=0
    ###Sort and filter the junction data
    for transcript in transcript_exon_db:
        gene,chr,strand = transcript_gene_db[transcript]
        transcript_exon_db[transcript].sort()
        if strand == '-': transcript_exon_db[transcript].reverse()
        index=0
        ###Introduced a filter to remove transcripts from UCSC with no supporting Ensembl exons (distinct unknown transcript type)
        ###Since these are already exon regions, we exclude them from alt. splicing/promoter assignment.
        proceed_status = checkforEnsemblExons(transcript_exon_db[transcript])
        ###Loop through the exons in each transcript
        if proceed_status == 'yes':
            for (start_pos,ste,spe) in transcript_exon_db[transcript]:
                if (index+1) != len(transcript_exon_db[transcript]): ###Don't grab past the last exon in the transcript
                    start_pos2,ste2,spe2 = transcript_exon_db[transcript][index+1]
                    exon_junction = (ste.ExonRegionNumbers(),spe.ExonRegionNumbers()),(ste2.ExonRegionNumbers(),spe2.ExonRegionNumbers())
                    try: exon_junction_db[gene].append(exon_junction)
                    except KeyError: exon_junction_db[gene] = [exon_junction]
                    try: junction_transcript_db[gene,exon_junction].append(transcript)
                    except KeyError: junction_transcript_db[gene,exon_junction] = [transcript]
                index+=1
        else:
            #print transcript
            rt +=1
    print rt, "transcripts removed from analysis with no Ensembl exon evidence. Results in more informative splicing annotations downstream"
    print len(junction_transcript_db), 'length of junction_transcript_db'
    print len(exon_junction_db),'length of exon_junction_db'
    
    ###Stringent count, since it requires all region information for each exon and common splice events occur for just one region to another     
    ###example: (((8, 1), (8, 1)), ((9, 1), (9, 1))), (((8, 1), (8, 1)), ((9, 1), (9, 2)))
    putative_as_junction_db={}
    for gene in exon_junction_db:
        junctions = exon_junction_db[gene]
        junction_count={}
        for junction in exon_junction_db[gene]:
            try: junction_count[junction]+=1
            except KeyError: junction_count[junction]=1
        count_junction={}; count_list=[]
        for junction in junction_count:
            count = junction_count[junction]
            try: count_junction[count].append(junction)
            except KeyError: count_junction[count] = [junction]
            count_list.append(count)
        count_list = unique.unique(count_list); count_list.sort() ###number of unique counts
        if len(count_list)>1 and gene in gene_transcript_multiple: ###Otherwise, there is no variation in junction number between transcripts
            transcript_number = gene_transcript_multiple[gene]
            max_count = count_list[-1] ###most common junction - not alternatively spliced
            if max_count == transcript_number: ###Ensures we are grabbing everything but constitutive exons (max_count can include AS if no Ensembl constitutive).
                for count in count_junction:
                    if count != max_count:
                        junctions = count_junction[count]
                        junctions.sort()
                        try: putative_as_junction_db[gene]+=junctions
                        except KeyError: putative_as_junction_db[gene]=junctions
            else:
                try: putative_as_junction_db[gene]+=junctions
                except KeyError: putative_as_junction_db[gene]=junctions
        elif gene in gene_transcript_multiple:
            ###If there are multiple transcripts, descriminating these is difficult, just include all junctions for that gene
            try: putative_as_junction_db[gene]+=junctions
            except KeyError: putative_as_junction_db[gene]=junctions
    return exon_junction_db,putative_as_junction_db,exon_junction_db

def compareJunctions(putative_as_junction_db,exon_regions):
    ###Find splice events based on structure based evidence
    print len(putative_as_junction_db),'genes being examined for AS/alt-promoters in Ensembl'
    critical_exon_db={}; j=0; global add_to_for_terminal_exons; add_to_for_terminal_exons={}
    complex3prime_event_db={}; complex5prime_event_db={}; cassette_exon_record={}
    for gene in putative_as_junction_db:
        for j1 in putative_as_junction_db[gene]:
            for j2 in putative_as_junction_db[gene]: ### O^n squared query
                if j1 != j2:
                    temp_junctions = [j1,j2]; temp_junctions.sort(); junction1,junction2 = temp_junctions
                    splice_junctions=[]; critical_exon=[]; splice_type=''
                    e1a,e2a = junction1; e1b,e2b = junction2 ###((8, 2), (8, 2)) break down the exon into exon_block,region tubles.
                    e1a3,e1a5 = e1a; e2a3,e2a5 = e2a; e1b3,e1b5 = e1b; e2b3,e2b5 = e2b ###(8, 2) break down the exons into single tuples (designating 5' and 3' ends of the exon): 
                    e1a3_block,e1a3_reg = e1a3; e1a5_block,e1a5_reg = e1a5; e2a3_block,e2a3_reg = e2a3; e2a5_block,e2a5_reg = e1a5
                    e1b3_block,e1b3_reg = e1b3; e1b5_block,e1b5_reg = e1b5 ;e2b3_block,e2b3_reg = e2b3; e2b5_block,e2b5_reg = e1b5
                    splice_junctions = [(e1a5,e2a3),(e1b5,e2b3)] ###three junctions make up the cassette event, record the two evidenced by this comparison and agglomerate after all comps
                    ###IMPORTANT NOTE: The temp_junctions are sorted, but doesn't mean that splice_junctions is sorted correctly... must account for this
                    splice_junctions2 = customLSDeepCopy(splice_junctions); splice_junctions2.sort()
                    if splice_junctions2 != splice_junctions: ###Then the sorting is wrong and the down-stream method won't work
                        ###Must re-do the above assingments
                        junction2,junction1 = temp_junctions
                        e1a,e2a = junction1; e1b,e2b = junction2 ###((8, 2), (8, 2)) break down the exon into exon_block,region tubles.
                        e1a3,e1a5 = e1a; e2a3,e2a5 = e2a; e1b3,e1b5 = e1b; e2b3,e2b5 = e2b ###(8, 2) break down the exons into single tuples (designating 5' and 3' ends of the exon): 
                        e1a3_block,e1a3_reg = e1a3; e1a5_block,e1a5_reg = e1a5; e2a3_block,e2a3_reg = e2a3; e2a5_block,e2a5_reg = e1a5
                        e1b3_block,e1b3_reg = e1b3; e1b5_block,e1b5_reg = e1b5 ;e2b3_block,e2b3_reg = e2b3; e2b5_block,e2b5_reg = e1b5
                        splice_junctions = [(e1a5,e2a3),(e1b5,e2b3)]
                    if e1a5_block == e2a3_block or e1b5_block == e2b3_block: continue ###suggests splicing within a block... we won't deal with these
                    if e1a5 == e1b5:  ###If 5'exons in the junction are the same
                        ###make sure the difference isn't in the 5' side of the next splice junction (or exon end)
                        if e2a3 != e2b3: #(((5, 1), (5, 1)*), ((6, 1)*, (6, 1)))  --  (((5, 1), (5, 1)*), ((7, 1)*, (7, 1)))
                            if e2a3_block == e2b3_block: #[(((1, 1), (1, 1)*), ((2, 1)*, (2, 1))) ----(((1, 1), (1, 1)*), ((2, 3)*, (2, 3)))]
                                splice_type = "alt-3'"; critical_exon = pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon)
                                y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y] 
                            else:
                                critical_exon = [e2a3]; splice_type = 'cassette-exon' #[(((1, 1), (1, 1)*), ((2, 1)*, (2, 1))) ----(((1, 1), (1, 1)*), ((3, 1)*, (3, 1)))]    
                                try: add_to_for_terminal_exons[gene,e2a3_block].append(e2b3)
                                except KeyError: add_to_for_terminal_exons[gene,e2a3_block] = [e2b3]
                                try: cassette_exon_record[gene,e2a3_block].append(e1a5_block)
                                except KeyError: cassette_exon_record[gene,e2a3_block] = [e1a5_block]
                                y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                                
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y]
                                #print critical_exon,splice_type,splice_junctions
                    if splice_type =='' and e2a3 == e2b3: ###If 3'exons in the junction are the same
                        if e1a5 != e1b5:
                            if e1a5_block == e1b5_block: ###simple alt 5' splice site encountered
                                splice_type = "alt-5'"; critical_exon = pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon)#[(((1, 1), (1, 1)*), ((2, 1)*, (2, 1))) ----(((1, 1), (1, 3)*), ((2, 1)*, (2, 1)))]
                                y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y] 
                            else:
                                splice_type = 'cassette-exon'; critical_exon = [e1b5] #[(((1, 1), (1, 1)*), ((3, 1)*, (3, 1))) ----(((2, 1), (2, 1)*), ((3, 1)*, (3, 1)))]       
                                try: add_to_for_terminal_exons[gene,e1b5_block].append(e1a5)
                                except KeyError: add_to_for_terminal_exons[gene,e1b5_block] = [e1a5]
                                try: cassette_exon_record[gene,e1b5_block].append(e2b3_block)
                                except KeyError: cassette_exon_record[gene,e1b5_block] = [e2b3_block]                                
                                #if gene == 'ENSG00000128606' and critical_exon == [(4, 1)] : print junction1,junction2;kill
                                y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y]
                                #print critical_exon,splice_type,splice_junctions
                    if splice_type =='' and e2a3_block == e2b3_block and e1a5_block != e2a3_block and e1b5_block != e2b3_block: ###Begin looking at complex examples: If 3'exon blocks in the junction are the same
                        if e1a5_block == e1b5_block: #alt5'-alt3' [(((1, 1), (1, 1)*), ((2, 1)*, (2, 1))) ----(((1, 3), (1, 3)*), ((2, 3)*, (2, 3)))]
                            critical_exon = pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon); splice_type = "alt5'-alt3'"
                            if len(critical_exon)>0:
                                alt5_exon = [critical_exon[0]]; alt3_exon = [critical_exon[1]]
                                #print alt5_exon,critical_exon,critical_exon;kill
                                y = CriticalExonInfo(gene,alt5_exon,"alt-5'",splice_junctions)
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y]
                                y = CriticalExonInfo(gene,alt3_exon,"alt-3'",splice_junctions)
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y] 
                        else: #cassette-alt3' [(((1, 1), (1, 1)*), ((4, 1)*, (4, 1))) ----(((2, 1), (2, 3)*), ((4, 3)*, (4, 3)))]
                            critical_exon = pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon)
                            splice_type = "alt-3'"
                            y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                            try: critical_exon_db[gene].append(y)
                            except KeyError: critical_exon_db[gene] = [y]
                            if e1a5_block < e1b5_block:
                                critical_exon = [e1b5]
                                try: add_to_for_terminal_exons[gene,e1b5_block] = [e1a5]
                                except KeyError: add_to_for_terminal_exons[gene,e1b5_block].append(e1a5)
                                complex3prime_event_db[gene,e1b5_block] = e1a5
                                try: cassette_exon_record[gene,e1b5_block].append(e2b3_block)
                                except KeyError: cassette_exon_record[gene,e1b5_block] = [e2b3_block]
                            elif e1a5_block != e1b5_block:
                                critical_exon = [e1a5]
                                try: add_to_for_terminal_exons[gene,e1a5_block] = [e1b5]
                                except KeyError: add_to_for_terminal_exons[gene,e1a5_block].append(e1b5)
                                complex3prime_event_db[gene,e1a5_block] = e1b5
                                try: cassette_exon_record[gene,e1a5_block].append(e2a3_block)
                                except KeyError: cassette_exon_record[gene,e1a5_block] = [e2a3_block]
                            splice_type = "cassette-exon"
                            y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                            try: critical_exon_db[gene].append(y)
                            except KeyError: critical_exon_db[gene] = [y] 
                    if splice_type =='' and e1a5_block == e1b5_block and e1a5_block != e2a3_block and e1b5_block != e2b3_block:
                        #alt5'-cassette' [(((1, 1), (1, 1)*), ((4, 1)*, (4, 1))) ----(((1, 1), (1, 3)*), ((5, 1)*, (5, 1)))]
                        critical_exon = pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon)
                        splice_type = "alt-5'"
                        y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                        try: critical_exon_db[gene].append(y)
                        except KeyError: critical_exon_db[gene] = [y]
                            
                        if e2a3_block < e2b3_block:
                            critical_exon = [e2a3]
                            try: add_to_for_terminal_exons[gene,e2a3_block] = [e2b3]
                            except KeyError: add_to_for_terminal_exons[gene,e2a3_block].append(e2b3)
                            complex5prime_event_db[gene,e2a3_block] = e2b3
                            try: cassette_exon_record[gene,e2a3_block].append(e1a5_block)
                            except KeyError: cassette_exon_record[gene,e2a3_block] = [e1a5_block]
                        elif e2a3_block != e2b3_block:
                            critical_exon = [e2b3]
                            try: add_to_for_terminal_exons[gene,e2b3_block] = [e2a3]
                            except KeyError: add_to_for_terminal_exons[gene,e2b3_block].append(e2a3)
                            complex5prime_event_db[gene,e2b3_block] = e2a3
                            try: cassette_exon_record[gene,e2b3_block].append(e1b5_block)
                            except KeyError: cassette_exon_record[gene,e2b3_block] = [e1b5_block]
                        splice_type = "cassette-exon"
                        y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                        try: critical_exon_db[gene].append(y)
                        except KeyError: critical_exon_db[gene] = [y]
                    if splice_type =='' and e1a5_block<e1b5_block and e2b3_block>e2a3_block and e2a3_block>e1b5_block:
                        #mx-mix [(((1, 1), (1, 1)*), ((4, 1)*, (4, 1))) ----(((2, 1), (2, 1)*), ((5, 1)*, (5, 1)))]
                        critical_exon = [e2a3,e1b5]#; mx_event_db[gene,e2a3] = e1b5; mx_event_db[gene,e1b5] = e2a3
                        splice_type = 'cassette-exon'
                        #"""
                        try: add_to_for_terminal_exons[gene,e2a3_block].append(e2b3)
                        except KeyError: add_to_for_terminal_exons[gene,e2a3_block] = [e2b3]                        
                        try: add_to_for_terminal_exons[gene,e1b5_block].append(e1a5)
                        except KeyError: add_to_for_terminal_exons[gene,e1b5_block] = [e1a5]
                        #"""       
                        try: cassette_exon_record[gene,e2a3_block].append(e1a5_block)
                        except KeyError: cassette_exon_record[gene,e2a3_block] = [e1a5_block]
                        try: cassette_exon_record[gene,e1b5_block].append(e2b3_block)
                        except KeyError: cassette_exon_record[gene,e1b5_block] = [e2b3_block]   
                        #print 'mx-mx',critical_exon, splice_junctions
                        y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                        try: critical_exon_db[gene].append(y)
                        except KeyError: critical_exon_db[gene] = [y]
                        #print splice_type,critical_exon, gene 
                    if splice_type =='' and e1a5_block<e1b5_block and e2a3_block>e2b3_block:
                        #(((2, 1), (2, 1)), ((9, 6), (9, 9))) (((4, 2), (4, 4)), ((7, 1), (7, 1))) ###one junction inside another
                        splice_type = "cassette-exon"; critical_exon = [e1b5,e2b3]
                        #"""

                        try: add_to_for_terminal_exons[gene,e2b3_block].append(e2a3)
                        except KeyError: add_to_for_terminal_exons[gene,e2b3_block] = [e2a3]                        
                        try: add_to_for_terminal_exons[gene,e1b5_block].append(e1a5)
                        except KeyError: add_to_for_terminal_exons[gene,e1b5_block] = [e1a5]
                        
                        try: cassette_exon_record[gene,e2b3_block].append(e1b5_block)
                        except KeyError: cassette_exon_record[gene,e2b3_block] = [e1b5_block]
                        try: cassette_exon_record[gene,e1b5_block].append(e2b3_block)
                        except KeyError: cassette_exon_record[gene,e1b5_block] = [e2b3_block]
                        #"""
                        y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                        try: critical_exon_db[gene].append(y)
                        except KeyError: critical_exon_db[gene] = [y]      
    ###Determine unique splice events and improve the annotations
    print len(critical_exon_db), 'genes identified from Ensembl, with alternatively regulated junctions'
    cassette_exon_record = eliminate_redundant_dict_values(cassette_exon_record)
    """for i in cassette_exon_record:
        print i, cassette_exon_record[i]"""
    add_to_for_terminal_exons = eliminate_redundant_dict_values(add_to_for_terminal_exons)
    ###Store all region information in a dictionary for efficient recall
    global region_db; region_db={}
    for gene in exon_regions:
        for rd in exon_regions[gene]:
            try: region_db[gene,rd.ExonRegionNumbers()]=rd
            except AttributeError: print gene, rd;kill
            
    alternative_exon_db={}; critical_junction_db={}; critical_gene_junction_db={}
    for gene in critical_exon_db:
        critical_exon_junctions={}; critical_exon_splice_type={}
        for sd in critical_exon_db[gene]:
            for critical_exon in sd.CriticalExonRegion():
                try: critical_exon_junctions[critical_exon]+=sd.Junctions()
                except KeyError: critical_exon_junctions[critical_exon]=sd.Junctions()
                try: critical_exon_splice_type[critical_exon].append(sd.SpliceType())
                except KeyError: critical_exon_splice_type[critical_exon]=[sd.SpliceType()]
                for junction in sd.Junctions():
                    try: critical_junction_db[tuple(junction)].append(junction)
                    except KeyError: critical_junction_db[tuple(junction)]=[sd.SpliceType()]
        critical_exon_junctions = eliminate_redundant_dict_values(critical_exon_junctions)
        critical_exon_splice_type = eliminate_redundant_dict_values(critical_exon_splice_type)
        for critical_exon in critical_exon_junctions:
            cj = critical_exon_junctions[critical_exon]
            splice_events = critical_exon_splice_type[critical_exon]
            status = 'stop'
            #print splice_events,critical_exon
            if splice_events == ['cassette-exon'] or ((gene,critical_exon[0]) in complex3prime_event_db) or ((gene,critical_exon[0]) in complex5prime_event_db):
                exons_blocks_joined_to_critical = cassette_exon_record[gene,critical_exon[0]]
                cassette_status = check_exon_polarity(critical_exon[0],exons_blocks_joined_to_critical)
                if len(critical_exon_junctions[critical_exon])<3 or cassette_status == 'no': ###Thus, it is not supported by 3 independent junctions
                    if len(exons_blocks_joined_to_critical)<2 or cassette_status == 'no':
                        if cj[0][1] == cj[1][1]:
                            splice_events = ['alt-N-term']; second_critical_exon = add_to_for_terminal_exons[gene,critical_exon[0]]; status = 'add_another'
                        elif cj[0][0] == cj[1][0]:
                            splice_events = ['alt-C-term']; second_critical_exon = add_to_for_terminal_exons[gene,critical_exon[0]]; status = 'add_another'
                        else:
                            if critical_exon == cj[0][1]: splice_events = ['alt-C-term'] ###this should be the only alt-exon
                elif (gene,critical_exon[0]) in complex3prime_event_db:
                    #print '3prime',splice_events,critical_exon
                    if (gene,critical_exon[0]) in add_to_for_terminal_exons:
                        #print critical_exon,len(cassette_exon_record[gene,critical_exon[0]]),cassette_exon_record[gene,critical_exon[0]];kill
                        if len(exons_blocks_joined_to_critical)<2 or cassette_status == 'no':
                            second_critical_exon = add_to_for_terminal_exons[gene,critical_exon[0]]
                            splice_events = ['alt-N-term']; status = 'add_another'
                elif (gene,critical_exon[0]) in complex5prime_event_db:
                    #print '5prime',splice_events,critical_exon
                    if (gene,critical_exon[0]) in add_to_for_terminal_exons:
                        #print critical_exon,len(cassette_exon_record[gene,critical_exon[0]]),cassette_exon_record[gene,critical_exon[0]];kill
                        if len(exons_blocks_joined_to_critical)<2 or cassette_status == 'no':
                            second_critical_exon = add_to_for_terminal_exons[gene,critical_exon[0]]
                            splice_events = ['alt-C-term']; status = 'add_another'
                """if 'mx-mx' in splice_events and (gene,critical_exon) in mx_event_db:
                    ###if one exon is a true cassette exon, then the mx-mx is not valid
                    if (gene,critical_exon[0]) in add_to_for_terminal_exons:
                        second_critical_exon = add_to_for_terminal_exons[gene,critical_exon[0]]
                        #print gene,critical_exon,second_critical_exon;kill"""
            splice_events = string.join(splice_events,'|'); exon_junction_str_list=[]
            splice_events = string.replace(splice_events, 'cassette-exons','cassette-exon(s)')
            ###if the junctions comprising the splice event for an alt-cassette show evidence of multiple exons, annotate as such
            if "alt5'-cassette" in splice_events: #or "cassette-alt3'"
                for junction in cj:
                    splicing_annotations = critical_junction_db[tuple(junction)]
                    if 'cassette-exons' in splicing_annotations:
                        splice_events = string.replace(splice_events, "alt5'-cassette","alt5'-cassette(s)"); break
            if "cassette-alt3'" in splice_events: #or "cassette-alt3'"
                for junction in cj:
                    splicing_annotations = critical_junction_db[tuple(junction)]
                    if 'cassette-exons' in splicing_annotations:
                        splice_events = string.replace(splice_events, "cassette-alt3'","cassette(s)-alt3'"); break
            ###Currently, 'cassette-exon(s)' is redundant with "alt5'-cassette(s)" or "cassette(s)-alt3'", so simplify
            if "alt5'-cassette(s)" in splice_events and 'cassette-exon(s)' in splice_events:
                splice_events = string.replace(splice_events, 'cassette-exon(s)','')
            if "cassette(s)-alt3'" in splice_events and 'cassette-exon(s)' in splice_events:
                splice_events = string.replace(splice_events, 'cassette-exon(s)','')
            splice_events = string.replace(splice_events, '||','|')
            for j in cj:
                nj=[]
                for exon in j: e = 'E'+str(exon[0])+'.'+str(exon[1]); nj.append(e)
                try: critical_gene_junction_db[gene].append(nj)
                except KeyError: critical_gene_junction_db[gene] = [nj]
                nj = string.join(nj,'-')
                exon_junction_str_list.append(nj)
            exon_junction_str = string.join(exon_junction_str_list,'|')
            try: rd = region_db[gene,critical_exon] ###('ENSG00000213588', (26, 1)) and ('ENSG00000097007', (79, 1)) didn't work
            except KeyError:
                ###Occurs as a results of either exons or transcripts with sketchy or complex assignments
                null = []
            se = rd.AssociatedSplicingEvent()
            if len(se)>1:
                if splice_events not in se: se = se+'|'+ splice_events
            else: se = splice_events
            rd.setSpliceData(se,exon_junction_str)
                    
            #print critical_exon,se,exon_junction_str,gene
            if status == 'add_another':
                for critical_exon in second_critical_exon: 
                    rd = region_db[gene,critical_exon]
                    se = rd.AssociatedSplicingEvent()
                    if len(se)>1:
                        if splice_events not in se:
                            if se != 'cassette-exon': se = se+'|'+ splice_events
                    else: se = splice_events
                    rd.setSpliceData(se,exon_junction_str)
                    #print critical_exon, se, exon_junction_str, gene,'second'
            """
            ###create an index for easy searching of exon content in downstream modules
            critical_exon_block = critical_exon[0]
            if gene in alternative_exon_db:
                block_db = alternative_exon_db[gene]
                try: block_db[critical_exon_block].append(rd)
                except KeyError: block_db[critical_exon_block] = [rd]
            else:
                block_db = {}; block_db[critical_exon_block]=[rd]
                alternative_exon_db[gene]=block_db"""
            
    ###Since setSpliceData will update the existing instance, we can just re-roder the region db for easy searching in downstream modules
    ### (note: the commented out code above could be useful for exon-structure output)
    for gene in exon_regions:
        block_db = {}
        for rd in exon_regions[gene]:
            try: block_db[rd.ExonNumber()].append(rd)
            except KeyError: block_db[rd.ExonNumber()] = [rd]
        exon_regions[gene] = block_db ###Replace the existing list with a dictionary for faster look-ups
    return exon_regions,critical_gene_junction_db

def check_exon_polarity(critical_exon_block,exons_blocks_joined_to_critical):
    g=0;l=0
    for joined_exon_blocks in exons_blocks_joined_to_critical:
        if joined_exon_blocks>critical_exon_block: g+=1
        if joined_exon_blocks<critical_exon_block: l+=1
    if g>0 and l>0: cassette_status = 'yes'
    else: cassette_status = 'no'
    return cassette_status

def pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon):
    e1a5_block,e1a5_reg = e1a5
    e2a3_block,e2a3_reg = e2a3
    e1b5_block,e1b5_reg = e1b5
    e2b3_block,e2b3_reg = e2b3

    if e1a5_block == e1b5_block:                    
        if e1a5_reg < e1b5_reg: critical_exon += [e1b5]
        elif e1a5_reg != e1b5_reg: critical_exon += [e1a5]
    if e2a3_block == e2b3_block:
        if e2a3_reg < e2b3_reg: critical_exon += [e2a3]
        elif e2a3_reg != e2b3_reg: critical_exon += [e2b3]
    return critical_exon

def customDBDeepCopy(db):
    db2={}
    for i in db:
        for e in db[i]:
            try: db2[i].append(e)
            except KeyError: db2[i]=[e]
    return db2

def customLSDeepCopy(ls):
    ls2=[]
    for i in ls: ls2.append(i)
    return ls2

if __name__ == '__main__':
    ###KNOWN PROBLEMS: the junction analysis program calls exons as cassette-exons if there a new C-terminal exon occurs downstream of that exon in a different transcript (ENSG00000197991).
    species = 'Hs'
    #species = 'Mm'
    test = 'yes'
    data_type = 'ncRNA'
    data_type = 'mRNA'
    test_gene = ['ENSG00000201902']#,'ENSG00000154889','ENSG00000156026','ENSG00000148584','ENSG00000063176','ENSG00000126860'] #['ENSG00000105968']
    meta_test = ["ENSG00000215305","ENSG00000179676","ENSG00000170484","ENSG00000138180","ENSG00000100258","ENSG00000132170","ENSG00000105767","ENSG00000105865","ENSG00000108523","ENSG00000150045","ENSG00000156026"]
    #test_gene = meta_test
    #critical_exon_seq_db = import_sequence_data(gene_seq_filename,critical_exon_seq_db,species);kill
    import_dir = '/AltDatabase/ensembl/'+species
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for file_name in dir_list:    #loop through each file in the directory to output results
        dir_file = 'AltDatabase/ensembl/'+species+'/'+file_name
        if 'exon' in dir_file: exon_file = dir_file
        elif 'transcript' in dir_file: trans_file = dir_file
        elif 'gene' in dir_file: gene_seq_file = dir_file
        elif 'Exon_cDNA' in dir_file: exon_trans_file = dir_file
        elif 'Domain' in dir_file: domain_file = dir_file
    #"""
    exon_annotation_db,transcript_gene_db,gene_transcript,transcript_exon_db,intron_retention_db,ucsc_splicing_annot_db = getEnsExonStructureData(species,data_type)
    KILL

    #exon_db = customDBDeepCopy(exon_annotation_db)
    #"""
    exon_annotation_db2 = annotate_exons(exon_annotation_db)
    #kill
    exon_db2 = customDBDeepCopy(exon_annotation_db2) ##having problems with re-writting contents of this db when I don't want to
    
    exon_clusters,intron_clusters,exon_regions,intron_region_db = exon_clustering(exon_db2); exon_db2={}
    #"""
    exon_junction_db,putative_as_junction_db,exon_junction_db = processEnsExonStructureData(exon_annotation_db,exon_regions,transcript_gene_db,gene_transcript,transcript_exon_db,intron_retention_db)
    
    #ej = {};ej['ENSG00000124104'] = exon_junction_db['ENSG00000124104']
    #exon_regions,critical_gene_junction_db = compareJunctions(putative_as_junction_db,exon_regions)
    #exportSubGeneViewerData(exon_regions,critical_gene_junction_db,intron_region_db,intron_retention_db)        
    #ENSG00000149792 possible retained 3' intron... check out
    kill
    use_exon_data='no';get_splicing_factors = 'yes'
    rna_processing_ensembl = GO_parsing.parseAffyGO(use_exon_data,get_splicing_factors,species)
    ensembl_annot_file = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations_simple.txt'
    ensembl_annotation_db = getEnsemblAnnotations(ensembl_annot_file,rna_processing_ensembl)    
    
    ###Print ovelap statistics for exon-blocks
    x=0;y=0; m=0; l=0
    for key in exon_clusters:
        for exon_block in exon_clusters[key]:
            if len(exon_block[2])>1: y += 1; x += 1; m += len(exon_block[2]); l += len(exon_block[2])
            else: x += 1; m += 1
            #if x < 50:
                #print key[0], exon_block[2], len(exon_block[2]),x,y,m,l
    print 'x',x,'y',y,'m',m,'l',l
    """
    for gene in exon_regions:
        db = exon_regions[gene]
        for block in db:
            for rd in db[block]:
                try: print rd.AssociatedSplicingEvent(),rd.AssociatedSplicingJunctions();kill
                except AttributeError: continue"""