###ExonSeqModule
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
import copy
import time
import update

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".TXT" or entry[-3:] == ".fa":
            dir_list2.append(entry)
    return dir_list2
           
def importSplicingAnnotationDatabase(filename,array_type):
    global exon_db; fn=filepath(filename)
    print 'importing', filename
    exon_db={}; count = 0; x = 0
    for line in open(fn,'r').readlines():             
        probeset_data,null = string.split(line,'\n')  #remove endline
        if x == 0: x = 1
        else:
            try: probeset_id, exon_id, ensembl_gene_id, transcript_cluster_id, chromosome, strand, probeset_start, probeset_stop, affy_class, constitutitive_probeset, ens_exon_ids, ens_const_exons, exon_region, exon_region_start, exon_region_stop, splicing_event, splice_junctions = string.split(probeset_data,'\t')
            except Exception: t = string.split(probeset_data,'\t'); print len(t), t;kill
            #probe_data = AffyExonSTData(probeset_id,ensembl_gene_id,exon_id,ens_exon_ids,transcript_cluster_id, chromosome, strand, probeset_start, probeset_stop, affy_class, constitutitive_probeset, exon_region, splicing_event, splice_junctions)
            if ':' in probeset_id:
                id1,id2 = string.split(probeset_id,':')
                if id2[0]=='E' or id2[0]=='I': probeset_id = probeset_id ### This import method applies to 
                else: probeset_id = id2
            ### Store exon region start and stop since this corresponds to the critical-exon-region-seq_update.txt file sequence coordinates
            probe_data = AffyExonSTDataSimple(probeset_id,ensembl_gene_id,exon_region,ens_exon_ids,probeset_start,probeset_stop,strand,chromosome)
            exon_db[probeset_id] = probe_data
    print len(exon_db), 'probesets imported'
    return exon_db

class SplicingAnnotationData:
    def ArrayType(self):
        self._array_type = array_type
        return self._array_type
    def Probeset(self): return self._probeset
    def IncludeProbeset(self):
        include_probeset = 'yes'
        if self.ArrayType() == 'AltMouse':
            if  filter_probesets_by == 'exon': 
                if '-' in self.ExonID() or '|' in self.ExonID(): ###Therfore the probeset represents an exon-exon junction or multi-exon probeset
                    include_probeset = 'no'
            if  filter_probesets_by != 'exon': 
                if '|' in self.ExonID(): include_probeset = 'no'
        if self.ArrayType() == 'exon':
            if filter_probesets_by == 'core':
                if self.ProbesetClass() != 'core': include_probeset = 'no'
        return include_probeset
    def ExonID(self): return self._exonid
    def GeneID(self): return self._geneid
    def ExternalGeneID(self): return self._external_gene
    def ProbesetType(self):
        ###e.g. Exon, junction, constitutive(gene)
        return self._probeset_type
    def GeneStructure(self): return self._block_structure
    def SecondaryExonID(self): return self._block_exon_ids
    def Chromosome(self): return self._chromosome
    def Strand(self): return self._strand
    def ProbeStart(self): return int(self._start)
    def ProbeStop(self): return int(self._stop)
    def ProbesetClass(self): 
        ###e.g. core, extendended, full
        return self._probest_class
    def ExternalExonIDs(self): return self._external_exonids
    def ExternalExonIDList(self):
        external_exonid_list = string.split(self.ExternalExonIDs(),'|')
        return external_exonid_list
    def Constitutive(self): return self._constitutive_status
    def SecondaryGeneID(self): return self._secondary_geneid
    def SetExonSeq(self,seq): self._exon_seq = seq
    def ExonSeq(self): return string.upper(self._exon_seq)
    def SetJunctionSeq(self,seq): self._junction_seq = seq
    def JunctionSeq(self): return string.upper(self._junction_seq)
    def RecipricolProbesets(self): return self._junction_probesets
    def ExonRegionID(self): return self._exon_region
    def Report(self):
        output = self.Probeset() +'|'+ self.ExternalGeneID()
        return output
    def __repr__(self): return self.Report()

class AffyExonSTData(SplicingAnnotationData):
    def __init__(self,probeset_id,ensembl_gene_id,exon_id,ens_exon_ids,transcript_cluster_id, chromosome, strand, probeset_start, probeset_stop, affy_class, constitutitive_probeset, exon_region, splicing_event, splice_junctions):
        self._geneid = ensembl_gene_id; self._external_gene = ensembl_gene_id; self._exonid = exon_id
        self._constitutive_status = constitutitive_probeset; self._start = probeset_start; self._stop = probeset_stop
        self._external_exonids = ens_exon_ids; self._secondary_geneid = transcript_cluster_id; self._chromosome = chromosome
        self._probest_class = affy_class; self._probeset=probeset_id
        self._exon_region=exon_region;  self._splicing_event=splicing_event; self._splice_junctions=splice_junctions
        if self._exonid[0] == 'U': self._probeset_type = 'UTR'
        elif self._exonid[0] == 'E': self._probeset_type = 'exonic'
        elif self._exonid[0] == 'I': self._probeset_type = 'intronic'
        def ExonRegionID(self): return self._exon_region
        def SplicingEvent(self): return self._splicing_event
        def SpliceJunctions(self): return self._splice_junctions
        def Constitutive(self):
            if len(self._splicing_event)>0: return 'no' ###Over-ride affymetrix probeset file annotations if an exon is alternatively spliced
            else: return self._constitutive_status

class AffyExonSTDataSimple(SplicingAnnotationData):
    def __init__(self,probeset_id,ensembl_gene_id,exon_region,ens_exon_ids,probeset_start,probeset_stop,strand,chr):
        self._geneid = ensembl_gene_id; self._external_gene = ensembl_gene_id; self._exon_region = exon_region
        self._external_exonids = ens_exon_ids; self._probeset=probeset_id
        self._chromosome = chr
        try: self._start=int(probeset_start); self._stop = int(probeset_stop)
        except Exception: self._start=probeset_start; self._stop = probeset_stop
        self._strand = strand
    
class JunctionDataSimple(SplicingAnnotationData):
    def __init__(self,probeset_id,ensembl_gene_id,array_geneid,junction_probesets,critical_exons):
        self._geneid = ensembl_gene_id; self._junction_probesets = junction_probesets; self._exonid = critical_exons
        self._probeset=probeset_id; self._external_gene = array_geneid; self._external_exonids = ''
        
def importProbesetSeqeunces(filename,array_type,exon_db,chromosome,species):
    #output_file = 'input/'+species+'/filtered_probeset_sequence.txt'
    #fn=filepath(output_file)
    #datar = open(fn,'w')
                                  
    print 'importing', filename
    probeset_seq_db={}; probesets_parsed=0; probesets_added=0
    chromosome = 'chr'+str(chromosome)
    print "Begining generic fasta import of",filename
    fn=filepath(filename)
    sequence = ''; x = 0;count = 0
    for line in open(fn,'r').xreadlines():
        try: data, newline= string.split(line,'\n')
        except ValueError: continue
        try:
            if data[0] == '>':
                    try: 
                        try:
                            y = exon_db[probeset]
                            sequence = string.upper(sequence)
                            gene = y.GeneID(); y.SetExonSeq(sequence)
                            try: probeset_seq_db[gene].append(y)
                            except KeyError: probeset_seq_db[gene] = [y]
                            sequence = ''; t= string.split(data,';'); probesets_added+=1
                            probeset=t[0]; probeset_data = string.split(probeset,':'); probeset = probeset_data[-1]
                            chr = t[2]; chr = string.split(chr,'='); chr = chr[-1]; probesets_parsed +=1
                        except KeyError:
                            sequence = ''; t= string.split(data,';')
                            probeset=t[0]; probeset_data = string.split(probeset,':'); probeset = probeset_data[-1]
                            chr = t[2]; chr = string.split(chr,'='); chr = chr[-1]; probesets_parsed +=1
                            #if chr == chromosome: go = 'yes'
                            #else: go = 'no'
                    except UnboundLocalError: ###Occurs for the first entry
                        t= string.split(data,';')
                        probeset=t[0]; probeset_data = string.split(probeset,':'); probeset = probeset_data[-1]
                        chr = t[2]; chr = string.split(chr,'='); chr = chr[-1]; probesets_parsed +=1
                        #if chr == chromosome: go = 'yes'
                        #else: go = 'no'
            else: sequence = sequence + data
        except IndexError: continue
        
    try:
        y = exon_db[probeset]
        sequence = string.upper(sequence)
        gene = y.GeneID(); y.SetExonSeq(sequence)
        try: probeset_seq_db[gene].append(y)
        except KeyError: probeset_seq_db[gene] = [y]
    except KeyError: null=[]

    #datar.close()
    #exportAssociations(probeset_seq_db,species)
    print len(probeset_seq_db), probesets_added, probesets_parsed, len(exon_db)
    #print probeset_seq_db['ENSG00000156006']
    return probeset_seq_db

def exportAssociations(probeset_seq_db,species):
    output_file = 'AltAnalyze/'+species+'/SequenceData/'+species+'/filtered_probeset_sequence.txt'
    fn=filepath(output_file)
    data = open(fn,'w')
    for probeset in probeset_seq_db:
        seq = probeset_seq_db[probeset]
        data.write(probeset+'\t'+seq+'\n')
    data.close()

def getParametersAndExecute(probeset_seq_file,array_type,species):
    probeset_annotations_file = 'AltDatabase/'+species+'/exon/'+species+'_Ensembl_probesets.txt'
    ###Import probe-level associations
    exon_db = importSplicingAnnotationDatabase(probeset_annotations_file,array_type)
    start_time = time.time(); chromosome = 1
    probeset_seq_db = importProbesetSeqeunces(probeset_seq_file,array_type,exon_db,chromosome,species)
    if array_type == 'gene':
        probeset_annotations_file = string.replace(probeset_annotations_file,'exon','gene')
        gene_exon_db = importSplicingAnnotationDatabase(probeset_annotations_file,array_type)
        translation_db = exportAllExonToGeneArrayAssociations(exon_db,gene_exon_db)
        probeset_seq_db = convertExonProbesetSequencesToGene(probeset_seq_db,gene_exon_db,translation_db)
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "Analyses finished in %d seconds" % time_diff
    return probeset_seq_db

def convertExonProbesetSequencesToGene(probeset_seq_db,gene_exon_db,translation_db):
    probeset_gene_seq_db={}; probeset_db = {}
    print len(probeset_seq_db), 'gene entries with exon probeset sequence'
    for probeset in gene_exon_db:
        y = gene_exon_db[probeset]
        if probeset in translation_db:
            try:
                for z in probeset_seq_db[y.GeneID()]:
                    if z.Probeset() in translation_db[probeset]:
                        y.SetExonSeq(z.ExonSeq()) ### Use the exon array sequence if the exon regions are the same
                        try: probeset_gene_seq_db[y.GeneID()].append(y)
                        except KeyError: probeset_gene_seq_db[y.GeneID()] = [y]
                        probeset_db[probeset]=[]
            except KeyError: null=[]
    print len(probeset_db), 'gene probesets annotated with exon sequence'
    return probeset_gene_seq_db

def exportAllExonToGeneArrayAssociations(exon_db,gene_exon_db):
    ### Above function basically does this but some exon array probesets may not be sequence matched (probably not necessary)

    output_file = 'AltDatabase/'+species+'/gene/'+species+'_gene-exon_probesets.txt'
    fn=filepath(output_file)
    data = open(fn,'w')
    data.write('gene_probeset\texon_probeset\n')
    
    exon_db_gene = {}
    for probeset in exon_db:
        y = exon_db[probeset]
        try: exon_db_gene[y.GeneID()].append(y)
        except KeyError: exon_db_gene[y.GeneID()] = [y]
    #print gene_exon_db['10411047'].ExonRegionID(),exon_db['4355948'].ExonRegionID()
    translation_db={}
    for probeset in gene_exon_db:
        y = gene_exon_db[probeset]
        try:
            for z in exon_db_gene[y.GeneID()]:
                if z.ExonRegionID() == y.ExonRegionID():
                    try: translation_db[probeset].append(z.Probeset())
                    except KeyError: translation_db[probeset] = [z.Probeset()]
        except KeyError: null=[]
        
    ### Mop-ups (get imperfect relationships)
    for probeset in gene_exon_db:
        if probeset not in translation_db:
            y = gene_exon_db[probeset]
            try:
                for z in exon_db_gene[y.GeneID()]:
                    if (z.ExonRegionID()+'|') in (y.ExonRegionID()+'|') or (y.ExonRegionID()+'|') in (z.ExonRegionID()+'|'):
                        #print z.ExonRegionID(), y.ExonRegionID();kill
                        try: translation_db[probeset].append(z.Probeset())
                        except KeyError: translation_db[probeset] = [z.Probeset()]
            except KeyError: null=[]

    translation_db2={}
    for probeset in translation_db:
        y = gene_exon_db[probeset]
        #if probeset == '10411047': print '******',translation_db[probeset]
        for exon_probeset in translation_db[probeset]:
            z = exon_db[exon_probeset]
            coordinates = [y.ProbeStart(), y.ProbeStop(), z.ProbeStart(), z.ProbeStop()]
            coordinates.sort(); include = 'yes'
            ### Ensures that the coordinates overlap versus are separate
            if [y.ProbeStart(), y.ProbeStop()] == coordinates[:2]: include ='no'
            elif [y.ProbeStart(), y.ProbeStop()] == coordinates[-2:]: include = 'no'
            ### Include a probeset if the exon and gene probeset locations overlap or if there is only one associations
            if include == 'yes' or len(translation_db[probeset]) == 1:
                data.write(probeset+'\t'+exon_probeset+'\n')
                try: translation_db2[probeset].append(exon_probeset)
                except Exception: translation_db2[probeset] = [exon_probeset]
        if probeset not in translation_db2:
            ### If nothing got added, then just add the last match
            data.write(probeset+'\t'+exon_probeset+'\n')
            try: translation_db2[probeset].append(exon_probeset)
            except Exception: translation_db2[probeset] = [exon_probeset]
                
    data.close()
    print 'Out of all',len(gene_exon_db),'gene array probesets',len(translation_db2), 'had corresponding exon array probesets found.'
    return translation_db2

def getExonAnnotationsAndSequence(probeset_seq_file,array_type,species):
    probeset_annotations_file = 'AltDatabase/'+species+'/exon/'+species+'_Ensembl_probesets.txt'
    exon_db = importSplicingAnnotationDatabase(probeset_annotations_file,array_type)
    chromosome = 1
    ###By running this next function, we can update exon_db to include probeset sequence data
    array_type=''
    try: probeset_seq_db = importProbesetSeqeunces(probeset_seq_file,array_type,exon_db,chromosome,species)
    except IOError: null = []
    return exon_db

def import_ensembl_unigene(species):
    filename = 'AltDatabase/'+species+'/SequenceData/'+species+'_Ensembl-Unigene.txt'
    fn=filepath(filename); unigene_ensembl={}
    print 'importing', filename
    for line in open(fn,'r').xreadlines():        
        data, newline= string.split(line,'\n')
        ensembl,unigene = string.split(data,'\t')
        if len(unigene)>1 and len(ensembl)>1:
            try: unigene_ensembl[unigene].append(ensembl)
            except KeyError: unigene_ensembl[unigene] = [ensembl]
    print 'len(unigene_ensembl)',len(unigene_ensembl)
    return unigene_ensembl

def importEnsemblAnnotations(species):
    filename = 'AltDatabase/'+species+'/SequenceData/'+species+'_Ensembl-annotations.txt'
    fn=filepath(filename); symbol_ensembl={}
    print 'importing', filename
    for line in open(fn,'r').xreadlines():        
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t'); ensembl = t[0]
        try: symbol = t[2]
        except IndexError: symbol = ''
        if len(symbol)>1 and len(ensembl)>1:
            try: symbol_ensembl[symbol].append(ensembl)
            except KeyError: symbol_ensembl[symbol] = [ensembl]
    print 'len(symbol_ensembl)',len(symbol_ensembl)
    symbol_ensembl = eliminate_redundant_dict_values(symbol_ensembl)
    return symbol_ensembl

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def importmiRNATargetPredictionsAdvanced(species):
    filename = 'AltDatabase/'+species+'/SequenceData/miRBS-combined_gene-target-sequences.txt'
    print 'importing', filename; count = 0
    source_db={} ### collect stats on the input file sources
    fn=filepath(filename); ensembl_mirna_db={}; unigene_ids={}; x = 4000; z=0; nulls={}
    for line in open(fn,'r').xreadlines():
        line = string.replace(line,"'",''); line = string.replace(line,'"','')
        data, newline = string.split(line,'\n')
        microrna,ensembl,three_prime_utr_anchor_site_seq,sources = string.split(data,'\t')
        sources_list = string.split(sources,'|')
        for source in sources_list:
            try: source_db[source]+=1
            except Exception: source_db[source]=1
        y = MicroRNAData(ensembl,microrna,three_prime_utr_anchor_site_seq,sources)
        count+=1
        #y = microrna,three_prime_utr_anchor_site_seq,sources
        try: ensembl_mirna_db[ensembl].append(y)
        except KeyError: ensembl_mirna_db[ensembl] = [y]
    ensembl_mirna_db = eliminate_redundant_dict_values(ensembl_mirna_db)
    print count, "microRNA to target relationships imported"
    for source in source_db:
        print source,':',source_db[source],
    print ''
    return ensembl_mirna_db

def importmiRNATargetPredictions(species):
    unigene_ensembl = import_ensembl_unigene(species)
    symbol_ensembl = importEnsemblAnnotations(species)
    filename = 'AltDatabase/'+species+'/SequenceData//microRNA-target-annotated.txt'
    print 'importing', filename
    fn=filepath(filename); ensembl_mirna_db={}; unigene_ids={}; x = 4000; z=0; nulls={};k=0
    for line in open(fn,'r').xreadlines():
      if k == 0:k=1
      else:
        line = string.replace(line,"'",''); line = string.replace(line,'"','')
        data, newline = string.split(line,'\n')
        refseqid, unigene, symbol, symbol2, llid, llrepprotacc, microrna, rank, pictar_score, annotation, num_anchor_sites, three_prime_utr_anchor_site_seq = string.split(data,'\t')
        sequences = string.split(three_prime_utr_anchor_site_seq,' ')
        ensembls = []
        unigene_ids[unigene]=[]
        if unigene in unigene_ensembl: ensembls = unigene_ensembl[unigene]
        elif symbol in symbol_ensembl: ensembls = symbol_ensembl[symbol]
        if len(ensembls)>5: print len(ensembls)
        for ensembl in ensembls:
            for sequence in sequences:
                y = MicroRNAData(ensembl,microrna,sequence,'')
                #y = microrna,three_prime_utr_anchor_site_seq,''
                try: ensembl_mirna_db[ensembl].append(y)
                except KeyError: ensembl_mirna_db[ensembl] = [y]
        if len(ensembls)<1: nulls[unigene] = []
    ensembl_mirna_db = eliminate_redundant_dict_values(ensembl_mirna_db)
    print len(ensembl_mirna_db), "genes with associated microRNA data out of",len(unigene_ids)
    return ensembl_mirna_db

class MicroRNAData:
    def __init__(self, gene_id, microrna, sequence, source):
        self._gene_id = gene_id; self._microrna = microrna; self._sequence = sequence; self._source = source
        self.start_set = []
        self.end_set = []
        self.exon_sets=[]
    def GeneID(self): return self._gene_id
    def MicroRNA(self): return self._microrna
    def Sequence(self): return self._sequence
    def Source(self): return self._source
    def setEnsExons(self,ens_exons): self.exon_sets+=ens_exons
    def EnsExons(self):
        exonsets_unique = unique.unique(self.exon_sets)
        return string.join(exonsets_unique,'|')
    def setCoordinates(self,chr,strand,start,end):
        self.start_set.append(start)
        self.end_set.append(end)
        self.chr = chr
        self.strand = strand
    def Coordinates(self):
        x=0; coords=[]
        for i in self.start_set:
            coord = self.Chr()+':'+str(i)+'-'+str(self.end_set[x])
            coords.append(coord)
            x+=1
        coords = unique.unique(coords)
        coords = string.join(coords,'|') ###If multiple coordinates
        return coords
    def Chr(self): return self.chr
    def Strand(self): return self.strand
    def SummaryValues(self):
        output = self.GeneID()+'|'+self.MicroRNA()+'|'+self.Sequence()
        return output
    def __repr__(self): return self.SummaryValues()

def alignmiRNAData(array_type,mir_source,species,stringency,ensembl_mirna_db,splice_event_db):
    output_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_ensembl_microRNAs.txt'
    if mir_source == 'pictar': fn=filepath(output_file); data = open(fn,'w')
    print "Aligning microRNAs to probesets"
    added = {} ###not sure where
    probeset_miRNA_db={}
    for gene in ensembl_mirna_db:
        for y in ensembl_mirna_db[gene]:
            miRNA = y.MicroRNA(); miRNA_seq = y.Sequence(); sources = y.Source()
            if mir_source == 'pictar':
                if (gene,miRNA,miRNA_seq) not in added:
                    data.write(gene+'\t'+miRNA+'\t'+miRNA_seq+'\n'); added[(gene,miRNA,miRNA_seq)]=[]
            if gene in splice_event_db:
                for ed in splice_event_db[gene]:
                    probeset = ed.Probeset()
                    probeset_seq = ed.ExonSeq()
                    #exonid = ed.ExonID()
                    proceed = 'no'; hit = 'no'
                    if len(miRNA_seq)>0 and len(probeset_seq)>0:
                        miRNA_seqs = string.split(miRNA_seq,'|') ### Can be multiples
                        for mseq in miRNA_seqs:
                            if mseq in probeset_seq:
                                hit = 'yes'
                                if stringency == 'strict':
                                    if '|' in sources:  proceed='yes'
                                else: proceed='yes'
                        if proceed == 'yes' and hit == 'yes':
                            yc = getMicroRNAGenomicCoordinates(ed,y) ### Add's genomic coordinate information to these objects
                            try: probeset_miRNA_db[probeset].append(yc)
                            except KeyError: probeset_miRNA_db[probeset] = [yc]
    if mir_source == 'pictar': data.close()
    
    probeset_miRNA_db = eliminate_redundant_dict_values(probeset_miRNA_db)
    if stringency == 'lax': export_type = 'any'
    else: export_type = 'multiple'
    output_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_probeset_microRNAs_'+export_type+'.txt'
    coord_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_probeset_coord_microRNAs_'+export_type+'.txt'
    k=0
    fn=filepath(output_file)
    cfn=filepath(coord_file)
    data = open(fn,'w'); coord_data = open(cfn,'w')
    for probeset in probeset_miRNA_db:
        for y in probeset_miRNA_db[probeset]:
            miRNA = y.MicroRNA(); miRNA_seq = y.Sequence(); sources = y.Source()
            data.write(probeset+'\t'+miRNA+'\t'+miRNA_seq+'\t'+sources+'\n')
            coord_data.write(probeset+'\t'+miRNA+'\t'+miRNA_seq+'\t'+sources+'\t'+y.GeneID()+'\t'+y.Coordinates()+'\t'+y.EnsExons()+'\n')
        k+=1
    data.close(); coord_data.close()
    print k, 'entries written to', output_file

def getMicroRNAGenomicCoordinates(ed,y):
    ### This function is used for internal QC
    miRNA = y.MicroRNA(); miRNA_seq = y.Sequence(); sources = y.Source()
    ps_start = ed.ProbeStart(); ps_end = ed.ProbeStop(); strand = ed.Strand(); chr = ed.Chromosome()
    probeset_seq = ed.ExonSeq()
    mi_start = string.find(probeset_seq,miRNA_seq)
    mi_end = mi_start+len(miRNA_seq)
    if strand == '+':
        genomic_start = ps_start+mi_start
        genomic_end = ps_start+mi_end
    else:
        genomic_start = ps_end-mi_start
        genomic_end = ps_end-mi_end
    yc = copy.deepcopy(y) ### Multiple sites for the same miRNA can exist per gene
    yc.setCoordinates(chr,strand,genomic_start,genomic_end)
    yc.setEnsExons(ed.ExternalExonIDList())

    return yc

def checkProbesetSequenceFile(species):
    """ Check to see if the probeset sequence file is present and otherwise download AltAnalyze hosted version"""
    probeset_seq_file = getProbesetSequenceFile(species)
    if probeset_seq_file == None:
        dir = '/AltDatabase/'+species+'/'+array_type
        filename = update.getFileLocations(species,'exon_seq')
        filename = dir[1:]+'/'+ filename
        update.downloadCurrentVersion(filename,'exon','.fa')
        probeset_seq_file = getProbesetSequenceFile(species)
    return probeset_seq_file

def getProbesetSequenceFile(species):
    probeset_seq_file = None
    import_dir = '/AltDatabase'+'/'+species+'/exon'
    filedir = import_dir[1:]+'/'
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for input_file in dir_list:    #loop through each file in the directory to output results
        ### No probeset sequence file currently for gene arrays, so use the same exon region probeset sequence from the exon array
        if '.probeset.fa' in input_file: probeset_seq_file = filedir+input_file
    return probeset_seq_file

def runProgram(Species,Array_type,Process_microRNA_predictions,miR_source,Stringency):
    global species; species = Species; global array_type; array_type = Array_type
    global process_microRNA_predictions; process_microRNA_predictions = Process_microRNA_predictions
    global mir_source; mir_source = miR_source 
    global stringency; stringency = Stringency
    probeset_seq_file = checkProbesetSequenceFile(species)
    
    probeset_annotations_file = 'AltDatabase/'+species+'/exon/'+species+'_Ensembl_probesets.txt'
    splice_event_db = getParametersAndExecute(probeset_seq_file,array_type,species)

    if process_microRNA_predictions == 'yes':
        print 'stringency:',stringency 
        if mir_source == 'pictar':
            ensembl_mirna_db = importmiRNATargetPredictions(species)
        else:
            ensembl_mirna_db = importmiRNATargetPredictionsAdvanced(species)
        alignmiRNAData(array_type,mir_source,species,stringency,ensembl_mirna_db,splice_event_db)
        
if __name__ == '__main__':
    species = 'Hs'; array_type = 'exon'
    process_microRNA_predictions = 'yes'
    mir_source = 'multiple'; stringency = 'strict'
    print array_type
    runProgram(species,array_type,process_microRNA_predictions,mir_source,stringency)

    stringency = 'lax'
    runProgram(species,array_type,process_microRNA_predictions,mir_source,stringency); sys.exit()
    
    print "******Select Species*******"
    print "1) Human"
    print "2) Mouse"
    print "3) Rat"
    inp = sys.stdin.readline(); inp = inp.strip()
    if inp == '1': species = 'Hs'
    if inp == '2': species = 'Mm'
    if inp == '3': species = 'Rn'
    
    print "******Source Data*******"
    print "1) Multiple miR database sources"
    print "2) PicTar"
    inp = sys.stdin.readline(); inp = inp.strip()
    if inp == '1': mir_source = 'multiple'
    if inp == '2': mir_source = 'pictar'

    print "******Analysis Stringency*******"
    print "1) Atleast two overlapping probeset-miR binding site predictions"
    print "2) Any probeset-miR binding site predictions (lax)"
    inp = sys.stdin.readline(); inp = inp.strip()
    if inp == '1': stringency = 'strict'
    if inp == '2': stringency = 'lax'
    
