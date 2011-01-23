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
import copy
import time
import export
import EnsemblImport
import JunctionArrayEnsemblRules
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

def reformatExonFile(species,type):
    if type == 'exon':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
        export_path = 'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.txt'
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
            export_title = ['probeset_id','exon_id','ensembl_gene_id','transcript_cluster_id','chromosome','strand','probeset_start','probeset_stop']
            export_title +=['affy_class','constitutive_probeset','ens_exon_ids','ens_constitutive_status','exon_region','exon-region-start(s)','exon-region-stop(s)','splice_events','splice_junctions']
            export_title = string.join(export_title,'\t')+'\n'; export_data.write(export_title)
        else:
            gene, exonid, chr, strand, start, stop, constitutive_call, ens_exon_ids, splice_events, splice_junctions = t
            if constitutive_call == 'yes': ens_constitutive_status = '1'
            else: ens_constitutive_status = '0'
            export_values = [gene+':'+exonid, exonid, gene, '', chr, strand, start, stop, '', constitutive_call, ens_exon_ids, ens_constitutive_status]
            export_values+= [exonid, start, stop, splice_events, splice_junctions]
            export_values = string.join(export_values,'\t')+'\n'; export_data.write(export_values)
    export_data.close()

def importExonAnnotations(species,type):
    if type == 'exon':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
    else:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_junction.txt'

    fn=filepath(filename); x=0; exon_annotation_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            gene, exonid, chr, strand, start, stop, constitutive_call, ens_exon_ids, splice_events, splice_junctions = t
            ea = EnsemblImport.ExonAnnotationsSimple(gene, ens_exon_ids, constitutive_call, exonid, splice_events, splice_junctions)
            if type == 'junction_coordinates': 
                exon1_start,exon1_stop = string.split(start,'|')
                exon2_start,exon2_stop = string.split(stop,'|')
                if strand == '-':
                    exon1_stop,exon1_start = exon1_start,exon1_stop
                    exon2_stop,exon2_start = exon2_start,exon2_stop
                try: exon_annotation_db[chr,int(exon1_stop),int(exon2_start)].append(ea)
                except KeyError: exon_annotation_db[chr,int(exon1_stop),int(exon2_start)]=[ea]              
            else:
                try: exon_annotation_db[gene].append(ea)
                except KeyError: exon_annotation_db[gene]=[ea]
    return exon_annotation_db

def exportKnownJunctionComparisons(species):
    gene_junction_db = JunctionArrayEnsemblRules.importEnsemblUCSCAltJunctions(species,'standard')
    gene_intronjunction_db = JunctionArrayEnsemblRules.importEnsemblUCSCAltJunctions(species,'_intronic')
    for i in gene_intronjunction_db: gene_junction_db[i]=[]
    
    junction_export = 'AltDatabase/' + species + '/RNASeq/'+ species + '_junction_comps.txt'
    fn=filepath(junction_export); data = open(fn,'w')
    print "Exporting",junction_export
    title = 'gene'+'\t'+'critical_exon'+'\t'+'exclusion_junction_region'+'\t'+'inclusion_junction_region'+'\t'+'exclusion_probeset'+'\t'+'inclusion_probeset'+'\t'+'data_source'+'\n'
    data.write(title); temp_list=[]
    for (gene,critical_exon,incl_junction,excl_junction) in gene_junction_db:
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
    
    ### Import the full Ensembl exon sequence (not just the probeset region) for miRNA binding site analysis   
    analysis_type = 'get_sequence'; array_type = 'RNASeq'
    dir = 'AltDatabase/'+species+'/SequenceData/chr/'+species; gene_seq_filename = dir+'_gene-seq-2000_flank.fa'
    ensembl_exon_db = EnsemblImport.import_sequence_data(gene_seq_filename,ensembl_exon_db,species,analysis_type)
    
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
    
    junction_annotation_db = importExonAnnotations(species,'junction')
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

class JunctionData:
    def __init__(self,chr,strand,exon1_stop,exon2_start,reads,junction_id,condition):
        self.chr = chr; self.strand = strand; self._chr = chr
        self.exon2_start = exon2_start
        self.reads = reads; self.junction_id = junction_id
        self.condition = condition
    def Chr(self): return self.chr
    def Strand(self): return self.strand
    def Exon1Stop(self): return self.exon1_stop
    def Exon2Start(self): return self.exon2_start
    def Reads(self): return self.reads
    def JunctionID(self): return self.junction_id
    def Condition(self): return self.condition
    def __repr__(self): return "JunctionData values"
    
def importBEDFile(filename,species):
    fn=filepath(filename); x=0; neg_count=0; pos_count=0; junction_db={}
    condition = export.findFilename(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            format_description = data; x=1
            algorithm = 'Unknown'
            if 'TopHat' in format_description: algorithm = 'TopHat'
            if 'Splice Junctions' in format_description: algorithm = 'HMMSplicer'
        else:
            chr, exon1_start, exon2_stop, junction_id, reads, strand, null, null, null, null, lengths, null = t
            exon1_len,exon2_len=string.split(lengths,',')[:2]; exon1_len = int(exon1_len); exon2_len = int(exon2_len)
            exon1_start = int(exon1_start); exon2_stop = int(exon2_stop)
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
            if proceed == 'yes':
                if strand == '+': pos_count+=1
                else: neg_count+=1
                ji = JunctionData(chr,strand,exon1_stop,exon2_start,reads,junction_id,condition)
                try: junction_db[chr,exon1_stop,exon2_start].append(ji)
                except Exception: junction_db[chr,exon1_stop,exon2_start] = [ji]

    print len(junction_db),'junctions present in bed file, ('+str(pos_count),str(neg_count)+') by strand'
    return junction_db

def getEnsemblAssociations(species,data_type,test_status,force):
    ### Get UCSC associations (download databases if necessary)
    import UCSCImport
    mRNA_Type = 'mrna'; run_from_scratch = 'yes'
    export_all_associations = 'no' ### YES only for protein prediction analysis
    UCSCImport.runUCSCEnsemblAssociations(species,mRNA_Type,export_all_associations,run_from_scratch,force)
    
    null = EnsemblImport.getEnsemblAssociations(species,data_type,test_status); null=[]
    reformatExonFile(species,'exon'); reformatExonFile(species,'junction')
    exportKnownJunctionComparisons(species)
    getExonAndJunctionSequences(species)

def alignJunctionsToEnsembl(species,dir):
    ens_junction_coord_db = importExonAnnotations(species,'junction_coordinates')
    print len(ens_junction_coord_db),'Ensembl/UCSC junctions imported'
    junction_db = importBEDFile(dir,species)

    count=0    
    for key in junction_db:
        if key in ens_junction_coord_db: count +=1
    print count,  'junctions found in Ensembl/UCSC and',len(junction_db)-count,'are novel.'
    
if __name__ == '__main__':
    species = 'Hs'
    test_status = 'yes'
    data_type = 'ncRNA'
    data_type = 'mRNA'
    alignJunctionsToEnsembl(species,'canonical.bed'); sys.exit()
    getEnsemblAssociations(species,data_type,test_status,'yes'); sys.exit()
