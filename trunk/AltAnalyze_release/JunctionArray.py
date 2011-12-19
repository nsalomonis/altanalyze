###JunctionArray
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
import math
import reorder_arrays
import ExonArray
import EnsemblImport
import ExonArrayEnsemblRules
import JunctionArrayEnsemblRules
import export
import RNASeq

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def verifyFile(filename):
    fn=filepath(filename)
    try:
        for line in open(fn,'rU').xreadlines():break
    except Exception:
        import update; reload(update); server_folder = 'AltMouse'
        update.downloadCurrentVersion(filename,server_folder,'txt')

########### Recent code for dealing with comprehensive Affymetrix Junction Arrays

########### Begin Analyses ###########
class ExonAnnotationData:
    def Probeset(self): return self._probeset
    def ProbesetName(self): return self._psr
    def ExonClusterID(self): return self._exon_cluster_id
    def EnsemblGeneID(self):
        geneid = self._geneid
        if 'ENS' in self._geneid:
            if ',' in self._geneid:
                ens=[]
                ids = string.split(self._geneid,',')
                for id in ids:
                    if 'ENS' in id: ens.append(id)
                geneid = unique.unique(ens)[-1]
        else: geneid=''
        return geneid
    def Symbol(self):
        try: symbols = string.split(self._symbols,',')
        except Exception: print self._probeset,self._symbols;kill
        return symbols
    def TranscriptCluster(self): return self._transcript_cluster
    def ProbesetType(self):
        ###e.g. Exon, junction, constitutive(gene)
        return self._probeset_type
    def Chromosome(self):
        if len(self._chromosome_info)>0:
            null,chr = string.split(self._chromosome_info,'=chr')
            chromosome,null=string.split(chr,':')
        else: chromosome = 'not-assinged'
        return chromosome
    def Strand(self):
        if self._strand == '-': self._strand = '-1'
        else: self._strand = '1'
        return self._strand
    def ProbesetClass(self): 
        ###e.g. core, extendended, full
        #return self._probest_class
        return 'core'
    def ExternalExonClusterIDs(self): return self._exon_clusters
    def ExternalExonClusterIDList(self):
        external_exonid_list = string.split(self.ExternalExonClusterIDs(),'|')
        return external_exonid_list
    def Constitutive(self): return self._constitutive_status
    def Sequence(self): return string.lower(self._seq)
    def JunctionSequence(self): return string.replace(self.Sequence(),'|','')
    def JunctionSequences(self):
        seq1, seq2 = string.split(self.Sequence(),'|')
        return seq1, seq2
    def Report(self):
        output = self.Probeset()+'|'+self.ExternalGeneID()
        return output
    def __repr__(self): return self.Report()

class PSRAnnotation(ExonAnnotationData):
    def __init__(self,psr,probeset,ucsclink,transcript_cluster,strand,geneids,symbols,exon_clusters,constitutive,seq,probeset_type):
        self._transcript_cluster = transcript_cluster; self._geneid = geneids; self._exon_clusters = exon_clusters;
        self._constitutive_status = constitutive; self._symbols = symbols
        self._strand = strand; self._chromosome_info = ucsclink; self._probeset = probeset; self._psr = psr; self._seq = seq
        self._probeset_type = probeset_type

class EnsemblInformation:
    def __init__(self, chr, strand, gene, symbol, description):
        self._chr = chr; self._strand = strand; self._gene = gene; self._description = description
        self._symbol = symbol
    def GeneID(self): return self._gene
    def Chromosome(self): return self._chr
    def Strand(self): return self._strand
    def Description(self): return self._description
    def Symbol(self): return self._symbol
    def __repr__(self): return self.GeneID()
    
def importJunctionArrayAnnotations(species,array_type):

    import EnsemblImport
    ens_gene_chr_db = EnsemblImport.importEnsGeneData(species) ### retrieves chromosome and strand info for each gene
    ensembl_annotations = 'AltDatabase/ensembl/'+ species + '/'+species+ '_Ensembl-annotations_simple.txt'
    ensembl_annotation_db = importGeneric(ensembl_annotations)
    
    extraction_type = 'Ensembl'
    tc_ensembl_annotations = importJunctionArrayAnnotationMappings(array_type,species,extraction_type)
        
    ensembl_symbol_db={}; ensembl_gene_db={}
    for ens_geneid in ensembl_annotation_db:
        description, symbol = ensembl_annotation_db[ens_geneid]
        if ens_geneid in ens_gene_chr_db:
            chr,strand = ens_gene_chr_db[ens_geneid]
            ei = EnsemblInformation(chr,strand,ens_geneid,symbol,description)
            if len(symbol)>0:
                try: ensembl_symbol_db[symbol].append(ei)
                except KeyError: ensembl_symbol_db[symbol] =[ei]
                ensembl_gene_db[ens_geneid] = ei

    primary_gene_annotation_export = 'AltDatabase/'+species +'/'+ array_type +'/'+ array_type+ '_gene_annotations.txt'    
    ens_match=0; sym_match=0; ensembl_associations={}; gene_annotation_db={}         
    for transcript_cluster_id in tc_ensembl_annotations:
        ti = tc_ensembl_annotations[transcript_cluster_id]
        ti.Strand(); ti.Chromosome(); ti.Symbol(); ens_geneids=[]
        if ti.EnsemblGeneID() in ensembl_gene_db: ens_geneids = [ti.EnsemblGeneID()]; ens_match+=1
        else:
            for symbol in ti.Symbol():
                if symbol in ensembl_symbol_db:
                    for ei in ensembl_symbol_db[symbol]:
                        #print [ei.Chromosome(),ti.Chromosome(),ei.Strand(),ti.Strand()];kill
                        if ei.Chromosome() == ti.Chromosome() and ei.Strand() == ti.Strand():
                            ens_geneids.append(ei.GeneID()); sym_match+=1
        if len(ens_geneids)>0:
            for ens_geneid in ens_geneids:
                try: ensembl_associations[transcript_cluster_id].append(ens_geneid)
                except KeyError: ensembl_associations[transcript_cluster_id] = [ens_geneid]
                ei = ensembl_gene_db[ens_geneid]
                gene_annotation_db[transcript_cluster_id]=[ei.Description(),ens_geneid,ei.Symbol(),'']
    exportDB(primary_gene_annotation_export,gene_annotation_db)
    ensembl_associations = eliminate_redundant_dict_values(ensembl_associations)
    print ens_match, 'direct Ensembl-Ensembl gene mapping and', sym_match, 'indirect Symbol-chromosome mapping'
    print len(tc_ensembl_annotations)-len(ensembl_associations),'unmapped transcript clusters'
    """
    u=0 ### print transcript clusters without gene IDs
    for i in tc_ensembl_annotations:
        if i not in ensembl_associations:
            if u<15:
                print i, tc_ensembl_annotations[i].EnsemblGeneID(); u+=1
                """
    
    exportArrayIDEnsemblAssociations(ensembl_associations,species,array_type) ###Use these For LinkEST program
    return ensembl_associations

def pickShortestExonIDDiff(exon_to_exon):
    if '|' not in exon_to_exon:
        five_exon,three_exon=string.split(exon_to_exon,'_to_')
        return five_exon,three_exon
    else:
        exon_comps = string.split(exon_to_exon,'|'); diff_list=[]
        for exon_comp in exon_comps:
            five_exon,three_exon=string.split(exon_comp,'_to_')
            try: diff=abs(int(five_exon[5:])-int(three_exon[5:]))
            except Exception: print five_exon, three_exon;kill
            diff_list.append((diff,[five_exon,three_exon]))
        diff_list.sort()
        return diff_list[0][1]
    
def importJunctionArrayAnnotationMappings(array_type,species,extraction_type):
    print 'Importing junction array sequence mapping'
    export_dir = 'AltDatabase/'+species+'/'+array_type+'/'
    filename = export_dir+string.lower(species[0])+'jay.r2.annotation_map'
    fn=filepath(filename)

    if extraction_type == 'sequence':
        probeset_junctionseq_export = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'_critical-junction-seq.txt'
        fn2=filepath(probeset_junctionseq_export); dw = open(fn2,'w'); print "Exporting",probeset_junctionseq_export
    
    probeset_translation_db={}; x=0; tc=0; j=0; p=0; k=0; tc_db=(); transcript_cluster_count={}; transcript_cluster_count2={}
    global probeset_db; global junction_comp_db; junction_comp_db={}; global junction_alinging_probesets
    ps_db={}; jc_db={}; left_ec={}; right_ec={}; psr_ec={}; probeset_db={}; junction_alinging_probesets={}; nonconstitutive_junctions={}
    for line in open(fn,'r').xreadlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x<5: x+=1
        elif x>2:
            if len(t)==15: ###Transcript Cluster ID Lines
                probeset, probeset_name, ucsclink, transcript_cluster, genome_pos, strand, transcripts, geneids, symbols, descriptions, TR_count, JUC_count, PSR_count, EXs_count, ECs_count = t
                type = 'TranscriptCluster'; seq=''; exon_cluster=''; constitutive=''
                psri = PSRAnnotation(probeset_name,probeset,ucsclink,transcript_cluster,strand,geneids,symbols,exon_cluster,constitutive,seq,type)
            elif len(t)==28:###Junction ID Lines
                probeset, probeset_name, ucsclink, transcript_cluster, genome_pos, strand, transcripts, geneids, symbols, descriptions, TR_count, JUC_count, PSR_count, EXs_count, ECs_count, junction_number, original_seq, exon_to_exon, observed_speculative, strand, five_PSR, three_PSR, five_EC, three_EC, Rel_5EC, Rel_3EC, constitutive, blat_junction = t
                type = 'Junction'; exon_cluster = [five_EC,three_EC]
                if constitutive == 'alternative': nonconstitutive_junctions[probeset]=[]
                five_exon,three_exon = pickShortestExonIDDiff(exon_to_exon)
                junction_alinging_probesets[probeset] = [five_PSR,five_exon], [three_PSR,three_exon]; seq = blat_junction
                psri = PSRAnnotation(probeset_name,probeset,ucsclink,transcript_cluster,strand,geneids,symbols,exon_cluster,constitutive,seq,type)
            elif len(t)==24: ###Probeset ID Lines
                probeset, probeset_name, ucsclink, transcript_cluster, genome_pos, strand, transcripts, geneids, symbols, descriptions, TR_count, JUC_count, PSR_count, EXs_count, ECs_count, PSR_region, genome_pos2, strand, exon_cluster, constitutive, TR_hits, percent_TR_hits, location_5to3_percent,seq = t
                type = 'Exon'
                psri = PSRAnnotation(probeset_name,probeset,ucsclink,transcript_cluster,strand,geneids,symbols,exon_cluster,constitutive,seq,type)
            else:
                #if k<40 and len(t)>5: print len(t),t; k+=1
                type = 'null'
            ### Exon clusters are equivalent to exon blocks in this schema and can be matched between junctions and exons
            store = 'yes'
            if extraction_type == 'Ensembl':
                if type != 'TranscriptCluster': store = 'no'
            elif extraction_type == 'sequence':
                store = 'no'
                if type == 'Exon' or type == 'Junction':
                    transcript_cluster_count[psri.TranscriptCluster()]=[]
                    if psri.TranscriptCluster() in ensembl_associations:
                        ens_geneid = ensembl_associations[psri.TranscriptCluster()][0]
                        critical_junctions=''
                        if type == 'Junction':
                            dw.write(probeset+'\t'+psri.JunctionSequence()+'\t\t\n')
                            seq = psri.JunctionSequences()[0]; exon_id = probeset+'|5'
                            seq_data = ExonSeqData(exon_id,psri.TranscriptCluster(),psri.TranscriptCluster()+':'+exon_id,critical_junctions,seq)
                            try: probeset_db[ens_geneid].append(seq_data)
                            except Exception: probeset_db[ens_geneid] = [seq_data]  
                            
                            seq = psri.JunctionSequences()[1]; exon_id = probeset+'|3'
                            seq_data = ExonSeqData(exon_id,psri.TranscriptCluster(),psri.TranscriptCluster()+':'+exon_id,critical_junctions,seq)
                            try: probeset_db[ens_geneid].append(seq_data)
                            except Exception: probeset_db[ens_geneid] = [seq_data]
                            transcript_cluster_count2[psri.TranscriptCluster()]=[]
                        elif type == 'Exon':
                            dw.write(probeset+'\t'+psri.Sequence()+'\t\t\n')
                            seq = psri.Sequence(); exon_id = probeset
                            seq_data = ExonSeqData(exon_id,psri.TranscriptCluster(),psri.TranscriptCluster()+':'+exon_id,critical_junctions,seq)
                            try: probeset_db[ens_geneid].append(seq_data)
                            except Exception: probeset_db[ens_geneid] = [seq_data]
                            transcript_cluster_count2[psri.TranscriptCluster()]=[]
            if store == 'yes':
                try: probeset_db[probeset] = psri
                except Exception: null=[]
            if type == 'TranscriptCluster':
                tc+=1
            if type == 'Junction':
                j+=1
                if extraction_type == 'comparisons':
                    ### Store the left exon-cluster and right exon-cluster for each junction
                    try: left_ec[five_EC].append(probeset)
                    except KeyError: left_ec[five_EC]=[probeset]
                    try: right_ec[three_EC].append(probeset)
                    except KeyError: right_ec[three_EC]=[probeset]
            if type == 'Exon':
                p+=1
                if extraction_type == 'comparisons':
                    try: psr_ec[exon_cluster].append(probeset)
                    except KeyError: psr_ec[exon_cluster]=[probeset]
            """
            print 'psid',psid; print 'probeset',probeset; print 'ucsclink',ucsclink
            print 'transcript_cluster',transcript_cluster; print 'transcripts',transcripts
            print 'geneids',geneids; print 'symbols',symbols; print 'seq',seq; kill"""
            x+=1
    print 'TCs:',tc, 'Junctions:',j, 'Exons:',p, 'Total:',x

    #print 'JUC0900017373',probeset_db['JUC0900017373'].Sequence()
    #print 'JUC0900017385',probeset_db['JUC0900017385'].Sequence();kill

    if extraction_type == 'sequence':
        dw.close()
        print len(probeset_db),'Entries exported from Junction annotation file'
        return probeset_db
    if extraction_type == 'Ensembl':
        print len(probeset_db),'Entries exported from Junction annotation file'
        return probeset_db 
    if extraction_type == 'comparisons':
        global junction_inclusion_db; global ensembl_exon_db; global exon_gene_db
        junction_inclusion_db = JunctionArrayEnsemblRules.reimportJunctionComps(species,array_type,'original')
        ensembl_exon_db,exon_gene_db = JunctionArrayEnsemblRules.importAndReformatEnsemblJunctionAnnotations(species,array_type,nonconstitutive_junctions)
        
        global failed_db; failed_db={}
        global passed_db; passed_db={}
        print len(junction_inclusion_db)
        identifyCompetitiveJunctions(right_ec,"3'")
        identifyCompetitiveJunctions(left_ec,"5'")
        print 'len(passed_db)',len(passed_db),'len(failed_db)',len(failed_db)
        print 'len(junction_inclusion_db)',len(junction_inclusion_db)
        exportUpdatedJunctionComps(species,array_type)
    
def exportUpdatedJunctionComps(species,array_type):
    db_version = unique.getCurrentGeneDatabaseVersion() ### Only need this since we are exporting to root_dir for RNASeq
    if array_type == 'RNASeq': species,root_dir=species
    else: root_dir = ''
    probeset_junction_export = root_dir+'AltDatabase/'+db_version+'/'+ species + '/'+array_type+'/'+ species + '_junction_comps_updated.txt'; lines_exported=0
    
    if array_type == 'RNASeq':
        data,status = RNASeq.AppendOrWrite(probeset_junction_export) ### Creates a new file or appends if already existing (import is chromosome by chromosome)
    else:
        data = export.ExportFile(probeset_junction_export); status = 'not found'
    
    if array_type != 'RNASeq': print "Exporting",probeset_junction_export
    if status == 'not found':
        title = 'gene'+'\t'+'critical_exon'+'\t'+'exclusion_junction_region'+'\t'+'inclusion_junction_region'+'\t'+'exclusion_probeset'+'\t'+'inclusion_probeset'+'\t'+'data_source'+'\n'
        data.write(title)
    for i in junction_inclusion_db:
        critical_exons=[]
        for ji in junction_inclusion_db[i]:
            #value = string.join([ji.GeneID(),ji.CriticalExon(),ji.ExclusionJunction(),ji.InclusionJunction(),ji.ExclusionProbeset(),ji.InclusionProbeset(),ji.DataSource()],'\t')+'\n'
            ### Combine all critical exons for a probeset pair
            critical_exons.append(ji.CriticalExon())
        critical_exons = unique.unique(critical_exons); critical_exons = string.join(critical_exons,'|'); ji.setCriticalExons(critical_exons); lines_exported+=1
        data.write(ji.OutputLine())
    data.close()
    if array_type != 'RNASeq':
        print lines_exported,'for',probeset_junction_export
    
    
def identifyCompetitiveJunctions(exon_cluster_db,junction_type):
    """To identify critical exons (e.g., the alternatively spliced exon sequence for two alternative exon-junctions), this script:
    1) Finds pairs of junctions that contain the same 5' or 3' exon-cluster (genomic overlapping transcript exons)
    2) Determines which junction has exons that are closes in genomic space, between the pair of junctions (based on exon-cluster ID number or exon ID)
    3) Selects the non-common exon and stores the junction sequence for that exon
    4) Selects any exon probeset ID that is annotated as overlapping with the critical exon

    The greatest assumption with this method is that the critical exon is choosen based on the numerical ID in the exon-cluster or exon ID (when the exon-clusters
    between the two junctions are the same). For example looked at, this appears to be true (e.g., two exons that make up a junction have a difference of 1 in their ID),
    but this may not always be the case. Ideally, this method is more extensively tested by evaluating junction and exon sequences mapped to genomic coordinates
    and AltAnalyze exon block and region coordinates to verify the critical exon selection."""
    passed=0; failed=0; already_added=0
    if junction_type == "5'": index = 1
    else: index = 0
    for ec in exon_cluster_db:
        if len(exon_cluster_db[ec])>1:
            junction_comps={} ### Calculate all possible pairwise-junction comparisons
            for junction1 in exon_cluster_db[ec]:
                for junction2 in exon_cluster_db[ec]:
                    if junction1 != junction2: temp = [junction1,junction2]; temp.sort(); junction_comps[tuple(temp)]=[]

            for (junction1,junction2) in junction_comps:
                store_data = 'no'
                if (junction1,junction2) in junction_inclusion_db or (junction2,junction1) in junction_inclusion_db:
                    already_added+=1
                elif junction1 in ensembl_exon_db and junction2 in ensembl_exon_db: ### Thus, these are mapped to the genome
                    ed1 = ensembl_exon_db[junction1]; ed2 = ensembl_exon_db[junction2]
                    ensembl_gene_id = ed1.GeneID()
                    diff1 = ed1.JunctionDistance(); diff2 = ed2.JunctionDistance()
                    ### Using the ranked exon-cluster IDs
                    psri1 = probeset_db[junction1]; exon1a = psri1.ExternalExonClusterIDs()[0]; exon1b = psri1.ExternalExonClusterIDs()[-1]
                    psri2 = probeset_db[junction2]; exon2a = psri2.ExternalExonClusterIDs()[0]; exon2b = psri2.ExternalExonClusterIDs()[-1]
                    diffX1 = abs(int(exon1a[5:])-int(exon1b[5:])); diffX2 = abs(int(exon2a[5:])-int(exon2b[5:]))
                    junction1_exon_id = ed1.ExonID(); junction2_exon_id = ed2.ExonID()
                    if diffX1==0 or diffX2==0: null=[] ### splicing occurs within a single exon-cluster
                    elif diff1<diff2: ### Thus the first junction contains the critical exon
                        #critical_exon_seq = psri1.JunctionSequences()[index] ### if left most exon in junction is common, then choose the most proximal right exon as critical
                        incl_junction_probeset = junction1; excl_junction_probeset = junction2
                        incl_junction_id = junction1_exon_id; excl_junction_id = junction2_exon_id
                        incl_exon_probeset,incl_exon_id = junction_alinging_probesets[junction1][index]
                        store_data = 'yes'
                    elif diff2<diff1:
                        incl_junction_probeset = junction2; excl_junction_probeset = junction1
                        incl_junction_id = junction2_exon_id; excl_junction_id = junction1_exon_id
                        incl_exon_probeset,incl_exon_id = junction_alinging_probesets[junction2][index]
                        store_data = 'yes'
                    if store_data == 'yes':
                        critical_exon_id = string.split(incl_junction_id,'-')[index]; critical_exon_id = string.replace(critical_exon_id,'.','-')
                        if incl_exon_probeset in ensembl_exon_db:
                            if (excl_junction_probeset,incl_exon_probeset) in junction_inclusion_db or (incl_exon_probeset,excl_junction_probeset) in junction_inclusion_db:
                                already_added+=1
                            else:
                                critical_exon_id = ensembl_exon_db[incl_exon_probeset]
                                ji=JunctionArrayEnsemblRules.JunctionInformation(ensembl_gene_id,critical_exon_id,excl_junction_id,critical_exon_id,excl_junction_probeset,incl_exon_probeset,'Affymetrix')
                                try: junction_inclusion_db[excl_junction_probeset,incl_exon_probeset].append(ji)
                                except Exception: junction_inclusion_db[excl_junction_probeset,incl_exon_probeset] = [ji]
                                #value = string.join([ji.GeneID(),ji.CriticalExon(),ji.ExclusionJunction(),ji.InclusionJunction(),ji.ExclusionProbeset(),ji.InclusionProbeset(),ji.DataSource()],'\t')+'\n'
                                #print ji.OutputLine();kill
                                #print [[critical_exon_id,junction2,ed2.ExonID(), ed1.JunctionCoordinates(), ed2.JunctionCoordinates(), diff1,diff2]]

                        passed+=1
                        passed_db[junction1,junction2]=[]
                        ji=JunctionArrayEnsemblRules.JunctionInformation(ensembl_gene_id,critical_exon_id,excl_junction_id,incl_junction_id,excl_junction_probeset,incl_junction_probeset,'Affymetrix')
                        #print ensembl_gene_id,critical_exon_id,excl_junction_id,incl_junction_id,excl_junction_probeset,incl_junction_probeset;kill
                        #print [critical_exon_id,junction1,junction2,ed1.ExonID(),ed2.ExonID(), ed1.JunctionCoordinates(), ed2.JunctionCoordinates(), diff1,diff2]
                        try: junction_inclusion_db[exclusion_junction,inclusion_junction].append(ji)
                        except Exception: junction_inclusion_db[excl_junction_probeset,incl_junction_probeset] = [ji]

    print 'already_added:',already_added,'passed:',passed,'failed:',failed

def identifyJunctionComps(species,array_type):
    JunctionArrayEnsemblRules.getJunctionComparisonsFromExport(species,array_type)
    
    extraction_type = 'comparisons'
    tc_ensembl_annotations = importJunctionArrayAnnotationMappings(array_type,species,extraction_type)

    inferJunctionComps(species,array_type)
    
def filterForCriticalExons(species,array_type):
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_'+array_type+'_probesets.txt' 
    importForFiltering(species,array_type,filename,'exclude_junction_psrs')
    
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_probeset_microRNAs_any.txt' 
    importForFiltering(species,array_type,filename,'include_critical_exon_ids')

    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_probeset_microRNAs_multiple.txt' 
    importForFiltering(species,array_type,filename,'include_critical_exon_ids')

    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_domain_aligning_probesets.txt' 
    importForFiltering(species,array_type,filename,'include_critical_exon_ids')

    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_indirect_domain_aligning_probesets.txt' 
    importForFiltering(species,array_type,filename,'include_critical_exon_ids')

    filename = 'AltDatabase/'+species+'/'+array_type+'/exon/probeset-protein-annotations-exoncomp.txt' 
    importForFiltering(species,array_type,filename,'exclude_critical_exon_ids')

    filename = 'AltDatabase/'+species+'/'+array_type+'/exon/probeset-domain-annotations-exoncomp.txt' 
    importForFiltering(species,array_type,filename,'exclude_critical_exon_ids')
    
def importForFiltering(species,array_type,filename,export_type):
    fn=filepath(filename); dbase={}; x = 0
    print 'Filtering:',filename
    dbase['filename'] = filename
    ###Import expression data (non-log space)
    for line in open(fn,'r').xreadlines():             
        data = cleanUpLine(line); splitup = 'no'
        if x == 0: x=1; dbase['title'] = line; x+=1   ###Grab expression values
        if x !=0:
            key = string.split(data,'\t')[0]
            if ':' in key:
                old_key = key
                key = string.split(key,':')[1]
                line = string.replace(line,old_key,key)
            if '|' in key: ### Get rid of |5 or |3
                line = string.replace(line,key,key[:-2])
                if export_type == 'exclude_critical_exon_ids': splitup = 'yes'
            if splitup == 'no':
                try: dbase[key].append(line)
                except Exception: dbase[key] = [line]
    #print len(dbase)
    filterExistingFiles(species,array_type,dbase,export_type)
    
def importGenericAppend(filename,key_db):
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        key_db[t[0]] = t[1:]
    return key_db
   
def importGenericAppendDBList(filename,key_db):
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        try: key_db[t[0]].append(t[1])
        except KeyError:  key_db[t[0]] = [t[1]]
    return key_db

def combineExonJunctionAnnotations(species,array_type):
    ###Currently used for RNASeq databases to minimize the number of files supplied to the user
    collapseSequenceFiles(species,array_type)
    overRideJunctionEntriesWithExons(species,array_type)
    collapseDomainAlignmentFiles(species,array_type,species+'_Ensembl_domain_aligning_probesets.txt')
    collapseDomainAlignmentFiles(species,array_type,species+'_Ensembl_indirect_domain_aligning_probesets.txt')
    
def collapseDomainAlignmentFiles(species,array_type,filename):
    original_filename = 'AltDatabase/'+species+'/'+array_type+'/'+filename
    domain_db = importGenericAppendDBList(original_filename,{})
    
    filename = 'AltDatabase/'+species+'/'+array_type+'/junction/'+filename
    domain_db = importGenericAppendDBList(filename,domain_db); del domain_db['Probeset']
    header = 'Probeset\tInterPro-Description\n'
    exportGenericList(domain_db,original_filename,header)
    
def exportGenericList(db,filename,header):
    data_export = export.ExportFile(filename)
    if len(header)>0: data_export.write(header)
    print 'Re-writing',filename
    for key in db:
        for i in db[key]: data_export.write(string.join([key]+[i],'\t')+'\n')
    data_export.close()
    
def collapseSequenceFiles(species,array_type):
    original_filename = 'AltDatabase/'+species+'/'+array_type+'/SEQUENCE-protein-dbase_exoncomp.txt'
    seq_db = importGenericAppend(original_filename,{})
    
    filename = 'AltDatabase/'+species+'/'+array_type+'/exon/SEQUENCE-protein-dbase_exoncomp.txt'
    try: seq_db = importGenericAppend(filename,seq_db)
    except Exception: print 'SEQUENCE-protein-dbase_exoncomp.txt - exon version not found'
    
    filename = 'AltDatabase/'+species+'/'+array_type+'/junction/SEQUENCE-protein-dbase_exoncomp.txt'
    try: seq_db = importGenericAppend(filename,seq_db)
    except Exception: print 'SEQUENCE-protein-dbase_exoncomp.txt - junction version not found'
    exportGeneric(seq_db,original_filename,[])
    
def exportGeneric(db,filename,header):
    data_export = export.ExportFile(filename)
    if len(header)>0: data_export.write(header)
    print 'Re-writing',filename
    for key in db:
        data_export.write(string.join([key]+db[key],'\t')+'\n')
    data_export.close()
    
def overRideJunctionEntriesWithExons(species,array_type):
    filename1 = 'AltDatabase/'+species+'/'+array_type+'/exon/probeset-protein-annotations-exoncomp.txt'
    filename2 = 'AltDatabase/'+species+'/'+array_type+'/junction/probeset-protein-annotations-exoncomp.txt'
    overRideExistingEntries(filename1,filename2)
    
    filename1 = 'AltDatabase/'+species+'/'+array_type+'/exon/probeset-domain-annotations-exoncomp.txt' 
    filename2 = 'AltDatabase/'+species+'/'+array_type+'/junction/probeset-domain-annotations-exoncomp.txt'
    overRideExistingEntries(filename1,filename2)
 
def overRideExonEntriesWithJunctions(species,array_type):
    filename1 = 'AltDatabase/'+species+'/'+array_type+'/junction/'+species+'_Ensembl_domain_aligning_probesets.txt' 
    filename2 = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_domain_aligning_probesets-filtered.txt'
    overRideExistingEntries(filename1,filename2)
    
    filename1 = 'AltDatabase/'+species+'/'+array_type+'/junction/'+species+'_Ensembl_indirect_domain_aligning_probesets.txt'
    filename2 = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_indirect_domain_aligning_probesets-filtered.txt' 
    overRideExistingEntries(filename1,filename2)
    
    filename1 = 'AltDatabase/'+species+'/'+array_type+'/junction/probeset-protein-annotations-exoncomp.txt'
    filename2 = 'AltDatabase/'+species+'/'+array_type+'/exon/probeset-protein-annotations-exoncomp-filtered.txt'
    overRideExistingEntries(filename1,filename2)
    
    filename1 = 'AltDatabase/'+species+'/'+array_type+'/junction/probeset-domain-annotations-exoncomp.txt'
    filename2 = 'AltDatabase/'+species+'/'+array_type+'/exon/probeset-domain-annotations-exoncomp-filtered.txt' 
    overRideExistingEntries(filename1,filename2)
 
def overRideExistingEntries(file_include,file_exclude):
    ### Imports two files and over-rides entries in one with another
    
    ### These are the filtered entries to replace
    fn=filepath(file_include); dbase_include={}; x = 0
    for line in open(fn,'r').xreadlines():             
        data = cleanUpLine(line)
        key = string.split(data,'\t')[0]
        try: dbase_include[key].append(line)
        except Exception: dbase_include[key] = [line]
        x+=1 
    print x;title=''
    fn=filepath(file_exclude); dbase_exclude={}; x = 0
    for line in open(fn,'r').xreadlines():             
        data = cleanUpLine(line)
        if x == 0: x=1; title = line; x+=1 
        if x != 0:
            key = string.split(data,'\t')[0]
            try: dbase_exclude[key].append(line)
            except Exception: dbase_exclude[key] = [line]
            x+=1 
    print x
    count=0
    for key in dbase_exclude: count+=1
    print file_exclude, count
    
    count=0
    for key in dbase_include:
        dbase_exclude[key] = dbase_include[key]
        count+=1
    print file_exclude, count
    dbase_exclude = eliminate_redundant_dict_values(dbase_exclude)
    data_export = export.ExportFile(file_exclude)
    count=0
    print 'Re-writing',file_exclude,'with junction aligned entries.'
    try: data_export.write(title)
    except Exception: null=[] ### Occurs when no alternative isoforms present for this genome
    for key in dbase_exclude:
        for line in dbase_exclude[key]:
            data_export.write(line); count+=1
    data_export.close()
    print count

def clearObjectsFromMemory(db_to_clear):
    db_keys={}
    for key in db_to_clear: db_keys[key]=[]
    for key in db_keys:
        try: del db_to_clear[key]
        except Exception: 
            try:
                for i in key: del i ### For lists of tuples
            except Exception: del key ### For plain lists
            
class JunctionInformationSimple:
    def __init__(self,critical_exon,excl_junction,incl_junction,excl_probeset,incl_probeset):
        self._critical_exon = critical_exon; self.excl_junction = excl_junction; self.incl_junction = incl_junction
        self.excl_probeset = excl_probeset; self.incl_probeset = incl_probeset
        #self.critical_exon_sets = string.split(critical_exon,'|')
        self.critical_exon_sets = [critical_exon]
    def CriticalExon(self):
        ce = str(self._critical_exon)
        if '-' in ce: ce = string.replace(ce,'-','.')
        return ce
    def CriticalExonList(self):
        critical_exon_str = self.CriticalExon()
        critical_exons = string.split(critical_exon_str,'|')
        return critical_exons
    def setCriticalExons(self,critical_exons): self._critical_exon = critical_exons
    def setCriticalExonSets(self,critical_exon_sets): self.critical_exon_sets = critical_exon_sets
    def setInclusionProbeset(self,incl_probeset): self.incl_probeset = incl_probeset
    def setInclusionJunction(self,incl_junction): self.incl_junction = incl_junction
    def CriticalExonSets(self): return self.critical_exon_sets ### list of critical exons (can select any or all for functional analysis)
    def InclusionJunction(self): return self.incl_junction
    def ExclusionJunction(self): return self.excl_junction
    def InclusionProbeset(self): return self.incl_probeset
    def ExclusionProbeset(self): return self.excl_probeset
    def setInclusionLookup(self,incl_junction_probeset): self.incl_junction_probeset = incl_junction_probeset
    def InclusionLookup(self): return self.incl_junction_probeset
    def __repr__(self): return self.GeneID()
     
def getPutativeSpliceEvents(species,array_type,exon_db,agglomerate_inclusion_probesets,root_dir):
    alt_junction_db={}; critical_exon_db={}; critical_agglomerated={}; exon_inclusion_agglom={}; incl_junctions_agglom={}; exon_dbase={}
    exon_inclusion_db={}; comparisons=0

    ### Previously, JunctionArrayEnsemblRules.reimportJunctionComps (see above) used for import---> too memory intensive
    if array_type == 'junction': root_dir=''
    filename = root_dir+'AltDatabase/' + species + '/'+array_type+'/'+ species + '_junction_comps_updated.txt'
    fn=filepath(filename); junction_inclusion_db={}; x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line); junction_info=[]
        gene,critical_exon,excl_junction,incl_junction,excl_probeset,incl_probeset,source = string.split(data,'\t')
        """
        if gene == 'ENSG00000140464':
            a=0; b=0
            if excl_probeset in exon_db: a = 1
            if incl_probeset in exon_db: b = 1
            #print incl_probeset, a, b, excl_probeset, critical_exon
            """
        try:
            null=exon_db[excl_probeset] ### Exclusion needs to be present
            if incl_probeset in exon_db:
                ji = JunctionInformationSimple(critical_exon,excl_junction,incl_junction,excl_probeset,incl_probeset)
                junction_info.append(ji)
                #print [ji.InclusionProbeset(),ji.ExclusionProbeset()]
            if array_type == 'RNASeq':
                critical_exons = string.split(critical_exon,'|')
                for ce in critical_exons:
                    critical_exon_probeset = gene+':'+ce
                    ji=JunctionInformationSimple(ce,excl_junction,ce,excl_probeset,critical_exon_probeset)
                    junction_info.append(ji); ji.setInclusionLookup(incl_probeset) ### Use this ID to get protein and domain annotations
                    """
                    if gene == 'ENSG00000140464' and ce == 'E5.2':
                        a=0; b=0
                        if ji.ExclusionProbeset() in exon_db: a = 1
                        if ji.InclusionProbeset() in exon_db: b = 1
                        print [ji.InclusionProbeset()],a,b;kill
                        """
                    #print [ji.InclusionProbeset(),ji.ExclusionProbeset()];kill
            for ji in junction_info:
                try:
                    geneid=exon_db[ji.InclusionProbeset()].GeneID() ### This inclusion needs to be present 
                    if agglomerate_inclusion_probesets == 'yes':
                        exclProbeset = ji.ExclusionProbeset(); inclProbeset = ji.InclusionProbeset()
                        exon_inclusion_agglom[exclProbeset] = ji ### Just need one example
                        try: critical_exon_db[exclProbeset].append(ji.CriticalExon())
                        except Exception: critical_exon_db[exclProbeset]=[ji.CriticalExon()]
                        try: critical_agglomerated[exclProbeset]+=ji.CriticalExonList()
                        except Exception: critical_agglomerated[exclProbeset]=ji.CriticalExonList()
                        try: incl_junctions_agglom[exclProbeset].append(ji.InclusionJunction())
                        except Exception: incl_junctions_agglom[exclProbeset]=[ji.InclusionJunction()]
                        try: exon_inclusion_db[exclProbeset].append(inclProbeset)
                        except Exception: exon_inclusion_db[exclProbeset]=[inclProbeset]
                    else:
                        try: alt_junction_db[geneid].append(ji)
                        except Exception: alt_junction_db[geneid] = [ji]
                        comparisons+=1
                except KeyError: null=[]
        except KeyError: null=[]
        
    #print comparisons, "Junction comparisons in database"
    if agglomerate_inclusion_probesets == 'yes':
        alt_junction_agglom={}
        for excl in exon_inclusion_db:
            ji = exon_inclusion_agglom[excl]
            ed = exon_db[ji.InclusionProbeset()]; ed1 = ed
            geneid = ed.GeneID() ### If two genes are present for trans-splicing, over-ride with the one in the database
            critical_exon_sets = unique.unique(critical_exon_db[excl])
            incl_probesets = unique.unique(exon_inclusion_db[excl])
            exon_inclusion_db[excl] = incl_probesets
            critical_exons = unique.unique(critical_agglomerated[excl]); critical_exons.sort()
            incl_junctions = unique.unique(incl_junctions_agglom[excl]); incl_junctions.sort()
            ji.setCriticalExons(string.join(critical_exons,'|'))
            ji.setInclusionJunction(string.join(incl_junctions,'|'))
            ji.setInclusionProbeset(string.join(incl_probesets,'|'))
            ji.setCriticalExonSets(critical_exon_sets)
            ed1.setProbeset(string.replace(incl_probesets[0],'@',':')) ### Actually needs to be the first entry to match of re-import of a filtered list for exon_db (full not abbreviated)
            #if '|' in ji.InclusionProbeset(): print ji.InclusionProbeset(), string.replace(incl_probesets[0],'@',':');sys.exit()
            #print string.join(incl_probesets,'|'),ji.InclusionProbeset();kill 
            ### Create new agglomerated inclusion probeset entry
            #ed1.setProbeset(ji.InclusionProbeset()) ### Agglomerated probesets
            ed1.setDisplayExonID(string.join(incl_junctions,'|'))
            exon_db[ji.InclusionProbeset()] = ed1 ### Agglomerated probesets
            #if 'ENSMUSG00000032497:E23.1-E24.1' in ji.InclusionProbeset():
            #print ji.InclusionProbeset();sys.exit()
            #if '198878' in ji.InclusionProbeset(): print ji.InclusionProbeset(),excl
            try: alt_junction_db[geneid].append(ji)
            except Exception: alt_junction_db[geneid] = [ji]            

        del exon_inclusion_agglom
        critical_exon_db={}

        if array_type == 'RNASeq':
            ### Need to remove the @ from the IDs
            for e in exon_inclusion_db:
                incl_probesets=[]
                for i in exon_inclusion_db[e]:
                    incl_probesets.append(string.replace(i,'@',':'))
                exon_inclusion_db[e] = incl_probesets
            
    #clearObjectsFromMemory(junction_inclusion_db); junction_inclusion_db=[]
    critical_agglomerated=[];exon_inclusion_agglom={}; incl_junctions_agglom={}
    """ Not used for junction or RNASeq platforms
    if array_type == 'AltMouse':
        for probeset in array_id_db:
            try:
                geneid = exon_db[probeset].GeneID()  
                exons = exon_db[probeset].ExonID()
                exon_dbase[geneid,exons] = probeset
            except Exception: null=[]
    """
    #print '--------------------------------------------'
    ### Eliminate redundant entries
    objects_to_delete=[]
    for geneid in alt_junction_db:
        junction_temp_db={}; junction_temp_ls=[]
        for ji in alt_junction_db[geneid]: ### Redundant entries can be present
            id = ji.ExclusionProbeset(),ji.InclusionProbeset()
            if id in junction_temp_db: objects_to_delete.append(ji)
            else: junction_temp_db[id]=ji
        for i in junction_temp_db:
            ji = junction_temp_db[i]; junction_temp_ls.append(ji)
        alt_junction_db[geneid]=junction_temp_ls
    """
    for ji in alt_junction_db['ENSG00000140464']:
        print ji.ExclusionProbeset(), ji.InclusionProbeset(), ji.CriticalExon(), ji.ExclusionJunction(), ji.InclusionJunction()
    kill
    """
    clearObjectsFromMemory(objects_to_delete); objects_to_delete=[]
    return alt_junction_db,critical_exon_db,exon_dbase,exon_inclusion_db,exon_db

def getPutativeSpliceEventsOriginal(species,array_type,exon_db,agglomerate_inclusion_probesets,root_dir):
    junction_inclusion_db = JunctionArrayEnsemblRules.reimportJunctionComps((species,root_dir),array_type,'updated')
    alt_junction_db={}; critical_exon_db={}; critical_agglomerated={}; exon_inclusion_agglom={}; incl_junctions_agglom={}; exon_dbase={}
    exon_inclusion_db={}; comparisons=0
    
    for i in junction_inclusion_db:
        critical_exons=[]
        for ji in junction_inclusion_db[i]:
            #ji.GeneID(),ji.CriticalExon(),ji.ExclusionJunction(),ji.InclusionJunction(),ji.ExclusionProbeset(),ji.InclusionProbeset(),ji.DataSource()
            if agglomerate_inclusion_probesets == 'yes':
                if ji.InclusionProbeset() in exon_db and ji.ExclusionProbeset() in exon_db:
                    if array_type == 'RNASeq':
                        exclProbeset = ji.ExclusionProbeset(); inclProbeset=JunctionArrayEnsemblRules.formatID(ji.InclusionProbeset())
                    else: exclProbeset = ji.ExclusionProbeset(); inclProbeset = ji.InclusionProbeset()
                    exon_inclusion_agglom[exclProbeset] = ji ### Just need one example
                    try: critical_exon_db[exclProbeset].append(ji.CriticalExon())
                    except Exception: critical_exon_db[exclProbeset]=[ji.CriticalExon()]
                    try: critical_agglomerated[exclProbeset]+=ji.CriticalExonList()
                    except Exception: critical_agglomerated[exclProbeset]=ji.CriticalExonList()
                    try: incl_junctions_agglom[exclProbeset].append(ji.InclusionJunction())
                    except Exception: incl_junctions_agglom[exclProbeset]=[ji.InclusionJunction()]
                    try: exon_inclusion_db[exclProbeset].append(inclProbeset)
                    except Exception: exon_inclusion_db[exclProbeset]=[inclProbeset]
            else:
                try:
                    geneid = exon_db[ji.InclusionProbeset()].GeneID() ### If two genes are present for trans-splicing, over-ride with the one in the database
                    try: alt_junction_db[geneid].append(ji)
                    except Exception: alt_junction_db[geneid] = [ji]
                    comparisons+=1
                except Exception: geneid = ji.GeneID() ### If not in the local user datasets (don't think these genes need to be added)
                
    #print comparisons, "Junction comparisons in database"
    if agglomerate_inclusion_probesets == 'yes':
        alt_junction_agglom={}
        for excl in exon_inclusion_db:
            ji = exon_inclusion_agglom[excl]
            ed = exon_db[ji.InclusionProbeset()]; ed1 = ed
            geneid = ed.GeneID() ### If two genes are present for trans-splicing, over-ride with the one in the database
            critical_exon_sets = unique.unique(critical_exon_db[excl])
            incl_probesets = unique.unique(exon_inclusion_db[excl])
            exon_inclusion_db[excl] = incl_probesets
            critical_exons = unique.unique(critical_agglomerated[excl]); critical_exons.sort()
            incl_junctions = unique.unique(incl_junctions_agglom[excl]); incl_junctions.sort()
            ji.setCriticalExons(string.join(critical_exons,'|'))
            ji.setInclusionJunction(string.join(incl_junctions,'|'))
            ji.setInclusionProbeset(string.join(incl_probesets,'|'))
            ji.setCriticalExonSets(critical_exon_sets)
            ed1.setProbeset(string.replace(incl_probesets[0],'@',':')) ### Actually needs to be the first entry to match of re-import of a filtered list for exon_db (full not abbreviated)
            #if '|' in ji.InclusionProbeset(): print ji.InclusionProbeset(), string.replace(incl_probesets[0],'@',':');sys.exit()
            #print string.join(incl_probesets,'|'),ji.InclusionProbeset();kill 
            ### Create new agglomerated inclusion probeset entry
            #ed1.setProbeset(ji.InclusionProbeset()) ### Agglomerated probesets
            ed1.setDisplayExonID(string.join(incl_junctions,'|'))
            exon_db[ji.InclusionProbeset()] = ed1 ### Agglomerated probesets
            #if 'ENSMUSG00000032497:E23.1-E24.1' in ji.InclusionProbeset():
            #print ji.InclusionProbeset();sys.exit()
            #if '198878' in ji.InclusionProbeset(): print ji.InclusionProbeset(),excl
            try: alt_junction_db[geneid].append(ji)
            except Exception: alt_junction_db[geneid] = [ji]            

        del exon_inclusion_agglom
        critical_exon_db={}
        
    if array_type == 'RNASeq':
        ### Need to remove the @ from the IDs
        for e in exon_inclusion_db:
            incl_probesets=[]
            for i in exon_inclusion_db[e]:
                incl_probesets.append(string.replace(i,'@',':'))
            exon_inclusion_db[e] = incl_probesets

    ### Eliminate redundant entries
    objects_to_delete=[]
    for geneid in alt_junction_db:
        junction_temp_db={}; junction_temp_ls=[]
        for ji in alt_junction_db[geneid]: ### Redundant entries can be present
            id = ji.ExclusionProbeset(),ji.InclusionProbeset()
            if id in junction_temp_db: objects_to_delete.append(ji)
            else: junction_temp_db[id]=ji
        for i in junction_temp_db:
            ji = junction_temp_db[i]; junction_temp_ls.append(ji)
        alt_junction_db[geneid]=junction_temp_ls
                
    clearObjectsFromMemory(objects_to_delete); objects_to_delete=[]
    
    return alt_junction_db,critical_exon_db,exon_dbase,exon_inclusion_db,exon_db

def filterExistingFiles(species,array_type,db,export_type):
    """Remove probesets entries (including 5' and 3' junction exons) from the database that don't indicate possible critical exons"""
    export_exon_filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_'+array_type+'_probesets.txt'        
    ensembl_probeset_db = ExonArrayEnsemblRules.reimportEnsemblProbesetsForSeqExtraction(export_exon_filename,'probesets',{})
    
    critical_junction_db = {}; critical_probeset_db={}; crit1={}
    junction_inclusion_db = JunctionArrayEnsemblRules.reimportJunctionComps(species,array_type,'updated')
    for ids in junction_inclusion_db:
        for jd in junction_inclusion_db[ids]:
            critical_exon_id = jd.ParentCriticalExon()
            critical_id = jd.GeneID()+':'+jd.CriticalExon()
            critical_exon_ids = string.split(critical_exon_id,'|')
            critical_junction_db[jd.ExclusionProbeset(),jd.InclusionProbeset()]=critical_exon_ids,critical_id
            crit1[critical_id]=[]
    """
    for id in crit1:
        if 'ENSMUSG00000066842' in id: print id
    stop
    """
    
    #print len(crit1);
    crit2={}
        
    for (pX,probeset) in critical_junction_db:
        ###Keep only junction probesets that contain possible critical exons
        p1 = probeset+'|5'; p2 = probeset+'|3'
        c1s,critical_id = critical_junction_db[(pX,probeset)]; proceed = 'no'
        #print p1, p2, c1s, critical_id
        #for probeset in db: print [probeset];kill
        if probeset in ensembl_probeset_db and probeset in db:
            critical_probeset_db[probeset,critical_id]=db[probeset]
            crit2[probeset]=[]
        else:
            if p1 in ensembl_probeset_db and p1 in db:
                c2s = ensembl_probeset_db[p1]; p = p1
                c2s = string.split(c2s,'|')
                for c1 in c1s:
                    if c1 in c2s:
                        critical_probeset_db[p,critical_id]=db[p]
                        crit2[probeset]=[]                
            if p2 in ensembl_probeset_db and p2 in db:
                c2s = ensembl_probeset_db[p2]; p = p2
                c2s = string.split(c2s,'|')
                for c1 in c1s:
                    if c1 in c2s:
                        critical_probeset_db[p,critical_id]=db[p]
                        crit2[probeset]=[]

    for probeset in ensembl_probeset_db: ### For non-junction probesets
        if '|' not in probeset:
            if probeset in db: critical_probeset_db[probeset,probeset]=db[probeset]; crit2[probeset]=[]


    critical_probeset_db = eliminate_redundant_dict_values(critical_probeset_db)
    print len(crit2),'len(crit2)'
    x=0
    """
    for probeset in db:
        if probeset not in crit2:
            x+=1
            if x<20: print probeset """
    
    print len(critical_probeset_db),': length of filtered db', len(db), ': length of db'
    """
    for probeset in ensembl_probeset_db:
        ###Keep only probesets that contain possible critical exons
        if '|' in probeset:
            if probeset[:-2] in critical_junction_db and probeset in db:
                critical_probeset_db[probeset[:-2]]=db[probeset]
        elif probeset in db: critical_probeset_db[probeset]=db[probeset] """
    """
    for id in critical_probeset_db:
        if 'ENSMUSG00000066842' in id[1]: print id
    stop
    """
    
    if export_type == 'exclude_junction_psrs':
        critical_probeset_db['title'] = db['title']
        critical_probeset_db['filename'] = db['filename']
        exportFiltered(critical_probeset_db)
    else:
        for p in db:
            if '|' not in p: probeset = p
            else: probeset = p[:-2]
            if probeset not in crit2:
                ### Add back any junction probesets that do not have a critical exon component
                critical_probeset_db[probeset,probeset]=db[p]
                
        if export_type == 'exclude_critical_exon_ids':
            critical_probeset_db2={}
            for (p,cid) in critical_probeset_db:
                if ':' in cid or '|' in p:
                    critical_probeset_db2[p[:-2],p[:-2]] = critical_probeset_db[(p,cid)]
                else: critical_probeset_db2[p,p] = critical_probeset_db[(p,cid)]
            critical_probeset_db = critical_probeset_db2

        critical_probeset_db['title'] = db['title']
        critical_probeset_db['filename'] = db['filename']
        exportFiltered(critical_probeset_db)
    
########### Code originally designed for AltMouseA array database builds (adapted for use with Mouse and Human Junction Arrays)        
def filterExpressionData(filename1,filename2,pre_filtered_db,constitutive_db):
    fn2=filepath(filename2)
    probeset_translation_db={}
    ###Import probeset number/id relationships (note: forced to use numeric IDs for Plier/Exact analysis)
    if analysis_method != 'rma':
        for line in open(fn2,'r').xreadlines():             
            data = cleanUpLine(line)
            probeset_number,probeset_id = string.split(data,'\t')
            probeset_translation_db[probeset_number]=probeset_id

    fn=filepath(filename1)
    exp_dbase={}; d = 0; x = 0
    ###Import expression data (non-log space)
    try:
        for line in open(fn,'r').xreadlines():             
          data = cleanUpLine(line)
          if data[0] != '#' and x == 1:   ###Grab expression values
            tab_delimited_data = string.split(data,'\t')
            z = len(tab_delimited_data)
            probeset = tab_delimited_data[0]
            if analysis_method == 'rma': exp_vals = tab_delimited_data[1:]
            else: exp_vals = convertToLog2(tab_delimited_data[1:])
            ###Filter results based on whether a sufficient number of samples where detected as Present
            if probeset in pre_filtered_db:
                if probeset in probeset_translation_db: original_probeset_id = probeset_translation_db[probeset]
                else: original_probeset_id = probeset  ###When p-values are generated outside of Plier
                if original_probeset_id in constitutive_db:
                    percent_present = pre_filtered_db[probeset]
                    if percent_present > 0.99: exp_dbase[original_probeset_id] = exp_vals
                    #else: print percent_present,original_probeset_id; kill
                else: exp_dbase[original_probeset_id] = exp_vals
                
          elif data[0] != '#' and x == 0:  ###Grab labels
            array_names = []
            tab_delimited_data = string.split(data,'\t')
            for entry in tab_delimited_data: array_names.append(entry)
            x += 1
    except IOError: exp_dbase = exp_dbase
    print len(exp_dbase),"probesets imported with expression values"

            
    ###If the arrayid column header is missing, account for this
    if len(array_names) == z:
        array_names = array_names[1:]
    
    null,filename = string.split(filename1,'\\')
    filtered_exp_export = 'R_expression_raw_data\\'+filename[:-4]+'-filtered.txt'
    fn=filepath(filtered_exp_export); data = open(fn,'w'); title = 'probeset_id'
    for array in array_names: title = title +'\t'+ array
    data.write(title+'\n')
    for probeset in exp_dbase:
        exp_vals = probeset
        for exp_val in exp_dbase[probeset]:
            exp_vals = exp_vals +'\t'+ str(exp_val)
        data.write(exp_vals+'\n')
    data.close()
    #return filtered_exp_export

def convertToLog2(data_list):
    new_list=[]
    for item in data_list:
        new_list.append(math.log(float(item)+1,2))
    return new_list
        
def getAnnotations(filename,p,Species,Analysis_Method,constitutive_db):
    global species; species = Species
    global analysis_method; analysis_method = Analysis_Method
    array_type = 'AltMouse'
    
    filtered_junctions_list = ExonArray.getFilteredExons(filename,p)
    probe_id_translation_file = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'-probeset_translation.txt'
    
    filtered_exp_export_file = filterExpressionData(filename,probe_id_translation_file,filtered_junctions_list,constitutive_db)
    
    return filtered_exp_export_file

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def importGeneric(filename):
    verifyFile(filename)
    fn=filepath(filename); key_db = {}; x = 0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if len(t[1:]) == 1:
            try: key_db[t[0]].append(t[1])
            except KeyError: key_db[t[0]] = [t[1]]
        else: key_db[t[0]] = t[1:]
    return key_db

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def importAnnotateCriticalExonSequences(species,array_type):
    ensembl_associations = importArrayAnnotations(species,array_type)
    critical_exon_seq_db = importCriticalExonSeq('AltDatabase/'+species+'/'+array_type+'/'+ array_type+'_critical-exon-seq.txt',ensembl_associations)
    return critical_exon_seq_db

def importArrayAnnotations(species,array_type):
    primary_gene_annotation_file = 'AltDatabase/'+species +'/'+ array_type +'/'+ array_type+ '_gene_annotations.txt'
    ensembl_array_gene_annotation_file = 'AltDatabase/'+species+'/'+ array_type + '/'+array_type+ '-Ensembl.txt'
    ensembl_annotations = 'AltDatabase/ensembl/'+ species + '/'+species+ '_Ensembl-annotations_simple.txt'
    array_gene_annotations = importGeneric(primary_gene_annotation_file)
    ensembl_associations = importGeneric(ensembl_array_gene_annotation_file)
    ensembl_annotation_db = importGeneric(ensembl_annotations)
    
    ensembl_symbol_db={}
    for ens_geneid in ensembl_annotation_db:
        description, symbol = ensembl_annotation_db[ens_geneid]
        #print symbol;klll
        if len(symbol)>0:
            try: ensembl_symbol_db[symbol].append(ens_geneid)
            except KeyError: ensembl_symbol_db[symbol] =[ens_geneid]

    ### Update array Ensembl annotations        
    for array_geneid in array_gene_annotations:    
        t = array_gene_annotations[array_geneid]; description=t[0];entrez=t[1];symbol=t[2]
        if symbol in ensembl_symbol_db:
            ens_geneids = ensembl_symbol_db[symbol]
            for ens_geneid in ens_geneids:
                try: ensembl_associations[array_geneid].append(ens_geneid)
                except KeyError: ensembl_associations[array_geneid] = [ens_geneid]
    ensembl_associations = eliminate_redundant_dict_values(ensembl_associations)
    exportArrayIDEnsemblAssociations(ensembl_associations,species,array_type) ###Use these For LinkEST program
    return ensembl_associations

def exportDB(filename,db):
    fn=filepath(filename); data = open(fn,'w')
    for key in db:
        values = string.join([key]+db[key],'\t')+'\n'; data.write(values)
    data.close()

def exportFiltered(db):
    filename = db['filename']; title = db['title']
    filename = string.replace(filename,'.txt','-filtered.txt')
    print 'Writing',filename
    del db['filename']; del db['title']
    fn=filepath(filename); data = open(fn,'w'); data.write(title)
    for (old,new) in db:
        for line in db[(old,new)]: ### Replace the old ID with the new one
            if old not in line and '|' in old:
                old = old[:-2]
            line = string.replace(line,old,new); data.write(line)
    data.close()
    
def exportArrayIDEnsemblAssociations(ensembl_associations,species,array_type):   
    annotation_db_filename = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'-Ensembl_relationships.txt'
    fn=filepath(annotation_db_filename); data = open(fn,'w')
    title = ['ArrayGeneID','Ensembl']; title = string.join(title,'\t')+'\n'
    data.write(title)

    for array_geneid in ensembl_associations:
        for ens_geneid in ensembl_associations[array_geneid]:
            values = [array_geneid,ens_geneid]; values = string.join(values,'\t')+'\n'; data.write(values)
    data.close()

def exportCriticalExonLocations(species,array_type,critical_exon_seq_db):   
    location_db_filename = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'_critical_exon_locations.txt'
    fn=filepath(location_db_filename); data = open(fn,'w')
    title = ['Affygene','ExonID','Ensembl','start','stop','gene-start','gene-stop','ExonSeq']; title = string.join(title,'\t')+'\n'
    data.write(title)

    for ens_geneid in critical_exon_seq_db:
        for cd in critical_exon_seq_db[ens_geneid]:
            try:
                values = [cd.ArrayGeneID(),cd.ExonID(),ens_geneid,cd.ExonStart(),cd.ExonStop(),cd.GeneStart(), cd.GeneStop(), cd.ExonSeq()]
                values = string.join(values,'\t')+'\n'
                data.write(values)
            except AttributeError: null = []
    data.close()
    
class ExonSeqData:
    def __init__(self,exon,array_geneid,probeset_id,critical_junctions,critical_exon_seq):
        self._exon = exon; self._array_geneid = array_geneid; self._critical_junctions = critical_junctions
        self._critical_exon_seq = critical_exon_seq; self._probeset_id = probeset_id
    def ProbesetID(self): return self._probeset_id
    def ArrayGeneID(self): return self._array_geneid
    def ExonID(self): return self._exon
    def CriticalJunctions(self): return self._critical_junctions
    def ExonSeq(self): return string.upper(self._critical_exon_seq)
    def setExonStart(self,exon_start): self._exon_start = exon_start
    def setExonStop(self,exon_stop): self._exon_stop = exon_stop
    def setGeneStart(self,gene_start): self._gene_start = gene_start
    def setGeneStop(self,gene_stop): self._gene_stop = gene_stop
    def ExonStart(self): return str(self._exon_start)
    def ExonStop(self): return str(self._exon_stop)
    def GeneStart(self): return str(self._gene_start)
    def GeneStop(self): return str(self._gene_stop)
    
def importCriticalExonSeq(filename,ensembl_associations):
    fn=filepath(filename); key_db = {}; x = 0
    verifyFile(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x = 1
        else:
            arraygeneid_exon,critical_junctions,critical_exon_seq = string.split(data,'\t')
            if len(critical_exon_seq)>5:
                array_geneid, exon = string.split(arraygeneid_exon,':')
                if array_geneid in ensembl_associations:
                    ens_geneids = ensembl_associations[array_geneid]
                    for ens_geneid in ens_geneids:
                        seq_data = ExonSeqData(exon,array_geneid,arraygeneid_exon,critical_junctions,critical_exon_seq)
                        try: key_db[ens_geneid].append(seq_data)
                        except KeyError: key_db[ens_geneid] = [seq_data]
    return key_db

def updateCriticalExonSequences(array_type, filename,ensembl_probeset_db):
    exon_seq_db_filename = filename[:-4]+'_updated.txt'
    fn=filepath(exon_seq_db_filename); data = open(fn,'w')
    
    critical_exon_seq_db={}
    for ens_gene in ensembl_probeset_db:
        for probe_data in ensembl_probeset_db[ens_gene]:
            exon_id,((probe_start,probe_stop,probeset_id,exon_class,transcript_clust),ed) = probe_data
            try: critical_exon_seq_db[probeset_id] = ed.ExonSeq()
            except AttributeError: null=[] ### Occurs when no sequence data is associated with exon (probesets without exon associations)

    ensembl_probeset_db=[]; key_db = {}; x = 0
    if array_type == 'AltMouse':      
        fn1=filepath(filename)
        verifyFile(filename)
        for line in open(fn1,'rU').xreadlines():
            line_data = cleanUpLine(line)
            if x == 0: x = 1; data.write(line)
            else:
                arraygeneid_exon,critical_junctions,critical_exon_seq = string.split(line_data,'\t')
                if arraygeneid_exon in critical_exon_seq_db:
                    critical_exon_seq = critical_exon_seq_db[arraygeneid_exon]
                    values = [arraygeneid_exon,critical_junctions,critical_exon_seq]
                    values = string.join(values,'\t')+'\n'
                    data.write(values)
                else: data.write(line)
    elif array_type == 'junction':
        ### We don't need any of the additional information used for AltMouse arrays
        for probeset in critical_exon_seq_db:
            critical_exon_seq = critical_exon_seq_db[probeset]
            if ':' in probeset:
                probeset = string.split(probeset,':')[1]
            values = [probeset,'',critical_exon_seq]
            values = string.join(values,'\t')+'\n'
            data.write(values)  
    data.close()
    print exon_seq_db_filename, 'exported....'

def inferJunctionComps(species,array_type):
    if len(array_type) == 3:
        ### This indicates that the ensembl_probeset_db is already included
        array_type,ensembl_probeset_db,root_dir = array_type
        comps_type = ''
    else:
        export_exon_filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_probesets.txt'  
        ensembl_probeset_db = ExonArrayEnsemblRules.reimportEnsemblProbesetsForSeqExtraction(export_exon_filename,'junction-regions',{})
        comps_type = 'updated'; root_dir = ''

    if array_type != 'RNASeq':
        print "Import junction probeset region IDs for",species
        print "Preparing region IDs for analysis of possible reciprocal junctions"
    putative_as_junction_db={}; probeset_juntion_db={}; common_exon_blocks_exon={}; common_exon_blocks_intron={}; count=0
    for gene in ensembl_probeset_db:
        for (probeset,regionid) in ensembl_probeset_db[gene]:
            regionids = string.split(regionid,'|')
            for regionid in regionids:
                if '-' in regionid:
                    if 'I' in regionid: exons_type = 'exon-intron'
                    else: exons_type = 'exons'
                    exon_5prime_original, exon_3prime_original = string.split(regionid,'-')
                    
                    exon_5prime = string.split(exon_5prime_original,'.')
                    if '_' in exon_5prime[1]: exon_5prime[1] = float(string.replace(exon_5prime[1],'_','.'))
                    else: exon_5prime[1] = int(exon_5prime[1])
                    e1a3 = (int(exon_5prime[0][1:]),int(exon_5prime[1])) ### The first is an int for the region - since it hybs early
                    e1a5 = (int(exon_5prime[0][1:]),exon_5prime[1]) 
                    e1 = e1a3, e1a5

                    exon_3prime = string.split(exon_3prime_original,'.')
                    if '_' in exon_3prime[1]: exon_3prime[1] = float(string.replace(exon_3prime[1],'_','.'))
                    else:
                        try: exon_3prime[1] = int(exon_3prime[1])
                        except Exception: print exon_3prime;kill
                    e2a3 = (int(exon_3prime[0][1:]),exon_3prime[1])
                    e2a5 = (int(exon_3prime[0][1:]),int(exon_3prime[1])) ### The second is an int for the region - since it hybs late
                    e2 = e2a3, e2a5

                    if exons_type == 'exons':
                        count+=1
                        try: putative_as_junction_db[gene].append((e1,e2))
                        except Exception: putative_as_junction_db[gene] = [(e1,e2)]
                        
                        ### This matches the recorded junction ID from EnsemblImport.compareJunctions()
                        try: probeset_juntion_db[gene,(e1a5,e2a3)].append(probeset)
                        except Exception: probeset_juntion_db[gene,(e1a5,e2a3)] = [probeset]
                        
                        ### Defines exon-intron and exon-exon reciprical junctions based on shared exon blocks
                        block = e1a3[0]; side = 'left'
                        try: common_exon_blocks_exon[side,gene,block].append([regionid,probeset])
                        except KeyError: common_exon_blocks_exon[side,gene,block] = [[regionid,probeset]]
                        block = e2a3[0]; side = 'right'
                        try: common_exon_blocks_exon[side,gene,block].append([regionid,probeset])
                        except KeyError: common_exon_blocks_exon[side,gene,block] = [[regionid,probeset]]
                    else:
                        ### Defines exon-intron and exon-exon reciprical junctions based on shared exon blocks
                        if 'I' in exon_5prime: block = e2a3[0]; side = 'right'; critical_intron = exon_5prime_original
                        else: block = e1a3[0]; side = 'left'; critical_intron = exon_3prime_original
                        try: common_exon_blocks_intron[side,gene,block].append([regionid,probeset,critical_intron])
                        except KeyError: common_exon_blocks_intron[side,gene,block] = [[regionid,probeset,critical_intron]]
        
    if array_type != 'RNASeq':               
        print count, 'probed junctions being compared to identify putative reciprocal junction comparisons'
    critical_exon_db, critical_gene_junction_db = EnsemblImport.compareJunctions(species,putative_as_junction_db,{})
    if array_type != 'RNASeq':
        print len(critical_exon_db),'genes with alternative reciprocal junctions pairs found'
    
    global junction_inclusion_db; count=0; redundant=0; junction_annotations={}; critical_exon_annotations={}
    junction_inclusion_db = JunctionArrayEnsemblRules.reimportJunctionComps(species,array_type,(comps_type,ensembl_probeset_db))
    for gene in critical_exon_db:
        for sd in critical_exon_db[gene]:
            junction_pairs = getJunctionPairs(sd.Junctions())
            """
            if len(junction_pairs)>1 and len(sd.CriticalExonRegion())>1:
                print
                .Junctions()
                print sd.CriticalExonRegion();kill"""
            for (junction1,junction2) in junction_pairs:
                critical_exon = sd.CriticalExonRegion()
                excl_junction,incl_junction = determineExclIncl(junction1,junction2,critical_exon)
                incl_junction_probeset = probeset_juntion_db[gene,incl_junction][0]
                excl_junction_probeset = probeset_juntion_db[gene,excl_junction][0]
                source = 'Inferred'
                incl_junction=formatJunctions(incl_junction)
                excl_junction=formatJunctions(excl_junction)
                critical_exon=string.replace(formatJunctions(critical_exon),'-','|'); count+=1
                ji=JunctionArrayEnsemblRules.JunctionInformation(gene,critical_exon,excl_junction,incl_junction,excl_junction_probeset,incl_junction_probeset,source)
                if (excl_junction_probeset,incl_junction_probeset) not in junction_inclusion_db:
                    try: junction_inclusion_db[excl_junction_probeset,incl_junction_probeset].append(ji)
                    except KeyError: junction_inclusion_db[excl_junction_probeset,incl_junction_probeset] = [ji]
                    junction_str = string.join([excl_junction,incl_junction],'|')
                    #splice_event_str = string.join(sd.SpliceType(),'|')
                    try: junction_annotations[ji.InclusionProbeset()].append((junction_str,sd.SpliceType()))
                    except KeyError: junction_annotations[ji.InclusionProbeset()] = [(junction_str,sd.SpliceType())]
                    try: junction_annotations[ji.ExclusionProbeset()].append((junction_str,sd.SpliceType()))
                    except KeyError: junction_annotations[ji.ExclusionProbeset()] = [(junction_str,sd.SpliceType())]
                    critical_exons = string.split(critical_exon,'|')
                    for critical_exon in critical_exons:
                        try: critical_exon_annotations[gene+':'+critical_exon].append((junction_str,sd.SpliceType()))
                        except KeyError: critical_exon_annotations[gene+':'+critical_exon] = [(junction_str,sd.SpliceType())]
                else: redundant+=1
                
    if array_type != 'RNASeq':
        print count, 'Inferred junctions identified with',redundant, 'redundant.'
    
    ### Compare exon and intron blocks for intron alinging junctions
    junction_inclusion_db = annotateNovelIntronSplicingEvents(common_exon_blocks_intron,common_exon_blocks_exon,junction_inclusion_db)

    if len(root_dir)>0: exportUpdatedJunctionComps((species,root_dir),array_type)
    else: exportUpdatedJunctionComps(species,array_type)
    
    clearObjectsFromMemory(junction_inclusion_db); junction_inclusion_db=[]
    if array_type == 'RNASeq':    
        ### return these annotations for RNASeq analyses
        return junction_annotations,critical_exon_annotations

def annotateNovelIntronSplicingEvents(common_exon_blocks_intron,common_exon_blocks_exon,junction_inclusion_db):
    ### Add exon-intron, exon-exon reciprical junctions determined based on common block exon (same side of the junction)
    new_intron_events=0
    for key in common_exon_blocks_intron:
        (side,gene,block) = key; source='Inferred-Intron'
        if key in common_exon_blocks_exon:
            for (excl_junction,excl_junction_probeset) in common_exon_blocks_exon[key]:
                for (incl_junction,incl_junction_probeset,critical_intron) in common_exon_blocks_intron[key]:
                    ji=JunctionArrayEnsemblRules.JunctionInformation(gene,critical_intron,excl_junction,incl_junction,excl_junction_probeset,incl_junction_probeset,source)
                    if (excl_junction_probeset,incl_junction_probeset) not in junction_inclusion_db:
                        try: junction_inclusion_db[excl_junction_probeset,incl_junction_probeset].append(ji)
                        except Exception: junction_inclusion_db[excl_junction_probeset,incl_junction_probeset] = [ji]
                        new_intron_events+=1
    #print new_intron_events, 'novel intron-splicing events added to database'
    """
    ### While the below code seemed like a good idea, the current state of RNA-seq alignment tools produced a rediculous amount of intron-intron junctions (usually in the same intron)
    ### Without supporting data (e.g., other junctions bridging these intron junction to a validated exon), we must assume these juncitons are not associated with the alinging gene
    new_intron_events=0 ### Compare Intron blocks to each other
    for key in common_exon_blocks_intron:
        (side,gene,block) = key; source='Inferred-Intron'
        for (excl_junction,excl_junction_probeset,critical_intron1) in common_exon_blocks_intron[key]:
            for (incl_junction,incl_junction_probeset,critical_intron2) in common_exon_blocks_intron[key]:
                    if (excl_junction,excl_junction_probeset) != (incl_junction,incl_junction_probeset): ### If comparing entries in the same list, don't compare an single entry to itself
                        ji=JunctionArrayEnsemblRules.JunctionInformation(gene,critical_intron1+'|'+critical_intron2,excl_junction,incl_junction,excl_junction_probeset,incl_junction_probeset,source)
                        if (excl_junction_probeset,incl_junction_probeset) not in junction_inclusion_db and (incl_junction_probeset,excl_junction_probeset) not in junction_inclusion_db:
                            try: junction_inclusion_db[excl_junction_probeset,incl_junction_probeset].append(ji)
                            except Exception: junction_inclusion_db[excl_junction_probeset,incl_junction_probeset] = [ji]
                            new_intron_events+=1
    """
    #print new_intron_events, 'novel intron-splicing events added to database'
    return junction_inclusion_db

def determineExclIncl(junction1,junction2,critical_exons):
    #((3, 2), (6, 1))
    for critical_exon in critical_exons:
        if critical_exon in junction1: incl_junction = junction1; excl_junction = junction2
        if critical_exon in junction2: incl_junction = junction2; excl_junction = junction1
    return excl_junction,incl_junction

def formatJunctions(junction):
    #((3, 2), (6, 1))
    exons_to_join=[]
    for i in junction:
        exons_to_join.append('E'+str(i[0])+'.'+string.replace(str(i[1]),'.','_'))
    junction_str = string.join(exons_to_join,'-')
    return junction_str
    
def getJunctionPairs(junctions):
    ### Although the pairs of junctions (exclusion 1st, inclusion 2nd) are given, need to separate out the pairs to report reciprical junctions
    # (((3, 2), (6, 1)), ((4, 2), (6, 1)), ((3, 2), (6, 1)), ((4, 2), (6, 1))))
    count = 0; pairs=[]; pairs_db={}
    for junction in junctions:
        count +=1; pairs.append(junction)
        if count==2: pairs_db[tuple(pairs)]=[]; count = 0; pairs=[]
    return pairs_db

def getJunctionExonLocations(species,array_type):
    global ensembl_associations
    ensembl_associations = importJunctionArrayAnnotations(species,array_type)
    extraction_type = 'sequence'
    exon_seq_db=importJunctionArrayAnnotationMappings(array_type,species,extraction_type)
    getLocations(species,array_type,exon_seq_db)
    
def identifyCriticalExonLocations(species,array_type):
    critical_exon_seq_db = importAnnotateCriticalExonSequences(species,array_type)
    getLocations(species,array_type,critical_exon_seq_db)
    
def getLocations(species,array_type,critical_exon_seq_db):
    analysis_type = 'get_locations'
    dir = 'AltDatabase/'+species+'/SequenceData/chr/'+species; gene_seq_filename = dir+'_gene-seq-2000_flank.fa'
    critical_exon_seq_db = EnsemblImport.import_sequence_data(gene_seq_filename,critical_exon_seq_db,species,analysis_type)
    exportCriticalExonLocations(species,array_type,critical_exon_seq_db)

def reAnnotateCriticalExonSequences(species,array_type):
    export_exon_filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_'+array_type+'_probesets.txt'        
    ensembl_probeset_db = ExonArrayEnsemblRules.reimportEnsemblProbesetsForSeqExtraction(export_exon_filename,'null',{})
    
    analysis_type = 'get_sequence'
    dir = 'AltDatabase/'+species+'/SequenceData/chr/'+species; gene_seq_filename = dir+'_gene-seq-2000_flank.fa'
    ensembl_probeset_db = EnsemblImport.import_sequence_data(gene_seq_filename,ensembl_probeset_db,species,analysis_type)

    critical_exon_file = 'AltDatabase/'+species+'/'+ array_type + '/' + array_type+'_critical-exon-seq.txt'
    updateCriticalExonSequences(array_type, critical_exon_file, ensembl_probeset_db)
    
if __name__ == '__main__':
    """Module has methods for annotating Junction associated critical exon sequences with up-to-date genome coordinates and analysis options for
    junciton arrays from AnalyzeExpressionDatasets"""
    m = 'Mm'; h = 'Hs'; species = h; array_type = 'RNASeq' ###In theory, could be another type of junction or combination array
    combineExonJunctionAnnotations(species,array_type);sys.exit()
    filterForCriticalExons(species,array_type)
    overRideExonEntriesWithJunctions(species,array_type);sys.exit()
    
    #inferJunctionComps(species,array_type); sys.exit()
    
    identifyJunctionComps(species,array_type);sys.exit()
    filterForCriticalExons(species,array_type);sys.exit()
    reAnnotateCriticalExonSequences(species,array_type)

    #getJunctionExonLocations(species,array_type)
    sys.exit()
    
    import_dir = '/AltDatabase/exon/'+species; expr_file_dir = 'R_expression_raw_data\exp.altmouse_es-eb.dabg.rma.txt'
    dagb_p = 1; Analysis_Method = 'rma'

    #identifyCriticalExonLocations(species,array_type)
    #JunctionArrayEnsemblRules.getAnnotations(species,array_type)
    
    ### Only needs to be run once, to update the original 
    #reAnnotateCriticalExonSequences(species,array_type); sys.exit() 
    
    
    #getAnnotations(expr_file_dir,dagb_p,Species,Analysis_Method)
    
