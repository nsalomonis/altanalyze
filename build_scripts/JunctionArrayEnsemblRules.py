###JunctionArrayEnsemblRules
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
try: from build_scripts import ExonArrayEnsemblRules
except Exception: pass
try: from build_scripts import EnsemblImport
except Exception: pass
import shutil
try: from build_scripts import JunctionArray
except Exception: pass
import update

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list
        
def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def eliminateRedundant(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

############# Affymetrix NextGen Junction Array Code

def importEnsemblUCSCAltJunctions(species,type):
    if type == 'standard':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_alternative_junctions.txt'
    else: 
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_alternative_junctions'+type+'.txt'
    fn=filepath(filename); x = 0; gene_junction_db={}
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        if x == 0: x=1
        else:
            gene,critical_exon,junction1,junction2,splice_event = string.split(data,'\t')
            if critical_exon in junction1: incl_junction = junction1; excl_junction = junction2 
            else: incl_junction = junction2; excl_junction = junction1
            gene_junction_db[gene,critical_exon,incl_junction,excl_junction]=[]
            #try: gene_junction_db[gene,incl_junction,excl_junction].append(critical_exon)
            #except KeyError: gene_junction_db[gene,incl_junction,excl_junction] = [critical_exon]

    print len(gene_junction_db), 'alternative junction-pairs identified from Ensembl and UCSC'
    return gene_junction_db

def getJunctionComparisonsFromExport(species,array_type):
    type = 'standard'
    gene_junction_db = importEnsemblUCSCAltJunctions(species,type)
    
    ### Retrieve probesets with exon-junctions associated - these are critical exons
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_'+array_type+'_probesets.txt'
    gene_probeset_db = ExonArrayEnsemblRules.reimportEnsemblProbesetsForSeqExtraction(filename,'junctions',{})
    left={}; right={}; gene_db={}; gene_exon_db={}; nonjunction_aligning={}
    for gene in gene_probeset_db:
        for (probe_data,ed) in gene_probeset_db[gene]:
            probeset, strand, probeset_start, probeset_stop = probe_data
            region_id = string.replace(ed.RegionNumber(),'-','.')
            original_region_id = region_id
            region_ids = string.split(region_id,'|')
            gene_db[probeset[:-2]]=gene
            #ed.AssociatedSplicingJunctions()
            r_starts=string.split(ed.ExonStart(),'|'); r_stops=string.split(ed.ExonStop(),'|')
            for region_id in region_ids:
                if '|5' in probeset:
                    try: left[probeset[:-2]].append(region_id)
                    except Exception: left[probeset[:-2]]=[region_id]
                    if strand == '+': ### If the junction probesets DO NOT align to the region coordinates, then the probeset maps to a junction outside the database
                        if probeset_stop not in r_stops: nonjunction_aligning[probeset[:-2]] = original_region_id+'_'+probeset_stop,'left'
                    elif probeset_start not in r_starts: nonjunction_aligning[probeset[:-2]] = original_region_id+'_'+probeset_start,'left'
                elif '|3' in probeset:
                    try: right[probeset[:-2]].append(region_id)
                    except Exception: right[probeset[:-2]]=[region_id]
                    if strand == '+':
                        if probeset_start not in r_starts: nonjunction_aligning[probeset[:-2]] = original_region_id+'_'+probeset_start,'right'
                    elif probeset_stop not in r_stops: nonjunction_aligning[probeset[:-2]] = original_region_id+'_'+probeset_stop,'right'
                else:
                    if '_' in region_id: print killer
                    try: gene_exon_db[gene,region_id].append(probeset)
                    except Exception: gene_exon_db[gene,region_id] = [probeset]

    print 'len(nonjunction_aligning)',len(nonjunction_aligning)
    gene_exon_db = eliminateRedundant(gene_exon_db)            
    junction_db={} ### Get the exon-region IDs for an exon-junction
    for probeset in left:
        gene = gene_db[probeset]
        if probeset in right:
            for region1 in left[probeset]:
                for region2 in right[probeset]:
                    junction = region1+'-'+region2
                    try: junction_db[gene,junction].append(probeset)
                    except Exception: junction_db[gene,junction] = [probeset]

    probeset_junction_export = 'AltDatabase/' + species + '/'+array_type+'/'+ species + '_junction_comps.txt'
    
    fn=filepath(probeset_junction_export); data = open(fn,'w')
    print "Exporting",probeset_junction_export
    title = 'gene'+'\t'+'critical_exon'+'\t'+'exclusion_junction_region'+'\t'+'inclusion_junction_region'+'\t'+'exclusion_probeset'+'\t'+'inclusion_probeset'+'\t'+'data_source'+'\n'
    data.write(title); temp_list=[]
    
    for (gene,critical_exon,incl_junction,excl_junction) in gene_junction_db:
        if (gene,incl_junction) in junction_db:
            incl_junction_probesets = junction_db[gene,incl_junction]
            if (gene,excl_junction) in junction_db:
                excl_junction_probesets = junction_db[gene,excl_junction]
                for incl_junction_probeset in incl_junction_probesets:
                    for excl_junction_probeset in excl_junction_probesets:
                        try:
                            for incl_exon_probeset in gene_exon_db[gene,critical_exon]:
                                if incl_junction_probeset in nonjunction_aligning or excl_junction_probeset in nonjunction_aligning: null=[]
                                else: ### Ensure the probeset DOES map to the annotated junctions
                                    temp_list.append(string.join([gene,critical_exon,excl_junction,critical_exon,excl_junction_probeset,incl_exon_probeset,'AltAnalyze'],'\t')+'\n')
                        except Exception: null=[]
                        if incl_junction_probeset in nonjunction_aligning:
                            new_region_id, side = nonjunction_aligning[incl_junction_probeset]
                            incl_junction = renameJunction(incl_junction,side,new_region_id)
                        if excl_junction_probeset in nonjunction_aligning:
                            new_region_id, side = nonjunction_aligning[excl_junction_probeset]
                            excl_junction = renameJunction(excl_junction,side,new_region_id)
                        if excl_junction_probeset!=incl_junction_probeset:
                            temp_list.append(string.join([gene,critical_exon,excl_junction,incl_junction,excl_junction_probeset,incl_junction_probeset,'AltAnalyze'],'\t')+'\n')
    temp_list = unique.unique(temp_list)
    for i in temp_list: data.write(i)
    data.close()
    print 'Number of compared junctions exported', len(temp_list)
    

def renameJunction(junction_id,side,new_region_id):
    try: l,r = string.split(junction_id,'-')
    except Exception: print junction_id;kill
    if side == 'left': junction_id = new_region_id+'-'+r
    else: junction_id = l+'-'+new_region_id
    return junction_id
    
class JunctionInformation:
    def __init__(self,gene,critical_exon,excl_junction,incl_junction,excl_probeset,incl_probeset,source):
        self._gene = gene; self._critical_exon = critical_exon; self.excl_junction = excl_junction; self.incl_junction = incl_junction
        self.excl_probeset = excl_probeset; self.incl_probeset = incl_probeset; self.source = source
        self.critical_exon_sets = [critical_exon]
    def GeneID(self): return str(self._gene)
    def CriticalExon(self):
        ce = str(self._critical_exon)
        if '-' in ce: ce = string.replace(ce,'-','.')
        return ce
    def CriticalExonList(self):
        critical_exon_str = self.CriticalExon()
        critical_exons = string.split(critical_exon_str,'|')
        return critical_exons
    def ParentCriticalExon(self):
        ce = str(self.CriticalExon())
        if '_' in ce:
            ces = string.split(ce,'|')
            ces2=[]
            for ce in ces:
                if '_' in ce:
                    ce = string.split(ce,'_')[0]
                    ces2.append(ce)
            ce = string.join(ces2,'|')
        return ce
    def setCriticalExons(self,critical_exons): self._critical_exon = critical_exons
    def setCriticalExonSets(self,critical_exon_sets): self.critical_exon_sets = critical_exon_sets
    def CriticalExonSets(self): return self.critical_exon_sets ### list of critical exons (can select any or all for functional analysis)
    def setInclusionJunction(self,incl_junction): self.incl_junction = incl_junction
    def InclusionJunction(self): return self.FormatJunction(self.incl_junction,self.InclusionProbeset())
    def ExclusionJunction(self): return self.FormatJunction(self.excl_junction,self.ExclusionProbeset())
    def FormatJunction(self,junction,probeset):
        ### Used for RNASeq when trans-splicing occurs
        probeset_split = string.split(probeset,':') ### Indicates trans-splicing - two gene IDs listed
        if len(probeset_split)>2:
            junction = probeset
        elif 'E0.1' in junction:
            junction = 'U'+junction[1:]
        return str(junction)
    def setInclusionProbeset(self,incl_probeset): self.incl_probeset = incl_probeset
    def InclusionProbeset(self): return self.FormatProbeset(self.incl_probeset)
    def ExclusionProbeset(self): return self.FormatProbeset(self.excl_probeset)
    def FormatProbeset(self,probeset):
        if ':' in probeset:
            probeset = string.split(probeset,':')[1]
        elif '@' in probeset:
            probeset = reformatID(probeset)
        return str(probeset)
    def DataSource(self): return str(self.source)
    def OutputLine(self):
        value = string.join([self.GeneID(),self.CriticalExon(),self.ExclusionJunction(),self.InclusionJunction(),self.ExclusionProbeset(),self.InclusionProbeset(),self.DataSource()],'\t')+'\n'
        return value
    def __repr__(self): return self.GeneID()
    
def formatID(id):
    ### JunctionArray methods handle IDs with ":" different than those that lack this
    return string.replace(id,':','@')

def reformatID(id):
    ### JunctionArray methods handle IDs with ":" different than those that lack this
    return string.replace(id,'@',':')

def reimportJunctionComps(species,array_type,file_type):
    if len(species[0])>1:
        species, root_dir = species
    else: species = species
    if len(file_type) == 2:
        file_type, filter_db = file_type; filter='yes'
    else: filter_db={}; filter='no'
        
    
    filename = 'AltDatabase/' + species + '/'+array_type+'/'+ species + '_junction_comps.txt'
    if file_type == 'updated':
        filename = string.replace(filename,'.txt','_updated.txt')
        if array_type == 'RNASeq': filename = root_dir+filename
    fn=filepath(filename); junction_inclusion_db={}; x=0
    for line in open(fn,'rU').xreadlines():
        if x==0: x=1
        else:
            data = cleanUpLine(line)
            try: gene,critical_exon,excl_junction,incl_junction,excl_junction_probeset,incl_junction_probeset,source = string.split(data,'\t')
            except Exception: print data;kill
            proceed = 'yes'
            if filter == 'yes':
                if gene not in filter_db: proceed = 'no'
            if proceed == 'yes':
                if array_type == 'RNASeq':
                    ### Reformat IDs, as junction arrays interpret ':' as a separator to remove the gene ID
                    excl_junction_probeset = formatID(excl_junction_probeset)
                    incl_junction_probeset = formatID(incl_junction_probeset)
                    if file_type == 'updated':
                        ### Add exclusion-junction versus inclusion-exon
                        critical_exons = string.split(critical_exon,'|')
                        for ce in critical_exons:
                            critical_exon_probeset = formatID(gene+':'+ce)
                            ji=JunctionInformation(gene,critical_exon,excl_junction,ce,excl_junction_probeset,critical_exon_probeset,source)
                            junction_inclusion_db[excl_junction_probeset,critical_exon_probeset] = [ji]
                ji=JunctionInformation(gene,critical_exon,excl_junction,incl_junction,excl_junction_probeset,incl_junction_probeset,source)
                if ':' in excl_junction_probeset:
                    tc,excl_junction_probeset=string.split(excl_junction_probeset,':')
                    tc,incl_junction_probeset=string.split(incl_junction_probeset,':')
                try: junction_inclusion_db[excl_junction_probeset,incl_junction_probeset].append(ji)
                except Exception: junction_inclusion_db[excl_junction_probeset,incl_junction_probeset] = [ji]
    if array_type != 'RNASeq':
        print len(junction_inclusion_db),'reciprocol junctions imported'
    return junction_inclusion_db

def importAndReformatEnsemblJunctionAnnotations(species,array_type,nonconstitutive_junctions):
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_'+array_type+'_probesets.txt'
    export_filepath = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_probesets.txt'
    efn=filepath(export_filepath); export_data = open(efn,'w')
    
    fn=filepath(filename); x = 0; ensembl_exon_db={}; left={}; right={}; exon_gene_db={}; nonjunction_aligning={}
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        if x == 0: x=1; export_data.write(data+'\n')
        else:
            t = string.split(data,'\t')
            probeset, exon_id, ensembl_gene_id, transcript_cluster_id, chr, strand, probeset_start, probeset_stop, affy_class, constitutitive_probeset, ens_exon_ids, exon_annotations,regionid,r_start,r_stop,splice_event,splice_junctions = t
            if len(regionid)<1: regionid = exon_id; t[12] = exon_id     
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            tc,probeset=string.split(probeset,':'); regionid = string.replace(regionid,'-','.'); original_region_id = regionid
            r_starts=string.split(r_start,'|'); r_stops=string.split(r_stop,'|')
            ed = EnsemblImport.ExonStructureData(ensembl_gene_id, chr, strand, probeset_start, probeset_stop, constitutitive_probeset, ens_exon_ids, []); ed.reSetExonID(regionid)
            if '|5' in probeset:
                left[probeset[:-2]] = ed,t
                if strand == '+': ### If the junction probesets DO NOT align to the region coordinates, then the probeset maps to a junction outside the database
                    if probeset_stop not in r_stops: nonjunction_aligning[probeset[:-2]] = original_region_id+'_'+probeset_stop,'left'
                elif probeset_start not in r_starts: nonjunction_aligning[probeset[:-2]] = original_region_id+'_'+probeset_start,'left'
            elif '|3' in probeset:
                right[probeset[:-2]] = ed,t
                if strand == '+':
                    if probeset_start not in r_starts: nonjunction_aligning[probeset[:-2]] = original_region_id+'_'+probeset_start,'right'
                elif probeset_stop not in r_stops: nonjunction_aligning[probeset[:-2]] = original_region_id+'_'+probeset_stop,'right'
            else:
                t[0] = probeset
                ensembl_exon_db[probeset] = ed
                export_data.write(string.join(t,'\t')+'\n')
                regionids = string.split(regionid,'|')
                for regionid in regionids: exon_gene_db[ensembl_gene_id,regionid] = probeset
                
    for probeset in left:
        if probeset in right:
            l,pl = left[probeset]; r,pr = right[probeset]
            if l.Constitutive() != r.Constitutive(): l.setConstitutive('no') ### used to determine if a junciton is alternative or constitutive
            if probeset in nonconstitutive_junctions: l.setConstitutive('no')
            l.setJunctionCoordinates(l.ExonStart(),l.ExonStop(),r.ExonStart(),r.ExonStop())            
            ens_exon_idsl = pl[10]; ens_exon_idsr = pr[10]; exon_idl = pl[1]; exon_idr = pr[1]
            regionidl = pl[12]; regionidr = pr[12]; splice_junctionsl = pl[-1]; splice_junctionsr = pr[-1]
            exon_idl = string.replace(exon_idl,'-','.'); exon_idr = string.replace(exon_idr,'-','.')
            regionidl_block = string.split(regionidl,'-')[0]; regionidr_block = string.split(regionidr,'-')[0]
                            
            if regionidl_block != regionidr_block: ### Otherwise, the junction is probing a single exon block and thus is not informative     
                regionidl = string.replace(regionidl,'-','.'); regionidr = string.replace(regionidr,'-','.')
                exon_id = exon_idl+'-'+exon_idr; regionid = regionidl+'-'+regionidr
                
                if probeset in nonjunction_aligning:
                    new_region_id, side = nonjunction_aligning[probeset]
                    regionid = renameJunction(regionid,side,new_region_id)
                
                l.reSetExonID(regionid); ensembl_exon_db[probeset] = l
                
                splice_junctionsl+=splice_junctionsr
                ens_exon_idsl = string.split(ens_exon_idsl,'|'); ens_exon_idsr = string.split(ens_exon_idsr,'|')
                ens_exon_ids=string.join(unique.unique(ens_exon_idsl+ens_exon_idsr),'|')
                pl[10] = ens_exon_ids; pl[12] = regionid; pl[1] = exon_id; pl[-1] = splice_junctionsl
                pl[13] = l.ExonStart()+'|'+l.ExonStop(); pl[14] = r.ExonStart()+'|'+r.ExonStop()
                strand = pl[5]
                if strand == '+':
                    pl[6] = l.ExonStop(); pl[7] = r.ExonStart() ### juncstion splice-sites
                else:
                    pl[6] = l.ExonStart(); pl[7] = r.ExonStop() ### juncstion splice-sites

                pl[0] = probeset; pl[9] = l.Constitutive()
            
                pl = string.join(pl,'\t')+'\n'
                export_data.write(pl)

    export_data.close()
    return ensembl_exon_db,exon_gene_db

################## AltMouse and Generic Junction Array Analysis

def importCriticalExonLocations(species,array_type,ensembl_exon_db,force):
    ###ensembl_exon_db[(geneid,chr,strand)] = [[E5,exon_info]] #exon_info = (exon_start,exon_stop,exon_id,exon_annot)
    ###ensembl_probeset_db[geneid,chr,strand].append(probeset_data) #probeset_data = [start,stop,probeset_id,exon_class,transcript_cluster_id]
    gene_info_db = {}
    for (ens_geneid,chr,strand) in ensembl_exon_db: gene_info_db[ens_geneid] = chr,strand
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'_critical_exon_locations.txt'
    array_ensembl={}

    ###Get the most recent gene-symbol annotations (applicable with a new Ensembl build for the same genomic build)
    ensembl_symbol_db = getEnsemblAnnotations(species)
    primary_gene_annotation_file = 'AltDatabase/'+species +'/'+ array_type +'/'+ array_type+ '_gene_annotations.txt'
    update.verifyFile(primary_gene_annotation_file,array_type)
    array_gene_annotations = JunctionArray.importGeneric(primary_gene_annotation_file)
            
    for array_geneid in array_gene_annotations:    
        t = array_gene_annotations[array_geneid]; description=t[0];entrez=t[1];symbol=t[2]
        if symbol in ensembl_symbol_db and len(symbol)>0 and len(array_geneid)>0:
            ens_geneid = ensembl_symbol_db[symbol]
            if len(ens_geneid)>0: array_ensembl[array_geneid]= ens_geneid
          
    update.verifyFile(filename,array_type)  
    ensembl_probeset_db = importJunctionLocationData(filename,array_ensembl,gene_info_db,test)
       
    print len(ensembl_probeset_db), "Genes inlcuded in",array_type,"location database"
    return ensembl_probeset_db

def importJunctionLocationData(filename,array_ensembl,gene_info_db,test):
    fn=filepath(filename); key_db = {}; x = 0; y = 0; ensembl_probeset_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x = 1
        else:
            k=0
            array_geneid,exonid,ens_geneid,start,stop,gene_start,gene_stop,exon_seq = string.split(data,'\t')
            probeset_id = array_geneid+':'+exonid
            if test == 'yes':
                if array_geneid in test_cluster: k=1
            else: k = 1
            if k==1:
                if array_geneid in array_ensembl: ens_geneid = array_ensembl[array_geneid]; y+=1 ### Over-ride any outdated assocations
                if ens_geneid in gene_info_db:
                    chr,strand = gene_info_db[ens_geneid]
                    probeset_data = [start,stop,probeset_id,'core',array_geneid]
                    try: ensembl_probeset_db[ens_geneid,'chr'+chr,strand].append(probeset_data)
                    except KeyError: ensembl_probeset_db[ens_geneid,'chr'+chr,strand] = [probeset_data]
    print y,"ArrayID to Ensembl genes re-annotated..."
    return ensembl_probeset_db

def getEnsemblAnnotations(species):
    filename = 'AltDatabase/ensembl/'+ species + '/'+species+ '_Ensembl-annotations_simple.txt'
    ensembl_annotation_db = {}
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        ensembl_gene_id,description,symbol = string.split(data,'\t')
        ensembl_annotation_db[symbol] = ensembl_gene_id
    return ensembl_annotation_db

def getAnnotations(Species,array_type,reannotate_exon_seq,force):
    """Annotate Affymetrix exon array data using files Ensembl data (sync'ed to genome release)."""
    global species; species = Species; global test; global test_cluster
    test = 'no'; test_cluster = ['TC0701360']; data_type = 'mRNA'

    global ensembl_exon_db; global ensembl_exon_db; global exon_clusters; global exon_region_db
    ensembl_exon_db,ensembl_annot_db,exon_clusters,intron_clusters,exon_region_db,intron_retention_db,ucsc_splicing_annot_db,ens_transcript_db = EnsemblImport.getEnsemblAssociations(species,data_type,test)
    ensembl_probeset_db = importCriticalExonLocations(species,array_type,ensembl_exon_db,force) ###Get Pre-computed genomic locations for critical exons
    ensembl_probeset_db = ExonArrayEnsemblRules.annotateExons(ensembl_probeset_db,exon_clusters,ensembl_exon_db,exon_region_db,intron_retention_db,intron_clusters,ucsc_splicing_annot_db); constitutive_gene_db={}
    ExonArrayEnsemblRules.exportEnsemblLinkedProbesets(array_type,ensembl_probeset_db,species)
    print "\nCritical exon data exported coordinates, exon associations and splicing annotations exported..."
    
    ### Change filenames to reflect junction array type
    export_filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_probesets.txt'; ef=filepath(export_filename)
    export_replacement = string.replace(export_filename,'_probe','_'+array_type+'_probe')
    er=filepath(export_replacement); shutil.copyfile(ef,er); os.remove(ef) ### Copy file to a new name

    ### Export full exon seqeunce for probesets/critical exons to replace the original incomplete sequence (used for miRNA analyses)
    if reannotate_exon_seq == 'yes':
        JunctionArray.reAnnotateCriticalExonSequences(species,array_type)
    
def annotateJunctionIDsAsExon(species,array_type):
    from build_scripts import ExonSeqModule
    probeset_annotations_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_junction_probesets-filtered.txt'
    if array_type == 'RNASeq':
        probeset_annotations_file = string.replace(probeset_annotations_file,'junction_probesets-filtered','exons')
    junction_exon_db = ExonSeqModule.importSplicingAnnotationDatabase(probeset_annotations_file,array_type)
    probeset_annotations_file = 'AltDatabase/'+species+'/exon/'+species+'_Ensembl_probesets.txt'
    exon_db = ExonSeqModule.importSplicingAnnotationDatabase(probeset_annotations_file,array_type)
    
    ### Extract unique exon regions from Exon Array annotations
    multiple_exon_regions={}; unique_exon_regions={}
    for probeset in exon_db:
        y = exon_db[probeset]
        geneid = y.GeneID()
        if '|' in y.ExonRegionID():
            exonids = string.split(y.ExonRegionID(),'|')
            for exonid in exonids: multiple_exon_regions[geneid,exonid] = y
        else:
            unique_exon_regions[geneid,y.ExonRegionID()] = y
    ### Add missing exons to unique
    for uid in multiple_exon_regions:
        if uid not in unique_exon_regions: unique_exon_regions[uid]=multiple_exon_regions[uid]

    """
        for i in unique_exon_regions:
            if 'ENSMUSG00000066842' in i:
                print i
    stop
    """
    
    ### Extract unique exon regions from Junction Array annotation
    junction_to_exonids={}
    for probeset in junction_exon_db:
        if 'ENSMUSG00000066842' in probeset: print probeset
        y = junction_exon_db[probeset]
        geneid = y.GeneID()
        if '|' in y.ExonRegionID():
            exonids = string.split(y.ExonRegionID(),'|')
            if probeset == 'ENSMUSG00000066842|E60.1': print [[exonids]]
            for exonid in exonids:
                if (geneid,exonid) in unique_exon_regions:
                    y = unique_exon_regions[geneid,exonid]
                    if probeset == 'ENSMUSG00000066842:E60.1': print [y.Probeset()]
                    junction_to_exonids[probeset] = y.Probeset()
        else:
            if (geneid,string.replace(y.ExonRegionID(),'.','-')) in unique_exon_regions:
                #if ':' in probeset: print [probeset,y.ExonRegionID()];kill
                y = unique_exon_regions[geneid,string.replace(y.ExonRegionID(),'.','-')]
                junction_to_exonids[probeset] = y.Probeset()
                
    output_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_'+array_type+'-exon_probesets.txt'
    fn=filepath(output_file); data = open(fn,'w')
    data.write(array_type+'_probeset\texon_probeset\n')
    
    for probeset in junction_to_exonids:
        exon_probeset = junction_to_exonids[probeset]
        data.write(probeset+'\t'+exon_probeset+'\n')    
    data.close()
    
if __name__ == '__main__':
    m = 'Mm'; h = 'Hs'
    Species = m
    array_type = 'RNASeq' ###In theory, could be another type of junciton or combination array

    annotateJunctionIDsAsExon(Species,array_type); sys.exit()
    
    #reimportJunctionComps(Species,array_type,'original');kill
    #JunctionArray.getJunctionExonLocations(Species,array_type)
    """
    ### Get UCSC associations (download databases if necessary)
    from build_scripts import UCSCImport
    mRNA_Type = 'mrna'; run_from_scratch = 'yes'; force='no'
    export_all_associations = 'no' ### YES only for protein prediction analysis
    UCSCImport.runUCSCEnsemblAssociations(Species,mRNA_Type,export_all_associations,run_from_scratch,force)
    """
    getAnnotations(Species,array_type,'yes','no')
    JunctionArray.identifyJunctionComps(Species,array_type)
    #importAndReformatEnsemblJunctionAnnotations(Species,array_type)
    
    
