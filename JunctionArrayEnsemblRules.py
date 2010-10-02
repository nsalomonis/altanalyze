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

import sys, string
import os.path
import unique
import ExonArrayEnsemblRules
import EnsemblImport
import shutil
import JunctionArray

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

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

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
            
    try:
        if force == 'no': array_gene_annotations = JunctionArray.importGeneric(primary_gene_annotation_file)
        else: initiate_exception
    except Exception:
        import update; reload(update)
        update.downloadCurrentVersion(primary_gene_annotation_file,array_type,'txt')
        array_gene_annotations = JunctionArray.importGeneric(primary_gene_annotation_file)
            
    for array_geneid in array_gene_annotations:    
        t = array_gene_annotations[array_geneid]; description=t[0];entrez=t[1];symbol=t[2]
        if symbol in ensembl_symbol_db and len(symbol)>0 and len(array_geneid)>0:
            ens_geneid = ensembl_symbol_db[symbol]
            if len(ens_geneid)>0: array_ensembl[array_geneid]= ens_geneid
            
    try:
        if force == 'no': ensembl_probeset_db = importAltMouseLocationData(filename,array_ensembl,gene_info_db,test)
        else: initiate_exception
    except Exception:
        import update; reload(update)
        update.downloadCurrentVersion(filename,array_type,'txt')
        ensembl_probeset_db = importAltMouseLocationData(filename,array_ensembl,gene_info_db,test)
       
    print len(ensembl_probeset_db), "Genes inlcuded in",array_type,"location database"
    return ensembl_probeset_db

def importAltMouseLocationData(filename,array_ensembl,gene_info_db,test):
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

def getAnnotations(Species,array_type,force):
    """Annotate Affymetrix exon array data using files Ensembl data (sync'ed to genome release)."""
    global species; species = Species; global test; global test_cluster
    test = 'no'; test_cluster = ['G7055024']; data_type = 'mRNA'

    global ensembl_exon_db; global ensembl_exon_db; global exon_clusters; global exon_region_db
    ensembl_exon_db,ensembl_annot_db,exon_clusters,intron_clusters,exon_region_db,intron_retention_db,ucsc_splicing_annot_db,ens_transcript_db = EnsemblImport.getEnsemblAssociations(species,data_type,test)
    ensembl_probeset_db = importCriticalExonLocations(species,array_type,ensembl_exon_db,force) ###Get Pre-computed genomic locations for critical exons
    ensembl_probeset_db = ExonArrayEnsemblRules.annotateExons(ensembl_probeset_db,exon_clusters,ensembl_exon_db,exon_region_db,intron_retention_db,intron_clusters,ucsc_splicing_annot_db); constitutive_gene_db={}
    ExonArrayEnsemblRules.exportEnsemblLinkedProbesets(ensembl_probeset_db,species)
    print "\nCritical exon data exported coordinates, exon associations and splicing annotations exported..."
    
    ### Change filenames to reflect junction array type
    export_filename = 'AltDatabase/'+species+'/exon/'+species+'_Ensembl_probesets.txt'; ef=filepath(export_filename)
    export_replacement = string.replace(export_filename,'_probe','_'+array_type+'_probe')
    export_replacement = string.replace(export_replacement,'exon',array_type)
    er=filepath(export_replacement); shutil.copyfile(ef,er); os.remove(ef) ### Copy file to a new name

    ### Export full exon seqeunce for probesets/critical exons to replace the original incomplete sequence (used for miRNA analyses)
    #JunctionArray.reAnnotateCriticalExonSequences(species,array_type)
    
if __name__ == '__main__':
    m = 'Mm'; h = 'Hs'
    Species = m
    array_type = 'AltMouse' ###In theory, could be another type of junciton or combination array

    getAnnotations(Species,array_type,'no')
    
    