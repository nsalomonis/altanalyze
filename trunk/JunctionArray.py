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
        
def filterExpressionData(filename1,filename2,pre_filtered_db,constitutive_db):
    fn2=filepath(filename2)
    probeset_translation_db={}
    ###Import probeset number/id relationships (note: forced to use numeric IDs for Plier/Exact analysis)
    if analysis_method != 'rma':
        for line in open(fn2,'r').readlines():             
            data = cleanUpLine(line)
            probeset_number,probeset_id = string.split(data,'\t')
            probeset_translation_db[probeset_number]=probeset_id

    fn=filepath(filename1)
    exp_dbase={}; d = 0; x = 0
    ###Import expression data (non-log space)
    try:
        for line in open(fn,'r').readlines():             
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
    for line in open(fn,'rU').readlines():
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
    for line in open(fn,'rU').readlines():
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

def updateCriticalExonSequences(filename,ensembl_probeset_db):
    exon_seq_db_filename = filename[:-4]+'_updated.txt'
    fn=filepath(exon_seq_db_filename); data = open(fn,'w')
    
    critical_exon_seq_db={}
    for ens_gene in ensembl_probeset_db:
        for probe_data in ensembl_probeset_db[ens_gene]:
            exon_id,((probe_start,probe_stop,probeset_id,exon_class,transcript_clust),ed) = probe_data
            try: critical_exon_seq_db[probeset_id] = ed.ExonSeq()
            except AttributeError: null=[] ### Occurs when no sequence data is associated with exon (probesets without exon associations)
            
    fn1=filepath(filename); key_db = {}; x = 0
    verifyFile(filename)
    for line in open(fn1,'rU').readlines():
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
    data.close()
    print exon_seq_db_filename, 'exported....'
    
def identifyCriticalExonLocations(species,array_type):
    critical_exon_seq_db = importAnnotateCriticalExonSequences(species,array_type)
    dir = 'AltDatabase/ensembl/'+species+'/'; gene_seq_filename = dir+species+'_gene-seq-2000_flank'
    analysis_type = 'get_locations'  
    critical_exon_seq_db = EnsemblImport.import_sequence_data(gene_seq_filename,critical_exon_seq_db,species,analysis_type)
    exportCriticalExonLocations(species,array_type,critical_exon_seq_db)

def reAnnotateCriticalExonSequences(species,array_type):
    export_exon_filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_'+array_type+'_probesets.txt'        
    ensembl_probeset_db = ExonArrayEnsemblRules.reimportEnsemblProbesetsForSeqExtraction(export_exon_filename,'null',{})
    
    analysis_type = 'get_sequence'
    dir = 'AltDatabase/ensembl/'+species+'/'; gene_seq_filename = dir+species+'_gene-seq-2000_flank'
    ensembl_probeset_db = EnsemblImport.import_sequence_data(gene_seq_filename,ensembl_probeset_db,species,analysis_type)

    critical_exon_file = 'AltDatabase/'+species+'/'+ array_type + '/' + array_type+'_critical-exon-seq.txt'
    updateCriticalExonSequences(critical_exon_file, ensembl_probeset_db)
    
if __name__ == '__main__':
    """Module has methods for annotating Junction associated critical exon sequences with up-to-date genome coordinates and analysis options for
    junciton arrays from AnalyzeExpressionDatasets"""
    m = 'Mm'; h = 'Hs'; species = m; array_type = 'AltMouse' ###In theory, could be another type of junciton or combination array

    import_dir = '/AltDatabase/exon/'+species; expr_file_dir = 'R_expression_raw_data\exp.altmouse_es-eb.dabg.rma.txt'
    dagb_p = 1; Analysis_Method = 'rma'

    identifyCriticalExonLocations(species,array_type)
    JunctionArrayEnsemblRules.getAnnotations(species,array_type)
    
    ### Only needs to be run once, to update the original 
    #reAnnotateCriticalExonSequences(species,array_type); sys.exit() 
    
    
    #getAnnotations(expr_file_dir,dagb_p,Species,Analysis_Method)
