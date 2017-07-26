import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os.path
import unique
from stats_scripts import statistics
import copy
import time
from build_scripts import ExonSeqModule
import update

dirfile = unique

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".TXT" or entry[-3:] == ".fa":
            dir_list2.append(entry)
    return dir_list2

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

#################### Begin Analyzing Datasets #################### 
def importProbesetSeqeunces(filename,exon_db,species):
    print 'importing', filename
    fn=filepath(filename)
    probeset_seq_db={}; x = 0;count = 0
    for line in open(fn,'r').xreadlines():
        data, newline = string.split(line,'\n'); t = string.split(data,'\t')
        probeset = t[0]; sequence = t[-1]; sequence = string.upper(sequence)
        try:
            y = exon_db[probeset]; gene = y.GeneID(); y.SetExonSeq(sequence)
            try: probeset_seq_db[gene].append(y)
            except KeyError: probeset_seq_db[gene] = [y]
        except KeyError: null=[] ### Occurs if there is no Ensembl for the critical exon or the sequence is too short to analyze
    print len(probeset_seq_db), "length of gene - probeset sequence database"
    return probeset_seq_db
        
def importSplicingAnnotationDatabaseAndSequence(species,array_type,biotype):
    array_ens_db={}
    if array_type == 'AltMouse':
        filename = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'-Ensembl_relationships.txt'
        update.verifyFile(filename,array_type) ### Will force download if missing
        fn=filepath(filename); x = 0
        for line in open(fn,'r').xreadlines():
            data, newline = string.split(line,'\n'); t = string.split(data,'\t')
            if x==0: x=1
            else: 
                array_gene,ens_gene = t
                try: array_ens_db[array_gene].append(ens_gene)
                except KeyError: array_ens_db[array_gene]=[ens_gene]

    filename = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'_critical-junction-seq.txt'         
    fn=filepath(filename); probeset_seq_db={}; x = 0
    for line in open(fn,'r').xreadlines():
        data, newline = string.split(line,'\n'); t = string.split(data,'\t')
        if x==0: x=1
        else: 
            probeset,probeset_seq,junction_seq = t; junction_seq=string.replace(junction_seq,'|','')
            probeset_seq_db[probeset] = probeset_seq,junction_seq
            
    ###Import reciprocol junctions, so we can compare these directly instead of hits to nulls and combine with sequence data
    ###This short-cuts what we did in two function in ExonSeqModule with exon level data
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'_junction-comparisons.txt'
    fn=filepath(filename); probeset_gene_seq_db={}; x = 0
    for line in open(fn,'r').xreadlines():
        data, newline = string.split(line,'\n'); t = string.split(data,'\t')
        if x==0: x=1
        else: 
            array_gene,probeset1,probeset2,critical_exons = t #; critical_exons = string.split(critical_exons,'|')
            probesets = [probeset1,probeset2]
            if array_type == 'junction' or array_type == 'RNASeq': array_ens_db[array_gene]=[array_gene]
            if array_gene in array_ens_db:
                ensembl_gene_ids = array_ens_db[array_gene]
                for probeset_id in probesets:
                    if probeset_id in probeset_seq_db:
                        probeset_seq,junction_seq = probeset_seq_db[probeset_id]
                        if biotype == 'gene':
                            for ensembl_gene_id in ensembl_gene_ids:
                                probe_data = ExonSeqModule.JunctionDataSimple(probeset_id,ensembl_gene_id,array_gene,probesets,critical_exons)
                                probe_data.SetExonSeq(probeset_seq)
                                probe_data.SetJunctionSeq(junction_seq)
                                try: probeset_gene_seq_db[ensembl_gene_id].append(probe_data)
                                except KeyError: probeset_gene_seq_db[ensembl_gene_id] = [probe_data]
                        else: ### Used for probeset annotations downstream of sequence alignment in LinkEST, analagous to exon_db for exon analyses
                            probe_data = ExonSeqModule.JunctionDataSimple(probeset_id,ensembl_gene_ids,array_gene,probesets,critical_exons)
                            probe_data.SetExonSeq(probeset_seq)
                            probe_data.SetJunctionSeq(junction_seq)                            
                            probeset_gene_seq_db[probeset_id] = probe_data                
    print len(probeset_gene_seq_db),"genes with probeset sequence associated"
    return probeset_gene_seq_db

def getParametersAndExecute(probeset_seq_file,array_type,species,data_type):
    if data_type == 'critical-exons':
        if array_type == 'RNASeq': probeset_annotations_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_exons.txt'
        else: probeset_annotations_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_'+array_type+'_probesets.txt'
        ###Import probe-level associations
        exon_db = ExonSeqModule.importSplicingAnnotationDatabase(probeset_annotations_file,array_type)
        start_time = time.time()
        probeset_seq_db = importProbesetSeqeunces(probeset_seq_file,exon_db,species)  ###Do this locally with a function that works on tab-delimited as opposed to fasta sequences (exon array)
        end_time = time.time(); time_diff = int(end_time-start_time)
    elif data_type == 'junctions':
        start_time = time.time(); biotype = 'gene' ### Indicates whether to store information at the level of genes or probesets
        probeset_seq_db = importSplicingAnnotationDatabaseAndSequence(species,array_type,biotype)
        end_time = time.time(); time_diff = int(end_time-start_time)
    print "Analyses finished in %d seconds" % time_diff
    return probeset_seq_db

def runProgram(Species,Array_type,mir_source,stringency,Force):
    global species; global array_type; global force
    process_microRNA_predictions = 'yes'

    species = Species; array_type = Array_type; force = Force
    
    import_dir = '/AltDatabase/'+species+'/'+array_type
    filedir = import_dir[1:]+'/'
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    probeset_seq_file=''
    for input_file in dir_list:    #loop through each file in the directory to  results
        if 'critical-exon-seq_updated' in input_file: probeset_seq_file = filedir+input_file
        elif 'critical-exon-seq' in input_file: probeset_seq_file2 = filedir+input_file
    if len(probeset_seq_file)==0: probeset_seq_file=probeset_seq_file2
        
    data_type = 'critical-exons'
    try: splice_event_db = getParametersAndExecute(probeset_seq_file,array_type,species,data_type)
    except UnboundLocalError:
        probeset_seq_file = 'AltDatabase/'+species+'/'+array_type+'/'+array_type+'_critical-exon-seq_updated.txt'
        update.downloadCurrentVersion(probeset_seq_file,array_type,'txt')
        splice_event_db = getParametersAndExecute(probeset_seq_file,array_type,species,data_type)
        
    if process_microRNA_predictions == 'yes':
        print 'stringency:',stringency
        try:
            ensembl_mirna_db = ExonSeqModule.importmiRNATargetPredictionsAdvanced(species)
            ExonSeqModule.alignmiRNAData(array_type,mir_source,species,stringency,ensembl_mirna_db,splice_event_db)
        except Exception: pass
        
if __name__ == '__main__':
    species = 'Mm'; array_type = 'junction'
    process_microRNA_predictions = 'yes'
    mir_source = 'multiple'
    force = 'yes'
    
    runProgram(species,array_type,mir_source,'lax',force)
    runProgram(species,array_type,mir_source,'strict',force)
