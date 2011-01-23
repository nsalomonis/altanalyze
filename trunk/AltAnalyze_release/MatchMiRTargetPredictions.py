###MatchTargetPredictions
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
import export

############ File Import Functions #############

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".all" or entry[-5:] == ".data" or entry[-3:] == ".fa":
            dir_list2.append(entry)
    return dir_list2

def returnDirectories(sub_dir):
    dir=os.path.dirname(dirfile.__file__)
    dir_list = os.listdir(dir + sub_dir)
    ###Below code used to prevent FILE names from being included
    dir_list2 = []
    for entry in dir_list:
        if "." not in entry: dir_list2.append(entry)
    return dir_list2

class GrabFiles:
    def setdirectory(self,value): self.data = value
    def display(self): print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        files = getDirectoryFiles(self.data,search_term)
        if len(files)<1: print 'files not found'
        return files
    def returndirectory(self):
        dir_list = getAllDirectoryFiles(self.data)
        return dir_list

def getEnsemblAnnotations(filename,symbol_ensembl_db):
    fn=filepath(filename)
    print 'importing', filename
    verifyFile(filename,species) ### Makes sure file is local and if not downloads.
    for line in open(fn,'r').xreadlines():        
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t'); ensembl=t[0]
        try: symbol = t[2]
        except IndexError: symbol = ''
        if len(symbol)>1 and len(ensembl)>1:
            try: symbol_ensembl_db[string.upper(symbol)].append(ensembl)
            except KeyError: symbol_ensembl_db[string.upper(symbol)] = [ensembl]
    return symbol_ensembl_db
    
def getCurrentEnsembls(symbol_ensembl_old):
    ### Get Ensembl gene annotaitons for the current version of Ensembl only
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations_simple.txt'
    symbol_ensembl_current = getEnsemblAnnotations(filename,{})
    
    all_current_ensembls={}
    for symbol in symbol_ensembl_current:
        if len(symbol_ensembl_current[symbol][0])<3: print symbol_ensembl_current[symbol][0]
        all_current_ensembls[symbol_ensembl_current[symbol][0]] = []
    
    ### Augment with UniProt Symbols and Aliases
    uniprot_symbol_ensembl = importUniProtAnnotations()
    
    for symbol in uniprot_symbol_ensembl:
        ensembls = uniprot_symbol_ensembl[symbol]
        if symbol not in symbol_ensembl_current:
            for ensembl in ensembls:
                if ensembl in all_current_ensembls:
                    try: symbol_ensembl_current[symbol].append(ensembl)
                    except Exception: symbol_ensembl_current[symbol] = [ensembl]
    for symbol in symbol_ensembl_old:
        ensembls = symbol_ensembl_old[symbol]
        if symbol not in symbol_ensembl_current:
            for ensembl in ensembls:
                if ensembl in all_current_ensembls:
                    try: symbol_ensembl_current[symbol].append(ensembl)
                    except Exception: symbol_ensembl_current[symbol] = [ensembl]
    return symbol_ensembl_current

def processEnsemblAnnotations():
    global redundant_ensembl_by_build; redundant_ensembl_by_build={}; symbol_ensembl={}
    
    ###This is done for archived gene annotations from Ensembl         
    filename = 'AltDatabase/miRBS/'+species+'/'+species+'_Ensembl-annotations_simple.txt'
    symbol_ensembl = getEnsemblAnnotations(filename,symbol_ensembl)
    
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations_simple.txt'
    symbol_ensembl_current = getCurrentEnsembls(symbol_ensembl)
    
    for symbol in symbol_ensembl_current:
        for ensembl in symbol_ensembl_current[symbol]:
            try: symbol_ensembl[symbol].append(ensembl)
            except Exception: symbol_ensembl[symbol] = [ensembl]
                        
    for symbol in symbol_ensembl:
      if len(symbol_ensembl[symbol])>1: ###Thus there are more than one Ensembl gene IDs that correspond... probably due to versioning (which we've seen)
        for ensembl in symbol_ensembl[symbol]:
            for ensembl2 in symbol_ensembl[symbol]:
                if ensembl != ensembl2:
                    try: redundant_ensembl_by_build[ensembl].append(ensembl2)
                    except KeyError: redundant_ensembl_by_build[ensembl] = [ensembl2]
    print 'len(symbol_ensembl)',len(symbol_ensembl)

    for transcript in ens_gene_to_transcript:
        if len(ens_gene_to_transcript[transcript])>1:
            for ensembl in ens_gene_to_transcript[transcript]:
                for ensembl2 in ens_gene_to_transcript[transcript]:
                    if ensembl != ensembl2:
                        try: redundant_ensembl_by_build[ensembl].append(ensembl2)
                        except KeyError: redundant_ensembl_by_build[ensembl] = [ensembl2]
    redundant_ensembl_by_build = eliminateRedundant(redundant_ensembl_by_build)
    return symbol_ensembl,symbol_ensembl_current

def getAllDirectoryFiles(import_dir):
    all_files = []
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        data_dir = import_dir[1:]+'/'+data
        all_files.append(data_dir)
    return all_files

def getDirectoryFiles(import_dir,search_term):
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    matches=[]
    for data in dir_list:    #loop through each file in the directory to output results
        data_dir = import_dir[1:]+'/'+data
        if search_term in data_dir: matches.append(data_dir)
    return matches

class MicroRNATargetData:
    def __init__(self,gene,gene_symbol,mir,mir_sequences,source):
        self._geneid = gene; self._symbol = gene_symbol; self._mir = mir; self._source = source; self._mir_sequences = mir_sequences
    def MicroRNA(self): return self._mir
    def GeneID(self): return self._geneid
    def Symbol(self): return self._symbol
    def Source(self): return self._source
    def Sequences(self): return self._mir_sequences
    def Report(self):
        output = self.GeneID()
        return output
    def __repr__(self): return self.Report()

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def mirHitImport():
    try: 
        filename = 'AltDatabase/miRBS/'+species+'/'+'hits.txt'
        print 'parsing', filename
        fn=filepath(filename); x=0
        for line in open(fn,'rU').xreadlines():         
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0: x=1
            else:
                mir = t[0]; fold = t[1]; mir_hit_db[mir] = fold
    except IOError: mir_hit_db={}

def importExpressionData():
    global upregulated; global downregulated
    try:
        filename = 'AltDatabase/miRBS/'+species+'/'+'expression-data.txt'
        print 'parsing', filename
        fn=filepath(filename); x=0; upregulated={}; downregulated={}
        for line in open(fn,'rU').xreadlines():         
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0: x=1
            else:
                gene = t[0]; log_fold_change = float(t[-2]); ttest = float(t[-1])
                if ttest<0.05 and log_fold_change>1: upregulated[gene] = log_fold_change,ttest
                if ttest<0.05 and log_fold_change<-1: downregulated[gene] = log_fold_change,ttest
    except IOError: upregulated={}; downregulated={}
            
def pictarImport(parse_sequences,type,added):
    """Annotations originally from the file: ng1536-S3.xls, posted as supplementary data at:
    http://www.nature.com/ng/journal/v37/n5/suppinfo/ng1536_S1.html. The file being parsed here has been pre-matched to Ensembl IDs
    using the ExonModule of LinkEST, for human."""
    mir_sequences=[]
    if species == 'Mm': filename = 'AltDatabase/miRBS/'+species+'/'+'pictar-target-annotated.txt'; tax = '10090'
    if species == 'Rn': filename = 'AltDatabase/miRBS/'+'Mm'+'/'+'pictar-target-annotated.txt'; tax = '10116'
    #if species == 'Hs': filename = 'AltDatabase/miRBS/'+species+'/'+'pictar-conserved-targets-2005.txt'; tax = '9606'
    if type == 'pre-computed':
        if species == 'Hs': filename = 'AltDatabase/miRBS/'+species+'/'+'pictar-conserved-targets-2005.txt'; tax = '9606'
    else:
        if species == 'Hs': filename = 'AltDatabase/miRBS/'+'Mm'+'/'+'pictar-target-annotated.txt'; tax = '9606'
 
    print 'parsing', filename; count=0
    print 'len(symbol_ensembl)', len(symbol_ensembl)
    verifyFile(filename,species) ### Makes sure file is local and if not downloads.
    fn=filepath(filename); x=1
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            if species == 'Hs':
                if type == 'pre-computed':
                    ensembl_geneid, mir, mir_sequences = t; ensembl_geneids = [ensembl_geneid]
                else:
                    symbol=string.upper(t[2]);mir=t[6];mir_sequences=t[11]
                    if symbol in symbol_ensembl and len(symbol)>0: ensembl_geneids=symbol_ensembl[symbol]
                    else: ensembl_geneids=['']                    
            if species == 'Mm':
                mm_symbol=string.upper(t[3]);mir=t[6];mir_sequences=t[11]; mir = string.replace(mir,'hsa','mmu')
                if mm_symbol in symbol_ensembl and len(mm_symbol)>0: ensembl_geneids=symbol_ensembl[mm_symbol]
                else: ensembl_geneids=['']
            if species == 'Rn':
                mm_symbol=string.upper(t[3]);mir=t[6];mir_sequences=t[11]; mir = string.replace(mir,'hsa','rno')
                if mm_symbol in symbol_ensembl and len(mm_symbol)>0: ensembl_geneids=symbol_ensembl[mm_symbol]
                else: ensembl_geneids=['']
            for ensembl_geneid in ensembl_geneids:
                if len(ensembl_geneid)>1 and (ensembl_geneid,mir) not in added:
                    if parse_sequences == 'yes':
                        if (mir,ensembl_geneid) in combined_results:
                            combined_results[(mir,ensembl_geneid)].append(string.upper(mir_sequences)); count+=1
                    else:
                        #if count < 800 and '-125b' in mir: print ensembl_geneid, mir, mm_symbol; count+=1
                        #elif count>799: kill
                        y = MicroRNATargetData(ensembl_geneid,'',mir,mir_sequences,'pictar'); count+=1
                        try: microRNA_target_db[mir].append(y)
                        except KeyError: microRNA_target_db[mir] = [y]
                        added[(ensembl_geneid,mir)]=[]
                    
    print count, 'miRNA-target relationships added for PicTar'
    return added

def sangerImport(parse_sequences):
    """"Sanger center (miRBase) sequence was provided as a custom (requested) dump of their v5 target predictions
    (http://microrna.sanger.ac.uk/targets/v5/), containing Ensembl gene IDs, microRNA names, and putative target
    sequences, specific for either mouse or human. Mouse was requested in late 2005 whereas human in late 2007.
    These same annotation files, missing the actual target sequence but containing an ENS transcript and coordinate
    locations for that build (allowing seqeunce extraction with the appropriate Ensembl build) exist at:
    http://microrna.sanger.ac.uk/cgi-bin/targets/v5/download.pl"""
    
    if species == 'Hs': filename = 'AltDatabase/miRBS/'+species+'/'+'mirbase-v5_homo_sapiens.mirna.txt'; prefix = 'hsa-'
    if species == 'Rn': filename = 'AltDatabase/miRBS/'+species+'/'+'sanger_miR_target_predictions.txt'; prefix = 'rno-'
    if species == 'Mm': filename = 'AltDatabase/miRBS/'+species+'/'+'sanger_miR_target_predictions.txt'; prefix = 'mmu-'    

    print 'parsing', filename; count=0
    fn=filepath(filename); x=1; mir_sequences=[]
    verifyFile(filename,species) ### Makes sure file is local and if not downloads.
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            if species == 'Hs':
                try: mir = t[1]; ens_transcript = t[2]; ensembl_geneid = t[17]; mir_sequences = string.upper(t[14])
                except IndexError: print line;kill
            elif species == 'Mm':
                ens_transcript,mir,mir_sequences = t
                if ens_transcript in ens_gene_to_transcript: ensembl_geneids = ens_gene_to_transcript[ens_transcript]; ensembl_geneid = ensembl_geneids[0]
                else: ensembl_geneids=[]
            elif species == 'Rn':
                ensembl_geneid,mir,mir_sequences = t
                mir_sequences = string.lower(mir_sequences); mir = string.replace(mir,'hsa','rno'); mir = string.replace(mir,'mmu','rno')
                ensembl_geneids=[ensembl_geneid]
            geneid_ls=[]
            if ensembl_geneid in redundant_ensembl_by_build: ###Thus there are redundant geneids
                geneid_ls = redundant_ensembl_by_build[ensembl_geneid]+[ensembl_geneid]
            else: geneid_ls = [ensembl_geneid]
            if species == 'Hs':
                if ens_transcript in ens_gene_to_transcript: geneid_ls+= ens_gene_to_transcript[ens_transcript]
            geneid_ls = unique.unique(geneid_ls)
            if len(geneid_ls) == 1 and geneid_ls[0]=='': null =[] ###not a valid gene
            elif prefix in mir:
                for ensembl_geneid in geneid_ls:
                    if parse_sequences == 'yes':
                        if (mir,ensembl_geneid) in combined_results:
                            mir_sequences = string.replace(mir_sequences,'-',''); mir_sequences = string.replace(mir_sequences,'=',''); count+=1
                            combined_results[(mir,ensembl_geneid)].append(string.upper(mir_sequences))
                    else:
                        if prefix in mir:
                            y = MicroRNATargetData(ensembl_geneid,'',mir,mir_sequences,'mirbase'); count+=1
                            try: microRNA_target_db[mir].append(y)
                            except KeyError: microRNA_target_db[mir] = [y]
    print count, 'miRNA-target relationships added for mirbase'

def downloadFile(file_type):
    import UI; import update
    file_location_defaults = UI.importDefaultFileLocations()
    try:
        fld = file_location_defaults[file_type]
        url = fld.Location()
    except Exception:
        for fl in fld:
            if species in fl.Species(): url = fl.Location()
    if 'Target' in file_type: output_dir = 'AltDatabase/miRBS/'
    else: output_dir = 'AltDatabase/miRBS/'+species + '/'
    gz_filepath, status = update.download(url,output_dir,'')
    if status == 'not-removed':
        try: os.remove(gz_filepath) ### Not sure why this works now and not before
        except Exception: status = status
    
    filepath = string.replace(gz_filepath,'.zip','.txt')
    filepath = string.replace(filepath,'.gz','.txt')
    filepath = string.replace(filepath,'.txt.txt','.txt')
    
    return filepath
                        
def TargetScanImport(parse_sequences,force):
    """The TargetScan data is currently extracted from a cross-species conserved family file. This file only contains
    gene symbol, microRNA name and 3'UTR seed locations."""
    if species == 'Mm': tax = '10090'; prefix = 'mmu-'
    if species == 'Hs': tax = '9606'; prefix = 'hsa-'
    if species == 'Rn': tax = '10116'; prefix = 'rno-'

    global l
    if force == 'yes':    
        ### Then download the latest annotations and sequences
        target_scan_target_file = downloadFile('TargetScanGenes')
        target_scan_sequence_file = downloadFile('TargetScanSequences')

        ### Cross-species TargetScan file with UTR seqeunces for all genes with reported targets in the conserved family file
        ### Although this file includes valid sequence data that appears to match up to the target file, the target file
        ### appears to only list the seed seqeunce location (UTR start and stop) and not the full binding sequence and thus
        ### is not ammenable to probe set alignment.
        print 'parsing', target_scan_sequence_file
        fn=filepath(target_scan_sequence_file); x=0; target_scan_gene_utr_seq={}
        for line in open(fn,'rU').xreadlines():         
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0: x=1
            else:
                symbol = string.upper(t[2]); tax_id = t[3]; utr_seq = t[4]
                if tax_id == tax:
                    utr_seq_no_gaps = string.replace(utr_seq,'-','')
                    utr_seq_no_gaps = string.replace(utr_seq_no_gaps,'U','T')
                    if symbol in symbol_ensembl_current and len(utr_seq_no_gaps)>0:
                        target_scan_gene_utr_seq[symbol] = utr_seq_no_gaps
        print 'UTR sequence for',len(target_scan_gene_utr_seq),'TargetScan genes stored in memory.'
    else:    
        ### Archival file
        target_scan_target_file = 'AltDatabase/miRBS/target-scan-Conserved_Family_Conserved_Targets_Info.txt'
        target_scan_gene_utr_seq={}
        
    mir_sequences = []; count=0
    print 'parsing', target_scan_target_file
    verifyFile(target_scan_target_file,species) ### Makes sure file is local and if not downloads.
    fn=filepath(target_scan_target_file); x=0; k=[]; l=[]
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            if force == 'yes':
                mir = t[0]; geneid = t[1]; gene_symbol = string.upper(t[2]); taxid = t[3]; utr_start = int(t[4]); utr_end  = int(t[5])
            else:
                ### Old format
                mir = t[0]; gene_symbol = string.upper(t[1]); taxid = t[2]; utr_start = t[3]; utr_end = t[4]
            if '/' in mir:
                mir_list=[]
                mirs = string.split(mir,'/')
                for mirid in mirs[1:]:
                    mirid = 'miR-'+mirid
                    mir_list.append(mirid)
                mir_list.append(mirs[0])
            else: mir_list = [mir]
            if taxid == tax: ###human
                #target_scan_gene_utr_seq[symbol] = utr_seq_no_gaps
                if gene_symbol in symbol_ensembl_current: ensembl_geneids = symbol_ensembl_current[gene_symbol]; proceed = 'yes'; k.append(gene_symbol)
                else: proceed = 'no'; l.append(gene_symbol)
                if gene_symbol in target_scan_gene_utr_seq:
                    ### TargetScan provides the core, while processed miRs are typically 22nt - seems to approximate other databases better
                    adj_start = utr_start-15
                    if adj_start < 0: adj_start=0
                    mir_sequences = target_scan_gene_utr_seq[gene_symbol][adj_start:utr_end+1]
                    #if string.lower(gene_symbol) == 'tns3' and mir == 'miR-182': print mir,gene_symbol,taxid,utr_start,utr_end,mir_sequences
                else: mir_sequences=[]
                ###Already multiple geneids associated with each symbol so don't need to worry about renundancy
                if proceed == 'yes':
                    for ensembl_geneid in ensembl_geneids:
                        for mir in mir_list:
                            #if ensembl_geneid == 'ENSG00000137815' and mir == 'miR-214': print mir,gene_symbol,taxid,utr_start,utr_end,mir_sequences,target_scan_gene_utr_seq[gene_symbol];sys.exit()
                            if parse_sequences == 'yes':
                                if (prefix+mir,ensembl_geneid) in combined_results:
                                    combined_results[(prefix+mir,ensembl_geneid)].append(mir_sequences); count+=1
                            else:
                                #if ensembl_geneid == 'ENSMUSG00000029467': print mir
                                y = MicroRNATargetData(ensembl_geneid,gene_symbol,mir_sequences,prefix+mir,'TargetScan')
                                count+=1
                                try: microRNA_target_db[prefix+mir].append(y)
                                except KeyError: microRNA_target_db[prefix+mir] = [y]
    k = unique.unique(k); l = unique.unique(l)
    print 'ensembls-found:',len(k),', not found:',len(l)
    print l[:10]
    print count, 'miRNA-target relationships added for TargetScan'
    
def mirandaImport(parse_sequences,force):
    """Miranda data is avaialble from two file types from different websites. The first is human-centric with multi-species
    target alignment information organized by Ensembl gene ID (http://cbio.mskcc.org/research/sander/data/miRNA2003/mammalian/index.html).
    A larger set of associations was also pulled from species specific files (http://www.microrna.org/microrna/getDownloads.do),
    where gene symbol was related to Ensembl gene. Both files provided target microRNA sequence."""
    
    if force == 'yes':
        ### Then download the latest annotations and sequences
        filename = downloadFile('miRanda')

        print 'parsing', filename; count=0; null_count=[]
        fn=filepath(filename); x=1; mir_sequences=[]
        verifyFile(filename,species) ### Makes sure file is local and if not downloads.
        for line in open(fn,'rU').xreadlines():         
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0: x=1
            else:
                symbol = string.upper(t[3]); mir = t[1]; entrez_gene = t[2]; mir_sequences = string.upper(t[8])
                mir_sequences = string.replace(mir_sequences,'-',''); mir_sequences = string.replace(mir_sequences,'=','')
                mir_sequences = string.replace(mir_sequences,'U','T')
                ensembl_gene_ids = []
                if symbol in symbol_ensembl_current:
                    ensembl_gene_ids = symbol_ensembl_current[symbol]
                else: ensembl_gene_ids=[]; null_count.append(symbol)

                if len(ensembl_gene_ids) > 0:
                    for ensembl_geneid in ensembl_gene_ids:
                        if parse_sequences == 'yes':
                            if (mir,ensembl_geneid) in combined_results:
                                combined_results[(mir,ensembl_geneid)].append(string.upper(mir_sequences)); count+=1
                        else:
                            y = MicroRNATargetData(ensembl_geneid,'',mir,mir_sequences,'miranda'); count+=1
                            try: microRNA_target_db[mir].append(y)
                            except KeyError: microRNA_target_db[mir] = [y]
        print count, 'miRNA-target relationships added for miRanda'
        null_count = unique.unique(null_count)
        print len(null_count), 'missing symbols',null_count[:10]
    else:
        filename = 'AltDatabase/miRBS/miranda-by_gene140.txt'
        print 'len(ens_gene_to_transcript)',len(ens_gene_to_transcript)
        print 'parsing', filename; count=0; bad_line=0
        fn=filepath(filename); x=1 ###since the first line contains data
        if species == 'Hs': mir_type = 'hsa'; ens_transcript_prefix = 'ENST'; trans_index = 0
        if species == 'Mm': mir_type = 'mmu'; ens_transcript_prefix = 'ENSMUST'; trans_index = 2
        if species == 'Rn': mir_type = 'rno'; ens_transcript_prefix = 'ENSRNOT'; trans_index = 3
        verifyFile(filename,species) ### Makes sure file is local and if not downloads.
        for line in open(fn,'rU').xreadlines():         
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0: x=1
            else:
                try:
                    ens_transcript = t[trans_index]; mirs = t[8]; mir_seq_db = {}
                    if parse_sequences == 'yes':
                        mir_sequences = t[9:]
                        for mir_seq_data in mir_sequences:
                            mir_seq_data = string.split(mir_seq_data,"5'_")
                            if mir_type in mir_seq_data[0]: ###require that the microRNA must be human
                                mir_info = string.split(mir_seq_data[0],":"); microRNA = mir_info[1]; mir_info = string.split(microRNA," "); microRNA = mir_info[0]
                                for info in mir_seq_data[1:]:
                                    if ens_transcript_prefix in info: ###again, ensure that it's matched to a human sequence (human mirs can be matched to mouse sequence)
                                        mir_seq = string.split(info,"_3':"); mir_seq = mir_seq[0]; mir_seq=string.replace(mir_seq,'-',''); mir_seq = string.upper(mir_seq)
                                        try: mir_seq_db[microRNA].append(mir_seq)
                                        except KeyError: mir_seq_db[microRNA] = [mir_seq]
                        mir_seq_db = eliminateRedundant(mir_seq_db)
                    multi_species_mir_list = string.split(mirs,' '); multi_species_mir_list2=[]
                    for mir in multi_species_mir_list:
                        if mir_type in mir:  multi_species_mir_list2.append(mir)
                        else: null = ''
                    ens_transcript = string.replace(ens_transcript,'_P','')
                    try:
                        geneid = ens_gene_to_transcript[ens_transcript][0]
                        if len(geneid)<1: geneid = t[1]; geneid_data = string.split(geneid,' '); geneid = geneid_data[0]
                    except KeyError: geneid = t[1]; geneid_data = string.split(geneid,' '); geneid = geneid_data[0]
    
                    if len(geneid)>0:
                        if geneid in redundant_ensembl_by_build: ###Thus there are redundant geneids
                            geneid_ls = redundant_ensembl_by_build[geneid]+[geneid]
                        else: geneid_ls = [geneid]
                        for geneid in geneid_ls:
                            for mir in multi_species_mir_list2:
                                if parse_sequences == 'yes':
                                    if mir in mir_seq_db:
                                        if (mir,geneid) in combined_results:
                                            mir_sequences = mir_seq_db[mir]
                                            for mir_sequence in mir_sequences:
                                                combined_results[(mir,geneid)].append(string.upper(mir_sequence)); count+=1
                                else:
                                    mir_sequences = []
                                    y = MicroRNATargetData(geneid,'',mir,mir_sequences,'miranda'); count+=1    
                                    try: microRNA_target_db[mir].append(y)
                                    except KeyError: microRNA_target_db[mir] = [y]
                except IndexError: bad_line+=1
    
        if species == 'Mm':
            filename = 'AltDatabase/miRBS/'+species+'/'+'mouse_predictions_4dwnload.txt'
            print 'parsing', filename; count2=0; added={}; added2={}
            fn=filepath(filename); x=1 ###since the first line contains data
            if species == 'Hs': mir_type = 'hsa'; ens_transcript_prefix = 'ENST'; trans_index = 0
            if species == 'Mm': mir_type = 'mmu'; ens_transcript_prefix = 'ENSMUST'; trans_index = 2
            if species == 'Rn': mir_type = 'rno'; ens_transcript_prefix = 'ENSRNOT'; trans_index = 2
            verifyFile(filename,species) ### Makes sure file is local and if not downloads.
            for line in open(fn,'rU').xreadlines():         
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if x==0: x=1
                else:
                    symbol = string.upper(t[2]); mir = t[4]; mir_sequence = t[7]; mir_sequence = string.replace(string.upper(mir_sequence),'U','T')
                    if symbol in symbol_ensembl and len(symbol)>0: geneid_ls = symbol_ensembl[symbol]
                    else: geneid_ls = []
                    mir_sequence=string.replace(mir_sequence,'-',''); mir_sequence=string.replace(mir_sequence,'=','')
                    for ensembl_geneid in geneid_ls:
                        if parse_sequences == 'yes':
                            if (ensembl_geneid,mir,mir_sequence) not in added2:
                                combined_results[(mir,ensembl_geneid)].append(string.upper(string.upper(mir_sequence))); count2+=1
                                added2[(ensembl_geneid,mir,mir_sequence)]=[]
                        else:
                            if (ensembl_geneid,mir) not in added:
                                y = MicroRNATargetData(ensembl_geneid,'',mir,[mir_sequence],'miranda'); count2+=1
                                try: microRNA_target_db[mir].append(y)
                                except KeyError: microRNA_target_db[mir] = [y]
                                added[(ensembl_geneid,mir)]=[]
            print count2, 'miRNA-target relationships added for miranda2'
                                
        print bad_line, "bad lines"
        print count, 'miRNA-target relationships added for miranda'
    
def importEnsTranscriptAssociations(ens_gene_to_transcript,type):
    ###This function is used to extract out EnsExon to EnsTranscript relationships to find out directly
    ###which probesets associate with which transcripts and then which proteins
    if type == 'current':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    else: filename = 'AltDatabase/miRBS/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    print 'parsing', filename
    fn=filepath(filename)
    verifyFile(filename,species)
    for line in open(fn,'rU').readlines():
        data,null = string.split(line,'\n')
        t = string.split(data,'\t')
        ens_geneid,chr,strand,exon_start,exon_end,ens_exonid,constitutive_exon,ens_transcript_id = t
        try: ens_gene_to_transcript[ens_transcript_id].append(ens_geneid)
        except KeyError:ens_gene_to_transcript[ens_transcript_id] = [ens_geneid]
    ens_gene_to_transcript = eliminateRedundant(ens_gene_to_transcript)
    print 'len(ens_gene_to_transcript)',len(ens_gene_to_transcript)
    return ens_gene_to_transcript

def findMirTargetOverlaps():
    print 'Number of microRNAs to be analyzed:',len(mir_hit_db)
    print "Number of microRNAs which have predicted targets:", len(microRNA_target_db)
    combined_output = 'AltDatabase/'+species+'/SequenceData/'+'miRBS-combined_gene-targets.txt'
    data1 = export.ExportFile(combined_output)
    mir_gene_mult_hits={}
    for mir in microRNA_target_db:
        if mir in mir_hit_db:
            output_file = 'AltDatabase/'+species+'/SequenceData/'+mir+'_gene-targets.txt'; output_file = string.replace(output_file,'*','-')
            data = export.ExportFile(output_file); delete = {}
            data = export.ExportFile(output_file); delete = {}
            title = ['TargetID','Evidence']; title = string.join(title,'\t')+'\n'; data.write(title)
        matches=[]; source_mir_db = {}; hit=0
        source_db={}
        for tg in microRNA_target_db[mir]:
            try:
                try: source_db[tg.GeneID()].append(tg.Source())
                except KeyError: source_db[tg.GeneID()] = [tg.Source()]
            except TypeError: print tg.Source(),tg.GeneID();kill
        source_db = eliminateRedundant(source_db) 
        for gene in source_db:
            for source in source_db[gene]: source_mir_db[source]=[]
            sources = string.join(source_db[gene],'|')
            y = MicroRNATargetData(gene,'',mir,'',sources)
            if mir in mir_hit_db:
                try: mir_gene_mult_hits[mir].append(y)
                except KeyError: mir_gene_mult_hits[mir] = [y]
                else: delete[mir] = []
            values = [gene,sources]; values = string.join(values,'\t')+'\n'
            values2 = [mir,gene,sources]; values2 = string.join(values2,'\t')+'\n'
            if len(source_db[gene])>1:
                matches.append(gene)
                if mir in mir_hit_db: data.write(values); hit+=1
            data1.write(values2)
        if mir in mir_hit_db:
            print len(source_db),'genes associated with',mir+'.',len(source_mir_db),'sources linked to gene. Targets with more than one associated Dbase:', len(matches)
        if mir in mir_hit_db:
            data.close()
            if len(source_mir_db)<4 or hit<5:
                try: del mir_gene_mult_hits[mir]
                except KeyError: null = []
                os.remove(fn)
    data1.close()
    for mir in mir_gene_mult_hits:
        mir_fold = mir_hit_db[mir]
        if mir_fold>0: mir_direction = 'up'
        else: mir_direction = 'down'
        output_file = 'AltDatabase/regulated/'+mir+'_regualted-gene-targets.txt'; output_file = string.replace(output_file,'*','-')
        fn=filepath(output_file);data = open(fn,'w')
        title = ['TargetID','Evidence','CSd40_vs_ESCs-log_fold','CSd40_vs_ESCs-ttest']; title = string.join(title,'\t')+'\n'; data.write(title)
        hit = 0
        if mir in mir_gene_mult_hits:
            for y in mir_gene_mult_hits[mir]:
                gene = y.GeneID(); sources = y.Source()
                if mir_direction == 'up': ###thus look for down-regulated targets
                    if gene in downregulated:
                        hit +=1
                        log_fold,ttest = downregulated[gene]
                        values = [gene,sources,str(log_fold),str(ttest)]; values = string.join(values,'\t')+'\n';data.write(values)
                if mir_direction == 'down': ###thus look for up-regulated targets
                    if gene in upregulated:
                        log_fold,ttest = upregulated[gene]
                        values = [gene,sources,str(log_fold),str(ttest)]; values = string.join(values,'\t')+'\n';data.write(values)
            
            print hit,'regualted target genes associated with',mir
        data.close()
        if hit<5: os.remove(fn)
        
def eliminateRedundant(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def importUniProtAnnotations():
    #g = GrabFiles(); g.setdirectory('/AltDatabase/miRBS/'+species); filenames = g.searchdirectory('UniProt')
    filename = 'AltDatabase/uniprot/'+species+'/uniprot_sequence.txt'
    fn=filepath(filename); uniprot_symbol_ensembl={}
    for line in open(fn,'rU').readlines():             
        data,null = string.split(line,'\n')  #remove endline
        t = string.split(data,'\t')  #remove endline
        symbols = t[3]; ###ASPN; ORFNames=RP11-77D6.3-002, hCG_1784540
        symbols = string.split(symbols,';'); primary_symbol = symbols[0]
        if len(symbols)>1:
            null,symbols = string.split(symbols[1],'=')
            symbols = string.split(symbols,', '); symbols = symbols+[primary_symbol]
        else: symbols = [primary_symbol]
        ensembl_ids = t[4]; ensembl_ids = string.split(ensembl_ids,",")
        for symbol in symbols:
            if len(symbol)>0:
                if len(ensembl_ids)>0: uniprot_symbol_ensembl[string.upper(symbol)] = ensembl_ids #; print symbol,ensembl_ids;kil
    return uniprot_symbol_ensembl

def exportCombinedMirResultSequences():
    combined_input = 'AltDatabase/'+species+'/SequenceData/'+'miRBS-combined_gene-targets.txt'
    parse_sequences = 'yes'; global combined_results; combined_results={}
    fn=filepath(combined_input); combined_results2={} ###the second is for storing source information
    for line in open(fn,'rU').xreadlines():
        data,null = string.split(line,'\n')  #remove endline
        t = string.split(data,'\t')  #remove endline
        mir, gene, sources = t
        combined_results[mir,gene] = []
        combined_results2[mir,gene] = sources
    microRNA_target_db = {}; mir_hit_db = {}
    mirandaImport(parse_sequences,'yes')
    TargetScanImport(parse_sequences,'yes')
    try: sangerImport(parse_sequences)
    except Exception: null=[]
    added = pictarImport(parse_sequences,'pre-computed',{})
    added = pictarImport(parse_sequences,'symbol-based',added)
    output_file = 'AltDatabase/'+species+'/SequenceData/'+'miRBS-combined_gene-target-sequences.txt'
    combined_results = eliminateRedundant(combined_results); fn=filepath(output_file);data = open(fn,'w')
    for (mir,gene) in combined_results:
        sequences = combined_results[(mir,gene)]
        sources = combined_results2[(mir,gene)]
        sequences = string.join(sequences,'|')
        sequences = string.replace(sequences,' ','|') ### Seems to occur with PicTar
        data.write(mir+'\t'+gene+'\t'+sequences+'\t'+sources+'\n')
    data.close()

def runProgram(Species,Force,Only_add_sequence_to_previous_results):
    global species; global only_add_sequence_to_previous_results; global symbol_ensembl; global force
    global ens_gene_to_transcript; global microRNA_target_db; global mir_hit_db; global parse_sequences
    global symbol_ensembl_current; combined_results={}
    species = Species; compare_to_user_data = 'no'; force = Force
    only_add_sequence_to_previous_results = Only_add_sequence_to_previous_results
    
    ens_gene_to_transcript = importEnsTranscriptAssociations({},'archive')
    ens_gene_to_transcript = importEnsTranscriptAssociations(ens_gene_to_transcript,'current')
    symbol_ensembl,symbol_ensembl_current = processEnsemblAnnotations()

    microRNA_target_db = {}; mir_hit_db = {}
    if only_add_sequence_to_previous_results != 'yes':
        parse_sequences = 'no'
        try: del symbol_ensembl['']
        except KeyError: null=[]
        
        try: sangerImport(parse_sequences)
        except Exception: null=[]        
        TargetScanImport(parse_sequences,'yes')
        mirandaImport(parse_sequences,'yes')

        added = pictarImport(parse_sequences,'pre-computed',{})
        added = pictarImport(parse_sequences,'symbol-based',added)
        if compare_to_user_data == 'yes': importExpressionData();mirHitImport()
        findMirTargetOverlaps(); exportCombinedMirResultSequences()
    else: exportCombinedMirResultSequences()

def verifyFile(filename,species):
    fn=filepath(filename)
    try:
        for line in open(fn,'rU').xreadlines():break
    except Exception:
        if species in filename: server_folder = species ### Folder equals species unless it is a universal file
        else: server_folder = 'all'
        import update; reload(update)
        update.downloadCurrentVersion(filename,server_folder,'txt')
    
if __name__ == '__main__':
    a = 'Mm'; b = 'Hs'; c = 'Rn'
    species = b; force = 'no'
    runProgram(species,force,'no'); sys.exit()