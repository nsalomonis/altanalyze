import sys, string
import os.path
import unique
dirfile = unique

############ File Import Functions #############

def read_directory(sub_dir):
    dir=os.path.dirname(dirfile.__file__)
    dir_list = os.listdir(dir + sub_dir)
    #add in code to prevent folder names from being included
    dir_list2 = [] 
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".all" or entry[-5:] == ".data" or entry[-3:] == ".fa":
            dir_list2.append(entry)
    return dir_list2

def filepath(filename):
    dir=os.path.dirname(dirfile.__file__)       #directory file is input as a variable under the main            
    fn=os.path.join(dir,filename)
    return fn

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

def importEnsemblAnnotations():
    global redundant_ensembl_by_build; redundant_ensembl_by_build={}
    filename = 'input/'+species+'/'+species+'_Ensembl-annotations_simple.txt'
    fn=filepath(filename); symbol_ensembl={}
    print 'importing', filename
    for line in open(fn,'r').xreadlines():        
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t'); ensembl=t[0]
        try: symbol = t[2]
        except IndexError: symbol = ''
        if len(symbol)>1 and len(ensembl)>1:
            try: symbol_ensembl[symbol].append(ensembl)
            except KeyError: symbol_ensembl[symbol] = [ensembl]
    for symbol in symbol_ensembl:
      if len(symbol_ensembl[symbol])>1: ###Thus there are more than one Ensembl gene IDs that correspond... probably due to versioning (which we've seen)
        for ensembl in symbol_ensembl[symbol]:
            for ensembl2 in symbol_ensembl[symbol]:
                if ensembl != ensembl2:
                    try: redundant_ensembl_by_build[ensembl].append(ensembl2)
                    except KeyError: redundant_ensembl_by_build[ensembl] = [ensembl2]
    print 'len(symbol_ensembl)',len(symbol_ensembl)
    return symbol_ensembl

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
        filename = 'input/'+species+'/'+'hits.txt'
        print 'parsing', filename
        fn=filepath(filename); x=0
        for line in open(fn,'r').xreadlines():         
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0: x=1
            else:
                mir = t[0]; fold = t[1]
                mir_hit_db[mir] = fold
    except IOError: mir_hit_db={}

def importExpressionData():
    global upregulated; global downregulated
    try:
        filename = 'input/'+species+'/'+'expression-data.txt'
        print 'parsing', filename
        fn=filepath(filename); x=0; upregulated={}; downregulated={}
        for line in open(fn,'r').xreadlines():         
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0: x=1
            else:
                gene = t[0]; log_fold_change = float(t[-2]); ttest = float(t[-1])
                if ttest<0.05 and log_fold_change>1: upregulated[gene] = log_fold_change,ttest
                if ttest<0.05 and log_fold_change<-1: downregulated[gene] = log_fold_change,ttest
    except IOError: upregulated={}; downregulated={}
            
def pictarImport(parse_sequences):
    """Annotations originally from the file: ng1536-S3.xls, posted as supplementary data at:
    http://www.nature.com/ng/journal/v37/n5/suppinfo/ng1536_S1.html. The file being parsed here has been pre-matched to Ensembl IDs
    using the ExonModule of LinkEST, for human."""
    mir_sequences=[]
    if species == 'Mm': filename = 'input/'+species+'/'+'pictar-target-annotated.txt'; tax = '10090'; prefix = 'hsa-'
    if species == 'Hs': filename = 'input/'+species+'/'+'pictar-conserved-targets-2005.txt'; tax = '9606'; prefix = 'mmu-'
    
    print 'parsing', filename; count=0
    fn=filepath(filename); x=1; added={}
    for line in open(fn,'r').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            if species == 'Hs': ensembl_geneid, mir, mir_sequences = t; ensembl_geneids = [ensembl_geneid]
            if species == 'Mm':
                mm_symbol=t[3];mir=t[6];mir_sequences=t[11]; mir = string.replace(mir,'hsa','mmu')
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

def sangerImport(parse_sequences):
    """"Sanger center (miRBase) sequence was provided as a custom (requested) dump of their v5 target predictions
    (http://microrna.sanger.ac.uk/targets/v5/), containing Ensembl gene IDs, microRNA names, and putative target
    sequences, specific for either mouse or human. Mouse was requested in late 2005 whereas human in late 2007.
    These same annotation files, missing the actual target sequence but containing an ENS transcript and coordinate
    locations for that build (allowing seqeunce extraction with the appropriate Ensembl build) exist at:
    http://microrna.sanger.ac.uk/cgi-bin/targets/v5/download.pl"""
    
    if species == 'Hs': filename = 'input/'+species+'/'+'mirbase-v5_homo_sapiens.mirna.txt'; prefix = 'hsa-'
    if species == 'Mm': filename = 'input/'+species+'/'+'sanger_miR_target_predictions.txt'; prefix = 'mmu-'    

    print 'parsing', filename; count=0
    fn=filepath(filename); x=1; mir_sequences=[]
    
    for line in open(fn,'r').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            if species == 'Hs':
                try: mir = t[1]; ensembl_geneid = t[17]; mir_sequences = string.upper(t[14])
                except IndexError: print line;kill
            elif species == 'Mm':
                ens_transcript,mir,mir_sequences = t
                if ens_transcript in ens_gene_to_transcript: ensembl_geneid = ens_gene_to_transcript[ens_transcript][0]
                else: ensembl_geneid=''
            if ensembl_geneid in redundant_ensembl_by_build: ###Thus there are redundant geneids
                geneid_ls = redundant_ensembl_by_build[ensembl_geneid]+[ensembl_geneid]
            else: geneid_ls = [ensembl_geneid]
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

def TargetScanImport():
    """The TargetScan data is currently extracted from a cross-species conserved family file. This file only contains
    gene symbol, microRNA name and 3'UTR seed locations."""
    if species == 'Mm': tax = '10090'; prefix = 'mmu-'
    if species == 'Hs': tax = '9606'; prefix = 'hsa-'
    global l
    mir_sequences = []; count=0
    uniprot_symbol_ensembl = importUniProtAnnotations()
    for symbol in uniprot_symbol_ensembl:
        ensembls = uniprot_symbol_ensembl[symbol]
        if symbol not in symbol_ensembl:
            symbol_ensembl[symbol] = ensembls
    if species != 'Hs':
        symbol_ensembl2={}
        for symbol in symbol_ensembl:
            ensembls = symbol_ensembl[symbol]
            symbol = string.upper(symbol)
            symbol_ensembl2[symbol] = ensembls
    else: symbol_ensembl2 = symbol_ensembl

    """
    ### Cross-species TargetScan file with UTR seqeunces for all genes with reported targets in the conserved family file
    ### Although this file includes valid sequence data that appears to match up to the target file, the target file
    ### appears to only list the seed seqeunce location (UTR start and stop) and not the full binding sequence and thus
    ### is not ammenable to probe set alignment.
    filename = 'input/TargetScan-UTR_Sequences.txt'
    print 'parsing', filename
    fn=filepath(filename); x=0; target_scan_gene_utr_seq={}
    for line in open(fn,'r').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            symbol = t[1]; tax_id = t[3]; utr_seq = t[4]
            if tax_id == tax:
                utr_seq_no_gaps = string.replace(utr_seq,'-','')
                target_scan_gene_utr_seq[symbol] = utr_seq_no_gaps
                if symbol == 'EGR1':
                    print utr_seq; print utr_seq[410:420];kill"""
        
    filename = 'input/target-scan-Conserved_Family_Conserved_Targets_Info.txt'
    print 'parsing', filename
    fn=filepath(filename); x=0; k=[]; l=[]
    for line in open(fn,'r').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            mir = t[0]; gene_symbol = t[1]; taxid = t[2]; utr_start = t[3]; utr_start = t[3]
            if '/' in mir:
                mir_list=[]
                mirs = string.split(mir,'/')
                for mirid in mirs[1:]:
                    mirid = 'miR-'+mirid
                    mir_list.append(mirid)
                mir_list.append(mirs[0])
            else: mir_list = [mir]
            if taxid == tax: ###human
                target_scan_gene_utr_seq[symbol] = utr_seq_no_gaps
                if gene_symbol in symbol_ensembl2: ensembl = symbol_ensembl2[gene_symbol]; proceed = 'yes'; k.append(gene_symbol)
                else: proceed = 'no'; l.append(gene_symbol)
                ###Already multiple geneids associated with each symbol so don't need to worry about renundancy
                if proceed == 'yes':
                    for geneid in ensembl:
                        for mir in mir_list:
                            #if geneid == 'ENSMUSG00000029467': print mir
                            y = MicroRNATargetData(geneid,gene_symbol,mir_sequences,prefix+mir,'TargetScan')
                            count+=1
                            try: microRNA_target_db[prefix+mir].append(y)
                            except KeyError: microRNA_target_db[prefix+mir] = [y]
    k = unique.unique(k); l = unique.unique(l)
    print 'ensembls-found:',len(k),', not found:',len(l)
    print count, 'miRNA-target relationships added for TargetScan'
    
def mirandaImport(parse_sequences):
    """Miranda data is avaialble from two file types from different websites. The first is human-centric with multi-species
    target alignment information organized by Ensembl gene ID (http://cbio.mskcc.org/research/sander/data/miRNA2003/mammalian/index.html).
    A larger set of associations was also pulled from species specific files (http://www.microrna.org/microrna/getDownloads.do),
    where gene symbol was related to Ensembl gene. Both files provided target microRNA sequence."""
    
    filename = 'input/miranda-by_gene140.txt'
    print 'parsing', filename; count=0; bad_line=0
    fn=filepath(filename); x=1 ###since the first line contains data
    if species == 'Hs': mir_type = 'hsa'; ens_transcript_prefix = 'ENST'; trans_index = 0
    if species == 'Mm': mir_type = 'mmu'; ens_transcript_prefix = 'ENSMUST'; trans_index = 2
    for line in open(fn,'r').xreadlines():         
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
        filename = 'input/'+species+'/'+'mouse_predictions_4dwnload.txt'
        print 'parsing', filename; count2=0; added={}; added2={}
        fn=filepath(filename); x=1 ###since the first line contains data
        if species == 'Hs': mir_type = 'hsa'; ens_transcript_prefix = 'ENST'; trans_index = 0
        if species == 'Mm': mir_type = 'mmu'; ens_transcript_prefix = 'ENSMUST'; trans_index = 2
        for line in open(fn,'r').xreadlines():         
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0: x=1
            else:
                symbol = t[2]; mir = t[4]; mir_sequence = t[7]; mir_sequence = string.replace(string.upper(mir_sequence),'U','T')
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
    
def importEnsTranscriptAssociations():
    ###This function is used to extract out EnsExon to EnsTranscript relationships to find out directly
    ###which probesets associate with which transcripts and then which proteins
    filename = 'input/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    print 'parsing', filename
    fn=filepath(filename); global ens_gene_to_transcript;  ens_gene_to_transcript={}
    for line in open(fn,'r').readlines():
        data,null = string.split(line,'\n')
        t = string.split(data,'\t')
        ens_geneid,chr,strand,exon_start,exon_end,ens_exonid,constitutive_exon,ens_transcript_id = t
        try: ens_gene_to_transcript[ens_transcript_id].append(ens_geneid)
        except KeyError:ens_gene_to_transcript[ens_transcript_id] = [ens_geneid]
    ens_gene_to_transcript = eliminateRedundant(ens_gene_to_transcript)    
    return ens_gene_to_transcript

def importArchivedEnsTranscriptAssociations(ens_gene_to_transcript,symbol_ensembl):
    try:
        ###This function is used to extract out EnsExon to EnsTranscript relationships to find out directly
        ###which probesets associate with which transcripts and then which proteins
        filename = 'input/'+species+'/'+species+'_Ensembl_transcript.txt'
        print 'parsing', filename
        fn=filepath(filename)
        for line in open(fn,'r').readlines():
            data,null = string.split(line,'\n');t = string.split(data,'\t')
            ens_geneid, ens_transcript_id, symbol, null, uniprot, entrezgene, ref_prot  = t
            #if ' ' not in ens_transcript_id: print ens_transcript_id, ens_geneid, symbol;kill
            try: ens_gene_to_transcript[ens_transcript_id].append(ens_geneid)
            except KeyError:ens_gene_to_transcript[ens_transcript_id] = [ens_geneid]
            try: symbol_ensembl[symbol].append(ens_geneid)
            except KeyError:symbol_ensembl[symbol] = [ens_geneid]
        ens_gene_to_transcript = eliminateRedundant(ens_gene_to_transcript)
        symbol_ensembl = eliminateRedundant(symbol_ensembl)
    except IOError: null=[]
    return ens_gene_to_transcript,symbol_ensembl

def findMirTargetOverlaps():
    print 'Number of microRNAs to be analyzed:',len(mir_hit_db)
    print "Number of microRNAs which have predicted targets:", len(microRNA_target_db)
    combined_output = 'output/'+species+'/'+'combined_gene-targets.txt'
    fn1=filepath(combined_output);data1 = open(fn1,'w')
    mir_gene_mult_hits={}
    for mir in microRNA_target_db:
        if mir in mir_hit_db:
            output_file = 'output/'+species+'/'+mir+'_gene-targets.txt'; output_file = string.replace(output_file,'*','-')
            fn=filepath(output_file);data = open(fn,'w'); delete = {}
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
        output_file = 'output/regulated/'+mir+'_regualted-gene-targets.txt'; output_file = string.replace(output_file,'*','-')
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
    g = GrabFiles(); g.setdirectory('/input/'+species); filenames = g.searchdirectory('UniProt')
    uniprot_symbol_ensembl={}
    fn=filepath(filenames[0]); x=1
    for line in open(fn,'r').readlines():             
        data,null = string.split(line,'\n')  #remove endline
        t = string.split(data,'\t')  #remove endline
        symbols = t[2]; ###ASPN; ORFNames=RP11-77D6.3-002, hCG_1784540
        symbols = string.split(symbols,';')
        primary_symbol = symbols[0]
        if len(symbols)>1:
            null,symbols = string.split(symbols[1],'=')
            symbols = string.split(symbols,', ')
            symbols = symbols+[primary_symbol]
        else: symbols = [primary_symbol]
        ensembl_ids = t[3]; ensembl_ids = string.split(ensembl_ids,",")
        for symbol in symbols:
            if len(symbol)>0:
                if len(ensembl_ids)>0: uniprot_symbol_ensembl[symbol] = ensembl_ids #; print symbol,ensembl_ids;kil
    return uniprot_symbol_ensembl

def exportCombinedMirResultSequences():
    combined_input = 'output/'+species+'/'+'combined_gene-targets.txt'
    global combined_results
    fn=filepath(combined_input); combined_results={}; combined_results2={} ###the second is for storing source information
    for line in open(fn,'r').xreadlines():
        data,null = string.split(line,'\n')  #remove endline
        t = string.split(data,'\t')  #remove endline
        mir, gene, sources = t
        combined_results[mir,gene] = []
        combined_results2[mir,gene] = sources
    microRNA_target_db = {}; mir_hit_db = {}
    mirandaImport(parse_sequences)
    sangerImport(parse_sequences)
    pictarImport(parse_sequences)
    
    output_file = 'output/'+species+'/'+'combined_gene-target-sequences.txt'
    combined_results = eliminateRedundant(combined_results)
    fn=filepath(output_file);data = open(fn,'w')
    for (mir,gene) in combined_results:
        sequences = combined_results[(mir,gene)]
        sources = combined_results2[(mir,gene)]
        sequences = string.join(sequences,'|')
        data.write(mir+'\t'+gene+'\t'+sequences+'\t'+sources+'\n')
    data.close()

if __name__ == '__main__':
    a = 'Mm'; b = 'Hs'
    species = b
    only_add_sequence_to_previous_results = 'no'
    dirfile = unique

    symbol_ensembl = importEnsemblAnnotations()
    ens_gene_to_transcript = importEnsTranscriptAssociations()
    ens_gene_to_transcript,symbol_ensembl = importArchivedEnsTranscriptAssociations(ens_gene_to_transcript,symbol_ensembl)    
    if only_add_sequence_to_previous_results != 'yes':
        parse_sequences = 'no'
        microRNA_target_db = {}; mir_hit_db = {}
        try: del symbol_ensembl['']
        except KeyError: null=[]
        
        TargetScanImport()
        mirandaImport(parse_sequences)
        sangerImport(parse_sequences)
        pictarImport(parse_sequences)

        ###Not required if you just want to export all associations
        importExpressionData()
        mirHitImport()
        
        findMirTargetOverlaps()
    else:
        symbol_ensembl = importEnsemblAnnotations()
        microRNA_target_db = {}; mir_hit_db = {}
        parse_sequences = 'yes'
        exportCombinedMirResultSequences()
    