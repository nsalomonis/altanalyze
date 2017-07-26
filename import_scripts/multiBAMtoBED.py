### hierarchical_clustering.py
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

""" Batch script for extracting many junction.bed and building exon.bed files from
an input set of BAM files in a directory. Requires a reference text file containing
exon regions (currently provided from AltAnalyze - see ReferenceExonCoordinates
folder). Can produce only junction.bed files, only a combined exon reference or only
exon.bed files optionally. Can run using a single processor or multiple simultaneous
processes (--m flag)."""

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import export
import time
import shutil
import unique
import subprocess
from import_scripts import BAMtoJunctionBED
from import_scripts import BAMtoExonBED
import getopt
import traceback
    
################# General data import methods #################

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def cleanUpLine(line):
    data = string.replace(line,'\n','')
    data = string.replace(data,'\c','')
    data = string.replace(data,'\r','')
    data = string.replace(data,'"','')
    return data
           
def getFolders(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Only get folder names
    for entry in dir_list:
        if '.' not in entry: dir_list2.append(entry)
    return dir_list2

def getFiles(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Only get folder names
    for entry in dir_list:
        if '.' in entry: dir_list2.append(entry)
    return dir_list2

def getFolders(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Only get folder names
    for entry in dir_list:
        if '.' not in entry: dir_list2.append(entry)
    return dir_list2

def parallelBAMProcessing(directory,refExonCoordinateFile,bed_reference_dir,analysisType=[],useMultiProcessing=False,MLP=None,root=None):
    paths_to_run=[]
    errors=[]
    if '.bam' in directory:
        ### Allow a single BAM file to be specifically analyzed (e.g., bsub operation)
        bam_file = directory
        bam_file = string.replace(directory,'\\','/')
        directory = string.join(string.split(directory,'/')[:-1],'/')
    else:
        bam_file = None
        
    outputExonCoordinateRefBEDfile = str(bed_reference_dir)
    bed_reference_dir = string.replace(bed_reference_dir,'\\','/')
    ### Check if the BAM files are located in the target folder (not in subdirectories)
    files = getFiles(directory)
    for file in files:
        if '.bam' in file and '.bai' not in file:
            source_file = directory+'/'+file
            source_file = filepath(source_file)
            output_filename = string.replace(file,'.bam','')
            output_filename = string.replace(output_filename,'=','_')
            destination_file = directory+'/'+output_filename+'__exon.bed'
            destination_file = filepath(destination_file)
            paths_to_run.append((source_file,refExonCoordinateFile,bed_reference_dir,destination_file))

    ### Otherwise, check subdirectories for BAM files
    folders = getFolders(directory)
    if len(paths_to_run)==0:
        for top_level in folders: 
            try:
                files = getFiles(directory+'/'+top_level)
                for file in files:
                    if '.bam' in file and '.bai' not in file:
                        source_file = directory+'/'+file
                        source_file = filepath(source_file)
                        destination_file = directory+'/'+top_level+'__exon.bed'
                        destination_file = filepath(destination_file)
                        paths_to_run.append((source_file,refExonCoordinateFile,bed_reference_dir,destination_file))
            except Exception: pass
    
    ### If a single BAM file is indicated
    if bam_file != None:
        output_filename = string.replace(bam_file,'.bam','')
        output_filename = string.replace(output_filename,'=','_')
        destination_file = output_filename+'__exon.bed'
        paths_to_run = [(bam_file,refExonCoordinateFile,bed_reference_dir,destination_file)]

    if 'reference' in analysisType and len(analysisType)==1:
        augmentExonReferences(directory,refExonCoordinateFile,outputExonCoordinateRefBEDfile)
        sys.exit()
        
    if useMultiProcessing:
        pool_size = MLP.cpu_count()
        if len(paths_to_run)<pool_size:
            pool_size = len(paths_to_run)
        print 'Using %d processes' % pool_size
        if len(paths_to_run) > pool_size:
            pool_size = len(paths_to_run)
        
        if len(analysisType) == 0 or 'junction' in analysisType:
            print 'Extracting junction alignments from BAM files...',
            pool = MLP.Pool(processes=pool_size)
            try: results = pool.map(runBAMtoJunctionBED, paths_to_run) ### worker jobs initiated in tandem
            except ValueError:
                print_out = '\WARNING!!! No Index found for the BAM files (.bam.bai). Sort and Index using Samtools prior to loading in AltAnalyze'
                print traceback.format_exc()
                if root!=None:
                    import UI
                    UI.WarningWindow(print_out,'Exit');sys.exit()
                    
            try:pool.close(); pool.join(); pool = None
            except Exception: pass
            print_out=None
            for sample,missing in results:
                if len(missing)>1:
                    print_out = '\nWarning!!! %s chromosomes not found in: %s (PySam platform-specific error)' % (string.join(missing,', '),sample)
                    
            if root!=None and print_out!=None:
                try:
                    import UI
                    UI.WarningWindow(print_out,'Continue')
                except Exception: pass  
                    
            print len(paths_to_run), 'BAM files','processed'
        if len(analysisType) == 0 or 'reference' in analysisType:
            #print 'Building exon reference coordinates from Ensembl/UCSC and all junctions...',
            augmentExonReferences(directory,refExonCoordinateFile,outputExonCoordinateRefBEDfile)
            #print 'completed'

        print 'Extracting exon alignments from BAM files...',
        if len(analysisType) == 0 or 'exon' in analysisType:
            pool = MLP.Pool(processes=pool_size)
            results = pool.map(runBAMtoExonBED, paths_to_run) ### worker jobs initiated in tandem
            try:pool.close(); pool.join(); pool = None
            except Exception: pass
            print len(paths_to_run), 'BAM files','processed'

    else:
        if len(analysisType) == 0 or 'junction' in analysisType:
            for i in paths_to_run:
                runBAMtoJunctionBED(i)
        if len(analysisType) == 0 or 'reference' in analysisType:
            augmentExonReferences(directory,refExonCoordinateFile,outputExonCoordinateRefBEDfile)
        if len(analysisType) == 0 or 'exon' in analysisType:
            for i in paths_to_run:
                runBAMtoExonBED(i)

def runBAMtoJunctionBED(paths_to_run):
    bamfile_dir,refExonCoordinateFile,bed_reference_dir,output_bedfile_path = paths_to_run
    output_bedfile_path = string.replace(bamfile_dir,'.bam','__junction.bed')
    #if os.path.exists(output_bedfile_path) == False: ### Only run if the file doesn't exist
    results = BAMtoJunctionBED.parseJunctionEntries(bamfile_dir,multi=True,ReferenceDir=refExonCoordinateFile)
    #else: print output_bedfile_path, 'already exists.'
    return results
    
def runBAMtoExonBED(paths_to_run):
    bamfile_dir,refExonCoordinateFile,bed_reference_dir,output_bedfile_path = paths_to_run
    if os.path.exists(output_bedfile_path) == False: ### Only run if the file doesn't exist
        BAMtoExonBED.parseExonReferences(bamfile_dir,bed_reference_dir,multi=True,intronRetentionOnly=False)
    else:
        print output_bedfile_path, 'already exists... re-writing'
        BAMtoExonBED.parseExonReferences(bamfile_dir,bed_reference_dir,multi=True,intronRetentionOnly=False)

def getChrFormat(directory):
    ### Determine if the chromosomes have 'chr' or nothing
    files = getFiles(directory)
    chr_status=True
    for file in files:
        firstLine=True
        if 'junction' in file and '.bed' in file:
            for line in open(directory+'/'+file,'rU').xreadlines():
                if firstLine: firstLine=False
                else:
                    t = string.split(line)
                    chr = t[0]
                    if 'chr' not in chr:
                        chr_status = False
                        break
            break
    return chr_status
            
def augmentExonReferences(directory,refExonCoordinateFile,outputExonCoordinateRefBEDfile):
    print 'Building reference bed file from all junction.bed files'
    splicesite_db={} ### reference splice-site database (we only want to add novel splice-sites to our reference)
    real_splicesites={}
    introns={}
    novel_db={}
    reference_toplevel = string.join(string.split(outputExonCoordinateRefBEDfile,'/')[:-1],'/')
    try: os.mkdir(reference_toplevel) ### If the bed folder doesn't exist
    except Exception: pass
    chr_status = getChrFormat(directory)

    o = open (outputExonCoordinateRefBEDfile,"w")
    #refExonCoordinateFile = '/Users/saljh8/Desktop/Code/AltAnalyze/AltDatabase/EnsMart72/ensembl/Mm/Mm_Ensembl_exon.txt'
    reference_rows=0
    if '.gtf' in refExonCoordinateFile: firstLine = False
    else: firstLine = True
    for line in open(refExonCoordinateFile,'rU').xreadlines():
        if firstLine: firstLine=False
        else:
            line = line.rstrip('\n')
            reference_rows+=1
            t = string.split(line,'\t'); #'gene', 'exon-id', 'chromosome', 'strand', 'exon-region-start(s)', 'exon-region-stop(s)', 'constitutive_call', 'ens_exon_ids', 'splice_events', 'splice_junctions'
            geneID, exon, chr, strand, start, stop = t[:6]
            if chr_status == False:
                chr = string.replace(chr,'chr','')
            o.write(string.join([chr,start,stop,geneID+':'+exon,'',strand],'\t')+'\n')
            start = int(start); stop = int(stop)
            #geneID = string.split(exon,':')[0]
            splicesite_db[chr,start]=geneID
            splicesite_db[chr,stop]=geneID
            if 'I' in exon:
                try: introns[geneID].append([start,stop])
                except Exception: introns[geneID] = [[start,stop]]  
            
    files = getFiles(directory)
    for file in files:
        firstLine=True
        if 'junction' in file and '.bed' in file:
            for line in open(directory+'/'+file,'rU').xreadlines():
                if firstLine: firstLine=False
                else:
                    line = line.rstrip('\n')
                    t = string.split(line,'\t'); #'12', '6998470', '6998522', 'ENSG00000111671:E1.1_ENSE00001754003', '0', '-'  
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
                    seq_length = abs(float(exon1_stop-exon2_start)) ### Junction distance
                    key = chr,exon1_stop,exon2_start
                    if (chr,exon1_stop) not in splicesite_db: ### record the splice site and position of the max read
                        if (chr,exon2_start) in splicesite_db: ### only include splice sites where one site is known
                            geneID = splicesite_db[(chr,exon2_start)]
                            novel_db[chr,exon1_stop,strand] = exon1_start,geneID,5
                        real_splicesites[chr,exon2_start]=None
                    elif (chr,exon2_start) not in splicesite_db: ### record the splice site and position of the max read
                        if (chr,exon1_stop) in splicesite_db: ### only include splice sites where one site is known
                            #if 121652702 ==exon2_start:
                            #print chr, exon1_start,exon1_stop,exon2_start,exon2_stop, strand;sys.exit()
                            geneID = splicesite_db[(chr,exon1_stop)]
                            novel_db[chr,exon2_start,strand] = exon2_stop,geneID,3
                        real_splicesites[chr,exon1_stop]=None
                    else:
                        real_splicesites[chr,exon1_stop]=None
                        real_splicesites[chr,exon2_start]=None
                    
    print len(novel_db), 'novel splice sites and', len(real_splicesites), 'known splice sites.'
    
    gene_organized={}
    for (chr,pos1,strand) in novel_db:
        pos2,geneID,type = novel_db[(chr,pos1,strand)]
        try: gene_organized[chr,geneID,strand].append([pos1,pos2,type])
        except Exception: gene_organized[chr,geneID,strand] = [[pos1,pos2,type]]
    
    def intronCheck(geneID,coords):
        ### see if the coordinates are within a given intron
        try:
            for ic in introns[geneID]:
                if withinQuery(ic,coords):
                    return True
        except Exception:
            pass
        
    def withinQuery(ls1,ls2):
        imax = max(ls1)
        imin = min(ls1)
        qmax = max(ls2)
        qmin = min(ls2)
        if qmin >= imin and qmax <= imax:
            return True
        else:
            return False
    
    ### Compare the novel splice site locations in each gene
    added=[]
    for (chr,geneID,strand) in gene_organized:
        gene_organized[(chr,geneID,strand)].sort()
        if strand == '-':
            gene_organized[(chr,geneID,strand)].reverse()
        i=0
        set = gene_organized[(chr,geneID,strand)]
        for (pos1,pos2,type) in set:
            k = [pos1,pos2]
            annotation='novel'
            if i==0 and type == 3:
                if len(set)>1:
                    if set[i+1][-1]==5:
                        l = [set[i+1][0],pos1]
                        if (max(l)-min(l))<300 and intronCheck(geneID,l):
                            k=l
                            #print chr,k
                            annotation='novel-paired'
            elif type == 5:
                if set[i-1][-1]==3:
                    l = [set[i-1][0],pos1]
                    if (max(l)-min(l))<300 and intronCheck(geneID,l):
                        k=l
                        #print chr,k
                        annotation='novel-paired'
            k.sort(); i+=1
            if k not in added:
                values = string.join([chr,str(k[0]),str(k[1]),geneID+':'+annotation,'',strand],'\t')+'\n'
                added.append(k)
                o.write(values)
    o.close()

if __name__ == '__main__':
    import multiprocessing as mlp
    refExonCoordinateFile = ''
    outputExonCoordinateRefBEDfile = ''
    #bam_dir = "H9.102.2.6.bam"
    #outputExonCoordinateRefBEDfile = 'H9.102.2.6__exon.bed'
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a directory containing BAM files as input in the command-line"
        print "Example: python multiBAMtoBED.py --i /Users/me/BAMfiles --g /Users/me/ReferenceExonCoordinates/Hs_Ensembl_exon_hg19.txt --r /Users/me/ExonBEDRef/Hs_Ensembl_exon-cancer_hg19.bed --a exon --a junction --a reference"
        print "Example: python multiBAMtoBED.py --i /Users/me/BAMfiles --a junction"
        sys.exit()
    else:
        analysisType = []
        useMultiProcessing=False
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','g=','r=','a=','m='])
        for opt, arg in options:
            if opt == '--i': bam_dir=arg
            elif opt == '--g': refExonCoordinateFile=arg
            elif opt == '--r': outputExonCoordinateRefBEDfile=arg
            elif opt == '--a': analysisType.append(arg) ### options are: all, junction, exon, reference
            elif opt == '--m': ### Run each BAM file on a different processor
                if arg == 'yes': useMultiProcessing=True
                elif arg == 'True': useMultiProcessing=True
                else: useMultiProcessing=False
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
        if len(analysisType) == 0 or 'all' in analysisType:
            analysisType = ['exon','junction','reference']
            try:
                refExonCoordinateFile = refExonCoordinateFile
                outputExonCoordinateRefBEDfile = outputExonCoordinateRefBEDfile
            except Exception:
                print 'Please provide a exon coordinate text file using the option --g and a output coordinate file path (--r) to generate exon.bed files'
                analysisType = ['junction']
                refExonCoordinateFile = ''
                outputExonCoordinateRefBEDfile = ''

    try: bam_dir = bam_dir
    except Exception: print 'You must specify a directory of BAM files or a single bam file with --i';sys.exit()
    try: refExonCoordinateFile = refExonCoordinateFile
    except Exception: print 'You must specify a AltAnalyze exon coordinate text file with --g';sys.exit()
    try: outputExonCoordinateRefBEDfile = outputExonCoordinateRefBEDfile
    except Exception: print 'You must specify an output path for the exon.bed reference file location with --r (e.g., --r /users/Hs_exon.bed)';sys.exit()
    
    parallelBAMProcessing(bam_dir,refExonCoordinateFile,outputExonCoordinateRefBEDfile,analysisType=analysisType,useMultiProcessing=useMultiProcessing,MLP=mlp)
    