###BAMtoJunctionBED
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

"""This script can be run on its own to extract a single BAM file at a time or
indirectly by multiBAMtoBED.py to extract junction.bed files (Tophat format)
from many BAM files in a single directory at once. Currently uses the Tophat
predicted Strand notation opt('XS') for each read. This can be substituted with
strand notations from other aligners (check with the software authors)."""

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import pysam
import copy,getopt
import time
import traceback
try: import export
except Exception: pass
try: import unique
except Exception: pass

try:
    import TabProxies
    import ctabix
    import csamtools
    import cvcf
except Exception:
    try:
        if os.name != 'posix': print traceback.format_exc()
    except Exception: pass

def getSpliceSites(cigarList,X):
    cummulative=0
    coordinates=[]
    for (code,seqlen) in cigarList:
        if code == 0:
            cummulative+=seqlen
        if code == 3:
            #if strand == '-':
            five_prime_ss = str(X+cummulative)
            cummulative+=seqlen  ### add the intron length
            three_prime_ss = str(X+cummulative+1) ### 3' exon start (prior exon splice-site + intron length)
            coordinates.append([five_prime_ss,three_prime_ss])
            up_to_intron_dist = cummulative
    return coordinates, up_to_intron_dist

def writeJunctionBedFile(junction_db,jid,o):
    strandStatus = True
    for (chr,jc,tophat_strand) in junction_db:
        if tophat_strand==None:
            strandStatus = False
        break
    
    if strandStatus== False: ### If no strand information in the bam file filter and add known strand data
        junction_db2={}
        for (chr,jc,tophat_strand) in junction_db:
            original_chr = chr
            if 'chr' not in chr:
                chr = 'chr'+chr
            for j in jc:
                try:
                    strand = splicesite_db[chr,j]
                    junction_db2[(original_chr,jc,strand)]=junction_db[(original_chr,jc,tophat_strand)]
                except Exception: pass
        junction_db = junction_db2
    for (chr,jc,tophat_strand) in junction_db:
        x_ls=[]; y_ls=[]; dist_ls=[]
        read_count = str(len(junction_db[(chr,jc,tophat_strand)]))
        for (X,Y,dist) in junction_db[(chr,jc,tophat_strand)]:
            x_ls.append(X); y_ls.append(Y); dist_ls.append(dist)
        outlier_start = min(x_ls); outlier_end = max(y_ls); dist = str(max(dist_ls))
        exon_lengths = outlier_start
        exon_lengths = str(int(jc[0])-outlier_start)+','+str(outlier_end-int(jc[1])+1)
        junction_id = 'JUNC'+str(jid)+':'+jc[0]+'-'+jc[1] ### store the unique junction coordinates in the name
        output_list = [chr,str(outlier_start),str(outlier_end),junction_id,read_count,tophat_strand,str(outlier_start),str(outlier_end),'255,0,0\t2',exon_lengths,'0,'+dist]
        o.write(string.join(output_list,'\t')+'\n')
          
def writeIsoformFile(isoform_junctions,o):
    for coord in isoform_junctions:
        isoform_junctions[coord] = unique.unique(isoform_junctions[coord])
        if '+' in coord:
            print coord, isoform_junctions[coord] 
    
    if '+' in coord:
        sys.exit()
    
def verifyFileLength(filename):
    count = 0
    try:
        fn=unique.filepath(filename)
        for line in open(fn,'rU').xreadlines():
            count+=1
            if count>9: break
    except Exception: null=[]
    return count

def retreiveAllKnownSpliceSites(returnExonRetention=False,DesignatedSpecies=None,path=None):
    ### Uses a priori strand information when none present
    import export, unique
    chromosomes_found={}
    try: parent_dir = export.findParentDir(bam_file)
    except Exception: parent_dir = export.findParentDir(path)
    species = None
    for file in os.listdir(parent_dir):
        if 'AltAnalyze_report' in file and '.log' in file:
            log_file = unique.filepath(parent_dir+'/'+file)
            log_contents = open(log_file, "rU")
            species_tag = '	species: '
            for line in log_contents:
                line = line.rstrip()
                if species_tag in line:
                    species = string.split(line,species_tag)[1]
    if species == None:
        try: species = IndicatedSpecies
        except Exception: species = DesignatedSpecies
    
    splicesite_db={}
    gene_coord_db={}
    try:
        if ExonReference==None:
            exon_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
            length = verifyFileLength(exon_dir)
    except Exception:
        #print traceback.format_exc();sys.exit()
        length = 0
    if length==0:
        exon_dir = ExonReference
    refExonCoordinateFile = unique.filepath(exon_dir)
    firstLine=True
    for line in open(refExonCoordinateFile,'rU').xreadlines():
        if firstLine: firstLine=False
        else:
            line = line.rstrip('\n')
            t = string.split(line,'\t'); #'gene', 'exon-id', 'chromosome', 'strand', 'exon-region-start(s)', 'exon-region-stop(s)', 'constitutive_call', 'ens_exon_ids', 'splice_events', 'splice_junctions'
            geneID, exon, chr, strand, start, stop = t[:6]
            spliceEvent = t[-2]
            #start = int(start); stop = int(stop)
            #geneID = string.split(exon,':')[0]
            try:
                gene_coord_db[geneID,chr].append(int(start))
                gene_coord_db[geneID,chr].append(int(stop))
            except Exception:
                gene_coord_db[geneID,chr] = [int(start)]
                gene_coord_db[geneID,chr].append(int(stop))
            if returnExonRetention:
                if 'exclusion' in spliceEvent or 'exclusion' in spliceEvent:
                    splicesite_db[geneID+':'+exon]=[]
            else:
                splicesite_db[chr,start]=strand
                splicesite_db[chr,stop]=strand
                if len(chr)<5 or ('GL0' not in chr and 'GL' not in chr and 'JH' not in chr and 'MG' not in chr):
                    chromosomes_found[string.replace(chr,'chr','')] = []
    for i in gene_coord_db:
        gene_coord_db[i].sort()
        gene_coord_db[i] = [gene_coord_db[i][0],gene_coord_db[i][-1]]
    return splicesite_db,chromosomes_found,gene_coord_db

def parseJunctionEntries(bam_dir,multi=False, Species=None, ReferenceDir=None):
    global bam_file
    global splicesite_db
    global IndicatedSpecies
    global ExonReference
    IndicatedSpecies = Species
    ExonReference = ReferenceDir
    bam_file = bam_dir
    splicesite_db={}; chromosomes_found={}

    start = time.time()
    try: import collections; junction_db=collections.OrderedDict()
    except Exception:
        try: import ordereddict; junction_db = ordereddict.OrderedDict()
        except Exception: junction_db={}
    original_junction_db = copy.deepcopy(junction_db)
    
    bamf = pysam.Samfile(bam_dir, "rb" )
    ### Is there are indexed .bai for the BAM? Check.
    try:
        for entry in bamf.fetch():
            codes = map(lambda x: x[0],entry.cigar)
            break
    except Exception:
        ### Make BAM Index
        if multi == False:
            print 'Building BAM index file for', bam_dir
        bam_dir = str(bam_dir)
        #On Windows, this indexing step will fail if the __init__ pysam file line 51 is not set to - catch_stdout = False
        pysam.index(bam_dir)
        bamf = pysam.Samfile(bam_dir, "rb" )

    chromosome = False
    barcode_pairs={}
    bam_reads=0
    count=0
    jid = 1
    prior_jc_start=0
    import Bio; from Bio.Seq import Seq
    l1 = None; l2=None
    o = open (string.replace(bam_dir,'.bam','.export2.txt'),"w")
    spacer='TGGT'
    for entry in bamf.fetch():
        #if entry.query_name == 'M03558:141:GW181002:1:2103:13361:6440':
        if spacer in entry.seq:
            if entry.seq.index(spacer) == 14:
                viral_barcode = entry.seq[:48]
                try:
                    mate = bamf.mate(entry)
                    mate_seq = Seq(mate.seq)
                    cell_barcode=str(mate_seq.reverse_complement())[:16]
                    if (viral_barcode,cell_barcode) not in barcode_pairs:
                        o.write(viral_barcode+'\t'+cell_barcode+'\n')
                    barcode_pairs[viral_barcode,cell_barcode]=[]
                    if 'ATAGCGGGAACATGTGGTCATGGTACTGACGTTGACACGTACGTCATA' == viral_barcode:
                        print entry.query_name, cell_barcode, mate_seq
                except:
                    pass
                #print viral_barcode, mate.seq;sys.exit()
                    

        count+=1
        #if count==100: sys.exit()
    bamf.close()
    o.close()

if __name__ == "__main__":
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a BAM file as input in the command-line"
        print "Example: python BAMtoJunctionBED.py --i /Users/me/sample1.bam"
        sys.exit()
    else:
        Species = None
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','species='])
        for opt, arg in options:
            if opt == '--i': bam_dir=arg ### full path of a BAM file
            elif opt == '--species': Species=arg ### species for STAR analysis to get strand
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
            
    try: parseJunctionEntries(bam_dir,Species=Species)
    except ZeroDivisionError:
        print [sys.argv[1:]],'error'; error

""" Benchmarking notes: On a 2017 MacBook Pro with 16GB of RAM and a local 7GB BAM file (solid drive), 9 minutes (526s) to complete writing a junction.bed.
To simply search through the file without looking at the CIGAR, the script takes close to 5 minutes (303s)"""
