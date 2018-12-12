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
try:
    from pysam import libctabixproxies
except:
    print traceback.format_exc()
    pass

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

def exportIndexes(input_dir):
    import unique
    bam_dirs = unique.read_directory(input_dir)
    print 'Building BAM index files',
    for file in bam_dirs:
        if string.lower(file[-4:]) == '.bam':
            bam_dir = input_dir+'/'+file
            bamf = pysam.Samfile(bam_dir, "rb" )
            ### Is there an indexed .bai for the BAM? Check.
            try:
                for entry in bamf.fetch():
                    codes = map(lambda x: x[0],entry.cigar)
                    break
            except Exception:
                ### Make BAM Indexv lciv9df8scivx 
                print '.',
                bam_dir = str(bam_dir)
                #On Windows, this indexing step will fail if the __init__ pysam file line 51 is not set to - catch_stdout = False
                pysam.index(bam_dir)
                bamf = pysam.Samfile(bam_dir, "rb" )        

def parseJunctionEntries(bam_dir,multi=False, Species=None, ReferenceDir=None):
    global bam_file
    global splicesite_db
    global IndicatedSpecies
    global ExonReference
    IndicatedSpecies = Species
    ExonReference = ReferenceDir
    bam_file = bam_dir
    try: splicesite_db,chromosomes_found, gene_coord_db = retreiveAllKnownSpliceSites()
    except Exception:
        print traceback.format_exc()
        splicesite_db={}; chromosomes_found={}

    start = time.time()
    try: import collections; junction_db=collections.OrderedDict()
    except Exception:
        try: import ordereddict; junction_db = ordereddict.OrderedDict()
        except Exception: junction_db={}
    original_junction_db = copy.deepcopy(junction_db)
    
    bamf = pysam.Samfile(bam_dir, "rb" )
    ### Is there an indexed .bai for the BAM? Check.
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
    chromosomes={}
    bam_reads=0
    count=0
    jid = 1
    prior_jc_start=0
    l1 = None; l2=None
    o = open (string.replace(bam_dir,'.bam','__junction.bed'),"w")
    o.write('track name=junctions description="TopHat junctions"\n')
    export_isoform_models = False
    if export_isoform_models:
        io = open (string.replace(bam_dir,'.bam','__isoforms.txt'),"w")
        isoform_junctions = copy.deepcopy(junction_db)
    outlier_start = 0; outlier_end = 0; read_count = 0; c=0
    for entry in bamf.fetch():
      bam_reads+=1
      try: cigarstring = entry.cigarstring
      except Exception:
          codes = map(lambda x: x[0],entry.cigar)
          if 3 in codes: cigarstring = 'N'
          else: cigarstring = None
    
      if cigarstring != None:
        if 'N' in cigarstring: ### Hence a junction
            if prior_jc_start == 0: pass
            elif (entry.pos-prior_jc_start) > 5000 or bamf.getrname( entry.rname ) != chromosome: ### New chr or far from prior reads
                writeJunctionBedFile(junction_db,jid,o)
                #writeIsoformFile(isoform_junctions,io)
                junction_db = copy.deepcopy(original_junction_db) ### Re-set this object
                jid+=1
            
            chromosome = bamf.getrname( entry.rname )
            chromosomes[chromosome]=[] ### keep track
            X=entry.pos
            #if entry.query_name == 'SRR791044.33673569':
            #print chromosome, entry.pos, entry.reference_length, entry.alen, entry.query_name
            Y=entry.pos+entry.alen
            prior_jc_start = X

            try: tophat_strand = entry.opt('XS') ### TopHat knows which sequences are likely real splice sites so it assigns a real strand to the read
            except Exception:
                #if multi == False:  print 'No TopHat strand information';sys.exit()
                tophat_strand = None
            coordinates,up_to_intron_dist = getSpliceSites(entry.cigar,X)
            #if count > 100: sys.exit()
            #print entry.query_name,X, Y, entry.cigarstring, entry.cigar, tophat_strand
            for (five_prime_ss,three_prime_ss) in coordinates:
                jc = five_prime_ss,three_prime_ss
                #print X, Y, jc, entry.cigarstring, entry.cigar
                try: junction_db[chromosome,jc,tophat_strand].append([X,Y,up_to_intron_dist])
                except Exception: junction_db[chromosome,jc,tophat_strand] = [[X,Y,up_to_intron_dist]]
                
            if export_isoform_models:
                try:
                    mate = bamf.mate(entry) #https://groups.google.com/forum/#!topic/pysam-user-group/9HM6nx_f2CI
    
                    if 'N' in mate.cigarstring:
                        mate_coordinates,mate_up_to_intron_dist = getSpliceSites(mate.cigar,mate.pos)
                    else: mate_coordinates=[]
                except Exception: mate_coordinates=[]
                #print coordinates,mate_coordinates
                junctions = map(lambda x: tuple(x),coordinates)
                if len(mate_coordinates)>0:
                    try:
                        isoform_junctions[chromosome,tuple(junctions),tophat_strand].append(mate_coordinates)
                    except Exception:
                        isoform_junctions[chromosome,tuple(junctions),tophat_strand] = [mate_coordinates]
                else:
                    if (chromosome,tuple(junctions),tophat_strand) not in isoform_junctions:
                        isoform_junctions[chromosome,tuple(junctions),tophat_strand] = []
                
            count+=1
    writeJunctionBedFile(junction_db,jid,o) ### One last read-out
    if multi == False:
        print bam_reads, count, time.time()-start, 'seconds required to parse the BAM file'
    o.close()
    bamf.close()
    
    missing_chromosomes=[]
    for chr in chromosomes_found:
        if chr not in chromosomes:
            chr = string.replace(chr,'chr','')
            if chr not in chromosomes_found:
                if chr != 'M' and chr != 'MT':
                    missing_chromosomes.append(chr)
    #missing_chromosomes = ['A','B','C','D']
    try: bam_file = export.findFilename(bam_file)
    except Exception: pass
    return bam_file, missing_chromosomes
    

if __name__ == "__main__":
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a BAM file as input in the command-line"
        print "Example: python BAMtoJunctionBED.py --i /Users/me/sample1.bam"
        sys.exit()
    else:
        Species = None
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','species=','r='])
        for opt, arg in options:
            if opt == '--i': bam_dir=arg ### full path of a BAM file
            elif opt == '--species': Species=arg ### species for STAR analysis to get strand
            elif opt == '--r': reference_dir=arg ### An exon.bed reference file (created by AltAnalyze from junctions, multiBAMtoBED.py or other) - required for STAR to get strand if XS field is empty
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
            
    try: parseJunctionEntries(bam_dir,Species=Species,ReferenceDir=reference_dir)
    except ZeroDivisionError:
        print [sys.argv[1:]],'error'; error

""" Benchmarking notes: On a 2017 MacBook Pro with 16GB of RAM and a local 7GB BAM file (solid drive), 9 minutes (526s) to complete writing a junction.bed.
To simply search through the file without looking at the CIGAR, the script takes close to 5 minutes (303s)"""
