###BAMtoExonBED
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
indirectly by multiBAMtoBED.py to extract exon.bed files (Tophat format)
from many BAM files in a single directory at once. Requires an exon.bed reference
file for exon coordinates (genomic bins for which to sum unique read counts).
Excludes junction reads within each interval"""

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import pysam
import copy
import time
import getopt

def findGeneVariants(species,symbols,bam_dir,variants=None):
    global insertion_db
    insertion_db={}
    print symbols
    print bam_dir
    if len(symbols)>0:
        ### Search for genes not for coordinates
        search_locations = geneCoordinates(species,symbols)
    else:
        ### Search for coordinates and not genes
        search_locations = variantCoordinates(variants)
        
    ### Discover the variants
    variant_db = findVariants(bam_dir,search_locations)
    
    variant_filtered_db={}
    for var in variant_db:
        #print var, variant_db[var]
        if variant_db[var]>3:
            #print var,variant_db[var]
            variant_filtered_db[var] = variant_db[var]
    
    ### Quantify the variants versus background
    pileupAnalysis(bam_dir,variant_filtered_db)

def variantCoordinates(variants):
    search_locations=[]
    contents = open(variants, "rU")
    for line in contents:
        line = line.rstrip()
        chr,start,end,symbol = string.split(line,'\t')
        if 'chr' not in chr: chr = 'chr'+chr
        strand = 'NA'
        search_locations.append([chr,strand,start,end,symbol])
    return search_locations

def geneCoordinates(species,symbols):
    genes=[]
    from build_scripts import EnsemblImport
    ensembl_annotation_db = EnsemblImport.reimportEnsemblAnnotations(species,symbolKey=True)
    for symbol in symbols:
        if symbol in ensembl_annotation_db:
            ens_geneid = ensembl_annotation_db[symbol]
            genes.append((ens_geneid,symbol))
        else:
            print symbol, 'not found'
    
    ### Get gene genomic locations
    gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,'RNASeq','key_by_array')
    search_locations=[]
    for (gene,symbol) in genes:
        chr,strand,start,end = gene_location_db[gene]
        #if symbol == 'SRSF10': chr = 'chr1'; strand = '-'; start = '24295573'; end = '24306953'
        if len(chr)>6: print symbol, 'bad chromosomal reference:',chr
        else:
            search_locations.append([chr,strand,start,end,symbol])
        
    return search_locations
    
def findVariants(bam_dir,search_locations,multi=False):
    start_time = time.time()    
    bamfile = pysam.Samfile(bam_dir, "rb" )
    output_bed_rows=0
    #https://www.biostars.org/p/76119/
    variant_db={}
    reference_rows=0
    o = open (string.replace(bam_dir,'.bam','__variant.txt'),"w")
    for (chr,strand,start,stop,symbol) in search_locations: ### read each line one-at-a-time rather than loading all in memory
            read_count=0
            reference_rows+=1
            stop=int(stop)+100 ### buffer for single variants
            start=int(start)-100 ### buffer for single variants
            for alignedread in bamfile.fetch(chr, int(start),int(stop)):
                md = alignedread.opt('MD')
                omd = md
                codes = map(lambda x: x[0],alignedread.cigar)
                cigarstring = alignedread.cigarstring
                #print symbol,cigarstring
                if 1 in codes and alignedread.pos:
                    ### Thus an insertion is present
                    cigarstring = alignedread.cigarstring
                    chr = bamfile.getrname(alignedread.rname)
                    pos = alignedread.pos
                    def getInsertions(cigarList,X):
                        cummulative=0
                        coordinates=[]
                        for (code,seqlen) in cigarList:
                            if code == 0 or code == 3:
                                cummulative+=seqlen
                            if code == 1:
                                coordinates.append(X+cummulative)
                        return coordinates
                    coordinates = getInsertions(alignedread.cigar,pos)
                    """
                    print pos
                    print coordinates
                    print alignedread.seq
                    print codes
                    print alignedread.cigar
                    print cigarstring
                    
                    print md;sys.exit()
                    """
                    for pos in coordinates:
                        try: variant_db[chr,pos,symbol]+=1
                        except Exception: variant_db[chr,pos,symbol] = 1
                        insertion_db[chr,pos]=[]
                    continue
                    
                try:
                    int(md) ### If an integer, no mismatches or deletions present
                    continue
                except Exception:
                    #print chr, int(start),int(stop)
                    #print alignedread.get_reference_sequence()
                    #print alignedread.seq
                    md = string.replace(md,'C','A')
                    md = string.replace(md,'G','A')
                    md = string.replace(md,'T','A')
                    md = string.split(md,'A')
                    pos = alignedread.pos
                    chr = bamfile.getrname(alignedread.rname)
                    #if omd == '34^GA16': print md, pos
                    for i in md[:-1]:
                        try:
                            pos+=int(i)+1
                        except Exception:
                            if i == '':
                                pos+=+1
                            elif '^' in i: ### position is equal to the last position
                                pos+=int(string.split(i,'^')[0])+1
                                #pass
                        #if 'CGGATCC' in alignedread.seq: print string.split(alignedread.seq,'CGGATCC')[1],[pos]
                        try: variant_db[chr,pos,symbol]+=1
                        except Exception: variant_db[chr,pos,symbol] = 1
                        
                    #codes = map(lambda x: x[0],alignedread.cigar)
            output_bed_rows+=1
    o.close()
    bamfile.close()
    if multi==False:
        print time.time()-start_time, 'seconds to assign reads for %d entries from %d reference entries' % (output_bed_rows,reference_rows)
    #print variant_db;sys.exit()
    return variant_db

def pileupAnalysis(bam_dir,search_locations,multi=False):
    start_time = time.time()    
    bamfile = pysam.Samfile(bam_dir, "rb" )
    reference_rows=0
    output_bed_rows=0
    #https://www.biostars.org/p/76119/
    variant_db={}
    o = open (string.replace(bam_dir,'.bam','__variant.txt'),"w")
    entries = ['chr','position','rare-allele frq','type','depth','gene','variant_info','alt_frq']
    o.write(string.join(entries,'\t')+'\n')
    #print 'Analyzing',len(search_locations),'variants'
    for (chr,pos,symbol) in search_locations: ### read each line one-at-a-time rather than loading all in memory
        pos = int(pos)
        read_count=0
        reference_rows+=1
        nucleotide_frequency={}
        for pileupcolumn in bamfile.pileup(chr,pos,pos+1):
            # Skip columns outside desired range
            #print pos, pileupcolumn.pos, pileupcolumn.cigarstring, pileupcolumn.alignment.pos
            if pileupcolumn.pos == (pos-1):
                for pileupread in pileupcolumn.pileups:
                    try: nt = pileupread.alignment.query_sequence[pileupread.query_position]
                    except Exception,e:
                        if 'D' in pileupread.alignment.cigarstring:
                            nt = 'del'
                        else:
                            nt = 'ins'
                    try: nucleotide_frequency[nt]+=1
                    except Exception: nucleotide_frequency[nt]=1
        nt_freq_list=[]
        nt_freq_list_tuple=[]
        for nt in nucleotide_frequency:
            nt_freq_list.append(nucleotide_frequency[nt])
            nt_freq_list_tuple.append([nucleotide_frequency[nt],nt])
        s = sum(nt_freq_list)
        nt_freq_list.sort()
        nt_freq_list_tuple.sort()
        
        try:
            frq = float(search_locations[chr,pos,symbol])/s ### This fixes that (number of insertions from before)
        except Exception: frq = '1.000000'; print symbol, pos, nucleotide_frequency, search_locations[chr,pos,symbol]
        if (chr,pos) in insertion_db:
            #print 'insertion', chr, pos
            call = 'insertion'
            ### For insertions if the inserted base matches the reference base, incorrect freq will be reported
        elif 'del' in nucleotide_frequency:
            #frq = float(nt_freq_list[-2])/s
            call = 'del'
        else:
            #frq = float(nt_freq_list[-2])/s
            call = 'mismatch'
        if len(nt_freq_list)>1 or call == 'insertion':
            if frq>0.01:
                frq = str(frq)[:4]
                most_frequent_frq,most_frequent_nt = nt_freq_list_tuple[-1]
                try:
                    second_most_frequent_frq,second_most_frequent_nt = nt_freq_list_tuple[-2]
                    alt_frq = str(float(second_most_frequent_frq)/most_frequent_frq)
                except Exception:
                    second_most_frequent_frq = 'NA'; second_most_frequent_nt='NA'
                    alt_frq = 'NA'
                
                variant_info = most_frequent_nt+'('+str(most_frequent_frq)+')|'+second_most_frequent_nt+'('+str(second_most_frequent_frq)+')'
                entries = [chr,str(pos),str(frq),call,str(s),symbol,variant_info,alt_frq]
                o.write(string.join(entries,'\t')+'\n')
                output_bed_rows+=1
    
    o.close()
    bamfile.close()
    if multi==False:
        print time.time()-start_time, 'seconds to assign reads for %d entries from %d reference entries' % (output_bed_rows,reference_rows)

if __name__ == "__main__":
    #bam_dir = "H9.102.2.6.bam"
    #reference_dir = 'H9.102.2.6__exon.bed'
    ################  Comand-line arguments ################
    symbols=[]
    variantFile = None
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a BAM file as input in the command-line"
        print "Example: python BAMtoExonBED.py --i /Users/me/sample1.bam --r /Users/me/Hs_exon-cancer_hg19.bed"
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','species=','g=','v='])
        for opt, arg in options:
            if opt == '--i': bam_dir=arg ### A single BAM file location (full path)
            elif opt == '--species': species=arg 
            elif opt == '--g': symbols.append(arg)
            elif opt == '--v': variantFile = arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

    findGeneVariants(species,symbols,bam_dir,variants=variantFile)
