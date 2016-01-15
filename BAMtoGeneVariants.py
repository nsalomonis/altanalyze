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

import pysam
import string,os,sys,copy
import time
import getopt

def findGeneVariants(species,symbols,bam_dir):
    search_locations = geneCoordinates(species,symbols)
    variant_db = findVariants(bam_dir,search_locations)
    
    variant_filtered_db={}
    for var in variant_db:
        #print var, variant_db[var]
        if variant_db[var]>8:
            #print var,variant_db[var]
            variant_filtered_db[var] = variant_db[var]
            
    pileupAnalysis(bam_dir,variant_filtered_db)

def geneCoordinates(species,symbols):
    genes=[]
    import EnsemblImport
    ensembl_annotation_db = EnsemblImport.reimportEnsemblAnnotations(species,symbolKey=True)
    for symbol in symbols:
        ens_geneid = ensembl_annotation_db[symbol]
        genes.append((ens_geneid,symbol))
    
    ### Get gene genomic locations
    gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,'RNASeq','key_by_array')
    search_locations=[]
    for (gene,symbol) in genes:
        chr,strand,start,end = gene_location_db[gene]
        if symbol == 'SRSF10': chr = 'chr1'; strand = '-'; start = '24295573'; end = '24306953'
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
            for alignedread in bamfile.fetch(chr, int(start),int(stop)):
                md = alignedread.opt('MD')
                omd = md
                try:
                    int(md) ### If an integer, no mismatches or deletions/insertion present
                    continue
                except Exception:
                    #print chr, int(start),int(stop)
                    #print alignedread.get_reference_sequence()
                    #print alignedread.seq;sys.exit()
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
    entries = ['chr','position','rare-allele frq','type','depth','gene']
    o.write(string.join(entries,'\t')+'\n')
    for (chr,pos,symbol) in search_locations: ### read each line one-at-a-time rather than loading all in memory
        pos = int(pos)
        read_count=0
        reference_rows+=1
        nucleotide_frequency={}
        for pileupcolumn in bamfile.pileup(chr,pos,pos+1):
            # Skip columns outside desired range
            if pileupcolumn.pos == (pos-1):
                for pileupread in pileupcolumn.pileups:
                    try: nt = pileupread.alignment.query_sequence[pileupread.query_position]
                    except Exception,e:
                        if 'D' in pileupread.alignment.cigarstring:
                            nt = 'del'
                        else:
                            continue
                    try: nucleotide_frequency[nt]+=1
                    except Exception: nucleotide_frequency[nt]=1

        nt_freq_list=[]
        for nt in nucleotide_frequency:
            nt_freq_list.append(nucleotide_frequency[nt])
        s = sum(nt_freq_list)
        nt_freq_list.sort()

        if len(nt_freq_list)>1:
            frq = float(nt_freq_list[-2])/s
            if frq>0.01:
                if 'del' in nucleotide_frequency: call = 'del'
                else: call = 'mismatch'
                frq = str(frq)[:4]
                entries = [chr,str(pos),str(frq),call,str(s),symbol]
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
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a BAM file as input in the command-line"
        print "Example: python BAMtoExonBED.py --i /Users/me/sample1.bam --r /Users/me/Hs_exon-cancer_hg19.bed"
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','species=','g='])
        for opt, arg in options:
            if opt == '--i': bam_dir=arg ### A single BAM file location (full path)
            elif opt == '--species': species=arg 
            elif opt == '--g': symbols.append(arg)
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

    findGeneVariants(species,symbols,bam_dir)
