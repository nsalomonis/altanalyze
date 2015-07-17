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

def parseExonReferences(bam_dir,reference_exon_bed,multi=False):
    start_time = time.time()    
    bamfile = pysam.Samfile(bam_dir, "rb" )
    reference_rows=0
    output_bed_rows=0
    o = open (string.replace(bam_dir,'.bam','__exon.bed'),"w")
    for line in open(reference_exon_bed,'rU').xreadlines(): ### read each line one-at-a-time rather than loading all in memory
        line = line.rstrip('\n')
        reference_rows+=1
        ref_entries = string.split(line,'\t'); #'12', '6998470', '6998522', 'ENSG00000111671:E1.1_ENSE00001754003', '0', '-'
        chr,start,stop,exon,null,strand = ref_entries[:6]
        read_count=0
        try:
            #if exon == 'ENSMUSG00000001472:E17.1':
            #chr = '12'; start = '6998470'; stop = '6998522'
            for alignedread in bamfile.fetch(chr, int(start),int(stop)):
                proceed = True
                if alignedread.cigarstring == None: pass
                else:
                    ### Exclude junction reads ("N")
                    if 'N' in alignedread.cigarstring:
                        X=int(alignedread.pos)
                        Y=int(alignedread.pos+alignedread.alen)
                        start= int(start)
                        stop = int(stop)
                        proceed = False
                        a = [X,Y]; a.sort()
                        b = [X,Y,start,stop]; b.sort()
                        if a[0]==b[1] or a[1]==b[2]: ### Hence, the read starts or ends in that interval
                            proceed = True
                        if proceed == False:
                            ### Also search for cases were part of the read is contained within the exon
                            import BAMtoJunctionBED
                            coordinates,up_to_intron_dist = BAMtoJunctionBED.getSpliceSites(alignedread.cigar,X)
                            for (five_prime_ss,three_prime_ss) in coordinates:
                                five_prime_ss,three_prime_ss=int(five_prime_ss),int(three_prime_ss)
                                if five_prime_ss==start or three_prime_ss==start or five_prime_ss==stop or three_prime_ss==stop:
                                    proceed = True
                                    #print five_prime_ss, three_prime_ss, start, stop;sys.exit()
                if proceed: read_count+=1
            entries = [chr,str(start),str(stop),exon,null,strand,str(read_count),'0',str(int(stop)-int(start)),'0']
            o.write(string.join(entries,'\t')+'\n')
            output_bed_rows+=1
        except Exception,e:
            ### Occurs also due to non-chromosome contigs in the annotation file
            if 'bamfile without index' in e:
                print 'Please ensure an index exists for the bam file:',bam_dir;sys.exit()
    o.close()
    bamfile.close()
    if multi==False:
        print time.time()-start_time, 'seconds to assign reads for %d entries from %d reference entries' % (output_bed_rows,reference_rows)

if __name__ == "__main__":
    #bam_dir = "H9.102.2.6.bam"
    #reference_dir = 'H9.102.2.6__exon.bed'
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a BAM file as input in the command-line"
        print "Example: python BAMtoExonBED.py --i /Users/me/sample1.bam --r /Users/me/Hs_exon-cancer_hg19.bed"
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','g=','r='])
        for opt, arg in options:
            if opt == '--i': bam_dir=arg ### A single BAM file location (full path)
            elif opt == '--r': reference_dir=arg ### An exon.bed reference file (created by AltAnalyze from junctions, multiBAMtoBED.py or other)
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

    parseExonReferences(bam_dir,reference_dir)

