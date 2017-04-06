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
import traceback

def AppendOrWrite(export_path):
    status = verifyFile(export_path)
    if status == 'not found':
        export_data = open(export_path,'w') ### Write this new file
    else:
        export_data = open(export_path,'a') ### Appends to existing file
    return export_data

def verifyFile(filename):
    status = 'not found'
    try:
        for line in open(filename,'rU').xreadlines(): status = 'found';break
    except Exception: status = 'not found'
    return status

class IntronRetenionReads():
    def __init__(self):
        self.spanning=0
        self.containing=0
    def setCombinedPositions(self,combined_pos):
        self.combined_pos = combined_pos
    def setIntronSpanningRead(self,read_pos):
        self.read_pos = read_pos
        self.spanning=1
    def setIntronMateRead(self):
        self.containing=1
    def IntronSpanningRead(self):
        return self.read_pos
    def CombinedPositions(self):
        return self.combined_pos
    def For_and_Rev_Present(self, ):
        if (self.containing+self.spanning)==2:
            return True
        else:
            return False

def parseExonReferences(bam_dir,reference_exon_bed,multi=False,intronRetentionOnly=False, MateSearch=False):
    start_time = time.time()    
    bamfile = pysam.Samfile(bam_dir, "rb" )
    reference_rows=0
    output_bed_rows=0
    o = open (string.replace(bam_dir,'.bam','__exon.bed'),"w")
    #io = AppendOrWrite(string.replace(bam_dir,'.bam','__junction.bed'))
    io = open (string.replace(bam_dir,'.bam','__intronJunction.bed'),"w")
    intron_count=0
    for line in open(reference_exon_bed,'rU').xreadlines(): ### read each line one-at-a-time rather than loading all in memory
        line = line.rstrip('\n')
        reference_rows+=1
        #if reference_rows==6000: break
        ref_entries = string.split(line,'\t'); #'12', '6998470', '6998522', 'ENSG00000111671:E1.1_ENSE00001754003', '0', '-'
        chr,start,stop,exon,null,strand = ref_entries[:6]
        read_count=0;
        five_intron_junction_count=0
        three_intron_junction_count=0
        intronJunction={}
        try:    
            #if exon == 'ENSMUSG00000001472:E17.1':
            #chr = '12'; start = '6998470'; stop = '6998522'
            if intronRetentionOnly:
                if ':E' in exon:
                    continue
            for alignedread in bamfile.fetch(chr, int(start),int(stop)):
                proceed = True
                try: cigarstring = alignedread.cigarstring
                except Exception:
                    codes = map(lambda x: x[0],alignedread.cigar)
                    if 3 in codes: cigarstring = 'N'
                    else: cigarstring = None
                
                try: read_strand = alignedread.opt('XS') ### TopHat/STAR knows which sequences are likely real splice sites so it assigns a real strand to the read
                except Exception,e:
                    #if multi == False:  print 'No TopHat strand information';sys.exit()
                    read_strand = None ### TopHat doesn't predict strand for many reads                
                
                if read_strand==None or read_strand==strand: ### Tries to ensure the propper strand reads are considered (if strand read info available)
                    if cigarstring == None: pass
                    else:
                        ### Exclude junction reads ("N")
                        if 'N' in cigarstring:
                            if intronRetentionOnly==False:
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
                        else:
                            ### Below code is for more accurate estimation of intron retention
                            #"""
                            try:
                                if 'I' in exon:
                                    X=int(alignedread.pos)
                                    Y=int(alignedread.pos+alignedread.alen)
                                    start= int(start)
                                    stop = int(stop)
                                    read_pos = [X,Y]; read_pos.sort()
                                    combined_pos = [X,Y,start,stop]; combined_pos.sort()
                                    if MateSearch==False:
                                        if read_pos[0]==combined_pos[0] or read_pos[1]==combined_pos[-1]:
                                            ### Hence, the read starts or ends OUTSIDE of that interval (Store the overlap read coordinates)
                                            if alignedread.query_name in intronJunction:
                                                ir = intronJunction[alignedread.query_name]
                                                ir.setIntronSpanningRead(read_pos)
                                                ir.setCombinedPositions(combined_pos)
                                            else:
                                                ir = IntronRetenionReads()
                                                ir.setIntronSpanningRead(read_pos)
                                                ir.setCombinedPositions(combined_pos)
                                                intronJunction[alignedread.query_name] = ir
                                        else:
                                            intron_boundaries = [start,stop]; intron_boundaries.sort()
                                            if intron_boundaries[0]==combined_pos[0] and intron_boundaries[-1]==combined_pos[-1]: ### Hence, the read occurs entirely within the intron
                                                ### Store the "MATE" information (intron contained read)
                                                if alignedread.query_name in intronJunction:
                                                    ir = intronJunction[alignedread.query_name]
                                                    ir.setIntronMateRead()
                                                else:
                                                    ir = IntronRetenionReads()
                                                    ir.setIntronMateRead()
                                                    intronJunction[alignedread.query_name] = ir
                                        if alignedread.query_name in intronJunction:
                                            ir = intronJunction[alignedread.query_name]
                                            found = ir.For_and_Rev_Present()
                                            if found:
                                                combined_pos = ir.CombinedPositions()
                                                read_pos = ir.IntronSpanningRead()
                                                if read_pos[0]==combined_pos[0]:
                                                    five_intron_junction_count+=1 ### intron junction read that spans the 5' intron-exon
                                                    #print alignedread.query_name, exon, start, stop, X, Y, strand;sys.exit()
                                                elif read_pos[1]==combined_pos[-1]:
                                                    three_intron_junction_count+=1 ### intron junction read that spans the 3' intron-exon
                                    else:
                                        mate = bamfile.mate(alignedread) ### looup the paired-end mate for this read
                                        try: cigarstring = mate.cigarstring
                                        except Exception:
                                            codes = map(lambda x: x[0],mate.cigar)
                                            if 3 in codes: cigarstring = 'N'
                                            else: cigarstring = None
                                        if 'N' not in cigarstring:
                                            RX=int(mate.pos)
                                            RY=int(mate.pos+mate.alen)
                                            intron_boundaries = [start,stop]; intron_boundaries.sort()
                                            combined_pos2 = [RX,RY,start,stop]; combined_pos2.sort()
                                            if intron_boundaries[0]==combined_pos2[0] and intron_boundaries[-1]==combined_pos2[-1]:
                                                if read_pos[0]==intron_boundaries[0]:
                                                    five_intron_junction_count+=1 ### intron junction read that spans the 5' intron-exon
                                                    #print alignedread.query_name,exon, start, stop, X, Y, RX, RY, strand, read_strand;sys.exit()
                                                elif read_pos[1]==intron_boundaries[-1]:
                                                    three_intron_junction_count+=1 ### intron junction read that spans the 3' intron-exon
                                        
                            except Exception,e: ### Usually an unmapped read
                                #print traceback.format_exc();sys.exit()
                                pass
                             #"""
                    if proceed: read_count+=1
            entries = [chr,str(start),str(stop),exon,null,strand,str(read_count),'0',str(int(stop)-int(start)),'0']
            o.write(string.join(entries,'\t')+'\n')
            output_bed_rows+=1
            if 'I' in exon and five_intron_junction_count>0 and three_intron_junction_count>0:
                intron_count+=1
            #"""
            if 'I' in exon and five_intron_junction_count>0 and three_intron_junction_count>0:
                if strand=='-': increment = -1
                else: increment = -1
                total_count = five_intron_junction_count+three_intron_junction_count
                outlier_start = start-10+increment; outlier_end = start+10+increment
                junction_id = exon+'-'+str(start)
                exon_lengths = '10,10'; dist = '0,0'
                entries = [chr,str(outlier_start),str(outlier_end),junction_id,str(five_intron_junction_count),strand,str(outlier_start),str(outlier_end),'255,0,0\t2',exon_lengths,'0,'+dist]
                io.write(string.join(entries,'\t')+'\n')
                
                ### 3' junction
                if strand=='-': increment = 0
                else: increment = 0
                outlier_start = stop-10+increment; outlier_end = stop+10+increment
                junction_id = exon+'-'+str(stop)
                exon_lengths = '10,10'; dist = '0,0'
                entries = [chr,str(outlier_start),str(outlier_end),junction_id,str(three_intron_junction_count),strand,str(outlier_start),str(outlier_end),'255,0,0\t2',exon_lengths,'0,'+dist]
                io.write(string.join(entries,'\t')+'\n')
            #"""
        except Exception,e:
            #print e;sys.exit()
            ### Occurs also due to non-chromosome contigs in the annotation file
            if 'bamfile without index' in e:
                print 'Please ensure an index exists for the bam file:',bam_dir;sys.exit()
    o.close()
    io.close()
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

    parseExonReferences(bam_dir,reference_dir,intronRetentionOnly=False)

