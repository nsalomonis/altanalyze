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
import copy,math
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

def parseExonReferences(bam_dir,reference_exon_bed,multi=False,intronRetentionOnly=False, MateSearch=False, species=None):
    start_time = time.time()    
    bamfile = pysam.Samfile(bam_dir, "rb" )
    reference_rows=0
    output_bed_rows=0
    exportCoordinates=True
    try:
        from import_scripts import BAMtoJunctionBED
        retainedIntrons, chrm_found, gene_coord_db = BAMtoJunctionBED.retreiveAllKnownSpliceSites(returnExonRetention=True,DesignatedSpecies=species,path=bam_dir)
    except Exception:
        #print traceback.format_exc();sys.exit()
        retainedIntrons={}
    if intronRetentionOnly==False:
        o = open (string.replace(bam_dir,'.bam','__exon.bed'),"w")
    #io = AppendOrWrite(string.replace(bam_dir,'.bam','__junction.bed'))
    io = open (string.replace(bam_dir,'.bam','__intronJunction.bed'),"w")
    if exportCoordinates:
        eo = open (reference_exon_bed[:-4]+'__minimumIntronIntervals.bed',"w")
    intron_count=0
    
    quick_look = 0
    paired = False ### Indicates if paired-end reads exist in the file
    for entry in bamfile.fetch():
        example_chromosome = bamfile.getrname(entry.rname)
        try: mate = bamfile.mate(entry); paired=True
        except Exception: pass
        quick_look+=1
        if quick_look>20: ### Only examine the first 20 reads
            break

    ### Import the gene data and remove introns that overlap with exons on the opposite or same strand
    exonData_db={}
    exon_sorted_list=[]
    for line in open(reference_exon_bed,'rU').xreadlines(): ### read each line one-at-a-time rather than loading all in memory
        line = line.rstrip('\n')
        reference_rows+=1
        #if reference_rows==6000: break
        ref_entries = string.split(line,'\t'); #'12', '6998470', '6998522', 'ENSG00000111671:E1.1_ENSE00001754003', '0', '-'
        chr,start,stop,exon,null,strand = ref_entries[:6]
        start = int(start)
        stop = int(stop)
        if 'novel' not in exon:
            exon_sorted_list.append([chr,start,stop,exon,strand])

    exon_sorted_list.sort()
    exon_sorted_filtered=[]
    exon_index=0
    
    def checkOverlap(region, query_regions):
        ### Do not trust intron retention estimates if the intron is within an exon
        chr,start,stop,exon,strand = region
        gene = string.split(exon,':')[0]
        overlap = False
        for (chr,start2,stop2,exon2,strand2) in query_regions:
            if ':E' in exon2 and gene not in exon2:
                x = [start,stop,start2,stop2]
                x.sort()
                """
                if exon2 == 'ENSG00000155657:E363.1' and exon == 'ENSG00000237298:I10.1':
                    if x[1]==start and x[2]==stop: print 'kill'
                    print x, start, stop;sys.exit()
                """
                if x[:2] == [start,stop] or x[-2:] == [start,stop]:
                    ### intron NOT overlapping with an exon
                    #print region; print chr,start2,stop2,exon2,strand2;sys.exit()
                    pass
                else: ### Partial or complete exon overlap with an intron in a second gene
                    if x[0]==start and x[-1]==stop and (start2-start)>100 and (stop-stop2)>100:
                        ### Classic small intronic RNA example
                        pass
                    elif x[1]==start and x[2]==stop:
                        ### The exon spans the entire intron
                        overlap = True
                    elif (stop2-start2)<50 and (stop-start)>400:
                        ### Considered a minor insignificant overlap that can't account for paired-end mapping
                        pass
                    
                    elif ((x[1]-x[0])<50 or (x[3]-x[2])<50 or (x[2]-x[1])<50) and (stop-start)>400:
                        ### Considered a minor insignificant overlap that can't account for paired-end mapping
                        ### Minor partial overlap of either side of the intron with an exon
                        if (stop2-start)>50 or (stop-start2)>50: ### Then it is more than a minor overlap
                             overlap = True
                    else:
                        overlap = True
                    """
                    if exon == 'ENSG00000230724:I7.1' and exon2 == 'ENSG00000255229:E1.1':
                        print (x[1]-x[0]), (x[3]-x[2]), (x[2]-x[1]), (stop-start);sys.exit()
                    """
        return overlap

    exonOverlappingIntrons=0
    totalIntrons=0
    for region in exon_sorted_list:
        chr,start,stop,exon,strand = region
        if ':I' in exon or exon in retainedIntrons:
            try:
                totalIntrons+=1
                query_regions =  exon_sorted_list[exon_index-10:exon_index]+exon_sorted_list[exon_index+1:exon_index+10]
                overlap = checkOverlap(region,query_regions)
                """
                if 'ENSG00000262880:I3.1' == exon:
                    for i in query_regions: print i
                    print chr,start,stop,exon,strand
                    print overlap;sys.exit()
                    """
                if overlap == True:
                    """
                    if exonOverlappingIntrons>10000000:
                        for i in query_regions: print i
                        print chr,start,stop,exon,strand
                        print overlap;sys.exit()
                    """
                    exon_index+=1
                    exonOverlappingIntrons+=1
                    continue
            except Exception: pass
        exon_index+=1
        #try: exonData_db[gene].append([chr,start,stop,exon,strand])
        #except Exception: exonData_db[gene]=[[chr,start,stop,exon,strand]]
        exon_sorted_filtered.append(region)

    #print exonOverlappingIntrons, 'introns overlapping with exons in distinct genes out of',totalIntrons;sys.exit()
    for (chr,start,stop,exon,strand) in exon_sorted_filtered:
        read_count=0;
        five_intron_junction_count=0
        three_intron_junction_count=0
        intronJunction={}
        exportIntervals=[]
        #if exon != 'ENSG00000167107:I10.1': continue
        if 'chr' in example_chromosome and 'chr' not in chr: ### Ensures the reference chromosome names match the query
            chr = 'chr'+chr
        try:    
            #if exon == 'ENSMUSG00000001472:E17.1':
            #chr = '12'; start = '6998470'; stop = '6998522'
            if intronRetentionOnly:
                if ':E' in exon and exon not in retainedIntrons:
                    continue
            if ':I' in exon or exon in retainedIntrons:
                INTRON = True
            else:
                INTRON = False
            start,stop=int(start),int(stop)
            regionLen = abs(start-stop)
            interval_read_count=0
            if exportCoordinates:
                if regionLen<700:
                    exportIntervals = [[start-50,stop+50]] ### Buffer intron into the exon
                else:
                    exportIntervals = [[start-50,start+350],[stop-350,stop+50]]
                    #print exportIntervals, start, stop;sys.exit()
                for interval in exportIntervals:
                    interval = map(str,interval)
                    eo.write(string.join([chr]+interval+[exon,'',strand],'\t')+'\n')
                #chr = '*'
            for alignedread in bamfile.fetch(chr, start,stop):
                proceed = True
                interval_read_count+=1
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
                                proceed = False
                                a = [X,Y]; a.sort()
                                b = [X,Y,start,stop]; b.sort()
                                if a[0]==b[1] or a[1]==b[2]: ### Hence, the read starts or ends in that interval
                                    proceed = True
                                if proceed == False:
                                    ### Also search for cases were part of the read is contained within the exon
                                    from import_scripts import BAMtoJunctionBED
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
                                if INTRON:
                                    X=int(alignedread.pos)
                                    Y=int(alignedread.pos+alignedread.alen)
                                    read_pos = [X,Y]; read_pos.sort()
                                    combined_pos = [X,Y,start,stop]; combined_pos.sort()
                                    if MateSearch==False:
                                        if read_pos[0]==combined_pos[0] or read_pos[1]==combined_pos[-1]:
                                            ### Hence, the read starts or ends OUTSIDE of that interval (Store the overlap read coordinates)
                                            if alignedread.query_name in intronJunction: ### occurs when the other mate has been stored to this dictionary
                                                ir = intronJunction[alignedread.query_name]
                                                ir.setIntronSpanningRead(read_pos)
                                                ir.setCombinedPositions(combined_pos)
                                            else:
                                                ir = IntronRetenionReads() ### first time the read-name added to this dictionary
                                                ir.setIntronSpanningRead(read_pos)
                                                ir.setCombinedPositions(combined_pos)
                                                intronJunction[alignedread.query_name] = ir
                                                if paired == False:
                                                    ir.setIntronMateRead() ### For single-end FASTQ (not accurate)
                                        else:
                                            intron_boundaries = [start,stop]; intron_boundaries.sort()
                                            if intron_boundaries[0]==combined_pos[0] and intron_boundaries[-1]==combined_pos[-1]: ### Hence, the read occurs entirely within the intron
                                                ### Store the "MATE" information (intron contained read)
                                                if alignedread.query_name in intronJunction:
                                                    ir = intronJunction[alignedread.query_name]
                                                    ir.setIntronMateRead()
                                                    found = ir.For_and_Rev_Present()
                                                else:
                                                    ir = IntronRetenionReads()
                                                    ir.setIntronMateRead()
                                                    intronJunction[alignedread.query_name] = ir
                                        if alignedread.query_name in intronJunction:
                                            ir = intronJunction[alignedread.query_name]
                                            found = ir.For_and_Rev_Present()
                                            if found: ### if the intron is less 500 (CAN CAUSE ISSUES IF READS BLEED OVER ON BOTH SIDES OF THE INTRON)
                                                combined_pos = ir.CombinedPositions()
                                                read_pos = ir.IntronSpanningRead()
                                                if read_pos[0]==combined_pos[0]:
                                                    five_intron_junction_count+=1 ### intron junction read that spans the 5' intron-exon
                                                    #print alignedread.query_name, exon, start, stop, X, Y, strand, read_pos, combined_pos
                                                elif read_pos[1]==combined_pos[-1]:
                                                    three_intron_junction_count+=1 ### intron junction read that spans the 3' intron-exon
                                            elif regionLen<500:
                                                combined_pos = ir.CombinedPositions()
                                                read_pos = ir.IntronSpanningRead()
                                                intron_read_overlap = combined_pos[2]-combined_pos[1]
                                                #print intron_read_overlap
                                                if intron_read_overlap>25:
                                                    if read_pos[0]==combined_pos[0]:
                                                        five_intron_junction_count+=1 ### intron junction read that spans the 5' intron-exon
                                                        #print alignedread.query_name, exon, start, stop, X, Y, strand, read_pos, combined_pos
                                                    elif read_pos[1]==combined_pos[-1]:
                                                        #print read_pos, combined_pos, intron_read_overlap
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
            if intronRetentionOnly==False:
                entries = [chr,str(start),str(stop),exon,null,strand,str(read_count),'0',str(int(stop)-int(start)),'0']
                o.write(string.join(entries,'\t')+'\n')
                output_bed_rows+=1
            #"""
            if INTRON and five_intron_junction_count>4 and three_intron_junction_count>4:
                interval_read_count = interval_read_count/2 ### if paired-end reads
                #print interval_read_count, five_intron_junction_count, three_intron_junction_count
                #print abs((math.log(five_intron_junction_count,2)-math.log(three_intron_junction_count,2)))
                if abs((math.log(five_intron_junction_count,2)-math.log(three_intron_junction_count,2)))<2: ### if > 4 fold difference
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
                    intron_count+=1
                    output_bed_rows+=1
                    #if output_bed_rows==1000: break
            #"""
        except Exception,e:
            #print e;sys.exit()
            ### Occurs also due to non-chromosome contigs in the annotation file
            if 'bamfile without index' in e:
                print 'Please ensure an index exists for the bam file:',bam_dir;sys.exit()
    try: o.close()
    except Exception: pass
    try: eo.close()
    except Exception: pass
    try: io.close()
    except Exception: pass
    bamfile.close()
    if multi==False:
        print time.time()-start_time, 'seconds to assign reads for %d entries from %d reference entries' % (output_bed_rows,reference_rows)

if __name__ == "__main__":
    #bam_dir = "H9.102.2.6.bam"
    #reference_dir = 'H9.102.2.6__exon.bed'
    ################  Comand-line arguments ################
    species = None
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a BAM file as input in the command-line"
        print "Example: python BAMtoExonBED.py --i /Users/me/sample1.bam --r /Users/me/Hs_exon-cancer_hg19.bed"
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','g=','r=','s=','species='])
        for opt, arg in options:
            if opt == '--i': bam_dir=arg ### A single BAM file location (full path)
            elif opt == '--s': species = arg
            elif opt == '--species': species = arg
            elif opt == '--r': reference_dir=arg ### An exon.bed reference file (created by AltAnalyze from junctions, multiBAMtoBED.py or other)
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

    parseExonReferences(bam_dir,reference_dir,intronRetentionOnly=True,species=species)

