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

import pysam
import string,os,sys,copy,getopt
import time

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
            
def parseJunctionEntries(bam_dir,multi=False):
    start = time.time()
    try: import collections; junction_db=collections.OrderedDict()
    except Exception:
        try: import ordereddict; junction_db = ordereddict.OrderedDict()
        except Exception: junction_db={}
    original_junction_db = copy.deepcopy(junction_db)
    
    bamf = pysam.Samfile(bam_dir, "rb" )
    chromosome = False
    count=0
    jid = 1
    prior_jc_start=0
    o = open (string.replace(bam_dir,'.bam','__junction.bed'),"w")
    o.write('track name=junctions description="TopHat junctions"\n')
    outlier_start = 0; outlier_end = 0; read_count = 0
    for entry in bamf.fetch():
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
                junction_db = copy.deepcopy(original_junction_db) ### Re-set this object
                jid+=1
            chromosome = bamf.getrname( entry.rname ) 
            X=entry.pos
            Y=entry.pos+entry.alen
            prior_jc_start = X
            if entry.is_reverse:
                strand = '-' ### This is the strand the seq aligns to but not necessarily the REAL strand the mRNA aligns to (see XS below)
            else:                
                strand = '+'
            try: tophat_strand = entry.opt('XS') ### TopHat knows which sequences are likely real splice sites so it assigns a real strand to the read
            except Exception:
                if multi == False:
                    print 'No TopHat strand information';sys.exit()
            coordinates,up_to_intron_dist = getSpliceSites(entry.cigar,X)
            for (five_prime_ss,three_prime_ss) in coordinates:
                jc = five_prime_ss,three_prime_ss
                #print X, Y, jc, entry.cigarstring, entry.cigar
                try: junction_db[chromosome,jc,tophat_strand].append([X,Y,up_to_intron_dist])
                except Exception: junction_db[chromosome,jc,tophat_strand] = [[X,Y,up_to_intron_dist]]
            count+=1
    writeJunctionBedFile(junction_db,jid,o) ### One last read-out
    if multi == False:
        print time.time()-start, 'seconds required to parse the BAM file'
    o.close()
    bamf.close()

if __name__ == "__main__":
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a BAM file as input in the command-line"
        print "Example: python BAMtoJunctionBED.py --i /Users/me/sample1.bam"
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i='])
        for opt, arg in options:
            if opt == '--i': bam_dir=arg ### full path of a BAM file
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
            
    try: parseJunctionEntries(bam_dir)
    except ZeroDivisionError:
        print [sys.argv[1:]],'error'; error

