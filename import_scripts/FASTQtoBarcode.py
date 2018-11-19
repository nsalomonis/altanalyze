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
import Bio; from Bio.Seq import Seq

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def parseFASTQFile(fn):
    count=0
    spacer='TGGT'
    global_count=0
    read2_viral_barcode={}
    read1_cellular_barcode={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        global_count+=1
        if count == 0:
            read_id = string.split(data,' ')[0][1:]
            count+=1
        elif count == 1:
            sequence = data
            count+=1
        else:
            count+=1
        if count == 4:
            count = 0
            if 'R2' in fn:
                if spacer in sequence:
                    if sequence.index(spacer) == 14:
                        viral_barcode = sequence[:48]
                        read2_viral_barcode[read_id]=viral_barcode
                else: ### Reverse complement
                    sequence = Seq(sequence)
                    sequence=str(sequence.reverse_complement())
                    if spacer in sequence:
                        if sequence.index(spacer) == 14:
                            viral_barcode = sequence[:48]
                            read2_viral_barcode[read_id]=viral_barcode    
            if 'R1' in fn:
                if 'TTTTT' in sequence:
                    cell_barcode = sequence[:16]
                    read1_cellular_barcode[read_id]=cell_barcode
                elif 'AAAAA' in sequence:  ### Reverse complement
                    sequence = Seq(sequence)
                    cell_barcode=str(sequence.reverse_complement())[:16]
                    read1_cellular_barcode[read_id]=cell_barcode
    if 'R2' in fn:
        return read2_viral_barcode
    else:
        return read1_cellular_barcode
                    
def outputPairs(fastq_dir,read1_cellular_barcode,read2_viral_barcode):
    outdir = fastq_dir+'.viral_barcodes.txt'
    o = open (outdir,"w")
    unique_pairs={}
    for uid in read2_viral_barcode:
        if uid in read1_cellular_barcode:
            cellular = read1_cellular_barcode[uid]
            viral = read2_viral_barcode[uid]
            if (viral,cellular) not in unique_pairs:
                o.write(viral+'\t'+cellular+'\n')
                unique_pairs[(viral,cellular)]=[]
    o.close()
    
if __name__ == '__main__':
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a SAM file as input in the command-line"
        print "Example: python BAMtoJunctionBED.py --i /Users/me/sample1.fastq"
        sys.exit()
    else:
        Species = None
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i='])
        for opt, arg in options:
            if opt == '--i': fastq_dir=arg ### full path of a BAM file
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    if 'R1' in fastq_dir:
        r1 = fastq_dir
        r2 = string.replace(fastq_dir,'R1','R2')
    else:
        r1 = string.replace(fastq_dir,'R2','R1')
        r2= fastq_dir
    read2_viral_barcode = parseFASTQFile(r2)
    read1_cellular_barcode = parseFASTQFile(r1)
    
    outputPairs(fastq_dir,read1_cellular_barcode,read2_viral_barcode)