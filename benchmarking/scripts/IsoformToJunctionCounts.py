#!/usr/bin/env python
import os,sys,string,inspect
### Import python modules from an upstream directory
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
parentdir = os.path.dirname(parentdir)
sys.path.insert(0,parentdir)
import unique
import export

"""Goal: Import de novo Cufflinks transcripts and expression to infer dPSI for junctions"""

def verifyFile(filename):
    status = False
    try:
        fn=unique.filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = True;break
    except Exception: status = False
    return status

def getJunctions(species,input_dir,fpkm_tracking_dir,gtf_dir):

    ensembl_chr_db={}
    gene_intron_db={}
    ### Get the chr for each Ensembl gene to carefully match to Cufflinks Symbol
    exon_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
    fn = unique.filepath(exon_dir)
    header = True
    for line in open(fn,'rU').xreadlines():
        line = line.rstrip('\n')
        t = string.split(line,'\t')
        if header:
            header = False
        else:
            ensembl = t[0]
            exon = t[1]
            chr = t[2]
            strand = t[3]
            start =int(t[4])
            end = int(t[5])
            ensembl_chr_db[ensembl]=chr,strand
            if 'I' in exon:
                intronID = ensembl+':'+exon
                coords = [start,end]
                coords.sort()
                coords.append(intronID)
                try: gene_intron_db[ensembl].append(coords)
                except Exception: gene_intron_db[ensembl] = [coords]
        
    symbol_ensembl_db={}
    altanalyze_gene_annot_dir = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations.txt'
    fn = unique.filepath(altanalyze_gene_annot_dir)
    for line in open(fn,'rU').xreadlines():
        line = line.rstrip('\n')
        ensemblID,symbol,description,null = string.split(line,'\t')
        if 'Ensembl Gene ID' != ensemblID:
            chr,strand = ensembl_chr_db[ensemblID]
            symbol_ensembl_db[symbol]=ensemblID,chr,strand
            
    ### Import gene, isoform and expression information
    isoform_fpkm={}
    isoform_annotation={}
    header = True
    for line in open(fpkm_tracking_dir,'rU').xreadlines():
        line = line.rstrip('\n')
        t = string.split(line,'\t')
        if header:
            gi = t.index('gene_id')
            si = t.index('gene_short_name')
            ti = t.index('tracking_id')
            fi = t.index('FPKM')
            header = False
        else:
            geneID = t[gi]
            symbol = t[si]
            transcriptID = t[ti]
            fpkm = float(t[fi])
            isoform_fpkm[transcriptID]=fpkm
            isoform_annotation[transcriptID]=symbol

    transcript_junctions={}
    transcript_exons = {}
    for line in open(gtf_dir,'rU').xreadlines():
        line = line.rstrip('\n')
        t = string.split(line,'\t')
        chr = t[0]
        type = t[2]
        start = int(t[3])
        stop = int(t[4])
        strand = t[6]
        annotations = t[8]
        annotations = string.split(annotations,'"')
        geneID = annotations[1]
        transcriptID = annotations[3]
        if strand == '-':
           start, stop = stop, start
        if type == 'exon':
            try: transcript_exons[transcriptID,strand].append([chr,start,stop])
            except Exception: transcript_exons[transcriptID,strand] = [[chr,start,stop]]
    
    parent = findParentDir(fpkm_tracking_dir)
    output_file = input_dir+'/'+parent+'__junction.bed'
    io = export.ExportFile(output_file)
    io.write('track name=junctions description="TopHat junctions"\n')
    
    unmatched_gene={}
    intron_retention=0
    for (transcriptID,strand) in transcript_exons:
        exons = transcript_exons[transcriptID,strand]
        numExons = len(exons)
        fpkm = isoform_fpkm[transcriptID]
        symbol = isoform_annotation[transcriptID]
        chr,start,stop = exons[0]
        if len(symbol)>0 and symbol in symbol_ensembl_db:
            ensemblID, EnsChr,EnsStrand = symbol_ensembl_db[symbol]
            strand = EnsStrand ### Likely more accurate
            if strand == '-':
                exons.reverse()
            if chr == EnsChr: ### Double check this is the correct gene symbol for that gene ID
                i=0
                while i<(numExons-1):
                    chr,start,stop = exons[i]
                    chr,next_start,next_stop = exons[i+1]
                    junction = stop,next_start
                    """
                    if stop == 63886987 or start == 63886987 or next_start == 63886987 or next_stop == 63886987:
                        print junction, start, stop, next_start, next_stop"""
                    if len(symbol)>0 and symbol in symbol_ensembl_db:
                        if chr == EnsChr: ### Double check this is the correct gene symbol for that gene ID
                            try: transcript_junctions[ensemblID,chr,junction,strand]+=fpkm
                            except Exception: transcript_junctions[ensemblID,chr,junction,strand]=fpkm
                    i+=1

                ### Identify retained introns in transcripts
                fpkm  = int(fpkm*10)
                increment=-1
                for exon in exons:
                    chr,start,end = exon
                    if ensemblID in gene_intron_db and fpkm>0:
                        for (intron5,intron3,intronID) in gene_intron_db[ensemblID]:
                            start_end = [start,end]; start_end.sort(); start,end = start_end
                            coords = [start,end]+[intron5,intron3]
                            coords.sort()
                            if coords[0] == start and coords[-1] == end:
                                """
                                if ensemblID == 'ENSG00000167522':
                                    if strand == '+': print [start,end], [intron5,intron3];sys.exit()
                                """
                                start = intron5
                                stop = intron3
                                if strand=='-': increment = -1
                                else: increment = -1
                                outlier_start = start-10+increment; outlier_end = start+10+increment
                                junction_id = intronID+'-'+str(start)
                                exon_lengths = '10,10'; dist = '0,0'
                                entries = [chr,str(outlier_start),str(outlier_end),junction_id,str(fpkm),strand,str(outlier_start),str(outlier_end),'255,0,0\t2',exon_lengths,'0,'+dist]
                                io.write(string.join(entries,'\t')+'\n')
                                
                                ### 3' junction
                                if strand=='-': increment = 0
                                else: increment = 0
                                outlier_start = stop-10+increment; outlier_end = stop+10+increment
                                junction_id = intronID+'-'+str(stop)
                                exon_lengths = '10,10'; dist = '0,0'
                                entries = [chr,str(outlier_start),str(outlier_end),junction_id,str(fpkm),strand,str(outlier_start),str(outlier_end),'255,0,0\t2',exon_lengths,'0,'+dist]
                                io.write(string.join(entries,'\t')+'\n')            
                                intron_retention +=1
            else:
                unmatched_gene[symbol]=None
        else:
            unmatched_gene[symbol]=None
                
    print len(unmatched_gene), 'Unmatched gene symbols'
    
    num_junctions=0
    for (ensemblID,chr,junction,strand) in transcript_junctions:
        fpkm = int(transcript_junctions[ensemblID,chr,junction,strand]*10)
        if fpkm>0:
            fpkm = str(fpkm)
            junction = list(junction)
            junction.sort()
            junction_id = ensemblID+':'+str(junction[0])+'-'+str(junction[1])
            start = int(junction[0])
            end = int(junction[1])
            if strand == '-':
                alt_start = start
                alt_end = end-1
            else:
                alt_start = start
                alt_end = end-1       
            alt_start = str(alt_start)
            alt_end = str(alt_end)
            #exon_lengths = '10,10'; dist = '0,0'
            exon_lengths = '0,0'; dist = '0,0'
            entries = [chr,alt_start,alt_end,junction_id,fpkm,strand,alt_start,alt_end,'255,0,0\t2',exon_lengths,'0,'+dist]
            io.write(string.join(entries,'\t')+'\n')
            num_junctions+=1
                        
    io.close()
    print intron_retention, 'transcripts with intron retention'
    print num_junctions, 'junctions exported to:',output_file
    
def findParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    return string.split(filename,'/')[-2]

if __name__ == '__main__':
    import multiprocessing as mlp
    import getopt
    species = 'Hs'
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []
        useMultiProcessing=False
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i='])
        for opt, arg in options:
            if opt == '--i': input_dir=arg
            elif opt == '--species': species=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    dir_list = unique.read_directory(input_dir)

    for file in dir_list:
        if 'isoforms.fpkm_tracking' in file:
            fpkm_tracking_dir = input_dir+'/'+file
        elif 'transcripts.gtf' in file:
            gtf_dir = input_dir+'/'+file
            
    getJunctions(species,input_dir,fpkm_tracking_dir,gtf_dir)
