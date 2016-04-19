#!/usr/bin/env python
import numpy as np
import pylab as pl
import sys,string
import os
from misopy import index_gff as ifg
import subprocess
import multiprocessing
import time
import unique
import traceback
samplelis=[]
samp=[]
sample_read={}
PSIJunctions=[]
new_lis=[]
lis=[]
added_gff_entries={}

def reimportFeatures(featureFile):
    gene_event_db={}
    featureFile = unique.filepath(featureFile)
    head=0
    for line in open(featureFile,'rU').xreadlines():
     #for k in range(len(a['AltAnalyze_ID'])):
     if head ==0: head=1
     else:
	line = line.rstrip('\n')
	event=string.split(line,'\t')[0] #example event: ENSMUSG00000025915:E17.2-E17.5=chr1:9885753-9886047
	event = string.replace(event,':','__')
	event_split=string.split(event,'__')
	for i in range(len(event_split)):
	    if "ENS" in event_split[i] or '00000' in event_split[i]:
		if '-' in event_split[i]:
		    ji=string.split(event_split[i],'-')
		    gene=ji[1]
		else:
		    gene=event_split[i]
		if gene in gene_event_db:
		    gene_event_db[gene].append(event)
		else:
		    gene_event_db[gene]=[event]
    return gene_event_db

def writegene(chrom,start_l,end_l,strand,uid):
    #start_l = str(int(start_l)-1000)
    #end_l = str(int(end_l)+2000)
    if 'M' not in chrom:
	export_in.write(chrom+'\t'+'SE'+'\t'+'gene'+'\t'+start_l+'\t'+end_l+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'ID='+uid+';'+'Name='+uid+';'+'\n')
	export_in.write(chrom+'\t'+'SE'+'\t'+'mRNA'+'\t'+start_l+'\t'+end_l+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'ID='+uid+'.A;'+'Parent='+uid+';'+'\n')

def manualWriteExon(chrom,start_l,end_l,strand,uid):
    if 'M' not in chrom:
	if strand== '-': i1=-5; i2=5
	else: i1=5; i2=-5 
	export_in.write(chrom+'\t'+'SE'+'\t'+'exon'+'\t'+start_l+'\t'+str(int(float(start_l)))+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'ID='+uid+';'+'Parent='+uid+'.A;'+'\n')
	export_in.write(chrom+'\t'+'SE'+'\t'+'exon'+'\t'+end_l+'\t'+str(int(float(end_l)))+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'ID='+uid+';'+'Parent='+uid+'.A;'+'\n')

def importReciprocalJunctions(inputpsi,PSIJunctions):
    ### Also include other predicted splicing events
    alt_dir = string.split(inputpsi,'AlternativeOutput')[0]+'AlternativeOutput'
    files = unique.read_directory(alt_dir)
    added=0
    already_added=0
    for file in files:
	if 'ASPIRE-exon-inclusion-results' in file or 'linearregres-exon-inclusion-results' in file:
            alt_exon_path = alt_dir+'/'+file
	    header=True
	    for line in open(alt_exon_path,'rU').xreadlines():
		line = line.rstrip(os.linesep)
		if header: header=False
		else:
		    t=string.split(line,'\t')
		    inclusion_junction = t[8]
		    exclusion_junction = t[10]
		    pair = inclusion_junction+' '+exclusion_junction
		    pair = string.replace(pair,':','__')
		    if pair in PSIJunctions:
			already_added+=1
		    else:
			PSIJunctions.append(pair)
			added+=1
    return PSIJunctions
    
def importPSIJunctions(fname):
    fname = unique.filepath(fname)
    header=True
    for line in open(fname,'rU').xreadlines():
        line = line.rstrip(os.linesep)
	if header: header = False
	else:
	    t=string.split(line,'\t')
	    try:
		### Re-order these to have the exclusion be listed first
		j1a,j1b = string.split(t[2],'-')
		j2a,j2b = string.split(t[3],'-')
		j1a = string.split(j1a,':')[1]
		j2a = string.split(j2a,':')[1]
		j1a = int(float(string.split(j1a,'.')[0][1:]))
		j1b = int(float(string.split(j1b,'.')[0][1:]))
		j2a = int(float(string.split(j2a,'.')[0][1:]))
		j2b = int(float(string.split(j2b,'.')[0][1:]))
		#print [j1a,j2a,j1b,j2b], t[2], t[3]
		if j1a>j2a or j1b<j2b:
		    val = t[2]+' '+t[3]
		else:
		    val=t[3]+' '+t[2]
	    except Exception:
		#print traceback.format_exc();sys.exit()
		val=t[2]+' '+t[3]
	    if '-' not in t[2]:
		val = t[3]+' '+t[2]
	    val = string.replace(val,":","__")
	    lis.append(val)
	    #print t[0]
    return lis

def write_index(gene_n,new_lis):
    
    #print gene_n
    #print new_lis
    #print gene_n, new_lis
    valid=0
    ch=string.split(gene_n," ")
    if '-' in ch[1]: ch1=ch[1]
    else: ch1=ch[0]

    gene_mo=string.split(ch1,'__')
    gene_mo1=gene_mo[1].replace('E',"")
    gene_mo1=gene_mo1.replace('I',"")
    f=gene_mo1.split('-')
    c=f[0]
    
    k=len(gene_mo)
    gen_mo2=gene_mo[k-1].replace('E',"")
    gen_mo2=gen_mo2.replace('I',"")
    d=gen_mo2
    
    if ('_' in c):
	c1=string.split(c,'_')
	c=c1[0]
    if ('_' in d):
	c1=string.split(d,'_')
	d=c1[0]
    
    count=0
    for j in range(len(new_lis)):
       count=0
       if ch1 in new_lis[j]:
    	    c1=string.split(new_lis[j],"=")
	    chro=string.split(c1[1],'__')
	    chrom=chro[0]
	    lengt=string.split(chro[1],'-')
	    start_l=int(lengt[0])
	    end_l=int(lengt[1])
	    if start_l>end_l:
		a='-'
		b=start_l
		start_l=end_l
		end_l=b
	    else: a='+'
            count=1
	    break
       else: continue
    
    if count==1:
      #print start_l,end_l
      for j in range(len(new_lis)):
	q=string.split(new_lis[j],"=")
	if '-' in q[0]:
	    continue
	else:
	    cho=string.split(q[0],'__')
	    ge=cho[1].replace('E',"")
	    #print ge
	    if ('I' not in ge) and (';' not in ge) and ('-' not in ge):
             #if (ge >=c and ge <= d):
		#print g
                chro=string.split(q[1],'__')
		chrom=chro[0]
		lengt=string.split(chro[1],'-')
		start_e=int(lengt[0])
		end_e=int(lengt[1])
		valid=0
		##sp=string.split(lis1[i],':')
		#if sp[1]==gene_mo[1]:
	        if(start_e >= start_l and end_e <= end_l):
		    valid=1 
	        elif(start_e< start_l and end_e == start_l):
		    start_l=start_e
		    #if start_e==end_e: start_e=int(end_e)-10
		    
		    valid=1
		elif(end_e>end_l and start_e ==end_l):
		    end_l=end_e
		    #if start_e==end_e: end_e=int(start_l)+10
		    valid=1

		if valid==1:
		    #print chrom,start_l,end_l,a,ch1
		    if 'M' not in chrom:
			key = q[0],ch1
			if key not in added_gff_entries: 
			    export_in.write(chrom+'\t'+'SE'+'\t'+'exon'+'\t'+str(start_e)+'\t'+str(end_e)+'\t'+'.'+'\t'+a+'\t'+'.'+'\t'+'ID='+q[0]+';'+'Parent='+ch1+'.A;'+'\n') 
			    added_gff_entries[key]=[]
	    else:
		continue
		#export_in.write(chrom+'\t'+'SE'+'\t'+'exon'+'\t'+start_l+'\t'+end_l+'\t'+'.'+'\t'+a+'\t'+'.'+'\t'+'ID='+q[0]+'Parent='+ch[1]+'.A;'+'\n')

    start_l=int(start_l)-300 ### increase the window size
    end_l=int(end_l)+300 ### increase the window size
    return chrom,str(start_l),str(end_l),a,ch1

def Indexing(counts,inputpsi,output):
    
    gene_feature_db=reimportFeatures(counts)
    PSIJunctions=importPSIJunctions(inputpsi)
    PSIJunctions=importReciprocalJunctions(inputpsi,PSIJunctions)
    for i in range(len(PSIJunctions)):
	#if 'ENSMUSG00000066007:E5.1-ENSMUSG00000063245:E3.1' in PSIJunctions[i]:
        event_split=string.split(PSIJunctions[i],"__")
	gene = event_split[0]
	
        if gene in gene_feature_db:
	    try: chrom,start_l,end_l,a,ch1=write_index(PSIJunctions[i],gene_feature_db[gene])
	    except Exception: pass ### Not handling trans-splicing well

	for k in range(len(event_split)):
		if "ENS" in event_split[k] or '00000' in event_split[k]:
		    if '-' in event_split[k]:
			ji=string.split(event_split[k],'-')
			geneID=ji[1]
		    else:
			geneID=event_split[k]
		    if geneID != event_split[0]:
		        if geneID in gene_feature_db:
		           try: chrom,start_l,end_l,a,ch1= write_index(PSIJunctions[i],gene_feature_db[geneID])
			   except Exception: pass ### Not handling trans-splicing well
	try: writegene(chrom,start_l,end_l,a,ch1) #('chrX', '38600355', '38606652', '+', 'ENSMUSG00000000355__E1.2-E5.1')
        except Exception: pass
    
    for gene in gene_feature_db:
	for event in gene_feature_db[gene]:
	    if '100.1' in event:
		event, position_data = string.split(event,'=')
		chr, pos = string.split(position_data,'__')
		start,stop = string.split(pos,'-')
		if float(stop)>float(start): strand = '+'
		else: strand = '-'
		geneID = string.split(event,'__')[0] ### just use the gene ID
		#chrom,start_l,end_l,a,ch1= write_index(event+' '+string.replace(event,'100','99'),gene_feature_db[geneID])
		manualWriteExon(chr,start,stop,strand,geneID)
		writegene(chr,start,stop,strand,geneID)
	
    time.sleep(2)
    export_in.close()
    newout=findParentDir(output)+'sashimi_index/'
    try: ifg.index_gff(output,newout)
    except Exception:
	print traceback.format_exc();sys.exit()
	print('error in indexing')

def findParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1
    return filename[:x]
        
def remoteIndexing(species,fl):
    global export_in
    try:
	countinp = fl.CountsFile()
	root_dir = fl.RootDir()
    except Exception:
	root_dir = fl
        search_dir = root_dir+'/ExpressionInput'
        files = unique.read_directory(search_dir)
        for file in files:
            if 'counts.' in file and 'steady-state.txt' not in file:
                    countinp = search_dir+'/'+file
		    
    inputpsi = root_dir+'/AltResults/AlternativeOutput/'+species+'_RNASeq_top_alt_junctions-PSI.txt'
    outputdir=findParentDir(inputpsi)
    output=outputdir+"events_sashimi.gff"
    export_in=open(output,'w')
    
    ### Sometimes only junctions are in the count file so create a new file with detected junctions and all exons
    featuresEvaluated = extractFeatures(species,countinp)

    Indexing(featuresEvaluated,inputpsi,output)
  
def extractFeatures(species,countinp):
    import export
    ExonsPresent=False
    lastgene = None
    lastend = None
    genes_detected={}
    count=0
    first_last_exons = {} ### Make a fake junction comprised of the first and last exon
    if 'counts.' in countinp:
	feature_file = string.replace(countinp,'counts.','features.')
	fe = export.ExportFile(feature_file)
	firstLine = True
	for line in open(countinp,'rU').xreadlines():
	    if firstLine: firstLine=False
	    else:
		feature_info = string.split(line,'\t')[0]
		fe.write(feature_info+'\n')
		junction_annotation = string.split(feature_info,'=')[0]
		if '-' in junction_annotation:
		    geneid = string.split(junction_annotation,':')[0]
		    genes_detected[geneid]=[]
		if ExonsPresent == False:
		    exon = string.split(feature_info,'=')[0]
		    if '-' not in exon:
			ExonsPresent = True
			
	### Add exon-info if necessary

	exons_file = unique.filepath('AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt')
        firstLine = True
	for line in open(exons_file,'rU').xreadlines():
	    if firstLine: firstLine=False
	    else:
		line = line.rstrip('\n')
		t = string.split(line,'\t')
		gene,exon,chr,strand,start,end = t[:6]
		if gene!=lastgene:
		    if len(genes_detected)==0 or gene in genes_detected: ### restrict to detected genes
			first_last_exons[gene] = [(chr,start)]
		    if len(genes_detected)==0 or lastgene in genes_detected: ### restrict to detected genes
			try: first_last_exons[lastgene].append(lastend)
			except Exception:
			    pass ### occurs for the first gene	
		if ExonsPresent == False:
		    fe.write(gene+':'+exon+'='+chr+':'+start+'-'+end+'\n')
		lastgene = gene; lastend = end
	if len(genes_detected)==0 or lastgene in genes_detected:
	    first_last_exons[lastgene].append(lastend)
	
	### Add a fake junction for the whole gene
	for gene in first_last_exons:
	    (chr,start),end = first_last_exons[gene]
	    fe.write(gene+':E1.1-E100.1'+'='+chr+':'+start+'-'+end+'\n')
	fe.close()
    return feature_file	

def obtainTopGeneResults():
    pass

if __name__ == '__main__':
    remoteIndexing('Mm','/Volumes/salomonis1/projects/Grimes/Grimes_Single_cell_bams/')
    sys.exit()
    #"""
    countinp = os.path.abspath(os.path.expanduser(sys.argv[1]))
    inputpsi = os.path.abspath(os.path.expanduser(sys.argv[2]))
    outputdir=findParentDir(inputpsi)
    output=outputdir+"events_sashimi.gff"
    export_in=open(output,'w')
    Indexing(countinp,inputpsi,output)
    #"""