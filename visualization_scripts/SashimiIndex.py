#!/usr/bin/env python
import numpy as np
import pylab as pl
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

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
added_gff_entries={}

class PositionData:
    def __init__(self,position_str):
	self.chr,interval = string.split(position_str,'__') #chr16__18810389-18807423
	self.start,self.end = string.split(interval,'-')
    def Chr(self): return self.chr
    def Start(self): return int(self.start)
    def End(self): return int(self.end)

def reimportFeatures(featureFile):
    """ Import the exon and gene coordinates """
    gene_event_db={}
    featureFile = unique.filepath(featureFile)
    head=0
    for line in open(featureFile,'rU').xreadlines():
     #for k in range(len(strand['AltAnalyze_ID'])):
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
		featureID,position = string.split(event,'=') ### store the feature (exon or junction) position and ID separately
		pd = PositionData(position)
		if gene in gene_event_db:
		    feature_db = gene_event_db[gene]
		    feature_db[featureID] = pd
		else:
		    feature_db = {featureID:pd}
		    gene_event_db[gene]=feature_db
    return gene_event_db

def writegene(chromosome,junction_start,junction_end,strand,uid):
    #junction_start = str(int(junction_start)-1000)
    #junction_end = str(int(junction_end)+2000)
    temp = [int(junction_start),int(junction_end)]
    temp.sort()
    junction_start,junction_end = str(temp[0]),str(temp[1])
    if 'M' not in chromosome:
	gff_export_obj.write(chromosome+'\t'+'SE'+'\t'+'gene'+'\t'+junction_start+'\t'+junction_end+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'ID='+uid+';'+'Name='+uid+';'+'\n')
	gff_export_obj.write(chromosome+'\t'+'SE'+'\t'+'mRNA'+'\t'+junction_start+'\t'+junction_end+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'ID='+uid+'.STRAND;'+'Parent='+uid+';'+'\n')

def manualWriteExon(chromosome,junction_start,junction_end,strand,uid):
    if 'M' not in chromosome:
	if strand== '-': i1=-5; i2=5
	else: i1=5; i2=-5 
	gff_export_obj.write(chromosome+'\t'+'SE'+'\t'+'exon'+'\t'+junction_start+'\t'+str(int(float(junction_start)))+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'ID='+uid+';'+'Parent='+uid+'.STRAND;'+'\n')
	gff_export_obj.write(chromosome+'\t'+'SE'+'\t'+'exon'+'\t'+junction_end+'\t'+str(int(float(junction_end)))+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'ID='+uid+';'+'Parent='+uid+'.STRAND;'+'\n')

def importReciprocalJunctions(PSIFileDir,PSIJunctions):
    ### Also include other predicted splicing events from ASPIRE or LinearRegression
    alt_dir = string.split(PSIFileDir,'AlternativeOutput')[0]+'AlternativeOutput'
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
		    inclusion_junction = string.replace(t[8],':','__')
		    exclusion_junction = string.replace(t[10],':','__')
		    pair = inclusion_junction,exclusion_junction
		    if pair in PSIJunctions:
			already_added+=1
		    else:
			PSIJunctions.append(pair)
			added+=1
    return PSIJunctions
    
def importPSIJunctions(fname):
    All_PSI_Reciprocol_Junctions=[]
    fname = unique.filepath(fname)
    header=True
    for line in open(fname,'rU').xreadlines():
        line = line.rstrip(os.linesep)
	if header: header = False
	else:
	    t=string.split(line,'\t')
	    junction1 = t[2]
	    junction2 = t[3]
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
		event1 = string.replace(junction2,":","__") ### first listed junction
		event2 = string.replace(junction2,":","__") ### second listed junction
		if j1a>j2a or j1b<j2b:
		    event_pair = event1,event2
		else:
		    event_pair=event2,event1
	    except Exception:
		#print traceback.format_exc();sys.exit()
		event_pair=event1,event2
	    if '-' not in event1:
		event_pair = event2,event1
	    All_PSI_Reciprocol_Junctions.append(event_pair)
    return All_PSI_Reciprocol_Junctions

def checkForNovelJunctions(exon):
    ### if a novel junction, make this a non-novel for indexing
    exon_alt = string.replace(exon,'__',':') #revert back to :
    if '_' in exon_alt:
	exon = string.split(exon_alt,'_')[0] ### Only report the non-novel prefix exonID
	exon = string.replace(exon,':','__')
    if 'I' in exon:
	exon = string.replace(exon,'I','E')
	### We should also increment the exon ID by 1, to skip the intron
	### e.g., ENSMUSG00000000934__I6.1 => ENSMUSG00000000934__E7.1
	try:
	    gene,block_region = string.split(exon,'__E')
	    block,region = string.split(block_region,'.')
	    exon = gene+'__E'+str(int(block)+1)+'.'+region
	except Exception:
	    block,region = string.split(exon[1:],'.')
	    exon = 'E'+str(int(block)+1)+'.'+region
	    
    if 'U0' in exon:
	exon = string.replace(exon,'U0','E1')
    return exon

def exportToGff(incl_junction,excl_junction,feature_positions,gene):
    """ Write the exon and gene coordinates to strand gff file for any alternatively regulated and flanking exons """
    
    proceed = False
    if '-' in excl_junction: selected_event=excl_junction ### should be the exclusion junction
    else: selected_event=incl_junction

    ### For the first and last exon, determine their positions then later loop through the exons in between
    gene_prefix = gene+'__'
    pds = feature_positions[selected_event]
    e1,e2 = string.split(selected_event,'-')
    e1 = checkForNovelJunctions(e1)
    e2 = checkForNovelJunctions(e2)
    chr = pds.Chr()
    if pds.Start() > pds.End(): strand = '-'
    else: strand = '+'
    try: pd1 = feature_positions[e1]
    except Exception:
	e1 = string.split(e1,'.')[0]+'.1' ### occurs with IDS such as ENSMUSG00000002107__E26.27, where .27 is not in the database (not sure why)
	try: pd1 = feature_positions[e1]
	except Exception: ### Occurs with gene entry, e.g., ENSMUSG00000028180__E1.1-E100.1
	    pd1 = PositionData(chr+'__'+str(pds.Start())+'-'+str(pds.Start())) #chr3:157544964-157545122
	    selected_event = gene  ### Again, set this for gene entry coordinates
    try: pd2 = feature_positions[gene_prefix+e2]
    except Exception:
	e2 = string.split(e2,'.')[0]+'.1' ### occurs with IDS such as ENSMUSG00000002107__E26.27, where .27 is not in the database (not sure why)
	try: pd2 = feature_positions[gene_prefix+e2]
	except Exception: ### Occurs with gene entry, e.g., ENSMUSG00000028180__E1.1-E100.1
	    pd2 = PositionData(chr+'__'+str(pds.End())+'-'+str(pds.End())) #chr3:157544964-157545122
	    selected_event = gene ### Again, set this for gene entry coordinates
    #print pd.Start(),pd.End(), pd1.Start(),pd1.End(), pd2.Start(), pd2.End()
    first_start = pd1.Start()
    last_end = pd2.End()
    
    ### Now, loop through all gene features and only select the ones that are in between the spliced exons
    for exonID in feature_positions:
	proceed = False
	if '-' not in exonID: ### Hence strand junction
	    geneID, block_region = string.split(exonID,'__')
	    exon_region = string.replace(block_region,'E','')
	    if ('I' not in exon_region) and (';' not in exon_region) and ('-' not in exon_region):
		pd = feature_positions[exonID]
		exon_start = pd.Start()
		exon_end = pd.End()
	        if(exon_start >= first_start and exon_end <= last_end):
		    proceed = True
	        if(exon_start <= first_start and exon_end >= last_end):
		    proceed = True

		if proceed:
		    export_string = chr+'\t'+'SE'+'\t'+'exon'+'\t'+str(exon_start)+'\t'+str(exon_end)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'ID='+exonID+';'+'Parent='+selected_event+'.STRAND;'+'\n'
		    if 'M' not in chr:
			key = exonID,selected_event
			#if key not in added_gff_entries: 
			gff_export_obj.write(export_string) 
			added_gff_entries[key]=[]
    if strand == '-':
	junction_exon_start=int(pd2.Start())-300 ### increase the window size
	junction_exon_end=int(pd1.End())+300 ### increase the window size
    else:
	junction_exon_start=int(pd1.Start())-300 ### increase the window size
	junction_exon_end=int(pd2.End())+300 ### increase the window size
    return chr,str(junction_exon_start),str(junction_exon_end),strand,selected_event

def Indexing(counts,PSIFileDir,output):
    """ Indexing is the process of creating strand machine readable binary for SashimiPlot """
    
    gene_feature_db=reimportFeatures(counts) ### import the exon and gene coordinates
    PSIJunctions=importPSIJunctions(PSIFileDir) ### Retreive all PSI junction pairs
    PSIJunctions=importReciprocalJunctions(PSIFileDir,PSIJunctions) ### Include other Junctions from ASPIRE/LinRegress
    exported = False
    for (junction1,junction2) in PSIJunctions:
	#if 'ENSMUSG00000066007:E5.1-ENSMUSG00000063245:E3.1' in PSIJunctions[event_pair]:
	geneID = string.split(junction1,'__')[0]
        if geneID in gene_feature_db:
	    try:
		chr,junction_start,junction_end,strand,selected_event=exportToGff(junction1,junction2,gene_feature_db[geneID],geneID)
		exported=True
	    except Exception: pass ### Not handling trans-splicing well
	if exported:
	    ### For each exon entry exported for this event, we need a gene entry with the boundary exon positions exported (returned from the above function)
	    try: writegene(chr,junction_start,junction_end,strand,selected_event) #('chrX', '38600355', '38606652', '+', 'ENSMUSG00000000355__E1.2-E5.1')
	    except Exception: pass

    ### Export coordinate information for the entire gene
    for geneID in gene_feature_db:
	feature_db = gene_feature_db[geneID]
	for event in feature_db:
	    if 'E1.1-E100.1' in event:
		pd = feature_db[event]
		chr = pd.Chr()
		start,stop = pd.Start(), pd.End()
		if start>stop: strand = '-'
		else: strand = '+'
		try:
		    chr,junction_start,junction_end,strand,selected_event=exportToGff('null',event,gene_feature_db[geneID],geneID)
		    writegene(chr,junction_start,junction_end,strand,selected_event) #('chrX', '38600355', '38606652', '+', 'ENSMUSG00000000355__E1.2-E5.1')
		except Exception:
		    pass

    time.sleep(2)
    gff_export_obj.close()
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
    """ Begin building strand gff and index files for SashimiPlot based on the AltAnalyze database
    exon, junction and gene annotations """
    
    global gff_export_obj
    try:
	### When fl is strand dataset information object
	countsFileDir = fl.CountsFile() ### Counts file containing exon and junction positions
	root_dir = fl.RootDir() ### Root folder location
    except Exception:
	### STRAND proper object may not be supplied with this information. Use the root directory alone to infer these
	root_dir = fl
	search_dir = root_dir+'/ExpressionInput'
	files = unique.read_directory(search_dir) ### all files in ExpressionInput
	for file in files:
	    if 'counts.' in file and 'steady-state.txt' not in file:
		countsFileDir = search_dir+'/'+file ### counts file with exon positions

    PSIFileDir = root_dir+'/AltResults/AlternativeOutput/'+species+'_RNASeq_top_alt_junctions-PSI.txt'
    OutputDir=findParentDir(PSIFileDir)
    output=OutputDir+"events_sashimi.gff"
    gff_export_obj=open(output,'w')
    
    ### Sometimes only junctions are in the count file so create strand new file with detected junctions and all exons
    ### This information and associated featrues is extracted from the counts file
    featuresEvaluated = extractFeatures(species,countsFileDir)

    ### Compile and export the coordinates to gff format and index these coordinates for fast retreival by Miso
    Indexing(featuresEvaluated,PSIFileDir,output)
  
def extractFeatures(species,countsFileDir):
    import export
    ExonsPresent=False
    lastgene = None
    lastend = None
    genes_detected={}
    count=0
    first_last_exons = {} ### Make strand fake junction comprised of the first and last exon
    if 'counts.' in countsFileDir:
	### The feature_file contains only ExonID or Gene IDs and associated coordinates
	feature_file = string.replace(countsFileDir,'counts.','features.')
	fe = export.ExportFile(feature_file)
	firstLine = True
	for line in open(countsFileDir,'rU').xreadlines():
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
			first_last_exons[gene,strand] = [(chr,start)]
		    if len(genes_detected)==0 or lastgene in genes_detected: ### restrict to detected genes
			try: first_last_exons[lastgene,laststrand].append(lastend)
			except Exception:
			    pass ### occurs for the first gene	
		if ExonsPresent == False:
		    fe.write(gene+':'+exon+'='+chr+':'+start+'-'+end+'\n')
		lastgene = gene; lastend = end; laststrand = strand
	if len(genes_detected)==0 or lastgene in genes_detected:
	    first_last_exons[lastgene,laststrand].append(lastend)
	
	### Add strand fake junction for the whole gene
	for (gene,strand) in first_last_exons:
	    (chr,start),end = first_last_exons[gene,strand]
	    if strand == '-':
		start,end = end,start # Need to encode strand in this annotation, do this by strand orienting the positions
	    fe.write(gene+':E1.1-E100.1'+'='+chr+':'+start+'-'+end+'\n')
	fe.close()
    return feature_file	### return the location of the exon and gene coordinates file

def obtainTopGeneResults():
    pass

if __name__ == '__main__':
    Species = 'Hs'
    root = '/Volumes/salomonis2-1/Ichi_data/bams/Insertion_analysis/'
    #root = '/Users/saljh8/Desktop/Grimes/GEC14074/'
    remoteIndexing(Species,root)
    sys.exit()
    #"""
    countsFileDir = os.path.abspath(os.path.expanduser(sys.argv[1]))
    PSIFileDir = os.path.abspath(os.path.expanduser(sys.argv[2]))
    OutputDir=findParentDir(PSIFileDir)
    output=OutputDir+"events_sashimi.gff"
    gff_export_obj=open(output,'w')
    Indexing(countsFileDir,PSIFileDir,output)
    #"""