#!/usr/bin/env python

import numpy as np
import pylab as pl
import sys,string
import os
import misopy
from misopy.sashimi_plot import sashimi_plot as ssp
import subprocess
import multiprocessing
import time
import subprocess
import random
import argparse
import math
import unique
import traceback
count_sum_array_db={}
sampleReadDepth={}
geneSymbol_db={}

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def update_plot_settings(bamdir,group_psi_values,sample_headers):
    ### This functions writes out the sample orders, colors and sequence coverage for each BAM files for SashimiPlot
    bams=[]
    sample_colors=[]
    sample_coverage=[]
    colors = ['red','blue','green','grey','orange','purple','yellow','peach','pink','violet','magenta','navy']
    colors = colors*300
    color_index=0

    for group in group_psi_values:
	for index in group_psi_values[group]:
	    g=sample_headers[index].replace('.bed','.bam')
	    bams.append('"'+g+'"')
	    sample_colors.append('"'+colors[color_index]+'"')
	    sample_coverage.append(str(int(sampleReadDepth[index])))
	color_index+=1 ### reset for the new group
    bams = string.join(bams,',')
    sample_colors = string.join(sample_colors,',')
    sample_coverage = string.join(sample_coverage,',')
    
    export_pl=open(unique.filepath('Config/sashimi_plot_settings.txt'),'w')
    export_pl.write('[data]\n')
    export_pl.write('bam_prefix = '+bamdir+'\n')
    export_pl.write('bam_files =['+bams+']\n')

    export_pl.write('\n')
    export_pl.write('[plotting]')
    export_pl.write('\n') 
    export_pl.write('fig_width = 7 \nfig_height = 7 \nintron_scale = 30 \nexon_scale = 4 \nlogged = False\n')
    export_pl.write('font_size = 6 \nbar_posteriors = False \nnyticks = 4 \nnxticks = 4 \n')
    export_pl.write('show_ylabel = False \nshow_xlabel = True \nshow_posteriors = False \nnumber_junctions = True \n')
    export_pl.write('resolution = .5 \nposterior_bins = 40 \ngene_posterior_ratio = 5 \n')
    export_pl.write('colors =['+sample_colors+']\n')
    export_pl.write('coverages =['+sample_coverage+']\n')
    
    export_pl.write('bar_color = "b" \nbf_thresholds = [0, 1, 2, 5, 10, 20]')
    export_pl.close()
        
def sashmi_plot_list(bamdir,eventsToVisualizeFilename,PSIFilename):
    splicing_events=[]
    
    ### Import the splicing events to visualize from an external text file (multiple formats supported)
    type = None
    expandedSearch = False
    firstLine = True
    for line in open(eventsToVisualizeFilename,'rU').xreadlines():
	line = cleanUpLine(line)
	t = string.split(line,'\t')
	if firstLine:
	    if 'junctionID-1' in t:
		j1i = t.index('junctionID-1')
		j2i = t.index('junctionID-2')
		type='ASPIRE'
		expandedSearch = True
	    if 'ANOVA' in t:
		type='PSI'
	    elif 'independent confirmation' in t:
		type='confirmed'
		expandedSearch = True
	    elif 'ANOVA' in eventsToVisualizeFilename:
		type = 'ANOVA'
	    firstLine=False
	if '|' in t[0]:
	    type = 'ANOVA'
	if ' ' in t[0] and ':' in t[0]:
	    splicing_events.append(t[0])
	elif type=='ASPIRE':
	    splicing_events.append(t[j1i] +' '+ t[j2i])
	    splicing_events.append(t[j2i] +' '+ t[j1i])
	elif type=='ANOVA':
	    try:
		a,b = string.split(t[0],'|')
		a = string.split(a,':')
		a = string.join(a[1:],':')
		splicing_events.append(a +' '+ b)
		splicing_events.append(b +' '+ a)
	    except Exception: pass
	elif type=='PSI':
	    try:
		j1,j2 = string.split(t[0],'|')
		a,b,c = string.split(j1,':')
		j1 = b+':'+c
		splicing_events.append(j1 +' '+ j2)
		splicing_events.append(j2 +' '+ j1)
	    except Exception:
		#print traceback.format_exc();sys.exit()
		pass
	elif type=='confirmed':
	    try:
		event_pair1 = string.split(t[1],'|')[0]
		a,b,c,d = string.split(event_pair1,'-')
		splicing_events.append(a+'-'+b +' '+ c+'-'+d)
		splicing_events.append(c+'-'+d +' '+ a+'-'+b)
	    except Exception: pass

    spliced_junctions=[] ### Alternatively, compare to just one of the junctions
    for splicing_event in splicing_events:
	j1,j2 = string.split(splicing_event,' ')
	spliced_junctions.append(j1)
	spliced_junctions.append(j2)
    
    if len(splicing_events)==0:
	forceNoCompatibleEventsInFile
    
    print 'Exporting plots',
    
    groups_file = 'None'
    dir_list = unique.read_directory(root_dir+'/ExpressionInput')
    for file in dir_list:
         if 'groups.' in file:
            groups_file = root_dir+'ExpressionInput/'+file
    if groups_file != None:
	try:
	    import ExpressionBuilder
            sample_group_db = ExpressionBuilder.simplerGroupImport(groups_file)
	    groups=[]
	    for sample in sample_group_db:
		if sample_group_db[sample] not in groups:
		    groups.append(sample_group_db[sample]) ### create an ordered list of unique group
	except Exception:
	    groups = ['None']
            #print traceback.format_exc()
            pass
	    
    firstLine = True
    setting = unique.filepath("Config/sashimi_plot_settings.txt")
    psi_parent_dir=findParentDir(PSIFilename)
    index_dir=psi_parent_dir+"trial_index/"

    import collections
    analyzed_events=[]
    for line in open(PSIFilename,'rU').xreadlines():
	line = cleanUpLine(line)
	t = string.split(line,'\t')
	if firstLine:
	    sample_headers = t[11:]
	    index=0
	    sample_group_index={}
	    for s in sample_headers:
		group = sample_group_db[s]
		sample_group_index[index]=group
		sampleReadDepth[index]=count_sum_array_db[s]
		index+=1
	    firstLine = False
	else:
	    splicing_event = val=t[2]+' '+t[3]
	    if t[2] in analyzed_events and t[3] in analyzed_events:
		continue
	    #if 'ENSG00000095794' not in splicing_event: continue
	    #else: print splicing_event
	    if ":U" in splicing_event or "-U" in splicing_event:
		continue
	    else:
		### First check to see if the full splicing event matches the entry
		### If not (and not a PSI regulation hits list), look for an individual junction match
		if splicing_event in splicing_events or (expandedSearch and (t[2] in spliced_junctions or t[3] in spliced_junctions)):
		    geneID = string.split(t[2],':')[0]
		    symbol = t[0]
		    geneSymbol_db[geneID]=symbol
		    index=0
		    analyzed_events.append(t[2])
		    analyzed_events.append(t[3])
		    import collections
		    initial_group_psi_values={}
		    try: group_psi_values = collections.OrderedDict()
		    except Exception:
			try:
			    import ordereddict
			    group_psi_values = ordereddict.OrderedDict()
			except Exception:
			    group_psi_values={}			
		    for i in t[11:]: ### Value PSI range in the input file
			try: group = sample_group_index[index]
			except Exception: group=None
			try:
			    try: initial_group_psi_values[group].append([float(i),index])
			    except Exception: initial_group_psi_values[group] = [[float(i),index]]
			except Exception:
			    #print traceback.format_exc();sys.exit()
			    pass ### Ignore the NULL values
			index+=1
		    ### limit the number of events reported and sort based on the PSI values in each group
		    if 'None' in groups and len(groups)==1:
			initial_group_psi_values['None'].sort()
			group_size = len(initial_group_psi_values['None'])/2
			filtered_group_index1 = map(lambda x: x[1], initial_group_psi_values['None'][:group_size])
			filtered_group_index2 = map(lambda x: x[1], initial_group_psi_values['None'][group_size:])
			group_psi_values['low']=filtered_group_index1
			group_psi_values['high']=filtered_group_index2
		    else:
			for group in groups:
			    if group in initial_group_psi_values:
				initial_group_psi_values[group].sort()
				filtered_group_indexes = map(lambda x: x[1], initial_group_psi_values[group][:5])
				group_psi_values[group]=filtered_group_indexes
		    update_plot_settings(bamdir,group_psi_values,sample_headers)
	
		    try: formatted_splice_event=string.replace(t[3],':','__')
		    except Exception: pass

		    try: os.makedirs(outputdir)
		    except Exception: pass
		    
		    #print '********',[formatted_splice_event],[index_dir],outputdir
		    #print symbol,formatted_splice_event,
		    try:
			ssp.plot_event(formatted_splice_event,index_dir,setting,outputdir)
			#print 'success'
		    except Exception:
			pass
			#print '^^^^^^^^^^^^',[formatted_splice_event],[index_dir],outputdir;sys.exit()
			#print traceback.format_exc()
			#print 'failed'

		    try: formatted_splice_event=string.replace(t[2],':','__')
		    except Exception: pass
		    #print symbol,formatted_splice_event,
		    try:
			ssp.plot_event(formatted_splice_event,index_dir,setting,outputdir)
			#print 'success'
		    except Exception:
			pass
			#print '^^^^^^^^^^^^',[formatted_splice_event],[index_dir],outputdir;sys.exit()
			#print traceback.format_exc()
			#print 'failed'
		    
    for filename in os.listdir(outputdir):
	newname=string.split(filename,'/')
	#print newname[0]
	if newname[0] in geneSymbol_db:
	    new_path = geneSymbol_db[newname[0]]+'-'+filename
	    #new_path = string.replace()
	    os.rename(filename,new_path)
	else:
	    continue

def findParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1
    return filename[:x]      

def Sashimiplottting(bamdir,countsin,PSIFilename,eventsToVisualizeFilename):
    PSIFilename = unique.filepath(PSIFilename)

    header=True
    junction_max=[]
    countsin = unique.filepath(countsin)
    count_sum_array=[]
    for line in open(countsin,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if header:
            samples = t[1:]
            header=False
            count_sum_array=[0]*len(samples)
        else:
            values = map(float,t[1:])
            count_sum_array = [sum(value) for value in zip(*[count_sum_array,values])]

    index=0
    for sample in samples:
	count_sum_array_db[sample] = count_sum_array[index]
	index+=1

    eventsToVisualizeFilename = unique.filepath(eventsToVisualizeFilename)

    sashmi_plot_list(bamdir,eventsToVisualizeFilename,PSIFilename)

def remoteSashimiPlot(species,fl,bamdir,eventsToVisualizeFilename):
    global PSIFilename
    global outputdir
    global root_dir
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
    
    PSIFilename = root_dir+'/AltResults/AlternativeOutput/'+species+'_RNASeq_top_alt_junctions-PSI.txt'
    
    import ExpressionBuilder
    dir_list = unique.read_directory(root_dir+'/ExpressionInput')
    for file in dir_list:
        if 'exp.' in file and 'steady-state' not in file:
            exp_file = root_dir+'/ExpressionInput/'+file
    global sample_group_db
    sample_group_db = ExpressionBuilder.simplerGroupImport(exp_file)
    
    #outputdir=findParentDir(PSIFilename)+"sashimiplots"
    outputdir = root_dir+'/ExonPlots'
    outputdir = root_dir+'/SashimiPlots'
    try: os.mkdir(unique.filepath(outputdir))
    except Exception: pass

    Sashimiplottting(bamdir,countinp,PSIFilename,eventsToVisualizeFilename)

    for filename in os.listdir(outputdir):
	if '.pdf' in filename or '.png' in filename:
	    newname=string.split(filename,'__')
	    if newname[0] in geneSymbol_db:
		new_filename = str(filename)
		if '__' in filename:
		    new_filename = string.split(filename,'__')[1]
		elif '\\' in filename:
		    new_filename = string.split(filename,'\\')[1]
		elif '/' in filename:
		    new_filename = string.split(filename,'/')[1]
	        nnname=geneSymbol_db[newname[0]]+'-SashimiPlot_'+new_filename
		os.rename(os.path.join(outputdir, filename), os.path.join(outputdir,nnname))
	    else:
		continue

if __name__ == '__main__':
    root_dir = '/Volumes/SEQ-DATA/Jared/'
    eventsToVisualizeFilename = '/Volumes/SEQ-DATA/SingleCell-Churko/Filtered/Unsupervised-AllExons/NewVersion/AltResults/AlternativeOutput/exp.Hs_RNASeq_top_alt_junctions-PSI-clust-ANOVA.txt'
    #eventsToVisualizeFilename = '/Volumes/SEQ-DATA/SingleCell-Churko/Filtered/Unsupervised-AllExons/NewVersion/AltResults/Clustering/Combined-junction-exon-evidence.txt'
    #eventsToVisualizeFilename = '/Volumes/salomonis1/projects/Bex1-RIP/Input/AltAnalyze_new/AltResults/Clustering/top50/Combined-junction-exon-evidence.txt'
    #eventsToVisualizeFilename = '/Volumes/salomonis2/2015-08-05_TRAF6_KO_RNA_Seq_pair_end/bams/AltResults/AlternativeOutput/Mm_RNASeq_TRAF6 KO_vs_wt.ExpCutoff-5.0_average-ASPIRE-exon-inclusion-results.txt'
    bamdir = '/Volumes/SEQ-DATA/SingleCell-Churko/Filtered/Unsupervised-AllExons/NewVersion/'
    remoteSashimiPlot('Hs',root_dir,bamdir,eventsToVisualizeFilename)
    sys.exit()
