#!/usr/bin/env python

import numpy as np
import pylab as pl
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

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
import anydbm
import dbhash
count_sum_array_db={}
sampleReadDepth={}

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
      
def importSplicingEventsToVisualize(eventsToVisualizeFilename):
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
	else:
	   splicing_events.append(t[0])
    splicing_events = unique.unique(splicing_events)
    return splicing_events,expandedSearch

def sashmi_plot_list(bamdir,eventsToVisualizeFilename,PSIFilename,events=None):
    try:
	import gene_associations
	gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
	from import_scripts import OBO_import; symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)	
    except Exception:
	symbol_to_gene={}


    if events==None:
        splicing_events,expandedSearch = importSplicingEventsToVisualize(eventsToVisualizeFilename)
    else:
        ### Replace any ":" from the input events
        #for i in range(len(events)): events[i] = string.replace(events[i],':','__')
        expandedSearch = True
        
        for i in range(len(events)):
            gene = string.split(events[i],'__')[0]
            if gene in gene_to_symbol:
                symbol = gene_to_symbol[gene][0]
            elif 'ENS' not in gene or 'G0000' in gene:
                if gene in symbol_to_gene:
                    ensID = symbol_to_gene[gene][0]
                    symbol = gene
                    events[i] = ensID ### translate this ID to an Ensembl gene ID for propper SashimiPlot lookup
        splicing_events = events ### optionally get from supplied variable

    if len(splicing_events)==0:
        print eventsToVisualizeFilename
        forceNoCompatibleEventsInFile
    
    print 'Exporting plots',
    
    ### Determine Groups for Coloring
    groups_file = 'None'
    dir_list = unique.read_directory(root_dir+'/ExpressionInput')

    for file in dir_list:
         if 'groups.' in file:
            groups_file = root_dir+'/ExpressionInput/'+file

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

    processed_events = formatAndSubmitSplicingEventsToSashimiPlot(PSIFilename, bamdir, splicing_events, sample_group_db, groups, False)
    mopup_events = getMopUpEvents(splicing_events, processed_events)

    ### Do the same for supplied gene queries or junctions that didn't map above using the gene expression values as a guide
    #print len(splicing_events),len(processed_events),len(mopup_events)
    processed_events = formatAndSubmitSplicingEventsToSashimiPlot(steady_state_exp_file,bamdir,mopup_events,sample_group_db,groups,expandedSearch)
    if len(processed_events)>0:
        mopup_events = getMopUpEvents(mopup_events, processed_events)
        processed_events = formatAndSubmitSplicingEventsToSashimiPlot(PSIFilename, bamdir, mopup_events, sample_group_db, groups, True)
    return gene_to_symbol

def getMopUpEvents(splicing_events, processed_events):
    mopup_events = []
    for event in splicing_events:
        add = True
        if event in processed_events:
            add = False
        if ' ' in event:
            try:
                j1, j2 = string.split(event, ' ')
                if j1 in processed_events:
                    add = False
                if j2 in processed_events:
                    add = False
            except Exception:
                pass
        if add:
            mopup_events.append(event)
    return mopup_events

def reorderEvents(events):
    splicing_events = events
    index = 0
    for e in events:
        j1o, j2o = string.split(e, ' ')
        gene, j1 = string.split(j1o, ':')
        gene, j2 = string.split(j2o, ':')
        if '-' in j1 and '-' in j2:
            j1a, j1b = string.split(j1, '-')
            j2a, j2b = string.split(j2, '-')
            j1a_block, j1a_region = string.split(j1a[1:], '.')
            j2a_block, j2a_region = string.split(j2a[1:], '.')
            j1b_block, j1b_region = string.split(j1b[1:], '.')
            j2b_block, j2b_region = string.split(j2b[1:], '.')
            if int(j1b_block) < int(j2b_block) and int(j1a_block) < int(j2a_block):
                pass ### Occurs for complex cassette exon splicing events but matches SashimiIndex's selection for exclusion
            elif int(j1b_block) > int(j2b_block):
                new_e = j2o + ' ' + j1o
                splicing_events[index] = new_e
            elif int(j1a_block) < int(j2a_block):
                new_e = j2o + ' ' + j1o
                splicing_events[index] = new_e
            elif int(j1b_region) > int(j2b_region):
                new_e = j2o + ' ' + j1o
                splicing_events[index] = new_e
            elif int(j1a_region) < int(j2a_region):
                new_e = j2o + ' ' + j1o
                splicing_events[index] = new_e

        index += 1
    return splicing_events

def formatAndSubmitSplicingEventsToSashimiPlot(filename,bamdir,splicing_events,sample_group_db,groups,expandedSearch):
    ### Begin exporting parameters and events for SashimiPlot visualization
    firstLine = True
    setting = unique.filepath("Config/sashimi_plot_settings.txt")
    psi_parent_dir=findParentDir(filename)
    if 'PSI' not in filename:
        index_dir=string.split(psi_parent_dir,'ExpressionInput')[0]+"AltResults/AlternativeOutput/sashimi_index/"
    else:
        index_dir=psi_parent_dir+"sashimi_index/"

    spliced_junctions=[] ### Alternatively, compare to just one of the junctions
    for splicing_event in splicing_events:
        try:
            j1,j2 = string.split(splicing_event,' ')
            spliced_junctions.append(j1)
            spliced_junctions.append(j2)    
        except Exception:
            spliced_junctions.append(splicing_event) ### single gene ID or junction

    if 'PSI' not in filename:
        splicing_events_db = {}
        for event in splicing_events:
            event = string.replace(event,':','__')
            if ' ' in event:
                event = string.split(event,' ')[-1]
            gene = string.split(event,"__")[0]
            try: splicing_events_db[gene].append(event)
            except Exception: splicing_events_db[gene] = [event]
        splicing_events = splicing_events_db
    
    import collections
    analyzed_junctions=[]
    processed_events=[]
    #restrictToTheseGroups=['WT','R636S_Het','R636S_Homo'] #Meg HSCP-1 , Myelocyte Mono
    restrictToTheseGroups = None
    for line in open(filename,'rU').xreadlines():
        line = cleanUpLine(line)
        t = string.split(line,'\t')
        if firstLine:
            if 'PSI' in filename:
                sampleIndexBegin = 11
                sample_headers = t[sampleIndexBegin:]
            else:
                sampleIndexBegin = 1
                sample_headers = t[sampleIndexBegin:]
                if '.bed' not in sample_headers[0]: ### Add .bed if removed manually
                    sample_headers = map(lambda s: s+'.bed',sample_headers)
            index=0
            sample_group_index={}
            for s in sample_headers:
                group = sample_group_db[s]
                sample_group_index[index]=group
                try: sampleReadDepth[index]=count_sum_array_db[s]
                except Exception: sampleReadDepth[index]=count_sum_array_db[s]
                index+=1
            firstLine = False
        else:
            if 'PSI' in filename:
                splicing_event = val=t[2]+' '+t[3]
                j1=t[2]
                j2=t[3]
                if t[2] in analyzed_junctions and t[3] in analyzed_junctions:
                    continue
            else:
                splicing_event = t[0] ### The gene ID
                j1 = t[0]
                j2 = t[0]
            if ":U" in splicing_event or "-U" in splicing_event:
                continue
            else:
                ### First check to see if the full splicing event matches the entry
                ### If not (and not a PSI regulation hits list), look for an individual junction match
                if splicing_event in splicing_events or (expandedSearch and (j1 in spliced_junctions or j2 in spliced_junctions)):
                    if splicing_event in processed_events: continue
                    if j2 in processed_events: continue
                    if j1 in processed_events: continue
                    processed_events.append(splicing_event)
                    processed_events.append(j1)
                    processed_events.append(j2)
                    #print processed_events, splicing_event
                    if 'PSI' in filename:
                        geneID = string.split(t[2],':')[0]
                        symbol = t[0]
                        analyzed_junctions.append(t[2])
                        analyzed_junctions.append(t[3])
                    else: ### For exp.dataset-steady-state.txt files
                        geneID = splicing_event
                        events = splicing_events[geneID]
                    index=0

                    import collections
                    initial_group_psi_values={}
                    try: group_psi_values = collections.OrderedDict()
                    except Exception:
                        try:
                            import ordereddict
                            group_psi_values = ordereddict.OrderedDict()
                        except Exception:
                            group_psi_values={}			
                    for i in t[sampleIndexBegin:]: ### Value PSI range in the input file
                        try: group = sample_group_index[index]
                        except Exception: group=None
                        try:
                            try: initial_group_psi_values[group].append([float(i),index])
                            except Exception: initial_group_psi_values[group] = [[float(i),index]]
                        except Exception:
                            #print traceback.format_exc();sys.exit()
                            pass ### Ignore the NULL values
                        index+=1
                    if restrictToTheseGroups !=None: ### Exclude unwanted groups
                        initial_group_psi_values2={}
			groups2 = collections.OrderedDict()
                        for group in groups:
			    if group in initial_group_psi_values:
				if group in restrictToTheseGroups:
				    initial_group_psi_values2[group]=initial_group_psi_values[group]
				    groups2[group]=[]
                        initial_group_psi_values = initial_group_psi_values2
			groups = groups2
                    ### limit the number of events reported and sort based on the PSI values in each group
                    if 'None' in groups and len(groups)==1:
                        initial_group_psi_values['None'].sort()
                        group_size = len(initial_group_psi_values['None'])/2
                        filtered_group_index1 = map(lambda x: x[1], initial_group_psi_values['None'][:group_size])
                        filtered_group_index2 = map(lambda x: x[1], initial_group_psi_values['None'][group_size:])
                        group_psi_values['low']=filtered_group_index1
                        group_psi_values['high']=filtered_group_index2
                    else:
                        gn=0
                        for group in groups:
			    #print group
			    gn+=1
			    #if gn>4: break
                            if group in initial_group_psi_values:
                                initial_group_psi_values[group].sort()
				if len(groups)>7:
				    filtered_group_indexes = map(lambda x: x[1], initial_group_psi_values[group][:1])
				elif len(groups)>5:
				    filtered_group_indexes = map(lambda x: x[1], initial_group_psi_values[group][:2])
				elif len(groups)>3:
				    filtered_group_indexes = map(lambda x: x[1], initial_group_psi_values[group][:4])
				else:
				    filtered_group_indexes = map(lambda x: x[1], initial_group_psi_values[group][:5])
                                group_psi_values[group]=filtered_group_indexes
                    try: update_plot_settings(bamdir,group_psi_values,sample_headers)
		    except Exception:
			print 'Cannot update the settings file. Likely permissions issue.'

		    try:
			reordered = reorderEvents([t[2] + ' ' + t[3]])
			reordered = string.split(reordered[0], ' ')
		    except Exception:
			reordered = [t[2] + ' ' + t[3]]
			reordered = string.split(reordered[0], ' ')
		    #print reordered
                    if 'PSI' in filename:
                        try: formatted_splice_event = string.replace(reordered[1], ':', '__')
                        except Exception: pass
                        ### Submit the query
                        try: ssp.plot_event(formatted_splice_event,index_dir,setting,outputdir); success = True
                        except Exception:
                            success = False
                            #print traceback.format_exc()
                    
                    else:
                        for event in events:
                            try:
                                ssp.plot_event(event,index_dir,setting,outputdir)
                                #print 'success' #formatted_splice_event='ENSMUSG00000000355__E4.1-E5.1'
                            except Exception: ### If it fails, output the gene-level plot
                                try: ssp.plot_event(geneID,index_dir,setting,outputdir); success = True
                                except Exception:
                                    success = False
                                    #print traceback.format_exc()
		    """
                    ### Second attempt
                    if 'PSI' in filename and success==False: ### Only relevant when parsing the junction pairs but not genes
                        try: formatted_splice_event=string.replace(reordered[0],':','__')
                        except Exception: pass
                        try: ssp.plot_event(formatted_splice_event,index_dir,setting,outputdir); # print 'success'
                        except Exception: pass
		    """
    return processed_events

def findParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1
    return filename[:x]      

def Sashimiplottting(bamdir,countsin,PSIFilename,eventsToVisualizeFilename,events=None):
    PSIFilename = unique.filepath(PSIFilename)
    header=True
    junction_max=[]
    countsin = unique.filepath(countsin)
    count_sum_array=[]
    count=0
    for line in open(countsin,'rU').xreadlines():
	data = cleanUpLine(line)
	t = string.split(data,'\t')
	if header:
	    samples = []
	    for s in t[1:]:
		if '.bed' not in s: s+='.bed'
		samples.append(s)
	    header=False
	    count_sum_array=[0]*len(samples)
	else:
	    values = map(float,t[1:])
	    count_sum_array = [sum(value) for value in zip(*[count_sum_array,values])]
	    count+=1
	    if count >30000 and 'salomonis' in bamdir: break

    index=0
    for sample in samples:
        count_sum_array_db[sample] = count_sum_array[index]
        index+=1

    if events==None:
        #print 'Preparing Sashimi-Input:',eventsToVisualizeFilename
        eventsToVisualizeFilename = unique.filepath(eventsToVisualizeFilename)

    gene_to_symbol=sashmi_plot_list(bamdir,eventsToVisualizeFilename,PSIFilename,events=events)
    return gene_to_symbol

def remoteSashimiPlot(Species,fl,bamdir,eventsToVisualizeFilename,events=None,show=False):
    global PSIFilename
    global outputdir
    global root_dir
    global steady_state_exp_file
    global species
    species = Species
    
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

    ### Export BAM file indexes
    from import_scripts import BAMtoJunctionBED
    try: BAMtoJunctionBED.exportIndexes(root_dir)
    except:
        print 'BAM file indexing failed...'
        print traceback.format_exc()
        
    PSIFilename = root_dir+'/AltResults/AlternativeOutput/'+species+'_RNASeq_top_alt_junctions-PSI.txt'
    
    import ExpressionBuilder
    dir_list = unique.read_directory(root_dir+'/ExpressionInput')
    for file in dir_list:
        if 'exp.' in file and 'steady-state' not in file:
            exp_file = root_dir+'/ExpressionInput/'+file
        elif 'exp.' in file and 'steady-state' in file:
            steady_state_exp_file = root_dir+'/ExpressionInput/'+file
    global sample_group_db
    sample_group_db = ExpressionBuilder.simplerGroupImport(exp_file)
    
    #outputdir=findParentDir(PSIFilename)+"sashimiplots"
    outputdir = root_dir+'/ExonPlots'
    outputdir = root_dir+'/SashimiPlots'
    try: os.mkdir(unique.filepath(outputdir))
    except Exception: pass
    
    if show:
        s = open(outputdir+'/show.txt','w')
        s.write('TRUE'); s.close()
    else:
        s = open(outputdir+'/show.txt','w')
        s.write('FALSE'); s.close()

    geneSymbol_db=Sashimiplottting(bamdir,countinp,PSIFilename,eventsToVisualizeFilename,events=events)
    for filename in os.listdir(outputdir):
        if '.pdf' in filename or '.png' in filename:
            fn = string.replace(filename,'.pdf','')
            fn = string.replace(fn,'.png','')
            newname=string.split(fn,'__')
            if newname[0] in geneSymbol_db:
                new_filename = str(filename)
                if '__' in filename:
                    new_filename = string.split(filename,'__')[1]
                elif '\\' in filename:
                    new_filename = string.split(filename,'\\')[1]
                elif '/' in filename:
                    new_filename = string.split(filename,'/')[1]
                nnname=geneSymbol_db[newname[0]][0]+'-SashimiPlot_'+new_filename
                try: os.rename(os.path.join(outputdir, filename), os.path.join(outputdir,nnname))
                except Exception:
                    if 'already exists' in traceback.format_exc():
                        ### File already exists, delete the new one
                        try: os.remove(os.path.join(outputdir,nnname))
                        except Exception: pass
                        ### Now right the new one
                        try: os.rename(os.path.join(outputdir, filename), os.path.join(outputdir,nnname))
                        except Exception: pass
                    pass
            else:
                continue
    print ''

def justConvertFilenames(species,outputdir):
    import gene_associations
    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    from import_scripts import OBO_import; symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
    
    for filename in os.listdir(outputdir):
        if '.pdf' in filename or '.png' in filename:
            fn = string.replace(filename,'.pdf','')
            fn = string.replace(fn,'.png','')
            newname=string.split(fn,'__')

            if newname[0] in gene_to_symbol:
                new_filename = str(filename)
                if '__' in filename:
                    new_filename = string.split(filename,'__')[1]
                elif '\\' in filename:
                    new_filename = string.split(filename,'\\')[1]
                elif '/' in filename:
                    new_filename = string.split(filename,'/')[1]
                nnname=gene_to_symbol[newname[0]][0]+'-SashimiPlot_'+new_filename
                try: os.rename(os.path.join(outputdir, filename), os.path.join(outputdir,nnname))
                except Exception: pass
            else:
                continue
            
if __name__ == '__main__':
    root_dir = '/Volumes/salomonis2/Bruce_conklin_data/Pig-RNASeq/bams/'
    events = ['ENSSSCG00000024564']
    events = None
    eventsToVisualizeFilename = None
    eventsToVisualizeFilename = '/Volumes/salomonis2/Bruce_conklin_data/Pig-RNASeq/events.txt'
    bamdir = root_dir
    remoteSashimiPlot('Ss', root_dir, bamdir, eventsToVisualizeFilename, events=events, show=False)
    sys.exit()
