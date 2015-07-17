#!/usr/bin/env python
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
from pyPdf import PdfFileReader, PdfFileWriter
import argparse
import math
import unique
count_sum_array=[];
dem=dict()
sample_read=dict()
gene_label=[]
gene_sym=dict()

def genelist(fname):
    fname = unique.filepath(fname)
    lis=[]
    for line in open(fname,'rU').xreadlines():
        line = cleanUpLine(line)
        t=string.split(line,'\t')
	gene=string.split(t[2],':')
	val=t[2]+' '+t[3]
	lis.append(val)
	
        if gene[0] in gene_sym:
		continue	
	else:
		gene_sym[gene[0]]=t[0]
	       
	
        #print t[0]
    return lis,gene_sym

def sample(fname):
    fname = unique.filepath(fname)
    head=0
    samplelis=[]
    for line in open(fname,'rU').xreadlines():
        line = cleanUpLine(line)
	if head ==0:
	    t=string.split(line,'\t')
	    #print t
	    for p in range(11,len(t)):
		samplelis.append(t[p])
	    head=1
        else:
	    break;
    return samplelis

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def update_plot_settings(bamdir,list1,list2,samp):
    export_pl=open(unique.filepath('Config/sashimi_plot_settings.txt'),'w')
    export_pl.write('[data]')
    export_pl.write('\n')
    export_pl.write('bam_prefix = '+bamdir+'\n')
    export_pl.write('bam_files =[')
  
    for i in range(len(list1)):
        g=samp[list1[i]].replace('.bed','.bam')
	print i
        if i==len(list1)-1 and len(list2)==0:
            export_pl.write('"'+g+'"]')
        else:
            export_pl.write('"'+g+'",')      
    for j in range(len(list2)):
	print j
        g=samp[list2[j]].replace('.bed','.bam')
        export_pl.write('"'+g+'"')
        if j==len(list2)-1:
            export_pl.write(']')
        else:
            export_pl.write(',')
    
	
    export_pl.write('\n')
    export_pl.write('[plotting]')
    export_pl.write('\n') 
    export_pl.write('fig_width = 7 \nfig_height = 7 \nintron_scale = 30 \nexon_scale = 4 \nlogged = False\n')
    export_pl.write('font_size = 6 \nbar_posteriors = False \nnyticks = 4 \nnxticks = 4 \n')
    export_pl.write('show_ylabel = False \nshow_xlabel = True \nshow_posteriors = False \nnumber_junctions = True \n')
    export_pl.write('resolution = .5 \nposterior_bins = 40 \ngene_posterior_ratio = 5 \n')
    export_pl.write('colors =[')
    for i in range(len(list1)):
        export_pl.write('"'+'red'+'"')
        if i==len(list1)-1 and len(list2)==0:
            export_pl.write(']')
        else:
            export_pl.write(',')
    for j in range(len(list2)):
        export_pl.write('"'+'blue'+'"')
        if j==len(list2)-1:
            export_pl.write(']')
        else:
            export_pl.write(',')
    export_pl.write('\n')       
    export_pl.write('coverages =[')
    for i in range(len(list1)):
        
        e=sample_read[samp[list1[i]]]
        export_pl.write(str(int(e)))
        if i==len(list1)-1 and len(list2)==0:
            export_pl.write(']')
        else:
            export_pl.write(',')
    for j in range(len(list2)):
        e=sample_read[samp[list2[j]]]
        export_pl.write(str(int(e)))
        if j==len(list2)-1:
            export_pl.write(']')
        else:
            export_pl.write(',')
    export_pl.write('\n')
    export_pl.write('bar_color = "b" \nbf_thresholds = [0, 1, 2, 5, 10, 20]')
    export_pl.close()
        
def sashmi_plot_list(bamdir,fname,gene_label,lines,samp,gene_sym):
  line=fname.readlines()
  fname.close()
  for li in line:
    if ":U" in li or "-U" in li:
	
	continue
    else:
    
     li=cleanUpLine(li)
     print li
     
    #dem[0]=['ENSG00000132424:I10.1 ENSG00000132424:E10.1-E11.1','ENSG00000146147:E10.3-E11.1 ENSG00000146147:E9.3-E15.1']
     de=string.split(li,'\t')
     dem[0]=de
     print dem[0]
     for key in dem:
      for i in range(len(dem[key])):
	list1=[]
        list2=[]
        try:
	    k=gene_label.index(dem[key][i])
	    flag=1
	    lt=cleanUpLine(lines[k])
	
	    t=string.split(lt,'\t')
	    print t
	
	    t=t[11:]
	    print t
	
    
            #list3=[]
            #ind=[]
        
	    for x in range(len(t)):
                print x,t[x]
                if(t[x]!=''):
                    if float(t[x]) < 0.8:
                        list1.append(x)
			#print x
                        #print 'list1:'+str(x)
                    else:
                        list2.append(x)
			#print x
                       # print str(x)
                 
                else:
                    continue
	
	    if len(list1)>5:
                list1=list1[1:5]
	    if len(list2)>5:
                list2=list2[1:5]
	    print len(list1),len(list2)
	except Exception:
	    
	    for ij in range(len(samp)):
		list1.append(ij)
	
	update_plot_settings(bamdir,list1,list2,samp)
	
	a=string.split(dem[key][i]," ")
	if '-' in a[1]:
                
                ch1=a[1]
                f=string.split(a[0],':')
	else:
        	ch1=a[0]
                f=string.split(a[1],':')
        event=findParentDir(inputpsi)
        event=event+"trial_index/"
        setting =unique.filepath("Config/sashimi_plot_settings.txt")
        name=ch1
	#outputdir=findParentDir(inputpsi)+"sashimiplots"
	try: os.makedirs(outputdir)
	except Exception: pass
	
    try:
	ssp.plot_event(ch1,event,setting,outputdir)
    except Exception:
	 print "error2"
	 continue
  #outputdir=findParentDir(inputpsi)+"sashimiplots" 
  for filename in os.listdir(outputdir):
    newname=string.split(filename,'/')
    #print newname[0]
    if newname[0] in gene_sym:
	new_path = gene_sym[newname[0]]+'-'+filename
	#new_path = string.replace()
	os.rename(filename,new_path)
    else:
	continue

  
def findParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1
    return filename[:x]      

def Sashimiplottting(bamdir,countsin,inputpsi,genelis):
    inputpsi = unique.filepath(inputpsi)
    text_file = open(inputpsi,'rU')
    lines = text_file.readlines()
   
    text_file.close()
    samp=sample(inputpsi)
    gene_label,gene_sym=genelist(inputpsi)

    header=True
    junction_max=[]
    countsin = unique.filepath(countsin)
    for line in open(countsin,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if header:
            samples = t[1:]
            header=False
            exon_sum_array=[0]*len(samples)
            count_sum_array=[0]*len(samples)
        else:
            values = map(float,t[1:])
            count_sum_array = [sum(value) for value in zip(*[count_sum_array,values])]

        for i in range(len(samp)):
	   sample_read[samp[i]]=count_sum_array[i]
          #print samp[i],sample_read[samp[i]]

    genelis = unique.filepath(genelis)
    gene_file=open(genelis,'rU')

    sashmi_plot_list(bamdir,gene_file,gene_label,lines,samp,gene_sym)

def remoteSashimiPlot(species,fl,bamdir,genelis):
    global inputpsi
    global outputdir
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

    #outputdir=findParentDir(inputpsi)+"sashimiplots"
    outputdir = root_dir+'/ExonPlots'
    outputdir = root_dir+'/SashimiPlots'
    try: os.mkdir(unique.filepath(outputdir))
    except Exception: pass
    Sashimiplottting(bamdir,countinp,inputpsi,genelis)

    gene_label,gene_sym=genelist(inputpsi)
    for filename in os.listdir(outputdir):
	if '.pdf' in filename:
	    newname=string.split(filename,':')
	    if newname[0] in gene_sym:
		new_filename = str(filename)
		if ':' in filename:
		    new_filename = string.split(filename,':')[1]
		elif '\\' in filename:
		    new_filename = string.split(filename,'\\')[1]
		elif '/' in filename:
		    new_filename = string.split(filename,'/')[1]
	        nnname=gene_sym[newname[0]]+'-SashimiPlot_'+new_filename
		os.rename(os.path.join(outputdir, filename), os.path.join(outputdir,nnname))
	    else:
		continue

if __name__ == '__main__':
    root_dir = '/Users/saljh8/Desktop/Grimes/GEC14074/'
    genelis = '/Users/saljh8/Desktop/Grimes/GEC14074/AltResults/AlternativeOutput/test.txt'
    remoteSashimiPlot('Mm',root_dir,root_dir,genelis)
    sys.exit()
