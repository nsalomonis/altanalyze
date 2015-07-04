#!/usr/bin/env python
import numpy as np
import pylab as pl
import sys,string
import os
from misopy import index_gff as ifg
import subprocess
import multiprocessing
import time
samplelis=[]
samp=[]
sample_read=dict()
index_read=dict()
gene_label=[]
new_lis=[]
lis=[]


def indexdic(fname):
    head=0
    for line in open(fname,'rU').xreadlines():
#for k in range(len(a['AltAnalyze_ID'])):
     
     if head ==0:
	head=1
	continue
     else:
	a=string.split(line,'\t')

    #p=a['AltAnalyze_ID'][k]
	p=a[0]
	j=string.split(p,':')
    #print j[0]
	for i in range(len(j)):
	    if "ENS" in j[i]:
		if '-' in j[i]:
		    ji=string.split(j[i],'-')
		    jj=ji[1]
		else:
		    jj=j[i]
	    #print jj,'first check'
		if jj in index_read:
			index_read[jj].append(p)
		else:
			index_read[jj]=[p,]
    return index_read

def writegene(chrom,start_l,end_l,a,ch1):
    export_in.write(chrom+'\t'+'SE'+'\t'+'gene'+'\t'+start_l+'\t'+end_l+'\t'+'.'+'\t'+a+'\t'+'.'+'\t'+'ID='+ch1+';'+'Name='+ch1+';'+'\n')
    export_in.write(chrom+'\t'+'SE'+'\t'+'mRNA'+'\t'+start_l+'\t'+end_l+'\t'+'.'+'\t'+a+'\t'+'.'+'\t'+'ID='+ch1+'.A;'+'Parent='+ch1+';'+'\n')
    


def genelist(fname):
    for line in open(fname,'rU').xreadlines():
        line = line.rstrip(os.linesep)
        t=string.split(line,'\t')
	val=t[2]+' '+t[3]
        lis.append(val)
        #print t[0]
    return lis

def write_index(gene_n,new_lis):
    #print gene_n
    #print new_lis
    valid=0
    ch=string.split(gene_n," ")
    if '-' in ch[1]:
	ch1=ch[1]
    else:
	ch1=ch[0]
    gene_mo=string.split(ch1,':')
    gene_mo1=gene_mo[1].replace('E',"")
    gene_mo1=gene_mo1.replace('I',"")
    f=gene_mo1.split('-')
    c=f[0]
    
    k=len(gene_mo)
    gen_mo2=gene_mo[k-1].replace('E',"")
    gen_mo2=gen_mo2.replace('I',"")
    d=gen_mo2
    
    #print c,d
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
	             chro=string.split(c1[1],':')
	             chrom=chro[0]
		     lengt=string.split(chro[1],'-')
	             start_l=lengt[0]
		     
	             end_l=lengt[1]
		     if start_l>end_l:
			    a='-'
			    b=start_l
			    start_l=end_l
			    end_l=b
		     else:
			    a='+'
		     
		     
                     count=1
		     break
		    
		     
       else:
	    continue
    
    if count==1:
      #print start_l,end_l
      for j in range(len(new_lis)):
	q=string.split(new_lis[j],"=")
	#print q
	if '-' in q[0]:
	    continue
	else:
	    cho=string.split(q[0],':')
	    
	    ge=cho[1].replace('E',"")
	    #print ge
	    if ('I' not in ge) and (';' not in ge) and ('-' not in ge):
	     
             #if (ge >=c and ge <= d):
		#print g
                chro=string.split(q[1],':')
		chrom=chro[0]
		lengt=string.split(chro[1],'-')
		start_e=lengt[0]
		end_e=lengt[1]
		valid=0
	##sp=string.split(lis1[i],':')
	#if sp[1]==gene_mo[1]:
	        if(start_e >= start_l and end_e <= end_l):
		    valid=1 
	        elif(start_e< start_l and end_e == start_l):
		    start_l=start_e
		    #if start_e==end_e:
			#start_e=int(end_e)-10
		    
		    valid=1
		elif(end_e>end_l and start_e ==end_l):
		    end_l=end_e
		    #if start_e==end_e:
			#end_e=int(start_l)+10
		    valid=1
		if valid==1:
		       #print chrom,start_l,end_l,a,ch1
		       export_in.write(chrom+'\t'+'SE'+'\t'+'exon'+'\t'+str(start_e)+'\t'+str(end_e)+'\t'+'.'+'\t'+a+'\t'+'.'+'\t'+'ID='+q[0]+';'+'Parent='+ch1+'.A;'+'\n') 
	    else:
		continue
		#export_in.write(chrom+'\t'+'SE'+'\t'+'exon'+'\t'+start_l+'\t'+end_l+'\t'+'.'+'\t'+a+'\t'+'.'+'\t'+'ID='+q[0]+'Parent='+ch[1]+'.A;'+'\n')
    #try:
   
    return chrom,start_l,end_l,a,ch1
    

def Indexing(counts,inputpsi,output):
    
    head=0
    
    index_read=indexdic(counts)
    #print len(index_read)
#for k in range(len(a['AltAnalyze_ID'])):
 
    gene_label=genelist(inputpsi)
    for i in range(len(gene_label)):
  #if 'ENSMUSG00000066007:E5.1-ENSMUSG00000063245:E3.1' in gene_label[i]:
       try:
        if head==0:
	    head=1
	    continue
        p=string.split(gene_label[i],":")
        #print gene_label[i]
     
        if p[0] in index_read:
	    chrom,start_l,end_l,a,ch1=write_index(gene_label[i],index_read[p[0]])
	for k in range(len(p)):
    	#print p[k]
		if "ENS" in p[k]:
		    if '-' in p[k]:
			ji=string.split(p[k],'-')
			jj=ji[1]
		    else:
			jj=p[k]
	    
		    if jj != p[0]:
		        if jj in index_read:

		           chrom,start_l,end_l,a,ch1= write_index(gene_label[i],index_read[jj])
  

        
	writegene(chrom,start_l,end_l,a,ch1)
       except Exception:
		continue
    
    time.sleep(2)
    export_in.close()
    newout=findParentDir(output)+'trial_index/'
    try:
	ifg.index_gff(output,newout)
    except Exception:
	print('error in indexing')


  
def findParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1
    return filename[:x]      
        

if __name__ == '__main__':
  countinp = os.path.abspath(os.path.expanduser(sys.argv[1]))
  inputpsi = os.path.abspath(os.path.expanduser(sys.argv[2]))
  outputdir=findParentDir(inputpsi)
  output=outputdir+"events_trial.gff"
  export_in=open(output,'w')
  Indexing(countinp,inputpsi,output)

