#!/usr/bin/env python

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import numpy as np
import pylab as pl
import os.path
from collections import defaultdict
from sklearn.cluster import KMeans

try:from stats_scripts import statistics
except Exception: import statistics
import random
import UI
import export; reload(export)
import re
from stats_scripts import fishers_exact_test
import traceback
import warnings
import math
import export

def FishersExactTest(r,n,R,N):
    z=0.0
    """
    N is the total number of genes measured (Ensembl linked from denom) (total number of ) (number of exons evaluated)
    R is the total number of genes meeting the criterion (Ensembl linked from input) (number of exonic/intronic regions overlaping with any CLIP peeks)
    n is the total number of genes in this specific MAPP (Ensembl denom in MAPP) (number of exonic/intronic regions associated with the SF)
    r is the number of genes meeting the criterion in this MAPP (Ensembl input in MAPP) (number of exonic/intronic regions with peeks overlapping with the SF)
    
    With these values, we must create a 2x2 contingency table for a Fisher's Exact Test
    that reports:
     
    +---+---+    a is the # of IDs in the term regulated
    | a | b |    b is the # of IDs in the term not-regulated 
    +---+---+    c is the # of IDs not-in-term and regulated
    | c | d |    d is the # of IDs not-in-term and not-regulated
    +---+---+
    
    If we know r=20, R=80, n=437 and N=14480
    +----+-----+    
    | 20 | 417 |  437   
    +----+-----+    
    | 65 |13978| 14043
    +----+-----+
      85  14395  14480
    """
    if (R-N) == 0: return 0
    elif r==0 and n == 0: return 0
    else:
        try:
            #try:
            z = (r - n*(R/N))/math.sqrt(n*(R/N)*(1-(R/N))*(1-((n-1)/(N-1))))
          
            #except ZeroDivisionError: print 'r,_n,R,N: ', r,_n,R,N;kill
        except Exception: print (r - n*(R/N)), n*(R/N),(1-(R/N)),(1-((n-1)/(N-1))),r,n,N,R;kill
    a = r; b = n-r; c=R-r; d=N-R-b
    table = [[int(a),int(b)], [int(c),int(d)]]

    """
    print a,b; print c,d
    import fishers_exact_test; table = [[a,b], [c,d]]
    ft = fishers_exact_test.FishersExactTest(table)
    print ft.probability_of_table(table); print ft.two_tail_p()
    print ft.right_tail_p(); print ft.left_tail_p()
    """
    
    try: ### Scipy version - cuts down rutime by ~1/3rd the time
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides import warnings
            oddsratio, pvalue = stats.fisher_exact(table)
       # print pvalue
        return pvalue,z
    except Exception:
        ft = fishers_exact_test.FishersExactTest(table)
       # print ft.two_tail_p()
        return ft.two_tail_p(),z
    

            
def header_file(fname, delimiter=None,Expand="no"):
    head=0
    header=[]
    newheader=[]
    with open(fname, 'rU') as fin:
        for line in fin:
            #print line
            line = line.rstrip(os.linesep)
            header=string.split(line,'\t')
        
            if Expand=="yes":
                if head==0:
                    
                    for i in range(1,len(header)):
                        iq=header[i]
                        #iq=string.split(header[i],".")[0]
                        newheader.append(iq)
                    head=1
                else:break
            else:
                if len(header)<3 or Expand=="no":
                    if header[0] not in newheader:
                        newheader.append(header[0]) 
    #print len(newheader)           
    return newheader

def Enrichment(Inputfile,mutdict,mutfile,Expand,header):
    import collections
    import mappfinder
    X=defaultdict(list)
    prev=""
    head=0
    group=defaultdict(list)
    enrichdict=defaultdict(float)
    mut=export.findFilename(mutfile)
    dire=export.findParentDir(Inputfile)
    output_dir = dire+'MutationEnrichment'
    export.createExportFolder(output_dir)

    exportnam=output_dir+'/Enrichment_Results.txt'
    export_enrich=open(exportnam,"w")
    exportnam=output_dir+'/Enrichment_tophits.txt'
    export_hit=open(exportnam,"w")
   
    export_enrich.write("Mutations"+"\t"+"Cluster"+"\t"+"r"+"\t"+"R"+"\t"+"n"+"\t"+"Sensitivity"+"\t"+"Specificity"+"\t"+"z-score"+"\t"+"Fisher exact test"+"\t"+"adjp value"+"\n")
    if Expand=="yes":
        header2=header_file(Inputfile,Expand="yes")
        
        for line in open(Inputfile,'rU').xreadlines():
            if head >0:
                line=line.rstrip('\r\n')
                q= string.split(line,'\t')
                for i in range(1,len(q)):
                    if q[i]==str(1):
                        #group[q[0]].append(header2[i-1])
                        group[header2[i-1]].append(q[0])
           
            else:
                head+=1
                continue
    else:
        for line in open(Inputfile,'rU').xreadlines():
            line=line.rstrip('\r\n')
            line=string.split(line,'\t')
            #for i in range(1,len(line)):
            group[line[2]].append(line[0])
   
    total_Scores={}
    for kiy in mutdict:
        if kiy =="MDP":
            print mutdict[kiy]
        groupdict={}
        remaining=[]
        remaining=list(set(header) - set(mutdict[kiy]))
        groupdict[1]=mutdict[kiy]
        groupdict[2]=remaining
       # export_enrich1.write(kiy)
        for key2 in group:
           
            
            r=float(len(list(set(group[key2])))-len(list(set(group[key2]) - set(mutdict[kiy]))))
            n=float(len(group[key2]))
            R=float(len(set(mutdict[kiy])))
            N=float(len(header))
        
            if r==0 or R==1.0:
                print kiy,key2,r,n,R,N
                pval=float(1)
                z=float(0)
                null_z = 0.000
                zsd = mappfinder.ZScoreData(key2,r,R,z,null_z,n)
                zsd.SetP(pval)
            else:
                try: z = Zscore(r,n,N,R)
                except : z = 0.0000
                ### Calculate a Z-score assuming zero matching entries
                try: null_z = Zscore(0,n,N,R)
                except Exception: null_z = 0.000
               
                
                try:
                    pval = mappfinder.FishersExactTest(r,n,R,N)
                    zsd = mappfinder.ZScoreData(key2,r,R,z,null_z,n)
                    zsd.SetP(pval)
                except Exception:
                    pval=1.0
                    zsd = mappfinder.ZScoreData(key2,r,R,z,null_z,n)
                    zsd.SetP(pval)
                    #pass
                
          
            if kiy in total_Scores:
                    signature_db = total_Scores[kiy]
                    signature_db[key2]=zsd ### Necessary format for the permutation function
            else:
                    signature_db={key2:zsd}
                    total_Scores[kiy] = signature_db
    sorted_results=[]
    mutlabels={}
    for kiy in total_Scores:
        
        signature_db = total_Scores[kiy]
        ### Updates the adjusted p-value instances
        mappfinder.adjustPermuteStats(signature_db)
        for signature in signature_db:
            zsd = signature_db[signature]
           
            results = [kiy,signature,zsd.Changed(),zsd.Measured(),zsd.InPathway(),str(float(zsd.PercentChanged())/100.0),str(float(float(zsd.Changed())/float(zsd.InPathway()))), zsd.ZScore(), zsd.PermuteP(), zsd.AdjP()] #string.join(zsd.AssociatedIDs(),'|')
            sorted_results.append([signature,float(zsd.PermuteP()),results])
    sorted_results.sort() ### Sort by p-value
    prev=""
    for (sig,p,values) in sorted_results:
        if sig!=prev:
            flag=True
            export_hit.write(string.join(values,'\t')+'\n')
        if flag:
            if (float(values[5])>=0.5 and float(values[6])>=0.5) or float(values[5])>=0.6 :
                mutlabels[values[1]]=values[0]
                flag=False
                export_hit.write(string.join(values,'\t')+'\n')
        export_enrich.write(string.join(values,'\t')+'\n')
        prev=sig
    if len(sorted_results)==0:
            export_enrich.write(string.join([splicing_factor,'NONE','NONE','NONE','NONE','NONE','NONE'],'\t')+'\n')
    export_enrich.close()
    #print mutlabels
    return mutlabels
            
        
def findsiggenepermut(mutfile):
    
    samplelist=[]
    mutdict=defaultdict(list)
    head=0
    #File with all the sample names
    for exp1 in open(mutfile,"rU").xreadlines():
        #print exp1
        lin=exp1.rstrip('\r\n')
        lin=string.split(lin,"\t")
        if len(lin)>3:
            if head==0:
                for i in lin[1:]:
            
                    samplelist.append(i)
                head=1
                continue
            else:
                
                for j in range(1,len(lin)):
                    if lin[j]==str(1):
                        mutdict[lin[0]].append(samplelist[j-1])
                        
        else:
            mutdict[lin[2]].append(lin[0])
   
    return  mutdict

def Zscore(r,n,N,R):
    """where N is the total number of events measured: 
    R is the total number of events meeting the criterion:
    n is the total number of events in this specific reference gene-set: 
    r is the number of events meeting the criterion in the examined reference gene-set: """
    N=float(N) ### This bring all other values into float space
    z = (r - n*(R/N))/math.sqrt(n*(R/N)*(1-(R/N))*(1-((n-1)/(N-1))))
    return z               
                
if __name__ == '__main__':

    import getopt
  
    mutdict=defaultdict(list)
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['Inputfile=','Reference=','Expand='])
        for opt, arg in options:
            if opt == '--Inputfile': Inputfile=arg
            elif opt == '--Reference':Reference=arg
            elif opt =='--Expand': Expand=arg
            
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    mutfile=Reference
         
    header=header_file(mutfile)
    mutdict=findsiggenepermut(mutfile)
    Enrichment(Inputfile,mutdict,mutfile,Expand,header)

        
        
        