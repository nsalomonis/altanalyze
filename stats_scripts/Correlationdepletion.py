#!/usr/bin/env python

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import numpy as np
import pylab as pl
import os.path
import scipy
from collections import defaultdict
from sklearn.cluster import KMeans
from import_scripts import sampleIndexSelection
import export


def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue
            
def header_file(fname, delimiter=None):
    head=0
    header=[]
    with open(fname, 'rU') as fin:
        for line in fin:
            if head==0:
                line = line.rstrip(os.linesep)
                line=string.split(line,'\t')
                for i in line:
                    if ":" in i:
                        i=string.split(i,":")
                        header.append(i[1])
                    else:
                        header.append(i)
                
                del header[:1]
                head=1
            else:break
    return header

def DepleteSplicingevents(commonkeys,keylabel,count,InputFile):
    eventlist=[]
    #exportname=keylabel[:-4]+'correlationSelected_0.3.txt'
    #export_res=open(exportname,"w")
    exportdep=InputFile[:-4]+'cor_depleted.txt'
    export_res1=open(exportdep,"w")
    #export_res.write("splicingevent"+"\t"+"key"+"\n")
    for event in commonkeys:
        if commonkeys[event]==count:
            eventlist.append(event)
     #       export_res.write(event+"\t"+str(1)+"\n")
    head=0
    for line in open(keylabel,'rU').xreadlines():
        if head==0:
            export_res1.write(line)
            head=1
            continue
        else:
            line1=line.rstrip('\r\n')
            q= string.split(line1,'\t')
            if q[0] in eventlist:
                export_res1.write(line)
    return exportdep
            
def FindCorrelations(filename,PSIfile,name):
    X=defaultdict(list)
    prev=""
    head=0
    for line in open(filename,'rU').xreadlines():
        if head >0:
            val=[]
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            for i in range(1,len(q)):
                val.append(float(q[i]))
            flag=0
            for i in range(len(name)):
                key1=q[0]+"_vs"
                key2="vs_"+q[0]+".txt"
            
                if key1 in name[i] or key2 in name[i]:
                    flag=1
            
            if flag==1:
                if q[0]==prev:
                    X[prev].append(val)
                else:
                    prev=q[0]
                    X[prev].append(val)
        else:
            #print line
            head+=1
            continue
    
    head=0
    matrix=[]
    eventnames=[]
    Y={}
    eventkeys=defaultdict(list)
    commonkeys=defaultdict(int)
    count=len(X)
    for key in X:
        #print key
        X[key]=np.array(X[key])
        #print X[key].shape
        mat=[]
        mat=zip(*X[key])
        mat=np.array(mat)
        #print mat.shape
        mat=np.mean(mat,axis=1)
        
        Y[key]=np.array(mat)
   
       
    counter=defaultdict(int)
    for line in open(PSIfile,'rU').xreadlines():
        if head >0:
            for key in Y:
                
                list1=[]
                list2=[]
                mean=0.0
                
                line=line.rstrip('\r\n')
                #print line
                q= string.split(line,'\t')
                eventnames.append(q[0])
                #print len(Y[key])
                for i in range(1,len(q)):
                    try:
                        
                        list1.append(float(q[i]))
                        list2.append(float(Y[key][i-1]))
                    except Exception:
                            continue
                #print len(list1),len(list2)
                rho=scipy.stats.pearsonr(list1,list2)
                if abs(rho[0])<0.3:
                    
                    commonkeys[q[0]]+=1
                    counter[key]+=1
                else:
                    eventkeys[key].append([q[0],rho[0]])
        
        else:
           # print line
            head+=1
                
            continue
    #for key in counter:
        #print counter[key]
     #   export_key=open(PSIfile[:-4]+str(key)+'.txt',"w")
      #  for i,j in eventkeys[key]:
       #     export_key.write(i+"\t"+str(j)+"\n")
        
    return commonkeys,count

if __name__ == '__main__':

    import getopt
  

    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['Guidefile=','PSIfile='])
        for opt, arg in options:
            if opt == '--Guidefile': Guidefile=arg
            elif opt =='--PSIfile':PSIfile=arg
           
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
       
#filename="/Users/meenakshi/Documents/leucegene/ICGS/Clustering-exp.Hs_RNASeq_top_alt_junctions367-Leucegene-75p_no149-Guide1 TRAK1&ENSG00000182606&I1.1_42075542-E2.1__E-hierarchical_cosine_correlation.txt"          
#PSIfile="/Users/meenakshi/Documents/leucegene/ExpressionInput/exp.Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-367-Leucegene-75p-unique-filtered-filtered.txt"
#keylabel="/Users/meenakshi/Documents/leucegene/ExpressionInput/exp.round2_glmfilteredKmeans_label.txt"
    header=header_file(Guidefile)
    output_file=PSIfile[:-4]+"-filtered.txt"
    sampleIndexSelection.filterFile(PSIfile,output_file,header)
    commonkeys,count=FindCorrelations(Guidefile,output_file)
    DepleteSplicingevents(commonkeys,output_file,count)