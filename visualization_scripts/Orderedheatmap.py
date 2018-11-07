#!/usr/bin/env python
from __future__ import print_function
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import numpy as np
import pylab as pl
import sys,string
import os
import os.path
import scipy
import operator
from collections import OrderedDict
from collections import defaultdict
from operator import itemgetter
from visualization_scripts import clustering; reload(clustering)
import export

def Classify(filename,Mutlabels={},dire="",flag=True):
    count=0
    start=1
    orderdict=OrderedDict()
    countdict=OrderedDict()
  
    countlst=[]
    
    Y=[]
    head=0
    rownames=[]
    colnames=[]
    q=[]
    Z=[]
    if dire!="":
        output_dir = dire+'Results'
        export.createExportFolder(output_dir)
        if flag:
            output_file=output_dir+"/Consolidated-Increasing"+".txt"
        else:
            output_file=output_dir+"/Consolidated-Decreasing"+".txt"
    else:
        output_file=filename[:-4]+"-ordered.txt"
    export_object = open(output_file,'w')
    for line in open(filename,'rU').xreadlines():
        if head >0:
            val=[]
            counter2=0
            val2=[]
            me=0.0
            
            line=line.rstrip('\r\n')
            
            q= string.split(line,'\t')
           # rownames.append(q[0])
            if q[0]=="":
                continue
            orderdict[q[0]]=[q[0],]
            for i in range(start,len(q)):
                try:
                    val2.append(float(q[i]))
                    try:
                        orderdict[q[0]].append(float(q[i]))
                    except Exception:
                        orderdict[q[0]]=[float(q[i]),]
                    try:
                        countdict[i].append(float(q[i]))
                    except Exception:
                        countdict[i]=[float(q[i]),]
                except Exception:
                    continue
            
            
            count+=1
        else:
            #export_object.write(line)
            head=1
            line=line.rstrip('\r\n')
            
            q= string.split(line,'\t')
            header=q
            continue
    
    for i in countdict:
       
        countlst.append(sum(countdict[i]))
    #print countlst
    
    B=sorted(range(len(countlst)),key=lambda x:countlst[x],reverse=flag)
    C=sorted(range(len(countlst)),key=lambda x:B[x])
    
    qu=0
    for i in orderdict.keys():
        Y.append(orderdict[i])
        qu+=1
        #print Y
    
    for i in range(0,len(C)):
        jk= C.index(i)+1
        #print jk
        #print Y[jk]
        Y=sorted(Y,key=itemgetter(jk))
        
        #orderdict=OrderedDict(sorted(orderdict,key=itemgetter(jk)))
        #colnames.append(header[C.index(i)+1])
    
    Y=np.array(Y)
    Y=zip(*Y)
    Y=np.array(Y)
    Z.append(Y[0,:])
    for i in range(0,len(C)):
        jk= C.index(i)+1
        Z.append(Y[jk,:])
    Z=np.array(Z)
    q= Z.shape

    export_object.write("uid")
  
        
    for i in range(q[1]):
        export_object.write("\t"+Z[0][i])
    export_object.write("\n")
    for ij in range(1,q[0]):
        jk= C.index(ij-1)+1
        if header[jk] in Mutlabels:
           export_object.write(Mutlabels[header[jk]]) 
        else:
            export_object.write(header[jk])
        for jq in range(0,q[1]):
            export_object.write("\t"+str(Z[ij][jq]))
        export_object.write("\n")
    export_object.close()
        
    graphic_links=[]
    row_method = None
    column_method=None
    column_metric='cosine'
    row_metric='cosine'
    color_gradient = 'yellow_black_blue'
    transpose=False
    graphic_links = clustering.runHCexplicit(output_file,graphic_links, row_method, row_metric, column_method, column_metric, color_gradient, transpose, display=False, Normalize=False)
        
if __name__ == '__main__':

    import getopt
    group=[]
    grplst=[]
    name=[]
    matrix={}
    compared_groups={}

    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1: sys.exit()
    else:
        analysisType = []
    options, remainder = getopt.getopt(sys.argv[1:],'', ['Guidedir='])
    for opt, arg in options:
        if opt == '--Guidedir':
            Guidedir=arg
         
    Classify(Guidedir)