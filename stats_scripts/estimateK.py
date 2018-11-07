#!/usr/bin/env python

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import os,string
import numpy as np
from sklearn.preprocessing import scale
from numpy import linalg as LA
import scipy
def estimateK(inputfile):
    header=[]
    X=[]
    head=0
    counter=0
    hgv={}
    hgvgenes=[]
    diclst={}
    for line in open(inputfile,'rU').xreadlines():
        if head==0:
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            header=q
            head=1
            continue
        else:
            val=[]
            
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            #header.append(q[0])
            for i in range(1,len(q)):
                try:
                    val.append(float(q[i]))
                   
                except Exception:
                    continue
           
                
            counter+=1
            
             #   break

        X.append(val)
    #X=zip(*X)
    X=np.array(X)
    n=float(X.shape[0])
    p=float(X.shape[1])
    print n
    print p
    
    
    X=scale(X)
    Xt=np.transpose(X)
    
    muTW=float((np.sqrt(n-1))+float(np.sqrt(p)))**2.0

    sigmaTW=(float(np.sqrt(n - 1.0)) + float(np.sqrt(p))) * (1.0/float(np.sqrt(n - 1)) + 1.0/float(np.sqrt(p)))**(1.0/3.0)

    sigmaHat=np.dot(Xt,X)
   
    bd = 3.273 * sigmaTW + muTW
    print bd
    w,v = LA.eig(sigmaHat)
    w=w.tolist()

    k=0
    for i in range(len(w)):
        try:
            if w[i]>bd:
                k=k+1
        except Exception:
            if w[i].real>bd:
                k=k+1
    print k
    return k
    
inputfile="/Volumes/Pass/Immune-complete/ExpressionInput/OncoInputs/NMFInput-Round1.txt"
estimateK(inputfile)
#inputfile="/Volumes/Pass/Singlecellbest/Pollen_upd/ExpressionInput/SamplePrediction/input-CORRELATED-FEATURES.txt"
#estimateK(inputfile)
#inputfile="/Volumes/Pass/Singlecellbest/Usoskin_upd/ExpressionInput/SamplePrediction/input-CORRELATED-FEATURES.txt"
#estimateK(inputfile)
#inputfile="/Volumes/Pass/Singlecellbest/Zeisel_upd/ExpressionInput/SamplePrediction/input-CORRELATED-FEATURES.txt"
#estimateK(inputfile)
    
    
    