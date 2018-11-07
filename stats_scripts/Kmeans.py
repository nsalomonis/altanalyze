#!/usr/bin/env python

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import numpy as np
import pylab as pl
import os.path
from collections import defaultdict
from sklearn.cluster import KMeans
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
                header=string.split(line,'\t')
                head=1
            else:break
    return header

def KmeansAnalysis(filename,header,InputFile,turn):
    X=defaultdict(list)
    prev=""
    head=0
    for line in open(filename,'rU').xreadlines():
        if head >1:
            val=[]
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            for i in range(2,len(q)):
                val.append(float(q[i]))  
            if q[1]==prev:
                X[prev].append(val)
            else:
                prev=q[1]
                X[prev].append(val)
        else:
            head+=1
            continue
    
    for key in X:
        print key
        X[key]=np.array(X[key])
        print X[key].shape
        mat=[]
        dire= export.findParentDir(export.findParentDir(InputFile)[:-1])
        output_dir = dire+'SVMOutputs'
        if os.path.exists(output_dir)==False:
            export.createExportFolder(output_dir)
   
        exportname=output_dir+'/round'+str(turn)+'Kmeans_result.txt'
        #exportname=filename[:-4]+key+'.txt'
        export_results=open(exportname,"w")
        mat=zip(*X[key])
        mat=np.array(mat)
        print mat.shape
        kmeans = KMeans(n_clusters=2, random_state=0).fit(mat)
        
        y=kmeans.labels_
        #cent=kmeans.cluster_centers_
        y=y.tolist()
        total=len(y)
        cent_1=y.count(0)
        cent_2=y.count(1)
        print cent_1,cent_2
        export_results.write("uid"+"\t"+"group"+"\n")
        if cent_1<cent_2:
            count=2
            for j in y:
                if j==0:
                    export_results.write(header[count]+"\t"+"1"+"\n")
                else:
                    export_results.write(header[count]+"\t"+"0"+"\n")
                count+=1
        else:
            count=2
            for j in y:
                if j==1:
                    export_results.write(header[count]+"\t"+"1"+"\n")
                else:
                    export_results.write(header[count]+"\t"+"0"+"\n")
                count+=1
                
if __name__ == '__main__':

    import getopt
  

    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['Guidefile='])
        for opt, arg in options:
            if opt == '--Guidefile': Guidefile=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
                
    #Guidefile="/Users/meenakshi/Documents/leucegene/ICGS/Round2_cor_0.6_280default/Clustering-exp.round2_insignificantU2like-Guide1 DDX5&ENSG00000108654&E3.4-E3.9__ENSG0000010-hierarchical_cosine_correlation.txt"          
    header=header_file(Guidefile)
    KmeansAnalysis(Guidefile)

