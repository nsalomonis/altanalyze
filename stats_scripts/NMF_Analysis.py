#!/usr/bin/env python

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
    import numpy as np
    from sklearn.cluster import KMeans
    import nimfa
    from sklearn.decomposition import NMF
import os.path
from collections import defaultdict
import traceback
import export
from visualization_scripts import Orderedheatmap
#import statistics

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filterRows(input_file,output_file,filterDB=None,logData=False):
    orderlst={}
    counter=[]
    export_object = open(output_file,'w')
    firstLine = True
    Flag=0

    for line in open(input_file,'rU').xreadlines(): ### Original expression file (source IDs)
        #for i in filterDB:
            flag1=0
            data = cleanUpLine(line)     
            values = string.split(data,'\t')
            if firstLine:
                firstLine = False
                if Flag==0:
                    export_object.write(line)
            else:
                #print values[0], filterDB
                #sys.exit()
                uid = values[0]
                if uid in filterDB:
                    counter=[index for index, value in enumerate(filterDB) if value == uid]
                    for it in range(0,len(counter)):
                        orderlst[counter[it]]=line
    try:
        for i in range(0,len(orderlst)):
            export_object.write(orderlst[i])
    except Exception:
        print i,filterDB[i]
    export_object.close()
    print 'Filtered rows printed to:',output_file

def FilterGuideGeneFile(Guidefile,Guidefile_block,expressionInputFile,iteration,platform,uniqueIDs,symbolIDs):
    """ Filters the original input expression file for Guide3 genes/events. Needed
    Since NMF only can deal with positive values [Guide3 has negative values]"""
    
    root_dir = export.findParentDir(expressionInputFile)[:-1]
    if 'ExpressionInput' in root_dir:
        root_dir = export.findParentDir(root_dir)
    
    if 'Clustering' in Guidefile:
        count=1
        flag=True
        rank_Count=0
        prev=0
    else:
        count=0
    val=[]
    head=0
    for line in open(Guidefile_block,'rU').xreadlines():
        if head >count:
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            #val.append(q[0])
            if flag:
                if int(q[1])==prev:
                    continue
                else:
                    rank_Count+=1
                    prev=int(q[1])
        else:
            head+=1
            continue
    head=0
    for line in open(Guidefile,'rU').xreadlines():
        line=line.rstrip('\r\n')
        q= string.split(line,'\t')
        n=len(q)
        if head >count:
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            uid = q[0]
            if uid not in uniqueIDs:
                if uid in symbolIDs:
                    uid = symbolIDs[uid]
            val.append(uid)
            if platform != "PSI" and head==2:
                rank_Count=rank_Count+int(q[1])
                print rank_Count
            head=head+1
        else:
            head+=1
            if platform != "PSI" and q[0]=="column_clusters-flat":
                    rank_Count=int(q[n-1])
            continue

    output_dir = root_dir+'/NMF-SVM'
    if os.path.exists(output_dir)==False:
        export.createExportFolder(output_dir)
    
    output_file = output_dir+'/NMFInput-Round'+str(iteration)+'.txt'
    filterRows(expressionInputFile,output_file,filterDB=val)
    
    return output_file,rank_Count


def NMFAnalysis(expressionInputFile,NMFinputDir,Rank,platform,iteration=0,strategy="conservative"):

    root_dir = export.findParentDir(NMFinputDir)[:-1]
    if 'ExpressionInput' in root_dir:
        root_dir = export.findParentDir(root_dir)
    if 'NMF-SVM' in root_dir:
        root_dir = export.findParentDir(root_dir)
        
    export.findFilename(NMFinputDir)
        
    X=[]
    header=[]
    head=0
    exportnam=root_dir+'/NMF-SVM/NMF/round'+str(iteration)+'NMFsnmf_versionr'+str(Rank)+'.txt'
    export_res=export.ExportFile(exportnam)
    exportnam_bin=root_dir+'/NMF-SVM/NMF/round'+str(iteration)+'NMFsnmf_binary'+str(Rank)+'.txt'
    export_res1=export.ExportFile(exportnam_bin)
    exportnam_bint=root_dir+'/NMF-SVM/NMF/round'+str(iteration)+'NMFsnmf_binary_t_'+str(Rank)+'.txt'
    export_res5=export.ExportFile(exportnam_bint)
    MF_input = root_dir+'/NMF-SVM/ExpressionInput/exp.NMF-MarkerFinder.txt'
    export.customFileCopy(expressionInputFile,root_dir+'/NMF-SVM/ExpressionInput/exp.NMF-MarkerFinder.txt')
    export_res4=open(string.replace(MF_input,'exp.','groups.'),"w")
    export_res7=open(string.replace(MF_input,'exp.','comps.'),"w")
    exportnam2=root_dir+'/NMF-SVM/SubtypeAnalyses/round'+str(iteration)+'Metadata'+str(Rank)+'.txt'
    export_res2=export.ExportFile(exportnam2)
    exportnam3=root_dir+'/NMF-SVM/SubtypeAnalyses/round'+str(iteration)+'Annotation'+str(Rank)+'.txt'
    export_res3=export.ExportFile(exportnam3)
    if 'Clustering' in NMFinputDir:
        count=1
        start=2
    else:
        count=0
        start=1
    #print Rank
    for line in open(NMFinputDir,'rU').xreadlines():
        line=line.rstrip('\r\n')
        q= string.split(line,'\t')
        if head >count:
            val=[]
            val2=[]
            me=0.0
            
            for i in range(start,len(q)):
                try:
                    val2.append(float(q[i]))
                except Exception:
                    continue
            me=np.median(val2)
            for i in range(start,len(q)):
                try:
                    val.append(float(q[i]))
                except Exception:
                    val.append(float(me))
            #if q[1]==prev:
            X.append(val)
          
        else:
            export_res1.write(line)
            export_res.write(line)
            export_res1.write("\n")
            #export_res4.write(line)
            #export_res4.write("\n")
            export_res.write("\n")
            header=q
            head+=1
            continue   
    group=defaultdict(list)
        
    sh=[]
    X=np.array(X)
    #print X.shape
    mat=[]
    #mat=X
    mat=zip(*X)
    mat=np.array(mat)
    #print mat.shape
    #model = NMF(n_components=15, init='random', random_state=0)
    #W = model.fit_transform(mat)
    nmf = nimfa.Snmf(mat,seed="nndsvd", rank=int(Rank), max_iter=20,n_run=1,track_factor=False,theta=0.95)
    nmf_fit = nmf()
    W = nmf_fit.basis()
    W=np.array(W)
    #np.savetxt("basismatrix2.txt",W,delimiter="\t")
    H=nmf_fit.coef()
    H=np.array(H)
   # np.savetxt("coefficientmatrix2.txt",H,delimiter="\t")
    #print W.shape
    sh=W.shape
    export_res3.write("uid\tUID\tUID\n")
    if int(Rank)==2:
        par=1
    else:
        par=2
    #for i in range(sh[1]):
    #    val=W[:,i]
    #    me=np.mean(val)
    #    st=np.std(val)
    #    export_res2.write(header[i+1])
    #    for j in range(sh[0]):
    #        if float(W[i][j])>=float(me+(par*st)):
    #          
    #            export_res2.write("\t"+str(1))
    #        else:
    #            export_res2.write("\t"+str(0))
    #       
    #    export_res2.write("\n")
    if platform != 'PSI':
        sh=W.shape
        Z=[]
        export_res5.write("uid")
        export_res2.write("uid")
        for i in range(sh[1]):
            
            export_res5.write("\t"+'V'+str(i))
            export_res2.write("\t"+'V'+str(i))
            export_res3.write('V'+str(i)+"\t"+"Covariate"+"\t"+str(1)+"\n")
            
        export_res5.write("\n")
        export_res2.write("\n")
        export_res3.write("\n")
        for i in range(sh[0]):
            new_val=[]
            val=W[i,:]
            export_res2.write(header[i+1])
            export_res5.write(header[i+1])
            export_res4.write(header[i+1])
            flag=True
            for j in range(sh[1]):
                if W[i][j]==max(val) and flag:
                    export_res5.write("\t"+str(1))
                    export_res2.write("\t"+str(1))
                    new_val.append(1)
                    export_res4.write("\t"+str(j+1)+"\t"+'V'+str(j))
                    flag=False
                else:
                    export_res5.write("\t"+str(0))
                    export_res2.write("\t"+str(0))
                    new_val.append(0)
                
            Z.append(new_val)
            export_res5.write("\n")
            export_res2.write("\n")
            export_res4.write("\n")
        W=zip(*W)
        W=np.array(W)
        sh=W.shape
        Z=zip(*Z)
        Z=np.array(Z)
        for i in range(sh[0]):
            export_res.write('V'+str(i))
            export_res1.write('V'+str(i))
            for j in range(sh[1]):
                export_res.write("\t"+str(W[i][j]))
                export_res1.write("\t"+str(Z[i][j]))
            export_res.write("\n")
            export_res1.write("\n")
            
        export_res.close()
        export_res1.close()
        export_res2.close()
        export_res5.close()
        Orderedheatmap.Classify(exportnam_bint)    
        
        return exportnam,exportnam_bin,exportnam2,exportnam3
    
    else:
        W=zip(*W)
        W=np.array(W)
        sh=W.shape
        Z=[]
        for i in range(sh[0]):
            new_val=[]
            val=W[i,:]
            num=sum(i > 0.10 for i in val)
            if num >40 or num <3:
                compstd=True
            else:
                compstd=False
            me=np.mean(val)
            st=np.std(val)
            #print 'V'+str(i)
            export_res.write('V'+str(i))
            export_res1.write('V'+str(i))
           
            for j in range(sh[1]):
                
                if compstd:   
                    if float(W[i][j])>=float(me+(par*st)):
                    
                        export_res1.write("\t"+str(1))
                        new_val.append(1)
                    else:
                        export_res1.write("\t"+str(0))
                        new_val.append(0)
                else:
                    if float(W[i][j])>0.1:
                    
                        export_res1.write("\t"+str(1))
                        new_val.append(1)
                    else:
                        export_res1.write("\t"+str(0))
                        new_val.append(0)
                export_res.write("\t"+str(W[i][j]))
                
            Z.append(new_val)
            export_res.write("\n")
            export_res1.write("\n")
       # Z=zip(*Z)
        Z=np.array(Z)
        sh=Z.shape
        Z_new=[]
        val1=[]
        Z1=[]
        dellst=[]
        export_res2.write("uid")
        export_res5.write("uid")
        for i in range(sh[0]):
            indices=[]
            val1=Z[i,:]
            sum1=sum(val1)
            flag=False
            indices=[index for index, value in enumerate(val1) if value == 1]
            for j in range(sh[0]):
                val2=[]
                
                if i!=j:
                    val2=Z[j,:]
                    
                    sum2=sum([val2[x] for x in indices])
                    summ2=sum(val2)
                    try:
                        if float(sum2)/float(sum1)>0.5:
                            if summ2>sum1:
                                flag=True
                                #print str(i)
                    except Exception:
                        continue
            if flag==False:
    
                Z1.append(val1)
                export_res2.write("\t"+'V'+str(i))
                export_res5.write("\t"+'V'+str(i))
                export_res3.write('V'+str(i)+"\t"+"Covariate"+"\t"+str(1)+"\n")
        
        export_res2.write("\n")
        export_res5.write("\n")
        Z1=np.array(Z1)
        Z=Z1
        Z=zip(*Z)
        Z=np.array(Z)
        sh=Z.shape
            
        for i in range(sh[0]):
            val1=Z[i,:]
            #print sum(val1)
            #if sum(val)>2: 
            if sum(val1)>2:
                val=[0 if x==1 else x for x in val1]
            else:
                val=val1
            me=np.mean(val)
            st=np.std(val)
            export_res2.write(header[i+1])
            export_res5.write(header[i+1])
            for j in range(sh[1]):
                if strategy=="conservative":
                    export_res2.write("\t"+str(val1[j]))
                    export_res5.write("\t"+str(val1[j]))
                else:
                   export_res2.write("\t"+str(val[j]))
                   export_res5.write("\t"+str(val[j])) 
            export_res2.write("\n")
            export_res5.write("\n")
            Z_new.append(val)
        Z_new=zip(*Z_new)
        Z_new=np.array(Z_new)
        
        sh=Z_new.shape

        export_res5.close()
        Orderedheatmap.Classify(exportnam_bint)    
        if strategy=="conservative":
            return exportnam,exportnam_bin,exportnam2,exportnam3
        else:
            return exportnam,exportnam_bin,exportnam2,exportnam3
    
if __name__ == '__main__':

    import getopt
  
    mutdict=defaultdict(list)
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['Guidefile=','Rank=','PSI='])
        for opt, arg in options:
            if opt == '--Guidefile': Guidefile=arg
            elif opt == '--Rank':Rank=arg
            elif opt == '--PSI':PSI=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
   
    inputfile=Guidefile
 
    Rank=Rank
    if Rank>1:
       NMFAnalysis(inputfile,Rank,platform="RNASeq")
    else:
        pass
   