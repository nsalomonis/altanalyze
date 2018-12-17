#!/usr/bin/env python

#Copyright 2017 Cincinnati Children's Hospital Medical Center, Research Foundation
#Author Meenakshi Venkatasubramanian - altanalyze@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

""" ICGS-NMF Module (Combatible with ICGS2 and splice-ICGS)
https://github.com/venkatmi/oncosplice
Steps applied in this workflow:
1 - Run splice-ICGS (Feature Selection)
2 - Block identification (Rank analysis)
3 - NMF Analysis (Initial subtype identification)
4 - Filter Event Annotation
5 - Meta data analysis (differential expression)
6 - Expand clusters (SVM sample classification)
7 - Mutation enrichment (MAF or VCF - optional)
8 - Correlation depletion (excluded biological confounding signatures)
"""

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import traceback
import sys, string, os
import RNASeq
import numpy as np
from stats_scripts import RNASeq_blockIdentification
from stats_scripts import NMF_Analysis; reload(NMF_Analysis)
from stats_scripts import filterEventAnnotation
from stats_scripts import metaDataAnalysis
from stats_scripts import ExpandSampleClusters; reload(ExpandSampleClusters)
from import_scripts import sampleIndexSelection
from stats_scripts import Correlationdepletion
import UI
import multiprocessing as mlp
import export
upd_guides=[]
import operator
from collections import OrderedDict
from collections import defaultdict
from stats_scripts import Kmeans
from stats_scripts import MutationEnrichment_adj as ME
from visualization_scripts import Orderedheatmap
from visualization_scripts import clustering; reload(clustering)
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.neighbors import KDTree
import community
import collections
from scipy.stats import variation
import networkx as nx
from sklearn.preprocessing import scale
from numpy import linalg as LA
import scipy
import warnings

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
    try:
        n=float(X.shape[0])
        p=float(X.shape[1])
    except: ### dimension error - assume k=30
        return 15
        
    X=scale(X)
    Xt=np.transpose(X)
    
    muTW=float((np.sqrt(n-1))+float(np.sqrt(p)))**2.0

    sigmaTW=(float(np.sqrt(n - 1.0)) + float(np.sqrt(p))) * (1.0/float(np.sqrt(n - 1)) + 1.0/float(np.sqrt(p)))**(1.0/3.0)

    sigmaHat=np.dot(Xt,X)
   
    bd = 3.273 * sigmaTW + muTW
    
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
    return k

def caldist(X,i,keys,keylist):
        D=[]
        Xxd=[]
        newlist=[]
    #for i in range(len(visited)):
        #Xd=np.array(X[i])
        #Xd=Xd.reshape(1, -1)
        for ii in keys:
            if ii==i: continue
            newlist.append(ii)
            Xxd.append(X[ii].tolist())
        
        Xxd=np.array(Xxd)
        Xd=X[i]
        
        #Xd=Xxd
        #Xxd=Xxd.tolist()
        Xd=Xd.reshape(1, -1)
        D=pairwise_distances(Xd,Xxd,metric='euclidean').tolist()
        
        for q in range(len(np.argsort(D)[0])):
            if newlist[q] in keylist:
                continue
            else:
                key1=newlist[q]
                break
        return key1

def hgvfinder(inputfile):
    header=[]
    X=[]
    head=0
    counter=0
    hgv={}
    hgvgenes=[]
   
    for line in open(inputfile,'rU').xreadlines():
        if head==0:
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            count=len(q)-1
            if count >20000:
                umap=True
            else:
                umap=False
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
            coun=len(set(val))
            qq=q[0].lower()
          
            if (qq.startswith("rpl") or qq.startswith("rps") or qq.startswith("mt-") or qq.startswith("ig")) and umap:
                continue
            else:
             if coun >5:
                disp=float(np.var(val))/float(np.mean(val))
                #me=float(np.mean(val))
                hgv[q[0]]=disp
                
            counter+=1
            #if counter%500==0: print counter,
             #   break
    #with open('hgv_0.1.txt', 'w') as f:
    #    for item in hgv:
    #        f.write(str(item)+"\t"+str(hgv[item]))
    #        f.write("\n")
    #

    hgv= sorted(hgv.items(), key=operator.itemgetter(1),reverse=True)
    counter=0

    for item,item2 in hgv:
        if counter<500: 
           
            hgvgenes.append(item)
            counter+=1
        
    output_file=inputfile[:-4]+'-filtered.txt'
    #copy sample index selection file-mv
    sampleIndexSelection.filterRows(inputfile,output_file,hgvgenes)
    return output_file,count

def UMAPsampling(inputfile):
    """ This function performs downsampling of the input data using UMAP for non-linear
    Uniform Manifold Approximation and Projection dimensionality reduction used to identify
    initial distribution of cells, then Louvain clustering using the minimum resolution to
    identify discrete initial clusters. """
    
    header=[]
    X=[]
    head=0
    counter=0
    hgv={}
    hgvgenes=[]
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
    
    X=zip(*X)
    X=np.array(X)
  
    n=X.shape[0]
    sampmark=[]    
    nn=X.shape[0]
    nm=X.shape[1]
    try: import umap
    except:
        from visualization_scripts.umap_learn import umap
    print "running UMAP"  
    model=umap.UMAP()

    from annoy import AnnoyIndex
    t=AnnoyIndex(nm)
    for i in range(nn):
      try:  t.add_item(i,X[i])
      except Exception: print i
    t.build(100)
   # t.save('test.ann')
    #u=AnnoyIndex(nm)
    diclst={}
    #u.load('test.ann')
    #n1=25
    print "creating graphs"
    for i in range(nn):
    #ind = tree.query([Xtemp[i]],k=10,return_distance=False,dualtree=True)
        ind=t.get_nns_by_item(i,10)
        diclst[i]=ind

    print "creating graphs"
    G=nx.from_dict_of_lists(diclst)
   # nx.write_adjlist(G,"test.adjlist")
    #G=nx.read_adjlist("test.adjlist")
    dendrogram= community.generate_dendrogram(G)
    #for level in range(len(dendrogram) - 1):
    level=0
    pr= community.partition_at_level(dendrogram,level)
    commun={}
    comval={}
    for key1 in pr:
        try: commun[pr[key1]].append(key1)
        except Exception: commun[pr[key1]]=[key1,]
        try: comval[pr[key1]].append(X[int(key1)])
        except Exception: comval[pr[key1]]=[X[int(key1)],]
    
    print "Finding medians"
    comindices=[]
   
    for key1 in comval:
        k=10000/len(comval)
        if k<1: k=1
        k2=len(comval[key1])
        matri=np.array(comval[key1])
        matri=np.array(matri)
   
    #n=matri.shape[0]
        D=pairwise_distances(matri,metric='euclidean').tolist()
        D=np.array(D)
    
        dist=np.mean(D,0)
        
        if k2<k:
            k=k2
       
        count=0
        for i in np.argsort(dist):
            if count<k:
                comindices.append(commun[key1][i])
                count=count+1
    sampmark=[]

    for key1 in comindices:
        #if count<2500:
            #print key1
        key=int(key1)
        sampmark.append(header[key+1])
    return sampmark
    
def PageRankSampling(inputfile):
    """ Google PageRank algorithm from networkX for graph-based link analysis """
    header=[]
    X=[]
    head=0
    counter=0
    hgv={}
    hgvgenes=[]
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
        X.append(val)
    
    X=zip(*X)
    X=np.array(X)
    n=X.shape[0]
    sampmark1=[]
    
    for iq in range(0,n,20000):
        jj=2500
        if iq+24999>n:
            j=n-iq
        else:
            j=24999
        jj=int(float(j+1)/4.0)
        jj=2500
        #if jj<2500 and n<3000:
            #jj=n

        Xtemp=X[iq:iq+j,]
        nn=Xtemp.shape[0]
        nm=Xtemp.shape[1]
 
        diclst={}
        from annoy import AnnoyIndex
        t=AnnoyIndex(nm)
        for i in range(nn):
            t.add_item(i,Xtemp[i])
            
        t.build(100)
        t.save('test.ann')
        u=AnnoyIndex(nm)
        u.load('test.ann')
        tree = KDTree(X, leaf_size=10, metric='euclidean')
        #n1=25
        for i in range(nn):
            #ind = tree.query([Xtemp[i]],k=10,return_distance=False,dualtree=True)
            ind=u.get_nns_by_item(i,10)
            diclst[i]=ind
           # diclst[i]=ind.tolist()[0]
        
        print "creating graphs"
        G=nx.from_dict_of_lists(diclst)
        #nx.write_adjlist(G,"test.adjlist")
        #G=nx.read_adjlist("test.adjlist")
        print "computing page rank"
        pr=nx.pagerank(G)
        pr= sorted(pr.items(), key=operator.itemgetter(1),reverse=True)
      
        count=0
        pr1=OrderedDict()
        for (key1,key2) in pr:
            if count<jj:
                #print key1
                key1=iq+int(key1)
                pr1[key1,key2]=[]
                
                #print header[key1-1]
                sampmark1.append(key1)
            count+=1
        #with open('pangranresults_0.1.txt', 'w') as f:
        #for (key1,key2) in pr:
        #f.write(str(key1)+"\t"+str(key2)+"\n")
        #f.write("\n")
        samp=[]
        sampdict={}
        sampmark=[]
        for (key1,key2) in pr1:
            if len(samp)<len(pr1):
                if key1 not in samp:
                    sampdict[key1]=[]
                    neighbours=list(G.adj[key1])
                    
                    samp.append(key1)
                    for ii in range(len(neighbours)):
                        if neighbours[ii] not in samp and neighbours[ii] in sampmark1:    
                            sampdict[key1].append(neighbours[ii])
                            samp.append(neighbours[ii])
                else:
                    dup=[]
                    for key in sampdict:
                        if key1 in sampdict[key]:
                            neighbours=list(G.adj[key1])
                            
                            for ii in range(len(neighbours)):
                                if neighbours[ii] not in samp and neighbours[ii] in sampmark1:    
                                    sampdict[key].append(neighbours[ii])
                                    samp.append(neighbours[ii])         
        key=pr[0][0]
        keylist=[]
        keylist.append(key)
        while len(keylist) <len(sampdict):
            key=caldist(X,key,sampdict,keylist)
            keylist.append(key)
        for keys in range(len(keylist)):
            sampmark.append(header[keylist[keys]+1])
            for i in range(len(sampdict[keylist[keys]])):
                sampmark.append(header[sampdict[keylist[keys]][i]+1])

        #with open('pangranresults_0.1.txt', 'w') as f:
        #for item in range(len(sampmark)):
        #f.write(str(sampmark[item])+"\n")
        #f.write("\n")

    samptemp=[]
    for i in range(len(header)):
       if header[i] in sampmark:
            samptemp.append(header[i])
    
    sampmark=samptemp
    if len(sampmark)>3000 and n<5000:
        output_file=inputfile[:-4]+'-filtered.txt'
        sampleIndexSelection.filterFile(inputfile,output_file,sampmark)
        sampmark=sampling(output_file)
        return sampmark
    else:
        
        return sampmark
  

def filterPSIValues(filename):
    fn = filepath(filename)
    firstRow=True
    header = True
    rows=0
    filtered=0
    new_file = filename[:-4]+'-75p.txt'
    ea = export.ExportFile(new_file)

    for line in open(fn,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            header = False
            eventindex=t.index('EventAnnotation')
            t = [t[1]]+t[eventindex+1:]
            header_length = len(t)-1
            minimum_values_present = int(header_length)-1
            not_detected = header_length-minimum_values_present
            new_line = string.join(t,'\t')+'\n'
            ea.write(new_line)
        else:
            t = [t[1]]+t[eventindex+1:]
            missing_values_at_the_end = (header_length+1)-len(t)
            missing = missing_values_at_the_end+t.count('')
            if missing<not_detected:
                new_line = string.join(t,'\t')+'\n'
                ea.write(new_line)
                filtered+=1
        rows+=1
    ea.close()
    return newfile
    
def header_list(EventAnnot):
    head=0
    header=[]
    with open(EventAnnot, 'rU') as fin:
        for line in fin:
            if head==0:
                line = line.rstrip(os.linesep)
                line=string.split(line,'\t')
                startpos=line.index('EventAnnotation')
                header.append('UID')
                for i in range(startpos+1,len(line)):
                        header.append(line[i])
                head=1
            else:break
    return header

def grpDict(grplst):
    head=0
    header={}
    with open(grplst, 'rU') as fin:
        for line in fin:
                line = line.rstrip(os.linesep)
                line=string.split(line,'\t')
                #for i in range(len(line)):
                try:header[line[2]].append(line[0])
                except Exception: header[line[2]]=[line[0],]
    return header

def FindTopUniqueEvents(Guidefile,psi,Guidedir):
    head=0
    guidekeys=[]
    exportnam=os.path.join(Guidedir,"SplicingeventCount1.txt")
    export_class=open(exportnam,"a")

    tempkeys={}
    global upd_guides
    global train
    omitcluster=0
    
    unique_clusters={}
    for line in open(Guidefile,'rU').xreadlines():
        if head==0:
            line1=line.rstrip('\r\n')
            q= string.split(line1,'\t')
            head=1
            try:
                uid=q.index('UID')
                adjp=q.index('rawp')
                dpsi=q.index('dPSI')
                Clusterid=q.index('UpdatedClusterID')
                cutoff=0.1
                continue
            except Exception:
                uid=q.index('Symbol')
                adjp=q.index('rawp')
                dpsi=q.index('LogFold')
                Clusterid=q.index('Symbol')
                cutoff=0.58
        else:
            line1=line.rstrip('\r\n')
            q= string.split(line1,'\t')
            if abs(float(q[dpsi]))>cutoff and float(q[adjp])<0.01:
                try:
                    tempkeys[q[Clusterid]].append([q[uid],float(q[adjp]),q[adjp+1]])
                except KeyError:
                    tempkeys[q[Clusterid]]=[[q[uid],float(q[adjp]),q[adjp+1]],]
    for i in tempkeys:
        if len(tempkeys[i])>1:
            tempkeys[i].sort(key=operator.itemgetter(1),reverse=False)
            try:
                unique_clusters[0].append(tempkeys[i][0])
            except KeyError:
                unique_clusters[0]=[tempkeys[i][0],]
        else:
            try:
                unique_clusters[0].append(tempkeys[i][0])
            except KeyError:
                unique_clusters[0]=[tempkeys[i][0],]
    try:
        if len(unique_clusters[0])>1:     
            unique_clusters[0].sort(key=operator.itemgetter(1))
            if len(unique_clusters[0])>10:
                guidekeys=unique_clusters[0][0:150]
                for i in range(0,len(guidekeys)):
                    upd_guides.append(guidekeys[i][0])
            else:
                        omitcluster=1
        else:
            omitcluster=1
        export_class.write(psi+"\t"+str(len(unique_clusters[0]))+"\n")
    except Exception:
        omitcluster=1
    return omitcluster

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def MergeResults(dire):
    file_index={}
    count=0
    for filename in os.listdir(dire):

        if ("Results_max" in filename or "Kmeans" in filename)  and "._" not in filename and "ordered" not in filename and "max_t" not in filename:
            file_index[filename]=count
            count+=1

    keylist={}
    heads={}
    for filename in os.listdir(dire):
        if "Results_max" in filename or "Kmeans" in filename:
          if  "._" not in filename and "ordered" not in filename and "max_t" not in filename:
            Guidefile=os.path.join(dire, filename)
            head=0
            for line in open(Guidefile,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                header=[]
                if head==0:
                    head=1
                    for i in range(1,len(t)):
                        header.append(t[i])
                    heads[filename]=header
                        
                    continue
                else:
                    
                    val=[]
                    key=t[0]
                    for i in range(1,len(t)):
                        val.append(t[i])
                    if key not in keylist:
                        keylist[key]=[[file_index[filename],val],]
                    else:
                        keylist[key].append([file_index[filename],val])
    exportnam=os.path.join(dire,"MergedResult.txt")
    export_class=open(exportnam,"w")
    export_class.write("uid")
    
    for filename in file_index:
        export_class.write("\t")
        export_class.write(string.join(heads[filename],"\t"))
    export_class.write("\n")

    for key in keylist:
        export_class.write(key)
        
        for filename in file_index:
            
            for val1,val2 in keylist[key]:
                if file_index[filename]==val1:
                    export_class.write("\t")
                    export_class.write(string.join(val2,"\t"))
               
                    break
        export_class.write("\n")  
    return exportnam

def DetermineClusterFitness(allgenesfile,markerfile,filterfile,BinarizedOutput,rho_cutoff):
    """ Determines whether a cluster has mutiple unique genes and hence should be used for SVM (AKA cluster fitness) """
    header=True
    genes=[]
    nametemp=[]
    for line in open(BinarizedOutput,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            header = False
        else:
            val=[]
            for i in range(1,len(t)):
                val.append(float(t[i]))
            if sum(val)>2:
                nametemp.append(t[0])
            
    header=False
    genes=[]
    for line in open(filterfile,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            header = False
        else:
            genes.append(t[0])
    allgenes={}
    header=True
    name=[]
    for line in open(allgenesfile,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        uid = t[0] ### source ID not converted symbol
        rho = t[2]
        cluster = t[4]
        if header:
            header = False
        else:
            if float(rho)>0.3:
                allgenes[uid]=cluster

    header=True
    markerdict={}
    counter=1
    group=[]
    name=[]
    common_geneIDs=0
    marker_count=0
    for line in open(markerfile,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        uid = t[0]
        rho = t[2]
        cluster = t[4]
        marker_count+=1
        if header:
            header = False
        else:
            if uid in genes:
                common_geneIDs+=1
                #if rho_cutoff>0.4:rho_cutoff=0.4
                rho_cutoff=0.3
                #print rho_cutoff
                #rho_cutoff=0.2
                if float(rho)>rho_cutoff and cluster == allgenes[uid]:
                    try: markerdict[cluster].append([uid,float(rho)])
                    except Exception: markerdict[cluster]=[[uid,float(rho)]]
    if common_geneIDs<marker_count:
        print 'WARNING... only',common_geneIDs, 'out of', marker_count, 'gene IDs matched after conversion.'
        
    for key in markerdict:
        countr=1
        if len(markerdict[key])>=1 and key in nametemp:
            name.append(key+"_vs_Others.txt")
            group.append(counter)
      
            for i,j in markerdict[key] :
                #if countr<30:
                    upd_guides.append(i)
                    countr+=1
            counter+=1

    return upd_guides,name,group

def sortFile(allgenesfile,rho_cutoff):
    markergenes={}
    val=[]
    header=True
    for line in open(allgenesfile,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            header = False
        else:
            values=[]
            for i in range(len(t)):
                if i ==0:
                    values.append(t[i])
                if i ==2:
                    values.append(float(t[2]))
                if i==4:
                    values.append(t[4])
            
            val.append(values)
    val = sorted(val, key = operator.itemgetter(1),reverse=True)
    val = sorted(val, key = operator.itemgetter(2))
    count=0
    prev="NA"
    markerlst=[]
    markergrps={}
    for i in range(len(val)):
     
        if val[i][2]==prev:
            if count<60 and val[i][1]>=0.1: #rho_cutoff
                try:markergrps[val[i][2]].append(val[i][0])
                except Exception:markergrps[val[i][2]]=[val[i][0],]
                markerlst.append(val[i][0])
                count=count+1
                prev=val[i][2]
            else:
                prev=val[i][2]
                continue
        else:
            count=0
            if val[i][1]>=0.1:
                try:markergrps[val[i][2]].append(val[i][0])
                except Exception:markergrps[val[i][2]]=[val[i][0],]
                markerlst.append(val[i][0])
                count=count+1
                prev=val[i][2]
                
    return markergrps,markerlst

def generateMarkerheatmap(processedInputExpFile,output_file,NMFSVM_centroid_cluster_dir,groupsdict,markergrps,header1,outputDir,root_dir,species,uniqueIDs):
    """ Produces a final MarkerFinder result from ICGS-NMF """
    
    matrix={}
    header=True
    samples=[]
    samples2=[]
    samples_all=[]
    samples2_all=[]
    genes=[]
    genes2=[]
    exportnam2=root_dir+'/ICGS-NMF/FinalGroups.txt'
    export_class2=open(exportnam2,"w")
    
    for line in open(NMFSVM_centroid_cluster_dir,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        sampleOrder=[]
        if header:
            for i in range(len(t)):
                if ":" in t[i]:
                    val=string.split(t[i],":")[1]
                    gr=string.split(val,"-")[1]
                    gr=string.split(gr,"_")[0]
                    gr=gr.replace("V","")
                    #sampleOrder.append(string.split(val,"-")[1])
                    sampleOrder.append(gr)
            break
        
    header=True
    samp=[]
    for line in open(processedInputExpFile,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            for i in range(1,len(t)):
                samp.append(t[i])
            header=False
            continue
        else:
            for i in range(1,len(t)):
                matrix[t[0],samp[i-1]]=t[i]
                
    for i in range(len(sampleOrder)):
       
        for j in range(len(groupsdict[sampleOrder[i]])):
            export_class2.write(groupsdict[sampleOrder[i]][j]+"\t"+str(i+1)+"\t"+sampleOrder[i]+"\n")
            if groupsdict[sampleOrder[i]][j] in header1:
                samples.append(sampleOrder[i]+":"+groupsdict[sampleOrder[i]][j])
                samples2.append(groupsdict[sampleOrder[i]][j])
    for i in range(len(sampleOrder)):
        for j in range(len(markergrps[sampleOrder[i]])):
            uid = markergrps[sampleOrder[i]][j]
            genes.append(uid)
            if uid in uniqueIDs:
                symbol = uniqueIDs[uid]
            else:
                symbol = uid
            genes2.append(sampleOrder[i]+":"+uid+' '+symbol)
            
    Outfile = outputDir+'/'+'MarkerFinder-subsampled-ordered.txt'
    exportnam=open(Outfile,"w")     
    exportnam.write("uid"+"\t"+"row_clusters-flat")
   
    for i in range(len(samples)):
        exportnam.write("\t"+samples[i])
    exportnam.write("\n")
    exportnam.write("column_clusters-flat"+"\t")
    for i in range(len(samples)):
        exportnam.write("\t"+"NA")
    exportnam.write("\n")
    
    for i in range(len(genes)):
        exportnam.write(genes2[i]+"\t"+"NA")
        for j in range(len(samples)):
            exportnam.write("\t"+matrix[genes[i],samples2[j]])
        exportnam.write("\n")
    exportnam.close()

    export_class2.close()
    graphic_links=[]
    row_method=None
    column_method=None
    column_metric='euclidean'
    row_metric='correlation'
    color_gradient = 'yellow_black_blue'
    transpose=False
    import UI
    Species=species
    platform="RNASeq"
    Vendor=""
    gsp = UI.GeneSelectionParameters(Species,platform,Vendor)
    gsp.setPathwaySelect('None Selected')
    gsp.setGeneSelection('')
    gsp.setOntologyID('')
    gsp.setGeneSet('None Selected')
    gsp.setJustShowTheseIDs('')                
    gsp.setTranspose(False)
    gsp.setNormalize('median')
    gsp.setGeneSelection('')
    #gsp.setClusterGOElite('GeneOntology')
    gsp.setClusterGOElite('BioMarkers')
    graphic_links = clustering.runHCexplicit(Outfile,graphic_links, row_method, row_metric, column_method,column_metric,color_gradient, gsp, display=False, Normalize=True,contrast=5)
    
    if len(samp)>len(header1):
        Outfile1 = outputDir+'/'+'MarkerFinder-Allsamples-ordered.txt'
        exportnam1=open(Outfile1,"w")
        for i in range(len(sampleOrder)):
            for j in range(len(groupsdict[sampleOrder[i]])):
                samples_all.append(sampleOrder[i]+":"+groupsdict[sampleOrder[i]][j])
                samples2_all.append(groupsdict[sampleOrder[i]][j])
        exportnam1.write("uid"+"\t"+"row_clusters-flat")
        for i in range(len(samples_all)):
            exportnam1.write("\t"+samples_all[i])
        exportnam1.write("\n")
        exportnam1.write("column_clusters-flat"+"\t")
        for i in range(len(samples_all)):
            exportnam1.write("\t"+"NA")
        exportnam1.write("\n")

    
        for i in range(len(genes)):
            exportnam1.write(genes2[i]+"\t"+"NA")
            for j in range(len(samples_all)):
                
                exportnam1.write("\t"+matrix[genes[i],samples2_all[j]])
            exportnam1.write("\n")
        exportnam1.close()
        graphic_links = clustering.runHCexplicit(Outfile1,graphic_links, row_method, row_metric, column_method,column_metric,color_gradient, gsp, display=False, Normalize=True,contrast=5)
        status = 'subsampled'
    else:
        status = 'not-subsampled'

    return status, graphic_links

def callICGS(processedInputExpFile,species,rho_cutoff,dynamicCorrelation,platform,gsp):
    
    #Run ICGS recursively to dynamically identify the best rho cutoff
    graphic_links3,n = RNASeq.singleCellRNASeqWorkflow(species,platform,processedInputExpFile,mlp,exp_threshold=0, rpkm_threshold=0, parameters=gsp)
    if n>5000 and dynamicCorrelation:
            rho_cutoff=rho_cutoff+0.1
            gsp.setRhoCutoff(rho_cutoff)
            print 'Increasing the Pearson rho threshold to:',rho_cutoff
            graphic_links3,n,rho_cutoff=callICGS(processedInputExpFile,species,rho_cutoff,dynamicCorrelation,platform,gsp)
    return graphic_links3,n,rho_cutoff

def getAllSourceIDs(fileame,species):
    unique_ids={}
    symbol_ids={}
    IDtype='Symbol'
    count=0
    typeCount = 0
    for line in open(fileame,'rU').xreadlines():
        data = cleanUpLine(line)
        uid = string.split(data,'\t')[0]
        unique_ids[uid]=''
        if count<100:
            if 'ENS' in uid:
                typeCount+=1
                IDtype='Ensembl'
            else:
                try:
                    int(uid)
                    typeCount+=1
                    IDtype='EntrezGene'
                except Exception:
                    pass
        count+=1

    ### Check to see if these IDs are Ensembl IDs or EntrezGene
    if typeCount>50: ### If over half of the IDs are EntrezGene or Ensembl
        count=0
        try:
            import gene_associations
            gene_annotations = gene_associations.importGeneData(species,IDtype)
        except:
            gene_annotations={}
            
        for uid in gene_annotations:
            if uid in unique_ids:
                unique_ids[uid]=gene_annotations[uid].Symbol() #### Convert to Symbol
                if 'LRG_' not in uid:
                    symbol_ids[gene_annotations[uid].Symbol()]=uid
                    count+=1

        print count, IDtype, 'IDs with corresponding gene symbols out of', len(unique_ids)
    return unique_ids, symbol_ids

def CompleteICGSWorkflow(root_dir,processedInputExpFile,EventAnnot,iteration,rho_cutoff,dynamicCorrelation,platform,species,scaling,gsp):
    """ Run the entire ICGS-NMF workflow, recursively """
    
    originalExpFile = EventAnnot
    
    ### Store a list of all valid original IDs (for ID conversions)
    uniqueIDs, symbolIDs = getAllSourceIDs(processedInputExpFile,species)
    
    if platform=='PSI':
        ### For splice-ICGS, the method performs signature depletion (removes correlated events from the prior round) on the Event Annotation file
        FilteredEventAnnot=filterEventAnnotation.FilterFile(processedInputExpFile,EventAnnot,iteration)
        graphic_links3 = RNASeq.singleCellRNASeqWorkflow(species, 'exons', processedInputExpFile,mlp,exp_threshold=0, rpkm_threshold=0, parameters=gsp)
    else:
        ### For single-cell RNA-Seq - run ICGS recursively to dynamically identify the best rho cutoff
        graphic_links3,n,rho_cutoff=callICGS(processedInputExpFile,species,rho_cutoff,dynamicCorrelation,platform,gsp)
    Guidefile=graphic_links3[-1][-1]
    Guidefile=Guidefile[:-4]+'.txt'

    #Guidefile="/Volumes/Pass/ICGS2_testrun/ExpressionInput/amplify/DataPlots/Clustering-exp.input-Guide3 AREG GZMA BTG1 CCL5 TMSB4X ITGA2B UBE2C IRF-hierarchical_euclidean_correlation.txt"
    #rho_cutoff=0.2
    try:
        print "Running block identification for rank analyses - Round"+str(iteration)
        try:
            RNASeq_blockIdentification.correlateClusteredGenesParameters(Guidefile,rho_cutoff=0.4,hits_cutoff=4,hits_to_report=50,ReDefinedClusterBlocks=True,filter=True) 
            Guidefile_block=Guidefile[:-4]+'-BlockIDs.txt'
        except Exception:
            Guidefile_block=Guidefile
            
        ### Filters the original expression file for the guide3 genes [returns a filename similar to NMFInput-Round1.txt]
        NMFinput,Rank=NMF_Analysis.FilterGuideGeneFile(Guidefile,Guidefile_block,processedInputExpFile,iteration,platform,uniqueIDs,symbolIDs)
        #NMFinput="/Volumes/Pass/ICGS2_testrun/ExpressionInput/ICGS-interim/NMFInput-Round1.txt"
        try: k = int(gsp.K())
        except: k = None; #print traceback.format_exc()
        
        if k==None:
            Rank=estimateK(NMFinput)
            Rank=Rank*2
        else:
            Rank=k
        print "The number target number of clusters (k/rank) is:",k
        
        if Rank>1:
            if Rank>2 and platform=='PSI':
                Rank=30
            if Rank<5 and platform!='PSI':
                Rank=10
            print "Running NMF analyses for dimension reduction using "+str(Rank)+" ranks - Round"+str(iteration)
            
            ### This function prepares files for differential expression analsyis (MetaDataAnalysis), MarkerFinder
            filteredInputExpFile = string.replace(processedInputExpFile,'exp.','filteredExp.')
            if '-OutliersRemoved' in Guidefile:
                filteredInputExpFile = string.replace(filteredInputExpFile,'.txt','-OutliersRemoved.txt')
            NMFResult,BinarizedOutput,Metadata,Annotation=NMF_Analysis.NMFAnalysis(filteredInputExpFile,NMFinput,Rank,platform,iteration)
            
            if platform == 'PSI':
                print "Identifying cluster-specific differential splicing events"
                findmarkers=False
            else:
                print 'Identifying cell-population specific genes'
                findmarkers=True
                
            if findmarkers:
                import markerFinder
                ### Default path for the NMF clustered groups for MarkerFinder analysis
                input_exp_file=root_dir+'/NMF-SVM/ExpressionInput/exp.NMF-MarkerFinder.txt'
                logTransform = False
                    ### Work around when performing this analysis on an alternative exon input cluster file
                group_exp_file = input_exp_file
                fl = UI.ExpressionFileLocationData(input_exp_file,'','',''); fl.setOutputDir(root_dir)
                fl.setSpecies(species); fl.setVendor("3'array")
                rpkm_threshold = 0.00
                fl.setRPKMThreshold(rpkm_threshold)
                fl.setCorrelationDirection('up')
                compendiumType = 'protein_coding'
                genesToReport = 60
                correlateAll = True
                markerFinder.analyzeData(input_exp_file,species,platform,compendiumType,geneToReport=genesToReport,correlateAll=correlateAll,AdditionalParameters=fl,logTransform=logTransform)
                print 'MarkerFinder analysis complete'
                
                #markerfile="/Volumes/Pass/Final_scicgs/ExpressionOutput/MarkerFinder/MarkerGenes_correlations-ReplicateBased.txt"
                allgenesfile = root_dir+'/NMF-SVM/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt'
                markerfile = root_dir+'/NMF-SVM/ExpressionOutput/MarkerFinder/MarkerGenes_correlations-ReplicateBased.txt'
                guides=[]
                ### See if any unique genes are found in a cluster before using it for SVM
                guides,name,group=DetermineClusterFitness(allgenesfile,markerfile,input_exp_file,BinarizedOutput,rho_cutoff)
                counter=group
            else:
                if platform=="PSI":
                    rootdir,CovariateQuery=metaDataAnalysis.remoteAnalysis(species,FilteredEventAnnot,Metadata,'PSI',0.1,use_adjusted_p,0.05,Annotation)
                else:
                    rootdir,CovariateQuery=metaDataAnalysis.remoteAnalysis(species,processedInputExpFile,Metadata,'RNASeq',0.58,use_adjusted_p,0.05,Annotation)
                counter=1
                Guidedir=rootdir+CovariateQuery
                PSIdir=rootdir+'ExpressionProfiles'
    
                global upd_guides
                upd_guides=[]
                name=[]
                group=[]
                
                for filename in os.listdir(Guidedir):
                    if filename.startswith("PSI."):
                        Guidefile=os.path.join(Guidedir, filename)
                        psi=string.replace(filename,"PSI.","")
                    if filename.startswith("GE."):
                        Guidefile=os.path.join(Guidedir, filename)
                        psi=string.replace(filename,"GE.","")
                        PSIfile=os.path.join(PSIdir, psi)
                        omitcluster=FindTopUniqueEvents(Guidefile,psi,Guidedir)
                       
                        if omitcluster==0:
                            group.append(counter)
                            name.append(psi)
                            counter+=1
                upd_guides=[x for x in upd_guides if x != ""]
            upd_guides=guides
            upd_guides=list(set(upd_guides))
           
            scaling=True
            grplst=[]
            
            ############ Perform SVM classification to assign individual cells to valid-NMF clusters #############
            ### The below analysis is performed on the down-sampled expression file
            if counter>2:
                output_dir = root_dir+'/NMF-SVM'
                if os.path.exists(output_dir)==False:
                    export.createExportFolder(output_dir)
                
                #output_file = output_dir+'/SVMInput-Round'+str(iteration)+'.txt'
                #ExpandSampleClusters.filterRows(processedInputExpFile,output_file,filterDB=upd_guides,logData=False)
                
                if scaling:
                    output_fil=EventAnnot
                    output_file=output_dir+'/SVMInput-Round'+str(iteration)+'.txt'
                    #output_file1 = "/Users/meenakshi/Documents/Singlecellbest/exp.exp.CD34+.v5-log2_filtered.txt"
                    ExpandSampleClusters.filterRows(EventAnnot,output_file,filterDB=upd_guides,logData=False)
                else:
                    output_file = output_dir+'/SVMInput-Round'+str(iteration)+'.txt' 
                    ExpandSampleClusters.filterRows(processedInputExpFile,output_file,filterDB=upd_guides,logData=False)
                
                header=ExpandSampleClusters.header_file(output_file)
                
                print "Running SVM prediction for improved subtypes - Round"+str(iteration)
                ### Create teh training data for SVM
                train,null=ExpandSampleClusters.TrainDataGeneration(output_file,BinarizedOutput,name,scaling,exports=False,rootDir=root_dir)
            
                ### Determine the medoids (use medoids for SVM but centroids for clustering)
                grplst.append(group)
                
                Expand=False ### If Expand == True, use all down-sampled cells for classification rather than medoids (similar cellHarmony)
                if Expand==True: 
                    grplst=[]
                    group=ExpandSampleClusters.Findgroups(BinarizedOutput,name)
                    grplst.append(group)
                    
                ### Perform SVM
                ExpandSampleClusters.Classify(header,train,output_file,grplst,name,iteration,platform,output_dir,root_dir)
                
                ### Create a groups file for the downsampled (or original) file
                groupsfile = string.replace(originalExpFile,'exp.','groups.')
                groupsfile_downsampled = string.replace(processedInputExpFile,'exp.','groups.')
                finalgrpfile=root_dir+"/ICGS-NMF/FinalGroups.txt"
                if groupsfile_downsampled == groupsfile:
                    pass
                else:
                    export.customFileCopy(finalgrpfile,groupsfile_downsampled)
                export.customFileCopy(finalgrpfile,groupsfile[:-4]+'-ICGS.txt')
                export.customFileCopy(finalgrpfile,groupsfile[:-4]+'-markers.txt')
                from shutil import copyfile
                ### Don't overwrite the original groups
                updated_expfile = originalExpFile[:-4]+'-ICGS.txt'
                copyfile(originalExpFile, updated_expfile)
                if groupsfile_downsampled == groupsfile:
                    processedInputExpFile = updated_expfile
                groupsfile=groupsfile[:-4]+'-ICGS.txt'
                
                ### Identify markers for the our final un-ordered clusters (clustering will need to be run after this)
                markerFinder.analyzeData(processedInputExpFile,species,platform,compendiumType,geneToReport=genesToReport,correlateAll=correlateAll,AdditionalParameters=fl,logTransform=logTransform)
                allgenesfile=root_dir+"/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt"
                markergrps,markerlst=sortFile(allgenesfile,rho_cutoff)
                ### To plot the heatmap, use the MarkerFinder genes (function pulls those genes out)
                ExpandSampleClusters.filterRows(EventAnnot,processedInputExpFile[:-4]+'-markers.txt',filterDB=markerlst,logData=False) ### the processedInputExpFile is overwritten
                
                groupsdict=grpDict(groupsfile)
                matrix, column_header, row_header, dataset_name, group_db = clustering.importData(updated_expfile,geneFilter=markerlst)
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
    
                    #matrix = map(np.array, zip(*matrix)) ### coverts these to tuples
                    #column_header, row_header = row_header, column_header
                    finalOutputDir=root_dir+"/ICGS-NMF/"
                    #clustering.tSNE(np.array(matrix),column_header,dataset_name,group_db,display=False,showLabels=False,species=species,reimportModelScores=False)
                    try:
                        clustering.runUMAP(np.array(matrix),column_header,dataset_name,group_db,display=False,
                            showLabels=False,species=species,reimportModelScores=False,rootDir=root_dir,finalOutputDir=finalOutputDir)
                    except:
                        print traceback.format_exc()
        
                #clustering.tSNE(processedInputExpFile,group_db=groupsdict,display=True,showLabels=False,row_header=None,colorByGene=None,species=None,reimportModelScores=False)
                ##MV need to code 
                #Orderedfile,groupsdict=FindcentroidGroups(filtered,groupfile)
                SVMBinOutput=root_dir+"/NMF-SVM/SVMOutputs/round1SVC_Results_max.txt"
                SVMBinOutput_t=root_dir+"/NMF-SVM/SVMOutputs/round1SVC_Results_max_t.txt"
                import csv
                from itertools import izip
                a = izip(*csv.reader(open(SVMBinOutput,"rb"),delimiter='\t'))
                csv.writer(open(SVMBinOutput_t, "wb"),delimiter='\t').writerows(a)
                scaling=False ### will calculate centroids rather than medoids
                centroids,centroid_heatmap_input=ExpandSampleClusters.TrainDataGeneration(processedInputExpFile[:-4]+'-markers.txt',SVMBinOutput_t,name,scaling,exports=True,rootDir=root_dir)
                scaling=True
                
                graphic_links=[]
                row_method = "hopach"
                column_method="hopach"
                column_metric='cosine'
                row_metric='correlation'
                color_gradient = 'yellow_black_blue'
                transpose=False
                graphic_links = clustering.runHCexplicit(centroid_heatmap_input,graphic_links, row_method, row_metric, column_method,column_metric,color_gradient, transpose, display=False, Normalize=True)
                NMFSVM_centroid_cluster_dir=graphic_links[0][1][:-4]+'.txt'
                outputDir = root_dir+"/NMF-SVM/SVMOutputs"
                header=ExpandSampleClusters.header_file(NMFinput)
                status,graphic_links2=generateMarkerheatmap(processedInputExpFile[:-4]+'-markers.txt',output_file,NMFSVM_centroid_cluster_dir,groupsdict,markergrps,header,outputDir,root_dir,species,uniqueIDs)
                import shutil
                if status=='not-subsampled':
                    NMFSVM_centroid_cluster_dir=graphic_links2[0][1][:-4]
                    shutil.copy(NMFSVM_centroid_cluster_dir+'.txt',root_dir+"/ICGS-NMF/FinalMarkerHeatmap.txt")
                    shutil.copy(NMFSVM_centroid_cluster_dir+'.png',root_dir+"/ICGS-NMF/FinalMarkerHeatmap.png")
                    shutil.copy(NMFSVM_centroid_cluster_dir+'.pdf',root_dir+"/ICGS-NMF/FinalMarkerHeatmap.pdf")
                    shutil.copy(allgenesfile,root_dir+"/ICGS-NMF/MarkerGenes.txt")
                else:
                    NMFSVM_centroid_cluster_dir=graphic_links2[0][1][:-4]
                    NMFSVM_centroid_cluster_dir1=graphic_links2[1][1][:-4]
                    shutil.copy(NMFSVM_centroid_cluster_dir+'.txt',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_sampled.txt")
                    shutil.copy(NMFSVM_centroid_cluster_dir+'.png',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_sampled.png")
                    shutil.copy(NMFSVM_centroid_cluster_dir+'.pdf',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_sampled.pdf")
                    shutil.copy(NMFSVM_centroid_cluster_dir1+'.txt',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_all.txt")
                    shutil.copy(NMFSVM_centroid_cluster_dir1+'.png',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_all.png")
                    shutil.copy(NMFSVM_centroid_cluster_dir1+'.pdf',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_all.pdf")
                    shutil.copy(allgenesfile,root_dir+"/ICGS-NMF/MarkerGenes.txt")
                
                ### write final groups ordered
                #exportGroups(root_dir+"/ICGS-NMF/FinalMarkerHeatmap.txt",root_dir+"/ICGS-NMF/FinalGroups.txt",platform)
                
                if scaling:
                    flag=False
                    return flag,processedInputExpFile,EventAnnot,graphic_links3+graphic_links2
    
                header=Correlationdepletion.header_file(NMFResult)
                
                output_file=output_dir+'/DepletionInput-Round'+str(iteration)+".txt"
                sampleIndexSelection.filterFile(processedInputExpFile[:-4]+'-markers.txt',output_file,header)
                print "Running Correlation Depletion - Round"+str(iteration)
                commonkeys,count=Correlationdepletion.FindCorrelations(NMFResult,output_file,name)
                Depleted=Correlationdepletion.DepleteSplicingevents(commonkeys,output_file,count,processedInputExpFile)
                processedInputExpFile=Deplete
                flag=True
                
            else:
                try:
                    print "Running K-means analyses instead of NMF - Round"+str(iteration)
                    header=[]
                    header=Kmeans.header_file(Guidefile_block)
                    Kmeans.KmeansAnalysis(Guidefile_block,header,processedInputExpFile,iteration)
                    flag=False
                except Exception:
                    flag=False
                
        else:
            if Rank==1:
                try:
                    print "Running K-means analyses instead of NMF - Round"+str(iteration)
                    header=[]
                    header=Kmeans.header_file(Guidefile_block)
                    Kmeans.KmeansAnalysis(Guidefile_block,header,processedInputExpFile,iteration)
            
                    flag=False
                except Exception:
                    flag=False
            else:
                flag=False
         
        return flag,processedInputExpFile,EventAnnot
    except:
        print traceback.format_exc()
        print 'WARNING!!!! Error encountered in the NMF ICGS analysis... See the above report.'
        flag=False
        return flag,processedInputExpFile,EventAnnot,graphic_links3

def exportGroups(cluster_file,outdir,platform):
    lineNum=1
    for line in open(cluster_file,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if lineNum==1: names = t[2:]; lineNum+=1
        elif lineNum==2: clusters = t[2:]; lineNum+=1
        else: break

    out_obj = export.ExportFile(outdir)
    for name in names:
        cluster = clusters[names.index(name)]
        if platform == 'RNASeq':
            if 'junction_quantification' not in name and '.bed' not in name:
                name = name+'.bed'
            elif 'junction_quantification.txt' not in name and '.txt' not in name and '.bed' not in name:
                name = name+'.txt'
        if ':' in name:
            group,name = string.split(name,':')
            if cluster=='NA': cluster = group
        out_obj.write(name+'\t'+cluster+'\t'+cluster+'\n')
    out_obj.close()
    
def runICGS_NMF(inputExpFile,scaling,platform,species,gsp,enrichmentInput='',dynamicCorrelation=True):
    """ Runs downsampling analysis and prepares files for ICGS. After running ICGS, peform enrichment analyses """
    
    print 'Running ICGS-NMF'
    
    ### Find the parent dir of the output directory (expression file from the GUI will be stored in the output dir [ExpressionInput])
    root_dir = export.findParentDir(inputExpFile)[:-1]
    if 'ExpressionInput' in inputExpFile:
        root_dir = export.findParentDir(root_dir)
    exp_file_name = export.findFilename(inputExpFile)
    
    ### Assign the expression filename (for single-cell RNA-Seq rather than splicing)
    if 'exp.' not in exp_file_name:
        exp_file_name = 'exp.' + exp_file_name
    
    ########## Perform Downsampling for large datasets ##########
    
    ### Use dispersion (variance by mean) to define initial variable genes
    inputExpFileVariableGenesDir,n=hgvfinder(inputExpFile) ### returns filtered expression file with 500 variable genes
    if n>2500 and scaling:
        if n>25000: ### For extreemly large datasets, UMAP is used as a preliminary downsampling before pagerank
            inputExpFileScaled=inputExpFile[:-4]+'-UMAP-downsampled.txt'
            ### UMAP and Louvain clustering for down-sampling from >25,000 to 10,000 cells
            sampmark=UMAPsampling(inputExpFileVariableGenesDir) ### returns list of UMAP downsampled cells
            ### Filer the original expression file using these downsampled cells
            sampleIndexSelection.filterFile(inputExpFile,inputExpFileScaled,sampmark)
            ### Use dispersion (variance by mean) to define post-UMAP/Louvain selected cell variable genes
            inputExpFileVariableGenesDir,n=hgvfinder(inputExpFileScaled) ### returns filtered expression file with 500 variable genes
            ### Run PageRank on the UMAP/Louvain/dispersion downsampled dataset
            sampmark=PageRankSampling(inputExpFileVariableGenesDir)
        else:
            ### Directly run PageRank on the initial dispersion based dataset
            sampmark=PageRankSampling(inputExpFileVariableGenesDir)
        
        ### Write out final downsampled results to a new file
        output_dir = root_dir+'/ExpressionInput'
        try: export.createExportFolder(output_dir)
        except: pass ### Already exists
        processedInputExpFile = root_dir+'/ExpressionInput/'+exp_file_name[:-4]+'-PageRank-downsampled.txt' ### down-sampled file
        sampleIndexSelection.filterFile(inputExpFile,processedInputExpFile,sampmark)
    else:
        output_dir = root_dir+'/ExpressionInput'
        try: export.createExportFolder(output_dir)
        except: pass ### Already exists
        if platform == 'PSI':
            ### The PSI file name by default is not informative
            processedInputExpFile=output_dir+"/exp.spliceICGS-input.txt"
            export.customFileCopy(inputExpFile,processedInputExpFile)
        elif 'ExpressionInput' not in inputExpFile:
            processedInputExpFile = root_dir+'/'+exp_file_name
            export.customFileCopy(inputExpFile,processedInputExpFile)
        else: processedInputExpFile = inputExpFile
  
    flag=True
    iteration=1 ### Always equal to 1 for scRNA-Seq but can increment for splice-ICGS

    ### Recursively run ICGS with NMF
    flag,processedInputExpFile,inputExpFile,graphic_links3=CompleteICGSWorkflow(root_dir,processedInputExpFile,
            inputExpFile,iteration,gsp.RhoCutoff(),dynamicCorrelation,platform,species,scaling,gsp)
  
    if platform == 'PSI':
        output_dir = root_dir+'/SVMOutputs'
        Combinedres=MergeResults(output_dir)
    
    mutlabels={}
    if enrichmentInput!='':
        print "Running Mutation Enrichment Analyses"
        Expand="yes"
        mutdict=defaultdict(list)
        header=ME.header_file(enrichmentInput)
      
        mutdict=ME.findsiggenepermut(enrichmentInput)
      
        mutlabels=ME.Enrichment(Combinedres,mutdict,enrichmentInput,Expand,header)
    if platform == 'PSI':
        print "Generating the final consolidated results"
        Orderedheatmap.Classify(Combinedres,mutlabels,dire)
        Orderedheatmap.Classify(Combinedres,mutlabels,dire,False)

    print "successfully completed"
    return graphic_links3
    
if __name__ == '__main__':
    """
    processedInputExpFile="/Volumes/Pass/ICGS2_testrun/ExpressionInput/exp.input.txt"
    
    matrix, column_header, row_header, dataset_name, group_db =clustering.importData(processedInputExpFile)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
    
        #matrix = map(np.array, zip(*matrix)) ### coverts these to tuples
        #column_header, row_header = row_header, column_header
        directory=export.findParentDir(export.findParentDir(processedInputExpFile)[:-1])+"ICGS-NMF/"
        #clustering.tSNE(np.array(matrix),column_header,dataset_name,group_db,display=False,showLabels=False,species="Hs",reimportModelScores=False)
        clustering.umap(np.array(matrix),column_header,dataset_name,group_db,display=False,showLabels=False,species="Hs",reimportModelScores=False,directory=directory)
    sys.exit()
    """
    import getopt
    rho_cutoff=0.2
    dynamicCorrelation="optimize" 
    Mutationref=""
    platform="RNASeq"
    scaling=True
    species="Hs"
    row_method = 'hopach'
    column_method = 'hopach'
    row_metric = 'correlation'
    column_metric = 'euclidean'
    color_gradient = 'yellow_black_blue'
    contrast=3
    vendor = "RNASeq"
    GeneSelection = ''
    PathwaySelection = ''
    GeneSetSelection = 'None Selected'
    excludeCellCycle = False
    
    restrictBy = 'protein_coding'
    #restrictBy = 'None'
    featurestoEvaluate = 'Genes'
    ExpressionCutoff = 0
    CountsCutoff = 0
    if platform=="PSI":
        FoldDiff = 1.2
    else:
        FoldDiff=4.0
        ExpressionCutoff = 1
    SamplesDiffering = 4
    JustShowTheseIDs=''
    removeOutliers = False
    PathwaySelection=[]
    array_type="RNASeq"
    rho_cutoff=rho_cutoff

    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['Input=','rho=','dynamicCorrelation=','Mutationref=','platform=','scaling=','species=','ExpressionCutoff=','CountsCutoff=','FoldDiff=','SamplesDiffering=','removeOutliers=','featurestoEvaluate=','restrictBy=','excludeCellCycle=','column_metric=','column_method=','row_method=','row_metric'])
        for opt, arg in options:
            #if opt == '--processedInputExpFile': processedInputExpFile=arg
            if opt=='--Input':EventAnnot=arg # input file
            if opt=='--rho':rho_cutoff=arg # rho cutoff
            if opt=='--dynamicCorrelation':
                if string.lower(dynamicCorrelation) == 'true' or string.lower(dynamicCorrelation) == 'optimize':
                    dynamicCorrelation=True #constant using the provided correlation,iteratively optimize correlation cutoff"
                else:
                    dynamicCorrelation=False
            if opt=='--Mutationref':Mutationref=arg #reference file provided for enrichment (format groups file)
            if opt=='--platform':platform=arg
            if opt=='--scaling':scaling=arg # True to scale for large datasets, False run with all samples
            if opt=='--species':species=arg
            if opt=='--ExpressionCutoff':ExpressionCutoff=arg
            if opt=='--CountsCutoff':CountsCutoff=arg
            if opt=='--FoldDiff':FoldDiff=arg
            if opt=='--SamplesDiffering':SamplesDiffering=arg
            if opt=='--removeOutliers':removeOutliers=arg
            if opt=='--featurestoEvaluate':featurestoEvaluate=arg
            if opt=='--restrictBy':restrictBy=arg
            if opt=='--column_metric':column_metric=arg
            if opt=='--column_method':column_method=arg
            if opt=='--row_method':row_method=arg
            if opt=='--row_metric':row_metric=arg
            
    gsp = UI.GeneSelectionParameters(species,array_type,vendor)
    gsp.setGeneSet(GeneSetSelection)
    gsp.setPathwaySelect(PathwaySelection)
    gsp.setGeneSelection(GeneSelection)
    gsp.setJustShowTheseIDs(JustShowTheseIDs)
    gsp.setNormalize('median')
    gsp.setSampleDiscoveryParameters(ExpressionCutoff,CountsCutoff,FoldDiff,SamplesDiffering,removeOutliers,featurestoEvaluate,restrictBy,excludeCellCycle,column_metric,column_method,rho_cutoff) 
    
    runICGS_NMF(EventAnnot,scaling,dynamicCorrelation,platform,species,Mutationref,gsp)
