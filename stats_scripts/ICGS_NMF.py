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
warnings.filterwarnings('ignore')

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
                community=True
            else:
                community=False
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
          
            if (qq.startswith("rpl") or qq.startswith("rps") or qq.startswith("mt-") or qq.startswith("ig")) and community:
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

def community_sampling(inputfile):
    """ This function performs downsampling of the input data using networkx to identify
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
    

    from annoy import AnnoyIndex
    t=AnnoyIndex(nm,metric="euclidean")
    for i in range(nn):
      try:  t.add_item(i,X[i])
      except Exception: print i
    t.build(100)
    ### t.save('ICGS.ann')
    ### u=AnnoyIndex(nm,metric="euclidean")
    diclst={}
    #### u.load('ICGS.ann')
    #n1=25
    print "creating graphs"
    for i in range(nn):
    #ind = tree.query([Xtemp[i]],k=10,return_distance=False,dualtree=True)
        ind=t.get_nns_by_item(i,10)
        diclst[i]=ind

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
    
def PageRankSampling(inputfile,downsample_cutoff):
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
        jj=downsample_cutoff
        if iq+24999>n:
            j=n-iq
        else:
            j=24999
        jj=int(float(j+1)/4.0)
        jj=downsample_cutoff
        #if jj<downsample_cutoff and n<3000:
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
        
        t.save('ICGS.ann')
       
        u=AnnoyIndex(nm)
        u.load('ICGS.ann')
       
        #tree = KDTree(X, leaf_size=10, metric='euclidean')
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
    if len(sampmark)>downsample_cutoff:
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
    if (common_geneIDs+2)<marker_count:
        print 'WARNING... only',common_geneIDs, 'out of', marker_count, 'gene IDs matched after conversion.'
        
    for key in markerdict:
        countr=1
        if len(markerdict[key])>=2 and key in nametemp:
            name.append(key+"_vs_Others.txt")
            group.append(counter)
      
            for i,j in markerdict[key] :
                #if countr<30:
                    upd_guides.append(i)
                    countr+=1
            counter+=1

    return upd_guides,name,group

def sortFile(allgenesfile,rho_cutoff,name):
    
    markergenes={}
    val=[]
    header=True
    namelst=[]
    for i in range(len(name)):
        s=string.split(name[i],"_")[0]
        namelst.append(s)
    
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
                    
                    if "V" in t[4] and t[4] in namelst:		
                        t[4]=string.replace(t[4],"V","")
                        values.append(t[4])
                        
                    else:
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
    samples3=[]
    samples_all=[]
    samples2_all=[]
    groups_list=[]
    groups_list_all=[]
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
                samples.append(groupsdict[sampleOrder[i]][j])
                groups_list.append(sampleOrder[i])
                samples2.append(groupsdict[sampleOrder[i]][j])
                samples3.append(sampleOrder[i]+':'+groupsdict[sampleOrder[i]][j])
    for i in range(len(sampleOrder)):
        for j in range(len(markergrps[sampleOrder[i]])):
            uid = markergrps[sampleOrder[i]][j]
            genes.append(uid)
            if uid in uniqueIDs:
                symbol = uniqueIDs[uid]
            else:
                symbol = uid
            genes2.append((sampleOrder[i],uid))
            
    MF_subsampled_export = outputDir+'/'+'MarkerFinder-subsampled-ordered.txt'
    exportnam=open(MF_subsampled_export,"w")     
    
    exportnam.write(string.join(['UID','row_clusters-flat']+samples3,'\t')+'\n')
    exportnam.write(string.join(['column_clusters-flat','']+groups_list,'\t')+'\n')
    i=0

    for i in range(len(genes)):
        exportnam.write(genes2[i][1]+"\t"+genes2[i][0])
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
    graphic_links = clustering.runHCexplicit(MF_subsampled_export,graphic_links, row_method, row_metric, column_method,column_metric,color_gradient, gsp, display=False, Normalize=True,contrast=5)
    graphic_links[-1][0] = MF_subsampled_export
    
    if len(samp)>len(header1):
        MF_all_export = outputDir+'/'+'MarkerFinder-Allsamples-ordered.txt'
        all_cells_export=open(MF_all_export,"w")
        
        for i in range(len(sampleOrder)):
            for j in range(len(groupsdict[sampleOrder[i]])):
                samples_all.append(sampleOrder[i]+":"+groupsdict[sampleOrder[i]][j])
                groups_list_all.append(sampleOrder[i])
                samples2_all.append(groupsdict[sampleOrder[i]][j])
                
        all_cells_export.write(string.join(['UID','row_clusters-flat']+samples_all,'\t')+'\n')
        all_cells_export.write(string.join(['column_clusters-flat','']+groups_list_all,'\t')+'\n')
    
        for i in range(len(genes)):
            all_cells_export.write(genes2[i][1]+"\t"+genes2[i][0])
            for j in range(len(samples_all)):
                
                all_cells_export.write("\t"+matrix[genes[i],samples2_all[j]])
            all_cells_export.write("\n")
        all_cells_export.close()
        graphic_links = clustering.runHCexplicit(MF_all_export,graphic_links, row_method, row_metric, column_method,column_metric,color_gradient, gsp, display=False, Normalize=True,contrast=5)
        graphic_links[-1][0] = MF_all_export
        status = 'subsampled'
    else:
        status = 'not-subsampled'

    return status, graphic_links

def callICGS(processedInputExpFile,species,rho_cutoff,dynamicCorrelation,platform,gsp):
    
    #Run ICGS recursively to dynamically identify the best rho cutoff
    graphic_links3,n = RNASeq.singleCellRNASeqWorkflow(species,platform,processedInputExpFile,mlp,dynamicCorrelation, rpkm_threshold=0, parameters=gsp)
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
        graphic_links3 = RNASeq.singleCellRNASeqWorkflow(species, 'exons', processedInputExpFile,mlp, rpkm_threshold=0, parameters=gsp)
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
            k=estimateK(NMFinput)
            Rank=k*2
            if Rank>2 and platform=='PSI':
                Rank=30
            if Rank<5 and platform!='PSI':
                Rank=10
            ### This function prepares files for differential expression analsyis (MetaDataAnalysis), MarkerFinder
            filteredInputExpFile = string.replace(processedInputExpFile,'exp.','filteredExp.')
            
            if '-OutliersRemoved' in Guidefile:
                filteredInputExpFile = string.replace(filteredInputExpFile,'.txt','-OutliersRemoved.txt')

            try: NMFResult,BinarizedOutput,Metadata,Annotation=NMF_Analysis.NMFAnalysis(filteredInputExpFile,NMFinput,Rank,platform,iteration)
            except:
                try:
                    Rank=k*1.5
                    NMFResult,BinarizedOutput,Metadata,Annotation=NMF_Analysis.NMFAnalysis(filteredInputExpFile,NMFinput,Rank,platform,iteration)
                except:
                    Rank=k
                    NMFResult,BinarizedOutput,Metadata,Annotation=NMF_Analysis.NMFAnalysis(filteredInputExpFile,NMFinput,Rank,platform,iteration)
        else:
            Rank=k
            
            print "Running NMF analyses for dimension reduction using "+str(Rank)+" ranks - Round"+str(iteration)
            print "The number target number of clusters (k/rank) is:",k
            filteredInputExpFile = string.replace(processedInputExpFile,'exp.','filteredExp.')
            
            if '-OutliersRemoved' in Guidefile:
                filteredInputExpFile = string.replace(filteredInputExpFile,'.txt','-OutliersRemoved.txt')
            try:
                NMFResult,BinarizedOutput,Metadata,Annotation=NMF_Analysis.NMFAnalysis(filteredInputExpFile,NMFinput,Rank,platform,iteration)
            except Exception:
                "Exception, choose a lower k value."
        if Rank>1:

            
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
                counter=len(group)
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
                markergrps,markerlst=sortFile(allgenesfile,rho_cutoff,name)
                if len(markergrps)!=len(name):
                    allgenesfile1 = root_dir+'/NMF-SVM/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt'
                    markergrps,markerlst=sortFile(allgenesfile1,rho_cutoff,name)
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
                    NMFSVM_centroid_cluster_graphics_dir=graphic_links2[0][1][:-4]
                    NMFSVM_centroid_cluster_dir=graphic_links2[0][0][:-4]
                    shutil.copy(NMFSVM_centroid_cluster_dir+'.txt',root_dir+"/ICGS-NMF/FinalMarkerHeatmap.txt")
                    shutil.copy(NMFSVM_centroid_cluster_graphics_dir+'.png',root_dir+"/ICGS-NMF/FinalMarkerHeatmap.png")
                    shutil.copy(NMFSVM_centroid_cluster_graphics_dir+'.pdf',root_dir+"/ICGS-NMF/FinalMarkerHeatmap.pdf")
                    shutil.copy(allgenesfile,root_dir+"/ICGS-NMF/MarkerGenes.txt")
                else:
                    NMFSVM_centroid_cluster_graphics_dir=graphic_links2[0][1][:-4]
                    NMFSVM_centroid_cluster_dir=graphic_links2[0][0][:-4]
                    NMFSVM_centroid_cluster_graphics_dir2=graphic_links2[1][1][:-4]
                    NMFSVM_centroid_cluster_dir2=graphic_links2[1][0][:-4]
                    
                    NMFSVM_centroid_cluster_dir=graphic_links2[0][0][:-4]
                    NMFSVM_centroid_cluster_dir1=graphic_links2[1][0][:-4]
                    shutil.copy(NMFSVM_centroid_cluster_dir+'.txt',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_sampled.txt")
                    shutil.copy(NMFSVM_centroid_cluster_graphics_dir+'.png',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_sampled.png")
                    shutil.copy(NMFSVM_centroid_cluster_graphics_dir+'.pdf',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_sampled.pdf")
                    shutil.copy(NMFSVM_centroid_cluster_dir2+'.txt',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_all.txt")
                    shutil.copy(NMFSVM_centroid_cluster_graphics_dir2+'.png',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_all.png")
                    shutil.copy(NMFSVM_centroid_cluster_graphics_dir2+'.pdf',root_dir+"/ICGS-NMF/FinalMarkerHeatmap_all.pdf")
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
                if iteration<2:
                    gsp.setK(k)
                    iteration=1
                    flag,processedInputExpFile,inputExpFile,graphic_links3=CompleteICGSWorkflow(root_dir,processedInputExpFile,inputExpFile,iteration,gsp.RhoCutoff(),dynamicCorrelation,platform,species,scaling,gsp)
                else:
                    "No groups found!!! Re-analyze the data with a small k"
                #try:
                 #   print "Running K-means analyses instead of NMF - Round"+str(iteration)
                  #  print "Extremely sparse data!! choose a small k"
                  #  header=[]
                   # header=Kmeans.header_file(Guidefile_block)
                    #Kmeans.KmeansAnalysis(Guidefile_block,header,processedInputExpFile,iteration)
                    #flag=False
                #except Exception:
                 #   flag=False
                
        else:
            if Rank==1:
                try:
                    print "Running K-means analyses instead of NMF - Round"+str(iteration)
                    print "Extremely sparse data!! choose a small k"
                    header=[]
                    header=Kmeans.header_file(Guidefile_block)
                    Kmeans.KmeansAnalysis(Guidefile_block,header,processedInputExpFile,iteration)
            
                    flag=False
                except Exception:
                    flag=False
            else:
                flag=False
         
        return flag,processedInputExpFile,EventAnnot,graphic_links3
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
    """ Export the filtered expression file then run downsampling analysis and prepares files for ICGS. After running ICGS, peform enrichment analyses """
    
    try: downsample_cutoff = gsp.DownSample()
    except: downsample_cutoff = 2500
    print 'DownSample threshold =',downsample_cutoff, 'cells'
    print 'Filtering the expression dataset (be patient).',
    print_out, inputExpFile = RNASeq.singleCellRNASeqWorkflow(species,platform,inputExpFile,mlp,rpkm_threshold=0,parameters=gsp,reportOnly=True)
    
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
    
    if n>downsample_cutoff and scaling:
       
        if n>25000: ### For extreemly large datasets, Louvain is used as a preliminary downsampling before pagerank
            print 'Performing Community Clustering...'
            inputExpFileScaled=inputExpFile[:-4]+'-Louvain-downsampled.txt'
            ### Louvain clustering for down-sampling from >25,000 to 10,000 cells
            sampmark=community_sampling(inputExpFileVariableGenesDir) ### returns list of Louvain downsampled cells
            ### Filer the original expression file using these downsampled cells
            sampleIndexSelection.filterFile(inputExpFile,inputExpFileScaled,sampmark)
            ### Use dispersion (variance by mean) to define post-Louvain selected cell variable genes
            inputExpFileVariableGenesDir,n=hgvfinder(inputExpFileScaled) ### returns filtered expression file with 500 variable genes
            ### Run PageRank on the Louvain/dispersion downsampled dataset
            sampmark=PageRankSampling(inputExpFileVariableGenesDir,downsample_cutoff)
        else:
            ### Directly run PageRank on the initial dispersion based dataset
            sampmark=PageRankSampling(inputExpFileVariableGenesDir,downsample_cutoff)
       
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
    processedInputExpFile='/Users/saljh8/Desktop/DemoData/ICGS-Mm/ExpressionInput/exp.Bladder-10X_P4_3-VarGenes-ICGS-markers.txt'
    output_file='/Users/saljh8/Desktop/DemoData/ICGS-Mm//NMF-SVM/SVMInput-Round1.txt'
    NMFSVM_centroid_cluster_dir='/Users/saljh8/Desktop/DemoData/ICGS-Mm/NMF-SVM/centroids//DataPlots/Clustering-exp.MF-hierarchical_cosine_correlation.txt'
    groupsdict={'1': ['GTAACTGCAGCATACT', 'GCTGCTTTCAAAGTAG', 'TAAGAGAAGCTATGCT', 'CGGAGTCAGTACGCCC', 'CTAGCCTAGCTGCGAA', 'GAACATCAGTGTTAGA', 'TTGGAACGTTAGGGTG', 'TGGCTGGGTCGCTTCT', 'ACTTTCAAGATGGGTC', 'CGTGAGCGTGCACTTA', 'CACACCTAGCCCAATT', 'GTGCAGCTCAGTTAGC', 'CCGTGGAGTGAGGGTT', 'CTGCGGAGTATGCTTG', 'CACACCTGTTTGACAC', 'CTGTGCTCAATGGACG', 'AAGACCTAGGTGGGTT', 'ACAGCTATCATACGGT', 'CTGCCTACAGCCACCA', 'GCGACCAAGGGATCTG', 'CTTGGCTCACGAAACG', 'AGCGGTCTCTGAGTGT', 'GTACTCCTCGCCTGTT', 'TTTACTGCAGGACGTA', 'GAATAAGGTCGGCATC', 'AGTGAGGTCAGCATGT', 'GCATGATCATATGGTC', 'TGCGGGTTCTTTACGT', 'TTGGCAACATGATCCA', 'GTTCTCGGTGCTTCTC', 'TAAGCGTAGAAACCGC', 'CCGGTAGAGGCAGTCA', 'CGTAGGCCACGGACAA', 'GGGCATCTCCGTTGTC', 'ATCGAGTTCTCCAGGG', 'TAGGCATTCTATCCTA', 'GTTTCTAAGGACTGGT', 'TCCACACCAGCTGCTG', 'CTGCTGTTCGTTACGA', 'AACTCTTAGTTGCAGG', 'TTTACTGCACACGCTG', 'GTTAAGCCACCAGGTC', 'AGGTCCGAGGACACCA', 'CACATAGCAGATGGCA', 'CGGACTGTCTTGTATC', 'TGGACGCTCTGCAGTA'], '0': ['GTCACAATCTACTCAT', 'TCTTCGGGTTTAGGAA', 'GAAATGACAATCCAAC', 'TACGGATAGGTACTCT', 'AGGGTGACAGTCGTGC', 'GTAACGTTCAGGCGAA', 'GTCACAATCACTTCAT', 'AGGGAGTCAATCGGTT', 'GCGAGAATCCCTCTTT', 'CCCAGTTCACATGACT'], '3': ['GTGTGCGTCAGTTCGA', 'TTGCCGTCACGCCAGT', 'AAGACCTAGATCCGAG', 'CAGCATAAGATGGGTC'], '2': ['ACGGAGATCGGAATCT', 'AGACGTTTCTTACCGC', 'GGACAGACATTCACTT', 'CATTCGCTCTCGTATT', 'AGTGTCACAGTATCTG', 'GTCGGGTAGGCCGAAT', 'GTCATTTGTCTTGCGG', 'CAGCTGGTCCTTTACA', 'AACTCTTTCATAACCG', 'CGGGTCATCGTGGACC', 'TAGCCGGAGGACTGGT', 'CTCGAGGAGGCGCTCT', 'CAACCAAAGACGCACA', 'CTCTGGTAGTGTACCT', 'AACTCCCGTCGGGTCT', 'CACTCCACACATCCAA', 'TACTTGTGTGTTCTTT', 'GACGGCTAGTCGTTTG', 'GGTGCGTAGAGGTACC', 'TGAGAGGAGACATAAC', 'GACGTTAGTAGGCTGA', 'AGAGCTTTCTGTTGAG', 'CTCTAATAGTGAATTG', 'CCTTCGATCTCCGGTT', 'ATTTCTGGTCAGGACA', 'GTACGTACACACCGAC', 'GGCTCGAAGCGTGAGT', 'GAAGCAGCACCGATAT', 'CTCGAGGAGGGTTTCT', 'GTGCAGCGTGGTACAG', 'ATCACGACAACAACCT', 'TGGGAAGAGTAGGCCA', 'AACCGCGTCCAACCAA'], '5': ['AAGTCTGAGATAGTCA', 'CACCAGGCAAGGCTCC', 'ACATACGCAGCTCCGA', 'TTGTAGGCATCCGGGT', 'CACATTTTCTGGTATG', 'GATCGCGTCACCCTCA', 'GAACCTACATTAGGCT', 'GTGGGTCAGATGTGGC', 'CTCATTAGTGAAAGAG', 'CGATTGATCTAACTTC', 'TCTTTCCTCCGCATAA'], '4': ['AGTCTTTAGCTAGTGG', 'AGTGAGGAGTGTACTC', 'CCTACACAGGGCATGT', 'AGGGATGAGTCGCCGT', 'TGAGAGGAGCGATATA', 'TTGGCAAGTCCGTTAA', 'TGGACGCTCCGATATG', 'GGCGTGTCATTGGGCC', 'CAGAATCTCTCAAGTG', 'GACTAACCATCGACGC', 'GTATCTTTCGCCATAA', 'GCAAACTAGAGCCCAA', 'AAAGTAGAGATGCCAG', 'GTCTCGTCAAACTGTC', 'ACTGAACCATCCGGGT', 'ACGGGCTGTCTGATCA', 'CTCGTACTCAATACCG', 'CCGGTAGCAGTGGGAT'], '7': ['GTAGTCAAGACGACGT', 'ATTATCCAGTAAGTAC', 'TTTGTCAGTTGCGTTA', 'GATTCAGCAGTGGAGT', 'TTGAACGTCTCTTATG', 'CTCTAATCACGGCGTT', 'TTCTACAAGGCAGTCA', 'GACGTTATCAGGTTCA', 'ATTTCTGCATGTTGAC', 'CTGAAACGTAGTACCT', 'TCAACGACAATGGTCT', 'AAGGAGCGTGCAACTT', 'GAGGTGATCCAGTATG', 'GCGAGAACACATGGGA', 'ACGAGCCAGATAGCAT', 'CACCTTGAGGACCACA', 'GGAATAATCATTCACT', 'CGTCCATAGTCCTCCT', 'CCAGCGATCGCTGATA', 'ACCAGTACACGTCTCT', 'CTAGAGTGTCATATCG', 'GTATCTTGTAGAGGAA', 'CGAGAAGAGGACTGGT', 'AGCTTGAAGACTAAGT'], '6': ['GCGACCACACCAACCG', 'TGAGCATAGTGGAGTC', 'CCGTTCAGTTCCCTTG']}
    markergrps={'1': ['Krt15', 'Igfbp2', 'Lmo1', 'Krt5', 'Gsta3', 'Acaa1b', 'Gstm1', 'Cbr2', 'Ly6d', 'Fxyd3', '2200002D01Rik', 'Gsto1', 'Klf5', 'Gsta4', 'Wfdc2', 'Aqp3', 'Perp', 'Foxq1', 'Uba52', 'Mgst3', '1190003J15Rik', 'Sfn', 'Spint2', 'Aldh3a1', 'Sdc1', 'Rab25', 'Krt19', 'Avpi1', 'Krt7', 'Ezr', 'Mal', 'Sprr1a', 'Gsdmc2', 'Fau', 'Krt4', 'Ctse', 'Krt18', 'Gclc', 'Hilpda', 'Cyb5', 'S100a6', 'Lgals3', 'Cldn4', 'Krt8', 'Akr1b8', 'Dnmt3l', 'Cldn7', 'Tnfaip8', 'Cox4i1', 'Fmo5', 'Mgst2', 'Mgst1', '1810011O10Rik', 'Nqo1', 'Upk1a', 'Kcnk1', 'Foxa1', 'S100a14', 'Ptprf', 'Net1'], '0': ['Cck', 'Dpep1', 'Apoe', 'Naalad2', 'Efemp1', 'Lmo4', 'Hsd11b1', 'Cfl2', 'Aldh2', 'Hjurp', 'Wsb1', 'Tmem158', 'Pdia6', 'Psmb3', 'Mgat2', 'Pex5l', 'Arsb', 'Cyba', 'Fam114a1', 'Triobp', 'Sf3b2', 'Grn', 'Samhd1', 'Mll5', '1110007C09Rik', 'Sfxn3', 'Sphk2', 'Tln1', 'Ext2', 'Sarnp', 'Lrpap1', 'Topors', 'Scyl1', 'Fuca1', 'Aida', 'Bscl2', 'Ppp1r10', 'Med28', 'Rap2a', 'St13', 'Fos', 'Fam76a', 'Wtap', 'Pdlim5', 'Ypel3', 'Myd88', 'Sbds', 'Gng10', 'Sgce', 'Angptl2', 'H2afv', 'Klhdc2', 'Rnh1', 'Larp7', 'Ggnbp2', 'Pcna', 'Sqrdl', 'Sdf4', 'Fam46a', 'Osbpl9'], '3': ['Cd52', 'Coro1a', 'Laptm5', 'Cytip', 'Srgn', 'H2-Aa', 'Ctss', 'H2-Eb1', 'H2-DMa', 'Cd74', 'H2-Ab1', 'Il1b', 'Arl4c', 'Psmb8', 'Arhgdib', 'Ckb', 'Vasp', 'Stk17b', 'Efhd2', 'Cotl1', 'H2afz', 'H2-D1', 'Syngr2', 'Tnfaip3', 'Ifngr1', 'Tgfb1', 'Sema4a', 'Wnk1', 'Actr3', 'Sh3bgrl3', 'Pim1', 'D16Ertd472e', 'Hmgb2', 'Fam107b', 'H2afy', 'Actg1', 'Ehd1', 'Plscr1', 'Arpc5', 'Eif4e', 'Psap', 'H3f3a', 'Tob2', 'Npm1', 'Atp6v0e', 'Rnaseh2c', 'Dek', 'G3bp1', 'Ube2l3', 'Arpc3', 'Psma5', 'Pole4', 'Picalm', 'Cxcl16', 'Hpcal1', 'Cmpk1', 'Bax', 'Clic1', 'Stx7', 'Cct4'], '2': ['Car3', 'Tnc', 'Tagln', 'Myl9', 'Acta2', 'Hhip', 'Gpx3', 'Rbp4', 'Col6a3', 'Gas6', 'Cpz', 'Spon1', 'Dkk2', 'Tgfbi', 'Col4a1', 'Aldh1a2', 'Actg2', 'Ppic', 'Thbs2', 'Cxcl14', '3632451O06Rik', 'Col12a1', 'Col11a1', 'Mylk', 'Tpm2', 'Srpx', 'Smoc2', 'Mmp2', 'Crispld2', 'Wfdc1', 'Leprel2', 'Tmem119', 'Prelp', 'Mdk', 'Creb3l1', 'Spon2', 'Emilin1', 'Col4a2', 'Sparc', 'Kcnk2', 'Mfap2', 'Col6a2', 'Ctsl', 'Mmp14', 'Mfap4', 'Csrp1', 'Grem2', 'Vcam1', 'Cald1', 'Cd63', 'Fth1', 'Cpxm1', 'Ctsk', 'Mgp', 'Il4ra', 'Col6a1', 'Chpf', 'F2r', 'Gadd45b', 'Avpr1a'], '5': ['Psca', 'Upk1b', 'Sytl2', 'Upk3a', 'Spint1', 'Sprr2b', 'Cxadr', 'Upk2', 'Oit1', 'Snhg11', 'Trim29', 'Nipal1', 'Tmprss13', 'Sprr2g', 'Serinc2', 'Rab27b', 'Dsg2', 'Mpzl2', 'Far1', 'Tmem123', 'Klhdc3', 'Capn5', 'Atp6v0a1', 'Fdps', 'Ppargc1a', 'Rnf4', 'Ikzf2', 'Senp2', 'Ppfibp2', 'Ccrn4l', 'Gstm4', 'Efnb2', 'Wasl', 'D15Ertd621e', 'Sgpl1', 'Pttg1ip', 'Socs7', 'Sart1', 'Hsd17b12', 'Pnpla2', 'H2-K1', 'Trim25', 'Cd82', 'Mapk1', 'Srrm2', 'Rnf40', 'Nsun2', 'Smarca4', 'Erdr1', 'Pnn', 'Ergic1', 'Cd55', 'Rtn4', 'D17Wsu92e', 'Ctsh', 'Wdr1', 'Ugcg', 'Actn4', 'Plat', 'Dnajb6'], '4': ['Scd1', 'Ftl1', 'Junb', 'Gfpt2', 'Gda', 'Egr1', 'Galnt1', 'Tnxb', 'Has1', 'Fosb', 'Fxyd1', 'Akap12', 'Sepp1', 'Ddr2', 'Ifrd1', 'Hspb8', 'Gpx1', 'Tpm3', 'Emp1', 'Zfp36', 'Sema3c', 'Sphk1', 'Tmsb10', 'Gpr153', 'Hn1', 'Lhfp', 'Hspg2', 'Crip1', 'Socs3', 'Cxcl1', 'Gsk3a', 'Adck4', 'Ltbp4', 'Gja1', 'Klf4', 'Mapre1', 'Uap1', 'Fmo1', 'Errfi1', 'Lsm7', 'Tppp3', 'Akap2', 'Calu', 'Csrnp1', 'Fkbp8', 'Jund', 'Lpgat1', 'Crtap', 'Ptrf', 'Cebpb', 'Pxn', 'Ces2g', 'Pros1', 'Clic4', 'Adamts1', 'Add3', 'Fkbp1a', 'Eif5a', 'Ifit1', 'Lpar1'], '7': ['Clec3b', 'Cd34', 'Lum', 'Rbp1', 'Pi16', 'Col15a1', 'Htra3', 'Mfap5', 'Scara5', 'Gsn', 'Pdlim2', 'Sparcl1', 'Prss23', 'Fbln1', 'Entpd2', 'Krtdap', 'Igfbp6', 'Npy1r', 'Dmkn', 'Gpc3', 'Ifitm1', 'Itm2a', 'Anxa2', 'Nsg1', 'Sdpr', 'Sfrp1', 'Osr1', 'Nbl1', 'Chodl', 'Igfbp4', 'S100a10', 'Pmp22', 'Prkcdbp', 'Asz1', 'Col14a1', 'Ccdc80', 'Inmt', 'Plxdc2', 'Abi3bp', 'Nrp2', 'Pbx1', 'Bmp4', 'Tuba1a', 'Lbp', 'Cfb', 'Col8a1', 'Ifi27l2a', 'Timp3', 'Cfh', 'Dcn', 'Cygb', 'Slc25a4', 'Tcf21', 'Bglap2', 'Cntn4', 'Pcolce', 'Laptm4a', 'Cyr61', 'Meis2', 'Ifitm2'], '6': ['Tm4sf1', 'Cav2', 'Ddit4', 'Bcam', 'Gng11', 'Tinagl1', 'Utrn', 'Tns1', 'Crip2', 'Ripk1', '0610010K14Rik', 'Klf13', 'Calm1', 'Mzt2', 'Ly6c1', 'Tsc22d1', 'Isca1', 'Psme1', 'Mex3c', 'Limd1', 'Isg15', 'Nrp1', 'Sbno2', 'Serinc3', 'Kcne4', 'Ccnd2', 'Ednrb', 'Arrdc3', 'Vps26b', 'Prnp', 'Tmem66', 'Pdcd10', 'Ankrd11', 'Mapk3', 'Btg1', 'Ssbp3', 'Rab5c', 'Dhx15', 'Sap18', 'Glul', 'Fam162a', 'Trpm7', 'Ilk', 'Rap1a', 'Zbtb7a', 'Pttg1', 'Dennd5b', 'Ap2m1', 'Pten', 'Fosl2', 'Ccdc85b', 'Tpst2', 'Grb2', 'Eif4h', 'Eif3e', 'Jkamp', 'Rere', 'Cdc26', 'Psmd3', 'Mcl1']}
    header1=['UID', 'AAAGTAGAGATGCCAG', 'AACCGCGTCCAACCAA', 'AACTCCCGTCGGGTCT', 'AACTCTTAGTTGCAGG', 'AACTCTTTCATAACCG', 'AAGACCTAGATCCGAG', 'AAGACCTAGGTGGGTT', 'AAGGAGCGTGCAACTT', 'AAGTCTGAGATAGTCA', 'ACAGCTATCATACGGT', 'ACATACGCAGCTCCGA', 'ACCAGTACACGTCTCT', 'ACGAGCCAGATAGCAT', 'ACGGAGATCGGAATCT', 'ACGGGCTGTCTGATCA', 'ACTGAACCATCCGGGT', 'ACTTTCAAGATGGGTC', 'AGACGTTTCTTACCGC', 'AGAGCTTTCTGTTGAG', 'AGCGGTCTCTGAGTGT', 'AGCTTGAAGACTAAGT', 'AGGGAGTCAATCGGTT', 'AGGGATGAGTCGCCGT', 'AGGGTGACAGTCGTGC', 'AGGTCCGAGGACACCA', 'AGTCTTTAGCTAGTGG', 'AGTGAGGAGTGTACTC', 'AGTGAGGTCAGCATGT', 'AGTGTCACAGTATCTG', 'ATCACGACAACAACCT', 'ATCGAGTTCTCCAGGG', 'ATTATCCAGTAAGTAC', 'ATTTCTGCATGTTGAC', 'ATTTCTGGTCAGGACA', 'CAACCAAAGACGCACA', 'CACACCTAGCCCAATT', 'CACACCTGTTTGACAC', 'CACATAGCAGATGGCA', 'CACATTTTCTGGTATG', 'CACCAGGCAAGGCTCC', 'CACCTTGAGGACCACA', 'CACTCCACACATCCAA', 'CAGAATCTCTCAAGTG', 'CAGCATAAGATGGGTC', 'CAGCTGGTCCTTTACA', 'CATTCGCTCTCGTATT', 'CCAGCGATCGCTGATA', 'CCCAGTTCACATGACT', 'CCGGTAGAGGCAGTCA', 'CCGGTAGCAGTGGGAT', 'CCGTGGAGTGAGGGTT', 'CCGTTCAGTTCCCTTG', 'CCTACACAGGGCATGT', 'CCTTCGATCTCCGGTT', 'CGAGAAGAGGACTGGT', 'CGATTGATCTAACTTC', 'CGGACTGTCTTGTATC', 'CGGAGTCAGTACGCCC', 'CGGGTCATCGTGGACC', 'CGTAGGCCACGGACAA', 'CGTCCATAGTCCTCCT', 'CGTGAGCGTGCACTTA', 'CTAGAGTGTCATATCG', 'CTAGCCTAGCTGCGAA', 'CTCATTAGTGAAAGAG', 'CTCGAGGAGGCGCTCT', 'CTCGAGGAGGGTTTCT', 'CTCGTACTCAATACCG', 'CTCTAATAGTGAATTG', 'CTCTAATCACGGCGTT', 'CTCTGGTAGTGTACCT', 'CTGAAACGTAGTACCT', 'CTGCCTACAGCCACCA', 'CTGCGGAGTATGCTTG', 'CTGCTGTTCGTTACGA', 'CTGTGCTCAATGGACG', 'CTTGGCTCACGAAACG', 'GAAATGACAATCCAAC', 'GAACATCAGTGTTAGA', 'GAACCTACATTAGGCT', 'GAAGCAGCACCGATAT', 'GAATAAGGTCGGCATC', 'GACGGCTAGTCGTTTG', 'GACGTTAGTAGGCTGA', 'GACGTTATCAGGTTCA', 'GACTAACCATCGACGC', 'GAGGTGATCCAGTATG', 'GATCGCGTCACCCTCA', 'GATTCAGCAGTGGAGT', 'GCAAACTAGAGCCCAA', 'GCATGATCATATGGTC', 'GCGACCAAGGGATCTG', 'GCGACCACACCAACCG', 'GCGAGAACACATGGGA', 'GCGAGAATCCCTCTTT', 'GCTGCTTTCAAAGTAG', 'GGAATAATCATTCACT', 'GGACAGACATTCACTT', 'GGCGTGTCATTGGGCC', 'GGCTCGAAGCGTGAGT', 'GGGCATCTCCGTTGTC', 'GGTGCGTAGAGGTACC', 'GTAACGTTCAGGCGAA', 'GTAACTGCAGCATACT', 'GTACGTACACACCGAC', 'GTACTCCTCGCCTGTT', 'GTAGTCAAGACGACGT', 'GTATCTTGTAGAGGAA', 'GTATCTTTCGCCATAA', 'GTCACAATCACTTCAT', 'GTCACAATCTACTCAT', 'GTCATTTGTCTTGCGG', 'GTCGGGTAGGCCGAAT', 'GTCTCGTCAAACTGTC', 'GTGCAGCGTGGTACAG', 'GTGCAGCTCAGTTAGC', 'GTGGGTCAGATGTGGC', 'GTGTGCGTCAGTTCGA', 'GTTAAGCCACCAGGTC', 'GTTCTCGGTGCTTCTC', 'GTTTCTAAGGACTGGT', 'TAAGAGAAGCTATGCT', 'TAAGCGTAGAAACCGC', 'TACGGATAGGTACTCT', 'TACTTGTGTGTTCTTT', 'TAGCCGGAGGACTGGT', 'TAGGCATTCTATCCTA', 'TCAACGACAATGGTCT', 'TCCACACCAGCTGCTG', 'TCTTCGGGTTTAGGAA', 'TCTTTCCTCCGCATAA', 'TGAGAGGAGACATAAC', 'TGAGAGGAGCGATATA', 'TGAGCATAGTGGAGTC', 'TGCGGGTTCTTTACGT', 'TGGACGCTCCGATATG', 'TGGACGCTCTGCAGTA', 'TGGCTGGGTCGCTTCT', 'TGGGAAGAGTAGGCCA', 'TTCTACAAGGCAGTCA', 'TTGAACGTCTCTTATG', 'TTGCCGTCACGCCAGT', 'TTGGAACGTTAGGGTG', 'TTGGCAACATGATCCA', 'TTGGCAAGTCCGTTAA', 'TTGTAGGCATCCGGGT', 'TTTACTGCACACGCTG', 'TTTACTGCAGGACGTA', 'TTTGTCAGTTGCGTTA']
    outputDir='/Users/saljh8/Desktop/DemoData/ICGS-Mm//NMF-SVM/SVMOutputs'
    root_dir='/Users/saljh8/Desktop/DemoData/ICGS-Mm/'
    species='Mm'
    uniqueIDs={'Cygb': '', 'Strn3': '', 'Nsa2': '', 'Atp13a3': '', 'Nampt': '', 'Ehf': '', 'Tspan3': '', 'Bckdk': '', 'Crk': '', 'Tspan8': '', 'Syt8': '', 'Pkd1': '', 'Zfhx3': '', 'Gpx1': '', 'Gpx3': '', 'Gpx4': '', 'Gpx7': '', 'Gpx8': '', 'Bcl3': '', 'Avpr1a': '', 'Syf2': '', 'Cirbp': '', 'Rnf128': '', 'Cstf3': '', 'Uqcrq': '', 'Hspd1': '', 'Akt2': '', 'Scara5': '', 'Uqcrh': '', 'Myadm': '', 'Uqcrb': '', 'Wdr89': '', 'Hsd11b1': '', 'Cebpb': '', 'Cebpa': '', 'Ezr': '', 'Taf1d': '', 'Nid1': '', 'Hsdl2': '', 'Sap18': '', 'Apoc1': '', 'Ccar1': '', 'Pitpna': '', 'Rab27b': '', 'Rbm8a': '', 'Tmed5': '', 'Rtf1': '', 'Npdc1': '', 'Btg2': '', 'Btg1': '', 'Scpep1': '', 'Ndrg1': '', 'Ndrg2': '', 'Wasl': '', 'Ppap2b': '', 'Ppp1r2': '', 'H13': '', 'Krtcap2': '', 'Cpz': '', 'Pole4': '', 'Cdv3': '', 'Lama2': '', 'Lama4': '', 'Lama5': '', 'Gng10': '', 'Katna1': '', 'Gng12': '', 'Hint1': '', 'Slc25a5': '', 'Amfr': '', 'Col16a1': '', 'Pgk1': '', 'Cxcl1': '', 'Stub1': '', 'D4Wsu53e': '', 'Golim4': '', 'St7': '', 'Ccrn4l': '', 'Arhgef25': '', 'Reep3': '', 'Rap2b': '', 'Spint1': '', 'Mrpl14': '', 'Psmd8': '', 'Mrpl17': '', 'Ier5': '', 'Psmd7': '', 'Junb': '', 'Psmd1': '', 'B2m': '', 'Psmd3': '', 'Ier3': '', 'Ncl': '', 'Ifrd1': '', 'Lmna': '', 'Dbi': '', 'Phb2': '', 'Sorcs2': '', 'Ermp1': '', 'Cmas': '', 'Gstp2': '', 'Gstp1': '', 'Mcl1': '', 'Hes1': '', 'Tmem208': '', 'Pttg1': '', 'Lamtor2': '', 'Atp5f1': '', 'Serpinb6a': '', 'Wbp5': '', 'Mcfd2': '', 'Imp3': '', 'Acat1': '', 'Tmem176a': '', 'Tmem176b': '', 'Sars': '', 'Tppp3': '', 'Cd52': '', 'Soat1': '', 'Cd55': '', 'Ahnak': '', '0610010K14Rik': '', 'Emp2': '', 'Emp3': '', 'Emp1': '', 'Ifi27l2a': '', 'Ggh': '', 'Mmp14': '', 'Smap1': '', 'Dbnl': '', '2310036O22Rik': '', 'Fam114a1': '', 'Tpr': '', 'Dbp': '', 'Aldoa': '', 'Rlim': '', 'Cox6c': '', 'Lrrfip2': '', 'Akt1': '', 'Cope': '', 'Mgll': '', 'Nsmce1': '', 'Igtp': '', 'Gcc1': '', 'Ifitm3': '', 'Ano1': '', 'Nr1d1': '', 'Nr1d2': '', 'Irf1': '', 'Irf6': '', 'Nme2': '', 'Cstb': '', 'Irf9': '', 'Irf8': '', 'Ppargc1a': '', 'Sfpq': '', 'Clk1': '', 'Clk4': '', 'Tmem159': '', 'Herpud1': '', 'Pgrmc1': '', 'Mettl9': '', 'Tpp1': '', 'Tspan4': '', 'Rfk': '', 'Blvrb': '', 'Fuca1': '', 'Ing2': '', 'Ddit4': '', 'F830016B08Rik': '', 'Il11ra1': '', 'Zfp36l1': '', 'Gtf2h5': '', 'Cdc42se1': '', 'Cst3': '', 'Ncor1': '', 'Usmg5': '', 'Wdr1': '', 'Scarf2': '', 'Marcks': '', 'Uba1': '', 'Pcolce': '', 'Flna': '', 'Flnb': '', 'Lypd3': '', 'Creg1': '', 'Sep15': '', 'Ndufb8': '', 'Clint1': '', 'Ptgfrn': '', 'Tgm2': '', 'Zfp706': '', 'Sltm': '', 'Gstm1': '', 'Gstm2': '', 'Gstm4': '', 'Abhd12': '', 'Ilk': '', 'Cct8': '', 'Cct2': '', 'Cct7': '', 'Cct4': '', 'Cct5': '', 'Prom2': '', 'Lamp1': '', 'Arpc2': '', 'Fech': '', 'Dcun1d5': '', 'Tmem66': '', 'Tmem64': '', 'Raly': '', 'Rala': '', 'Hdac1': '', 'Hdac2': '', 'Litaf': '', 'Pros1': '', 'Notch1': '', 'Sfxn3': '', 'Vamp3': '', 'Clns1a': '', 'Fus': '', 'Rarg': '', 'Sphk2': '', 'Vamp8': '', 'Fbn1': '', 'Zcchc24': '', 'Rnaset2b': '', 'Set': '', 'Nbl1': '', 'Psme2': '', 'Ndufb11': '', 'Psme1': '', 'Ikzf2': '', 'H2afv': '', 'Myl12b': '', 'Myl12a': '', 'H2afy': '', 'H2afz': '', '2200002D01Rik': '', 'Fam84a': '', 'Tnfrsf12a': '', 'Cox7c': '', 'Atp5g3': '', 'Atp5g2': '', 'Cnn3': '', 'Cnn2': '', 'Hadhb': '', 'Spry2': '', 'Akr1b8': '', 'Tmem30a': '', 'Skp1a': '', 'Klf13': '', 'Cyb5b': '', 'Bmpr1a': '', 'Fam107b': '', 'Wls': '', 'Snrpd2': '', 'Dnlz': '', 'Snrpd1': '', 'Cd44': '', 'Cd47': '', 'Rap1b': '', 'Rap1a': '', 'Serp1': '', 'Phc2': '', 'Npy1r': '', 'Pcmtd1': '', 'Ncoa4': '', 'Dnajb4': '', 'Sumo1': '', 'Klf10': '', 'Sumo3': '', 'Sumo2': '', 'Dkk2': '', 'Dkk3': '', 'Ier3ip1': '', 'Chmp3': '', 'Rnf7': '', 'Rnf4': '', 'Net1': '', 'Podn': '', 'H2-D1': '', 'Ythdc1': '', 'Cd2ap': '', 'F2r': '', 'Col27a1': '', 'Rnh1': '', 'Odc1': '', 'Hjurp': '', 'Tubb6': '', 'Tubb5': '', 'Mtpn': '', 'Tmem167': '', 'Cdc26': '', 'Pdap1': '', 'Ext2': '', 'Dgat2': '', 'Gsk3a': '', 'Gsk3b': '', 'Ucp2': '', 'Ghitm': '', 'Isca1': '', 'Sdc1': '', 'Sdc2': '', 'Sdc4': '', 'Itgav': '', 'Nav1': '', 'BC031181': '', 'Gng5': '', 'Cnot4': '', 'Gtf2e2': '', 'Srrm2': '', 'Phf12': '', 'Phf13': '', 'Kdelc2': '', 'Bgn': '', 'Tmem120a': '', 'Csrp1': '', 'Limd1': '', 'Csrp2': '', 'Hnrnpul2': '', 'Ube2j1': '', 'Ralgds': '', 'Nop56': '', 'Nsg1': '', 'Fkbp7': '', 'Cdk4': '', 'Sarnp': '', 'Fkbp2': '', 'Fkbp3': '', 'Mt1': '', 'Trip12': '', 'Pea15a': '', 'Fkbp8': '', 'Fkbp9': '', 'Nkd1': '', 'Wasf2': '', 'Pbrm1': '', 'Sqstm1': '', 'Snrpb2': '', 'Ltbp3': '', 'Gltp': '', 'Ap2m1': '', 'Hbegf': '', 'Nipbl': '', 'Snrnp27': '', 'Pid1': '', 'C2': '', 'Ppp1r14b': '', 'M6pr': '', 'Edf1': '', 'Pde12': '', 'Pappa': '', 'Pdgfrb': '', 'Pdgfra': '', 'Ftl1': '', 'Prkar1a': '', 'Cox7a2': '', 'Ndufab1': '', 'Actn4': '', 'Capn1': '', 'Sat1': '', 'Ptp4a2': '', 'Gng11': '', 'Pltp': '', 'Npm1': '', 'Kcne4': '', 'Galnt2': '', 'Dstn': '', 'Xrn2': '', 'Cnih': '', 'Cd9': '', 'Srsf10': '', 'Med28': '', 'Rassf1': '', 'Prpf4b': '', 'Pmepa1': '', 'B3gat3': '', 'Dnm1': '', 'Dnm2': '', 'Tmem222': '', 'Ednrb': '', 'Brwd1': '', 'Ict1': '', 'Sphk1': '', 'Col11a1': '', 'Hmgb2': '', 'Hmgb1': '', 'Sparcl1': '', 'Ramp2': '', 'Tiparp': '', 'Ndufs7': '', 'Leprot': '', 'Bdnf': '', 'Uchl1': '', 'Gabarapl1': '', 'Mtch1': '', 'Gabarapl2': '', '3230401D17Rik': '', 'Eif1ax': '', 'Tmbim6': '', 'Sprr2g': '', 'Tmbim1': '', 'Sprr2b': '', 'UID': '', 'Cd34': '', 'Ech1': '', 'Tmem55b': '', 'Cited2': '', 'Plin2': '', 'Arl3': '', 'Gaa': '', 'Leprel2': '', 'Nfic': '', 'Nfib': '', 'Nfia': '', 'Ehd4': '', 'Cyr61': '', 'Ehd2': '', 'Ehd1': '', 'Trf': '', 'Scand1': '', 'Atxn7': '', 'Abl2': '', 'Atxn1': '', 'Slc6a6': '', 'Plscr1': '', 'Aurkaip1': '', 'Nme1': '', 'Atp5k': '', 'Atp5j': '', 'Csnk1d': '', 'Atp5h': '', 'Atp5o': '', 'Sra1': '', 'Atp5l': '', 'Atp5b': '', 'Atp5e': '', 'Atp5d': '', 'Bcap31': '', 'Plac8': '', 'Txlna': '', 'Zfand5': '', 'Zfand3': '', 'Ubxn4': '', 'Ubxn1': '', 'Timm8b': '', 'Cpeb2': '', 'Ly6c1': '', 'Spint2': '', 'Ift20': '', 'Emilin1': '', 'Lman2': '', 'Rdx': '', 'Banf1': '', 'Bambi': '', 'Rap2a': '', 'Has3': '', 'Elf1': '', 'Psmd9': '', 'Elf3': '', 'Yipf5': '', 'Yipf4': '', 'Gsn': '', 'Hprt': '', 'Spag9': '', 'Hk2': '', 'Wnt2': '', 'Fosb': '', 'Wnt4': '', 'Psmd4': '', 'Sfr1': '', 'Magoh': '', 'F11r': '', 'Hist1h1c': '', 'Sod2': '', 'Sod1': '', 'Gbp2': '', 'Gbp7': '', 'Rock2': '', 'Jund': '', 'Ddost': '', 'Ier2': '', 'Gpc6': '', 'Slc1a5': '', 'Lpar1': '', 'Dad1': '', 'Cox7a2l': '', 'Tgfb1': '', 'Samd4b': '', 'Socs3': '', 'Tmprss2': '', 'Arid5a': '', 'Cbr2': '', 'Col4a5': '', 'Vcam1': '', 'Col4a1': '', 'Arl4d': '', 'Uqcrc2': '', 'Arl4c': '', 'Nr1h2': '', 'D17Wsu104e': '', 'Sdpr': '', 'Hnrpdl': '', 'Tm4sf1': '', 'H6pd': '', 'Tacstd2': '', 'Uqcrfs1': '', 'Zmym2': '', 'Prickle1': '', 'Zmym5': '', 'Tubb4b': '', 'Vcan': '', 'Top1': '', 'Isg15': '', 'Tomm20': '', 'Gpc3': '', 'Tgfbi': '', 'Igf2r': '', 'Dpm3': '', 'Slc50a1': '', 'Acvrl1': '', 'Napa': '', 'Twsg1': '', 'Metap2': '', 'Id2': '', 'Id3': '', 'Cotl1': '', 'Id1': '', 'Dtnbp1': '', 'Galnt1': '', 'Nrp2': '', 'Procr': '', 'Slc9a3r1': '', 'Nr2f2': '', 'Shisa5': '', 'Lims1': '', 'Luzp1': '', 'Camk1': '', 'Ufc1': '', 'Psmg4': '', 'Krtdap': '', 'Gsta4': '', 'Gsta3': '', 'Ppp1r15a': '', 'Ppp1r15b': '', 'Pxn': '', 'Hist2h2aa1': '', 'Tmem205': '', 'Tsn': '', 'Fdps': '', 'Pabpn1': '', 'Copz2': '', 'Cyp1b1': '', 'Rasl11b': '', 'Rasl11a': '', 'Cmpk1': '', 'Ctbp1': '', 'Pon2': '', 'Iscu': '', 'Fnta': '', 'Kcnk1': '', 'Serpine1': '', 'Dap': '', 'Serpine2': '', 'Sec61g': '', 'Sept7': '', 'Rbm39': '', 'Sept2': '', 'Nap1l1': '', 'Sgce': '', 'Tecr': '', 'Col1a2': '', 'Col1a1': '', 'Thy1': '', 'Ptprs': '', '0610009D07Rik': '', 'Dpep1': '', 'Ptprf': '', 'Eci1': '', 'Ipo5': '', 'Ptprk': '', 'Prdm1': '', 'Dpysl2': '', 'Arl6ip5': '', 'Arl6ip1': '', 'Pten': '', 'Grn': '', 'Chpf': '', 'Zfp266': '', 'H2-T23': '', 'Prmt1': '', 'Dhx40': '', 'H2-Eb1': '', '1110008P14Rik': '', 'Calm2': '', 'Calm3': '', 'Calm1': '', 'AW112010': '', 'Aga': '', 'Clec11a': '', 'Basp1': '', 'Abca1': '', 'Sprr1a': '', 'Wtap': '', 'Pdcl3': '', 'Cox5b': '', 'Btbd1': '', 'Qsox1': '', 'Vapb': '', 'Vapa': '', 'Tmem109': '', 'Frmd6': '', 'Myd88': '', 'Ostc': '', 'Yipf3': '', 'Pmp22': '', 'Sub1': '', 'H3f3b': '', 'H3f3a': '', 'Serbp1': '', 'Eif4e': '', 'Tns1': '', 'Tagln': '', 'Eif4b': '', 'Eif4g2': '', 'Eif4g1': '', 'Asz1': '', 'Bola3': '', 'Bola2': '', 'Ubn2': '', 'Ubn1': '', 'Gnai2': '', 'Naa50': '', 'Bhlhe40': '', 'Glul': '', 'Akap12': '', 'Akap13': '', 'Acaa1b': '', 'B4galt1': '', 'Pgls': '', '2010107E04Rik': '', 'Mif': '', 'Icam1': '', 'Anxa11': '', 'Aard': '', 'Efhd2': '', 'Mlf2': '', 'Ntan1': '', 'Tpm4': '', 'Tpm3': '', 'Tpm2': '', 'Tpm1': '', 'Sfn': '', 'Efemp1': '', 'Dmkn': '', 'Mmp2': '', 'Rbm7': '', 'Pnp2': '', 'Krt23': '', 'Tmem33': '', 'Rarres2': '', 'Fam46a': '', 'Klhdc2': '', 'Klhdc3': '', 'Ddah2': '', 'Xbp1': '', 'Dnaja2': '', 'Capza2': '', 'Pgd': '', 'Dnaja1': '', 'Tpst2': '', 'Mrpl42': '', 'Serinc1': '', 'Serinc3': '', 'Sf1': '', 'Wdr26': '', 'Psma6': '', 'Sec11a': '', 'Apoe': '', 'Tomm7': '', 'Dnmt3l': '', 'Sdcbp': '', 'Arpp19': '', 'Ilkap': '', 'Mrpl52': '', 'Eif1b': '', 'Mlec': '', 'Grhl3': '', 'Ssr2': '', 'Zfp36': '', 'Mgp': '', 'Srp14': '', 'Cytip': '', 'Cd248': '', 'Papss1': '', 'Scamp5': '', 'Antxr1': '', 'Gstt1': '', 'Polr2e': '', 'Polr2f': '', 'Shh': '', 'Polr2l': '', 'Polr2i': '', 'Polr2k': '', 'Pbx1': '', 'Ptgs2': '', 'Ptgs1': '', 'Dock9': '', 'Cd24a': '', 'Cfh': '', 'Cfb': '', 'Cdk11b': '', 'Serpina3g': '', 'Rbm25': '', 'Golgb1': '', 'Skil': '', 'Etnk1': '', 'Perp': '', 'Impdh2': '', 'Lum': '', 'D16Ertd472e': '', 'Hsd17b12': '', 'Gspt1': '', 'Hsd17b11': '', 'Vasn': '', 'Stmn2': '', 'Ube2d3': '', 'Fkbp1a': '', 'Baiap2': '', 'Vasp': '', 'Slc25a11': '', 'Nudc': '', 'Sh3bgrl': '', 'Per1': '', 'Utrn': '', 'Fblim1': '', 'Ptges': '', 'Grem2': '', 'Fkbp14': '', 'Lamc1': '', 'Fkbp10': '', 'Fkbp11': '', 'App': '', 'Pvrl2': '', 'Wsb1': '', 'Ptov1': '', 'Minos1': '', 'Grcc10': '', 'Ifitm2': '', 'Ifitm1': '', 'Prdx4': '', 'Topors': '', 'Scd2': '', 'Scd1': '', 'Hoxa10': '', 'Purb': '', 'Ldb1': '', 'Cycs': '', 'Prdx2': '', 'Ephx1': '', 'Slbp': '', 'Arpc3': '', 'Pdia6': '', 'Arpc1a': '', 'Pdia4': '', 'Pdia3': '', 'Arpc5': '', 'Arpc4': '', 'Tmem119': '', 'Itpr1': '', 'Dusp11': '', 'Col4a2': '', 'Myh10': '', 'Sgpl1': '', 'Arsb': '', 'Sft2d1': '', 'Meis2': '', 'Srp9': '', 'Nufip2': '', 'Gps1': '', 'U2af1': '', 'Txnl1': '', 'U2af2': '', 'Paqr5': '', 'Lgals9': '', 'Paqr6': '', 'Tnfsf9': '', 'Bud31': '', 'Tmed10': '', 'Cul1': '', 'Calr': '', 'Vcl': '', 'Calu': '', 'Fstl1': '', 'Vcp': '', 'Capn5': '', 'Ptp4a1': '', 'Capn2': '', 'Tm9sf3': '', 'Tm9sf2': '', 'Dab2': '', 'Rrbp1': '', 'Zfp131': '', 'Praf2': '', 'Cadm3': '', 'H2-K1': '', 'Rp9': '', 'Phb': '', 'Mustn1': '', 'Mfap5': '', 'Mfap4': '', 'Mfap2': '', 'Cdh1': '', 'Prelid1': '', 'Plxdc2': '', 'Plk3': '', 'Plk2': '', 'Cdc42ep5': '', 'Tatdn2': '', 'Col6a1': '', 'Col6a3': '', 'Col6a2': '', 'Timm17b': '', 'Timm17a': '', 'Myh9': '', 'Hspa1a': '', 'Glrx3': '', 'Rras': '', 'Ap3s1': '', 'Snrnp70': '', 'Eif5b': '', 'Eif5a': '', 'Pdlim5': '', 'Tmem29': '', 'Peli1': '', 'Pdlim1': '', 'Pdlim2': '', 'Ppfibp2': '', 'Txn2': '', 'Prdm6': '', 'Txn1': '', 'Polr1d': '', 'Gpaa1': '', 'Prdm2': '', 'Ski': '', 'Tgif1': '', 'Rab5a': '', 'Pnrc1': '', 'Rab5c': '', 'Paics': '', 'Ddr2': '', 'Thsd4': '', 'Lrpap1': '', 'Creb3l1': '', 'Plec': '', 'Rgs2': '', 'Serf2': '', 'Psma2': '', 'Psma3': '', 'Psma1': '', 'Ablim1': '', 'Psma7': '', 'Psma4': '', 'Psma5': '', 'Mrpl48': '', 'Bcl7c': '', 'Cdk2ap1': '', 'Cdk2ap2': '', 'Apoa1bp': '', 'Nolc1': '', 'Oaz2': '', 'Gltscr2': '', 'Oaz1': '', 'Krt15': '', 'Pdk4': '', 'Dsg2': '', 'Krt19': '', 'Krt18': '', 'Atp5c1': '', 'Ctdnep1': '', 'Cd164': '', 'Meis1': '', 'Fam110c': '', 'Fosl1': '', 'Crem': '', 'Fosl2': '', 'Lgals3': '', 'Tmsb4x': '', 'Lgals1': '', 'Dcn': '', 'Pigyl': '', 'Ptgr1': '', 'Nrd1': '', 'Tnfrsf21': '', 'Golga7': '', 'Plbd2': '', 'Atxn7l3b': '', 'Csnk2b': '', '1190003J15Rik': '', 'Nedd8': '', 'Ppp2r1a': '', 'Srpr': '', 'Ckb': '', 'Tcp11l2': '', 'Nedd4': '', 'Sec62': '', 'Sec63': '', 'Rbp1': '', 'Stxbp2': '', 'Rbp4': '', 'Gfpt1': '', 'Gfpt2': '', 'Hdgf': '', 'Ankrd12': '', 'Spcs1': '', 'Ywhaz': '', 'Cisd1': '', 'Spcs2': '', 'Ywhaq': '', 'Celf2': '', 'Ube2e3': '', 'Ywhah': '', 'Colec12': '', 'Slc25a25': '', 'Ywhab': '', 'Hras1': '', 'Ywhae': '', 'Rpp21': '', 'Yif1b': '', 'Ece1': '', 'Pi16': '', 'Etfa': '', 'Gja1': '', 'Etfb': '', 'Las1l': '', 'Golph3': '', 'Xdh': '', 'Hibadh': '', 'Vdac2': '', 'Vdac1': '', 'Ankrd11': '', 'Rrp1': '', 'Aes': '', 'Rcn3': '', 'Emg1': '', 'Rcn1': '', 'Esyt1': '', 'Phip': '', 'Prrx2': '', 'Smarca2': '', 'Prkd3': '', 'Tmem123': '', 'Tmem158': '', 'Sidt2': '', 'Etf1': '', 'Oaf': '', 'Trex1': '', 'Hn1': '', 'Oat': '', 'Igfbp6': '', 'Igfbp7': '', 'Igfbp4': '', 'Ppig': '', 'Phf5a': '', 'Rnasek': '', 'Ppic': '', 'Ppib': '', 'Ppia': '', 'Lsmd1': '', 'Igf1r': '', 'Khdrbs1': '', 'Nrp1': '', 'D15Ertd621e': '', 'Grb2': '', 'Hmg20b': '', 'Ubb': '', 'Ubc': '', 'Thoc7': '', 'Plin3': '', 'Ubl5': '', 'Zfp935': '', 'Rnase4': '', 'Ubl3': '', 'Wnk1': '', 'C1d': '', 'Nop10': '', 'C1s': '', 'Il1b': '', 'Tmem63a': '', 'Tax1bp3': '', 'Ociad1': '', 'Chd4': '', 'Aida': '', 'Erf': '', 'Rer1': '', 'Cox6a1': '', 'Cdc42': '', 'Acsl3': '', 'Psap': '', 'Hnrnpu': '', 'Tmed1': '', 'Erp44': '', 'Hnrnpk': '', 'Hnrnpl': '', 'Hnrnpm': '', 'Fyttd1': '', 'Hnrnpc': '', 'Hnrnpd': '', 'Pde5a': '', 'Dph3': '', 'Spop': '', 'Tbrg1': '', 'Stk25': '', 'Gtf2a2': '', 'D17Wsu92e': '', 'Marcksl1': '', 'Thrap3': '', 'AI462493': '', 'Pum2': '', 'Pum1': '', 'Mdh2': '', 'Mdh1': '', 'Tnks2': '', 'Pbxip1': '', 'Fzd1': '', 'Dnajc1': '', 'Snx18': '', '2900097C17Rik': '', 'Glo1': '', 'Myc': '', 'Cxcl12': '', 'Cxcl14': '', 'Cxcl16': '', 'Osr1': '', 'Cnpy3': '', 'Cnpy2': '', 'Psmb8': '', 'Psmb7': '', 'Psmb6': '', 'Psmb5': '', 'Psmb4': '', 'Psmb3': '', 'Psmb2': '', 'Psmb1': '', 'Adh1': '', 'Smg1': '', 'Tnfrsf1a': '', '4932438A13Rik': '', 'Srpx': '', 'Trib1': '', 'Pgam1': '', 'Srp72': '', 'Pabpc1': '', 'Mrps24': '', 'Mrps25': '', 'Pabpc4': '', 'Slc7a11': '', 'Kcnk2': '', 'Mrps21': '', 'Llph': '', 'Tnrc6b': '', 'Zfp423': '', 'Plp2': '', 'Gtf2i': '', 'Luc7l3': '', 'Luc7l2': '', 'Fer1l4': '', 'Rab4b': '', 'Sepw1': '', 'Egr2': '', 'Egr1': '', 'Me1': '', 'Wfdc1': '', 'Smoc2': '', 'Wfdc2': '', 'Hsbp1': '', 'Ndufa6': '', 'Ndufa7': '', 'Ndufa4': '', 'Ndufa5': '', 'Ndufa2': '', 'Ndufa3': '', 'Ndufa1': '', 'Atp6ap1': '', 'Sec13': '', 'Ndufa8': '', 'Coro1b': '', 'Coro1a': '', 'Pycard': '', 'Srxn1': '', 'Adk': '', 'Samm50': '', 'Huwe1': '', 'Pink1': '', 'Bag3': '', 'Bag1': '', 'Bag6': '', 'H2-DMa': '', 'Chmp4b': '', 'Tspan13': '', 'Vdac3': '', '9530068E07Rik': '', 'Tnc': '', 'Serinc2': '', 'Tbc1d16': '', 'Ccnl1': '', 'Ddx17': '', 'Ccnl2': '', 'Tmem9': '', 'Ndufa11': '', 'Ndufa12': '', 'Ndufa13': '', 'Thbd': '', 'Tmem5': '', 'Srsf5': '', 'Atp6v0b': '', 'Srsf7': '', 'Srsf6': '', 'Srsf3': '', 'Srsf2': '', 'Bglap2': '', 'Esd': '', 'Abcd3': '', 'Rcn2': '', 'Rheb': '', 'Pomp': '', 'Slc38a2': '', 'Nfkbia': '', 'Nfkbiz': '', 'Ndufs2': '', 'Ndufs3': '', 'Ndufs4': '', 'Ndufs6': '', 'Yif1a': '', 'Ndufs8': '', 'Psenen': '', 'Degs1': '', 'Degs2': '', 'Sugt1': '', 'Hcfc1r1': '', 'Fmnl2': '', 'Bmyc': '', 'Dynlrb1': '', 'Ahsa1': '', 'Atp1b3': '', 'Atg3': '', 'Arfgap3': '', 'Jak1': '', 'Swi5': '', 'Zfp36l2': '', 'Gsdmc2': '', 'Copz1': '', 'Hpcal1': '', 'Csde1': '', 'Lima1': '', 'Rbbp7': '', 'Rbbp6': '', 'Rbbp4': '', 'Ccni': '', 'Tob1': '', 'Tob2': '', 'Hnrnpa3': '', 'Hnrnpa0': '', 'Siva1': '', 'Klf6': '', 'Klf5': '', 'Klf4': '', 'Klf2': '', 'Mbd1': '', 'Cox6b2': '', 'Islr': '', 'Klf9': '', 'Snx3': '', 'Snx1': '', 'Zc3h11a': '', 'Asah1': '', 'Arap1': '', 'Pex5l': '', 'Gnb2l1': '', 'Clta': '', 'Rnf40': '', 'Picalm': '', 'Romo1': '', 'Notum': '', 'Higd2a': '', 'Rabac1': '', 'Msrb2': '', 'Evi5': '', 'Hnrnpab': '', 'Ddb1': '', 'Chka': '', 'Impact': '', 'Tsen34': '', 'Pim3': '', 'Pim1': '', 'Dst': '', 'Zfp703': '', 'Dsp': '', 'Ap2b1': '', 'Dnajb9': '', 'Dnajb1': '', 'Dnajb6': '', 'Cdr2l': '', 'Ptplb': '', 'Sec24d': '', 'Sec24a': '', '1810055G02Rik': '', 'Capg': '', 'Fgfr1': '', 'Dynll1': '', 'Dynll2': '', 'Prrc2c': '', 'Prrc2b': '', 'Prrc2a': '', 'Ccdc80': '', 'Pdcd6ip': '', 'Tgfb1i1': '', 'Txnl4a': '', 'Pde4b': '', 'Psmc3': '', 'Psmc4': '', 'Psmc5': '', 'Ints6': '', 'Numa1': '', 'Hsp90b1': '', 'Srsf11': '', 'Atp6v0a1': '', 'Atp2b1': '', 'Diap1': '', 'Nkain4': '', 'Canx': '', 'Lrp10': '', 'Ube2h': '', 'Dpy30': '', 'Atp5a1': '', 'Casp4': '', 'Mrps33': '', 'Vti1b': '', 'Oit1': '', 'Sepp1': '', 'Actg1': '', 'Actg2': '', 'Hif1a': '', 'Usp9x': '', 'Ifi27l1': '', 'Tead1': '', 'Mtdh': '', 'Crispld2': '', 'Slx1b': '', 'Rnf10': '', 'Rnf11': '', 'Azin1': '', '3632451O06Rik': '', 'C1qbp': '', 'Prss23': '', 'Crip1': '', 'Ewsr1': '', 'Fam134a': '', 'Crip2': '', 'Mdk': '', '2810417H13Rik': '', 'Ahcyl2': '', 'Sec61b': '', 'Fxyd1': '', 'Adam23': '', 'Ddrgk1': '', 'Ilf2': '', 'Ifit3': '', 'Mapk8ip1': '', 'Ifit1': '', 'Slc35f5': '', 'Errfi1': '', 'Ak3': '', 'Dclk1': '', 'Ak1': '', 'Adprh': '', 'Vim': '', 'Mzt2': '', 'Tmem38a': '', 'Vat1': '', 'Csnk1a1': '', 'Ddx24': '', 'Tsc22d1': '', 'Tsc22d3': '', 'Tsc22d4': '', 'Hilpda': '', 'Ssr1': '', 'Sparc': '', 'Ssr3': '', 'Ssr4': '', 'Mxra7': '', 'Tacc1': '', 'Mxra8': '', 'Bsg': '', 'Tacc2': '', 'Cox16': '', 'Cox17': '', 'Fam102b': '', 'Lpp': '', 'Slc12a6': '', 'Kdm6b': '', 'Copb1': '', 'Gls': '', 'Hhip': '', 'Birc6': '', 'Osgin1': '', 'Alyref': '', 'Eef1b2': '', 'Aprt': '', 'Dctpp1': '', 'Gypc': '', 'Rere': '', 'Nfat5': '', 'Rerg': '', 'Lgmn': '', 'Ccndbp1': '', 'Higd1a': '', 'Vps37b': '', 'Acaa2': '', 'Tuba4a': '', 'Cdc37': '', 'Trmt112': '', '1600029D21Rik': '', 'Anxa2': '', 'Anxa5': '', 'Sh3bgrl3': '', 'Anxa7': '', 'Anxa6': '', 'Lsp1': '', 'Nenf': '', 'Ovol1': '', 'Atf5': '', 'Atf4': '', 'Lrrc58': '', 'Lrrc59': '', 'Smarca4': '', 'Atf3': '', 'Tm2d1': '', 'Tm2d2': '', 'Nipal1': '', 'Nr4a2': '', 'Nr4a1': '', 'Bax': '', 'Ergic3': '', 'Bad': '', 'Pcbp2': '', 'Pcbp1': '', 'Gdi2': '', 'Cmtm3': '', 'Erh': '', 'Cmtm7': '', 'Lbp': '', 'Col8a1': '', 'Ran': '', 'H2-Ab1': '', 'Gpnmb': '', 'Il33': '', 'Eny2': '', 'Stt3b': '', 'Stt3a': '', 'Plagl1': '', 'Has1': '', 'Sbds': '', 'Mesdc2': '', 'Fn1': '', 'Tgfbr2': '', 'Rtn4': '', 'Asph': '', 'Rtn3': '', 'Tpi1': '', 'Timm13': '', 'Psca': '', 'Pa2g4': '', '2700094K13Rik': '', 'Aqp3': '', 'Ap2s1': '', 'Ddx3x': '', 'Syngr2': '', 'Prr13': '', 'Aebp1': '', 'Mbd3': '', 'Tcf12': '', 'Syne1': '', 'Fndc3b': '', 'Ybx1': '', '1700025G04Rik': '', 'Pla2g16': '', 'Kazald1': '', 'Taldo1': '', 'Nono': '', 'Clic4': '', 'H1f0': '', 'Hexim1': '', 'Fxyd6': '', 'Sirt2': '', 'Fxyd4': '', 'Fxyd5': '', 'Fxyd3': '', 'Hspb1': '', 'Phf21a': '', 'Hspb8': '', 'Clptm1l': '', 'Wdr61': '', 'Myof': '', 'Abhd16a': '', 'Hectd1': '', 'Upk3a': '', 'Upk2': '', 'Acox1': '', 'Ergic2': '', 'Paip2': '', 'Pdpn': '', 'Ergic1': '', 'Derl1': '', 'Phpt1': '', 'Dazap1': '', 'Dazap2': '', 'Lrp1': '', 'Sertad1': '', 'Atp2a2': '', 'Clic1': '', 'Acin1': '', 'Gabarap': '', 'Cd200': '', 'Aplp2': '', 'Ociad2': '', 'Sertad2': '', 'Add3': '', 'Acp1': '', 'Add1': '', 'Itm2c': '', 'Itm2b': '', 'Itm2a': '', 'Taf10': '', 'Pik3r1': '', 'Sema4a': '', 'Tssc4': '', 'Cops7a': '', 'Matn2': '', 'Klhl21': '', 'Tmco1': '', 'Rab6a': '', 'Srebf1': '', 'Ube2v1': '', 'Numbl': '', 'Ndufc1': '', 'Ndufc2': '', 'Fmo1': '', 'Fmo5': '', 'Hmox1': '', 'Shfm1': '', 'Mef2a': '', 'Krt7': '', 'Krt5': '', 'Krt4': '', 'Krt8': '', 'Actr3': '', 'Actr2': '', 'Fam168b': '', 'Ddx3y': '', 'Adamtsl5': '', 'Gas7': '', 'Adam10': '', 'Cd151': '', 'Lhfp': '', 'Tax1bp1': '', 'Nucks1': '', 'Srgn': '', 'Prelp': '', 'Pttg1ip': '', 'Ddx39b': '', 'Ptgis': '', 'Chmp2b': '', 'Chmp2a': '', 'Cish': '', 'Fmod': '', 'Dhx15': '', 'Hebp1': '', '1700094D03Rik': '', 'Ash1l': '', 'Nhp2l1': '', 'Aff4': '', 'Pnpla2': '', 'Abcf1': '', 'Nsun2': '', 'Adrm1': '', 'Iqgap1': '', 'Tceb2': '', 'Tceb1': '', 'Tln1': '', 'Tspan31': '', 'Cox6b1': '', 'Matr3': '', 'Ccdc85b': '', 'Ube2l3': '', 'Msi2': '', 'Mbnl1': '', 'Mbnl2': '', 'Sox4': '', 'Pcna': '', 'Cttn': '', 'Glis2': '', 'C1ra': '', 'Hipk3': '', 'Pcnp': '', 'Zbtb7a': '', 'Fos': '', 'Lsm4': '', 'Eef1g': '', 'Eef1d': '', 'Lsm7': '', 'Mylk': '', 'Fam103a1': '', 'Guk1': '', 'Cbx6': '', 'Ndufb10': '', 'Cbx3': '', 'Iigp1': '', 'Inmt': '', 'Parl': '', 'Hspa12b': '', 'Zmat2': '', 'Hnrnpf': '', 'Slc43a3': '', 'Ccl2': '', 'Eif1a': '', 'Cpxm1': '', 'H2-Q4': '', 'Tmf1': '', 'Myl6': '', 'Tm7sf3': '', 'Myl9': '', 'Calcrl': '', 'Gorasp2': '', 'Mgst2': '', 'Trappc6b': '', 'Arhgdib': '', 'St13': '', 'Stk17b': '', 'Vps29': '', 'Adamts1': '', 'Adamts2': '', 'Mt2': '', 'Kif5b': '', 'Ifi35': '', 'Wac': '', 'Mdm2': '', 'Pigp': '', 'Laptm4a': '', 'Hspa5': '', 'Cyb5': '', 'Rcan1': '', 'Hspa9': '', 'Hspa8': '', 'Churc1': '', 'Sec61a1': '', 'Dnajb11': '', '2410015M20Rik': '', 'Parva': '', 'Ctdsp2': '', 'Emd': '', 'Mia3': '', 'Cox5a': '', 'Naca': '', 'Cuta': '', 'Car3': '', 'Rab18': '', 'Cyba': '', 'Tpd52l2': '', 'Got2': '', 'Rab10': '', 'Ptpn1': '', 'Got1': '', 'Rab14': '', 'Hint2': '', 'Htra3': '', 'Cuedc2': '', 'Vma21': '', 'Il6': '', 'Naalad2': '', 'Snf8': '', 'Ctnnd1': '', 'Drap1': '', 'Acta2': '', '2810403A07Rik': '', 'Copa': '', 'Rnaseh2c': '', 'Il4ra': '', 'Atp5g1': '', 'Bclaf1': '', 'Pcf11': '', 'Ddx5': '', 'Ddx6': '', 'Cd302': '', 'Anp32b': '', 'H2afj': '', 'Anp32a': '', 'Rbms1': '', 'Rnf19b': '', 'Mfap1a': '', 'Pebp1': '', 'Adipor2': '', 'Lrrfip1': '', 'Ndufb9': '', 'Adipor1': '', 'Ndufb7': '', 'Ndufb6': '', 'Ndufb5': '', 'Ndufb4': '', 'Ndufb3': '', 'Ndufb2': '', 'Serping1': '', 'Fbxo30': '', 'Aamp': '', 'Sik1': '', 'Sbsn': '', 'Syncrip': '', 'Ces2g': '', 'Vkorc1': '', 'Nrarp': '', 'Son': '', 'Tapbp': '', 'Mrps16': '', 'Mrps14': '', 'Bst2': '', 'Manf': '', 'Mylip': '', 'Zfp207': '', 'Xpnpep1': '', 'Ctsh': '', 'Ctsk': '', 'Ctsl': '', 'Arid5b': '', 'Sgk1': '', 'Ctsa': '', 'Ctsb': '', 'Ctsd': '', 'Ctse': '', 'Ctsz': '', 'Rexo2': '', 'Ctss': '', 'Gns': '', 'Ganab': '', 'Col14a1': '', 'Mapk3': '', 'Mapk1': '', 'Prkcdbp': '', 'Eprs': '', 'Cldn4': '', 'Ckap4': '', 'Cldn7': '', 'Ptges3': '', 'Rad23a': '', 'Sepn1': '', 'Lman1': '', 'Cltc': '', 'Cltb': '', 'Fbln2': '', 'Fbln1': '', 'Cnrip1': '', 'Maged2': '', 'Fbln5': '', 'Socs7': '', 'Socs1': '', 'Ptma': '', 'Gnas': '', 'Dhrs7': '', 'Idh1': '', 'Pkn2': '', 'Dhrs3': '', 'Laptm5': '', 'Ptms': '', 'Atp1a2': '', 'Atp1a1': '', 'Snrpb': '', 'Snrpc': '', 'Dnajc3': '', 'Snrpf': '', 'Snrpg': '', 'Snrpe': '', 'Os9': '', 'Uqcr10': '', 'Uqcr11': '', 'Suds3': '', 'Gna13': '', 'Gna11': '', 'S100a11': '', 'S100a10': '', 'S100a13': '', 'S100a14': '', 'S100a16': '', 'Slc3a2': '', 'Cdkn1a': '', 'Cxcl10': '', 'Txnip': '', 'Ugp2': '', 'Wfdc15b': '', 'Slc14a1': '', 'Ppp1r10': '', 'Cox4i1': '', 'Alkbh5': '', 'Nucb1': '', 'Nucb2': '', 'Tomm6': '', 'Sbno2': '', 'Entpd5': '', 'Tomm5': '', 'Entpd2': '', 'Ccm2': '', 'Msn': '', 'Mll5': '', 'Map1lc3a': '', 'Nol7': '', 'Map1lc3b': '', 'Phlda1': '', 'Phlda3': '', 'Efnb2': '', 'Bri3': '', 'Prdx6': '', 'Prdx5': '', 'Nupr1': '', 'Prdx3': '', 'Comt': '', 'Prdx1': '', 'Tgoln1': '', 'Pitx2': '', 'Pitx1': '', 'Eif4h': '', 'Loxl1': '', 'Loxl2': '', 'Ebpl': '', 'Hspa1b': '', 'Puf60': '', 'P4ha1': '', 'Eif4ebp1': '', 'Ifngr1': '', 'Abi3bp': '', 'Htra1': '', 'Cyc1': '', 'Atp6ap2': '', 'Use1': '', 'Sypl': '', 'Hsp90aa1': '', 'Egln2': '', 'Eif3i': '', 'Abhd14a': '', 'Upk1a': '', 'Upk1b': '', 'Lrrc8a': '', 'Pofut2': '', 'Zfp750': '', 'Sqrdl': '', 'Max': '', 'Fam120a': '', 'Taf9': '', 'Ostf1': '', 'Map2k2': '', 'Mal': '', 'Dtx3': '', 'Sf3b1': '', 'Park7': '', 'Hdlbp': '', 'Eif5': '', 'Glud1': '', 'Cast': '', 'Ssu72': '', 'Eif3g': '', 'Stx5a': '', 'Rbm42': '', 'Cyb5r3': '', 'Larp7': '', 'Eif6': '', 'Anapc16': '', 'Anapc11': '', 'Anapc13': '', 'Tmsb10': '', 'Mast4': '', 'Ppdpf': '', 'Morf4l2': '', 'Morf4l1': '', 'Fbl': '', 'Chic2': '', 'Selm': '', 'Rbm5': '', 'Rbm3': '', 'Cnih4': '', 'Cand1': '', 'Selk': '', 'Adck4': '', 'Olfml3': '', 'Tmcc3': '', 'Foxa1': '', 'Serpinf1': '', 'Ptrf': '', 'Atpif1': '', 'Kdm5b': '', 'S100a1': '', 'Rnf187': '', 'Csrnp1': '', 'Jkamp': '', 'Slc39a1': '', 'Ggnbp2': '', 'Slc39a7': '', 'Edem1': '', 'Scp2': '', 'Tinagl1': '', 'Gngt2': '', 'Surf1': '', 'Trpm7': '', 'Surf4': '', 'Cxadr': '', 'Cops8': '', 'Kdelr1': '', 'Atp5j2': '', 'Anapc7': '', 'Cops6': '', 'Cops4': '', 'Hnrnpa2b1': '', 'Fundc2': '', 'Dennd5b': '', 'Bcam': '', 'Hexa': '', 'Hexb': '', 'Ifnar2': '', 'C4b': '', 'Dcaf8': '', 'Palld': '', 'Ube2l6': '', 'Pfn1': '', 'Mycbp2': '', 'Serpinh1': '', 'Spty2d1': '', 'Stat6': '', 'Rab2a': '', 'Stat3': '', 'Stat2': '', 'Ldha': '', 'Fndc1': '', 'Rpia': '', 'Pmm1': '', 'Sept10': '', 'Ldhb': '', 'Cldn23': '', 'Dync1i2': '', 'Cldn25': '', 'Ahr': '', 'Mgst3': '', 'Mgst1': '', 'Palm': '', 'Rbms3': '', 'Rbms2': '', 'Jag1': '', 'Cald1': '', 'Ube4b': '', 'Bnip3l': '', 'Tprgl': '', '0610007P14Rik': '', 'Tex264': '', 'Srd5a1': '', 'Tubb2a': '', 'Bloc1s1': '', 'Mknk2': '', 'Pnp': '', 'Zyx': '', 'Pnn': '', 'Eif3h': '', 'Tnfsf12': '', 'Eif3k': '', 'Eif3l': '', 'Eif3m': '', 'Ypel3': '', 'Sf3b5': '', 'Eif3a': '', 'Eif3b': '', 'Eif3c': '', 'H2-Aa': '', 'Eif3e': '', 'Eif3f': '', 'Sf3b2': '', 'Atp8b1': '', 'Cfl2': '', 'Cfl1': '', 'S100a6': '', 'Timm23': '', 'F3': '', 'Timp3': '', 'Timp2': '', 'Timp1': '', 'Efna5': '', 'Senp2': '', 'Efna1': '', 'Rdh10': '', 'Sh3glb1': '', 'Tcf25': '', 'Tcf21': '', 'Ccdc124': '', 'Capns1': '', 'Pla2g7': '', 'Lpgat1': '', 'Uqcrc1': '', 'Gpi1': '', 'Selenbp1': '', 'Hspg2': '', 'Rbpms': '', 'Emid1': '', 'Usf2': '', 'Prnp': '', 'Crtap': '', 'Nptn': '', 'Irf2bp2': '', 'Gapdh': '', 'Lgals3bp': '', 'Rab11b': '', 'Ugcg': '', 'Rab34': '', 'Stx7': '', 'Pcyox1': '', 'Stra13': '', 'Tspo': '', 'Tnfaip8': '', 'Coq10b': '', 'Tnfaip3': '', 'Tnfaip2': '', 'Ctnnb1': '', 'Avpi1': '', 'Chmp1a': '', 'Tmod3': '', 'Tslp': '', 'Lamb1': '', 'Lamb2': '', 'Ssb': '', 'Gnl3': '', 'Atox1': '', 'Aldh3a1': '', 'Thbs2': '', 'Cat': '', 'Exosc7': '', 'Sema3c': '', 'Ccdc23': '', 'Enpp2': '', 'Clec3b': '', 'Arpc1b': '', 'Pja1': '', 'Ssbp3': '', 'Atp6v1f': '', 'Atp6v1a': '', 'Ralbp1': '', 'Rsrc2': '', 'Tmem219': '', 'Gpr153': '', 'Arl4a': '', 'Arf3': '', 'Arf1': '', 'Arf6': '', 'Arf5': '', 'Arf4': '', 'Fdx1': '', 'Mex3c': '', 'Tmed7': '', 'Rpn1': '', 'Rpn2': '', 'Tmed3': '', 'Tmed2': '', 'Ghr': '', 'Ppp4c': '', 'Dguok': '', 'Tmed9': '', 'Nqo1': '', 'Pdcd6': '', 'Pdcd4': '', 'Reep5': '', 'Ndufv3': '', 'Rab11a': '', 'Ppp1cb': '', 'Ppp1cc': '', 'Ppp1ca': '', 'Commd3': '', 'Commd6': '', 'Cnbp': '', 'Tkt': '', 'Osbpl9': '', 'Atp6v1g1': '', 'Rbx1': '', 'Txnrd1': '', 'Tsg101': '', 'Snhg11': '', 'Oxct1': '', 'Dync1h1': '', 'Ndufa4l2': '', 'Atrn': '', 'Bin1': '', 'Mmp23': '', 'Spon2': '', 'Spon1': '', 'Ndel1': '', 'Pdcd10': '', 'Sar1a': '', 'Sar1b': '', 'Lsm6': '', 'Fscn1': '', 'Bscl2': '', 'Lsm1': '', 'Igfbp2': '', 'Col12a1': '', 'Ranbp1': '', 'Tpst1': '', 'Tmem50a': '', 'Far1': '', 'Mapre1': '', 'Tcf4': '', 'Tcf3': '', 'Erp29': '', 'Nfil3': '', 'Zbtb20': '', 'Dusp6': '', 'Dusp1': '', 'Pam': '', 'Nadk': '', 'Rabggtb': '', 'Mfge8': '', 'Ugdh': '', 'Chodl': '', 'Metrnl': '', 'Uba52': '', 'Ms4a4d': '', 'Col15a1': '', 'Fam76a': '', 'Rragc': '', 'Ppa1': '', 'Mat2a': '', 'Eif4a1': '', 'Eif4a2': '', 'Hp1bp3': '', 'Col18a1': '', '1110007C09Rik': '', 'Por': '', 'Ripk1': '', 'Hnrnph1': '', 'Ripk4': '', 'Dlgap4': '', 'Lasp1': '', 'Cdh11': '', 'Ngfrap1': '', 'Wbp11': '', 'Actb': '', 'Anapc5': '', 'Ubr5': '', 'Ubr4': '', 'Rab25': '', 'C1qtnf2': '', 'C1qtnf1': '', 'Rhov': '', 'Rhoj': '', 'Tcp1': '', 'Fis1': '', 'Myeov2': '', 'Rhoc': '', 'Rhob': '', 'Rhoa': '', 'Alpl': '', 'Bzw1': '', 'Gnai3': '', 'Eef2': '', 'Eif2s2': '', 'Midn': '', 'Dlg1': '', 'Cyp2f2': '', 'Fam195b': '', 'Efemp2': '', 'Pamr1': '', 'Slc2a1': '', 'Tmem59': '', 'Cd63': '', 'Ctnna1': '', 'Foxq1': '', 'Slc25a3': '', 'Ube2r2': '', 'Slc25a4': '', 'Npc2': '', 'Mafk': '', 'Ost4': '', 'Fth1': '', 'Maff': '', 'Gas6': '', 'Areg': '', 'Actn1': '', 'Rac1': '', 'Angptl2': '', 'Ccdc12': '', '1810037I17Rik': '', 'Mrfap1': '', '1110004F10Rik': '', 'Mrc2': '', 'Mrpl33': '', 'Mrpl30': '', 'Vps26b': '', 'Fhl1': '', 'Anpep': '', 'Eid1': '', 'Gata6': '', 'Camk2n1': '', 'Gata3': '', 'Ly6e': '', 'Ly6d': '', 'Ly6a': '', 'Nfe2l2': '', 'Nfe2l1': '', 'Sri': '', 'Srm': '', 'Axl': '', 'Api5': '', 'Hmgn3': '', 'Hmgn2': '', 'Hmgn1': '', 'Akr1a1': '', 'Vmp1': '', 'Ccnd3': '', 'Ccnd2': '', 'Tnxb': '', '1810009A15Rik': '', 'Sdha': '', 'Sdhb': '', 'Acadl': '', 'Sdhd': '', 'Tra2b': '', 'Cd74': '', 'Smad7': '', 'Jun': '', 'Smad3': '', 'Tmprss13': '', 'Jup': '', 'Tceal8': '', 'Moxd1': '', 'Gem': '', 'Sfrp2': '', 'Sfrp1': '', 'Fam162a': '', 'Atp6v0e': '', 'Denr': '', 'Cetn3': '', 'Cetn2': '', 'Mprip': '', 'Cox8a': '', 'BC056474': '', 'Atp6v0d1': '', 'Nop58': '', 'Ang': '', 'Sytl2': '', 'Elovl1': '', 'Dctn2': '', 'Dctn3': '', 'Ube2i': '', 'Hadh': '', 'Ube2k': '', 'Triobp': '', 'Ube2n': '', 'Plaur': '', 'Ube2b': '', 'Gnb1': '', 'Gnb2': '', 'Slc38a10': '', 'Arrdc3': '', 'Ube2s': '', '2700060E02Rik': '', 'Tuba1c': '', 'Tuba1b': '', 'Tuba1a': '', 'Mid1ip1': '', 'Aimp1': '', 'Atxn10': '', 'Carhsp1': '', 'Trappc2l': '', '1500012F01Rik': '', 'Tgfbr1': '', 'Akap2': '', 'Akap9': '', 'Lamp2': '', 'Brd4': '', 'Thbs1': '', 'Brd2': '', 'Ugt2b34': '', 'Sdf4': '', 'Sdf2': '', 'Sys1': '', 'Strap': '', 'Ets2': '', 'Pafah1b1': '', 'Pafah1b2': '', 'Capzb': '', 'P4hb': '', 'Polr2m': '', 'Erdr1': '', 'Hspe1': '', 'Rsu1': '', 'Zfand2b': '', 'Phf23': '', 'Zfand2a': '', 'Itgb1': '', 'Itgb5': '', 'Arih1': '', 'Vps25': '', 'Arl8b': '', 'Wdr92': '', 'Vps28': '', 'Pigk': '', 'Fau': '', 'Hsp90ab1': '', 'Rnd3': '', 'Pigt': '', '1810011O10Rik': '', 'Caprin1': '', 'Col3a1': '', 'Gsto1': '', 'Mpzl2': '', 'Mpzl1': '', 'Cryab': '', 'Arhgdia': '', 'Acaa1a': '', 'Tram1': '', 'Anxa1': '', 'Epcam': '', 'Cav2': '', 'Cav1': '', 'Rad23b': '', 'Serpina3n': '', 'Btf3': '', 'Ppp2ca': '', 'Maged1': '', 'Cox7b': '', 'Mgat1': '', 'Mgat2': '', 'G3bp1': '', 'G3bp2': '', 'Tor1aip2': '', 'Tor1aip1': '', 'Tbcb': '', 'Tbca': '', 'Txndc9': '', 'Runx1': '', 'Txndc5': '', 'Bcl10': '', 'Tmem43': '', 'Foxp1': '', 'Scyl1': '', 'Col5a2': '', 'Prkca': '', 'Col5a1': '', 'Kdelr2': '', 'Kdelr3': '', 'Lmbrd1': '', 'Fcgrt': '', '2810428I15Rik': '', 'Rab1': '', 'Rab11fip1': '', 'Cd81': '', 'Cd82': '', 'Ltbp4': '', 'Cck': '', 'Trp53i13': '', 'Pfdn1': '', 'Mrpl27': '', 'Pfdn2': '', 'Pfdn5': '', 'Mrpl20': '', 'Mrpl23': '', 'Gadd45a': '', 'Gadd45b': '', 'Vsig2': '', 'Anxa4': '', 'Dynlt3': '', 'Gadd45g': '', 'Srrm1': '', 'Sec31a': '', 'Rwdd1': '', 'Erlec1': '', 'Epn1': '', 'Ecm1': '', 'Tmem214': '', '8430408G22Rik': '', 'Ptbp1': '', 'Plat': '', 'Chchd2': '', 'Chchd1': '', 'Lmo1': '', 'Arid1a': '', 'Lmo4': '', 'Epha2': '', 'Lmo7': '', 'Prkcsh': '', 'Tagln2': '', 'Tmem45a': '', 'Tbc1d10a': '', 'Ppm1a': '', 'Gclc': '', 'Aldh2': '', 'Dek': '', 'Nhp2': '', 'Aldh1a2': '', 'Irgm1': '', 'Samhd1': '', 'Trim29': '', 'Hs3st1': '', 'Ivl': '', 'Trim25': '', 'Txndc17': '', 'Atp6v1e1': '', 'Tef': '', 'Itpkc': '', 'Fabp5': '', 'Cntn4': '', 'Sart1': '', 'Ppp1r14a': '', 'Itih5': '', 'Ryk': '', 'Ttc3': '', 'D8Ertd738e': '', 'Uap1': '', 'Tmem14c': '', 'Pcbd2': '', 'Atl3': '', 'Rsl1d1': '', 'Hsd17b10': '', 'Ppp2r4': '', 'Stmn1': '', 'Bmp4': '', 'Bmp2': '', 'Bmp1': '', 'Myo1c': '', '1110008F13Rik': '', 'Tmem147': '', 'Gda': '', 'Tmem140': ''}
    

    generateMarkerheatmap(processedInputExpFile,output_file,NMFSVM_centroid_cluster_dir,groupsdict,markergrps,header1,outputDir,root_dir,species,uniqueIDs)
    sys.exit()
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
