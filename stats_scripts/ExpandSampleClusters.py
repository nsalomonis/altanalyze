#!/usr/bin/env python

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import numpy as np
import pylab as pl
import os.path
import scipy
from import_scripts import sampleIndexSelection
import matplotlib.pyplot as plt
import export

from sklearn import datasets, linear_model
from sklearn.preprocessing import StandardScaler

from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import pairwise_distances

from visualization_scripts import Orderedheatmap

#from sklearn import cross_validation

from sklearn import svm
from sklearn.multiclass import OneVsOneClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC
from sklearn import linear_model
import operator
from collections import OrderedDict
from collections import defaultdict
upd_guides=[]

#upd_guides.append("uid")
def FindTopUniqueEvents(Guidefile,psi,Guidedir):
    head=0
    guidekeys=[]
    exportnam=os.path.join(Guidedir,"SplicingeventCount1.txt")
    export_class=open(exportnam,"w")
    #commonkeys=[]
    tempkeys={}
    global upd_guides
    global train
    omitcluster=0
    
    unique_clusters={}

    for line in open(Guidefile,'rU').xreadlines():
        if head==0:
            head=1
            continue
        else:
            line1=line.rstrip('\r\n')
            q= string.split(line1,'\t')
            
            if abs(float(q[8]))>0.15:
                try:
                    tempkeys[q[2]].append([q[0],float(q[10]),q[11]])
                except KeyError:
                    tempkeys[q[2]]=[[q[0],float(q[10]),q[11]],]
    for i in tempkeys:
       
        if len(tempkeys[i])>1:
            #print tempkeys[i]
            tempkeys[i].sort(key=operator.itemgetter(1),reverse=False)
            #print tempkeys[i][0]
            
            try:
                unique_clusters[0].append(tempkeys[i][0])
            except KeyError:
                unique_clusters[0]=[tempkeys[i][0],]
          
        else:
            try:
                unique_clusters[0].append(tempkeys[i][0])
            except KeyError:
                unique_clusters[0]=[tempkeys[i][0],]
          
    unique_clusters[0].sort(key=operator.itemgetter(1))
  
    if len(unique_clusters[0])>100:
        guidekeys=unique_clusters[0]
        for i in range(0,len(guidekeys)):
            
            #upd_guides[i]=[upd_guides[i][3],upd_guides[i][4]]
            upd_guides.append(guidekeys[i][0])
    else:
        omitcluster=1
    export_class.write(psi+"\t"+str(len(unique_clusters[0]))+"\n")

    return omitcluster
    #return upd_guides,train

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filterRows(input_file,output_file,filterDB=None,logData=False,firstLine=True):
    orderlst={}
    counter=[]
    export_object = open(output_file,'w')
    row_count=0
    for line in open(input_file,'rU').xreadlines():
        #for i in filterDB:
            flag1=0
            row_count+=1
            data = cleanUpLine(line)
            line = data+'\n'
            values = string.split(data,'\t')
            if row_count == 1:
                initial_length = len(values)
            if len(values)!=initial_length:
                print values
            if firstLine:
                firstLine = False
                export_object.write(line)
            else:
                if values[0] in filterDB:
                    counter=[index for index, value in enumerate(filterDB) if value == values[0]]
                    for it in range(0,len(counter)):
                        orderlst[counter[it]]=line
    for i in range(0,len(orderlst)):
            try: export_object.write(orderlst[i])
            except Exception:
                print i,filterDB[i]
                continue
    export_object.close()
    print 'Filtered rows printed to:',output_file
    
def filterRows_data(input_file,output_file,filterDB=None,logData=False):
    filteredevents=[]
    tempevents=[]
    orderlst={}
    counter=[]
    export_object = open(output_file,'w')
    firstLine = True
    Flag=0;

    for i in filterDB:
        event=string.split(i,"|")[0]
        tempevents.append(event)
    for line in open(input_file,'rU').xreadlines():
        #for i in filterDB:
            flag1=0
            data = cleanUpLine(line)
            values = string.split(data,'\t')
            event=string.split(values[0],"|")[0]
                
            if firstLine:
                firstLine = False
                if Flag==0:
                    export_object.write(line)
            else:
                if event in tempevents:
                    counter=[index for index, value in enumerate(tempevents) if value == event]
                    #print counter
                    filteredevents.append(event)
                    for it in range(0,len(counter)):
                        orderlst[counter[it]]=line
                    if logData:
                        line = string.join([values[0]]+map(str,(map(lambda x: math.log(float(x)+1,2),values[1:]))),'\t')+'\n'
   
    for i in range(0,len(tempevents)):
        if i in orderlst:
            export_object.write(orderlst[i])
            if "\n" not in orderlst[i]:
                #print i
                export_object.write("\n") 
    export_object.close()
    tempevents2=[]
    #print 'Filtered rows printed to:',output_file
    for i in range(len(tempevents)):
        if tempevents[i] in filteredevents:
            tempevents2.append(tempevents[i])
    
   # print len(tempevents2) 
    return tempevents2

def findParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1
    return filename[:x]

def Classify(header,Xobs,output_file,grplst,name,turn,platform,output_dir,root_dir):
    count=0
    start=1
    Y=[]
    head=0
    for line in open(output_file,'rU').xreadlines():
        if head >count:
            val=[]
            counter2=0
            val2=[]
            me=0.0
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
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
            Y.append(val)
        
        else:
            head+=1
            continue

    Xobs=zip(*Xobs)

    Xobs=np.array(Xobs)
    Xobs=zip(*Xobs)

    Xobs=np.array(Xobs)
    X=grplst
    X=zip(*X)
    X=np.array(X)
    #print X
    Y=zip(*Y)
    Y=np.array(Y)

    #np.savetxt("/Volumes/MyPassport/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Leucegene/July-2017/PSI/ExpressionProfiles/DataPlots/complete_KNN.txt",q)
    #if platform=="PSI":
    #else:

    output_dir = output_dir+'/SVMOutputs'
    output_dir2 = root_dir+'/ICGS-NMF'
    if os.path.exists(output_dir)==False:
        export.createExportFolder(output_dir)
    if os.path.exists(output_dir2)==False:
        export.createExportFolder(output_dir2)
    #exportnam=output_dir+'/round'+str(turn)+'SVC_test_50cor.txt'
    #export_class=open(exportnam,"w")
    exportnam1=output_dir+'/round'+str(turn)+'SVC_decision_func.txt'
    export_class1=open(exportnam1,"w")
    
    if platform=="PSI":
        exportnam2=output_dir+'/round'+str(turn)+'SVC_Results.txt'
        export_class2=open(exportnam2,"w")
    else:
        exportnam2=output_dir2+'/FinalGroups.txt'
        export_class2=open(exportnam2,"w")
        
    exportnam3=output_dir+'/round'+str(turn)+'SVC_Results_max.txt'
    export_class3=open(exportnam3,"w")
    #export_class2.write("uid"+"\t"+"group"+"\t"+"class"+"\n")
    regr = LinearSVC()
    regr.fit(Xobs,X[:,0])
    q=regr.predict(Y)
    #print q
    count=1
    ordersamp={}
    order=[]
    for i in q:
        gr=string.split(name[int(i)-1],"_")[0]
        gr=gr.replace("V","")
        
        #export_class2.write(header[count]+"\t"+str(i)+"\t"+name[int(i)-1]+"\n")
       # export_class2.write(header[count]+"\t"+str(i)+"\t"+gr+"\n")
        ordersamp[header[count]]=[name[int(i)-1],str(i)]
        count+=1
    #print len(X[:,0])
    if len(X[:,0])>2:
        prob_=regr.fit(Xobs,X[:,0]).decision_function(Y)
        #k=list(prob_)

        export_class1.write("uid")
        #export_class2.write("uid")
        export_class3.write("uid")
        for ni in name:
            export_class1.write("\t"+"R"+str(turn)+"-"+ni)
            #export_class2.write("\t"+"R"+str(turn)+"-"+ni)
            export_class3.write("\t"+"R"+str(turn)+"-"+ni)
        export_class1.write("\n")
        #export_class2.write("\n")
        export_class3.write("\n")
        #print prob_
        for iq in range(0,len(header)-1):
            export_class1.write(header[iq+1])
            #export_class2.write(header[iq+1])
            export_class3.write(header[iq+1])
            for jq in range(0,len(name)):
                export_class1.write("\t"+str(prob_[iq][jq]))
                if prob_[iq][jq]==max(prob_[iq,:]):
                    #print ordersamp[header[iq+1]],name[jq]
                    if ordersamp[header[iq+1]][0]==name[jq]: 
                        order.append([header[iq+1],name[jq],prob_[iq][jq],ordersamp[header[iq+1]][1]])
                    export_class3.write("\t"+str(1))
                else:
                    export_class3.write("\t"+str(0))
                
            export_class1.write("\n")
            #export_class2.write("\n")
            export_class3.write("\n")
        export_class1.close()
        export_class3.close()
    else:
        if platfrm=="PSI":
            prob_=regr.fit(Xobs,X[:,0]).decision_function(Y)
            #k=list(prob_)
            export_class1.write("uid"+"\t")
            export_class2.write("uid"+"\t")
            export_class1.write("group")
            export_class2.write("round"+str(turn)+"-V1"+"\t"+"round"+str(turn)+"-V2"+"\n")
            #for ni in name:
            #   export_class1.write("\t"+ni)
            #   export_class2.write("\t"+ni)
            export_class1.write("\n")
            export_class2.write("\n")
            #print prob_
            #export_class1.write(header[1])
            #export_class2.write(header[1])
            for iq in range(0,len(header)-1):
                export_class1.write(header[iq+1])
                export_class2.write(header[iq+1])
                #for jq in range(0,len(X[:,0])):
                export_class1.write("\t"+str(prob_[iq]))
                if prob_[iq]>0.5:                    
                    export_class2.write("\t"+str(1)+"\t"+str(0))            
                else:
                    if prob_[iq]<-0.5:  
                        export_class2.write("\t"+str(0)+"\t"+str(1))
                    else:
                        export_class2.write("\t"+str(0)+"\t"+str(0))
                export_class1.write("\n")
                export_class2.write("\n")
    order = sorted(order, key = operator.itemgetter(2),reverse=True)
    order = sorted(order, key = operator.itemgetter(1))
    for i in range(len(order)):
        #export_class2.write(order[i][0]+"\t"+order[i][3]+"\t"+order[i][1]+"\n")
        gr=string.split(order[i][1],"_")[0]
        gr=gr.replace("V","")
        
        #export_class2.write(header[count]+"\t"+str(i)+"\t"+name[int(i)-1]+"\n")
        export_class2.write(order[i][0]+"\t"+order[i][3]+"\t"+gr+"\n")
    
    export_class2.close()
    if platform=="PSI":
        Orderedheatmap.Classify(exportnam2)
    else:
        Orderedheatmap.Classify(exportnam3)
       
def header_file(fname, delimiter=None):
    head=0
    header=[]
    new_head=[]
    with open(fname, 'rU') as fin:
        for line in fin:
            if head==0:
                line = line.rstrip(os.linesep)
                header=string.split(line,'\t')
                for i in header:
                    if ":" in i:
                        
                        i=string.split(i,":")[1]
                    new_head.append(i)     
                head=1
            else:break   
    return new_head

def avg(array):
    total = sum(map(float, array))
    average = total/len(array)
    return average

def Findgroups(NMF_annot,name):
    head=0
    groups=[]
    counter=1
    line=[]
    for exp1 in open(NMF_annot,"rU").xreadlines():
        lin=exp1.rstrip('\r\n')
        lin=string.split(lin,"\t")
        mapping={}
        if head==0:
            head=1
            continue
        else:
            line.append(lin)
            
    for j in range(0,len(name)):
        #print name[j]
            for i in range(len(line)):
                lin=line[i]
                key=lin[0]
                key1=key+"_vs"
                key2="vs_"+key+".txt"
                if key1 in name[j] or key2 in name[j]:
                    tot=0
                    for q in range(1,len(lin)):
                        
                        if lin[q]=='1':
                            tot=tot+1
                            groups.append(counter)
            
            counter=counter+1
    groups=np.asarray(groups)
    return groups

def TrainDataGeneration(output_file,NMF_annot,name,scaling=False,exports=False,rootDir='',calculatecentroid=True):
    head=0
    groups=[1,2]
    centroid_heatmap_input=''
    matrix=defaultdict(list)
    compared_groups={}
    for exp1 in open(NMF_annot,"rU").xreadlines():
        lin=exp1.rstrip('\r\n')
        lin=string.split(lin,"\t")
        mapping={}
        if head==0:
            header=lin
            head=1
            continue
        else:
            for i in range(1,len(lin)):
                if lin[i]=='1':
                    try:mapping[1].append(header[i])
                    except Exception: mapping[1]=[header[i]]
                else:
                    try:mapping[0].append(header[i])
                    except Exception: mapping[0]=[header[i]]
            head2=0
            #print len(mapping[1]),len(mapping[0])
            #print lin[0]
            eventname=[]
            for exp2 in open(output_file,"rU").xreadlines():
                    lin2=exp2.rstrip('\r\n')
                    lin2=string.split(lin2,"\t")
                    if head2==0:
                        group_db={}
                        index=0
                        try:
                            if len(mapping[1])>0 and len(mapping[0])>0:
                                for i in lin2[1:]:
                                    if i in mapping[1]:
                                        
                                        try: group_db[1].append(index)
                                        except Exception: group_db[1] = [index]
                                    else:
                                        try: group_db[2].append(index)
                                        except Exception: group_db[2] = [index]
                                    index+=1
                        except Exception:
                            break
                                
                        #print len(group_db[1])   
                        head2=1
                        continue
                    else:
                        key = lin2[0]
                        lin2=lin2[1:]
                        grouped_floats=[]
                        associated_groups=[]
            
                        ### string values
                        gvalues_list=[]
                        for i in group_db[1]:
                                try:
                                    x=float(lin2[i])
                                    gvalues_list.append(x)
                                    
                                except Exception:
                                    #try: gvalues_list.append('') ### Thus are missing values
                                    #except Exception: pass
                                    pass
                        if calculatecentroid:
                            try:
                            
                                matrix[lin[0]].append(avg(gvalues_list))
                                eventname.append(key)
                            except Exception:
                                matrix[lin[0]].append(float(0))
                                eventname.append(key)
                        else:
                        #matrix[lin[0]].append(gvalues_list)
                            try:
                                matrix[lin[0]].append(gvalues_list)
                                eventname.append(key)
                            except Exception:
                                matrix[lin[0]].append(float(0))
                                eventname.append(key)
    
    #export_class=open(exportnam,"w")
    #export_class.write('uid')
    #for i in range(len(eventname)):
    #export_class.write('\t'+eventname[i])
    #export_class.write('\n')

    keylist=[]
    matri=[]
    train=[]
    for j in range(0,len(name)):
        mediod=[]
        for key in matrix:
            """
            temp_ls=[]
            for ls in matrix[key]:
                temp_ls.append(np.median(ls))
            print key,sum(temp_ls)
            if sum(temp_ls)==0 and scaling==False:
                continue
            """
            if exports:
                key1=string.split(key,"-")[1]
                key2=string.split(key,"-")[1]
            else:
                key1=key+"_vs"
                key2="vs_"+key+".txt"
       
            if key1 in name[j] or key2 in name[j]:
                keylist.append(key)
                #print name[j]
                if calculatecentroid:
                    train.append(matrix[key])
                else:
                #mediod.append(matrix[key])
                    matri=zip(*matrix[key])
        #print len(matri)
        
      
        k=0
        if calculatecentroid==False:
            matri=np.array(matri)
            #print matri.shape
            n=matri.shape[0]
            D=pairwise_distances(matri,metric='euclidean').tolist()
            D=np.array(D)
       
            dist=np.mean(D,0)
            for i in np.argsort(dist):
                if k<1:
                   
                    train.append(np.array(matri[i]))
                    k=k+1
         
        #Xst=10000000
        #temp=[]
        #for i in range(n):
        #    Xd=np.array(matri[i])
        #    Xsd=Xd.reshape(1, -1)
        #   
        #    Xxd=pairwise_distances(Xsd,matri,metric='euclidean')
        #    #print Xxd
        #    dist=np.mean(Xxd)
        #    #print dist
        #    if dist<Xst:
        #        Xst=dist
        #        temp=Xd
        #train.append(temp)
        
    if exports:
        train1=zip(*train)
        centroid_heatmap_input=rootDir+'NMF-SVM/centroids/exp.MF.txt'
        centroid_heatmap_groups=rootDir+'NMF-SVM/centroids/groups.MF.txt'
        centroid_heatmap_obj = export.ExportFile(centroid_heatmap_input)
        centroid_heatmap_groups_obj = export.ExportFile(centroid_heatmap_groups)
        centroid_heatmap_obj.write("uid")
        train1=zip(*train)
        count=1

        for k in range(len(keylist)):
            centroid_heatmap_obj.write("\t"+keylist[k])
            centroid_heatmap_groups_obj.write(keylist[k]+"\t"+str(k)+"\t"+str(k)+"\n")
            count=count+1
        centroid_heatmap_obj.write("\n")
     
        train=zip(*train)
        for i in range(len(train)):
            ls=[]
            ls.append(eventname[i])
            for j in range(len(train[i])):
                ls.append(train[i][j])
            s=sum(ls[1:])
            #print eventname[i], s
            if s==0 and scaling == False:
                pass
            else:
                centroid_heatmap_obj.write(string.join(map(str,ls),'\t')+'\n')
        centroid_heatmap_obj.close()
        centroid_heatmap_groups_obj.close()
        #for i in range(len(matrix[key])):
        #export_class.write('\t'+str(matrix[key][i]))
        #export_class.write('\n')
        train=zip(*train)
    scaling=False
    train=np.array(train)
 
    if scaling==True:
        trainlst=[None]*len(train[0])      
        for i in range(len(train)):
            for j in range(len(train[i])):
                try:trainlst[j]=trainlst[j]+train[i][j]
                except Exception:trainlst[j]=train[i][j]
        trainlst=zip(*trainlst)
        train=np.array(trainlst)
    return train, centroid_heatmap_input

if __name__ == '__main__':
    import getopt
    group=[]
    grplst=[]
    name=[]
    matrix={}
    compared_groups={}
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['Guidedir=','PSIdir=','PSI=','NMF_annot='])
        for opt, arg in options:
            if opt == '--Guidedir': Guidedir=arg
            elif opt =='--PSIdir':PSIdir=arg
            elif opt =='--PSI':PSI=arg
            elif opt =='--NMF_annot':NMF_annot=arg
           
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    #commonkeys=[]
    counter=1
    #filename="/Users/meenakshi/Documents/leucegene/ICGS/Clustering-exp.Hs_RNASeq_top_alt_junctions367-Leucegene-75p_no149-Guide1 TRAK1&ENSG00000182606&I1.1_42075542-E2.1__E-hierarchical_cosine_correlation.txt"          
    #PSIfile="/Users/meenakshi/Documents/leucegene/ExpressionInput/exp.Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-367-Leucegene-75p-unique-filtered-filtered.txt"
    #keylabel="/Users/meenakshi/Documents/leucegene/ExpressionInput/exp.round2_glmfilteredKmeans_label.txt"
    for filename in os.listdir(Guidedir):
        if filename.startswith("PSI."):
            Guidefile=os.path.join(Guidedir, filename)
            psi=string.replace(filename,"PSI.","")
            PSIfile=os.path.join(PSIdir, psi)
            print Guidefile,PSIfile
  
            #output_file=PSIfile[:-4]+"-filtered.txt"
            #sampleIndexSelection.filterFile(PSIfile,output_file,header)
            omitcluster=FindTopUniqueEvents(Guidefile,psi,Guidedir)
            print omitcluster
            if omitcluster==0:
                group.append(counter)
                name.append(psi)
                counter+=1
        
    output_file=PSI[:-4]+"-filtered.txt"  

    print len(upd_guides)
    filterRows(PSI,output_file,filterDB=upd_guides,logData=False)
    header=header_file(output_file)
    
    train=TrainDataGeneration(output_file,NMF_annot,name)
    grplst.append(group)
    print grplst
    Classify(header,train,output_file,grplst,name)