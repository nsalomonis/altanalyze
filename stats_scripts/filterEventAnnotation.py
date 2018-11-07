import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import numpy as np

import os.path
from collections import defaultdict
import export

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
    export_object = export.ExportFile(output_file)
    firstLine = True
    Flag=0;
    print len(filterDB)
    #for i in filterDB:
    for line in open(input_file,'rU').xreadlines():
        #for i in filterDB:
            flag1=0
            
            data = cleanUpLine(line)
           
            values = string.split(data,'\t')
          
            
            if firstLine:
                firstLine = False
                k=values.index('UID')
                if Flag==0:
                    export_object.write(line)
            else:
                
                if values[k] in filterDB:
                    counter=[index for index, value in enumerate(filterDB) if value == values[k]]
                        #print counter
                    for it in range(0,len(counter)):
                        orderlst[counter[it]]=line
                   
                        #export_object.write(line)
                        #firstLine=True
                       # Flag=1;
               
                    
                #else:
                   # max_val = max(map(float,values[1:]))
                #min_val = min(map(float,values[1:]))
                #if max_val>0.1:

                     #   export_object.write(line)
    try:
        for i in range(0,len(orderlst)):
            export_object.write(orderlst[i])
    except Exception:
        print i,filterDB[i]
    
    
         
    export_object.close()
    print 'Filtered rows printed to:',output_file

def FilterFile(Guidefile,PSI,turn=0):
    if 'Clustering' in Guidefile:
        count=1
    
    else:
        count=0
    val=[]
    head=0
    for line in open(Guidefile,'rU').xreadlines():
        if head >count:
            
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            val.append(q[0])  
        else:
            head+=1
            continue
   
    dire = export.findParentDir(export.findParentDir(Guidefile)[:-1])
   
    output_dir = dire+'SubtypeAnalyses-Results'
    if os.path.exists(output_dir)==False:
        export.createExportFolder(output_dir)
    
    #output_file = output_dir+'/round'+str(turn)+'/'+export.findFilename(PSI)+'-filtered.txt'
    output_file = output_dir+'/round'+str(turn)+'/'+export.findFilename(PSI)[:-4]+'-filtered.txt'
    filterRows(PSI,output_file,filterDB=val)
    
    return output_file

if __name__ == '__main__':

    import getopt

    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['PSIfile=','PSIEvent='])
        for opt, arg in options:
            if opt == '--PSIfile': PSIfile=arg
            elif opt == '--PSIEvent':PSIEvent=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    #mutfile="/Users/meenakshi/Desktop/Leucegene-data1/Mutation_Annotations.txt"          
    #Guidefile="/Users/meenakshi/Documents/leucegene/ICGS/Round2_cor_0.6_280default/Clustering-exp.round2_insignificantU2like-Guide1 DDX5&ENSG00000108654&E3.4-E3.9__ENSG0000010-hierarchical_cosine_correlation.txt"          
    inputfile=FilterFile(PSIfile,PSIEvent)