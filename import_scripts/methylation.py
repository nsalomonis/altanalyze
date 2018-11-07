###ExonArray
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

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

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os.path
import unique
import time
import export

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    #add in code to prevent folder names from being included
    dir_list2 = []
    for file in dir_list:
        if '.txt' in file: dir_list2.append(file)
    return dir_list2

################# Begin Analysis

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def importAnnotations(filename):
    firstLine = True
    fn = filepath(filename)
    rows = 0
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line);
        tab_delimited_data = string.split(data,'\t')
        if rows > 10: sys.exit()
        print tab_delimited_data#;sys.exit()
        rows+=1

def correlateMethylationData(filename,betaLow=0.4,betaHigh=0.6,counts=-1):
    ### Takes a filtered pre-processed beta-value file as input
    firstLine = True
    rows=0; filtered=0
    for line in open(filename,'rU').xreadlines():             
        data = cleanUpLine(line);
        t = string.split(data,'\t')
        if firstLine:
            header = t
            if len(t)>5 and 'Illumina_name' in header:
                delimiter = -50
                annot_export_object.write(string.join([t[0]]+t[delimiter:],'\t')+'\n')
            else:
                delimiter = len(header)
            headers = t[1:delimiter]
            firstLine = False
            export_object.write(string.join([t[0]]+headers,'\t')+'\n')
        else:
            probeID = t[0]
            #try: beta_values = map(float,t[1:50])
            beta_values = map(lambda x: conFloat(x,t[1:delimiter]),t[1:delimiter])
            if '' in beta_values:
                print beta_values;sys.exit()
            high = sum(betaHighCount(x,betaHigh) for x in beta_values)
            low = sum(betaLowCount(x,betaLow) for x in beta_values)
    

def importMethylationData(filename,betaLow=0.4,betaHigh=0.6,counts=-1, filter=None):
    annot_file = filepath('AltDatabase/ucsc/Hs/Illumina_methylation_genes.txt')
    export_object = open(filename[:-4]+'-filtered.txt','w')
    
    print filename[:-4]+'-filtered.txt', counts
    firstLine = True
    rows=0; filtered=0
    for line in open(filename,'rU').xreadlines():             
        data = cleanUpLine(line);
        t = string.split(data,'\t')
        #export_object.write(string.join(t,'\t')+'\n')
        #"""
        if firstLine:
            header = t
            if len(t)>5 and 'Illumina_name' in header:
                delimiter = -50
                annot_export_object = open(annot_file,'w')
                annot_export_object.write(string.join([t[0]]+t[delimiter:],'\t')+'\n')
            else:
                delimiter = len(header)
            headers = t[1:delimiter]
            firstLine = False
            export_object.write(string.join([t[0]]+headers,'\t')+'\n')
        else:
            probeID = t[0]
            #try: beta_values = map(float,t[1:50])
            beta_values = map(lambda x: conFloat(x,t[1:delimiter]),t[1:delimiter])
            if '' in beta_values:
                print beta_values;sys.exit()
            high = sum(betaHighCount(x,betaHigh) for x in beta_values)
            low = sum(betaLowCount(x,betaLow) for x in beta_values)
            #if rows<50: print high, low, max(beta_values), min(beta_values)
            #else:sys.exit()
            #export_object.write(string.join(t[:delimiter])+'\n')
            if high>=counts and low>=counts:
            #if (high-low) > 0.2:
                #if rows<50: print 1
                if filter!=None:
                    if probeID in filter: proceed=True; probeID = str(filter[probeID])+':'+probeID
                    else: proceed = False
                else: proceed = True
                if proceed:
                    filtered+=1
                    export_object.write(string.join([probeID]+map(str,beta_values),'\t')+'\n')
            if 'Illumina_name' in header:
                annot_export_object.write(string.join([t[0]]+t[delimiter:],'\t')+'\n')
            rows+=1
        #"""
    export_object.close()
    if delimiter == '-50':
        annot_export_object.close()
    print filtered, rows

def conFloat(x,betaValues):
    try: x = float(x)
    except Exception:  x=None
    if x== None or x == 0:
        floats=[]
        for i in betaValues:
            if i=='': pass
            elif float(i)==0: pass
            else: floats.append(float(i))
        try: return min(floats)
        except Exception: print betaValues;sys.exit()
    else:
        return x

def betaHighCount(x,betaHigh):
    if x>betaHigh:
        return 1
    else: return 0
    
def betaLowCount(x,betaLow):
    if x<betaLow:
        return 1
    else: return 0

def getIDsFromFile(filename):
    filterIDs = {}
    fn = filepath(filename)
    for line in open(fn,'rU').xreadlines():      
        data = cleanUpLine(line);
        t = string.split(data,'\t')
        filterIDs[string.lower(t[0])]=[]
    return filterIDs

def getRegionType(filename,featureType=None,chromosome=None,filterIDs=None):
    if filterIDs !=None:
        filterIDs = getIDsFromFile(filterIDs)
        
    firstLine = True
    fn = filepath(filename)
    count=0; filter_db={}
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line);
        t = string.split(data,',')
        if firstLine:
            if len(t[2]) >0:
                header = t
                firstLine=False
                chr_ind = header.index('CHR')
                pos_ind = header.index('Coordinate_36')
                tss_ind = header.index('UCSC_RefGene_Group')
                gene_name = header.index('UCSC_RefGene_Name')
        else:
            probeID = t[0]
            count+=1
            try: gene_names = string.split(t[gene_name],';')
            except Exception: gene_names = []
            try:
                if chromosome != None:
                    if t[chr_ind] == chromosome:
                        if filterIDs !=None:
                            for gene in gene_names:
                                if string.lower(gene) in filterIDs:
                                    filter_db[probeID]=t[pos_ind]
                        else:
                            filter_db[probeID]=t[pos_ind]
                if 'promoter' in string.lower(featureType):
                    if 'TSS' in t[tss_ind]:
                        if filterIDs !=None:
                            for gene in gene_names:
                                if string.lower(gene) in filterIDs:
                                    filter_db[probeID]=t[pos_ind]
                        else:
                            filter_db[probeID]=t[pos_ind]
                if 'mir' in string.lower(featureType) or 'micro' in string.lower(featureType):
                    if 'mir' in string.lower(t[gene_name]) or 'let' in string.lower(t[gene_name]):
                        if filterIDs !=None:
                            for gene in gene_names:
                                if string.lower(gene) in filterIDs:
                                    filter_db[probeID]=t[pos_ind]
                        else:
                            filter_db[probeID]=t[pos_ind]
                if filterIDs !=None:
                            for gene in gene_names:
                                if string.lower(gene) in filterIDs:
                                    filter_db[probeID]=t[pos_ind]
            except Exception:
                pass
    
    print len(filter_db), 'probes remaining'
    return filter_db

if __name__ == '__main__':
    import getopt
    featureType = 'promoter'
    featureType = 'all'
    Species = 'Hs'
    filter_db=None
    chromosome=None
    numRegulated = -1
    analysis = 'filter'
    filterIDs = None
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a methylation beta-value file as input in the command-line"
        print "Example: python methylation.py --i /Users/me/sample1.txt --g /Users/me/human.gtf"
        sys.exit()
    else:
        analysisType = []
        useMultiProcessing=False
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','a=','t=','r=','c=','f='])
        for opt, arg in options:
            if opt == '--i': input_file=arg
            elif opt == '--a': analysis=arg
            elif opt == '--t': featureType=arg
            elif opt == '--r': numRegulated=int(arg)
            elif opt == '--c': chromosome=arg
            elif opt == '--f': filterIDs=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    
    if analysis == 'filter':
        filename = 'AltDatabase/ucsc/Hs/wgEncodeHaibMethyl450CpgIslandDetails.txt'
        #input_file = '/Volumes/SEQ-DATA/PCBC/Methylation/Methylome70allBValues_aronowAnnotations.txt'
        if featureType!= 'all' or chromosome != None or filterIDs!=None:
            filter_db = getRegionType(filename,featureType=featureType,chromosome=chromosome,filterIDs=filterIDs)
        importMethylationData(input_file,filter = filter_db,counts=numRegulated); sys.exit()
        #importAnnotations(methylation_file);sys.exit()
    if analysis == 'correlate':
        ### Performs all pairwise correlations between probes corresponding to a gene
        correlateMethylationData(input_file)
        