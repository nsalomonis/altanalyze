### ADT.py
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

import math
import sys,string,os
import statistics
import export
import copy

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def replicateFACS(adt_cptt,adt_annotations,groups_file,steps):

    export_results = steps[:-4]+'-associations.txt'
    eos = export.ExportFile(export_results)
    
    import collections
    FACS_rules = collections.OrderedDict()
    for line in open(steps,'rU').xreadlines():
        data = cleanUpLine(line)
        protocol,ADT,direction,distribution, regional_global = string.split(data,'\t')
        FACS_rules[protocol,ADT]=direction,distribution,regional_global
          
    adt_names={} ### ADT name lookup table
    firstRow = True
    for line in open(adt_annotations,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            firstRow = False
            header = t[1:]
        else:
            adt_names[t[0]]=t[1]

    clusters={} ### ICGS, cellHarmony cluster or capture notations
    for line in open(groups_file,'rU').xreadlines():
        data = cleanUpLine(line)
        cell_barcode,group_number,cluster = string.split(data,'\t')
        clusters[cell_barcode]=cluster
           
    adt_values = {}
    firstRow = True
    for line in open(adt_cptt,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            firstRow = False
            cell_barcodes = t[1:]
        else:
            adt = t[0][:-5]
            adt = adt_names[adt]
            values = map(float,t[1:])
            adt_values[adt] = values
    
    excluded_barcodes={}
    alt1_excluded_barcodes={}
    alt2_excluded_barcodes={}
    occurance1 = True
    occurance2 = True
    for (protocol,ADT) in FACS_rules:
        direction,distribution,regional_global = FACS_rules[protocol,ADT]
        if ADT in adt_values:
            values = adt_values[ADT]
            if regional_global == 'regional':
                ### Compute statistics for only the filtered barcodes - similar to normal FACS
                index = 0
                values_regional=[]
                cell_barcodes_regional=[]
                for barcode in cell_barcodes:
                    if barcode not in excluded_barcodes_prime:
                        cell_barcodes_regional.append(barcode)
                        values_regional.append(values[index])
                    index+=1
                values = values_regional
            else:
                cell_barcodes_regional = cell_barcodes
            lower25,median_val,upper75th,int_qrt_range = statistics.iqr(list(values))
            avg = statistics.avg(list(values))
            stdev = statistics.stdev(list(values))
            #if ADT == 'CD3': print lower25,median_val,upper75th,int_qrt_range;sys.exit()
            if direction == '<':
                threshold = upper75th
            else:
                threshold = lower25
            if protocol == 'primary':
                excluded_barcodes = percentageMatching(ADT,threshold,values,cell_barcodes_regional,clusters,direction,excluded_barcodes)
                excluded_barcodes_prime = excluded_barcodes
            if protocol == 'alt1':
                if occurance1:
                    alt1_excluded_barcodes = copy.deepcopy(excluded_barcodes_prime)
                    occurance1 = False
                alt1_excluded_barcodes = percentageMatching(ADT,threshold,values,cell_barcodes_regional,clusters,direction,alt1_excluded_barcodes)
            if protocol == 'alt2':
                if occurance2:
                    alt2_excluded_barcodes = copy.deepcopy(excluded_barcodes_prime)
                    occurance2 = False
                alt2_excluded_barcodes = percentageMatching(ADT,threshold,values,cell_barcodes_regional,clusters,direction,alt2_excluded_barcodes)
    
    ### Find the filtered barcodes from the different strategies
    final_selected_barcodes={}
    filtered_clusters={}
    all_clusters={}
    for barcode in cell_barcodes:
        if barcode in clusters:
            cluster = clusters[barcode]
            try: all_clusters[cluster]+=1
            except: all_clusters[cluster]=1
            if barcode not in excluded_barcodes_prime:
                passed=False
                if len(alt1_excluded_barcodes)>0 and barcode not in alt1_excluded_barcodes:
                    final_selected_barcodes[barcode]=[]
                    passed=True
                if len(alt2_excluded_barcodes)>0 and barcode not in alt2_excluded_barcodes:
                    final_selected_barcodes[barcode]=[]
                    passed=True
                if passed:
                    try: filtered_clusters[cluster]+=1
                    except: filtered_clusters[cluster]=1
    for cluster in filtered_clusters:
        print cluster, direction,100*filtered_clusters[cluster]/(all_clusters[cluster]*1.000), filtered_clusters[cluster],all_clusters[cluster]
    for barcode in final_selected_barcodes:
        eos.write(barcode+'\n')
    eos.close()
    
def percentageMatching(ADT,threshold,values,cell_barcodes,clusters,direction,excluded_barcodes):
    index = 0
    filered_cell_barcodes={}
    filtered_clusters={}
    all_clusters={}
    #print ADT,threshold, max(values),direction
    for val in values:
        cell_barcode = cell_barcodes[index]
        if cell_barcode in clusters:
            cluster = clusters[cell_barcode]
            try: all_clusters[cluster]+=1
            except: all_clusters[cluster]=1
        if cell_barcode not in excluded_barcodes:
                if direction==">":
                    if val>threshold:
                        filered_cell_barcodes[cell_barcode]=[]
                        try: filtered_clusters[cluster]+=1
                        except: filtered_clusters[cluster]=1
                    else:
                        excluded_barcodes[cell_barcode]=[]
                elif direction=="<":
                    if val<threshold:
                        filered_cell_barcodes[cell_barcode]=[]
                        try: filtered_clusters[cluster]+=1
                        except: filtered_clusters[cluster]=1
                    else:
                        excluded_barcodes[cell_barcode]=[]
        index+=1
    #print adt,filtered_clusters,all_clusters
    for cluster in filtered_clusters:
        #pass
        print ADT,cluster, direction,100*filtered_clusters[cluster]/(all_clusters[cluster]*1.000), filtered_clusters[cluster],all_clusters[cluster],len(excluded_barcodes)
    return excluded_barcodes

if __name__ == '__main__':
    
    adt_cptt='/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/CPTT/exp.Biolegend-ADT-1.txt'
    adt_annotations='/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/feature_reference.txt'
    groups_file='/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/CPTT/AltAnalyze/cellHarmony/cellHarmony/QueryGroups.cellHarmony_captures-filtered.txt'
    steps='/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Isolation-Strategy/HSCP.txt'
    
    ################  Comand-line arguments ################
    import getopt

    #python stats_scripts/ADT.py --i /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/CPTT/exp.Biolegend-ADT-1.txt --a /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/feature_reference.txt --g /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/CPTT/AltAnalyze/cellHarmony/cellHarmony/QueryGroups.cellHarmony_captures-filtered.txt --s /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Isolation-Strategy/HSCP.txt
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print 'WARNING!!!! Too commands supplied.'
        
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','a=','g=','s=',])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i':
                adt_cptt = arg
            if opt == '--a':
                adt_annotations = arg
            if opt == '--g':
                groups_file = arg
            if opt == '--s':
                steps = arg
                
    replicateFACS(adt_cptt,adt_annotations,groups_file,steps)