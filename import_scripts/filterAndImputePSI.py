#!/usr/bin/env python

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import export
import numpy as np
import os.path
from numpy import corrcoef, sum, log, arange
from scipy.stats.stats import pearsonr
import traceback
from stats_scripts import statistics
import unique

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def missingValueMedian(ls):
    ls2=[]
    for i in ls:
        try: ls2.append(float(i))
        except: pass
    return statistics.median(ls2)

def filterPSIValues(filename,impute=False,cutoff=75,returnValues=False):
    fn = filepath(filename)
    firstRow=True
          
    header = True
    rows=0
    filtered=0
    cutoff_name = str(int(cutoff))
    cutoff = cutoff/100.000
    new_file = filename[:-4]+'-%sp.txt' % cutoff_name
    new_file_clust = new_file[:-4]+'-clustID.txt'
    ea = export.ExportFile(new_file)
    eac = export.ExportFile(new_file_clust)
    added=[]
    PSI_db={}
    
    for line in open(fn,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            new_t=[]
            for x in t:
                if '.counts' in x:
                    x = string.split(x,'.junction_quantification.txt')[0]
                new_t.append(x)
            t = new_t
            header = False
            t = [t[8]]+t[11:]
            header_length = len(t)-1
            minimum_values_present = int(cutoff*int(header_length))
            not_detected = header_length-minimum_values_present
            new_header = string.join(t,'\t')+'\n'
            ea.write(new_header)
        else:
            cID = t[7]
            t = [t[8]]+t[11:]
            missing_values_at_the_end = (header_length+1)-len(t)
            missing = missing_values_at_the_end+t.count('')
            uid = t[0]
            if missing<not_detected:
                #if cID not in added:
                added.append(cID)
                if impute:
                    avg = missingValueMedian(t)
                    t = [str(avg) if x=='' else x for x in t+missing_values_at_the_end*['']]
                    if '' in t: print t;sys.exit()
                    t[0] = string.replace(t[0],':','__')
                    t[0] = string.replace(t[0],'|','&')
                    t[0] = string.split(t[0],'&')[0]
                new_line = string.join(t,'\t')+'\n'
                if returnValues:
                    PSI_db[uid]=new_line
                ea.write(new_line)
                eac.write(uid+'\t'+cID+'\n')
                filtered+=1
        rows+=1
    print rows, filtered
    ea.close()
    eac.close()
    #removeRedundantCluster(new_file,new_file_clust)
    return PSI_db,new_header

def removeRedundantCluster(filename,clusterID_file):
    from scipy import stats
    import ExpressionBuilder
    sort_col=0
    export_count=0
    ### Sort the filtered PSI model by gene name
    ExpressionBuilder.exportSorted(filename, sort_col, excludeHeader=True)
    new_file = filename[:-4]+'-unique.txt'
    ea = export.ExportFile(new_file)
    
    event_clusterID_db={}
    for line in open(clusterID_file,'rU').xreadlines():
        data = line.rstrip()
        eventID,clusterID = string.split(data,'\t')
        event_clusterID_db[eventID]=clusterID

        def compareEvents(events_to_compare,export_count):
                ### This is where we compare the events and write out the unique entries
                if len(events_to_compare)==1:
                    ea.write(events_to_compare[0][-1])
                    export_count+=1
                else:
                    exclude={}
                    compared={}
                    for event1 in events_to_compare:
                        if event1[0] not in exclude:
                            ea.write(event1[-1])
                            exclude[event1[0]]=[]
                            export_count+=1
                        for event2 in events_to_compare:
                            if event2[0] not in exclude:
                                if event1[0] != event2[0] and (event1[0],event2[0]) not in compared:
                                    uid1,values1,line1 = event1
                                    uid2,values2,line2 = event2
                                    coefr=numpy.ma.corrcoef(values1,values2)
                                    #rho,p = stats.pearsonr(values1,values2)
                                    rho = coefr[0][1]
                                    if rho>0.9 or rho<-0.9:
                                        exclude[event2[0]]=[]
                                    compared[event1[0],event2[0]]=[]
                                    compared[event2[0],event1[0]]=[]
                                    
                    for event in events_to_compare:
                        if event[0] not in exclude:
                            ea.write(event[-1]) ### write out the line
                            exclude.append(event[0])
                            export_count+=1
                return export_count

    header = True
    rows=0
    filtered=0
    prior_cID = 0
    events_to_compare=[]
    for line in open(filename,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            ea.write(line)
            header_row = t
            header=False
        else:
            uid = t[0]
            cID = event_clusterID_db[uid]
            empty_offset = len(header_row)-len(t)
            t+=['']*empty_offset
            values = ['0.000101' if x=='' else x for x in t[1:]]
            values = map(float,values)
            values = numpy.ma.masked_values(values,0.000101)
            if prior_cID==0: prior_cID = cID ### Occurs for the first entry
            if cID == prior_cID:
                ### Replace empty values with 0
                events_to_compare.append((uid,values,line))
            else:
                export_count = compareEvents(events_to_compare,export_count)
                events_to_compare=[(uid,values,line)]
            prior_cID = cID

    if len(events_to_compare)>0: ### If the laster cluster set not written out yet
        export_count = compareEvents(events_to_compare,export_count)

    ea.close()
    print export_count,'Non-redundant splice-events exported'
    
if __name__ == '__main__':
    import getopt
    impute=False
    cutoff = 75
    returnValues = False
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        #python import_scripts/filterAndImputePSI.py --i /Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Krithika/Cancer-PSI-Counts/counts/AltAnalyze/AltResults/AlternativeOutput-min-10-reads/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt --percentExp 50 --impute yes
        print "Insufficient arguments";sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','percentExp=','impute=','returnValues='])
        for opt, arg in options:
            if opt == '--i': filename=arg
            elif opt == '--percentExp': cutoff=float(arg)
            elif opt == '--returnValues': returnValues = False
            elif opt == '--impute':
                if string.lower(arg) == 'true' or string.lower(arg) == 'yes':
                    impute = True
                else:
                    impute = False
            
    filterPSIValues(filename,impute=impute,cutoff=cutoff,returnValues=returnValues)
