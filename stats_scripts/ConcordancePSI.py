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
import UI

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data
                      
def compareEventLists(folder,minimumOverlap=10):
    import collections
    event_db = collections.OrderedDict()
    groups_list=['']
    files = UI.read_directory(folder)
    file_headers = {}
    for file in files:
        if '.txt' in file and 'PSI.' in file or '__' in file:
            ls={}
            event_db[file[:-4]]=ls
            groups_list.append(file[:-4])
            fn = folder+'/'+file
            firstLine = True
            for line in open(fn,'rU').xreadlines():
                data = line.rstrip()
                t = string.split(data,'\t')
                if firstLine:
                    file_headers[file[:-4]] = t ### Store the headers
                    cid = t.index('ClusterID')
                    try: event_index = t.index('Event-Direction')
                    except:
                        try: event_index = t.index('Inclusion-Junction') ### legacy
                        except:
                            print [file]
                            print 'Event-Direction error'; sys.exit()
                    firstLine= False
                    continue
                uid = t[0]
                uid = string.split(uid,'|')[0]
                #uid = t[cid]
                if 'U2AF1-l' in file or 'U2AF1-E' in file:
                    if t[2] == "inclusion":
                        ls[(uid,t[event_index])]=t ### Keep the event data for output
                else:
                    ls[(uid,t[event_index])]=t ### Keep the event data for output

    def convertEvents(events):
        opposite_events=[]
        for (event,direction) in events:
            if direction == 'exclusion':
                direction = 'inclusion'
            else:
                direction = 'exclusion'
            opposite_events.append((event,direction))
        return opposite_events
        
    ea1 = export.ExportFile(folder+'/overlaps-same-direction.txt')
    ea2 = export.ExportFile(folder+'/overlaps-opposite-direction.txt')
    ea3 = export.ExportFile(folder+'/concordance.txt')
    #ea4 = export.ExportFile(folder+'/overlap-same-direction-events.txt')
    ea1.write(string.join(groups_list,'\t')+'\n')
    ea2.write(string.join(groups_list,'\t')+'\n')
    ea3.write(string.join(groups_list,'\t')+'\n')
    
    comparison_db={}
    best_hits={}
    for comparison1 in event_db:
        events1 = event_db[comparison1]
        hits1=[comparison1]
        hits2=[comparison1]
        hits3=[comparison1]
        best_hits[comparison1]=[]
        for comparison2 in event_db:
            events2 = event_db[comparison2]
            events3 = convertEvents(events2)
            overlapping_events = list(set(events1).intersection(events2))                  
            overlap = len(overlapping_events)
            inverse_overlap = len(set(events1).intersection(events3)) ### Get opposite events
            ### Calculate ratios based on the size of the smaller set
            min_events1 = min([len(events1),len(events2)]) 
            min_events2 = min([len(events1),len(events3)])
            denom = overlap+inverse_overlap
            if denom == 0: denom = 0.00001
            #comparison_db[comparison1,comparison2]=overlap
            if min_events1 == 0: min_events1 = 1
            if (overlap+inverse_overlap)<minimumOverlap:
                hits1.append('0.5')
                hits2.append('0.5')
                hits3.append('0.5|0.5')
            else:
                hits1.append(str((1.00*overlap)/min_events1))
                hits2.append(str((1.00*inverse_overlap)/min_events1))
                hits3.append(str(1.00*overlap/denom)+'|'+str(1.00*inverse_overlap/denom)+':'+str(overlap+inverse_overlap))
                if 'Leu' not in comparison2:
                    comp_name = string.split(comparison2,'_vs')[0]
                    best_hits[comparison1].append([abs(1.00*overlap/denom),'cor',comp_name])
                    best_hits[comparison1].append([abs(1.00*inverse_overlap/denom),'anti',comp_name])
            if comparison1 != comparison2:
                if len(overlapping_events)>0:
                    #ea4.write(string.join(['UID',comparison1]+file_headers[comparison1]+[comparison2]+file_headers[comparison2],'\t')+'\n')
                    pass
                overlapping_events.sort()
                for event in overlapping_events:
                    vals = string.join([event[0],comparison1]+event_db[comparison1][event]+[comparison2]+event_db[comparison2][event],'\t')
                    #ea4.write(vals+'\n')
                    pass
        ea1.write(string.join(hits1,'\t')+'\n')
        ea2.write(string.join(hits2,'\t')+'\n')
        ea3.write(string.join(hits3,'\t')+'\n')
    ea1.close()
    ea2.close()
    ea3.close()
    #ea4.close()
    for comparison in best_hits:
        best_hits[comparison].sort()
        best_hits[comparison].reverse()
        hits = best_hits[comparison][:10]
        hits2=[]
        for (score,dir,comp) in hits:
            h = str(score)[:4]+'|'+dir+'|'+comp
            hits2.append(h)
        print comparison,'\t',string.join(hits2,', ')

def summarizePSIresults(folder, gene_class_file,dPSI_cutoff=0.1):
    if gene_class_file != None:
        TFs = simpleListImport(gene_class_file)
    else:
        TFs = None
        
    ### Import PSI results and report number of impacted TFs
    files = UI.read_directory(folder)
    eo = export.ExportFile(folder + '/CommonEvents.txt')
    events_direction = {}
    
    for file in files:
        TFs_in_file = []
        filename = folder + '/' + file
        if '.txt' in file and 'PSI.' in file:
            condition = string.split(file[4:],'-')[0]
            header = True
            count = 0
            for line in open(filename, 'rU').xreadlines():
                if header:
                    data = cleanUpLine(line)
                    t = string.split(data, '\t')
                    uidd = t.index('UID')
                    cid = t.index('ClusterID')
                    dpsid = t.index('dPSI')
                    adjpid = t.index('adjp')
                    try: event_index = t.index('Event-Direction')
                    except:
                        try: event_index = t.index('Inclusion-Junction') ### legacy
                        except: print file, 'Event-Direction error';sys.exit()
                    header = False
                else:
                    data = cleanUpLine(line)
                    t = string.split(data, '\t')
                    uid = t[uidd]
                    uid = string.split(uid,'|')[0]
                    symbol = string.split(uid, ':')[0]
                    dPSI = abs(float(t[dpsid]))
                    direction = t[event_index]
    
                    if TFs != None:
                        if symbol in TFs and symbol not in TFs_in_file and dPSI > dPSI_cutoff:
                            TFs_in_file.append(symbol)
                            if symbol not in all_TFs:
                                all_TFs.append(symbol)
                            count += 1
                            print file+"\t"+t[-5]+"\t"+t[-4]+"\t"+t[0]
                            ##print file, count, len(all_TFs), string.join(TFs_in_file, ',')
                    else:
                        if uid in events_direction:
                            direction_db = events_direction[uid]
                            try: direction_db[direction].append(condition)
                            except: direction_db[direction] = [condition]
                        else:
                            direction_db={}
                            direction_db[direction] = [condition]
                            events_direction[uid] = direction_db
                            
    eo.write(string.join(["event","incl","incl_count","excl","excl_count"],'\t')+'\n')
    for event in events_direction:
        direction_db = events_direction[event]
        try:
            excl = string.join(direction_db['exclusion'],"|")
            excl_count = len(direction_db['exclusion'])
        except:
            excl = ''
            excl_count = 0
        try:
            incl = string.join(direction_db['inclusion'],"|")
            incl_count = len(direction_db['inclusion'])
        except:
            incl = ''
            incl_count = 0
        eo.write(string.join([event,incl,str(incl_count),excl,str(excl_count)],'\t')+'\n')
    eo.close()

if __name__ == '__main__':
    import getopt
    minimumOverlap=10
    dPSI_cutoff=01
    concordance = True
    gene_class_file=None
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
      print "Insufficient arguments";sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','concordance=','overlap=','genes=','concordance='])
        for opt, arg in options:
            if opt == '--i': folder=arg
            elif opt == '--overlap': minimumOverlap=float(arg)
            elif opt == '--genes': gene_class_file = arg
            elif opt == '--dPSI': dPSI_cutoff = float(arg)
            elif opt == '--concordance':
                if string.lower(arg) == 'true' or string.lower(arg) == 'yes':
                    concordance = True
                else:
                    concordance = False
                    
    if concordance:
        compareEventLists(folder,minimumOverlap=10)
    else:
        summarizePSIresults(folder,gene_class_file,dPSI_cutoff)
