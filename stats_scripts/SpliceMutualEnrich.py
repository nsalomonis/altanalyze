#!/usr/bin/env python
import numpy as np
import sys,string
import os
import os.path
from collections import defaultdict

try:from stats_scripts import statistics
except Exception: import statistics
import export; reload(export)
import re
from stats_scripts import fishers_exact_test
import traceback
import warnings
import math
import export

def importDenominator(denom_dir):
    denominator_events={}
    firstRow=True
    for line in open(denom_dir,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if firstRow:
            uid_index = t.index('UID')
            firstRow=False
        else:
            uid = t[uid_index]
            denominator_events[uid]=[]
    return denominator_events

def importEvents(folder,denom={}):
    ### For example, knockdown signatures
    import collections
    unique_events = {}
    import UI
    files = UI.read_directory(folder)
    comparison_events={}
    for file in files:
        if '.txt' in file and 'PSI.' in file:
            fn = folder+'/'+file
            firstLine = True
            comparison = file[:-4]
            comparison_events[comparison,'inclusion']=[]
            comparison_events[comparison,'exclusion']=[]
            for line in open(fn,'rU').xreadlines():
                data = line.rstrip()
                t = string.split(data,'\t')
                if firstLine:
                    try: event_index = t.index('Event-Direction')
                    except:
                        try: event_index = t.index('Inclusion-Junction') ### legacy
                        except: print file, 'Event-Direction error';sys.exit()
                    firstLine= False
                    continue
                event = t[0]
                #event = string.split(event,'|')[0]
                unique_events[event]=[]
                event_dictionary = comparison_events[comparison,t[event_index]]
                if len(denom)>0:
                    if event in denom:
                        event_dictionary.append(event)
                else:
                    event_dictionary.append(event)
    return unique_events,comparison_events
                
def performMutualEnrichment(unique_inp_events,event_inp_dictionary,unique_ref_events,event_ref_dictionary):
    N = len(unique_inp_events)
    N = 88000
    for (comparison,direction) in event_inp_dictionary:
        if direction == 'inclusion': alt_direction = 'exclusion'
        else: alt_direction = 'inclusion'
        comparison_events1 = event_inp_dictionary[(comparison,direction)]
        comparison_events2 = event_inp_dictionary[(comparison,alt_direction)]
        for (reference_comp,ref_direction) in event_ref_dictionary:
            if direction == ref_direction and direction == 'inclusion':
                if ref_direction == 'inclusion': alt_ref_direction = 'exclusion'
                else: alt_ref_direction = 'inclusion'
                ref_events1 = event_ref_dictionary[(reference_comp,ref_direction)]
                ref_events2 = event_ref_dictionary[(reference_comp,alt_ref_direction)]
                concordant1 = len(list(set(comparison_events1) & set(ref_events1)))
                concordant2 = len(list(set(comparison_events2) & set(ref_events2)))
                r1 = concordant1+concordant2
                n = len(ref_events1)+len(ref_events2)
                R = len(comparison_events1)+len(comparison_events2)
    
                disconcordant1 = len(list(set(comparison_events1) & set(ref_events2)))
                disconcordant2 = len(list(set(comparison_events2) & set(ref_events1)))
    
                r2 = disconcordant1+disconcordant2
                #n = r1+r2
                
                try: z_concordant = Zscore(r1,n,N,R)
                except ZeroDivisionError: z_concordant = 0.0000
                
                try: z_discordant = Zscore(r2,n,N,R)
                except ZeroDivisionError: z_discordant = 0.0000
                
                try: null_z = Zscore(0,n,N,R)
                except ZeroDivisionError: null_z = 0.000
                ### Calculate a Fischer's Exact P-value
                import mappfinder
                pval1 = mappfinder.FishersExactTest(r1,n,R,N)
                pval2 = mappfinder.FishersExactTest(r2,n,R,N)
                ### Store these data in an object
                #zsd = mappfinder.ZScoreData(signature,r,n,z,null_z,n)
                #zsd.SetP(pval)
                
                print comparison+'\t'+reference_comp+'\t'+ref_direction+'\t'+str(z_concordant)+'\t'+str(z_discordant)+'\t'+str(r2)+'\t'+str(n)+'\t'+str(pval1)+'\t'+str(pval2)
            
def Zscore(r,n,N,R):
    """where N is the total number of events measured: 
    R is the total number of events meeting the criterion:
    n is the total number of events in this specific reference gene-set: 
    r is the number of events meeting the criterion in the examined reference gene-set: """
    N=float(N) ### This bring all other values into float space
    z = (r - n*(R/N))/math.sqrt(n*(R/N)*(1-(R/N))*(1-((n-1)/(N-1))))
    return z               
                
if __name__ == '__main__':

    import getopt
  
    mutdict=defaultdict(list)
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','r=','d='])
        for opt, arg in options:
            if opt == '--i': input_directory=arg
            elif opt == '--r':reference_directory=arg
            elif opt == '--d': denominator_directory=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
                
    denominator_events = importDenominator(denominator_directory)
    unique_ref_events,event_ref_dictionary = importEvents(reference_directory,denom=denominator_events)
    #unique_inp_events,event_inp_dictionary = importEvents(input_directory,denom=unique_ref_events)
    unique_inp_events,event_inp_dictionary = importEvents(input_directory)
    performMutualEnrichment(unique_inp_events,event_inp_dictionary,unique_ref_events,event_ref_dictionary)
