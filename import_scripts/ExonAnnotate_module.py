###ExonAnalyze_module
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
import unique
import copy

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list
            
def identifyPutativeSpliceEvents(exon_db,constituitive_probeset_db,array_id_db,agglomerate_inclusion_probesets,onlyAnalyzeJunctions):
    exon_dbase = {}; probeset_comparison_db = {}; x = 0; y = 0
    ### Grab all probesets where we can identify a potential exon inclusion/exclusion event
    if len(array_id_db) == 0: array_id_db = exon_db ### Used when exporting all comparitive junction data
        
    for probeset in array_id_db:
        if probeset in exon_db:
            affygene = exon_db[probeset].GeneID() #exon_db[probeset] = affygene,exons,ensembl,block_exon_ids,block_structure,comparison_info    
            exons = exon_db[probeset].ExonID() #get rid of last pipe
            if probeset not in constituitive_probeset_db:
              #thus, there is a 'gene' probeset for that gene, but we don't want to look at the gene probesets
              if '|' not in exons: #get rid of any block exons or ambiguities) 
                try: x += 1; probeset_comparison_db[affygene].append(exons)
                except KeyError: x += 1; probeset_comparison_db[affygene] = [exons]
            exon_dbase[affygene,exons] = probeset

    print "Number of putative probeset comparisons:",x
    
    probe_level_db = {}
    for affygene in probeset_comparison_db:
        for exon_probeset1 in probeset_comparison_db[affygene]:
            for exon_probeset2 in probeset_comparison_db[affygene]:
                if exon_probeset1 != exon_probeset2:
                    if '-' in exon_probeset1: #get both pair-wise possibilities with this, to grab junctions
                        e1a,e1b = string.split(exon_probeset1,'-')
                        e1 = e1a,e1b
                        try:
                            e2a,e2b = string.split(exon_probeset2,'-')
                            e2 = e2a,e2b
                        except ValueError: e2 = exon_probeset2
                        try: probe_level_db[affygene,e1].append(e2)
                        except KeyError: probe_level_db[affygene,e1] = [e2]
                    else: ### Required when exon_probeset1 is a single exon rather than a junction
                        if '-' in exon_probeset2:
                            e2a,e2b = string.split(exon_probeset2,'-')
                            e2 = e2a,e2b
                            e1 = exon_probeset1
                            try: probe_level_db[affygene,e2].append(e1)
                            except KeyError: probe_level_db[affygene,e2] = [e1]
    #print "Looking for exon events defined by probeset exon associations"
    alt_junction_db,critical_exon_db = independently_rank_analyze_junction_sets(probe_level_db,onlyAnalyzeJunctions)
    #print "Associations Built\n"

    ### Rearange alt_junction_db and agglomerate data for inclusion probesets
    exon_inclusion_db={}; exon_inclusion_event_db={}; alt_junction_db_collapsed={}
    if agglomerate_inclusion_probesets == 'yes':
        for affygene in alt_junction_db:
            alt_junction_db[affygene].sort() ### Should be no need to sort later if we do this
            for event in alt_junction_db[affygene]:
                ### event = [('ei', 'E16-E17'), ('ex', 'E16-E18')]
                event1 = event[0][0]; exon_set1 = event[0][1]; exon_set2 = event[1][1]            
                probeset1 = exon_dbase[affygene,exon_set1]; probeset2 = exon_dbase[affygene,exon_set2]
                if event1 == 'ei':
                    ###First generate the original fold values for export summary, then the adjusted
                    try: exon_inclusion_db[probeset2].append(probeset1)
                    except KeyError: exon_inclusion_db[probeset2] = [probeset1]
                    try: exon_inclusion_event_db[(affygene, probeset2, event[1])].append(event)
                    except KeyError: exon_inclusion_event_db[(affygene, probeset2, event[1])] = [event]
                else: ### Store all the missing mutual exclusive splicing events
                    try: alt_junction_db_collapsed[affygene].append(event)
                    except KeyError: alt_junction_db_collapsed[affygene] = [event]
        
        ###Create a new alt_junction_db with merged inclusion events
        for key in exon_inclusion_event_db:
            affygene = key[0]; excl_probeset=key[1]; excl_event = key[2]
            ###Collect critical exon information from each inclusion exon-set to agglomerate and delete old entries
            new_critical_exon_list=[]; incl_exon_sets=[]
            for event in exon_inclusion_event_db[key]:
                incl_exon_set = event[0][1]; incl_exon_sets.append(incl_exon_set) ### Don't sort since this will throw off probeset relationships: incl_exon_sets.sort()
                if len(exon_inclusion_event_db[key])>1:  ###If the original list of events > 1
                    critical_exon_list = critical_exon_db[affygene,tuple(event)][1]
                    for exon in critical_exon_list: new_critical_exon_list.append(exon)
                    #del critical_exon_db[affygene,tuple(event)]
            new_critical_exon_list = unique.unique(new_critical_exon_list); new_critical_exon_list.sort()
            new_critical_exon_list = [1,new_critical_exon_list]
            incl_exon_sets_str = string.join(incl_exon_sets,'|') ### New inclusion exon group
            event = [('ei',incl_exon_sets_str),excl_event] ### Store new inclusion exon group
            try: alt_junction_db_collapsed[affygene].append(event)
            except KeyError: alt_junction_db_collapsed[affygene] = [event]
            ###Replace exon_dbase entries with new combined probeset IDs
            incl_probesets = exon_inclusion_db[excl_probeset]
            incl_probesets_str = string.join(incl_probesets,'|')
            if len(incl_exon_sets)>1: ###Often there will be only a single inclusion probeset
                """for exons in incl_exon_sets:
                    key = affygene,exons
                    try: del exon_dbase[key] ###delete individual inclusion exons and replace with a single inclusion agglomerate
                    except KeyError: continue ###Can occur more than once, if an exon participates in more than one splicing event
                """
                exon_dbase[affygene,incl_exon_sets_str] = incl_probesets_str
                critical_exon_db[affygene,tuple(event)] = new_critical_exon_list
                ###Create a new probeset entry in exon_db for the agglomerated probesets
                new_block_exon_ids=[] #exon_db[probeset] = affygene,exons,ensembl,block_exon_ids,block_structure
                for probeset in incl_probesets:
                    edat = exon_db[probeset]; ensembl = edat.ExternalGeneID(); block_exon_ids = edat.SecondaryExonID(); block_structure = edat.GeneStructure()
                    new_block_exon_ids.append(block_exon_ids)
                new_block_exon_ids = string.join(new_block_exon_ids,'')
                edat = exon_db[incl_probesets[0]]; edat1 = edat; edat1.setDisplayExonID(incl_exon_sets_str) #; edat1.setExonID(edat.ExonID()) ### Use the first inclusion probeset instance for storing all instance data
                edat1.setSecondaryExonID(new_block_exon_ids); edat1.setProbeset(incl_probesets[0])
                exon_db[incl_probesets_str] = edat1
        print "Length of original splice event database:",len(alt_junction_db)
        print "Length of agglomerated splice event database:",len(alt_junction_db_collapsed)
        alt_junction_db = alt_junction_db_collapsed  ### Replace with agglomerated database
        ### End Rearangement
        
    return alt_junction_db,critical_exon_db,exon_dbase,exon_inclusion_db,exon_db

def independently_rank_analyze_junction_sets(probe_level_db,onlyAnalyzeJunctions):
    ### The below code is used to identify sets of junctions and junction and exon sets anti-correlated with each other
    ### independently storing the critical exons involved
    #probe_level_db[affygene,exons1].append(exons2)
    x = 0
    critical_exon_db = {}
    alt_junction_db = {}
    probe_level_db = eliminate_redundant_dict_values(probe_level_db)
    for key in probe_level_db:
        critical_exon_list = []
        affygene = key[0]
        exon_pair1 = key[1]
        e1a = int(exon_pair1[0][1:])
        e1b = int(exon_pair1[1][1:])
        for exon_pair2 in probe_level_db[key]: #exon_pair2 could be a single exon
            s = 0 #moved this down!!!!!!!!!
            if exon_pair2[0] == 'E': # thus, exon_pair2 is actually a single exon
                e2 = int(exon_pair2[1:])
                s=1
            else:
                e2a = int(exon_pair2[0][1:])
                e2b = int(exon_pair2[1][1:])
            if s==0:
                e1_pair = e1a,e1b
                e2_pair = e2a,e2b
            if s==1: # thus, exon_pair2 is actually a single exon  
                e1_pair = e1a,e1b
                if e1a < e2 and e1b > e2 and onlyAnalyzeJunctions == 'no': # e.g. E3-E5 vs. E4
                    e1 = 'ex','E'+str(e1a)+'-'+'E'+str(e1b); e2x = 'ei','E'+str(e2) 
                    critical_exons = [e1,e2x]
                    critical_exons.sort(); critical_exon_list.append(critical_exons)
                    critical_exon_db[affygene,tuple(critical_exons)] = [1,['E'+str(e2)]] ###The 1 indicates that the exon can be called up or down, since it is an ei or ex event vs. mx
            ### Note: everything except for the last one should have two instances added to the database
            elif (e1b == e2b and e1a > e2a): # e.g. E2-E3 vs. E1-E3
                    e1 = 'ei','E'+str(e1a)+'-'+'E'+str(e1b); e2 = 'ex','E'+str(e2a)+'-'+'E'+str(e2b)
                    critical_exons = [e1,e2]
                    critical_exons.sort();critical_exon_list.append(critical_exons)
                    critical_exon_db[affygene,tuple(critical_exons)] = [1,['E'+str(e1a)]]
                    #print affygene, exon_pair1,e1a,e1b,'----',exon_pair2,e2a,e2b
            elif (e1b == e2b and e1a < e2a): # e.g. E1-E3 vs. E2-E3
                    e1 = 'ex','E'+str(e1a)+'-'+'E'+str(e1b); e2 = 'ei','E'+str(e2a)+'-'+'E'+str(e2b)
                    critical_exons = [e1,e2]
                    critical_exons.sort();critical_exon_list.append(critical_exons)
                    critical_exon_db[affygene,tuple(critical_exons)] = [1,['E'+str(e2a)]]
            elif (e1a == e2a and e1b < e2b): # e.g. E2-E3 vs. E2-E4
                    e1 = 'ei','E'+str(e1a)+'-'+'E'+str(e1b); e2 = 'ex','E'+str(e2a)+'-'+'E'+str(e2b)
                    critical_exons = [e1,e2]
                    critical_exons.sort();critical_exon_list.append(critical_exons)
                    critical_exon_db[affygene,tuple(critical_exons)] = [1,['E'+str(e1b)]]
            elif (e1a == e2a and e1b > e2b): # e.g. E2-E4 vs. E2-E3
                    e1 = 'ex','E'+str(e1a)+'-'+'E'+str(e1b); e2 = 'ei','E'+str(e2a)+'-'+'E'+str(e2b)
                    critical_exons = [e1,e2]
                    critical_exons.sort();critical_exon_list.append(critical_exons)
                    critical_exon_db[affygene,tuple(critical_exons)] = [1,['E'+str(e2b)]]
            elif (e1a < e2a and e1b > e2a) and (e1a < e2b and e1b > e2b): # e.g. E2-E6 vs. E3-E5
                    e1 = 'ex','E'+str(e1a)+'-'+'E'+str(e1b); e2 = 'ei','E'+str(e2a)+'-'+'E'+str(e2b)
                    critical_exons = [e1,e2]
                    critical_exons.sort();critical_exon_list.append(critical_exons)
                    critical_exon_db[affygene,tuple(critical_exons)] = [1,['E'+str(e2a),'E'+str(e2b)]]
            elif (e1a > e2a and e1b < e2a) and (e1a > e2b and e1b < e2b): # e.g. E3-E5 vs. E2-E6
                    e1 = 'ei','E'+str(e1a)+'-'+'E'+str(e1b); e2 = 'ex','E'+str(e2a)+'-'+'E'+str(e2b)
                    critical_exons = [e1,e2]
                    critical_exons.sort();critical_exon_list.append(critical_exons)
                    critical_exon_db[affygene,tuple(critical_exons)] = [1,['E'+str(e1a),'E'+str(e1b)]]
            elif (e1a < e2a and e1b > e2a): # e.g. E2-E6 vs. E3-E8
                    e1 = 'mx','E'+str(e1a)+'-'+'E'+str(e1b); e2 = 'mx','E'+str(e2a)+'-'+'E'+str(e2b)
                    critical_exons = [e1,e2]
                    critical_exon_list.append(critical_exons)
                    critical_exon_db[affygene,tuple(critical_exons)] = [2,['E'+str(e1b),'E'+str(e2a)]]
            elif (e1a < e2b and e1b > e2b): # e.g. E2-E6 vs. E1-E3
                    e1 = 'mx','E'+str(e1a)+'-'+'E'+str(e1b); e2 = 'mx','E'+str(e2a)+'-'+'E'+str(e2b)
                    critical_exons = [e1,e2]
                    critical_exon_list.append(critical_exons)
                    critical_exon_db[affygene,tuple(critical_exons)] = [2,['E'+str(e1a),'E'+str(e2b)]]
            if len(critical_exon_list)>0:
              for entry in critical_exon_list:
                try:
                    alt_junction_db[affygene].append(entry)
                except KeyError:
                    alt_junction_db[affygene] = [entry]
    alt_junction_db = eliminate_redundant_dict_values(alt_junction_db)

    return alt_junction_db, critical_exon_db

def exportJunctionComparisons(alt_junction_db,critical_exon_db,exon_dbase):
    competitive_junction_export = 'AltDatabase\Mm\AltMouse\AltMouse_junction-comparisons.txt'
    fn=filepath(competitive_junction_export)
    data = open(fn,'w')
    title = ['Affygene','probeset1','probeset2','critical-exons']; title = string.join(title,'\t')+'\n'; data.write(title)
    for affygene in alt_junction_db:
        alt_junction_db[affygene].sort() ### Should be no need to sort later if we do this
        for event in alt_junction_db[affygene]:
            ### event = [('ei', 'E16-E17'), ('ex', 'E16-E18')]
            exon_set1 = event[0][1]; exon_set2 = event[1][1]            
            probeset1 = exon_dbase[affygene,exon_set1]; probeset2 = exon_dbase[affygene,exon_set2]
            critical_exon_list = critical_exon_db[affygene,tuple(event)];critical_exon_list = critical_exon_list[1]
            critical_exon_list = string.join(critical_exon_list,'|')
            export_data = string.join([affygene]+[probeset1,probeset2,critical_exon_list],'\t')+'\n'
            data.write(export_data)
    data.close()

def annotate_splice_event(exons1,exons2,block_structure):
    #1(E13|E12)-2(E11)-3(E10)-4(E9|E8)-5(E7|E6|E5)-6(E4|E3|E2|E1)
    splice_event = ''; evidence = 'clear'
    string.replace(block_structure,')','|')
    block_list = string.split(block_structure,'-')
    #[1(E13|E12|,2(E11|,3(E10|,4(E9|E8|,5(E7|E6|E5|,6(E4|E3|E2|E1|]
    ###Perform a membership query
    try: exon1a,exon1b = string.split(exons1,'-') ###***
    except ValueError: exon1a = exons1; exon1b = exons1; evidence = 'unclear'
    try: exon2a,exon2b = string.split(exons2,'-')  
    except ValueError: exon2a = exons2; exon2b = exons2; evidence = 'unclear'
    a = '';b = '';block_a='';block_b=''
    if exon1a == exon2a: a = 'same'
    if exon1b == exon2b: b = 'same'
    ex1a_m = exon_membership(exon1a,block_list);ex2a_m = exon_membership(exon2a,block_list)
    ex1b_m = exon_membership(exon1b,block_list);ex2b_m = exon_membership(exon2b,block_list)
    #print ex1a_m, ex2a_m,ex1b_m,ex2b_m;dog
    if ex1a_m == ex2a_m: block_a = 'same'
    if ex1b_m == ex2b_m: block_b = 'same'
    ### Correct for strand differences
    strand = "+"
    if ex1a_m > ex1b_m:  #strand therefore is negative
        strand = "-"
    if (abs(ex1a_m - ex2a_m) == 1) or (abs(ex1b_m - ex2b_m) == 1): alternative_exons = 'one'
    else: alternative_exons = 'multiple'
    if (ex1a_m == -1) or (ex2a_m == -1) or (ex1b_m == -1) or (ex2b_m == -1): splice_event = "retained_intron"
    elif block_a == 'same' and b == 'same': splice_event = "alt5'"
    elif block_b == 'same' and a == 'same': splice_event = "alt3'"
    elif (block_a == 'same' and block_b != 'same'):
        if a == 'same':
            if alternative_exons == 'one': splice_event = "cassette-exon"
            else: splice_event = "cassette-exons"
        else:
            if alternative_exons == 'one': splice_event = "alt5'-cassette-exon"
            else: splice_event = "alt5'-cassette-exons"
    elif (block_b == 'same' and block_a != 'same'):
        if b == 'same':
            if alternative_exons == 'one': splice_event = "cassette-exon"
            else: splice_event = "cassette-exons"
        else:
            if alternative_exons == 'one': splice_event = "cassette-exon-alt3'"
            else: splice_event = "cassette-exons-alt3'"
    else:
        if alternative_exons == 'one': splice_event = "alt5'-cassette-exon-alt3'"
        else: splice_event = "alt5'-cassette-exons-alt3'"
    if evidence == 'unclear':
        ###If the first probeset is a junction and the second is an exon, are junction exons 2 blocks way
        if (abs(ex1a_m - ex2a_m) == 1) and (abs(ex1b_m - ex2b_m) == 1): splice_event = "cassette-exon"
        elif (block_a == 'same' and block_b != 'same'):
            if alternative_exons == 'one': splice_event = "alt5'"
            else: splice_event = "alt5'-cassette-exons"
        elif (block_a != 'same' and block_b == 'same'):
            if alternative_exons == 'one': splice_event = "alt3'"
            else: splice_event = "cassette-exons-alt3'"
        else: splice_event = "unclear"
        
    if strand == "-":
        if splice_event == "alt5'": splice_event = "alt3'"
        elif splice_event == "alt3'": splice_event = "alt5'"
        elif splice_event == "alt5'-cassette-exon": splice_event = "cassette-exon-alt3'"
        elif splice_event == "alt5'-cassette-exons": splice_event = "cassette-exons-alt3'"
        elif splice_event == "cassette-exons-alt3'": splice_event = "alt5'-cassette-exons"
        elif splice_event == "cassette-exon-alt3'": splice_event = "alt5'-cassette-exon"
    #print splice_event
    return splice_event

def exon_membership(exon,block_structure):
    i=0; x = -1
    exon_temp1 = exon+'|'; exon_temp2 = exon+')'
    for exon_block in block_structure:
        if exon_temp1 in exon_block or exon_temp2 in exon_block:
            x = i
        i += 1
    return x

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

if __name__ == '__main__':
    exons1 = 'E9-E12'
    exons2 = 'E11-E15'
    block_structure = '1(E1)-2(E2)-3(E3|E4)-4(E5)-5(E6)-6(E7|E8|E9|E10|E11)-7(E12)-8(E13|E14)-9(E15)-10(E16|E17)-11(E18)-12(E19|E20)-13(E21|E22)-14(E23|E24)'
    a = annotate_splice_event(exons1,exons2,block_structure)
    print a
    #alt_junction_db,critical_exon_db,exon_dbase,exon_inclusion_db = identifyPutativeSpliceEvents(exon_db,constituitive_probeset_db,agglomerate_inclusion_probesets)