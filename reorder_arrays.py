###reorder_arrays
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

import sys, string
import os.path
import unique
from stats_scripts import statistics
import math

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def reorderArrayHeaders(data_headers,array_order,comp_group_list,array_linker_db):
    ###array_order gives the final level order sorted, followed by the original index order as a tuple                   
    data_headers2 = {}; array_linker_db2 = {}; ranked_array_headers = []; raw_data_comps={}; group_name_db = {}

    for x in array_order:
        y = x[1]  ### this is the new first index
        group = x[2]; group_name = x[3]
        group_name_db[group] = group_name            
        ### for example y = 5, therefore the data[row_id][5] entry is now the first
        try: data_headers2[group].append(data_headers[y])
        except KeyError: data_headers2[group]= [data_headers[y]]

    raw_data_comp_headers = {}
    for comp in comp_group_list:
        temp_raw = []
        group1 = int(comp[0]);group2 = int(comp[1])
        comp = str(comp[0]),str(comp[1])
        g1_headers = data_headers2[group1]
        g2_headers = data_headers2[group2]
        g1_name = group_name_db[group1]
        g2_name = group_name_db[group2]
        for header in g2_headers: temp_raw.append(g2_name+':'+header)
        for header in g1_headers: temp_raw.append(g1_name+':'+header)
        raw_data_comp_headers[comp] = temp_raw
        
    for array_name in array_linker_db: array_linker_db2[array_linker_db[array_name]]=array_name
    ###Determine the number of arrays in each group for f-test analysis
    group_count={}
    for x in array_order:
        original_index = x[1]; group = x[2]; group_name = x[3]
        array_name = array_linker_db2[original_index]; ranked_array_headers.append(group_name+':'+array_name)
        try: group_count[group] += 1
        except KeyError: group_count[group] = 1
    group_count_list=[]; group_count_list2=[]
    for group_number in group_count:
        count = group_count[group_number]
        group_count_list.append((group_number,count))
    group_count_list.sort()
    #print group_count_list
    for (group_number,count) in group_count_list: group_count_list2.append(count)
            
    #return expbuilder_value_db,group_count_list2,ranked_array_headers,raw_data_comps,raw_data_comp_headers
    return group_count_list2,raw_data_comp_headers

def filterBlanks(data_list):
    data_list_new=[]
    for i in data_list:
        if i=='':pass
        else: data_list_new.append(i)
    return data_list_new

def reorder(data,data_headers,array_order,comp_group_list,probeset_db,include_raw_data,array_type,norm,fl,logvalues=True,blanksPresent=False):
    ###array_order gives the final level order sorted, followed by the original index order as a tuple                   
    expbuilder_value_db = {}; group_name_db = {}; summary_filtering_stats = {}; pval_summary_db= {}
    replicates = 'yes'
    
    stat_result_names = ['avg-','log_fold-','fold-','rawp-','adjp-']
    group_summary_result_names = ['avg-']
    
    ### Define expression variables
    try: probability_statistic = fl.ProbabilityStatistic()
    except Exception: probability_statistic = 'unpaired t-test'
    try: gene_exp_threshold = math.log(fl.GeneExpThreshold(),2)
    except Exception: gene_exp_threshold = 0
    try: gene_rpkm_threshold = float(fl.RPKMThreshold())
    except Exception: gene_rpkm_threshold = 0
    try: FDR_statistic = fl.FDRStatistic()
    except Exception: FDR_statistic = 'Benjamini-Hochberg'
    calculateAsNonLog=True
    if blanksPresent:
        calculateAsNonLog=False
    
    ### Begin processing sample expression values according to the organized groups
    for row_id in data:
        try: gene = probeset_db[row_id][0]
        except TypeError: gene = '' #not needed if not altsplice data
        data_headers2 = {} #reset each time
        grouped_ordered_array_list = {}
        for x in array_order:
            y = x[1]  #this is the new first index
            group = x[2]
            group_name = x[3]
            group_name_db[group] = group_name
            #for example y = 5, therefore the data[row_id][5] entry is now the first
            try:
                try: new_item = data[row_id][y]
                except IndexError: print row_id,data[row_id],len(data[row_id]),y,len(array_order),array_order;kill
                if logvalues==False and calculateAsNonLog and array_type == 'RNASeq':
                    new_item = math.pow(2,new_item)
            except TypeError: new_item = ''  #this is for a spacer added in the above function
            try: grouped_ordered_array_list[group].append(new_item)
            except KeyError: grouped_ordered_array_list[group] = [new_item]
            try: data_headers2[group].append(data_headers[y])
            except KeyError: data_headers2[group]= [data_headers[y]]
        #perform statistics on each group comparison - comp_group_list: [(1,2),(3,4)]
        stat_results = {}
        group_summary_results = {}
        for comp in comp_group_list:
            group1 = int(comp[0])
            group2 = int(comp[1])
            group1_name = group_name_db[group1]
            group2_name = group_name_db[group2]
            groups_name = group1_name + "_vs_" + group2_name
            data_list1 = grouped_ordered_array_list[group1] 
            data_list2 = grouped_ordered_array_list[group2] #baseline expression
            if blanksPresent: ### Allows for empty cells
                data_list1 = filterBlanks(data_list1)
                data_list2 = filterBlanks(data_list2)
            try: avg1 = statistics.avg(data_list1)
            except Exception: avg1 = ''
            try: avg2 = statistics.avg(data_list2)
            except Exception: avg2=''
            try:
                if (logvalues == False and array_type != 'RNASeq') or (logvalues==False and calculateAsNonLog):
                    fold = avg1/avg2
                    log_fold = math.log(fold,2)
                    if fold<1: fold = -1.0/fold
                else:
                    log_fold = avg1 - avg2
                    fold = statistics.log_fold_conversion(log_fold) 
            except Exception:
                log_fold=''; fold=''
            try:
                #t,df,tails = statistics.ttest(data_list1,data_list2,2,3) #unpaired student ttest, calls p_value function
                #t = abs(t); df = round(df); p = str(statistics.t_probability(t,df))
                p = statistics.runComparisonStatistic(data_list1,data_list2,probability_statistic)
            except Exception: p = 1; sg = 1; N1=0; N2=0
            comp = group1,group2
            if array_type == 'RNASeq': ### Also non-log but treated differently
                if 'RPKM' == norm: adj = 0
                else: adj = 1
                if calculateAsNonLog == False:
                    try: avg1 = math.pow(2,avg1)-adj; avg2 = math.pow(2,avg2)-adj
                    except Exception: avg1=''; avg2=''
                if 'RPKM' == norm:
                    if avg1 < gene_rpkm_threshold and avg2 < gene_rpkm_threshold:
                        log_fold = 'Insufficient Expression'
                        fold = 'Insufficient Expression'
                else:
                    if avg1 < gene_exp_threshold and avg2 < gene_exp_threshold:
                        log_fold = 'Insufficient Expression'
                        fold = 'Insufficient Expression'
                    #if row_id=='ENSG00000085514':
                    #if fold=='Insufficient Expression':
                    #print [norm, avg1, avg2, fold, comp, gene_exp_threshold, gene_rpkm_threshold, row_id]
                    #5.96999111075 7.72930768675 Insufficient Expression (3, 1) 1.0 ENSG00000085514
            if gene_rpkm_threshold!=0 and calculateAsNonLog: ### Any other data
                a1 = nonLogAvg(data_list1)
                a2 = nonLogAvg(data_list2)
                #print [a1,a2,gene_rpkm_threshold]
                if a1<gene_rpkm_threshold and a2<gene_rpkm_threshold:
                    log_fold = 'Insufficient Expression'
                    fold = 'Insufficient Expression'
                #print log_fold;kill
            try:
                gs = statistics.GroupStats(log_fold,fold,p)
                stat_results[comp] = groups_name,gs,group2_name
                if probability_statistic == 'moderated t-test':
                    gs.setAdditionalStats(data_list1,data_list2) ### Assuming equal variance
                if probability_statistic == 'moderated Welch-test':
                    gs.setAdditionalWelchStats(data_list1,data_list2) ### Assuming unequal variance
            except Exception:
                null=[]; replicates = 'no' ### Occurs when not enough replicates
                #print comp, len(stat_results); kill_program
            group_summary_results[group1] = group1_name,[avg1]
            group_summary_results[group2] = group2_name,[avg2]

        ### Replaces the below method to get the largest possible comparison fold and ftest p-value
        grouped_exp_data = []; avg_exp_data = []
        for group in grouped_ordered_array_list:
            data_list = grouped_ordered_array_list[group]
            if blanksPresent: ### Allows for empty cells
                data_list = filterBlanks(data_list)
            if len(data_list)>0: grouped_exp_data.append(data_list)
            try: avg = statistics.avg(data_list); avg_exp_data.append(avg)
            except Exception:
                avg = ''
                #print row_id, group, data_list;kill
        try: avg_exp_data.sort(); max_fold = avg_exp_data[-1]-avg_exp_data[0]
        except Exception: max_fold = 'NA'
        try: ftestp = statistics.OneWayANOVA(grouped_exp_data)
        except Exception: ftestp = 1
        gs = statistics.GroupStats(max_fold,0,ftestp)
        summary_filtering_stats[row_id] = gs
        
        stat_result_list = []
        for entry in stat_results:
            data_tuple = entry,stat_results[entry]
            stat_result_list.append(data_tuple)
        stat_result_list.sort()
        
        grouped_ordered_array_list2 = []
        for group in grouped_ordered_array_list:
            data_tuple = group,grouped_ordered_array_list[group]
            grouped_ordered_array_list2.append(data_tuple)
        grouped_ordered_array_list2.sort() #now the list is sorted by group number
        
        ###for each rowid, add in the reordered data, and new statistics for each group and for each comparison
        for entry in grouped_ordered_array_list2:
            group_number = entry[0]
            original_data_values = entry[1]
            if include_raw_data == 'yes': ###optionally exclude the raw values
                for value in original_data_values:
                    if array_type == 'RNASeq':
                        if norm == 'RPKM': adj = 0
                        else: adj = 1
                        if calculateAsNonLog == False:
                            value = math.pow(2,value)-adj
                    try: expbuilder_value_db[row_id].append(value)
                    except KeyError: expbuilder_value_db[row_id] = [value]
            if group_number in group_summary_results:
                group_summary_data = group_summary_results[group_number][1] #the group name is listed as the first entry
                for value in group_summary_data:
                    try: expbuilder_value_db[row_id].append(value)
                    except KeyError: expbuilder_value_db[row_id] = [value]
            for info in stat_result_list:
                if info[0][0] == group_number: #comp,(groups_name,[avg1,log_fold,fold,ttest])
                    comp = info[0]; gs = info[1][1]
                    expbuilder_value_db[row_id].append(gs.LogFold())
                    expbuilder_value_db[row_id].append(gs.Fold())
                    expbuilder_value_db[row_id].append(gs.Pval())
                    ### Create a placeholder and store the position of the adjusted p-value to be calculated
                    expbuilder_value_db[row_id].append('') 
                    gs.SetAdjPIndex(len(expbuilder_value_db[row_id])-1)
                    gs.SetPvalIndex(len(expbuilder_value_db[row_id])-2)
                    pval_summary_db[(row_id,comp)] = gs

    ###do the same for the headers, but at the dataset level (redundant processes)
    array_fold_headers = []; data_headers3 = []
    try:
        for group in data_headers2:
            data_tuple = group,data_headers2[group]  #e.g. 1, ['X030910_25_hl.CEL', 'X030910_29R_hl.CEL', 'X030910_45_hl.CEL'])
            data_headers3.append(data_tuple)
        data_headers3.sort()
    except UnboundLocalError:
        print data_headers,'\n',array_order,'\n',comp_group_list,'\n'; kill_program
    
    for entry in data_headers3:
        x = 0 #indicates the times through a loop
        y = 0 #indicates the times through a loop
        group_number = entry[0]
        original_data_values = entry[1]
        if include_raw_data == 'yes': ###optionally exclude the raw values
            for value in original_data_values:
                array_fold_headers.append(value)
        if group_number in group_summary_results:
            group_name = group_summary_results[group_number][0]
            group_summary_data = group_summary_results[group_number][1]
            for value in group_summary_data:
                combined_name = group_summary_result_names[x] + group_name  #group_summary_result_names = ['avg-']
                array_fold_headers.append(combined_name)
                x += 1 #increment the loop index

        for info in stat_result_list:
            if info[0][0] == group_number:  #comp,(groups_name,[avg1,log_fold,fold,ttest],group2_name)
                groups_name = info[1][0]
                only_add_these = stat_result_names[1:]
                for value in only_add_these:
                    new_name = value + groups_name
                    array_fold_headers.append(new_name)

    ###For the raw_data only export we need the headers for the different groups (data_headers2) and group names (group_name_db)       
    raw_data_comp_headers = {}
    for comp in comp_group_list:
        temp_raw = []
        group1 = int(comp[0]);group2 = int(comp[1])
        comp = str(comp[0]),str(comp[1])
        g1_headers = data_headers2[group1]
        g2_headers = data_headers2[group2]
        g1_name = group_name_db[group1]
        g2_name = group_name_db[group2]
        for header in g2_headers: temp_raw.append(g2_name+':'+header)
        for header in g1_headers: temp_raw.append(g1_name+':'+header)
        raw_data_comp_headers[comp] = temp_raw

    ###Calculate adjusted ftest p-values using BH95 sorted method
    statistics.adjustPermuteStats(summary_filtering_stats)
    
    ### Calculate adjusted p-values for all p-values using BH95 sorted method
    round=0
    for info in comp_group_list:
        compid = int(info[0]),int(info[1]); pval_db={}
        for (rowid,comp) in pval_summary_db:
            if comp == compid:
                gs = pval_summary_db[(rowid,comp)]
                pval_db[rowid] = gs

        if 'moderated' in probability_statistic and replicates == 'yes':
            ### Moderates the original reported test p-value prior to adjusting
            try: statistics.moderateTestStats(pval_db,probability_statistic)
            except Exception:
                if round == 0:
                    if replicates == 'yes':
                        print 'Moderated test failed due to issue with mpmpath or out-of-range values\n   ... using unmoderated unpaired test instead!'
                null=[] ### Occurs when not enough replicates
            round+=1
            
        if FDR_statistic == 'Benjamini-Hochberg':
            statistics.adjustPermuteStats(pval_db)
        else:
            ### Calculate a qvalue (https://github.com/nfusi/qvalue)
            import numpy; from stats_scripts import qvalue; pvals = []; keys = []
            for key in pval_db: pvals.append(pval_db[key].Pval()); keys.append(key)
            pvals = numpy.array(pvals)
            pvals = qvalue.estimate(pvals)
            for i in range(len(pvals)): pval_db[keys[i]].SetAdjP(pvals[i])
            
        for rowid in pval_db:
            gs = pval_db[rowid]
            expbuilder_value_db[rowid][gs.AdjIndex()] = gs.AdjP() ### set the place holder to the calculated value
            if 'moderated' in probability_statistic:
                expbuilder_value_db[rowid][gs.RawIndex()] = gs.Pval() ### Replace the non-moderated with a moderated p-value
                
    pval_summary_db=[]            
    ###Finished re-ordering lists and adding statistics to expbuilder_value_db
    return expbuilder_value_db, array_fold_headers, summary_filtering_stats, raw_data_comp_headers

def nonLogAvg(data_list):
    return statistics.avg(map(lambda x: math.pow(2,x)-1,data_list))

if __name__ == '__main__':
    print array_cluster_final