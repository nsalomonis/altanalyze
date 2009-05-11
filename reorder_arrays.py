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
import statistics
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
            
    #return data2,group_count_list2,ranked_array_headers,raw_data_comps,raw_data_comp_headers
    return group_count_list2,raw_data_comp_headers
   
def reorderArraysOnly(data,array_order,comp_group_list):
    ###array_order gives the final level order sorted, followed by the original index order as a tuple                   
    data2 = {}
    for row_id in data:
        grouped_ordered_array_list = {}
        data_headers2 = {} #reset each time
        for x in array_order:
            y = x[1]  ### this is the new first index
            group = x[2]   
            ### for example y = 5, therefore the data[row_id][5] entry is now the first
            try: new_item = data[row_id][y]
            except TypeError: print y,x,array_order; kill
            try: data2[row_id].append(new_item)
            except KeyError: data2[row_id] = [new_item]
            ###Used for comparision analysis
            try: grouped_ordered_array_list[group].append(new_item)
            except KeyError: grouped_ordered_array_list[group] = [new_item]
        ###*******Include a database with the raw values saved for permuteAltAnalyze*******
        for info in comp_group_list:
            group1 = int(info[0]); group2 = int(info[1]); comp = str(info[0]),str(info[1])
            g1_data = grouped_ordered_array_list[group1]
            g2_data = grouped_ordered_array_list[group2]
            print g1_data
            data = comparision_export_db[comp]
            values = [row_id]+g1_data+g2_data; values = string.join(values,'\t')+'\n'
            #raw_data_comps[row_id,comp] = temp_raw
            data.write(values)

def reorder(data,data_headers,array_order,comp_group_list,probeset_db,include_raw_data):
    ###array_order gives the final level order sorted, followed by the original index order as a tuple                   
    data2 = {}; stat_summary_data = {}; group_name_db = {}; summary_filtering_stats = {}; raw_data_comps = {}
    
    stat_result_names = ['avg-','log_fold-','fold-','ttest-']
    group_summary_result_names = ['avg-']
    
    for row_id in data:
        try: affygene = probeset_db[row_id][0]
        except TypeError: affygene = '' #not needed if not altsplice data
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
                except IndexError: print row_id, len(data[row_id]),y,len(array_order),array_order;kill
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
            avg1 = statistics.avg(data_list1)
            try: avg2 = statistics.avg(data_list2)
            except ValueError: print data_list2,row_id
            #add std-error later
            log_fold = avg1 - avg2
            fold = statistics.log_fold_conversion(log_fold)
            try:
                t,df,tails = statistics.ttest(data_list1,data_list2,2,3) #unpaired student ttest, calls p_value function
                t = abs(t)
                df = round(df) #Excel doesn't recognize fractions in a DF
                ttest = '=tdist('+str(t)+','+str(df)+','+str(tails)+')'
                p = str(statistics.t_probability(t,df))
            except Exception: p = 1
            comp = group1,group2
            try:
                stat_results[comp] = groups_name,[avg2,log_fold,fold,p],group2_name #this structure is a little weird, since we will use this data for two differnent output files
            except TypeError:
                print comp, len(stat_results), dog
            group_summary_results[group1] = group1_name,[avg1]
            group_summary_results[group2] = group2_name,[avg2]

        stat_result_list = []
        p_list = []
        fold_list = []
        for entry in stat_results:
            data_tuple = entry,stat_results[entry]
            stat_result_list.append(data_tuple)
            ### add in script to grab the smallest p and larges fold for all comparisons
            log_fold = abs(stat_results[entry][1][1])
            p = stat_results[entry][1][3]
            p_list.append(float(p))
            fold_list.append(float(log_fold))
        stat_result_list.sort()
        p_list.sort()
        p_list.sort(); fold_list.sort(); fold_list.reverse()
        summary_filtering_stats[row_id] = p_list[0],fold_list[0]

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
                    try: data2[row_id].append(value)
                    except KeyError: data2[row_id] = [value]
            if group_number in group_summary_results:
                group_summary_data = group_summary_results[group_number][1] #the group name is listed as the first entry
                for value in group_summary_data:
                    try: data2[row_id].append(value)
                    except KeyError: data2[row_id] = [value]
            for info in stat_result_list:
                if info[0][0] == group_number: #comp,(groups_name,[avg1,log_fold,fold,ttest])
                    comp = info[0]
                    only_add_these = info[1][1][1:] #don't add the first value which is used in the second output file(avg)
                    for value in only_add_these:
                        try: data2[row_id].append(value)
                        except KeyError: data2[row_id] = [value]
                    ###Create a second database storing just the statistical information
                    groups_name = info[1][0]
                    stat_summary_data[row_id,comp] = info[1][1] #use this to output multiple files for input into altanalyze later
                    
        ###*******Include a database with the raw values saved for permuteAltAnalyze*******
        for info in comp_group_list:
            temp_raw = []
            group1 = int(info[0]); group2 = int(info[1]); comp = str(info[0]),str(info[1])
            g1_data = grouped_ordered_array_list[group1]
            g2_data = grouped_ordered_array_list[group2]
            for value in g2_data: temp_raw.append(value)
            for value in g1_data: temp_raw.append(value)
            raw_data_comps[row_id,comp] = temp_raw

    ###do the same for the headers, but at the dataset level (redundant processes)
    data_headers4 = []
    data_headers3 = []
    stat_summary_data_names = {}
    try:
        for group in data_headers2:
            data_tuple = group,data_headers2[group]  #e.g. 1, ['X030910_25_hl.CEL', 'X030910_29R_hl.CEL', 'X030910_45_hl.CEL'])
            data_headers3.append(data_tuple)
        data_headers3.sort()
    except UnboundLocalError:
        print data_headers,'\n'
        print array_order,'\n'
        print comp_group_list,'\n'
        kill
    
    for entry in data_headers3:
        x = 0 #indicates the times through a loop
        y = 0 #indicates the times through a loop
        group_number = entry[0]
        original_data_values = entry[1]
        if include_raw_data == 'yes': ###optionally exclude the raw values
            for value in original_data_values:
                data_headers4.append(value)
        if group_number in group_summary_results:
            group_name = group_summary_results[group_number][0]
            group_summary_data = group_summary_results[group_number][1]
            for value in group_summary_data:
                combined_name = group_summary_result_names[x] + group_name  #group_summary_result_names = ['avg-']
                data_headers4.append(combined_name)
                x += 1 #increment the loop index
        temp = []
        for info in stat_result_list:
            #print group_number, info[0]   
            if info[0][0] == group_number:  #comp,(groups_name,[avg1,log_fold,fold,ttest],group2_name)
                groups_name = info[1][0]
                group2_name = info[1][2]
                comp = info[0]
                comp = str(comp[0]),str(comp[1])
                baseline_avg_name = stat_result_names[0] + group2_name
                temp.append(baseline_avg_name) #add avg
                only_add_these = stat_result_names[1:]
                for value in only_add_these:
                    new_name = value + groups_name
                    data_headers4.append(new_name)
                    temp.append(new_name)
                stat_summary_data_names[comp] = temp

    ###For the raw_data only export we need the headers for the different groups (data_headers2), the comparisons (last stat_result_list entry is fine) and group names (group_name_db)       

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
        
    #print data_headers4, stat_summary_data_names
    """for entry in data2:
        print entry
        print data2[entry],dog """

    ###Finished re-ordering lists and adding statistics to data2 and stat_summary_data
    return data2, data_headers4, stat_summary_data, stat_summary_data_names, summary_filtering_stats, raw_data_comp_headers, raw_data_comps

def reorder_dabg(data,data_headers,array_order,comp_group_list,probeset_db):
    ###array_order gives the final level order sorted, followed by the original index order as a tuple                   
    data2 = {}
    stat_summary_data = {}
    group_name_db = {}
    summary_filtering_stats = {}
    raw_data_comps = {}
    
    stat_result_names = ['avg-','log_fold-','fold-','ttest-']
    group_summary_result_names = ['avg-']
    
    for row_id in data:
        try:
            affygene = probeset_db[row_id][0]
        except TypeError:
            affygene = '' #not needed if not altsplice data
        data_headers2 = {} #reset each time
        grouped_ordered_array_list = {}
        ###Grab each array in the order provided in array_order for each probeset
        for x in array_order:
            y = x[1]  #this is the new first index
            group = x[2]
            group_name = x[3]
            group_name_db[group] = group_name
            #for example y = 5, therefore the data[row_id][5] entry is now the first
            new_item = data[row_id][y]
            ###Add the dabg p-value to a new array for each group
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
            avg1 = statistics.avg(data_list1)
            try: avg2 = statistics.avg(data_list2)
            except ValueError: print data_list2,row_id
            #add std-error later
            log_fold = avg1 - avg2
            fold = statistics.log_fold_conversion(log_fold)
            t,df,tails = statistics.ttest(data_list1,data_list2,2,3) #unpaired student ttest, calls p_value function
            t = abs(t)
            df = round(df) #Excel doesn't recognize fractions in a DF
            ttest = '=tdist('+str(t)+','+str(df)+','+str(tails)+')'
            p = str(statistics.t_probability(t,df))
            comp = group1,group2
            try:
                stat_results[comp] = groups_name,[avg2,log_fold,fold,p],group2_name #this structure is a little weird, since we will use this data for two differnent output files
            except TypeError:
                print comp, len(stat_results), dog
            group_summary_results[group1] = group1_name,[avg1]
            group_summary_results[group2] = group2_name,[avg2]

        stat_result_list = []
        p_list = []
        fold_list = []
        for entry in stat_results:
            data_tuple = entry,stat_results[entry]
            stat_result_list.append(data_tuple)
            ### add in script to grab the smallest p and larges fold for all comparisons
            log_fold = abs(stat_results[entry][1][1])
            p = stat_results[entry][1][3]
            p_list.append(float(p))
            fold_list.append(float(log_fold))
        stat_result_list.sort()
        p_list.sort()
        p_list.sort(); fold_list.sort(); fold_list.reverse()
        summary_filtering_stats[row_id] = p_list[0],fold_list[0]

        grouped_ordered_array_list2 = []
        for group in grouped_ordered_array_list:
            data_tuple = group,grouped_ordered_array_list[group]
            grouped_ordered_array_list2.append(data_tuple)
        grouped_ordered_array_list2.sort() #now the list is sorted by group number
        
        ###for each rowid, add in the reordered data, and new statistics for each group and for each comparison
        for entry in grouped_ordered_array_list2:
            group_number = entry[0]
            original_data_values = entry[1]
            for value in original_data_values:
                try:
                    data2[row_id].append(value)
                except KeyError:
                    data2[row_id] = [value]
            if group_number in group_summary_results:
                group_summary_data = group_summary_results[group_number][1] #the group name is listed as the first entry
                for value in group_summary_data:
                    data2[row_id].append(value)
            for info in stat_result_list:
                if info[0][0] == group_number: #comp,(groups_name,[avg1,log_fold,fold,ttest])
                    comp = info[0]
                    only_add_these = info[1][1][1:] #don't add the first value which is used in the second output file(avg)
                    for value in only_add_these:
                        data2[row_id].append(value)
                    ###Create a second database storing just the statistical information
                    groups_name = info[1][0]
                    stat_summary_data[row_id,comp] = info[1][1] #use this to output multiple files for input into altanalyze later
                    
        ###*******Include a database with the raw values saved for permuteAltAnalyze*******
        for info in comp_group_list:
            temp_raw = []
            group1 = int(info[0])
            group2 = int(info[1])
            comp = str(info[0]),str(info[1])
            g1_data = grouped_ordered_array_list[group1]
            g2_data = grouped_ordered_array_list[group2]
            for value in g2_data:
                temp_raw.append(value)
            for value in g1_data:
                temp_raw.append(value)
            raw_data_comps[row_id,comp] = temp_raw

    ###do the same for the headers, but at the dataset level (redundant processes)
    data_headers4 = []
    data_headers3 = []
    stat_summary_data_names = {}
    try:
        for group in data_headers2:
            data_tuple = group,data_headers2[group]  #e.g. 1, ['X030910_25_hl.CEL', 'X030910_29R_hl.CEL', 'X030910_45_hl.CEL'])
            data_headers3.append(data_tuple)
        data_headers3.sort()
    except UnboundLocalError:
        print data_headers,'\n'
        print array_order,'\n'
        print comp_group_list,'\n'
        kill
    
    for entry in data_headers3:
        x = 0 #indicates the times through a loop
        y = 0 #indicates the times through a loop
        group_number = entry[0]
        original_data_values = entry[1]
        for value in original_data_values:
            data_headers4.append(value)
        if group_number in group_summary_results:
            group_name = group_summary_results[group_number][0]
            group_summary_data = group_summary_results[group_number][1]
            for value in group_summary_data:
                combined_name = group_summary_result_names[x] + group_name  #group_summary_result_names = ['avg-']
                data_headers4.append(combined_name)
                x += 1 #increment the loop index
        temp = []
        for info in stat_result_list:
            #print group_number, info[0]   
            if info[0][0] == group_number:  #comp,(groups_name,[avg1,log_fold,fold,ttest],group2_name)
                groups_name = info[1][0]
                group2_name = info[1][2]
                comp = info[0]
                comp = str(comp[0]),str(comp[1])
                baseline_avg_name = stat_result_names[0] + group2_name
                temp.append(baseline_avg_name) #add avg
                only_add_these = stat_result_names[1:]
                for value in only_add_these:
                    new_name = value + groups_name
                    data_headers4.append(new_name)
                    temp.append(new_name)
                stat_summary_data_names[comp] = temp

    ###For the raw_data only export we need the headers for the different groups (data_headers2), the comparisons (last stat_result_list entry is fine) and group names (group_name_db)       

    raw_data_comp_headers = {}
    for comp in comp_group_list:
        temp_raw = []
        group1 = int(comp[0])
        group2 = int(comp[1])
        comp = str(comp[0]),str(comp[1])
        g1_headers = data_headers2[group1]
        g2_headers = data_headers2[group2]
        g1_name = group_name_db[group1]
        g2_name = group_name_db[group2]
        for header in g2_headers:
            temp_raw.append(g2_name+':'+header)
        for header in g1_headers:
            temp_raw.append(g1_name+':'+header)
        raw_data_comp_headers[comp] = temp_raw
        
    #print data_headers4, stat_summary_data_names
    """for entry in data2:
        print entry
        print data2[entry],dog """

    ###Finished re-ordering lists and adding statistics to data2 and stat_summary_data
    return data2, data_headers4, stat_summary_data, stat_summary_data_names, summary_filtering_stats, raw_data_comp_headers, raw_data_comps
              
    
if __name__ == '__main__':
    print array_cluster_final