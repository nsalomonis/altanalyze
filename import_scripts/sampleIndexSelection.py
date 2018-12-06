import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os
import math
import traceback
try: from stats_scripts import statistics
except Exception: pass

def makeTestFile():
    all_data = [['name','harold','bob','frank','sally','kim','jim'],
        ['a','0','0','1','2','0','5'],['b','0','0','1','2','0','5'],
        ['c','0','0','1','2','0','5'],['d','0','0','1','2','0','5']]
    
    input_file = 'test.txt'
    export_object = open(input_file,'w')
    for i in all_data:
        export_object.write(string.join(i,'\t')+'\n')
    export_object.close()
    return input_file

def filterFile(input_file,output_file,filter_names,force=False,calculateCentroids=False,comparisons=[]):
    if calculateCentroids:
        filter_names,group_index_db=filter_names
        
    export_object = open(output_file,'w')
    firstLine = True
    for line in open(input_file,'rU').xreadlines():
        data = cleanUpLine(line)
        if '.csv' in input_file:
            values = string.split(data,',')
        else:
            values = string.split(data,'\t')
        if firstLine:
            uid_index = 0
            if data[0]!='#':
                if force:
                    filter_names2=[]
                    for f in filter_names:
                        if f in values: filter_names2.append(f)
                    filter_names = filter_names2
                try: sample_index_list = map(lambda x: values.index(x), filter_names)
                except:
                    ### If ":" in header name
                    if ':' in line:
                        values2=[]
                        for x in values:
                            if ':' in x:
                                x=string.split(x,':')[1]
                            values2.append(x)
                        values = values2
                        sample_index_list = map(lambda x: values.index(x), filter_names)
                    elif ':' in filter_names:
                        filter_names = map(lambda x: string.split(filter_names,':')[1], filter_names)
                    else:
                        for x in filter_names:
                            if x not in values:
                                print x,
                            print 'are missing';kill
                        
                firstLine = False
                header = values
            if 'PSI_EventAnnotation' in input_file:
                uid_index = values.index('UID')
            if calculateCentroids:
                if len(comparisons)>0:
                    export_object.write(string.join(['UID']+map(lambda x: x[0]+'-fold',comparisons),'\t')+'\n') ### Use the numerator group name                  
                else:
                    clusters = map(str,group_index_db)
                    export_object.write(string.join([values[uid_index]]+clusters,'\t')+'\n')
                continue ### skip the below code
        try: filtered_values = map(lambda x: values[x], sample_index_list) ### simple and fast way to reorganize the samples
        except Exception:
            print 'aaaaaa'
            print traceback.format_exc()
            print len(values), len(sample_index_list)
            print input_file, len(filter_names)
            for i in filter_names:
                if i not in header:
                    print i, 'not found'
            sys.exit()
            
            sys.exit()
            ### For PSI files with missing values at the end of each line, often
            if len(header) != len(values):
                diff = len(header)-len(values)
                values+=diff*['']
            filtered_values = map(lambda x: values[x], sample_index_list) ### simple and fast way to reorganize the samples
            #print values[0]; print sample_index_list; print values; print len(values); print len(prior_values);kill
        prior_values=values
        ######################## Begin Centroid Calculation ########################
        if calculateCentroids:
            mean_matrix=[]
            means={}
            for cluster in group_index_db:
                #### group_index_db[cluster] is all of the indeces for samples in a noted group, cluster is the actual cluster name (not number)
                try: mean=statistics.avg(map(lambda x: float(filtered_values[x]), group_index_db[cluster]))
                except Exception: mean = map(lambda x: filtered_values[uid][x], group_index_db[cluster]) ### Only one value
                means[cluster]=mean
                mean_matrix.append(str(mean))
            filtered_values = mean_matrix
            if len(comparisons)>0:
                fold_matrix=[]
                for (group2, group1) in comparisons:
                    fold = means[group2]-means[group1]
                    fold_matrix.append(str(fold))
                filtered_values = fold_matrix
        ########################  End Centroid Calculation  ######################## 
        export_object.write(string.join([values[uid_index]]+filtered_values,'\t')+'\n')
    export_object.close()
    print 'Filtered columns printed to:',output_file
    return output_file

def filterRows(input_file,output_file,filterDB=None,logData=False,exclude=False):
    export_object = open(output_file,'w')
    firstLine = True
    uid_index=0
    for line in open(input_file,'rU').xreadlines():
        data = cleanUpLine(line)
        values = string.split(data,'\t')
        if 'PSI_EventAnnotation' in input_file and firstLine:
            uid_index = values.index('UID')
        if firstLine:
            try:uid_index=values.index('UID')
            except Exception:
                try: uid_index=values.index('uid')
                except Exception: uid_index = 0
            firstLine = False
            export_object.write(line)
        else:
            if filterDB!=None:
                if values[uid_index] in filterDB:
                    if logData:
                        line = string.join([values[0]]+map(str,(map(lambda x: math.log(float(x)+1,2),values[1:]))),'\t')+'\n'
                    if exclude==False:
                        export_object.write(line)
                elif exclude: ### Only write out the entries NOT in the filter list
                    export_object.write(line)
            else:
                max_val = max(map(float,values[1:]))
                #min_val = min(map(float,values[1:]))
                #if max_val>0.1:
                if max_val < 0.1:
                    export_object.write(line)
    export_object.close()
    print 'Filtered rows printed to:',output_file
    
def getFiltersFromHeatmap(filter_file):
    import collections
    alt_filter_list=None
    group_index_db = collections.OrderedDict()
    index=0
    for line in open(filter_file,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if t[1] == 'row_clusters-flat':
            filter_list = string.split(data,'\t')[2:]
            if ':' in data:
                alt_filter_list = map(lambda x: string.split(x,":")[1],string.split(data,'\t')[2:])
        elif t[0] == 'column_clusters-flat':
            cluster_list = string.split(data,'\t')[2:]
            if 'NA' in cluster_list: ### When MarkerFinder notated groups
                sample_names = map(lambda x: string.split(x,":")[1],filter_list)
                cluster_list = map(lambda x: string.split(x,":")[0],filter_list)
                filter_list = sample_names
            elif alt_filter_list != None: ### When undesired groups notated in the sample names
                filter_list = alt_filter_list
            index=0
            for sample in filter_list:
                cluster=cluster_list[index]
                try: group_index_db[cluster].append(index)
                except Exception: group_index_db[cluster] = [index]
                index+=1
    return filter_list, group_index_db
       
def getComparisons(filter_file):
    """Import group comparisons when calculating fold changes"""
    groups={}
    for line in open(filter_file,'rU').xreadlines():
        data = cleanUpLine(line)
        sample,group_num,group_name = string.split(data,'\t')
        groups[group_num]=group_name
        
    comparisons=[]
    comparison_file = string.replace(filter_file,'groups.','comps.')
    for line in open(comparison_file,'rU').xreadlines():
        data = cleanUpLine(line)
        group2,group1 = string.split(data,'\t')
        group2 = groups[group2]
        group1 = groups[group1]
        comparisons.append([group2,group1])
    return comparisons

def getFilters(filter_file,calculateCentroids=False):
    """Import sample list for filtering and optionally sample to groups """
    filter_list=[]
    if calculateCentroids:
        import collections
        group_index_db = collections.OrderedDict()
    
    index=0
    for line in open(filter_file,'rU').xreadlines():
        data = cleanUpLine(line)
        sample = string.split(data,'\t')[0]
        filter_list.append(sample)
        if calculateCentroids:
            if 'row_clusters-flat' in data:
                forceHeatmapError
            sample,group_num,group_name = string.split(data,'\t')
            try: group_index_db[group_name].append(index)
            except Exception: group_index_db[group_name] = [index]
        index+=1
    if calculateCentroids:
        return filter_list,group_index_db
    else:
        return filter_list
    
""" Filter a dataset based on number of genes with expression above the indicated threshold """

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    #https://stackoverflow.com/questions/36598136/remove-all-hex-characters-from-string-in-python
    data = data.decode('utf8').encode('ascii', errors='ignore') ### get rid of bad quotes
    return data

def statisticallyFilterFile(input_file,output_file,threshold,minGeneCutoff=499):
    if 'exp.' in input_file:
        counts_file = string.replace(input_file,'exp.','geneCount.')
    else:
        counts_file = input_file[:-4]+'-geneCount.txt'
    sample_expressed_genes={}
    header=True
    junction_max=[]
    count_sum_array=[]
    for line in open(input_file,'rU').xreadlines():
        data = cleanUpLine(line)
        if '.csv' in input_file:
            t = string.split(data,',')
        else:
            t = string.split(data,'\t')
        if header:
            header_len = len(t)
            full_header = t
            samples = t[1:]
            header=False
            count_sum_array=[0]*len(samples)
        else:
            if len(t)==(header_len+1):
                ### Correct header with a missing UID column
                samples = full_header
                count_sum_array=[0]*len(samples)
                print 'fixing bad header'
            try: values = map(float,t[1:])
            except Exception:
                if 'NA' in t[1:]:
                    tn = [0 if x=='NA' else x for x in t[1:]] ### Replace NAs
                    values = map(float,tn)
                else:
                    tn = [0 if x=='' else x for x in t[1:]] ### Replace NAs
                    values = map(float,tn)       
                
            binarized_values = []
            for v in values:
                if v>threshold: binarized_values.append(1)
                else: binarized_values.append(0)
            count_sum_array = [sum(value) for value in zip(*[count_sum_array,binarized_values])]
            
    index=0
    distribution=[]
    count_sum_array_db={}
    samples_to_retain =[]
    samples_to_exclude = []
    for sample in samples:
        count_sum_array_db[sample] = count_sum_array[index]
        distribution.append(count_sum_array[index])
        index+=1
    from stats_scripts import statistics
    distribution.sort()
    avg = int(statistics.avg(distribution))
    stdev = int(statistics.stdev(distribution))
    min_exp = int(min(distribution))
    cutoff = avg - (stdev*2)
    dev = 2
    print 'The average number of genes expressed above %s is %s, (SD is %s, min is %s)' % (threshold,avg,stdev,min_exp)
    if cutoff<0:
        if (stdev-avg)>0:
            cutoff = avg - (stdev/2); dev = 0.5
            print cutoff, 'genes expressed selected as a default cutoff to include cells (2-stdev away)'
        else:
            cutoff = avg - stdev; dev = 1
            print cutoff, 'genes expressed selected as a default cutoff to include cells (1-stdev away)'
    if min_exp>cutoff:
        cutoff = avg - stdev; dev = 1
    
    print 'Using a default cutoff of >=500 genes per cell expressed/cell'

    import export
    eo = export.ExportFile(counts_file)
    eo.write('Sample\tGenes Expressed(threshold:'+str(threshold)+')\n')
    for sample in samples: ### keep the original order
        if count_sum_array_db[sample]>minGeneCutoff:
            samples_to_retain.append(sample)
        else:
            samples_to_exclude.append(sample)
        eo.write(sample+'\t'+str(count_sum_array_db[sample])+'\n')

    if len(samples_to_retain)<4: ### Don't remove any if too few samples
        samples_to_retain+=samples_to_exclude
    else:
        print len(samples_to_exclude), 'samples removed (< %d genes expressed)' % minGeneCutoff
    eo.close()
    print 'Exporting the filtered expression file to:'
    print output_file
    filterFile(input_file,output_file,samples_to_retain)

if __name__ == '__main__':
    ################  Comand-line arguments ################

    import getopt
    filter_rows=False
    filter_file=None
    force=False
    exclude = False
    calculateCentroids=False
    geneCountFilter=False
    expressionCutoff=1
    returnComparisons=False
    comparisons=[]
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        filter_names = ['test-1','test-2','test-3']
        input_file = makeTestFile()
        
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','f=','r=','median=','medoid=', 'fold=', 'folds=',
                            'centroid=','force=','minGeneCutoff=','expressionCutoff=','geneCountFilter='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i': input_file=arg
            elif opt == '--f': filter_file=arg
            elif opt == '--median' or opt=='--medoid' or opt=='--centroid': calculateCentroids = True
            elif opt == '--fold': returnComparisons = True
            elif opt == '--r':
                if arg == 'exclude':
                    filter_rows=True
                    exclude=True
                else:
                    filter_rows=True
            elif opt == '--force': force=True
            elif opt == '--geneCountFilter': geneCountFilter=True
            elif opt == '--expressionCutoff': expressionCutoff=float(arg)
            elif opt == '--minGeneCutoff': minGeneCutoff=int(arg)
            
    output_file = input_file[:-4]+'-filtered.txt'
    if geneCountFilter:
        statisticallyFilterFile(input_file,input_file[:-4]+'-OutlierRemoved.txt',1,minGeneCutoff=199); sys.exit()
    if filter_rows:
        filter_names = getFilters(filter_file)
        filterRows(input_file,output_file,filterDB=filter_names,logData=False,exclude=exclude)
    elif calculateCentroids:
        output_file = input_file[:-4]+'-mean.txt'
        if returnComparisons:
            comparisons = getComparisons(filter_file)
            output_file = input_file[:-4]+'-fold.txt'
        try: filter_names,group_index_db = getFilters(filter_file,calculateCentroids=calculateCentroids)
        except Exception:
            print traceback.format_exc()
            filter_names,group_index_db = getFiltersFromHeatmap(filter_file)
        filterFile(input_file,output_file,(filter_names,group_index_db),force=force,calculateCentroids=calculateCentroids,comparisons=comparisons)
    else:
        filter_names = getFilters(filter_file)
        filterFile(input_file,output_file,filter_names,force=force)

