import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os
import math

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

def filterFile(input_file,output_file,filter_names,force=False):
    export_object = open(output_file,'w')
    firstLine = True
    for line in open(input_file,'rU').xreadlines():
        data = cleanUpLine(line)
        if '.csv' in input_file:
            values = string.split(data,',')
        else:
            values = string.split(data,'\t')
        if firstLine:
            if data[0]!='#':
                if force:
                    filter_names2=[]
                    for f in filter_names:
                        if f in values: filter_names2.append(f)
                    filter_names = filter_names2
                sample_index_list = map(lambda x: values.index(x), filter_names)
                firstLine = False   
                header = values
        try: filtered_values = map(lambda x: values[x], sample_index_list) ### simple and fast way to reorganize the samples
        except Exception:
            ### For PSI files with missing values at the end of each line, often
            if len(header) != len(values):
                diff = len(header)-len(values)
                values+=diff*['']
            filtered_values = map(lambda x: values[x], sample_index_list) ### simple and fast way to reorganize the samples
            #print values[0]; print sample_index_list; print values; print len(values); print len(prior_values);kill
        prior_values=values
        export_object.write(string.join([values[0]]+filtered_values,'\t')+'\n')
    export_object.close()
    print 'Filtered columns printed to:',output_file
    return output_file

def filterRows(input_file,output_file,filterDB=None,logData=False,exclude=False):
    export_object = open(output_file,'w')
    firstLine = True
    for line in open(input_file,'rU').xreadlines():
        data = cleanUpLine(line)
        values = string.split(data,'\t')
        if firstLine:
            firstLine = False
            export_object.write(line)
        else:
            if filterDB!=None:
                if values[0] in filterDB:
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
    
def getFilters(filter_file):
    filter_list=[]
    for line in open(filter_file,'rU').xreadlines():
        data = cleanUpLine(line)
        sample = string.split(data,'\t')[0]
        filter_list.append(sample)
    return filter_list

"""" Filter a dataset based on number of genes with expression above the indicated threshold"""

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def statisticallyFilterFile(input_file,output_file,threshold):
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
            samples = t[1:]
            header=False
            count_sum_array=[0]*len(samples)
        else:
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
    cutoff = 499
    import export
    eo = export.ExportFile(counts_file)
    eo.write('Sample\tGenes Expressed(threshold:'+str(threshold)+')\n')
    for sample in samples: ### keep the original order
        if count_sum_array_db[sample]>cutoff:
            samples_to_retain.append(sample)
        else:
            samples_to_exclude.append(sample)
        eo.write(sample+'\t'+str(count_sum_array_db[sample])+'\n')

    if len(samples_to_retain)<4: ### Don't remove any if too few samples
        samples_to_retain+=samples_to_exclude
    else:
        print len(samples_to_exclude), 'samples removed (< 500 genes expressed)' # (%s)' % (dev,string.join(samples_to_exclude,', '))
    eo.close()
    print 'Exporting the filtered expression file to:'
    print output_file
    filterFile(input_file,output_file,samples_to_retain)

def combineDropSeq(input_dir):
    import unique
    files = unique.read_directory(input_dir)
    combinedGeneExpression={}
    for input_file in files: #:70895507-70895600
        header=True
        if '.txt' in input_file:
            for line in open(input_dir+'/'+input_file,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if header:
                    header_row = line
                    samples = t[1:]
                    header=False
                else:
                    values = map(float,t[1:])
                    gene = t[0]
                    if gene in combinedGeneExpression:
                        prior_values = combinedGeneExpression[gene]
                        count_sum_array = [sum(value) for value in zip(*[prior_values,values])]
                    else:
                        count_sum_array = values
                    combinedGeneExpression[gene] = count_sum_array

    input_file = input_dir+'/test.txt'
    export_object = open(input_file,'w')
    export_object.write(string.join(['UID']+samples,'\t')+'\n')
    for gene in combinedGeneExpression:
        values = string.join(map(str,[gene]+combinedGeneExpression[gene]),'\t')
        export_object.write(values+'\n')
    export_object.close()

if __name__ == '__main__':
    ################  Comand-line arguments ################
    #statisticallyFilterFile('/Volumes/salomonis2/Driscoll/DriscollPREPOST/run1964_PRE_normalized.txt','/Volumes/salomonis2/Driscoll/DriscollPREPOST/run1964_PRE_normalized2.txt',1); sys.exit()
    import getopt
    filter_rows=False
    filter_file=None
    force=False
    exclude = False
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        filter_names = ['bob','sally','jim']
        input_file = makeTestFile()
        
        #Filtering samples in a datasets
        #python SampleSelect.py --i /Users/saljh8/Desktop/C4-hESC/ExpressionInput/exp.C4.txt --f /Users/saljh8/Desktop/C4-hESC/ExpressionInput/groups.C4.txt
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','f=','r=','force='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i': input_file=arg
            elif opt == '--f': filter_file=arg
            elif opt == '--r':
                if arg == 'exclude':
                    filter_rows=True
                    exclude=True
                else:
                    filter_rows=True
            elif opt == '--force': force=True
            
    output_file = input_file[:-4]+'-filtered.txt'
    if filter_rows:
        filter_names = getFilters(filter_file)
        filterRows(input_file,output_file,filterDB=filter_names,logData=False,exclude=exclude)
    elif filter_file ==None:
        combineDropSeq(input_file)
    else:
        filter_names = getFilters(filter_file)
        filterFile(input_file,output_file,filter_names,force=force)

