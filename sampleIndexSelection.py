import sys,string
import os

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

def filterFile(input_file,output_file,filter_names):
    export_object = open(output_file,'w')
    firstLine = True
    for line in open(input_file,'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstLine:
            if data[0]!='#':
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

def filterRows(input_file,output_file,filterDB=None):
    export_object = open(output_file,'w')
    firstLine = True
    for line in open(input_file,'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstLine:
            firstLine = False
            export_object.write(line)
        else:
            if filterDB!=None:
                if values[0] in filterDB:
                    export_object.write(line)
            else:
                max_val = max(map(float,values[1:]))
                #min_val = min(map(float,values[1:]))
                #if max_val>5:
                if max_val < 0.1:
                    export_object.write(line)
    export_object.close()
    print 'Filtered rows printed to:',output_file
    
def getFilters(filter_file):
    filter_list=[]
    for line in open(filter_file,'rU').xreadlines():
        data = line.rstrip()
        sample = string.split(data,'\t')[0]
        filter_list.append(sample)
    return filter_list

if __name__ == '__main__':
    ################  Comand-line arguments ################
    import getopt
    filter_rows=False
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        filter_names = ['bob','sally','jim']
        input_file = makeTestFile()
        
        #Filtering samples in a datasets
        #python SampleSelect.py --i /Users/saljh8/Desktop/C4-hESC/ExpressionInput/exp.C4.txt --f /Users/saljh8/Desktop/C4-hESC/ExpressionInput/groups.C4.txt
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','f=','r='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i': input_file=arg
            elif opt == '--f': filter_file=arg
            elif opt == '--r': filter_rows=True
            
    output_file = input_file[:-4]+'-output.txt'
    if filter_rows:
        filterRows(input_file,output_file)
    else:
        filter_names = getFilters(filter_file)
        filterFile(input_file,output_file,filter_names)

