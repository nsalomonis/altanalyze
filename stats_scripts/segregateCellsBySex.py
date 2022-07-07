import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import export
import unique
import traceback
import numpy as np

""" Determine and segregate single-cell based on sex associated gene expression """

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def computeOnSexTranscripts(fn):
    """ Import a flat single-cell expression matrix """
    female = ['TSIX','XIST','ENSG00000229807','ENSMUSG00000086503','ENSG00000270641','ENSMUSG00000085715']
    male = ['EIF2S3Y','DDX3Y','UTY','SRY','ENSMUSG00000069049','ENSMUSG00000069049', 'ENSG00000067048', 'ENSMUSG00000069045','ENSMUSP00000070012', 'ENSG00000183878','ENSG00000184895','ENSG00000067048']
    informative_genes = male+female
    female_data = []
    male_data = []
    firstRow=True
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            header = t[1:]
            firstRow = False
        else:
            gene = t[0]
            gene_upper = string.upper(gene)
            if gene_upper in informative_genes:
                values = map(float,t[1:])
                if gene_upper in male:
                    male_data.append(values)
                else:
                    female_data.append(values)
    print len(male_data)
    print len(female_data)
    male_data = np.array(male_data)
    female_data = np.array(female_data)
    male_data = np.average(male_data, axis=0)
    female_data = np.average(female_data, axis=0)
    ratios = [i - j for i, j in zip(female_data, male_data)]
    
    female_count=0
    male_count=0
    uknown=0
    i=0
    females=[]
    males=[]
    unknowns=[]
    for ratio in ratios:
        if ratio>0:
            female_count+=1
            females.append(header[i])
        elif ratio<0:
            male_count+=1
            males.append(header[i])
        else:
            uknown+=1
            unknowns.append(header[i])
        i+=1
    
    male_ratio = 100*(male_count*1.00)/((male_count*1.00)+(female_count*1.00)+uknown*1.00)
    female_ratio = 100*(female_count*1.00)/((male_count*1.00)+(female_count*1.00)+uknown*1.00)
    uknown_ratio = 100*(uknown*1.00)/((male_count*1.00)+(female_count*1.00)+uknown*1.00)

    print 'Male cell percentage:',male_ratio
    print 'Female cell percentage:',female_ratio
    print 'Unknown cell percentage:',uknown_ratio
    
    from import_scripts import sampleIndexSelection
    sampleIndexSelection.filterFile(fn,fn[:-4]+'-female.txt',['row_clusters-flat']+females,force=True)
    sampleIndexSelection.filterFile(fn,fn[:-4]+'-male.txt',['row_clusters-flat']+males,force=True)
    
    export_object = open(fn[:-4]+'-groups.txt','w')
    def writeData(header,group):
        export_object.write(header+'\t'+group+'\n')
    map(lambda x: writeData(x,'male'),males)
    map(lambda x: writeData(x,'female'),females)
    map(lambda x: writeData(x,'unknown'),unknowns)
    export_object.close()
                        
if __name__ == '__main__':
    ################  Comand-line arguments ################
    import getopt

    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print 'WARNING!!!! Too commands supplied.'
        
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i':
                inputFile = arg
                
    computeOnSexTranscripts(inputFile)

