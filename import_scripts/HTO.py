### ADT.py
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

import math
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import export
import copy
from import_scripts import ChromiumProcessing
from import_scripts import sampleIndexSelection
from import_scripts import CountsNormalize

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def HTO(filename,ratio=0.6,DataType='log'):
    firstRow = True
    hashed={}
    export_results = filename[:-4]+'-count-assigned-singlets.txt'
    eos = export.ExportFile(export_results)
    multiplet_results = filename[:-4]+'-count-assigned-multiplets.txt'
    eo = export.ExportFile(multiplet_results)	
    count = 0; singlets = 0; multiplets = 0
    print 'ratio:',ratio
    multiplet_list=[]
    for line in open(filename,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t'); count+=1
        if firstRow:
            firstRow = False
            header = t[1:]
            count_sum_array=[0]*len(header)
        else:
            cell = t[0]
            values = map(float,t[1:])
            if DataType == 'log':
                values = map(lambda x: math.pow(2,x),values)
            sum_values = sum(values)
            index=0
            singlet=None
            multiplet=[]
            for val in values:
                if val>0:
                    if (val/sum_values)>0.33:
                        multiplet.append(header[index])
                    if (val/sum_values)>ratio:
                        singlet = header[index]
                index+=1
            if len(multiplet)>1:
                eo.write(cell+'\t'+string.join(multiplet,'|')+'\n'); multiplets+=1
                multiplet_list.append(cell)
            elif singlet!=None:
                eos.write(cell+'\t'+singlet+'\n'); singlets+=1
            
            count_sum_array = [sum(value) for value in zip(*[count_sum_array,values])]
    eos.close()
    eo.close()
    print singlets,'singlets with one cell assignment out of',count
    print multiplets,'multiplets out of',count


    def calculateCPTT(val,barcode_sum):
        if val==0:
            return 0
        else:
            return 10000.00*val/barcode_sum

    export_results = filename[:-4]+'-norm-assigned-singlets.txt'
    eos = export.ExportFile(export_results)
    multiplet_results = filename[:-4]+'-norm-assigned-multiplets.txt'
    eo = export.ExportFile(multiplet_results)	
    count = 0; singlets = 0; multiplets = 0
    firstRow=True
    for line in open(filename,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t'); count+=1
        if firstRow:
            firstRow = False
            barcodes = t[1:]
        else:
            cell = t[0]
            values = map(float,t[1:])
            if DataType == 'log':
                values = map(lambda x: math.pow(2,x),values)
            sum_values = sum(values)
            index=0
            singlet=None
            multiplet=[]
            
            index2=0
            cptt_values = []
            for barcode in barcodes:
                try: barcode_sum = count_sum_array[index2]
                except:
                    print [gene], len(count_sum_array), len(values),len(barcodes);sys.exit()
                val = values[index2]
                cptt_val = calculateCPTT(val,barcode_sum)
                cptt_values.append(cptt_val)
                index2+=1
                
            sum_values = sum(cptt_values)
            for val in cptt_values:
                if val>0:
                    if (val/sum_values)>0.3:
                        multiplet.append(header[index])
                    if (val/sum_values)>ratio:
                        singlet = header[index]
                index+=1
            
            if len(multiplet)>1 or cell in multiplet_list:
                eo.write(cell+'\t'+string.join(multiplet,'|')+'\n'); multiplets+=1
            elif singlet!=None:
                eos.write(cell+'\t'+singlet+'\n'); singlets+=1
    eos.close()
    eo.close()
    print singlets,'singlets with one cell assignment out of',count
    print multiplets,'multiplets out of',count

if __name__ == '__main__':    
	################  Comand-line arguments ################
	import getopt

	ratio = 0.6
	counts_path = None
	h5 = None
	#python stats_scripts/ADT.py --i /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/CPTT/exp.Biolegend-ADT-1.txt --a /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/feature_reference.txt --g /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Biolegend/CPTT/AltAnalyze/cellHarmony/cellHarmony/QueryGroups.cellHarmony_captures-filtered.txt --s /Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Isolation-Strategy/HSCP.txt
	if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
		print 'WARNING!!!! Too commands supplied.'
		
	else:
		options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','HTO=','ratio='])
		#print sys.argv[1:]
		for opt, arg in options:
			if opt == '--i':
				h5 = arg
			if opt == '--HTO':
				counts_path = arg
			if opt == '--ratio': ###
				ratio = float(arg)
	
	filter_file=None
	genome = 'hg19'
	dataset_name = '10X_filtered'
	if counts_path == None:
		cptt_path = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Mm-100k-CITE-Seq/Kyle-Sorted/KPORT-HTO/KFport2_matrix_CPTT.txt'
		cptt_path = ChromiumProcessing.import10XSparseMatrix(h5,genome,dataset_name,geneIDs = False, minReads=1000, maxCells=25000)
		counts_path = string.replace(cptt_path,'_CPTT','')
		output_file = counts_path[:-4]+'-HTO.txt'
		sampleIndexSelection.filterRows(counts_path,output_file,logData=False,stringMatch='-ADT')
	else:
		output_file = counts_path
	input_file = sampleIndexSelection.transposeMatrix(output_file)
	#input_file = CountsNormalize.normalizeDropSeqCountsMemoryEfficient(input_file)
	HTO(input_file,ratio=ratio,DataType='non-log')
	