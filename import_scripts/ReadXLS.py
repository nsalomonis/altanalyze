import sys, string
import os.path
import unique
import copy
import time
import math
import export
from xlrd import open_workbook
import traceback

sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

""" Methods for reading Excel files """

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def makeUnique(item):
    db1={}; list1=[]; k=0
    for i in item:
        try: db1[i]=[]
        except TypeError: db1[tuple(i)]=[]; k=1
    for i in db1:
        if k==0: list1.append(i)
        else: list1.append(list(i))
    list1.sort()
    return list1

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def importExcelFile(filename, column=False, top=None):
    """ Iterate through each tab to get relevant data """
    from xlrd import open_workbook
    worksheet_data={}
    results = False
    wb = open_workbook(filename)
    for s in wb.sheets():
        worksheet = s.name
        rows=[]
        count=0
        for row in range(s.nrows):
            count+=1
            values = []
            if column!=False: ### Return results from a single designated column number
                if top!=None:
                    if count>top:
                        continue
                try: values.append(str(s.cell(row,column-1).value))
                except Exception: pass
            else:
                for col in range(s.ncols):
                    if top!=None:
                        if count>top:
                            continue
                    try: values.append(str(s.cell(row,col).value))
                    except Exception: pass
            rows.append(values)
        worksheet_data[worksheet]=rows
    return worksheet_data

def processWorksheetMarkers(worksheet_data):
    """ Compile the ID to category results """""
    compiled_results=[]
    for worksheet in worksheet_data:
        firstRow=True
        for value in worksheet_data[worksheet]:
            if firstRow:
                category = value
                firstRow = False
            else:
                compiled_results.append([category,value])
    return compiled_results

def exportCompiledResults(output,compiled_results):
    export_object = open(output,'w')
    for (category,value) in compiled_results:
        try: export_object.write(category[0]+'\t'+value[0]+'\n')
        except: pass
    export_object.close()
     
def computeSignEnrichmentScore(reference,query,useMaxLength=False,randomize=False):
    """ Given the total number of features in the reference signature and
    the total in the compared query (e.g., folds, dPSI values), compare the
    sign of the feature in the reference and query to compute an enrichment
    score """
    ref_len = len(reference)
    query_len = len(query)
    downScore = 0
    upScore = 0
    for feature in query:
        if feature in reference:
            ref_sign = reference[feature]            
            if randomize:
                signs = ['-','+']
                query_sign = random.shuffle(signs)[0]
            else:
                query_sign = query[feature]
            if ref_sign == '-':
                if query_sign == '-':
                    downScore+=1
                else:
                    downScore-=1
            elif ref_sign == '+':
                if query_sign == '+':
                    upScore+=1
                else:
                    upScore-=1
                    
    if useMaxLength:
        ref_len = max(ref_len,query_len)
    upScore = upScore/(ref_len*1.000)
    downScore = downScore/(ref_len*1.000)
    Score1 = (upScore+downScore)*0.5
    
    if Score1>0: sign = 1
    else: sign = -1
    
    Score1T = sign* math.log(1+(100*abs(Score1)),10)
    return Score1T

def computeSignPval(score_db,random_db):
    from scipy.stats import norm
    stats_db = {}
    for (ref,query) in score_db:
        Score1T = score_db[(ref,query)]
        RScore_stdev = numpy.std(random_db[(ref,query)])
        z = (Score1T)/RScore_stdev
        p =  2* (1-norm.cdf(z))
        stats_db[ref,query] = z,p
        
    return stats_db
    
if __name__ == '__main__':
    import getopt
    top = None
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['xls=','i=','column=','output=','o=','top='])
        for opt, arg in options:
            if opt == '--xls' or opt == '--i': xls=arg
            if opt == '--output' or opt == '--o': output=arg
            if opt == '--top': top=int(arg)
            if opt == '--column':
                try: column = int(arg)
                except: column = arg
    
    worksheet_data = importExcelFile(xls,column=column,top=top)
    compiled_results = processWorksheetMarkers(worksheet_data)
    exportCompiledResults(output,compiled_results)
    