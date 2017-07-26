###GO_parsing
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

import math
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os.path
import unique

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv":
            dir_list2.append(entry)
    return dir_list2

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def parse_affymetrix_annotations(filename):
    temp_affy_db = {}
    x=0
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():             
        probeset_data,null = string.split(line,'\n')  #remove endline
        affy_data = string.split(probeset_data[1:-1],'","')  #remove endline
        if x==0:
            if probeset_data[0] == '#':
                continue
            x +=1
            affy_headers = affy_data
        else:
            x +=1
            probesets = affy_data[0]
            temp_affy_db[probesets] = affy_data[1:]
    for header in affy_headers:
        x = 0; eg = ''; gs = ''
        while x < len(affy_headers):
            if 'rocess' in affy_headers[x]: gb = x - 1
            if 'omponent' in affy_headers[x]: gc = x - 1
            if 'olecular' in affy_headers[x]: gm = x - 1
            if 'athway' in affy_headers[x]: gp = x - 1
            if 'Gene Symbol' in affy_headers[x]: gs = x - 1
            if 'Ensembl' in affy_headers[x]: eg = x - 1
            x += 1
        ###Below code used if human exon array parsed
        global analyze_human_exon_data
        analyze_human_exon_data = 'no'
        if eg == '':
            x = 0
            while x < len(affy_headers):
                if 'mrna_assignment' in affy_headers[x]:
                    eg = x - 1
                    analyze_human_exon_data = 'yes'
                x+=1
    for probeset in temp_affy_db:
        affy_data = temp_affy_db[probeset]
        try:
            go_bio = affy_data[gb]
        except IndexError:
            ###Occurs due to a new line error
            continue
        go_com = affy_data[gc]
        go_mol = affy_data[gm]
        genmapp = affy_data[gp]
        if gs == '': symbol = ''
        else: symbol = affy_data[gs]
        if analyze_human_exon_data == 'no':
            ensembl = affy_data[eg]
        else:
            ensembl_data = affy_data[eg]
            ensembl=''
            try:
                if 'gene:ENSMUSG' in ensembl_data:
                    ensembl_data = string.split(ensembl_data,'gene:ENSMUSG')
                    ensembl_data = string.split(ensembl_data[1],' ')
                    ensembl = 'ENSMUSG'+ ensembl_data[0]
                if 'gene:ENSG' in ensembl_data:
                    ensembl_data = string.split(ensembl_data,'gene:ENSG')
                    ensembl_data = string.split(ensembl_data[1],' ')
                    ensembl = 'ENSG'+ ensembl_data[0]
            except IndexError:
                continue
        goa=[]
      
        goa = merge_go_annoations(go_bio,goa)
        goa = merge_go_annoations(go_com,goa)
        goa = merge_go_annoations(go_mol,goa)
        goa = merge_go_annoations(genmapp,goa)

        goa=unique.unique(goa); goa.sort(); 
        goa = string.join(goa,'')
        try:
            ensembl = string.split(ensembl,' /// ')
        except ValueError:
            ensembl = [ensembl]
        for ensembl_id in ensembl:
            if len(goa)>10:
                go_annotations[ensembl_id] = goa, symbol

def merge_go_annoations(go_category,goa):
    dd = ' // '
    td = ' /// '
    if td in go_category:
        go_split = string.split(go_category,td)
        for entry in go_split:
          if analyze_human_exon_data == 'no':
            try:
                a,b,null = string.split(entry,dd )
                entry = b+dd
                goa.append(entry)
            except ValueError:
                ###occurs with GenMAPP entries
                a,null = string.split(entry,dd )
                entry = a+dd
                goa.append(entry)
          else:
            try:
                f,a,b,null = string.split(entry,dd )
                entry = b+dd
                goa.append(entry)
            except ValueError:
                ###occurs with GenMAPP entries
                f,null,a = string.split(entry,dd )
                entry = a+dd
                goa.append(entry)
    else:
        if go_category != '---':
            if dd in go_category:
              if analyze_human_exon_data == 'no':
                try:
                    a,null = string.split(go_category,dd)
                    entry = a+dd
                except ValueError:
                    a,b,null = string.split(go_category,dd )
                    entry = b+dd
                goa.append(entry)
              else:
                try:
                    f,null,a = string.split(go_category,dd)
                    entry = a+dd
                except ValueError:
                    f,a,b,null = string.split(go_category,dd )
                    entry = b+dd
                goa.append(entry)        
            else:
                goa.append(go_category)
    return goa
            
def getDirectoryFiles(dir,species):
    dir_list = read_directory('/AltDatabase/'+species+'/exon')  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = dir[1:]+'/'+data
        if 'MoEx-1_0-st-transcript-annot.csv' in affy_data_dir and species == 'Mm':
            transcript_annotation_file = affy_data_dir
        elif 'HuEx-1_0-st-transcript-annot.csv' in affy_data_dir and species == 'Hs':
            transcript_annotation_file = affy_data_dir
        elif 'annot.hg' in affy_data_dir:
            probeset_transcript_file = affy_data_dir
    return transcript_annotation_file

def parseAffyGO(use_exon_data,get_splicing_factors,species):
    print "Parsing out Affymetrix GO annotations"
    global go_annotations; go_annotations={}
    import_dir = '/AltDatabase/affymetrix/'+species+'/'
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for affy_data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = import_dir[1:]+affy_data
        if use_exon_data == 'yes':
            affy_data_dir = getDirectoryFiles('/AltDatabase/'+species+'/exon',species)
        try: parse_affymetrix_annotations(affy_data_dir)
        except Exception: null=[]
    print len(go_annotations),"Ensembl gene annotations parsed"
    if get_splicing_factors == 'yes':
        go = go_annotations; mRNA_processing_ensembl=[]
        for entry in go_annotations:
            if 'splicing' in go[entry][0] or ('mRNA' in go[entry][0] and 'processing' in go[entry][0]):
                mRNA_processing_ensembl.append(entry)
        print len(mRNA_processing_ensembl),"mRNA processing/splicing regulators identified"
        return mRNA_processing_ensembl
    else:
        return go_annotations

if __name__ == '__main__':
    go = parseAffyGO('yes','no','Hs')