###SubGeneViewerExport
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
import os.path
import unique
import export
dirfile = unique

############ File Import Functions #############
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    #add in code to prevent folder names from being included
    dir_list2 = [] 
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".all" or entry[-5:] == ".data" or entry[-3:] == ".fa":
            dir_list2.append(entry)
    return dir_list2

def returnDirectories(sub_dir):
    dir=os.path.dirname(dirfile.__file__)
    dir_list = os.listdir(dir + sub_dir)
    ###Below code used to prevent FILE names from being included
    dir_list2 = []
    for entry in dir_list:
        if "." not in entry: dir_list2.append(entry)
    return dir_list2

class GrabFiles:
    def setdirectory(self,value): self.data = value
    def display(self): print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        files = getDirectoryFiles(self.data,search_term)
        if len(files)<1: print 'files not found'
        return files
    def returndirectory(self):
        dir_list = getAllDirectoryFiles(self.data)
        return dir_list

def getAllDirectoryFiles(import_dir):
    all_files = []
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        data_dir = import_dir[1:]+'/'+data
        all_files.append(data_dir)
    return all_files

def getDirectoryFiles(import_dir,search_term):
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    matches=[]
    for data in dir_list:    #loop through each file in the directory to output results
        data_dir = import_dir[1:]+'/'+data
        if search_term in data_dir: matches.append(data_dir)
    return matches

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

############### Main Program ###############

def importAnnotationData(filename):
    fn=filepath(filename); x=1
    global gene_symbol_db; gene_symbol_db={}
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            gene = t[0]
            try: symbol = t[1]
            except IndexError: symbol = ''
            if len(symbol)>0: gene_symbol_db[gene] = symbol
            
def importGeneData(filename,data_type):
    fn=filepath(filename); x=0; gene_db={}
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:x=1
        else:
            proceed = 'yes'
            if data_type == 'junction': gene, region5, region3 = t; value_str = region5+':'+region3
            if data_type == 'feature':
                probeset, gene, feature, region = t; value_str = region,feature+':'+region+':'+probeset ###access region data later
                #if (gene,region) not in region_db: region_db[gene,region] = feature,probeset  ### Needed for processed structure table (see two lines down)
                try: region_db[gene,region].append((feature,probeset))  ### Needed for processed structure table (see two lines down)
                except KeyError: region_db[gene,region] = [(feature,probeset)]
                try: region_count_db[(gene,region)]+=1
                except KeyError: region_count_db[(gene,region)]=1
                ###have to add in when parsing structure probeset values for nulls (equal to 0)
            if data_type == 'structure':
                gene, exon, type, block, region, const, start, annot = t; region_id = exon
                if len(annot)<1: annot = '---'
                if (gene,exon) in region_db:
                    probeset_data = region_db[(gene,exon)]
                    for (feature,probeset) in probeset_data:
                        count = str(region_count_db[(gene,exon)]) ###here, feature is the label (reversed below)
                        value_str = feature+':'+exon+':'+probeset+':'+type+':'+count+':'+const+':'+start+':'+annot
                        if gene in gene_symbol_db: ###Only incorporate gene data with a gene symbol, since Cytoscape currently requires this
                            try: gene_db[gene].append(value_str)
                            except KeyError: gene_db[gene] = [value_str]
                            proceed = 'no'
                else: ### Occurs when no probeset is present: E.g. the imaginary first and last UTR region if doesn't exit
                    feature = exon ###feature contains the region information, exon is the label used in Cytoscape
                    exon,null = string.split(exon,'.')
                    probeset = '0'
                    count = '1'
                    null_value_str = exon,exon+':'+feature+':'+probeset ###This is how Alex has it... to display the label without the '.1' first
                    try: feature_db[gene].append(null_value_str)
                    except KeyError: feature_db[gene] = [null_value_str]
                    value_str = exon+':'+feature+':'+probeset+':'+type+':'+count+':'+const+':'+start+':'+annot
                if gene in structure_region_db:
                    order_db = structure_region_db[gene]
                    order_db[exon] = block
                else:
                    order_db = {}
                    order_db[exon] = block
                    structure_region_db[gene] = order_db
            if gene in gene_symbol_db and proceed == 'yes': ###Only incorporate gene data with a gene symbol, since Cytoscape currently requires this
                try: gene_db[gene].append(value_str)
                except KeyError: gene_db[gene] = [value_str]
    return gene_db
            
def exportData(gene_db,data_type,species):
    export_file = 'AltDatabase/ensembl/SubGeneViewer/'+species+'/Xport_sgv_'+data_type+'.csv'
    if data_type == 'feature': title = 'gene'+'\t'+'symbol'+'\t'+'sgv_feature'+'\n'
    if data_type == 'structure': title = 'gene'+'\t'+'symbol'+'\t'+'sgv_structure'+'\n'
    if data_type == 'splice': title = 'gene'+'\t'+'symbol'+'\t'+'sgv_splice'+'\n'
    data = export.createExportFile(export_file,'AltDatabase/ensembl/SubGeneViewer/'+species)
    #fn=filepath(export_file); data = open(fn,'w')
    data.write(title)
    for gene in gene_db:
        try:
            symbol = gene_symbol_db[gene]
            value_str_list = gene_db[gene]
            value_str = string.join(value_str_list,',')
            values = string.join([gene,symbol,value_str],'\t')+'\n'; data.write(values)
        except KeyError: null = []
    data.close()
    print "exported to",export_file

def customLSDeepCopy(ls):
    ls2=[]
    for i in ls: ls2.append(i)
    return ls2

def reorganizeData(species):
    global region_db; global region_count_db; global structure_region_db; global feature_db
    region_db={}; region_count_db={}; structure_region_db={}
    import_dir = '/AltDatabase/ensembl/'+species
    g = GrabFiles(); g.setdirectory(import_dir)
    exon_struct_file = g.searchdirectory('exon-structure')
    feature_file = g.searchdirectory('feature-data')
    junction_file = g.searchdirectory('junction-data')
    annot_file = g.searchdirectory('Ensembl-annotations.')

    importAnnotationData(annot_file[0])
    ### Run the files through the same function which has options for different pieces of data. Feature data is processed a bit differently
    ### since fake probeset data is supplied for intron and UTR features not probed for
    splice_db = importGeneData(junction_file[0],'junction')
    feature_db = importGeneData(feature_file[0],'feature')
    structure_db = importGeneData(exon_struct_file[0],'structure')
    
    for gene in feature_db:
        order_db = structure_region_db[gene]
        temp_list0 = []; temp_list = []; rank = 1
        for (region,value_str) in feature_db[gene]:
            ###First, we have to get the existing order... this is important because when we sort, it screw up ranking within an intron with many probesets
            temp_list0.append((rank,region,value_str)); rank+=1
        for (rank,region,value_str) in temp_list0:
            try: block_number = order_db[region]
            except KeyError: print gene, region, order_db;kill
            temp_list.append((int(block_number),rank,value_str)) ###Combine the original ranking plus the ranking included from taking into account regions not covered by probesets
        temp_list.sort()
        temp_list2 = []
        for (block,rank,value_str) in temp_list:
            temp_list2.append(value_str)
        feature_db[gene] = temp_list2
        
    exportData(splice_db,'splice',species)
    exportData(structure_db,'structure',species)
    exportData(feature_db,'feature',species)
    
if __name__ == '__main__':
    dirfile = unique
    species = 'Hs'
    reorganizeData(species)
