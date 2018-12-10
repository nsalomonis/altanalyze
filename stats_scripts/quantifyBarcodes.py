### quantify Viral Barcodes
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

import sys,string,os,copy
import export
import unique
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

command_args = string.join(sys.argv,' ')
if len(sys.argv[1:])>0 and '--' in command_args: commandLine=True
else: commandLine=False

display_label_names = True

import traceback

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def importViralBarcodeReferences(barcode1,barcode2):
    ### Derive all possible viral barcode combination sequences truncated to 38nt (not ideal but necessary with our read length)
    b1_ls=[]
    spacer='TGGT'
    for line in open(barcode1,'rU').xreadlines():
        b1 = cleanUpLine(line)
        b1_ls.append(b1)
    reference_38mers={}
    for line in open(barcode2,'rU').xreadlines():
        b2 = cleanUpLine(line)
        for b1 in b1_ls:
            #reference_38mers[b1+spacer+b2[:20]]=[]
            reference_38mers[b1+spacer+b2]=[]
    return reference_38mers
        
def processBarcodes(viral_barcode_file,cell_cluster_file,reference_38mers):
    eo = export.ExportFile(viral_barcode_file[:-4]+'-cleaned.txt')
    ### Import a file with the sample names in the groups file in the correct order
    viral_barcodes={}
    repair={}
    short={}
    cluster_header=[]
    
    cell_clusters={}
    for line in open(cell_cluster_file,'rU').xreadlines():
        data = cleanUpLine(line)
        cell, cluster, cluster_name = string.split(data,'\t')
        cell_clusters[cell]=cluster_name
        if cluster_name not in cluster_header:
            cluster_header.append(cluster_name)
        
    cells_with_virus={}
    for line in open(viral_barcode_file,'rU').xreadlines():
        data = cleanUpLine(line)
        cellular, viral = string.split(data,'\t')
        if cellular in cell_clusters:
            try:
                if viral not in cells_with_virus[cellular]:
                    cells_with_virus[cellular].append(viral)
            except Exception: cells_with_virus[cellular]=[viral]
            if len(viral)<48:
            #if len(viral)<38:
                if viral not in repair:
                    repair[viral]=[cellular]
                else:
                    if cellular not in repair[viral]:
                        repair[viral].append(cellular)
            else:
                #short[viral[:35]]=viral
                try:
                    if cellular not in viral_barcodes[viral]:
                        viral_barcodes[viral].append(cellular)
                except Exception: viral_barcodes[viral] = [cellular]

    ### Repair the short sequences
    for viral_short in repair:
        cellular_barcodes = repair[viral_short]
        if viral_short[:35] in short:
            viral = short[viral_short[:35]]
            for cellular in cellular_barcodes:
                try:
                    if cellular not in viral_barcodes[viral]:
                        viral_barcodes[viral].append(cellular)
                except Exception: viral_barcodes[viral] = [cellular]
    print len(viral_barcodes),'unique viral barcodes present'

    #print cells_with_virus['ACGCCGATCTGTTGAG']
    #print cells_with_virus['CAGAATCCAAACTGCT']
    #sys.exit()
    
    valid_barcodes = 0
    for viral in viral_barcodes:
        if viral in reference_38mers:
            valid_barcodes+=1
    print valid_barcodes, 'unique valid viral barcodes present'
    #"""
    ### If the viral barcodes have frequent errors - associate the error with the reference in a cell-specific manner
    ### Only one virus for cell should be present unless it is a doublet
    print len(cells_with_virus), 'cells with viral barcodes'
    doublet_cell={}
    mismatch_to_match={}
    cells_with_valid_barcodes=0
    viral_barcodes_overide={}
    cellular_barcodes_overide={}
    for cellular in cells_with_virus:
        cell_5prime={}
        cell_3prime={}
        ref_sequences=[]
        if len(cells_with_virus[cellular])>1:
            for i in cells_with_virus[cellular]:
                try: cell_5prime[i[:10]].append(i)
                except Exception: cell_5prime[i[:10]]=[i]
                try: cell_3prime[i[-10:]].append(i)
                except Exception: cell_3prime[i[-10:]]=[i]
                if i in reference_38mers:
                    ref_sequences.append(i)
            if len(ref_sequences)>0:
                cells_with_valid_barcodes+=1 ### Determine how many cells have valid viral barcodes
            cell_5prime_ls=[]
            cell_3prime_ls=[]
            for i in cell_5prime:
                cell_5prime_ls.append([len(cell_5prime[i]),i])
            for i in cell_3prime:
                cell_3prime_ls.append([len(cell_3prime[i]),i])
            cell_5prime_ls.sort(); cell_3prime_ls.sort()
            
            for seq in ref_sequences:
                if cell_5prime_ls[-1][1] in seq and cell_3prime_ls[-1][1] in seq:
                    ref_seq = seq
            try: viral_barcodes_overide[ref_seq].append(cellular)
            except: viral_barcodes_overide[ref_seq]=[cellular]
            cellular_barcodes_overide[cellular]=[ref_seq]
            for y in cell_5prime[cell_5prime_ls[-1][1]]:
                mismatch_to_match[y] = ref_seq
            for y in cell_3prime[cell_3prime_ls[-1][1]]:
                mismatch_to_match[y] = ref_seq
                
        else:
            for i in cells_with_virus[cellular]:
                if i in reference_38mers:
                    cells_with_valid_barcodes+=1 ### Determine how many cells have valid viral barcodes
                try: viral_barcodes_overide[i].append(cellular)
                except: viral_barcodes_overide[i]=[cellular]

    viral_barcodes = viral_barcodes_overide
    cells_with_virus = cellular_barcodes_overide
    
    ### Update the viral_barcodes dictionary
    viral_barcodes2={}; cells_with_virus2={}
    for v in viral_barcodes:
        cell_barcodes = viral_barcodes[v]
        proceed = False
        if v in mismatch_to_match:
            v = mismatch_to_match[v]
            proceed = True
        elif v in reference_38mers:
            proceed = True
        if proceed:
            if v in viral_barcodes2:
                for c in cell_barcodes:
                    if c not in viral_barcodes2:
                        viral_barcodes2[v].append(c)
            else:
                viral_barcodes2[v] = cell_barcodes


    
    print cells_with_valid_barcodes, 'cells with valid viral barcodes.'
    viral_barcodes = viral_barcodes2
    ### Update the cells_with_virus dictionary
    for v in viral_barcodes:
        cell_barcodes = viral_barcodes[v]
        for c in cell_barcodes:
            if c in cells_with_virus2:
                if v not in cells_with_virus2[c]:
                    cells_with_virus2[c].append(v)
            else:
                cells_with_virus2[c]=[v]
    cells_with_virus = cells_with_virus2
    
    for c in cells_with_virus:
        if len(cells_with_virus[c])>1:
            doublet_cell[c]=[]
    print len(doublet_cell),'doublets'
    #print cells_with_virus['ACGCCGATCTGTTGAG']
    #print cells_with_virus['CAGAATCCAAACTGCT']
    #sys.exit()
    
    print len(cells_with_virus),'updated cells with virus'
    print len(viral_barcodes),'updated unique viral barcodes'
    #"""
             
    #reference_38mers={}
    
    multi_cell_mapping=0
    unique_cells={}
    multiMappingFinal={}
    import collections
    import unique
    event_db = collections.OrderedDict()
    for cluster in cluster_header:
        event_db[cluster]='0'
    k_value = 1
    import unique
    cluster_hits_counts={}
    cluster_pairs={}
    custom=[]
    for viral in viral_barcodes:
        clusters=[]
        k=len(unique.unique(viral_barcodes[viral]))
        if k>k_value:
            proceed=True
            if len(reference_38mers)>0:
                if viral in reference_38mers:
                    proceed = True
                else: proceed = False
            if proceed:
                viral_cluster_db = copy.deepcopy(event_db) ### copy this
                multi_cell_mapping+=1
                cell_tracker=[]
                multilin=[]
                for cell in viral_barcodes[viral]:
                        #if cell not in doublet_cell:
                        cell_tracker.append(cell)
                        try: unique_cells[cell].append(viral)
                        except: unique_cells[cell] = [viral]
                        if cell in cell_clusters:
                            cluster = cell_clusters[cell]
                            if 'Multi-Lin c4-Mast' in cluster:
                                multilin.append(cell)
                            viral_cluster_db[cluster]='1'
                            clusters.append(cluster)
                c1= unique.unique(clusters)
                #if c1 == ['Multi-Lin c4-Mast']:
                #if c1 == ['MultiLin','MEP','Myelo-1'] or  c1 == ['MultiLin','MEP','Myelo-2'] or  c1 == ['MultiLin','MEP','Myelo-4']:
                #if 'Multi-Lin c4-Mast' in c1 and ('ERP-primed' not in c1 and 'MEP' not in c1 and 'MKP-primed' not in c1 and 'MKP' not in c1 and 'ERP' not in c1) and 'Monocyte' not in c1 and 'e-Mono' not in c1 and ('Gran' in c1 or 'Myelo-1' in c1 or 'Myelo-2' in c1 and 'Myelo-3' in c1 and 'Myelo-4' in c1):
                #if 'Multi-Lin' in c1 and ('e-Mono' in c1 or 'Monocyte' in c1) and ('ERP-primed' in c1 or 'MEP' in c1 or 'MKP-primed' in c1 or 'MKP' in c1) and ('Gran' in c1 or 'Myelo-4' in c1 or 'Myelo-1' in c1 or 'Myelo-2' in c1 or 'Myelo-3' in c1):
                if 'Multi-Lin c4-Mast' in c1:
                    for cell in multilin:
                        print string.join(c1,'|')+'\t'+cell+'\t'+viral
                    custom+=viral_barcodes[viral]
                    
                multiMappingFinal[viral]=viral_cluster_db

        ### Count the number of cluster pairs to make a weighted network
        for c1 in clusters:
            for c2 in clusters:
                if c1 != c2:
                    try:
                        cx = cluster_pairs[c1]
                        try: cx[c2]+=1
                        except: cx[c2]=1
                    except:
                        cx={}
                        cx[c2]=1
                        cluster_pairs[c1] = cx
        clusters = string.join(unique.unique(clusters),'|')
        try: cluster_hits_counts[clusters]+=1
        except Exception: cluster_hits_counts[clusters]=1
    sys.exit()
    print custom
        
    for cluster in cluster_pairs:
        cluster_counts=[]
        cx = cluster_pairs[cluster]
        for c2 in cx:
            count=cx[c2]
            cluster_counts.append([count,c2])
        cluster_counts.sort()
        cluster_counts.reverse()
        #print cluster, cluster_counts
    print len(multiMappingFinal)

    final_ranked_cluster_hits=[]
    for clusters in cluster_hits_counts:
        final_ranked_cluster_hits.append([cluster_hits_counts[clusters],clusters])
    final_ranked_cluster_hits.sort()
    final_ranked_cluster_hits.reverse()
    for (counts,clusters) in final_ranked_cluster_hits:
        print counts, clusters
    
    eo.write(string.join(['UID']+cluster_header,'\t')+'\n')
    for viral_barcode in multiMappingFinal:
        cluster_db = multiMappingFinal[viral_barcode]
        hits=[]
        for cluster in cluster_db:
            hits.append(cluster_db[cluster])
        eo.write(string.join([viral_barcode]+hits,'\t')+'\n')
    eo.close()
    
    eo = export.ExportFile(viral_barcode_file[:-4]+'-cells-'+str(k_value)+'.txt')
    for cell in unique_cells:
        #eo.write(cell+'\t1\t1\t'+str(len(unique_cells[cell]))+'\t'+string.join(unique_cells[cell],'|')+'\n')
        eo.write(cell+'\t1\t1\t\n')
    eo.close()
    #print multi_cell_mapping
    #print len(unique_cells)
    
        

if __name__ == '__main__':
    
    cellBarcodes = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Viral-tracking/3-prime/cellular-viral_und-det.txt'
    cellClusters = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Viral-tracking/3-prime/groups.cellHarmony-Celexa5prime.txt'
    #cellClusters = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Viral-tracking/3-prime/groups.cellHarmony-Celexa5prime-MultiMerge.txt'
    #cellClusters = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Viral-tracking/3-prime/groups.cellHarmony-Celexa5prime-Baso.txt'
    barcode1 = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Viral-tracking/3-prime/14mer.txt'
    barcode2 = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/Viral-tracking/3-prime/30mer.txt'
    references = importViralBarcodeReferences(barcode1,barcode2)
    #references={}
    processBarcodes(cellBarcodes,cellClusters,references);sys.exit()
    import getopt
    filter_rows=False
    filter_file=None
    genome = 'hg19'
    dataset_name = '10X_filtered'
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Insufficient options provided";sys.exit()
        #Filtering samples in a datasets
        #python 10XProcessing.py --i /Users/test/10X/outs/filtered_gene_bc_matrices/ --g hg19 --n My10XExperiment
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','g=','n='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i': matrices_dir=arg
            elif opt == '--g': genome=arg
            elif opt == '--n': dataset_name=arg
     