import sys,string,os,shutil
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
from scipy import sparse, io
import numpy
import LineageProfilerIterate
import cluster_corr
import export
from import_scripts import ChromiumProcessing
import traceback

""" cellHarmony without alignment """

def cellHarmony(species,platform,query_exp_file,exp_output,
            customMarkers=False,useMulti=False,fl=None,customLabels=None):
    """ Prepare pre-aligned result files in a pre-defined format for cellHarmony post-aligment
    differential and visualization analyses """
    
    customLabels = fl.Labels()  
    reference_exp_file = customMarkers ### pre-formatted from Seurat or other outputs
        
    export_directory = os.path.abspath(os.path.join(query_exp_file, os.pardir))
    if 'ExpressionInput' in query_exp_file:
        ### Change to the root directory above ExpressionINput
        export_directory = os.path.abspath(os.path.join(export_directory, os.pardir))
    dataset_name = string.replace(string.split(query_exp_file,'/')[-1][:-4],'exp.','')
    try: os.mkdir(export_directory+'/cellHarmony/')
    except: pass
    try: os.mkdir(export_directory+'/cellHarmony/CellClassification/')
    except: pass
    try: os.mkdir(export_directory+'/cellHarmony/OtherFiles/')
    except: pass
    
    ### Get the query and reference cells, dataset names
    refererence_cells, query_cells, reference_dataset, query_dataset = importCelltoClusterAnnotations(customLabels)  ### Get the reference and query cells in their respective order
    
    ### copy and re-name the input expression file to the output cellHarmony directory
    if len(reference_dataset)>0 and len(query_dataset)>0:
        target_exp_dir = export_directory+'/cellHarmony/exp.'+reference_dataset+'__'+query_dataset+'-AllCells.txt'
    else:
        target_exp_dir = export_directory+'/cellHarmony/exp.cellHarmony-reference__Query-AllCells.txt'
        reference_dataset = 'cellHarmony-reference'
    shutil.copy(query_exp_file,target_exp_dir) 
    
    ### filter and export the heatmap with just reference cells
    cell_cluster_order = simpleHeaderImport(reference_exp_file)
    filtered_reference_cells=[]
    filtered_query_cells_db={}
    filtered_query_cells=[]
    representative_refcluster_cell = {}
    for cell_id in cell_cluster_order:
        if cell_id in refererence_cells:
            filtered_reference_cells.append(cell_id)
            cluster_label = refererence_cells[cell_id].Label()
            ### Identifies where to place each query cell
            try: representative_refcluster_cell[cluster_label].append(cell_id)
            except: representative_refcluster_cell[cluster_label] = [cell_id]
        elif cell_id in query_cells:
            filtered_query_cells_db[cell_id]=query_cells[cell_id]
            filtered_query_cells.append(cell_id)
    
    #reference_output_file = export.findParentDir(reference_exp_file)+'/'+reference_dataset+'.txt'
    reference_output_file = export_directory+'/cellHarmony/OtherFiles/'+reference_dataset+'.txt'
    reference_output_file2 = export_directory+'/cellHarmony/exp.'+reference_dataset+'__'+query_dataset+'-Reference.txt'
    query_output_file =export_directory+'/'+query_dataset+'.txt'
    ### Write out separate refernece and query files
    from import_scripts import sampleIndexSelection
    sampleIndexSelection.filterFile(reference_exp_file,reference_output_file,['row_clusters-flat']+filtered_reference_cells,force=True)
    sampleIndexSelection.filterFile(target_exp_dir,query_output_file,filtered_query_cells,force=True)
    shutil.copy(reference_output_file,reference_output_file2)
    
    ### export the CellClassification file
    output_classification_file = export_directory+'/cellHarmony/CellClassification/CellClassification.txt'
    exportCellClassifications(output_classification_file,filtered_query_cells_db,filtered_query_cells,representative_refcluster_cell)
    labels_file = export_directory+'/labels.txt'
    exportLabels(labels_file,filtered_reference_cells,refererence_cells)
    fl.setLabels(labels_file)
    
    print 'Files formatted for cellHarmony... running differential expression analyses'
    try:
        print reference_output_file
        print query_output_file
        print output_classification_file
        LineageProfilerIterate.harmonizeClassifiedSamples(species, reference_output_file, query_output_file, output_classification_file,fl=fl)
    except:
        print '\nFAILED TO COMPLETE THE FULL CELLHARMONY ANALYSIS (SEE LOG FILE)...'
        print traceback.format_exc()
    
    return True

def exportCellClassifications(output_file,query_cells,filtered_query_cells,representative_refcluster_cell):
    """ Match the Louvain cellHarmony export format for the classification file """
    
    header = 'Query Barcode\tRef Barcode\tCorrelation\tQuery Partition\tRef Partition\tLabel\n'
    
    o = open(output_file,'w')
    o.write(header)
    for query_barcode in filtered_query_cells:
        CI = query_cells[query_barcode]
        cluster_number = CI.ClusterNumber()
        label = CI.Label()
        ref_barcode = representative_refcluster_cell[label][-1]
        values = [query_barcode,ref_barcode,'1.0',cluster_number,cluster_number,label]
        o.write(string.join(values,'\t')+'\n')
    o.close()

def exportLabels(labels_file,filtered_reference_cells,refererence_cells):
    l = open(labels_file,'w')
    for cell_id in filtered_reference_cells:
        CI = refererence_cells[cell_id]
        cluster_number = CI.ClusterNumber()
        label = CI.Label()
        values = [cell_id,cluster_number,label]
        l.write(string.join(values,'\t')+'\n')
    l.close()
    
def simpleHeaderImport(filename):
    for line in open(filename,'rU').xreadlines():
        data = cleanUpLine(line)
        if '\t' in data:
            t = string.split(data,'\t')
        else:
            t = string.split(data,',')
        header = t[2:]
        header2 = []
        for h in header:
            if ":" in h:
                h = string.split(h,':')[-1]
            header2.append(h)
        break
    return header2

class CellInfo:
    def __init__(self,cell_id, cluster_number, dataset_name, dataset_type, label):
        self.cell_id = cell_id; self.cluster_number = cluster_number; self.dataset_name = dataset_name
        self.dataset_type = dataset_type; self.label = label
    def CellID(self): return self.cell_id
    def ClusterNumber(self): return self.cluster_number
    def DatasetName(self): return self.dataset_name
    def DataSetType(self): return self.dataset_type
    def Label(self): return self.label
    def __repr__(self):
        return self.CellID()+'|'+self.Label()+'|'+self.DataSetType()
    
def importCelltoClusterAnnotations(filename):
    firstRow = True
    refererence_cells={}
    query_cells={}
    for line in open(filename,'rU').xreadlines():
        data = cleanUpLine(line)
        if '\t' in data:
            t = string.split(data,'\t')
        else:
            t = string.split(data,',')
        if firstRow:
            ci = t.index('cell_id')
            cn = t.index('cluster_number')
            try: cm = t.index('cluster_name')
            except: cm = False
            dn = t.index('dataset_name')
            dt = t.index('dataset_type')
            firstRow = False
        else:
            cell_id = t[ci]
            cluster_number = t[cn]
            dataset_name = t[dn]
            dataset_type = t[dt]
            if cm != False:
                cluster_name = t[cm]
                label = cluster_name + '_c'+cluster_number
            else:
                label = 'c'+cluster_number
            
            if string.lower(dataset_type)[0] == 'r':
                dataset_type = 'Reference'
                reference_dataset = dataset_name
                CI = CellInfo(cell_id, cluster_number, dataset_name, dataset_type, label)
                refererence_cells[cell_id]=CI
            else:
                dataset_type = 'Query'
                query_dataset = dataset_name
                CI = CellInfo(cell_id, cluster_number, dataset_name, dataset_type, label)             
                query_cells[cell_id]=CI
    return refererence_cells, query_cells, reference_dataset, query_dataset

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

if __name__ == '__main__':
    platform = 'RNASeq'
    cellHarmony(genome,platform,args.query_h5,None,
            customMarkers=args.reference_h5,useMulti=False,fl=None,customLabels=labels)
    