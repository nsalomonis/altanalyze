### WikiPathways_webservice.py
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

import base64
import export
import time
import unique
import traceback
try:
    import gene_associations
    """
    import suds
    from suds.client import Client
    wsdl = 'http://www.wikipathways.org/wpi/webservice/webservice.php?wsdl'
    client = Client(wsdl) """
    from wikipathways_api_client import WikipathwaysApiClient
    wikipathways_api_client_instance = WikipathwaysApiClient()
except Exception:
    #print traceback.format_exc()
    None ### Occurs when functions in this module are resued in different modules
            
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

class PathwayData:
    def __init__(self,wpname):
        self._wpname = wpname
    def WPName(self): return self._wpname
    def setSourceIDs(self,id_db): self.source_ids = id_db
    def SourceIDs(self): return self.source_ids
    def Report(self):
        output = self.WPName()
        return output
    def __repr__(self): return self.Report()

def getPathwayAs(pathway_db,species_code,mod):
    begin_time = time.time()
    for wpid in pathway_db:
        #print [wpid],'pathway_db',len(pathway_db)
        file_type = 'gpml'
        #file_type = 'application/gpml+xml'
        processor_time = str(time.clock())
        #try: export.deleteFolder('BuildDBs/WPs') ### clear any remaining pathway files
        #except Exception: pass
        #wp_id_data = client.service.getPathwayAs(fileType = file_type,pwId = wpid, revision = 0)
        kwargs = {
            'identifier': 'WP2062',
            'version': 0,
            'file_format': 'application/gpml+xml'}
        #wp_id_data = wikipathways_api_client_instance.get_pathway_as(**kwargs)
        wp_id_data = wikipathways_api_client_instance.get_pathway_as(file_format = file_type,identifier = wpid, version = 0)
        #wp_id_data = base64.b64decode(wp_id_data)
        gpml_path = filepath('BuildDBs/WPs/'+processor_time+'/'+wpid+'.gpml')
        #print gpml_path
        outfile = export.ExportFile(gpml_path)
        outfile.write(wp_id_data); outfile.close()
        gene_system_list = string.split(wp_id_data,'\n')
        parent_path = export.findParentDir(gpml_path)
        pathway_db = gene_associations.getGPMLGraphData(parent_path,species_code,mod) ### get GPML data back
        
        #os.remove(gpml_path) ### Only store the file temporarily
        try: export.deleteFolder('BuildDBs/WPs/'+processor_time) ### clear any remaining pathway files
        except Exception: pass
        
    end_time = time.time(); time_diff = float(end_time-begin_time)
    """
    try: print "WikiPathways data imported in %d seconds" % time_diff
    except Exception: null=None ### Occurs when transitioning back from the Official Database download window (not sure why) -- TclError: can't invoke "update" command
    """
    return pathway_db

def getHexadecimalColorRanges(fold_db,analysis_type):
    all_folds=[]
    folds_to_gene={}
    for gene in fold_db:
        fold = fold_db[gene]
        all_folds.append(fold)
        try: folds_to_gene[fold].append(gene)
        except Exception: folds_to_gene[fold] = [gene]
    all_folds.sort() ### Sorted range of folds
    if analysis_type == 'Lineage':
        vmax = max(all_folds) ### replaces the old method of getting centered colors
    else:
        try: vmax,vmin=getColorRange(all_folds) ### We want centered colors for this (unlike with Lineage analysis)
        except Exception: vmax,vmin = 0,0
    ### Normalize these values from 0 to 1
    norm_folds_db={}; norm_folds=[]; color_all_yellow = False
    for f in all_folds:
        if analysis_type != 'Lineage':
            try: n=(f-vmin)/(vmax-vmin) ### normalized
            except Exception: n = 1; color_all_yellow = True
        else:
            n=(f-1.0)/(vmax-1.0) ### customized -> 1 is the lowest
        norm_folds_db[n]=f
        norm_folds.append(n)
        
    ### Calculate the color tuple for R, G and then B (blue to red)
    ### 00,00,255 (blue) -> 255,00,00 (red)
    ### 207, 207, 207 (grey) -> 255, 64, 64 (indian red 2) -> http://web.njit.edu/~kevin/rgb.txt.html
    r0 = 207; g0 = 207; b0 = 207
    rn = 255; gn = 64; bn = 64
    
    if color_all_yellow:
        for gene in fold_db:
            fold_db[gene] = 'FFFF00'
    else:
        if analysis_type != 'Lineage':
            ### blue -> grey
            r0 = 0; g0 = 191; b0 = 255
            rn = 207; gn = 207; bn = 207
            ### grey -> red
            r20 = 207; g20 = 207; b20 = 207
            r2n = 255; g2n = 0; b2n = 0
    
        gene_colors_hex = {}
        for n in norm_folds:
            ri=int(r0+(n*(rn-r0)))
            gi=int(g0+(n*(gn-g0)))
            bi=int(b0+(n*(bn-b0)))
            rgb=ri,gi,bi ###blue to grey for non-lineage analyses
            #rgb = (255,0,0) 
            if analysis_type != 'Lineage':
                r2i=int(r20+(n*(r2n-r20)))
                g2i=int(g20+(n*(g2n-g20)))
                b2i=int(b20+(n*(b2n-b20)))
                rgb2=r2i,g2i,b2i ### grey->red
            f = norm_folds_db[n] ### get the original fold
            genes = folds_to_gene[f] ## look up the gene(s) for that fold
            if f<=1 and analysis_type == 'Lineage': ### only show positive z-scores with color
                rgb = (207, 207, 207)
            if f>0 and analysis_type == 'Genes':
                rgb = rgb2
            hex = '#%02x%02x%02x' % rgb
            #print f,n,rgb,hex
            for gene in genes:
                fold_db[gene] = hex[1:]
    return fold_db

def getColorRange(x):
    """ Determines the range of colors, centered at zero, for normalizing cmap """
    vmax=max(x)
    vmin=min(x)
    vmax = max([vmax,abs(vmin)])
    vmin = -1*vmax
    return vmax,vmin
    
def getGraphIDAssociations(id_color_db,pathway_db,key_by):
    graphID_pathway_db={}
    for pathway in pathway_db:
        wpi = pathway_db[pathway] ### all data for the pathway is stored in this single object - wpi.Pathway() is the pathway name
        graphID_pathway_db[pathway,wpi.Pathway()]=db={} ### add a new dictionary key (pathway) and initialize a new dictionary inside
        for gi in wpi.PathwayGeneData():
            if key_by == 'Label':
                if gi.Label() in id_color_db:
                    hex_color = id_color_db[gi.Label()]
                    graphID_pathway_db[pathway,wpi.Pathway()][gi.GraphID()] = hex_color ### set the key,value of the child dictionary
                elif len(wpi.PathwayGeneData())<277: ### more than this and parser.pxi thows an error - memory associated?
                    ### add it as white node
                    graphID_pathway_db[pathway,wpi.Pathway()][gi.GraphID()] = 'FFFFFF'
            else:
                try:
                    for mod_id in gi.ModID():
                        if mod_id in id_color_db:
                            hex_color = id_color_db[mod_id]
                            graphID_pathway_db[pathway,wpi.Pathway()][gi.GraphID()] = hex_color ### set the key,value of the child dictionary
                except Exception: None ### No MOD translation for this ID
    return graphID_pathway_db
            
def viewLineageProfilerResults(filename,graphic_links):
    global graphic_link
    graphic_link=graphic_links ### This is a list of tuples containing name and file location
    
    ### Log any potential problems
    log_file = filepath('webservice.log')
    log_report = open(log_file,'w')
    
    root_dir = export.findParentDir(filename)
    root_dir = string.replace(root_dir,'ExpressionOutput/Clustering','DataPlots')
    if 'DataPlots' not in root_dir: ### Occurs when directly supplying an input matrix by the user
        root_dir += '/DataPlots/'
        try: os.mkdir(root_dir) ### May need to create this directory
        except Exception: None
    id_db,column_headers = importDataSimple(filename,'LineageProfiler')
    log_report.write('LineageProfiler input ID file imported successfully\n')
    pathway_db={}
    pathway_db['WP2062'] = PathwayData('TissueFateMap')
    ### MOD and species are not particularly important for Lineage analysis
    pathway_db = getPathwayAs(pathway_db,'Hs','Ensembl')
    log_report.write('Pathway data imported from GPML files obtained from webservice\n')
    i=0
    group_id_db={} ### store the results separately for each sample
    ### When analyzing z-scores, you can have multiple samples you wish to visualize results for (not so for regulated gene lists)
    for biological_group in column_headers:
        group_id_db[biological_group]=db={}
        for gene in id_db:
            group_id_db[biological_group][gene] = id_db[gene][i] ### get the index value of that biological group (z-score change)
        i+=1
    for biological_group in group_id_db:
        group_specific = group_id_db[biological_group]
        analysis_type = 'Lineage'
        id_color_db = getHexadecimalColorRanges(group_specific,analysis_type) ### example "id_db" is key:tissue, value:z-score
        graphID_db = getGraphIDAssociations(id_color_db,pathway_db,'Label')
        file_type = 'png' ### svg, pdf, png
        getColoredPathway(root_dir,graphID_db,file_type,'-'+biological_group)
        file_type = 'pdf' ### svg, pdf, png
        getColoredPathway(root_dir,graphID_db,file_type,'-'+biological_group)
        
    log_report.write('Pathways colored and images saved to disk. Exiting webservice.\n')
    log_report.close()
    return graphic_link
    
def visualizePathwayAssociations(filename,species,mod_type,wpid,imageExport=True):
    ### Log any potential problems
    log_file = filepath('webservice.log')
    log_report = open(log_file,'w')
    if wpid == None:
        force_invalid_pathway
        
    global mod
    global species_code
    global graphic_link
    graphic_link={}
    mod = mod_type
    species_code = species
    root_dir = export.findParentDir(filename)
    criterion_name = export.findFilename(filename)[:-4]
    log_report.write('Filename: %s and WPID %s\n' % (filename,wpid))
    if 'GO-Elite/input' in root_dir:
        root_dir = string.replace(root_dir,'GO-Elite/input','WikiPathways')
    else:
        root_dir+='WikiPathways/'
    analysis_type = 'Genes'
    id_db,column_headers = importDataSimple(filename,'GO-Elite')
    log_report.write('GO-Elite input ID file imported successfully\n')
    log_report.write('%d IDs imported\n' % len(id_db))
    pathway_db={}
    pathway_db[wpid] = PathwayData(None) ### only need to analyze object (method allows for analysis of any number)
    pathway_db = getPathwayAs(pathway_db,species_code,mod)
    log_report.write('Pathway data imported from GPML files obtained from webservice\n')
    id_color_db = getHexadecimalColorRanges(id_db,analysis_type) ### example id_db" is key:gene, value:fold
    graphID_db = getGraphIDAssociations(id_color_db,pathway_db,'MOD')
    if imageExport != 'png':
        file_type = 'pdf' ### svg, pdf, png
        getColoredPathway(root_dir,graphID_db,file_type,'-'+criterion_name,WPID=wpid)
    if imageExport != 'pdf':
        file_type = 'png' ### svg, pdf, png
        getColoredPathway(root_dir,graphID_db,file_type,'-'+criterion_name,WPID=wpid)
    log_report.write('Pathways colored and image data returned. Exiting webservice.\n')
    log_report.close()
    return graphic_link

def getColoredPathway(root_dir,graphID_db,file_type,dataset_name,WPID=None):
    for (wpid,name) in graphID_db:
        ### Example: graphId="ffffff90"; wpid = "WP2062"; color = "0000ff"
        if WPID==wpid or WPID==None:
            graphID_list = []
            hex_color_list = []
            for graphID in graphID_db[(wpid,name)]:
                graphID_list.append(graphID)
                hex_color_list.append(graphID_db[(wpid,name)][graphID]) ### order is thus the same for both
            #hex_color_list = ["0000ff"]*11
            #print len(graphID_list),graphID_list
            #print len(hex_color_list),hex_color_list
            #print file_type
            if len(graphID_list)==0:
                continue
            ### revision = 0 is the most current version
            #file = client.service.getColoredPathway(pwId=wpid,revision=0,graphId=graphID_list,color=hex_color_list,fileType=file_type)

            file = wikipathways_api_client_instance.get_colored_pathway(identifier=wpid,version=0,
                                element_identifiers=graphID_list,colors=hex_color_list,file_format=file_type)
            #file = base64.b64decode(file) ### decode this file
            name = string.replace(name,':','-')
            name = string.replace(name,'/','-')
            name = string.replace(name,'\\','-') ### will otherwise create a sub-directory
            output_filename = root_dir+wpid+'_'+name+dataset_name+'.'+file_type
            outfile = export.ExportFile(output_filename)
            if file_type == 'png':
                if wpid == 'WP2062': ### This is the LineageMapp
                    graphic_link.append(('LineageProfiler'+dataset_name,output_filename))
                else:
                    graphic_link['WP'] = output_filename
            outfile.write(file); outfile.close()
            #http://au.answers.yahoo.com/question/index?qid=20111029100100AAqxS8l
            #http://stackoverflow.com/questions/2374427/python-2-x-write-binary-output-to-stdout
        
def importDataComplex(filename,input_type,MOD=None,Species=None):
    ### If multiple system codes and mods exist
    None
    
def importDataSimple(filename,input_type,MOD=None,Species=None):
    id_db={}
    fn = filepath(filename)
    x=0
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if data[0]=='#' and x==0: x=0
        elif x==0:
            column_headers = t[1:]
            if input_type != 'LineageProfiler':
                try: column_headers = t[2] ### exclude the ID, system code and p-value column headers
                except Exception: column_headers = 'NA'
            x=1
        else:
            if x==1 and input_type != 'LineageProfiler':
                ### get system conversions
                system_code = t[1]
                import GO_Elite
                from import_scripts import OBO_import
                system_codes,source_types,mod_types = GO_Elite.getSourceData()
                source_data = system_codes[system_code]
                try:
                    Mod=mod ### global established in upstream functions
                    speciescode = species_code
                except Exception:
                    Mod=MOD
                    speciescode = Species
                if source_data == Mod:
                    source_is_mod = True
                else:
                    source_is_mod = False
                    mod_source = Mod+'-'+source_data+'.txt'
                    gene_to_source_id = gene_associations.getGeneToUid(speciescode,('hide',mod_source))
                    source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
            if input_type != 'LineageProfiler':
                if source_is_mod == True:
                    try: id_db[t[0]] = float(t[2])
                    except Exception: id_db[t[0]] = 'NA'
                elif t[0] in source_to_gene:
                    mod_ids = source_to_gene[t[0]]
                    for mod_id in mod_ids:
                        try: value = t[2]
                        except Exception: value = 'NA'
                        if value == '+': value = 1
                        elif value == '-': value = -1
                        try: id_db[mod_id] = float(value) ### If multiple Ensembl IDs in dataset, only record the last associated fold change
                        except Exception: id_db[mod_id] = 'NA'
                        break
            else:
                id_db[t[0]]= map(float,t[1:]) ### Applies to LineageProfiler
            x+=1
    #print len(id_db),column_headers
    return id_db,column_headers

def getColoredPathwayTest():
    fileType = 'png' ### svg, pdf
    graphId="ffffff90"; wpid = "WP2062"; color = "0000ff"
    graphId=["ffffff90","ffffffe5"]
    color = ["0000ff","0000ff"]
    ### revision = 0 is the most current version
    #file = client.service.getColoredPathway(pwId=wpid,revision=0,graphId=graphId,color=color,fileType=fileType)
    kwargs = {
        'identifier': 'WP2062',
        'version': 0,
        'element_identifiers': ["ffffff90","ffffffe5"],
        'colors': ["#0000FF","#0000FF"],
        'file_format': 'image/svg+xml'}
    file = wikipathways_api_client_instance.get_colored_pathway(identifier=wpid,version=0,element_identifiers=graphId,colors=color,file_format=fileType)
    #file = base64.b64decode(file) ### decode this file
    outfile = export.ExportFile(wpid+'.png')
    outfile.write(file); outfile.close()

def getAllSpeciesPathways(species_full):
    #import GO_Elite
    #species_names = GO_Elite.remoteSpeciesData()
    #species_full = string.replace(species_names[species_code],'_',' ')
    #species_full = 'Mycobacterium tuberculosis'; pathway_db = {}
    #pathways_all = client.service.listPathways(organism = species_full)
    pathways_all = wikipathways_api_client_instance.list_pathways(organism = species_full)
    pathway_db={}
    for pathway in pathways_all:
        #wpid = pathway[0]; wpname = pathway[2]
        wpid = pathway['identifier']
        wpname = pathway['name']
        pathway_db[wpid] = PathwayData(wpname)
    return pathway_db

if __name__ == '__main__':
    pathway_db = getAllSpeciesPathways('Homo sapiens');
    for i in pathway_db:
        print i
    
    getPathwayAs(pathway_db,'','');sys.exit()
    getColoredPathwayTest();sys.exit()
    filename = "/Users/saljh8/Desktop/PCBC_MetaData_Comparisons/AltAnalyzeExon/Methylation_Variance/GO-Elite_adjp-2fold/regulated/GE.poor_vs_good-fold2.0_adjp0.05.txt"
    visualizePathwayAssociations(filename,'Hs','Ensembl','WP2857')
    #viewLineageProfilerResults(filename,[]); sys.exit()
    filename = "/Users/nsalomonis/Desktop/code/AltAnalyze/datasets/3'Array/Merrill/GO-Elite/input/GE.ko_vs_wt.txt"
    pathway_db = getAllSpeciesPathways('Homo sapiens')
    for i in pathway_db:
        print i, pathway_db[i].WPName(), len(pathway_db)
        

