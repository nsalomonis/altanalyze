###AltAnalyze
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
#import pkg_resources
#import distutils
from stats_scripts import statistics
import sys, string
import os.path
import unique
import update
import UI
import copy
import export; reload(export)
import ExpressionBuilder; reload(ExpressionBuilder)
from build_scripts import ExonAnalyze_module; reload(ExonAnalyze_module)
from import_scripts import ExonAnnotate_module; reload(ExonAnnotate_module)
from import_scripts import ResultsExport_module
from build_scripts import FeatureAlignment
import GO_Elite
import time
import webbrowser
import random
import traceback
import shutil

try:
    import multiprocessing as mlp
except Exception:
    mlp=None
    print 'Note: Multiprocessing not supported for this verison python.'

try:
    from scipy import stats
except Exception:
    pass ### scipy is not required but is used as a faster implementation of Fisher Exact Test when present

try:
    from PIL import Image as PIL_Image
    try: import ImageTk
    except Exception: from PIL import ImageTk
    import PIL._imaging
    import PIL._imagingft
except Exception:
    print traceback.format_exc()
    pass #print 'Python Imaging Library not installed... using default PNG viewer'
    
use_Tkinter = 'no'
debug_mode = 'no'
analysis_start_time = time.time()

command_args = string.join(sys.argv,' ')
if len(sys.argv[1:])>1 and '-' in command_args and '--GUI' not in command_args:
    runningCommandLine = True
else:
    runningCommandLine = False

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    dir_list2 = [] #add in code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv" or entry[-4:] == ".TXT":
            dir_list2.append(entry)
    return dir_list2

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

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

def returnLargeGlobalVars():
    ### Prints all large global variables retained in memory (taking up space)
    all = [var for var in globals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all:
        try:
            if len(globals()[var])>500:
                print var, len(globals()[var])
        except Exception: null=[]

def clearObjectsFromMemory(db_to_clear):
    db_keys={}
    try:
        for key in db_to_clear: db_keys[key]=[]
    except Exception:
        for key in db_to_clear: del key ### if key is a list
        
    for key in db_keys:
        try: del db_to_clear[key]
        except Exception: 
            try:
                for i in key: del i ### For lists of tuples
            except Exception: del key ### For plain lists

def importGeneric(filename):
    fn=filepath(filename); key_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        key_db[t[0]] = t[1:]
    return key_db

def importGenericFiltered(filename,filter_db):
    fn=filepath(filename); key_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        key = t[0]
        if key in filter_db: key_db[key] = t[1:]
    return key_db

def importGenericFilteredDBList(filename,filter_db):
    fn=filepath(filename); key_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        try:
            null=filter_db[t[0]]
            try: key_db[t[0]].append(t[1])
            except KeyError:  key_db[t[0]] = [t[1]]
        except Exception: null=[]
    return key_db

def importGenericDBList(filename):
    fn=filepath(filename); key_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        try: key_db[t[0]].append(t[1])
        except KeyError:  key_db[t[0]] = [t[1]]
    return key_db

def importExternalDBList(filename):
    fn=filepath(filename); key_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        try: key_db[t[0]].append(t[1:])
        except Exception:  key_db[t[0]] = [t[1:]]
    return key_db

def FindDir(dir,term):
    dir_list = unique.read_directory(dir)
    dir_list2=[]
    dir_list.sort()
    for i in dir_list:
        if term == i: dir_list2.append(i)
    if len(dir_list2)==0:
        for i in dir_list:
            if term in i: dir_list2.append(i)        
    dir_list2.sort(); dir_list2.reverse()
    if len(dir_list2)>0:  return dir_list2[0]
    else: return ''

def openFile(file_dir):
    if runningCommandLine:
        pass
    elif os.name == 'nt':
        try: os.startfile('"'+file_dir+'"')
        except Exception:  os.system('open "'+file_dir+'"')
    elif 'darwin' in sys.platform: os.system('open "'+file_dir+'"')
    elif 'linux' in sys.platform: os.system('xdg-open "'+file_dir+'"')

def openCytoscape(parent_dir,application_dir,application_name):
    cytoscape_dir = FindDir(parent_dir,application_dir); cytoscape_dir = filepath(parent_dir+'/'+cytoscape_dir)
    app_dir = FindDir(cytoscape_dir,application_name)
    app_dir = cytoscape_dir+'/'+app_dir
    if 'linux' in sys.platform:
        app_dir = app_dir
        app_dir2 = cytoscape_dir+'/Cytoscape'
        try: createCytoscapeDesktop(cytoscape_dir)
        except Exception: null=[]
        dir_list = unique.read_directory('/usr/bin/') ### Check to see that JAVA is installed
        if 'java' not in dir_list: print 'Java not referenced in "usr/bin/. If not installed,\nplease install and re-try opening Cytoscape'
        try:
            jar_path = cytoscape_dir+'/cytoscape.jar'
            main_path = cytoscape_dir+'/cytoscape.CyMain'
            plugins_path = cytoscape_dir+'/plugins'
            os.system('java -Dswing.aatext=true -Xss5M -Xmx512M -jar '+jar_path+' '+main_path+'  -p '+plugins_path+' &')
            print 'Cytoscape jar opened:',jar_path
        except Exception:
            print 'OS command to open Java failed.'      
            try:
                try: openFile(app_dir2); print 'Cytoscape opened:',app_dir2
                except Exception:
                    os.chmod(app_dir,0777)
                    openFile(app_dir2)
            except Exception:
                try: openFile(app_dir)
                except Exception:
                    os.chmod(app_dir,0777)
                    openFile(app_dir)
    else:
        try: openFile(app_dir)
        except Exception:
            os.chmod(app_dir,0777)
            openFile(app_dir)

def createCytoscapeDesktop(cytoscape_dir):

    cyto_ds_output = cytoscape_dir+'/Cytoscape.desktop'
    data = export.ExportFile(cyto_ds_output)
    
    cytoscape_desktop = cytoscape_dir+'/Cytoscape'; #cytoscape_desktop = '/hd3/home/nsalomonis/Cytoscape_v2.6.1/Cytoscape'
    cytoscape_png = cytoscape_dir+ '/.install4j/Cytoscape.png'; #cytoscape_png = '/hd3/home/nsalomonis/Cytoscape_v2.6.1/.install4j/Cytoscape.png'
    data.write('[Desktop Entry]'+'\n')
    data.write('Type=Application'+'\n')
    data.write('Name=Cytoscape'+'\n')
    data.write('Exec=/bin/sh "'+cytoscape_desktop+'"'+'\n')
    data.write('Icon='+cytoscape_png+'\n')
    data.write('Categories=Application;'+'\n')
    data.close()

########### Parse Input Annotations ###########
def ProbesetCalls(array_type,probeset_class,splice_event,constitutive_call,external_exonid):
    include_probeset = 'yes'

    if array_type == 'AltMouse':
        exonid = splice_event
        if filter_probesets_by == 'exon':
            if '-' in exonid or '|' in exonid: ###Therfore the probeset represents an exon-exon junction or multi-exon probeset
                include_probeset = 'no'
        if filter_probesets_by != 'exon': 
            if '|' in exonid: include_probeset = 'no'
        if constitutive_call == 'yes': include_probeset = 'yes'
    else:
        if avg_all_for_ss == 'yes' and (probeset_class == 'core' or len(external_exonid)>2): constitutive_call = 'yes'
        #if len(splice_event)>2 and constitutive_call == 'yes' and avg_all_for_ss == 'no':  constitutive_call = 'no'
        if constitutive_call == 'no' and len(splice_event)<2 and len(external_exonid)<2:  ###otherwise these are interesting probesets to keep
            if filter_probesets_by != 'full':
                if filter_probesets_by == 'extended':
                    if probeset_class == 'full': include_probeset = 'no'
                elif filter_probesets_by == 'core':
                    if probeset_class != 'core': include_probeset = 'no'
    return include_probeset,constitutive_call

def EvidenceOfAltSplicing(slicing_annot):
    splice_annotations = ["ntron","xon","strangeSplice","Prime","3","5","C-term"]; as_call = 0
    splice_annotations2 = ["ntron","assette","strangeSplice","Prime","3","5"]
    for annot in splice_annotations:
        if annot in slicing_annot: as_call = 1

    if as_call == 1:
        if "C-term" in slicing_annot and ("N-" in slicing_annot or "Promoter" in slicing_annot):
            as_call = 0
            for annot in splice_annotations2:
                if annot in slicing_annot: as_call = 1
        elif "bleed" in slicing_annot and ("N-" in slicing_annot or "Promoter" in slicing_annot):
            as_call = 0
            for annot in splice_annotations2:
                if annot in slicing_annot: as_call = 1
    return as_call

########### Begin Analyses ###########
class SplicingAnnotationData:
    def ArrayType(self):
        self._array_type = array_type
        return self._array_type
    def Probeset(self): return self._probeset
    def setProbeset(self,probeset): self._probeset = probeset
    def ExonID(self): return self._exonid
    def setDisplayExonID(self,exonid): self._exonid = exonid
    def GeneID(self): return self._geneid
    def Symbol(self):
        symbol = ''
        if self.GeneID() in annotate_db:
            y = annotate_db[self.GeneID()]
            symbol = y.Symbol()
        return symbol
    def ExternalGeneID(self): return self._external_gene
    def ProbesetType(self):
        ###e.g. Exon, junction, constitutive(gene)
        return self._probeset_type
    def GeneStructure(self): return self._block_structure
    def SecondaryExonID(self): return self._block_exon_ids
    def setSecondaryExonID(self,ids): self._block_exon_ids = ids
    def setLocationData(self, chromosome, strand, probeset_start, probeset_stop):
        self._chromosome = chromosome; self._strand = strand
        self._start = probeset_start; self._stop = probeset_stop
    def LocationSummary(self):
        location = self.Chromosome()+':'+self.ProbeStart()+'-'+self.ProbeStop()+'('+self.Strand()+')'
        return location
    def Chromosome(self): return self._chromosome
    def Strand(self): return self._strand
    def ProbeStart(self): return self._start
    def ProbeStop(self): return self._stop
    def ProbesetClass(self): 
        ###e.g. core, extendended, full
        return self._probest_class
    def ExternalExonIDs(self): return self._external_exonids
    def ExternalExonIDList(self):
        external_exonid_list = string.split(self.ExternalExonIDs(),'|')
        return external_exonid_list
    def Constitutive(self): return self._constitutive_status
    def setTranscriptCluster(self,secondary_geneid): self._secondary_geneid = secondary_geneid
    def setNovelExon(self,novel_exon): self._novel_exon = novel_exon
    def NovelExon(self): return self._novel_exon
    def SecondaryGeneID(self): return self._secondary_geneid
    def setExonRegionID(self,exon_region): self._exon_region = exon_region
    def ExonRegionID(self): return self._exon_region
    def SplicingEvent(self):
        splice_event = self._splicing_event
        if len(splice_event)!=0:
            if splice_event[0] == '|': splice_event = splice_event[1:]
        return splice_event
    def SplicingCall(self): return self._splicing_call
    def SpliceJunctions(self): return self._splice_junctions
    def Delete(self): del self
    def Report(self):
        output = self.ArrayType() +'|'+ self.ExonID() +'|'+ self.ExternalGeneID()
        return output
    def __repr__(self): return self.Report()

class AltMouseData(SplicingAnnotationData):
    def __init__(self,affygene,exons,ensembl,block_exon_ids,block_structure,probe_type_call):
        self._geneid = affygene; self._external_gene = ensembl; self._exonid = exons; self._secondary_geneid = ensembl
        self._probeset_type = probe_type_call; self._block_structure = block_structure; self._block_exon_ids = block_exon_ids
        self._external_exonids = 'NA';
        self._constitutive_status = 'no'
        self._splicing_event = ''
        self._secondary_geneid = 'NA'
        self._exon_region = ''
        if self._probeset_type == 'gene': self._constitutive_status = 'yes'
        else: self._constitutive_status = 'no'

class AffyExonSTData(SplicingAnnotationData):
    def __init__(self,ensembl_gene_id,exon_id,ens_exon_ids, constitutive_call_probeset, exon_region, splicing_event, splice_junctions, splicing_call):
        self._geneid = ensembl_gene_id; self._external_gene = ensembl_gene_id; self._exonid = exon_id
        self._constitutive_status = constitutive_call_probeset#; self._start = probeset_start; self._stop = probeset_stop
        self._external_exonids = ens_exon_ids; #self._secondary_geneid = transcript_cluster_id#; self._chromosome = chromosome; self._strand = strand
        self._exon_region=exon_region;  self._splicing_event=splicing_event; self._splice_junctions=splice_junctions; self._splicing_call = splicing_call
        if self._exonid[0] == 'U': self._probeset_type = 'UTR'
        elif self._exonid[0] == 'E': self._probeset_type = 'exonic'
        elif self._exonid[0] == 'I': self._probeset_type = 'intronic'

class AffyExonSTDataAbbreviated(SplicingAnnotationData):
    def __init__(self,ensembl_gene_id,exon_id,splicing_call):
        self._geneid = ensembl_gene_id; self._exonid = exon_id; self._splicing_call = splicing_call

def importSplicingAnnotations(array_type,Species,probeset_type,avg_ss_for_all,root_dir):
    global filter_probesets_by; filter_probesets_by = probeset_type
    global species; species = Species; global avg_all_for_ss; avg_all_for_ss = avg_ss_for_all; global exon_db; exon_db={}
    global summary_data_db; summary_data_db={}; global remove_intronic_junctions; remove_intronic_junctions = 'no'
    if array_type == 'RNASeq':
        probeset_annotations_file = root_dir+'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_junctions.txt'
    else: probeset_annotations_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_probesets.txt'
    filtered_arrayids={};filter_status='no'
    constitutive_probeset_db,exon_db,genes_being_analyzed = importSplicingAnnotationDatabase(probeset_annotations_file,array_type,filtered_arrayids,filter_status)
    return exon_db, constitutive_probeset_db

def importSplicingAnnotationDatabase(filename,array_type,filtered_arrayids,filter_status):
    begin_time = time.time()
    probesets_included_by_new_evidence = 0; export_exon_regions = 'yes'
    if 'fake' in array_type: array_type = string.replace(array_type,'-fake',''); original_arraytype = 'RNASeq'
    else: original_arraytype = array_type
    if filter_status == 'no': global gene_transcript_cluster_db; gene_transcript_cluster_db={}; gene_transcript_cluster_db2={}; global last_exon_region_db; last_exon_region_db = {}
    else: new_exon_db={}
    fn=filepath(filename)
    last_gene = ' '; last_exon_region = ''
    constitutive_probeset_db = {}; constitutive_gene = {}
    count = 0; x = 0; constitutive_original = {}
    #if filter_status == 'yes': exon_db = {}
    if array_type == 'AltMouse':
        for line in open(fn,'rU').xreadlines():
            probeset_data = cleanUpLine(line)  #remove endline
            probeset,affygene,exons,transcript_num,transcripts,probe_type_call,ensembl,block_exon_ids,block_structure,comparison_info = string.split(probeset_data,'\t')
            ###note: currently exclude comparison_info since not applicable for existing analyses
            if x == 0: x = 1
            else:
                if exons[-1] == '|': exons = exons[0:-1]
                if affygene[-1] == '|': affygene = affygene[0:-1]; constitutive_gene[affygene]=[]
                if probe_type_call == 'gene': constitutive_call = 'yes' #looked through the probe annotations and the gene seems to be the most consistent constitutive feature
                else: constitutive_call = 'no'
                include_call,constitutive_call = ProbesetCalls(array_type,'',exons,constitutive_call,'')
                if include_call == 'yes':
                    probe_data = AltMouseData(affygene,exons,ensembl,block_exon_ids,block_structure,probe_type_call)  #this used to just have affygene,exon in the values (1/17/05)
                    exon_db[probeset] = probe_data
                    if filter_status == 'yes':  new_exon_db[probeset] = probe_data
                if constitutive_call == 'yes': constitutive_probeset_db[probeset] = affygene
        genes_being_analyzed = constitutive_gene
    else:
        for line in open(fn,'rU').xreadlines():             
            probeset_data = cleanUpLine(line)  #remove endline
            if x == 0: x = 1
            else:
                try: probeset_id, exon_id, ensembl_gene_id, transcript_cluster_id, chromosome, strand, probeset_start, probeset_stop, affy_class, constitutive_call_probeset, external_exonid, ens_const_exons, exon_region, exon_region_start, exon_region_stop, splicing_event, splice_junctions = string.split(probeset_data,'\t')
                except Exception: print probeset_data;force_error
                if affy_class == 'free': affy_class = 'full' ### Don't know what the difference is
                include_call,constitutive_call = ProbesetCalls(array_type,affy_class,splicing_event,constitutive_call_probeset,external_exonid)
                #if 'ENSG00000163904:E11.5' in probeset_id: print probeset_data
                    #print array_type,affy_class,splicing_event,constitutive_call_probeset,external_exonid,constitutive_call,include_call;kill
                if array_type == 'junction' and '.' not in exon_id: exon_id = string.replace(exon_id,'-','.'); exon_region = string.replace(exon_region,'-','.')
                if ensembl_gene_id != last_gene: new_gene = 'yes'
                else: new_gene = 'no'
                if filter_status == 'no' and new_gene == 'yes':
                    if '.' in exon_id: ### Exclude junctions
                        if '-' not in last_exon_region and 'E' in last_exon_region: last_exon_region_db[last_gene] = last_exon_region
                    else: last_exon_region_db[last_gene] = last_exon_region
                last_gene = ensembl_gene_id
                if len(exon_region)>1: last_exon_region = exon_region ### some probeset not linked to an exon region
                ###Record the transcript clusters assoicated with each gene to annotate the results later on
                if constitutive_call_probeset!=constitutive_call: probesets_included_by_new_evidence +=1#; print probeset_id,[splicing_event],[constitutive_call_probeset];kill
                proceed = 'no'; as_call = 0
                if array_type == 'RNASeq' or array_type == 'junction': include_call = 'yes' ### Constitutive expression is not needed
                if remove_intronic_junctions == 'yes':
                    if 'E' not in exon_id: include_call = 'no' ### Remove junctions that only have splice-sites within an intron or UTR
                if include_call == 'yes' or constitutive_call == 'yes':
                    #if proceed == 'yes':
                    as_call = EvidenceOfAltSplicing(splicing_event)
                    if filter_status == 'no':
                        probe_data = AffyExonSTDataAbbreviated(ensembl_gene_id, exon_id, as_call)
                        if array_type != 'RNASeq':
                            probe_data.setTranscriptCluster(transcript_cluster_id)
                            try:
                                if export_exon_regions == 'yes':
                                    probe_data.setExonRegionID(exon_region)
                            except Exception: null=[]
                    else:
                        probe_data = AffyExonSTData(ensembl_gene_id, exon_id, external_exonid, constitutive_call, exon_region, splicing_event, splice_junctions, as_call)
                        probe_data.setLocationData(chromosome, strand, probeset_start, probeset_stop)
                        if array_type != 'RNASeq':
                            probe_data.setTranscriptCluster(transcript_cluster_id)
                        else:
                            probe_data.setNovelExon(affy_class)
                    if filter_status == 'yes':
                        try: ### saves memory
                            null = filtered_arrayids[probeset_id]
                            new_exon_db[probeset_id] = probe_data
                        except KeyError: null = []
                    else: exon_db[probeset_id] = probe_data
                    if constitutive_call == 'yes' and filter_status == 'no': ###only perform function when initially running
                        constitutive_probeset_db[probeset_id] = ensembl_gene_id
                        try: constitutive_gene[ensembl_gene_id].append(probeset_id)
                        except Exception: constitutive_gene[ensembl_gene_id] = [probeset_id]
                        ###Only consider transcript clusters that make up the constitutive portion of the gene or that are alternatively regulated
                        if array_type != 'RNASeq':
                            try: gene_transcript_cluster_db[ensembl_gene_id].append(transcript_cluster_id)
                            except KeyError: gene_transcript_cluster_db[ensembl_gene_id] = [transcript_cluster_id]
                if constitutive_call_probeset == 'yes' and filter_status == 'no': ###only perform function when initially running
                    try: constitutive_original[ensembl_gene_id].append(probeset_id)
                    except KeyError: constitutive_original[ensembl_gene_id] = [probeset_id]
                    if array_type != 'RNASeq':
                        try: gene_transcript_cluster_db2[ensembl_gene_id].append(transcript_cluster_id)
                        except KeyError: gene_transcript_cluster_db2[ensembl_gene_id] = [transcript_cluster_id]
        ###If no constitutive probesets for a gene as a result of additional filtering (removing all probesets associated with a splice event), add these back
        original_probesets_add = 0; genes_being_analyzed = {}
        for gene in constitutive_gene: genes_being_analyzed[gene]=[]
        
        for gene in constitutive_original:
            if gene not in constitutive_gene:
                genes_being_analyzed[gene] = [gene]
                constitutive_gene[gene]=[]
                original_probesets_add +=1
                gene_transcript_cluster_db[gene] = gene_transcript_cluster_db2[gene]
                for probeset in constitutive_original[gene]: constitutive_probeset_db[probeset] = gene
        
        #if array_type == 'junction' or array_type == 'RNASeq':
        ### Added the below in 1.16!!!
        ### If no constitutive probesets for a gene assigned, assign all gene probesets
        for probeset in exon_db:
            gene = exon_db[probeset].GeneID()
            proceed = 'no'
            exonid = exon_db[probeset].ExonID()
            ### Rather than add all probesets, still filter based on whether the probeset is in an annotated exon
            if 'E' in exonid and 'I' not in exonid and '_' not in exonid: proceed = 'yes'
            if proceed == 'yes':
                if gene not in constitutive_gene:
                    constitutive_probeset_db[probeset] = gene
                    genes_being_analyzed[gene] = [gene]
                    ### DO NOT ADD TO constitutive_gene SINCE WE WANT ALL mRNA ALIGNING EXONS/JUNCTIONS TO BE ADDED!!!!
                    #constitutive_gene[gene]=[]
                    
        gene_transcript_cluster_db = eliminate_redundant_dict_values(gene_transcript_cluster_db)

    #if affygene == 'ENSMUSG00000023089': print [abs(fold_change_log)],[log_fold_cutoff];kill

    if array_type == 'RNASeq':
        import RNASeq
        try: last_exon_region_db = RNASeq.importExonAnnotations(species,'distal-exon','')
        except Exception: null=[]

    constitutive_original=[]; constitutive_gene=[]
    #clearObjectsFromMemory(exon_db); constitutive_probeset_db=[];genes_being_analyzed=[] ### used to evaluate how much memory objects are taking up
    #print 'remove_intronic_junctions:',remove_intronic_junctions
    #print constitutive_gene['ENSMUSG00000031170'];kill ### Determine if avg_ss_for_all is working
    if original_arraytype == 'RNASeq': id_name = 'exon/junction IDs'
    else: id_name = 'array IDs'
    print len(exon_db),id_name,'stored as instances of SplicingAnnotationData in memory'
    #print len(constitutive_probeset_db),'array IDs stored as constititive'
    #print probesets_included_by_new_evidence, 'array IDs were re-annotated as NOT constitutive based on mRNA evidence'
    if array_type != 'AltMouse': print original_probesets_add, 'genes not viewed as constitutive as a result of filtering',id_name,'based on splicing evidence, added back'
    end_time = time.time(); time_diff = int(end_time-begin_time)
    #print filename,"import finished in %d seconds" % time_diff
    if filter_status == 'yes': return new_exon_db
    else:
        summary_data_db['gene_assayed'] = len(genes_being_analyzed)
        try: exportDenominatorGenes(genes_being_analyzed)
        except Exception: null=[]
        return constitutive_probeset_db,exon_db,genes_being_analyzed

def exportDenominatorGenes(genes_being_analyzed):
    goelite_output = root_dir+'GO-Elite/denominator/AS.denominator.txt'
    goelite_data = export.ExportFile(goelite_output)
    systemcode = 'En'
    goelite_data.write("GeneID\tSystemCode\n")
    for gene in genes_being_analyzed:
        if array_type == 'AltMouse':
            try: gene = annotate_db[gene].ExternalGeneID()
            except KeyError: null = []
        goelite_data.write(gene+'\t'+systemcode+'\n')
    try: goelite_data.close()
    except Exception: null=[]

def performExpressionAnalysis(filename,constitutive_probeset_db,exon_db,annotate_db,dataset_name):
    #if analysis_method == 'splicing-index': returnLargeGlobalVars();kill ### used to ensure all large global vars from the reciprocal junction analysis have been cleared from memory
    #returnLargeGlobalVars()
    
    """import list of expression values for arrayids and calculates statistics"""
    global fold_dbase; global original_conditions; global normalization_method
    stats_dbase = {}; fold_dbase={}; ex_db={}; si_db=[]; bad_row_import = {}; count=0
    global array_group_name_db; array_group_name_db = {}
    global array_group_db; array_group_db = {}; 
    global array_raw_group_values; array_raw_group_values = {}; global original_array_names; original_array_names=[]
    global max_replicates; global equal_replicates; global array_group_list
    array_index_list = []  ###Use this list for permutation analysis
    fn=filepath(filename); line_num = 1
    for line in open(fn,'rU').xreadlines():
      data = cleanUpLine(line); t = string.split(data,'\t'); probeset = t[0]
      if t[0]== '#': null=[] ### Don't import line
      elif line_num == 1:
          line_num += 1 #makes this value null for the next loop of actual array data
          ###Below ocucrs if the data is raw opposed to precomputed
          if ':' in t[1]:
              array_group_list = []; x=0 ###gives us an original index value for each entry in the group
              for entry in t[1:]:
                  original_array_names.append(entry)
                  aa = string.split(entry,':')
                  try: array_group,array_name = aa
                  except Exception: array_name = string.join(aa[1:],':'); array_group = aa[0]
                  try:
                      array_group_db[array_group].append(x)
                      array_group_name_db[array_group].append(array_name)
                  except KeyError:
                      array_group_db[array_group] = [x]
                      array_group_name_db[array_group] = [array_name]
                      ### below only occurs with a new group addition
                      array_group_list.append(array_group) #use this to generate comparisons in the below linked function
                  x += 1
          else:
              #try: print data_type
              #except Exception,exception:
              #print exception
              #print traceback.format_exc()
              print_out = 'The AltAnalyze filtered expression file "'+filename+'" is not propperly formatted.\n Review formatting requirements if this file was created by another application.\n'
              print_out += "\nFirst line\n"+line
              try: UI.WarningWindow(print_out,'Exit'); print print_out
              except Exception: print print_out
              badExit()
      else:
          #if probeset in exon_db:
          #if exon_db[probeset].GeneID() == 'ENSG00000139970':
          ###Use the index values from above to assign each expression value to a new database
          temp_group_array = {}
          line_num+=1
          for group in array_group_db:
              if count == 0: array_index_list.append(array_group_db[group])
              for array_index in array_group_db[group]:
                  try: exp_val = float(t[array_index+1])
                  except Exception:
                    if 'Gene_ID' not in line: bad_row_import[probeset]=line; exp_val = 1
                  ###appended is the numerical expression value for each array in the group (temporary array)
                  try: temp_group_array[group].append(exp_val)  #add 1 since probeset is the first column
                  except KeyError: temp_group_array[group] = [exp_val]
          if count == 0: array_index_list.sort(); count = 1
          ####store the group database within the probeset database entry
          try:
              null = exon_db[probeset] ###To conserve memory, don't store any probesets not used for downstream analyses (e.g. not linked to mRNAs)
              #if 'ENSG00000139970' in probeset:
              #print [max_exp]
              #print t[1:];kill
              #max_exp = max(map(float, t[1:]))
              #if len(array_raw_group_values)>10000: break
              #if max_exp>math.log(70,2):
              array_raw_group_values[probeset] = temp_group_array
          except KeyError:
            #print probeset
            pass
              
    print len(array_raw_group_values), 'sequence identifiers imported out of', line_num-1
    if len(bad_row_import)>0:
        print len(bad_row_import), "Rows with an unexplained import error processed and deleted."
        print "Example row:"; x=0
        for i in bad_row_import:
            if x==0: print bad_row_import[i]
            try: del array_raw_group_values[i]
            except Exception: null=[]
            x+=1

    ### If no gene expression reporting probesets were imported, update constitutive_probeset_db to include all mRNA aligning probesets
    cs_genedb={}; missing_genedb={}; addback_genedb={}; rnaseq_cs_gene_db={}
    for probeset in constitutive_probeset_db:
        gene = constitutive_probeset_db[probeset]
        #if gene == 'ENSG00000185008': print [probeset]
        try:
            null=array_raw_group_values[probeset]; cs_genedb[gene]=[]
            if gene == probeset: rnaseq_cs_gene_db[gene]=[] ### If RPKM normalization used, use the gene expression values already calculated
        except Exception: missing_genedb[gene]=[] ### Collect possible that are missing from constitutive database (verify next)
        
    for gene in missing_genedb:
        try: null=cs_genedb[gene]
        except Exception: addback_genedb[gene]=[]
    for probeset in array_raw_group_values:
        try:
            gene = exon_db[probeset].GeneID()
            try:
                null=addback_genedb[gene]
                if 'I' not in probeset and 'U' not in probeset: ### No intron or UTR containing should be used for constitutive expression
                    null=string.split(probeset,':')
                    if len(null)<3: ### No trans-gene junctions should be used for constitutive expression
                        constitutive_probeset_db[probeset]=gene
            except Exception: null=[]
        except Exception: null=[]

    for probeset in constitutive_probeset_db:
        gene = constitutive_probeset_db[probeset]
        #if gene == 'ENSG00000185008': print [[probeset]]
        
    ### Only examine values for associated exons when determining RNASeq constitutive expression (when exon data is present)
    normalization_method = 'raw'
    if array_type == 'RNASeq':
        junction_count=0; constitutive_probeset_db2={}
        for uid in constitutive_probeset_db:
            if '-' in uid: junction_count+=1
        if len(rnaseq_cs_gene_db)>0: ### If filtered RPKM gene-level expression data present, use this instead (and only this)
            normalization_method = 'RPKM'
            constitutive_probeset_db={} ### Re-set this database
            for gene in rnaseq_cs_gene_db:
                constitutive_probeset_db[gene]=gene
        elif junction_count !=0 and len(constitutive_probeset_db) != junction_count:
            ### occurs when there is a mix of junction and exon IDs
            for uid in constitutive_probeset_db:
                if '-' not in uid: constitutive_probeset_db2[uid] = constitutive_probeset_db[uid]
            constitutive_probeset_db = constitutive_probeset_db2; constitutive_probeset_db2=[]
        
    """
    for probeset in constitutive_probeset_db:
        gene = constitutive_probeset_db[probeset]
        if gene == 'ENSG00000185008': print [probeset]
        """

    ###Build all putative splicing events
    global alt_junction_db; global exon_dbase; global critical_exon_db; critical_exon_db={}
    if array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null'):
        ### Applies to reciprocal junction analyses only
        if array_type == 'AltMouse':
            alt_junction_db,critical_exon_db,exon_dbase,exon_inclusion_db,exon_db = ExonAnnotate_module.identifyPutativeSpliceEvents(exon_db,constitutive_probeset_db,array_raw_group_values,agglomerate_inclusion_probesets,onlyAnalyzeJunctions)
            print 'Number of Genes with Examined Splice Events:',len(alt_junction_db)
        elif (array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null':
            from build_scripts import JunctionArray
            alt_junction_db,critical_exon_db,exon_dbase,exon_inclusion_db,exon_db = JunctionArray.getPutativeSpliceEvents(species,array_type,exon_db,agglomerate_inclusion_probesets,root_dir)
            print 'Number of Genes with Examined Splice Events:',len(alt_junction_db)
        #alt_junction_db=[]; critical_exon_db=[]; exon_dbase=[]; exon_inclusion_db=[]

        if agglomerate_inclusion_probesets == 'yes':
            array_raw_group_values = agglomerateInclusionProbesets(array_raw_group_values,exon_inclusion_db)
            exon_inclusion_db=[]
        
        ### For datasets with high memory requirements (RNASeq), filter the current and new databases
        ### Begin this function after agglomeration to ensure agglomerated probesets are considered
        reciprocal_probesets = {}
        if array_type == 'junction' or array_type == 'RNASeq':
            for affygene in alt_junction_db:
                for event in alt_junction_db[affygene]:
                   reciprocal_probesets[event.InclusionProbeset()]=[]
                   reciprocal_probesets[event.ExclusionProbeset()]=[]
            not_evalutated={}
            for probeset in array_raw_group_values:
                try: null=reciprocal_probesets[probeset]
                except Exception:
                    ### Don't remove constitutive probesets
                    try: null=constitutive_probeset_db[probeset]
                    except Exception: not_evalutated[probeset]=[]
            #print 'Removing',len(not_evalutated),'exon/junction IDs not evaulated for splicing'
            for probeset in not_evalutated:
                del array_raw_group_values[probeset]
                
    ###Check to see if we have precomputed expression data or raw to be analyzed
    x=0; y=0; array_raw_group_values2={}; probesets_to_delete=[] ### Record deleted probesets
    if len(array_raw_group_values)==0:
        print_out = "No genes were considered 'Expressed' based on your input options. Check to make sure that the right species database is indicated and that the right data format has been selected (e.g., non-log versus log expression)."
        try: UI.WarningWindow(print_out,'Exit')
        except Exception: print print_out; print "Exiting program"
        badExit()
                
    elif len(array_raw_group_values)>0:
        ###array_group_list should already be unique and correctly sorted (see above)
        for probeset in array_raw_group_values:
            data_lists=[]
            for group_name in array_group_list:
                data_list = array_raw_group_values[probeset][group_name] ###nested database entry access - baseline expression
                if global_addition_factor > 0: data_list = addGlobalFudgeFactor(data_list,'log')
                data_lists.append(data_list)
            if len(array_group_list)==2:
                data_list1 = data_lists[0]; data_list2 = data_lists[-1]; avg1 = statistics.avg(data_list1); avg2 = statistics.avg(data_list2)
                log_fold = avg2 - avg1
                try:
                    #t,df,tails = statistics.ttest(data_list1,data_list2,2,3) #unpaired student ttest, calls p_value function
                    #t = abs(t); df = round(df) #Excel doesn't recognize fractions in a DF
                    #p = statistics.t_probability(t,df)
                    p = statistics.runComparisonStatistic(data_list1,data_list2,probability_statistic)
                    if p == -1:
                        if len(data_list1)>1 and len(data_list2)>1:
                            print_out = "The probability statistic selected ("+probability_statistic+") is not compatible with the\nexperimental design. Please consider an alternative statistic or correct the problem.\nExiting AltAnalyze."
                            try: UI.WarningWindow(print_out,'Exit')
                            except Exception: print print_out; print "Exiting program"
                            badExit()
                        else: p = 1
                except Exception: p = 1
                fold_dbase[probeset] = [0]; fold_dbase[probeset].append(log_fold)
                stats_dbase[probeset]=[avg1]; stats_dbase[probeset].append(p)
                ###replace entries with the two lists for later permutation analysis
                if p == -1: ### should by p == 1: Not sure why this filter was here, but mistakenly removes probesets where there is just one array for each group
                    del fold_dbase[probeset];del stats_dbase[probeset]; probesets_to_delete.append(probeset); x += 1
                    if x == 1: print 'Bad data detected...', data_list1, data_list2
                elif (avg1 < expression_threshold and avg2 < expression_threshold and p > p_threshold) and array_type != 'RNASeq': ### Inserted a filtering option to exclude small variance, low expreession probesets
                    del fold_dbase[probeset];del stats_dbase[probeset]; probesets_to_delete.append(probeset); x += 1
                else: array_raw_group_values2[probeset] = [data_list1,data_list2]
            else: ###Non-junction analysis can handle more than 2 groups
                index=0
                for data_list in data_lists:
                    try: array_raw_group_values2[probeset].append(data_list)
                    except KeyError: array_raw_group_values2[probeset] = [data_list]
                    if len(array_group_list)>2: ### Thus, there is some variance for this probeset
                        ### Create a complete stats_dbase containing all fold changes
                        if index==0:
                            avg_baseline = statistics.avg(data_list); stats_dbase[probeset] = [avg_baseline]
                        else:
                            avg_exp = statistics.avg(data_list)
                            log_fold = avg_exp - avg_baseline
                            try: fold_dbase[probeset].append(log_fold)
                            except KeyError: fold_dbase[probeset] = [0,log_fold]
                    index+=1

    if array_type == 'RNASeq': id_name = 'exon/junction IDs'
    else: id_name = 'array IDs'
    
    array_raw_group_values = array_raw_group_values2; array_raw_group_values2=[]
    print x, id_name,"excluded prior to  analysis... predicted not detected"            

    global original_avg_const_exp_db; global original_fold_dbase
    global avg_const_exp_db; global permute_lists; global midas_db
    if len(array_raw_group_values)>0:
        adj_fold_dbase, nonlog_NI_db, conditions, gene_db, constitutive_gene_db, constitutive_fold_change, original_avg_const_exp_db = constitutive_exp_normalization(fold_dbase,stats_dbase,exon_db,constitutive_probeset_db)
        stats_dbase=[] ### No longer needed after this point
        original_fold_dbase = fold_dbase; avg_const_exp_db = {}; permute_lists = []; y = 0; original_conditions = conditions; max_replicates,equal_replicates = maxReplicates()
        gene_expression_diff_db = constitutive_expression_changes(constitutive_fold_change,annotate_db) ###Add in constitutive fold change filter to assess gene expression for ASPIRE
        while conditions > y:
            avg_const_exp_db = constitutive_exp_normalization_raw(gene_db,constitutive_gene_db,array_raw_group_values,exon_db,y,avg_const_exp_db); y+=1
            #print len(avg_const_exp_db),constitutive_gene_db['ENSMUSG00000054850']
        ###Export Analysis Results for external splicing analysis (e.g. MiDAS format)            
        if run_MiDAS == 'yes' and normalization_method != 'RPKM': ### RPKM has negative values which will crash MiDAS
            status = ResultsExport_module.exportTransitResults(array_group_list,array_raw_group_values,array_group_name_db,avg_const_exp_db,adj_fold_dbase,exon_db,dataset_name,apt_location)
            print "Finished exporting input data for MiDAS analysis"
            try: midas_db = ResultsExport_module.importMidasOutput(dataset_name)
            except Exception: midas_db = {} ### Occurs if there are not enough samples to calculate a MiDAS p-value
        else: midas_db = {}

        ###Provides all pairwise permuted group comparisons
        if array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null'): 
            permute_lists = statistics.permute_arrays(array_index_list)
            ### Now remove probesets from the analysis that were used to evaluate gene expression
            for probeset in constitutive_probeset_db:
                try: null = reciprocal_probesets[probeset]
                except Exception:
                    try: del array_raw_group_values[probeset]
                    except Exception: null=[]
            not_evalutated=[]; reciprocal_probesets=[]
        constitutive_probeset_db=[]
        
        ### Above, all conditions were examined when more than 2 are present... change this so that only the most extreeem are analyzed further
        if len(array_group_list)>2 and analysis_method == 'splicing-index' and (array_type == 'exon' or array_type == 'gene' or explicit_data_type != 'null'): ### USED FOR MULTIPLE COMPARISONS
            print 'Calculating splicing-index values for multiple group comparisons (please be patient)...',

            """
            if len(midas_db)==0:
                print_out = 'Warning!!! MiDAS failed to run for multiple groups. Please make\nsure there are biological replicates present for your groups.\nAltAnalyze requires replicates for multi-group (more than two) analyses.'
                try: UI.WarningWindow(print_out,'Exit')
                except Exception: print print_out; print "Exiting program"
                badExit()"""

            if filter_for_AS == 'yes':
                for probeset in exon_db:
                    as_call = exon_db[probeset].SplicingCall()
                    if as_call == 0:
                        try: del nonlog_NI_db[probeset]
                        except KeyError: null=[]
                
            if export_NI_values == 'yes':
                export_exon_regions = 'yes'
                ### Currently, we don't deal with raw adjusted expression values, just group, so just export the values for each group
                summary_output = root_dir+'AltResults/RawSpliceData/'+species+'/'+analysis_method+'/'+dataset_name[:-1]+'.txt'
                print "Exporting all normalized intensities to:\n"+summary_output
                adjoutput = export.ExportFile(summary_output)
                title = string.join(['Gene\tExonID\tprobesetID']+original_array_names,'\t')+'\n'; adjoutput.write(title)


            ### Pick which data lists have the most extreem values using the NI_dbase (adjusted folds for each condition)
            original_increment = int(len(nonlog_NI_db)/20); increment = original_increment; interaction = 0
            for probeset in nonlog_NI_db:
                if interaction == increment: increment+=original_increment; print '*',
                interaction +=1
                
                geneid = exon_db[probeset].GeneID(); ed = exon_db[probeset]
                index=0; NI_list=[] ### Add the group_name to each adj fold value
                for NI in nonlog_NI_db[probeset]:
                    NI_list.append((NI,index)); index+=1 ### setup to sort for the extreeme adj folds and get associated group_name using the index
                    
                raw_exp_vals = array_raw_group_values[probeset]
                adj_exp_lists={} ### Store the adjusted expression values for each group
                if geneid in avg_const_exp_db:
                    k=0; gi=0; adj_exp_vals = []
                    for exp_list in raw_exp_vals:
                        for exp in exp_list:
                            adj_exp_val = exp-avg_const_exp_db[geneid][k]
                            try: adj_exp_lists[gi].append(adj_exp_val)
                            except Exception: adj_exp_lists[gi] = [adj_exp_val]
                            if export_NI_values == 'yes': adj_exp_vals.append(str(adj_exp_val))
                            k+=1
                        gi+=1
                    if export_NI_values == 'yes':
                        #print geneid+'-'+probeset, adj_exp_val, [ed.ExonID()];kill
                        if export_exon_regions == 'yes':
                            try: ### Thid will only work if ExonRegionID is stored in the abreviated AffyExonSTData object - useful in comparing results between arrays (exon-region centric)
                                if (array_type == 'exon' or array_type == 'gene') or '-' not in ed.ExonID(): ### only include exon entries not junctions
                                    exon_regions = string.split(ed.ExonRegionID(),'|')
                                    for er in exon_regions:
                                        if len(er)>0: er = er
                                        else:
                                            try: er = ed.ExonID()
                                            except Exception: er = 'NA'
                                        ev = string.join([geneid+'\t'+er+'\t'+probeset]+adj_exp_vals,'\t')+'\n'
                                        if len(filtered_probeset_db)>0:
                                            if probeset in filtered_probeset_db: adjoutput.write(ev) ### This is used when we want to restrict to only probesets known to already by changed
                                        else: adjoutput.write(ev)
                            except Exception:
                                ev = string.join([geneid+'\t'+'NA'+'\t'+probeset]+adj_exp_vals,'\t')+'\n'; adjoutput.write(ev)
    
                NI_list.sort()
                examine_pairwise_comparisons = 'yes'
                
                if examine_pairwise_comparisons == 'yes':
                    k1=0; k2=0; filtered_NI_comps = []
                    NI_list_rev = list(NI_list); NI_list_rev.reverse()
                    NI1,index1 = NI_list[k1]; NI2,index2 = NI_list_rev[k2]; abs_SI = abs(math.log(NI1/NI2,2))
                    
                    if abs_SI<alt_exon_logfold_cutoff:
                        ### Indicates that no valid matches were identified - hence, exit loop and return an NI_list with no variance
                        NI_list = [NI_list[0],NI_list[0]]
                    else:
                        ### Indicates that no valid matches were identified - hence, exit loop and return an NI_list with no variance
                        constit_exp1 = original_avg_const_exp_db[geneid][index1]
                        constit_exp2 = original_avg_const_exp_db[geneid][index2]
                        ge_fold = constit_exp2-constit_exp1
                        #print 'original',abs_SI,k1,k2, ge_fold, constit_exp1, constit_exp2
                        if abs(ge_fold) < log_fold_cutoff: filtered_NI_comps.append([abs_SI,k1,k2])
                        else:
                            for i1 in NI_list:
                                k2=0
                                for i2 in NI_list_rev:
                                    NI1,index1 = i1; NI2,index2 = i2; abs_SI = abs(math.log(NI1/NI2,2))
                                    #constit_exp1 = original_avg_const_exp_db[geneid][index1]
                                    #constit_exp2 = original_avg_const_exp_db[geneid][index2]
                                    #ge_fold = constit_exp2-constit_exp1
                                    #if abs(ge_fold) < log_fold_cutoff: filtered_NI_comps.append([abs_SI,k1,k2])
                                    #print k1,k2, i1, i2, abs_SI, abs(ge_fold), log_fold_cutoff, alt_exon_logfold_cutoff
                                    if abs_SI<alt_exon_logfold_cutoff: break
                                    else:
                                        constit_exp1 = original_avg_const_exp_db[geneid][index1]
                                        constit_exp2 = original_avg_const_exp_db[geneid][index2]
                                        ge_fold = constit_exp2-constit_exp1
                                        if abs(ge_fold) < log_fold_cutoff:
                                            filtered_NI_comps.append([abs_SI,k1,k2])
                                            #if k1 == 49 or k1 == 50 or k1 == 51: print probeset, abs_SI, k1, k2, abs(ge_fold),log_fold_cutoff, index1, index2, NI1, NI2, constit_exp1,constit_exp2
                                    k2+=1
                                k1+=1
                    
                    if len(filtered_NI_comps)>0:
                        #print filtered_NI_comps
                        #print NI_list_rev
                        #print probeset,geneid
                        #print len(filtered_NI_comps)
                        #print original_avg_const_exp_db[geneid]
                        filtered_NI_comps.sort()
                        si,k1,k2 = filtered_NI_comps[-1]
                        NI_list = [NI_list[k1],NI_list_rev[k2]]
                        """
                        NI1,index1 = NI_list[0]; NI2,index2 = NI_list[-1]
                        constit_exp1 = original_avg_const_exp_db[geneid][index1]
                        constit_exp2 = original_avg_const_exp_db[geneid][index2]
                        ge_fold = constit_exp2-constit_exp1
                        print probeset, si, ge_fold, NI_list"""
                        #print k1,k2;sys.exit()
  
                index1 = NI_list[0][1]; index2 = NI_list[-1][1]
                nonlog_NI_db[probeset] = [NI_list[0][0],NI_list[-1][0]] ### Update the values of this dictionary
                data_list1 = array_raw_group_values[probeset][index1]; data_list2 = array_raw_group_values[probeset][index2]
                avg1 = statistics.avg(data_list1); avg2 = statistics.avg(data_list2); log_fold = avg2 - avg1
                group_name1 = array_group_list[index1]; group_name2 = array_group_list[index2]
                try:
                    #t,df,tails = statistics.ttest(data_list1,data_list2,2,3) #unpaired student ttest, calls p_value function
                    #t = abs(t); df = round(df); ttest_exp_p = statistics.t_probability(t,df)
                    ttest_exp_p = statistics.runComparisonStatistic(data_list1,data_list2,probability_statistic)
                except Exception: ttest_exp_p = 1                                
                fold_dbase[probeset] = [0]; fold_dbase[probeset].append(log_fold)
                if ttest_exp_p == -1: del fold_dbase[probeset]; probesets_to_delete.append(probeset); x += 1
                elif avg1 < expression_threshold and avg2 < expression_threshold and (ttest_exp_p > p_threshold and ttest_exp_p != 1): ### Inserted a filtering option to exclude small variance, low expreession probesets
                    del fold_dbase[probeset]; probesets_to_delete.append(probeset); x += 1
                else:
                    constit_exp1 = original_avg_const_exp_db[geneid][index1]
                    constit_exp2 = original_avg_const_exp_db[geneid][index2]
                    ge_fold = constit_exp2-constit_exp1
                    normInt1 = (avg1-constit_exp1); normInt2 = (avg2-constit_exp2)
                    adj_fold = normInt2 - normInt1
                    splicing_index = -1*adj_fold; abs_splicing_index = abs(splicing_index)
                    #print probeset, splicing_index, ge_fold, index1, index2
                    #normIntList1 = adj_exp_lists[index1]; normIntList2 = adj_exp_lists[index2]
                    all_nI=[]
                    for g_index in adj_exp_lists: all_nI.append(adj_exp_lists[g_index])
                    try: normIntensityP = statistics.OneWayANOVA(all_nI) #[normIntList1,normIntList2] ### This stays an ANOVA independent of the algorithm choosen since groups number > 2
                    except Exception: normIntensityP = 'NA'
                    
                    if (normInt1*normInt2)<0: opposite_SI_log_mean = 'yes'
                    else: opposite_SI_log_mean = 'no'
          
                    abs_log_ratio = abs(ge_fold)
                    if probeset in midas_db:
                        try: midas_p = float(midas_db[probeset])
                        except ValueError: midas_p = 'NA'
                    else: midas_p = 'NA'
                    
                    #if 'ENSG00000059588' in geneid: print probeset, splicing_index, constit_exp1, constit_exp2, ge_fold,group_name2+'_vs_'+group_name1, index1, index2  
                    if abs_splicing_index>alt_exon_logfold_cutoff and (midas_p < p_threshold or midas_p == 'NA'): #and abs_log_ratio>1 and ttest_exp_p<0.05: ###and ge_threshold_count==2
                        exonid = ed.ExonID(); critical_exon_list = [1,[exonid]]
                        ped = ProbesetExpressionData(avg1, avg2, log_fold, adj_fold, ttest_exp_p, group_name2+'_vs_'+group_name1)
                        sid = ExonData(splicing_index,probeset,critical_exon_list,geneid,normInt1,normInt2,normIntensityP,opposite_SI_log_mean)
                        sid.setConstitutiveExpression(constit_exp1); sid.setConstitutiveFold(ge_fold); sid.setProbesetExpressionData(ped) 
                        si_db.append((splicing_index,sid))
                    else:
                        ### Also record the data for probesets that are excluded... Used by DomainGraph
                        eed = ExcludedExonData(splicing_index,geneid,normIntensityP)
                        ex_db[probeset] = eed

            if array_type == 'RNASeq': id_name = 'exon/junction IDs'
            else: id_name = 'array IDs'
            print len(si_db),id_name,"with evidence of Alternative expression"

            original_fold_dbase = fold_dbase; si_db.sort()
            summary_data_db['denominator_exp_events']=len(nonlog_NI_db)
            del avg_const_exp_db; del gene_db; del constitutive_gene_db; gene_expression_diff_db={}
            if export_NI_values == 'yes': adjoutput.close()

        ### Above, all conditions were examined when more than 2 are present... change this so that only the most extreeem are analyzed further
        elif len(array_group_list)>2 and (array_type == 'junction' or array_type == 'RNASeq' or array_type == 'AltMouse'): ### USED FOR MULTIPLE COMPARISONS
            excluded_probeset_db={}
            group_sizes = []; original_array_indices = permute_lists[0] ###p[0] is the original organization of the group samples prior to permutation
            for group in original_array_indices: group_sizes.append(len(group))
            
            if analysis_method == 'linearregres': ### For linear regression, these scores are non-long
                original_array_raw_group_values = copy.deepcopy(array_raw_group_values)
                for probeset in array_raw_group_values:
                    ls_concatenated=[]
                    for group in array_raw_group_values[probeset]: ls_concatenated+=group
                    ls_concatenated = statistics.log_fold_conversion_fraction(ls_concatenated)
                    array_raw_group_values[probeset] = ls_concatenated
                pos1=0; pos2=0; positions=[]
                for group in group_sizes:
                    if pos1 == 0:  pos2 = group; positions.append((pos1,pos2))
                    else: pos2 = pos1+group; positions.append((pos1,pos2))
                    pos1 = pos2
                    
            if export_NI_values == 'yes':
                export_exon_regions = 'yes'
                ### Currently, we don't deal with raw adjusted expression values, just group, so just export the values for each group
                summary_output = root_dir+'AltResults/RawSpliceData/'+species+'/'+analysis_method+'/'+dataset_name[:-1]+'.txt'
                print "Exporting all normalized intensities to:\n"+summary_output
                adjoutput = export.ExportFile(summary_output)
                title = string.join(['gene\tprobesets\tExonRegion']+original_array_names,'\t')+'\n'; adjoutput.write(title)
                    
            events_examined= 0; denominator_events=0; fold_dbase=[]; adj_fold_dbase=[]; scores_examined=0
            splice_event_list=[]; splice_event_list_mx=[]; splice_event_list_non_mx=[]; event_mx_temp = []; permute_p_values={}; probeset_comp_db={}#use this to exclude duplicate mx events
            for geneid in alt_junction_db:
                affygene = geneid
                for event in alt_junction_db[geneid]:
                    if array_type == 'AltMouse':
                        #event = [('ei', 'E16-E17'), ('ex', 'E16-E18')] 
                        #critical_exon_db[affygene,tuple(critical_exons)] = [1,'E'+str(e1a),'E'+str(e2b)] --- affygene,tuple(event) == key, 1 indicates both are either up or down together
                        event_call = event[0][0] + '-' + event[1][0]
                        exon_set1 = event[0][1]; exon_set2 = event[1][1]
                        probeset1 = exon_dbase[affygene,exon_set1]
                        probeset2 = exon_dbase[affygene,exon_set2]
                        critical_exon_list = critical_exon_db[affygene,tuple(event)]
                    if array_type == 'junction' or array_type == 'RNASeq':
                        event_call = 'ei-ex' ### Below objects from JunctionArrayEnsemblRules - class JunctionInformation
                        probeset1 = event.InclusionProbeset(); probeset2 = event.ExclusionProbeset()
                        exon_set1 = event.InclusionJunction(); exon_set2 = event.ExclusionJunction()
                        try: novel_event = event.NovelEvent()
                        except Exception: novel_event = 'known'
                        critical_exon_list = [1,event.CriticalExonSets()]
                    key,jd = formatJunctionData([probeset1,probeset2],geneid,critical_exon_list[1])
                    if array_type == 'junction' or array_type == 'RNASeq':
                        try: jd.setSymbol(annotate_db[geneid].Symbol())
                        except Exception: null=[]
                    #if '|' in probeset1: print probeset1, key,jd.InclusionDisplay();kill
                    probeset_comp_db[key] = jd ### This is used for the permutation analysis and domain/mirBS import
                    dI_scores=[]
                    if probeset1 in nonlog_NI_db and probeset2 in nonlog_NI_db and probeset1 in array_raw_group_values and probeset2 in array_raw_group_values:
                        events_examined+=1
                        if analysis_method == 'ASPIRE':
                            index1=0; NI_list1=[]; NI_list2=[] ### Add the group_name to each adj fold value
                            for NI in nonlog_NI_db[probeset1]: NI_list1.append(NI)
                            for NI in nonlog_NI_db[probeset2]: NI_list2.append(NI)
                            for NI1_g1 in NI_list1: 
                                NI2_g1 = NI_list2[index1]; index2=0
                                for NI1_g2 in NI_list1:
                                    try: NI2_g2 = NI_list2[index2]
                                    except Exception: print index1, index2, NI_list1, NI_list2;kill
                                    if index1 != index2:
                                        b1 = NI1_g1; e1 = NI1_g2
                                        b2 = NI2_g1; e2 = NI2_g2  
                                        try:
                                            dI = statistics.aspire_stringent(b1,e1,b2,e2); Rin = b1/e1; Rex = b2/e2
                                            if (Rin>1 and Rex<1) or (Rin<1 and Rex>1):
                                                if dI<0: i1,i2 = index2,index1 ### all scores should indicate upregulation
                                                else: i1,i2=index1,index2
                                                dI_scores.append((abs(dI),i1,i2))
                                        except Exception:
                                            #if array_type != 'RNASeq': ### RNASeq has counts of zero and one that can cause the same result between groups and probesets
                                            #print probeset1, probeset2, b1, e1, b2, e2, index1, index2, events_examined;kill
                                            ### Exception - Occurs for RNA-Seq but can occur for array data under extreemly rare circumstances (Rex=Rin even when different b1,e1 and b2,ed values)
                                            null=[]
                                    index2+=1
                                index1+=1
                            dI_scores.sort()
                        if analysis_method == 'linearregres':        
                            log_fold,i1,i2 = getAllPossibleLinearRegressionScores(probeset1,probeset2,positions,group_sizes)
                            dI_scores.append((log_fold,i1,i2))
                            raw_exp_vals1 = original_array_raw_group_values[probeset1]; raw_exp_vals2 = original_array_raw_group_values[probeset2]
                        else: raw_exp_vals1 = array_raw_group_values[probeset1]; raw_exp_vals2 = array_raw_group_values[probeset2]
                        adj_exp_lists1={}; adj_exp_lists2={} ### Store the adjusted expression values for each group
                        if geneid in avg_const_exp_db:
                            gi=0; l=0; adj_exp_vals = []; anova_test=[]
                            for exp_list in raw_exp_vals1:
                                k=0; anova_group=[]
                                for exp in exp_list:
                                    adj_exp_val1 = exp-avg_const_exp_db[geneid][l]
                                    try: adj_exp_lists1[gi].append(adj_exp_val1)
                                    except Exception: adj_exp_lists1[gi] = [adj_exp_val1]
                                    adj_exp_val2 = raw_exp_vals2[gi][k]-avg_const_exp_db[geneid][l]
                                    try: adj_exp_lists2[gi].append(adj_exp_val2)
                                    except Exception: adj_exp_lists2[gi] = [adj_exp_val2]
                                    anova_group.append(adj_exp_val2-adj_exp_val1)
                                    if export_NI_values == 'yes':
                                            #if analysis_method == 'ASPIRE':
                                            adj_exp_vals.append(str(adj_exp_val2-adj_exp_val1))
                                            ### BELOW CODE PRODUCES THE SAME RESULT!!!!
                                            """folds1 = statistics.log_fold_conversion_fraction([exp])
                                            folds2 = statistics.log_fold_conversion_fraction([raw_exp_vals2[gi][k]])
                                            lr_score = statistics.convert_to_log_fold(statistics.simpleLinRegress(folds1,folds2))
                                            adj_exp_vals.append(str(lr_score))"""
                                    k+=1; l+=0
                                gi+=1; anova_test.append(anova_group)
                            if export_NI_values == 'yes':
                                if export_exon_regions == 'yes':
                                    exon_regions = string.join(critical_exon_list[1],'|')
                                    exon_regions = string.split(exon_regions,'|')
                                    for er in exon_regions:
                                        ev = string.join([geneid+'\t'+probeset1+'-'+probeset2+'\t'+er]+adj_exp_vals,'\t')+'\n'
                                        if len(filtered_probeset_db)>0:
                                            if probeset1 in filtered_probeset_db and probeset2 in filtered_probeset_db:
                                                adjoutput.write(ev) ### This is used when we want to restrict to only probesets known to already by changed
                                        else: adjoutput.write(ev)
                            try: anovaNIp = statistics.OneWayANOVA(anova_test) ### This stays an ANOVA independent of the algorithm choosen since groups number > 2
                            except Exception: anovaNIp='NA'
                            
                        if len(dI_scores)>0 and geneid in avg_const_exp_db:
                          dI,index1,index2 = dI_scores[-1]; count=0
                          probesets = [probeset1, probeset2]; index=0
     
                          key,jd = formatJunctionData([probeset1,probeset2],affygene,critical_exon_list[1])
                          if array_type == 'junction' or array_type == 'RNASeq':
                              try: jd.setSymbol(annotate_db[affygene].Symbol())
                              except Exception:null=[]
                          probeset_comp_db[key] = jd ### This is used for the permutation analysis and domain/mirBS import
                            
                          if max_replicates >2 or equal_replicates==2: permute_p_values[(probeset1,probeset2)] = [anovaNIp, 'NA', 'NA', 'NA']
                          index=0
                          for probeset in probesets:  
                            if analysis_method == 'linearregres':
                                data_list1 = original_array_raw_group_values[probeset][index1]; data_list2 = original_array_raw_group_values[probeset][index2]
                            else: data_list1 = array_raw_group_values[probeset][index1]; data_list2 = array_raw_group_values[probeset][index2]
                            baseline_exp = statistics.avg(data_list1); experimental_exp = statistics.avg(data_list2); fold_change = experimental_exp - baseline_exp
                            group_name1 = array_group_list[index1]; group_name2 = array_group_list[index2]
                            try:
                                ttest_exp_p = statistics.runComparisonStatistic(data_list1,data_list2,probability_statistic)
                            except Exception: ttest_exp_p = 'NA'
                            if ttest_exp_p==1: ttest_exp_p = 'NA'
                            if index == 0:
                                try: adj_fold = statistics.avg(adj_exp_lists1[index2]) - statistics.avg(adj_exp_lists1[index1])
                                except Exception:
                                    print raw_exp_vals1,raw_exp_vals2, avg_const_exp_db[geneid]
                                    print probeset,probesets,adj_exp_lists1,adj_exp_lists2,index1,index2;kill
                                ped1 = ProbesetExpressionData(baseline_exp, experimental_exp, fold_change, adj_fold, ttest_exp_p, group_name2+'_vs_'+group_name1)
                            else:
                                adj_fold = statistics.avg(adj_exp_lists2[index2]) - statistics.avg(adj_exp_lists2[index1])
                                ped2 = ProbesetExpressionData(baseline_exp, experimental_exp, fold_change, adj_fold, ttest_exp_p, group_name2+'_vs_'+group_name1)
                                                            
                            constit_exp1 = original_avg_const_exp_db[geneid][index1]
                            constit_exp2 = original_avg_const_exp_db[geneid][index2]
                            ge_fold = constit_exp2-constit_exp1
                            index+=1

                          try:
                            pp1 = statistics.runComparisonStatistic(adj_exp_lists1[index1], adj_exp_lists1[index2],probability_statistic)
                            pp2 = statistics.runComparisonStatistic(adj_exp_lists2[index1], adj_exp_lists2[index2],probability_statistic)
                          except Exception:  pp1 = 'NA'; pp2 = 'NA'
                          if analysis_method == 'ASPIRE' and len(dI_scores)>0:
                              p1 = JunctionExpressionData(adj_exp_lists1[index1], adj_exp_lists1[index2], pp1, ped1)
                              p2 = JunctionExpressionData(adj_exp_lists2[index1], adj_exp_lists2[index2], pp2, ped2)
                              ### ANOVA p-replaces the below p-value
                              """try: baseline_scores, exp_scores, pairwiseNIp = calculateAllASPIREScores(p1,p2)
                              except Exception: baseline_scores = [0]; exp_scores=[dI]; pairwiseNIp = 0 """
                          #if pairwiseNIp == 'NA': pairwiseNIp = 0 ### probably comment out
                              
                        if len(dI_scores)>0:
                            scores_examined+=1
                            if probeset in midas_db:
                                try: midas_p = float(midas_db[probeset])
                                except ValueError: midas_p = 'NA'
                            else: midas_p = 'NA'
                            if dI>alt_exon_logfold_cutoff and (anovaNIp < p_threshold or perform_permutation_analysis == 'yes' or anovaNIp == 'NA' or anovaNIp == 1): #and abs_log_ratio>1 and ttest_exp_p<0.05: ###and ge_threshold_count==2
                                #print [dI, probeset1,probeset2, anovaNIp, alt_exon_logfold_cutoff];kill
                                ejd = ExonJunctionData(dI,probeset1,probeset2,pp1,pp2,'upregulated',event_call,critical_exon_list,affygene,ped1,ped2)
                                ejd.setConstitutiveFold(ge_fold); ejd.setConstitutiveExpression(constit_exp1)
                                if array_type == 'RNASeq':
                                    ejd.setNovelEvent(novel_event)
                                splice_event_list.append((dI,ejd))
                            else: excluded_probeset_db[affygene+':'+critical_exon_list[1][0]] = probeset1, affygene, dI, 'NA', anovaNIp
                            
            statistics.adjustPermuteStats(permute_p_values)
            ex_db = splice_event_list, probeset_comp_db, permute_p_values, excluded_probeset_db
            original_fold_dbase = fold_dbase; original_avg_const_exp_db=[]; nonlog_NI_db = []; fold_dbase=[]
            summary_data_db['denominator_exp_events']=events_examined
            del avg_const_exp_db; del gene_db; del constitutive_gene_db; gene_expression_diff_db={}
            if export_NI_values == 'yes': adjoutput.close()
            print len(splice_event_list), 'alternative exons out of %s exon events examined' % events_examined
    fold_dbase=[]; original_fold_dbase=[]; exon_db=[]; constitutive_gene_db=[]; addback_genedb=[]
    gene_db=[]; missing_genedb=[]
    """
    print 'local vars'
    all = [var for var in locals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all:
        try:
            if len(locals()[var])>500: print var, len(locals()[var])
        except Exception: null=[]
    """
    return conditions,adj_fold_dbase,nonlog_NI_db,dataset_name,gene_expression_diff_db,midas_db,ex_db,si_db

class ProbesetExpressionData:
    def __init__(self, baseline_exp, experimental_exp, fold_change, adj_fold, ttest_raw_exp, annotation):
        self.baseline_exp = baseline_exp; self.experimental_exp = experimental_exp
        self.fold_change = fold_change; self.adj_fold = adj_fold
        self.ttest_raw_exp = ttest_raw_exp; self.annotation = annotation
    def BaselineExp(self): return str(self.baseline_exp)
    def ExperimentalExp(self): return str(self.experimental_exp)
    def FoldChange(self): return str(self.fold_change)
    def AdjFold(self): return str(self.adj_fold)
    def ExpPval(self): return str(self.ttest_raw_exp)
    def Annotation(self): return self.annotation
    def __repr__(self): return self.BaselineExp()+'|'+FoldChange()
    
def agglomerateInclusionProbesets(array_raw_group_values,exon_inclusion_db):
    ###Combine expression profiles for inclusion probesets that correspond to the same splice event
    for excl_probeset in exon_inclusion_db:
        inclusion_event_profiles = []
        if len(exon_inclusion_db[excl_probeset])>1:
            for incl_probeset in exon_inclusion_db[excl_probeset]:
                if incl_probeset in array_raw_group_values and excl_probeset in array_raw_group_values:
                    array_group_values = array_raw_group_values[incl_probeset]
                    inclusion_event_profiles.append(array_group_values)
                    #del array_raw_group_values[incl_probeset]  ###Remove un-agglomerated original entry
        if len(inclusion_event_profiles) > 0: ###Thus, some probesets for this splice event in input file
            combined_event_profile = combine_profiles(inclusion_event_profiles)
            ###Combine inclusion probesets into a single ID (identical manner to that in ExonAnnotate_module.identifyPutativeSpliceEvents
            incl_probesets = exon_inclusion_db[excl_probeset]
            incl_probesets_str = string.join(incl_probesets,'|')
            array_raw_group_values[incl_probesets_str] = combined_event_profile
    return array_raw_group_values

def combine_profiles(profile_list):
    profile_group_sizes={}
    for db in profile_list:
        for key in db: profile_group_sizes[key] = len(db[key])
        break

    new_profile_db={}
    for key in profile_group_sizes:
        x = profile_group_sizes[key] ###number of elements in list for key
        new_val_list=[]; i = 0
        while i<x:
            temp_val_list=[]
            for db in profile_list:
                if key in db: val = db[key][i]; temp_val_list.append(val)
            i+=1; val_avg = statistics.avg(temp_val_list); new_val_list.append(val_avg)
        new_profile_db[key] = new_val_list
    return new_profile_db

def constitutive_exp_normalization(fold_db,stats_dbase,exon_db,constitutive_probeset_db):
    """For every expression value, normalize to the expression of the constitutive gene features for that condition,
    then store those ratios (probeset_exp/avg_constitutive_exp) and regenerate expression values relative only to the
    baseline avg_constitutive_exp, for all conditions, to normalize out gene expression changes"""
    #print "\nParameters:"
    #print "Factor_out_expression_changes:",factor_out_expression_changes
    #print "Only_include_constitutive_containing_genes:",only_include_constitutive_containing_genes
    #print "\nAdjusting probeset average intensity values to factor out condition specific expression changes for optimal splicing descrimination"
    gene_db = {}; constitutive_gene_db = {}
    ### organize everything by gene
    for probeset in fold_db: conditions = len(fold_db[probeset]); break
    remove_diff_exp_genes = remove_transcriptional_regulated_genes
    if conditions > 2: remove_diff_exp_genes = 'no'
    
    for probeset in exon_db:
      affygene = exon_db[probeset].GeneID() #exon_db[probeset] = affygene,exons,ensembl,block_exon_ids,block_structure,comparison_info
      if probeset in fold_db:
        try: gene_db[affygene].append(probeset)
        except KeyError: gene_db[affygene] = [probeset]
        if probeset in constitutive_probeset_db and (only_include_constitutive_containing_genes == 'yes' or factor_out_expression_changes == 'no'):
            #the second conditional is used to exlcude constitutive data if we wish to use all probesets for
            #background normalization rather than just the designated 'gene' probesets.
            if probeset in stats_dbase:
                try: constitutive_gene_db[affygene].append(probeset)
                except KeyError: constitutive_gene_db[affygene] = [probeset]

    if len(constitutive_gene_db)>0:
        ###This is blank when there are no constitutive and the above condition is implemented
        gene_db2 = constitutive_gene_db
    else: gene_db2 = gene_db
    avg_const_exp_db = {}

    for affygene in gene_db2:
        probeset_list = gene_db2[affygene]
        x = 0
        while x < conditions:
            ### average all exp values for constitutive probesets for each condition
            exp_list=[]
            for probeset in probeset_list:
                probe_fold_val = fold_db[probeset][x]
                baseline_exp = stats_dbase[probeset][0]
                exp_val = probe_fold_val + baseline_exp
                exp_list.append(exp_val)
            avg_const_exp = statistics.avg(exp_list)
            try: avg_const_exp_db[affygene].append(avg_const_exp)
            except KeyError: avg_const_exp_db[affygene] = [avg_const_exp]
            x += 1
    
    adj_fold_dbase={}; nonlog_NI_db={}; constitutive_fold_change={}
    for affygene in avg_const_exp_db:   ###If we only wish to include propper constitutive probes, this will ensure we only examine those genes and probesets that are constitutive
        probeset_list = gene_db[affygene]
        x = 0
        while x < conditions:
            exp_list=[]
            for probeset in probeset_list:
                expr_to_subtract = avg_const_exp_db[affygene][x]
                baseline_const_exp = avg_const_exp_db[affygene][0]
                probe_fold_val = fold_db[probeset][x]
                baseline_exp = stats_dbase[probeset][0]
                exp_val = probe_fold_val + baseline_exp
    
                exp_val_non_log = statistics.log_fold_conversion_fraction(exp_val)
                expr_to_subtract_non_log = statistics.log_fold_conversion_fraction(expr_to_subtract)
                baseline_const_exp_non_log = statistics.log_fold_conversion_fraction(baseline_const_exp)
                if factor_out_expression_changes == 'yes':
                    exp_splice_valff = exp_val_non_log/expr_to_subtract_non_log
                else: #if no, then we just normalize to the baseline constitutive expression in order to keep gene expression effects (useful if you don't trust constitutive feature expression levels)
                    exp_splice_valff = exp_val_non_log/baseline_const_exp_non_log
                
                constitutive_fold_diff = expr_to_subtract_non_log/baseline_const_exp_non_log
                ###To calculate adjusted expression, we need to get the fold change in the constitutive avg (expr_to_subtract/baseline_const_exp) and divide the experimental expression
                ###By this fold change.
                ge_adj_exp_non_log = exp_val_non_log/constitutive_fold_diff #gives a GE adjusted expression
                
                try: ge_adj_exp = math.log(ge_adj_exp_non_log,2)
                except ValueError: print probeset,ge_adj_exp_non_log,constitutive_fold_diff,exp_val_non_log,exp_val,baseline_exp, probe_fold_val, dog
                adj_probe_fold_val = ge_adj_exp - baseline_exp
                ### Here we normalize probeset expression to avg-constitutive expression by dividing probe signal by avg const.prove sig (should be < 1)
                ### refered to as steady-state normalization
                if array_type != 'AltMouse' or (probeset not in constitutive_probeset_db):
                    """Can't use constitutive gene features since these have no variance for pearson analysis
                    Python will approximate numbers to a small decimal point range. If the first fold value is
                    zero, often, zero will be close to but not exactly zero. Correct below """
                    try:
                        adj_fold_dbase[probeset].append(adj_probe_fold_val)
                    except KeyError:
                        if abs(adj_probe_fold_val - 0) < 0.0000001: #make zero == exactly to zero
                            adj_probe_fold_val = 0
                        adj_fold_dbase[probeset] = [adj_probe_fold_val]

                try: nonlog_NI_db[probeset].append(exp_splice_valff) ###ratio of junction exp relative to gene expression at that time-point              
                except KeyError: nonlog_NI_db[probeset] = [exp_splice_valff]
                n = 0
                #if expr_to_subtract_non_log != baseline_const_exp_non_log: ###otherwise this is the first value in the expression array
                if x!=0: ###previous expression can produce errors when multiple group averages have identical values
                    fold_change = expr_to_subtract_non_log/baseline_const_exp_non_log
                    fold_change_log = math.log(fold_change,2)
                    constitutive_fold_change[affygene] = fold_change_log
                    ### If we want to remove any genes from the analysis with large transcriptional changes
                    ### that may lead to false positive splicing calls (different probeset kinetics)
                    if remove_diff_exp_genes == 'yes':
                        if abs(fold_change_log) > log_fold_cutoff:
                            del constitutive_fold_change[affygene]
                            try: del adj_fold_dbase[probeset]
                            except KeyError: n = 1
                            try: del nonlog_NI_db[probeset]
                            except KeyError: n = 1
                """elif expr_to_subtract_non_log == baseline_const_exp_non_log:  ###This doesn't make sense, since n can't equal 1 if the conditional is false (check this code again later 11/23/07)
                    if n == 1:
                        del adj_fold_dbase[probeset]
                        del nonlog_NI_db[probeset]"""
            x += 1

    print "Intensity normalization complete..."

    if factor_out_expression_changes == 'no':
        adj_fold_dbase = fold_db #don't change expression values
    print len(constitutive_fold_change), "genes undergoing analysis for alternative splicing/transcription"
    summary_data_db['denominator_exp_genes']=len(constitutive_fold_change)
    """    
    mir_gene_count = 0
    for gene in constitutive_fold_change:
        if gene in gene_microRNA_denom: mir_gene_count+=1
    print mir_gene_count, "Genes with predicted microRNA binding sites undergoing analysis for alternative splicing/transcription"
    """

    global gene_analyzed; gene_analyzed = len(constitutive_gene_db)
    return adj_fold_dbase, nonlog_NI_db, conditions, gene_db, constitutive_gene_db,constitutive_fold_change, avg_const_exp_db

class TranscriptionData:
    def __init__(self, constitutive_fold, rna_processing_annotation):
        self._constitutive_fold = constitutive_fold; self._rna_processing_annotation = rna_processing_annotation
    def ConstitutiveFold(self): return self._constitutive_fold
    def ConstitutiveFoldStr(self): return str(self._constitutive_fold)
    def RNAProcessing(self): return self._rna_processing_annotation
    def __repr__(self): return self.ConstitutiveFoldStr()+'|'+RNAProcessing()
        
def constitutive_expression_changes(constitutive_fold_change,annotate_db):
    ###Add in constitutive fold change filter to assess gene expression for ASPIRE
    gene_expression_diff_db = {}
    for affygene in constitutive_fold_change:
        constitutive_fold = constitutive_fold_change[affygene]; rna_processing_annotation=''
        if affygene in annotate_db:
            if len(annotate_db[affygene].RNAProcessing()) > 4: rna_processing_annotation = annotate_db[affygene].RNAProcessing()
        ###Add in evaluation of RNA-processing/binding factor
        td = TranscriptionData(constitutive_fold,rna_processing_annotation)
        gene_expression_diff_db[affygene] = td
    return gene_expression_diff_db

def constitutive_exp_normalization_raw(gene_db,constitutive_gene_db,array_raw_group_values,exon_db,y,avg_const_exp_db):
    """normalize expression for raw expression data (only for non-baseline data)"""
    #avg_true_const_exp_db[affygene] = [avg_const_exp]
    temp_avg_const_exp_db={}
    
    for probeset in array_raw_group_values:
        conditions = len(array_raw_group_values[probeset][y]); break #number of raw expresson values to normalize

    for affygene in gene_db:
        ###This is blank when there are no constitutive or the above condition is implemented
        if affygene in constitutive_gene_db:
            probeset_list = constitutive_gene_db[affygene]  
            z = 1
        else:  ###so we can analyze splicing independent of gene expression even if no 'gene' feature is present
            probeset_list = gene_db[affygene]
            z = 0
        x = 0
        while x < conditions:
            ### average all exp values for constitutive probesets for each conditionF
            exp_list=[]
            for probeset in probeset_list:
                try: exp_val = array_raw_group_values[probeset][y][x] ### try statement is used for constitutive probes that were deleted due to filtering in performExpressionAnalysis
                except KeyError: continue
                exp_list.append(exp_val)
            try: avg_const_exp = statistics.avg(exp_list) 
            except Exception: avg_const_exp = 'null'
            if only_include_constitutive_containing_genes == 'yes' and avg_const_exp != 'null': 
                if z == 1:
                    try: avg_const_exp_db[affygene].append(avg_const_exp)
                    except KeyError: avg_const_exp_db[affygene] = [avg_const_exp]
                    try: temp_avg_const_exp_db[affygene].append(avg_const_exp)
                    except KeyError: temp_avg_const_exp_db[affygene] = [avg_const_exp]
            elif avg_const_exp != 'null': ###***
                try: avg_const_exp_db[affygene].append(avg_const_exp)
                except KeyError: avg_const_exp_db[affygene] = [avg_const_exp]
                try: temp_avg_const_exp_db[affygene].append(avg_const_exp)
                except KeyError: temp_avg_const_exp_db[affygene] = [avg_const_exp]
            x += 1        
        
    if analysis_method == 'ANOVA':
        global normalized_raw_exp_ratios; normalized_raw_exp_ratios = {}
        for affygene in gene_db:
            probeset_list = gene_db[affygene]
            for probeset in probeset_list:
                while x < group_size:
                    new_ratios = [] ### Calculate expression ratios relative to constitutive expression
                    exp_val = array_raw_group_values[probeset][y][x]
                    const_exp_val = temp_avg_const_exp_db[affygene][x]
                    ###Since the above dictionary is agglomerating all constitutive expression values for permutation,
                    ###we need an unbiased way to grab just those relevant const. exp. vals. (hence the temp dictionary)
                    #non_log_exp_val = statistics.log_fold_conversion_fraction(exp_val)
                    #non_log_const_exp_val = statistics.log_fold_conversion_fraction(const_exp_val)
                    #non_log_exp_ratio = non_log_exp_val/non_log_const_exp_val
                    log_exp_ratio = exp_val - const_exp_val
                    try: normalized_raw_exp_ratios[probeset].append(log_exp_ratio)
                    except KeyError: normalized_raw_exp_ratios[probeset] = [log_exp_ratio]
    return avg_const_exp_db

######### Z Score Analyses #######

class ZScoreData:
    def __init__(self,element,changed,measured,zscore,null_z,gene_symbols):
        self._element = element; self._changed = changed; self._measured = measured
        self._zscore = zscore; self._null_z = null_z; self._gene_symbols = gene_symbols
    def ElementID(self): return self._element
    def Changed(self): return str(self._changed)
    def Measured(self): return str(self._measured)
    def AssociatedWithElement(self): return str(self._gene_symbols)
    def ZScore(self): return str(self._zscore)
    def SetP(self,p): self._permute_p = p
    def PermuteP(self): return str(self._permute_p)
    def SetAdjP(self,adjp): self._adj_p = adjp
    def AdjP(self): return str(self._adj_p)
    def PercentChanged(self):
        try: pc = float(self.Changed())/float(self.Measured())*100
        except Exception: pc = 0
        return str(pc)
    def NullZ(self): return self._null_z
    def Report(self):
        output = self.ElementID()
        return output
    def __repr__(self): return self.Report()

class FDRStats(ZScoreData):
    def __init__(self,p): self._permute_p = p
    def AdjP(self): return str(self._adj_p)
        
def countGenesForElement(permute_input_list,probeset_to_gene,probeset_element_db):
    element_gene_db={}
    for probeset in permute_input_list:
        try:
            element_list = probeset_element_db[probeset]
            gene = probeset_to_gene[probeset]
            for element in element_list:
                try: element_gene_db[element].append(gene)
                except KeyError: element_gene_db[element] = [gene]
        except KeyError: null=[]

    ### Count the number of unique genes per element
    for element in element_gene_db:
        t = {}
        for i in element_gene_db[element]: t[i]=[]
        element_gene_db[element] = len(t)
    return element_gene_db

def formatGeneSymbolHits(geneid_list):
    symbol_list=[]
    for geneid in geneid_list:
        symbol = ''
        if geneid in annotate_db: symbol = annotate_db[geneid].Symbol()
        if len(symbol)<1: symbol = geneid
        symbol_list.append(symbol)
    symbol_str = string.join(symbol_list,', ')
    return symbol_str

def zscore(r,n,N,R):
    z = (r - n*(R/N))/math.sqrt(n*(R/N)*(1-(R/N))*(1-((n-1)/(N-1)))) #z = statistics.zscore(r,n,N,R)
    return z

def calculateZScores(hit_count_db,denom_count_db,total_gene_denom_count,total_gene_hit_count,element_type):
    N = float(total_gene_denom_count)  ###Genes examined
    R = float(total_gene_hit_count)  ###AS genes
    for element in denom_count_db:
        element_denom_gene_count = denom_count_db[element]
        n = float(element_denom_gene_count) ###all genes associated with element
        if element in hit_count_db:
            element_hit_gene_count = len(hit_count_db[element])
            gene_symbols = formatGeneSymbolHits(hit_count_db[element])
            r = float(element_hit_gene_count) ###regulated genes associated with element
        else: r = 0; gene_symbols = ''
        
        try: z = zscore(r,n,N,R)
        except Exception: z = 0; #print 'error:',element,r,n,N,R; kill
        try: null_z = zscore(0,n,N,R)
        except Exception: null_z = 0; #print 'error:',element,r,n,N,R; kill
        
        zsd = ZScoreData(element,r,n,z,null_z,gene_symbols)
        if element_type == 'domain': original_domain_z_score_data[element] = zsd
        elif element_type == 'microRNA': original_microRNA_z_score_data[element] = zsd
        permuted_z_scores[element] = [z]
        if perform_element_permutation_analysis == 'no':
            ### The below is an alternative to the permute t-statistic that is more effecient
            p = FishersExactTest(r,n,R,N)
            zsd.SetP(p)
    return N,R

######### Begin Permutation Analysis #######
def calculatePermuteZScores(permute_element_inputs,element_denominator_gene_count,N,R):
    ###Make this code as efficient as possible
    for element_input_gene_count in permute_element_inputs:  
        for element in element_input_gene_count:
            r = element_input_gene_count[element]
            n = element_denominator_gene_count[element]
            try: z = statistics.zscore(r,n,N,R)
            except Exception: z = 0
            permuted_z_scores[element].append(abs(z))
            #if element == '0005488':
            #a.append(r)
                
def calculatePermuteStats(original_element_z_score_data):
    for element in original_element_z_score_data:
        zsd = original_element_z_score_data[element]
        z = abs(permuted_z_scores[element][0])
        permute_scores = permuted_z_scores[element][1:] ###Exclude the true value
        nullz = zsd.NullZ()
        if abs(nullz) == z: ###Only add the nullz values if they can count towards the p-value (if equal to the original z)
            null_z_to_add = permutations - len(permute_scores)
            permute_scores+=[abs(nullz)]*null_z_to_add ###Add null_z's in proportion to the amount of times there were not genes found for that element
        if len(permute_scores)>0:
            p = permute_p(permute_scores,z)
        else: p = 1
        #if p>1: p=1
        zsd.SetP(p)

def FishersExactTest(r,n,R,N):
    a = r; b = n-r; c=R-r; d=N-R-b
    table = [[int(a),int(b)], [int(c),int(d)]]
    
    try: ### Scipy version - cuts down rutime by ~1/3rd the time
        oddsratio, pvalue = stats.fisher_exact(table)
        return pvalue
    except Exception:
        ft = fishers_exact_test.FishersExactTest(table)
        return ft.two_tail_p()
        
def adjustPermuteStats(original_element_z_score_data):
    #1. Sort ascending the original input p value vector.  Call this spval.  Keep the original indecies so you can sort back.
    #2. Define a new vector called tmp.  tmp= spval.  tmp will contain the BH p values.
    #3. m is the length of tmp (also spval)
    #4. i=m-1
    #5  tmp[ i ]=min(tmp[i+1], min((m/i)*spval[ i ],1)) - second to last, last, last/second to last
    #6. i=m-2
    #7  tmp[ i ]=min(tmp[i+1], min((m/i)*spval[ i ],1))
    #8  repeat step 7 for m-3, m-4,... until i=1
    #9. sort tmp back to the original order of the input p values.
    
    spval=[]
    for element in original_element_z_score_data:
        zsd = original_element_z_score_data[element]
        p = float(zsd.PermuteP())
        spval.append([p,element])
        
    spval.sort(); tmp = spval; m = len(spval); i=m-2; x=0 ###Step 1-4
    
    while i > -1:
        tmp[i]=min(tmp[i+1][0], min((float(m)/(i+1))*spval[i][0],1)),tmp[i][1]; i -= 1
        
    for (adjp,element) in tmp:
        zsd = original_element_z_score_data[element]
        zsd.SetAdjP(adjp)
    spval=[]
    
def permute_p(null_list,true_value):
    y = 0; z = 0; x = permutations
    for value in null_list:
        if value >= true_value: y += 1
    #if true_value > 8: global a; a = null_list; print true_value,y,x;kill
    return (float(y)/float(x))  ###Multiply probabilty x2?

######### End Permutation Analysis #######
def exportZScoreData(original_element_z_score_data,element_type):
    element_output = root_dir+'AltResults/AlternativeOutput/' + dataset_name + analysis_method+'-'+element_type+'-zscores.txt'
    data = export.ExportFile(element_output)

    headers = [element_type+'-Name','Number Changed','Number Measured','Percent Changed', 'Zscore','PermuteP','AdjP','Changed GeneSymbols']
    headers = string.join(headers,'\t')+'\n'
    data.write(headers); sort_results=[]
    #print "Results for",len(original_element_z_score_data),"elements exported to",element_output
    for element in original_element_z_score_data:
        zsd=original_element_z_score_data[element]
        try: results = [zsd.Changed(), zsd.Measured(), zsd.PercentChanged(), zsd.ZScore(), zsd.PermuteP(), zsd.AdjP(), zsd.AssociatedWithElement()]
        except AttributeError: print element,len(permuted_z_scores[element]);kill
        results = [element] + results
        results = string.join(results,'\t') + '\n'
        sort_results.append([float(zsd.PermuteP()),-1/float(zsd.Measured()),results])

    sort_results.sort()
    for values in sort_results:
        results = values[2]
        data.write(results)
    data.close()

def getInputsForPermutationAnalysis(exon_db):
    ### Filter fold_dbase, which is the proper denominator
    probeset_to_gene = {}; denominator_list = []
    for probeset in exon_db:
        proceed = 'no'
        if filter_for_AS == 'yes':
            as_call = exon_db[probeset].SplicingCall()
            if as_call == 1: proceed = 'yes'
        else: proceed = 'yes'
        if proceed == 'yes':
            gene = exon_db[probeset].GeneID()
            probeset_to_gene[probeset] = gene
            denominator_list.append(probeset)
    return probeset_to_gene,denominator_list

def getJunctionSplicingAnnotations(regulated_exon_junction_db):
        filter_status = 'yes'
        ###########  Import critical exon annotation for junctions, build through the exon array analysis pipeline - link back to probesets
        filtered_arrayids={}; critical_probeset_annotation_db={}
        if array_type == 'RNASeq' and explicit_data_type == 'null':
            critical_exon_annotation_file = root_dir+'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_exons.txt'
        elif array_type == 'RNASeq' and explicit_data_type != 'null':
            critical_exon_annotation_file = root_dir+'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_junctions.txt'
        else:
            critical_exon_annotation_file = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_"+array_type+"_probesets.txt"
            critical_exon_annotation_file = filename=getFilteredFilename(critical_exon_annotation_file)
        for uid in regulated_exon_junction_db:
            gene = regulated_exon_junction_db[uid].GeneID()
            critical_exons = regulated_exon_junction_db[uid].CriticalExons()
            """### It appears that each critical exon for junction arrays can be a concatenation of multiple exons, making this unnecessary
            if len(critical_exons)>1 and array_type == 'junction':
                critical_exons_joined = string.join(critical_exons,'|')
                filtered_arrayids[gene+':'+critical_exon].append(uid)"""
            for critical_exon in critical_exons:
                try:
                    try: filtered_arrayids[gene+':'+critical_exon].append(uid)
                    except TypeError: print gene, critical_exon, uid;kill
                except KeyError: filtered_arrayids[gene+':'+critical_exon]=[uid]
        critical_exon_annotation_db = importSplicingAnnotationDatabase(critical_exon_annotation_file,'exon-fake',filtered_arrayids,filter_status);null=[] ###The file is in exon centric format, so designate array_type as exon
        for key in critical_exon_annotation_db:
            ced = critical_exon_annotation_db[key]
            for junction_probesets in filtered_arrayids[key]:
                try: critical_probeset_annotation_db[junction_probesets].append(ced) ###use for splicing and Exon annotations
                except KeyError: critical_probeset_annotation_db[junction_probesets] = [ced]
        for junction_probesets in critical_probeset_annotation_db:
            if len(critical_probeset_annotation_db[junction_probesets])>1: ###Thus multiple exons associated, must combine annotations
                exon_ids=[]; external_exonids=[]; exon_regions=[]; splicing_events=[]
                for ed in critical_probeset_annotation_db[junction_probesets]:
                    ensembl_gene_id = ed.GeneID();  transcript_cluster_id = ed.ExternalGeneID()
                    exon_ids.append(ed.ExonID()); external_exonids.append(ed.ExternalExonIDs());  exon_regions.append(ed.ExonRegionID()); se = string.split(ed.SplicingEvent(),'|')
                    for i in se: splicing_events.append(i)
                splicing_events = unique.unique(splicing_events) ###remove duplicate entries
                exon_id = string.join(exon_ids,'|'); external_exonid = string.join(external_exonids,'|'); exon_region = string.join(exon_regions,'|'); splicing_event = string.join(splicing_events,'|')
                probe_data = AffyExonSTData(ensembl_gene_id, exon_id, external_exonid, '', exon_region, splicing_event, '','')
                if array_type != 'RNASeq': probe_data.setTranscriptCluster(transcript_cluster_id)
                critical_probeset_annotation_db[junction_probesets] = probe_data
            else:
                critical_probeset_annotation_db[junction_probesets] = critical_probeset_annotation_db[junction_probesets][0]
        return critical_probeset_annotation_db

def determineExternalType(external_probeset_db):
    external_probeset_db2={}
    if 'TC' in external_probeset_db:
        temp_index={}; i=0; type = 'JETTA'
        for name in external_probeset_db['TC'][0]: temp_index[i]=i; i+=1
        if 'PS:norm_expr_fold_change' in temp_index: NI_fold_index = temp_index['PS:norm_expr_fold_change']
        if 'MADS:pv_1over2' in temp_index: MADS_p1_index = temp_index['MADS:pv_1over2']
        if 'MADS:pv_2over1' in temp_index: MADS_p2_index = temp_index['MADS:pv_2over1']
        if 'TC:expr_fold_change' in temp_index: MADS_p2_index = temp_index['MADS:pv_2over1']
        if 'PsId' in temp_index: ps_index = temp_index['PsId']
        for tc in external_probeset_db:
            for list in external_probeset_db[tc]:
                try: NI_fold = float(list[NI_fold_index])
                except Exception: NI_fold = 1
                try: MADSp1 = float(list[MADS_p1_index])
                except Exception: MADSp1 = 1
                try: MADSp2 = float(list[MADS_p2_index])
                except Exception: MADSp1 = 1
                if MADSp1<MADSp2: pval = MADSp1
                else: pval = MADSp2
                probeset = list[ps_index]
                external_probeset_db2[probeset] = NI_fold,pval
    else:
        type = 'generic'
        a = []; b = []
        for id in external_probeset_db:
            #print external_probeset_db[id]
            try: a.append(abs(float(external_probeset_db[id][0][0])))
            except Exception: null=[]
            try: b.append(abs(float(external_probeset_db[id][0][1])))
            except Exception: null=[]
        a.sort(); b.sort(); pval_index = None; score_index = None
        if len(a)>0:
            if max(a) > 1: score_index = 0
            else: pval_index = 0
        if len(b)>0:
            if max(b) > 1: score_index = 1
            else: pval_index = 1
        for id in external_probeset_db:
            if score_index != None: score = external_probeset_db[id][0][score_index]
            else: score = 1
            if pval_index != None: pval = external_probeset_db[id][0][pval_index]
            else: pval = 1
            external_probeset_db2[id] = score,pval
    return external_probeset_db2, type
                
def importExternalProbesetData(dataset_dir):
    excluded_probeset_db={}; splice_event_list=[]; p_value_call={}; permute_p_values={}; gene_expression_diff_db={}
    analyzed_probeset_db = {}
    external_probeset_db = importExternalDBList(dataset_dir)
    external_probeset_db, ext_type = determineExternalType(external_probeset_db)

    for probeset in exon_db: analyzed_probeset_db[probeset] = []
    ### Used to restrict the analysis to a pre-selected set of probesets (e.g. those that have a specifc splicing pattern)
    if len(filtered_probeset_db)>0:
        temp_db={}
        for probeset in analyzed_probeset_db: temp_db[probeset]=[]
        for probeset in temp_db:
            try: filtered_probeset_db[probeset]
            except KeyError: del analyzed_probeset_db[probeset]

    ### Used to restrict the analysis to a pre-selected set of probesets (e.g. those that have a specifc splicing annotation)
    if filter_for_AS == 'yes':
        for probeset in exon_db:
            as_call = exon_db[probeset].SplicingCall()
            if as_call == 0:
                try: del analyzed_probeset_db[probeset]
                except KeyError: null=[]
                
    for probeset in analyzed_probeset_db:
        ed = exon_db[probeset]; geneid = ed.GeneID()
        td = TranscriptionData('',''); gene_expression_diff_db[geneid] = td
        if probeset in external_probeset_db:
            exonid = ed.ExonID(); critical_exon_list = [1,[exonid]]
            splicing_index,normIntensityP = external_probeset_db[probeset]
            group1_ratios=[]; group2_ratios=[];exp_log_ratio=''; ttest_exp_p='';normIntensityP='';opposite_SI_log_mean=''
            sid = ExonData(splicing_index,probeset,critical_exon_list,geneid,group1_ratios,group2_ratios,normIntensityP,opposite_SI_log_mean)
            splice_event_list.append((splicing_index,sid))
        else:
            ### Also record the data for probesets that are excluded... Used by DomainGraph
            eed = ExcludedExonData(0,geneid,'NA')
            excluded_probeset_db[probeset] = eed        
    print len(splice_event_list), 'pre-filtered external results imported...\n'
    return splice_event_list, p_value_call, permute_p_values, excluded_probeset_db, gene_expression_diff_db

def splicingAnalysisAlgorithms(nonlog_NI_db,fold_dbase,dataset_name,gene_expression_diff_db,exon_db,ex_db,si_db,dataset_dir):
    protein_exon_feature_db={}; global regulated_exon_junction_db; global critical_exon_annotation_db; global probeset_comp_db; probeset_comp_db={}
    if original_conditions == 2: print "Beginning to run", analysis_method, "algorithm on",dataset_name[0:-1],"data"
    if run_from_scratch == 'Annotate External Results':
        splice_event_list, p_value_call, permute_p_values, excluded_probeset_db, gene_expression_diff_db = importExternalProbesetData(dataset_dir)
    elif analysis_method == 'ASPIRE' or analysis_method == 'linearregres':
        original_exon_db = exon_db
        if original_conditions > 2:
            splice_event_list, probeset_comp_db, permute_p_values, excluded_probeset_db = ex_db
            splice_event_list, p_value_call, permute_p_values, exon_db, regulated_exon_junction_db = furtherProcessJunctionScores(splice_event_list, probeset_comp_db, permute_p_values)
        else: 
            splice_event_list, probeset_comp_db, permute_p_values, excluded_probeset_db = analyzeJunctionSplicing(nonlog_NI_db)
            splice_event_list, p_value_call, permute_p_values, exon_db, regulated_exon_junction_db = furtherProcessJunctionScores(splice_event_list, probeset_comp_db, permute_p_values)
    elif analysis_method == 'splicing-index':
        regulated_exon_junction_db = {}
        if original_conditions > 2:
            excluded_probeset_db = ex_db; splice_event_list = si_db;
            clearObjectsFromMemory(ex_db); clearObjectsFromMemory(si_db)
            ex_db=[]; si_db=[]; permute_p_values={}; p_value_call=''
        else: splice_event_list, p_value_call, permute_p_values, excluded_probeset_db = analyzeSplicingIndex(fold_dbase)
    elif analysis_method == 'FIRMA':
        regulated_exon_junction_db = {}
        splice_event_list, p_value_call, permute_p_values, excluded_probeset_db = FIRMAanalysis(fold_dbase)
                
    global permuted_z_scores; permuted_z_scores={}; global original_domain_z_score_data; original_domain_z_score_data={}
    global original_microRNA_z_score_data; original_microRNA_z_score_data={}
    nonlog_NI_db=[] ### Clear memory of this large dictionary
    try: clearObjectsFromMemory(original_avg_const_exp_db); clearObjectsFromMemory(array_raw_group_values)
    except Exception: null=[]
    try: clearObjectsFromMemory(avg_const_exp_db)
    except Exception: null=[]
    try: clearObjectsFromMemory(alt_junction_db)
    except Exception: null=[]
    try: clearObjectsFromMemory(fold_dbase); fold_dbase=[]
    except Exception: null=[]

    microRNA_full_exon_db,microRNA_count_db,gene_microRNA_denom = ExonAnalyze_module.importmicroRNADataExon(species,array_type,exon_db,microRNA_prediction_method,explicit_data_type,root_dir)
    #print "MicroRNA data imported"

    if use_direct_domain_alignments_only == 'yes':
        protein_ft_db_len,domain_associated_genes = importProbesetAligningDomains(exon_db,'gene')
    else: protein_ft_db_len,domain_associated_genes = importProbesetProteinCompDomains(exon_db,'gene','exoncomp')

    if perform_element_permutation_analysis == 'yes':
        probeset_to_gene,denominator_list = getInputsForPermutationAnalysis(exon_db)

    if array_type == 'gene' or array_type == 'junction' or array_type == 'RNASeq':
        exon_gene_array_translation_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_'+array_type+'-exon_probesets.txt'
        try: exon_array_translation_db = importGeneric(exon_gene_array_translation_file)
        except Exception: exon_array_translation_db={} ### Not present for all species
        
    exon_hits={}; clearObjectsFromMemory(probeset_comp_db); probeset_comp_db=[]
    ###Run analyses in the ExonAnalyze_module module to assess functional changes
    for (score,ed) in splice_event_list:
        geneid = ed.GeneID()
        if analysis_method == 'ASPIRE' or 'linearregres' in analysis_method:
            pl = string.split(ed.Probeset1(),'|'); probeset1 = pl[0] ### When agglomerated, this is important
            uid = (probeset1,ed.Probeset2())
        else: uid = ed.Probeset1()
        gene_exon = geneid,uid; exon_hits[gene_exon] = ed
        #print probeset1,ed.Probeset1(),ed.Probeset2(),gene_exon,ed.CriticalExons()

    dataset_name_original = analysis_method+'-'+dataset_name[8:-1]
    global functional_attribute_db; global protein_features
    
    ### Possibly Block-out code for DomainGraph export             
    ###########  Re-import the exon_db for significant entries with full annotaitons
    exon_db={}; filtered_arrayids={}; filter_status='yes' ###Use this as a means to save memory (import multiple times - only storing different types relevant information)
    for (score,entry) in splice_event_list:
        try: probeset = original_exon_db[entry.Probeset1()].Probeset()
        except Exception: probeset = entry.Probeset1()
        pl = string.split(probeset,'|'); probeset = pl[0]; filtered_arrayids[probeset] = [] ### When agglomerated, this is important
        if array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null'):
            try: probeset = entry.Probeset2(); filtered_arrayids[probeset] = []
            except AttributeError: null =[] ###occurs when running Splicing 
    exon_db = importSplicingAnnotationDatabase(probeset_annotations_file,array_type,filtered_arrayids,filter_status);null=[] ###replace existing exon_db (probeset_annotations_file should be a global)

    ###domain_gene_changed_count_db is the number of genes for each domain that are found for regulated probesets
    if array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null'):
        if use_direct_domain_alignments_only == 'yes':
            protein_features,domain_gene_changed_count_db,functional_attribute_db = importProbesetAligningDomains(regulated_exon_junction_db,'probeset')
        else: protein_features,domain_gene_changed_count_db,functional_attribute_db = importProbesetProteinCompDomains(regulated_exon_junction_db,'probeset','exoncomp')
    else:
        if use_direct_domain_alignments_only == 'yes':
            protein_features,domain_gene_changed_count_db,functional_attribute_db = importProbesetAligningDomains(exon_db,'probeset')
        else: protein_features,domain_gene_changed_count_db,functional_attribute_db = importProbesetProteinCompDomains(exon_db,'probeset','exoncomp')

    filtered_microRNA_exon_db = ExonAnalyze_module.filterMicroRNAProbesetAssociations(microRNA_full_exon_db,exon_hits)
    microRNA_full_exon_db=[]
    ###add microRNA data to functional_attribute_db
    microRNA_hit_gene_count_db = {}; all_microRNA_gene_hits={}; microRNA_attribute_db={}; probeset_mirBS_db={}
    for (affygene,uid) in filtered_microRNA_exon_db: ###example ('G7091354', 'E20|') [('hsa-miR-130a', 'Pbxip1'), ('hsa-miR-130a', 'Pbxip1'
        ###3-1-08
        miR_list = []
        microRNA_symbol_list = filtered_microRNA_exon_db[(affygene,uid)]
        for mir_key in microRNA_symbol_list:
            microRNA,gene_symbol,miR_seq, miR_sources = mir_key
            #if 'ENS' in microRNA: print microRNA; kill  ### bug in some miRNA annotations introduced in the build process
            specific_microRNA_tuple = (microRNA,'~')
            try: microRNA_hit_gene_count_db[microRNA].append(affygene)
            except KeyError: microRNA_hit_gene_count_db[microRNA] = [affygene]
            ###Create a database with the same structure as "protein_exon_feature_db"(below) for over-representation analysis (direction specific), after linking up splice direction data
            try: microRNA_attribute_db[(affygene,uid)].append(specific_microRNA_tuple)
            except KeyError: microRNA_attribute_db[(affygene,uid)] = [specific_microRNA_tuple]
            miR_data = microRNA+':'+miR_sources
            miR_list.append(miR_data) ###Add miR information to the record
            function_type = ('miR-sequence: ' +'('+miR_data+')'+miR_seq,'~') ###Add miR sequence information to the sequence field of the report
            try: functional_attribute_db[(affygene,uid)].append(function_type)
            except KeyError: functional_attribute_db[(affygene,uid)]=[function_type]
            #print (affygene,uid), [function_type];kill
            if perform_element_permutation_analysis == 'yes':
                try: probeset_mirBS_db[uid].append(microRNA)
                except KeyError: probeset_mirBS_db[uid] = [microRNA]
                
        miR_str = string.join(miR_list,','); miR_str = '('+miR_str+')'
        function_type = ('microRNA-target'+miR_str,'~')
        try: functional_attribute_db[(affygene,uid)].append(function_type)
        except KeyError: functional_attribute_db[(affygene,uid)]=[function_type]
        all_microRNA_gene_hits[affygene] = []
            
    ###Replace the gene list for each microRNA hit with count data            
    microRNA_hit_gene_count_db = eliminate_redundant_dict_values(microRNA_hit_gene_count_db)
            
    ###Combines any additional feature alignment info identified from 'ExonAnalyze_module.characterizeProteinLevelExonChanges' (e.g. from Ensembl or junction-based queries rather than exon specific) and combines
    ###this with this database of (Gene,Exon)=[(functional element 1,'~'),(functional element 2,'~')] for downstream result file annotatations
    domain_hit_gene_count_db = {}; all_domain_gene_hits = {}; probeset_domain_db={}
    for entry in protein_features:
        gene,uid = entry
        for data_tuple in protein_features[entry]:
            domain,call = data_tuple
            try: protein_exon_feature_db[entry].append(data_tuple)
            except KeyError: protein_exon_feature_db[entry] = [data_tuple]
            try: domain_hit_gene_count_db[domain].append(gene)
            except KeyError: domain_hit_gene_count_db[domain] = [gene]
            all_domain_gene_hits[gene]=[]
            if perform_element_permutation_analysis == 'yes':
                try: probeset_domain_db[uid].append(domain)
                except KeyError: probeset_domain_db[uid] = [domain]
    protein_features=[]; domain_gene_changed_count_db=[]
    ###Replace the gene list for each microRNA hit with count data
    domain_hit_gene_count_db = eliminate_redundant_dict_values(domain_hit_gene_count_db)

    ############   Perform Element Over-Representation Analysis ############
    """Domain/FT Fishers-Exact test: with "protein_exon_feature_db" (transformed to "domain_hit_gene_count_db") we can analyze over-representation of domain/features WITHOUT taking into account exon-inclusion or exclusion
    Do this using: "domain_associated_genes", which contains domain tuple ('Tyr_pkinase', 'IPR001245') as a key and count in unique genes as the value in addition to
    Number of genes linked to splice events "regulated" (SI and Midas p<0.05), number of genes with constitutive probesets

    MicroRNA Fishers-Exact test: "filtered_microRNA_exon_db" contains gene/exon to microRNA data. For each microRNA, count the representation in spliced genes microRNA (unique gene count - make this from the mentioned file)
    Do this using: "microRNA_count_db"""

    domain_gene_counts = {} ### Get unique gene counts for each domain
    for domain in domain_associated_genes:
        domain_gene_counts[domain] = len(domain_associated_genes[domain])
                
    total_microRNA_gene_hit_count = len(all_microRNA_gene_hits)
    total_microRNA_gene_denom_count = len(gene_microRNA_denom)
    Nm,Rm = calculateZScores(microRNA_hit_gene_count_db,microRNA_count_db,total_microRNA_gene_denom_count,total_microRNA_gene_hit_count,'microRNA')

    gene_microRNA_denom =[]
    
    summary_data_db['miRNA_gene_denom'] = total_microRNA_gene_denom_count
    summary_data_db['miRNA_gene_hits'] = total_microRNA_gene_hit_count
    summary_data_db['alt_events']=len(splice_event_list)
    
    total_domain_gene_hit_count = len(all_domain_gene_hits)
    total_domain_gene_denom_count = protein_ft_db_len ###genes connected to domain annotations
    Nd,Rd = calculateZScores(domain_hit_gene_count_db,domain_gene_counts,total_domain_gene_denom_count,total_domain_gene_hit_count,'domain')
        
    microRNA_hit_gene_counts={}; gene_to_miR_db={} ### Get unique gene counts for each miR and the converse
    for microRNA in microRNA_hit_gene_count_db:
        microRNA_hit_gene_counts[microRNA] = len(microRNA_hit_gene_count_db[microRNA])
        for gene in microRNA_hit_gene_count_db[microRNA]:
            try: gene_to_miR_db[gene].append(microRNA)
            except KeyError: gene_to_miR_db[gene] = [microRNA]
    gene_to_miR_db = eliminate_redundant_dict_values(gene_to_miR_db)

    if perform_element_permutation_analysis == 'yes':
        ###Begin Domain/microRNA Permute Analysis
        input_count = len(splice_event_list)  ### Number of probesets or probeset pairs (junction array) alternatively regulated
        original_increment = int(permutations/20); increment = original_increment
        start_time = time.time(); print 'Permuting the Domain/miRBS analysis %d times' % permutations
        x=0; permute_domain_inputs=[]; permute_miR_inputs=[]
        while x<permutations:
            if x == increment: increment+=original_increment; print '*',
            permute_input_list = random.sample(denominator_list,input_count); x+=1
            permute_domain_input_gene_counts = countGenesForElement(permute_input_list,probeset_to_gene,probeset_domain_db)
            permute_domain_inputs.append(permute_domain_input_gene_counts)
            permute_miR_input_gene_counts = countGenesForElement(permute_input_list,probeset_to_gene,probeset_mirBS_db)
            permute_miR_inputs.append(permute_miR_input_gene_counts)
        calculatePermuteZScores(permute_domain_inputs,domain_gene_counts,Nd,Rd)
        calculatePermuteZScores(permute_miR_inputs,microRNA_hit_gene_counts,Nm,Rm)
        calculatePermuteStats(original_domain_z_score_data)
        calculatePermuteStats(original_microRNA_z_score_data)
        adjustPermuteStats(original_domain_z_score_data)
        adjustPermuteStats(original_microRNA_z_score_data)
        exportZScoreData(original_domain_z_score_data,'ft-domain')
        exportZScoreData(original_microRNA_z_score_data,'microRNA')
        end_time = time.time(); time_diff = int(end_time-start_time)
        print "Enrichment p-values for Domains/miRBS calculated in %d seconds" % time_diff

    denominator_list=[]
    try: clearObjectsFromMemory(original_microRNA_z_score_data)
    except Exception: null=[]
    microRNA_hit_gene_count_db={}; microRNA_hit_gene_counts={};
    clearObjectsFromMemory(permuted_z_scores); permuted_z_scores=[]; original_domain_z_score_data=[]
    
    if (array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null')) and analysis_method != 'splicing-index':
        critical_probeset_annotation_db = getJunctionSplicingAnnotations(regulated_exon_junction_db)
        probeset_aligning_db = importProbesetAligningDomains(regulated_exon_junction_db,'perfect_match')
        
    else: probeset_aligning_db = importProbesetAligningDomains(exon_db,'perfect_match')
    
    ############   Export exon/junction level results ############    
    splice_event_db={}; protein_length_list=[]; aspire_gene_results={}
    critical_gene_exons={}; unique_exon_event_db={}; comparison_count={}; direct_domain_gene_alignments={}
    functional_attribute_db2={}; protein_exon_feature_db2={}; microRNA_exon_feature_db2={}
    external_exon_annot={}; gene_exon_region={}; gene_smallest_p={}; gene_splice_event_score={}; alternatively_reg_tc={}
    aspire_output = root_dir+'AltResults/AlternativeOutput/' + dataset_name + analysis_method+'-exon-inclusion-results.txt'
    data = export.ExportFile(aspire_output)
    goelite_output = root_dir+'GO-Elite/AltExon/AS.'+ dataset_name + analysis_method+'.txt'
    goelite_data = export.ExportFile(goelite_output); gcn=0

    #print 'LENGTH OF THE GENE ANNOTATION DATABASE',len(annotate_db)    
    if array_type != 'AltMouse':
        DG_output = root_dir+'AltResults/DomainGraph/' + dataset_name + analysis_method+'-DomainGraph.txt'
        DG_data = export.ExportFile(DG_output)
        
        ### Write out only the inclusion hits to a subdir
        SRFinder_inclusion = root_dir+'GO-Elite/exon/' + dataset_name + analysis_method+'-inclusion.txt'
        SRFinder_in_data = export.ExportFile(SRFinder_inclusion)
        SRFinder_in_data.write('probeset\tSystemCode\tdeltaI\tp-value\n')
    
        ### Write out only the exclusion hits to a subdir
        SRFinder_exclusion = root_dir+'GO-Elite/exon/' + dataset_name + analysis_method+'-exclusion.txt'
        SRFinder_ex_data = export.ExportFile(SRFinder_exclusion)
        SRFinder_ex_data.write('probeset\tSystemCode\tdeltaI\tp-value\n')

        ### Write out only the denominator set to a subdir
        SRFinder_denom = root_dir+'GO-Elite/exon_denominator/' + species+'-'+array_type+'.txt'
        SRFinder_denom_data = export.ExportFile(SRFinder_denom)
        SRFinder_denom_data.write('probeset\tSystemCode\n')
        
        ens_version = unique.getCurrentGeneDatabaseVersion()
        ProcessedSpliceData_output = string.replace(DG_output,'DomainGraph','ProcessedSpliceData') ### This is the same as the DG export but without converting the probeset IDs for non-exon arrays
        ProcessedSpliceData_data = export.ExportFile(ProcessedSpliceData_output)
        if ens_version == '':
            try:
                elite_db_versions = UI.returnDirectoriesNoReplace('/AltDatabase')
                if len(elite_db_versions)>0: ens_version = elite_db_versions[0]
            except Exception: null=[]
        ens_version = string.replace(ens_version,'EnsMart','ENS_')
        DG_data.write(ens_version+"\n")
        DG_data.write("Probeset\tGeneID\tRegulation call\tSI\tSI p-value\tMiDAS p-value\n")
        ProcessedSpliceData_data.write("ExonID(s)\tGeneID\tRegulation call\t"+analysis_method+"\t"+analysis_method+" p-value\tMiDAS p-value\n")
             
    if analysis_method == 'ASPIRE' or 'linearregres' in analysis_method:
        if perform_permutation_analysis == 'yes': p_value_type = 'permutation-values'
        else: p_value_type = 'FDR-'+p_value_call
        if array_type == 'AltMouse': gene_name = 'AffyGene'; extra_transcript_annotation = 'block_structure'; extra_exon_annotation = 'splice_event_description'
        if array_type == 'junction' or array_type == 'RNASeq':
            gene_name = 'Ensembl'; extra_transcript_annotation = 'transcript cluster ID'; extra_exon_annotation = 'distal exon-region-ID'
            goelite_data.write("GeneID\tSystemCode\tscore\tp-value\tSymbol\tExonIDs\n")
        if array_type == 'RNASeq':
            id1='junctionID-1'; id2='junctionID-2'; loc_column='exon/junction locations'
            extra_transcript_annotation = 'Known/Novel Feature'
        else: id1='probeset1'; id2='probeset2'; loc_column='probeset locations'
        title = [gene_name,analysis_method,'symbol','description','exons1','exons2','regulation_call','event_call',id1,'norm-p1',id2,'norm-p2','fold1','fold2']
        title +=['adj-fold1' ,'adj-fold2' ,extra_transcript_annotation,'critical_up_exons','critical_down_exons','functional_prediction','uniprot-ens_feature_predictions']
        title +=['peptide_predictions','exp1','exp2','ens_overlapping_domains','constitutive_baseline_exp',p_value_call,p_value_type,'permutation-false-positives']
        title +=['gene-expression-change', extra_exon_annotation ,'ExternalExonIDs','ExonRegionID','SplicingEvent','ExonAnnotationScore','large_splicing_diff',loc_column]
    else:
        goelite_data.write("GeneID\tSystemCode\tSI\tSI p-value\tMiDAS p-value\tSymbol\tExonID\n")
            
        if analysis_method == 'splicing-index':
            NIpval = 'SI_rawp'; splicing_score = 'Splicing-Index'; lowestp = 'lowest_p (MIDAS or SI)'; AdjPcolumn = 'Deviation-Value'; #AdjPcolumn = 'SI_adjp'
        else:
            NIpval = 'FIRMA_rawp'; splicing_score = 'FIRMA_fold'; lowestp = 'lowest_p (MIDAS or FIRMA)'; AdjPcolumn = 'Deviation-Value'; #AdjPcolumn = 'FIRMA_adjp'
        if array_type == 'RNASeq':
            id1='junctionID'; pval_column='junction p-value'; loc_column='junction location'
        else: id1='probeset'; pval_column='probeset p-value'; loc_column='probeset location'
        
        if array_type == 'RNASeq': secondary_ID_title = 'Known/Novel Feature'
        else: secondary_ID_title = 'alternative gene ID'
        
        title= ['Ensembl',splicing_score,'symbol','description','exons','regulation_call',id1,pval_column,lowestp,'midas p-value','fold','adjfold']
        title+=['up_exons','down_exons','functional_prediction','uniprot-ens_feature_predictions','peptide_predictions','ens_overlapping_domains','baseline_probeset_exp']
        title+=['constitutive_baseline_exp',NIpval,AdjPcolumn,'gene-expression-change']
        title+=[secondary_ID_title, 'ensembl exons', 'consitutive exon', 'exon-region-ID', 'exon annotations','distal exon-region-ID',loc_column]    
    title = string.join(title,'\t') + '\n'
    try:
        if original_conditions>2: title = string.replace(title,'regulation_call','conditions_compared')
    except Exception: null=[]
    data.write(title)

    ### Calculate adjusted normalized intensity p-values
    fdr_exon_stats={}
    if analysis_method != 'ASPIRE' and 'linearregres' not in analysis_method:    
        for (score,entry) in splice_event_list: ### These are all "significant entries"
            fds = FDRStats(entry.TTestNormalizedRatios())
            fdr_exon_stats[entry.Probeset1()] = fds
        for probeset in excluded_probeset_db: ### These are all "non-significant entries"
            fds = FDRStats(excluded_probeset_db[probeset].TTestNormalizedRatios())
            fdr_exon_stats[probeset] = fds        
    try: adjustPermuteStats(fdr_exon_stats)
    except Exception: null=[]
    
    ### Calculate score average and stdev for each gene to alter get a Deviation Value
    gene_deviation_db={}
    for (score,entry) in splice_event_list:
        dI = entry.Score(); geneID = entry.GeneID()
        try: gene_deviation_db[geneID].append(dI)
        except Exception: gene_deviation_db[geneID] = [dI]
    for i in excluded_probeset_db:
        entry = excluded_probeset_db[i]
        try: dI = entry.Score(); geneID = entry.GeneID()
        except Exception: geneID = entry[1]; dI = entry[-1]
        try: gene_deviation_db[geneID].append(dI)
        except Exception: None ### Don't include genes with no hits
    for geneID in gene_deviation_db:
        try:
            avg_dI=statistics.avg(gene_deviation_db[geneID])
            stdev_dI=statistics.stdev(gene_deviation_db[geneID])
            gene_deviation_db[geneID] = avg_dI,stdev_dI
        except Exception: 
            gene_deviation_db[geneID] = 'NA','NA'
    
    event_count = 0
    for (score,entry) in splice_event_list:
        event_count += 1
        dI = entry.Score(); probeset1 = entry.Probeset1(); regulation_call = entry.RegulationCall(); event_call = entry.EventCall();critical_exon_list = entry.CriticalExonTuple()
        probeset1_display = probeset1; selected_probeset = probeset1
        if agglomerate_inclusion_probesets == 'yes':
            if array_type == 'AltMouse':
                exons1 = original_exon_db[probeset1].ExonID()
                try: probeset1 = original_exon_db[probeset1].Probeset()
                except Exception: null=[]
            else:
                probeset1 = probeset1; exons1 = original_exon_db[probeset1].ExonID()
            try: selected_probeset = original_exon_db[probeset1].Probeset()
            except Exception: selected_probeset = probeset1
        else:
            try: exons1 = exon_db[probeset1].ExonID()
            except Exception:
                print probeset1, len(exon_db)
                for i in exon_db: print i; break
                kill
        critical_probeset_list = [selected_probeset]
        affygene = entry.GeneID()

        ### Calculate deviation value for each exon
        avg_dI,stdev_dI = gene_deviation_db[affygene]
        try: DV = deviation(dI,avg_dI,stdev_dI) ### Note: the dI values are always in log2 space, independent of platform
        except Exception: DV = 'NA'
            
        if affygene in annotate_db: description = annotate_db[affygene].Description(); symbol = annotate_db[affygene].Symbol()
        else: description = ''; symbol = ''
        ped1 = entry.ProbesetExprData1(); adjfold1 = ped1.AdjFold(); exp1 = ped1.BaselineExp(); fold1 = ped1.FoldChange(); rawp1 = ped1.ExpPval()

        ### Get Constitutive expression values
        baseline_const_exp = entry.ConstitutiveExpression()  ### For multiple group comparisosn
            
        #if affygene in gene_expression_diff_db: mean_fold_change = gene_expression_diff_db[affygene].ConstitutiveFoldStr()
        try: mean_fold_change = str(entry.ConstitutiveFold()) ### For multi-condition analyses, the gene expression is dependent on the conditions compared
        except Exception: mean_fold_change = gene_expression_diff_db[affygene].ConstitutiveFoldStr()
        
        if analysis_method == 'ASPIRE' or 'linearregres' in analysis_method:
            probeset2 = entry.Probeset2(); exons2 = exon_db[probeset2].ExonID(); rawp1 = str(entry.TTestNormalizedRatios()); rawp2 = str(entry.TTestNormalizedRatios2()); critical_probeset_list.append(probeset2)
            ped2 = entry.ProbesetExprData2(); adjfold2 = ped2.AdjFold(); exp2 = ped2.BaselineExp(); fold2 = ped2.FoldChange()
            try: location_summary=original_exon_db[selected_probeset].LocationSummary()+'|'+original_exon_db[probeset2].LocationSummary()
            except Exception:
                try: location_summary=exon_db[selected_probeset].LocationSummary()+'|'+exon_db[probeset2].LocationSummary()
                except Exception: location_summary=''    
            if array_type == 'AltMouse':
                extra_transcript_annotation = exon_db[probeset1].GeneStructure()
            else:
                try: extra_exon_annotation = last_exon_region_db[affygene]
                except KeyError: extra_exon_annotation = ''
                try:
                    tc1 = original_exon_db[probeset1].SecondaryGeneID()
                    tc2 = original_exon_db[probeset2].SecondaryGeneID() ### Transcript Cluster
                    probeset_tc = makeUnique([tc1,tc2])
                    extra_transcript_annotation = string.join(probeset_tc,'|')
                    try: alternatively_reg_tc[affygene] += probeset_tc
                    except KeyError: alternatively_reg_tc[affygene] = probeset_tc
                except Exception: extra_transcript_annotation=''
                if array_type == 'RNASeq':
                    try: extra_transcript_annotation = entry.NovelEvent() ### Instead of secondary gene ID, list known vs. novel reciprocal junction annotation
                    except Exception: None
            
            exp_list = [float(exp1),float(exp2),float(exp1)+float(fold1),float(exp2)+float(fold2)]; exp_list.sort();  exp_list.reverse()
            probeset_tuple = (probeset1,probeset2)
        else:
            try: exp_list = [float(exp1),float(exp1)+float(fold1)]; exp_list.sort();  exp_list.reverse()
            except Exception: exp_list = ['']
            probeset_tuple = (probeset1)
        highest_exp = exp_list[0]
        
        ###Use permuted p-value or lowest expression junction p-value based on the situtation
        ###This p-value is used to filter out aspire events for further analyses
        if len(p_value_call)>0:
            if probeset_tuple in permute_p_values:
                lowest_raw_p, pos_permute, total_permute, false_pos = permute_p_values[probeset_tuple]
            else: lowest_raw_p = "NA"; pos_permute = "NA"; total_permute = "NA"; false_pos = "NA"
        else:
            if analysis_method == 'ASPIRE' or 'linearregres' in analysis_method: raw_p_list = [entry.TTestNormalizedRatios(),entry.TTestNormalizedRatios2()]  #raw_p_list = [float(rawp1),float(rawp2)]; raw_p_list.sort()
            else:
                try: raw_p_list = [float(entry.TTestNormalizedRatios())]  ###Could also be rawp1, but this is more appropriate
                except Exception: raw_p_list = [1] ### Occurs when p='NA'
            raw_p_list.sort()
            lowest_raw_p = raw_p_list[0]; pos_permute = "NA"; total_permute = "NA"; false_pos = "NA"
        
        if perform_permutation_analysis == 'yes':
            p_value_extra = str(pos_permute)+' out of '+str(total_permute)
        else: p_value_extra = str(pos_permute)
        up_exons = ''; down_exons = ''; up_exon_list = []; down_exon_list = []; gene_exon_list=[]
        exon_data = critical_exon_list
        variable = exon_data[0]
        if variable == 1 and regulation_call == 'upregulated':
            for exon in exon_data[1]:
                up_exons = up_exons + exon + ',';up_exon_list.append(exon)
                key = affygene,exon+'|'; gene_exon_list.append(key)
        elif variable == 1 and regulation_call == 'downregulated':
            for exon in exon_data[1]:
                down_exons = down_exons + exon + ',';down_exon_list.append(exon)
                key = affygene,exon+'|';gene_exon_list.append(key)
        else:
            try: exon1 = exon_data[1][0]; exon2 = exon_data[1][1]
            except Exception: print exon_data;kill
            if adjfold1 > 0:
                up_exons = up_exons + exon1 + ',';down_exons = down_exons + exon2 + ','
                up_exon_list.append(exon1); down_exon_list.append(exon2)
                key = affygene,exon1+'|'; gene_exon_list.append(key);key = affygene,exon2+'|'; gene_exon_list.append(key)
            else:
                up_exons = up_exons + exon2 + ',';down_exons = down_exons + exon1 + ','
                up_exon_list.append(exon2); down_exon_list.append(exon1)
                key = affygene,exon1+'|'; gene_exon_list.append(key); key = affygene,exon2+'|'; gene_exon_list.append(key)
        up_exons = up_exons[0:-1];down_exons = down_exons[0:-1]

        try: ### Get comparisons group annotation data for multigroup comparison analyses
            if original_conditions>2:
                try: regulation_call = ped1.Annotation()
                except Exception: null=[]
        except Exception: null=[]
        
        ###Format functional results based on exon level fold change
        null = []
        #global a; a = exon_hits; global b; b=microRNA_attribute_db; kill
        """if 'G7100684@J934332_RC@j_at' in critical_probeset_list:
            print probeset1, probeset2, gene, critical_probeset_list, 'blah'
            if ('G7100684', ('G7100684@J934333_RC@j_at', 'G7100684@J934332_RC@j_at')) in functional_attribute_db:
                print functional_attribute_db[('G7100684', ('G7100684@J934333_RC@j_at', 'G7100684@J934332_RC@j_at'))];blah
            blah"""

        new_functional_attribute_str, functional_attribute_list2, seq_attribute_str,protein_length_list = format_exon_functional_attributes(affygene,critical_probeset_list,functional_attribute_db,up_exon_list,down_exon_list,protein_length_list)
        new_uniprot_exon_feature_str, uniprot_exon_feature_list, null, null = format_exon_functional_attributes(affygene,critical_probeset_list,protein_exon_feature_db,up_exon_list,down_exon_list,null)

        null, microRNA_exon_feature_list, null, null = format_exon_functional_attributes(affygene,critical_probeset_list,microRNA_attribute_db,up_exon_list,down_exon_list,null)
        if len(new_functional_attribute_str) == 0: new_functional_attribute_str = ' '
        if len(new_uniprot_exon_feature_str) == 0: new_uniprot_exon_feature_str = ' '
        if len(seq_attribute_str) > 12000: seq_attribute_str = 'The sequence is too long to report for spreadsheet analysis'
        ### Add entries to a database to quantify the number of reciprocal isoforms regulated
        reciprocal_isoform_data = [len(critical_exon_list[1]),critical_exon_list[1],event_call,regulation_call]
        try: float((lowest_raw_p))
        except ValueError: lowest_raw_p=0
        if (float((lowest_raw_p))<=p_threshold or false_pos < 2) or lowest_raw_p == 1 or lowest_raw_p == 'NA':
            try: unique_exon_event_db[affygene].append(reciprocal_isoform_data)
            except KeyError: unique_exon_event_db[affygene] = [reciprocal_isoform_data]
            
        ### Add functional attribute information to a new database
        for item in uniprot_exon_feature_list:
            attribute = item[0]
            exon = item[1]
            if (float((lowest_raw_p))<=p_threshold or false_pos < 2) or lowest_raw_p == 1 or lowest_raw_p == 'NA':
              try: protein_exon_feature_db2[affygene,attribute].append(exon)
              except KeyError: protein_exon_feature_db2[affygene,attribute]=[exon]
        ### Add functional attribute information to a new database
        """Database not used for exon/junction data export but for over-representation analysis (direction specific)"""
        for item in microRNA_exon_feature_list:
            attribute = item[0]
            exon = item[1]
            if (float((lowest_raw_p))<=p_threshold or false_pos < 2) or lowest_raw_p == 1 or lowest_raw_p == 'NA':
              try: microRNA_exon_feature_db2[affygene,attribute].append(exon)
              except KeyError: microRNA_exon_feature_db2[affygene,attribute]=[exon]              
        ### Add functional attribute information to a new database
        for item in functional_attribute_list2:
            attribute = item[0]
            exon = item[1]
            if (float((lowest_raw_p))<=p_threshold or false_pos < 2) or lowest_raw_p == 1 or lowest_raw_p == 'NA':
              try: functional_attribute_db2[affygene,attribute].append(exon)
              except KeyError: functional_attribute_db2[affygene,attribute]=[exon]
                    
        try:
            abs_fold = abs(float(mean_fold_change)); fold_direction = 'down'; fold1_direction = 'down'; fold2_direction = 'down'
            large_splicing_diff1 = 0; large_splicing_diff2 = 0; large_splicing_diff = 'null'; opposite_splicing_pattern = 'no'
            if float(mean_fold_change)>0: fold_direction = 'up'
            if float(fold1)>0: fold1_direction = 'up'
            if fold1_direction != fold_direction:
                if float(fold1)>float(mean_fold_change): large_splicing_diff1 = float(fold1)-float(mean_fold_change)
        except Exception:
            fold_direction = ''; large_splicing_diff = ''; opposite_splicing_pattern = ''
          
        if analysis_method != 'ASPIRE' and 'linearregres' not in analysis_method: ed = exon_db[probeset1]
        else:
            try: ed = critical_probeset_annotation_db[selected_probeset,probeset2]
            except KeyError:
                try: ed = exon_db[selected_probeset] ###not useful data here, but the objects need to exist
                except IOError: ed = original_exon_db[probeset1]
        ucsc_splice_annotations = ["retainedIntron","cassetteExon","strangeSplice","altFivePrime","altThreePrime","altPromoter","bleedingExon"]
        custom_annotations = ["alt-3'","alt-5'","alt-C-term","alt-N-term","cassette-exon","cassette-exon","exon-region-exclusion","intron-retention","mutually-exclusive-exon","trans-splicing"]

        custom_exon_annotations_found='no'; ucsc_annotations_found = 'no'; exon_annot_score=0
        if len(ed.SplicingEvent())>0: 
            for annotation in ucsc_splice_annotations:
                if annotation in ed.SplicingEvent(): ucsc_annotations_found = 'yes'
            for annotation in custom_annotations:
                if annotation in ed.SplicingEvent(): custom_exon_annotations_found = 'yes'

        if custom_exon_annotations_found == 'yes' and ucsc_annotations_found == 'no': exon_annot_score = 3
        elif ucsc_annotations_found == 'yes' and custom_exon_annotations_found == 'no': exon_annot_score = 4
        elif ucsc_annotations_found == 'yes' and custom_exon_annotations_found == 'yes': exon_annot_score = 5
        else: exon_annot_score = 2 
        try: gene_splice_event_score[affygene].append(exon_annot_score) ###store for gene level results
        except KeyError: gene_splice_event_score[affygene] = [exon_annot_score]
        try: gene_exon_region[affygene].append(ed.ExonRegionID()) ###store for gene level results
        except KeyError: gene_exon_region[affygene] = [ed.ExonRegionID()]          
                    
        if analysis_method == 'ASPIRE' or 'linearregres' in analysis_method:
            if float(fold2)>0: fold2_direction = 'up'
            if fold2_direction != fold_direction:
                if float(fold2)>float(mean_fold_change):
                    large_splicing_diff2 = float(fold2)-float(mean_fold_change)
                    if abs(large_splicing_diff2) > large_splicing_diff1: large_splicing_diff = str(large_splicing_diff2)
                    else: large_splicing_diff = str(large_splicing_diff1)
            if fold1_direction != fold2_direction and abs(float(fold1))>0.4 and abs(float(fold2))>0.4 and abs(float(mean_fold_change))< max([float(fold2),float(fold1)]):
                opposite_splicing_pattern = 'yes'                                                                                         

            ### Annotate splicing events based on exon_strucuture data
            if array_type == 'AltMouse':
                extra_exon_annotation = ExonAnnotate_module.annotate_splice_event(exons1,exons2,extra_transcript_annotation)
                try: splice_event_db[extra_exon_annotation] += 1
                except KeyError: splice_event_db[extra_exon_annotation] = 1

            try:
                direct_domain_alignments = probeset_aligning_db[selected_probeset,probeset2]
                try: direct_domain_gene_alignments[affygene]+=', '+direct_domain_alignments
                except KeyError: direct_domain_gene_alignments[affygene]=direct_domain_alignments
            except KeyError: direct_domain_alignments = ' '
            
            splicing_event = ed.SplicingEvent()
            if array_type == 'RNASeq':
                splicing_event = checkForTransSplicing(probeset1_display,splicing_event)
                splicing_event = checkForTransSplicing(probeset2,splicing_event)
                exp1 = covertLogExpressionToNonLog(exp1)
                exp2 = covertLogExpressionToNonLog(exp2)
                baseline_const_exp = covertLogExpressionToNonLog(baseline_const_exp)
                fold1 = covertLogFoldToNonLog(fold1)
                fold2 = covertLogFoldToNonLog(fold2)

                adjfold1 = covertLogFoldToNonLog(adjfold1)
                adjfold2 = covertLogFoldToNonLog(adjfold2)
                mean_fold_change = covertLogFoldToNonLog(mean_fold_change)

            ### Annotate splicing events based on pre-computed and existing annotations  
            values= [affygene,dI,symbol,fs(description),exons1,exons2,regulation_call,event_call,probeset1_display,rawp1,probeset2,rawp2,fold1,fold2,adjfold1,adjfold2]
            values+=[extra_transcript_annotation,up_exons,down_exons,fs(new_functional_attribute_str),fs(new_uniprot_exon_feature_str),fs(seq_attribute_str),exp1,exp2,fs(direct_domain_alignments)]
            values+=[str(baseline_const_exp),str(lowest_raw_p),p_value_extra,str(false_pos),mean_fold_change,extra_exon_annotation]
            values+=[ed.ExternalExonIDs(),ed.ExonRegionID(),splicing_event,str(exon_annot_score),large_splicing_diff,location_summary]
            exon_sets = abs(float(dI)),regulation_call,event_call,exons1,exons2,''

            ### Export significant reciprocol junction pairs and scores
            values_ps = [probeset1+'|'+probeset2,affygene,'changed',dI,'NA',str(lowest_raw_p)]; values_ps = string.join(values_ps,'\t')+'\n'
            try: ProcessedSpliceData_data.write(values_ps)
            except Exception: None
            
            values_ge = [affygene,'En',dI,str(lowest_raw_p),symbol,probeset1_display+' | '+probeset2]; values_ge = string.join(values_ge,'\t')+'\n'
            if array_type == 'junction' or array_type == 'RNASeq':   ### Only applies to reciprocal junction sensitive platforms (but not currently AltMouse)
                goelite_data.write(values_ge)

            if array_type == 'junction' or array_type == 'RNASeq': ### Only applies to reciprocal junction sensitive platforms (but not currently AltMouse)
                try: exon_probeset = exon_array_translation_db[affygene+':'+exon_data[1][0]][0]; probeset1 = exon_probeset; gcn+=1
                except Exception: probeset1 = None #probeset1 = affygene+':'+exon_data[1][0]
            try:
                null = int(probeset1) ### Must be an int to work in DomainGraph
                values_dg = [probeset1,affygene,'changed',dI,'NA',str(lowest_raw_p)]; values_dg = string.join(values_dg,'\t')+'\n'
                if array_type == 'junction' or array_type == 'RNASeq':
                    DG_data.write(values_dg)
                    values_srf = string.join([probeset1,'Ae',dI,str(lowest_raw_p)],'\t')+'\n'
                    if float(dI)>0:
                        SRFinder_ex_data.write(values_srf)
                    elif float(dI)<0:
                        SRFinder_in_data.write(values_srf)
            except Exception: null=[]
            
        else:
            si_pvalue = lowest_raw_p
            if si_pvalue == 1: si_pvalue = 'NA'
            if probeset1 in midas_db:
                midas_p = str(midas_db[probeset1])
                if float(midas_p)<lowest_raw_p: lowest_raw_p = float(midas_p) ###This is the lowest and SI-pvalue
            else: midas_p = ''
            
            ###Determine what type of exon-annotations are present to assign a confidence score
            if affygene in annotate_db: ###Determine the transcript clusters used to comprise a splice event (genes and exon specific)
                try:
                    gene_tc = annotate_db[affygene].TranscriptClusterIDs()
                    try: probeset_tc = [ed.SecondaryGeneID()]
                    except Exception: probeset_tc = [affygene]
                    for transcript_cluster in gene_tc: probeset_tc.append(transcript_cluster)
                    probeset_tc = makeUnique(probeset_tc)
                except Exception: probeset_tc = ''; gene_tc=''
            else:
                try:
                    try: probeset_tc = [ed.SecondaryGeneID()]
                    except Exception: probeset_tc = [affygene]
                    probeset_tc = makeUnique(probeset_tc)
                except Exception: probeset_tc = '';  gene_tc=''
            cluster_number = len(probeset_tc)
            try: alternatively_reg_tc[affygene] += probeset_tc
            except KeyError: alternatively_reg_tc[affygene] = probeset_tc

            try: last_exon_region = last_exon_region_db[affygene]
            except KeyError: last_exon_region = ''
            if cluster_number>1: exon_annot_score = 1
            
            direct_domain_alignments = ' '
            if array_type == 'exon' or array_type == 'gene' or explicit_data_type != 'null':
                try:
                    direct_domain_alignments = probeset_aligning_db[probeset1]
                    try: direct_domain_gene_alignments[affygene]+=', '+direct_domain_alignments
                    except KeyError: direct_domain_gene_alignments[affygene]=direct_domain_alignments
                except KeyError: direct_domain_alignments = ' '
            else:
                try: direct_domain_alignments = probeset_aligning_db[affygene+':'+exons1]
                except KeyError: direct_domain_alignments = ''
                
            if array_type == 'RNASeq':
                exp1 = covertLogExpressionToNonLog(exp1)
                baseline_const_exp = covertLogExpressionToNonLog(baseline_const_exp)
                fold1 = covertLogFoldToNonLog(fold1)

                adjfold1 = covertLogFoldToNonLog(adjfold1)
                mean_fold_change = covertLogFoldToNonLog(mean_fold_change)
                
            try: adj_SIp=fdr_exon_stats[probeset1].AdjP()
            except Exception: adj_SIp = 'NA'
            try: secondary_geneid = ed.SecondaryGeneID()
            except Exception: secondary_geneid = affygene
            if array_type == 'RNASeq':
                secondary_geneid = ed.NovelExon()
            
            ### Write Splicing Index results
            values= [affygene,dI,symbol,fs(description),exons1,regulation_call,probeset1,rawp1,str(lowest_raw_p),midas_p,fold1,adjfold1]
            values+=[up_exons,down_exons,fs(new_functional_attribute_str),fs(new_uniprot_exon_feature_str),fs(seq_attribute_str),fs(direct_domain_alignments),exp1]
            values+=[str(baseline_const_exp),str(si_pvalue),DV,mean_fold_change,secondary_geneid, ed.ExternalExonIDs()]
            values+=[ed.Constitutive(),ed.ExonRegionID(),ed.SplicingEvent(),last_exon_region,ed.LocationSummary()] #str(exon_annot_score)
            if probeset1 in filtered_probeset_db: values += filtered_probeset_db[probeset1]
            exon_sets = abs(float(dI)),regulation_call,event_call,exons1,exons1,midas_p
            probeset = probeset1 ### store original ID (gets converted below)
            
            ### Write DomainGraph results
            try: midas_p = str(midas_db[probeset1])
            except KeyError: midas_p = 'NA'

            ### Export significant exon/junction IDs and scores
            values_ps = [probeset1,affygene,'changed',dI,'NA',str(lowest_raw_p)]; values_ps = string.join(values_ps,'\t')+'\n'
            try: ProcessedSpliceData_data.write(values_ps)
            except Exception: None
            
            if array_type == 'gene' or array_type == 'junction' or array_type == 'RNASeq':
                if (array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null':
                    try: exon_probeset = exon_array_translation_db[affygene+':'+exon_data[1][0]][0]; probeset1 = exon_probeset; gcn+=1
                    except Exception: probeset1 = None ### don't write out a line
                else:
                    try: exon_probeset = exon_array_translation_db[probeset1][0]; probeset1 = exon_probeset; gcn+=1
                    except Exception: probeset1=None; #null=[]; #print gcn, probeset1;kill - force an error - new in version 2.0.8
            try:
                null = int(probeset1)
                values_dg = [probeset1,affygene,'changed',dI,str(si_pvalue),midas_p]; values_dg = string.join(values_dg,'\t')+'\n'
                DG_data.write(values_dg)
                values_srf = string.join([probeset1,'Ae',dI,str(lowest_raw_p)],'\t')+'\n'
                if float(dI)>0:
                    SRFinder_ex_data.write(values_srf)
                elif float(dI)<0:
                    SRFinder_in_data.write(values_srf)
            except Exception: null=[]
            values_ge = [affygene,'En',dI,str(si_pvalue),midas_p,symbol,probeset]; values_ge = string.join(values_ge,'\t')+'\n'
            goelite_data.write(values_ge)
            
        if len(ed.SplicingEvent())>2:
            try: external_exon_annot[affygene].append(ed.SplicingEvent())
            except KeyError: external_exon_annot[affygene] = [ed.SplicingEvent()]
            
        try: values = string.join(values,'\t')+'\n'
        except Exception: print values;kill
        data.write(values)
        ###Process data for gene level reports
        if float((lowest_raw_p))<=p_threshold or false_pos < 2 or lowest_raw_p == 1:
          try: comparison_count[affygene] += 1
          except KeyError: comparison_count[affygene] = 1
          try: aspire_gene_results[affygene].append(exon_sets)
          except KeyError: aspire_gene_results[affygene] = [exon_sets]
          for exon in up_exon_list:
            exon_info = exon,'upregulated'
            try: critical_gene_exons[affygene].append(exon_info)
            except KeyError: critical_gene_exons[affygene] = [exon_info]
          for exon in down_exon_list:
            exon_info = exon,'downregulated'
            try: critical_gene_exons[affygene].append(exon_info)
            except KeyError: critical_gene_exons[affygene] = [exon_info]
    data.close()
    print event_count, analysis_method, "results written to:", aspire_output,'\n'
    try: clearObjectsFromMemory(original_exon_db)
    except Exception: null=[]
    exon_array_translation_db=[]; original_exon_db=[]; probeset_to_gene=[]
    
    ### Finish writing the DomainGraph export file with non-significant probesets
    if array_type != 'AltMouse':
        for probeset in excluded_probeset_db:
            eed = excluded_probeset_db[probeset]

            try: midas_p = str(midas_db[probeset])
            except KeyError: midas_p = 'NA'
            
            ### Export significant exon/junction IDs and scores
            try: values_ps = [probeset,eed.GeneID(),'UC',eed.Score(),str(eed.TTestNormalizedRatios()),midas_p]
            except Exception: excl_probeset, geneid, score, rawp, pvalue = eed; values_ps = [probeset,geneid,'UC', str(score), str(rawp), str(pvalue)]
            values_ps = string.join(values_ps,'\t')+'\n'; ProcessedSpliceData_data.write(values_ps)
            
            ### Write DomainGraph results
            if array_type == 'gene' or array_type == 'junction' or array_type == 'RNASeq':
                try: exon_probeset = exon_array_translation_db[probeset][0]; probeset = exon_probeset; gcn+=1
                except Exception: probeset=None; # null=[] - force an error - new in version 2.0.8
                
            try: values_dg = [probeset,eed.GeneID(),'UC',eed.Score(),str(eed.TTestNormalizedRatios()),midas_p]
            except Exception:
                try:
                    excl_probeset, geneid, score, rawp, pvalue = eed
                    if ':' in probeset: probeset = excl_probeset ### Example: ENSMUSG00000029213:E2.1, make this just the numeric exclusion probeset - Not sure if DG handles non-numeric
                    values_dg = [probeset,geneid,'UC', str(score), str(rawp), str(pvalue)]
                except Exception: None
            try:
                null=int(probeset)
                values_dg = string.join(values_dg,'\t')+'\n'; DG_data.write(values_dg)
            except Exception: null=[]
        
        if array_type == 'gene' or array_type == 'junction' or array_type == 'RNASeq':
            for id in exon_array_translation_db:
                SRFinder_denom_data.write(exon_array_translation_db[id]+'\tAe\n')
        else:
            for probeset in original_exon_db:
                SRFinder_denom_data.write(probeset+'\tAe\n')
        DG_data.close()
        SRFinder_in_data.close()
        SRFinder_ex_data.close()
        SRFinder_denom_data.close()

    for affygene in direct_domain_gene_alignments:
        domains = string.split(direct_domain_gene_alignments[affygene],', ')
        domains = unique.unique(domains); domains = string.join(domains,', ')
        direct_domain_gene_alignments[affygene] = domains
        
    ### functional_attribute_db2 will be reorganized so save the database with another. Use this  
    functional_attribute_db = functional_attribute_db2
    functional_attribute_db2 = reorganize_attribute_entries(functional_attribute_db2,'no')
    external_exon_annot = eliminate_redundant_dict_values(external_exon_annot)

    protein_exon_feature_db = protein_exon_feature_db2
    protein_exon_feature_db2 = reorganize_attribute_entries(protein_exon_feature_db2,'no')
 
    ############   Export Gene Data ############        
    up_splice_val_genes = 0; down_dI_genes = 0; diff_exp_spliced_genes = 0; diff_spliced_rna_factor = 0
    ddI = 0; udI = 0

    summary_data_db['direct_domain_genes']=len(direct_domain_gene_alignments)
    summary_data_db['alt_genes']=len(aspire_gene_results)
    
    critical_gene_exons = eliminate_redundant_dict_values(critical_gene_exons)
    aspire_output_gene = root_dir+'AltResults/AlternativeOutput/' + dataset_name + analysis_method + '-exon-inclusion-GENE-results.txt'
    data = export.ExportFile(aspire_output_gene)
    if array_type == 'AltMouse': goelite_data.write("GeneID\tSystemCode\n")
    
    title = ['AffyGene','max_dI','midas-p (corresponding)','symbol','external gene ID','description','regulation_call','event_call']
    title +=['number_of_comparisons','num_effected_exons','up_exons','down_exons','functional_attribute','uniprot-ens_exon_features','direct_domain_alignments']
    title +=['pathways','mean_fold_change','exon-annotations','exon-region IDs','alternative gene ID','splice-annotation score']
    title = string.join(title,'\t')+'\n'
    data.write(title)
    for affygene in aspire_gene_results:
        if affygene in annotate_db:
            description = annotate_db[affygene].Description()
            symbol = annotate_db[affygene].Symbol()
            ensembl = annotate_db[affygene].ExternalGeneID()
            if array_type != 'AltMouse' and array_type != 'RNASeq': transcript_clusters = alternatively_reg_tc[affygene]; transcript_clusters = makeUnique(transcript_clusters); transcript_clusters = string.join(transcript_clusters,'|')
            else: transcript_clusters = affygene
            rna_processing_factor = annotate_db[affygene].RNAProcessing()
        else: description='';symbol='';ensembl=affygene;rna_processing_factor=''; transcript_clusters=''
        if ensembl in go_annotations: wpgo = go_annotations[ensembl]; goa = wpgo.Combined()
        else: goa = ''

        if array_type == 'AltMouse':
            if len(ensembl) >0: goelite_data.write(ensembl+'\tL\n')
            
        try: gene_splice_event_score[affygene].sort(); top_se_score = str(gene_splice_event_score[affygene][-1])
        except KeyError: top_se_score = 'NA'
        try: gene_regions = gene_exon_region[affygene]; gene_regions = makeUnique(gene_regions); gene_regions = string.join(gene_regions,'|')
        except KeyError: gene_regions = 'NA'
        if analysis_method == 'ASPIRE' or analysis_method == 'linearregres': number_of_comparisons = str(comparison_count[affygene])
        else: number_of_comparisons = 'NA'
        results_list = aspire_gene_results[affygene]
        results_list.sort(); results_list.reverse()
        max_dI = str(results_list[0][0])
        regulation_call = results_list[0][1]
        event_call = results_list[0][2]
        midas_p = results_list[0][-1]
        num_critical_exons = str(len(critical_gene_exons[affygene]))
        try: direct_domain_annots = direct_domain_gene_alignments[affygene]
        except KeyError: direct_domain_annots = ' '
        down_exons = ''; up_exons = ''; down_list=[]; up_list=[]
        for exon_info in critical_gene_exons[affygene]:
            exon = exon_info[0]; call = exon_info[1]
            if call == 'downregulated':
                down_exons = down_exons + exon + ','
                down_list.append(exon)
                ddI += 1
            if call == 'upregulated':
                up_exons = up_exons + exon + ','
                up_list.append(exon)
                udI += 1
        down_exons = down_exons[0:-1]
        up_exons = up_exons[0:-1]
        up_exons = add_a_space(up_exons); down_exons = add_a_space(down_exons)
        functional_annotation =''
        if affygene in functional_attribute_db2:
            number_of_functional_attributes = str(len(functional_attribute_db2[affygene]))
            attribute_list = functional_attribute_db2[affygene]
            attribute_list.sort()
            for attribute_exon_info in attribute_list:
                exon_attribute = attribute_exon_info[0]
                exon_list = attribute_exon_info[1]
                functional_annotation = functional_annotation + exon_attribute
                exons = '('
                for exon in exon_list: exons = exons + exon + ','
                exons = exons[0:-1] + '),'
                if add_exons_to_annotations == 'yes': functional_annotation = functional_annotation + exons
                else: functional_annotation = functional_annotation + ','
        functional_annotation = functional_annotation[0:-1]
        uniprot_exon_annotation = ''
        if affygene in protein_exon_feature_db2:
            number_of_functional_attributes = str(len(protein_exon_feature_db2[affygene]))
            attribute_list = protein_exon_feature_db2[affygene]; attribute_list.sort()
            for attribute_exon_info in attribute_list:
                exon_attribute = attribute_exon_info[0]
                exon_list = attribute_exon_info[1]
                uniprot_exon_annotation = uniprot_exon_annotation + exon_attribute
                exons = '('
                for exon in exon_list: exons = exons + exon + ','
                exons = exons[0:-1] + '),'
                if add_exons_to_annotations == 'yes': uniprot_exon_annotation = uniprot_exon_annotation + exons
                else: uniprot_exon_annotation = uniprot_exon_annotation + ','
        uniprot_exon_annotation = uniprot_exon_annotation[0:-1]
        if len(uniprot_exon_annotation) == 0: uniprot_exon_annotation = ' '
        if len(functional_annotation) == 0: functional_annotation = ' '
        if affygene in gene_expression_diff_db:
            mean_fold_change = gene_expression_diff_db[affygene].ConstitutiveFoldStr()
            try:
                if abs(float(mean_fold_change)) > log_fold_cutoff: diff_exp_spliced_genes += 1
            except Exception: diff_exp_spliced_genes = diff_exp_spliced_genes
        else: mean_fold_change = 'NC'
        if len(rna_processing_factor) > 2: diff_spliced_rna_factor +=1
        ###Add annotations for where in the gene structure these exons are (according to Ensembl)
        if affygene in external_exon_annot: external_gene_annot = string.join(external_exon_annot[affygene],', ')
        else: external_gene_annot = ''
        
        if array_type == 'RNASeq':
            mean_fold_change = covertLogFoldToNonLog(mean_fold_change)
                
        values =[affygene,max_dI,midas_p,symbol,ensembl,fs(description),regulation_call,event_call,number_of_comparisons]
        values+=[num_critical_exons,up_exons,down_exons,functional_annotation]
        values+=[fs(uniprot_exon_annotation),fs(direct_domain_annots),fs(goa),mean_fold_change,external_gene_annot,gene_regions,transcript_clusters,top_se_score]
        values = string.join(values,'\t')+'\n'
        data.write(values)
        ### Use results for summary statistics
        if len(up_list)>len(down_list): up_splice_val_genes +=1
        else: down_dI_genes +=1
    data.close()
    print "Gene-level results written"
    ###yes here indicates that although the truncation events will initially be filtered out, later they will be added
    ###back in without the non-truncation annotations....if there is no second database (in this case functional_attribute_db again)
    ###IF WE WANT TO FILTER OUT NON-NMD ENTRIES WHEN NMD IS PRESENT (FOR A GENE) MUST INCLUDE functional_attribute_db AS THE SECOND VARIABLE!!!!
    ###Currently, yes does nothing
    functional_annotation_db, null = grab_summary_dataset_annotations(functional_attribute_db,'','yes')

    upregulated_genes = 0; downregulated_genes = 0
    ###Calculate the number of upregulated and downregulated genes
    for affygene in gene_expression_diff_db:
        fold_val = gene_expression_diff_db[affygene].ConstitutiveFold()
        try:
            if float(fold_val) > log_fold_cutoff: upregulated_genes += 1
            elif abs(float(fold_val)) > log_fold_cutoff: downregulated_genes += 1
        except Exception: null=[]

    upregulated_rna_factor = 0; downregulated_rna_factor = 0
    ###Calculate the total number of putative RNA-processing/binding factors differentially regulated
    for affygene in gene_expression_diff_db:
        gene_fold = gene_expression_diff_db[affygene].ConstitutiveFold()
        rna_processing_factor = gene_expression_diff_db[affygene].RNAProcessing()
        if len(rna_processing_factor) > 1:
            if gene_fold>log_fold_cutoff: upregulated_rna_factor += 1
            elif abs(gene_fold)>log_fold_cutoff: downregulated_rna_factor += 1
  
    ###Generate three files for downstream functional summary
    ### functional_annotation_db2 is output to the same function as functional_annotation_db, ranked_uniprot_list_all to get all ranked uniprot annotations,
    ### and ranked_uniprot_list_coding_only to get only coding ranked uniprot annotations
    functional_annotation_db2, ranked_uniprot_list_all = grab_summary_dataset_annotations(protein_exon_feature_db,'','') #functional_attribute_db
    null, ranked_uniprot_list_coding_only = grab_summary_dataset_annotations(protein_exon_feature_db,functional_attribute_db,'') #functional_attribute_db

    functional_attribute_db=[]; protein_exon_feature_db=[]  
    ###Sumarize changes in avg protein length for each splice event
    up_protein_list=[];down_protein_list=[]; protein_length_fold_diff=[]
    for [down_protein,up_protein] in protein_length_list:
        up_protein = float(up_protein); down_protein = float(down_protein)
        down_protein_list.append(down_protein); up_protein_list.append(up_protein)
        if up_protein > 10 and down_protein > 10:
            fold_change = up_protein/down_protein; protein_length_fold_diff.append(fold_change)
    median_fold_diff = statistics.median(protein_length_fold_diff)
    try: down_avg=int(statistics.avg(down_protein_list)); up_avg=int(statistics.avg(up_protein_list))
    except Exception: down_avg=0; up_avg=0
    try:
        try:
            down_std=int(statistics.stdev(down_protein_list)); up_std=int(statistics.stdev(up_protein_list))
        except ValueError: ###If 'null' is returned fro stdev
            down_std = 0;up_std = 0
    except Exception:
        down_std = 0;up_std = 0
    if len(down_protein_list)>1 and len(up_protein_list)>1:
        try:
            #t,df,tails = statistics.ttest(down_protein_list,up_protein_list,2,3)
            #t = abs(t);df = round(df)
            #print 'ttest t:',t,'df:',df
            #p = str(statistics.t_probability(t,df))
            p = str(statistics.runComparisonStatistic(down_protein_list,up_protein_list,probability_statistic))
            #print dataset_name,p
        except Exception: p = 'NA'
        if p == 1: p = 'NA'
    else: p = 'NA'
    
    ###Calculate unique reciprocal isoforms for exon-inclusion, exclusion and mutual-exclusive events
    unique_exon_inclusion_count=0;unique_exon_exclusion_count=0;unique_mutual_exclusive_count=0;
    unique_exon_event_db = eliminate_redundant_dict_values(unique_exon_event_db)
    for affygene in unique_exon_event_db:
        isoform_entries = unique_exon_event_db[affygene]
        possibly_redundant=[]; non_redundant=[]; check_for_redundant=[]
        for entry in isoform_entries:
            if entry[0] == 1:  ### If there is only one regulated exon
                possibly_redundant.append(entry)
            else:
                non_redundant.append(entry)
                critical_exon_list = entry[1]
                for exon in critical_exon_list:
                    check_for_redundant.append(exon)
        for entry in possibly_redundant:
            exon = entry[1][0]
            if exon not in check_for_redundant:
                non_redundant.append(entry)
        for entry in non_redundant:
            if entry[2] == 'ei-ex':
                if entry[3] == 'upregulated': unique_exon_inclusion_count += 1
                else: unique_exon_exclusion_count += 1
            else: unique_mutual_exclusive_count += 1
    udI = unique_exon_inclusion_count; ddI = unique_exon_exclusion_count; mx = unique_mutual_exclusive_count

    ###Add splice event information to the functional_annotation_db
    for splice_event in splice_event_db:count = splice_event_db[splice_event]; functional_annotation_db.append((splice_event,count))
    
    if analysis_method == 'splicing-index' or analysis_method == 'FIRMA': udI='NA'; ddI='NA'
    
    summary_results_db[dataset_name[0:-1]] =  udI,ddI,mx,up_splice_val_genes,down_dI_genes,(up_splice_val_genes + down_dI_genes),upregulated_genes, downregulated_genes, diff_exp_spliced_genes, upregulated_rna_factor,downregulated_rna_factor,diff_spliced_rna_factor,down_avg,down_std,up_avg,up_std,p,median_fold_diff,functional_annotation_db
    result_list = exportComparisonSummary(dataset_name,summary_data_db,'log')
    
    ###Re-set this variable (useful for testing purposes)
    clearObjectsFromMemory(gene_expression_diff_db)
    clearObjectsFromMemory(splice_event_list); clearObjectsFromMemory(si_db); si_db=[]
    clearObjectsFromMemory(fdr_exon_stats)
    try: clearObjectsFromMemory(excluded_probeset_db); clearObjectsFromMemory(ex_db); ex_db=[]
    except Exception: ex_db=[]
    clearObjectsFromMemory(exon_db)
    #clearObjectsFromMemory(annotate_db)
    critical_probeset_annotation_db=[]; gene_expression_diff_db=[]; domain_associated_genes=[]; permute_p_values=[]
    permute_miR_inputs=[]; seq_attribute_str=[]; microRNA_count_db=[]; excluded_probeset_db=[]; fdr_exon_stats=[]
    splice_event_list=[]; critical_exon_db_len=len(critical_exon_db)#; critical_exon_db=[] deleting here will cause a global instance problem

    all_domain_gene_hits=[]; gene_splice_event_score=[]; unique_exon_event_db=[]; probeset_aligning_db=[]; ranked_uniprot_list_all=[];
    filtered_microRNA_exon_db=[]; permute_domain_inputs=[]; functional_annotation_db2=[]; functional_attribute_db2=[]; protein_length_list=[];
    ranked_uniprot_list_coding_only=[]; miR_str=[]; permute_input_list=[]; microRNA_exon_feature_db2=[]; alternatively_reg_tc=[];
    direct_domain_gene_alignments=[]; aspire_gene_results=[]; domain_gene_counts=[]; functional_annotation=[]; protein_exon_feature_db2=[];
    microRNA_attribute_db=[]; probeset_mirBS_db=[]; exon_hits=[]; critical_gene_exons=[]; gene_exon_region=[]; exon_db=[]; external_exon_annot=[];
    values=[]; down_protein_list=[]; functional_annotation_db=[]; protein_length_fold_diff=[]; comparison_count=[]; filtered_arrayids=[];
    domain_hit_gene_count_db=[]; up_protein_list=[]; probeset_domain_db=[]
    
    try: goelite_data.close()
    except Exception: null=[]
    """
    print 'local vars'
    all = [var for var in locals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all:
            try:
                if len(locals()[var])>500: print var, len(locals()[var])
            except Exception: null=[]
    """       
    return summary_results_db, summary_results_db2, aspire_output, aspire_output_gene, critical_exon_db_len
    
def deviation(dI,avg_dI,stdev_dI):
    dI = covertLogFoldToNonLogFloat(dI)
    avg_dI = covertLogFoldToNonLogFloat(avg_dI)
    stdev_dI = covertLogFoldToNonLogFloat(stdev_dI)
    return str(abs((dI-avg_dI)/stdev_dI))

def covertLogExpressionToNonLog(log_val):
    if normalization_method == 'RPKM':
        nonlog_val = (math.pow(2,float(log_val)))
    else:
        nonlog_val = (math.pow(2,float(log_val)))-1
    return str(nonlog_val)
    
def covertLogFoldToNonLog(log_val):
    try:
        if float(log_val)<0: nonlog_val = (-1/math.pow(2,(float(log_val))))
        else: nonlog_val = (math.pow(2,float(log_val)))
    except Exception: nonlog_val = log_val
    return str(nonlog_val)
    
def covertLogFoldToNonLogFloat(log_val):
    if float(log_val)<0: nonlog_val = (-1/math.pow(2,(float(log_val))))
    else: nonlog_val = (math.pow(2,float(log_val)))
    return nonlog_val
    
def checkForTransSplicing(uid,splicing_event):
    pl = string.split(uid,':')
    if len(pl)>2:
        if pl[0] not in pl[1]: ### Two different genes
            if len(splicing_event)>0: splicing_event+= '|trans-splicing'
            else: splicing_event = '|trans-splicing'
    return splicing_event

def fs(text):
    ### Formats a text entry to prevent delimiting a comma
    return '"'+text+'"'

def analyzeSplicingIndex(fold_dbase):
    """The Splicing Index (SI) represents the log ratio of the exon intensities between the two tissues after normalization
    to the gene intensities in each sample: SIi = log2((e1i/g1j)/(e2i/g2j)), for the i-th exon of the j-th gene in tissue
    type 1 or 2. The splicing indices are then subjected to a t-test to probe for differential inclusion of the exon into the gene.

    In order to determine if the change in isoform expression was statistically significant, a simple two-tailed t-test was carried
    out on the isoform ratios by grouping the 10 samples from either "tumor" or "normal" tissue.

    The method ultimately producing the highest proportion of true positives was to retain only: a) exons with a DABG p-value < 0.05,
    b) genes with a signal > 70, c) exons with a log ratio between tissues (i.e., the gene-level normalized fold change) > 0.5,
    d) Splicing Index p-values < 0.005 and e) Core exons.
    
    Gardina PJ, Clark TA, Shimada B, Staples MK, Yang Q, Veitch J, Schweitzer A, Awad T, Sugnet C, Dee S, Davies C, Williams A, Turpaz Y.
    Alternative splicing and differential gene expression in colon cancer detected by a whole genome exon array.
    BMC Genomics. 2006 Dec 27;7:325. PMID: 17192196
    """
    
    ### Used to restrict the analysis to a pre-selected set of probesets (e.g. those that have a specifc splicing pattern)
    if len(filtered_probeset_db)>0:
        temp_db={}
        for probeset in fold_dbase: temp_db[probeset]=[]
        for probeset in temp_db:
            try: filtered_probeset_db[probeset]
            except KeyError: del fold_dbase[probeset]

    ### Used to restrict the analysis to a pre-selected set of probesets (e.g. those that have a specifc splicing annotation)
    if filter_for_AS == 'yes':
        proceed = 0
        for probeset in exon_db:
            as_call = exon_db[probeset].SplicingCall()
            if as_call == 0:
                try: del fold_dbase[probeset]
                except KeyError: null=[]
            
    ### Used to the export relative individual adjusted probesets fold changes used for splicing index values   
    if export_NI_values == 'yes':
        summary_output = root_dir+'AltResults/RawSpliceData/'+species+'/'+analysis_method+'/'+dataset_name[:-1]+'.txt'
        data = export.ExportFile(summary_output)
        title = string.join(['gene\tExonID\tprobesets']+original_array_names,'\t')+'\n'; data.write(title)

    print 'Calculating splicing-index values (please be patient)...',
    if array_type == 'RNASeq': id_name = 'exon/junction IDs'
    else: id_name = 'array IDs'
    print len(fold_dbase),id_name,'beging examined'
    ###original_avg_const_exp_db contains constitutive mean expression values per group: G6953871 [7.71, 7.66]
    ###array_raw_group_values: Raw expression values in list of groups: G7072464@J935416_RC@j_at ([1.79, 2.16, 2.22], [1.68, 2.24, 1.97, 1.92, 2.12])
    ###avg_const_exp_db contains the raw constitutive expression values in a single list
    splicing_index_hash=[]; excluded_probeset_db={}; denominator_probesets=0; interaction = 0

    original_increment = int(len(exon_db)/20); increment = original_increment
    for probeset in exon_db:
        ed = exon_db[probeset]
        #include_probeset = ed.IncludeProbeset()
        if interaction == increment: increment+=original_increment; print '*',
        interaction +=1
        include_probeset = 'yes'  ###Moved this filter to import of the probeset relationship file
        ###Examines user input parameters for inclusion of probeset types in the analysis
        if include_probeset == 'yes':
            geneid = ed.GeneID()
            if probeset in fold_dbase and geneid in original_avg_const_exp_db:  ###used to search for array_raw_group_values, but when filtered by expression changes, need to filter by adj_fold_dbase
                denominator_probesets+=1
                ###Includes probesets with a calculated constitutive expression value for each gene and expression data for that probeset
                group_index = 0; si_interim_group_db={}; si_interim_group_str_db={}; ge_threshold_count=0; value_count = 0
                for group_values in array_raw_group_values[probeset]:
                    """gene_expression_value = math.pow(2,original_avg_const_exp_db[geneid][group_index])
                    ###Check to see if gene expression is > threshod for both conditions
                    if gene_expression_value>gene_expression_threshold:ge_threshold_count+=1"""
                    value_index = 0; ratio_hash=[]; ratio_str_hash=[]
                    for value in group_values:  ###Calculate normalized ratio's for each condition and save raw values for later permutation
                        #exp_val = math.pow(2,value);ge_val = math.pow(2,avg_const_exp_db[geneid][value_count]) ###To calculate a ttest we need the raw constitutive expression values, these are not in group list form but are all in a single list so keep count.
                        exp_val = value;ge_val = avg_const_exp_db[geneid][value_count]
                        exp_ratio = exp_val-ge_val; ratio_hash.append(exp_ratio); ratio_str_hash.append(str(exp_ratio))
                        value_index +=1; value_count +=1
                    si_interim_group_db[group_index] = ratio_hash
                    si_interim_group_str_db[group_index] = ratio_str_hash
                    group_index+=1
                group1_ratios = si_interim_group_db[0]; group2_ratios = si_interim_group_db[1]
                group1_mean_ratio = statistics.avg(group1_ratios); group2_mean_ratio = statistics.avg(group2_ratios)
                if export_NI_values == 'yes':
                    try: er = ed.ExonID()
                    except Exception: er = 'NA'
                    ev = string.join([geneid+'\t'+er+'\t'+probeset]+si_interim_group_str_db[0]+si_interim_group_str_db[1],'\t')+'\n'; data.write(ev)
                #if ((math.log(group1_mean_ratio,2))*(math.log(group2_mean_ratio,2)))<0: opposite_SI_log_mean = 'yes'
                if (group1_mean_ratio*group2_mean_ratio)<0: opposite_SI_log_mean = 'yes'
                else: opposite_SI_log_mean = 'no'
                try:
                    if calculate_normIntensity_p == 'yes':
                        try:
                            normIntensityP = statistics.runComparisonStatistic(group1_ratios,group2_ratios,probability_statistic)
                        except Exception: normIntensityP = 'NA' ### Occurs when analyzing two groups with no variance
                    else: normIntensityP = 'NA' ### Set to an always signficant value
                    if normIntensityP == 1: normIntensityP = 'NA'
                    splicing_index = group1_mean_ratio-group2_mean_ratio; abs_splicing_index = abs(splicing_index)

                    #if probeset == '3061323': print abs_splicing_index,normIntensityP,ed.ExonID(),group1_mean_ratio,group2_mean_ratio,math.log(group1_mean_ratio,2),math.log(group2_mean_ratio,2),((math.log(group1_mean_ratio,2))*(math.log(group2_mean_ratio,2))),opposite_SI_log_mean; kill
                    if probeset in midas_db:
                        try: midas_p = float(midas_db[probeset])
                        except ValueError:
                            midas_p = 0
                            #if abs_splicing_index>1 and normIntensityP < 0.05: print probeset,normIntensityP, abs_splicing_index;kill
                    else: midas_p = 0
                    #print ed.GeneID(),ed.ExonID(),probeset,splicing_index,normIntensityP,midas_p,group1_ratios,group2_ratios

                    if abs_splicing_index>alt_exon_logfold_cutoff and (normIntensityP < p_threshold or normIntensityP == 'NA' or normIntensityP == 1) and midas_p < p_threshold:
                        exonid = ed.ExonID(); critical_exon_list = [1,[exonid]]
                        constit_exp1 = original_avg_const_exp_db[geneid][0]
                        constit_exp2 = original_avg_const_exp_db[geneid][1]
                        ge_fold=constit_exp2-constit_exp1
                        
                        ### Re-define all of the pairwise values now that the two Splicing-Index groups to report have been determined
                        data_list1 = array_raw_group_values[probeset][0]; data_list2 = array_raw_group_values[probeset][1]
                        baseline_exp = statistics.avg(data_list1); experimental_exp = statistics.avg(data_list2); fold_change = experimental_exp - baseline_exp
                        try:
                            ttest_exp_p = statistics.runComparisonStatistic(data_list1,data_list2,probability_statistic)
                        except Exception: ttest_exp_p = 1                                
                        normInt1 = (baseline_exp-constit_exp1); normInt2 = (experimental_exp-constit_exp2); adj_fold = normInt2 - normInt1
                        ped = ProbesetExpressionData(baseline_exp, experimental_exp, fold_change, adj_fold, ttest_exp_p, '')
    
                        sid = ExonData(splicing_index,probeset,critical_exon_list,geneid,group1_ratios,group2_ratios,normIntensityP,opposite_SI_log_mean)
                        sid.setConstitutiveExpression(constit_exp1); sid.setConstitutiveFold(ge_fold); sid.setProbesetExpressionData(ped)
                        splicing_index_hash.append((splicing_index,sid))
                    else:
                        ### Also record the data for probesets that are excluded... Used by DomainGraph
                        eed = ExcludedExonData(splicing_index,geneid,normIntensityP)
                        excluded_probeset_db[probeset] = eed
                except Exception:
                    null = [] ###If this occurs, then most likely, the exon and constitutive probeset are the same

    print 'Splicing Index analysis complete'
    if export_NI_values == 'yes': data.close()
    splicing_index_hash.sort(); splicing_index_hash.reverse()
    print len(splicing_index_hash),id_name,"with evidence of Alternative expression"
    p_value_call=''; permute_p_values = {}; summary_data_db['denominator_exp_events']=denominator_probesets
    return splicing_index_hash,p_value_call,permute_p_values, excluded_probeset_db

def importResiduals(filename,probe_probeset_db):
    fn=filepath(filename); key_db = {}; x=0; prior_uid = ''; uid_gene_db={}
    for line in open(fn,'rU').xreadlines():
        if x == 0 and line[0] == '#': null=[]
        elif x == 0: x+=1
        else:
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            uid = t[0]; uid,probe = string.split(uid,'-')
            try:
                probeset = probe_probeset_db[probe]; residuals = t[1:]
                if uid == prior_uid:
                    try: uid_gene_db[probeset].append(residuals) ### Don't need to keep track of the probe ID
                    except KeyError: uid_gene_db[probeset] = [residuals]
                else: ### Hence, we have finished storing all residual data for that gene
                    if len(uid_gene_db)>0: calculateFIRMAScores(uid_gene_db); uid_gene_db={}
                    try: uid_gene_db[probeset].append(residuals) ### Don't need to keep track of the probe ID
                    except KeyError: uid_gene_db[probeset] = [residuals]
                    prior_uid = uid
            except Exception: null=[]
    ### For the last gene imported
    if len(uid_gene_db)>0: calculateFIRMAScores(uid_gene_db)

def calculateFIRMAScores(uid_gene_db):
    probeset_residuals={}; all_gene_residuals=[]; total_probes=0
    for probeset in uid_gene_db:
        residuals_list = uid_gene_db[probeset]; sample_db={}; total_probes+=len(residuals_list)
        ### For all probes in a probeset, calculate the median residual for each sample
        for residuals in residuals_list:
            index=0
            for residual in residuals:
                try: sample_db[index].append(float(residual))
                except KeyError: sample_db[index] = [float(residual)]
                all_gene_residuals.append(float(residual))
                index+=1
        for index in sample_db:
            median_residual = statistics.median(sample_db[index])
            sample_db[index] = median_residual
        probeset_residuals[probeset] = sample_db

    ### Calculate the Median absolute deviation
    """http://en.wikipedia.org/wiki/Absolute_deviation
    The median absolute deviation (also MAD) is the median absolute deviation from the median. It is a robust estimator of dispersion.
    For the example {2, 2, 3, 4, 14}: 3 is the median, so the absolute deviations from the median are {1, 1, 0, 1, 11} (or reordered as
    {0, 1, 1, 1, 11}) with a median absolute deviation of 1, in this case unaffected by the value of the outlier 14.
    Here, the global gene median will be expressed as res_gene_median.
    """

    res_gene_median = statistics.median(all_gene_residuals); subtracted_residuals=[]
    for residual in all_gene_residuals: subtracted_residuals.append(abs(res_gene_median-residual))
    gene_MAD = statistics.median(subtracted_residuals)
    #if '3263614' in probeset_residuals: print len(all_gene_residuals),all_gene_residuals
    for probeset in probeset_residuals:
        sample_db = probeset_residuals[probeset]
        for index in sample_db:
            median_residual = sample_db[index]
            try:
                firma_score = median_residual/gene_MAD
                sample_db[index] = firma_score
            except Exception: null=[]
            #if probeset == '3263614': print index, median_residual, firma_score, gene_MAD
        firma_scores[probeset] = sample_db

def importProbeToProbesets(fold_dbase):
    #print "Importing probe-to-probeset annotations (please be patient)..."
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_probeset-probes.txt'
    probeset_to_include={}
    gene2examine={}
    ### Although we want to restrict the analysis to probesets in fold_dbase, we don't want to effect the FIRMA model - filter later
    for probeset in fold_dbase:
        try: ed = exon_db[probeset]; gene2examine[ed.GeneID()]=[]
        except Exception: null=[]
    for gene in original_avg_const_exp_db: gene2examine[gene]=[]
    for probeset in exon_db:
        ed = exon_db[probeset]; geneid = ed.GeneID()
        if geneid in gene2examine:
            gene2examine[geneid].append(probeset) ### Store these so we can break things up
            probeset_to_include[probeset]=[]

    probeset_probe_db = importGenericFilteredDBList(filename,probeset_to_include)
    ### Get Residuals filename and verify it's presence
    #print "Importing comparison residuals..."
    filename_objects = string.split(dataset_name[:-1],'.p'); filename = filename_objects[0]+'.txt'
    if len(array_group_list)==2:
        filename = import_dir = root_dir+'AltExpression/FIRMA/residuals/'+array_type+'/'+species+'/'+filename
    else: filename = import_dir = root_dir+'AltExpression/FIRMA/FullDatasets/'+array_type+'/'+species+'/'+filename
    status = verifyFile(filename)
    if status != 'found':
        print_out = 'The residual file:'; print_out+= filename
        print_out+= 'was not found in the default location.\nPlease make re-run the analysis from the Beginning.'
        try: UI.WarningWindow(print_out,'Exit')
        except Exception: print print_out
        print traceback.format_exc(); badExit()

    print "Calculating FIRMA scores..."
    input_count = len(gene2examine)  ### Number of probesets or probeset pairs (junction array) alternatively regulated
    original_increment = int(input_count/20); increment = original_increment
    start_time = time.time(); x=0
            
    probe_probeset_db={}; gene_count=0; total_gene_count = 0; max_gene_count=3000; round = 1
    for gene in gene2examine:
        gene_count+=1; total_gene_count+=1; x+=1
        #if x == increment: increment+=original_increment; print '*',
        for probeset in gene2examine[gene]:
            for probe in probeset_probe_db[probeset]: probe_probeset_db[probe] = probeset
        if gene_count == max_gene_count:
            ### Import residuals and calculate primary sample/probeset FIRMA scores
            importResiduals(filename,probe_probeset_db)
            #print max_gene_count*round,"genes"
            print '*',
            gene_count=0; probe_probeset_db={}; round+=1 ### Reset these variables and re-run

    probeset_probe_db={}            
    ### Analyze residuals for the remaining probesets (< max_gene_count)
    importResiduals(filename,probe_probeset_db)
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "FIRMA scores calculted for",total_gene_count, "genes in %d seconds" % time_diff

def FIRMAanalysis(fold_dbase):
    """The FIRMA method calculates a score for each probeset and for each samples within a group of arrays, independent
    of group membership. However, in AltAnalyze, these analyses are performed dependent on group. The FIRMA score is calculated
    by obtaining residual values (residuals is a variable for each probe that can't be explained by the GC content or intensity
    of that probe) from APT, for all probes corresponding to a metaprobeset (Ensembl gene in AltAnalyze). These probe residuals
    are imported and the ratio of the median residual per probeset per sample divided by the absolute standard deviation of the
    median of all probes for all samples for that gene."""
    
    ### Used to restrict the analysis to a pre-selected set of probesets (e.g. those that have a specifc splicing pattern)
    if len(filtered_probeset_db)>0:
        temp_db={}
        for probeset in fold_dbase: temp_db[probeset]=[]
        for probeset in temp_db:
            try: filtered_probeset_db[probeset]
            except KeyError: del fold_dbase[probeset]

    ### Used to restrict the analysis to a pre-selected set of probesets (e.g. those that have a specifc splicing annotation)
    if filter_for_AS == 'yes':
        proceed = 0
        for probeset in exon_db:
            as_call = exon_db[probeset].SplicingCall()
            if as_call == 0:
                try: del fold_dbase[probeset]
                except KeyError: null=[]

    #print 'Beginning FIRMA analysis (please be patient)...'

    ### Used to the export relative individual adjusted probesets fold changes used for splicing index values   
    if export_NI_values == 'yes':
        sample_names_ordered = [] ### note: Can't use original_array_names since the order is potentially different (FIRMA stores sample data as indeces within dictionary keys)
        for group_name in array_group_list: ### THIS LIST IS USED TO MAINTAIN CONSISTENT GROUP ORDERING DURING ANALYSIS
            for sample_name in array_group_name_db[group_name]: sample_names_ordered.append(sample_name)
            
        summary_output = root_dir+'AltResults/RawSpliceData/'+species+'/'+analysis_method+'/'+dataset_name[:-1]+'.txt'
        data = export.ExportFile(summary_output)
        title = string.join(['gene-probesets']+sample_names_ordered,'\t')+'\n'; data.write(title)
    
    ### Import probes for probesets to be analyzed
    global firma_scores; firma_scores = {}
    importProbeToProbesets(fold_dbase)
    print 'FIRMA scores obtained for',len(firma_scores),'probests.'
    
    ### Group sample scores for each probeset and calculate statistics
    firma_hash=[]; excluded_probeset_db={}; denominator_probesets=0; interaction = 0
    original_increment = int(len(firma_scores)/20); increment = original_increment
    
    for probeset in firma_scores:
      if probeset in fold_dbase: ### Filter based on expression
        ed = exon_db[probeset]; geneid = ed.GeneID()
        if interaction == increment: increment+=original_increment; print '*',
        interaction +=1; denominator_probesets+=1
        sample_db = firma_scores[probeset]
        ###Use the index values from performExpressionAnalysis to assign each expression value to a new database
        firma_group_array = {}
        for group_name in array_group_db:
            for array_index in array_group_db[group_name]:
                firma_score = sample_db[array_index]
                try: firma_group_array[group_name].append(firma_score)
                except KeyError: firma_group_array[group_name] = [firma_score]

        ###array_group_list should already be unique and correctly sorted (see above)
        firma_lists=[]; index=0
        for group_name in array_group_list:
            firma_list = firma_group_array[group_name]
            if len(array_group_list)>2: firma_list = statistics.avg(firma_list), firma_list, index
            firma_lists.append(firma_list); index+=1
            
        if export_NI_values == 'yes': ### DO THIS HERE SINCE firma_lists IS SORTED BELOW!!!!
            try: er = ed.ExonID()
            except Exception: er = 'NA'
            export_list = [geneid+'\t'+er+'\t'+probeset]; export_list2=[]
            for firma_ls in firma_lists:
                if len(array_group_list)>2: firma_ls =firma_ls[1] ### See above modification of firma_list object for multiple group anlaysis
                export_list+=firma_ls
            for i in export_list: export_list2.append(str(i))
            ev = string.join(export_list2,'\t')+'\n'; data.write(ev)

        if len(array_group_list)==2:
            firma_list1 = firma_lists[0]; firma_list2 = firma_lists[-1]; firma_avg1 = statistics.avg(firma_list1); firma_avg2 = statistics.avg(firma_list2)
            index1=0; index2=1 ### Only two groups, thus only two indeces
        else: ### The below code deals with identifying the comparisons which yeild the greatest FIRMA difference
            firma_lists.sort(); index1=firma_lists[0][-1]; index2 = firma_lists[-1][-1]
            firma_list1 = firma_lists[0][1]; firma_list2 = firma_lists[-1][1]; firma_avg1 = firma_lists[0][0]; firma_avg2 = firma_lists[-1][0]
        if calculate_normIntensity_p == 'yes':
            try:
                normIntensityP = statistics.runComparisonStatistic(firma_list1,firma_list2,probability_statistic)
            except Exception: normIntensityP = 'NA' ### Occurs when analyzing two groups with no variance
        else: normIntensityP = 'NA'
        if normIntensityP == 1: normIntensityP = 'NA'
        firma_fold_change = firma_avg2 - firma_avg1
        firma_fold_change = -1*firma_fold_change   ### Make this equivalent to Splicing Index fold which is also relative to experimental not control
        if (firma_avg2*firma_avg1)<0: opposite_FIRMA_scores = 'yes'
        else: opposite_FIRMA_scores = 'no'
        
        if probeset in midas_db:
            try: midas_p = float(midas_db[probeset])
            except ValueError: midas_p = 0
        else: midas_p = 0
        #if probeset == '3263614': print firma_fold_change, normIntensityP, midas_p,'\n',firma_list1, firma_list2, [p_threshold];kill
        if abs(firma_fold_change)>alt_exon_logfold_cutoff and (normIntensityP < p_threshold or normIntensityP == 'NA') and midas_p < p_threshold:
            exonid = ed.ExonID(); critical_exon_list = [1,[exonid]]
            #gene_expression_values = original_avg_const_exp_db[geneid]
            constit_exp1 = original_avg_const_exp_db[geneid][index1]
            constit_exp2 = original_avg_const_exp_db[geneid][index2]
            ge_fold = constit_exp2-constit_exp1

            ### Re-define all of the pairwise values now that the two FIRMA groups to report have been determined
            data_list1 = array_raw_group_values[probeset][index1]; data_list2 = array_raw_group_values[probeset][index2]
            
            baseline_exp = statistics.avg(data_list1); experimental_exp = statistics.avg(data_list2); fold_change = experimental_exp - baseline_exp
            group_name1 = array_group_list[index1]; group_name2 = array_group_list[index2]
            try:
                ttest_exp_p = statistics.runComparisonStatistic(data_list1,data_list2,probability_statistic)
            except Exception: ttest_exp_p = 1                                
            normInt1 = (baseline_exp-constit_exp1); normInt2 = (experimental_exp-constit_exp2); adj_fold = normInt2 - normInt1
            ped = ProbesetExpressionData(baseline_exp, experimental_exp, fold_change, adj_fold, ttest_exp_p, group_name2+'_vs_'+group_name1)
                
            fid = ExonData(firma_fold_change,probeset,critical_exon_list,geneid,data_list1,data_list2,normIntensityP,opposite_FIRMA_scores)
            fid.setConstitutiveExpression(constit_exp1); fid.setConstitutiveFold(ge_fold); fid.setProbesetExpressionData(ped)

            firma_hash.append((firma_fold_change,fid))
            #print [[[probeset,firma_fold_change,normIntensityP,p_threshold]]]
        else:
            ### Also record the data for probesets that are excluded... Used by DomainGraph
            eed = ExcludedExonData(firma_fold_change,geneid,normIntensityP)
            excluded_probeset_db[probeset] = eed

    print 'FIRMA analysis complete'
    if export_NI_values == 'yes': data.close()
    firma_hash.sort(); firma_hash.reverse()
    print len(firma_hash),"Probesets with evidence of Alternative expression out of",len(excluded_probeset_db)+len(firma_hash)
    p_value_call=''; permute_p_values = {}; summary_data_db['denominator_exp_events']=denominator_probesets
    return firma_hash,p_value_call,permute_p_values, excluded_probeset_db

def getFilteredFilename(filename):
    if array_type == 'junction':
        filename = string.replace(filename,'.txt','-filtered.txt')
    return filename

def getExonVersionFilename(filename):
    original_filename = filename
    if array_type == 'junction' or array_type == 'RNASeq':
        if explicit_data_type != 'null':
            filename = string.replace(filename,array_type,array_type+'/'+explicit_data_type)
            ### Make sure the file exists, otherwise, use the original
            file_status = verifyFile(filename)
            #print [[filename,file_status]]
            if file_status != 'found': filename = original_filename
    return filename

def importProbesetAligningDomains(exon_db,report_type):
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_domain_aligning_probesets.txt'
    filename=getFilteredFilename(filename)
    probeset_aligning_db = importGenericDBList(filename)
    filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_indirect_domain_aligning_probesets.txt'
    filename=getFilteredFilename(filename)
    probeset_indirect_aligning_db = importGenericDBList(filename)

    if array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null'):
        new_exon_db={}; splicing_call_db={}
        for probeset_pair in exon_db:
            ### For junction analyses exon_db is really regulated_exon_junction_db, containing the inclusion,exclusion probeset tuple and an object as values
            ed = exon_db[probeset_pair]; geneid = ed.GeneID(); critical_exons = ed.CriticalExons()
            for exon in critical_exons:
                new_key = geneid+':'+exon
                try: new_exon_db[new_key].append(probeset_pair)
                except KeyError: new_exon_db[new_key] = [probeset_pair]
                try: splicing_call_db[new_key].append(ed.SplicingCall())
                except KeyError: splicing_call_db[new_key] = [ed.SplicingCall()]
        for key in new_exon_db:
            probeset_pairs = new_exon_db[key]; probeset_pair = probeset_pairs[0] ### grab one of the probeset pairs
            ed = exon_db[probeset_pair]; geneid = ed.GeneID()
            jd = SimpleJunctionData(geneid,'','','',probeset_pairs) ### use only those necessary fields for this function (probeset pairs will be called as CriticalExons)
            splicing_call_db[key].sort(); splicing_call = splicing_call_db[key][-1]; jd.setSplicingCall(splicing_call) ### Bug from 1.15 to have key be new_key?
            new_exon_db[key] = jd
        exon_db = new_exon_db

    gene_protein_ft_db={};domain_gene_count_db={};protein_functional_attribute_db={}; probeset_aligning_db2={}
    splicing_call_db=[]; new_exon_db=[] ### Clear memory

    for probeset in exon_db:
        #if probeset == '107650':
        #if probeset in probeset_aligning_db: print probeset_aligning_db[probeset];kill
        if probeset in probeset_aligning_db:
            proceed = 'no'
            if filter_for_AS == 'yes':
                as_call = exon_db[probeset].SplicingCall()
                if as_call == 1: proceed = 'yes'
            else: proceed = 'yes'
            gene = exon_db[probeset].GeneID()
            new_domain_list=[]; new_domain_list2=[]
            if report_type == 'gene' and proceed == 'yes':
                for domain in probeset_aligning_db[probeset]:
                    try: domain_gene_count_db[domain].append(gene)
                    except KeyError: domain_gene_count_db[domain] = [gene]
                    try: gene_protein_ft_db[gene].append(domain)
                    except KeyError: gene_protein_ft_db[gene]=[domain]
            elif proceed == 'yes':
                if array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null'):
                    probeset_list = exon_db[probeset].CriticalExons()
                else: probeset_list = [probeset]
                for id in probeset_list:
                    for domain in probeset_aligning_db[probeset]:
                        new_domain_list.append('(direct)'+domain)
                        new_domain_list2.append((domain,'+'))
                    new_domain_list = unique.unique(new_domain_list)
                    new_domain_list_str = string.join(new_domain_list,', ')
                    gene_protein_ft_db[gene,id] = new_domain_list2
                    probeset_aligning_db2[id] = new_domain_list_str
                    
    #print exon_db['107650']
    for probeset in exon_db:
        if probeset in probeset_indirect_aligning_db:
            proceed = 'no'
            if filter_for_AS == 'yes':
                as_call = exon_db[probeset].SplicingCall()
                if as_call == 1: proceed = 'yes'
            else: proceed = 'yes'            
            gene = exon_db[probeset].GeneID()
            new_domain_list=[]; new_domain_list2=[]
            if report_type == 'gene' and proceed == 'yes':
                for domain in probeset_indirect_aligning_db[probeset]:
                    try: domain_gene_count_db[domain].append(gene)
                    except KeyError: domain_gene_count_db[domain] = [gene]
                    try: gene_protein_ft_db[gene].append(domain)
                    except KeyError: gene_protein_ft_db[gene]=[domain]
            elif proceed == 'yes':
                if array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null'):
                    probeset_list = exon_db[probeset].CriticalExons()
                else: probeset_list = [probeset]
                for id in probeset_list:
                    for domain in probeset_indirect_aligning_db[probeset]:
                        new_domain_list.append('(indirect)'+domain)
                        new_domain_list2.append((domain,'-'))
                    new_domain_list = unique.unique(new_domain_list)
                    new_domain_list_str = string.join(new_domain_list,', ')
                    gene_protein_ft_db[gene,id] = new_domain_list2
                    probeset_aligning_db2[id] = new_domain_list_str

    domain_gene_count_db = eliminate_redundant_dict_values(domain_gene_count_db)
    gene_protein_ft_db = eliminate_redundant_dict_values(gene_protein_ft_db)

    if analysis_method == 'ASPIRE' or analysis_method == 'linearregres':
        clearObjectsFromMemory(exon_db);exon_db=[]
        try: clearObjectsFromMemory(new_exon_db)
        except Exception: null=[]
    probeset_indirect_aligning_db=[]; probeset_aligning_db=[]
    if report_type == 'perfect_match':
        gene_protein_ft_db=[];domain_gene_count_db=[];protein_functional_attribute_db=[]
        return probeset_aligning_db2
    elif report_type == 'probeset':
        probeset_aligning_db2=[]
        return gene_protein_ft_db,domain_gene_count_db,protein_functional_attribute_db  
    else:
        probeset_aligning_db2=[]; protein_functional_attribute_db=[]; probeset_aligning_db2=[]
        len_gene_protein_ft_db = len(gene_protein_ft_db); gene_protein_ft_db=[]
        return len_gene_protein_ft_db,domain_gene_count_db

def importProbesetProteinCompDomains(exon_db,report_type,comp_type):
    filename = 'AltDatabase/'+species+'/'+array_type+'/probeset-domain-annotations-'+comp_type+'.txt'
    if (array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type != 'null': filename=getFilteredFilename(filename)
    filename=getExonVersionFilename(filename)
    probeset_aligning_db = importGeneric(filename)
    filename = 'AltDatabase/'+species+'/'+array_type+'/probeset-protein-annotations-'+comp_type+'.txt'
    if (array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type != 'null': filename=getFilteredFilename(filename)
    filename=getExonVersionFilename(filename)
    
    gene_protein_ft_db={};domain_gene_count_db={}
    for probeset in exon_db:
        initial_proceed = 'no'; original_probeset = probeset
        if probeset in probeset_aligning_db: initial_proceed = 'yes'
        elif array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null'):
            ### For junction analyses exon_db is really regulated_exon_junction_db, containing the inclusion,exclusion probeset tuple and an object as values
            if '|' in probeset[0]: probeset1 = string.split(probeset[0],'|')[0]; probeset = probeset1,probeset[1]
            try: alternate_probeset_id = exon_db[probeset].InclusionLookup(); probeset = alternate_probeset_id,probeset[1]
            except Exception: null=[]
            probeset_joined = string.join(probeset,'|')
            #print [probeset_joined],[probeset]
            if probeset_joined in probeset_aligning_db: initial_proceed = 'yes'; probeset = probeset_joined
            elif probeset[0] in probeset_aligning_db: initial_proceed = 'yes'; probeset = probeset[0]
            elif probeset[1] in probeset_aligning_db: initial_proceed = 'yes'; probeset = probeset[1]
            #else: for i in probeset_aligning_db: print [i];kill
        if initial_proceed == 'yes':
            proceed = 'no'
            if filter_for_AS == 'yes':
                as_call = exon_db[original_probeset].SplicingCall()
                if as_call == 1: proceed = 'yes'
            else: proceed = 'yes'
            new_domain_list = []
            gene = exon_db[original_probeset].GeneID()
            if report_type == 'gene' and proceed == 'yes':
                for domain_data in probeset_aligning_db[probeset]:
                    try:
                        domain,call = string.split(domain_data,'|')
                    except Exception:
                        values = string.split(domain_data,'|')
                        domain = values[0]; call = values[-1] ### occurs when a | exists in the annotations from UniProt
                    try: domain_gene_count_db[domain].append(gene)
                    except KeyError: domain_gene_count_db[domain] = [gene]
                    try: gene_protein_ft_db[gene].append(domain)
                    except KeyError: gene_protein_ft_db[gene]=[domain]
            elif proceed == 'yes':
                for domain_data in probeset_aligning_db[probeset]:
                    try: domain,call = string.split(domain_data,'|')
                    except Exception:
                        values = string.split(domain_data,'|')
                        domain = values[0]; call = values[-1]
                    new_domain_list.append((domain,call))
                    #new_domain_list = string.join(new_domain_list,', ')
                gene_protein_ft_db[gene,original_probeset] = new_domain_list

    domain_gene_count_db = eliminate_redundant_dict_values(domain_gene_count_db)
    
    probeset_aligning_db=[] ### Clear memory
    probeset_aligning_protein_db = importGeneric(filename)
    
    probeset_pairs={} ### Store all possible probeset pairs as single probesets for protein-protein associations
    for probeset in exon_db:
        if len(probeset)==2:
            for p in probeset: probeset_pairs[p] = probeset
            
    if report_type == 'probeset':
        ### Below code was re-written to be more memory efficient by not storing all data in probeset-domain-annotations-*comp*.txt via generic import
        protein_functional_attribute_db={}; probeset_protein_associations={}; protein_db={}
        for probeset in exon_db:
            initial_proceed = 'no'; original_probeset = probeset
            if probeset in probeset_aligning_protein_db: initial_proceed = 'yes'
            elif array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null'):
                if '|' in probeset[0]: probeset1 = string.split(probeset[0],'|')[0]; probeset = probeset1,probeset[1]
                try: alternate_probeset_id = exon_db[probeset].InclusionLookup(); probeset = alternate_probeset_id,probeset[1]
                except Exception: null=[]                
                probeset_joined = string.join(probeset,'|')
                #print [probeset_joined],[probeset]
                if probeset_joined in probeset_aligning_protein_db: initial_proceed = 'yes'; probeset = probeset_joined
                elif probeset[0] in probeset_aligning_protein_db: initial_proceed = 'yes'; probeset = probeset[0]
                elif probeset[1] in probeset_aligning_protein_db: initial_proceed = 'yes'; probeset = probeset[1]
                #else: for i in probeset_aligning_db: print [i];kill
            if initial_proceed == 'yes':
                protein_data_list=probeset_aligning_protein_db[probeset]
                new_protein_list = []
                gene = exon_db[original_probeset].GeneID()
                for protein_data in protein_data_list:
                    protein_info,call = string.split(protein_data,'|')
                    if 'AA:' in protein_info:
                        protein_info_r = string.replace(protein_info,')','*')
                        protein_info_r = string.replace(protein_info_r,'(','*')
                        protein_info_r = string.split(protein_info_r,'*')
                        null_protein = protein_info_r[1]; hit_protein = protein_info_r[3]
                        probeset_protein_associations[original_probeset] = null_protein,hit_protein,call
                        protein_db[null_protein] = []; protein_db[hit_protein] = []
                    new_protein_list.append((protein_info,call))
                    #new_protein_list = string.join(new_domain_list,', ')
                protein_functional_attribute_db[gene,original_probeset] = new_protein_list

        filename = 'AltDatabase/'+species+'/'+array_type+'/SEQUENCE-protein-dbase_'+comp_type+'.txt'
        filename=getExonVersionFilename(filename)
        protein_seq_db = importGenericFiltered(filename,protein_db)
        for key in protein_functional_attribute_db:
            gene,probeset = key
            try:
                null_protein,hit_protein,call = probeset_protein_associations[probeset]
                null_seq = protein_seq_db[null_protein][0]; hit_seq = protein_seq_db[hit_protein][0]
                seq_attr = 'sequence: ' +'('+null_protein+')'+null_seq +' -> '+'('+hit_protein+')'+hit_seq
                protein_functional_attribute_db[key].append((seq_attr,call))
            except KeyError: null=[]
        protein_seq_db=[]; probeset_aligning_protein_db=[]
        return gene_protein_ft_db,domain_gene_count_db,protein_functional_attribute_db  
    else:
        probeset_aligning_protein_db=[]; len_gene_protein_ft_db = len(gene_protein_ft_db); gene_protein_ft_db=[]
        return len_gene_protein_ft_db,domain_gene_count_db      
        
class SimpleJunctionData:
    def __init__(self, geneid, probeset1, probeset2, probeset1_display, critical_exon_list):
        self._geneid = geneid; self._probeset1 = probeset1; self._probeset2 = probeset2
        self._probeset1_display = probeset1_display; self._critical_exon_list = critical_exon_list
    def GeneID(self): return self._geneid
    def Probeset1(self): return self._probeset1
    def Probeset2(self): return self._probeset2
    def InclusionDisplay(self): return self._probeset1_display
    def CriticalExons(self): return self._critical_exon_list
    def setSplicingCall(self,splicing_call):
        #self._splicing_call = EvidenceOfAltSplicing(slicing_annot)
        self._splicing_call = splicing_call
    def setSymbol(self,symbol): self.symbol = symbol
    def Symbol(self): return self.symbol
    def SplicingCall(self): return self._splicing_call
    def setInclusionLookup(self,incl_junction_probeset): self.incl_junction_probeset = incl_junction_probeset
    def InclusionLookup(self): return self.incl_junction_probeset
    
def formatJunctionData(probesets,affygene,critical_exon_list):
    if '|' in probesets[0]: ### Only return the first inclusion probeset (agglomerated probesets)
        incl_list = string.split(probesets[0],'|')
        incl_probeset = incl_list[0]; excl_probeset = probesets[1]
    else: incl_probeset = probesets[0]; excl_probeset = probesets[1]
    jd = SimpleJunctionData(affygene,incl_probeset,excl_probeset,probesets[0],critical_exon_list)
    key = incl_probeset,excl_probeset
    return key,jd

class JunctionExpressionData:
    def __init__(self, baseline_norm_exp, exper_norm_exp, pval, ped):
        self.baseline_norm_exp = baseline_norm_exp; self.exper_norm_exp = exper_norm_exp; self.pval = pval; self.ped = ped
    def ConNI(self):
        ls=[]
        for i in self.logConNI():
            ls.append(math.pow(2,i))
        return ls
    def ExpNI(self):
        ls=[]
        for i in self.logExpNI():
            ls.append(math.pow(2,i))
        return ls
    def ConNIAvg(self): return math.pow(2,statistics.avg(self.logConNI()))
    def ExpNIAvg(self): return math.pow(2,statistics.avg(self.logExpNI()))
    def logConNI(self): return self.baseline_norm_exp
    def logExpNI(self): return self.exper_norm_exp
    def Pval(self): return self.pval
    def ProbesetExprData(self): return self.ped
    def __repr__(self): return self.ConNI()+'|'+self.ExpNI()

    
def calculateAllASPIREScores(p1,p2):
    b1o = p1.ConNIAvg(); b2o = p2.ConNIAvg()
    e1o = p1.ExpNIAvg(); e2o = p2.ExpNIAvg(); original_score = statistics.aspire_stringent(b1o,e1o,b2o,e2o)
    
    index=0; baseline_scores=[] ### Loop through each control ratio and compare to control ratio mean
    for b1 in p1.ConNI():
        b2 = p2.ConNI()[index]
        score = statistics.aspire_stringent(b2,e2o,b1,e1o); index+=1
        baseline_scores.append(score)
        
    index=0; exp_scores=[] ### Loop through each experimental ratio and compare to control ratio mean
    for e1 in p1.ExpNI():
        e2 = p2.ExpNI()[index]
        score = statistics.aspire_stringent(b1o,e1,b2o,e2); index+=1
        exp_scores.append(score)
    
    try:
        aspireP = statistics.runComparisonStatistic(baseline_scores,exp_scores,probability_statistic)
    except Exception: aspireP = 'NA' ### Occurs when analyzing two groups with no variance
    if aspireP == 1: aspireP = 'NA'
    """
    if aspireP<0.05 and oscore>0.2 and statistics.avg(exp_scores)<0:
        index=0
        for e1 in p1.ExpNI():
            e2 = p2.ExpNI()[index]
            score = statistics.aspire_stringent(b1,e1,b2,e2)
            print p1.ExpNI(), p2.ExpNI(); print e1, e2
            print e1o,e2o; print b1, b2; print score, original_score
            print exp_scores, statistics.avg(exp_scores); kill"""

    return baseline_scores, exp_scores, aspireP

def stringListConvert(ls):
    ls2=[]
    for i in ls: ls2.append(str(i))
    return ls2
    
def analyzeJunctionSplicing(nonlog_NI_db):
    group_sizes = []; original_array_indices = permute_lists[0] ###p[0] is the original organization of the group samples prior to permutation
    for group in original_array_indices: group_sizes.append(len(group))
        
    ### Used to restrict the analysis to a pre-selected set of probesets (e.g. those that have a specifc splicing pattern)
    if len(filtered_probeset_db)>0:
        temp_db={}
        for probeset in nonlog_NI_db: temp_db[probeset]=[]
        for probeset in temp_db:
            try: filtered_probeset_db[probeset]
            except KeyError: del nonlog_NI_db[probeset]

    ### Used to the export relative individual adjusted probesets fold changes used for splicing index values   
    if export_NI_values == 'yes':
        global NIdata_export
        summary_output = root_dir+'AltResults/RawSpliceData/'+species+'/'+analysis_method+'/'+dataset_name[:-1]+'.txt'
        NIdata_export = export.ExportFile(summary_output)
        title = string.join(['inclusion-probeset','exclusion-probeset']+original_array_names,'\t')+'\n'; NIdata_export.write(title)

    ### Calculate a probeset p-value adjusted for constitutive expression levels (taken from splicing index method)
    xl=0

    probeset_normIntensity_db={}
    for probeset in array_raw_group_values:
            ed = exon_db[probeset]; geneid = ed.GeneID(); xl+=1
            #if geneid in alt_junction_db and geneid in original_avg_const_exp_db: ### Don't want this filter since it causes problems for Trans-splicing
            group_index = 0; si_interim_group_db={}; ge_threshold_count=0; value_count = 0
            ### Prepare normalized expression lists for recipricol-junction algorithms
            if geneid in avg_const_exp_db:
                for group_values in array_raw_group_values[probeset]:
                    value_index = 0; ratio_hash=[]
                    for value in group_values:  ###Calculate normalized ratio's for each condition and save raw values for later permutation
                        exp_val = value;ge_val = avg_const_exp_db[geneid][value_count]; exp_ratio = exp_val-ge_val
                        ratio_hash.append(exp_ratio); value_index +=1; value_count +=1
                    si_interim_group_db[group_index] = ratio_hash
                    group_index+=1
                group1_ratios = si_interim_group_db[0]; group2_ratios = si_interim_group_db[1]
                
                ### Calculate and store simple expression summary stats
                data_list1 = array_raw_group_values[probeset][0]; data_list2 = array_raw_group_values[probeset][1]
                baseline_exp = statistics.avg(data_list1); experimental_exp = statistics.avg(data_list2); fold_change = experimental_exp - baseline_exp
                #group_name1 = array_group_list[0]; group_name2 = array_group_list[1]
                try:
                    ttest_exp_p = statistics.runComparisonStatistic(data_list1,data_list2,probability_statistic)
                except Exception: ttest_exp_p = 'NA'
                if ttest_exp_p == 1: ttest_exp_p = 'NA'
                adj_fold = statistics.avg(group2_ratios) - statistics.avg(group1_ratios)
                ped = ProbesetExpressionData(baseline_exp, experimental_exp, fold_change, adj_fold, ttest_exp_p, '')

                try:
                    try:
                        normIntensityP = statistics.runComparisonStatistic(group1_ratios,group2_ratios,probability_statistic)
                    except Exception:
                        #print group1_ratios,group2_ratios,array_raw_group_values[probeset],avg_const_exp_db[geneid];kill
                        normIntensityP = 'NA' ###occurs for constitutive probesets
                except Exception: normIntensityP = 0
                if normIntensityP == 1: normIntensityP = 'NA'
                ji = JunctionExpressionData(group1_ratios, group2_ratios, normIntensityP, ped)
                probeset_normIntensity_db[probeset]=ji ### store and access this below
                #if probeset == 'G6899622@J916374@j_at': print normIntensityP,group1_ratios,group2_ratios;kill
                ###Concatenate the two raw expression groups into a single list for permutation analysis
                ls_concatenated = []
                for group in array_raw_group_values[probeset]:
                    for entry in group: ls_concatenated.append(entry)
                if analysis_method == 'linearregres': ###Convert out of log space
                    ls_concatenated = statistics.log_fold_conversion_fraction(ls_concatenated)
                array_raw_group_values[probeset] = ls_concatenated
            
    s = 0; t = 0; y = ''; denominator_events=0; excluded_probeset_db = {}
    splice_event_list=[]; splice_event_list_mx=[]; splice_event_list_non_mx=[]; event_mx_temp = []; permute_p_values={} #use this to exclude duplicate mx events
           
    for affygene in alt_junction_db:
        if affygene in original_avg_const_exp_db:
            constit_exp1 = original_avg_const_exp_db[affygene][0]
            constit_exp2 = original_avg_const_exp_db[affygene][1]
            ge_fold=constit_exp2-constit_exp1
            for event in alt_junction_db[affygene]:
                if array_type == 'AltMouse':
                    #event = [('ei', 'E16-E17'), ('ex', 'E16-E18')] 
                    #critical_exon_db[affygene,tuple(critical_exons)] = [1,'E'+str(e1a),'E'+str(e2b)] --- affygene,tuple(event) == key, 1 indicates both are either up or down together
                    event_call = event[0][0] + '-' + event[1][0]
                    exon_set1 = event[0][1]; exon_set2 = event[1][1]
                    probeset1 = exon_dbase[affygene,exon_set1]
                    probeset2 = exon_dbase[affygene,exon_set2]
                    critical_exon_list = critical_exon_db[affygene,tuple(event)]
                if array_type == 'junction' or array_type == 'RNASeq':
                    event_call = 'ei-ex' ### Below objects from JunctionArrayEnsemblRules - class JunctionInformation
                    probeset1 = event.InclusionProbeset(); probeset2 = event.ExclusionProbeset()
                    exon_set1 = event.InclusionJunction(); exon_set2 = event.ExclusionJunction()
                    try: novel_event = event.NovelEvent()
                    except Exception: novel_event = 'known'
                    critical_exon_list = [1,event.CriticalExonSets()]
                key,jd = formatJunctionData([probeset1,probeset2],affygene,critical_exon_list[1])
                if array_type == 'junction' or array_type == 'RNASeq':
                    try: jd.setSymbol(annotate_db[affygene].Symbol())
                    except Exception:null=[]
                #if '|' in probeset1: print probeset1, key,jd.InclusionDisplay();kill
                probeset_comp_db[key] = jd ### This is used for the permutation analysis and domain/mirBS import
                #print probeset1,probeset2, critical_exon_list,event_call,exon_set1,exon_set2;kill
                if probeset1 in nonlog_NI_db and probeset2 in nonlog_NI_db:
                    denominator_events+=1
                    try: p1 = probeset_normIntensity_db[probeset1]; p2 = probeset_normIntensity_db[probeset2]
                    except Exception:
                        print probeset1, probeset2
                        p1 = probeset_normIntensity_db[probeset1]
                        p2 = probeset_normIntensity_db[probeset2]
                #if '|' in probeset1: print  
                    pp1 = p1.Pval(); pp2 = p2.Pval()
                    
                    baseline_ratio1 = p1.ConNIAvg()
                    experimental_ratio1 = p1.ExpNIAvg()
                    baseline_ratio2 = p2.ConNIAvg()
                    experimental_ratio2 = p2.ExpNIAvg()
                    
                    ped1 = p1.ProbesetExprData()
                    ped2 = p2.ProbesetExprData()
                                    
                    Rin = ''; Rex = ''
                    r = 0 ###Variable used to determine if we should take the absolute value of dI for mutually exlcusive events
                    if event_call == 'ei-ex': #means probeset1 is an exon inclusion and probeset2 is an exon exclusion
                        Rin = baseline_ratio1/experimental_ratio1 # Rin=A/C
                        Rex = baseline_ratio2/experimental_ratio2 # Rin=B/D
                        I1=baseline_ratio1/(baseline_ratio1+baseline_ratio2)
                        I2=experimental_ratio1/(experimental_ratio1+experimental_ratio2)
                    ###When Rex is larger, the exp_ratio for exclusion is decreased in comparison to baseline.
                    ###Thus, increased inclusion (when Rin is small, inclusion is big)
                    if (Rin>1 and Rex<1): y = 'downregulated'
                    elif (Rin<1 and Rex>1): y = 'upregulated'
                    elif (Rex<Rin): y = 'downregulated'
                    else: y = 'upregulated'
                    temp_list = []

                    if event_call == 'mx-mx':
                        temp_list.append(exon_set1); temp_list.append(exon_set2);temp_list.sort()
                        if (affygene,temp_list) not in event_mx_temp: #use this logic to prevent mx entries being added more than once
                            event_mx_temp.append((affygene,temp_list))
                            ###Arbitrarily choose which exon-set will be Rin or Rex, does matter for mutually exclusive events
                            Rin = baseline_ratio1/experimental_ratio1 # Rin=A/C
                            Rex = baseline_ratio2/experimental_ratio2 # Rin=B/D
                            I1=baseline_ratio1/(baseline_ratio1+baseline_ratio2)
                            I2=experimental_ratio1/(experimental_ratio1+experimental_ratio2)
                            y = 'mutually-exclusive'; r = 1
                    if analysis_method == 'ASPIRE' and Rex != '':
                        #if affygene == 'ENSMUSG00000000126': print Rin, Rex, probeset1, probeset2
                        if (Rin>1 and Rex<1) or (Rin<1 and Rex>1):
                            s +=1
                            in1=((Rex-1.0)*Rin)/(Rex-Rin); in2=(Rex-1.0)/(Rex-Rin)
                            dI = ((in2-in1)+(I2-I1))/2.0 #modified to give propper exon inclusion
                            dI = dI*(-1) ### Reverse the fold to make equivalent to splicing-index and FIRMA scores
                            try: baseline_scores, exp_scores, aspireP = calculateAllASPIREScores(p1,p2)
                            except Exception: baseline_scores = [0]; exp_scores=[dI]; aspireP = 0
                            if export_NI_values == 'yes':
                                baseline_scores = stringListConvert(baseline_scores); exp_scores = stringListConvert(exp_scores)
                                ev = string.join([probeset1,probeset2]+baseline_scores+exp_scores,'\t')+'\n'; NIdata_export.write(ev)
                            if max_replicates >2 or equal_replicates==2:
                                permute_p_values[(probeset1,probeset2)] = [aspireP, 'NA', 'NA', 'NA']
                            if r == 1: dI = abs(dI)  ###Occurs when event is mutually exclusive
                            #if abs(dI)>alt_exon_logfold_cutoff: print [dI],pp1,pp2,aspireP;kill
                            
                            #print [affygene,dI,pp1,pp2,aspireP,event.CriticalExonSets(),probeset1,probeset2,alt_exon_logfold_cutoff,p_threshold]
                            if ((pp1<p_threshold or pp2<p_threshold) or pp1==1 or pp1=='NA') and abs(dI) > alt_exon_logfold_cutoff: ###Require that the splice event have a constitutive corrected p less than the user defined threshold
                                ejd = ExonJunctionData(dI,probeset1,probeset2,pp1,pp2,y,event_call,critical_exon_list,affygene,ped1,ped2)
                                """if probeset1 == 'ENSMUSG00000033335:E16.1-E17.1' and probeset2 == 'ENSMUSG00000033335:E16.1-E19.1':
                                    print [dI,pp1,pp2,p_threshold,alt_exon_logfold_cutoff]
                                    print baseline_scores, exp_scores, [aspireP]#;sys.exit()"""
                                ejd.setConstitutiveExpression(constit_exp1); ejd.setConstitutiveFold(ge_fold)
                                if perform_permutation_analysis == 'yes': splice_event_list.append((dI,ejd))
                                elif aspireP < permute_p_threshold or aspireP=='NA': splice_event_list.append((dI,ejd))
                                #if abs(dI)>.2: print probeset1, probeset2, critical_exon_list, [exon_set1], [exon_set2]
                                #if dI>.2 and aspireP<0.05: print baseline_scores,exp_scores,aspireP, statistics.avg(exp_scores), dI
                            elif array_type == 'junction' or array_type == 'RNASeq':
                                excluded_probeset_db[affygene+':'+event.CriticalExonSets()[0]] = probeset1, affygene, dI, 'NA', aspireP
                            if array_type == 'RNASeq':
                                try: ejd.setNovelEvent(novel_event)
                                except Exception: None
                    if analysis_method == 'linearregres' and Rex != '':
                            s+=1
                            log_fold,linregressP,rsqrd_status = getLinearRegressionScores(probeset1,probeset2,group_sizes)
                            log_fold = log_fold ### Reverse the fold to make equivalent to splicing-index and FIRMA scores
                            if max_replicates >2 or equal_replicates==2: permute_p_values[(probeset1,probeset2)] = [linregressP, 'NA', 'NA', 'NA']
                            if rsqrd_status == 'proceed':
                                if ((pp1<p_threshold or pp2<p_threshold) or pp1==1 or pp1=='NA') and abs(log_fold) > alt_exon_logfold_cutoff: ###Require that the splice event have a constitutive corrected p less than the user defined threshold
                                    ejd = ExonJunctionData(log_fold,probeset1,probeset2,pp1,pp2,y,event_call,critical_exon_list,affygene,ped1,ped2)
                                    ejd.setConstitutiveExpression(constit_exp1); ejd.setConstitutiveFold(ge_fold)
                                    if perform_permutation_analysis == 'yes': splice_event_list.append((log_fold,ejd))
                                    elif linregressP < permute_p_threshold: splice_event_list.append((log_fold,ejd))
                                    #if probeset1 ==  'G6990053@762121_762232_at' and probeset2 == 'G6990053@J926254@j_at':
                                    #print event_call, critical_exon_list,affygene, Rin, Rex, y, temp_list;kill
                                elif array_type == 'junction' or array_type == 'RNASeq':
                                    excluded_probeset_db[affygene+':'+event.CriticalExonSets()[0]] = probeset1, affygene, log_fold, 'NA', linregressP
                                if array_type == 'RNASeq':
                                    try: ejd.setNovelEvent(novel_event)
                                    except Exception: None
                    else: t +=1

    clearObjectsFromMemory(probeset_normIntensity_db)
    probeset_normIntensity_db={}; ### Potentially large memory object containing summary stats for all probesets
    statistics.adjustPermuteStats(permute_p_values)
    summary_data_db['denominator_exp_events']=denominator_events    
    print "Number of exon-events analyzed:", s
    print "Number of exon-events excluded:", t
    return splice_event_list, probeset_comp_db, permute_p_values, excluded_probeset_db

def maxReplicates():
    replicates=0; greater_than_two=0; greater_than_one=0; group_sizes=[]
    for probeset in array_raw_group_values:
        for group_values in array_raw_group_values[probeset]:
            try:
                replicates+=len(group_values); group_sizes.append(len(group_values))
                if len(group_values)>2: greater_than_two+=1
                elif len(group_values)>1: greater_than_one+=1
            except Exception: replicates+=len(array_raw_group_values[probeset]); break
        break
    group_sizes = unique.unique(group_sizes)
    if len(group_sizes) == 1: equal_replicates = group_sizes[0]
    else: equal_replicates = 0
    max_replicates = replicates/float(original_conditions)
    if max_replicates<2.01:
        if greater_than_two>0 and greater_than_one>0: max_replicates=3
    return max_replicates, equal_replicates

def furtherProcessJunctionScores(splice_event_list, probeset_comp_db, permute_p_values):
    splice_event_list.sort(); splice_event_list.reverse()
    print "filtered %s scores:" % analysis_method, len(splice_event_list)
    if perform_permutation_analysis == 'yes':
        ###*********BEGIN PERMUTATION ANALYSIS*********
        if max_replicates >2 or equal_replicates==2:
            splice_event_list, p_value_call, permute_p_values = permuteSplicingScores(splice_event_list)
        else:
            print "WARNING...Not enough replicates to perform permutation analysis."
            p_value_call=''; permute_p_values = {}
    else:
        if max_replicates >2 or equal_replicates==2:
            if probability_statistic == 'unpaired t-test':
                p_value_call=analysis_method+'-OneWayAnova'
            else:
                p_value_call=analysis_method+'-'+probability_statistic
        else:
            if probability_statistic == 'unpaired t-test':
                p_value_call='OneWayAnova'; permute_p_values = {}
            else:
                p_value_call=probability_statistic; permute_p_values = {}

    print len(splice_event_list), 'alternative events after subsequent filtering (optional)'
    ### Get ExonJunction annotaitons            
    junction_splicing_annot_db = getJunctionSplicingAnnotations(probeset_comp_db)
    
    regulated_exon_junction_db={}; new_splice_event_list=[]
    if filter_for_AS == 'yes': print "Filtering for evidence of Alternative Splicing"
    for (fold,ejd) in splice_event_list:
        proceed = 'no'
        if filter_for_AS == 'yes':
            try:
                ja = junction_splicing_annot_db[ejd.Probeset1(),ejd.Probeset2()]; splicing_call = ja.SplicingCall()
                if splicing_call == 1: proceed = 'yes'
            except KeyError: proceed = 'no'
        else: proceed = 'yes'
        if proceed == 'yes':
            key,jd = formatJunctionData([ejd.Probeset1(),ejd.Probeset2()],ejd.GeneID(),ejd.CriticalExons())
            regulated_exon_junction_db[key] = jd ### This is used for the permutation analysis and domain/mirBS import
            new_splice_event_list.append((fold,ejd))
            ### Add junction probeset lookup for reciprocal junctions composed of an exonid (not in protein database currently)
            if array_type == 'RNASeq' and '-' not in key[0]: ### Thus, it is an exon compared to a junction
                events = alt_junction_db[ejd.GeneID()]
                for ji in events:
                    if (ji.InclusionProbeset(),ji.ExclusionProbeset()) == key:
                        jd.setInclusionLookup(ji.InclusionLookup()) ### This is the source junction from which the exon ID comes from
                        probeset_comp_db[ji.InclusionLookup(),ji.ExclusionProbeset()]=jd
                        #print ji.InclusionProbeset(),ji.ExclusionProbeset(),'       ',ji.InclusionLookup()
    if filter_for_AS == 'yes': print len(new_splice_event_list), "remaining after filtering for evidence of Alternative splicing"
    filtered_exon_db = {}
    for junctions in probeset_comp_db:
            rj = probeset_comp_db[junctions] ### Add splicing annotations to the AltMouse junction DBs (needed for permutation analysis statistics and filtering)
            try: ja = junction_splicing_annot_db[junctions]; splicing_call = ja.SplicingCall(); rj.setSplicingCall(ja.SplicingCall())
            except KeyError: rj.setSplicingCall(0)
            if filter_for_AS == 'yes': filtered_exon_db[junctions] = rj
    for junctions in regulated_exon_junction_db:
            rj = regulated_exon_junction_db[junctions]
            try: ja = junction_splicing_annot_db[junctions]; rj.setSplicingCall(ja.SplicingCall())
            except KeyError: rj.setSplicingCall(0)
    if filter_for_AS == 'yes': probeset_comp_db = filtered_exon_db

    try: clearObjectsFromMemory(alt_junction_db)
    except Exception: null=[]
    
    return new_splice_event_list, p_value_call, permute_p_values, probeset_comp_db, regulated_exon_junction_db

class SplicingScoreData:
    def Method(self):
        ###e.g. ASPIRE
        return self._method
    def Score(self): return str(self._score)
    def Probeset1(self): return self._probeset1
    def Probeset2(self): return self._probeset2
    def RegulationCall(self): return self._regulation_call
    def GeneID(self): return self._geneid
    def CriticalExons(self): return self._critical_exon_list[1]
    def CriticalExonTuple(self): return self._critical_exon_list
    def TTestNormalizedRatios(self): return self._normIntensityP
    def TTestNormalizedRatios2(self): return self._normIntensityP2
    def setConstitutiveFold(self,exp_log_ratio): self._exp_log_ratio = exp_log_ratio
    def ConstitutiveFold(self): return str(self._exp_log_ratio)
    def setConstitutiveExpression(self,const_baseline): self.const_baseline = const_baseline
    def ConstitutiveExpression(self): return str(self.const_baseline)
    def setProbesetExpressionData(self,ped): self.ped1 = ped
    def ProbesetExprData1(self): return self.ped1
    def ProbesetExprData2(self): return self.ped2
    def setNovelEvent(self,novel_event): self._novel_event = novel_event
    def NovelEvent(self): return self._novel_event
    def EventCall(self):
        ###e.g. Exon inclusion (ei) Exon exclusion (ex), ei-ex, reported in that direction
        return self._event_call
    def Report(self):
        output = self.Method() +'|'+ self.GeneID() +'|'+ string.join(self.CriticalExons(),'|')
        return output
    def __repr__(self): return self.Report()

class ExonJunctionData(SplicingScoreData):
    def __init__(self,score,probeset1,probeset2,probeset1_p,probeset2_p,regulation_call,event_call,critical_exon_list,affygene,ped1,ped2):
        self._score = score; self._probeset1 = probeset1; self._probeset2 = probeset2; self._regulation_call = regulation_call
        self._event_call = event_call; self._critical_exon_list = critical_exon_list; self._geneid = affygene
        self._method = analysis_method; self._normIntensityP = probeset1_p; self._normIntensityP2 = probeset2_p
        self.ped1 = ped1; self.ped2=ped2
        
class ExonData(SplicingScoreData):
    def __init__(self,splicing_index,probeset,critical_exon_list,geneid,group1_ratios,group2_ratios,normIntensityP,opposite_SI_log_mean):
        self._score = splicing_index; self._probeset1 = probeset; self._opposite_SI_log_mean = opposite_SI_log_mean
        self._critical_exon_list = critical_exon_list; self._geneid = geneid
        self._baseline_ratio1 = group1_ratios; self._experimental_ratio1 = group2_ratios
        self._normIntensityP = normIntensityP
        self._method = analysis_method; self._event_call = 'exon-inclusion'
        if splicing_index > 0: regulation_call = 'downregulated'  ###Since baseline is the numerator ratio
        else: regulation_call = 'upregulated'
        self._regulation_call = regulation_call
    def OppositeSIRatios(self): return self._opposite_SI_log_mean

class ExcludedExonData(ExonData):
    def __init__(self,splicing_index,geneid,normIntensityP):
        self._score = splicing_index; self._geneid = geneid; self._normIntensityP = normIntensityP

def getAllPossibleLinearRegressionScores(probeset1,probeset2,positions,group_sizes):
    ### Get Raw expression values for the two probests
    p1_exp = array_raw_group_values[probeset1]
    p2_exp = array_raw_group_values[probeset2]
        
    all_possible_scores=[]; index1=0 ### Perform all possible pairwise comparisons between groups (not sure how this will work for 10+ groups)
    for (pos1a,pos2a) in positions:
        index2=0
        for (pos1b,pos2b) in positions:
            if pos1a != pos1b:
                p1_g1 = p1_exp[pos1a:pos2a]; p1_g2 = p1_exp[pos1b:pos2b]
                p2_g1 = p2_exp[pos1a:pos2a]; p2_g2 = p2_exp[pos1b:pos2b]
                #log_fold, linregressP, rsqrd = getAllLinearRegressionScores(probeset1,probeset2,p1_g1,p2_g1,p1_g2,p2_g2,len(group_sizes)) ### Used to calculate a pairwise group pvalue
                log_fold, rsqrd = performLinearRegression(p1_g1,p2_g1,p1_g2,p2_g2)
                if log_fold<0: i1,i2 = index2,index1 ### all scores should indicate upregulation
                else: i1,i2=index1,index2
                all_possible_scores.append((abs(log_fold),i1,i2))
            index2+=1
        index1+=1
    all_possible_scores.sort()
    try: log_fold,index1,index2 = all_possible_scores[-1]
    except Exception: log_fold=0; index1=0; index2=0
    return log_fold, index1, index2
    
def getLinearRegressionScores(probeset1,probeset2,group_sizes):
    ### Get Raw expression values for the two probests
    p1_exp = array_raw_group_values[probeset1]
    p2_exp = array_raw_group_values[probeset2]
    try:
        p1_g1 = p1_exp[:group_sizes[0]]; p1_g2 = p1_exp[group_sizes[0]:]
        p2_g1 = p2_exp[:group_sizes[0]]; p2_g2 = p2_exp[group_sizes[0]:]
    except Exception:
        print probeset1,probeset2
        print p1_exp
        print p2_exp
        print group_sizes
        force_kill
    log_fold, linregressP, rsqrd = getAllLinearRegressionScores(probeset1,probeset2,p1_g1,p2_g1,p1_g2,p2_g2,2)
    return log_fold, linregressP, rsqrd

def getAllLinearRegressionScores(probeset1,probeset2,p1_g1,p2_g1,p1_g2,p2_g2,groups):
    log_fold, rsqrd = performLinearRegression(p1_g1,p2_g1,p1_g2,p2_g2)
    try:
        ### Repeat for each sample versus baselines to calculate a p-value
        index=0; group1_scores=[]
        for p1_g1_sample in p1_g1:
            p2_g1_sample = p2_g1[index]
            log_f, rs = performLinearRegression(p1_g1,p2_g1,[p1_g1_sample],[p2_g1_sample])
            group1_scores.append(log_f); index+=1

        index=0; group2_scores=[]
        for p1_g2_sample in p1_g2:
            p2_g2_sample = p2_g2[index]
            log_f, rs = performLinearRegression(p1_g1,p2_g1,[p1_g2_sample],[p2_g2_sample])
            group2_scores.append(log_f); index+=1
        
        try:
            linregressP = statistics.runComparisonStatistic(group1_scores,group2_scores,probability_statistic)
        except Exception:
            linregressP = 0; group1_scores = [0]; group2_scores = [log_fold]
        if linregressP == 1: linregressP = 0
    except Exception:
        linregressP = 0; group1_scores = [0]; group2_scores = [log_fold]

    if export_NI_values == 'yes' and groups==2:
        group1_scores = stringListConvert(group1_scores)
        group2_scores = stringListConvert(group2_scores)
        ev = string.join([probeset1,probeset2]+group1_scores+group2_scores,'\t')+'\n'; NIdata_export.write(ev)
    return log_fold, linregressP, rsqrd

def performLinearRegression(p1_g1,p2_g1,p1_g2,p2_g2):
    return_rsqrd = 'no'
    if use_R == 'yes': ###Uses the RLM algorithm
        #print "Performing Linear Regression analysis using rlm."
        g1_slope = statistics.LinearRegression(p1_g1,p2_g1,return_rsqrd)
        g2_slope = statistics.LinearRegression(p1_g2,p2_g2,return_rsqrd)
    else: ###Uses a basic least squared method
        #print "Performing Linear Regression analysis using python specific methods."
        g1_slope = statistics.simpleLinRegress(p1_g1,p2_g1)
        g2_slope = statistics.simpleLinRegress(p1_g2,p2_g2)
    
    log_fold = statistics.convert_to_log_fold(g2_slope/g1_slope)
    rsqrd = 'proceed'
    #if g1_rsqrd > 0 and g2_rsqrd > 0: rsqrd = 'proceed'
    #else: rsqrd = 'hault'
    return log_fold, rsqrd

########### Permutation Analysis Functions ###########
def permuteLinearRegression(probeset1,probeset2,p):
    p1_exp = array_raw_group_values[probeset1]
    p2_exp = array_raw_group_values[probeset2]

    p1_g1, p1_g2 = permute_samples(p1_exp,p)
    p2_g1, p2_g2 = permute_samples(p2_exp,p)
    return_rsqrd = 'no'
    if use_R == 'yes': ###Uses the RLM algorithm
        g1_slope = statistics.LinearRegression(p1_g1,p2_g1,return_rsqrd)
        g2_slope = statistics.LinearRegression(p1_g2,p2_g2,return_rsqrd)
    else: ###Uses a basic least squared method
        g1_slope = statistics.simpleLinRegress(p1_g1,p2_g1)
        g2_slope = statistics.simpleLinRegress(p1_g2,p2_g2)
        
    log_fold = statistics.convert_to_log_fold(g2_slope/g1_slope)
    return log_fold

def permuteSplicingScores(splice_event_list):
    p_value_call = 'lowest_raw_p'
    permute_p_values = {}; splice_event_list2=[]
    if len(permute_lists) > 0:
        #tuple_data in splice_event_list = dI,probeset1,probeset2,y,event_call,critical_exon_list
        all_samples = []; a = 0
        for (score,x) in splice_event_list:
            ###NOTE: This reference dI differs slightly from the below calculated, since the values are calculated from raw relative ratios rather than the avg
            ###Solution: Use the first calculated dI as the reference
            score = score*(-1) ### Reverse the score to make equivalent to splicing-index and FIRMA scores
            ref_splice_val = score; probeset1 = x.Probeset1(); probeset2 = x.Probeset2(); affygene = x.GeneID()
            y = 0; p_splice_val_dist = []; count = 0; return_rsqrd = 'no'
            for p in permute_lists:            ###There are two lists in each entry
                count += 1
                permute = 'yes'
                if analysis_method == 'ASPIRE':
                    p_splice_val = permute_ASPIRE_filtered(affygene, probeset1,probeset2,p,y,ref_splice_val,x)
                elif analysis_method == 'linearregres':
                    slope_ratio = permuteLinearRegression(probeset1,probeset2,p)
                    p_splice_val = slope_ratio
                if p_splice_val != 'null': p_splice_val_dist.append(p_splice_val)
                y+=1
            p_splice_val_dist.sort()
            new_ref_splice_val = str(abs(ref_splice_val)); new_ref_splice_val = float(new_ref_splice_val[0:8]) #otherwise won't match up the scores correctly
            if analysis_method == 'linearregres':
                if ref_splice_val<0:
                    p_splice_val_dist2=[]
                    for val in p_splice_val_dist: p_splice_val_dist2.append(-1*val)
                    p_splice_val_dist=p_splice_val_dist2; p_splice_val_dist.reverse()
            p_val, pos_permute, total_permute, greater_than_true_permute = statistics.permute_p(p_splice_val_dist,new_ref_splice_val,len(permute_lists))
            #print p_val,ref_splice_val, pos_permute, total_permute, greater_than_true_permute,p_splice_val_dist[-3:];kill
            ###When two groups are of equal size, there will be 2 pos_permutes rather than 1
            if len(permute_lists[0][0]) == len(permute_lists[0][1]): greater_than_true_permute = (pos_permute/2) - 1 #size of the two groups are equal  
            else:greater_than_true_permute = (pos_permute) - 1
            if analysis_method == 'linearregres': greater_than_true_permute = (pos_permute) - 1 ###since this is a one sided test, unlike ASPIRE
                
            ###Below equation is fine if the population is large
            permute_p_values[(probeset1,probeset2)] = [p_val, pos_permute, total_permute, greater_than_true_permute]
            ###Remove non-significant linear regression results    
            if analysis_method == 'linearregres':
                if p_val <= permute_p_threshold or greater_than_true_permute < 2: splice_event_list2.append((score,x)) ###<= since many p=0.05
    print "Number of permutation p filtered splice event:",len(splice_event_list2)           
    if len(permute_p_values)>0: p_value_call = 'permuted_aspire_p-value'
    if analysis_method == 'linearregres': splice_event_list = splice_event_list2
    return splice_event_list, p_value_call, permute_p_values

def permute_ASPIRE_filtered(affygene,probeset1,probeset2,p,y,ref_splice_val,x):
    ### Get raw expression values for each permuted group for the two probesets
    b1,e1 = permute_dI(array_raw_group_values[probeset1],p)
    try: b2,e2 = permute_dI(array_raw_group_values[probeset2],p)
    except IndexError: print probeset2, array_raw_group_values[probeset2],p; kill
    ### Get the average constitutive expression values (averaged per-sample across probesets) for each permuted group
    try: bc,ec = permute_dI(avg_const_exp_db[affygene],p)
    except IndexError: print affygene, avg_const_exp_db[affygene],p; kill
    if factor_out_expression_changes == 'no':
        ec = bc
    ### Analyze the averaged ratio's of junction expression relative to permuted constitutive expression
    try: p_splice_val = abs(statistics.aspire_stringent(b1/bc,e1/ec,b2/bc,e2/ec))  ### This the permuted ASPIRE score
    except Exception: p_splice_val = 0
    #print p_splice_val, ref_splice_val, probeset1, probeset2, affygene; dog
    if y == 0:            ###The first permutation is always the real one
        ### Grab the absolute number with small number of decimal places
        try:
            new_ref_splice_val = str(p_splice_val); new_ref_splice_val = float(new_ref_splice_val[0:8])
            ref_splice_val = str(abs(ref_splice_val)); ref_splice_val = float(ref_splice_val[0:8]); y += 1
        except ValueError:
            ###Only get this error if your ref_splice_val is a null
            print y, probeset1, probeset2; print ref_splice_val, new_ref_splice_val, p
            print b1/bc,e1/ec,b2/bc,e2/ec; print (b1/bc)/(e1/ec), (b2/bc)/(e2/ec)
            print x[7],x[8],x[9],x[10]; kill
    return p_splice_val

def permute_samples(a,p):
    baseline = []; experimental = []
    for p_index in p[0]:
        baseline.append(a[p_index])  ###Append expression values for each permuted list
    for p_index in p[1]:
        experimental.append(a[p_index])
    return baseline, experimental
    
def permute_dI(all_samples,p):
    baseline, experimental = permute_samples(all_samples,p)
    #if get_non_log_avg == 'no':
    gb = statistics.avg(baseline); ge = statistics.avg(experimental)  ###Group avg baseline, group avg experimental value
    gb = statistics.log_fold_conversion_fraction(gb); ge = statistics.log_fold_conversion_fraction(ge)
    #else:
    #baseline = statistics.log_fold_conversion_fraction(baseline); experimental = statistics.log_fold_conversion_fraction(experimental)
    #gb = statistics.avg(baseline); ge = statistics.avg(experimental)  ###Group avg baseline, group avg experimental value      
    return gb,ge

def format_exon_functional_attributes(affygene,critical_probeset_list,functional_attribute_db,up_exon_list,down_exon_list,protein_length_list):
    ### Add functional attributes
    functional_attribute_list2=[]
    new_functional_attribute_str=''
    new_seq_attribute_str=''
    new_functional_attribute_list=[]
    if array_type == 'exon' or array_type == 'gene' or explicit_data_type != 'null': critical_probesets = critical_probeset_list[0]
    else: critical_probesets = tuple(critical_probeset_list)
    key = affygene,critical_probesets
    if key in functional_attribute_db:
        ###Grab exon IDs corresponding to the critical probesets
        if analysis_method == 'ASPIRE' or 'linearregres' in analysis_method:
            try: critical_exons = regulated_exon_junction_db[critical_probesets].CriticalExons() ###For junction arrays
            except Exception: print key, functional_attribute_db[key];kill
        else: critical_exons = [exon_db[critical_probesets].ExonID()] ###For exon arrays
        for exon in critical_exons:
            for entry in functional_attribute_db[key]:
                x = 0
                functional_attribute = entry[0]
                call = entry[1] # +, -, or ~
                if ('AA:' in functional_attribute) or ('ref' in functional_attribute):
                    x = 1
                if exon in up_exon_list:
                    ### design logic to determine whether up or down regulation promotes the functional change (e.g. NMD)
                    if 'ref' in functional_attribute:
                        new_functional_attribute = '(~)'+functional_attribute 
                        data_tuple = new_functional_attribute,exon
                    elif call == '+' or call == '~':
                        new_functional_attribute = '(+)'+functional_attribute 
                        data_tuple = new_functional_attribute,exon
                    elif call == '-':
                        new_functional_attribute = '(-)'+functional_attribute 
                        data_tuple = new_functional_attribute,exon
                    if 'AA:' in functional_attribute and '?' not in functional_attribute:
                        functional_attribute_temp = functional_attribute[3:]
                        if call == '+' or call == '~':
                            val1,val2 = string.split(functional_attribute_temp,'->')
                        else:
                            val2,val1 = string.split(functional_attribute_temp,'->')
                        val1,null = string.split(val1,'(')
                        val2,null = string.split(val2,'(')
                        protein_length_list.append([val1,val2])
                elif exon in down_exon_list:
                    if 'ref' in functional_attribute:
                        new_functional_attribute = '(~)'+functional_attribute 
                        data_tuple = new_functional_attribute,exon
                    elif call == '+' or call == '~':
                        new_functional_attribute = '(-)'+functional_attribute
                        data_tuple = new_functional_attribute,exon
                    elif call == '-':
                        new_functional_attribute = '(+)'+functional_attribute
                        data_tuple = new_functional_attribute,exon
                    if 'AA:' in functional_attribute and '?' not in functional_attribute:
                        functional_attribute_temp = functional_attribute[3:]
                        if call == '+' or call == '~':
                            val2,val1 = string.split(functional_attribute_temp,'->')
                        else:
                            val1,val2 = string.split(functional_attribute_temp,'->')
                        val1,null = string.split(val1,'(')
                        val2,null = string.split(val2,'(')
                        protein_length_list.append([val1,val2])
                if x == 0 or (exclude_protein_details != 'yes'):
                    try: new_functional_attribute_list.append(new_functional_attribute)
                    except UnboundLocalError:
                        print entry
                        print up_exon_list,down_exon_list
                        print exon, critical_exons
                        print critical_probesets, (key, affygene,critical_probesets)
                        for i in functional_attribute_db:
                            print i, functional_attribute_db[i]; kill
                ###remove protein sequence prediction_data
                if 'sequence' not in data_tuple[0]:
                    if x == 0 or exclude_protein_details == 'no':
                        functional_attribute_list2.append(data_tuple)
    ###Get rid of duplicates, but maintain non-alphabetical order
    new_functional_attribute_list2=[]
    for entry in new_functional_attribute_list:
        if entry not in new_functional_attribute_list2:
            new_functional_attribute_list2.append(entry)
    new_functional_attribute_list = new_functional_attribute_list2
    #new_functional_attribute_list = unique.unique(new_functional_attribute_list)
    #new_functional_attribute_list.sort()
    for entry in new_functional_attribute_list:
        if 'sequence' in entry: new_seq_attribute_str = new_seq_attribute_str + entry + ','
        else: new_functional_attribute_str = new_functional_attribute_str + entry + ','
    new_seq_attribute_str = new_seq_attribute_str[0:-1]
    new_functional_attribute_str = new_functional_attribute_str[0:-1]
    return new_functional_attribute_str, functional_attribute_list2, new_seq_attribute_str,protein_length_list

def grab_summary_dataset_annotations(functional_attribute_db,comparison_db,include_truncation_results_specifically):
    ###If a second filtering database present, filter the 1st database based on protein length changes

    fa_db={}; cp_db={} ###index the geneids for efficient recall in the next segment of code
    for (affygene,annotation) in functional_attribute_db:
        try: fa_db[affygene].append(annotation)
        except KeyError: fa_db[affygene]= [annotation]
    for (affygene,annotation) in comparison_db:
        try: cp_db[affygene].append(annotation)
        except KeyError: cp_db[affygene]= [annotation]

    functional_attribute_db_exclude = {}            
    for affygene in fa_db:
        if affygene in cp_db:
            for annotation2 in cp_db[affygene]:
                if ('trunc' in annotation2) or ('frag' in annotation2) or ('NMDs' in annotation2):
                    try: functional_attribute_db_exclude[affygene].append(annotation2)
                    except KeyError: functional_attribute_db_exclude[affygene] = [annotation2]
    functional_annotation_db = {}                    
    for (affygene,annotation) in functional_attribute_db:
        ### if we wish to filter the 1st database based on protein length changes
        if affygene not in functional_attribute_db_exclude: 
            try: functional_annotation_db[annotation] += 1
            except KeyError: functional_annotation_db[annotation] = 1
        elif include_truncation_results_specifically == 'yes':
            for annotation_val in functional_attribute_db_exclude[affygene]:
                try: functional_annotation_db[annotation_val] += 1
                except KeyError: functional_annotation_db[annotation_val] = 1
    annotation_list = []
    annotation_list_ranked = []
    for annotation in functional_annotation_db:
        if 'micro' not in annotation:
            count = functional_annotation_db[annotation]
            annotation_list.append((annotation,count))
            annotation_list_ranked.append((count,annotation))
    annotation_list_ranked.sort(); annotation_list_ranked.reverse()
    return annotation_list, annotation_list_ranked

def reorganize_attribute_entries(attribute_db1,build_attribute_direction_databases):
    attribute_db2 = {}; inclusion_attributes_hit_count={}; exclusion_attributes_hit_count={}
    genes_with_inclusion_attributes={}; genes_with_exclusion_attributes={}; 
    ###This database has unique gene, attribute information.  No attribute will now be represented more than once per gene
    for key in attribute_db1:
        ###Make gene the key and attribute (functional elements or protein information), along with the associated exons the values
        affygene = key[0];exon_attribute = key[1];exon_list = attribute_db1[key]
        exon_list = unique.unique(exon_list);exon_list.sort()
        attribute_exon_info = exon_attribute,exon_list #e.g. 5'UTR, [E1,E2,E3]
        try: attribute_db2[affygene].append(attribute_exon_info)
        except KeyError: attribute_db2[affygene] = [attribute_exon_info]
        ###Separate out attribute data by direction for over-representation analysis
        if build_attribute_direction_databases == 'yes':
            direction=exon_attribute[1:2];unique_gene_attribute=exon_attribute[3:]
            if direction == '+':
                try: inclusion_attributes_hit_count[unique_gene_attribute].append(affygene)
                except KeyError: inclusion_attributes_hit_count[unique_gene_attribute] = [affygene]
                genes_with_inclusion_attributes[affygene]=[]
            if direction == '-':
                try: exclusion_attributes_hit_count[unique_gene_attribute].append(affygene)
                except KeyError: exclusion_attributes_hit_count[unique_gene_attribute] = [affygene]
                genes_with_exclusion_attributes[affygene]=[]

    inclusion_attributes_hit_count = eliminate_redundant_dict_values(inclusion_attributes_hit_count)
    exclusion_attributes_hit_count = eliminate_redundant_dict_values(exclusion_attributes_hit_count)
    
    """for key in inclusion_attributes_hit_count:
        inclusion_attributes_hit_count[key] = len(inclusion_attributes_hit_count[key])
    for key in exclusion_attributes_hit_count:
        exclusion_attributes_hit_count[key] = len(exclusion_attributes_hit_count[key])"""

    if build_attribute_direction_databases == 'yes': return attribute_db2,inclusion_attributes_hit_count,genes_with_inclusion_attributes,exclusion_attributes_hit_count,genes_with_exclusion_attributes
    else: return attribute_db2

########### Misc. Functions ###########
def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def add_a_space(string):
    if len(string)<1:
        string = ' '
    return string

def convertToLog2(data_list):
    return map(lambda x: math.log(float(x), 2), data_list)
    
def addGlobalFudgeFactor(data_list,data_type):
    new_list = []
    if data_type == 'log':
        for item in data_list:
            new_item = statistics.log_fold_conversion_fraction(item)
            new_list.append(float(new_item) + global_addition_factor)
        new_list = convertToLog2(new_list)
    else:
        for item in data_list: new_list.append(float(item) + global_addition_factor)
    return new_list

def copyDirectoryPDFs(root_dir,AS='AS'):
    directories = ['AltResults/AlternativeOutputDirectoryDescription.pdf',
    'AltResultsDirectoryDescription.pdf',
    'ClusteringDirectoryDescription.pdf',
    'ExpressionInputDirectoryDescription.pdf',
    'ExpressionOutputDirectoryDescription.pdf',
    'GO-Elite/GO-Elite_resultsDirectoryDescription.pdf',
    'GO-EliteDirectoryDescription.pdf',
    'RootDirectoryDescription.pdf']
    import shutil
    for dir in directories:
        file = string.split(dir,'/')[-1]
        proceed=True
        if 'AltResult' in dir and AS!='AS': proceed=False
        if proceed:
            try: shutil.copyfile(filepath('Documentation/DirectoryDescription/'+file), filepath(root_dir+dir))
            except Exception: pass  

def restrictProbesets(dataset_name):
    ### Take a file with probesets and only perform the splicing-analysis on these (e.g. those already identified from a previous run with a specific pattern)
    ### Allows for propper denominator when calculating z-scores for microRNA and protein-domain ORA
    probeset_list_filename = import_dir = '/AltDatabaseNoVersion/filtering'; filtered_probeset_db={}
    if array_type == 'RNASeq': id_name = 'exon/junction IDs'
    else: id_name = 'array IDs'
    try:
        dir_list = read_directory(import_dir)
        fn_dir = filepath(import_dir[1:])
    except Exception: dir_list=[]; fn_dir=''
    if len(dir_list)>0:
        for file in dir_list:
            if file[:-4] in dataset_name:
                fn = fn_dir+'/'+file; fn = string.replace(fn,'AltDatabase','AltDatabaseNoVersion')
                filtered_probeset_db = importGeneric(fn)
                print len(filtered_probeset_db), id_name,"will be used to restrict analysis..."
    return filtered_probeset_db

def RunAltAnalyze():
  #print altanalyze_files
  #print '!!!!!starting to run alt-exon analysis'
  #returnLargeGlobalVars()
  global annotate_db; annotate_db={}; global splice_event_list; splice_event_list=[]; residuals_dirlist=[]
  global dataset_name; global constitutive_probeset_db; global exon_db; dir_list2=[]; import_dir2=''

  if array_type == 'AltMouse': import_dir = root_dir+'AltExpression/'+array_type
  elif array_type == 'exon':
      import_dir = root_dir+'AltExpression/ExonArray/'+species+'/'
  elif array_type == 'gene':
      import_dir = root_dir+'AltExpression/GeneArray/'+species+'/'
  elif array_type == 'junction':
      import_dir = root_dir+'AltExpression/JunctionArray/'+species+'/'
  else:
      import_dir = root_dir+'AltExpression/'+array_type+'/'+species+'/'
  
  #if analysis_method == 'ASPIRE' or analysis_method == 'linearregres' or analysis_method == 'splicing-index':    
  if array_type != 'AltMouse': gene_annotation_file = "AltDatabase/ensembl/"+species+"/"+species+"_Ensembl-annotations.txt"
  else: gene_annotation_file = "AltDatabase/"+species+"/"+array_type+"/"+array_type+"_gene_annotations.txt"
  annotate_db = ExonAnalyze_module.import_annotations(gene_annotation_file,array_type)

  ###Import probe-level associations    
  exon_db={}; filtered_arrayids={};filter_status='no'
  try: constitutive_probeset_db,exon_db,genes_being_analyzed = importSplicingAnnotationDatabase(probeset_annotations_file,array_type,filtered_arrayids,filter_status)
  except IOError:
        print_out = 'The annotation database: \n'+probeset_annotations_file+'\nwas not found. Ensure this file was not deleted and that the correct species has been selected.'
        try: UI.WarningWindow(print_out,'Exit'); print print_out
        except Exception: print print_out
        print traceback.format_exc()
        badExit()
  run=0
  ### Occurs when analyzing multiple conditions rather than performing a simple pair-wise comparison
  if run_from_scratch == 'Annotate External Results': import_dir = root_dir 
  elif analyze_all_conditions == 'all groups':
      import_dir = string.replace(import_dir,'AltExpression','AltExpression/FullDatasets')
      if array_type == 'AltMouse':
          import_dir = string.replace(import_dir,'FullDatasets/AltMouse','FullDatasets/AltMouse/Mm')
  elif analyze_all_conditions == 'both':
      import_dir2 = string.replace(import_dir,'AltExpression','AltExpression/FullDatasets')
      if array_type == 'AltMouse':
          import_dir2 = string.replace(import_dir2,'FullDatasets/AltMouse','FullDatasets/AltMouse/Mm')
      try: dir_list2 = read_directory(import_dir2)  #send a sub_directory to a function to identify all files in a directory
      except Exception:
        try:
            if array_type == 'exon': array_type_dir = 'ExonArray'
            elif array_type == 'gene': array_type_dir = 'GeneArray'
            elif array_type == 'junction': array_type_dir = 'GeneArray'
            else: array_type_dir = array_type
            import_dir2 = string.replace(import_dir2,'AltExpression/'+array_type_dir+'/'+species+'/','')
            import_dir2 = string.replace(import_dir2,'AltExpression/'+array_type_dir+'/','');
            dir_list2 = read_directory(import_dir2)
        except Exception:
            print_out = 'The expression files were not found. Please make\nsure you selected the correct species and array type.\n\nselected species: '+species+'\nselected array type: '+array_type+'\nselected directory:'+import_dir2
            try: UI.WarningWindow(print_out,'Exit'); print print_out
            except Exception: print print_out
            print traceback.format_exc()
            badExit()

  try: dir_list = read_directory(import_dir) #send a sub_directory to a function to identify all files in a directory
  except Exception:
        try:
            if array_type == 'exon': array_type_dir = 'ExonArray'
            elif array_type == 'gene': array_type_dir = 'GeneArray'
            elif array_type == 'junction': array_type_dir = 'JunctionArray'
            else: array_type_dir = array_type
            import_dir = string.replace(import_dir,'AltExpression/'+array_type_dir+'/'+species+'/','')
            import_dir = string.replace(import_dir,'AltExpression/'+array_type_dir+'/','');
            try: dir_list = read_directory(import_dir)
            except Exception:
                import_dir = root_dir
                dir_list = read_directory(root_dir) ### Occurs when reading in an AltAnalyze filtered file under certain conditions
        except Exception: 
            print_out = 'The expression files were not found. Please make\nsure you selected the correct species and array type.\n\nselected species: '+species+'\nselected array type: '+array_type+'\nselected directory:'+import_dir
            try: UI.WarningWindow(print_out,'Exit')
            except Exception: print print_out
            print traceback.format_exc()
            badExit()
  dir_list+=dir_list2

  ### Capture the corresponding files in the residual dir to make sure these files exist for all comparisons - won't if FIRMA was run on some files
  if analysis_method == 'FIRMA':
      try: 
          residual_dir = root_dir+'AltExpression/FIRMA/residuals/'+array_type+'/'+species+'/'
          residuals_dirlist = read_directory(residual_dir)
      except Exception: null=[]
      try:
          residual_dir = root_dir+'AltExpression/FIRMA/FullDatasets/'+array_type+'/'+species+'/'
          residuals_dirlist += read_directory(residual_dir)
      except Exception: null=[]
      dir_list_verified=[]
      for file in residuals_dirlist:
          for filename in dir_list:    
              if file[:-4] in filename: dir_list_verified.append(filename)
      dir_list = unique.unique(dir_list_verified)

  junction_biotype = 'no'
  if array_type == 'RNASeq':
      ### Check to see if user data includes junctions or just exons
      for probeset in exon_db:
          if '-' in probeset: junction_biotype = 'yes'; break

  if junction_biotype == 'no' and analysis_method != 'splicing-index' and array_type == 'RNASeq':
      dir_list=[] ### DON'T RUN ALTANALYZE WHEN JUST ANALYZING EXON DATA
      print 'No junction data to summarize... proceeding with exon analysis\n'
  elif len(dir_list)==0:
        print_out = 'No expression files available in the input directory:\n'+root_dir
        try: UI.WarningWindow(print_out,'Exit'); print print_out
        except Exception: print print_out
        badExit()

  dir_list = filterAltExpressionFiles(dir_list,altanalyze_files) ### Looks to see if the AltExpression files are for this run or from an older run
  for altanalyze_input in dir_list:    #loop through each file in the directory to output results
    ###Import probe-level associations
    if 'cel_files' in altanalyze_input:
        print_out = 'The AltExpression directory containing the necessary import file(s) is missing. Please verify the correct parameters and input directory were selected. If this error persists, contact us.'
        try: UI.WarningWindow(print_out,'Exit'); print print_out
        except Exception: print print_out
        badExit()
    if run>0: ### Only re-set these databases after the run when batch analysing multiple files
        exon_db={}; filtered_arrayids={};filter_status='no' ###Use this as a means to save memory (import multiple times - only storing different types relevant information)
        constitutive_probeset_db,exon_db,genes_being_analyzed = importSplicingAnnotationDatabase(probeset_annotations_file,array_type,filtered_arrayids,filter_status)

    if altanalyze_input in dir_list2: dataset_dir = import_dir2 +'/'+ altanalyze_input ### Then not a pairwise comparison
    else: dataset_dir = import_dir +'/'+ altanalyze_input
    dataset_name = altanalyze_input[:-4] + '-'
    print "Beginning to process",dataset_name[0:-1]
    
    ### If the user want's to restrict the analysis to preselected probesets (e.g., limma or FIRMA analysis selected)
    global filtered_probeset_db; filtered_probeset_db={}
    try: filtered_probeset_db = restrictProbesets(dataset_name)
    except Exception: null=[]

    if run_from_scratch != 'Annotate External Results':    
        ###Import expression data and stats and filter the expression data based on fold and p-value OR expression threshold
        try: conditions,adj_fold_dbase,nonlog_NI_db,dataset_name,gene_expression_diff_db,midas_db,ex_db,si_db = performExpressionAnalysis(dataset_dir,constitutive_probeset_db,exon_db,annotate_db,dataset_name)
        except IOError:
            #except Exception,exception:
            #print exception
            print traceback.format_exc()
            print_out = 'The AltAnalyze filtered expression file "'+dataset_name+'" is not propperly formatted. Review formatting requirements if this file was created by another application.'
            try: UI.WarningWindow(print_out,'Exit'); print print_out
            except Exception: print print_out
            badExit()
    else:
        conditions = 0; adj_fold_dbase={}; nonlog_NI_db={}; gene_expression_diff_db={}; ex_db={}; si_db={}
        defineEmptyExpressionVars(exon_db); adj_fold_dbase = original_fold_dbase
    ###Run Analysis
    summary_results_db, summary_results_db2, aspire_output, aspire_output_gene, number_events_analyzed = splicingAnalysisAlgorithms(nonlog_NI_db,adj_fold_dbase,dataset_name,gene_expression_diff_db,exon_db,ex_db,si_db,dataset_dir)
    aspire_output_list.append(aspire_output); aspire_output_gene_list.append(aspire_output_gene)
    
    try: clearObjectsFromMemory(exon_db); clearObjectsFromMemory(constitutive_probeset_db); constitutive_probeset_db=[]
    except Exception: null=[]
    try: clearObjectsFromMemory(last_exon_region_db);last_exon_region_db=[]
    except Exception: null=[]
    try: clearObjectsFromMemory(adj_fold_dbase);adj_fold_dbase=[]; clearObjectsFromMemory(nonlog_NI_db);nonlog_NI_db=[]
    except Exception: null=[]
    try: clearObjectsFromMemory(gene_expression_diff_db);gene_expression_diff_db=[]; clearObjectsFromMemory(midas_db);midas_db=[]
    except Exception: null=[]
    try: clearObjectsFromMemory(ex_db);ex_db=[]; clearObjectsFromMemory(si_db);si_db=[]
    except Exception: null=[]

    try: run+=1
    except Exception: run = 1
    
  if run>0: ###run = 0 if no filtered expression data present
      try: return summary_results_db, aspire_output_gene_list, number_events_analyzed
      except Exception:
        print_out = 'AltAnalyze was unable to find an expression dataset to analyze in:\n',import_dir,'\nor\n',import_dir2,'\nPlease re-run and select a valid input directory.'
        try: UI.WarningWindow(print_out,'Exit'); print print_out
        except Exception: print print_out
        badExit()
  else:
      try: clearObjectsFromMemory(exon_db); clearObjectsFromMemory(constitutive_probeset_db); constitutive_probeset_db=[]
      except Exception: null=[]
      try: clearObjectsFromMemory(last_exon_region_db);last_exon_region_db=[]
      except Exception: null=[]
      return None

def filterAltExpressionFiles(dir_list,current_files):
    dir_list2=[]
    try:
        if len(current_files) == 0: current_files = dir_list  ###if no filenames input
        for altanalzye_input in dir_list:    #loop through each file in the directory to output results
            if altanalzye_input in current_files:
                dir_list2.append(altanalzye_input)
        dir_list = dir_list2
    except Exception: dir_list = dir_list
    return dir_list

def defineEmptyExpressionVars(exon_db):
    global fold_dbase; fold_dbase={}; global original_fold_dbase; global critical_exon_db; critical_exon_db={}
    global midas_db; midas_db = {}; global max_replicates; global equal_replicates; max_replicates=0; equal_replicates=0
    for probeset in exon_db: fold_dbase[probeset]='',''
    original_fold_dbase = fold_dbase
    
def universalPrintFunction(print_items):
    log_report = open(log_file,'a')
    for item in print_items:
        if commandLineMode == 'no': ### Command-line has it's own log file write method (Logger)
            log_report.write(item+'\n')
        else: print item
    log_report.close()
    
class StatusWindow:
    def __init__(self,root,expr_var,alt_var,goelite_var,additional_var,exp_file_location_db):
            root.title('AltAnalyze version 2.1.2')
            statusVar = StringVar() ### Class method for Tkinter. Description: "Value holder for strings variables."
            self.root = root
            height = 450; width = 500
            if os.name != 'nt': height = 500; width = 600
            self.sf = PmwFreeze.ScrolledFrame(root,
                    labelpos = 'n', label_text = 'Results Status Window',
                    usehullsize = 1, hull_width = width, hull_height = height)
            self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
            self.frame = self.sf.interior()
            
            group = PmwFreeze.Group(self.sf.interior(),tag_text = 'Output')
            group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)
                
            Label(group.interior(),width=190,height=1000,justify=LEFT, bg='black', fg = 'white',anchor=NW,padx = 5,pady = 5, textvariable=statusVar).pack(fill=X,expand=Y)

            status = StringVarFile(statusVar,root) ### Likely captures the stdout
            sys.stdout = status
            for dataset in exp_file_location_db:
                fl = exp_file_location_db[dataset]; fl.setSTDOUT(sys.stdout)
            root.after(100, AltAnalyzeMain(expr_var, alt_var, goelite_var, additional_var, exp_file_location_db, root))
            try:
                root.protocol("WM_DELETE_WINDOW", self.deleteWindow)
                root.mainloop()
            except Exception: pass
    def deleteWindow(self):
        try: self.root.destroy()
        except Exception: pass
    def quit(self):
        try:
            self.root.quit()
            self.root.destroy()
        except Exception: pass
        sys.exit()

def exportComparisonSummary(dataset_name,summary_data_dbase,return_type):
    log_report = open(log_file,'a')
    result_list=[]
    for key in summary_data_dbase:
        if key != 'QC': ### The value is a list of strings
            summary_data_dbase[key] = str(summary_data_dbase[key])
    d = 'Dataset name: '+ dataset_name[:-1]; result_list.append(d+'\n')
    try:
        d = summary_data_dbase['gene_assayed']+':\tAll genes examined'; result_list.append(d)
        d = summary_data_dbase['denominator_exp_genes']+':\tExpressed genes examined for AS'; result_list.append(d)

        if explicit_data_type == 'exon-only': 
            d = summary_data_dbase['alt_events']+':\tAlternatively regulated probesets'; result_list.append(d)
            d = summary_data_dbase['denominator_exp_events']+':\tExpressed probesets examined'; result_list.append(d)  
        elif (array_type == 'AltMouse' or array_type == 'junction' or array_type == 'RNASeq') and (explicit_data_type == 'null' or return_type == 'print'):
            d = summary_data_dbase['alt_events']+':\tAlternatively regulated junction-pairs'; result_list.append(d)
            d = summary_data_dbase['denominator_exp_events']+':\tExpressed junction-pairs examined'; result_list.append(d)
        else: 
            d = summary_data_dbase['alt_events']+':\tAlternatively regulated probesets'; result_list.append(d)
            d = summary_data_dbase['denominator_exp_events']+':\tExpressed probesets examined'; result_list.append(d)
        d = summary_data_dbase['alt_genes']+':\tAlternatively regulated genes (ARGs)'; result_list.append(d)
        d = summary_data_dbase['direct_domain_genes']+':\tARGs - overlaping with domain/motifs'; result_list.append(d)
        d = summary_data_dbase['miRNA_gene_hits']+':\tARGs - overlaping with microRNA binding sites'; result_list.append(d)
    except Exception:
        pass
    
    result_list2=[]
    for d in result_list:
        if explicit_data_type == 'exon-only': d = string.replace(d,'probeset','exon')
        elif array_type == 'RNASeq': d = string.replace(d,'probeset','junction')
        result_list2.append(d)
    result_list = result_list2
    
    if return_type == 'log':
        for d in result_list: log_report.write(d+'\n')
        log_report.write('\n')
    log_report.close()
    return result_list
        
class SummaryResultsWindow:
    def __init__(self,tl,analysis_type,output_dir,dataset_name,output_type,summary_data_dbase):
        def showLink(event):
            try:
                idx = int(event.widget.tag_names(CURRENT)[1]) ### This is just the index provided below (e.g., str(0))
                #print [self.LINKS[idx]]
                if 'http://' in self.LINKS[idx]:
                    webbrowser.open(self.LINKS[idx])
                elif self.LINKS[idx][-1] == '/':
                    self.openSuppliedDirectory(self.LINKS[idx])
                else:
                   ### Instead of using this option to open a hyperlink (which is what it should do), we can open another Tk window
                    try: self.viewPNGFile(self.LINKS[idx])  ### ImageTK PNG viewer
                    except Exception:
                        try: self.ShowImageMPL(self.LINKS[idx]) ### MatPlotLib based dispaly
                        except Exception:
                            self.openPNGImage(self.LINKS[idx]) ### Native OS PNG viewer
                            #self.DisplayPlots(self.LINKS[idx]) ### GIF based display
            except Exception:
                null=[] ### anomalous error

        self.emergency_exit = False            
        self.LINKS = []
        self.tl = tl
        self.tl.title('AltAnalyze version 2.1.2')
        self.analysis_type = analysis_type

        filename = 'Config/icon.gif'
        fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(tl); can.pack(side='top'); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        use_scroll = 'yes'
    
        try: runGOElite = run_GOElite
        except Exception: runGOElite='decide_later'
  
        if 'QC' in summary_data_dbase:
            graphic_links = summary_data_dbase['QC'] ### contains hyperlinks to QC and Clustering plots
            if len(graphic_links)==0: del summary_data_dbase['QC'] ### This can be added if an analysis fails
        else:
            graphic_links = []
            
        label_text_str = 'AltAnalyze Result Summary'; height = 150; width = 500
        if analysis_type == 'AS' or 'QC' in summary_data_dbase: height = 330
        if analysis_type == 'AS' and 'QC' in summary_data_dbase: height = 330

        self.sf = PmwFreeze.ScrolledFrame(tl,
            labelpos = 'n', label_text = label_text_str,
            usehullsize = 1, hull_width = width, hull_height = height)
        self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
        self.frame = self.sf.interior()

        txt=Text(self.frame,bg='light gray',width=150, height=80)              
        txt.pack(expand=True, fill="both")
        #txt.insert(END, 'Primary Analysis Finished....\n')
        txt.insert(END, 'Results saved to:\n'+output_dir+'\n')
        
        f = Font(family="System", size=12, weight="bold")
        txt.tag_config("font", font=f)
        i=0
        
        copyDirectoryPDFs(output_dir,AS=analysis_type)
        if analysis_type == 'AS':
            txt.insert(END, '\n')
            result_list = exportComparisonSummary(dataset_name,summary_data_dbase,'print')        
            for d in result_list: txt.insert(END, d+'\n')
                
        if 'QC' in summary_data_dbase and len(graphic_links)>0:
            txt.insert(END, '\nQC and Expression Clustering Plots',"font")
            txt.insert(END, '\n\n  1) ')         
            for (name,file_dir) in graphic_links:
                txt.insert(END, name, ('link', str(i)))
                if len(graphic_links) > (i+1):
                    txt.insert(END, '\n  %s) ' % str(i+2))
                self.LINKS.append(file_dir)
                i+=1
            txt.insert(END, '\n\nView all primary plots in the folder ')
            txt.insert(END, 'DataPlots',('link', str(i))); i+=1
            self.LINKS.append(output_dir+'DataPlots/')
        else:
            url = 'http://altanalyze.readthedocs.io/en/latest/'
            self.LINKS=(url,'')
            txt.insert(END, '\nFor more information see the ')
            txt.insert(END, "AltAnalyze Online Help", ('link', str(0)))
            txt.insert(END, '\n\n')

        if runGOElite == 'run-immediately':
            txt.insert(END, '\n\nView all pathway enrichment results in the folder ')
            txt.insert(END, 'GO-Elite',('link', str(i))); i+=1
            self.LINKS.append(output_dir+'GO-Elite/')
            
        if analysis_type == 'AS':
            txt.insert(END, '\n\nView all splicing plots in the folder ')
            txt.insert(END, 'ExonPlots',('link', str(i))); i+=1
            try: self.LINKS.append(output_dir+'ExonPlots/')
            except Exception: pass
            
        txt.tag_config('link', foreground="blue", underline = 1)
        txt.tag_bind('link', '<Button-1>', showLink)
        txt.insert(END, '\n\n')

        open_results_folder = Button(tl, text = 'Results Folder', command = self.openDirectory)
        open_results_folder.pack(side = 'left', padx = 5, pady = 5);
        if analysis_type == 'AS':
            #self.dg_url = 'http://www.altanalyze.org/domaingraph.htm'
            self.dg_url = 'http://www.altanalyze.org/domaingraph.htm'
            dg_pdf_file = 'Documentation/domain_graph.pdf'; dg_pdf_file = filepath(dg_pdf_file); self.dg_pdf_file = dg_pdf_file
            text_button = Button(tl, text='Start DomainGraph in Cytoscape', command=self.SelectCytoscapeTopLevel)
            text_button.pack(side = 'right', padx = 5, pady = 5)
            self.output_dir = output_dir + "AltResults"
            self.whatNext_url = 'http://altanalyze.readthedocs.io/en/latest/' #http://www.altanalyze.org/what_next_altexon.htm'
            whatNext_pdf = 'Documentation/what_next_alt_exon.pdf'; whatNext_pdf = filepath(whatNext_pdf); self.whatNext_pdf = whatNext_pdf
            if output_type == 'parent': self.output_dir = output_dir ###Used for fake datasets
        else:
            if pathway_permutations == 'NA':
                self.output_dir = output_dir + "ExpressionOutput"
            else: self.output_dir = output_dir
            self.whatNext_url = 'http://altanalyze.readthedocs.io/en/latest/' #'http://www.altanalyze.org/what_next_expression.htm'
            whatNext_pdf = 'Documentation/what_next_GE.pdf'; whatNext_pdf = filepath(whatNext_pdf); self.whatNext_pdf = whatNext_pdf
        what_next = Button(tl, text='What Next?', command=self.whatNextlinkout)
        what_next.pack(side = 'right', padx = 5, pady = 5)
        quit_buttonTL = Button(tl,text='Close View', command=self.close)
        quit_buttonTL.pack(side = 'right', padx = 5, pady = 5)
        
        continue_to_next_win = Button(text = 'Continue', command = self.continue_win)
        continue_to_next_win.pack(side = 'right', padx = 10, pady = 10)
        quit_button = Button(root,text='Quit', command=self.quit)
        quit_button.pack(side = 'right', padx = 5, pady = 5)

        button_text = 'Help'; help_url = 'http://www.altanalyze.org/help_main.htm'; self.help_url = filepath(help_url)
        pdf_help_file = 'Documentation/AltAnalyze-Manual.pdf'; pdf_help_file = filepath(pdf_help_file); self.pdf_help_file = pdf_help_file
        help_button = Button(root, text=button_text, command=self.Helplinkout)
        help_button.pack(side = 'left', padx = 5, pady = 5)
        if self.emergency_exit == False:
            self.tl.protocol("WM_DELETE_WINDOW", self.tldeleteWindow)
            self.tl.mainloop() ###Needed to show graphic
        else:
            """ This shouldn't have to be called, but is when the topLevel window isn't closed first
            specifically if a PNG file is opened. the sys.exitfunc() should work but doesn't.
            work on this more later """
            #AltAnalyzeSetup('no')
            try: self._tls.quit(); self._tls.destroy()
            except Exception: None
            try:  self._tlx.quit(); self._tlx.destroy()
            except Exception: None
            try:  self._tlx.quit(); self._tlx.destroy()
            except Exception: None
            try:  self.tl.quit(); self.tl.destroy()
            except Exception: None
            try:  root.quit(); root.destroy()
            except Exception: None

            UI.getUpdatedParameters(array_type,species,'Process Expression file',output_dir)
            sys.exit() ### required when opening PNG files on Windows to continue (not sure why)
            #sys.exitfunc()
 
    def tldeleteWindow(self):
        try: self.tl.quit(); self.tl.destroy()
        except Exception: self.tl.destroy()
        
    def deleteTLWindow(self):
        self.emergency_exit = True
        try:  self._tls.quit(); self._tls.destroy()
        except Exception: None
        try:  self._tlx.quit(); self._tlx.destroy()
        except Exception: None
        self.tl.quit()
        self.tl.destroy()
        sys.exitfunc()
    def deleteWindow(self):
        self.emergency_exit = True
        try:  self._tls.quit(); self._tls.destroy()
        except Exception: None
        try:  self._tlx.quit(); self._tlx.destroy()
        except Exception: None
        try:
            self.tl.quit()
            self.tl.destroy()
        except Exception: None
        sys.exitfunc()
    def continue_win(self):
        self.emergency_exit = True
        try: self._tls.quit(); self._tls.destroy()
        except Exception: None
        try:  self._tlx.quit(); self._tlx.destroy()
        except Exception: None
        try: self.tl.quit(); self.tl.destroy()
        except Exception: pass
        root.quit()
        root.destroy()
        try: self.tl.grid_forget()
        except Exception: None
        try: root.grid_forget()
        except Exception: None
        sys.exitfunc()
    def openDirectory(self):
        if runningCommandLine:
            pass
        elif os.name == 'nt':
            try: os.startfile('"'+self.output_dir+'"')
            except Exception:  os.system('open "'+self.output_dir+'"')
        elif 'darwin' in sys.platform: os.system('open "'+self.output_dir+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+self.output_dir+'/"')
    def openSuppliedDirectory(self,dir):
        if runningCommandLine:
            pass
        elif os.name == 'nt':
            try: os.startfile('"'+self.output_dir+'"')
            except Exception:  os.system('open "'+dir+'"')
        elif 'darwin' in sys.platform: os.system('open "'+dir+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+dir+'/"')   
    def DGlinkout(self):
        try:
            altanalyze_path = filepath('') ### Find AltAnalye's path
            altanalyze_path = altanalyze_path[:-1]
        except Exception: null=[]
        if os.name == 'nt':
            parent_dir = 'C:/Program Files'; application_dir = 'Cytoscape_v'; application_name = 'Cytoscape.exe'
        elif 'darwin' in sys.platform:
            parent_dir = '/Applications'; application_dir = 'Cytoscape_v'; application_name = 'Cytoscape.app'
        elif 'linux' in sys.platform:
            parent_dir = '/opt'; application_dir = 'Cytoscape_v'; application_name = 'Cytoscape'
        try: openCytoscape(altanalyze_path,application_dir,application_name)
        except Exception: null=[]
        try: self._tls.destroy()
        except Exception: None
        try: ###Remove this cytoscape as the default
            file_location_defaults = UI.importDefaultFileLocations()
            del file_location_defaults['CytoscapeDir']
            UI.exportDefaultFileLocations(file_location_defaults)
        except Exception: null=[]
        self.GetHelpTopLevel(self.dg_url,self.dg_pdf_file)
    def Helplinkout(self): self.GetHelpTopLevel(self.help_url,self.pdf_help_file)
    def whatNextlinkout(self): self.GetHelpTopLevel(self.whatNext_url,self.whatNext_pdf)
    def ShowImageMPL(self,file_location):
        """ Visualization method using MatPlotLib """
        try:
            import matplotlib
            import matplotlib.pyplot as pylab
        except Exception:
            #print 'Graphical output mode disabled (requires matplotlib, numpy and scipy)'
            None
        fig = pylab.figure()
        pylab.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.00) ### Fill the plot area left to right
        ax = fig.add_subplot(111)
        ax.set_xticks([]) ### Hides ticks
        ax.set_yticks([])
        img= pylab.imread(file_location)
        imgplot = pylab.imshow(img)
        pylab.show()
    def viewPNGFile(self,png_file_dir):
        """ View PNG file within a PMW Tkinter frame """
        try: import ImageTk
        except Exception:
            from PIL import ImageTk
            from PIL import Image
        tlx = Toplevel(); self._tlx = tlx
        sf = PmwFreeze.ScrolledFrame(tlx, labelpos = 'n', label_text = '',
                usehullsize = 1, hull_width = 800, hull_height = 550)
        sf.pack(padx = 0, pady = 0, fill = 'both', expand = 1)
        frame = sf.interior()

        tlx.title(png_file_dir)
        
        img = ImageTk.PhotoImage(file=png_file_dir)
        can = Canvas(frame)
        can.pack(fill=BOTH, padx = 0, pady = 0)
        w = img.width()
        h = height=img.height()
        
        can.config(width=w, height=h)        
        can.create_image(2, 2, image=img, anchor=NW)
        tlx.mainloop()

    def openPNGImage(self,png_file_dir):
        if runningCommandLine:
            pass
        elif os.name == 'nt':
            try: os.startfile('"'+png_file_dir+'"')
            except Exception:  os.system('open "'+png_file_dir+'"')
        elif 'darwin' in sys.platform: os.system('open "'+png_file_dir+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+png_file_dir+'"')   
    def DisplayPlots(self,file_location):
        """ Native Tkinter method - Displays a gif file in a standard TopLevel window (nothing fancy) """
        tls = Toplevel(); self._tls = tls; nulls = '\t\t\t\t'; tls.title('AltAnalyze Plot Visualization')
        self.sf = PmwFreeze.ScrolledFrame(self._tls,
                labelpos = 'n', label_text = '', usehullsize = 1, hull_width = 520, hull_height = 500)
        self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
        self.frame = self.sf.interior()
        group = PmwFreeze.Group(self.sf.interior(),tag_text = file_location)
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)
        img = PhotoImage(file=filepath(file_location))
        can = Canvas(group.interior()); can.pack(side='left',padx = 10, pady = 20); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        tls.mainloop()

    def GetHelpTopLevel(self,url,pdf_file):
        try:
            config_db = UI.importConfigFile()
            ask_for_help = config_db['help'] ### hide_selection_option
        except Exception: ask_for_help = 'null'; config_db={}
        self.pdf_file = pdf_file; self.url = url
        if ask_for_help == 'null':
            message = ''; self.message = message; self.online_help = 'Online Documentation'; self.pdf_help = 'Local PDF File'
            tls = Toplevel(); self._tls = tls; nulls = '\t\t\t\t'; tls.title('Please select one of the options')
            self.sf = PmwFreeze.ScrolledFrame(self._tls,
                    labelpos = 'n', label_text = '', usehullsize = 1, hull_width = 320, hull_height = 200)
            self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
            self.frame = self.sf.interior()
            group = PmwFreeze.Group(self.sf.interior(),tag_text = 'Options')
            group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)
            filename = 'Config/icon.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
            can = Canvas(group.interior()); can.pack(side='left',padx = 10, pady = 20); can.config(width=img.width(), height=img.height())        
            can.create_image(2, 2, image=img, anchor=NW)
            l1 = Label(group.interior(), text=nulls);  l1.pack(side = 'bottom')
            text_button2 = Button(group.interior(), text=self.online_help, command=self.openOnlineHelp); text_button2.pack(side = 'top', padx = 5, pady = 5) 
            try: text_button = Button(group.interior(), text=self.pdf_help, command=self.openPDFHelp); text_button.pack(side = 'top', padx = 5, pady = 5)
            except Exception: text_button = Button(group.interior(), text=self.pdf_help, command=self.openPDFHelp); text_button.pack(side = 'top', padx = 5, pady = 5)
            text_button3 = Button(group.interior(), text='No Thanks', command=self.skipHelp); text_button3.pack(side = 'top', padx = 5, pady = 5) 
            c = Checkbutton(group.interior(), text = "Apply these settings each time", command=self.setHelpConfig); c.pack(side = 'bottom', padx = 5, pady = 0)
            tls.mainloop()
            try: tls.destroy()
            except Exception: None
        else:
            file_location_defaults = UI.importDefaultFileLocations()
            try:
                help_choice = file_location_defaults['HelpChoice'].Location()
                if help_choice == 'PDF': self.openPDFHelp()
                elif help_choice == 'http': self.openOnlineHelp()
                else: self.skip()
            except Exception: self.openPDFHelp() ### Open PDF if there's a problem
    def SelectCytoscapeTopLevel(self):
        try:
            config_db = UI.importConfigFile()
            cytoscape_type = config_db['cytoscape'] ### hide_selection_option
        except Exception: cytoscape_type = 'null'; config_db={}
        if cytoscape_type == 'null':
            message = ''; self.message = message
            tls = Toplevel(); self._tls = tls; nulls = '\t\t\t\t'; tls.title('Cytoscape Automatic Start Options')
            self.sf = PmwFreeze.ScrolledFrame(self._tls,
                    labelpos = 'n', label_text = '', usehullsize = 1, hull_width = 420, hull_height = 200)
            self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
            self.frame = self.sf.interior()
            group = PmwFreeze.Group(self.sf.interior(),tag_text = 'Options')
            group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)
            filename = 'Config/cyto-logo-smaller.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
            can = Canvas(group.interior()); can.pack(side='left',padx = 10, pady = 5); can.config(width=img.width(), height=img.height())        
            can.create_image(2, 2, image=img, anchor=NW)
            #"""
            self.local_cytoscape = 'AltAnalyze Bundled Version'; self.custom_cytoscape = 'Previously Installed Version'
            l1 = Label(group.interior(), text=nulls);  l1.pack(side = 'bottom')
            l3 = Label(group.interior(), text='Select version of Cytoscape to open:');  l3.pack(side = 'top', pady = 5)
            """
            self.local_cytoscape = '    No    '; self.custom_cytoscape = '   Yes   '
            l1 = Label(group.interior(), text=nulls);  l1.pack(side = 'bottom')
            l2 = Label(group.interior(), text='Note: Cytoscape can take up-to a minute to initalize', fg="red");  l2.pack(side = 'top', padx = 5, pady = 0)
            """
            text_button2 = Button(group.interior(), text=self.local_cytoscape, command=self.DGlinkout); text_button2.pack(padx = 5, pady = 5) 
            try: text_button = Button(group.interior(), text=self.custom_cytoscape, command=self.getPath); text_button.pack(padx = 5, pady = 5)
            except Exception: text_button = Button(group.interior(), text=self.custom_cytoscape, command=self.getPath); text_button.pack(padx = 5, pady = 5)
            l2 = Label(group.interior(), text='Note: Cytoscape can take up-to a minute to initalize', fg="blue");  l2.pack(side = 'bottom', padx = 5, pady = 0)
            c = Checkbutton(group.interior(), text = "Apply these settings each time and don't show again", command=self.setCytoscapeConfig); c.pack(side = 'bottom', padx = 5, pady = 0)
            #c2 = Checkbutton(group.interior(), text = "Open PDF of DomainGraph help rather than online help", command=self.setCytoscapeConfig); c2.pack(side = 'bottom', padx = 5, pady = 0)
            tls.mainloop()
            try: tls.destroy()
            except Exception: None
        else:
            file_location_defaults = UI.importDefaultFileLocations()
            try: cytoscape_app_dir = file_location_defaults['CytoscapeDir'].Location(); openFile(cytoscape_app_dir)
            except Exception:
                try: altanalyze_path = filepath(''); altanalyze_path = altanalyze_path[:-1]
                except Exception: altanalyze_path=''
                application_dir = 'Cytoscape_v'
                if os.name == 'nt': application_name = 'Cytoscape.exe'
                elif 'darwin' in sys.platform: application_name = 'Cytoscape.app'
                elif 'linux' in sys.platform: application_name = 'Cytoscape'
                try: openCytoscape(altanalyze_path,application_dir,application_name)
                except Exception: null=[]
        
    def setCytoscapeConfig(self):
        config_db={}; config_db['cytoscape'] = 'hide_selection_option'
        UI.exportConfigFile(config_db)

    def setHelpConfig(self):
        config_db={}; config_db['help'] = 'hide_selection_option'
        UI.exportConfigFile(config_db)
        
    def getPath(self):
        file_location_defaults = UI.importDefaultFileLocations()
        if os.name == 'nt': parent_dir = 'C:/Program Files'; application_dir = 'Cytoscape_v'; application_name = 'Cytoscape.exe'
        elif 'darwin' in sys.platform: parent_dir = '/Applications'; application_dir = 'Cytoscape_v'; application_name = 'Cytoscape.app'
        elif 'linux' in sys.platform: parent_dir = '/opt'; application_dir = 'Cytoscape_v'; application_name = 'Cytoscape'
        try:
            self.default_dir = file_location_defaults['CytoscapeDir'].Location()
            self.default_dir = string.replace(self.default_dir,'//','/')
            self.default_dir = string.replace(self.default_dir,'\\','/')
            self.default_dir = string.join(string.split(self.default_dir,'/')[:-1],'/')
        except Exception: 
            dir = FindDir(parent_dir,application_dir); dir = filepath(parent_dir+'/'+dir) 
            self.default_dir = filepath(parent_dir)
        try: dirPath = tkFileDialog.askdirectory(parent=self._tls,initialdir=self.default_dir)
        except Exception: 
            self.default_dir = ''
            try: dirPath = tkFileDialog.askdirectory(parent=self._tls,initialdir=self.default_dir)
            except Exception: 
                try: dirPath = tkFileDialog.askdirectory(parent=self._tls)
                except Exception: dirPath=''
        try:
            #print [dirPath],application_name
            app_dir = dirPath+'/'+application_name

            if 'linux' in sys.platform:
                try: createCytoscapeDesktop(cytoscape_dir)
                except Exception: null=[]
                dir_list = unique.read_directory('/usr/bin/') ### Check to see that JAVA is installed
                if 'java' not in dir_list: print 'Java not referenced in "usr/bin/. If not installed,\nplease install and re-try opening Cytoscape'
                try:
                    jar_path = dirPath+'/cytoscape.jar'
                    main_path = dirPath+'/cytoscape.CyMain'
                    plugins_path = dirPath+'/plugins'
                    os.system('java -Dswing.aatext=true -Xss5M -Xmx512M -jar '+jar_path+' '+main_path+'  -p '+plugins_path+' &')
                    print 'Cytoscape jar opened:',jar_path
                except Exception:
                    print 'OS command to open Java failed.'      
                    try: openFile(app_dir2); print 'Cytoscape opened:',app_dir2
                    except Exception: openFile(app_dir)
            else: openFile(app_dir)
    
            try: file_location_defaults['CytoscapeDir'].SetLocation(app_dir)
            except Exception:
                fl = UI.FileLocationData('', app_dir, 'all')
                file_location_defaults['CytoscapeDir'] = fl
            UI.exportDefaultFileLocations(file_location_defaults)
        except Exception: null=[]
        try: self._tls.destroy()
        except Exception: None
        self.GetHelpTopLevel(self.dg_url,self.dg_pdf_file)
    def openOnlineHelp(self):
        file_location_defaults = UI.importDefaultFileLocations()
        try:file_location_defaults['HelpChoice'].SetLocation('http')
        except Exception:
            fl = UI.FileLocationData('', 'http', 'all')
            file_location_defaults['HelpChoice'] = fl
        UI.exportDefaultFileLocations(file_location_defaults)
        webbrowser.open(self.url)
        #except Exception: null=[]
        try: self._tls.destroy()
        except Exception: None
    def skipHelp(self):
        file_location_defaults = UI.importDefaultFileLocations()
        try: file_location_defaults['HelpChoice'].SetLocation('skip')
        except Exception:
            fl = UI.FileLocationData('', 'skip', 'all')
            file_location_defaults['HelpChoice'] = fl
        UI.exportDefaultFileLocations(file_location_defaults)
        try: self._tls.destroy()
        except Exception: None
    def openPDFHelp(self):
        file_location_defaults = UI.importDefaultFileLocations()
        try:file_location_defaults['HelpChoice'].SetLocation('PDF')
        except Exception:
            fl = UI.FileLocationData('', 'PDF', 'all')
            file_location_defaults['HelpChoice'] = fl
        UI.exportDefaultFileLocations(file_location_defaults)
        if runningCommandLine:
            pass
        elif os.name == 'nt':
            try: os.startfile('"'+self.pdf_file+'"')
            except Exception:  os.system('open "'+self.pdf_file+'"')
        elif 'darwin' in sys.platform: os.system('open "'+self.pdf_file+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+self.pdf_file+'"')   
        try: self._tls.destroy()
        except Exception: None
    def quit(self):
        root.quit()
        root.destroy()
        sys.exit()
    def close(self):
        #self.tl.quit()   #### This was causing multiple errors in 2.0.7 - evaluate more!
        self.tl.destroy()
                
class StringVarFile:
    def __init__(self,stringVar,window):
        self.__newline = 0; self.__stringvar = stringVar; self.__window = window
    def write(self,s):
        try:
            log_report = open(log_file,'a')
            log_report.write(s); log_report.close() ### Variable to record each print statement
            new = self.__stringvar.get()
            for c in s:
                #if c == '\n': self.__newline = 1
                if c == '\k': self.__newline = 1### This should not be found and thus results in a continous feed rather than replacing a single line
                else:
                    if self.__newline: new = ""; self.__newline = 0
                    new = new+c
            self.set(new)
        except Exception: pass
    def set(self,s):
        try: self.__stringvar.set(s); self.__window.update()
        except Exception: pass
    def get(self):
        try:
            return self.__stringvar.get()
        except Exception: pass
    def flush(self):
        pass
            
def timestamp():
    import datetime
    today = str(datetime.date.today()); today = string.split(today,'-'); today = today[0]+''+today[1]+''+today[2]
    time_stamp = string.replace(time.ctime(),':','')
    time_stamp = string.replace(time_stamp,'  ',' ')
    time_stamp = string.split(time_stamp,' ') ###Use a time-stamp as the output dir (minus the day)
    time_stamp = today+'-'+time_stamp[3]
    return time_stamp

def callWXPython():
    import wx
    import AltAnalyzeViewer
    app = wx.App(False)
    AltAnalyzeViewer.remoteViewer(app)
    
def AltAnalyzeSetup(skip_intro):
    global apt_location; global root_dir;global log_file; global summary_data_db; summary_data_db={}; reload(UI)
    global probability_statistic; global commandLineMode; commandLineMode = 'no'
    if 'remoteViewer' == skip_intro:
        if os.name == 'nt':
            callWXPython()
        elif os.name == 'ntX':
            package_path = filepath('python')
            win_package_path = string.replace(package_path,'python','AltAnalyzeViewer.exe')
            import subprocess
            subprocess.call([win_package_path]);sys.exit()
        elif os.name == 'posix':
            package_path = filepath('python')
            #mac_package_path = string.replace(package_path,'python','AltAnalyze.app/Contents/MacOS/python')
            #os.system(mac_package_path+' RemoteViewer.py');sys.exit()
            mac_package_path = string.replace(package_path,'python','AltAnalyzeViewer.app/Contents/MacOS/AltAnalyzeViewer')
            import subprocess
            subprocess.call([mac_package_path]);sys.exit()
        """
        import threading
        import wx
        app = wx.PySimpleApp()
        t = threading.Thread(target=callWXPython)
        t.setDaemon(1)
        t.start()

        s = 1
        queue = mlp.Queue()
        proc = mlp.Process(target=callWXPython) ### passing sys.stdout unfortunately doesn't work to pass the Tk string
        proc.start()
        sys.exit()
        """
    reload(UI)
    expr_var, alt_var, additional_var, goelite_var, exp_file_location_db = UI.getUserParameters(skip_intro,Multi=mlp)
    """except Exception:
        if 'SystemExit' not in str(traceback.format_exc()):
            expr_var, alt_var, additional_var, goelite_var, exp_file_location_db = UI.getUserParameters('yes')
        else: sys.exit()"""
        
    for dataset in exp_file_location_db:
        fl = exp_file_location_db[dataset]
        apt_location = fl.APTLocation()
        root_dir = fl.RootDir()
        try: probability_statistic = fl.ProbabilityStatistic()
        except Exception: probability_statistic = 'unpaired t-test'
    
    time_stamp = timestamp()    
    log_file = filepath(root_dir+'AltAnalyze_report-'+time_stamp+'.log')
    log_report = open(log_file,'w'); log_report.close()
    if use_Tkinter == 'yes' and debug_mode == 'no':
        try:
            global root; root = Tk()
            StatusWindow(root,expr_var, alt_var, goelite_var, additional_var, exp_file_location_db)
            root.destroy()
        except Exception, exception:
            try:
                print traceback.format_exc()
                badExit()
            except Exception: sys.exit()
    else: AltAnalyzeMain(expr_var, alt_var, goelite_var, additional_var, exp_file_location_db,'')

def badExit():
        print "\n...exiting AltAnalyze due to unexpected error"
        try: 
            time_stamp = timestamp()    
            print_out = "Unknown error encountered during data processing.\nPlease see logfile in:\n\n"+log_file+"\nand report to altanalyze@gmail.com."
            try:
                if len(log_file)>0:
                    if commandLineMode == 'no':
                        if os.name == 'nt':
                            try: os.startfile('"'+log_file+'"')
                            except Exception:  os.system('open "'+log_file+'"')
                        elif 'darwin' in sys.platform: os.system('open "'+log_file+'"')
                        elif 'linux' in sys.platform: os.system('xdg-open "'+log_file+'"')
                if commandLineMode == 'no':
                    try: UI.WarningWindow(print_out,'Error Encountered!'); root.destroy()
                    except Exception: print print_out
            except Exception: sys.exit()
        except Exception: sys.exit()
        sys.exit()
        kill
            
def AltAnalyzeMain(expr_var,alt_var,goelite_var,additional_var,exp_file_location_db,root):
  ### Hard-coded defaults
  w = 'Agilent'; x = 'Affymetrix'; y = 'Ensembl'; z = 'any'; data_source = y; constitutive_source = z; manufacturer = x ### Constitutive source, is only really paid attention to if Ensembl, otherwise Affymetrix is used (even if default)
  ### Get default options for ExpressionBuilder and AltAnalyze    
  start_time = time.time()
  test_goelite = 'no'; test_results_pannel = 'no'

  global species; global array_type; global expression_data_format; global use_R; use_R = 'no'
  global analysis_method; global p_threshold; global filter_probeset_types
  global permute_p_threshold; global perform_permutation_analysis; global export_NI_values
  global run_MiDAS; global analyze_functional_attributes;  global microRNA_prediction_method
  global calculate_normIntensity_p; global pathway_permutations; global avg_all_for_ss; global analyze_all_conditions
  global remove_intronic_junctions 

  global agglomerate_inclusion_probesets; global expression_threshold; global factor_out_expression_changes
  global only_include_constitutive_containing_genes; global remove_transcriptional_regulated_genes; global add_exons_to_annotations
  global exclude_protein_details; global filter_for_AS; global use_direct_domain_alignments_only; global run_from_scratch
  global explicit_data_type; explicit_data_type = 'null'
  global altanalyze_files; altanalyze_files = []
  
  species,array_type,manufacturer,constitutive_source,dabg_p,raw_expression_threshold,avg_all_for_ss,expression_data_format,include_raw_data, run_from_scratch, perform_alt_analysis = expr_var
  analysis_method,p_threshold,filter_probeset_types,alt_exon_fold_variable,gene_expression_cutoff,remove_intronic_junctions,permute_p_threshold,perform_permutation_analysis, export_NI_values, analyze_all_conditions = alt_var
  calculate_normIntensity_p, run_MiDAS, use_direct_domain_alignments_only, microRNA_prediction_method, filter_for_AS, additional_algorithms = additional_var
  ge_fold_cutoffs,ge_pvalue_cutoffs,ge_ptype,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,pathway_permutations,mod,returnPathways = goelite_var
  
  original_remove_intronic_junctions = remove_intronic_junctions
  if run_from_scratch == 'Annotate External Results': analysis_method = 'external'
  if returnPathways == 'no' or returnPathways == 'None':
    returnPathways = None

  for dataset in exp_file_location_db:
    fl = exp_file_location_db[dataset]
    try: exon_exp_threshold = fl.ExonExpThreshold()
    except Exception: exon_exp_threshold = 'NA'
    try: gene_exp_threshold = fl.GeneExpThreshold()
    except Exception: gene_exp_threshold = 'NA'
    try: exon_rpkm_threshold = fl.ExonRPKMThreshold()
    except Exception: exon_rpkm_threshold = 'NA'
    try: rpkm_threshold = fl.RPKMThreshold() ### Gene-Level
    except Exception: rpkm_threshold = 'NA'
    fl.setJunctionExpThreshold(raw_expression_threshold) ### For RNA-Seq, this specifically applies to exon-junctions
    try: predictGroups = fl.predictGroups()
    except Exception: predictGroups = False
    try:
        if fl.excludeLowExpressionExons(): excludeLowExpExons = 'yes'
        else: excludeLowExpExons = 'no'
    except Exception: excludeLowExpExons = 'no'
    
    if test_goelite == 'yes': ### It can be difficult to get error warnings from GO-Elite, unless run here
      results_dir = filepath(fl.RootDir())
      elite_input_dirs = ['AltExonConfirmed','AltExon','regulated','upregulated','downregulated'] ### Run GO-Elite multiple times to ensure heatmaps are useful and to better organize results
      for elite_dir in elite_input_dirs:
           file_dirs = results_dir+'GO-Elite/'+elite_dir,results_dir+'GO-Elite/denominator',results_dir+'GO-Elite/'+elite_dir
           variables = species,mod,pathway_permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,returnPathways,file_dirs,root
           GO_Elite.remoteAnalysis(variables,'non-UI',Multi=mlp) 

  global perform_element_permutation_analysis; global permutations
  perform_element_permutation_analysis = 'yes'; permutations = 2000
  analyze_functional_attributes = 'yes' ### Do this by default (shouldn't substantially increase runtime)
  
  if run_from_scratch != 'Annotate External Results' and (array_type != "3'array" and array_type!='RNASeq'):
      if run_from_scratch !='Process AltAnalyze filtered':
          try: raw_expression_threshold = float(raw_expression_threshold)
          except Exception: raw_expression_threshold = 1
          if raw_expression_threshold<1:
              raw_expression_threshold = 1
              print "Expression threshold < 1, forcing to be a minimum of 1."
          try: dabg_p = float(dabg_p)
          except Exception: dabg_p = 0
          if dabg_p == 0 or dabg_p > 1:
              print "Invalid dabg-p value threshold entered,(",dabg_p,") setting to default of 0.05"
              dabg_p = 0.05

  if use_direct_domain_alignments_only == 'direct-alignment': use_direct_domain_alignments_only = 'yes'
  if run_from_scratch == 'Process CEL files': expression_data_format = 'log'
  print "Beginning AltAnalyze Analysis... Format:", expression_data_format

    
  if array_type == 'RNASeq': id_name = 'exon/junction IDs'
  else: id_name = 'array IDs'

  print_items=[]; #print [permute_p_threshold]; sys.exit()
  print_items.append("AltAnalyze version 2.1.2 - Expression Analysis Parameters Being Used...")
  print_items.append('\t'+'database'+': '+unique.getCurrentGeneDatabaseVersion())
  print_items.append('\t'+'species'+': '+species)
  print_items.append('\t'+'method'+': '+array_type)
  print_items.append('\t'+'manufacturer'+': '+manufacturer)
  print_items.append('\t'+'probability_statistic'+': '+probability_statistic)
  print_items.append('\t'+'constitutive_source'+': '+constitutive_source)
  print_items.append('\t'+'dabg_p'+': '+str(dabg_p))
  if array_type == 'RNASeq':
    print_items.append('\t'+'junction expression threshold'+': '+str(raw_expression_threshold))
    print_items.append('\t'+'exon_exp_threshold'+': '+str(exon_exp_threshold))
    print_items.append('\t'+'gene_exp_threshold'+': '+str(gene_exp_threshold))
    print_items.append('\t'+'exon_rpkm_threshold'+': '+str(exon_rpkm_threshold))
    print_items.append('\t'+'gene_rpkm_threshold'+': '+str(rpkm_threshold))
    print_items.append('\t'+'exclude low expressing exons for RPKM'+': '+excludeLowExpExons)
  else:
    print_items.append('\t'+'raw_expression_threshold'+': '+str(raw_expression_threshold))
  print_items.append('\t'+'avg_all_for_ss'+': '+avg_all_for_ss)
  print_items.append('\t'+'expression_data_format'+': '+expression_data_format)
  print_items.append('\t'+'include_raw_data'+': '+include_raw_data)
  print_items.append('\t'+'run_from_scratch'+': '+run_from_scratch)
  print_items.append('\t'+'perform_alt_analysis'+': '+perform_alt_analysis)
  if avg_all_for_ss == 'yes': cs_type = 'core'
  else: cs_type = 'constitutive'
  print_items.append('\t'+'calculate_gene_expression_using'+': '+cs_type)
  
  print_items.append("Alternative Exon Analysis Parameters Being Used..." )
  print_items.append('\t'+'analysis_method'+': '+analysis_method)
  print_items.append('\t'+'p_threshold'+': '+str(p_threshold))
  print_items.append('\t'+'filter_data_types'+': '+filter_probeset_types)
  print_items.append('\t'+'alt_exon_fold_variable'+': '+str(alt_exon_fold_variable))
  print_items.append('\t'+'gene_expression_cutoff'+': '+str(gene_expression_cutoff))
  print_items.append('\t'+'remove_intronic_junctions'+': '+remove_intronic_junctions)
  
  print_items.append('\t'+'avg_all_for_ss'+': '+avg_all_for_ss)
  print_items.append('\t'+'permute_p_threshold'+': '+str(permute_p_threshold))
  print_items.append('\t'+'perform_permutation_analysis'+': '+perform_permutation_analysis)
  print_items.append('\t'+'export_NI_values'+': '+export_NI_values)
  print_items.append('\t'+'run_MiDAS'+': '+run_MiDAS)
  print_items.append('\t'+'use_direct_domain_alignments_only'+': '+use_direct_domain_alignments_only)
  print_items.append('\t'+'microRNA_prediction_method'+': '+microRNA_prediction_method)
  print_items.append('\t'+'analyze_all_conditions'+': '+analyze_all_conditions)
  print_items.append('\t'+'filter_for_AS'+': '+filter_for_AS)
  
  if pathway_permutations == 'NA': run_GOElite = 'decide_later'
  else: run_GOElite = 'run-immediately'
  print_items.append('\t'+'run_GOElite'+': '+ run_GOElite)

  universalPrintFunction(print_items)

  if commandLineMode == 'yes': print 'Running command line mode:',commandLineMode

  summary_data_db['gene_assayed'] = 0
  summary_data_db['denominator_exp_genes']=0
  summary_data_db['alt_events'] = 0
  summary_data_db['denominator_exp_events'] = 0
  summary_data_db['alt_genes'] = 0
  summary_data_db['direct_domain_genes'] = 0
  summary_data_db['miRNA_gene_denom'] = 0
  summary_data_db['miRNA_gene_hits'] = 0
          
  if test_results_pannel == 'yes': ### It can be difficult to get error warnings from GO-Elite, unless run here
      graphic_links = []
      graphic_links.append(['test','Config/AltAnalyze_structure-RNASeq.jpg'])
      summary_data_db['QC']=graphic_links
      print_out = 'Analysis complete. AltAnalyze results\nexported to "AltResults/AlternativeOutput".'
      dataset = 'test'; results_dir=''
      print "Analysis Complete\n";
                
      if root !='' and root !=None:
          UI.InfoWindow(print_out,'Analysis Completed!')
          tl = Toplevel(); SummaryResultsWindow(tl,'GE',results_dir,dataset,'parent',summary_data_db)
      root.destroy(); sys.exit()
      
  global export_go_annotations; global aspire_output_list; global aspire_output_gene_list
  global filter_probesets_by; global global_addition_factor; global onlyAnalyzeJunctions
  global log_fold_cutoff; global aspire_cutoff; global annotation_system; global alt_exon_logfold_cutoff

  """dabg_p = 0.75; data_type = 'expression' ###used for expression analysis when dealing with AltMouse arrays
  a = "3'array"; b = "exon"; c = "AltMouse"; e = "custom"; array_type = c
  l = 'log'; n = 'non-log'; expression_data_format = l
  hs = 'Hs'; mm = 'Mm'; dr = 'Dr'; rn = 'Rn'; species = mm
  include_raw_data = 'yes'; expression_threshold = 70 ### Based on suggestion from BMC Genomics. 2006 Dec 27;7:325. PMID: 17192196, for hu-exon 1.0 st array
  avg_all_for_ss = 'no'  ###Default is 'no' since we don't want all probes averaged for the exon arrays"""
    
  ###### Run ExpressionBuilder ######
  """ExpressionBuilder is used to:
  (1) extract out gene expression values, provide gene annotations, and calculate summary gene statistics
  (2) filter probesets based DABG p-values and export to pair-wise comparison files
  (3) build array annotations files matched to gene structure features (e.g. exons, introns) using chromosomal coordinates
  options 1-2 are executed in remoteExpressionBuilder and option 3 is by running ExonArrayEnsembl rules"""

  try:
    additional_algorithm = additional_algorithms.Algorithm()
    additional_score = additional_algorithms.Score()
  except Exception: additional_algorithm = 'null'; additional_score = 'null'
  
  if analysis_method == 'FIRMA': analyze_metaprobesets = 'yes'
  elif additional_algorithm == 'FIRMA': analyze_metaprobesets = 'yes'
  else: analyze_metaprobesets = 'no'

  ### Check to see if this is a real or FAKE (used for demonstration purposes) dataset
  if run_from_scratch == 'Process CEL files' or 'Feature Extraction' in run_from_scratch: 
      for dataset in exp_file_location_db:
        if run_from_scratch == 'Process CEL files':
            fl = exp_file_location_db[dataset]
            pgf_file=fl.InputCDFFile()
            results_dir = filepath(fl.RootDir())
            if '_demo' in pgf_file: ### Thus we are running demo CEL files and want to quit immediately
              print_out = 'Analysis complete. AltAnalyze results\nexported to "AltResults/AlternativeOutput".'
              try:      
                  print "Analysis Complete\n";
                  if root !='' and root !=None:
                      UI.InfoWindow(print_out,'Analysis Completed!')
                      tl = Toplevel(); SummaryResultsWindow(tl,'AS',results_dir,dataset,'parent',summary_data_db)
              except Exception: null=[]
              skip_intro = 'yes'
              if pathway_permutations == 'NA' and run_from_scratch != 'Annotate External Results':
                  reload(UI)
                  UI.getUpdatedParameters(array_type,species,run_from_scratch,results_dir)
              try: AltAnalyzeSetup('no')
              except Exception: sys.exit()
      if 'CEL files' in run_from_scratch:
        from import_scripts import APT
        try:
            try:
                APT.probesetSummarize(exp_file_location_db,analyze_metaprobesets,filter_probeset_types,species,root)
                if analyze_metaprobesets == 'yes':
                    analyze_metaprobesets = 'no' ### Re-run the APT analysis to obtain probeset rather than gene-level results (only the residuals are needed from a metaprobeset run)
                    APT.probesetSummarize(exp_file_location_db,analyze_metaprobesets,filter_probeset_types,species,root)
            except Exception:
                    import platform
                    print "Trying to change APT binary access privileges"
                    for dataset in exp_file_location_db: ### Instance of the Class ExpressionFileLocationData
                        fl = exp_file_location_db[dataset]; apt_dir =fl.APTLocation()
                    if '/bin' in apt_dir: apt_file = apt_dir +'/apt-probeset-summarize' ### if the user selects an APT directory
                    elif os.name == 'nt': apt_file = apt_dir + '/PC/'+platform.architecture()[0]+'/apt-probeset-summarize.exe'
                    elif 'darwin' in sys.platform: apt_file = apt_dir + '/Mac/apt-probeset-summarize'
                    elif 'linux' in sys.platform:
                        if '32bit' in platform.architecture(): apt_file = apt_dir + '/Linux/32bit/apt-probeset-summarize'
                        elif '64bit' in platform.architecture(): apt_file = apt_dir + '/Linux/64bit/apt-probeset-summarize'
                        apt_file = filepath(apt_file)
                    os.chmod(apt_file,0777)
                    midas_dir = string.replace(apt_file,'apt-probeset-summarize','apt-midas')
                    os.chmod(midas_dir,0777)
                    APT.probesetSummarize(exp_file_location_db,analysis_method,filter_probeset_types,species,root)
        except Exception:
          print_out = 'AltAnalyze encountered an un-expected error while running Affymetrix\n'
          print_out += 'Power Tools (APT). Additional information may be found in the directory\n'
          print_out += '"ExpressionInput/APT" in the output directory. You may also encounter issues\n'
          print_out += 'if you are logged into an account with restricted priveledges.\n\n'
          print_out += 'If this issue can not be resolved, contact AltAnalyze help or run RMA outside\n'
          print_out += 'of AltAnalyze and import the results using the analysis option "expression file".\n'
          print traceback.format_exc()
          try:
              UI.WarningWindow(print_out,'Exit')
              root.destroy(); sys.exit()
          except Exception:
              print print_out; sys.exit()
      elif 'Feature Extraction' in run_from_scratch:
        from import_scripts import ProcessAgilentArrays
        try: ProcessAgilentArrays.agilentSummarize(exp_file_location_db)
        except Exception:
            print_out = 'Agilent array import and processing failed... see error log for details...'
            print traceback.format_exc()
            try:
                UI.WarningWindow(print_out,'Exit')
                root.destroy(); sys.exit()
            except Exception:
                print print_out; sys.exit()
        reload(ProcessAgilentArrays)

  if run_from_scratch == 'Process RNA-seq reads' or run_from_scratch == 'buildExonExportFiles':
      import RNASeq; reload(RNASeq); import RNASeq
      for dataset in exp_file_location_db: fl = exp_file_location_db[dataset]

      ### The below function aligns splice-junction coordinates to Ensembl exons from BED Files and
      ### exports AltAnalyze specific databases that are unique to this dataset to the output directory

      try: fastq_folder = fl.RunKallisto()
      except Exception: print traceback.format_exc()
      try: customFASTA = fl.CustomFASTA()
      except Exception: customFASTA = None
      processBEDfiles = True
      if len(fastq_folder)>0:
        ### Perform pseudoalignment with Kallisto on FASTQ files
        processBEDfiles=False
        try:
            RNASeq.runKallisto(species,dataset,root_dir,fastq_folder,mlp,returnSampleNames=False,customFASTA=customFASTA)
            biotypes = 'ran'
            dir_list = unique.read_directory(root_dir)
            ### If we are performing a splicing analysis
            if perform_alt_analysis != 'no' and perform_alt_analysis != 'expression':
                print '...Performing analyses on junction-RPKM versus Kallisto-TPM.'
                for file in dir_list:
                    if '.bam' in string.lower(file):
                        processBEDfiles=True
                    if '.bed' in string.lower(file):
                        processBEDfiles=True
                try: rpkm_threshold = fl.RPKMThreshold()
                except Exception: rpkm_threshold = []
                if isinstance(rpkm_threshold, int) ==False:
                    array_type = 'RNASeq'
                    fl.setArrayType(array_type)
                    fl.setBEDFileDir(root_dir)
                    fl.setRPKMThreshold(1.0)
                    fl.setExonExpThreshold(5.0)
                    fl.setGeneExpThreshold(200.0)
                    fl.setExonRPKMThreshold(0.5)
                    fl.setJunctionExpThreshold(5.0)
                    fl.setVendor('RNASeq')
                from import_scripts import BAMtoJunctionBED
                ### Export BAM file indexes
                try: BAMtoJunctionBED.exportIndexes(root_dir)
                except:
                    print 'BAM file indexing failed...'
                    print traceback.format_exc()
            else:
                print '...Performing analyses on Kallisto-TPM values directly.'
                array_type = "3'array"
                fl.setArrayType(array_type)
                vendor = 'other:Ensembl' ### Ensembl linked system name
                fl.setVendor(vendor)
        except Exception:
            print traceback.format_exc()
            biotypes='failed'
        
      if processBEDfiles:
          analyzeBAMs = False; bedFilesPresent = False
          dir_list = unique.read_directory(fl.BEDFileDir())
          for file in dir_list:
            if '.bam' in string.lower(file):
                analyzeBAMs=True
            if '.bed' in string.lower(file):
                bedFilesPresent=True

          if analyzeBAMs and bedFilesPresent==False:
            from import_scripts import multiBAMtoBED
            bam_dir = fl.BEDFileDir()
            refExonCoordinateFile = filepath('AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt')
            outputExonCoordinateRefBEDfile = bam_dir+'/BedRef/'+species+'_'+string.replace(dataset,'exp.','')
            analysisType = ['exon','junction','reference']
            #analysisType = ['junction']
            #print [fl.multiThreading()]
            multiBAMtoBED.parallelBAMProcessing(bam_dir,refExonCoordinateFile,outputExonCoordinateRefBEDfile,analysisType=analysisType,useMultiProcessing=fl.multiThreading(),MLP=mlp,root=root)
    
          biotypes = RNASeq.alignExonsAndJunctionsToEnsembl(species,exp_file_location_db,dataset,Multi=mlp)
      
      if biotypes == 'failed':
          print_out = 'No valid chromosomal positions in the input BED or BioScope files. Exiting AltAnalyze.'
          if len(fastq_folder)>0:
            if 'FTP' in traceback.format_exc():
                print_out = 'AltAnlayze was unable to retreive a transcript fasta sequence file from the Ensembl website. '
                print_out += 'Ensure you are connected to the internet and that the website http://ensembl.org is live.'
            else:
                print_out = 'An unexplained error was encountered with Kallisto analysis:\n'
                print_out += traceback.format_exc()
          try:
              UI.WarningWindow(print_out,'Exit')
              root.destroy(); sys.exit()
          except Exception:
              print print_out; sys.exit()
            
      reload(RNASeq)

      if root_dir in biotypes:
            print_out = 'Exon-level BED coordinate predictions exported to:\n'+biotypes
            print_out+= '\n\nAfter obtaining exon expression estimates, rename exon BED files to\n'
            print_out+= 'match the junction name (e.g., Sample1__exon.bed and Sample1__junction.bed)\n'
            print_out+= 'and re-run AltAnalyze (see tutorials at http://altanalyze.org for help).'
            UI.InfoWindow(print_out,'Export Complete')
            try: root.destroy(); sys.exit()
            except Exception: sys.exit()

  if predictGroups == True:
    expFile = fl.ExpFile()
    if array_type == 'RNASeq':
        exp_threshold=100; rpkm_threshold=10
    else:
        exp_threshold=200; rpkm_threshold=8
    RNASeq.singleCellRNASeqWorkflow(species, array_type, expFile, mlp, exp_threshold=exp_threshold, rpkm_threshold=rpkm_threshold)
    
  goelite_run = False
  if run_from_scratch == 'Process Expression file' or run_from_scratch == 'Process CEL files' or run_from_scratch == 'Process RNA-seq reads' or 'Feature Extraction' in run_from_scratch:
      if (fl.NormMatrix()=='quantile' or fl.NormMatrix()=='group') and 'Feature Extraction' not in run_from_scratch:
            from stats_scripts import NormalizeDataset
            try: NormalizeDataset.normalizeDataset(fl.ExpFile(),normalization=fl.NormMatrix(),platform=array_type)
            except Exception: print "Normalization failed for unknown reasons..."
      #"""
      status = ExpressionBuilder.remoteExpressionBuilder(species,array_type,
            dabg_p,raw_expression_threshold,avg_all_for_ss,expression_data_format,
            manufacturer,constitutive_source,data_source,include_raw_data,
            perform_alt_analysis,ge_fold_cutoffs,ge_pvalue_cutoffs,ge_ptype,
            exp_file_location_db,root)
      reload(ExpressionBuilder) ### Clears Memory
      #"""

      graphics=[]
      if fl.MarkerFinder() == 'yes':
            ### Identify putative condition-specific marker genees
            import markerFinder
            fl.setOutputDir(root_dir) ### This needs to be set here
            exp_file = fl.ExpFile()
            if array_type != "3'array": exp_file = string.replace(exp_file,'.txt','-steady-state.txt')
            markerFinder_inputs = [exp_file,fl.DatasetFile()] ### Output a replicate and non-replicate version
            markerFinder_inputs = [exp_file] ### Only considers the replicate and not mean analysis (recommended)
            for input_exp_file in markerFinder_inputs:
                ### This applies to an ExpressionOutput DATASET file compoosed of gene expression values (averages already present)
                try:
                    output_dir = markerFinder.getAverageExpressionValues(input_exp_file,array_type) ### Either way, make an average annotated file from the DATASET file
                except Exception:
                    print "Unknown MarkerFinder failure (possible filename issue or data incompatibility)..."
                    print traceback.format_exc()
                    continue
                if 'DATASET' in input_exp_file:
                    group_exp_file = string.replace(input_exp_file,'DATASET','AVERAGE')
                else:
                    group_exp_file = (input_exp_file,output_dir) ### still analyze the primary sample
                compendiumType = 'protein_coding'
                if expression_data_format == 'non-log': logTransform = True
                else: logTransform = False
                try: markerFinder.analyzeData(group_exp_file,species,array_type,compendiumType,AdditionalParameters=fl,logTransform=logTransform)
                except Exception: None
            ### Generate heatmaps (unclustered - order by markerFinder)
            try: graphics = markerFinder.generateMarkerHeatMaps(fl,array_type,graphics=graphics,Species=species)
            except Exception: print traceback.format_exc()

      remove_intronic_junctions = original_remove_intronic_junctions ### This var gets reset when running FilterDABG
      try:
           summary_data_db['QC'] = fl.GraphicLinks()+graphics ### provides links for displaying QC and clustering plots
      except Exception:
           null=[] ### Visualization support through matplotlib either not present or visualization options excluded
      #print '!!!!!finished expression builder'
      #returnLargeGlobalVars()
      expression_data_format = 'log' ### This variable is set from non-log in FilterDABG when present (version 1.16)
      try:
            parent_dir = fl.RootDir()+'/GO-Elite/regulated/'
            dir_list = read_directory(parent_dir)
            for file in dir_list:
                input_file_dir = parent_dir+'/'+file
                inputType = 'IDs'
                interactionDirs = ['WikiPathways','KEGG','BioGRID','TFTargets']
                output_dir = parent_dir
                degrees = 'direct'
                input_exp_file = input_file_dir
                gsp = UI.GeneSelectionParameters(species,array_type,manufacturer)
                gsp.setGeneSet('None Selected')
                gsp.setPathwaySelect('')
                gsp.setGeneSelection('')
                gsp.setOntologyID('')
                gsp.setIncludeExpIDs(True)
                UI.networkBuilder(input_file_dir,inputType,output_dir,interactionDirs,degrees,input_exp_file,gsp,'')
      except Exception:
        print traceback.format_exc()

                            
      if status == 'stop':
          ### See if the array and species are compatible with GO-Elite analysis
          system_codes = UI.getSystemInfo()
          go_elite_analysis_supported = 'yes'
          species_names = UI.getSpeciesInfo()
          for dataset in exp_file_location_db:
              fl = exp_file_location_db[dataset]; results_dir = filepath(fl.RootDir())

          ### Perform GO-Elite Analysis
          if pathway_permutations != 'NA':
              try:
                  print '\nBeginning to run GO-Elite analysis on alternative exon results'
                  elite_input_dirs = ['AltExonConfirmed','AltExon','regulated','upregulated','downregulated'] ### Run GO-Elite multiple times to ensure heatmaps are useful and to better organize results
                  for elite_dir in elite_input_dirs:
                        file_dirs = results_dir+'GO-Elite/'+elite_dir,results_dir+'GO-Elite/denominator',results_dir+'GO-Elite/'+elite_dir
                        input_dir = results_dir+'GO-Elite/'+elite_dir
                        variables = species,mod,pathway_permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,returnPathways,file_dirs,root
                        try: input_files = read_directory(input_dir) ### Are there any files to analyze?
                        except Exception: input_files=[]
                        if len(input_files)>0:
                            try: GO_Elite.remoteAnalysis(variables,'non-UI',Multi=mlp); goelite_run = True
                            except Exception,e:
                              print e
                              print "GO-Elite analysis failed"
                            try: GO_Elite.moveMAPPFinderFiles(file_dirs[0])
                            except Exception: print 'Input GO-Elite files could NOT be moved.'
                            try: GO_Elite.moveMAPPFinderFiles(file_dirs[1])
                            except Exception: print 'Input GO-Elite files could NOT be moved.'
              except Exception: pass
              if goelite_run == False:
                    print 'No GO-Elite input files to analyze (check your criterion).'
          print_out = 'Analysis complete. Gene expression\nsummary exported to "ExpressionOutput".'
          try:
              if use_Tkinter == 'yes':
                  print "Analysis Complete\n"; UI.InfoWindow(print_out,'Analysis Completed!')
                  tl = Toplevel(); SummaryResultsWindow(tl,'GE',results_dir,dataset,'parent',summary_data_db)
                  if pathway_permutations == 'NA' and run_from_scratch != 'Annotate External Results':                
                      if go_elite_analysis_supported == 'yes': 
                          UI.getUpdatedParameters(array_type,species,run_from_scratch,file_dirs)
                  try: AltAnalyzeSetup('no')
                  except Exception:
                      print traceback.format_exc()
                      sys.exit()
              else:  print '\n'+print_out; sys.exit()
          except Exception:
              #print 'Failed to report status through GUI.'
              sys.exit()
      else: altanalyze_files = status[1] ### These files are the comparison files to analyze
  elif run_from_scratch == 'update DBs':
      null=[] ###Add link to new module here (possibly)
      #updateDBs(species,array_type)
      sys.exit()
  if perform_alt_analysis != 'expression': ###Thus perform_alt_analysis = 'both' or 'alt' (default when skipping expression summary step)
      ###### Run AltAnalyze ######  
      global dataset_name; global summary_results_db; global summary_results_db2
      summary_results_db={}; summary_results_db2={}; aspire_output_list=[]; aspire_output_gene_list=[]

      onlyAnalyzeJunctions = 'no'; agglomerate_inclusion_probesets = 'no'; filter_probesets_by = 'NA'
      if array_type == 'AltMouse' or ((array_type == 'junction' or array_type == 'RNASeq') and explicit_data_type == 'null'):
          if filter_probeset_types == 'junctions-only':  onlyAnalyzeJunctions = 'yes'
          elif filter_probeset_types == 'combined-junctions': agglomerate_inclusion_probesets = 'yes'; onlyAnalyzeJunctions = 'yes'
          elif filter_probeset_types == 'exons-only': analysis_method = 'splicing-index'; filter_probesets_by = 'exon'
          if filter_probeset_types == 'combined-junctions' and array_type == 'junction' or array_type == 'RNASeq': filter_probesets_by = 'all'
      else: filter_probesets_by = filter_probeset_types

      c = 'Ensembl'; d = 'Entrez Gene'
      annotation_system = c
      expression_threshold = 0 ###This is different than the raw_expression_threshold (probably shouldn't filter so set to 0)

      if analysis_method == 'linearregres-rlm': analysis_method = 'linearregres';use_R = 'yes'

      if gene_expression_cutoff<1:
          gene_expression_cutoff = 2 ### A number less than one is invalid
          print "WARNING!!!! Invalid gene expression fold cutoff entered,\nusing the default value of 2, must be greater than 1."
      log_fold_cutoff = math.log(float(gene_expression_cutoff),2)
    
      if analysis_method != 'ASPIRE' and analysis_method != 'none' and analysis_method != 'MultiPath-PSI':
          if p_threshold <= 0 or p_threshold >1:
              p_threshold = 0.05 ### A number less than one is invalid
              print "WARNING!!!! Invalid alternative exon p-value threshold entered,\nusing the default value of 0.05."              
          if alt_exon_fold_variable<1: 
              alt_exon_fold_variable = 1 ### A number less than one is invalid
              print "WARNING!!!! Invalid alternative exon fold cutoff entered,\nusing the default value of 2, must be greater than 1."
          try: alt_exon_logfold_cutoff = math.log(float(alt_exon_fold_variable),2)
          except Exception: alt_exon_logfold_cutoff = 1
      else: alt_exon_logfold_cutoff = float(alt_exon_fold_variable)
      
      global_addition_factor = 0
      export_junction_comparisons = 'no' ### No longer accessed in this module - only in update mode through a different module
      
      factor_out_expression_changes = 'yes' ### Use 'no' if data is normalized already or no expression normalization for ASPIRE desired
      only_include_constitutive_containing_genes = 'yes'
      remove_transcriptional_regulated_genes = 'yes'
      add_exons_to_annotations = 'no'
      exclude_protein_details = 'no'

      if analysis_method == 'ASPIRE' or 'linearregres' in analysis_method: annotation_system = d
      if 'linear' in analysis_method: analysis_method = 'linearregres'
      if 'aspire' in analysis_method: analysis_method = 'ASPIRE'
      if array_type == 'AltMouse': species = 'Mm'
      #if export_NI_values == 'yes': remove_transcriptional_regulated_genes = 'no'

      ###Saves run-time while testing the software (global variable stored)
      #import_dir = '/AltDatabase/affymetrix/'+species
      #dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
      ### Get Ensembl-GO and pathway annotations from GO-Elite files
      universalPrintFunction(["Importing GO-Elite pathway/GO annotations"])
      global go_annotations; go_annotations={}
      from import_scripts import BuildAffymetrixAssociations
      go_annotations = BuildAffymetrixAssociations.getEnsemblAnnotationsFromGOElite(species)
      global probeset_annotations_file
      if array_type == 'RNASeq': probeset_annotations_file = root_dir+'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_junctions.txt'
      elif array_type == 'AltMouse': probeset_annotations_file = 'AltDatabase/'+species+'/'+array_type+'/'+'MASTER-probeset-transcript.txt'
      else: probeset_annotations_file = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_probesets.txt'

      #"""
      if analysis_method != 'none' and  analysis_method != 'MultiPath-PSI':
          analysis_summary = RunAltAnalyze() ### Only run if analysis methods is specified (only available for RNA-Seq and junction analyses)
      else: analysis_summary = None

      if analysis_summary != None:
          summary_results_db, aspire_output_gene_list, number_events_analyzed = analysis_summary
          summary_data_db2 = copy.deepcopy(summary_data_db)
          for i in summary_data_db2: del summary_data_db[i] ### If we reset the variable it violates it's global declaration... do this instead
          #universalPrintFunction(['Alternative Exon Results for Junction Comparisons:'])
          #for i in summary_data_db: universalPrintFunction([i+' '+ str(summary_data_db[i])])
          exportSummaryResults(summary_results_db,analysis_method,aspire_output_list,aspire_output_gene_list,annotate_db,array_type,number_events_analyzed,root_dir)
      else:
          ### Occurs for RNASeq when no junctions are present
          summary_data_db2={}
       
      if array_type == 'junction' or array_type == 'RNASeq':
          #Reanalyze junction array data separately for individual probests rather than recipricol junctions
          if array_type == 'junction': explicit_data_type = 'exon'
          elif array_type == 'RNASeq': explicit_data_type = 'junction'
          else: report_single_probeset_results = 'no'
          ### Obtain exon analysis defaults
          expr_defaults, alt_exon_defaults, functional_analysis_defaults, goelite_defaults = UI.importDefaults('exon',species)
          analysis_method, null, filter_probeset_types, null, null, alt_exon_fold_variable, null, null, null, null, null, null, null, calculate_normIntensity_p, null = alt_exon_defaults

          filter_probesets_by = filter_probeset_types
              
          if additional_algorithm == 'splicing-index' or additional_algorithm == 'FIRMA':
              analysis_method = additional_algorithm
              #print [analysis_method], [filter_probeset_types], [p_threshold], [alt_exon_fold_variable]
              try: alt_exon_logfold_cutoff = math.log(float(additional_score),2)
              except Exception: alt_exon_logfold_cutoff = 1
              agglomerate_inclusion_probesets = 'no'
              try:
                summary_results_db, aspire_output_gene_list, number_events_analyzed = RunAltAnalyze()
                exportSummaryResults(summary_results_db,analysis_method,aspire_output_list,aspire_output_gene_list,annotate_db,'exon',number_events_analyzed,root_dir)
                if len(summary_data_db2)==0: summary_data_db2 = summary_data_db; explicit_data_type = 'exon-only'
                #universalPrintFunction(['Alternative Exon Results for Individual Probeset Analyses:'])
                #for i in summary_data_db: universalPrintFunction([i+' '+ str(summary_data_db[i])])
              except Exception:
                print traceback.format_exc() 
                None
      #""" 
      ### Perform dPSI Analysis
      try:
        if 'counts.' in fl.CountsFile(): pass
        else:
            dir_list = read_directory(fl.RootDir()+'ExpressionInput')
            for file in dir_list:
                if 'exp.' in file and 'steady-state' not in file:
                    fl.setExpFile(fl.RootDir()+'ExpressionInput/'+file)
                    #print [fl.RootDir()+'ExpressionInput/'+file]
                if 'groups.' in file:
                    fl.setGroupsFile(search_dir+'/'+file)
      except Exception:
            search_dir = fl.RootDir()+'/ExpressionInput'
            files = unique.read_directory(fl.RootDir()+'/ExpressionInput')
            for file in files:
                if 'exp.' in file and 'steady-state.txt' not in file:
                    fl.setExpFile(search_dir+'/'+file)
                if 'groups.' in file:
                    fl.setGroupsFile(search_dir+'/'+file)
      try:
          #"""
          try:
             #"""
             graphic_links2,cluster_input_file=ExpressionBuilder.unbiasedComparisonSpliceProfiles(fl.RootDir(),
                    species,array_type,expFile=fl.CountsFile(),min_events=1,med_events=1)
             #"""
             from import_scripts import AugmentEventAnnotations
             psi_annotated = AugmentEventAnnotations.parse_junctionfiles(fl.RootDir()+'/AltResults/AlternativeOutput/',species,array_type) ### Added in 2.1.1 - adds cassette and domains annotations
          except Exception:
            print traceback.format_exc()
            pass
          #"""
          inputpsi = fl.RootDir()+'AltResults/AlternativeOutput/'+species+'_'+array_type+'_top_alt_junctions-PSI-clust.txt'
          
          useAdvancedMetaDataAnalysis = True
          ### Calculate ANOVA p-value stats based on groups
          if array_type !='gene' and array_type != 'exon':
                if useAdvancedMetaDataAnalysis:
                    from stats_scripts import metaDataAnalysis
                    if ge_ptype == 'adjp':
                        use_adjusted_pval = True
                    else:
                        use_adjusted_pval = False
                    try:
                        log_fold_cutoff = float(alt_exon_fold_variable)
                        if log_fold_cutoff == 0.1:
                            log_fold_cutoff = 0.1 ### For significant digits
                    except: log_fold_cutoff = 0.1
                    try:
                        if p_threshold <= 0 or p_threshold >1:
                            pvalThreshold = 0.05 ### A number less than one is invalid
                        else: pvalThreshold = p_threshold
                    except: pvalThreshold = ge_pvalue_cutoffs ### Use the gene expression p-value cutoff if NA
                    try:
                        graphics_alt = metaDataAnalysis.remoteAnalysis(species,psi_annotated,fl.GroupsFile(),
                                platform='PSI',log_fold_cutoff=0.1,use_adjusted_pval=use_adjusted_pval,
                                pvalThreshold=ge_pvalue_cutoffs)
                        try: summary_data_db['QC'] += graphics_alt
                        except Exception: summary_data_db['QC'] = graphics_alt
                        try: summary_data_db['QC'] += graphic_links2
                        except Exception: pass
                        
                    except Exception:
                        print traceback.format_exc()
                else:
                    matrix,compared_groups,original_data = statistics.matrixImport(inputpsi)
                    matrix_pvalues=statistics.runANOVA(inputpsi,matrix,compared_groups)
                    significantFilteredDir = statistics.returnANOVAFiltered(inputpsi,original_data,matrix_pvalues)
                    graphic_link1 = ExpressionBuilder.exportHeatmap(significantFilteredDir)
                    try: summary_data_db['QC']+=graphic_link1
                    except Exception: summary_data_db['QC']=graphic_link1
      except Exception:
        print traceback.format_exc()
      
      import RNASeq
      try:
        graphic_link = RNASeq.compareExonAndJunctionResults(species,array_type,summary_results_db,root_dir)
        try: summary_data_db['QC']+=graphic_link
        except Exception: summary_data_db['QC']=graphic_link
      except Exception:
        print traceback.format_exc()
    
      #""" 
      ### Export the top 15 spliced genes
      try:
            altresult_dir = fl.RootDir()+'/AltResults/'
            splicing_results_root = altresult_dir+'/Clustering/'
            dir_list = read_directory(splicing_results_root)
            gene_string=''
            altanalyze_results_folder = altresult_dir+'/RawSpliceData/'+species
            ### Lookup the raw expression dir
            expression_results_folder = string.replace(altresult_dir,'AltResults','ExpressionInput')
            expression_dir = UI.getValidExpFile(expression_results_folder)
      
            show_introns=False
            try: altresult_dir = UI.getValidSplicingScoreFile(altanalyze_results_folder)
            except Exception,e:
                  print traceback.format_exc()
            analysisType='plot'
            for file in dir_list:
                  if 'AltExonConfirmed' in file:
                      gene_dir = splicing_results_root+'/'+file
                      genes = UI.importGeneList(gene_dir,limit=50) ### list of gene IDs or symbols
                      gene_string = gene_string+','+genes
                      print 'Imported genes from',file,'\n'
                      analysisType='plot'
            for file in dir_list:
                  if 'Combined-junction-exon-evidence' in file and 'top' not in file:
                      gene_dir = splicing_results_root+'/'+file
                      try: isoform_dir = UI.exportJunctionList(gene_dir,limit=50) ### list of gene IDs or symbols
                      except Exception: print traceback.format_exc()

            UI.altExonViewer(species,array_type,expression_dir, gene_string, show_introns, analysisType, None); print 'completed'
            UI.altExonViewer(species,array_type,altresult_dir, gene_string, show_introns, analysisType, None); print 'completed'
      except Exception:
        #print traceback.format_exc()
        pass
        
      if array_type != 'exon' and array_type != 'gene':
            ### SashimiPlot Visualization
            try:
              expression_results_folder = fl.RootDir()+'/ExpressionInput/'
              expression_dir = UI.getValidExpFile(expression_results_folder)
              show_introns=False
              analysisType='plot'
              top_PSI_junction = inputpsi
              #isoform_dir2 = UI.exportJunctionList(top_PSI_junction,limit=50) ### list of gene IDs or symbols
              altoutput_dir = export.findParentDir(top_PSI_junction)
              isoform_dir2 = altoutput_dir+'/top'+str(50)+'/MultiPath-PSI.txt'
              gene_string = UI.importGeneList(isoform_dir2,limit=50)
              UI.altExonViewer(species,array_type,expression_dir, gene_string, show_introns, analysisType, None); print 'completed'
              UI.altExonViewer(species,array_type,altresult_dir, gene_string, show_introns, analysisType, None); print 'completed'
            except Exception:
              print traceback.format_exc()
            try:
                analyzeBAMs = False
                dir_list = unique.read_directory(fl.RootDir())
                for file in dir_list:
                  if '.bam' in string.lower(file):
                      analyzeBAMs=True
                if analyzeBAMs:
                  ### Create sashimi plot index
                  from visualization_scripts import SashimiIndex
                  SashimiIndex.remoteIndexing(species,fl)
                  from visualization_scripts import SashimiPlot
                  print 'Exporting Sashimi Plots for the top-predicted splicing events... be patient'
                  try: SashimiPlot.remoteSashimiPlot(species,fl,fl.RootDir(),isoform_dir) ### assuming the bam files are in the root-dir
                  except Exception: pass # print traceback.format_exc()
                  print 'completed'
                  try: SashimiPlot.remoteSashimiPlot(species,fl,fl.RootDir(),isoform_dir2) ### assuming the bam files are in the root-dir
                  except Exception: pass #print traceback.format_exc()
                  print 'completed'
                  ### Try again, in case the symbol conversion failed
                  SashimiPlot.justConvertFilenames(species,fl.RootDir()+'/SashimiPlots')
                else:
                  print 'No BAM files present in the root directory... skipping SashimiPlot analysis...'
            except Exception:
              print traceback.format_exc()
        
  try:
      clearObjectsFromMemory(exon_db); clearObjectsFromMemory(constitutive_probeset_db)
      clearObjectsFromMemory(go_annotations); clearObjectsFromMemory(original_microRNA_z_score_data)
      clearObjectsFromMemory(last_exon_region_db)
      """
      print 'local vars'
      all = [var for var in locals() if (var[:2], var[-2:]) != ("__", "__")]
      for var in all:
            try:
                if len(locals()[var])>500: print var, len(locals()[var])
            except Exception: null=[]
      """
  except Exception: null=[]

  #print '!!!!!finished'
  #returnLargeGlobalVars()
  
  end_time = time.time(); time_diff = int(end_time-start_time)
  universalPrintFunction(["Analyses finished in %d seconds" % time_diff])
  #universalPrintFunction(["Hit Enter/Return to exit AltAnalyze"])
  for dataset in exp_file_location_db:
       fl = exp_file_location_db[dataset]; results_dir = filepath(fl.RootDir())
        
  ### Perform GO-Elite Analysis
  if pathway_permutations != 'NA':
      goelite_run = False
      print '\nBeginning to run GO-Elite analysis on alternative exon results'
      elite_input_dirs = ['AltExonConfirmed','AltExon','regulated','upregulated','downregulated'] ### Run GO-Elite multiple times to ensure heatmaps are useful and to better organize results
      for elite_dir in elite_input_dirs:
          file_dirs = results_dir+'GO-Elite/'+elite_dir,results_dir+'GO-Elite/denominator',results_dir+'GO-Elite/'+elite_dir
          input_dir = results_dir+'GO-Elite/'+elite_dir 
          try: input_files = read_directory(input_dir) ### Are there any files to analyze?
          except Exception: input_files = []
          if len(input_files)>0:
              variables = species,mod,pathway_permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,returnPathways,file_dirs,root
              try: GO_Elite.remoteAnalysis(variables,'non-UI',Multi=mlp); goelite_run = True
              except Exception,e:
                print e
                print "GO-Elite analysis failed"
              try: GO_Elite.moveMAPPFinderFiles(file_dirs[0])
              except Exception: print 'Input GO-Elite files could NOT be moved.'
              try: GO_Elite.moveMAPPFinderFiles(file_dirs[1])
              except Exception: print 'Input GO-Elite files could NOT be moved.'
      if goelite_run == False:
        print 'No GO-Elite input files to analyze (check your criterion).'
  print_out = 'Analysis complete. AltAnalyze results\nexported to "AltResults/AlternativeOutput".'

  try:
      if root !='' and root !=None:
          print "Analysis Complete\n";
          UI.InfoWindow(print_out,'Analysis Completed!')
          try: dataset_name = dataset_name
          except: dataset_name = dataset
          #tl = Toplevel(); SummaryResultsWindow(tl,'AS',results_dir,dataset_name,'specific',summary_data_db2)
          tl = Toplevel(); SummaryResultsWindow(tl,'GE',results_dir,dataset_name,'specific',summary_data_db)
  except Exception:
    print traceback.format_exc()
    pass #print 'Failed to open GUI.'
  skip_intro = 'yes'
  if root !='' and root !=None:
      if pathway_permutations == 'NA' and run_from_scratch != 'Annotate External Results':
          try: UI.getUpdatedParameters(array_type,species,run_from_scratch,file_dirs)
          except Exception: pass
      try: AltAnalyzeSetup('no')
      except Exception: sys.exit()
          
def exportSummaryResults(summary_results_db,analysis_method,aspire_output_list,aspire_output_gene_list,annotate_db,array_type,number_events_analyzed,root_dir):       
      try:
          ResultsExport_module.outputSummaryResults(summary_results_db,'',analysis_method,root_dir)
          #ResultsExport_module.outputSummaryResults(summary_results_db2,'-uniprot_attributes',analysis_method)
          ResultsExport_module.compareAltAnalyzeResults(aspire_output_list,annotate_db,number_events_analyzed,'no',analysis_method,array_type,root_dir)
          ResultsExport_module.compareAltAnalyzeResults(aspire_output_gene_list,annotate_db,'','yes',analysis_method,array_type,root_dir) 
      except UnboundLocalError: print "...No results to summarize" ###Occurs if there is a problem parsing these files
          
def checkGOEliteProbesets(fn,species):
    ### Get all probesets in GO-Elite files
    mod_source = 'Ensembl'+'-'+'Affymetrix'
    import gene_associations
    try: ensembl_to_probeset_id = gene_associations.getGeneToUid(species,mod_source)
    except Exception: ensembl_to_probeset_id={}
    mod_source = 'EntrezGene'+'-'+'Affymetrix'
    try: entrez_to_probeset_id = gene_associations.getGeneToUid(species,mod_source)
    except Exception: entrez_to_probeset_id={}
    probeset_db={}
    for gene in ensembl_to_probeset_id:
        for probeset in ensembl_to_probeset_id[gene]: probeset_db[probeset]=[]
    for gene in entrez_to_probeset_id:
        for probeset in entrez_to_probeset_id[gene]: probeset_db[probeset]=[]

    ###Import an Affymetrix array annotation file (from http://www.affymetrix.com) and parse out annotations
    csv_probesets = {}; x=0; y=0
    fn=filepath(fn); status = 'no'
    for line in open(fn,'r').readlines():             
        probeset_data = string.replace(line,'\n','')  #remove endline
        probeset_data = string.replace(probeset_data,'---','')
        affy_data = string.split(probeset_data[1:-1],'","')
        if x==0 and line[0]!='#':
            x=1; affy_headers = affy_data
            for header in affy_headers:
                y = 0
                while y < len(affy_headers):
                    if 'Probe Set ID' in affy_headers[y] or 'probeset_id' in affy_headers[y]: ps = y
                    y+=1
        elif x == 1:
            try: probeset = affy_data[ps]; csv_probesets[probeset]=[]
            except Exception: null=[]
    for probeset in csv_probesets:
        if probeset in probeset_db: status = 'yes';break
    return status

class SpeciesData:
    def __init__(self, abrev, species, systems, taxid):
        self._abrev = abrev; self._species = species; self._systems = systems; self._taxid = taxid
    def SpeciesCode(self): return self._abrev
    def SpeciesName(self): return self._species
    def Systems(self): return self._systems
    def TaxID(self): return self._taxid
    def __repr__(self): return self.SpeciesCode()+'|'+SpeciesName
    
def getSpeciesInfo():
    ### Used by AltAnalyze
    UI.importSpeciesInfo(); species_names={}
    for species_full in species_codes:
        sc = species_codes[species_full]; abrev = sc.SpeciesCode()
        species_names[abrev] = species_full
    return species_codes,species_names

def importGOEliteSpeciesInfo():
    filename = 'Config/goelite_species.txt'; x=0
    fn=filepath(filename); species_codes={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        abrev,species,taxid,compatible_mods = string.split(data,'\t')
        if x==0: x=1
        else:
            compatible_mods = string.split(compatible_mods,'|')
            sd = SpeciesData(abrev,species,compatible_mods,taxid)
            species_codes[species] = sd
    return species_codes

def exportGOEliteSpeciesInfo(species_codes):
    fn=filepath('Config/goelite_species.txt'); data = open(fn,'w'); x=0
    header = string.join(['species_code','species_name','tax_id','compatible_algorithms'],'\t')+'\n'
    data.write(header)
    for species in species_codes:
        if 'other' not in species and 'all-' not in species:
            sd = species_codes[species]
            mods = string.join(sd.Systems(),'|')
            values = [sd.SpeciesCode(),sd.SpeciesName(),sd.TaxID(),mods]
            values = string.join(values,'\t')+'\n'
            data.write(values)
    data.close()

def TimeStamp():
    time_stamp = time.localtime()
    year = str(time_stamp[0]); month = str(time_stamp[1]); day = str(time_stamp[2])
    if len(month)<2: month = '0'+month
    if len(day)<2: day = '0'+day
    return year+month+day

def verifyFile(filename):
    status = 'not found'
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = 'found';break
    except Exception: status = 'not found'
    return status
        
def verifyFileLength(filename):
    count = 0
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            count+=1
            if count>9: break
    except Exception: null=[]
    return count

def verifyGroupFileFormat(filename):
    correct_format = False
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            if len(string.split(data,'\t'))==3:
                correct_format = True
            break
    except Exception: correct_format = False
    return correct_format

def parseExcludeGuides(excludeGuides):
    guides=[]
    for line in open(excludeGuides,'rU').readlines():
        data = cleanUpLine(line)
        if 'Guide' not in data:
            guides.append(data)
    return guides
        
def displayHelp():
    fn=filepath('Documentation/commandline.txt')
    print '\n################################################\nAltAnalyze Command-Line Help'
    for line in open(fn,'rU').readlines():
        print cleanUpLine(line)
    print '\n################################################ - END HELP'
    sys.exit()
    
def searchDirectory(directory,var):
    directory = unique.filepath(directory)

    files = unique.read_directory(directory)
    version = unique.getCurrentGeneDatabaseVersion()
    for file in files:
        if var in file:
            location = string.split(directory+'/'+file,version)[1][1:]
            return [location]
            break
    
###### Command Line Functions (AKA Headless Mode) ######
def commandLineRun():
    print 'Running commandline options'
    import getopt
    #/hd3/home/nsalomonis/normalization/mir1 - boxer
    #python AltAnalyze.py --species Mm --arraytype "3'array" --celdir "C:/CEL" --output "C:/CEL" --expname miR1_column --runGOElite yes --GEelitepval 0.01
    #python AltAnalyze.py --species Hs --arraytype "3'array" --FEdir "C:/FEfiles" --output "C:/FEfiles" --channel_to_extract "green/red ratio" --expname cancer --runGOElite yes --GEelitepval 0.01
    #python AltAnalyze.py --celdir "C:/CEL" --output "C:/CEL" --expname miR1_column
    #open ./AltAnalyze.app --celdir "/Users/nsalomonis/Desktop" --output "/Users/nsalomonis/Desktop" --expname test
    #python AltAnalyze.py --species Mm --arraytype "3'array" --expdir "C:/CEL/ExpressionInput/exp.miR1_column.txt" --output "C:/CEL" --runGOElite yes --GEelitepval 0.01
    #python AltAnalyze.py --species Mm --platform RNASeq --bedDir "/Users/nsalomonis/Desktop/code/AltAnalyze/datasets/BedFiles" --groupdir "/Users/nsalomonis/Desktop/code/AltAnalyze/datasets/BedFiles/ExpressionInput/groups.test.txt" --compdir "/Users/nsalomonis/Desktop/code/AltAnalyze/datasets/BedFiles/ExpressionInput/comps.test.txt" --output "/Users/nsalomonis/Desktop/code/AltAnalyze/datasets/BedFiles" --expname "test"
    #python AltAnalyze.py --species Mm --platform RNASeq --filterdir "/Users/nsalomonis/Desktop/code/AltAnalyze/datasets/BedFiles/" --output "/Users/nsalomonis/Desktop/code/AltAnalyze/datasets/BedFiles"
    #python AltAnalyze.py --expdir "/Users/nsalomonis/Desktop/Nathan/ExpressionInput/exp.test.txt" --exonMapFile "/Users/nsalomonis/Desktop/Nathan/hgu133_probe.txt" --species Hs --platform "3'array" --output "/Users/nsalomonis/Desktop/Nathan"
    #python AltAnalyze.py --species Hs --platform "3'array" --expname test --channelToExtract green --FEdir /Users/saljh8/Downloads/AgllentTest/ --output /Users/saljh8/Downloads/AgllentTest/
    global apt_location; global root_dir; global probability_statistic; global log_file; global summary_data_db; summary_data_db={}

    ###required
    marker_finder='no'
    manufacturer='Affymetrix'
    constitutive_source='Ensembl'
    ensembl_version = 'current'
    species_code = None
    species = None
    main_input_folder = None
    output_dir = None
    array_type = None
    input_annotation_file = None
    groups_file = None
    comps_file = None
    input_cdf_file = None
    exp_name = None
    run_GOElite = 'yes'
    visualize_qc_results = 'yes'
    run_lineage_profiler = 'yes'
    input_exp_file = ''
    cel_file_dir = ''
    input_stats_file = ''
    input_filtered_dir = ''
    external_annotation_dir = ''
    xhyb_remove = 'no'
    update_method = []
    update_dbs = 'no'
    analyze_all_conditions = 'no'
    return_all = 'no'
    additional_array_types = []
    remove_intronic_junctions = 'no'
    ignore_built_species = 'no'
    build_exon_bedfile = 'no'
    compendiumType = 'protein_coding'
    probability_statistic = 'unpaired t-test'
    specific_array_type = None
    additional_resources = [None]
    wpid = None
    mod = 'Ensembl'
    transpose = False
    input_file_dir = None
    denom_file_dir = None
    image_export = []
    selected_species = ['Hs','Mm','Rn'] ### These are the species that additional array types are currently supported
    selected_platforms = ['AltMouse','exon','gene','junction']
    returnPathways = 'no'
    compendiumPlatform = 'gene'
    exonMapFile = None
    platformType = None ### This option is used to store the orignal platform type
    perform_alt_analysis = 'no'
    mappedExonAnalysis = False ### Map the original IDs to the RNA-Seq exon database  (when True)
    microRNA_prediction_method = None
    pipelineAnalysis = True
    OntologyID=''
    PathwaySelection=''
    GeneSetSelection=''
    interactionDirs=[]
    inputType='ID list'
    Genes=''
    degrees='direct'
    includeExpIDs=True
    update_interactions=False
    data_type = 'raw expression'
    batch_effects = 'no'
    channel_to_extract = None
    normalization = False
    justShowTheseIDs = ''
    display=False
    accessoryAnalysis=''
    modelSize=None
    geneModel=False
    run_from_scratch = None
    systemToUse = None ### For other IDs
    custom_reference = False
    multiThreading = True
    genesToReport = 60
    correlateAll = True
    expression_data_format='log'
    runICGS=False
    IDtype=None
    runKallisto = False
    input_fastq_dir = ''
    ChromiumSparseMatrix=''
    perform_tests=False
    testType = "fast"
    inputTestData = "text"
    customFASTA = None
    filterFile = None
    PearsonThreshold = 0.1
    returnCentroids = False
    runCompleteWorkflow=True
    k=None
    original_arguments = sys.argv
    arguments=[]
    for arg in original_arguments:
        arg = string.replace(arg,'\xe2\x80\x9c','') ### These are non-standard forward quotes
        arg = string.replace(arg,'\xe2\x80\x9d','') ### These are non-standard reverse quotes
        arg = string.replace(arg,'\xe2\x80\x93','-') ### These are non-standard dashes
        arg = string.replace(arg,'\x96','-') ### These are non-standard dashes
        arg = string.replace(arg,'\x93','') ### These are non-standard forward quotes
        arg = string.replace(arg,'\x94','') ### These are non-standard reverse quotes
        arguments.append(arg)
    
    print '\nArguments input:',arguments,'\n'
    if '--help' in arguments[1:] or '--h' in arguments[1:]:
        try: displayHelp() ### Print out a help file and quit
        except Exception: print 'See: http://www.altanalyze.org for documentation and command-line help';sys.exit()
    
    if 'AltAnalyze' in arguments[1]:
        arguments = arguments[1:] ### Occurs on Ubuntu with the location of AltAnalyze being added to sys.argv (exclude this since no argument provided for this var)
    try:
        options, remainder = getopt.getopt(arguments[1:],'', ['species=', 'mod=','elitepval=', 'elitepermut=',
                                                         'method=','zscore=','pval=','num=',
                                                         'runGOElite=','denom=','output=','arraytype=',
                                                         'celdir=','expdir=','output=','statdir=',
                                                         'filterdir=','cdfdir=','csvdir=','expname=',
                                                         'dabgp=','rawexp=','avgallss=','logexp=',
                                                         'inclraw=','runalt=','altmethod=','altp=',
                                                         'probetype=','altscore=','GEcutoff=',
                                                         'exportnormexp=','calcNIp=','runMiDAS=',
                                                         'GEcutoff=','GEelitepval=','mirmethod=','ASfilter=',
                                                         'vendor=','GEelitefold=','update=','version=',
                                                         'analyzeAllGroups=','GEeliteptype=','force=',
                                                         'resources_to_analyze=', 'dataToAnalyze=','returnAll=',
                                                         'groupdir=','compdir=','annotatedir=','additionalScore=',
                                                         'additionalAlgorithm=','noxhyb=','platform=','bedDir=',
                                                         'altpermutep=','altpermute=','removeIntronOnlyJunctions=',
                                                         'normCounts=','buildExonExportFile=','groupStat=',
                                                         'compendiumPlatform=','rpkm=','exonExp=','specificArray=',
                                                         'ignoreBuiltSpecies=','ORAstat=','outputQCPlots=',
                                                         'runLineageProfiler=','input=','image=', 'wpid=',
                                                         'additional=','row_method=','column_method=',
                                                         'row_metric=','column_metric=','color_gradient=',
                                                         'transpose=','returnPathways=','compendiumType=',
                                                         'exonMapFile=','geneExp=','labels=','contrast=',
                                                         'plotType=','geneRPKM=','exonRPKM=','runMarkerFinder=',
                                                         'update_interactions=','includeExpIDs=','degrees=',
                                                         'genes=','inputType=','interactionDirs=','GeneSetSelection=',
                                                         'PathwaySelection=','OntologyID=','dataType=','combat=',
                                                         'channelToExtract=','showIntrons=','display=','join=',
                                                         'uniqueOnly=','accessoryAnalysis=','inputIDType=','outputIDType=',
                                                         'FEdir=','channelToExtract=','AltResultsDir=','geneFileDir=',
                                                         'AltResultsDir=','modelSize=','geneModel=','reference=',
                                                         'multiThreading=','multiProcessing=','genesToReport=',
                                                         'correlateAll=','normalization=','justShowTheseIDs=',
                                                         'direction=','analysisType=','algorithm=','rho=',
                                                         'clusterGOElite=','geneSetName=','runICGS=','IDtype=',
                                                         'CountsCutoff=','FoldDiff=','SamplesDiffering=','removeOutliers=',
                                                         'featurestoEvaluate=','restrictBy=','ExpressionCutoff=',
                                                         'excludeCellCycle=','runKallisto=','fastq_dir=','FDR=',
                                                         'reimportModelScores=','separateGenePlots=','ChromiumSparseMatrix=',
                                                         'test=','testType=','inputTestData=','customFASTA=','i=',
                                                         'excludeGuides=','cellHarmony=','BAM_dir=','filterFile=',
                                                         'correlationCutoff=','referenceType=','DE=','cellHarmonyMerge=',
                                                         'o=','dynamicCorrelation=','runCompleteWorkflow=','adjp=',
                                                         'fold=','performDiffExp=','centerMethod=', 'k='])
    except Exception:
        print traceback.format_exc()
        print "There is an error in the supplied command-line arguments (each flag requires an argument)"; sys.exit()

    for opt, arg in options:
        #print [opt, arg]
        if opt == '--species': species=arg
        elif opt == '--arraytype':
            if array_type != None: additional_array_types.append(arg)
            else: array_type=arg; platform = array_type
            if specific_array_type == None: specific_array_type = platform
        elif opt == '--test':
            try: perform_tests.append(arg)
            except Exception: perform_tests = [arg]
        elif opt == '--testType':
            testType = arg
        elif opt == '--inputTestData':
            inputTestData = arg
        elif opt == '--exonMapFile':
            perform_alt_analysis = 'yes' ### Perform alternative exon analysis
            exonMapFile = arg
        elif opt == '--specificArray': specific_array_type = arg  ### e.g., hGlue
        elif opt == '--celdir':
            arg = verifyPath(arg)
            cel_file_dir=arg
        elif opt == '--bedDir' or opt == '--BAM_dir':
            arg = verifyPath(arg)
            cel_file_dir=arg
        elif opt == '--ChromiumSparseMatrix':
            arg = verifyPath(arg)
            ChromiumSparseMatrix=arg
        elif opt == '--FEdir':
            arg = verifyPath(arg)
            cel_file_dir = arg
        elif opt == '--expdir':
            arg = verifyPath(arg)
            input_exp_file=arg
        elif opt == '--statdir':
            arg = verifyPath(arg)
            input_stats_file=arg
        elif opt == '--filterdir':
            arg = verifyPath(arg)
            input_filtered_dir=arg
        elif opt == '--groupdir':
            arg = verifyPath(arg)
            groups_file=arg
        elif opt == '--compdir':
            arg = verifyPath(arg)
            comps_file=arg
        elif opt == '--cdfdir':
            arg = verifyPath(arg)
            input_cdf_file=arg
        elif opt == '--csvdir':
            arg = verifyPath(arg)
            input_annotation_file=arg
        elif opt == '--expname': exp_name=arg
        elif opt == '--output' or opt == '--o':
            arg = verifyPath(arg)
            output_dir=arg
        elif opt == '--vendor': manufacturer=arg
        elif opt == '--runICGS': runICGS=True
        elif opt == '--IDtype': IDtype=arg
        elif opt == '--ignoreBuiltSpecies': ignore_built_species=arg
        elif opt == '--platform':
            if array_type != None: additional_array_types.append(arg)
            else: array_type=arg; platform = array_type
            if specific_array_type == None: specific_array_type = platform
        elif opt == '--update': update_dbs='yes'; update_method.append(arg)
        elif opt == '--version': ensembl_version = arg
        elif opt == '--compendiumPlatform': compendiumPlatform=arg ### platform for which the LineageProfiler compendium is built on
        elif opt == '--force': force=arg
        elif opt == '--input' or opt == '--i':
            arg = verifyPath(arg)
            input_file_dir=arg
            #input_exp_file=arg
            pipelineAnalysis = False ### If this option is entered, only perform the indicated analysis
        elif opt == '--image': image_export.append(arg)
        elif opt == '--wpid': wpid=arg
        elif opt == '--mod': mod=arg
        elif opt == '--runKallisto':
            if arg == 'yes' or string.lower(arg) == 'true':
                runKallisto = True
        elif opt == '--fastq_dir':
            input_fastq_dir = arg
        elif opt == '--customFASTA':
            customFASTA = arg
        elif opt == '--additional':
            if additional_resources[0] == None:
                additional_resources=[]
                additional_resources.append(arg)
            else:
                additional_resources.append(arg)
        elif opt == '--transpose':
            if arg == 'True': transpose = True
        elif opt == '--runLineageProfiler' or opt == '--cellHarmony': ###Variable declared here and later (independent analysis here or pipelined with other analyses later)
            run_lineage_profiler=arg
        elif opt == '--compendiumType': ### protein-coding, ncRNA, or exon
            compendiumType=arg
        elif opt == '--denom':
            denom_file_dir=arg ### Indicates that GO-Elite is run independent from AltAnalyze itself
        elif opt == '--accessoryAnalysis':
            accessoryAnalysis = arg
        elif opt == '--channelToExtract': channel_to_extract=arg
        elif opt == '--genesToReport': genesToReport = int(arg)
        elif opt == '--correlateAll': correlateAll = True
        elif opt == '--direction': direction = arg
        elif opt == '--logexp': expression_data_format=arg
        elif opt == '--geneRPKM': rpkm_threshold=arg
        elif opt == '--correlationCutoff': PearsonThreshold=float(arg)
        elif opt == '--DE':
            if string.lower(arg) == 'true':
                DE = True
            else:
                DE = False
        elif opt == '--referenceType':
            if string.lower(arg) == 'centroid': returnCentroids = True; CenterMethod='centroid'
            else: CenterMethod='median'; returnCentroids = True
        elif opt == '--multiThreading' or opt == '--multiProcessing':
            multiThreading=arg
            if multiThreading == 'yes': multiThreading = True
            elif 'rue' in multiThreading: multiThreading = True
            else: multiThreading = False

    if perform_tests != False:
        ### Requires the mouse RNASeq database
        ### python AltAnalyze.py --test --testType ICGS --inputTestData text
        ### python AltAnalyze.py --test --testType ICGS --inputTestData BAM
        ### python AltAnalyze.py --test --testType ICGS --inputTestData FASTQ
        ### python AltAnalyze.py --test --testType ICGS --inputTestData 10X
        count = verifyFileLength('AltDatabase/demo_data/ReadMe.txt')
        if count==0:
            file_location_defaults = UI.importDefaultFileLocations()
            goelite_url = file_location_defaults['goelite'].Location()
            fln,status = update.download(goelite_url+'TestData/demo_data.zip','AltDatabase/NoVersion','')
            if 'Internet' not in status: print "Demo data downloaded."
        
        if 'ICGS' in perform_tests:
            from tests.scripts import ICGS_test
            if runKallisto:
                inputTestData = "FASTQ"
            ICGS_test.runICGStest(testType=testType,inputData=inputTestData)
        sys.exit()
        
    if 'other' in manufacturer or 'Other' in manufacturer:
        ### For other IDs
        systemToUse = array_type
        if array_type == None:
            print 'Please indicate a ID type as --platform when setting vendor equal to "Other IDs"'; sys.exit()
        array_type = "3'array"
    if array_type == 'RNASeq': manufacturer = array_type
    if platformType == None: platformType = array_type
    if perform_alt_analysis == 'yes':
        if platform == "3'array":
            mappedExonAnalysis = True
            cel_file_dir = input_exp_file
            exp_name = export.findFilename(input_exp_file)
            exp_name = string.replace(exp_name,'.txt','')
            exp_name = string.replace(exp_name,'exp.','')
            input_exp_file = ''
            ### To perform alternative exon analyses for platforms without a dedicated database, must happing appropriate mapping info or array type data
            ### (will need to perform downstream testing for unsupported Affymetrix exon, gene and junction arrays)
            if exonMapFile == None and specific_array_type == None and cel_file_dir == '':
                print_out = "\nUnable to run!!! Please designate either a specific platfrom (e.g., --specificArray hgU133_2), select CEL files, or an "
                print_out += "exon-level mapping file location (--exonMapFile C:/mapping.txt) to perform alternative exon analyses for this platform."

        ### Will need to check here to see if the platform is supported (local or online files) OR wait until an error is encountered later
        
    """ Check to see if a database is already installed """
    try: current_species_dirs = unique.read_directory('/AltDatabase')
    except Exception: current_species_dirs=[]
    if len(current_species_dirs)==0 and update_dbs != 'yes':
        print "Please install a database before running AltAnalyze. Please note, AltAnalyze may need to install additional files later for RNASeq and LineageProfiler for some species, automatically. Make sure to list your platform as RNASeq if analyzing RNA-Seq data (--platform RNASeq)."
        print "Example:\n"
        print 'python AltAnalyze.py --species Hs --update Official --version EnsMart72';sys.exit()

    ######## Perform analyses independent from AltAnalyze database centric analyses that require additional parameters
    if len(image_export) > 0 or len(accessoryAnalysis)>0 or runICGS:
        if runICGS:
            #python AltAnalyze.py --runICGS yes --platform "RNASeq" --species Hs --column_method hopach --column_metric euclidean --rho 0.3 --ExpressionCutoff 1 --FoldDiff 4 --SamplesDiffering 3 --restrictBy protein_coding --excludeCellCycle conservative --removeOutliers yes --expdir /RNA-Seq/run1891_normalized.txt
            #python AltAnalyze.py --runICGS yes --expdir "/Users/saljh8/Desktop/demo/Myoblast/ExpressionInput/exp.myoblast.txt" --platform "3'array" --species Hs --GeneSetSelection BioMarkers --PathwaySelection Heart --column_method hopach --rho 0.4 --ExpressionCutoff 200 --justShowTheseIDs "NKX2-5 T TBX5" --FoldDiff 10 --SamplesDiffering 3 --excludeCellCycle conservative
            try: species = species
            except Exception: 'Please designate a species before continuing (e.g., --species Hs)'
            try: array_type = array_type
            except Exception: 'Please designate a species before continuing (e.g., --species Hs)'
            if len(cel_file_dir)>0: ### For BED files or BAM files
                if len(cel_file_dir) > 0: pass
                else: 'Please indicate a source folder (e.g., --bedDir /data/BAMFiles)'
            else:
                if len(input_exp_file) > 0: pass
                else: 'Please indicate a source folder or expression file (e.g., --expdir /dataset/singleCells.txt)'
            if array_type == 'Other' or 'Other' in array_type:
                if ':' in array_type:
                    array_type, IDtype = string.split(array_type)
                array_type == "3'array"
                if IDtype == None: IDtype = manufacturer
        
            row_method = 'hopach'
            column_method = 'hopach'
            row_metric = 'cosine'
            column_metric = 'cosine'
            color_gradient = 'yellow_black_blue'
            contrast=3
            vendor = manufacturer
            GeneSelection = ''
            PathwaySelection = ''
            GeneSetSelection = 'None Selected'
            excludeCellCycle = True
            rho_cutoff = 0.4
            restrictBy = None
            featurestoEvaluate = 'Genes'
            ExpressionCutoff = 1
            CountsCutoff = 0.9
            FoldDiff = 4
            SamplesDiffering = 4
            JustShowTheseIDs=''
            removeOutliers = False
            excludeGuides = None
            PathwaySelection=[]
            dynamicCorrelation=True
            runCompleteWorkflow=False
            if ChromiumSparseMatrix != '':
                rho_cutoff = 0.3
                column_metric = 'euclidean'
                restrictBy = 'protein_coding'
            for opt, arg in options: ### Accept user input for these hierarchical clustering variables
                if opt == '--row_method':
                    row_method=arg
                    if row_method == 'None': row_method = None
                elif opt == '--column_method':
                    column_method=arg
                    if column_method == 'None': column_method = None
                elif opt == '--row_metric': row_metric=arg
                elif opt == '--column_metric': column_metric=arg
                elif opt == '--color_gradient': color_gradient=arg
                elif opt == '--GeneSetSelection': GeneSetSelection=arg
                elif opt == '--PathwaySelection': PathwaySelection.append(arg)
                elif opt == '--genes': GeneSelection=arg
                elif opt == '--ExpressionCutoff': ExpressionCutoff=arg
                elif opt == '--normalization': normalization=arg
                elif opt == '--justShowTheseIDs': justShowTheseIDs=arg
                elif opt == '--rho': rho_cutoff=float(arg)
                elif opt == '--clusterGOElite':clusterGOElite=float(arg)
                elif opt == '--CountsCutoff':CountsCutoff=int(float(arg))
                elif opt == '--FoldDiff':FoldDiff=float(arg)
                elif opt == '--SamplesDiffering':SamplesDiffering=int(float(arg))
                elif opt == '--excludeGuides': excludeGuides=arg
                elif opt == '--dynamicCorrelation': dynamicCorrelation=arg
                elif opt == '--k': k=int(arg)
                elif opt == '--runCompleteWorkflow':
                    runCompleteWorkflow=arg
                    if string.lower(arg)=='false' or string.lower(arg)=='no':
                        runCompleteWorkflow = False
                    else:
                        runCompleteWorkflow = True
                elif opt == '--removeOutliers':
                    removeOutliers=arg
                    if removeOutliers=='yes' or removeOutliers=='True':
                        removeOutliers = True
                elif opt == '--featurestoEvaluate':featurestoEvaluate=arg
                elif opt == '--restrictBy':
                    if arg == 'None': restrictBy = None
                    else: restrictBy=arg
                elif opt == '--excludeCellCycle':
                    excludeCellCycle=arg
                    if excludeCellCycle == 'False' or excludeCellCycle == 'no': excludeCellCycle = False
                    elif excludeCellCycle == 'True' or excludeCellCycle == 'yes' or excludeCellCycle == 'conservative': excludeCellCycle = True
                elif opt == '--contrast':
                    try: contrast=float(arg)
                    except Exception: print '--contrast not a valid float';sys.exit()
                elif opt == '--vendor': vendor=arg
                elif opt == '--display':
                    if arg=='yes':
                        display=True
                    elif arg=='True':
                        display=True
                    else:
                        display=False
            if excludeGuides!=None:
                if '.txt' in excludeGuides:
                    try: excludeGuides = parseExcludeGuides(excludeGuides)
                    except Exception:
                        print 'Failure to parse input excludeGuides text file. Check to see if the correct file location is provided.'
                        sys.exit()
            if len(PathwaySelection)==0: PathwaySelection=''
            if len(GeneSetSelection)>0 or GeneSelection != '':
                gsp = UI.GeneSelectionParameters(species,array_type,vendor)
                gsp.setGeneSet(GeneSetSelection)
                gsp.setPathwaySelect(PathwaySelection)
                gsp.setGeneSelection(GeneSelection)
                gsp.setJustShowTheseIDs(JustShowTheseIDs)
                gsp.setNormalize('median')
                gsp.setExcludeGuides(excludeGuides)
                gsp.setK(k)
                gsp.setSampleDiscoveryParameters(ExpressionCutoff,CountsCutoff,FoldDiff,SamplesDiffering, dynamicCorrelation,  
                    removeOutliers,featurestoEvaluate,restrictBy,excludeCellCycle,column_metric,column_method,rho_cutoff) 

            import RNASeq
            mlp_instance = mlp
            
            if exp_name == None:
                exp_name = export.findFilename(input_exp_file)
                exp_name = string.replace(exp_name,'.txt','')
                exp_name = string.replace(exp_name,'exp.','')
            
            if cel_file_dir != '':
                expFile = output_dir + '/ExpressionInput/'+ 'exp.'+exp_name+'.txt' ### cel_file_dir will point to the input directory
            elif ChromiumSparseMatrix != '':
                expFile = output_dir + '/ExpressionInput/'+ 'exp.'+exp_name+'.txt'
            elif input_exp_file !='':
                if 'ExpressionInput' in input_exp_file: expFile = input_exp_file
                else:
                    ### Copy over expression file to ExpressionInput
                    expdir2 = string.replace(input_exp_file,'exp.','')
                    if output_dir == None:
                        ### Not suppplied, define relative to the input expression file
                        root_dir = export.findParentDir(input_exp_file)
                    else:
                        root_dir = output_dir
                    expFile = root_dir+'/ExpressionInput/exp.'+export.findFilename(expdir2)
                    export.copyFile(input_exp_file, expFile)

            ### Groups file if present to the output directory ExpressionInput folder
            if 'exp.' in input_exp_file:
                initial_groups = string.replace(input_exp_file,'exp.','groups.')
            else:
                initial_groups = export.findParentDir(input_exp_file)+'/groups.'+export.findFilename(input_exp_file)
            try:
                groups_file = string.replace(expFile,'exp.','groups.') ### destination file
                export.copyFile(initial_groups, groups_file)
                print 'Copied the groups file to ExpressionInput folder in the output directory'
            except Exception:
                print 'No groups file present in the input file folder.'

            global log_file
            try:
                if len(output_dir)>0:
                    root_dir = output_dir
                else:
                    forceError
            except Exception:
                try: root_dir = export.findParentDir(expFile)
                except Exception:
                    print 'Please include an output directory for the AltAnalyze results (e.g., --output /Data/Results)';sys.exit()
            root_dir = string.replace(root_dir,'/ExpressionInput','')
                    
            fl = UI.ExpressionFileLocationData('','','',''); fl.setFeatureNormalization('none')

            try: fl.setExpFile(expFile)
            except Exception:
                expFile = root_dir+'/ExpressionInput/exp.'+exp_name+'.txt'
                fl.setExpFile(expFile)

            fl.setArrayType(array_type)
            fl.setOutputDir(root_dir)
            fl.setMultiThreading(multiThreading)
            exp_file_location_db={}; exp_file_location_db[exp_name]=fl

            ### Assign variables needed to run Kallisto from FASTQ files
            if runKallisto and len(input_fastq_dir)==0:
                #python AltAnalyze.py --runICGS yes --platform "RNASeq" --species Mm --column_method hopach --rho 0.4 --ExpressionCutoff 1 --FoldDiff 4 --SamplesDiffering 1 --excludeCellCycle strict --output  /Users/saljh8/Desktop/Grimes/GEC14074 --expname test --fastq_dir /Users/saljh8/Desktop/Grimes/GEC14074 
                print 'Please include the flag "--fastq_dir" in the command-line arguments with an appropriate path';sys.exit()
            elif len(input_fastq_dir)>0:
                fl.setRunKallisto(input_fastq_dir)
                fl.setArrayType("3'array")
                fl.setMultiThreading(multiThreading)
                array_type = "3'array"
                if customFASTA!=None:
                    fl.setCustomFASTA(customFASTA)

            ### Assign variables needed to run BAMtoBED and/or BED file count analysis
            if len(cel_file_dir)>0 and array_type=='RNASeq':
                #python AltAnalyze.py --runICGS yes --platform "RNASeq" --species Mm --column_method hopach --rho 0.4 --ExpressionCutoff 1 --FoldDiff 4 --SamplesDiffering 1 --excludeCellCycle strict --output  /Users/saljh8/Desktop/Grimes/GEC14074 --expname test --bedDir /Users/saljh8/Desktop/Grimes/GEC14074 --multiProcessing no
                fl.setCELFileDir(cel_file_dir)
                fl.setMultiThreading(multiThreading)
                fl.setExonBedBuildStatus('no')
                fl.setFeatureNormalization('RPKM')
                fl.setArrayType(array_type)
                fl.setRootDir(root_dir)
                ### Import expression defaults
                expr_defaults, alt_exon_defaults, functional_analysis_defaults, goelite_defaults = UI.importDefaults(array_type,species)
                dabg_p, rpkm_threshold, gene_exp_threshold, exon_exp_threshold, exon_rpkm_threshold, expression_threshold, perform_alt_analysis, analyze_as_groups, expression_data_format, normalize_feature_exp, normalize_gene_data, avg_all_for_ss, include_raw_data, probability_statistic, FDR_statistic, batch_effects, marker_finder, visualize_qc_results, run_lineage_profiler, null = expr_defaults
                fl.setRPKMThreshold(rpkm_threshold)
                fl.setExonExpThreshold(exon_exp_threshold)
                fl.setGeneExpThreshold(gene_exp_threshold)
                fl.setExonRPKMThreshold(exon_rpkm_threshold)
                fl.setJunctionExpThreshold(expression_threshold)

            elif len(ChromiumSparseMatrix)>0:
                #python AltAnalyze.py --runICGS yes --platform "RNASeq" --species Mm --column_method hopach --rho 0.4 --ExpressionCutoff 1 --FoldDiff 4 --SamplesDiffering 1 --excludeCellCycle strict --output  /Users/saljh8/Desktop/Grimes/GEC14074 --expname test --ChromiumSparseMatrix /Users/saljh8/Desktop/Grimes/GEC14074 --multiProcessing no
                fl.setChromiumSparseMatrix(ChromiumSparseMatrix)
                fl.setMultiThreading(multiThreading)
                fl.setArrayType("3'array")
                array_type = "3'array"
                fl.setRootDir(root_dir)
            elif len(input_exp_file)>0:

                ### Dealing with an expression file which should not be treated as RNASeq workflow
                count = verifyFileLength(expFile[:-4]+'-steady-state.txt')
                if count<2:
                    if array_type != 'PSI' and array_type != 'exons':
                        array_type = "3'array" ### No steady-state file, must be an standard gene-level analysis
 
            time_stamp = timestamp()    
            log_file = filepath(root_dir+'/AltAnalyze_report-'+time_stamp+'.log')
            log_report = open(log_file,'w'); log_report.close()
            sys.stdout = Logger('')
            
            print "\nFull commandline:"
            try: print string.join(arguments,' ')
            except Exception: pass
            print ''

            count = verifyFileLength(expFile[:-4]+'-steady-state.txt')
            exonExpFile = str(expFile)
            if count>1:
                expFile = expFile[:-4]+'-steady-state.txt'
            elif array_type=='RNASeq' or len(ChromiumSparseMatrix)>0 or len(input_fastq_dir)>0:
                
                try:
                    ### Indicates that the steady-state file doesn't exist. The exp. may exist, be could be junction only so need to re-build from bed files here
                    values = species,exp_file_location_db,exp_name,mlp_instance
                    ### proceed to run the full discovery analysis here (Kallisto, BAM, BED, Chromium matrix)
                    UI.StatusWindow(values,'preProcessRNASeq') ### proceed to run the full discovery analysis here!!!
                    SteadyStateFile = expFile[:-4]+'-steady-state.txt'
                    status = verifyFile(SteadyStateFile)
                    if status == "found":
                        #fl.setExpFile(SteadyStateFile) ### This should not be over-written by a steady-state file
                        expFile = SteadyStateFile
                except Exception:
                    ### RNASeq is an official datatype that requires a steady-state file. However, for scRNA-Seq, usually the input is a text file or FASTQ which gets
                    ### changed to "3'array". We correct for this by excepting this error without doing anything else
                    #print traceback.format_exc();sys.exit()
                    pass

            if excludeCellCycle != False:
                print "Excluding Cell Cycle effects status:",excludeCellCycle

            ### Run ICGS through the GUI
            graphic_links = UI.RemotePredictSampleExpGroups(expFile, mlp_instance, gsp,(species,array_type)) ### proceed to run the full discovery analysis here!!!
                        
            ### Export Guide3 Groups automatically
            Guide3_results = graphic_links[-1][-1][:-4]+'.txt'
            new_groups_dir = RNASeq.exportGroupsFromClusters(Guide3_results,fl.ExpFile(),array_type,suffix='ICGS')
            
            exonExpFile,newExpFile,new_groups_dir  = UI.exportAdditionalICGSOutputs(expFile,Guide3_results,outputTSNE=True)
        
            fl.setExpFile(newExpFile) ### set this to the outlier removed version
            comps_file = string.replace(new_groups_dir,'groups.','comps.')
            fl.setGroupsFile(new_groups_dir)
            fl.setCompsFile(comps_file)
            exp_file_location_db[exp_name+'-ICGS'] = fl
            
            ### force MarkerFinder to be run
            input_exp_file = newExpFile  ### Point MarkerFinder to the new ICGS ordered copied expression file
            runMarkerFinder=False ### Not necessary for ICGS2 as MarkerFinder will already have been run
            if runMarkerFinder:
                update_method = ['markers']
            if runCompleteWorkflow == False:
                print 'ICGS run complete... halted prior to full differential comparison analysis'
                sys.exit()
            
        if 'WikiPathways' in image_export:
            #python AltAnalyze.py --input /Users/test/input/criterion1.txt --image WikiPathways --mod Ensembl --species Hs --wpid WP536
            if wpid==None:
                print 'Please provide a valid WikiPathways ID (e.g., WP1234)';sys.exit()
            if species==None:
                print 'Please provide a valid species ID for an installed database (to install: --update Official --species Hs --version EnsMart72Plus)';sys.exit()
            if input_file_dir==None:
                print 'Please provide a valid file location for your input IDs (also needs to inlcude system code and value column)';sys.exit()
            from visualization_scripts import WikiPathways_webservice
            try:
                print 'Attempting to output a WikiPathways colored image from user data'
                print 'mod:',mod
                print 'species_code:',species
                print 'wpid:',wpid
                print 'input GO-Elite ID file:',input_file_dir
                graphic_link = WikiPathways_webservice.visualizePathwayAssociations(input_file_dir,species,mod,wpid)
            except Exception,e:
                if 'force_no_matching_error' in traceback.format_exc():
                    print '\nUnable to run!!! None of the input IDs mapped to this pathway\n' 
                elif 'IndexError' in traceback.format_exc():
                    print '\nUnable to run!!! Input ID file does not have at least 3 columns, with the second column being system code\n'
                elif 'ValueError' in traceback.format_exc():
                    print '\nUnable to run!!! Input ID file error. Please check that you do not have extra rows with no data\n' 
                elif 'source_data' in traceback.format_exc():
                    print '\nUnable to run!!! Input ID file does not contain a valid system code\n'
                elif 'goelite' in traceback.format_exc():
                    print '\nUnable to run!!! A valid species database needs to first be installed. For example, run:'
                    print 'python AltAnalyze.py --update Official --species Hs --version EnsMart72\n'
                else:
                    print traceback.format_exc()
                    print '\nError generating the pathway "%s"' % wpid,'\n'
            try:
                printout = 'Finished exporting visualized pathway to:',graphic_link['WP']
                print printout,'\n'
            except Exception: None
            sys.exit()
        if 'FilterFile' in accessoryAnalysis:
            for opt, arg in options: ### Accept user input for these hierarchical clustering variables
                if opt == '--filterFile':
                    filterFile = arg
                if opt == '--input':
                    input_file = verifyPath(arg)
            output_file = input_file[:-4]+'-filtered.txt'
            from import_scripts import sampleIndexSelection
            filter_order = sampleIndexSelection.getFilters(filterFile)
            
            sampleIndexSelection.filterFile(input_file,output_file,filter_order)
            sys.exit()
        if 'MergeFiles' in accessoryAnalysis:
            #python AltAnalyze.py --accessoryAnalysis MergeFiles --input "C:\file1.txt" --input "C:\file2.txt" --output "C:\tables"
            files_to_merge=[]
            join_option='Intersection'
            uniqueOnly=False
            for opt, arg in options: ### Accept user input for these hierarchical clustering variables
                if opt == '--input':
                    arg = verifyPath(arg)
                    files_to_merge.append(arg)
                if opt == '--join': join_option = arg
                if opt == '--uniqueOnly': unique_only = arg
            if len(files_to_merge)<2:
                print 'Please designate two or more files to merge (--input)';sys.exit()
            UI.MergeFiles(files_to_merge, join_option, uniqueOnly, output_dir, None)
            sys.exit()

        if 'IDTranslation' in accessoryAnalysis:
            #python AltAnalyze.py --accessoryAnalysis IDTranslation --inputIDType Symbol --outputIDType RefSeq --input "C:\file1.txt" --species Hs
            inputIDType=None
            outputIDType=None
            for opt, arg in options: ### Accept user input for these hierarchical clustering variables
                if opt == '--inputIDType': inputIDType = arg
                if opt == '--outputIDType': outputIDType = arg
            if inputIDType==None or outputIDType==None:
                print 'Please designate an input ID type and and output ID type (--inputIDType Ensembl --outputIDType Symbol)'; sys.exit()
            if species==None:
                print "Please enter a valide species (--species)"; sys.exit()
            UI.IDconverter(input_file_dir, species, inputIDType, outputIDType, None)
            sys.exit()
            
        if 'hierarchical' in image_export:
            #python AltAnalyze.py --input "/Users/test/pluri.txt" --image hierarchical --row_method average --column_method single --row_metric cosine --column_metric euclidean --color_gradient red_white_blue --transpose False --PathwaySelection Apoptosis:WP254 --GeneSetSelection WikiPathways --species Hs --platform exon --display false
            if input_file_dir==None:
                print 'Please provide a valid file location for your input data matrix (must have an annotation row and an annotation column)';sys.exit()
            row_method = 'weighted'
            column_method = 'average'
            row_metric = 'cosine'
            column_metric = 'cosine'
            color_gradient = 'red_black_sky'
            contrast=2.5
            vendor = 'Affymetrix'
            GeneSelection = ''
            PathwaySelection = ''
            GeneSetSelection = 'None Selected'
            rho = None

            for opt, arg in options: ### Accept user input for these hierarchical clustering variables
                if opt == '--row_method':
                    row_method=arg
                    if row_method == 'None': row_method = None
                elif opt == '--column_method':
                    column_method=arg
                    if column_method == 'None': column_method = None
                elif opt == '--row_metric': row_metric=arg
                elif opt == '--column_metric': column_metric=arg
                elif opt == '--color_gradient': color_gradient=arg
                elif opt == '--GeneSetSelection': GeneSetSelection=arg
                elif opt == '--PathwaySelection': PathwaySelection=arg
                elif opt == '--genes': GeneSelection=arg
                elif opt == '--OntologyID': OntologyID=arg
                elif opt == '--normalization': normalization=arg
                elif opt == '--justShowTheseIDs': justShowTheseIDs=arg
                elif opt == '--rho': rho=arg
                elif opt == '--clusterGOElite':clusterGOElite=arg
                elif opt == '--contrast':
                    try: contrast=float(arg)
                    except Exception: print '--contrast not a valid float';sys.exit()
                elif opt == '--vendor': vendor=arg
                elif opt == '--display':
                    if arg=='yes':
                        display=True
                    elif arg=='True':
                        display=True
                    else:
                        display=False
            if len(GeneSetSelection)>0 or GeneSelection != '':
                gsp = UI.GeneSelectionParameters(species,array_type,vendor)
                gsp.setGeneSet(GeneSetSelection)
                gsp.setPathwaySelect(PathwaySelection)
                gsp.setGeneSelection(GeneSelection)
                gsp.setOntologyID(OntologyID)
                gsp.setTranspose(transpose)
                gsp.setNormalize(normalization)
                gsp.setJustShowTheseIDs(justShowTheseIDs)
                try: gsp.setClusterGOElite(clusterGOElite)
                except Exception: pass
                if rho!=None:
                    try:
                        float(rho)
                        gsp.setRhoCutoff(rho)
                    except Exception: print 'Must enter a valid Pearson correlation cutoff (float)'
                transpose = gsp ### this allows methods that don't transmit this object to also work
            if row_method == 'no': row_method = None
            if column_method == 'no': column_method = None
            if len(GeneSetSelection)>0:
                if species == None:
                    print "Please enter a valide species (--species)"; sys.exit()
            try:
                files = unique.read_directory(input_file_dir+'/')
                dir = input_file_dir
                for file in files:
                    filename = dir+'/'+file
                    UI.createHeatMap(filename, row_method, row_metric, column_method, column_metric, color_gradient, transpose, contrast, None, display=display)
            except Exception:
                UI.createHeatMap(input_file_dir, row_method, row_metric, column_method, column_metric, color_gradient, transpose, contrast, None, display=display)
            #from visualization_scripts import clustering; clustering.outputClusters([input_file_dir],[])
            sys.exit()
                        
        if 'PCA' in image_export or 't-SNE' in image_export or 'UMAP' in image_export or 'umap' in image_export:
            #AltAnalyze.py --input "/Users/nsalomonis/Desktop/folds.txt" --image PCA --plotType 3D --display True --labels yes
            #python AltAnalyze.py --input "/Users/nsalomonis/Desktop/log2_expression.txt" --image "t-SNE" --plotType 2D --display True --labels no --genes "ACTG2 ARHDIA KRT18 KRT8 ATP2B1 ARHGDIB" --species Hs --platform RNASeq --separateGenePlots True --zscore no
            #--algorithm "t-SNE"
            include_labels = 'yes'
            plotType = '2D'
            pca_algorithm = 'SVD'
            geneSetName = None
            zscore = True
            colorByGene=None
            separateGenePlots = False
            reimportModelScores = True
            if 't-SNE' in image_export:
                pca_algorithm = 't-SNE'
            if 'UMAP' in image_export or 'umap' in image_export:
                pca_algorithm = 'UMAP'
            for opt, arg in options: ### Accept user input for these hierarchical clustering variables
                #print opt,arg
                if opt == '--labels':
                    include_labels=arg
                    if include_labels == 'True' or include_labels == 'yes':
                        include_labels = 'yes'
                    else:
                        include_labels = 'no'
                if opt == '--plotType': plotType=arg
                if opt == '--algorithm': pca_algorithm=arg
                if opt == '--geneSetName': geneSetName=arg
                if opt == '--genes': colorByGene=arg
                if opt == '--reimportModelScores':
                    if arg == 'yes' or arg == 'True' or arg == 'true':
                        reimportModelScores = True
                    else:
                        reimportModelScores = False
                if opt == '--separateGenePlots':
                    if arg=='yes' or arg=='True' or arg == 'true':
                        separateGenePlots = True
                    else:
                        separateGenePlots = False
                if opt == '--zscore':
                    if arg=='yes' or arg=='True' or arg == 'true':
                        zscore=True
                    else:
                        zscore=False
                if opt == '--display':
                    if arg=='yes' or arg=='True' or arg == 'true':
                        display=True
            if input_file_dir==None:
                print 'Please provide a valid file location for your input data matrix (must have an annotation row and an annotation column)';sys.exit()
            UI.performPCA(input_file_dir, include_labels, pca_algorithm, transpose, None,
                          plotType=plotType, display=display, geneSetName=geneSetName, species=species, zscore=zscore, colorByGene=colorByGene, reimportModelScores=reimportModelScores, separateGenePlots=separateGenePlots)
            sys.exit()

        if 'VennDiagram' in image_export:
            # AltAnalyze.py --image "VennDiagram" --input "C:\file1.txt" --input "C:\file2.txt" --output "C:\graphs"
            files_to_merge=[]
            for opt, arg in options: ### Accept user input for these hierarchical clustering variables
                if opt == '--input':
                    arg = verifyPath(arg)
                    files_to_merge.append(arg)
                if opt == '--display':
                    if arg=='yes' or arg=='True' or arg == 'true':
                        display=True

            if len(files_to_merge)<2:
                print 'Please designate two or more files to compare (--input)';sys.exit()
            UI.vennDiagram(files_to_merge, output_dir, None, display=display)
            sys.exit()

        if 'AltExonViewer' in image_export:
            #python AltAnalyze.py --image AltExonViewer --AltResultsDir "C:\CP-hESC" --genes "ANXA7 FYN TCF3 NAV2 ETS2 MYLK ATP2A2" --species Hs --platform exon --dataType "splicing-index"
            genes=[]
            show_introns='no'
            geneFileDir=''
            analysisType='plot'
            for opt, arg in options: ### Accept user input for these hierarchical clustering variables
                if opt == '--genes':genes=arg
                elif opt == '--dataType': data_type = arg
                elif opt == '--showIntrons': show_introns = arg
                elif opt == '--AltResultsDir': altresult_dir = arg
                elif opt == '--geneFileDir': geneFileDir = arg
                elif opt == '--analysisType': analysisType=arg
                
            if altresult_dir == None:
                print 'Please include the location of the AltResults directory (--AltResultsDir)'; sys.exit()
            if len(genes)==0 and len(geneFileDir)==0:
                print "Please indicate the genes (--genes) or gene file location (--geneFileDir) for AltExonViewer";sys.exit()                
            if species == None:
                print "Please enter a valide species (--species)"; sys.exit()
            if array_type == None:
                print "Please enter a valide platform (--platform)"; sys.exit()
            if 'AltResults' not in altresult_dir:
                altresult_dir+='/AltResults/'

            if 'Sashimi' in analysisType:
                #python AltAnalyze.py --image AltExonViewer --AltResultsDir "/Users/saljh8/Desktop/Grimes/GEC14074/AltResults/" --genes "Dgat1 Dgat2 Tcf7l1" --species Mm --platform RNASeq --analysisType SashimiPlot
                analysisType = 'Sashimi-Plot'
                altresult_dir = string.split(altresult_dir,'AltResults')[0]
                if len(geneFileDir)>0: genes = geneFileDir
                geneFileDir=''
            elif 'raw' in data_type: ### Switch directories if expression
                altanalyze_results_folder = string.replace(altresult_dir,'AltResults','ExpressionInput')
                altresult_dir = UI.getValidExpFile(altanalyze_results_folder)
                if len(altresult_dir)==0:
                    print 'No valid expression input file (e.g., exp.MyExperiment.txt) found in',altanalyze_results_folder;sys.exit()
            else:
                altanalyze_results_folder = altresult_dir+'/RawSpliceData/'+species
                try: altresult_dir = UI.getValidSplicingScoreFile(altanalyze_results_folder)
                except Exception,e:
                    print "No files found in: "+altanalyze_results_folder; sys.exit()

            if len(geneFileDir)>0:
                try:
                    genes = UI.importGeneList(geneFileDir) ### list of gene IDs or symbols
                except Exception:
                    ### Can occur if a directory of files is selected
                    try:
                        files = unique.read_directory(geneFileDir+'/')
                        gene_string=''
                        for file in files:
                            if '.txt' in file:
                                filename = geneFileDir+'/'+file
                                genes = UI.importGeneList(filename) ### list of gene IDs or symbols
                                gene_string = gene_string+','+genes
                                print 'Imported genes from',file,'\n'
                                #print [altresult_dir];sys.exit()
                        UI.altExonViewer(species,platform,altresult_dir, gene_string, show_introns, analysisType, False)
                    except Exception: pass
                    sys.exit()
                        
            if len(genes)==0:
                print 'Please list one or more genes (--genes "ANXA7 FYN TCF3 NAV2 ETS2 MYLK ATP2A2")'; sys.exit()
            try: UI.altExonViewer(species,platform,altresult_dir, genes, show_introns, analysisType, False)
            except Exception:
                print traceback.format_exc() 
            sys.exit()
            
        if 'network' in image_export:
            #AltAnalyze.py --image network --species Hs --output "C:\GSE9440_RAW" --PathwaySelection Apoptosis:WP254 --GeneSetSelection WikiPathways
            for opt, arg in options: ### Accept user input for these hierarchical clustering variables
                if opt == '--update_interactions': update_interactions=arg
                elif opt == '--includeExpIDs': includeExpIDs=arg
                elif opt == '--degrees': degrees=arg
                elif opt == '--genes':
                    Genes=arg
                    inputType = 'IDs'
                elif opt == '--inputType': inputType=arg
                elif opt == '--interactionDirs': interactionDirs.append(arg)
                elif opt == '--GeneSetSelection': GeneSetSelection=arg
                elif opt == '--PathwaySelection': PathwaySelection=arg
                elif opt == '--OntologyID': OntologyID=arg
                elif opt == '--display': display=arg
            
            if update_interactions == 'yes': update_interactions = True
            else: update_interactions = False
            if input_file_dir == None: pass
            elif len(input_file_dir) == 0: input_file_dir = None
            if len(input_exp_file) == 0: input_exp_file = None
            if len(interactionDirs) == 0: interactionDirs=['WikiPathways']
            if interactionDirs == ['all']:
                interactionDirs = ['WikiPathways','KEGG','BioGRID','TFTargets','common-microRNATargets','all-microRNATargets','common-DrugBank','all-DrugBank']
            if interactionDirs == ['main']:
                interactionDirs = ['WikiPathways','KEGG','BioGRID','TFTargets']
            if interactionDirs == ['confident']:
                interactionDirs = ['WikiPathways','KEGG','TFTargets']

            if len(Genes) == 0: Genes = None
            if output_dir == None: pass
            elif len(output_dir) == 0: output_dir = None
            if len(GeneSetSelection) == 'None Selected': GeneSetSelection = None
            if includeExpIDs=='yes': includeExpIDs = True
            else: includeExpIDs = False

            gsp = UI.GeneSelectionParameters(species,array_type,manufacturer)
            gsp.setGeneSet(GeneSetSelection)
            gsp.setPathwaySelect(PathwaySelection)
            gsp.setGeneSelection(Genes)
            gsp.setOntologyID(OntologyID)
            gsp.setIncludeExpIDs(includeExpIDs)
            root = ''
            if species == None:
                print 'Please designate a species (--species).'; sys.exit()
            if output_dir == None:
                print 'Please designate an ouput directory (--output)'; sys.exit()
            if input_file_dir !=None:
                if '.txt' in input_file_dir or '.sif' in input_file_dir:
                    UI.networkBuilder(input_file_dir,inputType,output_dir,interactionDirs,degrees,input_exp_file,gsp,root)
                else:
                    parent_dir = input_file_dir
                    dir_list = read_directory(parent_dir)
                    for file in dir_list:
                        input_file_dir = parent_dir+'/'+file
                        try:
                            UI.networkBuilder(input_file_dir,inputType,output_dir,interactionDirs,degrees,input_exp_file,gsp,root)
                        except Exception:
                            print file, 'failed to produce network'
            else:
                UI.networkBuilder(None,inputType,output_dir,interactionDirs,degrees,input_exp_file,gsp,root)
            sys.exit()

    ########## Begin database dependent AltAnalyze workflows
    if ensembl_version != 'current' and 'markers' not in update_method:
        dbversion = string.replace(ensembl_version,'EnsMart','')
        UI.exportDBversion('EnsMart'+dbversion)
        gene_database = unique.getCurrentGeneDatabaseVersion()
        print 'Current database version:',gene_database
    if array_type == None and update_dbs != 'yes' and denom_file_dir == None:
        print "Please specify an array or data type (e.g., RNASeq, exon, gene, junction, AltMouse, 3'array)."; sys.exit()
    if 'archive' in update_method:
        ###
        print 'Archiving databases', ensembl_version
        try: archive_dir = 'ArchiveDBs/EnsMart'+ensembl_version+'/archive'; export.createDirPath(filepath(archive_dir))
        except Exception: null = [] ### directory already exists
        dirs = unique.read_directory('/ArchiveDBs/EnsMart'+ensembl_version)
        print len(dirs), dirs
        import shutil
        for species_dir in dirs:
            try:
                #print '/ArchiveDBs/EnsMart'+ensembl_version+'/'+species_dir+'/'+species_dir+'_RNASeq.zip'
                src = filepath('ArchiveDBs/EnsMart'+ensembl_version+'/'+species_dir+'/'+species_dir+'_RNASeq.zip')
                dstn = filepath('ArchiveDBs/EnsMart'+ensembl_version+'/archive/'+species_dir+'_RNASeq.zip')
                #export.copyFile(src, dstn)
                shutil.move(src, dstn)
                try:
                    srcj = string.replace(src,'RNASeq.','junction.'); dstnj = string.replace(dstn,'RNASeq.','junction.')
                    shutil.move(srcj, dstnj)
                except Exception: null=[]
                try:
                    src = string.replace(src,'_RNASeq.','.'); dstn = string.replace(dstn,'_RNASeq.','.')
                    shutil.move(src, dstn)
                except Exception:
                    print traceback.format_exc() 
                    pass
            except Exception:
                print traceback.format_exc() 
                pass
        sys.exit()
    if update_dbs == 'yes' and 'Official' not in update_method:
        if 'cleanup' in update_method:
            existing_species_dirs = unique.read_directory('/AltDatabase/ensembl')
            print 'Deleting EnsemblSQL directory for all species, ensembl version',ensembl_version
            for species in existing_species_dirs:
                export.deleteFolder('AltDatabase/ensembl/'+species+'/EnsemblSQL')
            existing_species_dirs = unique.read_directory('/AltDatabase')
            print 'Deleting SequenceData directory for all species, ensembl version',ensembl_version
            for species in existing_species_dirs:
                export.deleteFolder('AltDatabase/'+species+'/SequenceData')
            print 'Finished...exiting'
            sys.exit()
        if 'package' not in update_method and 'markers' not in update_method:
            ### Example: 
            ### python AltAnalyze.py --species all --arraytype all --update all --version 60
            ### tr -d \\r < AltAnalyze.py > AltAnalyze_new.py
            ### chmod +x AltAnalyze_new.py
            ### nohup ./AltAnalyze.py --update all --species Mm --arraytype gene --arraytype exon --version 60 2>&1 > nohup_v60_Mm.txt
        
            if array_type == 'all' and (species == 'Mm' or species == 'all'): array_type = ['AltMouse','exon','gene','junction','RNASeq']
            elif array_type == 'all' and (species == 'Hs' or species == 'Rn'): array_type = ['exon','gene','junction','RNASeq']
            else: array_type = [array_type]+additional_array_types
            if species == 'all' and 'RNASeq' not in array_type: species = selected_species ### just analyze the species for which multiple platforms are supported
            if species == 'selected': species = selected_species ### just analyze the species for which multiple platforms are supported
            elif species == 'all':
                all_supported_names = {}; all_species_names={}
                species_names = UI.getSpeciesInfo()
                for species in species_names: all_supported_names[species_names[species]]=species
                from build_scripts import EnsemblSQL
                child_dirs, ensembl_species, ensembl_versions = EnsemblSQL.getCurrentEnsemblSpecies('release-'+ensembl_version)
                for ens_species in ensembl_species:
                    ens_species = string.replace(ens_species,'_',' ')
                    if ens_species in all_supported_names:
                        all_species_names[all_supported_names[ens_species]]=[]
                
                del all_species_names['Hs']
                del all_species_names['Mm']
                del all_species_names['Rn']
                """
                del all_species_names['Go']
                del all_species_names['Bt']
                del all_species_names['Sc']
                
                del all_species_names['Ss']
                del all_species_names['Pv']
                del all_species_names['Pt']
                del all_species_names['La']
                del all_species_names['Tt']
                del all_species_names['Tr']
                del all_species_names['Ts']
                
                del all_species_names['Pb']
                del all_species_names['Pc']
                del all_species_names['Ec']
                
                del all_species_names['Tb']
                del all_species_names['Tg']
                del all_species_names['Dn']
                del all_species_names['Do']
                del all_species_names['Tn']
                del all_species_names['Dm']
                del all_species_names['Oc']
                del all_species_names['Og']
                
                del all_species_names['Fc']
                del all_species_names['Dr']
                del all_species_names['Me']
                
                del all_species_names['Cp']
                del all_species_names['Tt']
                del all_species_names['La']
                del all_species_names['Tr']
                del all_species_names['Ts']
                del all_species_names['Et'] ### No alternative isoforms?
                del all_species_names['Pc']
                del all_species_names['Tb']
                del all_species_names['Fc']
                del all_species_names['Sc']
                del all_species_names['Do']
                del all_species_names['Dn']
                del all_species_names['Og']
                del all_species_names['Ga']
                del all_species_names['Me']
                del all_species_names['Ml']
                del all_species_names['Mi']
                
                del all_species_names['St']
                del all_species_names['Sa']
                del all_species_names['Cs']
                del all_species_names['Vp']
                
                del all_species_names['Ch']
                del all_species_names['Ee']
                del all_species_names['Ac']"""
                sx=[]; all_species_names2=[] ### Ensure that the core selected species are run first
                for species in selected_species:
                    if species in all_species_names: sx.append(species)
                for species in all_species_names:
                    if species not in selected_species: all_species_names2.append(species)
                all_species_names = sx+all_species_names2
                species = all_species_names
            else: species = [species]
            
            update_uniprot='no'; update_ensembl='no'; update_probeset_to_ensembl='no'; update_domain='no'; update_miRs = 'no'; genomic_build = 'new'; update_miR_seq = 'yes'
            if 'all' in update_method:
                update_uniprot='yes'; update_ensembl='yes'; update_probeset_to_ensembl='yes'; update_domain='yes'; update_miRs = 'yes'
            if 'UniProt' in update_method: update_uniprot = 'yes'
            if 'Ensembl' in update_method: update_ensembl = 'yes'
            if 'Probeset' in update_method or 'ExonAnnotations' in update_method: update_probeset_to_ensembl = 'yes'
            if 'Domain' in update_method:
                update_domain = 'yes'
                try: from Bio import Entrez #test this
                except Exception: print 'The dependent module Bio is not installed or not accessible through the default python interpretter. Existing AltAnalyze.'; sys.exit()
            if 'miRBs' in update_method or 'miRBS' in update_method: update_miRs = 'yes'
            if 'NewGenomeBuild' in update_method: genomic_build = 'new'
            if 'current' in ensembl_version: print "Please specify an Ensembl version number (e.g., 60) before proceeding with the update.";sys.exit()
            try: force = force ### Variable is not declared otherwise
            except Exception: force = 'yes'; print 'force:',force
            existing_species_dirs={}
            update_all = 'no' ### We don't pass this as yes, in order to skip certain steps when multiple array types are analyzed (others are specified above)
            try: print "Updating AltDatabase the following array_types",string.join(array_type),"for the species",string.join(species)
            except Exception: print 'Please designate a valid platform/array_type (e.g., exon) and species code (e.g., Mm).'
            for specific_species in species:
                for platform_name in array_type:
                    if platform_name == 'AltMouse' and specific_species == 'Mm': proceed = 'yes'
                    elif platform_name == 'exon' or platform_name == 'gene':
                        from build_scripts import ExonArrayEnsemblRules
                        #### Check to see if the probeset.csv file is present
                        #try: probeset_transcript_file = ExonArrayEnsemblRules.getDirectoryFiles('/AltDatabase/'+specific_species+'/'+platform_name)
                        #except Exception: print "Affymetrix probeset.csv anotation file is not found. You must save this to",'/AltDatabase/'+specific_species+'/'+platform_name,'before updating (unzipped).'; sys.exit()
                        proceed = 'yes'
                    elif platform_name == 'junction' and (specific_species == 'Hs' or specific_species == 'Mm'): proceed = 'yes'
                    elif platform_name == 'RNASeq': proceed = 'yes'
                    else: proceed = 'no'
                    if proceed == 'yes':
                        print "Analyzing", specific_species, platform_name
                        if (platform_name != array_type[0]) and len(species)==1:
                            update_uniprot = 'no'; update_ensembl = 'no'; update_miR_seq = 'no' ### Don't need to do this twice in a row
                            print 'Skipping ensembl, uniprot and mir-sequence file import updates since already completed for this species',array_type,platform_name
                        if ignore_built_species == 'yes': ### Useful for when building all species for a new database build
                            try:
                                ### call this here to update with every species - if running multiple instances
                                existing_species_dirs = unique.read_directory('/AltDatabase/ensembl') 
                            except Exception: #ZeroDivisionError
                                try: os.mkdir(unique.filepath('AltDatabase/'))
                                except Exception: pass #already exists
                                existing_species_dirs = []
                        if specific_array_type != None and specific_array_type != platform_name: platform_name+='|'+specific_array_type ### For the hGlue vs. JAY arrays
                        if specific_species not in existing_species_dirs: ### Useful when running multiple instances of AltAnalyze to build all species
                            print 'update_ensembl',update_ensembl
                            print 'update_uniprot',update_uniprot
                            print 'update_probeset_to_ensembl',update_probeset_to_ensembl
                            print 'update_domain',update_domain
                            print 'update_miRs',update_miRs
                            update.executeParameters(specific_species,platform_name,force,genomic_build,update_uniprot,update_ensembl,update_probeset_to_ensembl,update_domain,update_miRs,update_all,update_miR_seq,ensembl_version)
                        else: print 'ignoring',specific_species 
            sys.exit()
    
    if 'package' in update_method:
        ### Example: python AltAnalyze.py --update package --species all --platform all --version 65
        
        if ensembl_version == 'current': print '\nPlease specify version of the database to package (e.g., --version 60).'; sys.exit()
        ensembl_version = 'EnsMart'+ensembl_version

        ### Get all possible species
        species_names = UI.getSpeciesInfo(); possible_species={}
        possible_species = species_names
        possible_arrays = ['exon','gene','junction','AltMouse','RNASeq']
        try:
            if species == 'all': possible_species = possible_species
            elif species == 'selected': possible_species = selected_species
            else: possible_species = [species]
        except Exception: species = possible_species
        if array_type == None or array_type == 'all': possible_arrays = possible_arrays
        else: possible_arrays = [array_type]+additional_array_types
        species_to_package={}
        
        dirs = unique.read_directory('/AltDatabase/'+ensembl_version)
        #print possible_arrays, possible_species; sys.exit()
        for species_code in dirs:
            if species_code in possible_species:
                array_types = unique.read_directory('/AltDatabase/'+ensembl_version+'/'+species_code)
                for arraytype in array_types:
                    if arraytype in possible_arrays:
                        if species_code in possible_species:
                            array_types = unique.read_directory('/AltDatabase/'+ensembl_version+'/'+species_code)
                            try: species_to_package[species_code].append(arraytype)
                            except Exception: species_to_package[species_code] = [arraytype]
        species_to_package = eliminate_redundant_dict_values(species_to_package)

        for species in species_to_package:
            files_to_copy =[species+'_Ensembl_domain_aligning_probesets.txt']
            files_to_copy+=[species+'_Ensembl_indirect_domain_aligning_probesets.txt']
            files_to_copy+=[species+'_Ensembl_probesets.txt']
            files_to_copy+=[species+'_Ensembl_exons.txt']
            #files_to_copy+=[species+'_Ensembl_junctions.txt']
            files_to_copy+=[species+'_exon_core.mps']
            files_to_copy+=[species+'_exon_extended.mps']
            files_to_copy+=[species+'_exon_full.mps']
            files_to_copy+=[species+'_gene_core.mps']
            files_to_copy+=[species+'_gene_extended.mps']
            files_to_copy+=[species+'_gene_full.mps']
            files_to_copy+=[species+'_gene-exon_probesets.txt']
            files_to_copy+=[species+'_probes_to_remove.txt']
            files_to_copy+=[species+'_probeset-probes.txt']
            files_to_copy+=[species+'_probeset_microRNAs_any.txt']
            files_to_copy+=[species+'_probeset_microRNAs_multiple.txt']
            files_to_copy+=['probeset-domain-annotations-exoncomp.txt']
            files_to_copy+=['probeset-protein-annotations-exoncomp.txt']
            #files_to_copy+=['probeset-protein-dbase_exoncomp.txt']
            files_to_copy+=['SEQUENCE-protein-dbase_exoncomp.txt']
            files_to_copy+=[species+'_Ensembl_junction_probesets.txt']
            files_to_copy+=[species+'_Ensembl_AltMouse_probesets.txt']
            files_to_copy+=[species+'_RNASeq-exon_probesets.txt']
            files_to_copy+=[species+'_junction-exon_probesets.txt']
            files_to_copy+=[species+'_junction_all.mps']
            files_to_copy+=['platform.txt'] ### Indicates the specific platform for an array type (e.g., HJAY for junction or hGlue for junction)
            files_to_copy+=[species+'_junction_comps_updated.txt']
            files_to_copy+=['MASTER-probeset-transcript.txt']
            files_to_copy+=['AltMouse-Ensembl.txt']
            files_to_copy+=['AltMouse_junction-comparisons.txt']
            files_to_copy+=['AltMouse_gene_annotations.txt']
            files_to_copy+=['AltMouse_annotations.txt']
            
            common_to_copy =['uniprot/'+species+'/custom_annotations.txt']
            common_to_copy+=['ensembl/'+species+'/'+species+'_Ensembl-annotations_simple.txt']
            common_to_copy+=['ensembl/'+species+'/'+species+'_Ensembl-annotations.txt']
            common_to_copy+=['ensembl/'+species+'/'+species+'_microRNA-Ensembl.txt']
            common_to_copy+=['ensembl/'+species+'/'+species+'_Ensembl_transcript-biotypes.txt']
            common_to_copy+=['ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt']
            common_to_copy+= searchDirectory("AltDatabase/ensembl/"+species+"/",'Ensembl_Protein')
            common_to_copy+= searchDirectory("AltDatabase/ensembl/"+species+"/",'ProteinFeatures')
            common_to_copy+= searchDirectory("AltDatabase/ensembl/"+species+"/",'ProteinCoordinates')
            common_to_copy+= searchDirectory("AltDatabase/uniprot/"+species+"/",'FeatureCoordinate')
            supported_arrays_present = 'no'
            for arraytype in selected_platforms:
                if arraytype in species_to_package[species]: supported_arrays_present = 'yes' #Hence a non-RNASeq platform is present
            if supported_arrays_present == 'yes':
                for file in common_to_copy:
                    ir = 'AltDatabase/'+ensembl_version+'/'
                    er = 'ArchiveDBs/'+ensembl_version+'/'+species+'/'+ensembl_version+'/'
                    export.copyFile(ir+file, er+file)
                    
            
            if 'RNASeq' in species_to_package[species]:
                common_to_copy+=['ensembl/'+species+'/'+species+'_Ensembl_junction.txt']
                common_to_copy+=['ensembl/'+species+'/'+species+'_Ensembl_exon.txt']
                for file in common_to_copy:
                    ir = 'AltDatabase/'+ensembl_version+'/'
                    er = 'ArchiveDBs/'+ensembl_version+'/'+species+'/'+ensembl_version+'/'
                    if species in selected_species:
                        er = 'ArchiveDBs/'+ensembl_version+'/'+species+'/RNASeq/'+ensembl_version+'/' ### This allows us to build the package archive in a separate directory for selected species, so separate but overlapping content can be packaged
                    export.copyFile(ir+file, er+file)

            for array_type in species_to_package[species]:
                ir = 'AltDatabase/'+ensembl_version+'/'+species+'/'+array_type+'/'
                er = 'ArchiveDBs/'+ensembl_version+'/'+species+'/'+ensembl_version+'/'+species+'/'+array_type+'/'
                if array_type == 'junction':
                    er = 'ArchiveDBs/'+ensembl_version+'/'+species+'/'+array_type+'/'
                if array_type == 'RNASeq' and species in selected_species:
                    er = 'ArchiveDBs/'+ensembl_version+'/'+species+'/RNASeq/'+ensembl_version+'/'+species+'/'+array_type+'/'
                for file in files_to_copy:
                    if array_type == 'RNASeq': file=string.replace(file,'_updated.txt','.txt')
                    filt_file = string.replace(file ,'.txt','-filtered.txt')
                    try: export.copyFile(ir+filt_file, er+filt_file); export_path = er+filt_file
                    except Exception:
                        try: export.copyFile(ir+file, er+file); export_path = er+file
                        except Exception: null = [] ### File not found in directory
                    if len(export_path)>0:
                        if 'AltMouse' in export_path or 'probes_' in export_path:
                            export.cleanFile(export_path)
                if array_type == 'junction':
                    subdir = '/exon/'
                    ir = 'AltDatabase/'+ensembl_version+'/'+species+'/'+array_type+subdir
                    er = 'ArchiveDBs/'+ensembl_version+'/'+species+'/'+array_type+subdir
                    for file in files_to_copy:
                        export_path=[]
                        filt_file = string.replace(file ,'.txt','-filtered.txt')
                        try: export.copyFile(ir+filt_file, er+filt_file); export_path = er+filt_file
                        except Exception:
                            try: export.copyFile(ir+file, er+file); export_path = er+file
                            except Exception: null = [] ### File not found in directory
                if array_type == 'RNASeq':
                    subdir = '/junction/'
                    ir = 'AltDatabase/'+ensembl_version+'/'+species+'/'+array_type+subdir
                    er = 'ArchiveDBs/'+ensembl_version+'/'+species+'/'+ensembl_version+'/'+species+'/'+array_type+subdir
                    if species in selected_species:
                        er = 'ArchiveDBs/'+ensembl_version+'/'+species+'/RNASeq/'+ensembl_version+'/'+species+'/'+array_type+subdir
                    for file in files_to_copy:
                        if 'SEQUENCE-protein-dbase' not in file and 'domain_aligning' not in file: ### This data is now combined into the main file
                            export_path=[]
                            filt_file = string.replace(file ,'.txt','-filtered.txt')
                            try: export.copyFile(ir+filt_file, er+filt_file); export_path = er+filt_file
                            except Exception:
                                try: export.copyFile(ir+file, er+file); export_path = er+file
                                except Exception: null = [] ### File not found in directory
                            
            if 'RNASeq' in species_to_package[species]:
                src = 'ArchiveDBs/'+ensembl_version+'/'+species+'/'+ensembl_version
                dst = 'ArchiveDBs/'+ensembl_version+'/'+species+'/'+species+'_RNASeq.zip'
                if species in selected_species:
                    src = 'ArchiveDBs/'+ensembl_version+'/'+species+'/RNASeq/'+ensembl_version
                update.zipDirectory(src); print 'Zipping',species, array_type, dst
                os.rename(src+'.zip', dst)
            if supported_arrays_present == 'yes':
                src = 'ArchiveDBs/'+ensembl_version+'/'+species+'/'+ensembl_version
                dst = 'ArchiveDBs/'+ensembl_version+'/'+species+'/'+species+'.zip'
                update.zipDirectory(src); print 'Zipping',species, array_type, dst
                os.rename(src+'.zip', dst)
            if 'junction' in species_to_package[species]:
                src = 'ArchiveDBs/'+ensembl_version+'/'+species+'/junction'
                dst = string.replace(src,'junction',species+'_junction.zip')
                update.zipDirectory(src); print 'Zipping',species+'_junction'
                os.rename(src+'.zip', dst)
        sys.exit()
        
    if 'markers' in update_method:
        if species == None or platform == None:
            print "WARNING! A species and platform (e.g., exon, junction, 3'array or RNASeq) must be defined to identify markers.";sys.exit()
        elif input_exp_file == '':
            print "WARNING! A input expression file must be supplied (e.g., ExpressionOutput/DATASET.YourExperimentName.txt) for this analysis.";sys.exit()
        else:
            #python AltAnalyze.py --update markers --platform gene --expdir "/home/socr/c/users2/salomoni/other/boxer/normalization/Mm_Gene-TissueAtlas/ExpressionInput/exp.meta.txt"
            #python AltAnalyze.py --update markers --platform gene --expdir "/home/socr/c/users2/salomoni/other/boxer/normalization/Mm_Gene-TissueAtlas/AltResults/RawSpliceData/Mm/splicing-index/meta.txt"
            #python AltAnalyze.py --update markers --platform "3'array" --expdir "/home/socr/c/users2/salomoni/other/boxer/normalization/U133/ExpressionOutput/DATASET-meta.txt"
            #python AltAnalyze.py --update markers --compendiumType ncRNA --platform "exon" --expdir "/home/socr/c/users2/salomoni/conklin/nsalomonis/normalization/Hs_Exon-TissueAtlas/ExpressionOutput/DATASET-meta.txt"
            #python AltAnalyze.py --update markers --platform RNASeq --species Mm --geneRPKM 1 --expdir /Users/saljh8/Desktop/Grimes/MergedRSEM/DN-Analysis/ExpressionInput/exp.DN.txt --genesToReport 200
            """The markerFinder module:
            1) takes an input ExpressionOutput file (DATASET.YourExperimentName.txt)
            2) extracts group average expression and saves to AVERAGE.YourExperimentName.txt to the ExpressionOutput directory
            3) re-imports AVERAGE.YourExperimentName.txt
            4) correlates the average expression of each gene to an idealized profile to derive a Pearson correlation coefficient
            5) identifies optimal markers based on these correlations for each tissue
            6) exports an expression file with just these marker genes and tissues

            This module can peform these analyses on protein coding or ncRNAs and can segregate the cell/tissue groups into clusters
            when a group notation is present in the sample name (e.g., 0~Heart, 0~Brain, 1~Stem Cell)"""
            
            import markerFinder
            if 'AltResults' in input_exp_file and 'Clustering' not in input_exp_file:
                ### This applies to a file compoosed of exon-level normalized intensities (calculae average group expression)
                markerFinder.getAverageExonExpression(species,platform,input_exp_file)
                if 'Raw' in input_exp_file:
                    group_exp_file = string.replace(input_exp_file,'Raw','AVERAGE')
                else:
                    group_exp_file = string.replace(input_exp_file,'FullDatasets','AVERAGE-FullDatasets')
                altexon_correlation_file = markerFinder.analyzeData(group_exp_file,species,platform,compendiumType,geneToReport=genesToReport,correlateAll=correlateAll,AdditionalParameters=fl)
                markerFinder.getExprValsForNICorrelations(platform,altexon_correlation_file,group_exp_file)
            else:
                ### This applies to an ExpressionOutput DATASET file compoosed of gene expression values (averages already present)
                import collections
                try: test_ordereddict=collections.OrderedDict()
                except Exception:
                    try: import ordereddict
                    except Exception:
                        ### This is needed to re-order the average file so that the groups are sequentially ordered when analyzing clustered groups (0~)
                        print 'Warning!!!! To run markerFinder correctly call python version 2.7x or greater (python 3.x not supported)'
                        print 'Requires ordereddict (also can install the library ordereddict). To call 2.7: /usr/bin/python2.7'
                        sys.exit()
                try:
                    output_dir = markerFinder.getAverageExpressionValues(input_exp_file,platform) ### Either way, make an average annotated file from the DATASET file
                    if 'DATASET' in input_exp_file:
                        group_exp_file = string.replace(input_exp_file,'DATASET','AVERAGE')
                    else:
                        group_exp_file = (input_exp_file,output_dir) ### still analyze the primary sample
                except Exception:
                    #print traceback.format_exc()
                    print 'No DATASET file present (used to obtain gene annotations)...'
                    ### Work around when performing this analysis on an alternative exon input cluster file
                    group_exp_file = input_exp_file
                fl = UI.ExpressionFileLocationData(input_exp_file,'','',''); fl.setOutputDir(export.findParentDir(export.findParentDir(input_exp_file)[:-1]))
                try: fl.setSpecies(species); fl.setVendor(vendor)
                except Exception: pass
                try:
                    rpkm_threshold = float(rpkm_threshold) ### If supplied, for any platform, use it
                    fl.setRPKMThreshold(rpkm_threshold)
                except Exception: pass
                if platform=='RNASeq':
                    try: rpkm_threshold = float(rpkm_threshold)
                    except Exception: rpkm_threshold = 1.0
                    fl.setRPKMThreshold(rpkm_threshold)

                try: correlationDirection = direction ### correlate to a positive or inverse negative in silico artificial pattern
                except Exception: correlationDirection = 'up'
                fl.setCorrelationDirection(correlationDirection)
        
                if expression_data_format == 'non-log': logTransform = True
                else: logTransform = False
                if 'topSplice' in input_exp_file:
                    markerFinder.filterRNASeqSpliceEvents(species,platform,fl,input_exp_file)
                    sys.exit()
                if 'stats.' in input_exp_file:
                    markerFinder.filterDetectionPvalues(species,platform,fl,input_exp_file)
                    sys.exit()
                else:
                    markerFinder.analyzeData(group_exp_file,species,platform,compendiumType,geneToReport=genesToReport,correlateAll=correlateAll,AdditionalParameters=fl,logTransform=logTransform)
                try: fl.setVendor(manufacturer)
                except Exception:
                    print '--vendor not indicated by user... assuming Affymetrix'
                    fl.setVendor('Affymetrix')
                try: markerFinder.generateMarkerHeatMaps(fl,array_type,convertNonLogToLog=logTransform,Species=species)
                except Exception: print traceback.format_exc() 
            print 'Cell/Tissue marker classification analysis finished';sys.exit()
    
    if 'EnsMart' in ensembl_version:
        UI.exportDBversion(ensembl_version)
    
    annotation_found = verifyFile(input_annotation_file)
    proceed = 'no'
    
    if 'Official' not in update_method and denom_file_dir == None: ### If running GO-Elite independent of AltAnalyze (see below GO_Elite call)
        try:
            time_stamp = timestamp()
            if len(cel_file_dir)>0:
                if output_dir == None:
                    output_dir = cel_file_dir
                    print "Setting output directory to the input path:", output_dir
            if output_dir == None and input_filtered_dir>0:
                output_dir = input_filtered_dir
            try:
                if '/' == output_dir[-1] or '\\' in output_dir[-2]: null=[]
                else: output_dir +='/'
            except:
                try: output_dir = export.findParentDir(input_file_dir)
                except:
                    output_dir = input_fastq_dir
            log_file = filepath(output_dir+'/AltAnalyze_report-'+time_stamp+'.log')
            log_report = open(log_file,'w'); log_report.close()
            sys.stdout = Logger('')
        except Exception,e:
            print traceback.format_exc() 
            print 'Please designate an output directory before proceeding (e.g., --output "C:\RNASeq)';sys.exit()

    try: print string.join(arguments,' ')
    except Exception: pass
            
    if mappedExonAnalysis:
        array_type = 'RNASeq' ### Although this is not the actual platform, the resulting data will be treated as RNA-Seq with parameters most suitable for arrays

    if len(external_annotation_dir)>0:
        run_from_scratch = 'Annotate External Results'
    if channel_to_extract != None:
        run_from_scratch = 'Process Feature Extraction files' ### Agilent Feature Extraction files as input for normalization
        manufacturer = 'Agilent'
        constitutive_source = 'Agilent'
        expression_threshold = 'NA'
        perform_alt_analysis = 'NA'
    if len(input_filtered_dir)>0:
        run_from_scratch ='Process AltAnalyze filtered'; proceed='yes'
    if len(input_exp_file)>0:
        run_from_scratch = 'Process Expression file'; proceed='yes'
        input_exp_file = string.replace(input_exp_file,'\\','/') ### Windows convention is \ rather than /, but works with /
        ief_list = string.split(input_exp_file,'/')
        if len(output_dir)>0: parent_dir = output_dir
        else: parent_dir = string.join(ief_list[:-1],'/')
        exp_name = ief_list[-1]


    if len(cel_file_dir)>0 or runKallisto == True:
        # python AltAnalyze.py --species Mm --platform RNASeq --runKallisto yes --expname test
        if exp_name == None:
            print "No experiment name defined. Please sumbit a name (e.g., --expname CancerComp) before proceeding."; sys.exit()
        else:
            dataset_name = 'exp.'+exp_name+'.txt'; exp_file_dir = filepath(output_dir+'/ExpressionInput/'+dataset_name)
        if runKallisto:
            run_from_scratch == 'Process RNA-seq reads'
        elif run_from_scratch!= 'Process Feature Extraction files':
            run_from_scratch = 'Process CEL files'; proceed='yes'
            if array_type == 'RNASeq': file_ext = '.BED'
            else: file_ext = '.CEL'
            try: cel_files,cel_files_fn = UI.identifyCELfiles(cel_file_dir,array_type,manufacturer)
            except Exception,e:
                print e
                if mappedExonAnalysis: pass
                else: print "No",file_ext,"files found in the directory:",cel_file_dir;sys.exit()
            if array_type != 'RNASeq': cel_file_list_dir = UI.exportCELFileList(cel_files_fn,cel_file_dir)

        if groups_file != None and comps_file != None:
            try: export.copyFile(groups_file, string.replace(exp_file_dir,'exp.','groups.'))
            except Exception: print 'Groups file already present in target location OR bad input path.'
            try: export.copyFile(comps_file, string.replace(exp_file_dir,'exp.','comps.'))
            except Exception: print 'Comparison file already present in target location OR bad input path.'
            groups_file = string.replace(exp_file_dir,'exp.','groups.')
            comps_file = string.replace(exp_file_dir,'exp.','comps.')
            if verifyGroupFileFormat(groups_file) == False:
                print "\nWarning! The format of your groups file is not correct. For details, see:\nhttp://altanalyze.readthedocs.io/en/latest/ManualGroupsCompsCreation\n"
                sys.exit()

        if array_type != 'RNASeq' and manufacturer!= 'Agilent':
            """Determine if Library and Annotations for the array exist, if not, download or prompt for selection"""
            try:
                ### For the HGLUE and HJAY arrays, this step is critical in order to have the commond-line AltAnalyze downloadthe appropriate junction database (determined from specific_array_type)
                specific_array_types,specific_array_type = UI.identifyArrayType(cel_files_fn)
                num_array_types = len(specific_array_types)
            except Exception:
                null=[]; num_array_types=1; specific_array_type=None
                if array_type == 'exon':
                    if species == 'Hs': specific_array_type = 'HuEx-1_0-st-v2'
                    if species == 'Mm': specific_array_type = 'MoEx-1_0-st-v2'
                    if species == 'Rn': specific_array_type = 'RaEx-1_0-st-v2'
                elif array_type == 'gene':
                    if species == 'Hs': specific_array_type = 'HuGene-1_0-st-v1'
                    if species == 'Mm': specific_array_type = 'MoGene-1_0-st-v1'
                    if species == 'Rn': specific_array_type = 'RaGene-1_0-st-v1'
                elif array_type == 'AltMouse': specific_array_type = 'altMouseA'
                """
                elif array_type == 'junction':
                    if species == 'Mm': specific_array_type = 'MJAY'
                    if species == 'Hs': specific_array_type = 'HJAY'
                """
            supproted_array_db = UI.importSupportedArrayInfo()

            if specific_array_type in supproted_array_db and input_cdf_file == None and input_annotation_file == None:
                sa = supproted_array_db[specific_array_type]; species = sa.Species(); array_type = sa.ArrayType()
                input_cdf_file, input_annotation_file, bgp_file, clf_file = UI.getAffyFilesRemote(specific_array_type,array_type,species)
            else: array_type = "3'array"

            cdf_found = verifyFile(input_cdf_file)
            annotation_found = verifyFile(input_annotation_file)
            if input_cdf_file == None:
                print [specific_array_type], 'not currently supported... Please provide CDF to AltAnalyze (commandline or GUI) or manually add to AltDatabase/affymetrix/LibraryFiles'; sys.exit()
                
            if cdf_found != "found":
                ### Copy valid Library files to a local AltAnalyze database directory
                input_cdf_file_lower = string.lower(input_cdf_file)
                if array_type == "3'array":
                    if '.cdf' in input_cdf_file_lower:
                        clf_file='';bgp_file=''; assinged = 'yes'
                        ###Thus the CDF or PDF file was confirmed, so copy it over to AltDatabase
                        icf_list = string.split(input_cdf_file,'/'); cdf_short = icf_list[-1]
                        destination_parent = 'AltDatabase/affymetrix/LibraryFiles/'
                        destination_parent = osfilepath(destination_parent+cdf_short)
                        info_list = input_cdf_file,destination_parent; UI.StatusWindow(info_list,'copy')
                    else: print "Valid CDF file not found. Exiting program.";sys.exit()          
                else:
                    if '.pgf' in input_cdf_file_lower:
                        ###Check to see if the clf and bgp files are present in this directory 
                        icf_list = string.split(input_cdf_file,'/'); parent_dir = string.join(icf_list[:-1],'/'); cdf_short = icf_list[-1]
                        clf_short = string.replace(cdf_short,'.pgf','.clf')
                        kil_short = string.replace(cdf_short,'.pgf','.kil') ### Only applies to the Glue array
                        if array_type == 'exon' or array_type == 'junction': bgp_short = string.replace(cdf_short,'.pgf','.antigenomic.bgp')
                        else: bgp_short = string.replace(cdf_short,'.pgf','.bgp')
                        dir_list = read_directory(parent_dir)
                        if clf_short in dir_list and bgp_short in dir_list:
                            pgf_file = input_cdf_file
                            clf_file = string.replace(pgf_file,'.pgf','.clf')
                            kil_file = string.replace(pgf_file,'.pgf','.kil') ### Only applies to the Glue array
                            if array_type == 'exon' or array_type == 'junction': bgp_file = string.replace(pgf_file,'.pgf','.antigenomic.bgp')
                            else: bgp_file = string.replace(pgf_file,'.pgf','.bgp')
                            assinged = 'yes'
                            ###Thus the CDF or PDF file was confirmed, so copy it over to AltDatabase
                            destination_parent = 'AltDatabase/affymetrix/LibraryFiles/'
                            info_list = input_cdf_file,osfilepath(destination_parent+cdf_short); UI.StatusWindow(info_list,'copy')
                            info_list = clf_file,osfilepath(destination_parent+clf_short); UI.StatusWindow(info_list,'copy')
                            info_list = bgp_file,osfilepath(destination_parent+bgp_short); UI.StatusWindow(info_list,'copy')
                            if 'Glue' in pgf_file:
                                info_list = kil_file,osfilepath(destination_parent+kil_short); UI.StatusWindow(info_list,'copy')
    if annotation_found != "found" and update_dbs == 'no' and array_type != 'RNASeq' and denom_file_dir == None and manufacturer != 'Agilent':
        ### Copy valid Annotation files to a local AltAnalyze database directory
        try:
            input_annotation_lower = string.lower(input_annotation_file)
            if '.csv' in input_annotation_lower:
                assinged = 'yes'
                ###Thus the CDF or PDF file was confirmed, so copy it over to AltDatabase
                icf_list = string.split(input_annotation_file,'/'); csv_short = icf_list[-1]
                destination_parent = 'AltDatabase/affymetrix/'+species+'/'
                info_list = input_annotation_file,filepath(destination_parent+csv_short); UI.StatusWindow(info_list,'copy')
        except Exception: print "No Affymetrix annotation file provided. AltAnalyze will use any .csv annotations files in AltDatabase/Affymetrix/"+species

    if array_type == 'PSI':
        array_type = "3'array"
        vendor = 'PSI'
    if 'Official' in update_method and species != None:
        proceed = 'yes'
    elif array_type != None and species != None:
        expr_defaults, alt_exon_defaults, functional_analysis_defaults, goelite_defaults = UI.importDefaults(array_type,species)
        ge_fold_cutoffs, ge_pvalue_cutoffs, ge_ptype, filter_method, z_threshold, p_val_threshold, change_threshold, ORA_algorithm, resources_to_analyze, goelite_permutations, mod, returnPathways, NA = goelite_defaults
        use_direct_domain_alignments_only,microRNA_prediction_method = functional_analysis_defaults
        analysis_method, additional_algorithms, filter_probeset_types, analyze_all_conditions, p_threshold, alt_exon_fold_variable, additional_score, permute_p_threshold, gene_expression_cutoff, remove_intronic_junctions, perform_permutation_analysis, export_NI_values, run_MiDAS, calculate_normIntensity_p, filter_for_AS = alt_exon_defaults
        dabg_p, rpkm_threshold, gene_exp_threshold, exon_exp_threshold, exon_rpkm_threshold, expression_threshold, perform_alt_analysis, analyze_as_groups, expression_data_format, normalize_feature_exp, normalize_gene_data, avg_all_for_ss, include_raw_data, probability_statistic, FDR_statistic, batch_effects, marker_finder, visualize_qc_results, run_lineage_profiler, null = expr_defaults
    elif denom_file_dir != None and species != None:
        proceed = 'yes' ### Only run GO-Elite
        expr_defaults, alt_exon_defaults, functional_analysis_defaults, goelite_defaults = UI.importDefaults('RNASeq',species) ### platform not relevant
        ge_fold_cutoffs, ge_pvalue_cutoffs, ge_ptype, filter_method, z_threshold, p_val_threshold, change_threshold, ORA_algorithm, resources_to_analyze, goelite_permutations, mod, returnPathways, NA = goelite_defaults
    else:
        print 'No species defined. Please include the species code (e.g., "--species Hs") and array type (e.g., "--arraytype exon") before proceeding.'
        print '\nAlso check the printed arguments above to see if there are formatting errors, such as bad quotes.'; sys.exit()     

    array_type_original = array_type
    #if array_type == 'gene': array_type = "3'array"
    for opt, arg in options:
        if opt == '--runGOElite': run_GOElite=arg
        elif opt == '--outputQCPlots': visualize_qc_results=arg
        elif opt == '--runLineageProfiler' or opt == '--cellHarmony' or opt == '--cellHarmonyMerge':
            if string.lower(arg) == 'yes' or string.lower(arg) == 'true':
                run_lineage_profiler = 'yes'
        elif opt == '--elitepermut': goelite_permutations=arg
        elif opt == '--method': filter_method=arg
        elif opt == '--zscore': z_threshold=arg
        elif opt == '--elitepval': p_val_threshold=arg
        elif opt == '--num': change_threshold=arg
        elif opt == '--dataToAnalyze': resources_to_analyze=arg
        elif opt == '--GEelitepval': ge_pvalue_cutoffs=arg
        elif opt == '--GEelitefold': ge_fold_cutoffs=arg
        elif opt == '--GEeliteptype': ge_ptype=arg
        elif opt == '--ORAstat': ORA_algorithm=arg
        elif opt == '--returnPathways': returnPathways=arg
        
        elif opt == '--FDR': FDR_statistic=arg
        elif opt == '--dabgp': dabg_p=arg
        elif opt == '--rawexp': expression_threshold=arg
        elif opt == '--geneRPKM': rpkm_threshold=arg
        elif opt == '--exonRPKM': exon_rpkm_threshold=arg
        elif opt == '--geneExp': gene_exp_threshold=arg
        elif opt == '--exonExp': exon_exp_threshold=arg
        elif opt == '--groupStat': probability_statistic=arg
        elif opt == '--avgallss': avg_all_for_ss=arg
        elif opt == '--logexp': expression_data_format=arg
        elif opt == '--inclraw': include_raw_data=arg
        elif opt == '--combat': batch_effects=arg
        elif opt == '--runalt': perform_alt_analysis=arg
        elif opt == '--altmethod': analysis_method=arg
        elif opt == '--altp': p_threshold=arg
        elif opt == '--probetype': filter_probeset_types=arg
        elif opt == '--altscore': alt_exon_fold_variable=arg
        elif opt == '--GEcutoff': gene_expression_cutoff=arg
        elif opt == '--removeIntronOnlyJunctions': remove_intronic_junctions=arg
        elif opt == '--normCounts': normalize_feature_exp=arg
        elif opt == '--normMatrix': normalize_gene_data=arg
        elif opt == '--altpermutep': permute_p_threshold=arg
        elif opt == '--altpermute': perform_permutation_analysis=arg
        elif opt == '--exportnormexp': export_NI_values=arg
        elif opt == '--buildExonExportFile': build_exon_bedfile = 'yes'
        elif opt == '--runMarkerFinder': marker_finder = arg
        elif opt == '--calcNIp': calculate_normIntensity_p=arg
        elif opt == '--runMiDAS': run_MiDAS=arg
        elif opt == '--analyzeAllGroups':
            analyze_all_conditions=arg
            if analyze_all_conditions == 'yes': analyze_all_conditions = 'all groups'
        elif opt == '--GEcutoff': use_direct_domain_alignments_only=arg
        elif opt == '--mirmethod': microRNA_prediction_method=arg
        elif opt == '--ASfilter': filter_for_AS=arg
        elif opt == '--noxhyb': xhyb_remove=arg
        elif opt == '--returnAll': return_all=arg
        elif opt == '--annotatedir': external_annotation_dir=arg
        elif opt == '--additionalScore': additional_score=arg
        elif opt == '--additionalAlgorithm': additional_algorithms=arg
        elif opt == '--modelSize':
            modelSize=arg
            try: modelSize = int(modelSize)
            except Exception: modelSize = None
        elif opt == '--geneModel':
            geneModel=arg # file location
            if geneModel == 'no' or 'alse' in geneModel:
                geneModel = False
        elif opt == '--reference':
            custom_reference = arg
        
    if run_from_scratch == 'Process Feature Extraction files': ### Agilent Feature Extraction files as input for normalization
        normalize_gene_data = 'quantile' ### required for Agilent
        proceed = 'yes'
    if returnPathways == 'no' or returnPathways == 'None':
        returnPathways = None
    if pipelineAnalysis == False:
        proceed = 'yes'
    if len(input_fastq_dir)>0:
        proceed = 'yes'
        run_from_scratch = 'Process RNA-seq reads'
        fl = UI.ExpressionFileLocationData('','','',''); fl.setFeatureNormalization('none')
        try: root_dir = root_dir
        except: root_dir = output_dir
        try: fl.setExpFile(expFile)
        except Exception:
            expFile = root_dir+'/ExpressionInput/exp.'+exp_name+'.txt'
            fl.setExpFile(expFile)

        fl.setArrayType(array_type)
        fl.setOutputDir(root_dir)
        fl.setMultiThreading(multiThreading)
        exp_file_location_db={}; exp_file_location_db[exp_name]=fl

        ### Assign variables needed to run Kallisto from FASTQ files
        if runKallisto and len(input_fastq_dir)==0:
            #python AltAnalyze.py --runICGS yes --platform "RNASeq" --species Mm --column_method hopach --rho 0.4 --ExpressionCutoff 1 --FoldDiff 4 --SamplesDiffering 1 --excludeCellCycle strict --output  /Users/saljh8/Desktop/Grimes/GEC14074 --expname test --fastq_dir /Users/saljh8/Desktop/Grimes/GEC14074 
            print 'Please include the flag "--fastq_dir" in the command-line arguments with an appropriate path';sys.exit()
        elif len(input_fastq_dir)>0:
            fl.setRunKallisto(input_fastq_dir)
            fl.setArrayType("3'array")
            array_type = "3'array"
            if customFASTA!=None:
                fl.setCustomFASTA(customFASTA)
    
    if proceed == 'yes':
        species_codes = UI.remoteSpeciesInfo()
        ### Update Ensembl Databases
        if 'Official' in update_method:
            file_location_defaults = UI.importDefaultFileLocations()
            db_versions_vendors,db_versions = UI.remoteOnlineDatabaseVersions()
            array_codes = UI.remoteArrayInfo()
            UI.getOnlineDBConfig(file_location_defaults,'')
            if len(species)==2:
                species_names = UI.getSpeciesInfo()
                species_full = species_names[species]
            else: species_full = species
            print 'Species name to update:',species_full
            db_version_list=[]
            for version in db_versions: db_version_list.append(version)
            db_version_list.sort(); db_version_list.reverse(); select_version = db_version_list[0]
            db_versions[select_version].sort()
            print 'Ensembl version',ensembl_version
            if ensembl_version != 'current':
                if len(ensembl_version) < 4: ensembl_version = 'EnsMart'+ensembl_version
                if ensembl_version not in db_versions:
                    try: UI.getOnlineEliteDatabase(file_location_defaults,ensembl_version,[species],'no',''); sys.exit()
                    except Exception:
                        ### This is only for database that aren't officially released yet for prototyping
                        print ensembl_version, 'is not a valid version of Ensembl, while',select_version, 'is.'; sys.exit()
                else: select_version = ensembl_version

            ### Export basic species information
            sc = species; db_version = ensembl_version
            if sc != None:
                for ad in db_versions_vendors[db_version]:
                    if ad.SpeciesCodes() == species_full:
                        for array_system in array_codes:
                            ac = array_codes[array_system]
                            compatible_species = ac.SpeciesCodes()
                            if ac.Manufacturer() in ad.Manufacturer() and ('expression' in ac.ArrayName() or 'RNASeq' in ac.ArrayName() or 'RNA-seq' in ac.ArrayName()):
                                if sc not in compatible_species: compatible_species.append(sc)
                            ac.setSpeciesCodes(compatible_species)
                UI.exportArrayInfo(array_codes)
                
            if species_full not in db_versions[select_version]:
                print db_versions[select_version]
                print species_full, ': This species is not available for this version %s of the Official database.' % select_version
            else:
                update_goelite_resources = 'no' ### This is handled separately below
                UI.getOnlineEliteDatabase(file_location_defaults,ensembl_version,[species],update_goelite_resources,'');
            
            ### Attempt to download additional Ontologies and GeneSets
            if additional_resources[0] != None: ### Indicates that the user requested the download of addition GO-Elite resources
                try:
                    from build_scripts import GeneSetDownloader
                    print 'Adding supplemental GeneSet and Ontology Collections'
                    if 'all' in additional_resources:
                        additionalResources = UI.importResourceList() ### Get's all additional possible resources
                    else: additionalResources = additional_resources
                    GeneSetDownloader.buildAccessoryPathwayDatabases([species],additionalResources,'yes')
                    print 'Finished adding additional analysis resources.'
                except Exception:
                    print 'Download error encountered for additional Ontologies and GeneSets...\nplease try again later.'
                status = UI.verifyLineageProfilerDatabases(species,'command-line')
                if status == False:
                    print 'Please note: LineageProfiler not currently supported for this species...'
                
            if array_type == 'junction' or array_type == 'RNASeq': ### Download junction databases
                try: UI.checkForLocalArraySupport(species,array_type,specific_array_type,'command-line')
                except Exception:
                    print 'Please install a valid gene database before proceeding.\n'
                    print 'For example: python AltAnalyze.py --species Hs --update Official --version EnsMart72';sys.exit()
                status = UI.verifyLineageProfilerDatabases(species,'command-line')
                print "Finished adding database"
            sys.exit()
        try:
            #print ge_fold_cutoffs,ge_pvalue_cutoffs, change_threshold, resources_to_analyze, goelite_permutations, p_val_threshold, z_threshold                    
            change_threshold = int(change_threshold)-1
            goelite_permutations = int(goelite_permutations);change_threshold = change_threshold
            p_val_threshold = float(p_val_threshold); z_threshold = float(z_threshold)
            if ORA_algorithm == 'Fisher Exact Test':
                goelite_permutations = 'FisherExactTest'
        except Exception,e:
            print e
            print 'One of the GO-Elite input values is inapporpriate. Please review and correct.';sys.exit()
        if run_GOElite == None or run_GOElite == 'no': goelite_permutations = 'NA' ### This haults GO-Elite from running
        else:
            if output_dir == None:
                print "\nPlease specify an output directory using the flag --output"; sys.exit()
        try: expression_threshold = float(expression_threshold)
        except Exception: expression_threshold = 1
        try: dabg_p = float(dabg_p)
        except Exception: dabg_p = 1 ### Occurs for RNASeq
        if microRNA_prediction_method == 'two or more': microRNA_prediction_method = 'multiple'
        else: microRNA_prediction_method = 'any'

        ### Run GO-Elite directly from user supplied input and denominator ID folders (outside of the normal workflows)
        if run_GOElite == 'yes' and pipelineAnalysis == False and '--runGOElite' in arguments:# and denom_file_dir != None:
            #python AltAnalyze.py --input "/Users/nsalomonis/Desktop/Mm_sample/input_list_small" --runGOElite yes --denom "/Users/nsalomonis/Desktop/Mm_sample/denominator" --mod Ensembl --species Mm
            """if denom_file_dir == None:
                print 'Please include a folder containing a valid denominator ID list for the input ID sets.'; sys.exit()"""
            try:
                if output_dir==None:
                    ### Set output to the same directory or parent if none selected
                    i = -1 ### 1 directory up
                    output_dir = string.join(string.split(input_file_dir,'/')[:i],'/')
                file_dirs = input_file_dir, denom_file_dir, output_dir
                import GO_Elite
                if ORA_algorithm == 'Fisher Exact Test':
                    goelite_permutations = 'FisherExactTest'
                goelite_var = species,mod,goelite_permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,returnPathways,file_dirs,''
                GO_Elite.remoteAnalysis(goelite_var,'non-UI',Multi=mlp)
                sys.exit()
            except Exception:
                print traceback.format_exc() 
                print "Unexpected error encountered. Please see log file."; sys.exit()
       
        if run_lineage_profiler == 'yes':
            status = UI.verifyLineageProfilerDatabases(species,'command-line')
            if status == False:
                print 'Please note: LineageProfiler not currently supported for this species...'

        if run_lineage_profiler == 'yes' and input_file_dir != None and pipelineAnalysis == False and ('--runLineageProfiler' in arguments or '--cellHarmony' in arguments or '--cellHarmonyMerge' in arguments):
            #python AltAnalyze.py --input "/Users/arrays/test.txt" --runLineageProfiler yes --vendor Affymetrix --platform "3'array" --species Mm --output "/Users/nsalomonis/Merrill"
            #python AltAnalyze.py --input "/Users/qPCR/samples.txt" --runLineageProfiler yes --geneModel "/Users/qPCR/models.txt" --reference "Users/qPCR/reference_profiles.txt"

            if array_type==None:
                print "Please include a platform name (e.g., --platform RNASeq)";sys.exit()
            if species==None:
                print "Please include a species name (e.g., --species Hs)";sys.exit()
            try:
                status = UI.verifyLineageProfilerDatabases(species,'command-line')
            except ValueError:
                ### Occurs due to if int(gene_database[-2:]) < 65: - ValueError: invalid literal for int() with base 10: ''
                    print '\nPlease install a valid gene database before proceeding.\n'
                    print 'For example: python AltAnalyze.py --species Hs --update Official --version EnsMart72\n';sys.exit()
            if status == False:
                print 'Please note: LineageProfiler not currently supported for this species...';sys.exit()
                
            try:
                FoldDiff=1.5
                performDiffExp=True
                pval = 0.05
                adjp = False
                for opt, arg in options: ### Accept user input for these hierarchical clustering variables
                    if opt == '--fold': FoldDiff=float(arg)
                    elif opt == '--pval': pval = float(arg)
                    elif opt == '--adjp': adjp = arg
                    elif opt == '--performDiffExp': performDiffExp = arg
                    elif opt == '--centerMethod': CenterMethod = arg
                fl = UI.ExpressionFileLocationData('','','','')
                fl.setSpecies(species)
                fl.setVendor(manufacturer)
                fl.setPlatformType(array_type)
                fl.setCompendiumType('protein_coding')
                if '--cellHarmony' in arguments:
                    fl.setClassificationAnalysis('cellHarmony')
                    fl.setPearsonThreshold(PearsonThreshold)
                    fl.setReturnCentroids(returnCentroids)
                    fl.setPeformDiffExpAnalysis(performDiffExp)
                    fl.setUseAdjPvalue(adjp)
                    fl.setPvalThreshold(pval)
                    fl.setFoldCutoff(FoldDiff)
                else:
                    fl.setClassificationAnalysis('LineageProfiler')
                #fl.setCompendiumType('AltExon')
                fl.setCompendiumPlatform(array_type)
                try: expr_input_dir
                except Exception: expr_input_dir = input_file_dir
                if '--cellHarmonyMerge' in arguments:
                    ICGS_files=[]
                    for opt, arg in options: ### Accept user input for these hierarchical clustering variables
                        if opt == '--input' or opt == '--i':
                            input_file = verifyPath(arg)
                            ICGS_files.append(input_file)
                    import LineageProfilerIterate
                    print 'center method =',CenterMethod
                    try: LineageProfilerIterate.createMetaICGSResults(ICGS_files,output_dir,CenterMethod =CenterMethod,species=species,PearsonThreshold=PearsonThreshold)
                    except:
                        LineageProfilerIterate.createMetaICGSResults(ICGS_files,output_dir,CenterMethod=CenterMethod)
                    sys.exit()
                try: CenterMethod=CenterMethod
                except: CenterMethod='centroid'
                UI.remoteLP(fl, expr_input_dir, manufacturer, custom_reference, geneModel, None, modelSize=modelSize, CenterMethod=CenterMethod) #,display=display
                #graphic_links = ExpressionBuilder.remoteLineageProfiler(fl,input_file_dir,array_type,species,manufacturer)
                print_out = 'Lineage profiles and images saved to the folder "DataPlots" in the input file folder.'
                print print_out
            except Exception:
                print traceback.format_exc() 
                print_out = 'Analysis error occured...\nplease see warning printouts.'
                print print_out
            sys.exit()
  
        if array_type == 'junction' or array_type == 'RNASeq': ### Download junction databases
            try: UI.checkForLocalArraySupport(species,array_type,specific_array_type,'command-line')
            except Exception:
                print 'Please install a valid gene database before proceeding.\n'
                print 'For example: python AltAnalyze.py --species Hs --update Official --version EnsMart72';sys.exit()
                    
        probeset_types = ['full','core','extended']
        if return_all == 'yes': ### Perform no alternative exon filtering when annotating existing FIRMA or MADS results
            dabg_p = 1; expression_threshold = 1; p_threshold = 1; alt_exon_fold_variable = 1
            gene_expression_cutoff = 10000; filter_probeset_types = 'full'; exon_exp_threshold = 1; rpkm_threshold = 0
            gene_exp_threshold = 1; exon_rpkm_threshold = 0
            if array_type == 'RNASeq':
                gene_exp_threshold = 0
        else:
            if array_type != "3'array":
                try:
                    p_threshold = float(p_threshold); alt_exon_fold_variable = float(alt_exon_fold_variable)
                    expression_threshold = float(expression_threshold); gene_expression_cutoff = float(gene_expression_cutoff)
                    dabg_p = float(dabg_p); additional_score = float(additional_score)
                    gene_expression_cutoff = float(gene_expression_cutoff)
                except Exception:
                    try: gene_expression_cutoff = float(gene_expression_cutoff)
                    except Exception: gene_expression_cutoff = 0
                try: rpkm_threshold = float(rpkm_threshold)
                except Exception: rpkm_threshold = -1
                try: exon_exp_threshold = float(exon_exp_threshold)
                except Exception: exon_exp_threshold = 0
                try: gene_exp_threshold = float(gene_exp_threshold)
                except Exception: gene_exp_threshold = 0
                try: exon_rpkm_threshold = float(exon_rpkm_threshold)
                except Exception: exon_rpkm_threshold = 0
                
                if filter_probeset_types not in probeset_types and array_type == 'exon':
                    print "Invalid probeset-type entered:",filter_probeset_types,'. Must be "full", "extended" or "core"'; sys.exit()
                elif array_type == 'gene' and filter_probeset_types == 'NA': filter_probeset_types = 'core'
                if dabg_p > 1 or dabg_p <= 0:
                    print "Invalid DABG p-value entered:",dabg_p,'. Must be > 0 and <= 1'; sys.exit()
                if expression_threshold <1:
                    print "Invalid expression threshold entered:",expression_threshold,'. Must be > 1'; sys.exit()
                if p_threshold > 1 or p_threshold <= 0:
                    print "Invalid alternative exon p-value entered:",p_threshold,'. Must be > 0 and <= 1'; sys.exit()
                if alt_exon_fold_variable < 1 and analysis_method != 'ASPIRE' and analysis_method != 'MultiPath-PSI':
                    print "Invalid alternative exon threshold entered:",alt_exon_fold_variable,'. Must be > 1'; sys.exit()
                if gene_expression_cutoff < 1:
                    print "Invalid gene expression threshold entered:",gene_expression_cutoff,'. Must be > 1'; sys.exit()
                if additional_score < 1:
                    print "Invalid additional score threshold entered:",additional_score,'. Must be > 1'; sys.exit()
                if array_type == 'RNASeq':
                    if rpkm_threshold < 0:
                        print "Invalid gene RPKM threshold entered:",rpkm_threshold,'. Must be >= 0'; sys.exit()
                    if exon_exp_threshold < 1:
                        print "Invalid exon expression threshold entered:",exon_exp_threshold,'. Must be > 1'; sys.exit()
                    if exon_rpkm_threshold < 0:
                        print "Invalid exon RPKM threshold entered:",exon_rpkm_threshold,'. Must be >= 0'; sys.exit()
                    if gene_exp_threshold < 1:
                        print "Invalid gene expression threshold entered:",gene_exp_threshold,'. Must be > 1'; sys.exit()
                        
        if 'FIRMA' in additional_algorithms and array_type == 'RNASeq':
            print 'FIRMA is not an available option for RNASeq... Changing this to splicing-index.'
            additional_algorithms = 'splicing-index'
            
        additional_algorithms = UI.AdditionalAlgorithms(additional_algorithms); additional_algorithms.setScore(additional_score)
        if array_type == 'RNASeq':
            try:
                if 'CEL' in run_from_scratch: run_from_scratch = 'Process RNA-seq reads'
                if build_exon_bedfile == 'yes': run_from_scratch = 'buildExonExportFiles'
                manufacturer = 'RNASeq'
            except Exception:
                ### When technically 3'array format
                array_type = "3'array"

        if run_from_scratch == 'Process AltAnalyze filtered': expression_data_format = 'log' ### This is switched to log no matter what, after initial import and analysis of CEL or BED files
                  
        ### These variables are modified from the defaults in the module UI as below
        excludeNonExpExons = True
        if avg_all_for_ss == 'yes': avg_all_for_ss = 'yes'
        elif 'all exon aligning' in avg_all_for_ss or 'known exons' in avg_all_for_ss or 'expressed exons' in avg_all_for_ss:
            if 'known exons' in avg_all_for_ss and array_type == 'RNASeq': excludeNonExpExons = False
            avg_all_for_ss = 'yes'
        else: avg_all_for_ss = 'no'
        
        if run_MiDAS == 'NA': run_MiDAS = 'no'
        
        if perform_alt_analysis == 'yes': perform_alt_analysis = 'yes'
        elif perform_alt_analysis == 'expression': perform_alt_analysis = 'expression'
        elif perform_alt_analysis == 'just expression': perform_alt_analysis = 'expression'
        elif perform_alt_analysis == 'no': perform_alt_analysis = 'expression'
        elif platform != "3'array": perform_alt_analysis = 'both'
        if systemToUse != None: array_type = systemToUse
        
        try: permute_p_threshold = float(permute_p_threshold)
        except Exception: permute_p_threshold = permute_p_threshold

        ### Store variables for AltAnalyzeMain
        expr_var = species,array_type,manufacturer,constitutive_source,dabg_p,expression_threshold,avg_all_for_ss,expression_data_format,include_raw_data,run_from_scratch,perform_alt_analysis
        alt_var = analysis_method,p_threshold,filter_probeset_types,alt_exon_fold_variable,gene_expression_cutoff,remove_intronic_junctions,permute_p_threshold,perform_permutation_analysis, export_NI_values, analyze_all_conditions
        additional_var = calculate_normIntensity_p, run_MiDAS, use_direct_domain_alignments_only, microRNA_prediction_method, filter_for_AS, additional_algorithms
        goelite_var = ge_fold_cutoffs,ge_pvalue_cutoffs,ge_ptype,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,goelite_permutations,mod,returnPathways
        
        if run_from_scratch == 'buildExonExportFiles':
            fl = UI.ExpressionFileLocationData('','','',''); fl.setExonBedBuildStatus('yes'); fl.setFeatureNormalization('none')
            fl.setCELFileDir(cel_file_dir); fl.setArrayType(array_type); fl.setOutputDir(output_dir)
            fl.setMultiThreading(multiThreading)
            exp_file_location_db={}; exp_file_location_db[dataset_name]=fl; parent_dir = output_dir
            perform_alt_analysis = 'expression'

        if run_from_scratch == 'Process Expression file':
            if len(input_exp_file)>0:
                if groups_file != None and comps_file != None:
                    if 'exp.' in input_exp_file: new_exp_file = input_exp_file
                    else:
                        new_exp_file = export.findParentDir(input_exp_file)+'exp.'+export.findFilename(input_exp_file)
                    if 'ExpressionInput' not in new_exp_file:
                        ### This expression file is not currently used (could make it the default after copying to this location)
                        if output_dir[-1] != '/' and output_dir[-1] != '\\':
                            output_dir += '/'
                        new_exp_file = output_dir+'ExpressionInput/'+export.findFilename(new_exp_file)
                        try: export.copyFile(input_exp_file, new_exp_file)
                        except Exception: print 'Expression file already present in target location.'
                    try: export.copyFile(groups_file, string.replace(new_exp_file,'exp.','groups.'))
                    except Exception: print 'Groups file already present in target location OR bad input path.'
                    try: export.copyFile(comps_file, string.replace(new_exp_file,'exp.','comps.'))
                    except Exception: print 'Comparison file already present in target location OR bad input path.'
                    groups_file = string.replace(new_exp_file,'exp.','groups.')
                    comps_file = string.replace(new_exp_file,'exp.','comps.')
                    input_exp_file = new_exp_file
                    if verifyGroupFileFormat(groups_file) == False:
                        print "\nWarning! The format of your groups file is not correct. For details, see:\nhttp://altanalyze.readthedocs.io/en/latest/ManualGroupsCompsCreation\n"
                        sys.exit()
                try:
                    cel_files, array_linker_db = ExpressionBuilder.getArrayHeaders(input_exp_file)
                    if len(input_stats_file)>1: ###Make sure the files have the same arrays and order first
                        cel_files2, array_linker_db2 = ExpressionBuilder.getArrayHeaders(input_stats_file)
                        if cel_files2 != cel_files:
                            print "The probe set p-value file:\n"+input_stats_file+"\ndoes not have the same array order as the\nexpression file. Correct before proceeding."; sys.exit()
                except Exception: print '\nWARNING...Expression file not found: "'+input_exp_file+'"\n\n'; sys.exit()
            
            exp_name = string.replace(exp_name,'exp.',''); dataset_name = exp_name; exp_name = string.replace(exp_name,'.txt','')
            groups_name = 'ExpressionInput/groups.'+dataset_name; comps_name = 'ExpressionInput/comps.'+dataset_name
            groups_file_dir = output_dir+'/'+groups_name; comps_file_dir = output_dir+'/'+comps_name
            groups_found = verifyFile(groups_file_dir)
            comps_found = verifyFile(comps_file_dir)
            if ((groups_found != 'found' or comps_found != 'found') and analyze_all_conditions != 'all groups') or (analyze_all_conditions == 'all groups' and groups_found != 'found'):
                files_exported = UI.predictGroupsAndComps(cel_files,output_dir,exp_name)
                if files_exported == 'yes': print "AltAnalyze inferred a groups and comps file from the CEL file names."
                elif run_lineage_profiler == 'yes' and input_file_dir != None and pipelineAnalysis == False and ('--runLineageProfiler' in arguments or '--cellHarmony' in arguments): pass 
                else: print '...groups and comps files not found. Create before running AltAnalyze in command line mode.';sys.exit()
            fl = UI.ExpressionFileLocationData(input_exp_file,input_stats_file,groups_file_dir,comps_file_dir)
            dataset_name = exp_name
            if analyze_all_conditions == "all groups":
                try: array_group_list,group_db = UI.importArrayGroupsSimple(groups_file_dir,cel_files)
                except Exception:
                    print '...groups and comps files not found. Create before running AltAnalyze in command line mode.';sys.exit()
                print len(group_db), 'groups found'
                if len(group_db) == 2: analyze_all_conditions = 'pairwise'
            exp_file_location_db={}; exp_file_location_db[exp_name]=fl
            
        elif run_from_scratch == 'Process CEL files' or run_from_scratch == 'Process RNA-seq reads' or run_from_scratch == 'Process Feature Extraction files':
            if groups_file != None and comps_file != None:
                try: shutil.copyfile(groups_file, string.replace(exp_file_dir,'exp.','groups.'))
                except Exception: print 'Groups file already present in target location OR bad input path.'
                try: shutil.copyfile(comps_file, string.replace(exp_file_dir,'exp.','comps.'))
                except Exception: print 'Comparison file already present in target location OR bad input path.'
            stats_file_dir = string.replace(exp_file_dir,'exp.','stats.')
            groups_file_dir = string.replace(exp_file_dir,'exp.','groups.')
            comps_file_dir = string.replace(exp_file_dir,'exp.','comps.')
            groups_found = verifyFile(groups_file_dir)
            comps_found = verifyFile(comps_file_dir)
            if ((groups_found != 'found' or comps_found != 'found') and analyze_all_conditions != 'all groups') or (analyze_all_conditions == 'all groups' and groups_found != 'found'):
                if mappedExonAnalysis: pass
                elif len(input_fastq_dir)<0: ### Don't check for FASTQ to allow for fast expression quantification even if groups not present
                    files_exported = UI.predictGroupsAndComps(cel_files,output_dir,exp_name)
                    if files_exported == 'yes': print "AltAnalyze inferred a groups and comps file from the CEL file names."
                    #else: print '...groups and comps files not found. Create before running AltAnalyze in command line mode.';sys.exit()
            fl = UI.ExpressionFileLocationData(exp_file_dir,stats_file_dir,groups_file_dir,comps_file_dir)
            exp_file_location_db={}; exp_file_location_db[dataset_name]=fl
            parent_dir = output_dir  ### interchangable terms (parent_dir used with expression file import)
            if analyze_all_conditions == "all groups":
                array_group_list,group_db = UI.importArrayGroupsSimple(groups_file_dir,cel_files)
                UI.exportGroups(exp_file_location_db,array_group_list)
                print len(group_db), 'groups found'
                if len(group_db) == 2: analyze_all_conditions = 'pairwise'
            try: fl.setRunKallisto(input_fastq_dir)
            except Exception: pass
        elif run_from_scratch == 'Process AltAnalyze filtered':
            if '.txt' in input_filtered_dir: ### Occurs if the user tries to load a specific file
                dirs = string.split(input_filtered_dir,'/')
                input_filtered_dir = string.join(dirs[:-1],'/')
            fl = UI.ExpressionFileLocationData('','','',''); dataset_name = 'filtered-exp_dir'
            dirs = string.split(input_filtered_dir,'AltExpression'); parent_dir = dirs[0]
            exp_file_location_db={}; exp_file_location_db[dataset_name]=fl
        
        for dataset in exp_file_location_db:
            fl = exp_file_location_db[dataset_name]
            file_location_defaults = UI.importDefaultFileLocations()
            apt_location = UI.getAPTLocations(file_location_defaults,run_from_scratch,run_MiDAS)
            fl.setAPTLocation(apt_location)
            if run_from_scratch == 'Process CEL files':
                if xhyb_remove == 'yes' and (array_type == 'gene' or array_type == 'junction'): xhyb_remove = 'no' ### This is set when the user mistakenly selects exon array, initially
                fl.setInputCDFFile(input_cdf_file); fl.setCLFFile(clf_file); fl.setBGPFile(bgp_file); fl.setXHybRemoval(xhyb_remove)
                fl.setCELFileDir(cel_file_dir); fl.setArrayType(array_type_original); fl.setOutputDir(output_dir)
            elif run_from_scratch == 'Process RNA-seq reads':
                fl.setCELFileDir(cel_file_dir); fl.setOutputDir(output_dir)
            elif run_from_scratch == 'Process Feature Extraction files':
                fl.setCELFileDir(cel_file_dir); fl.setOutputDir(output_dir)
            fl = exp_file_location_db[dataset]; fl.setRootDir(parent_dir)
            
            try: apt_location = fl.APTLocation()
            except Exception: apt_location = ''
            root_dir = fl.RootDir(); fl.setExonBedBuildStatus(build_exon_bedfile)
            fl.setMarkerFinder(marker_finder)
            fl.setFeatureNormalization(normalize_feature_exp)
            fl.setNormMatrix(normalize_gene_data)
            fl.setProbabilityStatistic(probability_statistic)
            fl.setProducePlots(visualize_qc_results)
            fl.setPerformLineageProfiler(run_lineage_profiler)
            fl.setCompendiumType(compendiumType)
            fl.setCompendiumPlatform(compendiumPlatform)
            fl.setVendor(manufacturer)
            try: fl.setFDRStatistic(FDR_statistic)
            except Exception: pass
            fl.setAnalysisMode('commandline')
            fl.setBatchEffectRemoval(batch_effects)
            fl.setChannelToExtract(channel_to_extract)
            fl.setMultiThreading(multiThreading)
            try: fl.setExcludeLowExpressionExons(excludeNonExpExons)
            except Exception: fl.setExcludeLowExpressionExons(True)
            if 'other' in manufacturer or 'Other' in manufacturer:
                ### For data without a primary array ID key
                manufacturer = "other:3'array"
                fl.setVendor(manufacturer)
        if array_type == 'RNASeq': ### Post version 2.0, add variables in fl rather than below            
            fl.setRPKMThreshold(rpkm_threshold)
            fl.setExonExpThreshold(exon_exp_threshold)
            fl.setGeneExpThreshold(gene_exp_threshold)
            fl.setExonRPKMThreshold(exon_rpkm_threshold)
            fl.setJunctionExpThreshold(expression_threshold)
            fl.setExonMapFile(exonMapFile)
            fl.setPlatformType(platformType)            
        
        ### Verify database presence
        try: dirs = unique.read_directory('/AltDatabase')
        except Exception: dirs=[]
        if species not in dirs:
            print '\n'+species,'species not yet installed. Please install before proceeding (e.g., "python AltAnalyze.py --update Official --species',species,'--version EnsMart75").'
        global commandLineMode; commandLineMode = 'yes'
        AltAnalyzeMain(expr_var, alt_var, goelite_var, additional_var, exp_file_location_db,None)
    else:
        print traceback.format_exc()
        print 'Insufficient Flags entered (requires --species and --output)'

def cleanUpCommandArguments():
    ### Needed on PC
    command_args = string.join(sys.argv,' ')
    arguments = string.split(command_args,' --')
    for argument in arguments:
        """
        argument_list = string.split(argument,' ')
        if len(argument_list)>2:
            filename = string.join(argument_list[1:],' ')
            argument = argument_list[0]+' '+string.replace(filename,' ','$$$')
        """
        argument_list = string.split(argument,' ')
        #argument = string.join(re.findall(r"\w",argument),'')
        if ':' in argument: ### Windows OS
            z = string.find(argument_list[1],':')
            if z!= -1 and z!=1: ### Hence, it is in the argument but not at the second position
                print 'Illegal parentheses found. Please re-type these and re-run.'; sys.exit()

def runCommandLineVersion():
    ### This code had to be moved to a separate function to prevent iterative runs upon AltAnalyze.py re-import
    command_args = string.join(sys.argv,' ')
    #try: cleanUpCommandArguments()
    #except Exception: null=[]
    #print 3,[sys.argv],
    if len(sys.argv[1:])>0 and '--' in command_args:
        if '--GUI' in command_args:
            ### Hard-restart of AltAnalyze while preserving the prior parameters
            command_arguments = string.split(command_args,' --')
            if len(command_arguments)>2:
                command_arguments = map(lambda x: string.split(x,' '),command_arguments)
                command_arguments = map(lambda (x,y): (x,string.replace(y,'__',' ')),command_arguments[2:])
                selected_parameters = [command_arguments[0][1]]
                user_variables={}
                for (o,v) in command_arguments: user_variables[o]=v
                AltAnalyzeSetup((selected_parameters,user_variables))
            else:
                AltAnalyzeSetup('no') ### a trick to get back to the main page of the GUI (if AltAnalyze has Tkinter conflict)
        try:
            commandLineRun()
        except Exception:
            print traceback.format_exc()

###### Determine Command Line versus GUI Control ######
command_args = string.join(sys.argv,' ')
if len(sys.argv[1:])>1 and '-' in command_args: null=[]
else:
    try:
        import Tkinter
        from Tkinter import *
        from visualization_scripts import PmwFreeze
        import tkFileDialog
        from tkFont import Font
        use_Tkinter = 'yes'
    except ImportError: use_Tkinter = 'yes'; print "\nPmw or Tkinter not found... Tkinter print out not available";

def testResultsPanel():
    from visualization_scripts import QC
    file = "/Users/nsalomonis/Desktop/code/AltAnalyze/datasets/3'Array/Merrill/ExpressionInput/exp.test.txt"
    #QC.outputArrayQC(file)
    global root; root = Tk()
    global pathway_permutations; pathway_permutations = 'NA'
    global log_file; log_file = 'null.txt'
    global array_type; global explicit_data_type
    global run_GOElite; run_GOElite = 'run-immediately'
    
    explicit_data_type = 'exon-only'
    array_type = 'RNASeq'
        
    fl = UI.ExpressionFileLocationData('','','','')
    graphic_links = []
    graphic_links.append(['PCA','PCA.png'])
    graphic_links.append(['HC','HC.png'])
    graphic_links.append(['PCA1','PCA.png'])
    graphic_links.append(['HC1','HC.png'])
    graphic_links.append(['PCA2','PCA.png'])
    graphic_links.append(['HC2','HC.png'])
    graphic_links.append(['PCA3','PCA.png'])
    graphic_links.append(['HC3','HC.png'])
    graphic_links.append(['PCA4','PCA.png'])
    graphic_links.append(['HC4','HC.png'])    
    summary_db={}
    summary_db['QC'] = graphic_links
    #summary_db={}
    fl.setGraphicLinks(graphic_links)
    
    summary_db['gene_assayed'] = 1
    summary_db['denominator_exp_genes'] = 1
    summary_db['alt_events'] = 1
    summary_db['denominator_exp_events'] = 1
    summary_db['alt_events'] = 1
    summary_db['denominator_exp_events'] = 1
    summary_db['alt_events'] = 1
    summary_db['denominator_exp_events'] = 1
    summary_db['alt_genes'] = 1
    summary_db['direct_domain_genes'] = 1
    summary_db['miRNA_gene_hits'] = 1
    #summary_db={}
    print_out = 'Analysis complete. AltAnalyze results\nexported to "AltResults/AlternativeOutput".'
    dataset = 'test'; results_dir=''
    print "Analysis Complete\n";
    if root !='' and root !=None:
        UI.InfoWindow(print_out,'Analysis Completed!')
        tl = Toplevel(); SummaryResultsWindow(tl,'GE',results_dir,dataset,'parent',summary_db)
    #sys.exit()
    
class Logger(object):
    def __init__(self,null): 
        self.terminal = sys.stdout
        self.log = open(log_file, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log = open(log_file, "a")
        self.log.write(message)
        self.log.close()

    def flush(self): pass
    
def verifyPath(filename):
    ### See if the file is in the current working directory
    new_filename = filename
    cwd = os.getcwd()
    try:
        files = unique.read_directory(cwd)
        if filename in files:
            new_filename = cwd+'/'+new_filename
    except Exception:
        pass
    try:
         ### For local AltAnalyze directories
        if os.path.isfile(cwd+'/'+filename):
            new_filename = cwd+'/'+filename
        else:
            files = unique.read_directory(cwd+'/'+filename)
            new_filename = cwd+'/'+filename
    except Exception:
            #print traceback.format_exc()
            pass

    return new_filename
    
def dependencyCheck():
    ### Make sure core dependencies for AltAnalyze are met and if not report back
    
    from pkgutil import iter_modules
    modules = set(x[1] for x in iter_modules())  ### all installed modules
    dependent_modules = ['string','csv','base64','getpass','requests']
    dependent_modules += ['warnings','sklearn','os','webbrowser']
    dependent_modules += ['scipy','numpy','matplotlib','igraph','pandas','patsy']
    dependent_modules += ['ImageTk','PIL','cairo','wx','fastcluster','pysam', 'Tkinter']
    dependent_modules += ['networkx','numba','umap','nimfa','lxml','annoy','llvmlite']
    print ''
    count=0
    for module in dependent_modules:
        if module not in modules:
            if 'ImageTk' != module and 'PIL' != module:
                print 'AltAnalyze depedency not met for:',module
            if 'fastcluster' == module:
                print '...Faster hierarchical cluster not supported without fastcluster'
            if 'pysam' == module:
                print '...BAM file access not supported without pysam'
            if 'scipy' == module:
                print '...Many required statistical routines not supported without scipy'
            if 'numpy' == module:
                print '...Many required statistical routines not supported without numpy'
            if 'matplotlib' == module:
                print '...Core graphical outputs not supported without matplotlib'
            if 'requests' == module:
                print '...Wikipathways visualization not supported without requests'
            if 'lxml' == module:
                print '...Wikipathways visualization not supported without lxml'
            if 'wx' == module:
                print '...The AltAnalyze Results Viewer requires wx'
            if 'ImageTk' == module or 'PIL' == module:
                if 'PIL' not in dependent_modules:
                    print 'AltAnalyze depedency not met for:',module
                    print '...Some graphical results displays require ImageTk and PIL'
            if 'Tkinter' == module:
                print '...AltAnalyze graphical user interface mode requires Tkinter'
            if 'igraph' == module or 'cairo' == module:
                print '...Network visualization requires igraph and cairo'
            if 'sklearn' == module:
                print '...t-SNE analysis requires sklearn'
            if 'pandas' == module or 'patsy' == module:
                print '...Combat batch effects correction requires pandas and patsy'
            count+=1
    if count>1:
        print '\nWARNING!!!! Some dependencies are not currently met.'
        print "This may impact AltAnalyze's performance\n"

def unpackConfigFiles():
    """ When pypi installs AltAnalyze in site-packages, a zip file for the Config
    and AltDatabase in the pypi installed AltAnalyze library directory. To allow
    for different flexible database versions to be written, AltDatabase and Config
    are written to the user home directory in the folder 'altanalyze'."""
    
    fn = filepath('Config/options.txt') ### See if a Config folder is already available
    fileExists = os.path.isfile(fn)
    if fileExists == False:
        import subprocess
        import shutil
        import site
        import zipfile
        from os.path import expanduser
        userdir = expanduser("~")
        
        try:
            os.mkdir(userdir+"/altanalyze")
        except:
            pass
        
        config_filepath = filepath("Config.zip")
        altdatabase_filepath = string.replace(config_filepath,'Config.zip','AltDatabase.zip')
        print '...Creating Config directory in:',userdir+"/altanalyze",
        with zipfile.ZipFile(config_filepath,"r") as zip_ref:
            zip_ref.extractall(userdir+"/altanalyze")
        with zipfile.ZipFile(altdatabase_filepath,"r") as zip_ref:
            zip_ref.extractall(userdir+"/altanalyze")
        print '...written'

if __name__ == '__main__':
    try: mlp.freeze_support()
    except Exception: pass

    try: unpackConfigFiles()
    except: pass
    #testResultsPanel()
    skip_intro = 'yes'; #sys.exit()
    #skip_intro = 'remoteViewer'
    runCommandLineVersion()
    
    dependencyCheck()
    if use_Tkinter == 'yes': AltAnalyzeSetup(skip_intro)

    """ To do list:
    1) RNA-Seq and LineageProfiler: threshold based RPKM expression filtering for binary absent present gene and exon calls
    3) SQLite for gene-set databases prior to clustering and network visualization
    5) (explored - not good) Optional algorithm type of PCA
    7) (partially) Integrate splicing factor enrichment analysis (separate module?)
    11) Update fields in summary combined alt.exon files (key by probeset)
    12) Check field names for junction, exon, RNA-Seq in summary alt.exon report
    14) Proper FDR p-value for alt.exon analyses (include all computed p-values)
    15) Add all major clustering and LineageProfiler options to UI along with stats filtering by default
    17) Support R check (and response that they need it) along with GUI gcrma, agilent array, hopach, combat
    18) Probe-level annotations from Ensembl (partial code in place) and probe-level RMA in R (or possibly APT) - google pgf for U133 array
    19) Update the software from the software
    
    Advantages of this tool kit:
    0) Easiest to use, hands down
    1) Established and novel functionality for transcriptome/proteomics analysis built in
    2) Independent and cooperative options for RNA-Seq and array analysis (splicing and gene expression)
    3) Superior functional analyses (TF-target, splicing-factor target, lineage markers, WikiPathway visualization)
    4) Options for different levels of users with different integration options (multiple statistical method options, option R support)
    5) Built in secondary analysis options for already processed data (graphing, clustering, biomarker discovery, pathway analysis, network visualization)
    6) Incorporates highly validated alternative exon identification methods, independent and jointly
    
    Primary Engineer Work:
    0) C-library calls and/or multithreading where applicable to improve peformance.
    1) MySQL or equivalent transition for all large database queries (e.g., HuEx 2.1 on-the-fly coordinate mapping).
    3) Isoform-domain network visualization and WP overlays.
    4) Webservice calls to in silico protein translation, domain prediction, splicing factor regulation.
    
    ### 2.0.9
    moncole integration
    generic and cell classification machine learning
    PCR primer design (gene centric after file selection)
    BAM->BED (local SAMTools)
    updated APT
    """