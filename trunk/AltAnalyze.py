import math
import statistics
import sys, string
import os.path
import unique
import UI
import export
import ExpressionBuilder
import ExonAnalyze_module
import ExonAnnotate_module
import ResultsExport_module
import FeatureAlignment
import time

dirfile = unique
py2app_adj = '/AltAnalyze.app/Contents/Resources/Python/site-packages.zip'

def filepath(filename):
    dir=os.path.dirname(dirfile.__file__)       #directory file is input as a variable under the main            
    fn=os.path.join(dir,filename)
    fn = string.replace(fn,py2app_adj,'')
    fn = string.replace(fn,'\\library.zip','') ###py2exe on some systems, searches for all files in the library file, eroneously
    return fn

def read_directory(sub_dir):
    dirfile = unique
    dir=os.path.dirname(dirfile.__file__)
    dir = string.replace(dir,py2app_adj,'')
    dir = string.replace(dir,'\\library.zip','')
    dir_list = os.listdir(dir + sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
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

def importGeneric(filename):
    fn=filepath(filename); key_db = {}; x = 0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x = 1
        else:
            t = string.split(data,'\t')
            key_db[t[0]] = t[1:]
    return key_db

########### Parse Input Annotations ###########
def parse_affymetrix_annotations(filename):
    temp_affy_db = {}; x=0
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():             
        probeset_data,null = string.split(line,'\n')  #remove endline
        affy_data = string.split(probeset_data[1:-1],'","')  #remove endline
        if x==0:
            x +=1; affy_headers = affy_data
        else:
            x +=1
            probesets = affy_data[0]
            temp_affy_db[probesets] = affy_data[1:]
    for header in affy_headers:
        x = 0
        while x < len(affy_headers):
            if 'rocess' in affy_headers[x]: gb = x - 1
            if 'omponent' in affy_headers[x]: gc = x - 1
            if 'unction' in affy_headers[x]: gm = x - 1
            if 'athway' in affy_headers[x]: gp = x - 1
            if 'Gene Symbol' in affy_headers[x]: gs = x - 1
            if annotation_system in affy_headers[x]: eg = x - 1
            x += 1
    for probeset in temp_affy_db:
        affy_data = temp_affy_db[probeset]
        go_bio = affy_data[gb]; go_com = affy_data[gc]; go_mol = affy_data[gm]; genmapp = affy_data[gp]; gene = affy_data[eg]; symbol = affy_data[gs]; goa=''
      
        goa = merge_go_annoations(go_bio,goa); goa = merge_go_annoations(go_com,goa)
        goa = merge_go_annoations(go_mol,goa); goa = merge_go_annoations(genmapp,goa)
        gene = string.split(gene,' /// ')
        for gene_id in gene:
            if len(goa)>10 and len(gene)<4: go_annotations[gene_id] = goa

def merge_go_annoations(go_category,goa):
    dd = ' // '; td = ' /// '
    if td in go_category:
        go_split = string.split(go_category,td)
        for entry in go_split:
            try:
                a,b,null = string.split(entry,dd)
                entry = b+dd
                goa = goa + entry
            except ValueError:
                a,null = string.split(entry,dd)
                entry = a + dd
                goa = goa + entry                
    else:
        if go_category != '---':
            if dd in go_category:
                try: a,null = string.split(go_category,dd); entry = a+dd
                except ValueError: a,b,null = string.split(go_category,dd); entry = b+dd
                goa = goa + entry        
            else: goa = goa + go_category
    return goa

def exportGOannotations(annotate_db):
    output_annotation_file = 'AltDatabase/export/GO-Entrez-Linked-Annotations.txt'
    fn=filepath(output_annotation_file); data = open(fn,'w')
    title = 'Affygene'+'\t'+'EntrezGene'+'\t'+'Symbol' +'\t'+ 'Description' +'\t'+ 'RNA-processing' +'\t'+ 'GO' +'\n'
    data.write(title)
    for affygene in annotate_db:
        a = annotate_db[affygene]
        description=a.Description(); symbol=a.Symbol(); entrez_id=a.ExternalGeneID(); rna_processing=a.RNAProcessing()
        if entrez_id in go_annotations:
            goa = go_annotations[entrez_id]
            values = affygene +'\t'+ entrez_id +'\t'+ symbol +'\t'+ description +'\t'+ rna_processing +'\t'+ goa +'\n'
            data.write(values)
    data.close()
    
class GeneAnnotationData:
    def __init__(self, geneid, description, symbol, external_geneid, rna_processing_annot):
        self._geneid = geneid; self._description = description; self._symbol = symbol
        self._rna_processing_annot = rna_processing_annot; self._external_geneid = external_geneid
    def GeneID(self): return self._geneid
    def Description(self): return self._description
    def Symbol(self): return self._symbol
    def ExternalGeneID(self): return self._external_geneid
    def TranscriptClusterIDs(self):
        if self.GeneID() in gene_transcript_cluster_db:
            transcript_cluster_list = gene_transcript_cluster_db[self.GeneID()]
            return transcript_cluster_list
        else:
            try: return transcript_cluster_list
            except UnboundLocalError: return [''] ### Occurs for SplicingIndex method with AltMouse data
    def RNAProcessing(self): return self._rna_processing_annot
    def Report(self):
        output = self.GeneID() +'|'+ self.Symbol() +'|'+ self.Description()
        return output
    def __repr__(self): return self.Report()
    
def import_annotations(filename,array_type):
    fn=filepath(filename); annotate_db = {}; x = 0
    if array_type != 'exon':
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            if x == 0: x = 1
            else:
                try: affygene, description, ll_id, symbol, rna_processing_annot = string.split(data,'\t')
                except ValueError: affygene, description, ll_id, symbol = string.split(data,'\t'); splicing_annotation = ''
                if '"' in description: null,description,null = string.split(description,'"')
                y = GeneAnnotationData(affygene, description, symbol, ll_id, rna_processing_annot)
                annotate_db[affygene] = y
    else:
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            rna_processing_annot=''
            try: ensembl, symbol, description, rna_processing_annot = string.split(data,'\t')
            except ValueError: ensembl, description, symbol = string.split(data,'\t')
            y = GeneAnnotationData(ensembl, description, symbol, ensembl, rna_processing_annot)
            annotate_db[ensembl] = y
    return annotate_db

def ProbesetCalls(array_type,probeset_class,splice_event,constitutive_call,external_exonid):
    include_probeset = 'yes'
    if array_type == 'exon':
        if len(splice_event)>2 and constitutive_call == 'yes':  constitutive_call = 'no'
        if constitutive_call == 'no' and len(splice_event)<2 and len(external_exonid)<2:  ###otherwise these are interesting probesets to keep
            if filter_probesets_by != 'full':
                if filter_probesets_by == 'extended':
                    if probeset_class == 'full': include_probeset = 'no'
                elif filter_probesets_by == 'core':
                    if probeset_class != 'core': include_probeset = 'no'
    elif array_type == 'AltMouse':
        exonid = splice_event
        if filter_probesets_by == 'exon':
            if '-' in exonid or '|' in exonid: ###Therfore the probeset represents an exon-exon junction or multi-exon probeset
                include_probeset = 'no'
        if filter_probesets_by != 'exon': 
            if '|' in exonid: include_probeset = 'no'
        if constitutive_call == 'yes': include_probeset = 'yes'
    return include_probeset,constitutive_call
                    
########### Begin Analyses ###########
class SplicingAnnotationData:
    def ArrayType(self):
        self._array_type = array_type
        return self._array_type
    def Probeset(self): return self._probeset
    def ExonID(self): return self._exonid
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
    def Chromosome(self): return self._chromosome
    def Strand(self): return self._stand
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
    def SecondaryGeneID(self): return self._secondary_geneid
    def ExonRegionID(self): return self._exon_region
    def SplicingEvent(self):
        splice_event = self._splicing_event
        if len(splice_event)!=0:
            if splice_event[0] == '|': splice_event = splice_event[1:]
        return splice_event
    def SpliceJunctions(self): return self._splice_junctions
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
        else: self._constitutive_status = 'yes'

class AffyExonSTData(SplicingAnnotationData):
    def __init__(self,ensembl_gene_id,exon_id,ens_exon_ids,transcript_cluster_id, affy_class, constituitive_call_probeset, exon_region, splicing_event, splice_junctions):
        self._geneid = ensembl_gene_id; self._external_gene = ensembl_gene_id; self._exonid = exon_id
        self._constitutive_status = constituitive_call_probeset#; self._start = probeset_start; self._stop = probeset_stop
        self._external_exonids = ens_exon_ids; self._secondary_geneid = transcript_cluster_id#; self._chromosome = chromosome; self._strand = strand
        self._probest_class = affy_class
        self._exon_region=exon_region;  self._splicing_event=splicing_event; self._splice_junctions=splice_junctions
        if self._exonid[0] == 'U': self._probeset_type = 'UTR'
        elif self._exonid[0] == 'E': self._probeset_type = 'exonic'
        elif self._exonid[0] == 'I': self._probeset_type = 'intronic'

class AffyExonSTDataAbbreviated(SplicingAnnotationData):
    def __init__(self,ensembl_gene_id,exon_id):self._geneid = ensembl_gene_id; self._exonid = exon_id
    
def importSplicingAnnotationDatabase(filename,array_type,filtered_arrayids,filter_status):
    begin_time = time.time()
    probesets_included_by_new_evidence = 0
    if filter_status == 'no': global gene_transcript_cluster_db; gene_transcript_cluster_db={}; gene_transcript_cluster_db2={}; global last_exon_region_db; last_exon_region_db = {}
    else: new_exon_db={}
    fn=filepath(filename)
    last_gene = ' '; last_exon_region = ''
    constituitive_probeset_db = {}; constituitive_gene = {}
    count = 0; x = 0; constitutive_original = {}
    #if filter_status == 'yes': exon_db = {}
    if array_type != 'exon':
        for line in open(fn,'rU').xreadlines():             
            probeset_data = cleanUpLine(line)  #remove endline
            probeset,affygene,exons,transcript_num,transcripts,probe_type_call,ensembl,block_exon_ids,block_structure,comparison_info = string.split(probeset_data,'\t')
            ###note: currently exclude comparison_info since not applicable for existing analyses
            if x == 0: x = 1
            else:
                if exons[-1] == '|': exons = exons[0:-1]
                if affygene[-1] == '|': affygene = affygene[0:-1]; constituitive_gene[affygene]=[]
                if probe_type_call == 'gene': constituitive_call = 'yes' #looked through the probe annotations and the gene seems to be the most consistent constituitive feature
                else: constituitive_call = 'no'
                include_call,constituitive_call = ProbesetCalls(array_type,'',exons,constituitive_call,'')
                if include_call == 'yes':
                    probe_data = AltMouseData(affygene,exons,ensembl,block_exon_ids,block_structure,probe_type_call)  #this used to just have affygene,exon in the values (1/17/05)
                    exon_db[probeset] = probe_data
                    if filter_status == 'yes':  new_exon_db[probeset] = probe_data
                if constituitive_call == 'yes': constituitive_probeset_db[probeset] = affygene
        genes_being_analyzed = constituitive_gene
    else:
        for line in open(fn,'rU').xreadlines():             
            probeset_data = cleanUpLine(line)  #remove endline
            if x == 0: x = 1
            else:
                probeset_id, exon_id, ensembl_gene_id, transcript_cluster_id, chromosome, strand, probeset_start, probeset_stop, affy_class, constituitive_call_probeset, external_exonid, ens_const_exons, exon_region, exon_region_start, exon_region_stop, splicing_event, splice_junctions = string.split(probeset_data,'\t')
                include_call,constituitive_call = ProbesetCalls(array_type,affy_class,splicing_event,constituitive_call_probeset,external_exonid)
                if ensembl_gene_id != last_gene: new_gene = 'yes'
                else: new_gene = 'no'
                if filter_status == 'no' and new_gene == 'yes':
                    last_exon_region_db[last_gene] = last_exon_region
                last_gene = ensembl_gene_id
                if len(exon_region)>1: last_exon_region = exon_region ### some probeset not linked to an exon region
                ###Record the transcript clusters assoicated with each gene to annotate the results later on
                if constituitive_call_probeset!=constituitive_call :probesets_included_by_new_evidence +=1#; print probeset_id,[splicing_event],[constituitive_call_probeset];kill
                if include_call == 'yes' or constituitive_call == 'yes':
                    if filter_status == 'no': probe_data = AffyExonSTDataAbbreviated(ensembl_gene_id, exon_id)
                    else: probe_data = AffyExonSTData(ensembl_gene_id, exon_id, external_exonid, transcript_cluster_id, affy_class, constituitive_call, exon_region, splicing_event, splice_junctions)
                    if filter_status == 'yes':
                        try: ### saves memory
                            null = filtered_arrayids[probeset_id]
                            new_exon_db[probeset_id] = probe_data
                        except KeyError: null = []
                    else: exon_db[probeset_id] = probe_data
                    if constituitive_call == 'yes' and filter_status == 'no': ###only perform function when initially running
                        constituitive_probeset_db[probeset_id] = ensembl_gene_id; constituitive_gene[ensembl_gene_id]=[]
                        ###Only consider transcript clusters that make up the constitutive portion of the gene or that are alternatively regulated
                        try: gene_transcript_cluster_db[ensembl_gene_id].append(transcript_cluster_id)
                        except KeyError: gene_transcript_cluster_db[ensembl_gene_id] = [transcript_cluster_id]
                if constituitive_call_probeset == 'yes' and filter_status == 'no': ###only perform function when initially running
                    try: constitutive_original[ensembl_gene_id].append(probeset_id)
                    except KeyError: constitutive_original[ensembl_gene_id] = [probeset_id]
                    try: gene_transcript_cluster_db2[ensembl_gene_id].append(transcript_cluster_id)
                    except KeyError: gene_transcript_cluster_db2[ensembl_gene_id] = [transcript_cluster_id]
        ###If no constitutive probesets for a gene as a result of additional filtering (removing all probesets associated with a splice event), add these back
        original_probesets_add = 0; genes_being_analyzed = constituitive_gene
        for gene in constitutive_original:
            if gene not in constituitive_gene:
                genes_being_analyzed[gene] = [gene]
                original_probesets_add +=1
                gene_transcript_cluster_db[gene] = gene_transcript_cluster_db2[gene]
                for probeset in constitutive_original[gene]: constituitive_probeset_db[probeset] = gene
                
        gene_transcript_cluster_db = eliminate_redundant_dict_values(gene_transcript_cluster_db)
    
    print len(exon_db),'array IDs stored as instances of SplicingAnnotationData in memory'
    print len(constituitive_probeset_db),'array IDs stored as constititive'
    print probesets_included_by_new_evidence, 'array IDs were re-annotated as NOT constitutive based on mRNA evidence'
    if array_type == 'exon': print original_probesets_add, 'genes not viewed as constitutive as a result of filtering probesets based on splicing evidence, added back'
    end_time = time.time(); time_diff = int(end_time-begin_time)
    print filename,"import finished in %d seconds" % time_diff
    if filter_status == 'yes': return new_exon_db
    else: return constituitive_probeset_db,exon_db,genes_being_analyzed

def expr_analysis(filename,constituitive_probeset_db,exon_db,annotate_db,dataset_name):
    """import list of expression values for arrayids and calculates statistics"""
    global fold_dbase; global stats_dbase
    stats_dbase = {}; array_id_db={}; fold_dbase={}
    array_group_db = {}; array_group_name_db = {}
    global array_raw_group_values; array_raw_group_values = {}; global original_array_names; original_array_names=[]
    fn=filepath(filename); first_line = 1
    for line in open(fn,'rU').xreadlines():
      data = cleanUpLine(line)
      data2 = string.split(data,'\t')
      probeset = data2[0]
      if first_line == 1:
          first_line = 0 #makes this value null for the next loop of actual array data
          ###Below ocucrs if the data is raw opposed to precomputed
          if ':' in data2[1]:
              data_type = 'raw'
              array_group_list = []
              x=0 ###gives us an original index value for each entry in the group
              for entry in data2[1:]:
                  original_array_names.append(entry)
                  array_group,array_name = string.split(entry,':')
                  try:
                      array_group_db[array_group].append(x)
                      array_group_name_db[array_group].append(array_name)
                  except KeyError:
                      array_group_db[array_group] = [x]
                      array_group_name_db[array_group] = [array_name]
                      ### below only occurs with a new group addition
                      array_group_list.append(array_group) #use this to generate comparisons in the below linked function
                  x += 1
          else: data_type = 'precomputed'
      elif data_type == 'precomputed':
          ### Import pre-computed statistics
          con_avg = float(data2[1]); rawp = float(data2[2])
          adjp = float(data2[3]); folds = data2[4:]; float_fold_list = []
          for fold in folds:
              try: float_fold_list.append(float(fold))
              except ValueError: continue #empty column
          fold_dbase[probeset] = float_fold_list
          stats_dbase[probeset]= con_avg, rawp, adjp
      elif data_type == 'raw':
          ###Use the index values from above to assign each expression value to a new database
          temp_group_array = {}; array_index_list = []  ###Use this list for permutation analysis
          for group in array_group_db:
              array_index_list.append(array_group_db[group])
              for array_index in array_group_db[group]:
                  exp_val = float(data2[array_index+1])
                  ###If non-log array data
                  if expression_data_format == 'non-log': exp_val = math.log(exp_val+1,2)
                  try:
                      ###appended is the numerical expression value for each array in the group (temporary array)
                      temp_group_array[group].append(exp_val)  #add 1 since probeset is the first column
                  except KeyError:
                      temp_group_array[group] = [exp_val]
          ####store the group database within the probeset database entry
          try:
              null = exon_db[probeset] ###To conserve memory, don't store any probesets not used for downstream analyses (e.g. not linked to mRNAs)
              array_raw_group_values[probeset] = temp_group_array
              array_index_list.sort();array_id_db[probeset] = probeset
          except KeyError:
              null = []

    ###Build all putative splicing events
    global alt_junction_db; global critical_exon_db; global exon_dbase; critical_exon_db={}
    if array_type != 'exon':
        alt_junction_db,critical_exon_db,exon_dbase,exon_inclusion_db,exon_db = ExonAnnotate_module.identifyPutativeSpliceEvents(exon_db,constituitive_probeset_db,array_id_db,agglomerate_inclusion_probesets,onlyAnalyzeJunctions)
        print 'Number of Genes with Examined Splice Events:',len(alt_junction_db)

    if agglomerate_inclusion_probesets == 'yes':
        array_raw_group_values = agglomerateInclusionProbesets(array_raw_group_values,exon_inclusion_db)

    ###Check to see if we have precomputed expression data or raw to be analyzed
    x=0; y=0; array_raw_group_values2 = {}; probesets_to_delete=[] ### Record deleted probesets
    if len(array_raw_group_values)>0:
      ###array_group_list should already be unique and correctly sorted (see above)
      for group_name in array_group_list:
        if y == 0: y+=1; group1_name = group_name
        else:
            group2_name = group_name
            for probeset in array_raw_group_values:
                data_list1 = array_raw_group_values[probeset][group1_name] #nested database entry access - baseline expression
                data_list2 = array_raw_group_values[probeset][group2_name]
                if get_non_log_avg == 'no':
                    if global_addition_factor > 0:
                        data_list1 = addGlobalFudgeFactor(data_list1,'log')
                        data_list2 = addGlobalFudgeFactor(data_list2,'log')
                    avg1 = statistics.avg(data_list1)  #control average
                    avg2 = statistics.avg(data_list2)
                else:
                    data_list1 = statistics.log_fold_conversion(data_list1)
                    data_list2 = statistics.log_fold_conversion(data_list2)
                    if global_addition_factor > 0:
                        data_list1 = addGlobalFudgeFactor(data_list1,'non-log')
                        data_list2 = addGlobalFudgeFactor(data_list2,'non-log')
                        avg1 = math.log(statistics.avg(data_list1),2)  #control average
                        avg2 = math.log(statistics.avg(data_list2),2)
                    else:
                        avg1 = math.log(statistics.avg(data_list1),2)  #control average
                        avg2 = math.log(statistics.avg(data_list2),2)
                    data_list1 = convertToLog2(data_list1) #convert the non-log fudge factor added values to log and store these
                    data_list2 = convertToLog2(data_list2)
                log_fold = avg2 - avg1
                fold = statistics.log_fold_conversion(log_fold)
                try:
                    t,df,tails = statistics.ttest(data_list1,data_list2,2,3) #unpaired student ttest, calls p_value function
                    t = abs(t); df = round(df) #Excel doesn't recognize fractions in a DF
                    p = statistics.t_probability(t,df)                
                except ZeroDivisionError: p = 1
                try: fold_dbase[probeset].append(log_fold)
                except KeyError: fold_dbase[probeset] = [0]; fold_dbase[probeset].append(log_fold)
                if probeset not in stats_dbase:
                    stats_dbase[probeset]=[avg1]; stats_dbase[probeset].append(p)
                if analysis_method != 'ANOVA' and export_splice_index_values != 'yes':  ###remove probesets where the two group means are identical... possibly remove
                    ###replace entries with the two lists for later permutation analysis
                    if p == 1:
                        #print fold_dbase[probeset], stats_dbase[probeset];kill
                        del fold_dbase[probeset];del stats_dbase[probeset]
                        probesets_to_delete.append(probeset); x += 1
                    elif avg1 < expression_threshold and avg2 < expression_threshold and p > p_threshold: ### Inserted a filtering option to exclude small variance, low expreession probesets
                        del fold_dbase[probeset];del stats_dbase[probeset]
                        probesets_to_delete.append(probeset); x += 1
                        #print probeset, avg1,avg2, p
                    else: array_raw_group_values2[probeset] = data_list1,data_list2
                else: ###Non-junction analysis can handle more than 2 groups
                    if probeset not in array_raw_group_values2:
                            array_raw_group_values2[probeset] = [data_list1]
                            array_raw_group_values2[probeset].append(data_list2)
                    else: array_raw_group_values2[probeset].append(data_list2)

    #print len(fold_dbase),len(stats_dbase),len(exon_db),len(constituitive_probeset_db);kill
    array_raw_group_values = array_raw_group_values2
    print x, "Probesets excluded prior to analysis... predicted not detected"            
    #if analysis_method == 'ASPIRE' or analysis_method == 'linearregres':
    global original_avg_const_exp_db; global original_fold_dbase
    adj_fold_dbase, relative_splicing_ratio, conditions, gene_db, constituitive_gene_db, constitutive_fold_change, original_avg_const_exp_db = constituitive_exp_normalization(fold_dbase,stats_dbase,exon_db,constituitive_probeset_db,factor_out_expression_changes,only_include_constitutive_containing_genes)

    original_fold_dbase = fold_dbase
    ###Add in constutive fold change filter to assess gene expression for ASPIRE
    gene_expression_diff_db = constitutive_expression_changes(constitutive_fold_change,annotate_db)
    ###Check to see if raw data for permutation is present for expression normalization
    if len(array_raw_group_values)>0:
        global avg_const_exp_db; global permute_lists
        avg_const_exp_db  = {}; permute_lists = []; y = 0
        while conditions > y:
            avg_const_exp_db = constituitive_exp_normalization_raw(gene_db,constituitive_gene_db,array_raw_group_values,exon_db,y,avg_const_exp_db)
            y+=1
        permute_lists = statistics.permute_arrays(array_index_list) ###Provides all pairwise permuted group comparisons
        ###Export Analysis Results for external splicing analysis (e.g. MiDAS format)
        if exportTransitResultsforAnalysis == 'yes':
            ResultsExport_module.exportTransitResults(array_group_list,array_raw_group_values,array_group_name_db,avg_const_exp_db,adj_fold_dbase,exon_db,dataset_name)
            print "Finished exporting input data for MiDAS analysis"; midas_db={}
        else:
            try: midas_db = ResultsExport_module.importMidasOutput(dataset_name)
            except IOError: midas_db={}

    return conditions,relative_splicing_ratio,stats_dbase,adj_fold_dbase,dataset_name,gene_expression_diff_db,midas_db

def agglomerateInclusionProbesets(array_raw_group_values,exon_inclusion_db):
    ###Combine expression profiles for inclusion probesets that correspond to the same splice event
    for excl_probeset in exon_inclusion_db:
        inclusion_event_profiles = []
        if len(exon_inclusion_db[excl_probeset])>1:
            for incl_probeset in exon_inclusion_db[excl_probeset]:
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

def constituitive_exp_normalization(fold_db,stats_dbase,exon_db,constituitive_probeset_db,factor_out_expression_changes,only_include_constitutive_containing_genes):
    """For every expression value, normalize to the expression of the constituitive gene features for that condition,
    then store those ratios (probeset_exp/avg_constituitive_exp) and regenerate expression values relative only to the
    baseline avg_constituitive_exp, for all conditions, to normalize out gene expression changes"""
    print "\nParameters:"
    print "Factor_out_expression_changes:",factor_out_expression_changes
    print "Only_include_constitutive_containing_genes:",only_include_constitutive_containing_genes
    print "\nAdjusting probeset average intensity values to factor out condition specific expression changes for"
    print "optimal splicing descrimination"
    gene_db = {}
    constituitive_gene_db = {}
    ### organize everything by gene
    x=0
    for probeset in fold_db:
        if x == 0: conditions = len(fold_db[probeset])
        x=1
        
    for probeset in exon_db:
      affygene = exon_db[probeset].GeneID() #exon_db[probeset] = affygene,exons,ensembl,block_exon_ids,block_structure,comparison_info
      if probeset in fold_db:
        try: gene_db[affygene].append(probeset)
        except KeyError: gene_db[affygene] = [probeset]
        if probeset in constituitive_probeset_db and (only_include_constitutive_containing_genes == 'yes' or factor_out_expression_changes == 'no'):
            #the second conditional is used to exlcude constitutive data if we wish to use all probesets for
            #background normalization rather than just the designated 'gene' probesets.
            if probeset in stats_dbase:
                baseline_exp = stats_dbase[probeset][0]
                try: constituitive_gene_db[affygene].append(probeset)
                except KeyError: constituitive_gene_db[affygene] = [probeset]

    if len(constituitive_gene_db)>0:
        ###This is blank when there are no constitutive and the above condition is implemented
        gene_db2 = constituitive_gene_db
    else: gene_db2 = gene_db
    avg_const_exp_db = {}               
    for affygene in gene_db2:
        probeset_list = gene_db2[affygene]
        x = 0
        while x < conditions:
            ### average all exp values for constituitive probesets for each condition
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
    
    adj_fold_dbase={}
    relative_splicing_ratio={}
    constitutive_fold_change={}
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
    
                exp_val_non_log = statistics.log_fold_conversion(exp_val)
                expr_to_subtract_non_log = statistics.log_fold_conversion(expr_to_subtract)
                baseline_const_exp_non_log = statistics.log_fold_conversion(baseline_const_exp)
                if factor_out_expression_changes == 'yes':
                    exp_splice_valff = exp_val_non_log/expr_to_subtract_non_log
                else: #if no, then we just normalize to the baseline constitutive expression in order to keep gene expression effects (useful if you don't trust constitutive feature expression levels)
                    exp_splice_valff = exp_val_non_log/baseline_const_exp_non_log
                
                constituitive_fold_diff = expr_to_subtract_non_log/baseline_const_exp_non_log
                ###To calculate adjusted expression, we need to get the fold change in the constituitive avg (expr_to_subtract/baseline_const_exp) and divide the experimental expression
                ###By this fold change.
                ge_adj_exp_non_log = exp_val_non_log/constituitive_fold_diff #gives a GE adjusted expression
                
                try: ge_adj_exp = math.log(ge_adj_exp_non_log,2)
                except ValueError: print probeset,ge_adj_exp_non_log,constituitive_fold_diff,exp_val_non_log,exp_val,baseline_exp, probe_fold_val, dog
                adj_probe_fold_val = ge_adj_exp - baseline_exp
                ### Here we normalize probeset expression to avg-constituitive expression by dividing probe signal by avg const.prove sig (should be < 1)
                ### refered to as steady-state normalization
                if (analysis_method == 'splicing-index') or (probeset not in constituitive_probeset_db):
                    """Can't use constituitive gene features since these have no variance for pearson analysis
                    Python will approximate numbers to a small decimal point range. If the first fold value is
                    zero, often, zero will be close to but not exactly zero. Correct below """
                    try:
                        adj_fold_dbase[probeset].append(adj_probe_fold_val)
                    except KeyError:
                        if abs(adj_probe_fold_val - 0) < 0.0000001: #make zero == exactly to zero
                            adj_probe_fold_val = 0
                        adj_fold_dbase[probeset] = [adj_probe_fold_val] 
                try: relative_splicing_ratio[probeset].append(exp_splice_valff) ###ratio of junction exp relative to gene expression at that time-point              
                except KeyError: relative_splicing_ratio[probeset] = [exp_splice_valff]
                n = 0
                #if expr_to_subtract_non_log != baseline_const_exp_non_log: ###otherwise this is the first value in the expression array
                if x!=0: ###previous expression can produce errors when multiple group averages have identical values
                    fold_change = expr_to_subtract_non_log/baseline_const_exp_non_log
                    fold_change_log = math.log(fold_change,2)
                    constitutive_fold_change[affygene] = fold_change_log
                    ### If we want to remove any genes from the analysis with large transcriptional changes
                    ### that may lead to false positive splicing calls (different probeset kinetics)
                    if remove_transcriptional_regulated_genes == 'yes':
                        if abs(fold_change_log) > log_fold_cutoff:
                            del constitutive_fold_change[affygene]
                            try: del adj_fold_dbase[probeset]
                            except KeyError: n = 1
                            try: del relative_splicing_ratio[probeset]
                            except KeyError: n = 1
                """elif expr_to_subtract_non_log == baseline_const_exp_non_log:  ###This doesn't make sense, since n can't equal 1 if the conditional is false (check this code again later 11/23/07)
                    if n == 1:
                        del adj_fold_dbase[probeset]
                        del relative_splicing_ratio[probeset]"""
            x += 1
    print "Intensity normalization complete..."

    if factor_out_expression_changes == 'no':
        adj_fold_dbase = fold_db #don't change expression values
    print len(constitutive_fold_change), "Genes undergoing analysis for alternative splicing/transcription"
    global gene_analyzed; gene_analyzed = len(constituitive_gene_db)
    return adj_fold_dbase, relative_splicing_ratio, conditions, gene_db, constituitive_gene_db,constitutive_fold_change, avg_const_exp_db

def constituitive_exp_normalization_raw(gene_db,constituitive_gene_db,array_raw_group_values,exon_db,y,avg_const_exp_db):
    """normalize expression for raw expression data (only for non-baseline data)"""
    #avg_true_const_exp_db[affygene] = [avg_const_exp]
    temp_avg_const_exp_db={}
    
    for probeset in array_raw_group_values:
        try:
            conditions = len(array_raw_group_values[probeset][y]) #number of raw expresson values to normalize
            break
        except IndexError: print array_raw_group_values[probeset], y;kill

    for affygene in gene_db:
        ###This is blank when there are no constitutive or the above condition is implemented
        if affygene in constituitive_gene_db:
            probeset_list = constituitive_gene_db[affygene]  
            z = 1
        else:  ###so we can analyze splicing independent of gene expression even if no 'gene' feature is present
            probeset_list = gene_db[affygene]
            z = 0
        x = 0
        while x < conditions:
            ### average all exp values for constituitive probesets for each conditionF
            exp_list=[]
            for probeset in probeset_list:
                try: exp_val = array_raw_group_values[probeset][y][x] ### try statement is used for constitutive probes that were deleted due to filtering in expr_analysis
                except KeyError: continue
                exp_list.append(exp_val)
            try: avg_const_exp = statistics.avg(exp_list) 
            except ZeroDivisionError: avg_const_exp = 'null'
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
                    #non_log_exp_val = statistics.log_fold_conversion(exp_val)
                    #non_log_const_exp_val = statistics.log_fold_conversion(const_exp_val)
                    #non_log_exp_ratio = non_log_exp_val/non_log_const_exp_val
                    log_exp_ratio = exp_val - const_exp_val
                    try: normalized_raw_exp_ratios[probeset].append(log_exp_ratio)
                    except KeyError: normalized_raw_exp_ratios[probeset] = [log_exp_ratio]
    return avg_const_exp_db

def formatGeneSymbolHits(geneid_list):
    symbol_list=[]
    for geneid in geneid_list:
        symbol = ''
        if geneid in annotate_db: symbol = annotate_db[geneid].Symbol()
        if len(symbol)<1: symbol = geneid
        symbol_list.append(symbol)
    symbol_str = string.join(symbol_list,', ')
    return symbol_str

def calculateZScores(hit_count_db,denom_count_db,total_gene_denom_count,total_gene_hit_count,z_results_db,gene_results_db,symbol_results_db):
    N = float(total_gene_denom_count)  ###Genes examined
    R = float(total_gene_hit_count)  ###AS genes

    for element in denom_count_db:
        element_denom_gene_count = denom_count_db[element]
        if element in hit_count_db:
            element_hit_gene_count = len(hit_count_db[element])
            gene_symbols = formatGeneSymbolHits(hit_count_db[element])
            n = float(element_denom_gene_count) ###all genes associated with element
            r = float(element_hit_gene_count) ###regulated genes associated with element
            #z = statistics.zscore(r,n,N,R)
            try: z = (r - n*(R/N))/math.sqrt(n*(R/N)*(1-(R/N))*(1-((n-1)/(N-1))))
            except ValueError: print 'error:',element,r,n,N,R
            try: z_results_db[element].append(str(z))
            except KeyError: z_results_db[element] = [str(z)]
            try: gene_results_db[element].append(str(r))
            except KeyError: gene_results_db[element] = [str(r)]
            try:
                try: symbol_results_db[element].append(gene_symbols)
                except KeyError: symbol_results_db[element] = [gene_symbols]
            except TypeError:
                print gene_symbols,symbol_results_db;kill
        else:
            try: z_results_db[element].append(str(0))
            except KeyError: z_results_db[element] = [str(0)]
            try: gene_results_db[element].append(str(0))
            except KeyError: gene_results_db[element] = [str(0)]
            try: symbol_results_db[element].append('')
            except KeyError: symbol_results_db[element] = ['']
            
    return z_results_db,gene_results_db,symbol_results_db

def exportZScoreData(results_db,headers,gene_results_db,denom_count_db,symbol_results_db,element_type):
    element_output = 'AltResults/AlternativeOutput/' + dataset_name + analysis_method+'-'+element_type+'-zscores.txt'
    data = export.createExportFile(element_output,'AltResults/AlternativeOutput')
    headers_z = []; headers_gene = []; headers_symbol = []
    for header in headers: headers_z.append(header+'-Zscore')
    for header in headers: headers_gene.append(header+'-gene-count')
    for header in headers: headers_symbol.append(header+'-gene-symbols')
    title = [element_type]+headers_z+headers_gene+headers_symbol+['all-denominator-gene-count']
    title = string.join(title,'\t')+'\n'; data.write(title)
    for element in results_db:
        element_denom_gene_count = str(denom_count_db[element]) ###This is the denominator for all (not inclusion or exclusion which we don't report)
        values = [element]+results_db[element]+gene_results_db[element]+symbol_results_db[element]+[element_denom_gene_count]
        values = string.join(values,'\t')+'\n'; data.write(values)
        
def splicing_analysis_algorithms(relative_splicing_ratio,fold_dbase,dataset_name,gene_expression_diff_db,exon_db):
    protein_exon_feature_db={}
    print "Begining to run", analysis_method, "algorithm on",dataset_name[0:-1],"data"
    if analysis_method == 'ASPIRE' or analysis_method == 'linearregres':
        splice_event_list, p_value_call, permute_p_values = analyzeJunctionSplicing(relative_splicing_ratio)
    if analysis_method == 'splicing-index':
        splice_event_list, p_value_call, permute_p_values = analyzeSplicingIndex(fold_dbase)
        
    exon_hits={}
    ###Run analyses in the ExonAnalyze_module module to assess functional changes
    for (score,ed) in splice_event_list:
        geneid = ed.GeneID()
        if array_type == 'exon': uid = ed.Probeset1()
        else: uid = (ed.Probeset1(),ed.Probeset2())
        gene_exon = geneid,uid; exon_hits[gene_exon] = ed
        
    dataset_name_original = analysis_method+'-'+dataset_name[8:-1]
    global functional_attribute_db; global protein_features
    
    filtered_microRNA_exon_db, exon_attribute_db, gene_top_attribute_db, functional_attribute_db, protein_features = ExonAnalyze_module.function_exon_analysis_main(array_type,exon_hits,dataset_name_original,microRNA_full_exon_db,protein_sequence_dbase,probeset_protein_db,protein_ft_db)
    
    ###add microRNA data to functional_attribute_db
    microRNA_hit_gene_count_db = {}; all_microRNA_gene_hits={}; microRNA_attribute_db={}
    for (affygene,uid) in filtered_microRNA_exon_db: ###example ('G7091354', 'E20|') [('hsa-miR-130a', 'Pbxip1'), ('hsa-miR-130a', 'Pbxip1'
        ###3-1-08
        miR_list = []
        microRNA_symbol_list = filtered_microRNA_exon_db[(affygene,uid)]
        for mir_key in microRNA_symbol_list:
            microRNA,gene_symbol,miR_seq, miR_sources = mir_key
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
        
        miR_str = string.join(miR_list,','); miR_str = '('+miR_str+')'
        function_type = ('microRNA-target'+miR_str,'~')
        try: functional_attribute_db[(affygene,uid)].append(function_type)
        except KeyError: functional_attribute_db[(affygene,uid)]=[function_type]
        all_microRNA_gene_hits[affygene] = []
            
    ###Replace the gene list for each microRNA hit with count data            
    microRNA_hit_gene_count_db = eliminate_redundant_dict_values(microRNA_hit_gene_count_db)
    """for microRNA in microRNA_hit_gene_count_db:
        microRNA_hit_gene_count_db[microRNA] = len(microRNA_hit_gene_count_db[microRNA])"""
            
    ###Combines any additional feature alignment info identified from 'ExonAnalyze_module.function_exon_analysis_main' (e.g. from Ensembl or junction-based queries rather than exon specific) and combines
    ###this with this database of (Gene,Exon)=[(functional element 1,'~'),(functional element 2,'~')] for downstream result file annotatations
    domain_hit_gene_count_db = {}; all_domain_gene_hits = {}
    for entry in protein_features:
        gene,uid = entry
        for data_tuple in protein_features[entry]:
            domain,call = data_tuple
            try:protein_exon_feature_db[entry].append(data_tuple)
            except KeyError:protein_exon_feature_db[entry] = [data_tuple]
            try: domain_hit_gene_count_db[domain].append(gene)
            except KeyError: domain_hit_gene_count_db[domain] = [gene]
            all_domain_gene_hits[gene]=[]
    ###Replace the gene list for each microRNA hit with count data
    domain_hit_gene_count_db = eliminate_redundant_dict_values(domain_hit_gene_count_db)
    """for domain in domain_hit_gene_count_db:
        domain_hit_gene_count_db[domain] = len(domain_hit_gene_count_db[domain])"""

    ############   Perform Element Over-Representation Analysis ############
    """Domain/FT Fishers-Exact test: with "protein_exon_feature_db" (transformed to "domain_hit_gene_count_db") we can analyze over-representation of domain/features WITHOUT taking into account exon-inclusion or exclusion
    Do this using: "domain_gene_counts", which contains domain tuple ('Tyr_pkinase', 'IPR001245') as a key and count in unique genes as the value in addition to
    Number of genes linked to splice events "regulated" (SI and Midas p<0.05), number of genes with constitutive probesets

    MicroRNA Fishers-Exact test: "filtered_microRNA_exon_db" contains gene/exon to microRNA data. For each microRNA, count the representation in spliced genes microRNA (unique gene count - make this from the mentioned file)
    Do this using: "microRNA_count_db"""
    
    total_microRNA_gene_hit_count = len(all_microRNA_gene_hits)
    total_microRNA_gene_denom_count = len(gene_microRNA_denom); microRNA_z_scores={}; microRNA_gene_count={}; microRNA_z_score_headers=[]; microRNA_symbol_results = {}
    microRNA_z_scores,microRNA_gene_count,microRNA_symbol_results = calculateZScores(microRNA_hit_gene_count_db,microRNA_count_db,total_microRNA_gene_denom_count,total_microRNA_gene_hit_count,microRNA_z_scores,microRNA_gene_count,microRNA_symbol_results)
    microRNA_z_score_headers+=['all']
    
    total_domain_gene_hit_count = len(all_domain_gene_hits);domain_z_scores={}; domain_z_score_headers=[]; domain_gene_count={}; domain_symbol_results={}
    total_domain_gene_denom_count = len(protein_ft_db) ###genes connected to domain annotations
    domain_z_scores,domain_gene_count,domain_symbol_results = calculateZScores(domain_hit_gene_count_db,domain_gene_counts,total_domain_gene_denom_count,total_domain_gene_hit_count,domain_z_scores,domain_gene_count,domain_symbol_results)
    domain_z_score_headers+=['all']

    ###########  Import critical exon annotation for junctions, build through the exon array analysis pipeline - link back to probesets
    exon_db={}; filtered_arrayids={}; filter_status='yes'; global critical_exon_annotation_db; critical_probeset_annotation_db={}
    if array_type != 'exon':
        critical_exon_annotation_file = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_"+array_type+"_probesets.txt"
        for key in exon_hits:
            gene,uid = key; probeset1,probeset2 = uid
            critical_exons = critical_exon_junction_db[uid]
            for critical_exon in critical_exons:
                try: filtered_arrayids[gene+':'+critical_exon].append(uid)
                except KeyError: filtered_arrayids[gene+':'+critical_exon]=[uid]
        critical_exon_annotation_db = importSplicingAnnotationDatabase(critical_exon_annotation_file,'exon',filtered_arrayids,filter_status);null=[] ###The file is in exon centric format, so designate array_type as exon
        for key in critical_exon_annotation_db:
            ced = critical_exon_annotation_db[key]
            for junction_probesets in filtered_arrayids[key]:
                try: critical_probeset_annotation_db[junction_probesets].append(ced) ###use for splicing and Exon annotations
                except KeyError: critical_probeset_annotation_db[junction_probesets] = [ced]
        for junction_probesets in critical_probeset_annotation_db:
            if len(critical_probeset_annotation_db[junction_probesets])>1: ###Thus multiple exons associated, must combine annotations
                exon_ids=[]; external_exonids=[]; exon_regions=[]; splicing_events=[]
                for ed in critical_probeset_annotation_db[junction_probesets]:
                    ensembl_gene_id = ed.GeneID();  transcript_cluster_id = ed.ExternalGeneID(); affy_class = ed.ProbesetClass()
                    exon_ids.append(ed.ExonID()); external_exonids.append(ed.ExternalExonIDs());  exon_regions.append(ed.ExonRegionID()); se = string.split(ed.SplicingEvent(),'|')
                    for i in se: splicing_events.append(i)
                splicing_events = unique.unique(splicing_events) ###remove duplicate entries
                exon_id = string.join(exon_ids,'|'); external_exonid = string.join(external_exonids,'|'); exon_region = string.join(exon_regions,'|'); splicing_event = string.join(splicing_events,'|')
                probe_data = AffyExonSTData(ensembl_gene_id, exon_id, external_exonid, transcript_cluster_id, affy_class, '', exon_region, splicing_event, '')
                critical_probeset_annotation_db[junction_probesets] = probe_data
            else:
                critical_probeset_annotation_db[junction_probesets] = critical_probeset_annotation_db[junction_probesets][0]

                
    ###########  Re-import the exon_db for significant entries with full annotaitons
    exon_db={}; filtered_arrayids={} ###Use this as a means to save memory (import multiple times - only storing different types relevant information)
    for (score,entry) in splice_event_list:
        probeset = entry.Probeset1(); filtered_arrayids[probeset] = []
        if array_type != 'exon':
            try: probeset = entry.Probeset2(); filtered_arrayids[probeset] = []
            except AttributeError: null =[] ###occurs when running Splicing 
    exon_db= importSplicingAnnotationDatabase(probeset_annotations_file,array_type,filtered_arrayids,filter_status);null=[] ###replace existing exon_db (probeset_annotations_file should be a global)

    ############   Export exon/junction level results ############    
    splice_event_db={}; protein_length_list=[]; aspire_gene_results={}
    critical_gene_exons={}; unique_exon_event_db={}; comparison_count={}
    functional_attribute_db2={}; protein_exon_feature_db2={}; microRNA_exon_feature_db2={}
    external_exon_annot={}; gene_exon_region={}; gene_smallest_p={}; gene_splice_event_score={}; alternatively_reg_tc={}
    aspire_output = 'AltResults/AlternativeOutput/' + dataset_name + analysis_method+'-exon-inclusion-results.txt'
    data = export.createExportFile(aspire_output,'AltResults/AlternativeOutput')
    if analysis_method != 'splicing-index':
        title = ['AffyGene','dI','symbol','description','exons1','exons2','regulation_call','event_call','probeset1','rawp1','probeset2','rawp2','fold1','fold2']
        title +=['adj-fold1' ,'adj-fold2' ,'block_structure','critical_up_exons','critical_down_exons','functional_prediction','uniprot-ens_feature_predictions']
        title +=['peptide_predictions','exp1','exp2','max_expression','constitutive_baseline_exp','p_value_call','permutation-values','permutation-false-positives']
        title +=['gene-expression-change','splice_event_description','ExternalExonIDs','ExonRegionID','SplicingEvent','ExonAnnotationScore','large_splicing_diff','opposite_splicing_pattern']
    else:
        title= ['Ensembl','dI','symbol','description','exons','regulation_call','event_call','probeset','lowest_p (MIDAS or SI)','midas p-value','fold','adjfold']
        title+=['up_exons','down_exons','functional_prediction','uniprot-ens_feature_predictions','peptide_predictions','exp','highest_exp']
        title+=['constitutive_baseline_exp','SI_p-value','probeset p-value','gene-expression-change']
        title+=['transcript cluster ID', 'ensembl exons', 'consitutive probeset', 'exon-region-ID', 'exon annotations','distal exon-region-ID','exon annotation score']    
    title = string.join(title,'\t') + '\n'
    data.write(title)

    event_count = 0
    for (score,entry) in splice_event_list:
        event_count += 1
        dI = entry.Score();probeset1 = entry.Probeset1(); regulation_call = entry.RegulationCall(); event_call = entry.EventCall();critical_exon_list = entry.CriticalExonTuple(); critical_probeset_list = [probeset1]
        affygene = exon_db[probeset1].GeneID(); exons1 = exon_db[probeset1].ExonID()
        if affygene in annotate_db: description = annotate_db[affygene].Description(); symbol = annotate_db[affygene].Symbol()
        else: description = ''; symbol = ''
        try:
            adjfold1 = str(fold_dbase[probeset1][1])
            exp1 = str(stats_dbase[probeset1][0]); fold1 = str(original_fold_dbase[probeset1][1])
        except KeyError:
            print [probeset1],len(stats_dbase),len(original_fold_dbase),len(fold_dbase),'\n',stats_dbase[probeset1],'\n',fold_dbase[probeset1],'\n',original_fold_dbase[probeset1];kill
        #rawp1=''; adjfold1=''; exp1=''; fold1=''
        baseline_const_exp = original_avg_const_exp_db[affygene][0]

        if analysis_method != 'splicing-index':   
            probeset2 = entry.Probeset2(); exons2 = exon_db[probeset2].ExonID(); rawp1 = str(entry.TTestNormalizedRatios()); rawp2 = str(entry.TTestNormalizedRatios2()); critical_probeset_list.append(probeset2)
            exp2 = str(stats_dbase[probeset2][0]); fold2 = str(original_fold_dbase[probeset2][1]); adjfold2 = str(fold_dbase[probeset2][1])
            block_structure = exon_db[probeset1].GeneStructure()
            exp_list = [float(exp1),float(exp2),float(exp1)+float(fold1),float(exp2)+float(fold2)]; exp_list.sort();  exp_list.reverse()
            probeset_tuple = (probeset1,probeset2)
        else:
            exp_list = [float(exp1),float(exp1)+float(fold1)]; exp_list.sort();  exp_list.reverse()
            probeset_tuple = (probeset1)
        highest_exp = exp_list[0]
        
        ###Use permuted p-value or lowest expression junction p-value based on the situtation
        ###This p-value is used to filter out aspire events for further analyses
        if p_value_call == 'permuted_aspire_p-value':
            if probeset_tuple in permute_p_values:
                lowest_raw_p, pos_permute, total_permute, false_pos = permute_p_values[probeset_tuple]
            else: lowest_raw_p = "NA"; pos_permute = "NA"; total_permute = "NA"; false_pos = "NA"
        else:
            if analysis_method != 'splicing-index': raw_p_list = [entry.TTestNormalizedRatios(),entry.TTestNormalizedRatios2()]  #raw_p_list = [float(rawp1),float(rawp2)]; raw_p_list.sort()
            else: raw_p_list = [float(entry.TTestNormalizedRatios())]  ###Could also be rawp1, but this is more appropriate
            raw_p_list.sort()
            lowest_raw_p = raw_p_list[0]; pos_permute = "NA"; total_permute = "NA"; false_pos = "NA"
            
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
        elif variable == 2:
            exon1 = exon_data[1][0]; exon2 = exon_data[1][1]
            #print exon_db[probeset1][0],fold_dbase[probeset1][1], exon_db[probeset1][1], fold_dbase[probeset2][1],exon_db[probeset1][2]
            if fold_dbase[probeset1][1] > 0:
                up_exons = up_exons + exon1 + ',';down_exons = down_exons + exon2 + ','
                up_exon_list.append(exon1); down_exon_list.append(exon2)
                key = affygene,exon1+'|';gene_exon_list.append(key);key = affygene,exon2+'|'; gene_exon_list.append(key)
            else:
                up_exons = up_exons + exon2 + ',';down_exons = down_exons + exon1 + ','
                up_exon_list.append(exon2); down_exon_list.append(exon1)
                key = affygene,exon1+'|'; gene_exon_list.append(key); key = affygene,exon2+'|'; gene_exon_list.append(key)
        up_exons = up_exons[0:-1];down_exons = down_exons[0:-1]
        
        ###Format functional results based on exon level fold change
        null = []
        #global a; a = exon_hits; global b; b=microRNA_attribute_db; kill
        new_functional_attribute_str, functional_attribute_list2, seq_attribute_str,protein_length_list = format_exon_functional_attributes(affygene,critical_probeset_list,functional_attribute_db,up_exon_list,down_exon_list,protein_length_list)
        new_uniprot_exon_feature_str, uniprot_exon_feature_list, null, null = format_exon_functional_attributes(affygene,critical_probeset_list,protein_exon_feature_db,up_exon_list,down_exon_list,null)
        null, microRNA_exon_feature_list, null, null = format_exon_functional_attributes(affygene,critical_probeset_list,microRNA_attribute_db,up_exon_list,down_exon_list,null)
        if len(new_functional_attribute_str) == 0: new_functional_attribute_str = ' '
        if len(new_uniprot_exon_feature_str) == 0: new_uniprot_exon_feature_str = ' '

        ### Add entries to a database to quantify the number of reciprocal isoforms regulated
        reciprocal_isoform_data = [len(critical_exon_list[1]),critical_exon_list[1],event_call,regulation_call]
        if float((lowest_raw_p)<=p_threshold or false_pos < 2):
            try: unique_exon_event_db[affygene].append(reciprocal_isoform_data)
            except KeyError: unique_exon_event_db[affygene] = [reciprocal_isoform_data]
            
        ### Add functional attribute information to a new database
        for item in uniprot_exon_feature_list:
            attribute = item[0]
            exon = item[1]
            if float((lowest_raw_p)<=p_threshold or false_pos < 2):
              try: protein_exon_feature_db2[affygene,attribute].append(exon)
              except KeyError: protein_exon_feature_db2[affygene,attribute]=[exon]
        ### Add functional attribute information to a new database
        """Database not used for exon/junction data export but for over-representation analysis (direction specific)"""
        for item in microRNA_exon_feature_list:
            attribute = item[0]
            exon = item[1]
            if float((lowest_raw_p)<=p_threshold or false_pos < 2):
              try: microRNA_exon_feature_db2[affygene,attribute].append(exon)
              except KeyError: microRNA_exon_feature_db2[affygene,attribute]=[exon]              
        ### Add functional attribute information to a new database
        for item in functional_attribute_list2:
            attribute = item[0]
            exon = item[1]
            if float((lowest_raw_p)<=p_threshold or false_pos < 2):
              try: functional_attribute_db2[affygene,attribute].append(exon)
              except KeyError: functional_attribute_db2[affygene,attribute]=[exon]
                    
        if affygene in gene_expression_diff_db: mean_fold_change = str(gene_expression_diff_db[affygene][0])
        else: mean_fold_change = 'NC'; print affygene; kill #ENSG00000205542

        abs_fold = abs(float(mean_fold_change)); fold_direction = 'down'; fold1_direction = 'down'; fold2_direction = 'down'
        large_splicing_diff1 = 0; large_splicing_diff2 = 0; large_splicing_diff = 'null'; opposite_splicing_pattern = 'no'
        if float(mean_fold_change)>0: fold_direction = 'up'
        if float(fold1)>0: fold1_direction = 'up'
        if fold1_direction != fold_direction:
            if float(fold1)>float(mean_fold_change): large_splicing_diff1 = float(fold1)-float(mean_fold_change)


        if array_type == 'exon': ed = exon_db[probeset1]
        else:
            try: ed = critical_probeset_annotation_db[probeset1,probeset2]
            except KeyError: ed = exon_db[probeset1] ###not useful data here, but the objects need to exist
        ucsc_splice_annotations = ["retainedIntron","cassetteExon","strangeSplice","altFivePrime","altThreePrime","altThreePrime","altPromoter","bleedingExon"]
        custom_annotations = ["alt-3'","alt-5'","alt-C-term","alt-N-term","cassette-exon","cassette-exon","exon-region-exclusion","intron-retention"]

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
                    
        if analysis_method != 'splicing-index':
            if float(fold2)>0: fold2_direction = 'up'
            if fold2_direction != fold_direction:
                if float(fold2)>float(mean_fold_change):
                    large_splicing_diff2 = float(fold2)-float(mean_fold_change)
                    if abs(large_splicing_diff2) > large_splicing_diff1: large_splicing_diff = str(large_splicing_diff2)
                    else: large_splicing_diff = str(large_splicing_diff1)
            if fold1_direction != fold2_direction and abs(float(fold1))>0.4 and abs(float(fold2))>0.4 and abs(float(mean_fold_change))< max([float(fold2),float(fold1)]):
                opposite_splicing_pattern = 'yes'                                                                                         

            ### Annotate splicing events based on exon_strucuture data
            splice_event = ExonAnnotate_module.annotate_splice_event(exons1,exons2,block_structure)
            try: splice_event_db[splice_event] += 1
            except KeyError: splice_event_db[splice_event] = 1

            ### Annotate splicing events based on pre-computed and existing annotations

            values= [affygene,dI,symbol,description,exons1,exons2,regulation_call,event_call,probeset1,rawp1,probeset2,rawp2,fold1,fold2,adjfold1,adjfold2]
            values+=[block_structure,up_exons,down_exons,new_functional_attribute_str,new_uniprot_exon_feature_str,seq_attribute_str,exp1,exp2,str(highest_exp)]
            values+=[str(baseline_const_exp),str(lowest_raw_p),str(pos_permute)+' out of '+str(total_permute),str(false_pos),mean_fold_change,splice_event]
            values+=[ed.ExternalExonIDs(),ed.ExonRegionID(),ed.SplicingEvent(),str(exon_annot_score),large_splicing_diff,opposite_splicing_pattern]
            exon_sets = abs(float(dI)),regulation_call,event_call,exons1,exons2,''            
        else:
            si_pvalue = lowest_raw_p; rawp1 = str(stats_dbase[probeset1][1])
            if probeset1 in midas_db:
                midas_p = str(midas_db[probeset1])
                if float(midas_p)<lowest_raw_p: lowest_raw_p = float(midas_p) ###This is the lowest and SI-pvalue
            else: midas_p = ''
            try:
                es = probeset_protein_db[probeset1]
                hit_exons = es.HitProteinExonIDs(); null_exons = es.NullProteinExonIDs()
            except KeyError: hit_exons =''; null_exons=''
            
            ###Determine what type of exon-annotations are present to assign a confidence score
            if affygene in annotate_db: ###Determine the transcript clusters used to comprise a splice event (genes and exon specific)
                gene_tc = annotate_db[affygene].TranscriptClusterIDs()
                probeset_tc = [ed.SecondaryGeneID()]
                for transcript_cluster in gene_tc: probeset_tc.append(transcript_cluster)
                probeset_tc = makeUnique(probeset_tc)
            cluster_number = len(probeset_tc)
            try: alternatively_reg_tc[affygene] += probeset_tc
            except KeyError: alternatively_reg_tc[affygene] = probeset_tc

            try: last_exon_region = last_exon_region_db[affygene]
            except KeyError: last_exon_region = ''
            if cluster_number>1: exon_annot_score = 1            
            
            values= [affygene,dI,symbol,description,exons1,regulation_call,event_call,probeset1,str(lowest_raw_p),midas_p,fold1,adjfold1]
            values+=[up_exons,down_exons,new_functional_attribute_str,new_uniprot_exon_feature_str,seq_attribute_str,exp1,str(highest_exp)]
            values+=[str(baseline_const_exp),str(si_pvalue),rawp1,mean_fold_change,ed.SecondaryGeneID(), ed.ExternalExonIDs()]
            values+=[ed.Constitutive(),ed.ExonRegionID(),ed.SplicingEvent(),last_exon_region,str(exon_annot_score)]

            exon_sets = abs(float(dI)),regulation_call,event_call,exons1,exons1,midas_p
            
        if len(ed.SplicingEvent())>2:
            try: external_exon_annot[affygene].append(ed.SplicingEvent())
            except KeyError: external_exon_annot[affygene] = [ed.SplicingEvent()]
            
        values = string.join(values,'\t')+'\n'
        data.write(values)
        ###Process data for gene level reports
        if float((lowest_raw_p)<=p_threshold or false_pos < 2):
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

    ### functional_attribute_db2 will be reorganized so save the database with another. Use this  
    functional_attribute_db = functional_attribute_db2
    functional_attribute_db2 = reorganize_attribute_entries(functional_attribute_db2,'no')
    external_exon_annot = eliminate_redundant_dict_values(external_exon_annot)
    protein_exon_feature_db = protein_exon_feature_db2
    microRNA_exon_feature_db = microRNA_exon_feature_db2
    protein_exon_feature_db2,incl_domain_hit_count,all_incl_domain_gene_hits,excl_domain_hit_count,all_excl_domain_gene_hits = reorganize_attribute_entries(protein_exon_feature_db2,'yes')
    microRNA_exon_feature_db2,incl_microRNA_hit_count,all_incl_microRNA_gene_hits,excl_microRNA_hit_count,all_excl_microRNA_gene_hits = reorganize_attribute_entries(microRNA_exon_feature_db,'yes')
    
    ############   Perform Element Over-Representation Analysis ############
    total_microRNA_gene_hit_count = len(all_incl_microRNA_gene_hits)
    microRNA_z_scores,microRNA_gene_count,microRNA_symbol_results = calculateZScores(incl_microRNA_hit_count,microRNA_count_db,total_microRNA_gene_denom_count,total_microRNA_gene_hit_count,microRNA_z_scores,microRNA_gene_count,microRNA_symbol_results)
    microRNA_z_score_headers+=['inclusion']

    total_microRNA_gene_hit_count = len(all_excl_microRNA_gene_hits)
    microRNA_z_scores,microRNA_gene_count,microRNA_symbol_results = calculateZScores(excl_microRNA_hit_count,microRNA_count_db,total_microRNA_gene_denom_count,total_microRNA_gene_hit_count,microRNA_z_scores,microRNA_gene_count,microRNA_symbol_results)
    microRNA_z_score_headers+=['exclusion']
    
    total_domain_gene_hit_count = len(all_incl_domain_gene_hits)
    domain_z_scores,domain_gene_count,domain_symbol_results = calculateZScores(incl_domain_hit_count,domain_gene_counts,total_domain_gene_denom_count,total_domain_gene_hit_count,domain_z_scores,domain_gene_count,domain_symbol_results)
    domain_z_score_headers+=['inclusion']

    total_domain_gene_hit_count = len(all_excl_domain_gene_hits)
    domain_z_scores,domain_gene_count,domain_symbol_results = calculateZScores(excl_domain_hit_count,domain_gene_counts,total_domain_gene_denom_count,total_domain_gene_hit_count,domain_z_scores,domain_gene_count,domain_symbol_results)
    domain_z_score_headers+=['exclusion']    

    exportZScoreData(microRNA_z_scores,microRNA_z_score_headers,microRNA_gene_count,microRNA_count_db,microRNA_symbol_results,'microRNA')
    exportZScoreData(domain_z_scores,domain_z_score_headers,domain_gene_count,domain_gene_counts,domain_symbol_results,'ft-domain')

    ############   Export Gene Data ############        
    up_splice_val_genes = 0; down_dI_genes = 0; diff_exp_spliced_genes = 0; diff_spliced_rna_factor = 0
    ddI = 0; udI = 0
    
    critical_gene_exons = eliminate_redundant_dict_values(critical_gene_exons)
    aspire_output_gene = 'AltResults/AlternativeOutput/' + dataset_name + analysis_method + '-exon-inclusion-GENE-results.txt'
    fn=filepath(aspire_output_gene)
    data = open(fn,'w')
    title = ['AffyGene','max_dI','midas-p (corresponding)','symbol','external gene ID','description','regulation_call','event_call']
    title +=['number_of_comparisons','num_critical_exons','up_exons','down_exons','functional_attribute','uniprot-ens_exon_features']
    title +=['go-annotations','mean_fold_change','exon-annotations','exon-region IDs','transcript-cluster-ids','splice-annotation score']
    title = string.join(title,'\t')+'\n'
    data.write(title)
    for affygene in aspire_gene_results:
        if affygene in annotate_db:
            description = annotate_db[affygene].Description()
            symbol = annotate_db[affygene].Symbol()
            ensembl = annotate_db[affygene].ExternalGeneID()
            if array_type == 'exon': transcript_clusters = alternatively_reg_tc[affygene]; transcript_clusters = makeUnique(transcript_clusters); transcript_clusters = string.join(transcript_clusters,'|')
            else: transcript_clusters = affygene
            rna_processing_factor = annotate_db[affygene].RNAProcessing()
        else: description='';symbol='';ensembl=affygene;rna_processing_factor=''; transcript_clusters=''
        if ensembl in go_annotations: goa = go_annotations[ensembl]
        else: goa = ''

        try: gene_splice_event_score[affygene].sort(); top_se_score = str(gene_splice_event_score[affygene][-1])
        except KeyError: top_se_score = 'NA'
        try: gene_regions = gene_exon_region[affygene]; gene_regions = makeUnique(gene_regions); gene_regions = string.join(gene_regions,'|')
        except KeyError: gene_regions = 'NA'
        number_of_comparisons = str(comparison_count[affygene])
        results_list = aspire_gene_results[affygene]
        results_list.sort(); results_list.reverse()
        max_dI = str(results_list[0][0])
        regulation_call = results_list[0][1]
        event_call = results_list[0][2]
        midas_p = results_list[0][-1]
        num_critical_exons = str(len(critical_gene_exons[affygene]))
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
            mean_fold_change = str(gene_expression_diff_db[affygene][0])
            try:
                if abs(float(mean_fold_change)) > log_fold_cutoff: diff_exp_spliced_genes += 1
            except TypeError:
                diff_exp_spliced_genes = diff_exp_spliced_genes
        else: mean_fold_change = 'NC'
        if len(rna_processing_factor) > 2: diff_spliced_rna_factor +=1
        ###Add annotations for where in the gene structure these exons are (according to Ensembl)
        if affygene in external_exon_annot: external_gene_annot = string.join(external_exon_annot[affygene],', ')
        else: external_gene_annot = ''
        values =[affygene,max_dI,midas_p,symbol,ensembl,description,regulation_call,event_call,number_of_comparisons]
        values+=[num_critical_exons,up_exons,down_exons,functional_annotation]
        values+=[uniprot_exon_annotation,goa,mean_fold_change,external_gene_annot,gene_regions,transcript_clusters,top_se_score]
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

    upregulated_genes = 0
    downregulated_genes = 0
    ###Calculate the number of upregulated and downregulated genes
    for affygene in gene_expression_diff_db:
        try: fold_val = gene_expression_diff_db[affygene][0]
        except TypeError: print gene_expression_diff_db[affygene]; kill
        if float(fold_val) > log_fold_cutoff: upregulated_genes += 1
        elif abs(float(fold_val)) > log_fold_cutoff: downregulated_genes += 1

    upregulated_rna_factor = 0; downregulated_rna_factor = 0
    ###Calculate the total number of putative RNA-processing/binding factors differentially regulated
    for affygene in gene_expression_diff_db:
        if len(gene_expression_diff_db[affygene][1]) > 1 and float(gene_expression_diff_db[affygene][0])>log_fold_cutoff:
            upregulated_rna_factor += 1
        elif len(gene_expression_diff_db[affygene][1]) > 1 and abs(float(gene_expression_diff_db[affygene][0]))>log_fold_cutoff:
            downregulated_rna_factor += 1
  
    ###Generate three files for downstream functional summary
    ### functional_annotation_db2 is output to the same function as functional_annotation_db, ranked_uniprot_list_all to get all ranked uniprot annotations,
    ### and ranked_uniprot_list_coding_only to get only coding ranked uniprot annotations
    functional_annotation_db2, ranked_uniprot_list_all = grab_summary_dataset_annotations(protein_exon_feature_db,'','') #functional_attribute_db
    null, ranked_uniprot_list_coding_only = grab_summary_dataset_annotations(protein_exon_feature_db,functional_attribute_db,'') #functional_attribute_db

    ###Sumarize changes in avg protein length for each splice event
    up_protein_list=[];down_protein_list=[]; protein_length_fold_diff=[]
    for [down_protein,up_protein] in protein_length_list:
        down_protein_list.append(down_protein); up_protein_list.append(up_protein)
        up_protein = float(up_protein); down_protein = float(down_protein)
        if up_protein > 10 and down_protein > 10:
            fold_change = up_protein/down_protein; protein_length_fold_diff.append(fold_change)
    median_fold_diff = statistics.median(protein_length_fold_diff)
    try: down_avg=int(statistics.avg(down_protein_list)); up_avg=int(statistics.avg(up_protein_list))
    except ZeroDivisionError: down_avg=0; up_avg=0
    try:
        try:
            down_std=int(statistics.stdev(down_protein_list)); up_std=int(statistics.stdev(up_protein_list))
        except ValueError: ###If 'null' is returned fro stdev
            down_std = 0;up_std = 0
    except ZeroDivisionError:
        down_std = 0;up_std = 0
    if len(down_protein_list)>1 and len(up_protein_list)>1:
        try:
            t,df,tails = statistics.ttest(down_protein_list,up_protein_list,2,3)
            t = abs(t);df = round(df)
            print 'ttest t:',t,'df:',df
            p = str(statistics.t_probability(t,df))
            print dataset_name,p
        except ZeroDivisionError: p = 'NA'
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
    
    summary_results_db[dataset_name[0:-1]] =  udI,ddI,mx,up_splice_val_genes,down_dI_genes,(up_splice_val_genes + down_dI_genes),upregulated_genes, downregulated_genes, diff_exp_spliced_genes, upregulated_rna_factor,downregulated_rna_factor,diff_spliced_rna_factor,down_avg,down_std,up_avg,up_std,p,median_fold_diff,functional_annotation_db

    ###Re-set this variable (useful for testing purposes)
    splice_event_list=[]
    print analysis_method,"results written to:", aspire_output,'\n'
    return summary_results_db, summary_results_db2, aspire_output, aspire_output_gene, len(critical_exon_db)

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

    ### Used to the export relative individual adjusted probesets fold changes used for splicing index values   
    if export_splice_index_values == 'yes':
        summary_output = 'AltResults/RawSpliceData/'+species+'/'+analysis_method+'/'+dataset_name[:-1]+'.txt'
        data = export.createExportFile(summary_output,'AltResults/RawSpliceData/'+species+'/'+analysis_method)
        title = string.join(['gene-probesets']+original_array_names,'\t')+'\n'; data.write(title)
    
    ###original_avg_const_exp_db contains constitutive mean expression values per group: G6953871 [7.71, 7.66]
    ###array_raw_group_values: Raw expression values in list of groups: G7072464@J935416_RC@j_at ([1.79, 2.16, 2.22], [1.68, 2.24, 1.97, 1.92, 2.12])
    ###avg_const_exp_db contains the raw constitutive expression values in a single list
    splicing_index_hash=[]
    for probeset in exon_db:
        ed = exon_db[probeset]
        #include_probeset = ed.IncludeProbeset()
        include_probeset = 'yes'  ###Moved this filter to import of the probeset relationship file
        ###Examines user input parameters for inclusion of probeset types in the analysis
        if include_probeset == 'yes':
            geneid = ed.GeneID()
            if probeset in fold_dbase and geneid in original_avg_const_exp_db:  ###used to search for array_raw_group_values, but when filtered by expression changes, need to filter by adj_fold_dbase
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
                if export_splice_index_values == 'yes':
                    ev = string.join([geneid+'-'+probeset]+si_interim_group_str_db[0]+si_interim_group_str_db[1],'\t')+'\n'; data.write(ev)
                #if ((math.log(group1_mean_ratio,2))*(math.log(group2_mean_ratio,2)))<0: opposite_SI_log_mean = 'yes'
                if (group1_mean_ratio*group2_mean_ratio)<0: opposite_SI_log_mean = 'yes'
                else: opposite_SI_log_mean = 'no'
                try:
                    group_ratio_p = ttestp(group1_ratios,group2_ratios,2,3)
                    #group_ratio_p = statistics.OneWayANOVA([group1_ratios,group2_ratios])
                    splicing_index = group1_mean_ratio-group2_mean_ratio; abs_splicing_index = abs(splicing_index)
                    exp_log_ratio = original_fold_dbase[probeset][1]; abs_log_ratio = abs(original_fold_dbase[probeset][1])
                    ttest_log_ratio = stats_dbase[probeset][1]
                    #if probeset == '2969916': print abs_splicing_index,group_ratio_p,ed.ExonID(),group1_mean_ratio,group2_mean_ratio,math.log(group1_mean_ratio,2),math.log(group2_mean_ratio,2),((math.log(group1_mean_ratio,2))*(math.log(group2_mean_ratio,2))),opposite_SI_log_mean; kill
                    if probeset in midas_db:
                        try: midas_p = float(midas_db[probeset])
                        except ValueError:
                            midas_p = 1
                            #if abs_splicing_index>1 and group_ratio_p < 0.05: print probeset,group_ratio_p, abs_splicing_index;kill
                    else: midas_p = 0
                    #if probeset == '3294457': print ed.GeneID(),ed.ExonID(),probeset,splicing_index,group_ratio_p,midas_p,group1_ratios,group2_ratios;kill

                    if abs_splicing_index>alt_exon_logfold_cutoff and group_ratio_p < p_threshold and midas_p < p_threshold: #and abs_log_ratio>1 and ttest_log_ratio<0.05: ###and ge_threshold_count==2
                        exonid = ed.ExonID(); critical_exon_list = [1,[exonid]]
                        gene_expression_values = original_avg_const_exp_db[geneid]
                        sid = ExonData(splicing_index,probeset,critical_exon_list,geneid,group1_ratios,group2_ratios,exp_log_ratio,ttest_log_ratio,gene_expression_values,group_ratio_p,opposite_SI_log_mean)
                        splicing_index_hash.append((splicing_index,sid))
                except ZeroDivisionError:
                    ###If this occurs, then most likely, the exon and constitutive probeset are the same
                    null = ''
    if export_splice_index_values == 'yes': data.close()
    splicing_index_hash.sort(); splicing_index_hash.reverse()
    print len(splicing_index_hash),"Probesets with evidence of Alternative Splicing"
    p_value_call=''; permute_p_values = {}
    return splicing_index_hash,p_value_call,permute_p_values

def ttestp(list1,list2,tails,variance):
    t,df,tails = statistics.ttest(list1,list2,tails,variance)
    p = statistics.t_probability(t,df)
    return p

def analyzeJunctionSplicing(relative_splicing_ratio):
    if analysis_method == 'linearregres':
        group_sizes = []; original_array_indices = permute_lists[0] ###p[0] is the original organization of the group samples prior to permutation
        for group in original_array_indices: group_sizes.append(len(group))
        
    ### Used to restrict the analysis to a pre-selected set of probesets (e.g. those that have a specifc splicing pattern)
    if len(filtered_probeset_db)>0:
        temp_db={}
        for probeset in relative_splicing_ratio: temp_db[probeset]=[]
        for probeset in temp_db:
            try: filtered_probeset_db[probeset]
            except KeyError: del relative_splicing_ratio[probeset]

    ### Used to the export relative individual adjusted probesets fold changes used for splicing index values   
    if export_splice_index_values == 'yes':
        summary_output = 'AltResults/RawSpliceData/'+species+'/'+analysis_method+'/'+dataset_name[:-1]+'.txt'
        data = export.createExportFile(summary_output,'AltResults/RawSpliceData/'+species+'/'+analysis_method)
        title = string.join(['probesets']+original_array_names,'\t')+'\n'; data.write(title)
        
    ### Calculate a probeset p-value adjusted for constitutive expression levels (taken from splicing index method)
    probeset_group_ratio_p={}
    for probeset in array_raw_group_values:
        ed = exon_db[probeset]; geneid = ed.GeneID()
        group_index = 0; si_interim_group_db={}; si_interim_group_str_db={}; ge_threshold_count=0; value_count = 0
        for group_values in array_raw_group_values[probeset]:
            value_index = 0; ratio_hash=[]; ratio_str_hash=[]
            for value in group_values:  ###Calculate normalized ratio's for each condition and save raw values for later permutation
                exp_val = value;ge_val = avg_const_exp_db[geneid][value_count]; exp_ratio = exp_val-ge_val
                ratio_hash.append(exp_ratio); ratio_str_hash.append(str(exp_ratio)); value_index +=1; value_count +=1
            si_interim_group_db[group_index] = ratio_hash
            si_interim_group_str_db[group_index] = ratio_str_hash; group_index+=1
        group1_ratios = si_interim_group_db[0]; group2_ratios = si_interim_group_db[1]
        if export_splice_index_values == 'yes':
            ev = string.join([probeset]+si_interim_group_str_db[0]+si_interim_group_str_db[1],'\t')+'\n'; data.write(ev)
        try: group_ratio_p = ttestp(group1_ratios,group2_ratios,2,3)
        except ZeroDivisionError: group_ratio_p = 1 ###occurs for constitutive probesets
        probeset_group_ratio_p[probeset]=group_ratio_p ### store and access this below
        #if probeset == 'G6899622@J916374@j_at': print group_ratio_p,group1_ratios,group2_ratios;kill
        ###Concatenate the two raw expression groups into a single list for permutation analysis
        ls_concatenated = []
        for group in array_raw_group_values[probeset]:
            for entry in group: ls_concatenated.append(entry)
        if analysis_method == 'linearregres': ###Convert out of log space
            ls_concatenated = statistics.log_fold_conversion(ls_concatenated)
        array_raw_group_values[probeset] = ls_concatenated
        
    s = 0; t = 0; y = ''
    splice_event_list=[]; splice_event_list_mx=[]; splice_event_list_non_mx=[]; event_mx_temp = []#use this to exclude duplicate mx events
    for affygene in alt_junction_db:
        for event in alt_junction_db[affygene]:
            #event = [('ei', 'E16-E17'), ('ex', 'E16-E18')] 
            #critical_exon_db[affygene,tuple(critical_exons)] = [1,'E'+str(e1a),'E'+str(e2b)] --- affygene,tuple(event) == key, 1 indicates both are either up or down together
            event_call = event[0][0] + '-' + event[1][0]
            exon_set1 = event[0][1]; exon_set2 = event[1][1]
            probeset1 = exon_dbase[affygene,exon_set1]
            probeset2 = exon_dbase[affygene,exon_set2]
            critical_exon_list = critical_exon_db[affygene,tuple(event)]
            #print probeset1,probeset2, critical_exon_list,event_call,exon_set1,exon_set2;kill
            if probeset1 in relative_splicing_ratio and probeset2 in relative_splicing_ratio:
                baseline_ratio1 = relative_splicing_ratio[probeset1][0]
                experimental_ratio1 = relative_splicing_ratio[probeset1][1]
                baseline_ratio2 = relative_splicing_ratio[probeset2][0]
                experimental_ratio2 = relative_splicing_ratio[probeset2][1]
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
                if (Rin<1 and Rex>1): y = 'upregulated'
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
                if (Rin>1 and Rex<1) or (Rin<1 and Rex>1):
                    if analysis_method == 'ASPIRE':
                        ###Calculate ASPIRE scores
                        s +=1
                        in1=((Rex-1.0)*Rin)/(Rex-Rin)
                        in2=(Rex-1.0)/(Rex-Rin)
                        ### dI = ((in1-in2)+(I1-I2))/2.0  #original equation
                        dI = ((in2-in1)+(I2-I1))/2.0 #modified to give propper exon inclusion
                        if r == 1: dI = abs(dI)  ###Occurs when event is mutually exclusive
                        pp1 = probeset_group_ratio_p[probeset1]; pp2 = probeset_group_ratio_p[probeset2]
                        if pp1<p_threshold or pp2<p_threshold: ###Require that the splice event have a constitutive corrected p less than the user defined threshold
                            ejd = ExonJunctionData(dI,probeset1,probeset2,pp1,pp2,y,event_call,critical_exon_list,affygene,baseline_ratio1,experimental_ratio1,baseline_ratio2,experimental_ratio2)
                            splice_event_list.append((dI,ejd))  
                    if analysis_method == 'linearregres':
                        s+=1
                        log_fold, rsqrd_status = performLinearRegression(probeset1,probeset2,group_sizes)
                        if rsqrd_status == 'proceed':
                            pp1 = probeset_group_ratio_p[probeset1]; pp2 = probeset_group_ratio_p[probeset2]
                            if pp1<p_threshold or pp2<p_threshold: ###Require that the splice event have a constitutive corrected p less than the user defined threshold
                                ejd = ExonJunctionData(log_fold,probeset1,probeset2,pp1,pp2,y,event_call,critical_exon_list,affygene,baseline_ratio1,experimental_ratio1,baseline_ratio2,experimental_ratio2)
                                splice_event_list.append((log_fold,ejd))                    
                else: t +=1
    
    print "All scores",len(splice_event_list)
    print "Number of exon-events analyzed:", s
    print "Number of exon-events excluded:", t
    
    splice_event_list2=[]
    for entry in splice_event_list:
        l=0
        if analysis_method == 'linearregres':
            if abs(entry[0])> alt_exon_logfold_cutoff: splice_event_list2.append(entry)
        if analysis_method == 'ASPIRE':
            if abs(entry[0])> aspire_cutoff: splice_event_list2.append(entry)
    
    splice_event_list = splice_event_list2; splice_event_list.sort(); splice_event_list.reverse()
    print "filtered %s scores:" % analysis_method, len(splice_event_list)
    if perform_permutation_analysis == 'yes':
        ###*********BEGIN PERMUTATION ANALYSIS*********
        splice_event_list, p_value_call, permute_p_values = permuteSplicingScores(splice_event_list)
    else:
        p_value_call=''; permute_p_values = {}
    return splice_event_list, p_value_call, permute_p_values


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
    def BaselineRatio1(self): return self._baseline_ratio1
    def BaselineRatio2(self): return self._baseline_ratio2
    def ExperimentalRatio1(self): return self._experimental_ratio1
    def ExperimentalRatio2(self): return self._experimental_ratio2
    def TTestNormalizedRatios(self): return self._group_ratio_p
    def TTestNormalizedRatios2(self): return self._group_ratio_p2
    def EventCall(self):
        ###e.g. Exon inclusion (ei) Exon exclusion (ex), ei-ex, reported in that direction
        return self._event_call
    def Report(self):
        output = self.Method() +'|'+ self.GeneID() +'|'+ self.CriticalExons()
        return output
    def __repr__(self): return self.Report()

class ExonJunctionData(SplicingScoreData):
    def __init__(self,score,probeset1,probeset2,probeset1_p,probeset2_p,regulation_call,event_call,critical_exon_list,affygene,baseline_ratio1,experimental_ratio1,baseline_ratio2,experimental_ratio2):
        self._score = score; self._probeset1 = probeset1; self._probeset2 = probeset2; self._regulation_call = regulation_call
        self._event_call = event_call; self._critical_exon_list = critical_exon_list; self._geneid = affygene
        self._baseline_ratio1 = baseline_ratio1; self._baseline_ratio2 = baseline_ratio2; self._experimental_ratio1 = experimental_ratio1
        self._experimental_ratio2 = experimental_ratio2; self._method = analysis_method; self._group_ratio_p = probeset1_p
        self._group_ratio_p2 = probeset2_p
        
class ExonData(SplicingScoreData):
    def __init__(self,splicing_index,probeset,critical_exon_list,geneid,group1_ratios,group2_ratios,exp_log_ratio,ttest_log_ratio,gene_expression_values,group_ratio_p,opposite_SI_log_mean):
        self._score = splicing_index; self._probeset1 = probeset; self._opposite_SI_log_mean = opposite_SI_log_mean
        self._critical_exon_list = critical_exon_list; self._geneid = geneid
        self._baseline_ratio1 = group1_ratios; self._experimental_ratio1 = group2_ratios
        self._gene_expression_values = gene_expression_values; self._exp_log_ratio = exp_log_ratio
        self._ttest_log_ratio = ttest_log_ratio; self._group_ratio_p = group_ratio_p
        self._method = analysis_method; self._event_call = 'exon-inclusion'
        if splicing_index > 0: regulation_call = 'downregulated'  ###Since baseline is the numerator ratio
        else: regulation_call = 'upregulated'
        self._regulation_call = regulation_call
    #def ExpressionLogRatio(self): return self._exp_log_ratio
    #def TTestLogRatio(self): return self._ttest_log_ratio
    #def GeneExpressionValues(self): return self._gene_expression_values
    def OppositeSIRatios(self): return self._opposite_SI_log_mean
  
    
def performLinearRegression(probeset1,probeset2,group_sizes):
    p1_exp = array_raw_group_values[probeset1]
    p2_exp = array_raw_group_values[probeset2]
    
    p1_g1 = p1_exp[:group_sizes[0]]; p1_g2 = p1_exp[group_sizes[0]:]
    p2_g1 = p2_exp[:group_sizes[0]]; p2_g2 = p2_exp[group_sizes[0]:]

    return_rsqrd = 'no'
    if use_R == 'yes': ###Uses the RLM algorithm
        g1_slope = statistics.LinearRegression(p1_g1,p2_g1,return_rsqrd)
        g2_slope = statistics.LinearRegression(p1_g2,p2_g2,return_rsqrd)
    else: ###Uses a basic least squared method
        g1_slope = statistics.simpleLinRegress(p1_g1,p2_g1)
        g2_slope = statistics.simpleLinRegress(p1_g2,p2_g2)
    
    log_fold = statistics.convert_to_log_fold(g2_slope/g1_slope)
    rsqrd = 'proceed'
    #if g1_rsqrd > 0 and g2_rsqrd > 0: rsqrd = 'proceed'
    #else: rsqrd = 'hault'
    return log_fold,rsqrd

########### Permutation Analysis Functions ###########
def permuteLinearRegression(probeset1,probeset2,p):
    p1_exp = array_raw_group_values[probeset1]
    p2_exp = array_raw_group_values[probeset2]

    p1_g1, p1_g2 = permute_samples(p1_exp,p)
    p2_g1, p2_g2 = permute_samples(p2_exp,p)
    return_rsqrd = 'no'
    g1_slope = statistics.LinearRegression(p1_g1,p2_g1,return_rsqrd)
    g2_slope = statistics.LinearRegression(p1_g2,p2_g2,return_rsqrd)
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
            permute_p_values[(probeset1,probeset2)] = p_val, pos_permute, total_permute, greater_than_true_permute
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
    p_splice_val = statistics.aspire_stringent(b1/bc,e1/ec,b2/bc,e2/ec)
    #print p_splice_val, ref_splice_val, probeset1, probeset2, affygene; dog
    if y == 0:            ###The first permutation is always the real one
        ### Grab the absolute number with small number of decimal places
        try:
            new_ref_splice_val = str(p_splice_val); new_ref_splice_val = float(new_ref_splice_val[0:8])
            ref_splice_val = str(abs(ref_splice_val)); ref_splice_val = float(ref_splice_val[0:8]); y += 1
        except ValueError:
            ###Only get this error if you ref_splice_val is a null
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
    if get_non_log_avg == 'no':
        gb = statistics.avg(baseline); ge = statistics.avg(experimental)  ###Group avg baseline, group avg experimental value
        gb = statistics.log_fold_conversion(gb); ge = statistics.log_fold_conversion(ge)
    else:
        baseline = statistics.log_fold_conversion(baseline); experimental = statistics.log_fold_conversion(experimental)
        gb = statistics.avg(baseline); ge = statistics.avg(experimental)  ###Group avg baseline, group avg experimental value      
    return gb,ge

def format_exon_functional_attributes(affygene,critical_probeset_list,functional_attribute_db,up_exon_list,down_exon_list,protein_length_list):
    ### Add functional attributes
    functional_attribute_list2=[]
    new_functional_attribute_str=''
    new_seq_attribute_str=''
    new_functional_attribute_list=[]
    if array_type == 'exon': critical_probesets = critical_probeset_list[0]
    else: critical_probesets = tuple(critical_probeset_list)
    key = affygene,critical_probesets
    if key in functional_attribute_db:
        ###Grab exon IDs corresponding to the critical probesets
        try: critical_exons = critical_exon_junction_db[critical_probesets] ###For junction arrays
        except KeyError: critical_exons = [exon_db[critical_probesets].ExonID()] ###For exon arrays
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
                    except UnboundLocalError: print entry;kill
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
    new_list=[]
    for item in data_list:
        new_list.append(math.log(float(item),2))
    return new_list
    
def addGlobalFudgeFactor(data_list,data_type):
    new_list = []
    if data_type == 'log':
        for item in data_list:
            new_item = statistics.log_fold_conversion(item)
            new_list.append(float(new_item) + global_addition_factor)
        new_list = convertToLog2(new_list)
    else:
        for item in data_list: new_list.append(float(item) + global_addition_factor)
    return new_list
    
def constitutive_expression_changes(constitutive_fold_change,annotate_db):
    ###Add in constutive fold change filter to assess gene expression for ASPIRE
    gene_expression_diff_db = {}
    rna_processing_regulated = {}
    for affygene in constitutive_fold_change:
        constitutive_fold = constitutive_fold_change[affygene]
        try: gene_expression_diff_db[affygene].append(constitutive_fold)
        except KeyError: gene_expression_diff_db[affygene] = [constitutive_fold]
        ###Add in evaluation of RNA-processing/binding factors
        #functional_annotation_db[affygene] = name, symbol,ll_id,splicing_annotation
        if affygene in annotate_db:
            if len(annotate_db[affygene].RNAProcessing()) > 4: rna_processing_regulated[affygene] = annotate_db[affygene].RNAProcessing()

    ###Select the largest fold change
    for affygene in gene_expression_diff_db:
        if len(gene_expression_diff_db[affygene]) > 1:
            fold_list = gene_expression_diff_db[affygene]
            fold_list.sort(); fold_list.reverse()
            if abs(fold_list[0]) > abs(fold_list[-1]): new_list = [fold_list[0]]
            else: new_list = [fold_list[-1]]
            if affygene in rna_processing_regulated:
                new_list.append(rna_processing_regulated[affygene])
                gene_expression_diff_db[affygene] = new_list
            else:
                new_list.append('')
                gene_expression_diff_db[affygene] = new_list
        else:
            new_list = gene_expression_diff_db[affygene]
            if affygene in rna_processing_regulated:
                new_list.append(rna_processing_regulated[affygene])
            else:
                new_list.append('')
                gene_expression_diff_db[affygene] = new_list
    return gene_expression_diff_db

def restrictProbesets():
    ### Take a file with probesets and only perform the splicing-analysis on these (e.g. those already identified from a previous run with a specific pattern)
    ### Allows for propper denominator when calculating z-scores for microRNA and protein-domain ORA
    probeset_list_filename = import_dir = 'AltDatabase/filtering/probesets.txt'
    filtered_probeset_db = importGeneric(probeset_list_filename)
    print len(filtered_probeset_db), "probesets will be used to restrict analysis..."
    return filtered_probeset_db

def RunAltAnalyze():

  if array_type == 'AltMouse': import_dir = '/AltExpression/'+array_type
  elif array_type == 'exon': import_dir = '/AltExpression/ExonArray/'+species+'/'
  
  if analysis_method == 'ASPIRE' or analysis_method == 'linearregres' or analysis_method == 'splicing-index':
    global microRNA_full_exon_db; global microRNA_count_db; global gene_microRNA_denom; global protein_exon_feature_db; global protein_sequence_dbase; global probeset_protein_db
    global protein_ft_db; global domain_gene_counts; global annotate_db; annotate_db={}; global splice_event_list; splice_event_list=[]; global critical_exon_junction_db
    global midas_db; global dataset_name; global exon_db; global constituitive_probeset_db
    
    if array_type == 'exon': gene_annotation_file = "AltDatabase/ensembl/"+species+"/"+species+"_Ensembl-annotations.txt"
    else: gene_annotation_file = "AltDatabase/"+species+"/"+array_type+"/"+array_type+"_gene_annotations.txt"
    annotate_db = import_annotations(gene_annotation_file,array_type)
    if export_go_annotations == 'yes': exportGOannotations(annotate_db)

    ###Import probe-level associations    
    exon_db={}; filtered_arrayids={};filter_status='no'
    constituitive_probeset_db,exon_db,genes_being_analyzed = importSplicingAnnotationDatabase(probeset_annotations_file,array_type,filtered_arrayids,filter_status)

    if array_type != 'exon': critical_exon_junction_db = ExonAnalyze_module.grabJunctionData(species,array_type,'probeset-pairs') ###used to translate probeset pairs into critical_exons
    else: critical_exon_junction_db = {}
      
    if analyze_functional_attributes == 'yes':
      exon_protein_sequence_file = "AltDatabase/"+species+"/"+array_type+"/"+"SEQUENCE-protein-dbase.txt"
      probeset_protein_db,protein_sequence_dbase = ExonAnalyze_module.importExonSequenceBuild(exon_protein_sequence_file,exon_db)
      protein_ft_db,domain_gene_counts = FeatureAlignment.grab_exon_level_feature_calls(probeset_protein_db,species,array_type,exon_db,genes_being_analyzed)
      ###Could use this directly for exon arrays by building up-front rather than on-the-fly.  Doesn't appear to be an analagous function for ensembl (that I wrote would have wrote before-NS)
      
      ###The 'domain_gene_counts' is not currently only for anything, but reports the number of domains for each gene on the array for ensembl and uniprot combined. Reports a tuple of (feature_type, specific_annotation) related to number of genes with it.
      start_time = time.time()
      print "Uniprot & Ensembl feature & domain data imported"
      print "MicroRNA data imported"
      microRNA_full_exon_db,microRNA_count_db,gene_microRNA_denom = ExonAnalyze_module.importmicroRNADataExon(species,array_type,exon_db,microRNA_prediction_method)
    else:
      probeset_protein_db={}; gene_microRNA_denom={}; domain_gene_counts={}
      protein_ft_db={}; protein_sequence_dbase ={};microRNA_count_db={}
      microRNA_full_exon_db = {}
  #"""

  run=1
  dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
  for altanalzye_input in dir_list:    #loop through each file in the directory to output results
    ###Import probe-level associations
    if run>1: ### Only re-set these databases after the run when batch analysing multiple files
        exon_db={}; filtered_arrayids={};filter_status='no' ###Use this as a means to save memory (import multiple times - only storing different types relevant information)
        constituitive_probeset_db,exon_db,genes_being_analyzed = importSplicingAnnotationDatabase(probeset_annotations_file,array_type,filtered_arrayids,filter_status)

    array_db = import_dir + "/"+ altanalzye_input
    array_db = array_db[1:] #not sure why, but the '\' needs to be there while reading initally but not while accessing the file late
    dataset_name = altanalzye_input[0:-4] + '-'
    print "Begining to process",dataset_name[0:-1]

    ###Import expression data and stats and filter the expression data based on fold and p-value OR expression threshold
    conditions,relative_splicing_ratio,stats_dbase,adj_fold_dbase,dataset_name,gene_expression_diff_db,midas_db = expr_analysis(array_db,constituitive_probeset_db,exon_db,annotate_db,dataset_name)
    ###Run Analysis
    if conditions == 2 and exportTransitResultsforAnalysis == 'no':
        summary_results_db, summary_results_db2, aspire_output, aspire_output_gene, number_events_analyzed = splicing_analysis_algorithms(relative_splicing_ratio,adj_fold_dbase,dataset_name,gene_expression_diff_db,exon_db)
        aspire_output_list.append(aspire_output); aspire_output_gene_list.append(aspire_output_gene)
    elif conditions != 2: print "Analysis not run...too many conditions to currently analyze"
    else: summary_results_db={}; number_events_analyzed = 0
    run+=1

  return summary_results_db, aspire_output_gene_list, number_events_analyzed

if __name__ == '__main__':
  ### Hard-coded defaults
  w = 'Agilent'; x = 'Affymetrix'; y = 'Ensembl'; z = 'any'; data_source = y; constitutive_source = z; manufacturer = x ### Constitutive source, is only really paid attention to if Ensembl, otherwise Affymetrix is used (even if default)
  use_R = 'no'
  
  ### Get default options for ExpressionBuilder and AltAnalyze
  expr_var, alt_var, additional_var, file_location_defaults = UI.getUserParameters()
  start_time = time.time()
  
  species,array_type,manufacturer,constitutive_source,dabg_p,raw_expression_threshold,avg_all_for_ss,expression_data_format,include_raw_data, run_from_scratch = expr_var
  analysis_method,p_threshold,filter_probeset_types,alt_exon_fold_variable,gene_expression_cutoff,permute_p_threshold,perform_permutation_analysis, export_splice_index_values = alt_var
  exportTransitResultsforAnalysis, analyze_functional_attributes, microRNA_prediction_method = additional_var

  print "Expression Analysis Parameters Being Used..."
  print '\t'+'species'+':',species
  print '\t'+'array_type'+':',array_type
  print '\t'+'manufacturer'+':',manufacturer
  print '\t'+'constitutive_source'+':',constitutive_source
  print '\t'+'dabg_p'+':',dabg_p
  print '\t'+'raw_expression_threshold'+':',raw_expression_threshold
  print '\t'+'avg_all_for_ss'+':',avg_all_for_ss
  print '\t'+'expression_data_format'+':',expression_data_format
  print '\t'+'include_raw_data'+':',include_raw_data
  print '\t'+'run_from_scratch'+':',run_from_scratch

  print "Alternative Exon Analysis Parameters Being Used..."  
  print '\t'+'analysis_method'+':',analysis_method
  print '\t'+'p_threshold'+':',p_threshold
  print '\t'+'filter_probeset_types'+':',filter_probeset_types
  print '\t'+'alt_exon_fold_variable'+':',alt_exon_fold_variable
  print '\t'+'gene_expression_cutoff'+':',gene_expression_cutoff
  print '\t'+'avg_all_for_ss'+':',avg_all_for_ss
  print '\t'+'permute_p_threshold'+':',permute_p_threshold
  print '\t'+'perform_permutation_analysis'+':',perform_permutation_analysis
  print '\t'+'export_splice_index_values'+':',export_splice_index_values
  print '\t'+'exportTransitResultsforAnalysis'+':',exportTransitResultsforAnalysis
  print '\t'+'analyze_functional_attributes'+':',analyze_functional_attributes
  print '\t'+'microRNA_prediction_method'+':',microRNA_prediction_method
      
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

  if run_from_scratch == 'raw input':
      ExpressionBuilder.remoteExpressionBuilder(species,array_type,dabg_p,raw_expression_threshold,avg_all_for_ss,expression_data_format,manufacturer,constitutive_source,data_source,include_raw_data)
  elif run_from_scratch == 'update DBs':
      null=[] ###Add link to new module
      #updateDBs(species,array_type)
   
  ###### Run AltAnalyze ######  
  global dataset_name; global summary_results_db; global summary_results_db2
  summary_results_db={}; summary_results_db2={}; aspire_output_list=[]; aspire_output_gene_list=[]

  onlyAnalyzeJunctions = 'no'; agglomerate_inclusion_probesets = 'no'; filter_probesets_by = 'NA'
  if array_type == 'AltMouse':
      if filter_probeset_types == 'junctions-only':  onlyAnalyzeJunctions = 'yes'
      elif filter_probeset_types == 'combined-junctions': agglomerate_inclusion_probesets = 'yes'; onlyAnalyzeJunctions = 'yes'
      elif filter_probeset_types == 'exons-only': analysis_method = 'splicing-index'; filter_probesets_by = 'exon'
  else: filter_probesets_by = filter_probeset_types

  c = 'Ensembl'; d = 'Entrez Gene'
  annotation_system = c
  expression_threshold = 0 ###This is different than the raw_expression_threshold (probably shouldn't filter so set to 0)

  if analysis_method == 'ASPIRE': aspire_cutoff = alt_exon_fold_variable
  if analysis_method == 'linearregres-rlm': analysis_method = 'linearregres';use_R = 'yes'

  if 'APT' in file_location_defaults: apt_location = file_location_defaults['APT'].Location()
  log_fold_cutoff = math.log(float(gene_expression_cutoff),2)
  alt_exon_logfold_cutoff = math.log(float(alt_exon_fold_variable),2)
      
  global_addition_factor = 0
  get_non_log_avg = 'no'
  export_junction_comparisons = 'no' ### No longer accessed in this module - only in update mode through a different module
  
  factor_out_expression_changes = 'yes' ### Use 'no' if data is normalized already or no expression normalization for ASPIRE desired
  only_include_constitutive_containing_genes = 'yes'
  remove_transcriptional_regulated_genes = 'yes'
  add_exons_to_annotations = 'no'
  exclude_protein_details = 'no'
  export_go_annotations = 'no'

  use_external_file_for_probeset_filtering = 'no'

  if analysis_method != 'splicing-index': annotation_system = d
  if array_type == 'AltMouse': species = 'Mm'
  if export_splice_index_values == 'yes': remove_transcriptional_regulated_genes = 'no'

  global filtered_probeset_db; filtered_probeset_db={}
  if use_external_file_for_probeset_filtering == 'yes': filtered_probeset_db = restrictProbesets()

  print "Parsing out Affymetrix GO annotations"
  global go_annotations; go_annotations={}
  ###Saves run-time while testing the software (global variable stored)
  import_dir = '/AltDatabase/affymetrix/'+species
  dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
  for affy_data in dir_list:    #loop through each file in the directory to output results
      affy_data_dir = 'AltDatabase/affymetrix/'+species+'/'+affy_data
      parse_affymetrix_annotations(affy_data_dir)  

  if array_type == 'exon': probeset_annotations_file = "AltDatabase/"+species+"/"+array_type+"/"+species+"_Ensembl_probesets.txt"
  else: probeset_annotations_file = "AltDatabase/"+species+"/"+array_type+"/"+"MASTER-probeset-transcript.txt"

  ### If exportTransitResultsforAnalysis == 'yes', modify the user options initially to export all probe set data for MiDAS
  if exportTransitResultsforAnalysis == 'yes':
      analyze_functional_attributes = 'no'; export_splice_index_values = 'no'; permute_p_threshold = 'no'
      if analysis_method == 'splicing-index': filter_probesets_by = 'full'
      else: onlyAnalyzeJunctions = 'no'; agglomerate_inclusion_probesets = 'no'
      print analysis_method, filter_probeset_types
      summary_results_db, aspire_output_gene_list, number_events_analyzed = RunAltAnalyze()
      
      ### Once MiDAS data is exported, reset the options to the original parameters and re-run
      analysis_method,p_threshold,filter_probeset_types,alt_exon_fold_variable,gene_expression_cutoff,permute_p_threshold,perform_permutation_analysis, export_splice_index_values = alt_var
      exportTransitResultsforAnalysis, analyze_functional_attributes, microRNA_prediction_method = additional_var
      exportTransitResultsforAnalysis='no'

      if array_type == 'AltMouse':
          if filter_probeset_types == 'junctions-only':  onlyAnalyzeJunctions = 'yes'
          elif filter_probeset_types == 'combined-junctions': agglomerate_inclusion_probesets = 'yes'; onlyAnalyzeJunctions = 'yes'
          elif filter_probeset_types == 'exons-only': analysis_method = 'splicing-index'; filter_probesets_by = 'exon'
      else: filter_probesets_by = filter_probeset_types

      a = r""+apt_location
      try:
          os.spawnl(os.P_NOWAIT, a) ###Open a APT command prompt
          print "\nAn Affymetrix Power Tools (APT) command prompt has been opened. Perform the following steps to build MiDAS statistics:"
      except OSError: print '\nBegin running the Affymetrix Power Tools (APT) Application and perform the following steps:'
      print '     1) Open each one of the new "command-*" files in the "AltResults/MIDAS" directory'
      print '     2) In the APT prompt, enter the first line in the "command" file (change direcotries to "AltResults/MIDAS")'
      print '     3) In the APT prompt, enter the first second in the "command" file (builds MIDAS p-value file to be read by AltAnalyze)'
      print '\nHit the Enter/Return key, when finished exporting MiDAS output to continue analysis'
      inp = sys.stdin.readline()
      
  summary_results_db, aspire_output_gene_list, number_events_analyzed = RunAltAnalyze()

  try:
    if exportTransitResultsforAnalysis == 'no':
        ResultsExport_module.output_summary_results(summary_results_db,'',analysis_method)
        #ResultsExport_module.output_summary_results(summary_results_db2,'-uniprot_attributes',analysis_method)
        ResultsExport_module.compare_ASPIRE_results(aspire_output_list,annotate_db,number_events_analyzed,'no',analysis_method,array_type)
        ResultsExport_module.compare_ASPIRE_results(aspire_output_gene_list,annotate_db,'','yes',analysis_method,array_type) 
  except UnboundLocalError: print "...No results to summarize" ###Occurs if there is a problem parsing these files
  
  end_time = time.time(); time_diff = int(end_time-start_time)
  print "Analyses finished in %d seconds" % time_diff
  print "Hit Enter/Return to exit AltAnalyze"
  inp = sys.stdin.readline(); sys.exit()
