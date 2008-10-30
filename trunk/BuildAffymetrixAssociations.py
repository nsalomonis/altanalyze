import sys, string
import os.path
import unique
import datetime
################# Parse directory files
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

def returnDirectories(sub_dir):
    dir=os.path.dirname(dirfile.__file__)
    dir = string.replace(dir,py2app_adj,'')
    dir = string.replace(dir,'\\library.zip','')
    dir_list = os.listdir(dir + sub_dir)
    ###Below code used to prevent folder names from being included
    dir_list2 = []
    for i in dir_list:
        if "." not in i: dir_list2.append(i)
    return dir_list2

def cleanUpLine(line):
    data = string.replace(line,'\n','')
    data = string.replace(data,'\c','')
    data = string.replace(data,'\r','')
    data = string.replace(data,'"','')
    return data

class GrabFiles:
    def setdirectory(self,value): self.data = value
    def display(self): print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        file = getDirectoryFiles(self.data,str(search_term))
        if len(file)<1: print search_term,'not found'
        return file
    def returndirectory(self):
        dir_list = read_directory(self.data)
        return dir_list
    
def getDirectoryFiles(import_dir, search_term):
    exact_file = ''
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = import_dir[1:]+'/'+data
        if search_term in affy_data_dir: exact_file = affy_data_dir
    return exact_file

################# Import and Annotate Data
class AffymetrixInformation:
    def __init__(self,probeset,symbol,ensembl,entrez,unigene,uniprot,description,goids,go_names,pathways):
        self._probeset = probeset; self._ensembl = ensembl; self._uniprot = uniprot; self._description = description
        self._entrez = entrez; self._symbol = symbol; self._unigene = unigene
        self._goids = goids; self._go_names = go_names; self._pathways = pathways
    def ArrayID(self): return self._probeset
    def Description(self): return self._description
    def Symbol(self): return self._symbol
    def Ensembl(self): return self._ensembl
    def EnsemblString(self):
        ens_str = string.join(self._ensembl,'|')
        return ens_str
    def Entrez(self): return self._entrez
    def EntrezString(self):
        entrez_str = string.join(self._entrez,'|')
        return entrez_str
    def Unigene(self): return self._unigene
    def UnigeneString(self):
        unigene_str = string.join(self._unigene,'|')
        return unigene_str
    def Uniprot(self): return self._uniprot
    def GOIDs(self): return self._goids
    def GOProcessIDs(self): return self._goids[0]
    def GOComponentIDs(self): return self._goids[1]
    def GOFunctionIDs(self): return self._goids[2]
    def GOProcessNames(self):
        go_names = string.join(self._go_names[0],' // ')
        return go_names
    def GOComponentNames(self):
        go_names = string.join(self._go_names[1],' // ')
        return go_names
    def GOFunctionNames(self):
        go_names = string.join(self._go_names[2],' // ')
        return go_names
    def GONames(self):
        go_names = self._go_names[0]+self._go_names[1]+self._go_names[2]
        return go_names
    def Pathways(self): return self._pathways
    def PathwayInfo(self):
        pathway_str = string.join(self.Pathways(),' // ')
        return pathway_str
    def GOPathwayInfo(self):
        pathway_str = string.join(self.GONames() + self.Pathways(),' // ')
        return pathway_str
    def resetEnsembl(self,ensembl): self._ensembl = ensembl
    def resetEntrez(self,entrez): self._entrez = entrez
    def ArrayValues(self):
        output = self.Symbol()+'|'+self.ArrayID()
        return output
    def __repr__(self): return self.ArrayValues()

class InferredEntrezInformation:
    def __init__(self,symbol,entrez,description):
        self._entrez = entrez; self._symbol = symbol; self._description = description
    def Entrez(self): return self._entrez
    def Symbol(self): return self._symbol
    def Description(self): return self._description
    def DataValues(self):
        output = self.Symbol()+'|'+self.Entrez()
        return output
    def __repr__(self): return self.DataValues()
    
def eliminate_redundant_dict_values(database):
    db1={}
    for key in database: list = unique.unique(database[key]); list.sort(); db1[key] = list
    return db1

######### Import New Data ##########
def parse_affymetrix_annotations(filename):
    ###Import an Affymetrix array annotation file (from http://www.affymetrix.com) and parse out annotations
    temp_affy_db = {}; x=0; y=0
    fn=filepath(filename)
    for line in open(fn,'r').readlines():             
        probeset_data,null = string.split(line,'\n')  #remove endline
        probeset_data = string.replace(probeset_data,'---','')
        affy_data = string.split(probeset_data[1:-1],'","')  #remove endline
        if x==0 and line[0]!='#':
            x=1; affy_headers = affy_data
            for header in affy_headers:
                y = 0
                while y < len(affy_headers):
                    if 'Probe Set ID' in affy_headers[y] or 'probeset_id' in affy_headers[y]: ps = y
                    if 'transcript_cluster_id' in affy_headers[y]: tc = y
                    if 'Gene Symbol' in affy_headers[y]: gs = y
                    if 'Ensembl' in affy_headers[y]: ens = y
                    if ('nigene' in affy_headers[y] or 'UniGene' in affy_headers[y]) and 'Cluster' not in affy_headers[y]: ug = y
                    if 'mrna_assignment' in affy_headers[y]: ma = y
                    if 'gene_assignment' in affy_headers[y]: ga = y
                    if 'Entrez' in affy_headers[y] or 'LocusLink' in affy_headers[y]: ll = y
                    if 'SwissProt' in affy_headers[y] or 'swissprot' in affy_headers[y]: sp = y
                    if 'Gene Title' in affy_headers[y]: gt = y
                    if 'rocess' in affy_headers[y]: bp = y
                    if 'omponent' in affy_headers[y]: cc = y
                    if 'unction' in affy_headers[y]: mf = y
                    if 'athway' in affy_headers[y]: gp = y
                    y += 1
        elif x == 1:
            ###If using the Affy 2.0 Annotation file structure, both probeset and transcript cluster IDs are present
            ###If transcript_cluster centric file (gene's only no probesets), then probeset = transcript_cluster
            try: 
                transcript_cluster = affy_data[tc]  ###Affy's unique Gene-ID
                probeset = affy_data[ps]
                if probeset != transcript_cluster: ###Occurs for transcript_cluster ID centered files
                    probesets = [probeset,transcript_cluster]
                else: probesets = [probeset]
                ps = tc; version = 2 ### Used to define where the non-UID data exists
            except UnboundLocalError: probesets = [affy_data[ps]]; version = 1
            
            uniprot = affy_data[sp]; unigene = affy_data[ug]; uniprot_list = string.split(uniprot,' /// ')
            symbol = ''; description = ''; pathway_data = affy_data[gp]
            for probeset in probesets:
                if version == 1: ###Applies to 3' biased arrays only (Conventional Format)
                    description = affy_data[gt]; symbol = affy_data[gs]; goa=''; entrez = affy_data[ll]
                    ensembl_list = []; ensembl = affy_data[ens]; ensembl_list = string.split(ensembl,' /// ')
                    entrez_list = string.split(entrez,' /// '); unigene_list = string.split(unigene,' /// ')
                    uniprot_list = string.split(uniprot,' /// ')
                    ###Process GO information if desired
                    if process_go == 'yes':
                        process = affy_data[bp]; component = affy_data[cc]; function = affy_data[mf]
                        process_goids, process_names = extractPathwayData(process,'GO',version)
                        component_goids, component_names = extractPathwayData(component,'GO',version)
                        function_goids, function_names = extractPathwayData(function,'GO',version)
                        goids = [process_goids,component_goids,function_goids]
                        go_names = [process_names,component_names,function_names]
                    else: goids=[]; go_names=[]
                    if extract_pathway_names == 'yes': null, pathways = extractPathwayData(pathway_data,'pathway',version)
                    else: pathways = []
                    ai = AffymetrixInformation(probeset, symbol, ensembl_list, entrez_list, unigene_list, uniprot_list, description, goids, go_names, pathways)
                    if len(entrez_list)<5: affy_annotation_db[probeset] = ai
                    #symbol_list = string.split(symbol,' /// '); description_list = string.split(description,' /// ')
                    if len(entrez_list)<2: ###Only store annotations for EntrezGene if there is only one listed ID, since the symbol, description and Entrez Gene are sorted alphabetically, not organized relative to each other (stupid)
                        iter = 0
                        for entrez in entrez_list:
                            #symbol = symbol_list[iter]; description = description_list[iter] ###grab the symbol that matches the EntrezGene entry
                            z = InferredEntrezInformation(symbol,entrez,description)
                            try: entrez_annotation_db[entrez] = z
                            except NameError: null=[]
                            iter += 1
                else: ### Applies to Exon, Transcript, whole geneome Gene arrays.
                    uniprot_list2 = []
                    for uniprot_id in uniprot_list:
                        if len(uniprot_id)>0:
                            try: a = int(uniprot_id[1]); uniprot_list2.append(uniprot_id)
                            except ValueError: null = []
                    uniprot_list = uniprot_list2
                    
                    mrna_associations = affy_data[ma]; ensembl_list=[]; descriptions=[]
                    ensembl_data = string.split(mrna_associations,' /// ')
                    for entry in ensembl_data:
                        annotations = string.split(entry,' // ')
                        #if probeset == '8148358': print annotations
                        for i in annotations:
                            if 'gene:ENS' in i:
                                description, ensembl_id = string.split(i,'gene:ENS'); descriptions.append((len(description),description))
                                ensembl_id = string.split(ensembl_id,' ')
                                ensembl_id = 'ENS'+ensembl_id[0]; ensembl_list.append(ensembl_id)

                    #if probeset == '8148358': print ensembl_list; kill
                    gene_assocs = affy_data[ga]; entrez_list=[]
                    entrez_data = string.split(gene_assocs,' /// ')
                    for entry in entrez_data:
                        try:
                            if len(entry)>0:
                                annotations = string.split(entry,' // ')
                                entrez_gene = int(annotations[-1]); entrez_list.append(str(entrez_gene))
                                symbol = annotations[1]; description = annotations[2]; descriptions.append((len(description),description))
                                #print entrez_gene,symbol, descriptions;kill
                                z = InferredEntrezInformation(symbol,entrez_gene,description)
                                try: entrez_annotation_db[str(entrez_gene)] = z ###create an inferred Entrez gene database
                                except NameError: null = []
                        except ValueError: null = []
                    gene_assocs = affy_data[ga]; gene_assocs = string.replace(gene_assocs,'---','')
                    unigene_data = string.split(unigene,' /// '); unigene_list = []
                    for entry in unigene_data:
                        if len(entry)>0:
                            annotations = string.split(entry,' // ')
                            try: null = int(annotations[-2][3:]); unigene_list.append(annotations[-2])
                            except ValueError: null = []
                    ###Only applies to optional GOID inclusion
                    if process_go == 'yes':
                        process = affy_data[bp]; component = affy_data[cc]; function = affy_data[mf]
                        process_goids, process_names = extractPathwayData(process,'GO',version)
                        component_goids, component_names = extractPathwayData(component,'GO',version)
                        function_goids, function_names = extractPathwayData(function,'GO',version)
                        goids = [process_goids,component_goids,function_goids]
                        go_names = [process_names,component_names,function_names]
                    else: goids=[]; go_names=[]
                    if extract_pathway_names == 'yes': null, pathways = extractPathwayData(pathway_data,[],version)
                    else: pathways = []
                    entrez_list=unique.unique(entrez_list); unigene_list=unique.unique(unigene_list); uniprot_list=unique.unique(uniprot_list); ensembl_list=unique.unique(ensembl_list)
                    descriptions2=[]
                    for i in descriptions: 
                        if 'cdna:known' not in i: descriptions2.append(i)
                    descriptions = descriptions2
                    if len(descriptions)>0:
                        descriptions.sort(); description = descriptions[-1][1]
                        if description[0] == ' ': description = description[1:] ### some entries begin with a blank
                    ai = AffymetrixInformation(probeset,symbol,ensembl_list,entrez_list,unigene_list,uniprot_list,description,goids,go_names,pathways)
                    if len(entrez_list)<5 and len(ensembl_list)<5: affy_annotation_db[probeset] = ai

def extractPathwayData(terms,type,version):
    goids = []; go_names = []
    buffer = ' /// '; small_buffer = ' // '
    go_entries = string.split(terms,buffer)
    for go_entry in go_entries:
        go_entry_info = string.split(go_entry,small_buffer)
        try:
            if len(go_entry_info)>1:
                if version == 1:
                    if type == 'GO': ### 6310 // DNA recombination // inferred from electronic annotation ///
                        goid, go_name, source = go_entry_info
                        while len(goid)< 7: goid = '0'+goid
                        goid = 'GO:'+goid
                    else: ### Calcium signaling pathway // KEGG ///
                        go_name, source = go_entry_info; goid = ''
                        if len(go_name)>1: go_name = source+'-'+go_name
                        else: go_name = ''
                if version == 2:
                    if type == 'GO': ### NM_153254 // GO:0006464 // protein modification process // inferred from electronic annotation /// 
                        accession, goid, go_name, source = go_entry_info
                    else: ### AF057061 // GenMAPP // Matrix_Metalloproteinases
                        accession, go_name, source = go_entry_info; goid = ''
                        if len(go_name)>1: go_name = source+'-'+go_name
                        else: go_name = ''
                goids.append(goid); go_names.append(go_name)
        except IndexError: goids = goids
    if extract_go_names != 'yes': go_names = [] ### Then don't store (save memory)
    goids = unique.unique(goids); go_names = unique.unique(go_names)
    return goids, go_names
          
def exportResultsSummary(dir_list):
    new_file = 'Databases/BuildDBs/BuiltDBs/'+species+'/Affymetrix_files_summarized.txt'
    today = str(datetime.date.today()); today = string.split(today,'-'); today = today[1]+'/'+today[2]+'/'+today[0]
    fn=filepath(new_file); data = open(fn,'w')
    for filename in dir_list: data.write(filename+'\t'+today+'\n')
    data.close()
    
def exportRelationshipDBs(species):
    new_file1 = 'Databases/BuildDBs/BuiltDBs/'+species+'/Ensembl-Affymetrix.txt'
    new_file2 = 'Databases/BuildDBs/BuiltDBs/'+species+'/EntrezGene-Affymetrix.txt'
    new_file3 = 'Databases/BuildDBs/BuiltDBs/'+species+'/EntrezGene.txt'
    today = str(datetime.date.today()); today = string.split(today,'-'); today = today[1]+'/'+today[2]+'/'+today[0]
    try: fn1=filepath(new_file1); data1 = open(fn1,'w')
    except IOError:
        try: new_dir = 'Databases/BuildDBs/BuiltDBs/'+species; fn = filepath(new_dir);os.mkdir(fn) ###Re-Create directory if deleted
        except OSError:
            new_dir = 'Databases/BuildDBs/BuiltDBs'; fn = filepath(new_dir);os.mkdir(fn)
            new_dir = 'Databases/BuildDBs/BuiltDBs/'+species; fn = filepath(new_dir);os.mkdir(fn)
        fn1=filepath(new_file1); data1 = open(fn1,'w')
    fn2=filepath(new_file2); data2 = open(fn2,'w')
    fn3=filepath(new_file3); data3 = open(fn3,'w')
    data3.write('ID'+'\t'+'Symbol'+'\t'+'Name'+'\t'+'Species'+'\t'+'Date'+'\t'+'Remarks'+'\n')
    if process_go == 'yes':
        new_file4 = 'Databases/BuildDBs/BuiltDBs/'+species+'/Ensembl-GeneOntology.txt'
        new_file5 = 'Databases/BuildDBs/BuiltDBs/'+species+'/EntrezGene-GeneOntology.txt'
        header = 'GeneID'+'\t'+'GOID'+'\n'
        fn4=filepath(new_file4); data4 = open(fn4,'w'); data4.write(header)
        fn5=filepath(new_file5); data5 = open(fn5,'w'); data5.write(header)
    for probeset in affy_annotation_db:
        ai = affy_annotation_db[probeset]
        ensembls = unique.unique(ai.Ensembl()); entrezs = unique.unique(ai.Entrez())
        for ensembl in ensembls:
            if len(ensembl)>0:
                data1.write(ensembl+'\t'+probeset+'\n')
                if process_go == 'yes':
                    goids = unique.unique(ai.GOIDs())
                    for goid in goids: data4.write(ensembl+'\t'+goid+'\n')
        for entrez in entrezs:
            if len(entrez)>0:
                data2.write(entrez+'\t'+probeset+'\n')
                if process_go == 'yes':
                    goids = unique.unique(ai.GOIDs())
                    for goid in goids: data5.write(entrez+'\t'+goid+'\n')
    for entrez in entrez_annotation_db:
        ea = entrez_annotation_db[entrez]
        if len(entrez)>0: data3.write(entrez+'\t'+ea.Symbol()+'\t'+ea.Description()+'\t'+species+'\t'+today+'\t'+''+'\n')
    data1.close(); data2.close(); data3.close()

def swapKeyValues(db):
    swapped={}
    for key in db:
        values = list(db[key]) ###If the value is not a list, make a list
        for value in values:
            try: swapped[value].append(key)
            except KeyError: swapped[value] = [key]
    swapped = eliminate_redundant_dict_values(swapped)
    return swapped

def integratePreviousAssociations():
    print 'Integrating associations from previous databases...'
    #print len(entrez_annotations),len(ens_to_uid), len(entrez_to_uid)
    for gene in entrez_annotations:
        if gene not in entrez_annotation_db:
            ### Add previous gene information to the new database
            y = entrez_annotations[gene]
            z = InferredEntrezInformation(y.Symbol(),gene,y.Description())
            entrez_annotation_db[gene] = z
    uid_to_ens = swapKeyValues(ens_to_uid); uid_to_entrez = swapKeyValues(entrez_to_uid)

    ###Add prior missing gene relationships for all probesets in the new database and that are missing    
    for uid in uid_to_ens:
        if uid in affy_annotation_db:
            y = affy_annotation_db[uid]
            ensembl_ids = uid_to_ens[uid]
            if y.Ensembl() == ['']: y.resetEnsembl(ensembl_ids)
        else:
            ensembl_ids = uid_to_ens[uid]; entrez_ids = []
            if uid in uid_to_entrez: entrez_ids = uid_to_entrez[uid]
            ai = AffymetrixInformation(uid, '', ensembl_ids, entrez_ids, [], [], '')
            affy_annotation_db[uid] = ai
    for uid in uid_to_entrez:
        if uid in affy_annotation_db:
            y = affy_annotation_db[uid]
            entrez_ids = uid_to_entrez[uid]
            if y.Entrez() == ['']: y.resetEntrez(entrez_ids)
        else:
            entrez_ids = uid_to_entrez[uid]; ensembl_ids = []
            if uid in uid_to_ens: ensembl_ids = uid_to_ens[uid]
            ai = AffymetrixInformation(uid, '', ensembl_ids, entrez_ids, [], [], '')
            affy_annotation_db[uid] = ai

def parseGene2GO(tax_id,species):
    import_dir = '/Databases/BuildDBs/Entrez/Gene2GO'
    g = GrabFiles(); g.setdirectory(import_dir)
    filename = g.searchdirectory('gene2go') ###Identify gene files corresponding to a particular MOD
    fn=filepath(filename); gene_go={}; x = 0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x == 0: x = 1 ###skip the first line
        else:
            taxid=t[0];entrez=t[1];goid=t[2]
            if taxid == tax_id:
                try: gene_go[entrez].append(goid)
                except KeyError: gene_go[entrez] = [goid]
    exportEntrezGO(gene_go,species)

def exportEntrezGO(gene_go,species):
    new_file = 'Databases/BuildDBs/BuiltDBs/'+species+'/EntrezGene-GeneOntology.txt'
    today = str(datetime.date.today()); today = string.split(today,'-'); today = today[1]+'/'+today[2]+'/'+today[0]
    y=0
    try: fn=filepath(new_file); data = open(fn,'w')
    except IOError:
        try: new_dir = 'Databases/BuildDBs/BuiltDBs/'+species; fn1 = filepath(new_dir);os.mkdir(fn1) ###Re-Create directory if deleted
        except OSError:
            new_dir = 'Databases/BuildDBs/BuiltDBs'; fn1 = filepath(new_dir);os.mkdir(fn1)
            new_dir = 'Databases/BuildDBs/BuiltDBs/'+species; fn1 = filepath(new_dir);os.mkdir(fn1)
        fn=filepath(new_file); data = open(fn,'w')
    data.write('EntrezGene ID'+'\t'+'GO ID'+'\n')
    for gene in gene_go:
        goids = unique.unique(gene_go[gene])
        for goid in goids: data.write(gene+'\t'+goid+'\n'); y+=1
    data.close()
    print 'Exported',y,'gene-GO relationships for species:',species
    print 
    
def extractAndIntegrateAffyData(species):
    global dirfile; dirfile = unique
    global affy_annotation_db; affy_annotation_db={}; global entrez_annotation_db; entrez_annotation_db = {}
    global ens_to_uid; global entrez_to_uid; global entrez_annotations
    
    import_dir = '/Databases/BuildDBs/Affymetrix/'+species
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for affy_data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = 'Databases/BuildDBs/Affymetrix/'+species+'/'+affy_data
        parse_affymetrix_annotations(affy_data_dir)
        
    import gene_associations
    if incorporate_previous_associations == 'yes':    
        ###dictionary gene to unique array ID
        mod_source1 = 'Ensembl'+'-'+'Affymetrix'; mod_source2 = 'EntrezGene'+'-'+'Affymetrix'
        ens_to_uid = gene_associations.getGeneToUid(species,mod_source1)
        entrez_to_uid = gene_associations.getGeneToUid(species,mod_source2)
        ### Gene IDs with annotation information    
        entrez_annotations = gene_associations.importGeneData(species,'EntrezGene')
        ### If we wish to combine old and new GO relationships - Unclear if this is a good idea
        """if process_go == 'yes':
            entrez_to_go = gene_associations.importGeneGOData(species,'EntrezGene','null')
            ens_to_go = gene_associations.importGeneGOData(species,'Ensembl','null')"""
        integratePreviousAssociations()
    exportRelationshipDBs(species)
    exportResultsSummary(dir_list)

def importAffymetrixAnnotations(dir,Species,Process_go,Extract_go_names,Extract_pathway_names):
    global species; global process_go; global extract_go_names; global extract_pathway_names
    global affy_annotation_db; affy_annotation_db={}

    species = Species; process_go = Process_go; extract_go_names = Extract_go_names; extract_pathway_names = Extract_pathway_names

    import_dir = '/'+dir+'/'+species
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for affy_data in dir_list:    #loop through each file in the directory to output results
        affy_data_dir = dir+'/'+species+'/'+affy_data
        parse_affymetrix_annotations(affy_data_dir)
    return affy_annotation_db

def runBuildAffymetrixAssociations():
    global incorporate_previous_associations; global species; global process_go; global extract_go_names
    global extract_pathway_names
    import_dir = '/Databases/BuildDBs/Affymetrix'
    dir_list = returnDirectories(import_dir)
    import_dir2 = '/Databases/BuildDBs/Entrez/Gene2GO'
    dir_list2 = returnDirectories(import_dir2)
    proceed = 'no'
    extract_go_names = 'no' ### Not currently used
    extract_pathway_names = 'no'

    if len(dir_list)>0 or len(dir_list2)>0:
        while proceed == 'no':
            print "\n*****Build Affymetrix Association Module version 1.0***** \nBuild GO-Elite Relational Database from Affymetrix or Unigene Annotation Files\n"
            print "Program Options"
            print "1) Build uid-gene tables and Entrez gene table from Affymetrix CSV files"
            print "2) Build gene-GO relationship tables for EntrezGene for any species"
            print "3) Quit"
            inp = sys.stdin.readline(); inp = inp.strip()
            if inp  == '1': parse_csv = 'yes';proceed = 'yes'
            elif inp == '2': parse_csv = 'no';proceed = 'yes'
            elif inp == '3': sys.exit()
            else: print "Sorry... that command is not an option\n"
            
        if parse_csv == 'yes':
            proceed = 'no'
            while proceed == 'no':
                print "Choose species directory containing Affymetrix CSV files:"
                x = 1
                for dir in dir_list: print str(x)+')',dir; x+=1
                inp = sys.stdin.readline(); inp = int(inp.strip())
                try: species = dir_list[int(inp)-1]; proceed = 'yes'
                except ValueError:
                    print "Sorry... that command is not an option\n"
                
            proceed = 'no'
            while proceed == 'no':
                print "Data to include in build"
                print "1) Include gene associations from gene and uid-gene files in current directories"
                print "2) Build new database exclusively from new Affymetrix CSV files"
                inp = sys.stdin.readline(); inp = inp.strip()
                if inp  == '1': incorporate_previous_associations = 'yes';proceed = 'yes'
                elif inp == '2': incorporate_previous_associations = 'no';proceed = 'yes'
                else: print "Sorry... that command is not an option\n"

                print "Extracting GO information"
                print "1) Skip this"
                print "2) Extract and infer gene-GeneOntology information from Affymetrix (not ideal - will create a new set of tables, rather than updating existing)"
                inp = sys.stdin.readline(); inp = inp.strip()
                if inp  == '1': process_go = 'no';proceed = 'yes'
                elif inp == '2': process_go = 'yes';proceed = 'yes'
                else: print "Sorry... that command is not an option\n"
                
            affy_annotation_db = extractAndIntegrateAffyData(species)
            print "New databases built...see the folder 'Databases\BuildDBs\BuiltDBs'"
        else:
            proceed = 'no'
            while proceed == 'no':
                print "Choose species to build gene-GO tables using EntrezGene:"
                print "note: you must have already downloaded the gene2go file from NCBI (see Gene2GO dir for URL)"
                print "1) Human"
                print "2) Mouse"
                print "3) Zebrafish"
                print "4) Rat"
                print "5) Other"
                inp = sys.stdin.readline(); inp = inp.strip()
                if inp  == '1': tax_id = '9606'; species = 'Hs';proceed = 'yes'
                elif inp == '2': tax_id = '10090'; species = 'Mm';proceed = 'yes'
                elif inp == '3': tax_id = '7955';species = 'Dr'; proceed = 'yes'
                elif inp == '4': tax_id = '10116';species = 'Rn'; proceed = 'yes'
                elif inp == '5': tax_id = '';proceed = 'yes'
                else: print "Sorry... that command is not an option\n"
            if tax_id == '':
                print "Enter taxid (see: http://www.ncbi.nlm.nih.gov/sites/entrez?db=taxonomy)"
                inp = sys.stdin.readline(); inp = inp.strip()
                tax_id = inp; species = inp
            try: parseGene2GO(tax_id,species)
            except IOError:
                print '\nThe file "gene2go.txt" was not found.\n'
                proceed = 'no'
                while proceed == 'no':
                    print "Would you like this program to fetch this file (~9MB)?"
                    print "1) Yes"
                    print "2) No"
                    inp = sys.stdin.readline(); inp = inp.strip()
                    if inp  == '1' or inp == 'y' or inp == 'Y': download_entrez_go = 'yes';proceed = 'yes'
                    elif inp  == '2' or inp == 'n' or inp == 'N': download_entrez_go = 'no';proceed = 'yes'
                    else: print "Sorry... that command is not an option\n"
                    if download_entrez_go == 'yes':
                        print 'downloading....'
                        url = 'http://conklinwolf.ucsf.edu/informatics/thesis/thesis_Salomonis_6-27-08.pdf'
                        download('ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz','Databases/BuildDBs/Entrez/Gene2GO/','txt')
                        #download(url,'Databases/BuildDBs/Entrez/Gene2GO')
    else: print "No Affymetrix folders or files present. Download these files to a new directory in the folder 'Databases\BuildDBs\Affymetrix'."
    print "Press any key to exit"
    inp = sys.stdin.readline()
if __name__ == '__main__':
    runBuildAffymetrixAssociations()