import sys, string
import os.path
import unique
import export
import gene_associations
import traceback
import time

################# Parse directory files

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    #add in code to prevent folder names from being included
    dir_list2 = []
    for file in dir_list:
        lf = string.lower(file)
        if '.txt' in lf or '.sif' in lf or '.tab' in lf: dir_list2.append(file)
    return dir_list2

################# Begin Analysis from parsing files

def getEnsemblGeneData(filename):
    fn=filepath(filename)
    global ensembl_symbol_db; ensembl_symbol_db={}; global symbol_ensembl_db; symbol_ensembl_db={}
    for line in open(fn,'rU').xreadlines():
        data,null = string.split(line,'\n')
        t = string.split(data,'\t')
        ensembl=t[0];symbol=t[1]
        ### Have to do this in order to get the WEIRD chromosomal assocaitions and the normal to the same genes
        try: symbol_ensembl_db[symbol].append(ensembl) 
        except Exception: symbol_ensembl_db[symbol] = [ensembl]
        try: symbol_ensembl_db[string.lower(symbol)].append(ensembl) 
        except Exception: symbol_ensembl_db[string.lower(symbol)] = [ensembl]
        try: symbol_ensembl_db[symbol.title()].append(ensembl) 
        except Exception: symbol_ensembl_db[symbol.title()] = [ensembl]
        ensembl_symbol_db[ensembl] = symbol

def getHMDBData(species):
    program_type,database_dir = unique.whatProgramIsThis()
    filename = database_dir+'/'+species+'/gene/HMDB.txt'

    x=0
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            try: hmdb_id,symbol,description,secondary_id,iupac,cas_number,chebi_id,pubchem_compound_id,Pathways,ProteinNames = t
            except Exception:
                ### Bad Tab introduced from HMDB
                hmdb_id = t[0]; symbol = t[1]; ProteinNames = t[-1]
            symbol_hmdb_db[symbol]=hmdb_id
            hmdb_symbol_db[hmdb_id] = symbol
    
            ProteinNames=string.split(ProteinNames,',')
            ### Add gene-metabolite interactions to databases
            for protein_name in ProteinNames:
                try:
                    for ensembl in symbol_ensembl_db[protein_name]: 
                        z = InteractionInformation(hmdb_id,ensembl,'HMDB','Metabolic')
                        interaction_annotation_dbase[ensembl,hmdb_id] = z ### This is the interaction direction that is appropriate
                        try: interaction_db[hmdb_id][ensembl]=1
                        except KeyError: db = {ensembl:1}; interaction_db[hmdb_id] = db ###weight of 1 (weights currently not-supported)
                        try: interaction_db[ensembl][hmdb_id]=1
                        except KeyError: db = {hmdb_id:1}; interaction_db[ensembl] = db ###weight of 1 (weights currently not-supported)
                except Exception: None
                 
def verifyFile(filename):
    status = 'not found'
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = 'found';break
    except Exception: status = 'not found'
    return status

def importInteractionDatabases(interactionDirs):
    """ Import multiple interaction format file types (designated by the user) """
    exclude=[]
    for file in interactionDirs:
        status = verifyFile(file)
        if status == 'not found':
            exclude.append(file)
    for i in exclude:
        interactionDirs.remove(i)
        
    for fn in interactionDirs:    #loop through each file in the directory to output results
        x=0; imported=0; stored=0
        file = export.findFilename(fn)
        count=0
        print "Parsing interactions from:",file
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            count+=1
            if x==0: x=1
            #elif 'PAZAR' in data or 'Amadeus' in data:x+=0
            else:
                obligatory = False
                imported+=1
                proceed = True
                source=''
                interaction_type = 'interaction'
                try:
                    symbol1,interaction_type, symbol2, ensembl1,ensembl2,source = t
                    ens_ls1=[ensembl1]; ens_ls2=[ensembl2]
                    if 'HMDB' in ensembl1:
                        ensembl1 = string.replace(ensembl1,' ','') ### HMDB ID sometimes proceeded by ' '
                        symbol_hmdb_db[symbol1]=ensembl1
                        hmdb_symbol_db[ensembl1] = symbol1
                        interaction_type = 'Metabolic'
                    if 'HMDB' in ensembl2:
                        ensembl2 = string.replace(ensembl2,' ','') ### HMDB ID sometimes proceeded by ' '
                        symbol_hmdb_db[symbol2]=ensembl2
                        hmdb_symbol_db[ensembl2] = symbol2
                        interaction_type = 'Metabolic'
                except Exception:
                    try:
                        ensembl1,ensembl2,symbol1,symbol2,interaction_type=t
                        if ensembl1 == '':
                            try:
                                ens_ls1 = symbol_ensembl_db[symbol1]
                                ens_ls2 = symbol_ensembl_db[symbol2]
                            except Exception: None
                    except Exception:
                        proceed = False
                if proceed: ### If the interaction data conformed to one of the two above types (typically two valid interacting gene IDs)
                    if (len(ens_ls1)>0 and len(ens_ls2)>0):
                        secondary_proceed = True
                        stored+=1
                        for ensembl1 in ens_ls1:
                            for ensembl2 in ens_ls2:
                                """
                                if (ensembl1,ensembl2) == ('ENSG00000111704','ENSG00000152284'):
                                    print t;sys.exit()
                                if (ensembl1,ensembl2) == ('ENSG00000152284','ENSG00000111704'):
                                    print t;sys.exit()
                                """
                                if 'WikiPathways' in file or 'KEGG' in file:
                                    if ensembl2 != ensembl1:
                                        if (ensembl2,ensembl1) in interaction_annotation_dbase:
                                            del interaction_annotation_dbase[(ensembl2,ensembl1)]
                                            ### Exclude redundant entries with fewer interaction details (e.g., arrow direction BIOGRID) - overwrite with the opposite gene arrangement below
                                        if (ensembl1,ensembl2) in interaction_annotation_dbase:
                                            if interaction_annotation_dbase[(ensembl1,ensembl2)].InteractionType() !='physical':
                                                secondary_proceed = False ### Don't overwrite a more informative annotation like transcriptional regulation or microRNA targeting
                                if 'DrugBank' in fn:
                                    source = 'DrugBank'
                                    interaction_type = 'drugInteraction'
                                    obligatory=True
                                    ensembl1, ensembl2 = ensembl2, ensembl1 ### switch the order of these (drugs reported as first ID and gene as the second)
                                
                                if secondary_proceed:
                                    z = InteractionInformation(ensembl1,ensembl2,source,interaction_type)
                                    interaction_annotation_dbase[ensembl1,ensembl2] = z
                                    #z = InteractionInformation(ensembl2,ensembl1,source,interaction_type)
                                    #interaction_annotation_dbase[ensembl2,ensembl1] = z
                                    try: interaction_db[ensembl1][ensembl2]=1
                                    except KeyError: db = {ensembl2:1}; interaction_db[ensembl1] = db ###weight of 1 (weights currently not-supported)
                                    try: interaction_db[ensembl2][ensembl1]=1
                                    except KeyError: db = {ensembl1:1}; interaction_db[ensembl2] = db ###weight of 1 (weights currently not-supported)
                                
                                if obligatory and source in obligatoryList: ### Include these in the final pathway if linked to any input node (e.g., miRNAs, drugs)
                                    try: obligatory_interactions[ensembl1][ensembl2]=1
                                    except KeyError: db = {ensembl2:1}; obligatory_interactions[ensembl1] = db ###weight of 1 (weights currentlynot-supported)
                                elif source in secondDegreeObligatoryCategories:
                                    try: second_degree_obligatory[ensembl1][ensembl2]=1
                                    except KeyError: db = {ensembl2:1}; second_degree_obligatory[ensembl1] = db ###weight of 1 (weights currently not-supported)
                                    
                else:
                    proceed = False
                    try:
                        ID1, null, ID2 = t
                        proceed = True
                    except Exception:
                        try:
                            ID1, ID2 = t
                            proceed = True
                        except Exception:
                            None
                            
                    if proceed:
                        if 'microRNATargets' in fn:
                            if 'mir' in ID2: prefix = 'MIR'
                            else: prefix = 'LET'
                            ID2='MIR'+string.split(ID2,'-')[2] ### Ensembl naming convention
                            source = 'microRNATargets'
                            interaction_type = 'microRNAInteraction'
                            obligatory=True
                        try: ID_ls1 = symbol_ensembl_db[ID1]
                        except Exception: ID_ls1 = [ID1]
                        try: ID_ls2 = symbol_ensembl_db[ID2]
                        except Exception: ID_ls2 = [ID2]
                        """if 'microRNATargets' in fn:
                            if '*' not in ID2: print ID_ls2;sys.exit()"""
                        addInteractions = True
                        for ID1 in ID_ls1:
                            for ID2 in ID_ls2:
                                z = InteractionInformation(ID2,ID1,source,interaction_type)
                                interaction_annotation_dbase[ID2,ID1] = z ### This is the interaction direction that is appropriate
                                try: interaction_db[ID1][ID2]=1
                                except KeyError: db = {ID2:1}; interaction_db[ID1] = db ###weight of 1 (weights currently supported)
                                try: interaction_db[ID2][ID1]=1
                                except KeyError: db = {ID1:1}; interaction_db[ID2] = db ###weight of 1 (weights currently supported)
                                    
                                if source in secondDegreeObligatoryCategories:
                                    try: second_degree_obligatory[ID1][ID2]=1
                                    except KeyError: db = {ID2:1}; second_degree_obligatory[ID1] = db ###weight of 1 (weights currently supported)
                                    
                                elif obligatory and source in obligatoryList: ### Include these in the final pathway if linked to any input node (e.g., miRNAs, drugs)
                                    try: obligatory_interactions[ID1][ID2]=1
                                    except KeyError: db = {ID2:1}; obligatory_interactions[ID1] = db ###weight of 1 (weights currently supported)
                           
    ### Evaluate the most promiscous interactors (e.g., UBC)
    remove_list=[]
    for ID in interaction_db:
        if len(interaction_db[ID])>2000:
            remove_list.append(ID)
            #print len(interaction_db[ID]),ensembl_symbol_db[ID]
    for ID in remove_list:
        #print 'removing', ID
        del interaction_db[ID]
    blackList[ID] = []

    print 'Imported interactions:',len(interaction_annotation_dbase)

class InteractionInformation:
    def __init__(self, ensembl1, ensembl2, source, interaction_type):
        self._ensembl1 = ensembl1; self._ensembl2 = ensembl2; self._source = source
        self._interaction_type = interaction_type
    def Ensembl1(self): return self._ensembl1
    def Ensembl2(self): return self._ensembl2
    def Source(self): return self._source
    def InteractionType(self): return self._interaction_type
    def Report(self):
        output = self.Ensembl1()+'|'+self.Ensembl2()
        return output
    def __repr__(self): return self.Report()

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def importqueryResults(species,dir_file,id_db):
    global query_db; query_db = {}
    query_interactions={} ### This is the final list of shown interactions
    
    if dir_file == None:
        fileRead = dir_file
    elif '.' in dir_file:
        fn=filepath(dir_file)
        fileRead = open(fn,'rU').xreadlines()
    else:
        fileRead = dir_file ### This is a list of IDs passed to this function rather than in a file
        
    if len(id_db)==0: ### Otherwise, already provided gene IDs to query
        translated=0
        try:
            x=0
            for line in fileRead:
                try:
                    data = cleanUpLine(line)
                    t = string.split(data,'\t')
                except Exception:
                    t = line
                if x==1: x = 1 ### no longer statement since the first row may be a valid ID(s)
                else:
                    id = t[0]
                    ensembl_ls1=[]
                    if id in ensembl_symbol_db:
                        symbol = ensembl_symbol_db[id]
                        query_db[id] = symbol
                        ensembl_ls1 = [id]
                        translated+=1
                    elif id in symbol_ensembl_db:
                        ensembl_ls1 = symbol_ensembl_db[id]
                        translated+=1
                        for ensembl in ensembl_ls1:
                            query_db[ensembl] = id
                    elif id in symbol_hmdb_db:
                            hmdb = symbol_hmdb_db[id]
                            query_db[hmdb] = id
                    elif id in hmdb_symbol_db:
                            symbol = hmdb_symbol_db[id]
                            query_db[id] = symbol
                    else:
                        query_db[id] = id ### Currently not dealt with
                        ensembl_ls1 = [id]
                    
                    ### If a SIF file add genes and interactions
                    if len(t)>1 and 'SIF' in inputDataType: ### Potentially SIF format
                        interaction_type = t[1]
                        try: id2 = t[2]
                        except Exception: id2 = t[1]; interaction_type = 'undetermined'
                        ensembl_ls2=[]
                        if id2 in ensembl_symbol_db:
                            symbol = ensembl_symbol_db[id2]
                            query_db[id2] = symbol
                            ensembl_ls2 = [id2]
                        elif id2 in symbol_ensembl_db:
                            ensembl_ls2 = symbol_ensembl_db[id2]
                            for ensembl in ensembl_ls2:
                                query_db[ensembl] = id2
                        elif id2 in symbol_hmdb_db:
                            hmdb = symbol_hmdb_db[id2]
                            query_db[hmdb] = id2
                        elif id2 in hmdb_symbol_db:
                            symbol = hmdb_symbol_db[id2]
                            query_db[id2] = symbol
                        else:
                            query_db[id2] = id2
                        for ensembl1 in ensembl_ls1:
                            for ensembl2 in ensembl_ls2:
                                try: query_interactions[ensembl1].append(ensembl2)
                                except Exception: query_interactions[ensembl1] = [ensembl2]
                                z = InteractionInformation(ensembl1,ensembl2,'custom',interaction_type)
                                interaction_annotation_dbase[ensembl1,ensembl2] = z          
        except Exception:
            print traceback.format_exc()
            print 'No valid directories or IDs provided. Exiting.'; kill
        if translated==0:
            from visualization_scripts import WikiPathways_webservice
            try: query_db = WikiPathways_webservice.importDataSimple(dir_file,None,MOD='Ensembl',Species=species)[0]
            except Exception: ### If metabolomics
                query_db = WikiPathways_webservice.importDataSimple(dir_file,None,MOD='HMDB',Species=species)[0]
            ### Translate the Ensembl IDs to symbols (where possible)
            for id in query_db:
                if id in ensembl_symbol_db:
                    symbol = ensembl_symbol_db[id]
                else:
                    symbol=id
                query_db[id] = symbol
    else:
        for id in id_db:
            if id_db[id]==None:
               try: id_db[id] = ensembl_symbol_db[id] ### Save symbol (done for imported pathway genes)
               except Exception: id_db[id]=id
        query_db = id_db ### Input gene IDs (not in a file)
        
    print 'Number of IDs from', dir_file, 'is', len(query_db)
    return query_db,query_interactions,dir_file

def associateQueryGenesWithInteractions(query_db,query_interactions,dir_file):
    suffix=''
    if dir_file!=None:
        if len(dir_file)!=0:
            suffix='-'+intNameShort+'_'+export.findFilename(dir_file)[:-4]
    if len(suffix)==0:
        try: suffix = '_'+FileName
        except Exception: None
    file_name = 'AltAnalyze-network'+suffix
    
    query_interactions_unique={}
    interacting_genes={}
    connections = 1
    primary=0
    secondary=0
    terciary=0
    for ensemblGene in query_db:
        if ensemblGene in interaction_db:
            for interacting_ensembl in interaction_db[ensemblGene]:
                if interacting_ensembl not in blackList:
                    ###Only allow direct interactions found in query
                    if interacting_ensembl in query_db:
                        try: query_interactions[ensemblGene].append(interacting_ensembl)
                        except KeyError: query_interactions[ensemblGene] = [interacting_ensembl]
                        try: query_interactions[interacting_ensembl].append(ensemblGene)
                        except KeyError: query_interactions[interacting_ensembl] = [ensemblGene]
                        primary+=1
                    if degrees == 2 or degrees == 'indirect':
                        try: interacting_genes[interacting_ensembl].append(ensemblGene)
                        except KeyError: interacting_genes[interacting_ensembl] = [ensemblGene]
                    elif degrees == 'allInteracting' or degrees == 'all possible':
                        try: query_interactions[ensemblGene].append(interacting_ensembl)
                        except KeyError: query_interactions[ensemblGene] = [interacting_ensembl]
                    if interacting_ensembl in secondaryQueryIDs: ### IDs in the expression file
                        secondary+=1 ### When indirect degrees selected, no additional power added by this (only for direct or shortest path)
                        try: query_interactions[ensemblGene].append(interacting_ensembl)
                        except KeyError: query_interactions[ensemblGene] = [interacting_ensembl]    
        if ensemblGene in second_degree_obligatory:
            for interacting_ensembl in second_degree_obligatory[ensemblGene]:
                try: interacting_genes[interacting_ensembl].append(ensemblGene)
                except KeyError: interacting_genes[interacting_ensembl] = [ensemblGene]

    ### Include indirect interactions to secondaryQueryIDs from the expression file
    if degrees == 2 or degrees == 'indirect':
        for ensemblGene in secondaryQueryIDs:
            if ensemblGene in interaction_db:
                for interacting_ensembl in interaction_db[ensemblGene]:
                    if interacting_ensembl not in blackList:
                        try:
                            interacting_genes[interacting_ensembl].append(ensemblGene)
                            terciary+=1#; print interacting_ensembl
                        except KeyError: None ### Only increase the interacting_genes count if the interacting partner is present from the primary query list
    #print primary,secondary,terciary
    
    ### Report the number of unique interacting genes
    for interacting_ensembl in interacting_genes:
        if len(interacting_genes[interacting_ensembl])==1:
            interacting_genes[interacting_ensembl] = 1
        else:
            unique_interactions = unique.unique(interacting_genes[interacting_ensembl])
            interacting_genes[interacting_ensembl] = len(unique_interactions)
    
    query_indirect_interactions={}; indirect_interacting_gene_list=[]; interacting_gene_list=[]; added=[] 
    if degrees=='shortestPath' or degrees=='shortest path': ### Typically identifying the single smallest path(s) between two nodes.
        query_indirect_interactions, indirect_interacting_gene_list, interacting_gene_list = evaluateShortestPath(query_db,interaction_db,10)
        
    else:
        if degrees==2 or degrees=='indirect' or len(secondDegreeObligatoryCategories)>0:
            for ensembl in interacting_genes:
                if interacting_genes[ensembl] > connections:
                    if ensembl in interaction_db: ### Only nodes removed due to promiscuity will not be found
                        for interacting_ensembl in interaction_db[ensembl]:
                            if interacting_ensembl in query_db or interacting_ensembl in secondaryQueryIDs:
                                try: query_indirect_interactions[interacting_ensembl].append(ensembl)
                                except KeyError: query_indirect_interactions[interacting_ensembl] = [ensembl]
                        ###Record the highest linked nodes
                        indirect_interacting_gene_list.append((interacting_genes[ensembl],ensembl)) 
        if len(obligatory_interactions)>0: ### Include always
            all_reported_genes = combineDBs(query_interactions,query_indirect_interactions) ### combinesDBs and returns a unique list of genes
            for ensemblGene in all_reported_genes: ###This only includes genes in the original input list
                if ensemblGene in obligatory_interactions:
                    for interacting_ensembl in obligatory_interactions[ensemblGene]:
                        #symbol = ensembl_symbol_db[ensemblGene]                    
                        try: query_interactions[ensemblGene].append(interacting_ensembl)
                        except KeyError: query_interactions[ensemblGene] = [interacting_ensembl]
    
    z = dict(query_interactions.items() + query_indirect_interactions.items())
    interaction_restricted_db={}
    for ensembl in z:
        interacting_nodes = z[ensembl]
        for node in interacting_nodes:
            if ensembl in interaction_restricted_db:
                db = interaction_restricted_db[ensembl]
                db[node] = 1
            else: interaction_restricted_db[ensembl] = {node:1}

            if node in interaction_restricted_db:
                db = interaction_restricted_db[node]
                db[ensembl] = 1
            else: interaction_restricted_db[node] = {ensembl:1}
            
    if degrees==2 or degrees=='indirect': ### get rid of non-specific interactions
        query_indirect_interactions, indirect_interacting_gene_list, interacting_gene_list = evaluateShortestPath(query_db,interaction_restricted_db,4)
        
    ###Record the highest linked nodes
    for ensembl in query_interactions:
        linked_nodes = len(unique.unique(query_interactions[ensembl]))
        interacting_gene_list.append((linked_nodes,ensembl))
    interacting_gene_list.sort(); interacting_gene_list.reverse()
    indirect_interacting_gene_list.sort();  indirect_interacting_gene_list.reverse()
    
    print "Length of query_interactions:",len(query_interactions)
    query_interactions_unique=[]
    for gene1 in query_interactions:
        for gene2 in query_interactions[gene1]:
            temp = []; temp.append(gene2); temp.append(gene1)#; temp.sort()
            if gene1 == gene2: interaction_type = 'self'
            else: interaction_type = 'distinct'
            temp.append(interaction_type); temp.reverse()
            query_interactions_unique.append(temp)
    for gene1 in query_indirect_interactions:
        for gene2 in query_indirect_interactions[gene1]:
            temp = []; temp.append(gene2); temp.append(gene1)#; temp.sort()
            if gene1 == gene2: interaction_type = 'self'
            else: interaction_type = 'indirect'
            temp.append(interaction_type); temp.reverse()
            query_interactions_unique.append(temp)
    query_interactions_unique = unique.unique(query_interactions_unique)
    query_interactions_unique.sort()
    

    ###Write out nodes linked to many other nodes
    new_file = outputDir+'/networks/'+file_name+ '-interactions_'+str(degrees)+'_degrees_summary.txt'
    data = export.ExportFile(new_file)
    for (linked_nodes,ensembl) in interacting_gene_list:
        try: symbol = query_db[ensembl]
        except KeyError: symbol = ensembl_symbol_db[ensembl]
        data.write(str(linked_nodes)+'\t'+ensembl+'\t'+symbol+'\t'+'direct'+'\n')
    for (linked_nodes,ensembl) in indirect_interacting_gene_list:
        try: symbol = query_db[ensembl]
        except KeyError:
            try: symbol = ensembl_symbol_db[ensembl]
            except KeyError: symbol = ensembl
            if 'HMDB' in symbol:
                try: symbol = hmdb_symbol_db[ensembl]
                except Exception: pass
        data.write(str(linked_nodes)+'\t'+ensembl+'\t'+symbol+'\t'+'indirect'+'\n')
    data.close()

    regulated_gene_db = query_db    
    sif_export,symbol_pair_unique = exportInteractionData(file_name,query_interactions_unique,regulated_gene_db)
    return sif_export,symbol_pair_unique

def combineDBs(db1,db2):
    ### combinesDBs and returns a unique list of genes
    new_db={}
    for i in db1:
        new_db[i]=[]
        for k in db1[i]:
            new_db[k]=[]
    for i in db2:
        new_db[i]=[]
        for k in db2[i]:
            new_db[k]=[]
    return new_db

def evaluateShortestPath(query_db,interaction_restricted_db,depth):
    interactions_found=0
    start_time = time.time()
    query_indirect_interactions={}; indirect_interacting_gene_list=[]; interacting_gene_list=[]; added=[] 
    print 'Performing shortest path analysis on %s IDs...' % len(query_db),
    for gene1 in query_db:
        for gene2 in query_db:
            if (gene1,gene2) not in added and (gene2,gene1) not in added: 
                if gene1 != gene2 and gene1 in interaction_restricted_db and gene2 in interaction_restricted_db:
                    try:
                        path = shortest_path(interaction_restricted_db,gene1,gene2,depth)
                        added.append((gene1,gene2))
                        i=1
                        while i<len(path): ### Add the relationship pairs
                            try: query_indirect_interactions[path[i-1]].append(path[i])
                            except Exception: query_indirect_interactions[path[i-1]]=[path[i]]
                            interactions_found+=1
                            i+=1
                    except Exception:
                        #tb = traceback.format_exc()
                        pass

        
    if len(query_indirect_interactions)==0:
            print 'None of the query genes interacting in the selected interaction databases...'; queryGeneError
    print interactions_found, 'interactions found in', time.time()-start_time, 'seconds'
    return query_indirect_interactions, indirect_interacting_gene_list, interacting_gene_list

def shortest_path(G, start, end, depth):
   #http://code.activestate.com/recipes/119466-dijkstras-algorithm-for-shortest-paths/
   import heapq
   
   def flatten(L):       # Flatten linked list of form [0,[1,[2,[]]]]
      while len(L) > 0:
         yield L[0]
         L = L[1]

   q = [(0, start, ())]  # Heap of (cost, path_head, path_rest).
   visited = set()       # Visited vertices.
   while True:
      (cost, v1, path) = heapq.heappop(q)
      if v1 not in visited and v1 in G:
         visited.add(v1)
         if v1 == end:
            final_path = list(flatten(path))[::-1] + [v1]
            if len(final_path)<depth:
                return final_path
            else:
                return None
         path = (v1, path)
         for (v2, cost2) in G[v1].iteritems():
            if v2 not in visited:
                heapq.heappush(q, (cost + cost2, v2, path))

def exportInteractionData(file_name,query_interactions_unique,regulated_gene_db):
    file_name = string.replace(file_name,':','-')
    new_file = outputDir+'/networks/'+file_name + '-interactions_'+str(degrees)+'.txt'
    sif_export = outputDir+'/networks/'+file_name + '-interactions_'+str(degrees)+'.sif'
    fn=filepath(new_file); fn2=filepath(sif_export)
    data = open(fn,'w'); data2 = open(fn2,'w')
    added = {} ### Don't add the same entry twice
    symbol_added={}; symbol_pair_unique={}
    for (interaction_type,gene1,gene2) in query_interactions_unique:
        try: symbol1 = query_db[gene1]
        except KeyError:
            try: symbol1 = ensembl_symbol_db[gene1]
            except KeyError: symbol1 = gene1
        if 'HMDB' in symbol1:
            symbol1 = hmdb_symbol_db[gene1]
        try: symbol2 = query_db[gene2]
        except KeyError:
            try: symbol2 = ensembl_symbol_db[gene2]
            except KeyError: symbol2 = gene2
        if 'HMDB' in symbol2:
            symbol2 = hmdb_symbol_db[gene2]
        gene_pair = ''; symbol_pair=''; direction = 'interactsWith'
        if (gene1,gene2) in interaction_annotation_dbase: gene_pair = gene1,gene2; symbol_pair = symbol1,symbol2
        elif (gene2,gene1) in interaction_annotation_dbase: gene_pair = gene2,gene1; symbol_pair = symbol2,symbol1
        else: print gene1, gene2, symbol1, symbol2; kill
        if len(gene_pair)>0:
            y = interaction_annotation_dbase[gene_pair]
            gene1,gene2 = gene_pair ### This is the proper order of the interaction
            symbol1,symbol2 = symbol_pair
            interaction_type = y.InteractionType()
            if interaction_type == 'drugInteraction':
                ### Switch their order
                gene1, gene2, symbol1, symbol2 = gene2, gene1, symbol2, symbol1
            direction = interaction_type
        if (gene_pair,direction) not in added:
            added[(gene_pair,direction)]=[]
            data.write(gene1+'\t'+gene2+'\t'+symbol1+'\t'+symbol2+'\t'+interaction_type+'\n')
            if len(symbol1)>1 and len(symbol2)>1 and (symbol_pair,direction) not in symbol_added:
                if symbol1 != symbol2:
                    data2.write(symbol1+'\t'+direction+'\t'+symbol2+'\n')
                    symbol_added[(symbol_pair,direction)]=[]
                    symbol_pair_unique[symbol_pair]=[]
    data.close(); data2.close()
    print "Interaction data exported"
    return sif_export,symbol_pair_unique

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def importInteractionData(interactionDirs):
    global interaction_db; interaction_db = {}
    global interaction_annotation_dbase; interaction_annotation_dbase = {}
    global obligatory_interactions; obligatory_interactions={}
    global second_degree_obligatory; second_degree_obligatory={}
    global blackList; blackList = {}
    ###Collect both Human and Mouse interactions (Mouse directly sorted in interaction_db
    importInteractionDatabases(interactionDirs)

def interactionPermuteTest(species,Degrees,inputType,inputDir,outputdir,interactionDirs,Genes=None,
                      geneSetType=None,PathwayFilter=None,OntologyID=None,directory=None,expressionFile=None,
                      obligatorySet=None,secondarySet=None,IncludeExpIDs=False):

    global degrees
    global outputDir
    global inputDataType
    global obligatoryList ### Add these if connected to anything
    global secondaryQueryIDs
    global secondDegreeObligatoryCategories ### Add if common to anything in the input - Indicates systems to apply this to
    global symbol_hmdb_db; symbol_hmdb_db={}; global hmdb_symbol_db; hmdb_symbol_db={} ### Create an annotation database for HMDB IDs
    global FileName
    secondaryQueryIDs = {}
    degrees = Degrees
    outputDir = outputdir
    inputDataType = inputType
    obligatoryList = obligatorySet
    secondDegreeObligatoryCategories=[]
    if obligatoryList == None:
        obligatoryList=[]
    if expressionFile == None:
        expressionFile = inputDir ### If it doesn't contain expression values, view as yellow nodes
    if secondarySet!= None and (degrees==1 or degrees=='direct'): ### If degrees == 2, this is redundant
        ### This currently adds alot of predictions - either make more stringent or currently exclude
        secondDegreeObligatoryCategories = secondarySet
    if PathwayFilter != None: FileName = PathwayFilter
    elif OntologyID != None: FileName = OntologyID
    elif Genes != None: FileName = Genes
    
    ### Import Ensembl-Symbol annotations
    getEnsemblGeneData('AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations.txt')
    
    ### Import interaction databases indicated in interactionDirs
    importInteractionData(interactionDirs)
    getHMDBData(species) ### overwrite the symbol annotation from any HMDB that comes from a WikiPathway or KEGG pathway that we also include (for consistent official annotation) 
    input_IDs = getGeneIDs(Genes)
    try: input_IDs = gene_associations.simpleGenePathwayImport(species,geneSetType,PathwayFilter,OntologyID,directory)
    except Exception: None
    permutations = 10000; p = 0
    secondaryQueryIDs = importqueryResults(species,expressionFile,{})[0]
    input_IDs,query_interactions,dir_file = importqueryResults(species,inputDir,input_IDs) ### Get the number of unique genes
    sif_file, original_symbol_pair_unique = associateQueryGenesWithInteractions(input_IDs,query_interactions,dir_file)
    #print len(original_symbol_pair_unique)
    ensembl_unique = map(lambda x: x, ensembl_symbol_db)
    
    interaction_lengths = []
    import random
    while p < permutations:
        random_inputs = random.sample(ensembl_unique,len(input_IDs))
        random_input_db={}
        #print len(random_inputs), len(input_IDs); sys.exit()
        for i in random_inputs: random_input_db[i]=i
        secondaryQueryIDs = importqueryResults(species,random_inputs,{})[0]
        input_IDs,query_interactions,dir_file = importqueryResults(species,inputDir,input_IDs)
        sif_file, symbol_pair_unique = associateQueryGenesWithInteractions(input_IDs,query_interactions,inputDir)
        #print len(symbol_pair_unique);sys.exit()
        interaction_lengths.append(len(symbol_pair_unique))
        p+=1

    interaction_lengths.sort(); interaction_lengths.reverse()
    y = len(original_symbol_pair_unique)
    print 'permuted length distribution:',interaction_lengths
    print 'original length:',y
    k=0
    for i in interaction_lengths:
        if i>=y: k+=1
    
    print 'p-value:',float(k)/float(permutations)

def buildInteractions(species,Degrees,inputType,inputDir,outputdir,interactionDirs,Genes=None,
                      geneSetType=None,PathwayFilter=None,OntologyID=None,directory=None,expressionFile=None,
                      obligatorySet=None,secondarySet=None,IncludeExpIDs=False):
    
    global degrees
    global outputDir
    global inputDataType
    global obligatoryList ### Add these if connected to anything
    global secondaryQueryIDs
    global secondDegreeObligatoryCategories ### Add if common to anything in the input - Indicates systems to apply this to
    global symbol_hmdb_db; symbol_hmdb_db={}; global hmdb_symbol_db; hmdb_symbol_db={} ### Create an annotation database for HMDB IDs
    global FileName
    global intNameShort
    secondaryQueryIDs = {}
    degrees = Degrees
    outputDir = outputdir
    inputDataType = inputType
    obligatoryList = obligatorySet
    secondDegreeObligatoryCategories=[]
    intNameShort=''
    if obligatoryList == None:
        obligatoryList=[]
    if expressionFile == None:
        expressionFile = inputDir ### If it doesn't contain expression values, view as yellow nodes
    if secondarySet!= None and (degrees==1 or degrees=='direct'): ### If degrees == 2, this is redundant
        ### This currently adds alot of predictions - either make more stringent or currently exclude
        secondDegreeObligatoryCategories = secondarySet
    if PathwayFilter != None:
        if len(PathwayFilter)==1:
            FileName = PathwayFilter[0]
        if isinstance(PathwayFilter, tuple) or isinstance(PathwayFilter, list):
            FileName = string.join(list(PathwayFilter),' ')
            FileName = string.replace(FileName,':','-')
        else:
            FileName = PathwayFilter
        if len(FileName)>40:
            FileName = FileName[:40]
    elif OntologyID != None: FileName = OntologyID
    elif Genes != None: FileName = Genes
    
    ### Import Ensembl-Symbol annotations
    getEnsemblGeneData('AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations.txt')
    if len(interactionDirs[0]) == 1: interactionDirs = [interactionDirs]    
    ### Import interaction databases indicated in interactionDirs
    for i in interactionDirs:
        print i
        i = export.findFilename(i)
        i=string.split(i,'-')[1]
        intNameShort+=i[0]

    importInteractionData(interactionDirs)
    getHMDBData(species) ### overwrite the symbol annotation from any HMDB that comes from a WikiPathway or KEGG pathway that we also include (for consistent official annotation) 
    
    input_IDs = getGeneIDs(Genes)
    try:
        if isinstance(PathwayFilter, tuple):
            for pathway in PathwayFilter:
                IDs = gene_associations.simpleGenePathwayImport(species,geneSetType,pathway,OntologyID,directory)
                for id in IDs:input_IDs[id]=None
        else:
            input_IDs = gene_associations.simpleGenePathwayImport(species,geneSetType,PathwayFilter,OntologyID,directory)
    except Exception: None
    if expressionFile == None or len(expressionFile)==0:
        expressionFile = exportSelectedIDs(input_IDs) ### create an expression file
    elif IncludeExpIDs: ### Prioritize selection of IDs for interactions WITH the primary query set (not among expression input IDs)
        secondaryQueryIDs = importqueryResults(species,expressionFile,{})[0]
    input_IDs,query_interactions,dir_file = importqueryResults(species,inputDir,input_IDs)
    sif_file,symbol_pair_unique = associateQueryGenesWithInteractions(input_IDs,query_interactions,dir_file)
    output_filename = exportGraphImage(species,sif_file,expressionFile)
    return output_filename

def exportSelectedIDs(input_IDs):
    expressionFile = outputDir+'/networks/IDList.txt'
    data = export.ExportFile(expressionFile)
    data.write('UID\tSystemCode\n')
    for id in input_IDs:
        if 'HMDB' in id:
            id = hmdb_symbol_db[id]
        data.write(id+'\tEn\n')
    data.close()
    return expressionFile
    
def exportGraphImage(species,sif_file,expressionFile):
    from visualization_scripts import clustering
    output_filename = clustering.buildGraphFromSIF('Ensembl',species,sif_file,expressionFile)
    return output_filename

def getGeneIDs(Genes):
    input_IDs={}
    if Genes == None: None
    elif len(Genes)>0:
        ### Get IDs from list of gene IDs
        Genes=string.replace(Genes,'|',',')
        Genes=string.replace(Genes,' ',',')
        if ',' in Genes: Genes = string.split(Genes,',')
        else: Genes = [Genes]
        for i in Genes:
            if len(i)>0:
                if i in symbol_ensembl_db:
                    for ensembl in symbol_ensembl_db[i]:
                        input_IDs[ensembl]=i ### Translate to Ensembl
                elif i in symbol_hmdb_db:
                    hmdb=symbol_hmdb_db[i]
                    symbol = hmdb_symbol_db[hmdb] ### Get the official symbol
                    input_IDs[hmdb]=symbol ### Translate to HMDB
                else:
                    try: input_IDs[i] = ensembl_symbol_db[i] ### If an input Ensembl ID
                    except Exception: input_IDs[i] = i ### Currently not dealt with
    return input_IDs

if __name__ == '__main__':
    Species = 'Hs'
    Degrees = 2
    inputType = 'IDs'
    inputDir=''
    inputDir='/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/Urine-AR-increased/met/networks/AltAnalyze-network_Met.inceased_AR_1.5fold_metabolite-interactions_shortest path.sif'
    inputDir='/Users/saljh8/Documents/1-dataAnalysis/PaulTang/ARVC_genes.txt'
    obligatorySet = []#['drugInteraction']#'microRNAInteraction'
    Genes = 'POU5F1,NANOG,TCF7L1,WNT1,CTNNB1,SOX2,TCF4,GSK3B'
    Genes = 'Glucose'; Degrees = 'shortestPath'; Degrees = 'indirect'; Degrees = 'all possible'
    Genes = ''; Degrees='indirect'
    interactionDirs = []
    Genes=''
    outputdir = filepath('AltAnalyze/test')
    outputdir = '/Users/saljh8/Desktop/Archived/Documents/1-manuscripts/Salomonis/SIDS-WikiPathways/Interactomics/'
    interaction_root = 'AltDatabase/goelite/'+Species+'/gene-interactions'
    files = read_directory('AltDatabase/goelite/'+Species+'/gene-interactions')
    rooot = '/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/CTOTC/AltAnalyze Based/GO-Elite/MarkerFinder/'
    expressionFile=None
    expressionFile = '/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/Urine-AR-increased/UrinProteomics_Kidney-All/GO-Elite/input/GE.AR_vs_STA-fold1.5_rawp0.05.txt'
    expressionFile = '/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/BKVN infection/GO-Elite/input/AR_vs_norm_adjp05.txt'
    expressionFile = '/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/Blood AR-BK/AR-STA/Batches/overlap/AR_vs_STA_p0.05_fold1_common.txt'
    expressionFile=None
    #files2 = read_directory(rooot)
    #inputType = 'SIF'
    for file in files:
        if 'micro' not in file and 'all-Drug' not in file and 'GRID' not in file and 'Drug' not in file and 'TF' not in file: # and 'TF' not in file and 'KEGG' not in file:
            interactionDirs.append(filepath(interaction_root+'/'+file))
    #"""
    inputDir='/Users/saljh8/Desktop/Archived/Documents/1-manuscripts/Salomonis/SIDS-WikiPathways/Interactomics/CoreGeneSet67/core_SIDS.txt'
    expressionFile = '/Users/saljh8/Desktop/Archived/Documents/1-manuscripts/Salomonis/SIDS-WikiPathways/Interactomics/Proteomics/proteomics_kinney.txt'
    
    interactionPermuteTest(Species,Degrees,inputType,inputDir,outputdir,interactionDirs,Genes=Genes,obligatorySet=obligatorySet,expressionFile=expressionFile, IncludeExpIDs=True)
    sys.exit()
    buildInteractions(Species,Degrees,inputType,inputDir,outputdir,interactionDirs,Genes=Genes,obligatorySet=obligatorySet,expressionFile=expressionFile, IncludeExpIDs=True)
    sys.exit()
    #"""
    #canonical Wnt signaling: GO:0060070
    # BioMarkers 'Pluripotent Stem Cells' 'gene-mapp'
    #inputDir = '/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/Diabetes-Blood/ACR/log2/MergedFiles-Symbol_ACR.txt'
    #inputDir = '/Users/nsalomonis/Desktop/dataAnalysis/SplicingFactors/RBM20_splicing_network.txt'; inputType = 'SIF'
    #inputDir = '/Users/nsalomonis/Documents/1-manuscripts/Salomonis/SIDS-WikiPathways/67_SIDS-genes.txt'
    #Genes=None
    #exportGraphImage(Species,'/Users/nsalomonis/Desktop/AltAnalyze/AltAnalyze/test/networks/AltAnalyze-network-interactions_1degrees.sif',inputDir);sys.exit()
    #buildInteractions(Species,Degrees,inputType,inputDir,outputdir,interactionDirs,Genes=None,obligatorySet=obligatorySet,geneSetType='BioMarkers',PathwayFilter='Pluripotent Stem Cells',directory='gene-mapp')
    buildInteractions(Species,Degrees,inputType,inputDir,outputdir,interactionDirs,Genes=Genes,obligatorySet=obligatorySet,expressionFile=expressionFile)

