###ExtractUniProtFunctAnnot
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
import copy

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".dat": dir_list2.append(entry)
    return dir_list2

def importEnsemblUniprot(filename):
    fn=filepath(filename); x = 0
    for line in open(fn,'r').xreadlines():
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            try:
                ensembl=t[0];uniprot=t[1]
                try: uniprot_ensembl_db[uniprot].append(ensembl)
                except KeyError: uniprot_ensembl_db[uniprot] = [ensembl]
            except Exception: null=[] ### Occurs when no file is located on the AltAnalyze server
    print len(uniprot_ensembl_db),"UniProt entries with Ensembl annotations"

def exportEnsemblUniprot(filename):
    import export
    export_data = export.ExportFile(filename)
    export_data.write(string.join(['ensembl','uniprot'],'\t')+'\n')
    for uniprot in uniprot_ensembl_db:
        for ensembl in uniprot_ensembl_db[uniprot]:
            export_data.write(string.join([ensembl,uniprot],'\t')+'\n')
    export_data.close()
    
def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def findSpeciesInUniProtFiles(force):
    ### Download all UniProt annotation files and grab all species names, TaxIDs and corresponding URLs
    
    import AltAnalyze
    ###Get species annotations from the GO-Elite config
    species_annot_db=AltAnalyze.importGOEliteSpeciesInfo(); tax_db={}
    for species_full in species_annot_db:
        taxid=species_annot_db[species_full].TaxID()
        tax_db[taxid]=species_full

    if force == 'yes':
        ### Should only need to be run if UniProt changes it's species to file associations or new species supported by Ensembl
        import export; import update
        filesearch = '_sprot_'
        all_swissprot = update.getFTPData('ftp.expasy.org','/databases/uniprot/current_release/knowledgebase/taxonomic_divisions',filesearch)
        for file in all_swissprot:
            gz_filepath, status = update.download(file,'uniprot_temp/','')        
            if status == 'not-removed':
                try: os.remove(gz_filepath) ### Not sure why this works now and not before
                except OSError: status = status
            
    species_uniprot_db={}; altanalyze_species_uniprot_db={}
    dir=read_directory('/uniprot_temp')
    for filename in dir:    
        fn=filepath('uniprot_temp/'+filename)
        for line in open(fn,'r').xreadlines():
            data = cleanUpLine(line)
            if data[0:2] == 'OX':
                taxid = string.split(data,'=')[1][:-1]
                if taxid in tax_db:
                    species_full = tax_db[taxid]
            elif data[0:2] == 'OS':
                species = data[5:]
                species = string.split(species,' ')[:2]
                species_full = string.join(species,' ')
            elif data[0] == '/':
                url = 'ftp.expasy.org/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/'+filename
                ss = string.split(species_full,' ')
                if len(ss)==2: ### Species name is in the format Homo sapiens - and '(' not in species_full and ')' not in species_full and '/' not in species_full
                    try: species_uniprot_db[species_full].append((taxid,'ftp://'+url+'.gz'))
                    except KeyError: species_uniprot_db[species_full] = [(taxid,'ftp://'+url+'.gz')]
                taxid = ''; species_full = ''
    from build_scripts import EnsemblImport
    species_uniprot_db = EnsemblImport.eliminate_redundant_dict_values(species_uniprot_db)
    ### Export all species to UniProt file relationships so this function needs to only be run once
    import export         
    up = export.ExportFile('Config/uniprot-species-file.txt')
    for species_full in species_uniprot_db:
        values = species_uniprot_db[species_full]
        if len(values)>1:
            found = 'no'
            for (taxid,url) in values:
                if taxid in tax_db:
                    if species_full == tax_db[taxid]: found='yes'; print 'ambiguity resolved:',species_full; break
                if found == 'yes': break
        else: (taxid,url) = values[0]
        up.write(string.join([species_full,taxid,url],'\t')+'\n')
    up.close()
    
def getUniProtURLsForAllSupportedSpecies():
    ### Import all UniProt supproted species and URLs
    species_uniprot_db={}
    fn=filepath('Config/uniprot-species-file.txt')
    for line in open(fn,'r').xreadlines():
        data = cleanUpLine(line)
        species_full,taxid,url = string.split(data,'\t')
        if 'Homo sapiens' not in species_full: ### There's a separate file for us humans (so egotistical!!!)
            species_uniprot_db[species_full] = taxid,url
        
    import AltAnalyze
    ###Get species annotations from the GO-Elite config
    species_annot_db=AltAnalyze.importGOEliteSpeciesInfo()
        
    ### Export all urls for currently supported species
    import UI
    file_location_defaults = UI.importDefaultFileLocations()
    
    location_db={}; species_added=[]
    for species_full in species_annot_db:
        if species_full in species_uniprot_db:
            taxid,url = species_uniprot_db[species_full]
            species_code = species_annot_db[species_full].SpeciesCode()
            try: location_db[url].append(species_code)
            except Exception: location_db[url] = [species_code]
            species_added.append(species_full)
            
    for species_full in species_annot_db:
        taxid = species_annot_db[species_full].TaxID()
        species_code = species_annot_db[species_full].SpeciesCode()
        if species_full not in species_added:
            for species_name in species_uniprot_db:
                tax,url = species_uniprot_db[species_name]
                if tax == taxid:
                    location_db[url].append(species_code)
                    print species_code
                
    for url in location_db:
        species = string.join(location_db[url],'|')
        fl = UI.FileLocationData('ftp', url, species)
        try: file_location_defaults['UniProt'].append(fl)
        except KeyError: file_location_defaults['UniProt'] = [fl]
    UI.exportDefaultFileLocations(file_location_defaults)
    
def import_uniprot_db(filename):
    fn=filepath(filename); global species_not_imported; species_not_imported=[]
    ac = '';sm='';id = '';sq = '';osd = ''; gn = '';dr = '';de = '';ft_string = ''; kw = ''; ft = []; ensembl = []; mgi = []; unigene = []; embl = []
    ft_call=''; rc=''; go=''; x = 0; y = 0; count = 0
    for line in open(fn,'r').xreadlines():
        data, newline= string.split(line,'\n'); 
        #if x<3: print data
        #else: kill
        if 'SRSF1_HUMAN' in data:
            count = 0
        if count == 1: print data
        if data[0:2] == 'ID': id += data[5:]
        elif "GO; GO:" in data: go += data[5:]
        elif data[0:2] == 'DE': de += data[5:]
        elif data[0:2] == 'KW': kw += data[5:] ### Keywords
        elif data[0:2] == 'AC': ac += data[5:]
        elif data[0:2] == 'OS': osd += data[5:]
        elif data[0:2] == 'RC': rc = rc + data[5:]
        elif data[0:2] == '  ': sq += data[5:]
        elif 'DR   Ensembl;' in data:
            null,dr= string.split(data,'Ensembl; '); dr = string.split(dr,'; '); ensembl+=dr[:-1]
        elif 'DR   MGI;' in data: null,dr,null= string.split(data,'; '); mgi.append(dr)
        elif 'DR   UniGene;' in data: null,dr,null= string.split(data,'; '); unigene.append(dr)
        elif 'DR   EMBL;' in data: null,dr,null,null,null= string.split(data,'; '); embl.append(dr)
        elif 'GN   Name=' in data:
            null,gn = string.split(data,'GN   Name='); gn = gn[0:-1]
        elif data[0:2] == 'FT':
            try:
                if len(ft_string) > 0 and data[5] == ' ': ft_string = ft_string + data[33:]
                elif len(ft_string) > 0 and data[5] != ' ': #if previous loop added data but the next ft line is a new piece of functional data
                    ft.append(ft_string) #append the previous value
                    ft_string = data[5:]
                else: ft_string = ft_string + data[5:]
            except IndexError:
                print ft_string;kill
        elif data[0:2] == 'CC': ###grab function description information
            if '-!-' in data: x=0;y=0
            if x == 1: ft_call = ft_call + data[8:]
            if y == 1: sm = sm + data[8:]
            ###if the CC entry is function, begin adding data
            if '-!- FUNCTION:' in data: ft_call = ft_call + data[19:];  x = 1
            if '-!- SIMILARITY' in data: sm = sm + data[21:]; y = 1            
        if data[0] == '/':
            ###Alternatively: if species_name in osd or 'trembl' in filename:
            if count == 1:
                count = 2
            if species_name == 'Mus musculus': alt_osd = 'mouse'
            else: alt_osd = 'alt_osd'
            try:
                if species_name in osd or alt_osd in osd: null=[]
            except TypeError: print species_name,osd,alt_osd;kill
                    
            if species_name in osd or alt_osd in osd:
              class_def,cellular_components = analyzeCommonProteinClassesAndCompartments(sm,kw,ft_call,ft_string,rc,de,go)
              ft_list2 = []; ac = string.split(ac,'; '); ac2=[]
              for i in ac:  i = string.replace(i,';',''); ac2.append(i)
              ac = ac2
              try: id = string.split(id,' ');id = id[0]
              except ValueError: id = id
              sq_str = ''; sq = string.split(sq,' ')
              for entry in sq: sq_str = sq_str + entry; sq = sq_str   
              ft.append(ft_string) #always need to add the current ft_string
              if len(ft_string) > 0: # or len(ft_string) == 0:
                for entry in ft:
                    entry = string.split(entry,'  ')
                    ft_list = []
                    for item in entry:
                        if len(item)>0:
                            try: item = int(item)
                            except ValueError:
                                if item[0] == ' ': item = item[1:]
                                if item[-1] == ' ': item = item[0:-1]
                                else: item = item
                            ft_list.append(item)
                    ft_list2.append(ft_list)
              if 'trembl' in filename: file_type = 'fragment'
              else: file_type = 'swissprot'
              alternate_ensembls=[]
              for secondary_ac in ac:
                  if secondary_ac in uniprot_ensembl_db:
                      for alt_ens in uniprot_ensembl_db[secondary_ac]: alternate_ensembls.append(alt_ens)
                  secondary_to_primary_db[secondary_ac] = id

              ### Update this database with new annotations from the file
              for ens in ensembl:
                for secondary_ac in ac:
                    try:
                        if ens not in uniprot_ensembl_db[secondary_ac]:
                            uniprot_ensembl_db[secondary_ac].append(ens)
                    except KeyError: uniprot_ensembl_db[secondary_ac]=[ens]

              ensembl += alternate_ensembls
              y = UniProtAnnotations(id,ac,sq,ft_list2,ensembl,gn,file_type,de,embl,unigene,mgi,ft_call,class_def,cellular_components)
              uniprot_db[id] = y
            else: species_not_imported.append(osd)
            ac = '';id = '';sq = '';osd = '';gn = '';dr = '';de = ''; ft_call=''; rc='';sm='';go=''; kw=''
            ft_string = '';ft = []; ensembl = []; mgi = []; unigene = []; embl = []
            
            x+=1
    print "Number of imported swissprot entries:", len(uniprot_db)

def analyzeCommonProteinClassesAndCompartments(sm,kw,ft_call,ft_string,rc,de,go):
    ### Used to assign "Common Protein Classes" annotations to Gene Expression summary file (ExpressionOutput folder)
    class_def=[]; annotation=[]; cellular_components = []
    if 'DNA-binding domain' in sm or 'Transcription' in go or 'Transcription regulation' in kw: class_def.append('transcription regulator')
    if 'protein kinase superfamily' in sm or 'Kinase' in go: class_def.append('kinase')
    if 'mRNA splicing' in kw or 'mRNA processing' in kw: class_def.append('splicing regulator')

    if 'G-protein coupled receptor' in sm or 'LU7TM' in sm:
        g_type = []
        if ('adenylate cyclase' in ft_call) or ('adenylyl cyclase'in ft_call):
                ###if both occur
                if (('stimulat' in ft_call) or ('activat' in ft_call)) and ('inhibit' in ft_call):
                    if 'inhibit aden' in ft_call: g_type.append('Gi')
                    if 'stimulate aden' in ft_call or 'activate aden' in ft_call: g_type.append('Gs')
                elif ('stimulat' in ft_call) or ('activat' in ft_call): g_type.append('Gs')
                elif ('inhibit' in ft_call): g_type.append('Gi')
        if ('cAMP' in ft_call):
            if ('stimulat' in ft_call) or ('activat' in ft_call): g_type.append('Gs')
            if ('inhibit' in ft_call): g_type.append('Gi')
        if ('G(s)' in ft_call): g_type.append('Gs')
        if ('G(i)' in ft_call): g_type.append('Gi')
        if ('pertussis' in ft_call and 'insensitive' not in ft_call): g_type.append('Gi')
        if ('G(i/0)' in ft_call) or ('G(i/o)' in ft_call): g_type.append('Gi')
        if ('G(o)' in ft_call): g_type.append('Go')
        if ('G(alpha)q' in ft_call): g_type.append('Gq')
        if ('G(11)' in ft_call): g_type.append('G11')
        if ('G(12)' in ft_call): g_type.append('G12')
        if ('G(13)' in ft_call): g_type.append('G13')
        if ('mobiliz' in ft_call and 'calcium' in ft_call and 'without formation' not in ft_call): g_type.append('Gq')
        if ('phosphatidyl' in ft_call and 'inositol' in ft_call) or ('G(q)' in ft_call) or ('phospholipase C' in ft_call):
                g_type.append('Gq')
        if ('inositol phos' in ft_call) or ('phosphoinositide' in ft_call) or ('PKC'in ft_call) or ('PLC' in ft_call):
            g_type.append('Gq')
        if ('intracellular' in ft_call and 'calcium' in ft_call) and 'nor induced' not in ft_call: g_type.append('Gq')
        if 'G-alpha-11' in ft_call: g_type.append('G11')
        if 'Orphan' in ft_call or 'orphan' in ft_call: g_type.append('orphan')
        if 'taste' in ft_call or 'Taste' in ft_call: g_type.append('taste')
        if 'vision' in ft_call or 'Vision' in ft_call: g_type.append('vision')
        if 'odorant' in ft_call or 'Odorant' in ft_call: g_type.append('oderant')
        if 'influx of extracellar calcium' in ft_call: g_type.append('Gq')
        if 'pheromone receptor' in ft_call or 'Pheromone receptor' in ft_call: g_type.append('pheromone')
        g_protein_list = unique.unique(g_type); g_protein_str = string.join(g_protein_list,'|')
        class_def.append('GPCR(%s)' % g_protein_str)
    elif 'receptor' in sm or 'Receptor' in go: class_def.append('receptor')
    if len(ft_string)>0: ### Add cellular component annotations
        if 'ecreted' in sm: k = 1; annotation.append('extracellular')
        if 'Extraceullar space' in sm: k = 1; annotation.append('extracellular')
        if 'ecreted' in go: k = 1; annotation.append('extracellular')
        if 'xtracellular' in go: k = 1; annotation.append('extracellular')
        if 'Membrane' in sm: k = 1; annotation.append('transmembrane')
        if 'TRANSMEM' in ft_string: k = 1; annotation.append('transmembrane')
        if 'integral to membrane' in go: k = 1; annotation.append('transmembrane')
        if 'Nucleus' in sm: k = 1; annotation.append('nucleus')
        if 'nucleus' in go: k = 1; annotation.append('nucleus')
        if 'Cytoplasm' in sm: k = 1; annotation.append('cytoplasm')
        if 'Mitochondrion' in sm: k = 1; annotation.append('mitochondrion')
        if 'SIGNAL' in ft_string: k = 1; annotation.append('signal')
        ###Generate probably secreted annotations
        if 'signal' in annotation and 'transmembrane' not in annotation:
            for entry in annotation:
                if entry != 'signal': cellular_components.append(entry)
            cellular_components.append('extracellular');annotation = cellular_components
        elif 'signal' in annotation and 'transmembrane' in annotation:
            for entry in annotation:
                if entry != 'signal': cellular_components.append(entry)
            annotation = cellular_components
            
    cellular_components = string.join(annotation,'|')
    class_def = string.join(class_def,'|')
    return class_def, cellular_components

class UniProtAnnotations:
    def __init__(self,primary_id,secondary_ids,sequence,ft_list,ensembl,name,file_type,description,embl,unigene,mgi,ft_call,class_def,cellular_components):
        self._primary_id = primary_id; self._sequence = sequence; self._name = name; self._secondary_ids = secondary_ids; self._class_def = class_def
        self._file_type = file_type; self._description = description; self._ensembl = ensembl; self._ft_list = ft_list
        self._embl = embl; self._unigene = unigene; self._mgi = mgi; self._ft_call = ft_call; self._cellular_components = cellular_components
    def PrimaryID(self): return self._primary_id
    def SecondaryIDs(self): return self._secondary_ids
    def Sequence(self): return self._sequence
    def Name(self): return self._name
    def FTList(self):
        new_FTList = [] ### Transform this set of feature information into objects
        exlcusion_list = ['CHAIN','VARIANT','CONFLICT','VAR_SEQ']
        for ft_entry in self._ft_list:
            try:
                if len(ft_entry)>3: feature, start, stop, description = ft_entry
                else: feature, start, stop = ft_entry; description = ''
                if feature not in exlcusion_list: ### Not informative annotations for AltAnalyze
                    dd = DomainData(feature,start,stop,description)
                    new_FTList.append(dd)
            except ValueError:
                new_FTList = new_FTList ### Occurs when no FT info present
        return new_FTList
    def FileType(self): return self._file_type
    def Description(self): return self._description
    def FunctionDescription(self):
        if 'yrighted' in self._ft_call: return ''
        else: return self._ft_call
    def CellularComponent(self): return self._cellular_components
    def ClassDefinition(self): return self._class_def
    def ReSetPrimarySecondaryType(self,primary_id,secondary_id,file_type):
        secondary_ids = [secondary_id]
        self._primary_id = primary_id
        self._secondary_ids = secondary_ids
        self._file_type = file_type
    def Ensembl(self):
        ens_list = unique.unique(self._ensembl)
        ens_str = string.join(ens_list,',')
        return ens_str
    def EnsemblList(self): return self._ensembl
    def EMBL(self):
        embl_str = string.join(self._embl,',')
        return embl_str
    def Unigene(self):
        unigene_str = string.join(self._unigene,',')
        return unigene_str
    def MGI(self):
        mgi_str = string.join(self._mgi,',')
        return mgi_str
    def DataValues(self):
        output = self.Name()+'|'+self.PrimaryID()
        return output
    def __repr__(self): return self.DataValues()

class DomainData:
    def __init__(self,feature,start,stop,description):
        self._feature = feature; self._start = start; self._stop = stop; self._description = description
    def Feature(self): return self._feature
    def Start(self): return self._start
    def Stop(self): return self._stop
    def Description(self): return self._description
    def DataValues(self):
        output = self.Feature()+'|'+self.Description()
        return output
    def __repr__(self): return self.DataValues()
    
def export():
    fasta_data = uniprot_fildir + 'uniprot_sequence.txt'
    fasta_data2 = uniprot_fildir + 'uniprot_feature_file.txt'
    fn=filepath(fasta_data); fn2=filepath(fasta_data2)
    data = open(fn,'w'); data2 = open(fn2,'w')

    custom_annotations = {}
    for id in uniprot_db:
        y = uniprot_db[id]; ac = ''; ac_list = y.SecondaryIDs(); sq = y.Sequence()
        ft_list = y.FTList(); ft_call = y.FunctionDescription(); ens_list = y.EnsemblList()
        ensembl = y.Ensembl(); mgi = y.MGI();embl = y.EMBL(); unigene = y.Unigene()
        gn = y.Name();de = y.Description()
        file_type = y.FileType(); ac = string.join(ac_list,',')
        
        if '-' in id: ac2 = id; id = ac; ac = ac2
        info = [id,ac,sq,gn,ensembl,de,file_type,unigene,mgi,embl]
        info = string.join(info,'\t')+'\n'
        data.write(info)

        for ens_gene in ens_list:
            if 'T0' not in ens_gene and 'P0' not in ens_gene: ### Exclude protein and transcript IDs
                custom_annot=string.join([ens_gene,y.CellularComponent(), y.ClassDefinition(),gn,de,id,ac,unigene],'\t')+'\n'
                if len(y.CellularComponent())>1 or len(y.ClassDefinition())>1: custom_annotations[ens_gene] = custom_annot
                                     
        if len(ft_list)>0:
            for dd in ft_list:  ### Export domain annotations
                try:
                    n = int(dd.Start()); n = int(dd.Stop()) ### Some now have ?. These are un-informative
                    info2 = string.join([id,ac,dd.Feature(),str(dd.Start()),str(dd.Stop()),dd.Description()],'\t') +'\n'
                    info2 = string.replace(info2,';,','')
                    data2.write(info2)
                except ValueError: null=[]
    data.close();data2.close()

    ###Export custom annotations for Ensembl genes
    output_file = uniprot_fildir + 'custom_annotations.txt'
    fn=filepath(output_file); data = open(fn,'w')
    for entry in custom_annotations: data.write(custom_annotations[entry])
    data.close()
    
def runExtractUniProt(species,species_full,uniprot_filename_url,trembl_filename_url,force):
    global uniprot_ensembl_db;uniprot_ensembl_db={}
    global uniprot_db;uniprot_db={}; global species_name; global uniprot_fildir
    global secondary_to_primary_db; secondary_to_primary_db={}
    import update; reload(update)
    
    species_name = species_full
    
    import UI; species_names = UI.getSpeciesInfo()
    species_full = species_names[species]
    species_full = string.replace(species_full,' ','_')

    uniprot_file = string.split(uniprot_filename_url,'/')[-1]; uniprot_file = string.replace(uniprot_file,'.gz','')
    trembl_file = string.split(trembl_filename_url,'/')[-1]; trembl_file = string.replace(trembl_file,'.gz','')
    uniprot_fildir = 'AltDatabase/uniprot/'+species+'/'
    uniprot_download_fildir = 'AltDatabase/uniprot/'
    uniprot_ens_file = species+'_Ensembl-UniProt.txt'; uniprot_ens_location = uniprot_fildir+uniprot_ens_file
    uniprot_location = uniprot_download_fildir+uniprot_file
    trembl_location = uniprot_download_fildir+trembl_file

    add_trembl_annotations = 'no' ### Currently we don't need these annotations    
    try: importEnsemblUniprot(uniprot_ens_location)
    except IOError:
        try:
            ### Download the data from the AltAnalyze website (if there)
            update.downloadCurrentVersion(uniprot_ens_location,species,'txt')
            importEnsemblUniprot(uniprot_ens_location)
        except Exception: null=[]
    try:
        uniprot_ens_location_built = string.replace(uniprot_ens_location,'UniProt','Uniprot-SWISSPROT')
        uniprot_ens_location_built = string.replace(uniprot_ens_location_built,'uniprot','Uniprot-SWISSPROT')
        importEnsemblUniprot(uniprot_ens_location_built)
    except Exception: null=[]
    
    ### Import UniProt annotations
    counts = update.verifyFile(uniprot_location,'counts')
    if force == 'no' or counts > 8: import_uniprot_db(uniprot_location)
    else:
        ### Directly download the data from UniProt
        gz_filepath, status = update.download(uniprot_filename_url,uniprot_download_fildir,'')

        if status == 'not-removed':
            try: os.remove(gz_filepath) ### Not sure why this works now and not before
            except OSError: status = status     
        import_uniprot_db(uniprot_location)
        
    if add_trembl_annotations == 'yes':
        ### Import TreMBL annotations
        try:
            if force == 'yes': uniprot_location += '!!!!!' ### Force an IOError
            import_uniprot_db(trembl_location)
        except IOError:
            ### Directly download the data from UniProt
            update.download(trembl_filename_url,uniprot_download_fildir,'')
            import_uniprot_db(trembl_location)        
    export()
    exportEnsemblUniprot(uniprot_ens_location)
    
if __name__ == '__main__':
    #findSpeciesInUniProtFiles('no'); sys.exit()
    getUniProtURLsForAllSupportedSpecies()
    sys.exit()
