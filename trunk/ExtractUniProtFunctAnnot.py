###ExtractUniProtFunctAnnot
#Copyright 2005-2008 J. Davide Gladstone Institutes, San Francisco California
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

import sys, string
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
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

def import_ensembl_uniprot_db(filename):
    fn=filepath(filename)
    for line in open(fn,'r').xreadlines():
        data, newline= string.split(line,'\n')
        t = string.split(data,'\t')
        ensembl,uniprot,null = t
        try: uniprot_ensembl_db[uniprot].append(ensembl)
        except KeyError: uniprot_ensembl_db[uniprot] = [ensembl]

def import_uniprot_db(filename):
    fn=filepath(filename); global species_not_imported; species_not_imported=[]
    ac = '';id = '';sq = '';osd = ''; gn = '';dr = '';de = '';ft_string = '';ft = []; ensembl = []; mgi = []; unigene = []; embl = []
    ft_call=''; rc=''; x = 0; y = 0
    for line in open(fn,'r').xreadlines():
        data, newline= string.split(line,'\n')
        #if x<3: print data
        #else: kill
        if data[0:2] == 'ID': id = id + data[5:]
        elif data[0:2] == 'DE': de += data[5:]
        elif data[0:2] == 'AC': ac += data[5:]
        elif data[0:2] == 'OS': osd += data[5:]
        elif data[0:2] == 'RC': rc = rc + data[5:]
        elif data[0:2] == '  ': sq += data[5:]
        elif 'DR   Ensembl;' in data: null,dr,null= string.split(data,'; '); ensembl.append(dr)
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
            ###if the CC entry is function, begin adding data
            if '-!- FUNCTION:' in data: ft_call = ft_call + data[19:];  x = 1              
        if data[0] == '/':
            ###Alternatively: if species_name in osd or 'trembl' in filename:
            if species_name == 'Mus musculus': alt_osd = 'mouse'
            else: alt_osd = 'alt_osd'
            try:
                if species_name in osd or alt_osd in osd: null=[]
            except TypeError: print species_name,osd,alt_osd;kill
                    
            if species_name in osd or alt_osd in osd:
              ft_list2 = []
              ac = string.split(ac,'; '); ac2=[]
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
              ensembl += alternate_ensembls
              y = UniProtAnnotations(id,ac,sq,ft_list2,ensembl,gn,file_type,de,embl,unigene,mgi,ft_call)
              uniprot_db[id] = y
            else: species_not_imported.append(osd)
            ac = '';id = '';sq = '';osd = '';gn = '';dr = '';de = ''; ft_call=''; rc=''
            ft_string = '';ft = []; ensembl = []; mgi = []; unigene = []; embl = []
            
            x+=1
    print "Number of imported swissprot entries:", len(uniprot_db)

class UniProtAnnotations:
    def __init__(self,primary_id,secondary_ids,sequence,ft_list,ensembl,name,file_type,description,embl,unigene,mgi,ft_call):
        self._primary_id = primary_id; self._sequence = sequence; self._name = name; self._secondary_ids = secondary_ids;
        self._file_type = file_type; self._description = description; self._ensembl = ensembl; self._ft_list = ft_list
        self._embl = embl; self._unigene = unigene; self._mgi = mgi; self._ft_call = ft_call
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
    def ReSetPrimarySecondaryType(self,primary_id,secondary_id,file_type):
        secondary_ids = [secondary_id]
        self._primary_id = primary_id
        self._secondary_ids = secondary_ids
        self._file_type = file_type
    def Ensembl(self):
        ens_list = unique.unique(self._ensembl)
        ens_str = string.join(ens_list,',')
        return ens_str
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
    fn=filepath(fasta_data)
    fn2=filepath(fasta_data2)
    data = open(fn,'w')
    data2 = open(fn2,'w')

    for id in uniprot_db:
        y = uniprot_db[id]
        ac = ''
        ac_list = y.SecondaryIDs()
        sq = y.Sequence()
        ft_list = y.FTList(); ft_call = y.FunctionDescription()
        ensembl = y.Ensembl()
        mgi = y.MGI();embl = y.EMBL(); unigene = y.Unigene()
        gn = y.Name();de = y.Description()
        file_type = y.FileType()
        ac = string.join(ac_list,',')
        
        if '-' in id: ac2 = id; id = ac; ac = ac2
        info = [id,ac,sq,gn,ensembl,de,file_type,unigene,mgi,embl]
        info = string.join(info,'\t')+'\n'
        data.write(info)
        
        if len(ft_list)>0:
            for dd in ft_list:  ### Export domain annotations
                try:
                    n = int(dd.Start()); n = int(dd.Stop()) ### Some now have ?. These are un-informative
                    info2 = string.join([id,ac,dd.Feature(),str(dd.Start()),str(dd.Stop()),dd.Description()],'\t') +'\n'
                    info2 = string.replace(info2,';,','')
                    data2.write(info2)
                except ValueError: null=[]
                
    data.close();data2.close()

def runExtractUniProt(species,species_full,uniprot_filename_url,trembl_filename_url,force):
    global uniprot_ensembl_db;uniprot_ensembl_db={}
    global uniprot_db;uniprot_db={}; global species_name; global uniprot_fildir
    global secondary_to_primary_db; secondary_to_primary_db={}
    
    species_name = species_full
    uniprot_file = string.split(uniprot_filename_url,'/')[-1]; uniprot_file = string.replace(uniprot_file,'.gz','')
    trembl_file = string.split(trembl_filename_url,'/')[-1]; trembl_file = string.replace(trembl_file,'.gz','')
    uniprot_fildir = 'AltDatabase/uniprot/'+species+'/'
    uniprot_ens_file = species+'_Ensembl-UniProt.txt'; uniprot_ens_location = uniprot_fildir+uniprot_ens_file
    uniprot_location = uniprot_fildir+uniprot_file
    trembl_location = uniprot_fildir+trembl_file
    try: import_ensembl_uniprot_db(uniprot_ens_location)
    except IOError:
        import update; reload(update)
        ### Download the data from the AltAnalyze website
        update.downloadCurrentVersion(uniprot_ens_file,uniprot_fildir,species,'txt')
        import_ensembl_uniprot_db(uniprot_ens_location)
    ### Import UniProt annotations
    try:
        if force == 'yes': uniprot_location += '!!!!!' ### Force an IOError
        import_uniprot_db(uniprot_location)
    except IOError:
        import update; reload(update)
        ### Directly download the data from UniProt
        update.download(uniprot_filename_url,uniprot_fildir,'')
        import_uniprot_db(uniprot_location)
    ### Import TreMBL annotations
    try:
        if force == 'yes': uniprot_location += '!!!!!' ### Force an IOError
        import_uniprot_db(trembl_location)
    except IOError:
        import update; reload(update)
        ### Directly download the data from UniProt
        update.download(trembl_filename_url,uniprot_fildir,'')
        import_uniprot_db(trembl_location)        
    export()
    
if __name__ == '__main__':
    sys.exit()