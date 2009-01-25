import sys, string
import os.path
import unique
import copy

dirfile = unique

def filepath(filename):
    dir=os.path.dirname(dirfile.__file__)       #directory file is input as a variable under the main            
    fn=os.path.join(dir,filename)
    return fn

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
    ac = '';id = '';sq = '';os = '';gn = '';dr = '';de = '';ft_string = '';ft = []; ensembl = []; mgi = []; unigene = []; embl = []
    ft_call=''; sm=''; rc=''; x = 0; y = 0
    for line in open(fn,'r').xreadlines():
        data, newline= string.split(line,'\n')
        #if x<3: print data
        #else: kill
        if data[0:2] == 'ID': id = id + data[5:]
        elif data[0:2] == 'DE': de += data[5:]
        elif data[0:2] == 'AC': ac += data[5:]
        elif data[0:2] == 'OS': os += data[5:]
        elif data[0:2] == 'RC': rc = rc + data[5:]
        elif data[0:2] == '  ': sq += data[5:]
        elif 'DR   Ensembl;' in data: null,dr,null= string.split(data,'; '); ensembl.append(dr)
        elif 'DR   MGI;' in data: null,dr,null= string.split(data,'; '); mgi.append(dr)
        elif 'DR   UniGene;' in data: null,dr,null= string.split(data,'; '); unigene.append(dr)
        elif 'DR   EMBL;' in data: null,dr,null,null,null= string.split(data,'; '); embl.append(dr)
        elif 'GN   Name=' in data:
            null,gn = string.split(data,'GN   Name=')
            gn = gn[0:-1]
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
            ###Alternatively: if species_type in os or 'trembl' in filename:
            if species_type == 'Mus musculus': alt_os = 'mouse'
            else: alt_os = 'alt_os'
            if species_type in os or alt_os in os:
              if 'G-protein coupled receptor' in sm: tissue_type,g_type = gpcrNLP(sm,ft_call,rc,de); gpcr_status = 'GPCR'
              else: tissue_type=[];g_type=[]; gpcr_status = ''
              ft_list2 = []
              ac = string.split(ac,'; ')
              ac2=[]
              for i in ac:  i = string.replace(i,';',''); ac2.append(i)
              ac = ac2
              try: id = string.split(id,' ');id = id[0]
              except ValueError: id = id
              sq_str = ''
              sq = string.split(sq,' ')
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
              y = UniProtAnnotations(id,ac,sq,ft_list2,ensembl,gn,file_type,de,embl,unigene,mgi,ft_call,tissue_type,g_type,gpcr_status)
              uniprot_db[id] = y
              ###Create a database of seconary IDs (really the primary used accession) to the UniProt official name (id) for linkin the splice-var data to annotations
            else: species_not_imported.append(os)
            ac = '';id = '';sq = '';os = '';gn = '';dr = '';de = ''; ft_call=''; sm=''; rc=''
            ft_string = '';ft = []; ensembl = []; mgi = []; unigene = []; embl = []; tissue_type=[];g_type=[]
            
            x+=1
    print "Number of imported swissprot entries:", len(uniprot_db)

def gpcrNLP(sm,ft_call,rc,de):
    g_type = []; tissue_type=[]
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
    ###Define Tissue type as neuronal if applicable
    if ('synapt' in ft_call) or ('brain' in ft_call) or ('neuro' in ft_call) or ('cortical' in ft_call):
        tissue_type.append('brain')
    if ('cerebellum' in ft_call) or ('neocortex' in ft_call) or ('olfactory' in ft_call):
        tissue_type.append('brain')
    if ('brain' in rc) or ('brain' in de): tissue_type.append('brain')
    return tissue_type,g_type

class UniProtAnnotations:
    def __init__(self,primary_id,secondary_ids,sequence,ft_list,ensembl,name,file_type,description,embl,unigene,mgi,ft_call,tissue_list,g_protein_list,gpcr_status):
        self._primary_id = primary_id; self._sequence = sequence; self._name = name; self._secondary_ids = secondary_ids; self._g_protein_list = g_protein_list
        self._file_type = file_type; self._description = description; self._ensembl = ensembl; self._ft_list = ft_list; self._tissue_list = tissue_list
        self._embl = embl; self._unigene = unigene; self._mgi = mgi; self._ft_call = ft_call; self._gpcr_status = gpcr_status
    def PrimaryID(self): return self._primary_id
    def SecondaryIDs(self): return self._secondary_ids
    def Sequence(self): return self._sequence
    def Name(self): return self._name
    def FTList(self): return self._ft_list
    def FileType(self): return self._file_type
    def Description(self): return self._description
    def FunctionDescription(self):
        if 'yrighted' in self._ft_call: return ''
        else: return self._ft_call
    def TissueTypes(self):
        tissue_list = unique.unique(self._tissue_type); tissue_str = string.join(tissue_list,',')
        return tissue_str
    def GProtein(self):
        g_protein_list = unique.unique(self._g_protein_list); g_protein_str = string.join(g_protein_list,',')
        return g_protein_str
    def GPCRStatus(self): return self._gpcr_status
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
    
def import_fasta(filename):
    print "Begining generic fasta import of",filename
    #>gi|10048413|ref|NM_020493.1| Mus musculus serum response factor (Srf), mRNA
    #>P48347-2 (P48347) Splice isoform 2 of P48347
    #AGCCAGCGCGTGGTCCCGGCCCCCTCCACCCGCGGTCTCGGCCGCGGCCAGCAGCCCCTGCCCCGCGGGG
    fn=filepath(filename)
    sequence = '|'; description = ''; x=0
    for line in open(fn,'r').xreadlines():
        try: data, newline= string.split(line,'\n')
        except ValueError: continue
        if data[0] == '>':
                #print len(sequence)
                if len(sequence) > 1:
                    if sp_id in secondary_to_primary_db:
                        primary_ac = secondary_to_primary_db[sp_id]
                        y = copy.deepcopy(uniprot_db[primary_ac])
                        sp_id = string.replace(sp_id,'>','')
                        variant_id = string.replace(variant_id,'>','')
                        y.ReSetPrimarySecondaryType(variant_id,sp_id,'variant')
                        uniprot_db[variant_id] = y; x+=1
                    sequence = '|';description = ''
                    data_list = string.split(data,' ')
                    variant_id = data_list[0]; sp_id = data_list[1][1:-1]
                    for entry in data_list[2:]: description = description + ' ' + entry
                else:                             #only used for the first entry, then the  above if is always used
                    data_list = string.split(data,' ')
                    variant_id = data_list[0]; sp_id = data_list[1][1:-1]
                    for entry in data_list[2:]: description = description + ' ' + entry
        try:
            if data[0] != '>': sequence = sequence + data
        except IndexError: continue

    print "Number of imported sequences:", x 

def export():
    fasta_data = species_type+'.uniprot_db.txt'
    fasta_data2 = species_type+'.feature_file.txt'
    fasta_data3 = species_type+'.GPCR.txt'
    fn=filepath(fasta_data)
    fn2=filepath(fasta_data2)
    fn3=filepath(fasta_data3)
    data = open(fn,'w')
    data2 = open(fn2,'w')
    data3 = open(fn3,'w')    

    for id in uniprot_db:
        y = uniprot_db[id]
        ac = ''; ac_list = y.SecondaryIDs();sq = y.Sequence(); ft_list = y.FTList()
        ensembl = y.Ensembl(); mgi = y.MGI();embl = y.EMBL(); unigene = y.Unigene()
        gn = y.Name(); file_type = y.FileType();de = y.Description()
        ac = string.join(ac_list,','); ft_call = y.FunctionDescription(); g_type = y.GProtein()
        gpcr_status = y.GPCRStatus()
        if '-' in id: ac2 = id; id = ac; ac = ac2
        info = [id,ac,sq,gn,ensembl,de,file_type,unigene,mgi,embl]
        info3 = [id,ac,gn,ensembl,de,file_type,unigene,mgi,embl,ft_call,g_type,gpcr_status]
        try: info = string.join(info,'\t')+'\n'
        except TypeError: print info;kill
        info3 = string.join(info3,'\t')+'\n'
        data.write(info); data3.write(info3)
        if len(ft_list)>0 and file_type!='variant':
            for ft_entry in ft_list:
                temp = id +'\t'+ ac
                for item in ft_entry: temp = temp +'\t'+ str(item)
                info2 = temp + '\n'
                info2 = string.replace(info2,';,','')
                data2.write(info2)
                
    data.close()
    data2.close()
    data3.close()
    
if __name__ == '__main__':
    #kill
    global uniprot_ensembl_db;uniprot_ensembl_db={}
    global uniprot_db;uniprot_db={}
    global secondary_to_primary_db; secondary_to_primary_db={}
    dirfile = unique
    species_type = 'Mus musculus'
    #species_type = 'Homo sapiens'

    if species_type == 'Mus musculus':
        input_file = 'Mm_Ensembl-UniProt.txt'; import_ensembl_uniprot_db(input_file)
        input_file = 'uniprot_sprot_rodents.dat'; import_uniprot_db(input_file)
        input_file = 'uniprot_trembl_rodents.dat'; import_uniprot_db(input_file)
        input_file = 'uniprot_sprot_varsplic.fasta.txt'; import_fasta(input_file)

    if species_type == 'Homo sapiens':
        input_file = 'Hs_Ensembl-UniProt.txt'; import_ensembl_uniprot_db(input_file)
        input_file = 'uniprot_sprot_human.dat'; import_uniprot_db(input_file)
        input_file = 'uniprot_trembl_human.dat'; import_uniprot_db(input_file)
        input_file = 'uniprot_sprot_varsplic.fasta.txt'; import_fasta(input_file)
    export()
