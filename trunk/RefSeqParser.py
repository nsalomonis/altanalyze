import sys, string
import os.path
import unique

dirfile = unique
py2app_adj = '/AltAnalyze.app/Contents/Resources/Python/site-packages.zip'

def filepath(filename):
    dir=os.path.dirname(dirfile.__file__)       #directory file is input as a variable under the main            
    fn=os.path.join(dir,filename)
    fn = string.replace(fn,py2app_adj,'')
    fn = string.replace(fn,'\\library.zip','') ###py2exe on some systems, searches for all files in the library file, eroneously
    return fn

def read_directory(sub_dir):
    dir=os.path.dirname(dirfile.__file__)
    dir = string.replace(dir,py2app_adj,'')
    dir = string.replace(dir,'\\library.zip','')
    dir_list = os.listdir(dir + sub_dir)
    return dir_list

def import_refseq(filename):
    fn=filepath(filename)
    refseq_mRNA_db = {}
    seq_begin = 0
    y = 0
    for line in open(fn,'r').readlines():
        data, newline= string.split(line,'\n')
        if data[0:5] == 'LOCUS':
            y += 1
            cds = ''
            seq = ''
            ac = data[12:]
            ac_ls= string.split(ac,' ')
            ac = ac_ls[0]
        elif data[5:8] == 'CDS':
            cds = cds + data[21:]
            try:
                cds_start,cds_end = string.split(cds,'..')
            except ValueError:
                if cds[0:4] == 'join':
                   cds,second_cds = string.split(cds[5:],',')
                   cds_start,cds_end = string.split(cds,'..')
            cds = int(cds_start),int(cds_end)
        elif seq_begin == 1:
            seq_temp = data[10:]
            seq_temp = string.split(seq_temp,' ')
            for sequence in seq_temp:
                seq = seq + sequence
        elif data[0:6] == 'ORIGIN':
            seq_begin = 1
        if data[0:2] == '//': #new entry
            if len(seq) > 0 and len(cds)>0 and len(ac)>0:
                #refseq_mRNA_db[ac] = cds,seq
                if  grab_sequence == "all":
                    retrieved_seq = seq
                elif grab_sequence == "3UTR":
                    retrieved_seq = seq[cds[1]:]
                elif grab_sequence == "5UTR":
                    retrieved_seq = seq[0:cds[0]]
                elif grab_sequence == "first_last_coding_30":
                    retrieved_seq = (seq[cds[0]-1:cds[0]-1+30],seq[cds[1]-30:cds[1]])
                if len(retrieved_seq) > 0:
                    refseq_mRNA_db[ac] = retrieved_seq
            ac = ''
            cds = ''
            seq = ''
            seq_begin = 0

    print "Number of imported refseq entries:", len(refseq_mRNA_db),'out of',y
    fasta_data = 'refseq_mRNA_converted_'+grab_sequence+'.txt'
    fn=filepath(fasta_data)
    data = open(fn,'w')
    for ac in refseq_mRNA_db:
        retrieved_seq = refseq_mRNA_db[ac]
        if len(retrieved_seq[0])>0: ###Thus if there are two sequences stored as a tuple or list
            if grab_sequence == "first_last_coding_30":
                retrieved_seq = retrieved_seq[0] +'\t'+ retrieved_seq[1]
        info = ac +'\t'+ str(len(retrieved_seq)) +'\t'+ retrieved_seq +'\n'
        data.write(info)
    data.close()
    print fasta_data, 'written'
    
if __name__ == '__main__':
    dirfile = unique
    
    #grab_sequence = "5UTR"
    #grab_sequence = "3UTR"
    grab_sequence = "all"
    #grab_sequence = "first_last_coding_30"
    input_file = 'mouse.rna.gbff'
    #input_file = 'human.rna.gbff'
    import_refseq(input_file)    
    
