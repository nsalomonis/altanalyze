###RefSeqParser
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

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def import_refseq(filename):
    fn=filepath(filename)
    refseq_mRNA_db = {}
    seq_begin = 0
    y = 0
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
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
    
    #grab_sequence = "5UTR"
    #grab_sequence = "3UTR"
    grab_sequence = "all"
    #grab_sequence = "first_last_coding_30"
    input_file = 'mouse.rna.gbff'
    #input_file = 'human.rna.gbff'
    import_refseq(input_file)    
    
