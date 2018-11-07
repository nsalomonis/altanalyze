from Bio import SeqIO
import sys, string
fasta_file = "/Users/saljh8/GitHub/altanalyze/AltDatabase/EnsMart72/Hs/SequenceData/Homo_sapiens.GRCh37.72.cdna.all.fa" # Input fasta file
result_file = "/Users/saljh8/GitHub/altanalyze/AltDatabase/EnsMart72/Hs/SequenceData/Homo_sapiens.GRCh37.72.cdna.all.filtered.fa" # Output fasta file
            
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        chr = string.split(seq.description,':')[3]
        try:
            float(chr)
            SeqIO.write([seq], f, "fasta")
        except: continue