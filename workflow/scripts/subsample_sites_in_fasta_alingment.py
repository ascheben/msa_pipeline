import sys
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def trim_fa(infile,max_len,outfa):
    infa = SeqIO.parse(infile, "fasta")
    # if iterating infa object, record will be skipped, thus use new FastaIterator
    seq_len = len(next(SeqIO.parse(infile, "fasta")))
    seq_cut = []
    if seq_len > max_len:
        # sample random indices
        for alignment in infa:
            # Get max_len seq from start of alignment
            subseq = alignment.seq[0:max_len]
            # Sampling random sites or blocks substantially slows down process
            record = SeqRecord(subseq,alignment.id,"","")
            seq_cut.append(record)
    else:
        seq_cut = infa
    SeqIO.write(seq_cut, outfa, "fasta")

trim_fa(snakemake.input[0],snakemake.config["maxNeutralSites"],snakemake.output[0])
