import sys
import itertools
from Bio import AlignIO

def maf2fa(inmaf,species,outfa):
    maf = AlignIO.parse(inmaf, "maf")
    # get list of species
    alndict = {}
    for line in species:
        line = line.strip()
        #line = line.split("\t")
        if line not in alndict:
            #spdict[line] = []
            alndict[line] = []
    num_sp = len(alndict)
    count = 0
    #print(spdict,csdict)
    for multiple_alignment in maf:
        spcheck = []
        aln_list = []
        # MAF format requires all sequences to be same length
        aln_len = len(multiple_alignment[0].seq)
        for record in multiple_alignment:
            #assume first record is ref
            chrom = str(record.id)
            species = chrom.split(".")[0]
            if species not in spcheck:
                spcheck.append(species)
                alndict[species].append(str(record.seq.upper()))
            else:
                # Replace duplicated entries with gaps
                alndict[species][-1] = "-" * aln_len
        for sp,seqs in alndict.items():
            if sp not in spcheck:
                alndict[sp].append("-" * aln_len)
    with open(outfa, 'w') as outfile:
        for sp,seq in alndict.items():
            outfile.write(">" + sp + "\n")
            seq = "".join(seq)
            outfile.write(seq + "\n")

maf2fa(snakemake.input[0],snakemake.params[0],snakemake.output[0])
