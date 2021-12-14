import sys
from Bio import AlignIO

def mafcoverage(inmaf,outbed):
    with open(outbed,'a') as out:
        for multiple_alignment in AlignIO.parse(inmaf, "maf"):
            # number species in alignment
            rows = len(multiple_alignment)
            # Note this will count positions with "-" or "N" characters as aligned
            record = multiple_alignment[0]
            chrom = ".".join(str(record.id).split(".")[1:])
            start = int(record.annotations["start"])
            size = int(record.annotations["size"])
            print(*[chrom,start,start+size,rows],sep="\t",file=out)

mafcoverage(snakemake.input[0], snakemake.output[0])
