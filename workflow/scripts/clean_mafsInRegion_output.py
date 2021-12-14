import sys
# Each target species occurs once in the output alignment
# Paralogous bases (when a species occurs more than once in an alignment block) are masked as 'N'
def clean_maf(inmaf,outmaf):
    with open(outmaf,'a') as out:
        with open(inmaf,'r') as maf:
            for line in maf:
                #line = line.strip()
                if line.startswith("#"):
                    print(line,file=out)
                elif not line.strip():
                    pass
                else:
                    line = line.strip()
                    line = line.split(" ")
                    if line[1].startswith("score="):
                        print("",file=out)
                        print(" ".join(line),file=out)
                    else:
                        # truncate seq to single base
                        seq = line[-1]
                        #print(seq)
                        if seq == "(null)":
                            seq = "N"
                        else:
                            seq = seq[0]
                        line[-1] = seq
                        line = " ".join(line)
                        print(line,file=out)

clean_maf(snakemake.input[0],snakemake.output[0]) 
