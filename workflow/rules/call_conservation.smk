# conver maf to fasta and run conservation tools 

from snakemake import available_cpu_count
from psutil import virtual_memory

available_mem_gb = lambda: '%dG' % (virtual_memory().available >> 30)
containerized: "docker://apscheben/msa_pipeline:latest"

SPECIES = config['species']
TREE = config['speciesTree']

rule call_conservation:
    input:
       'results/conservation/gerp.complete',
       'results/conservation/phylop.complete'

rule cds_fold:
    input:
       refFasta='data/{refname}.fa'.format(refname=config['refName'])
    output:
       txt=temp('results/tree/sitefold.txt'),
    params:
       gff='data/{ingff}'.format(ingff=config['refGFF'])
    conda:
      '../envs/biopython.yaml'
    log:
      'logs/cds_fold_log.txt'
    threads: 1
    script:
      '../scripts/cds_fold.py'

rule maf2bedfold:
    input:
       'results/tree/sitefold.txt'
    output:
       'results/tree/sitefold.bed'
    log:
      'logs/cds_fold_log.txt'
    threads: 1
    shell:
      """
      tr ' ' '\t' < {input} | awk '$3=="4"' | awk -v OFS="\t" '{{print $1,$2-1,$2}}' > {output}
      """

rule maf2fold:
    input:
       maf='results/roast/roast.maf',
       bed='results/tree/sitefold.bed'
    output:
       maf=temp('results/roast/roast.sitefold.maf')
    conda:
      '../envs/ucsc.yaml'
    log:
      'logs/maf2fold_log.txt'
    threads: 1
    shell:
      """
      mafsInRegion {input.bed} {output.maf} {input.maf} &> {log}
      """

rule clean_maf:
    input:
       maf='results/roast/roast.sitefold.maf'
    output:
       maf=temp('results/roast/roast.sitefold.clean.maf')
    threads: 1
    script:
      '../scripts/clean_mafsInRegion_output.py'

rule maf2fa:
    input:
       'results/roast/roast.maf' if not config['refGFF'] else 'results/roast/roast.sitefold.clean.maf'
    output:
       temp('results/tree/roast.fa')
    params:
       expand(['{species}','{refname}'],species=SPECIES,refname=config['refName'])
    log:
      'logs/maf2fa_log.txt'
    benchmark:
      'benchmark/maf2fa_bm.txt'
    conda:
      '../envs/biopython.yaml'
    threads: 4
    script:
      '../scripts/maf2mfa.py'

rule filtcols:
    input:
       'results/tree/roast.fa'
    output:
       temp('results/tree/roast_mincol.fa')
    log:
      'logs/filtcols_log.txt'
    benchmark:
      'benchmark/filtcols_bm.txt'
    conda:
      '../envs/biopython.yaml'
    threads: 4
    shell:
      """
      sed -i '/^>/! s/N/-/g'  {input} &>> {log}
      trimal -in {input} -out {output} -fasta -gt 0.9 &>> {log}
      """

rule trimcols:
    input:
       'results/tree/roast_mincol.fa'
    output:
       temp('results/tree/roast_mincol_maxlen.fa')
    params:
       maxSites=config["maxNeutralSites"]
    log:
      'logs/filtcols_log.txt'
    benchmark:
      'benchmark/filtcols_bm.txt'
    conda:
      '../envs/biopython.yaml'
    script:
      '../scripts/subsample_sites_in_fasta_alingment.py'

rule raxml:
    input:
       'results/tree/roast_mincol_maxlen.fa'
    output:
       'results/tree/raxml/RAxML_bestTree.roast_mincol_maxlen'
    params:
       workDir='results/tree',
       prefix='roast_mincol_maxlen'
    log:
      '../../logs/raxml_log.txt'
    benchmark:
      'benchmark/raxml_bm.txt'
    conda:
      '../envs/phast.yaml'
    threads: 4
    shell:
      """
      cd {params.workDir}
      raxmlHPC-PTHREADS -T {threads} -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s {params.prefix}.fa -n {params.prefix} &> {log}
      mv RAxML_* raxml/
      """

rule phylofit:
    input:
       tree='results/tree/raxml/RAxML_bestTree.roast_mincol_maxlen',
       alignment='results/tree/roast_mincol_maxlen.fa'
    output:
       'results/tree/neutral.tre'
    params:
       prefix='results/tree/neutral'
    log:
      'logs/phylofit_log.txt'
    benchmark:
      'benchmark/phylofit_bm.txt'
    conda:
      '../envs/phast.yaml'
    threads: 4
    shell:
      """
      phyloFit --out-root {params.prefix} --tree {input.tree}  --subst-mod REV --msa-format FASTA {input.alignment} &>> {log}
      tree_doctor -t {params.prefix}.mod > {output} 2>> {log}
      """ 

checkpoint rule maf_split:
    input:
       maf='results/roast/roast.maf',
       tree='results/tree/neutral.tre'
    output:
       dir=directory('results/conservation_raw')
    log:
      'logs/mafsplit_log.txt'
    benchmark:
      'benchmark/mafsplit_bm.txt'
    conda:
      '../envs/ucsc.yaml'
    threads: 1
    shell:
      """
      mkdir {output.dir}
      mafSplit -byTarget -useFullSequenceName splits.bed {output.dir}/ {input.maf} &> {log}
      """ 

rule gerp:
    input:
       maf='results/conservation_raw/{chr}.maf'
    output:
       temp=('results/conservation_raw/{chr}.maf.rates')
    params:
       refname=config['refName'],
       neutral='results/tree/neutral.tre'
    log:
      'logs/gerp_{chr}_log.txt'
    conda:
      '../envs/phast.yaml'
    threads: 1
    shell:
      """
      gerpcol -t {params.neutral} -f {input.maf} -e {params.refname} -j -z -x '.rates' &> {log}
      """

rule phylop:
    input:
       maf='results/conservation_raw/{chr}.maf'
    output:
       temp=('results/conservation_raw/{chr}.phylop.wig')
    params:
       neutral='results/tree/neutral.mod'
    log:
      'logs/phylop_{chr}_log.txt'
    conda:
      '../envs/phast.yaml'
    threads: 1
    shell:
      """
      phyloP --mode CONACC --method LRT --wig-scores {params.neutral} {input.maf} > {output} 2> {log}
      """

rule wig2bed:
    input:
       'results/conservation_raw/{chr}.phylop.wig'
    output:
       temp=('results/conservation_raw/{chr}.phylop.bed')
    conda:
      '../envs/bedtools.yaml'
    threads: 1
    shell:
      """
      wig2bed < {input} | cut -f1,2,3,5 > {output}
      """

rule add_phylop_header:
    input:
       'results/conservation_raw/{chr}.phylop.cov.bed'
    output:
       'results/phylop/{chr}.phylop.final.bed'
    threads: 1
    shell:
      """
      cat <(printf 'chr\tstart\tend\tphyloP_-log_pvalue\tTaxaAligned\n') {input} > {output}
      """

rule gerp2bed:
    input:
       'results/conservation_raw/{chr}.maf.rates'
    output:
       temp('results/conservation_raw/{chr}.maf.rates.bed')
    conda:
      '../envs/phast.yaml'
    threads: 1
    shell:
      """
      end_pos_bed=`wc -l {input}| cut -f1 -d' '`
      end_pos=$(($end_pos_bed - 1))
      paste <(seq 0 $end_pos) <(seq 1 $end_pos_bed) <(cat {input})| sed "s/^/{wildcards.chr}\t/" > {output}
      """

rule maf2cov:
    input:
       'results/conservation_raw/{chr}.maf'
    output:
       temp('results/conservation_raw/{chr}.maf.cov.bed')
    conda:
      '../envs/biopython.yaml'
    threads: 1
    script:
      '../scripts/maf2bedcov.py'

rule add_cov:
    input:
       cov='results/conservation_raw/{chr}.maf.cov.bed',
       rates='results/conservation_raw/{chr}.maf.rates.bed'
    output:
       temp('results/conservation_raw/{chr}.rates.cov.bed')
    conda:
      '../envs/bedtools.yaml'
    threads: 1
    shell:
      """
      bedtools intersect -loj -a {input.cov} -b {input.rates} | awk -v OFS='\t' '{{print $5,$6,$7,$8,$9,$4}}' > {output} 
      """
rule add_cov_phylop:
    input:
       cov='results/conservation_raw/{chr}.maf.cov.bed',
       bed='results/conservation_raw/{chr}.phylop.bed'
    output:
       temp('results/conservation_raw/{chr}.phylop.cov.bed')
    conda:
      '../envs/bedtools.yaml'
    threads: 1
    shell:
      """
      bedtools intersect -loj -a {input.cov} -b {input.bed} | awk -v OFS='\t' '{{print $5,$6,$7,$8,$4}}' > {output} 
      """

rule add_gerp_header:
    input:
      'results/conservation_raw/{chr}.rates.cov.bed' 
    output:
       'results/gerp/{chr}.rates.cov.final.bed'
    conda:
      '../envs/bedtools.yaml'
    threads: 1
    shell:
      """
      cat <(printf 'chr\tstart\tend\tGERP_ExpSubst\tGERP_RejSubstScore\tTaxaAligned\n') {input} > {output}
      """

def gerp_aggregate(wildcards):
    checkpoint_output = checkpoints.maf_split.get(**wildcards).output[0]
    return expand("results/gerp/{chr}.rates.cov.final.bed",
           chr=glob_wildcards(os.path.join(checkpoint_output, "{chr}.maf")).chr)

def phylop_aggregate(wildcards):
    checkpoint_output = checkpoints.maf_split.get(**wildcards).output[0]
    return expand("results/phylop/{chr}.phylop.final.bed",
           chr=glob_wildcards(os.path.join(checkpoint_output, "{chr}.maf")).chr)

rule gerp_clean:
    input:
       gerp_aggregate
    output:
       temp('results/conservation/gerp.complete')
    conda:
      '../envs/phast.yaml'
    threads: 1
    shell:
      """
      ls {input} > {output}
      """

rule phylop_clean:
    input:
       phylop_aggregate
    output:
       temp('results/conservation/phylop.complete')
    conda:
      '../envs/phast.yaml'
    threads: 1
    shell:
      """
      ls {input} > {output}
      """

