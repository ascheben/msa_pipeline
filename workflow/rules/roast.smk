# This series of commands will run multiz-tba "roast" program to create
# a combined MAF file that can be used as input to GERP++ or other analysis
# programs.

from snakemake import available_cpu_count
from psutil import virtual_memory

available_mem_gb = lambda: '%dG' % (virtual_memory().available >> 30)

SPECIES = config['species']
TREE = config['speciesTree']

if config["speciesTree"]:
    ruleorder: write_tree > make_tree
else:
    ruleorder: make_tree > write_tree

rule roast:
    input:
      roastMerge='results/roast/roast.maf'


rule write_tree:
    output:
       'results/tree/topology.nwk'
    params:
       tree=config["speciesTree"]
    threads: 1
    shell:
       """
       echo '{params.tree}' > {output} 
       """

rule make_tree:
    input:
       expand(['data/{species}.fa','data/{refname}.fa'],species=SPECIES,refname=config['refName'])
    output:
       mtree='results/tree/topology.nwk'
    log:
      'logs/topology_log.txt'
    benchmark:
      'benchmark/topology_bm.txt'
    conda:
      '../envs/mashtree.yaml'
    threads: 4
    shell:
      """ 
      mashtree --sketch-size 20000 --genomesize 1000000000 --numcpus {threads} --mindepth 0 {input} 2>>{log} > {output.mtree}
      """

rule clean_tree:
    input:
      top='results/tree/topology.nwk'
    output:
       mtree='results/tree/topology_clean.txt'
    log:
      'logs/topology_clean_log.txt'
    benchmark:
      'logs/topology_clean_bm.txt'
    conda:
      '../envs/phast.yaml'
    threads: 4
    shell:
      """
      tail -n 1 {input.top} | tail -n 1 | tree_doctor --no-branchlen -n - | sed -r 's/,/ /g; s/(.*);/\\1/' 2>{log} 1>{output.mtree}
      """

rule run_roast:
    #Create merged, multi-aligned MAF file from axtToMaf output
    input:
      mafs=rules.align.input,
      #toasts=expand('{refname}.{species}.toast2.maf',species=SPECIES,refname=config['refName']),
      toasts=expand('results/toast/{refname}.{species}.toast2.maf',species=SPECIES,refname=config['refName']),
      tree='results/tree/topology_clean.txt'
      #tree=config["speciesTree"] if config["speciesTree"] else 'results/tree/topology_clean.txt'
    output:
      'results/roast/roast.maf'
      #'results/roast/roast.maf'
    params:
      #toasts=expand('{refname}.{species}.toast2.maf',species=SPECIES,refname=config['refName']),
      #tree = TREE,
      inputDir='results/toast',
      # these should be in the config file
      roastPs=config["roastParams"],
      #roastMafFiles='*.toast2.maf',
      roastRef=config['refName'],
      outputDir='results/roast',
      #outputFile='{genomedir}/roastOutput/roast.output.X2.maf'
    benchmark:
      'benchmark/roast_bm.txt'
    log:
      '../../logs/roast.log'
    conda:
      '../envs/roast.yaml'
    threads: 1
    shell:
      """
      roast_tree=`cat {input.tree}` #&& \ 
      mkdir -p {params.outputDir} && \
      cd {params.inputDir} && \
      roast {params.roastPs}{params.roastRef} "$roast_tree" {input.toasts} ../roast/roast.maf &>{log}
      """
