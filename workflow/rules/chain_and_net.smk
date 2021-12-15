from snakemake import available_cpu_count
from psutil import virtual_memory

# Chaining is performed, we leave the maf file to be used as
# input to multiz-roast as a many-to-many alignment file.

available_mem_gb = lambda: '%dG' % (virtual_memory().available >> 30)
containerized: "docker://lynnjo/msa_pipeline:1.0.1"

SPECIES = config['species']

rule align:
    input:
      # the wildcards in chains will also handle definition of wildcards for fastas
      roastMafs=expand("results/toast/{refname}.{species}.toast2.maf",species=SPECIES,refname=config['refName'])


rule fa_to_2bit:
    # only do for species here.  Ref was done in last_align.smk when 
    # the index was created.  It made these rules simpler.
    input:
      speciesFasta='data/{species}.fa'
    output:
      speciesNib='results/genome/nib/{species}.2bit'
    conda:
      '../envs/ucsc.yaml'
    threads: 1
    shell:
      'faToTwoBit -long {input.speciesFasta} {output.speciesNib}'

rule axt_chain_prenet:
    input:
      'results/psl/{species}.psl'
    output:
      temp('results/psl/{refname}.{species}.preNet')
    log:
      'logs/{refname}.{species}.chainPreNet_log.txt'
    benchmark:
      'benchmark/{refname}.{species}_chainPreNet_bm.txt'
    conda:
      '../envs/ucsc.yaml'
    params:
      refTarget='data/{refname}.fa'.format(refname=config['refName']),
      speciesTarget='data/{species}.fa',
      pslFiles='results/psl/{species}.psl',
      refSizeFile='results/genome/{refname}.size'.format(refname=config['refName']),
      speciesSizeFile='results/genome/{species}.size',
    threads: 1
    shell:
      """
      axtChain -linearGap=loose -psl {params.pslFiles} \
      -faQ -faT {params.refTarget} {params.speciesTarget} /dev/stdout 2>>{log} | \
      chainPreNet /dev/stdin {params.refSizeFile} {params.speciesSizeFile} {output} 2>>{log}
      """

rule chain_net:
    input:
      rules.axt_chain_prenet.output
    output:
      targetNet=temp('results/psl/{refname}.{species}_refTarget.chainNet'),
      speciesNet=temp('results/psl/{refname}.{species}_speciesTarget.chainNet')
    params:
      refSizeFile='results/genome/{refname}.size'.format(refname=config['refName']),
      speciesSizeFile='results/genome/{species}.size',
    threads: 1
    log:
      'logs/{refname}.{species}_chainNet_log.txt'
    benchmark:
      'benchmark/{refname}.{species}-chainNet_bm.txt'
    conda:
      '../envs/ucsc.yaml'
    shell:
      #
      'chainNet {input} {params.refSizeFile} {params.speciesSizeFile} {output.targetNet} {output.speciesNet} 2> {log}' 

rule net_to_axt_to_maf:
    # Converting the *.net output from rule "netting" to an axt file, then running
    # axtToMaf to get a maf file.:
    input:
      twoBitFile=rules.fa_to_2bit.output,
      chainNetOutput=rules.chain_net.output.targetNet,
      chainPreNetFile=rules.axt_chain_prenet.output
    output:
      'results/toast/{refname}.{species}.toast2.maf'
    params:
      ref2BitFile='results/genome/nib/{refname}.2bit'.format(refname=config['refName']),
      refSizeFile='results/genome/{refname}.size'.format(refname=config['refName']),
      speciesSizeFile='results/genome/{species}.size',
      queryPrefix='{species}.',
      targetPrefix='{refname}.'
    threads: 1
    log:
      'logs/{refname}_{species}_netToAxt_log.txt'
    benchmark:
      'benchmark/{refname}_{species}_netToAxt_bm.txt'
    conda:
      '../envs/ucsc.yaml'
    shell:
      'netToAxt {input.chainNetOutput} {input.chainPreNetFile}  {params.ref2BitFile} {input.twoBitFile} /dev/stdout 2>>{log} | axtSort /dev/stdin /dev/stdout 2>>{log}| axtToMaf /dev/stdin {params.refSizeFile} {params.speciesSizeFile} {output} -qPrefix={params.queryPrefix} -tPrefix={params.targetPrefix} 2>>{log}'

