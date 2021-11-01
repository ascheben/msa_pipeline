from snakemake import available_cpu_count
from psutil import virtual_memory

# lastal parameters may now be specified in the config file by
# editing the "alignParams" parameter.

available_mem_gb = lambda: '%dG' % (virtual_memory().available >> 30)

SPECIES = config['species']

#print(config["splitFastaN"])

if config["splitFastaN"] > 1 and config["aligner"] == "last":
    ruleorder: last_align_split > last_align
    
else:
    ruleorder: last_align > last_align_split

# This is the final output, and the rule to call to run all rules in this file
rule last_all:
    input:
      # the wildcards in psls will also handle definition of wildcards for fastas
      psls=expand("results/psl/{species}.psl",species=SPECIES)


rule lastdb_index:
    input:
      fasta='data/{refname}.fa'
    output:
      # Create a fake file that tells us indexing the ref has finished.
      # file is actually created in the shell: directive with "touch"
      # This rule will also create the ref's 2bit file, which may be
      # used later in net_to_axt (but isn't used at the time of writing)
      'results/genome/{refname}lastdb_index.done',

    params:
      aligner=config['aligner'],
      indexBase='data/{refname}',
      refSizeFile='results/genome/{refname}.size',
    log:
      'logs/{refname}_lastdbIndex_log.txt'
    threads: 2
    benchmark:
      'benchmark/{refname}-lastdb.txt'
    conda:
      '../envs/last.yaml'
    shell:
      """ 
      echo "Running indexing" && \
      faSize -detailed {input.fasta} > {params.refSizeFile}
      if [ {params.aligner} = "last" ];then
        lastdb -R 10 -u YASS -c {params.indexBase} {input.fasta}
      elif [ {params.aligner} = "gsalign" ];then
        GSAlign index {input.fasta} {params.indexBase}
      fi
      touch {output}
      """

rule build_index:
    input:
      str(rules.lastdb_index.output).format(refname=config['refName']),
      fastaFile="data/{species}.fa"
    output:
      "results/genome/{species}.size"
    params:
      indexBase='data/{refname}'.format(refname=config['refName']),
      speciesSizeFile='results/genome/{species}.size',
      refNibDir='results/genome/nib',
      refFastaFile='data/{refname}.fa'.format(refname=config['refName']),
      refNib2Bit='results/genome/nib/{refname}.2bit'.format(refname=config['refName']),
    log:
      'logs/{species}_index_log.txt'
    benchmark:
      'benchmark/{species}-index_bm.txt'
    conda:
      '../envs/ucsc.yaml'
    threads: 1
    shell:
      # This shell will kick off for each fasta in the {genomedir}/fastas folder.  Each instance
      # of this rule gets 1 thread, but multiple lastal commands may be run, depending on the number of species
      # and the number of threads given on the command line.
      # NOTE: more threads means more memory used, and you could run out, so have to
      # temper the number of threads.
      # the file size from faSize is needed is the chain/net steps later as is the
      # ref .2bit file
      #
      """
      mkdir -p {params.refNibDir} && \
      faToTwoBit {params.refFastaFile} {params.refNib2Bit} && \
      faSize -detailed {input.fastaFile} > {params.speciesSizeFile}
      """

rule split_fasta:
    input:
      str(rules.lastdb_index.output).format(refname=config['refName']),
      speciesSizeFile='results/genome/{species}.size',
      fastaFile="data/{species}.fa" 
    output:
      splitDummy='results/genome/{species}_split/{species}.split'
    params:
      splitDir='results/genome/{species}_split',
      speciesPrefix='data/{species}',
      splitFastaN=config['splitFastaN']
    log:
      'logs/{species}_fasplit_log.txt'
    benchmark:
      'benchmark/{species}-split_bm.txt'
    conda:
      '../envs/split.yaml'
    #wildcard_constraints:
    #  N='\d+'
    threads: 1
    shell:
      """
      mkdir -p {params.splitDir} && \
      pyfasta split -n {params.splitFastaN} {input.fastaFile} &>{log}&& \
      mv {params.speciesPrefix}.[0-9]*.fa {params.splitDir} && \
      touch {output.splitDummy}
      """

rule last_align:
    input:
      str(rules.lastdb_index.output).format(refname=config['refName']),
      fastaFile="data/{species}.fa",
      speciesSizeFile='results/genome/{species}.size',
    output:
      name='results/psl/{species}.psl'
      #name='{genomedir}/psl/{species}.psl' if config['aligner'] == 'last' else '{genomedir}/psl/{species}.bam'
    params:
      indexBase='data/{refname}'.format(refname=config['refName']),
      refName=config['refName'],
      species='{species}',
      speciesPath='results/genome/{species}',
      speciesSizeFile='results/genome/{species}.size',
      lastParams=config['lastParams'],
      minimap2Params=config['minimap2Params'],
      gsalignParams=config['gsalignParams'],
      aligner=config['aligner'],
      lastSplitParams=config['lastSplit'],
      refNibDir='results/genome/nib',
      refFastaFile='data/{refname}.fa'.format(refname=config['refName']),
      refNib2Bit='results/genome/nib/{refname}.2bit'.format(refname=config['refName']),
    log:
      'logs/{species}_align_log.txt'
    benchmark:
      'benchmark/{species}-align_bm.txt'
    conda:
      '../envs/minimap2.yaml' if config['aligner'] == 'minimap2' else '../envs/last.yaml'
    threads: 1
    shell:
      # This shell will kick off for each fasta in the {genomedir}/fastas folder.  Each instance
      # of this rule gets 1 thread, but multiple lastal commands may be run, depending on the number of species
      # and the number of threads given on the command line.
      # NOTE: more threads means more memory used, and you could run out, so have to
      # temper the number of threads.
      # the file size from faSize is needed is the chain/net steps later as is the
      # ref .2bit file
      """
      if [ {params.aligner} = "minimap2" ];then
         echo "Running minimap2" && \
         #minimap2 -a -cx asm20 {params.refFastaFile} {input.fastaFile} | samtools sort | bamToPsl /dev/stdin {output}
         minimap2 {params.minimap2Params} {params.refFastaFile} {input.fastaFile} 2>>{log} | samtools sort | bamToPsl /dev/stdin {output.name} &>>{log}
      elif [ {params.aligner} = "last" ];then
         echo "Running lastal" && \
         lastal {params.lastParams} {params.indexBase} {input.fastaFile} {params.lastSplitParams} | maf-convert psl /dev/stdin 2>{log} 1>{output.name}
      elif [ {params.aligner} = "gsalign" ];then
         GSAlign {params.gsalignParams} -r {params.refFastaFile} -q {input.fastaFile} -o {params.speciesPath} -i {params.indexBase} &>>{log} && sed -i "s/^s qry\./s /" {params.speciesPath}.maf && sed -i "s/^s ref\./s /" {params.speciesPath}.maf && maf-convert psl {params.speciesPath}.maf 2>>{log} 1>{output.name}
      fi
      """

rule last_align_split:
    input:
      splitDummy='results/genome/{species}_split/{species}.split'
    output:
      'results/psl/{species}.psl'
    params:
      indexBase='data/{refname}'.format(refname=config['refName']),
      refName=config['refName'],
      splitDir='results/genome/{species}_split',
      speciesPath='results/genome/{species}',
      lastParams=config['lastParams'],
      lastSplitParams=config['lastSplit'],
    log:
      'logs/{species}_lastAlign_split_log.txt'
    benchmark:
      'benchmark/{species}-lastAlign_split_bm.txt'
    conda:
      '../envs/last.yaml'
    threads: config["splitFastaN"]
    shell:
      # This shell will kick off for each fasta in the {genomedir}/fastas folder.  Each instance
      # of this rule gets 1 thread, but multiple lastal commands may be run, depending on the number of species
      # and the number of threads given on the command line.
      # NOTE: more threads means more memory used, and you could run out, so have to
      # temper the number of threads.
      # the file size from faSize is needed is the chain/net steps later as is the
      # ref .2bit file
      #
      """
      echo "Running lastal on split files" && \
      ls {params.splitDir}/*fa| while read l; do echo "lastal {params.lastParams} {params.indexBase} $l {params.lastSplitParams} | maf-convert psl /dev/stdin| gzip > $l.psl.gz";done > {params.splitDir}/lastal.cmd && \
      ParaFly -c {params.splitDir}/lastal.cmd -CPU {threads} &>>{log} && \
      zcat {params.splitDir}/*.psl.gz | awk '$9!="++"' > {output} 
      """
