from snakemake import available_cpu_count
from psutil import virtual_memory

# lastal parameters may now be specified in the config file by
# editing the "alignParams" parameter.

available_mem_gb = lambda: '%dG' % (virtual_memory().available >> 30)
containerized: "docker://apscheben/msa_pipeline:latest"

SPECIES = config['species']

# split fasta index numbers
splitN = config["splitFastaN"]
padN = len(str(splitN))
Nlist = list(range(0, splitN))
padList = [str(item).zfill(padN) for item in Nlist]

# Prevent evaluation of lastParams as None
if not config["lastParams"]:
    LAST_PARAMS = ""
else:
    LAST_PARAMS = config["lastParams"]


if config["splitFastaN"] > 1 and config["aligner"] == "last":
    ruleorder: align_split > align_single
else:
    ruleorder: align_single > align_split

# This is the final output, and the rule to call to run all rules in this file
rule align_all:
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
      temp('results/genome/{refname}lastdb_index.done')

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
      '../envs/align.yaml'
    shell:
      """ 
      faSize -detailed {input.fasta} > {params.refSizeFile} 2>{log}
      if [ {params.aligner} = "last" ];then
        lastdb -R 10 -u YASS -c {params.indexBase} {input.fasta} &>>{log}
      elif [ {params.aligner} = "gsalign" ];then
        GSAlign index {input.fasta} {params.indexBase} &>>{log}
      fi
      touch {output} &>>{log}
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
      faSize -detailed {input.fastaFile} > {output}
      """

rule split_fasta:
    input:
      str(rules.lastdb_index.output).format(refname=config['refName']),
      speciesSizeFile='results/genome/{species}.size',
      fastaFile='data/{species}.fa' 
    output:
      splitFa=temp(expand('data/{{species}}.{index}.fa',index=padList)),
      flat=temp('data/{species}.fa.flat'),
      gdx=temp('data/{species}.fa.gdx'),
      splitDummy=temp('data/{species}.split')
    params:
      splitFastaN=config['splitFastaN']
    log:
      'logs/{species}_fasplit_log.txt'
    benchmark:
      'benchmark/{species}-split_bm.txt'
    conda:
      '../envs/split.yaml'
    threads: 1
    shell:
      """
      pyfasta split -n {params.splitFastaN} {input.fastaFile} &>{log} && \
      touch {output.splitDummy}
      """


rule align_single:
    input:
      str(rules.lastdb_index.output).format(refname=config['refName']),
      fastaFile="data/{species}.fa",
      speciesSizeFile='results/genome/{species}.size',
    output:
      name='results/psl/{species}.psl'
    params:
      indexBase='data/{refname}'.format(refname=config['refName']),
      refName=config['refName'],
      species='{species}',
      speciesPath='results/genome/{species}',
      speciesSizeFile='results/genome/{species}.size',
      lastParams=LAST_PARAMS,
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
      '../envs/align.yaml'
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
         minimap2 {params.minimap2Params} {params.refFastaFile} {input.fastaFile} 2>>{log} | samtools sort 2>>{log} | bamToPsl /dev/stdin {output.name} &>>{log}
      elif [ {params.aligner} = "last" ];then
         lastal {params.lastParams} {params.indexBase} {input.fastaFile} {params.lastSplitParams} | maf-convert psl /dev/stdin 2>{log} 1>{output.name}
      elif [ {params.aligner} = "gsalign" ];then
         GSAlign {params.gsalignParams} -r {params.refFastaFile} -q {input.fastaFile} -o {params.speciesPath} -i {params.indexBase} &>>{log} && sed -i "s/^s qry\./s /" {params.speciesPath}.maf && sed -i "s/^s ref\./s /" {params.speciesPath}.maf && maf-convert psl {params.speciesPath}.maf 2>>{log} 1>{output.name}
      fi
      """

rule align_split:
    input:
      splitFa=expand("data/{{species}}.{index}.fa",index=padList),
      splitDummy='data/{species}.split',
      speciesSizeFile='results/genome/{species}.size'
    output:
      psl='results/psl/{species}.psl',
      cmd=temp('results/genome/{species}.cmd'),
      cmdcomp=temp('results/genome/{species}.cmd.completed'),
      splitMaf=temp(expand('results/genome/{{species}}.{index}.maf',index=padList))
    params:
      indexBase='data/{refname}'.format(refname=config['refName']),
      refName=config['refName'],
      splitDir='results/genome/',
      speciesPath='results/genome/{species}',
      lastParams=LAST_PARAMS,
      lastSplitParams=config['lastSplit'],
      splitFa=expand('data/{{species}}.{index}.fa',index=padList)
    log:
      'logs/{species}_lastAlign_split_log.txt'
    benchmark:
      'benchmark/{species}_lastAlign_split_bm.txt'
    conda:
      '../envs/align.yaml'
    threads: config["splitFastaN"]
    shell:
      # This script will align split fasta files to the reference using parafly parallelization
      # It uses some potentially unsafe globbing and rm
      # These should be replaced with expand() inputs by eliminating 0 padding from split fasta names
      """
      ls {params.splitFa}| sed 's@.*/@@'| while read l; do echo "lastal {params.lastParams} {params.indexBase} data/$l > {params.splitDir}${{l%%.fa}}.maf";done > {params.speciesPath}.cmd && \
      ParaFly -c {params.speciesPath}.cmd -CPU {threads} &>>{log} && \
      cat {output.splitMaf} | sed '30,${{/^#/d;}}' | maf-sort /dev/stdin {params.lastSplitParams} | maf-convert psl /dev/stdin |awk '$9!="++"' > {output.psl}
      """
