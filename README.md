# MSA_pipeline

The msa_pipeline bundles existing tools to make multiple alignment of genomes easy.

## Quickstart

Set up the pipeline.

```
git clone https://github.com/ascheben/msa_pipeline.git
```

The pipeline is parametrised using a single config file that can be generated based on the provided template `./config/config.yaml`. The only required file input is a set of fasta files with the suffix `.fa` that should be placed in the directory `msa_pipeline/data`. The pipeline can then be executed in two steps as follows.

```
snakemake --use-conda -j 5 --configfile config/config.yaml -R align
snakemake --use-conda -j 1 --configfile config/config.yaml -R roast
```

The practical reason for splitting the pipeline into two steps is that the first pairwise alignment step is thread-intensive and the second multiple alignment step is memory-intensive and time-intensive. The two steps thus allow adjustment of computational resources, e.g. when submitting the pipeline as single scheduled jobs.


## Testing the pipeline

To test the pipeline before running on your own data, you can align some simulated mammalian sequences. This run should complete in <5 min on a desktop computer and uses 1 thread.

```
cd msa_pipeline
# download test mammalian fasta files (3Mb total) to the data directory
bash .test/mammals.sh ./data
# align sequences to reference and carry out chaining and netting
snakemake --use-conda -j 1 --configfile .test/mammals_config.yaml -R align
# generate multiple alignment
snakemake --use-conda -j 1 --configfile .test/mammals_config.yaml -R roast
```

The main multiple sequence alignment result is written to `./results/roast/roast.maf` and all other intermediate files are written to the directories in `./results`.

A larger set of whole plant genomes can also be used to test the pipeline and the parallelization efficiency. Set the number of threads with the `-j` flag; using 12 threads should allow the pipeline to run in <1h.

```
cd msa_pipeline
# download test Arabidopsis fasta files (700Mb total) to the data directory
bash .test/arabidopsis.sh ./data
# split fasta files and align to reference, then carry out chaining and netting
snakemake --use-conda -j 12 --configfile .test/arabidopsis_config.yaml -R align
# generate multiple alignment
snakemake --use-conda -j 1 --configfile .test/arabidopsis_config.yaml -R roast
```

### Modifying alignment parameters

Alignment parameters can be modified in the config file. For example, to run an alignment using the HOXD70 matrix with custom penalty setting the line could be changed to:

`lastParams: "-j 3 -u 1 -m50 -p HOXD70 -a 700 -b 20 -e 5000 -d 3500"`

## Additional information

### About

This repository is under development and is based on our previous repository [msa_pipeline](https://bitbucket.org/bucklerlab/msa_pipeline/).

### Requirements
   
Python 3 and snakemake. See [snakemake installation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).


