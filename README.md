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
mkdir data
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
   
Python 3, conda and snakemake. See [snakemake installation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Patching netToAxt for large, fragmented alignments

For diverged large plant genomes, it is common that large numbers of alignment chains are produced, which can lead to to `netToAxt` error: `chainId 273107809, can only handle up to 268435456`. Currently the only way to address this issue is to patch and recompile the utility and then to hardcode the path into the relevant snakemake rule net_to_axt_to_maf which can be found in `workflow/rules/chain_and_net.smk`. The script below shoud compile a patched netToAxt in `$(pwd)/bin`. 

```
wget http://hgdownload.cse.ucsc.edu/admin/exe/userApps.archive/userApps.v421.src.tgz
tar -xvzf userApps.v421.src.tgz
export MACHTYPE=x86_64
export BINDIR=$(pwd)/bin
export L="${LDFLAGS}"
mkdir -p "$BINDIR"
(cd userApps/kent/src/lib && make)
(cd userApps/kent/src/htslib && make)
(cd userApps/kent/src/jkOwnLib && make)
(cd userApps/kent/src/hg/lib && make)
(cd userApps/kent/src/hg/mouseStuff/netToAxt &&  sed -i 's/maxChainId (256\*1024\*1024)/maxChainId (256\*1024\*1024\*6)/' netToAxt.c && make)
chmod +x $BINDIR/netToAxt
```


