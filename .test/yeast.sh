#!/usr/bin/env bash

which wget 1&>/dev/null
WGET=`echo $?`
which curl 1&>/dev/null
CURL=`echo $?`
TOOL=""
PARAM=""
SNAKEFILE="workflow/Snakefile"
OUTDIR="data"
# Can be set to "--use-singularity"

if [ -z "$1" ]; then
    SINGULARITY=""
elif [ $1 = "--use-singularity" ]; then
    SINGULARITY="--use-singularity"   
else
    echo "Unrecognized positional argument. Only '--use-singularity' is accepted."
fi

if ! command -v snakemake &> /dev/null
then
    echo "snakemake not available. Exiting ..."
    exit 1
fi
if [ -f "$SNAKEFILE" ]; then
    echo "$SNAKEFILE found. Executing test analysis."
else 
    echo "$SNAKEFILE not found. This script must be executed from the directory /msa_pipeline"
    exit 1
fi

if [ "$WGET" == "1" ]; then
  TOOL="wget --no-check-certificate"
  PARAM="-O"
elif [ "$CURL" == "1" ]; then
  TOOL="curl"
  PARAM="-o"
fi

sleep 5s

## Download data

if [ -z "$TOOL" ]; then
  echo "Please install wget or curl."
else
  $TOOL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz $PARAM $OUTDIR/S288C.fa.gz ; gzip -f -d $OUTDIR/S288C.fa.gz
  $TOOL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/662/435/GCA_000662435.2_Sc_YJM993_v1/GCA_000662435.2_Sc_YJM993_v1_genomic.fna.gz $PARAM $OUTDIR/YJM993.fa.gz ; gzip -f -d $OUTDIR/YJM993.fa.gz
  $TOOL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/328/465/GCA_004328465.1_ASM432846v1/GCA_004328465.1_ASM432846v1_genomic.fna.gz $PARAM $OUTDIR/ySR128.fa.gz ; gzip -f -d $OUTDIR/ySR128.fa.gz
  $TOOL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/086/655/GCA_003086655.1_ASM308665v1/GCA_003086655.1_ASM308665v1_genomic.fna.gz $PARAM $OUTDIR/BY4742.fa.gz ; gzip -f -d $OUTDIR/BY4742.fa.gz
fi

if [[ -s "$OUTDIR/S288C.fa" ]] && [[ -s "$OUTDIR/YJM993.fa" ]] && [[ -s "$OUTDIR/ySR128.fa" ]] && [[ -s "$OUTDIR/BY4742.fa" ]]; then
    echo "Downloaded all test data."
else
    echo "Download failed!"
    exit 1
fi

## Run pipeline in three steps

snakemake --use-conda -p -j 4 $SINGULARITY --configfile .test/yeast_config.yaml -R align
snakemake --use-conda -j 4 $SINGULARITY --configfile .test/yeast_config.yaml -R roast
snakemake --use-conda -j 4 $SINGULARITY --configfile .test/yeast_config.yaml -R call_conservation

## Move results directory
TIMESTAMP=$(date -d "today" +"%Y%m%d%H%M")
mv results .test/yeast_results_$TIMESTAMP
echo "Completed tests analysis, check for results in .test/yeast_results_$TIMESTAMP"
