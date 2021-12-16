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
  $TOOL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_genomic.fna.gz $PARAM $OUTDIR/TAIR10.fa.gz ; gzip -f -d $OUTDIR/TAIR10.fa.gz
  $TOOL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_genomic.gff.gz $PARAM $OUTDIR/TAIR10.gff.gz; gzip -f -d $OUTDIR/TAIR10.gff.gz
  $TOOL https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/C24/C24.chr.all.v2.0.fasta.gz $PARAM $OUTDIR/C24.fa.gz ; gzip -f -d $OUTDIR/C24.fa.gz
  $TOOL https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/An-1/An-1.chr.all.v2.0.fasta.gz $PARAM $OUTDIR/An_1.fa.gz ; gzip -f -d $OUTDIR/An_1.fa.gz
  $TOOL https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Kyo/Kyo.chr.all.v2.0.fasta.gz $PARAM $OUTDIR/Kyo.fa.gz ; gzip -f -d $OUTDIR/Kyo.fa.gz
  $TOOL https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Ler/Ler.chr.all.v2.0.fasta.gz $PARAM $OUTDIR/Ler.fa.gz ; gzip -f -d $OUTDIR/Ler.fa.gz
  $TOOL https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Sha/Sha.chr.all.v2.0.fasta.gz $PARAM $OUTDIR/Sha.fa.gz ; gzip -f -d $OUTDIR/Sha.fa.gz
fi

if [[ -s "$OUTDIR/TAIR10.fa" ]] && [[ -s "$OUTDIR/TAIR10.gff" ]] && [[ -s "$OUTDIR/C24.fa" ]] && [[ -s "$OUTDIR/An_1.fa" ]] && [[ -s "$OUTDIR/Kyo.fa" ]] && [[ -s "$OUTDIR/Ler.fa" ]] && [[ -s "$OUTDIR/Sha.fa" ]]; then
    echo "Downloaded all test data."
else
    echo "Download failed!"
    exit 1
fi

## Run pipeline in three steps
snakemake --use-conda -p -j 4 $SINGULARITY --configfile .test/arabidopsis_config.yaml -R align
snakemake --use-conda -j 4 $SINGULARITY --configfile .test/arabidopsis_config.yaml -R roast
snakemake --use-conda -j 4 $SINGULARITY --configfile .test/arabidopsis_config.yaml -R call_conservation

## Move results directory
TIMESTAMP=$(date -d "today" +"%Y%m%d%H%M")
mv results .test/arabidopsis_results_$TIMESTAMP
echo "Completed tests analysis, check for results in .test/arabidopsis_results_$TIMESTAMP"
