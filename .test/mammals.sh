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
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simCow.chr6 $PARAM $OUTDIR/simCow.fa ; sed -i '/^>/s/.*chr/>chr/' $OUTDIR/simCow.fa 
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simDog.chr6 $PARAM $OUTDIR/simDog.fa ; sed -i '/^>/s/.*chr/>chr/' $OUTDIR/simDog.fa
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simHuman.chr6 $PARAM $OUTDIR/simHuman.fa ; sed -i '/^>/s/.*chr/>chr/' $OUTDIR/simHuman.fa
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simMouse.chr6 $PARAM $OUTDIR/simMouse.fa ; sed -i '/^>/s/.*chr/>chr/' $OUTDIR/simMouse.fa
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simRat.chr6 $PARAM $OUTDIR/simRat.fa ; sed -i '/^>/s/.*chr/>chr/' $OUTDIR/simRat.fa
fi

if [[ -s "$OUTDIR/simCow.fa" ]] && [[ -s "$OUTDIR/simDog.fa" ]] && [[ -s "$OUTDIR/simHuman.fa" ]] && [[ -s "$OUTDIR/simMouse.fa" ]]  && [[ -s "$OUTDIR/simRat.fa" ]]; then
    echo "Downloaded all test data."
else
    echo "Download failed!"
    exit 1
fi

## Run pipeline in three steps

snakemake --use-conda -p -j 4 $SINGULARITY --configfile .test/mammals_config.yaml -R align
snakemake --use-conda -j 4 $SINGULARITY --configfile .test/mammals_config.yaml -R roast
snakemake --use-conda -j 4 $SINGULARITY --configfile .test/mammals_config.yaml -R call_conservation

## Move results directory
TIMESTAMP=$(date -d "today" +"%Y%m%d%H%M")
mv results .test/mammals_results_$TIMESTAMP
echo "Completed tests analysis, check for results in .test/mammals_results_$TIMESTAMP"
