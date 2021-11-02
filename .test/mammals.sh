#!/usr/bin/env bash
which wget 1&>/dev/null
WGET=`echo $?`
which curl 1&>/dev/null
CURL=`echo $?`
TOOL=""
PARAM=""

if [ "$WGET" == "1" ]; then
  TOOL="wget"
  PARAM="-O"
elif [ "$CURL" == "1" ]; then
  TOOL="curl"
  PARAM="-o"
fi

if [ -n "$1" ];then
  if [ -d "$1" ]; then
    OUTDIR="$1"
    echo "Writing to $OUTDIR"
  else
    echo "Provided output directory is invalid. Writing to $PWD" 
    OUTDIR="."
  fi
else
  OUTDIR="."
  echo "Writing to $PWD"
fi

sleep 5s

if [ -z "$TOOL" ]; then
  echo "Please install wget or curl."
else
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simCow.chr6 $PARAM $OUTDIR/simCow.fa ; sed -i '/^>/s/.*chr/>chr/' $OUTDIR/simCow.fa 
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simDog.chr6 $PARAM $OUTDIR/simDog.fa ; sed -i '/^>/s/.*chr/>chr/' $OUTDIR/simDog.fa
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simHuman.chr6 $PARAM $OUTDIR/simHuman.fa ; sed -i '/^>/s/.*chr/>chr/' $OUTDIR/simHuman.fa
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simMouse.chr6 $PARAM $OUTDIR/simMouse.fa ; sed -i '/^>/s/.*chr/>chr/' $OUTDIR/simMouse.fa
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simRat.chr6 $PARAM $OUTDIR/simRat.fa ; sed -i '/^>/s/.*chr/>chr/' $OUTDIR/simRat.fa
fi
