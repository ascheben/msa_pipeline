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

if [ -z "$TOOL" ]; then
  echo "Please install wget or curl."
else
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simCow.chr6 $PARAM simCow.fa; sed -i '/^>/s/.*chr/>chr/' simCow.fa 
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simDog.chr6 $PARAM simDog.fa; sed -i '/^>/s/.*chr/>chr/' simDog.fa
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simHuman.chr6 $PARAM simHuman.fa; sed -i '/^>/s/.*chr/>chr/' simHuman.fa
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simMouse.chr6 $PARAM simMouse.fa; sed -i '/^>/s/.*chr/>chr/' simMouse.fa
  $TOOL https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simRat.chr6 $PARAM simRat.fa; sed -i '/^>/s/.*chr/>chr/' simRat.fa
fi
