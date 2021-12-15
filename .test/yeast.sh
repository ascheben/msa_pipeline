#!/usr/bin/env bash
which wget 1&>/dev/null
WGET=`echo $?`
which curl 1&>/dev/null
CURL=`echo $?`
TOOL=""
PARAM=""

if [ "$WGET" == "1" ]; then
  TOOL="wget --no-check-certificate"
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
  echo "You can also provide the relative or absolute path to an output directory as a positional argument to this script."
fi

sleep 5s

if [ -z "$TOOL" ]; then
  echo "Please install wget or curl."
else
  $TOOL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz $PARAM $OUTDIR/S288C.fa.gz ; gzip -d $OUTDIR/S288C.fa.gz
  $TOOL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/662/435/GCA_000662435.2_Sc_YJM993_v1/GCA_000662435.2_Sc_YJM993_v1_genomic.fna.gz $PARAM $OUTDIR/YJM993.fa.gz ; gzip -d $OUTDIR/YJM993.fa.gz
  $TOOL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/328/465/GCA_004328465.1_ASM432846v1/GCA_004328465.1_ASM432846v1_genomic.fna.gz $PARAM $OUTDIR/ySR128.fa.gz ; gzip -d $OUTDIR/ySR128.fa.gz
  $TOOL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/086/655/GCA_003086655.1_ASM308665v1/GCA_003086655.1_ASM308665v1_genomic.fna.gz $PARAM $OUTDIR/BY4742.fa.gz ; gzip -d $OUTDIR/BY4742.fa.gz
fi
