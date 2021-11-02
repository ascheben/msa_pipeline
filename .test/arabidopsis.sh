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
fi

sleep 5s

if [ -z "$TOOL" ]; then
  echo "Please install wget or curl."
else
  $TOOL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_genomic.fna.gz $PARAM $OUTDIR/TAIR10.fa.gz ; gzip -d $OUTDIR/TAIR10.fa.gz; sed -i '/^>/s/.*chromosome />chr/' $OUTDIR/TAIR10.fa; sed -i '/^>/s/ .*//' $OUTDIR/TAIR10.fa
  $TOOL https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/C24/C24.chr.all.v2.0.fasta.gz $PARAM $OUTDIR/C24.fa.gz ; gzip -d $OUTDIR/C24.fa.gz
  $TOOL https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/An-1/An-1.chr.all.v2.0.fasta.gz $PARAM $OUTDIR/An_1.fa.gz ; gzip -d $OUTDIR/An_1.fa.gz
  $TOOL https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Kyo/Kyo.chr.all.v2.0.fasta.gz $PARAM $OUTDIR/Kyo.fa.gz ; gzip -d $OUTDIR/Kyo.fa.gz
  $TOOL https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Ler/Ler.chr.all.v2.0.fasta.gz $PARAM $OUTDIR/Ler.fa.gz ; gzip -d $OUTDIR/Ler.fa.gz
  $TOOL https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Sha/Sha.chr.all.v2.0.fasta.gz $PARAM $OUTDIR/Sha.fa.gz ; gzip -d $OUTDIR/Sha.fa.gz
fi
