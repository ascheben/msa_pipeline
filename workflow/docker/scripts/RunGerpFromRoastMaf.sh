#!/bin/bash

# check args
array=( "$@" )
arrayLength=${#array[@]}
if [ ${arrayLength} -ge 6 ] ; then
        echo "Found correct number of arguments ..."
else
        echo "Error. Please check positional arguments. outgroup is optional and if present, should be the last parameter."
        echo "Usage: `basename $0` <maf> <gff3> <chr_name> <ref_name> <output_file> <threads> <comma separated outgroup> "
        echo "Example: `basename $0` alignment.maf reference.gff chr10 Zea_mays gerp.results.txt 6 Ecurvula,Othomaeum,Zjaponica"
        echo "Ensure all input files are in the working directory and all dependencies (mafSplit,trimal,raxml,gerpcol,phast,parallel) are installed."
        exit 1
fi

# Mount localMachine:/pathToDataDir to docker://tempFileDir/data

DATA_DIR=/tempFileDir/data

maf=$1
gff=$2
chr=$3
ref=$4
outfile=$5
threads=$6
outgroup=" "
if [ ${arrayLength} -eq 7 ] ; then
    outgroup=${array[6]}
fi

echo "RunGerpFromRoastMaf.sh: outgroup = $outgroup"
cd ${DATA_DIR}
/rungerp.sh $maf $gff $chr $ref $outfile $threads $outgroup 

