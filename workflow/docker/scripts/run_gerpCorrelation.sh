#!/bin/bash

# NOTE: I only want the last line echo'd to stdout.  Because the calling
# method in hyperopt will pull this data, and all I want it to pull is the correlation value

# Does this need to go into the docker?  If i is called from the docker, how would that
# return a result?
#  Docker would be needed if I want to run this in R and can't assume user has Rscript,
# but my docker would.

# for testing on cbsu, it is ok - Rscript exists and is on the path.

# check args
array=( "$@" )
arrayLength=${#array[@]}
if [ ${arrayLength} -eq 4 ] ; then
        echo "run_dockerGerp.sh: Found correct number of arguments ..." > /dev/stderr
else
        echo "Error. Please check positional arguments. outgroup is optional and if present, should be the last parameter."
        echo "Usage: run_gerpCorrelation.sh <PACMAD gerp file> <BOP gerp file> <maize-rice coordinates file> <outputFile> "
        exit 1
fi

# Set variables
PACMADfile=$1
BOPfile=$2
maizeRiceFile=$3
outputFile=$4


Rscript /GERP.cor.R $PACMADfile $BOPfile $maizeRiceFile  $outputFile

#echo "Rscript for GERP.cor.R"

result=$(tail -n 1 ${outputFile} | cut -f 3 )

echo $result
