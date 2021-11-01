#!/bin/bash

usage() { echo "Usage: $0 -d <FASTA_DIR> -r <REFERENCE_NAME> [--make-config-only] [-t <THREADS>] [-n <TREE>]" 1>&2;
          echo "Example 1: $0 -d /path/to/mammal/fastas -r GRCh38"
          echo "Example 2: $0 -d /path/to/plant/fastas -r TAIR10 -t 6 -n /path/to/tree.nwk"
          echo "Example 2: $0 -d /path/to/plant/fastas -r TAIR10 --make-config-only"

          exit 0; 
 }

WARN="\033[33mWARNING\033[0m"

for arg in "$@"; do
  shift
  case "$arg" in
    "--help") set -- "$@" "-h" ;;
    "--make-config-only") set -- "$@" "-p" ;;
    "--"*) usage ; exit 2 ;;
    *) set -- "$@" "$arg";;
  esac
done

while getopts "hpd:r:t:n:" o; do
    case "${o}" in
        d)
            d=${OPTARG}
            ;;
        r)
            r=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            # positive integer
            [[ $t =~ ^[0-9]+$ ]] || usage
            ;;
        n)
            n=${OPTARG}
            ;;
        p)
            p="yes"
            ;;
        h)
            usage
            ;;
        *)
            echo "Invalid options! Please check input"
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${d}" ] || [ -z "${r}" ]; then
    echo "Missing input for -d or -r argument(s)!"
    usage
fi

# Default threads is 1
if [ -z "${t}" ]; then
    t=1
fi


fastadir=${d}
refspecies=${r}
threads=${t}
tree=${n}

## ASSUMPTION IS THAT THE ANALYSIS IS RUN FROM msa_pipeline DIRECTORY

## GET CURRENT WORKING DIRECTORY
wd=$(pwd)

# check for snakefile
if [ ! $wd/Snakefile ]
    then
        echo "msa_pipeline snakefile not found in working dir: ${wd}! Exiting ..."
        exit 1
fi

if ! command -v singularity &> /dev/null
then
    echo "singularity not available. Exiting ..."
    exit 1
fi

## GET DIRECTORIES FOR FASTA FILES

# Warn if output dirs already exist
if [ ! -d $fastadir/data ]
    then
        mkdir $fastadir/data
    else
        echo -e "$WARN Directory $fastadir/data already exists - old results may be used by snakemake!"
fi
if [ ! -d $fastadir/analysis ]
    then
        mkdir $fastadir/analysis
    else
        echo -e "$WARN Directory $fastadir/analysis already exists - old results may be used by snakemake!"
fi
## Check input FASTA files

counter=0
for infile in ${fastadir}/*.fa
do
    if [ -f  "${infile}" ];then
        counter=$((counter+1))
        # count lines with lowercase characters
        masked_lines=`grep -v '>' ${infile} |tr -d [:upper:] | sed '/^$/d' | wc -l`
        if [ "${masked_lines}" -lt 1 ]
            then
                echo -e "$WARN ${infile} does not appear to be softmasked!"
                sleep 5s
        fi
    fi
done

if [ $counter -lt 1 ]; then
    echo "No FASTA files ending in .fa found. Exiting ..."
    exit 1
fi

## SET UP SINGULARITY FROM DOCKER IMAGE
## LATEST TAG CAN BE FOUND HERE - https://hub.docker.com/r/maizegenetics/msa_pipeline/tags?page=1&ordering=last_updated
if [ ! -f msa_pipeline.simg ]
    then
        echo "Singularity container msa_pipeline.simg not found in ${wd} . Building container ..."
        sleep 2s
        singularity build msa_pipeline.simg docker://maizegenetics/msa_pipeline:latest
    else
        echo "Using existing singularity container: msa_pipeline.simg"
fi

if [ -f $wd/config/config_$refspecies.yaml ];then
    echo -e "$WARN Existing snakemake config file will be used: $wd/config/config_$refspecies.yaml"
    sleep 2s
else

    if [ -n "$tree" ]; then
        ## IF TREE IS PROVIDED THEN USING THAT TREE, REMOVE BRANCH LENGTHS, REPLACE COMMAS WITH SPACE
        tail -n 1 $tree > $fastadir/data/topology_lastLine.nwk
        singularity run -B $wd:/msa_pipeline msa_pipeline.simg tree_doctor --no-branchlen -n $fastadir/data/topology_lastLine.nwk > $fastadir/data/topology_roast.nwk
        ## FORMATTING FOR THE CONFIG FILE
        sed -i -r 's/,/ /g; s/(.*);/\1/; 1s/^/speciesTree: "/; s/$/"/' $fastadir/data/topology_roast.nwk
    else
        ## GET THE TREE FROM ALL THE FASTA FILES
        singularity run -B $wd:/msa_pipeline msa_pipeline.simg mashtree.pl --tempdir tmp --sketch-size 20000 --genomesize 1000000000 --numcpus $threads --mindepth 0 $fastadir/*.fa > $fastadir/data/topology.nwk
        
        ## REMOVE BRANCH LENGTHS AND REPLACE COMMAS WITH SPACE
        tail -n 1 $fastadir/data/topology.nwk > $fastadir/data/topology_lastLine.nwk
        singularity run -B $wd:/msa_pipeline msa_pipeline.simg tree_doctor --no-branchlen -n $fastadir/data/topology_lastLine.nwk > $fastadir/data/topology_roast.nwk
        ## FORMATTING FOR THE CONFIG FILE
        sed -i -r 's/,/ /g; s/(.*);/\1/; 1s/^/speciesTree: "/; s/$/"/' $fastadir/data/topology_roast.nwk
    fi
    
    ## READ THE FILES IN THE FASTA DIR SPECIFIED BY USER AND STORE ALL THE FASTA FILES IN AN ARRAY
    array=($(ls $fastadir/*.fa | sed -e 's/\.fa$//' | sed 's@.*/@@'))
    
    ## TESTING THE ARRAY
    #printf '%s\n' "${array[@]}"
    
    ## REMOVE THE REF SPECIES FROM THE ARRAY
    new_array=()
    for value in "${array[@]}"; do
        [[ $value != $refspecies ]] && new_array+=($value)
    done
    
    ##TESTING SPECIES ARRAY
    #printf '%s\n' "${new_array[@]}"
    
    ## FORMATTING AND PRINTING SPECIES FOR THE CONFIG FILE
    for i in "${new_array[@]}";do
       echo "  - $i"
       # or do whatever with individual element of the array
    done > $fastadir/data/species.txt
    
    ## MAKING A COPY OF TEMPLATE CONFIG FILE
    cp $wd/config/initConfig.yaml $wd/config/config_$refspecies.yaml
    
    ## ADDING REFERENCE SPECIES, OTHER SPECIES AND TREE TO THE CONFIG FILE 
    sed -i "s/PLACEHOLDERREFERENCE/$refspecies/" $wd/config/config_$refspecies.yaml
    sed -i "s|PLACEHOLDERGENOMEDIR|$fastadir|" $wd/config/config_$refspecies.yaml
    sed -e "/PLACEHOLDERSPECIES/{r $fastadir/data/species.txt" -e 'd}' -i $wd/config/config_$refspecies.yaml
    cat $fastadir/data/topology_roast.nwk >> $wd/config/config_$refspecies.yaml

fi

if [ -n "$p" ]; then
    echo "Preparing pipeline complete."
    echo "Snakemake config file location: $wd/config/config_$refspecies.yaml"
else
    ## RUNNING THE FINAL STEPS FOR ALIGNMENT
    singularity run -B $wd:/msa_pipeline msa_pipeline.simg snakemake -j $threads -p -s $wd/Snakefile --configfile $wd/config/config_$refspecies.yaml --directory /msa_pipeline chain_finalize
    singularity run -B $wd:/msa_pipeline msa_pipeline.simg snakemake -j $threads -p -s $wd/Snakefile --configfile $wd/config/config_$refspecies.yaml --directory /msa_pipeline roast_gerp_finalize
    echo "Pipeline run complete."
    echo "Final multiple alignment written to: $fastadir/roast.output.X2.maf"
fi
