#!/bin/bash

# check args
array=( "$@" )
arrayLength=${#array[@]}

echo "number of args = $arrayLength"
if [ ${arrayLength} -ge 6 ] ; then
	echo "Starting GERP pipeline..."
else
	echo "Error. Please check positional arguments."
	echo "Usage: `basename $0` <maf> <gff3> <chr_name> <ref_name> <outputFile> <threads> <optional: comma separated outgroups>"
	echo "Example: `basename $0` alignment.maf reference.gff chr10 Zea_mays gerp.results.txt 6 Ecurvula,Othomaeum,Zjaponica "
        echo "Ensure all input files are in the working directory and all dependencies (mafSplit,trimal,raxml,gerpcol,phast,parallel) are installed."
	exit 1
fi

# Rename positional args
MAF=$1
# GFF must contain 'CDS' feature for CHR to get 4d sites
GFF=$2
# The chr will be extracted from the maf to calculate a neutral tree
CHR=$3
REF=$4
OUTFILE=$5
THREADS=$6
OUTGROUP=" "
if [ ${arrayLength} -eq 7 ] ; then
    OUTGROUP="-o ${array[6]}"
fi

echo "gff file = $GFF, chr = $CHR"
echo "outgroup = $OUTGROUP"
echo "Checking dependencies..."
# check dependencies
if ! type mafSplit &>/dev/null ; then
	echo "mafSplit not found"
	exit 1
fi
if ! type mafFilter &>/dev/null ; then
	echo "mafFilter not found"
	exit 1
fi
if ! type msa_view &>/dev/null ; then
	echo "phast msa_view not found"
	exit 1
fi
if ! type trimal &>/dev/null ; then
	echo "trimal not found"
	exit 1
fi
if ! type gerpcol &>/dev/null; then
	echo "gerp++ not found"
	exit 1
fi
if ! type parallel &>/dev/null; then
	echo "gnu parallel not found"
	exit 1
fi
echo "Dependencies ok. Running pipeline..."

# make outdir
var=`date +"%FORMAT_STRING"`
now=`date +"%m_%d_%Y_%I_%M_%S"`
mkdir gerp_${now}
# split maf into chr/scaff
mafSplit -byTarget dummy.bed -useFullSequenceName gerp_${now}/ ${MAF}
# The param -useFullSequenceName makes it easy to find mafs for each chr and scaff
# Note that chromosome in GFF needs to match contig name in MAF

echo "cding to gerp_now"
cd gerp_${now}/
echo "current directory = ${PWD}"

# Extract chr from GFF and rename to match MAF
awk -v chrom="${CHR}" '$1==chrom' ../${GFF} | sed "s/^${CHR}\t/${REF}.${CHR}\t/" > ${GFF%%.gff3}.${CHR}.gff3
# Remove alignments with few species or bases
mafFilter -minCol=3 -minRow=3 ${CHR}.maf > ${CHR}.filt.maf
#Get 4d codons from alignment
msa_view ${CHR}.filt.maf --4d --features ${GFF%%.gff3}.${CHR}.gff3 > 4d_codons.ss 
# clean up
rm -f ${GFF%%.gff3}.${CHR}.gff3
rm -f ${CHR}.filt.maf
# convert 4d-codons.ss to fasta
msa_view -i SS --out-format FASTA --randomize 4d_codons.ss | sed '/> /s/^> */>/'| sed '/^>/! s/\*/-/g'| sed '/^>/! s/N/-/g'> 4d_codons.fa
# filter and convert to acceptable phylip format for raxml
trimal -in 4d_codons.fa -out 4d_codons.filt.phy -phylip -gt 0.7
# run raxml to get tree
raxmlHPC-PTHREADS -T ${THREADS} -f a -m GTRGAMMA ${OUTGROUP} -p 12345 -x 12345 -# 100 -s 4d_codons.filt.phy -n mltree
#Extract just the 4d sites (in the 3rd codon positions)
msa_view 4d_codons.ss --in-format SS --out-format SS --tuple-size 1 > 4d_sites.ss
# Strip branch lengths 
tree_doctor --no-branchlen -n RAxML_bestTree.mltree > best_topology.tre
# Calculate neutral tree
phyloFit --tree best_topology.tre --subst-mod REV --msa-format SS 4d_sites.ss
# write out newick tree only
tree_doctor -t phyloFit.mod >neutral.tre

# Run Gerp on all mafs
find . -type f -name '*.maf' | while read alignment ; do
 echo "gerpcol -t neutral.tre -f ${alignment} -e ${REF} -j -z -x .gerp.rates"
done | parallel -j ${THREADS}
mv neutral.tre ../
for RATESFILE in *gerp.rates ; do
    lines=`cat ${RATESFILE}|wc -l`
    chrom=`echo ${RATESFILE} | sed 's/.maf.gerp.rates//' | sed 's/^chr//'`
    # add chrom and 1-based pos columns to rates files
    paste -d "\t" <(yes "${chrom}" | head -n ${lines}) <(seq 1 ${lines}) <(cat ${RATESFILE}) | tr -s ' '| tr ' ' '\t'  > ${RATESFILE%%rates}chrpos.rates
done
# Remove missing positions where GERP was not calculated
cat *.chrpos.rates | awk '$3 != 0' > ../${OUTFILE}
cd ..
rm -rf gerp_${now}
echo "...Finished! Temporary directory gerp_${now} has been deleted."
