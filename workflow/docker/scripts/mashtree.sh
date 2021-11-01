# The below commands can generate a newick tree to guide alignment

# install mashtree and phast 
#conda install -c bioconda mashtree
#conda install -c bioconda phast
mashtree --tempdir /path/to/tmp --sketch-size 20000 --genomesize 1000000000 --numcpus 12 --mindepth 0 *fa.gz > taxa.nwk
# remove branch lengths
tree_doctor --no-branchlen -n taxa.nwk > topology.nwk
