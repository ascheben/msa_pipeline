# MSA_pipeline

The msa_pipeline bundles existing tools to make multiple alignment of genomes easy.

## Quickstart

Set up the pipeline.

```
git clone https://<username>@bitbucket.org/bucklerlab/msa_pipeline.git
```

Then run the pipeline using our bash wrapper, providing a directory with FASTA files and a reference species name. Note that `singularity` must be installed to use the wrapper and that the pipeline must be executed from the msa_pipeline directory.
Note also that your fasta file names must NOT contain any periods other than the period preceding the "fa" suffix on the files.  Additional periods in the name create problems for Snakemake parsing.



```
cd msa_pipeline/
./msa_pipeline.sh
Usage: msa_pipeline.sh -d <FASTA_DIR> -r <REFERENCE_NAME> [--make-config-only] [-t <THREADS>] [-n <TREE>]
```

To test the pipeline before running on your own data, you can align some public yeast genomes. This run should complete in <5 min on a desktop computer and uses 1 thread by default.

```
mkdir test
wget -O test/BY4742.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/086/655/GCA_003086655.1_ASM308665v1/GCA_003086655.1_ASM308665v1_genomic.fna.gz
wget -O test/ySR128.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/328/465/GCA_004328465.1_ASM432846v1/GCA_004328465.1_ASM432846v1_genomic.fna.gz
wget -O test/S288C.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
wget -O test/YJM993.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/662/435/GCA_000662435.2_Sc_YJM993_v1/GCA_000662435.2_Sc_YJM993_v1_genomic.fna.gz
gzip -d test/*gz
./msa_pipeline.sh -d ${PWD}/test -r BY4742
```

The main multiple sequence alignment result is written to `./test/roast.output.X2.maf` and all other intermediate files are written to the new directories in `./test/`.

### Modifying alignment parameters

The `msa_pipeline.sh` wrapper can also be run with the flag `--make-config-only` to only produce a config file for snakemake containing the target fasta directory and species tree.  This file is written to `${wd}/config/config_${refspecies}.yaml`. LAST alignment parameters can be passed by modifying the "alignParams" entry.

For example, to run an alignment using the HOXD70 matrix with custom penalty setting the line could be changed to:

`alignParams: "-j 3 -u 1 -m50 -p HOXD70 -a 700 -b 20 -e 5000 -d 3500"`

Then the alignment steps of the pipeline can be executed, it will automatically use the config file based on the reference species name passed with the `-r` flag. Here we align the same yeast genomes but this time with our custom alignment parameters. 

```
mkdir test_custom
cp ${PWD}/test/*fa test_custom/
./msa_pipeline.sh -d ${PWD}/test_custom -r BY4742
```

**Note**: We recommend using Singularity to execute the snakemake pipeline. If using Docker, the config file parameter for `genomeDir:` must be a file path relative to the `msa_pipeline` directory. For example:

```

genomeDir: /workdir/msa_pipeline/test

```

Must be modified to:

```

genomeDir: test

```

### Calculating conservation scores

Using the multiple alignment, we can calculate a neutral tree based on fourfold degenerate codon positions and then use this to calibrate a conservation scoring method. Using the example yeast genomes aligned above, we can generate conservation scores across the genome as follows. 

```
wget -O test/S288C.gff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
gzip -d test/S288C.gff.gz
singularity run msa_pipeline.simg docker/scripts/rungerp.sh test/roast.output.X2.maf test/S288C.gff NC_001143.9 S288C test/yeast_S288C_gerp_rates.txt 1
```

## Additional information

### About

This repository was derived from Evan Rees' snakemake template repository. The repository holds code to align species to a reference and create a multi-alignment file from multiz-roast program. The repository also holds scripts to analyze the multi-sequence aligment file created from roast.  All programs/applications and scripts needed to run both the msa_pipeline and analysis on it are compiled into a single docker image.   

### Requirements
   
   Singularity, Python 3 and miniconda. See [miniconda installation guide][miniconda install].

### Detailed setup 
The msa_pipeline project is developed to be run within a docker or singularity container.  
Dependencies are managed via Docker and Conda.  The rules and example config files may be cloned from
 this repository:

```
git clone https://<username>@bitbucket.org/bucklerlab/msa_pipeline.git
```
Once basic snakemake structure with rules has been downloaded, pull the docker file from docker hub, where <latest tag> is the 
latest tag from the docker hub: 

``
docker pull maizegenetics/msa_pipeline:<latest tag>
``

The msa_pipeline directory is the top
level directory you will mount/bind for docker or singularity.


##### Repository layout:

| | | | |
| --- | --- | --- | --- |
| ├── | `.gitignore` | | Files to exclude from version control. IntelliJ can still track local changes to files excluded by `.gitignore`. Provided file ignores **everything** not mentioned in this table. |
| ├── | `LICENSE.md` | | Distribution license for your project. See [choosealicense.com][license].|
| ├── | `README.md` | | this! |
| ├── | `Snakefile` | | Main Snakemake workflow |
| ├── | `msa_pipeline.sh` | | Convenient wrapper to execute workflow |
| ├── | `config.yaml` | | Config file storing workflow parameters accessed by rules. Copy and conform to represent your species names |
| ├── | `envs` | | Rule-specific conda `environment.yaml` files |
| | └── | `snakemake.yaml` | Snakemake environment specification |
| ├── |`rules` | | Individual rule files in `.smk` format |
| ├── |`scripts` | | Scripts called by rules |


### Setup data directories for the msa_pipeline pipeline

    # copy the example config_maizeGrass.yaml file to a file name of your choice, and edit the "refName'
    # and 'species' fields as appropriate for your run.
        nano my_config.yaml
        ...

    The configuration file has parameter "genomeDir"  This is by default name "data".
     It is a folder which must be created under the msa_pipeline folder.  It is not created
     as part of the repository as the data to be stored here will be user specific and not
     stored in git hub. You may name the "data" directory anything you would like, but it
     must be created and the name must be consistent with the name you have in your configuration 
     file for the "genomeDir" parameter.
     
     The configuration file also has a parameter named "alignParams".  This is used to specify additional 
     parameters to the lastal alignment.  The paramenter in the example config file is a default for output
     (3 = gapped). Edit this field if you wish to run lastal with different or additional parameters.
     
     The configuration file  has a parameter named lastSplit.  There are 2 examples of parameters for
     lastSplit in the example conf_maizeGrass.yaml file.  The default is not to run lastSplit.  This may be
     changed by the user by changing the "lastSplit" parameter that is commented vs un-commented.
     
     Finally, the configuration file has a parameter named speciesTree.  This is a species tree required
     to run the "roast" program.  The species tree must be enclosed by double quotes, individually separated by a space or
     parenthesis, and should NOT include an ending semi-colon.  You may already have a tree for your species.  If not, there
      are programs available on the web for creating a species tree.  An example is the NCBI Taxonomy Browser. 
     Another program that generates a species tree is "mashtree".  An example tree in approriate format for the roast program is below:
     
     speciesTree: "((((Oropetium_thomaeum Eragrostis_curvula) (Cenchrus_americanus Setaria_italica)) (Saccharum_spontaneum Sorghum_bicolor)) Zea_mays)"
     
    For convenience, the program examples below assume you have chosen to keep the genomeDir name as shown in the 
    example config file. 
    
    Under the msa_pipeline directory, create a subdirectory named data and underneath that, a directory named "fastas".  
    Into msa_pipeline/data/fastas, copy your reference fasta file.  The file name must match the name of the reference 
    as specified in the configuration file's "refName" parameter, and it must have suffix ".fa".
    Example:  If the refName in the config file is "ZmB73v4", your fasta must be named "ZmB73v4.fa".
    This naming convention is to facilitate Snakemake file recognition.
    
    Into the msa_pipeline/data/fastas directory, copy the fasta file for every genome you have listed in
    the configuration file's species list.  These species show up under the "species:" tag in the config file.
    Replace the example "Athaliana" and "Sbicolor" with your species' names.  Use underscores
    if desired, but do not use spaces or periods, as those have meaning in the pipeline.  Add as
    many species as you wish to run.
    Copy your reference and species fastas to the data/fastas directory,
    using the fasta naming convention of <speciesname>.fa
    
    For example:  The refName in the example config file is "ZmB73v4", the species are "Athaliana" and "Sbicolor".
    Into directory msa_pipeline/data/fastas put the reference fasta file "ZmB73v4.fa".  Into this same directory
    you would add files Athaliana.fa and Sbicolor.fa.
  

#### Naming restraints
    The names of your reference and species must not contain periods (.).  The multiz-roast step uses the "." character 
    as a delimiter separating the species name from the chromosome/contig. If the species or reference name contains ".", 
    it will confuse roast.
    
### Run the msa_pipeline code

    You may (but don't have to ) first create the files that will be input for the roast program.  These can be run separately,
    by running the msa_pipeline code up through the chain_finalize rule.  Once all your species have successfully completed
    the alignment series of lastal/chain/netting, then you can run the roast_finalize rule to create the multi-sequence alignment.
    NOTE: You may instead run only the roast_gerp_finalize step.  This will run all the rules from aligning/chaining/netting
    that are needed to complete the roast step and then run the multiz-roast command.
    
    If you choose run the align/chain/net steps before the multiz step follow the instruction below:
    
    To run with your new config file use the --configfile parameter when invoking Snakemake
    example of running via Docker, with a  configfile called  msa_pipeline/config/myConfig.yaml:
    
    ```
     docker run --name <container name> --rm -v /<user prefix>/msa_pipeline:/analysis -t maizegenetics/msa_pipeline:<tag> snakemake -j <number of threads> -p -s /analysis/Snakefile --configfile /analysis/config/myConfig.yaml --directory /analysis chain_finalize
     ```
 
     Note in the example above, the -v command is used to mount the user's msa_pipeline directory to
     the /analysis directory in the docker container.  Docker will create the /analysis directory as needed
     inside the container.  See docker documentation for an explanation of the -v mount command:

    Once this step has finished successfully, you can run the roast step, using the same directory structure as
    above.  This allows Snakemake roast rule to find the files it needs as input.  You may also skip the docker 
    run above and run just this step.
   
    ```
    docker run --name <container name> --rm -v /<user prefix>/msa_pipeline:/analysis -t maizegenetics/msa_pipeline:<tag> snakemake -j <number of threads> -p -s /analysis/Snakefile --configfile /analysis/config/myConfig.yaml --directory /analysis roast_gerp_finalize
    ```
   
    If you are confident of your code, you may run just the roast_gerp_finalize step.  Snakemake will identify and run all rules 
    necessary to create files for the roast_gerp_finalize step.

    Singularity runs are similar to docker.  Many organizations require the user of Singularity over Docker due to permission and security issues. 
    It is simple to create a Singularity image from the msa_pipeline docker.  First, create the singularity container from the latest docker image.  The image here is 
    called "msa_pipeline.simg", and it is built from the docker image stored at docker hub in docker://maizegenetics/msa_pipeline:latest_tag where
    latest_tag is the latest tag number from the repository, e.g. docker://maizegenetics/msa_pipeline:0.05
    
    ``singularity build msa_pipeline.simg docker://maizegenetics/msa_pipeline:<latest_tag>
    ``

    You now have a singularity image.  To this image, you need to mount the folder that contains the msa_pipeline repository code
    to the /analysis directory inside the singularity container.  You then tell singularity to run "snakemake" as the command, you tell it how many
    threads it has available (-j 30 below), the location of the Snakemake file (-s /analysis/Snakefile), the name of the configuration file
    you are using (--configfile /analysis/config/myConfig.yaml), the directory internal to Singularity from which it should run
    (--directory /analysis)  and finally, the rule you wish to complete (chain_finalize).  Snakemake will run all rules necessary
    to create the output from the rule you specified.  Note that in these examples, the /analysis directory is a container
    specific directory where the msa_pipeline folder has been mounted.
    
    The -p option below instructs Singularity to print out the commands it is running.  This is useful to keep track of the program's progress.
   
    `` 
    singularity run -B /workdir/lcj34/gerpGrassAlign/msa_pipeline:/analysis msa_pipeline.simg snakemake -j 30 -p -s /analysis/Snakefile --configfile /analysis/config/myConfig.yaml --directory /analysis chain_finalize
    ``

    To perform a dry run of the msa_pipeline code change the -p option to -np (prints list of rules to run and workflow summary, but does not actually execute anything).
    This is very useful to verify you have your fasta files names correctly and stored in the correct place.
    
    You may also run snakemake using the --debug-dag option.  This will print the valu of wildcards and will
    assist with debugging if there are errors.
    
    Once you have successfully run the alignment steps, you can run multiz-roast as below:
   `` 
    singularity run -B /workdir/lcj34/gerpGrassAlign/msa_pipeline:/analysis msa_pipeline.simg snakemake -j 30 -p -s /analysis/Snakefile --configfile /analysis/config/myConfig.yaml --directory /analysis roast_gerp_finalize
    ``

    ### Other Snakemake Commands
    To list available target rules replace the "chain_finalize" directive in the singularity command above with --lt:

    To generate DAG of output files for visualization, change the "chain_finalize" directive in the singularity command above to:
        --dag | dot -Tpng > my_dag.png

    To generate a rulegraph of output Files for visualtionzation, change the "chain_finalize" directive in the singularity command above to:
        --rulegraph | dot -Tpng > my_rulegraph.png

### Analysis scripts
The docker associated with this pipeline also includes scripts that may be run to further analyze the data.  These scripts
are defined below.

#### GERP
The GERP script requires the MAF file output by the roast_gerp_finalize rule of this pipeline.  This script takes as input the
following:

* multi-sequence alignment file (MAF File output by multiz-roast)
* gff3 file for the reference species
* chromosome name as it appears in the MAF file 
* reference name as it appears in the MAF file 
* number of threads to use for processing


#### IntelliJ IDEA configuration

Install and enable required/optional plugins.

| plugin | description | required? |
| --- | --- | --- |
| [Python Community Edition] | Python editing /interpreter support | yes |
| [SnakeCharm] | SnakeMake language support | yes |
| [Bitbucket Linky] | Bitbucket integration | yes |
| [.ignore] | .gitignore formatting support | no |
| [Markdown Navigator] | Preview markdown files while editing | no |

Set the python binary from the newly installed `snakemake` environment as the SDK for a new/exiting project:

[Configure Snakemake support][snakemake setup]. 

#### Repository Specifics

[.ignore]: https://plugins.jetbrains.com/plugin/7495--ignore
[Python Community Edition]: https://plugins.jetbrains.com/plugin/7322-python-community-edition
[SnakeCharm]: https://plugins.jetbrains.com/plugin/11947-snakecharm/versions
[Markdown Navigator]: https://plugins.jetbrains.com/plugin/7896-markdown-navigator
[Bitbucket Linky]: https://plugins.jetbrains.com/plugin/8015-bitbucket-linky
[miniconda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
[cbsu-conda]: https://biohpc.cornell.edu/doc/software_install_slides_1.pdf
[medium-conda-overview]: https://medium.com/data-science-in-practice/saving-the-environment-with-anaconda-ad68e603d8c5
[Snakemake]: https://snakemake.readthedocs.io/en/stable/index.html
[snakemake setup]: https://github.com/JetBrains-Research/snakecharm/wiki#setup-snakemake-support
[Snakemake Tutorial]: https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html
[miniconda install]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
[license]: https://choosealicense.com/
[guidelines]: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility
