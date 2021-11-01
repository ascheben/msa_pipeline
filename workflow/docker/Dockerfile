# The docker is for use with the Buckler Lab msa_pipeline found here:
# https://bitbucket.org/bucklerlab/msa_pipeline/src/master/
# This file assumes that file openjdk-13.0.1_linux-x64_bin.tar.gz exists in
# this directory.

FROM continuumio/miniconda3

WORKDIR /

COPY ./scripts/rungerp.sh /rungerp.sh
COPY ./scripts/mashtree.sh /mashtree.sh
COPY ./scripts/RunGerpFromRoastMaf.sh /RunGerpFromRoastMaf.sh
COPY ./scripts/CreateRoastSpeciesTree.sh /CreateRoastSpeciesTree.sh
COPY ./scripts/GERP.cor.R /GERP.cor.R
COPY ./scripts/run_gerpCorrelation.sh /run_gerpCorrelation.sh

# make scripts executable
RUN chmod +x rungerp.sh
RUN chmod +x RunGerpFromRoastMaf.sh
RUN chmod +x mashtree.sh
RUN chmod +x CreateRoastSpeciesTree.sh
RUN chmod +x run_gerpCorrelation.sh


RUN apt-get update  # Ensure the package list is up to date
# libraries needed to compile last-1021 aligner
RUN apt-get -y install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev
RUN apt-get install -y vim
RUN apt-get install -y git
RUN apt-get install -y zip gzip tar
RUN apt-get install -y g++
RUN apt-get install -y rsync

# assumes the openjdk file exists in the same directory as this docker file.
RUN mkdir -p /programs

ADD ./openjdk-13.0.1_linux-x64_bin.tar.gz /programs

WORKDIR /programs
RUN ls

# add /programs (which contains axtChain, etc)
ENV PATH "$PATH:/programs"
ENV PATH "$PATH:/opt/conda/bin"

WORKDIR /programs

# get USCS programs - this is for chaining and netting, install into /programs
WORKDIR /programs
RUN rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./

RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda install -y python=3.6 psutil snakemake multiz last \
    && conda install -y trimal gerp raxml \
    && conda install -c bioconda/label/cf201901 mashtree \
    && conda install -y -c bioconda phast

# install R for running R script
RUN conda config --add channels r && conda install -y -c r r-base

# Add Java to the PATH
ENV PATH "$PATH:/programs/jdk-13.0.1/bin"

# URL for Kotlin install
ARG KOTLIN_URL="https://github.com/JetBrains/kotlin/releases/download/v1.3.61/kotlin-compiler-1.3.61.zip"

# installs kotlin
RUN cd /usr/lib && \
    wget -q $KOTLIN_URL && \
    unzip kotlin-compiler-*.zip && \
    rm kotlin-compiler-*.zip && \
    rm -f kotlinc/bin/*.bat

# add kotlin to the PATH
ENV PATH $PATH:/usr/lib/kotlinc/bin

WORKDIR /
