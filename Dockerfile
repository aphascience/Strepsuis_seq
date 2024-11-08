FROM ubuntu:24.04

####################### METADATA ##########################
LABEL Maintainer="Richard Ellis <richard.ellis@apha.gov.uk>"
LABEL base.image=ubuntu:24.04
LABEL software="Strepsuis_seq"
LABEL about.documentation="https://github.com/APHA-CSU/Strepsuis_seq"

#################### INSTALL BASICS #######################

RUN apt-get update && apt-get install --yes --no-install-recommends \
    wget \
    curl \
    unzip \
    bzip2 \
    gawk \
    gcc \
    make \
    python3-dev \
    python3-pip \
    zlib1g-dev \
    libncurses-dev \
    liblzma-dev \
    libbz2-dev \
    libcurl4-openssl-dev

################## INSTALL DEPENDANCIES ###################

RUN mkdir tools && cd tools

## install trimming tool (bbduk?)...
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_39.11.tar.gz --no-check-certificate && \
    tar -xvzf BBMap_39.11.tar.gz && rm -f BBMap_39.11.tar.gz && \
    export PATH=$PWD/bbmap:$PATH

## install bowtie
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-x86_64.zip && \
    unzip bowtie2-2.5.4-linux-x86_64.zip && \
    export PATH=$PWD/bowtie2-2.5.4-sra-linux-x86_64:$PATH

## install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar xjf samtools-1.21.tar.bz2 && rm -f samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && make && make install

## install srst2 (need to update to py3...)