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
    libcurl4-openssl-dev \
    bowtie2 \
    samtools

################## INSTALL DEPENDANCIES ###################

WORKDIR /home/tools

## install trimming tool (bbduk?)...
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_39.11.tar.gz --no-check-certificate && \
    tar -xvzf BBMap_39.11.tar.gz && rm -f BBMap_39.11.tar.gz && \
    export PATH=$PWD/bbmap:$PATH

## install srst2 (need to update to py3...)

WORKDIR /home
