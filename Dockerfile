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
    git \
    make
RUN apt-get install --yes --no-install-recommends \
    python3-dev \
    python3-pip \
    python3-pandas \
    python3-numpy \
    python3-scipy
RUN apt-get install --yes --no-install-recommends \
    zlib1g-dev \
    libncurses-dev \
    liblzma-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    bowtie2

################## INSTALL DEPENDENCIES ###################

WORKDIR /home/tools

## install trimming tool
# RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_39.11.tar.gz --no-check-certificate && \
#    tar -xvzf BBMap_39.11.tar.gz && rm -f BBMap_39.11.tar.gz

# ENV PATH="$PATH:/home/tools/bbmap"

## install samtools-0.1.18
RUN wget https://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2 --no-check-certificate && \
    tar -xjvf samtools-0.1.18.tar.bz2 && rm -f samtools-0.1.18.tar.bz2 && \
    cd samtools-01.18/ && make && cd ..

ENV PATH="$PATH:/home/tools/samtools-0.1.18/"

## install srst2 (updated for python3 compatibility)
RUN git clone https://github.com/APHA-CSU/srst2-py3.git

ENV PATH="$PATH:/home/tools/srst2-py3/scripts/"

WORKDIR /home
