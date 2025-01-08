#!/bin/bash

set -eo pipefail

# Input variables
ID=$1
raw1=$2
raw2=$3
adapters=$4

bbduk.sh in=$raw1 in2=$raw2 out="$ID"_trimmed_R1.fastq.gz out2="$ID"_trimmed_R2.fastq.gz ref=$adapters \
        ktrim=rl k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10 ftm=5
