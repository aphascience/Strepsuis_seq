#!/bin/bash

set -eo pipefail

# Input variables
ID=$1
raw1=$2
raw2=$3
adapters=$4

bbduk.sh in=$raw1 in2=$raw2 out="$ID"_trimmed_R1.fastq.gz out2="$ID"_trimmed_R2.fastq.gz ref=$adapters \
        ktrim=r qtrim=rl trimq=20 ftm=5 minlength=100
