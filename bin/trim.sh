#!/bin/bash

set -eo pipefail

# Input variables
raw1=$1
raw2=$2
adapters=$3

bbduk.sh in=$raw1 in2=$raw2 out=trimmed_R1.fastq.gz out2=trimmed_R2.fastq.gz ref=$adapters \
        ktrim=r qtrim=rl trimq=20 ftm=5 minlength=100
