#!/bin/bash

set -eo pipefail

index=$1
read1=$2
read2=$3

# map to cps2K
bowtie2 -x $index -1 $read1 -2 $read2 -a --very-sensitive-local --no-unal |
    samtools view -b1 - | samtools sort - -o sorted.bam

# create vcf for cps2K
samtools mpileup -Bf $index.fas sorted.bam -o cps2K.pileup 
