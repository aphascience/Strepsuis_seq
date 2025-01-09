#!/bin/bash

set -eo pipefail

index=$1
read1=$2
read2=$3
sample=$4

# map to cps2K
bowtie2 -x $index -1 $read1 -2 $read2 -a --very-sensitive-local --no-unal |
    samtools view -b1 - | samtools sort - -o sorted.bam

# create vcf for cps2K
samtools mpileup -Bf $index.fas sorted.bam -o "$sample"_cps2K.pileup 

base483=$(grep 483 Test_cps2K.pileup | awk '{print $5}' | grep -o -iE 'A|T|C|G|-' | sort | uniq -ic)

echo "$sample","$base483" > "$sample"_seroVerify.csv
