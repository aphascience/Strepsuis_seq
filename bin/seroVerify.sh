#!/bin/bash

set -eo pipefail

index=$1
read1=$2
read2=$3
sample=$4
serotype=$5

if [[ "$serotype" == "2\*" || "$serotype" == "1*" || "$serotype" == "1*?" || "$serotype" == "2" || "$serotype" == "2*" || "$serotype" == "2*?" ]]; then
    # map to cps2K
    bowtie2 -x $index -1 $read1 -2 $read2 -a --very-sensitive-local --no-unal |
        samtools view -b1 - | samtools sort - -o sorted.bam

    # create vcf for cps2K
    samtools mpileup -Bf $index.fas sorted.bam -o "$sample"_cps2K.pileup 

    verifyevidence=$(grep 483 "$sample"_cps2K.pileup | awk '{print $5}' | grep -o -iE 'A|T|C|G|-' | sort | uniq -ic)
    base483=($verifyevidence)

    echo "Sample,Pos483,Count" > "$sample"_seroVerify.csv
    echo "$sample","${base483[1]}",${base483[0]} >> "$sample"_seroVerify.csv

else
    echo "Verification not required"
fi
