#!/bin/bash

set -eo pipefail

index=$1
read1=$2
read2=$3
sample=$4
serotype=$5

if [[ "$serotype" == "1" || "$serotype" == "1*" || "$serotype" == "1*?" || "$serotype" == "2" || "$serotype" == "2*" || "$serotype" == "2*?" ]]; then
    # map to cps2K
    bowtie2 -x $index -1 $read1 -2 $read2 -a --very-sensitive-local --no-unal |
        ~/biotools/samtools-0.1.18/samtools view -S1 - | ~/biotools/samtools-0.1.18/samtools sort - sorted

    # create vcf for cps2K
    ~/biotools/samtools-0.1.18/samtools mpileup -DBuf $index.fas sorted.bam | ~/biotools/samtools-0.1.18/bcftools/bcftools view -c - > "$sample"_cps2K.vcf 

    verifyevidence=$(grep 483 "$sample"_cps2K.vcf | awk '{print $5}' | grep -o -iE 'A|T|C|G|-' | sort | uniq -ic)
    base483=($verifyevidence)

    echo "Sample,Pos483,Count" > "$sample"_seroVerify.csv
    echo "$sample","${base483[1]}",${base483[0]} >> "$sample"_seroVerify.csv

else
    echo "Verification not required"
fi
