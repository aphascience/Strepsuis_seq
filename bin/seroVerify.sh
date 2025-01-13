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

    # extract basecall and  evidence from vcf
    ref483=$(grep -w 483 "$sample"_cps2K.vcf | awk '{print $4}')
    alt483=$(grep -w 483 "$sample"_cps2K.vcf | awk '{print $5}')
    evidence=$(grep -w 483 "$sample"_cps2K.vcf | awk '{print $8}' | awk -F ';' '{print $5}' | awk -F '=' '{print $2}')
    if [ $alt483 = '.' ]; then
        base483=$ref483
    else
        base483=$alt483
    fi

    echo "Sample,Pos483,RefF,RefR,AltF,AltR" > "$sample"_seroVerify.csv
    echo "$sample","$base483","$evidence" >> "$sample"_seroVerify.csv

else
    echo "Verification not required"
fi
