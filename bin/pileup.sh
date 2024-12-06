#!/bin/bash

set -eo pipefail

sam_in=$1
ref=$2

samtools view $sam_in -huq 1 -o map.bam #| TODO: Fix piping of these commands
samtools sort map.bam -o sorted.bam
samtools mpileup -L 1000 -f $ref -Q 20 -q 1 sorted.bam -o output.pileup

rm map.bam
rm sorted.bam
