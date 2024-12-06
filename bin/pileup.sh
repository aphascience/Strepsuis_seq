#!/bin/bash

set -eo pipefail

sam_in=$1
#ref=$2

samtools view $sam_in -q 1 |
samtools sort - -o sorted.bam -O bam
#samtools mpileup............, '-L', '1000', '-f', fasta,
#					 '-Q', str(args.baseq), '-q', str(args.mapq), '-B', out_file_bam_sorted + '.bam'