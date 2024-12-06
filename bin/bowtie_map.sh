#!/bin/bash

set -eo pipefail

read1=$1
read2=$2
index=$3

bowtie2 -x $index -1 $read1 -2 $read2 -a -S map.sam \
    --very-sensitive-local \
    --no-unal
