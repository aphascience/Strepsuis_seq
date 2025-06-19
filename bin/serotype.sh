#!/bin/bash

set -eo pipefail

R1=$1
R2=$2
sample=$3
db=$4
def=$5

srst2.py --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 \
    --output serotype_$sample --mlst_db $db --mlst_definitions $def --mlst_delimiter "-" \
    --max_unaligned_overlap 75 --min_coverage 60 --max_divergence 20 --log
