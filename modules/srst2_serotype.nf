process SRST2_SEROTYPE {
    tag "$pairId"
    errorStrategy 'ignore'
    maxForks 2
    publishDir "${params.outdir}/serotype", mode: "copy"

    input:
        tuple val(pairId), file(R1), file(R2)

    output:
	    tuple val(pairId), file("serotype_${pairId}_*.txt"), emit: serotype_results
        tuple val(pairId), env(serotype), emit: srst2_sero
    
    script:
    """
    srst2.py --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 \
             --output serotype_$pairId --mlst_db ${params.sero_db} --mlst_definitions ${params.sero_def}\
             --mlst_delimiter "-" --max_unaligned_overlap 75 --min_coverage 60 --max_divergence 20 --log

    serotype=\$(awk 'FNR == 2 {print \$2}' *_serotype__results.txt)
    """
}