process SRST2_MLST {
    tag "$pairId"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/MLST", mode: "copy"

    input:
        tuple val(pairId), file(R1), file(R2)

    output:
	    tuple val(pairId), file("MLST_${pairId}*.txt")
	
    script:
    """
    srst2.py --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 --output MLST_$pairId\
             --mlst_db ${params.mlst_db} --mlst_definitions ${params.mlst_def} --mlst_delimiter "_"\
             --max_unaligned_overlap 75 --threads 2 --log
    """
}