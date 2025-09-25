process SRST2_VIRULENCE {
    tag "$pairId"
    errorStrategy 'ignore'
    maxForks 2
    publishDir "${params.outdir}/virulence", mode: "copy"

    input:
        tuple val(pairId), file(R1), file(R2)

    output:
	    tuple val(pairId), file("virulence_${pairId}*_genes_*.txt"), emit: virulence_genes
        tuple val(pairId), file("virulence_*_fullgenes_*.txt"), optional: true, emit: virulence_fullgenes
	
    script:
    """
    srst2.py --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 --output virulence_$pairId\
             --gene_db ${params.virulence_ref} --max_unaligned_overlap 75 --log\
             || true
    """
}