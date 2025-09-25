process VERIFY_SEROTYPE {
    tag "$pairId"
    publishDir "${params.outdir}/serotype", mode: "copy"
    input:
        tuple val(pairId), val(serotype), file(R1), file(R2)

    output:
        tuple val(pairId), file("*.csv"), optional: true

    script:
    """
    seroVerify.sh ${params.cps_ref} $R1 $R2 $pairId $serotype
    """
}
