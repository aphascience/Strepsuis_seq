#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process srst2_recN {
tag "$pairId"
    errorStrategy 'ignore'
    maxForks 1
    publishDir "${params.outdir}/recN", mode: "copy"

    input:
        tuple val(pairId), file(R1), file(R2)

    output:
	    tuple val(pairId), file("recN_${pairId}*_genes_*.txt"), emit: recN_genes
        tuple val(pairId), file("recN_*_fullgenes_*.txt"), optional: true, emit: recN_fullgenes
	
    script:
    """
    srst2.py --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 --output recN_$pairId\
             --gene_db ${params.recN_ref} --max_unaligned_overlap 75 --log
    """
}

process srst2_mlst {
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

process srst2_serotype {
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
    serotype.sh $R1 $R2 $pairId ${params.sero_db} ${params.sero_def}
    serotype=\$(awk 'FNR == 2 {print \$2}' *_serotype__results.txt)
    """
}

process srst2_virulence {
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

process verifySerotype {
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

process combineoutputs {
tag "combine"
    publishDir "${params.outdir}/CombinedResults", mode: "copy"

    input:
        file("recN.tsv")
        file("MLST.tsv")
        file("serotype.tsv")
        file("virulence.tsv")
        file("verify.csv")

    output:
        file("*.csv")

    script:
    """
        combine_outputs.py recN.tsv MLST.tsv serotype.tsv virulence.tsv verify.csv ${params.dataDir}
    """
}

workflow {
    Channel
        .fromFilePairs( params.reads, flat: true )
            .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
            .set { readPairs }

    //trim(readPairs)

    srst2_recN(readPairs)
    srst2_mlst(readPairs)
    srst2_serotype(readPairs)
    srst2_virulence(readPairs)
    srst2_serotype.out.srst2_sero
        .join( readPairs )
        .set { verify_reads }
    verifySerotype(verify_reads)
    srst2_recN.out.recN_genes
        .map { it[1] }
        .collectFile( name: "recN_table.tsv", sort: true, keepHeader: true, storeDir: "${params.outdir}" )
        .set { recN }
    srst2_mlst.out
        .map { it[1] }
        .collectFile( name: "mlst_table.tsv", sort: true, keepHeader: true, skip: 1, storeDir: "${params.outdir}" )
        .set { mlst }
    srst2_serotype.out.serotype_results
        .map { it[1] }
        .collectFile( name: "serotype_table.tsv", sort: true, keepHeader: true, skip: 1, storeDir: "${params.outdir}" )
        .set { serotype }
    srst2_virulence.out.virulence_genes
        .map { it[1] }
        .collectFile( name: "virulence_table.tsv", sort: true, keepHeader: true, skip: 1, storeDir: "${params.outdir}" )
        .set { virulence }
    verifySerotype.out
        .map { it[1] }
        .collectFile( name: "verification.csv", sort: true, keepHeader: true, skip: 1, storeDir: "${params.outdir}" )
        .set { verify }
    combineoutputs(recN, mlst, serotype, virulence, verify)
}