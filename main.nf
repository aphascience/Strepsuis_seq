#!/usr/bin/env nextflow

nextflow.enable.dsl=2


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