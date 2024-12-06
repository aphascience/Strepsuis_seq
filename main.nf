#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process trim {
    tag "pairId"
    input:
        tuple val(pairId), path(raw1), path(raw2)
    output:
        tuple val(pairId), path("trimmed_R1.fastq.gz"), path("trimmed_R2.fastq")
    """
    trim.sh raw1 raw2 ${params.adapters}
    """
}

process bowtie_map {
    tag "pairId"
    input:
        tuple val(pairId), path("trimmed_R1.fastq.gz"), path("trimmed_R2.fastq")
        path($projectDir/assets/Ss-recN)
    output:
        tuple val(pairId), path("map.sam")
    """
    bowtie_map.sh trimmed_R1.fastq.gz trimmed_R2.fastq.gz $projectDir/assets/Ss-recN
    """
}

process pileup {
    tag "pairId"
    input:
        tuple val(pairId), path("map.sam")
    
    
}

process srst2_recN {
tag { "recN_${pairId}" }
    publishDir "${params.outdir}/recN", mode: "copy"

    input:
        tuple pairId, file(reads)

    output:
	file("${pairId}_srst2*")	
	
    script:      
    geneDB = params.gene_db ? "--gene_db $gene_db" : ''
    mlstDB = params.mlst_db ? "--mlst_db $mlst_db" : ''
    mlstdef = params.mlst_db ? "--mlst_definitions $mlst_definitions" : ''
    mlstdelim = params.mlst_db ? "--mlst_delimiter $params.mlst_delimiter" : ''
    """
    srst2 --input_pe $reads --output ${pairId}_recN --min_coverage $params.min_gene_cov \
          --max_divergence $params.max_gene_divergence $mlstDB $mlstdef $mlstdelim $geneDB 
    """
}

process srst2_mlst {

}

process srst2_serotype {

}

process srst2_virulence {

}

process combineoutputs {

}

workflow {
    Channel
        .fromFilePairs( params.reads, flat: true )
            .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	        .set { readPairs }

    process trim(readPairs)

    process srst2_recN(trim.out)
    process srst2_mlst(trim.out)
    process srst2_serotype(trim.out)
    process srst2_virulence(trim.out)
    process combineoutputs(srst2_recN.out, srst2_mlst.out, srst2_serotype.out, srst2_virulence.out)
}