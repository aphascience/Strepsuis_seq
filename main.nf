#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* process trim {
    tag "$pairId"
    input:
        tuple val(pairId), path(raw1), path(raw2)
    output:
        tuple val(pairId), path("*_trimmed_R1.fastq.gz"), path("*_trimmed_R2.fastq.gz")
    script:
    """
    trim.sh $pairId $raw1 $raw2 ${params.adapters}
    """
}

process bowtie_map {
    tag "$pairId"
    input:
        tuple val(pairId), path("*_trimmed_R1.fastq.gz"), path("*_trimmed_R2.fastq")
        path($projectDir/assets/Ss-recN)
    output:
        tuple val(pairId), path("map.sam")
    """
    bowtie_map.sh trimmed_R1.fastq.gz trimmed_R2.fastq.gz $projectDir/assets/Ss-recN
    """
}

process pileup {
    tag "$pairId"
    input:
        tuple val(pairId), path("map.sam")
    output:
        tuple val(pairId), path(output.pileup)
    """
    pileup.sh map.sam $projectDir/assets/Ss-recN.fas
    """
    
}*/

process srst2_recN {
tag "$pairId"
    maxForks 1
    publishDir "${params.outdir}/recN", mode: "copy"

    input:
        tuple val(pairId), file(R1), file(R2)

    output:
	    tuple val(pairId), file("recN_${pairId}*_genes_*.txt"), emit: recN_genes
        tuple val(pairId), file("recN_*_fullgenes_*.txt"), optional: true, emit: recN_fullgenes
	
    script:
    """
    python3 srst2.py --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 --output recN_$pairId --gene_db $projectDir/assets/recN/Ss-recN.fas --log
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
    python3 srst2.py --threads 2 --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 --output MLST_$pairId --mlst_db $projectDir/assets/MLST/Ssuis_MLSTalleles.fas --mlst_definitions $projectDir/assets/MLST/Ssuis_Profiles.txt --mlst_delimiter "_" --log
    """
}

process srst2_serotype {
tag "$pairId"
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
    maxForks 2
    publishDir "${params.outdir}/virulence", mode: "copy"

    input:
        tuple val(pairId), file(R1), file(R2)

    output:
	    tuple val(pairId), file("virulence_${pairId}*_genes_*.txt"), emit: virulence_genes
        tuple val(pairId), file("virulence_*_fullgenes_*.txt"), optional: true, emit: virulence_fullgenes
	
    script:
    """
    python3 srst2.py --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 --output virulence_$pairId --gene_db $projectDir/assets/Virulence/Ssuis_virulence.fas --log
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
    seroVerify.sh $projectDir/assets/Serotype/cps2K $R1 $R2 $pairId $serotype
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
        combine_outputs.py recN.tsv MLST.tsv serotype.tsv virulence.tsv verify.csv
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