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
    maxForks 2
    publishDir "${params.outdir}/recN", mode: "copy"

    input:
        tuple val(pairId), file(R1), file(R2)

    output:
	    tuple val(pairId), file("recN_${pairId}*_genes_*.txt"), emit: recN_genes
        tuple val(pairId), file("recN_*_fullgenes_*.txt"), optional: true, emit: recN_fullgenes
	
    script:
    """
    python3 ~/srst2/scripts/srst2.py --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 --output recN_$pairId --gene_db $projectDir/assets/recN/Ss-recN.fas --log
    """
}

process srst2_mlst {
tag "$pairId"
    publishDir "${params.outdir}/MLST", mode: "copy"

    input:
        tuple val(pairId), file(R1), file(R2)

    output:
	    tuple val(pairId), file("MLST_${pairId}*.txt")
	
    script:
    """
    python3 ~/srst2/scripts/srst2.py --threads 2 --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 --output MLST_$pairId --mlst_db $projectDir/assets/MLST/Ssuis_MLSTalleles.fas --mlst_definitions $projectDir/assets/MLST/Ssuis_Profiles.txt --mlst_delimiter "_" --log
    """
}

process srst2_serotype {
tag "$pairId"
    maxForks 2
    publishDir "${params.outdir}/serotype", mode: "copy"

    input:
        tuple val(pairId), file(R1), file(R2)

    output:
	    tuple val(pairId), file("serotype_${pairId}_*.txt"), emit: serotype
	
    script:
    """
    python3 ~/srst2/scripts/srst2.py --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 --output serotype_$pairId --mlst_db $projectDir/assets/Serotype/Ssuis_serotype.fas --mlst_definitions $projectDir/assets/Serotype/Ssuis_serotype_def.txt --mlst_delimiter "-" --log
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
    python3 ~/srst2/scripts/srst2.py --input_pe $R1 $R2 --forward _S.*_R1_001 --reverse _S.*_R2_001 --output virulence_$pairId --gene_db $projectDir/assets/Virulence/Ssuis_virulence.fas --log
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

    output:
        file("*.csv")

    script:
    """
        combine_outputs.py recN.tsv MLST.tsv serotype.tsv virulence.tsv
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
    srst2_recN.out.recN_genes
        .map { it[1] }
        .collectFile( name: "recN_table.tsv", sort: true, keepHeader: true, storeDir: "${params.outdir}")
        .set { recN }
    srst2_mlst.out
        .map { it[1] }
        .collectFile(  name: "mlst_table.tsv", sort: true, keepHeader: true, skip: 1, storeDir: "${params.outdir}")
        .set { mlst }
    srst2_serotype.out.serotype
        .map { it[1] }
        .collectFile(  name: "serotype_table.tsv", sort: true, keepHeader: true, skip: 1, storeDir: "${params.outdir}")
        .set { serotype }
    srst2_virulence.out.virulence_genes
        .map { it[1] }
        .collectFile(  name: "virulence_table.tsv", sort: true, keepHeader: true, skip: 1, storeDir: "${params.outdir}")
        .set { virulence }
    combineoutputs(recN, mlst, serotype, virulence)
}