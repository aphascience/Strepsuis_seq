#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process trim {
    
}

process srst2_recN {

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