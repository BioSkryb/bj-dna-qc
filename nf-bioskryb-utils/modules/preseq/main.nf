nextflow.enable.dsl=2
params.timestamp = ""

include { PRESEQ_BAM2MR } from './bam2mr/main.nf'
include { PRESEQ_GC_EXTRAP } from './gcextrap/main.nf'

workflow PRESEQ_SUBWF{
    take:
        ch_bam
        ch_type
        ch_publish_dir
        ch_enable_publish
    main:
        PRESEQ_BAM2MR ( ch_bam,
                        ch_publish_dir,
                        ch_enable_publish
                      )

        PRESEQ_GC_EXTRAP ( 
                        PRESEQ_BAM2MR.out.mr,
                        ch_type,
                        ch_publish_dir,
                        ch_enable_publish
                     )
    
    
    emit:
        mr = PRESEQ_BAM2MR.out.mr
        coverage = PRESEQ_GC_EXTRAP.out.coverage
        preseqBam2mr_version = PRESEQ_BAM2MR.out.version
        preseqExtrap_version = PRESEQ_GC_EXTRAP.out.version
        
}