nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_FASTQ_MERGE {
    tag "${sample_name}"
    publishDir "${params.publish_dir}_${params.timestamp}/primary_analyses/output/${sample_name}", enabled:"$enable_publish"

    input:
    tuple val(sample_name), val(sample_number), val(lane_number), path(reads)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), path("*LMERGE*.fastq.gz"), emit: reads
    path("*LMERGE*.fastq.gz"), emit: merged_fastqs_to_compress
    
    script:
    
    """
    # echo ${sample_number.join('-')}
    # echo ${lane_number.join('-')}
    
    cat *R1_001* > ${sample_name}_S1_LMERGE_R1_001.fastq.gz &
    cat *R2_001* > ${sample_name}_S1_LMERGE_R2_001.fastq.gz
    wait
    
    """
}

workflow CUSTOM_FASTQ_MERGE_WF{
    take:
        ch_publish_dir
        ch_enable_publish
        ch_reads
        
    main:
        CUSTOM_FASTQ_MERGE ( 
                            ch_reads,
                            ch_publish_dir,
                            ch_enable_publish
                           )
                           
    emit:
        reads = CUSTOM_FASTQ_MERGE.out.reads
        merged_fastqs_to_compress = CUSTOM_FASTQ_MERGE.out.merged_fastqs_to_compress
    
}

workflow{
    
    ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
                            .map { tag, pair -> subtags = ( tag =~ /(.*)_(S\d+)_(L+\d+)/)[0]; [ subtags[1], subtags[2], subtags[3], pair ] }
                            .groupTuple()
                            .map { tag, sample, lane, pair -> [ tag, sample.flatten(), lane.flatten(), pair.flatten() ] }
                            
    CUSTOM_FASTQ_MERGE_WF(
                            params.publish_dir,
                            params.enable_publish,
                            ch_reads
                        )
}