nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_DATA_PROCESSING {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    tuple val(sample_name), path(bam), path(bai), path(dedup_metrics), path(metrics), path(preseq_cov)
    path fasta_ref
    val(publish_dir)
    val(enable_publish)


    output:
    // tuple val(sample_name), path("*_fastqqc_mqc.txt"), emit: metrics_tuple
    path "*_all_sentieonmetrics.txt", emit: metrics
    path("custom_processing_version.yml"), emit: version

    script:
    
    """
    #! /bin/bash
    set +u
    . /opt/sentieon/cloud_auth.sh no-op

    
    cat ${sample_name}.cov_sentieonmetrics.sample_interval_summary | tail -n+2 | sed 's|:|\t|g' | sed 's|-|\t|g' | cut -f1-4 > ${sample_name}_chromosome_read_proportions.tsv
    
    cat ${sample_name}.dedup.alignmentstat_sentieonmetrics.txt | grep "^PAIR" | cut -f 7 > ${sample_name}_pct_aligned
    
    cat ${sample_name}.dedup.alignmentstat_sentieonmetrics.txt | grep "^PAIR" | cut -f 6 > ${sample_name}_n_pe_trimmed_aligned
    
    
    
    samtools view -F 3854 ${bam} | awk '(\$7 ~ /=/ && \$9 < 10000) {print 1}' |wc -l > ${sample_name}_n_input_phi29_chimeras
    
    samtools view -c ${bam}  > ${sample_name}_preseq_read_count


    
    PRESEQ_COUNT=\$(cat ${sample_name}_preseq_cov.txt)
    PRESEQ_TOTAL=\$(cat ${sample_name}_preseq_read_count)
    
    N_PE_TRIMMED_ALIGNED_READS=\$(cat ${sample_name}_n_pe_trimmed_aligned)
    PHI29_CHIMERAS=\$(cat ${sample_name}_n_input_phi29_chimeras)
    
    
    parse_metrics_files.py \
      ${sample_name}.wgsmetricsalgo.sentieonmetrics.txt \
      ${sample_name}.dedup_sentieonmetrics.txt \
      ${sample_name}.dedup.alignmentstat_sentieonmetrics.txt \
      ${sample_name}_chromosome_read_proportions.tsv \
      ${sample_name}_pct_aligned \
      ${sample_name}.gcbias_summary.sentieonmetrics.txt \
      ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt \
      \$PRESEQ_COUNT \
      ${sample_name} \
      \$N_PE_TRIMMED_ALIGNED_READS \
      \$PHI29_CHIMERAS \
      \$PRESEQ_TOTAL
      
    mv ${sample_name}_all_metrics.tsv ${sample_name}_all_sentieonmetrics.txt

    echo custom_processing: v0.0.1 > custom_version.yml
    export SAMTOOLS_VER=\$(samtools --version 2>&1 |  sed -n -e '1p' | grep -Eo [0-9][.]*[0-9]*)
    echo samtools: \$SAMTOOLS_VER > samtools_version.yml
    cat custom_version.yml samtools_version.yml > custom_processing_version.yml
    """
}

workflow CUSTOM_DATA_PROCESSING_WF{
    take:
        ch_combine_outputs
        ch_reference
        ch_publish_dir
        ch_enable_publish
    main:
        CUSTOM_DATA_PROCESSING ( 
                                  ch_combine_outputs,
                                  ch_reference,
                                  ch_publish_dir,
                                  ch_enable_publish
                                )
    emit:
        metrics = CUSTOM_DATA_PROCESSING.out.metrics
        version = CUSTOM_DATA_PROCESSING.out.version
}

// workflow{
    
//     combine_outputs_b = ch_bam
//                               .join(ch_metrics
//                                               .join(ch_metrics_tuple
//                                                                     .join(ch_coverage)
//                                                     )
//                                     )
//     CUSTOM_DATA_PROCESSING_WF(
//                                 ch_combine_outputs,
//                                 params.reference,
//                                 params.publish_dir,
//                                 params.enable_publish
//                              )   
// }

