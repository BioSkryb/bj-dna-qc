nextflow.enable.dsl=2
params.timestamp = ""

process SENTIEON_DRIVER_DEDUP {
    tag "${sample_name}"
    publishDir "${params.publish_dir}_${params.timestamp}/secondary_analyses/alignment", enabled:"$enable_publish", pattern: "output/${sample_name}.bam*"

    input:
    tuple val(sample_name), path(bams), path(bais), path(locuscollector_score)
    path fasta_ref
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), path("output/${sample_name}.bam"), path("output/${sample_name}.bam.bai"), emit: bam
    tuple val(sample_name), path("*.dedup_sentieonmetrics.txt"), emit: metrics
    path("*.dedup_sentieonmetrics.txt"), emit: mqc_metrics
    path("output/${sample_name}.bam"), emit: bam_to_compress
    path("sentieon_driver_dedup_version.yml"), emit: version

    script:

    """
    
    set +u
    
    if [ \$LOCAL != "true" ]; then
        . /opt/sentieon/cloud_auth.sh no-op
    else
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
        echo \$SENTIEON_LICENSE
    fi
   
    mkdir output;

    sentieon driver -t $task.cpus -r ${fasta_ref}/genome.fa -i ${bams.join(" -i ")} \
           --algo Dedup --score_info ${locuscollector_score[0]} --output_dup_read_name --metrics ${sample_name}.dedup_sentieonmetrics.txt ${sample_name}.tmp_dup_qname.txt.gz
    
    
    sentieon driver -t $task.cpus -r ${fasta_ref}/genome.fa -i ${bams.join(" -i ")} \
           --algo Dedup --rmdup --dup_read_name ${sample_name}.tmp_dup_qname.txt.gz output/${sample_name}.bam
    
    export SENTIEON_VER="202308.01"
    echo Sentieon: \$SENTIEON_VER > sentieon_driver_dedup_version.yml
    """
}

workflow SENTIEON_DRIVER_DEDUP_WF{
    take:
        ch_combine_outputs
        ch_reference
        ch_publish_dir
        ch_enable_publish
        
    main:
        SENTIEON_DRIVER_DEDUP ( 
                                ch_combine_outputs,
                                ch_reference,
                                ch_publish_dir,
                                ch_enable_publish
                              )
    emit:
        bam = SENTIEON_DRIVER_DEDUP.out.bam
        metrics = SENTIEON_DRIVER_DEDUP.out.metrics
        mqc_metrics = SENTIEON_DRIVER_DEDUP.out.mqc_metrics
        bam_to_compress = SENTIEON_DRIVER_DEDUP.out.bam_to_compress
        version = SENTIEON_DRIVER_DEDUP.out.version
    
}

workflow{
    
    ch_bam = Channel.fromFilePairs(params.bam_dir + "/*{.bam,.bam.bai}")
    ch_locuscollector_score = Channel.fromPath(params.locus_collector)
    combine_outputs = ch_bam.join(ch_locuscollector_score)
    
    SENTIEON_DRIVER_DEDUP_WF ( 
                                combine_outputs,
                                params.reference,
                                params.publish_dir,
                                params.enable_publish
                             )
    
}
