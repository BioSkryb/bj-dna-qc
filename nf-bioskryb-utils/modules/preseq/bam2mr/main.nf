nextflow.enable.dsl=2
params.timestamp = ""

process PRESEQ_BAM2MR {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/metrics/${sample_name}/preseq", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(bam), path(bai)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), file("*.mr"), emit: mr
    path("preseq_bam2mr_version.yml"), emit: version

    script:

    """
    set +u
    . /opt/sentieon/cloud_auth.sh no-op

    touch ${sample_name}_to-mr_unpairables
    
    bam2mr -seg_len 100000 -o ${sample_name}.mr ${bam} 2> ${sample_name}_to-unpairables.mr

    export PRESEQ_BAM2MR_VER=\$(preseq 2>&1 | grep -Eo '[0-9].[0-9].[0-9]')
    echo Preseq: \$PRESEQ_BAM2MR_VER > preseq_bam2mr_version.yml
    """
}

workflow PRESEQ_BAM2MR_WF{
    take:
        ch_bam
        ch_publish_dir
        ch_enable_publish
    main:
        PRESEQ_BAM2MR ( ch_bam,
                        ch_publish_dir,
                        ch_enable_publish
                      )
    emit:
        version = PRESEQ_BAM2MR.out.version
        mr = PRESEQ_BAM2MR.out.mr
}

workflow{
    ch_bam = Channel.fromFilePairs(params.bam_dir + "/*{.bam,bam.bai}")
    PRESEQ_BAM2MR_WF(
                        ch_bam,
                        params.publish_dir,
                        params.enable_publish
                     )
}