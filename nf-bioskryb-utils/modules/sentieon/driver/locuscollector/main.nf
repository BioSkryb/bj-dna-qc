nextflow.enable.dsl=2
params.timestamp = ""

process SENTIEON_DRIVER_LOCUSCOLLECTOR {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(bams), path(bais)
    path fasta_ref
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), file("${sample_name}.locuscollector_score.gz*"), emit: locuscollector_score
    path("sentieon_driver_locuscollector_version.yml"), emit: version
    
    script:
    
    """
    set +u
        if [ \$NF_TEST != "true" ]; then
        . /opt/sentieon/cloud_auth.sh no-op
    else
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
        echo \$SENTIEON_LICENSE
    fi

    sentieon driver -t $task.cpus -r '${fasta_ref}/genome.fa' -i ${bams.join(" -i ")} \
       --algo LocusCollector --fun score_info '${sample_name}.locuscollector_score.gz'

    export SENTIEON_VER="202308.01"
    echo Sentieon: \$SENTIEON_VER > sentieon_driver_locuscollector_version.yml
    """
}

workflow SENTIEON_DRIVER_LOCUSCOLLECTOR_WF{
    take:
        ch_bam
        ch_reference
        ch_publish_dir
        ch_enable_publish
    
    main:
        SENTIEON_DRIVER_LOCUSCOLLECTOR ( 
                                            ch_bam,
                                            ch_reference,
                                            ch_publish_dir,
                                            ch_enable_publish
                                      )
    
    emit:
        locuscollector_score = SENTIEON_DRIVER_LOCUSCOLLECTOR.out.locuscollector_score
        version = SENTIEON_DRIVER_LOCUSCOLLECTOR.out.version
}

workflow{
    
    ch_bam = Channel.fromFilePairs(params.bam_dir + "/*{.bam,bam.bai}")
    
    SENTIEON_DRIVER_LOCUSCOLLECTOR_WF ( 
                            ch_bam,
                            params.reference,
                            params.publish_dir,
                            params.enable_publish
                         )
}
