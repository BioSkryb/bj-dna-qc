nextflow.enable.dsl=2
params.timestamp = ""

process FASTQC {
    tag "${sample_name}"
    publishDir "${params.publish_dir}_${params.timestamp}/primary_analyses/metrics/${sample_name}/fastqc/", enabled:params.publish_dir
    
    
    input:
    tuple val(sample_name), path(reads)
    val(publish_dir)
    val(enable_publish)

    
    output:
    path("${task.process.split(':')[1]}_*"), emit: report
    path("fastqc_version.yml"), emit: version
    
    script:
    """
    export _JAVA_OPTIONS="-XX:+ExitOnOutOfMemoryError -Xmx${task.memory.toGiga()-1}G"
    fastqc -f fastq --threads $task.cpus ${reads}
    for file in *_fastqc.{zip,html}; do echo \$file;mv {,${task.process.split(':')[1]}_}\$file;done

    export FASTQC_VER=\$(fastqc -v 2>&1 | grep -o 'v[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+')
    echo FastQC: \$FASTQC_VER > fastqc_version.yml
    """
}


workflow FASTQC_WF{
    take:
        ch_reads
        ch_publish_dir
        ch_enable_publish
    main:
        FASTQC ( ch_reads,
                         ch_publish_dir,
                         ch_enable_publish
                       )
    emit:
        results = FASTQC.out.report
        version = FASTQC.out.version
}

workflow{
    ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
    ch_reads.view()
    FASTQC_WF(
                        ch_reads,
                        params.publish_dir,
                        params.enable_publish
                     )
}
