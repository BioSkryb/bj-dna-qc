nextflow.enable.dsl=2
params.timestamp = ""


process QUALIMAP_BAMQC {
    tag "$sample_name"
    publishDir "${publish_dir}_${params.timestamp}/primary_analyses/metrics/${sample_name}/qualimap/", enabled:"$enable_publish"
    
    
    input:
    tuple val(sample_name), path(bam), path(bai)
    val(publish_dir)
    val(enable_publish)
    
    output:
    path("${task.process.split(':')[1]}_${sample_name}*"), emit: results
    path("qualimap_version.yml"), emit: version
    
    script:
    
    """
    unset DISPLAY
    mkdir tmp
    export _JAVA_OPTIONS="-Djava.io.tmpdir=./tmp" 
    #JAVA_MEM_SIZE=${task.memory.toGiga()}G
    #java_options="-Djava.awt.headless=true -Xmx\$JAVA_MEM_SIZE -XX:MaxPermSize=1024m"
    #java_options="-Djava.awt.headless=true -XX:MaxPermSize=1024m"
    #qualimap --java-mem-size=${task.memory.toGiga()}G
    qualimap \\
        --java-mem-size=${task.memory.toGiga()-1}G \\
        bamqc \\
        -bam ${bam} \\
        -p 'non-strand-specific' \\
        --collect-overlap-pairs \\
        -outdir ${task.process.split(':')[1]}_${sample_name} \\
        -nt $task.cpus
        
    echo QualiMap: \$(qualimap --help | sed  -n -e '4p' -e '/:7777/p' | grep -Eo [a-z]*.[0-9][0-9]*.*) > qualimap_version.yml
    """
}




workflow QUALIMAP_BAMQC_WF{
    take:
        ch_bam
        ch_publish_dir
        ch_enable_publish
    main:
        QUALIMAP_BAMQC ( ch_bam,
                         ch_publish_dir,
                         ch_enable_publish
                       )
    emit:
        results = QUALIMAP_BAMQC.out.results
        version = QUALIMAP_BAMQC.out.version
}

workflow {

    if (params.bam != "") {

        ch_bam_raw = Channel.fromFilePairs(params.bam, size: -1)
        ch_bam_raw
                  .map{ it -> it.flatten().collect() }
                  .set{ ch_bam }
    } else if(params.input_csv != "") {
        ch_bam = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                                .map { row -> [ row.sampleId, row.bam, row.bam + ".bai"  ] }
    }
    ch_bam.view()
    ch_bam.ifEmpty{ exit 1, "ERROR: No BAM files specified either via --bam or --input_csv" }

    QUALIMAP_BAMQC_WF(
                        ch_bam,
                        params.publish_dir,
                        params.enable_publish
                     )
}
