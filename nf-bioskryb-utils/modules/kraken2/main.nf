nextflow.enable.dsl=2
params.timestamp = ""

process KRAKEN2 {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/primary_analyses/metrics/${sample_name}/kraken2/", enabled:"$enable_publish"
    
    
    input:
    tuple val(sample_name), path(reads)
    path krakendb
    val(publish_dir)
    val(enable_publish)
    
    output:
    path("*unclassified*"), emit: unclassified
    path("*_kraken2_report.txt"), emit: report
    path("kraken2_version.yml"), emit: version
    
    script:
    """
    kraken2 \\
        --db ${krakendb} \\
        --threads $task.cpus \\
        --gzip-compressed \\
        --report ${sample_name}_kraken2_report.txt \\
        --output ${sample_name}_kraken2.out \\
        --unclassified-out ${sample_name}_unclassified_#.fastq \\
        --output ${sample_name}_kraken.out \\
        --paired ${reads}
        
    export KRAKEN2_VER=\$(kraken2 --version 2>&1 | sed  -n -e '1p' -e '/:7777/p' | grep -Eo [a-z]*[.]*[0-9][0-9]*.*)
    echo Kraken2: \$KRAKEN2_VER > kraken2_version.yml
    """
}

workflow KRAKEN2_WF{
    take:
        ch_reads
        ch_krakendb
        ch_publish_dir
        ch_enable_publish
    main:
        KRAKEN2 (
                    ch_reads,
                    ch_krakendb,
                    ch_publish_dir,
                    ch_enable_publish
                 )
    emit:
        unclassified = KRAKEN2.out.unclassified
        report = KRAKEN2.out.report
        version = KRAKEN2.out.version
}

workflow{
    ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
                            .map { tag, pair -> subtags = (tag =~ /(.*)_(S\d+)_(L0+\d+)/)[0]; [subtags[1], pair] }
    KRAKEN2_WF (
                    ch_reads,
                    params.krakendb,
                    params.publish_dir,
                    params.enable_publish
                  )
}

