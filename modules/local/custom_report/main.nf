params.publish_dir = ""
params.timestamp = ""

process CUSTOM_REPORT {
    tag 'report'
    publishDir "${params.publish_dir}_${params.timestamp}/secondary_analyses/metrics", enabled:"$enable_publish"

    
    input:
    path(files)
    path(input_csv)
    val(publish_dir)
    val(enable_publish)


    output:
    path "nf-preseq-pipeline*", emit: mqc
    path("custom_report_version.yml"), emit: version

    
    script:
    """
    python3 /scripts/custom_report.py -s '*sentieonmetrics.txt' -m '*MAPD' -q '*no_qc_fastp.json' -t '*sample_metadata.json' -l '*lorenzstats.tsv' -i '$input_csv' -o 'Summary'


    echo custom_report: v0.0.1 > custom_report_version.yml
    
    """
    
}

workflow CUSTOM_REPORT_WF{
    take:
        ch_metrics
        ch_input_csv
        ch_publish_dir
        ch_enable_publish
    main:
        CUSTOM_REPORT ( 
                        ch_metrics,
                        ch_input_csv,
                        ch_publish_dir,
                        ch_enable_publish
                      )
    emit:
        mqc = CUSTOM_REPORT.out.mqc
        version = CUSTOM_REPORT.out.version
}

workflow{
    
    ch_files = Channel.fromPath(params.metrics_file)
    ch_input_csv = Channel.fromPath("NO_FILE")
    
    CUSTOM_REPORT_WF ( ch_files, ch_input_csv, params.publish_dir, params.enable_publish )
    
}
