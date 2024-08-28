nextflow.enable.dsl=2
params.timestamp = ""

process FASTP {
    tag "${sample_name}"
    publishDir "${params.publish_dir}_${params.timestamp}/primary_analyses/metrics/${sample_name}/fastp/", enabled:"$enable_publish", pattern: "*fastp.json"

    input:
    tuple val(sample_name), path(reads)
    val(two_color_chemistry)
    val(qc_only_mode)
    val(adapter_sequence)
    val(adapter_sequence_r2)
    val(publish_dir)
    val(enable_publish)
    
    when:
    reads[0].size() > 10.KB

    output:
    tuple val(sample_name), path("*trim.fastq.gz"), emit: reads
    path("*fastp.json"), emit: report
    path("fastp_version.yml"), emit: version

    script:
  
    opts = []
    
    if (reads.size() == 2){ 
     
     if(!qc_only_mode){
        opts.add("--in1 ${reads[0]}")
        opts.add("--in2 ${reads[1]}") 
        opts.add("--out1 ${sample_name}_R1_trim.fastq.gz")
        opts.add("--out2 ${sample_name}_R2_trim.fastq.gz")
        opts.add("--json ${sample_name}_no_qc_fastp.json")
        opts.add("--html ${sample_name}_no_qc_fastp.html")
        if (two_color_chemistry){
            opts.add("--trim_poly_g")
            opts.add("--poly_g_min_len 10") 
            }
            
        if( adapter_sequence == "none" ){
        
            opts.add("--detect_adapter_for_pe")

            }else{
        
            opts.add("--adapter_sequence=${adapter_sequence}")
            opts.add("--adapter_sequence_r2=${adapter_sequence_r2}")
            
            }
        }

     else{
        opts.add("--in1 ${reads[0]}")
        opts.add("--in2 ${reads[1]}") 
        opts.add("-A")
        opts.add("-G")
        opts.add("-Q")
        opts.add("-L")
        opts.add("--out1 ${sample_name}_qc_only_R1_trim.fastq.gz")
        opts.add("--out2 ${sample_name}_qc_only_R2_trim.fastq.gz")
        opts.add("--json ${sample_name}_qc_only_fastp.json")
        opts.add("--html ${sample_name}_qc_only_fastp.html")
        
        if( adapter_sequence == "none" ){
        
            opts.add("--detect_adapter_for_pe")

        }else{
        
            opts.add("--adapter_sequence=${adapter_sequence}")
            opts.add("--adapter_sequence_r2=${adapter_sequence_r2}")
            
            }
    
        }
    }
    
   else{ 
   
     if(!qc_only_mode){
        opts.add("--in1 ${reads[0]}")
        opts.add("--out1 ${sample_name}_R1_trim.fastq.gz")
        opts.add("--json ${sample_name}_no_qc_fastp.json")
        opts.add("--html ${sample_name}_no_qc_fastp.html")
        if (two_color_chemistry){
            opts.add("--trim_poly_g")
            opts.add("--poly_g_min_len 10") 
            }
            
        if( adapter_sequence != "none" ){
        
            opts.add("--adapter_sequence=${adapter_sequence}")

            }
        }
        
     else{
        opts.add("--in1 ${reads[0]}")
        opts.add("-A")
        opts.add("-G")
        opts.add("-Q")
        opts.add("-L")
        opts.add("--out1 ${sample_name}_qc_only_R1_trim.fastq.gz")
        opts.add("--json ${sample_name}_qc_only_fastp.json")
        opts.add("--html ${sample_name}_qc_only_fastp.html")
        
        if( adapter_sequence != "none" ){
        
            opts.add("--adapter_sequence=${adapter_sequence}")

            }
    
        }
    
   }
    
    """
    df -h
    
    fastp \
        --thread $task.cpus \
        ${opts.join(' ')} \

        2> ${sample_name}_fastp.log
       
        
    export FASTP_VER=\$(fastp --version 2>&1 | sed -e "s/fastp //g")
    echo fastp: \$FASTP_VER > fastp_version.yml
    """
    
}

workflow FastpQCWF{
    take:
        ch_reads
        ch_two_colors_chemistry
        ch_adapter_sequence
        ch_adapter_sequence_r2
        ch_publish_dir
        ch_enable_publish

    
    main:
        FASTP ( 
                ch_reads, 
                ch_two_colors_chemistry,
                true,
                ch_adapter_sequence,
                ch_adapter_sequence_r2,
                ch_publish_dir,
                ch_enable_publish
              )
              
    emit:
        report = FASTP.out.report
    
}

workflow FastpNoQCWF{
    take:
        ch_reads
        ch_two_colors_chemistry
        ch_adapter_sequence
        ch_adapter_sequence_r2
        ch_publish_dir
        ch_enable_publish

    main:
        FASTP ( 
                ch_reads, 
                ch_two_colors_chemistry,
                false,
                ch_adapter_sequence,
                ch_adapter_sequence_r2,
                ch_publish_dir,
                ch_enable_publish
              )
    emit:
        reads = FASTP.out.reads
        version = FASTP.out.version
        report = FASTP.out.report
    
}
workflow FastpFull_WF{
    take:
        ch_reads
        ch_two_colors_chemistry
        ch_adapter_sequence
        ch_adapter_sequence_r2
        ch_publish_dir
        ch_enable_publish

    main:
        FASTP ( 
                ch_reads, 
                ch_two_colors_chemistry,
                false,
                ch_adapter_sequence,
                ch_adapter_sequence_r2,
                ch_publish_dir,
                ch_enable_publish
              )
    emit:
        reads = FASTP.out.reads
        version = FASTP.out.version
        report = FASTP.out.report
    
}

workflow{
    ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
                            .map { tag, pair -> subtags = (tag =~ /(.*)_(S\d+)_(L0+\d+)/)[0]; [subtags[1], pair] }
                            
    FastpWF(
            ch_reads,
            params.two_color_chemistry,
            params.adapter_sequence,
            params.adapter_sequence_r2,
            params.publish_dir,
            params.enable_publish
           )
}