nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_BAM_SUBSAMPLE {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/alignment/subsample", enabled:"$enable_publish"


    input:
    tuple val(sample_name), path(bam), path(bai), val(subsample)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val("${sample_name}_${subsample}"), path("*_${subsample}.bam"), path("*_${subsample}.bam.bai"), emit: bam
    tuple val("${sample_name}"), path("*_${subsample}.bam"), path("*_${subsample}.bam.bai"), emit: bam_with_samplename
    path("custom_bam_subsample_version.yml"), emit: version
    
    script:
    
    """
    #! /bin/bash
    set +u
    
    # Count properly aligned reads using best practices flags
    # -f 0x2: properly paired reads
    # -F excludes: secondary (0x100), QC-fail (0x200), duplicates (0x400), supplementary (0x800).
    export PROPERLY_ALIGNED=\$(samtools view -c -f 0x2 -F 0xF00 --threads $task.cpus ${bam})
    
    # Calculate fraction needed for subsampling
    export TARGET_READS=${subsample}
    frac=\$(awk -v t=\$TARGET_READS -v c=\$PROPERLY_ALIGNED 'BEGIN{p=t/c; if(p>1) p=1; printf("%.6f", p)}')
    
    echo "Target reads: \$TARGET_READS | Properly aligned reads: \$PROPERLY_ALIGNED | Fraction: \$frac"
    
    # Handle subsampling based on fraction
    if [ "\$frac" = "1.000000" ]; then
        # No subsampling needed - copy all properly aligned reads
        samtools view -@ $task.cpus -b -f 0x2 -F 0xF00 ${bam} -o ${sample_name}_${subsample}_unsorted.bam
    else
        # Subsample with reproducible seed
        # Remove leading "0." before concatenating with seed for reproducibility
        frac_digits=\${frac#0}
        samtools view -@ $task.cpus -b -f 0x2 -F 0xF00 -s 42\$frac_digits ${bam} -o ${sample_name}_${subsample}_unsorted.bam
    fi
    
    # Sort and index with optimal settings
    samtools sort -@ $task.cpus -l 6 -o ${sample_name}_${subsample}.bam ${sample_name}_${subsample}_unsorted.bam
    samtools index -@ $task.cpus ${sample_name}_${subsample}.bam
    
    # Clean up intermediate file
    rm ${sample_name}_${subsample}_unsorted.bam
    
    # Generate version file
    export SAMTOOLS_VER=\$(samtools --version 2>&1 | sed -n -e '1p' | grep -Eo '[0-9][.]*[0-9]*')
    echo "Samtools: \$SAMTOOLS_VER" > custom_bam_subsample_version.yml
    """
}

workflow CUSTOM_BAM_SUBSAMPLE_WF{
    take:
        ch_bam
        ch_publish_dir
        ch_enable_publish
        
    main:
        CUSTOM_BAM_SUBSAMPLE ( 
                                ch_bam,
                                ch_publish_dir,
                                ch_enable_publish
                             )
                           
    emit:
        bam = CUSTOM_BAM_SUBSAMPLE.out.bam
        version = CUSTOM_BAM_SUBSAMPLE.out.version
        
    
}

workflow{
    
    if (params.bam != "") {

        ch_bam_raw = Channel.fromFilePairs(params.bam, size: -1)
        ch_bam_raw
                  .map{ it -> it.flatten().collect() }
                  .set{ ch_bam }
    } else if(params.input_csv != "") {
        ch_bam = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                                .map { row -> [ row.biosampleName, row.bam, row.bam + ".bai"  ] }
    }
    ch_bam.view()
    ch_bam.ifEmpty{ exit 1, "ERROR: No BAM files specified either via --bam or --input_csv" }

    
    ch_subsample = Channel.fromList( [params.n_reads] )
    ch_bam_subsample = ch_bam.combine(ch_subsample)
                            
    CUSTOM_BAM_SUBSAMPLE_WF(
                            ch_bam_subsample,
                            params.publish_dir,
                            params.enable_publish
                           )
}
