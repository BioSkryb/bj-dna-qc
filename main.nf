nextflow.enable.dsl=2

include { printHeader; helpMessage } from './help' params ( params )
include { PRESEQ_WF } from './workflow/preseq.nf' params ( params )


if ( params.help ) {
    helpMessage()
    exit 0
}

params.reference     = params.genomes [ params.genome ] [ 'reference' ]
params.intervals     = params.genomes [ params.genome ] [ 'base_metrics_intervals' ]
params.ginko_ref_dir = params.genomes [ params.genome ] [ "${params.bin_size}" ] [ 'ginko_ref_dir' ]
params.krakendb      = params.genomes [ params.genome ] [ 'kraken2_db' ]

workflow {
    
    printHeader()
    
    // Check if publishDir specified
    
    if ( params.publish_dir == "") {
        exit 1, "ERROR: publish_dir is not defined.\nPlease add --publish_dir s3://<bucket_name>/<project_name> to specify where the pipeline outputs will be stored."
    }
    
    if ( params.reads ) {
        
        ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )    
        ch_reads.ifEmpty{ exit 1, "ERROR: cannot find any fastq files matching the pattern: ${params.reads}\nMake sure that the input file exists!" }

    } else if ( params.input_csv ) {
        
        ch_reads = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                            .map { row -> [ row.biosampleName, [ row.read1, row.read2 ] ] }
        ch_reads.ifEmpty{ exit 1, "ERROR: Input csv file is empty." }
    }
    
    

    ch_reads.view()
    ch_dummy_file = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true).collect()
    ch_dummy_file2 = Channel.fromPath("$projectDir/assets/dummy_file2.txt", checkIfExists: true).collect()
    
    ch_reference =  Channel.fromPath(params.reference)
    ch_intervals =  Channel.fromPath(params.intervals)
    
    ch_binref = Channel.fromPath(params.ginko_ref_dir + "variable_" + params.bin_size + "_" + params.ginko_readlen + "_bwa")
    ch_gcref = Channel.fromPath(params.ginko_ref_dir + "GC_variable_" + params.bin_size + "_" + params.ginko_readlen + "_bwa")
    ch_boundsref_file = Channel.fromPath(params.ginko_ref_dir + "bounds_variable_" + params.bin_size + "_" + params.ginko_readlen + "_bwa")

       
    PRESEQ_WF( 
                params.publish_dir,
                params.enable_publish,
                params.disable_publish,
                ch_reads,
                params.input_csv,
                ch_dummy_file,
                ch_dummy_file2,
                params.mode,
                params.multiqc_config,
                params.project,
                params.bin_size,
                params.min_ploidy,
                params.max_ploidy,
                params.min_bin_width,
                params.is_haplotype,
                ch_reference,
                ch_intervals,
                ch_binref,
                ch_gcref,
                ch_boundsref_file,
                params.seqtk_sample_seed
             )
    
    
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}
