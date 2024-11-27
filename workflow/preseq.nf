nextflow.enable.dsl=2
import groovy.json.JsonOutput
/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// MODULES
//

include { PUBLISH_INPUT_DATASET_WF } from '../nf-bioskryb-utils/modules/bioskryb/publish_input_dataset/main.nf' addParams(timestamp: params.timestamp)
include { CUSTOM_FASTQ_MERGE_WF } from '../nf-bioskryb-utils/modules/bioskryb/custom_fastq_merge/main.nf' addParams(timestamp: params.timestamp)
include { FastpNoQCWF } from '../nf-bioskryb-utils/modules/fastp/main.nf' addParams(timestamp: params.timestamp)
include { FastpQCWF } from '../nf-bioskryb-utils/modules/fastp/main.nf' addParams(timestamp: params.timestamp)
include { FASTQC } from '../nf-bioskryb-utils/modules/fastqc/main.nf' addParams(timestamp: params.timestamp)
include { SEQTK_WF } from '../nf-bioskryb-utils/modules/seqtk/sample/main.nf' addParams(timestamp: params.timestamp)
include { SENTIEON_BWA_WF } from '../nf-bioskryb-utils/modules/sentieon/bwa/mem/main.nf' addParams(timestamp: params.timestamp)
include { KRAKEN2_WF } from '../nf-bioskryb-utils/modules/kraken2/main.nf' addParams(timestamp: params.timestamp)
include { SENTIEON_DRIVER_LOCUSCOLLECTOR_WF } from '../nf-bioskryb-utils/modules/sentieon/driver/locuscollector/main.nf' addParams(timestamp: params.timestamp)
include { SENTIEON_DRIVER_DEDUP_WF } from '../nf-bioskryb-utils/modules/sentieon/driver/dedup/main.nf' addParams(timestamp: params.timestamp)
include { SENTIEON_DRIVER_METRICS_WF } from '../nf-bioskryb-utils/modules/sentieon/driver/metrics/main.nf' addParams(timestamp: params.timestamp)
include { QUALIMAP_BAMQC_WF } from '../nf-bioskryb-utils/modules/qualimap/bamqc/main.nf' addParams(timestamp: params.timestamp)
include { PRESEQ_BAM2MR_WF } from '../nf-bioskryb-utils/modules/preseq/bam2mr/main.nf' addParams(timestamp: params.timestamp)
include { PRESEQ_GC_EXTRAP_WF } from '../nf-bioskryb-utils/modules/preseq/gcextrap/main.nf' addParams(timestamp: params.timestamp)
include { CUSTOM_DATA_PROCESSING_WF } from '../modules/local/custom_data_processing/main.nf' addParams(timestamp: params.timestamp)
include { CUSTOM_CALCULATE_MAPD } from '../modules/local/custom_calculate_mapd/main.nf' addParams(timestamp: params.timestamp)
include { CUSTOM_REPORT_WF } from '../modules/local/custom_report/main.nf' addParams(timestamp: params.timestamp)
include { MULTIQC_WF } from '../nf-bioskryb-utils/modules/multiqc/main.nf' addParams(timestamp: params.timestamp)
include { REPORT_VERSIONS_WF } from '../nf-bioskryb-utils/modules/bioskryb/report_tool_versions/main.nf' addParams(timestamp: params.timestamp)
include { BAM_TO_BED_WF } from '../nf-bioskryb-utils/modules/ginkgo/bam_to_bed/main.nf' addParams(timestamp: params.timestamp)
include { GINKGO_BINUNSORT_WF } from '../nf-bioskryb-utils/modules/ginkgo/binunsort/main.nf' addParams(timestamp: params.timestamp)
include { GINKGO_SEGMENTATION_R_WF } from '../nf-bioskryb-utils/modules/ginkgo/segmentation_r/main.nf' addParams(timestamp: params.timestamp)
include { GINKGO_CNV_CALLER_WF } from '../nf-bioskryb-utils/modules/ginkgo/cnvcaller/main.nf' addParams(timestamp: params.timestamp)
include { GINKO_RDS_TO_FLAT } from '../nf-bioskryb-utils/modules/ginkgo/rds_to_flat/main.nf' addParams(timestamp: params.timestamp)
include { GINKO_PARSE_OUTPUTS } from '../nf-bioskryb-utils/modules/ginkgo/parse_ginko_outputs/main.nf' addParams(timestamp: params.timestamp)
include { PARSE_RDS_CNV_METRICS } from '../nf-bioskryb-utils/modules/ginkgo/parse_rds_cnv_metrics/main.nf' addParams(timestamp: params.timestamp)
include { BAM_LORENZ_COVERAGE_WF } from '../nf-bioskryb-utils/modules/bam_lorenz_coverage/main.nf' addParams(timestamp: params.timestamp)
include { COUNT_READS_FASTQ_WF } from '../nf-bioskryb-utils/modules/bioskryb/custom_read_counts/main.nf' addParams(timestamp: params.timestamp)


params.reference         = params.genomes [ params.genome ] [ 'reference' ]
params.intervals         = params.genomes [ params.genome ] [ 'base_metrics_intervals' ]
params.krakendb          = params.genomes [ params.genome ] [ 'kraken2_db' ]
params.blacklist_regions = params.genomes [ params.genome ] [ 'blacklist_regions' ]

workflow PRESEQ_WF {
    
    take:
        ch_publish_dir
        ch_enable_publish
        ch_disable_publish
        ch_reads
        ch_input_csv
        ch_dummy_file
        ch_dummy_file2
        ch_mode
        ch_multiqc_config
        ch_project_name
        ch_bin_size
        ch_min_ploidy
        ch_max_ploidy
        ch_min_bin_width
        ch_is_haplotype
        ch_reference
        ch_intervals
        ch_binref
        ch_gcref
        ch_boundsref_file
        ch_seqtk_sample_seed
        ch_min_reads
    
    main:
    
        if ( params.input_csv ) {
            
            PUBLISH_INPUT_DATASET_WF (
                                        ch_input_csv,
                                        ch_publish_dir,
                                        ch_enable_publish
                                    )
        } else {
            input_csv_dummy_file = file("input_csv_dummy_file.txt")
            ch_input_csv = Channel.fromPath( input_csv_dummy_file ).collect()
        }

        COUNT_READS_FASTQ_WF (
                                ch_reads,
                                params.publish_dir,
                                params.enable_publish
                            )
        COUNT_READS_FASTQ_WF.out.read_counts
        .map { sample_id, files, read_count -> 
            [sample_id, files, read_count.toInteger()]
        }
        .branch {
            small: it[2] < ch_min_reads
            large: it[2] >= ch_min_reads
        }
        .set { branched_reads }

        ch_reads_with_n_reads = branched_reads.large.map { biosampleName, reads, read_count ->
            return [ biosampleName, reads, params.n_reads ]
        }

        SEQTK_WF (
                        ch_reads_with_n_reads,
                        false,
                        params.read_length,
                        ch_seqtk_sample_seed,
                        ch_publish_dir,
                        ch_disable_publish
                     )
        ch_sample_metadata = SEQTK_WF.out.metadata
                     
        FastpNoQCWF ( 
                SEQTK_WF.out.reads, 
                params.two_color_chemistry,
                params.adapter_sequence,
                params.adapter_sequence_r2,
                ch_publish_dir,
                ch_enable_publish
              )
              
        ch_fastp_report = FastpNoQCWF.out.report
        
        FastpNoQCWF.out.reads.ifEmpty{ exit 1, "ERROR: subsample of fastq files resulted in fastq files < 10.KB in size." }
              
        ch_kraken2_report = Channel.empty()
        ch_kraken2_version = Channel.empty()

        if ( !params.skip_kraken ) {                         
    
            KRAKEN2_WF (
                        FastpNoQCWF.out.reads,
                        params.krakendb,
                        ch_publish_dir,
                        ch_enable_publish
                      )
                      
            ch_kraken2_report = KRAKEN2_WF.out.report
            ch_kraken2_version = KRAKEN2_WF.out.version
        }
        
        ch_fastqc_report = Channel.empty()
        ch_fastqc_version = Channel.empty()
        
        if ( !params.skip_fastqc && params.instrument != 'NovaSeq' ) {
            FASTQC (
                        FastpNoQCWF.out.reads,
                        ch_publish_dir,
                        ch_enable_publish
                   )
            ch_fastqc_report = FASTQC.out.report
            ch_fastqc_version = FASTQC.out.version
        }
        
        
        SENTIEON_BWA_WF ( 
                            FastpNoQCWF.out.reads,
                            params.reference,
                            ch_publish_dir,
                            ch_disable_publish
                         )
        
        SENTIEON_DRIVER_LOCUSCOLLECTOR_WF ( 
                                            SENTIEON_BWA_WF.out.bam,
                                            params.reference,
                                            ch_publish_dir,
                                            ch_disable_publish
                                      )
        
        combine_outputs_a = SENTIEON_BWA_WF.out.bam.join(SENTIEON_DRIVER_LOCUSCOLLECTOR_WF.out.locuscollector_score)
        
        SENTIEON_DRIVER_DEDUP_WF ( 
                                combine_outputs_a,
                                params.reference,
                                ch_publish_dir,
                                ch_enable_publish
                              )
        
        custom_output = SENTIEON_DRIVER_DEDUP_WF.out.bam.combine( ch_dummy_file )
        
        SENTIEON_DRIVER_DEDUP_WF.out.bam
            .collectFile( name: "bam_files.txt", newLine: true, sort: { it[0] }, storeDir: "${params.tmp_dir}" )
                { it[0] + "\t" + "${params.publish_dir}_${params.timestamp}/secondary_analyses/alignment/output/" + it[1].getName() }

        
        
        SENTIEON_DRIVER_METRICS_WF ( 
                                      custom_output,
                                      params.reference,
                                      params.intervals,
                                      ch_dummy_file2,
                                      ch_mode,
                                      ch_publish_dir,
                                      ch_enable_publish
                                   )
        
        ch_bam_lorenz_coverage_stats = Channel.empty()
        ch_bam_lorenz_coverage_version = Channel.empty()
        
        BAM_LORENZ_COVERAGE_WF (
                                SENTIEON_DRIVER_DEDUP_WF.out.bam,
                                ch_publish_dir,
                                ch_enable_publish
                            )
                            
        ch_bam_lorenz_coverage_stats = BAM_LORENZ_COVERAGE_WF.out.stats
        ch_bam_lorenz_coverage_version = BAM_LORENZ_COVERAGE_WF.out.version


        ch_qualimap_report = Channel.empty()
        ch_qualimap_version = Channel.empty()

        if ( !params.skip_qualimap ) {            
        
            QUALIMAP_BAMQC_WF ( 
                                SENTIEON_DRIVER_DEDUP_WF.out.bam,
                                ch_publish_dir,
                                ch_disable_publish
                              )
            ch_qualimap_report = QUALIMAP_BAMQC_WF.out.results
            ch_qualimap_version = QUALIMAP_BAMQC_WF.out.version
        }
        
        PRESEQ_BAM2MR_WF ( 
                            SENTIEON_DRIVER_DEDUP_WF.out.bam,
                            ch_publish_dir,
                            ch_disable_publish
                         )
        
        PRESEQ_GC_EXTRAP_WF ( 
                                PRESEQ_BAM2MR_WF.out.mr,
                                ch_publish_dir,
                                ch_disable_publish
                            )
        
        combine_outputs_b = SENTIEON_DRIVER_DEDUP_WF.out.bam
                                                    .join(SENTIEON_DRIVER_DEDUP_WF.out.metrics
                                                                                    .join(SENTIEON_DRIVER_METRICS_WF.out.metrics_tuple
                                                                                                                    .join(PRESEQ_GC_EXTRAP_WF.out.coverage)
                                                                                            )
                                                            )
        combine_outputs_b.view()
     
        
        
        // CNV CALLING = GINKGO
        ch_ginkgo_version = Channel.empty()
        ch_ginkgo_stats = Channel.empty()
        ch_bedtools_version = Channel.empty()
        if( !params.skip_ginkgo && (params.genome == "GRCh38" || params.genome == "GRCm39") ) {
            BAM_TO_BED_WF(
                            SENTIEON_DRIVER_DEDUP_WF.out.bam,
                            ch_bin_size,
                            ch_publish_dir,
                            ch_disable_publish
                         )
            ch_bedtools_version = BAM_TO_BED_WF.out.version

            GINKGO_BINUNSORT_WF (
                                 BAM_TO_BED_WF.out.bed_only,
                                 ch_binref.collect(),
                                 ch_bin_size,
                                 ch_publish_dir,
                                 ch_disable_publish
                               )
                               
            ch_mapped_files = GINKGO_BINUNSORT_WF.out.map{it -> it.last()}.collect()

            GINKGO_SEGMENTATION_R_WF (
                                      ch_mapped_files,
                                      ch_binref.collect(),
                                      ch_gcref.collect(),
                                      ch_boundsref_file.collect(),
                                      ch_min_ploidy,
                                      ch_max_ploidy,
                                      ch_min_bin_width,
                                      ch_bin_size,
                                      ch_is_haplotype,
                                      ch_publish_dir,
                                      ch_enable_publish
                                     )
            GINKGO_CNV_CALLER_WF (
                                    GINKGO_SEGMENTATION_R_WF.out.segcopy,
                                    ch_bin_size,
                                    ch_publish_dir,
                                    ch_enable_publish
                                 )
                                 
            GINKO_RDS_TO_FLAT(
                              GINKGO_SEGMENTATION_R_WF.out.RDS,
                              ch_bin_size,
                              ch_publish_dir,
                              ch_enable_publish
                            )
                            
            GINKO_PARSE_OUTPUTS(
                                GINKO_RDS_TO_FLAT.out.tsvs,
                                GINKGO_CNV_CALLER_WF.out.cnvs,
                                ch_binref.collect(),
                                ch_bin_size,
                                ch_publish_dir,
                                ch_enable_publish
                            )
                            
            PARSE_RDS_CNV_METRICS (
            
                            GINKGO_SEGMENTATION_R_WF.out.RDS,
                            ch_publish_dir,
                            ch_enable_publish
        
                            )
            ch_ginkgo_version = GINKGO_CNV_CALLER_WF.out.version
            ch_ginkgo_stats = PARSE_RDS_CNV_METRICS.out.collect()
        }
        
        CUSTOM_DATA_PROCESSING_WF ( 
                                  combine_outputs_b,
                                  params.reference,
                                  ch_publish_dir,
                                  ch_disable_publish
                                )

        ch_custom_calculate_mapd_metrics = Channel.empty()
        ch_custom_calculate_mapd_version = Channel.empty()
        
        if ( !params.skip_mapd ) {
            
            CUSTOM_CALCULATE_MAPD ( 
                                        SENTIEON_DRIVER_DEDUP_WF.out.bam,
                                        params.mapd_bin_size,
                                        params.blacklist_regions,
                                        params.reference,
                                        ch_publish_dir,
                                        ch_disable_publish
                                     )
            
            ch_custom_calculate_mapd_metrics = CUSTOM_CALCULATE_MAPD.out.metrics
            ch_custom_calculate_mapd_version = CUSTOM_CALCULATE_MAPD.out.version
        }
        
        
        CUSTOM_REPORT_WF ( 
                            CUSTOM_DATA_PROCESSING_WF.out.metrics
                            .collect()
                            .combine( ch_sample_metadata.collect().ifEmpty([]) )
                            .combine( ch_fastp_report.collect().ifEmpty([]) )
                            .combine( ch_bam_lorenz_coverage_stats.collect().ifEmpty([]) )
                            .combine( ch_custom_calculate_mapd_metrics.collect().ifEmpty([]) )
                            .combine( ch_ginkgo_stats.ifEmpty([]) ),
                            COUNT_READS_FASTQ_WF.out.combined_read_counts,
                            ch_min_reads,
                            ch_publish_dir,
                            ch_enable_publish
                        )
        

        ch_tool_versions = SEQTK_WF.out.version.take(1).ifEmpty([])
                            .combine(FastpNoQCWF.out.version.take(1))
                            .combine(ch_fastqc_version.take(1).ifEmpty([]))
                            .combine(SENTIEON_BWA_WF.out.version.take(1))
                            .combine(SENTIEON_DRIVER_LOCUSCOLLECTOR_WF.out.version.take(1))
                            .combine(SENTIEON_DRIVER_DEDUP_WF.out.version.take(1))
                            .combine(SENTIEON_DRIVER_METRICS_WF.out.version.take(1))
                            .combine(ch_qualimap_version.take(1).ifEmpty([]) )
                            .combine(PRESEQ_BAM2MR_WF.out.version.take(1))
                            .combine(PRESEQ_GC_EXTRAP_WF.out.version.take(1))
                            .combine(ch_kraken2_version.take(1).ifEmpty([]) )
                            .combine(ch_bam_lorenz_coverage_version.take(1).ifEmpty([]))
                            .combine((ch_ginkgo_version ?: Channel.empty()).ifEmpty([]))
                            .combine(ch_bedtools_version.take(1).ifEmpty([]))
                                                    
        REPORT_VERSIONS_WF(
                            ch_tool_versions,
                            ch_publish_dir,
                            ch_enable_publish
                          )

        collect_mqc = CUSTOM_REPORT_WF.out.mqc
                        .combine( ch_fastp_report.collect() )
                        .combine( ch_fastqc_report.collect().ifEmpty([]) )
                        .combine( ch_kraken2_report.collect().ifEmpty([]) )
                        .combine( ch_qualimap_report.collect().ifEmpty([]) )
                        .combine( REPORT_VERSIONS_WF.out.versions.collect().ifEmpty([]))

        params_meta = [
            session_id: workflow.sessionId,
            read_length: params.read_length,
            n_reads: params.n_reads,
            genome: params.genome,
            instrument: params.instrument,
            skip_fastqc: params.skip_fastqc,
            skip_kraken: params.skip_kraken,
            skip_qualimap: params.skip_qualimap,
            skip_CNV: params.skip_ginkgo
        ]

        if (!params.skip_ginkgo) {
            params_meta['bin_size'] = params.bin_size
            params_meta['min_ploidy'] = params.min_ploidy
            params_meta['max_ploidy'] = params.max_ploidy
        }
                                      
        MULTIQC_WF ( 
                  collect_mqc,
                  params_meta,
                  params.multiqc_config,
                  ch_publish_dir,
                  ch_enable_publish
                )
        
        MULTIQC_WF.out.report.ifEmpty{ exit 1, "ERROR: cannot generate any MULTIQC Report." }
  
}

workflow.onComplete {
    output                          = [:]
    output["pipeline_run_name"]     = workflow.runName
    output["pipeline_name"]         = workflow.manifest.name
    output["pipeline_version"]      = workflow.manifest.version
    output["pipeline_session_id"]   = workflow.sessionId
    output["output"]                = [:]
    output["output"]["bam"]         = [:]

    bam_outfile = file("$params.tmp_dir/bam_files.txt")
    bam_outfile_lines = bam_outfile.readLines()
    for ( bam_line : bam_outfile_lines ) {
        def (sample_name, bam_path) = bam_line.split('\t')
        output["output"]["bam"][sample_name] = [:]
        output["output"]["bam"][sample_name]["bam"] = bam_path
    }
    
    def output_json = JsonOutput.toJson(output)
    def output_json_pretty = JsonOutput.prettyPrint(output_json)
    File outputfile = new File("$params.tmp_dir/output.json")
    outputfile.write(output_json_pretty)
    println(output_json_pretty)
}

workflow.onError {
    
    output                          = [:]
    output["pipeline_run_name"]     = workflow.runName
    output["pipeline_name"]         = workflow.manifest.name
    output["pipeline_version"]      = workflow.manifest.version
    output["pipeline_session_id"]   = workflow.sessionId
    output["output"]                = [:]
    output["output"]["bam"]         = [:]
    output["output"]["vcf"]         = [:]
    
    def output_json = JsonOutput.toJson(output)
    def output_json_pretty = JsonOutput.prettyPrint(output_json)
    File outputfile = new File("$params.tmp_dir/output.json")
    outputfile.write(output_json_pretty)
    println(output_json_pretty)

    def subject = """\
        [nf-preseq-pipeline] FAILED: ${workflow.runName}
        """

    def msg = """\

        Pipeline execution summary 
        --------------------------------
        Script name       : ${workflow.scriptName ?: '-'}
        Script ID         : ${workflow.scriptId ?: '-'}
        Workflow session  : ${workflow.sessionId}
        Workflow repo     : ${workflow.repository ?: '-' }
        Workflow revision : ${workflow.repository ? "$workflow.revision ($workflow.commitId)" : '-'}
        Workflow profile  : ${workflow.profile ?: '-'}
        Workflow cmdline  : ${workflow.commandLine ?: '-'}
        Nextflow version  : ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})
        Error Report      : ${workflow.errorReport}
        """
        .stripIndent()
    
    log.info ( msg )
    
    if ( "${params.email_on_fail}" && workflow.exitStatus != 0 ) {
        sendMail(to: "${params.email_on_fail}", subject: subject, body: msg)
    }
}



