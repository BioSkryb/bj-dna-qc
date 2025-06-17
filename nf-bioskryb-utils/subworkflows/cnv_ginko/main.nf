nextflow.enable.dsl=2
 
include { BAM_TO_BED} from '../../modules/ginkgo/bam_to_bed/main.nf'
include { GINKGO_BINUNSORT } from '../../modules/ginkgo/binunsort/main.nf'
include { GINKGO_SEGMENTATION_R } from '../../modules/ginkgo/segmentation_r/main.nf'
include { GINKGO_CNV_CALLER } from '../../modules/ginkgo/cnvcaller/main.nf'
include { GINKO_RDS_TO_FLAT } from '../../modules/ginkgo/rds_to_flat/main.nf'
include { PARSE_RDS_CNV_METRICS } from '../../modules/ginkgo/parse_rds_cnv_metrics/main.nf'
include { GINKO_PARSE_OUTPUTS } from '../../modules/ginkgo/parse_ginko_outputs/main.nf'
 
workflow GINKO_WF{
     take:
        ch_bam
        ch_bin_size
        ch_binref
        ch_gcref
        ch_boundsref_file
        ch_min_ploidy
        ch_max_ploidy
        ch_min_bin_width
        ch_is_haplotype
        ch_publish_dir
        ch_enable_publish
        ch_disable_publish
        
     main:

        BAM_TO_BED (
                        ch_bam,
                        ch_bin_size,
                        ch_publish_dir,
                        ch_disable_publish
                     )

        GINKGO_BINUNSORT (
                             BAM_TO_BED.out.bed_only,
                             ch_binref.collect(),
                             ch_bin_size,
                             ch_publish_dir,
                             ch_disable_publish
                           )
                           
        ch_mapped_files = GINKGO_BINUNSORT.out.map{it -> it.last()}.collect()

        GINKGO_SEGMENTATION_R (
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
        GINKGO_CNV_CALLER (
                                GINKGO_SEGMENTATION_R.out.segcopy,
                                ch_bin_size,
                                ch_publish_dir,
                                ch_disable_publish
                            )
                             
        GINKO_RDS_TO_FLAT(
                          GINKGO_SEGMENTATION_R.out.RDS,
                          ch_bin_size,
                          ch_publish_dir,
                          ch_disable_publish
                        )
                        
        GINKO_PARSE_OUTPUTS(
                            GINKO_RDS_TO_FLAT.out.tsvs,
                            GINKGO_CNV_CALLER.out.cnvs,
                            ch_binref.collect(),
                            ch_bin_size,
                            ch_publish_dir,
                            ch_enable_publish
                           )
                           
        PARSE_RDS_CNV_METRICS ( 
                            GINKGO_SEGMENTATION_R.out.RDS,
                            ch_publish_dir,
                            ch_enable_publish)
                           
    emit:
        METRICS = PARSE_RDS_CNV_METRICS.out
        CNVS = GINKO_PARSE_OUTPUTS.out.tsvs
        graph = GINKGO_SEGMENTATION_R.out.jpeg
        ginkgo_version = GINKGO_CNV_CALLER.out.version
        ginkgo_stats = PARSE_RDS_CNV_METRICS.out.collect()
        bedtools_version = BAM_TO_BED.out.version

        
 }
  
 workflow{
    ch_reference = Channel.fromPath(params.reference)
    ch_bam = Channel.fromFilePairs(params.bam_dir + "/*" + params.sample_name + "*.{bam,bam.bai}")
    ch_binref = Channel.fromPath(params.ginko_ref_dir + "variable_" + params.bin_size + "_76_bwa")
    ch_gcref = Channel.fromPath(params.ginko_ref_dir + "GC_variable_" + params.bin_size + "_76_bwa")
    ch_boundsref_file = Channel.fromPath(params.ginko_ref_dir + "bounds_variable_" + params.bin_size + "_76_bwa")
    
    GINKO_WF(
        ch_bam,
        params.bin_size,
        ch_binref,
        ch_gcref,
        ch_boundsref_file,
        params.min_ploidy,
        params.max_ploidy,
        params.min_bin_width,
        params.is_haplotype,
        params.publish_dir,
        params.enable_publish,
        params.disable_publish
    )
    
 }