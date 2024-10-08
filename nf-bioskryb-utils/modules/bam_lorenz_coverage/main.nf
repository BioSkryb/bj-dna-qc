nextflow.enable.dsl=2
params.timestamp = ""

process BAM_LORENZ_COVERAGE {
  tag "${sample_name}"
  publishDir "${params.publish_dir}_${params.timestamp}/secondary_analyses/metrics/${sample_name}/bam_lorenz_coverage", enabled:"$enable_publish"

  input:
  tuple val(sample_name), path(bam), path(bai)
  val(publish_dir)
  val(enable_publish)

  output:
  path("*lorenz*"), emit: results
  path("*lorenzstats.tsv"), emit: stats
  path("bam_lorenz_coverage.yml"), emit: version

  script:
  """
  
  bam-lorenz-coverage -l ${sample_name}_lorenzcurve.tsv  -c ${sample_name}_lorenzcoverage.tsv -L ${sample_name}_lorenzcurve.svg -C ${sample_name}_lorenzcoverage.svg -s ${sample_name}_lorenzstats.tsv ${bam}

  export  LORENZ_VER=\$(bam-lorenz-coverage --version | awk '/version/{print \$3}')
  echo bam-lorenz-coverage: \$LORENZ_VER > bam_lorenz_coverage.yml
  
  """
}


workflow BAM_LORENZ_COVERAGE_WF{
    take:
        ch_bam
        ch_publish_dir
        ch_enable_publish
    main:
        BAM_LORENZ_COVERAGE ( ch_bam,
                         ch_publish_dir,
                         ch_enable_publish
                       )
    emit:
        results = BAM_LORENZ_COVERAGE.out.results
        stats = BAM_LORENZ_COVERAGE.out.stats
        version = BAM_LORENZ_COVERAGE.out.version
}

workflow{

    
     
    ch_bam_raw = Channel.fromFilePairs(params.bam, size: -1)
    ch_bam_raw
              .map{ it -> it.flatten().collect() }
              .set{ ch_bam }
              
              
    BAM_LORENZ_COVERAGE_WF(
                        ch_bam,
                        params.publish_dir,
                        params.enable_publish
                     )
}




