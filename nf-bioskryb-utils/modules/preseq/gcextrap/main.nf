nextflow.enable.dsl=2
params.timestamp = ""

process PRESEQ_GC_EXTRAP {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/metrics/${sample_name}/preseq", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(bam)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), file("*_preseq_cov.txt"), emit: coverage
    path("preseq_qc_extrap_version.yml"), emit: version

    script:

    """
    #! /bin/bash
    set +u
    . /opt/sentieon/cloud_auth.sh no-op

    echo "0" > ${sample_name}_preseq_cov.txt        
    {
      preseq gc_extrap -o ${sample_name}.curve ${bam} && tail -n 1 ${sample_name}.curve | cut -f 2 > ${sample_name}_preseq && PRESEQ_FLOAT=\$(cat ${sample_name}_preseq) && echo \${PRESEQ_FLOAT%.*} > ${sample_name}_preseq_cov.txt
    } || {
      echo "preseq failed; likely due to an inability to accurately model library complexity; using a dummy file with 0 count" &&        touch ${sample_name}.curve
    }

    export PRESEQ_BAM2MR_VER=\$(preseq 2>&1 | grep -Eo '[0-9].[0-9].[0-9]')
    echo Preseq: \$PRESEQ_BAM2MR_VER > preseq_qc_extrap_version.yml
    """
}

workflow PRESEQ_GC_EXTRAP_WF{
  take:
    ch_mr
    ch_publish_dir
    ch_enable_publish
  main:
    PRESEQ_GC_EXTRAP ( 
                        ch_mr,
                        ch_publish_dir,
                        ch_enable_publish
                     )
  emit:
    coverage = PRESEQ_GC_EXTRAP.out.coverage
    version = PRESEQ_GC_EXTRAP.out.version
}

workflow{
  ch_mr = Channel.fromPath(params.mr)
  PRESEQ_GC_EXTRAP_WF ( 
                        ch_mr,
                        params.publish_dir,
                        params.enable_publish
                     )
}