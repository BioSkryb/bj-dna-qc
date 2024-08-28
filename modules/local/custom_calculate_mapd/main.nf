nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_CALCULATE_MAPD {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"


    input:
    tuple val(sample_name), path(bam), path(_)
    val(bin_size)
    path(blacklist_regions)
    path fasta_ref
    val(publish_dir)
    val(enable_publish)


    output:
    path "*_MAPD", emit: metrics
    path("custom_calculate_mapd_version.yml"), emit: version

    
    script:
    """
    cut -f 1,2 ${fasta_ref}/genome.fa.fai > chrom.sizes
    bedtools makewindows -g chrom.sizes -w ${bin_size} > intervals

    bedtools intersect -v -a intervals -b ${blacklist_regions} > intervals_no_blacklist

    bedtools coverage -sorted -g chrom.sizes -bed -b ${bam[0]} -a intervals_no_blacklist -mean > ${sample_name}_coverage_by_intervals.tsv
    
    python3 /scripts/custom_calculate_mapd.py ${sample_name}_coverage_by_intervals.tsv > ${sample_name}_MAPD
    
    echo custom_calculate_map: v0.0.1 > custom_version.yml
    export BEDTOOLS_VER=\$(bedtools --version | sed -e "s/bedtools v//g")
    echo bedtools: \$BEDTOOLS_VER > bedtools_version.yml
    cat custom_version.yml bedtools_version.yml > custom_calculate_mapd_version.yml
    """
    
}
