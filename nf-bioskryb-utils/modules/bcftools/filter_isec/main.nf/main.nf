nextflow.enable.dsl=2
params.timestamp = ""

process BCFTOOLS_ISEC {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/mutationalcatalogs", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(input_vcf), path(input_vcf_index)
    path baseline_vcf
    path vcfeval_baseline_vcf_index
    path reference
    path dbsnp
    path dbsnp_index
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("${sample_name}.output.filteredv2.vcf.gz*"), emit: filtered_vcf
    path("${sample_name}.output.filteredv2.stats.txt"), emit: filtered_stats
    path("bcftools_version.yml"), emit: version

    script:
    """
    bcftools view --threads ${task.cpus} -i 'GT!="0/0"' ${input_vcf} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa | bcftools view --threads ${task.cpus} -i 'GT!="0/0"' | bcftools +fill-tags | bcftools view -Oz -o ${sample_name}_nonref.vcf.gz

    tabix -p vcf ${sample_name}_nonref.vcf.gz

    bcftools isec -C -O z -w1 ${sample_name}_nonref.vcf.gz ${baseline_vcf} > ${sample_name}_nonref_germline_filtered.vcf.gz

    tabix -p vcf ${sample_name}_nonref_germline_filtered.vcf.gz

    bcftools isec -C -O z -w1 ${sample_name}_nonref_germline_filtered.vcf.gz ${dbsnp} > ${sample_name}_nonref_germline_filtered_dbsnp_filtered.vcf.gz

    tabix -p vcf ${sample_name}_nonref_germline_filtered_dbsnp_filtered.vcf.gz

    bcftools filter -e 'INFO/DP<5 || QUAL<50' -Oz -o ${sample_name}.output.filteredv2.vcf.gz ${sample_name}_nonref_germline_filtered_dbsnp_filtered.vcf.gz

    tabix -p vcf ${sample_name}.output.filteredv2.vcf.gz

    bcftools stats ${sample_name}.output.filteredv2.vcf.gz > ${sample_name}.output.filteredv2.stats.txt

    # Check if there are any variants in the final VCF
    variant_count=\$(grep "number of records:" ${sample_name}.output.filteredv2.stats.txt | awk '{print \$NF}')

    if [ "\$variant_count" -eq 0 ]; then
        echo "Error: No variants found in the filtered VCF file"
        exit 1
    fi

    export BCFTOOLS_VER=\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    echo "bcftools: \$BCFTOOLS_VER" > bcftools_version.yml
    """
}