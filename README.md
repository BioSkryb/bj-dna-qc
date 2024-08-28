
# BJ-DNA-QC

The BioSkryb BJ-DNA-QC pipeline evaluates the quality of the single-cell library and it provides several qc metrics to assess the quality of the sequencing reads. 

**Input Options**

The input for the pipeline is fastq files. The input can be passed either directly as path to the input file directory or via a input.csv with a meta data.

- **Reads Directory Input**: Use the `--reads` parameter to specify the path to a directory containing the input files. By default, this parameter is set to `null`. For example, to use fastq files from a specific directory, you would use: 
`--reads 's3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/*R{1,2}_001.fastq.gz'`

- **CSV Metadata Input**: Alternatively, you can use the `--input_csv` parameter to specify a CSV file containing metadata. This parameter is also `null` by default. The CSV file should have 6 columns: `biosampleName`, `sampleId`, `reads`, `readLength`, `read1` and `read2`. 
The `biosampleName` column contains the name of the biosample, the `sampleId` column contains the sample name in Illumina specified name format, `reads` have the number of reads, `readLength` with length of the reads and `read1` and `read2` has the path to the input reads. For example:

```
biosampleName,sampleId,reads,readLength,read1,read2
DNAQC-test1-100reads,DNAQC-test1-100reads_S1_L001,100,76,s3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/DNAQC-test1-100reads_S1_L001_R1_001.fastq.gz,s3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/DNAQC-test1-100reads_S1_L001_R2_001.fastq.gz
DNAQC-test2-1000reads,DNAQC-test2-1000reads_S2_L001,1000,76,s3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/DNAQC-test2-1000reads_S2_L001_R1_001.fastq.gz,s3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/DNAQC-test2-1000reads_S2_L001_R2_001.fastq.gzDNAQC-test1-100reads_S1_L001
```

**Optional Modules**


This pipeline includes several optional modules. You can choose to include or exclude these modules by adjusting the following parameters:

- `--skip_kraken`: Set this to `true` to exclude the KRAKEN2 module. By default, it is set to `true`.
- `--skip_qualimap`: Set this to `true` to exclude the Qualimap module. By default, it is set to `true`.
- `--skip_fastqc`: Set this to `true` to exclude the fastqc module. By default, it is set to `true`.
- `--skip_ginkgo`: Set this to `true` to exclude the CNV - ginkgo module. By default, it is set to `true`.
- `--skip_mapd`: Set this to `true` to exclude the MAPD module. By default, it is set to `false`.

**Outputs**


The pipeline saves its output files in the designated "publish_dir" directory. The different QC metrics files are stored in the "secondary_analyses/metrics/<sample_name>/alignment_stats" subdirectory.


**LOCAL USAGE**

Pipeline uses dsl2 and have modularize the shared code into seperate utils repository (nf-biosrkyb-utils). The utils repository contains the modules/process independently written for each tool. The utils repository also contain subworkflows, and containers used across all the pipelines.

example-

** fastq input **

```
nextflow run main.nf \
    -r master -latest -profile batch_dev \
    --reads 's3://bioskryb-test-data/dnaseq/qc/smalltest/*R{1,2}_001.fastq.gz' \
    --publish_dir 's3://bioskryb-genomics-workflows-analysis-dev/analysis/test'
```


command options

```

Usage:
    nextflow run main.nf [options]

Script Options: see nextflow.config

    
    [required]

    --reads             FILE    Path to fastq files specified as a glob pattern
    OR
    --input_csv         FILE    Path to input csv file

    --publish_dir       DIR     Path to run output directory


    [optional]

    --genomes_base      STR     Path to the genomes
                                DEFAULT: s3://bioskryb-shared-data/genomes/

    --genome            STR     Reference genome to use. Available options - GRCh38
                                DEFAULT: GRCh38

    --timestamp         STR     User can specify timestamp otherwise uses runtime generated timestamp 

    --n_reads           VAL     Number of reads to sample for analysis eg. 2.5M == 5M paired reads
                                DEFAULT: 2000000

    --read_length       VAL     Desired read length for analysis and excess to be trimmed
                                DEFAULT: 75

    --email_on_fail     STR     Email to receive upon failure
                                DEFAULT: viren.amin@bioskryb.com

    --skip_kraken       STR     Skip KRAKEN2 module
                                DEFAULT: false

    --skip_qualimap     STR     Skip Qualimap module
                                DEFAULT: true

    --skip_fastqc       STR     Skip fastqc module
                                DEFAULT: false

    --skip_ginkgo        STR     Skip CNV - ginkgo module
                                DEFAULT: false

    --instrument        STR     Specify instrument. If 'NextSeq' or 'NovaSeq' set two_color_chemistry param true
                                DEFAULT: 

    --help              BOOL    Display help message

```


**nf-test**


The BioSkryb BJ-DNA-QC nextflow pipeline run is tested using the nf-test framework.

Installation:

nf-test has the same requirements as Nextflow and can be used on POSIX compatible systems like Linux or OS X. You can install nf-test using the following command:
```
wget -qO- https://code.askimed.com/install/nf-test | bash
```
It will create the nf-test executable file in the current directory. Optionally, move the nf-test file to a directory accessible by your $PATH variable.

Usage:

```
nf-test test
```

The nf-test for this repository is saved at tests/ folder.

```
    test("preseq_test") {
        when {
            params {
                // define parameters here. Example: 
                publish_dir = "${outputDir}/results"
                timestamp = "test"
            }
        }

        then {
            assertAll(
                
                // Check if the workflow was successful
                { assert workflow.success },

                // Verify existence of the multiqc report HTML file
                { 
                    assert new File("${outputDir}/results_test/multiqc/multiqc_report.html").exists()
                },

                // Check for a match in the all metrics MQC text file
                { 
                    assert snapshot(path("${outputDir}/results_test/secondary_analyses/metrics/nf-preseq-pipeline_all_metrics_mqc.txt"))
                        .match("all_metrics_mqc")
                },

                // Check for a match in the selected metrics MQC text file
                { 
                    assert snapshot(path("${outputDir}/results_test/secondary_analyses/metrics/nf-preseq-pipeline_selected_metrics_mqc.txt"))
                        .match("selected_metrics_mqc")
                },

                // Verify existence of the fastp JSON file
                { 
                    assert new File("${outputDir}/results_test/primary_analyses/metrics/H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64/fastp/H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64_no_qc_fastp.json").exists()
                },

                // Check for a match in the kraken report text file
                { 
                    assert snapshot(path("${outputDir}/results_test/primary_analyses/metrics/H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64/kraken2/H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64_kraken2_report.txt"))
                        .match("kraken_report")
                },

                // Verify existence of the qualimap HTML report
                { 
                    assert new File("${outputDir}/results_test/primary_analyses/metrics/H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64/qualimap/QUALIMAP_BAMQC_WF_H5L3TM-DNA-NA12878-23413-02F-100PG-V1-64/qualimapReport.html").exists()
                },

                // Check for a match in the Ginkgo segment summary text file
                { 
                    assert snapshot(path("${outputDir}/results_test/tertiary_analyses/cnv_ginkgo/AllSample-GinkgoSegmentSummary.txt"))
                        .match("ginkgo_summary")
                }
            )
            
        }

    }
```
