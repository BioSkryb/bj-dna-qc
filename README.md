
# BJ-DNA-QC

The BioSkryb BJ-DNA-QC pipeline evaluates the quality of the single-cell library and it provides several qc metrics to assess the quality of the sequencing reads.

# Running locally

Following are instructions for running BJ-DNA-QC in a local Ubuntu 18.04 server

## Install Java 11

```
sudo apt-get install default-jdk

java -version
#openjdk version "11.0.18" 2023-01-17
#OpenJDK Runtime Environment (build 11.0.18+10-post-Ubuntu-0ubuntu118.04.1)
#OpenJDK 64-Bit Server VM (build 11.0.18+10-post-Ubuntu-0ubuntu118.04.1, mixed mode, sharing)
```

## Install AWS CLI

```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

## Install Nextflow

```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

```

## Install Docker

```
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

```

## Download and Unzip Pipeline Repository
```
cd bj-dna-qc
```

## Sentieon License Setup

The Sentieon license is based on a lightweight floating license server process running on one node, and serving licenses though TCP to all other nodes. Normally this license server is running in a special non-computing node on the cluster periphery that has unrestricted access to the outside world through HTTPS, and serves the licenses to the rest of the nodes in the cluster by listening to a specific TCP port that needs to be open within the cluster. The license server needs to have external https access to validate the license, while the computing nodes do not need to have access to the internet.

Client will need to provide FQDN (hostname) or IP address of the machine that they would like to use to host the license server along with the port the license server will listen at to create a license file by Sentieon.

Sentieon license server supports connection through the proxy server. Set the standard `http_proxy` environment before starting the license server.

In order to run Sentieon software will need to start a Sentieon license server on eg. `b06x-pbs01.inet.xxxxxxxxx`; running the following command will setup the license server as a running daemon in your system:

```
export http_proxy=<proxy_server_name_and_port>
<SENTIEON_DIR>/bin/sentieon licsrvr --start --log <LOCATION_OF_LOG_FILE> <LICENSE_FILE> 
```
To run Sentieon software in the computing nodes, will need to set an environmental variable to tell Sentieon software the location of the license server and port. This can be added to bash profile or to the scripts that will drive the pipelines:
```
export SENTIEON_LICENSE=b06x-pbs01.inet.xxxxxxxxxxxx:xxxx
```

## Test Pipeline Execution

**Input Options**

The input for the pipeline can be passed via a input.csv with a meta data.

- **CSV Metadata Input**: The CSV file should have 4 columns: `biosampleName`, `reads`, `read1` and `read2`. 
The `biosampleName` column contains the name of the biosample, `reads` have the number of reads and `read1` and `read2` has the path to the input reads. For example:

```
biosampleName,reads,read1,read2
DNAQC-test1-100reads,100,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test1-100reads_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test1-100reads_S1_L001_R2_001.fastq.gz
DNAQC-test2-1000reads,1000,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test2-1000reads_S2_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test2-1000reads_S2_L001_R2_001.fastq.gzDNAQC-test1-100reads_S1_L001
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


**Command**

example-

** csv input **

```
nextflow run main.nf --input_csv $PWD/tests/data/inputs/input.csv
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
**Tool versions**

- `Seqtk: 1.3-r106`
- `fastp: 0.20.1`
- `FastQC: v0.11.9`
- `Sentieon: 202308.01`
- `QualiMap: v.2.2.2-dev`
- `Preseq: 2.0.3`
- `Kraken2: 2.1.3`
- `bam-lorenz-coverage: 2.3.0 GNU`
- `Ginkgo: 0.0.2`
- `bedtools: v2.28.0`

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