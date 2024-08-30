
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

The input for the pipeline is fastq files. The input can be passed either directly as path to the input file directory or via a input.csv with a meta data.

- **Reads Directory Input**: Use the `--reads` parameter to specify the path to a directory containing the input files. By default, this parameter is set to `null`. For example, to use fastq files from a specific directory, you would use: 
`--reads 's3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/*R{1,2}_001.fastq.gz'`

- **CSV Metadata Input**: Alternatively, you can use the `--input_csv` parameter to specify a CSV file containing metadata. This parameter is also `null` by default. The CSV file should have 6 columns: `biosampleName`, `sampleId`, `reads`, `readLength`, `read1` and `read2`. 
The `biosampleName` column contains the name of the biosample, the `sampleId` column contains the sample name in Illumina specified name format, `reads` have the number of reads, `readLength` with length of the reads and `read1` and `read2` has the path to the input reads. For example:

```
biosampleName,sampleId,reads,readLength,read1,read2
DNAQC-test1-100reads,DNAQC-test1-100reads_S1_L001,100,76,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test1-100reads_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test1-100reads_S1_L001_R2_001.fastq.gz
DNAQC-test2-1000reads,DNAQC-test2-1000reads_S2_L001,1000,76,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test2-1000reads_S2_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/DNAQC-test2-1000reads_S2_L001_R2_001.fastq.gzDNAQC-test1-100reads_S1_L001
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
nextflow run main.nf -c conf/test_csv_input.config
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
