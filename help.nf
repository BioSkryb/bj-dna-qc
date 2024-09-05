nextflow.enable.dsl=2

def printHeader() {
  
  log.info """\
  BJ-DNA-QC   P I P E L I N E
  ===================================
  glob fastq files    : ${ params.reads }
  csv fastq files     : ${ params.input_csv }
  publish_dir         : ${ params.publish_dir }
  timestamp           : ${ params.timestamp }
  genome              : ${ params.genome }
  n_reads             : ${ params.n_reads }
  skip_ginkgo         : ${ params.skip_ginkgo }
  \n
  """

}

def helpMessage() {

  yellow = "\033[0;33m"
  blue = "\033[0;34m"
  white = "\033[0m"
  red = "\033[0;31m"

  log.info """\
${blue}
    bj-dna-qc pipeline

    Usage:
        nextflow run main.nf [options]

    Script Options: see nextflow.config

${red}
        [required]
        --input_csv         FILE    Path to input csv file

${yellow}
        [optional]
        
        --genomes_base      STR     Path to the genomes
                                    DEFAULT: ${params.genomes_base}

        --genome            STR     Reference genome to use. Available options - GRCh38
                                    DEFAULT: ${params.genome}

        --publish_dir       DIR     Path to run output directory
                                    DEFAULT: ${params.publish_dir}
                                    
        --timestamp         STR     User can specify timestamp otherwise uses runtime generated timestamp 

        --n_reads           VAL     Number of reads to sample for analysis eg. 2.5M == 5M paired reads
                                    DEFAULT: ${params.n_reads}

        --read_length       VAL     Desired read length for analysis and excess to be trimmed
                                    DEFAULT: ${params.read_length}

        --email_on_fail     STR     Email to receive upon failure
                                    DEFAULT: ${params.email_on_fail}
                                    
        --skip_kraken       STR     Skip KRAKEN2 module
                                    DEFAULT: ${params.skip_kraken}
                                    
        --skip_qualimap     STR     Skip Qualimap module
                                    DEFAULT: ${params.skip_qualimap}
                                    
        --skip_mapd         STR     Skip MAPD module. MAPD is a measurement of the bin-to-bin variation in read coverage that is robust to the presence of CNVs, and is an indicator of the evenness of whole genome amplification (WGA)
                                    DEFAULT: ${params.skip_mapd}
        
        --skip_fastqc       STR     Skip fastqc module
                                    DEFAULT: ${params.skip_fastqc}
                                    
        --skip_ginkgo       STR     Skip CNV - ginkgo module
                                    DEFAULT: ${params.skip_ginkgo}
                                    
        --help              BOOL    Display help message
                                    
${white}
    """.stripIndent()
}

workflow{
  printHeader()
  helpMessage()
}
