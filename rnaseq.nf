#!/usr/bin/env nextflow

/*
*============================================
* Transcriptome Analysis Pipeline.
* https://github.com/palfalvi/rnaseq
*--------------------------------------------
*/


def helpMessage() {
    log.info """
    =======================================================
               Transcriptome mapping pipeline
             https://github.com/palfalvi/rnaseq
    =======================================================

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run palfalvi/rnaseq --transcriptome path/to/transcripts.fasta --read /path/to/reads/*R{1,2}.fastq.gz --mode 'salmon'

    Minimal arguments:
          --transcriptome                Transcript fasta file to map to. Not required for STAR mapping.
          --reads                        Path to raw reads in fastq or fastq.gz formats. For pair end reads, use {1,2}, e.g. /path/to/reads/*R{1,2}.fastq

    Options:
          --out                          Output diresctory. [results]
          --mode                         Mapping method to use. Accepted method: 'salmon', 'kallisto' and 'star'. Default is ['salmon']
          --single                       Required for single end read processing. [false]
          --skip_qc                      If specified, fastqc step is skipped and only mapping is performed. [false]
          --save_index                   Save index folder for later use. [false]
          --fastqc.cpus                  Number of threads to use for fastqc. [2]
          --mapping.cpus                 Number of threads to use for mapping. [20]
          --executor                     HPC executor, if available. As this workflow is optimized to NIBB-BIAS5 server, the default is ['pbspro']
          --save_index                   Saves index file for future reuse.
          --index                        External index file. Overrides index creations. If provided, transcriptome and genome options are deprecated.

    kallisto single end specific options:
          --fragment_length              Average fragment length for the sequencing. [300]
          --fragment_sd                  Fragment size SD for the sequencing. [1]

    STAR specific options:
          --genome                       Genome fasta file.
          --gtf                          GTF file.
          --sjdbGTFfeatureExon           sjdbGTFfeatureExon option from STAR. [exon]
    """.stripIndent()
}

params.help = false
if (params.help){
    helpMessage()
    exit 0
}


nextflow.enable.dsl=2

/*
* Include dsl2 modules
*/
include { SALMON } from './modules/SALMON.nf'
include { SALMONSE } from './modules/SALMONSE.nf'
include { KALLISTO } from './modules/KALLISTO.nf'
include { KALLISTOSE } from './modules/KALLISTOSE.nf'
include { STAR } from './modules/STAR.nf'
include { STARSE } from './modules/STARSE.nf'
include { star_idx } from './modules/star_idx.nf'
include { salmon_idx } from './modules/salmon_idx.nf'
include { kallisto_idx } from './modules/kallisto_idx.nf'
include { run_multiqc } from './modules/multiqc.nf'


workflow {

/*
* Check if reads or SRA are provided
*/

  if ( params.reads ) {
    // Local reads provided
    if ( params.single ) {
      // Single end reads are read as Path channels
      read_ch = Channel.fromPath( params.reads )
    } else {
      // Pair end reads are read as File Pair tuples
      read_pairs_ch = Channel.fromFilePairs( params.reads )
    }
  } else if ( params.sra ) {
    // SRA provided
    // Can we automatically detect if SRA is SE or PE?
    if ( params.single ) {
      // Single end SRA is provided,
      srain = Channel.fromSRA( params.sra )
      read_ch = srain[1]
    } else {
      read_pairs_ch = Channel.fromSRA( params.sra )
    }
  } else {
    error "No reads provided. Please specify reads with the --reads flag."
  }


/*
* Check if transcriptome/genome/gtf files are provided and read them in
* Also look for external index files
*/
	if (params.mode == 'salmon' || params.mode == 'kallisto') {
		if ( params.index ) {
      // Read in salmon or kallisto index
      idx = file( params.index )
      } else if (params.transcriptome) {
			  transcriptome = file( params.transcriptome )
        if (params.mode == 'salmon') {
          // Make salmon index
          salmon_idx(transcriptome)
          idx = salmon_idx.out
        } else if (params.mode == 'kallisto') {
          // Make kallisto index
          kallisto_idx(transcriptome)
          idx = kallisto_idx.out
          }
		} else {
	 		error "Transcriptome fasta file is not provided for ${params.mode} mapping. Please specify --transcriptome flag"
		}
	} else if (params.mode == 'star') {
		if ( params.index && params.gtf ) {
      // Read in index and gtf files
      idx = file( params.index )
      gtf = file( params.gtf )
    } else if ( params.genome && params.gtf ) {
      // Read in genome and gtf, make index
			genome = file( params.genome )
			gtf = file( params.gtf )
      star_idx(genome, gtf)
      idx = star_idx.out
    } else {
			error "Genome and/or GTF annotation file is not provided, but required for STAR mapping."
		}
	} else {
		error "Invalid mapping mode: ${params.mode}"
	}

/*
* Optional trimming
*/
  if ( params.trim ) {
    if ( params.single ) {
      // single end trim
    } else {
      // pair end trim
    }
  }

/*
* Main pipeline
*/
	if( params.mode == 'salmon' && !params.single ) {
	  	SALMON(idx, read_pairs_ch)
		  run_multiqc(SALMON.out, "$baseDir/${params.out}")
		}
	else if( params.mode == 'salmon' && params.single ) {
      SALMONSE(idx, read_ch)
      run_multiqc(SALMONSE.out, "$baseDir/${params.out}")
		}
	else if( params.mode == 'kallisto' && !params.single ) {
	  	KALLISTO(idx, read_pairs_ch)
	  	run_multiqc(KALLISTO.out, "$baseDir/${params.out}")
		}
	else if( params.mode == 'kallisto' && params.single ) {
      KALLISTOSE(idx, read_ch)
      run_multiqc(KALLISTOSE.out, "$baseDir/${params.out}")
		}
  else if( params.mode == 'star' && !params.single ) {
      STAR(read_pairs_ch, idx, gtf)
      run_multiqc(STAR.out, "$baseDir/${params.out}")
    }
  else if( params.mode == 'star' && params.single ) {
      STARSE(read_ch, idx, gtf)
      run_multiqc(STARSE.out, "$baseDir/${params.out}")
                }
	else {
		  error "Invalid mapping mode: ${params.mode}"
		}
}

workflow.onComplete {
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> Transcriptome mapping finished SUCCESSFULLY after $workflow.duration . Have fun with the data!"
    } else {
      log.info "[$workflow.complete] >> The script quit with ERROR after ${workflow.duration}. Please revise your code and resubmit jobs with the -resume option."
    }
}
