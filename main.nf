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
    nextflow run palfalvi/rnaseq --transcriptome path/to/transcripts.fasta --read /path/to/reads/*R{1,2}.fastq.gz

    Minimal arguments:
          --transcriptome                Transcript fasta file to map to. Not required for STAR mapping.
          --reads                        Path to raw reads in fastq or fastq.gz formats. For pair end reads, use {1,2}, e.g. /path/to/reads/*R{1,2}.fastq

    Options:
          --out                          Output diresctory. [results]
          --mode                         Mapping method to use. Accepted method: 'salmon', 'kallisto' and 'star'. Default is ['salmon']
          --single                       Required for single end read processing.
          --skip_trim                    If specified, fastp runs only quality check and the mapping is done on the raw input reads.
          --save_trimmed                 Saves trimmed fastq.gz files in 'trimmed' directory.
          --skip_qc                      If specified, fastp step (including adapter trimming) is omitted standard mapping is performed directly on input reads.
          --save_index                   Save transcriptome/genome index for later use.
          --index                        External index file. Overrides index creations. If provided, transcriptome and genome options are deprecated.
          --executor                     HPC executor, if available. As this workflow is optimized to NIBB-BIAS5 server, the default is ['pbspro']

    kallisto single end specific options:
          --fragment_length              Average fragment length for the sequencing. [300]
          --fragment_sd                  Fragment size SD for the sequencing. [1]

    STAR specific options:
          --genome                       Genome fasta file.
          --gtf                          GTF file.
          --sjdbGTFfeatureExon           sjdbGTFfeatureExon option from STAR. [exon]

    Computer allocation settings
          -profile                       Sets the running environment. Default is NIBB-BIAS5 PBSPro. 'cde' and 'local' are available to run on NIBB-CDE server or on local machine.
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

include { salmon_idx } from './modules/salmon_idx.nf'
include { salmon_quant } from './modules/salmon_quant.nf'
include { salmon_quantSE } from './modules/salmon_quantSE.nf'
include { kallisto_idx } from './modules/kallisto_idx.nf'
include { kallisto_quant } from './modules/kallisto_quant.nf'
include { kallisto_quantSE } from './modules/kallisto_quantSE.nf'
include { star_idx } from './modules/star_idx.nf'
include { star_align } from './modules/star_align.nf'
include { collect_star } from './modules/collect_star.nf'
include { star_alignSE } from './modules/star_alignSE.nf'
include { collect_starSE } from './modules/collect_starSE.nf'
include { run_fastp } from './modules/fastp.nf'
include { run_fastpSE } from './modules/fastpSE.nf'
include { run_fastp_qc } from './modules/fastp_qc.nf'
include { run_fastpSE_qc } from './modules/fastpSE_qc.nf'
include { run_multiqc } from './modules/multiqc.nf'


workflow {
  log.info """
  =======================================================
             Transcriptome mapping pipeline
           https://github.com/palfalvi/rnaseq
  =======================================================

  >> Running pipeline in $params.mode mode.
  """.stripIndent()


/*
* Check if reads or SRA are provided
*/

  if ( params.reads ) {
    // Local reads provided
    if ( params.single ) {
      // Single end reads are read as Path channels
      read_ch = Channel.fromPath( params.reads )
      log.info ">> Single end reads provided."
    } else {
      // Pair end reads are read as File Pair tuples
      read_pairs_ch = Channel.fromFilePairs( params.reads )
      log.info ">> Pair end reads provided."
    }
  } else if ( params.sra ) {
    // SRA provided
    error "SRA is currently not supported. Please download fastq files and provide to --reads"
    if ( params.single ) {
      // Single end SRA is provided,
      srain = Channel.fromSRA( params.sra )
      read_ch = srain
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

      log.info ">> Index file provided: $params.index"
      } else if (params.transcriptome) {
			  transcriptome = file( params.transcriptome )
        log.info ">> Transcriptome file provided: $params.transcriptome"
        log.info ">> Building $params.mode index ..."
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

      log.info ">> Index file provided: $params.index"
      log.info ">> GTF file provided: $params.gtf"
    } else if ( params.genome && params.gtf ) {
      // Read in genome and gtf, make index
			genome = file( params.genome )
			gtf = file( params.gtf )

      log.info ">> Genome file provided: $params.genome"
      log.info ">> GTF file provided: $params.gtf"
      log.info ">> Building $params.mode index ..."
      star_idx(genome, gtf)
      idx = star_idx.out
    } else {
			error "Genome and/or GTF annotation file is not provided, but required for STAR mapping."
		}
	} else {
		error "Invalid mapping mode: ${params.mode}"
	}

/*
* Main pipeline
*/
	if( params.mode == 'salmon' && !params.single ) {
    // salmon PE mode
    if ( params.skip_qc ) {
      // Don't run fastp, just quant on input reads.
      salmon_quant(idx, read_pairs_ch)
      run_multiqc(salmon_quant.out.collect(), "$baseDir/${params.out}")
    } else if ( params.skip_trim ) {
      // Run fastp for QC only, just quant on input reads.
      run_fastp_qc(read_pairs_ch)
      salmon_quant(idx, read_pairs_ch)
      run_multiqc(salmon_quant.out.collect().concat(run_fastp_qc.out.collect()), "$baseDir/${params.out}")
    } else {
      // Run fastp and quant on trimmed reads.
      run_fastp(read_pairs_ch)
      salmon_quant(idx, run_fastp.out.trimmed)
      run_multiqc(salmon_quant.out.collect(), "$baseDir/${params.out}")
    }
	}
	else if( params.mode == 'salmon' && params.single ) {
    // salmon SE mode
    if ( params.skip_qc ) {
      // Don't run fastp, just quant on input reads.
      salmon_quantSE(idx, read_ch)
      run_multiqc(salmon_quantSE.out.collect(), "$baseDir/${params.out}")
    } else if ( params.skip_trim ) {
      // Run fastp for QC only, just quant on input reads.
      run_fastpSE_qc(read_ch)
      salmon_quantSE(idx, read_ch)
      run_multiqc(salmon_quantSE.out.collect().concat(run_fastSE_qc.out.collect()), "$baseDir/${params.out}")
    } else {
      // Run fastp and quant on trimmed reads.
      run_fastpSE(read_ch)
      salmon_quantSE(idx, run_fastpSE.out.trimmed)
      run_multiqc(salmon_quantSE.out.collect(), "$baseDir/${params.out}")
    }
		}
	else if( params.mode == 'kallisto' && !params.single ) {
      // kallisto PE mode
      if ( params.skip_qc ) {
        // Don't run fastp, just quant on input reads.
        kallisto_quant(idx, read_pairs_ch)
        run_multiqc(kallisto_quant.out.collect(), "$baseDir/${params.out}")
      } else if ( params.skip_trim ) {
        // Run fastp for QC only, just quant on input reads.
        run_fastp_qc(read_pairs_ch)
        kallisto_quant(idx, read_pairs_ch)
        run_multiqc(kallisto_quant.out.collect().concat(run_fastp_qc.out.collect()), "$baseDir/${params.out}")
      } else {
        // Run fastp and quant on trimmed reads.
        run_fastp(read_pairs_ch)
        kallisto_quant(idx, run_fastp.out.trimmed)
        run_multiqc(kallisto_quant.out.collect(), "$baseDir/${params.out}")
      }
  	}
  else if( params.mode == 'kallisto' && params.single ) {
      // kallisto SE mode
      if ( params.skip_qc ) {
        // Don't run fastp, just quant on input reads.
        kallisto_quantSE(idx, read_ch)
        run_multiqc(kallisto_quantSE.out.collect(), "$baseDir/${params.out}")
      } else if ( params.skip_trim ) {
        // Run fastp for QC only, just quant on input reads.
        run_fastpSE_qc(read_ch)
        kallisto_quantSE(idx, read_ch)
        run_multiqc(kallisto_quantSE.out.collect().concat(run_fastpSE_qc.out.collect()), "$baseDir/${params.out}")
      } else {
        // Run fastp and quant on trimmed reads.
        run_fastpSE(read_ch)
        kallisto_quantSE(idx, run_fastpSE.out.trimmed)
        run_multiqc(kallisto_quantSE.out.collect(), "$baseDir/${params.out}")
      }
      // Run multiqc after salmon_quant finished.

  }
  else if( params.mode == 'star' && !params.single ) {
    // STAR PE module
    if ( params.skip_qc ) {
      // Don't run fastp, just map on reads
      star_align(idx, read_pairs_ch)
      collect_star(star_align.out, gtf)
      run_multiqc(collect_star.out.collect(), "$baseDir/${params.out}")
    } else if ( params.skip_trim ) {
      // Run fastp for QC only, just quant on input reads.
      run_fastp_qc(read_pairs_ch)
      star_align(idx, read_pairs_ch)
      collect_star(star_align.out, gtf)
      run_multiqc(collect_star.out.collect().concat(run_fast_qc.collect()), "$baseDir/${params.out}")
    } else {
      // Run fastp and align on trimmed reads
      run_fastp(read_pairs_ch)
      star_align(idx, run_fastp.out.trimmed)
      collect_star(star_align.out, gtf)
      run_multiqc(collect_star.out.collect(), "$baseDir/${params.out}")
    }
  }
  else if( params.mode == 'star' && params.single ) {
    // STAR SE module
    if ( params.skip_qc ) {
      // Don't run fastp, just map on reads
      star_alignSE(idx, read_ch)
      collect_starSE(star_alignSE.out, gtf)
      run_multiqc(collect_starSE.out.collect(), "$baseDir/${params.out}")
    } else if ( params.skip_trim ) {
      // Run fastp for QC only, just quant on input reads.
      run_fastpSE_qc(read_ch)
      star_alignSE(idx, read_ch)
      collect_starSE(star_alignSE.out, gtf)
      run_multiqc(collect_starSE.out.collect().concat(run_fastpSE_qc.out.collect()), "$baseDir/${params.out}")
    } else {
      // Run fastp and align on trimmed reads
      run_fastpSE(read_ch)
      star_alignSE(idx, run_fastpSE.out.trimmed)
      collect_starSE(star_alignSE.out, gtf)
      run_multiqc(collect_starSE.out.collect(), "$baseDir/${params.out}")
    }
  }
	else {
		  error "Invalid mapping mode: ${params.mode}"
	}
}

workflow.onComplete {
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> Transcriptome mapping finished SUCCESSFULLY after $workflow.duration ."
      log.info "[$workflow.complete] >> Your qunatification files are in $params.out/$params.mode"
      log.info "[$workflow.complete] >> You can find further help on https://github.com/palfalvi/rnaseq"
    } else {
      log.info "[$workflow.complete] >> The script quit with ERROR after ${workflow.duration}."
      log.info "[$workflow.complete] >> Please revise your code and resubmit jobs with the -resume option or reach out for help at https://github.com/palfalvi/rnaseq."
    }
}
