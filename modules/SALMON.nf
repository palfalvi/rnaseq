include { salmon_idx } from './salmon_idx.nf'
include { salmon_quant } from './salmon_quant.nf'
include { run_fastqc } from './fastqc.nf'

workflow SALMON {
        take:
                transcriptome 
                read_pairs_ch
                salmon_idx
        main:
		          if ( $salmon_idx ) {
		            salmon_idx(transcriptome)
                salmon_quant(salmon_idx.out, read_pairs_ch)
		          } else {
		            salmon_quant(salmon_idx, read_pairs_ch)
		          }
              
              run_fastqcSE(read_pairs_ch)

	      emit:
		          salmon_quant.out | concat(run_fastqc.out) | collect
}
