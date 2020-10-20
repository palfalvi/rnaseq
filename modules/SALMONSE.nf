include { salmon_idx } from './salmon_idx.nf'
include { salmon_quantSE } from './salmon_quantSE.nf'
include { run_fastqcSE } from './fastqcSE.nf'

workflow SALMONSE {
        take:
                transcriptome 
                read_ch
                salmon_idx
        main:
		          if ( $salmon_idx ) {
		            salmon_idx(transcriptome)
                salmon_quantSE(salmon_idx.out, read_ch)
		          } else {
		            salmon_quantSE(salmon_idx, read_ch)
		          }
              
              run_fastqcSE(read_ch)
                
	      emit:
	          	salmon_quantSE.out | concat(run_fastqcSE.out) | collect
}
