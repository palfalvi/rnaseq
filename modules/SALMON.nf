include { salmon_idx } from './salmon_idx.nf'
include { salmon_quant } from './salmon_quant.nf'
include { run_fastqc } from './fastqc.nf'

workflow SALMON {
        take:
                transcriptome 
                read_pairs_ch
                index
        main:
		         
		          salmon_idx(transcriptome)
		          index = salmon_idx.out
		          salmon_quant(index, read_pairs_ch)
              run_fastqcSE(read_pairs_ch)

	      emit:
		          salmon_quant.out | concat(run_fastqc.out) | collect
}
