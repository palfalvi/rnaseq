include { salmon_idx } from './salmon_idx.nf'
include { salmon_quantSE } from './salmon_quantSE.nf'
include { run_fastqcSE } from './fastqcSE.nf'

workflow SALMONSE {
        take:
                index 
                read_ch
        main:
                salmon_quantSE(index, read_ch)
              run_fastqcSE(read_ch)

	      emit:
	          	salmon_quantSE.out | concat(run_fastqcSE.out) | collect
}
