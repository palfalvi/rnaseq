include { salmon_idx } from './salmon_idx.nf'
include { salmon_quant } from './salmon_quant.nf'
include { run_fastqc } from './fastqc.nf'
include { run_fastp } from './fastp.nf'

workflow SALMON {
        take:
                index
                read_pairs_ch
        main:
              run_fastp(read_pairs_ch)
              salmon_quant(index, run_fastp.out.trimmed)
	      emit:
		          salmon_quant.out | concat(run_fastp.out.json) | collect
}
