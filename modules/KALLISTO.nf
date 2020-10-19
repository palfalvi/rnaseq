include { kallisto_idx } from './kallisto_idx.nf'
include { kallisto_quant } from './kallisto_quant.nf'
include { run_fastqc } from './fastqc.nf'

workflow KALLISTO {
        take:
                transcriptome
                read_pairs_ch
        main:
                kallisto_idx(transcriptome)
                kallisto_quant(kallisto_idx.out, read_pairs_ch)
                run_fastqc(read_pairs_ch)
	emit:
                kallisto_quant.out | concat(run_fastqc.out) | collect
}
