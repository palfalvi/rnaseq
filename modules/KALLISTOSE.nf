include { kallisto_idx } from './kallisto_idx.nf'
include { kallisto_quantSE } from './kallisto_quantSE.nf'
include { run_fastqcSE } from './fastqcSE.nf'

workflow KALLISTOSE {
        take:
                transcriptome
                read_ch
        main:
                kallisto_idx(transcriptome)
                kallisto_quantSE(kallisto_idx.out, read_ch)
                run_fastqcSE(read_ch)
	emit:
                kallisto_quantSE.out | concat(run_fastqcSE.out) | collect
}
