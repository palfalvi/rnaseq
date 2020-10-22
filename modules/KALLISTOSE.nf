include { kallisto_idx } from './kallisto_idx.nf'
include { kallisto_quantSE } from './kallisto_quantSE.nf'
include { run_fastpSE } from './fastpSE.nf'

workflow KALLISTOSE {
  take:
    index
    read_ch
  main:
    run_fastp(read_ch)
    kallisto_quantSE(index, run_fastpSE.out.trimmed)
	emit:
    kallisto_quantSE.out | concat(run_fastpSE.out.json) | collect
}
