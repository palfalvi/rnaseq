include { kallisto_idx } from './kallisto_idx.nf'
include { kallisto_quant } from './kallisto_quant.nf'
include { run_fastp } from './fastp.nf'

workflow KALLISTO {
  take:
    index
    read_pairs_ch
  main:
    run_fastp(read_pairs_ch)
    kallisto_quant(index, run_fastp.out.trimmed)
	emit:
    kallisto_quant.out | concat(run_fastp.out.json) | collect
}
