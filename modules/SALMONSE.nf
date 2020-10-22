include { salmon_idx } from './salmon_idx.nf'
include { salmon_quantSE } from './salmon_quantSE.nf'
include { run_fastqcSE } from './fastqcSE.nf'
include { run_fastpSE } from './fastpSE.nf'

workflow SALMONSE {
  take:
    index
    read_ch
  main:
    run_fastp(read_ch)
    salmon_quantSE(index, run_fastpSE.out.trimmed)
	emit:
	  salmon_quantSE.out | concat(run_fastpSE.out.json) | collect
}
