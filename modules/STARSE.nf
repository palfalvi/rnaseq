include { star_alignSE } from './star_alignSE.nf'
include { collect_starSE } from './collect_starSE.nf'
include { run_fastpSE } from './fastpSE.nf'

workflow STARSE {
  take:
    read_ch
    index
    gtf
  main:
    run_fastpSE(read_ch)
    star_alignSE(index, run_fastpSE.out.trimmed)
    collect_starSE(star_alignSE.out, gtf)
  emit:
    star_alignSE.out | concat(collect_starSE.out) | concat(run_fastpSE.out.json) | collect
}
