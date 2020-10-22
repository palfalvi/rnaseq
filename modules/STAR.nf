include { star_align } from './star_align.nf'
include { collect_star } from './collect_star.nf'
include { run_fastqc } from './fastqc.nf'

workflow STAR {
  take:
    read_pairs_ch
    index
    gtf
  main:
    run_fastp(read_pairs_ch)
    star_align(index, run_fastp.out.trimmed)
    collect_star(star_align.out, gtf)
  emit:
    star_align.out | concat(run_fastp.out.json) | collect
}
