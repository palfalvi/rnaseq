include { star_idx } from './star_idx.nf'
include { star_align } from './star_align.nf'
include { collect_star } from './collect_star.nf'
include { run_fastqc } from './fastqc.nf'

workflow STAR {
        take:
                genome
                gtf
                read_pairs_ch
        main:
                star_idx(genome, gtf)
                star_align(star_idx.out, read_pairs_ch)
                collect_star(star_align.out, gtf)
                run_fastqc(read_pairs_ch)
        emit:
                collect_star.out | concat(run_fastqc.out) | collect
}
