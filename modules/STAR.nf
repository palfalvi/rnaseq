include { star_index } from './star_idx.nf'
include { star_align } from './star_align.nf'
include { collect_star } from './collect_star.nf'
include { run_fastqc } from './fastqc.nf'

workflow STAR {
        take:
                read_pairs_ch
                index
                gtf
        main:
                star_align(index, read_pairs_ch)
                collect_star(star_align.out, gtf)
                run_fastqc(read_pairs_ch)
        emit:
                star_align.out | concat(run_fastqc.out) | collect
}
