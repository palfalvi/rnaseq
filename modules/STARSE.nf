include { star_index } from './star_idx.nf'
include { star_alignSE } from './star_alignSE.nf'
include { collect_star } from './collect_star.nf'
include { run_fastqcSE } from './fastqcSE.nf'

workflow STARSE {
        take:
                genome
                gtf
                read_ch
        main:
                star_index(genome, gtf)
                star_alignSE(star_index.out, read_ch)
                collect_star(star_alignSE.out, gtf)
                run_fastqcSE(read_ch)
        emit:
                star_alignSE.out | concat(collect_star.out) | concat(run_fastqcSE.out) | collect
}
