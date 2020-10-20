include { star_index } from './star_idx.nf'
include { star_alignSE } from './star_alignSE.nf'
include { collect_starSE } from './collect_starSE.nf'
include { run_fastqcSE } from './fastqcSE.nf'

workflow STARSE {
        take:
                genome
                gtf
                read_ch
        main:
                star_index(genome, gtf)
                star_alignSE(star_index.out, read_ch)
                collect_starSE(star_alignSE.out, gtf)
                run_fastqcSE(read_ch)
        emit:
                star_alignSE.out | concat(collect_starSE.out) | concat(run_fastqcSE.out) | collect
}
