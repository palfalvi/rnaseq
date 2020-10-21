include { star_alignSE } from './star_alignSE.nf'
include { collect_starSE } from './collect_starSE.nf'
include { run_fastqcSE } from './fastqcSE.nf'

workflow STARSE {
        take:
                read_ch
                index
                gtf
        main:
                star_alignSE(index, read_ch)
                collect_starSE(star_alignSE.out, gtf)
                run_fastqcSE(read_ch)
        emit:
                star_alignSE.out | concat(collect_starSE.out) | concat(run_fastqcSE.out) | collect
}
