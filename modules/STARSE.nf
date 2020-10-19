include { star_idx } from './star_idx.nf'
include { star_alignSE } from './star_alignSE.nf'
include { collect_star } from './collect_star.nf'
include { run_fastqcSE } from './fastqcSE.nf'

workflow STARSE {
        take:
                genome
                transcriptome
		transcriptome
                read_ch
        main:
                star_idx(genome, transcriptome)
                star_alignSE(star_idx.out, read_ch)
                collect_star(star_alignSE.out, transcriptome)
                run_fastqcSE(read_ch)
        emit:
                concat(collect_star.out) | concat(run_fastqcSE.out) | collect
}
