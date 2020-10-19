include { star_idx } from './star_idx.nf'
include { star_align } from './star_align.nf'
include { collect_star } from './collect_star.nf'
include { run_fastqc } from './fastqc.nf'

workflow STAR {
        take:
                genome
                transcriptome
		transcriptome
                read_pairs_ch
        main:
                star_idx(genome, transcriptome)
                star_align(star_idx.out, read_pairs_ch)
                collect_star(star_alignSE.out, transcriptome)
                run_fastqc(read_pairs_ch)
        emit:
                collect_star.out | concat(run_fastqc.out) | collect
}
