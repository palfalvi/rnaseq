process star_alignSE {
	tag "$reads.simpleName"

	publishDir "${params.out}/star", mode: 'copy'

	cpus "$params.mapping.cpus"
        conda 'star=2.7.6a'

        input:
                path genome_idx
                path reads
        output:
                path "$reads.simpleName"
        script:
                """
                echo "${reads.simpleName}"
		            """
}
