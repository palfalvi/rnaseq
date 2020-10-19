process star_alignSE {
	tag "$reads.simpleName"

	publishDir "${params.out}/star", mode: 'move'

	cpus params.mapping.cpus
        conda 'star=2.7.6a'

        input:
                path genome_idx
                path reads
        output:
                path "${read.simpleName}"
        script:
                """
                STAR --runThreadN $task.cpus --genomeDir $genome_idx --readFilesIn ${reads} --outFileNamePrefix ${reads.simpleName} --quantMode GeneCounts
		"""
}
