process star_align {
	tag "$sample_id"

	publishDir "${params.out}/star", mode: 'move'

	cpus params.mapping.cpus
        conda 'star=2.7.6a'

        input:
                path genome_idx
                tuple val(sample_id), file(reads)
        output:
                path sample_id
        script:
                """
                STAR --runThreadN $task.cpus --genomeDir $genome_idx --readFilesIn ${read[0]} ${read[1]} --outFileNamePrefix $sample_id --quantMode GeneCounts
		"""
}
