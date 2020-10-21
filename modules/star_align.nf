process star_align {
	tag "$sample_id"

	publishDir "${params.out}/star", mode: 'copy'

	cpus "$params.mapping.cpus"
        conda 'star=2.7.6a'

        input:
                path genome_idx
                tuple val(sample_id), file(reads)
        output:
                path sample_id
        script:
                """
                #mkdir $sample_id
                #STAR --runThreadN $task.cpus --genomeDir $genome_idx --readFilesIn $reads --readFilesCommand zcat --outFileNamePrefix ${sample_id}/${sample_id}_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
		            """
}
