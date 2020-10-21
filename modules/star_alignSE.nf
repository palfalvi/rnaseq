process star_alignSE {
	tag "$reads.simpleName"

	publishDir "${params.out}/star", mode: 'copy'

	cpus "$params.mapping.cpus"
        conda './conda-envs/star-env.yaml'

        input:
                path genome_idx
                path reads
        output:
                path "$reads.simpleName"
        script:
                """
                mkdir "${reads.simpleName}"
                STAR --runThreadN $task.cpus --genomeDir $genome_idx --readFilesIn $reads --readFilesCommand zcat --outFileNamePrefix ${reads.simpleName}/${reads.simpleName}_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
		            """
}
