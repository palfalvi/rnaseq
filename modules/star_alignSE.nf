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
                mkdir $reads.simpleName
                STAR --runThreadN $task.cpus --genomeDir $genome_idx --readFilesIn <(gunzip -c ${reads}) --outFileNamePrefix ${reads.simpleName}/${reads.simpleName}_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
		            """
}
