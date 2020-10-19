process star_alignSE {
	tag "$reads.simpleName"

	publishDir "${params.out}/star", mode: 'move'

	cpus params.mapping.cpus
        conda 'star=2.7.6a'

        input:
                path genome_idx
                path reads
        output:
                path "$reads.simpleName"

                if (reads.getExtension == 'fastq' || reads.getExtension == 'fq') {
                script:
                  """
                  STAR --runThreadN $task.cpus --genomeDir $genome_idx --readFilesIn ${reads} --outFileNamePrefix ${reads.simpleName} --quantMode GeneCounts
                  """
                } else {
                  script:
                  """
                  STAR --runThreadN $task.cpus --genomeDir $genome_idx --readFilesIn <(gunzip -c ${reads}) --outFileNamePrefix ${reads.simpleName} --quantMode GeneCounts --
                  """
                }

}
