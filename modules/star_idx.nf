

process star_idx {
        cpus params.mapping.cpus
        conda 'star=2.7.6a'

        input:
		path genome
                path transcriptome
        output:
                path "${genome.simpleName}_idx"
        script:
                """
		mkdir ${genome.simpleName}_idx
                STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir ${genome.baseName}_idx --genomeFastaFiles $genome --sjdbGTFfile $transcriptome --sjdbGTFfeatureExon "${params.sjdbGTFfeatureExon}"
                """
}
