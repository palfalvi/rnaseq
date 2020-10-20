process star_index {
        tag "$genome"
        cpus "$params.mapping.cpus"
        conda 'star=2.7.6a'

        input:
		path genome
                path gtf
        output:
                path "${genome.simpleName}_idx"
        script:
                """
		mkdir ${genome.simpleName}_idx
                STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir ${genome.simpleName}_idx --genomeFastaFiles $genome --sjdbGTFfile $gtf --sjdbGTFfeatureExon "${params.sjdbGTFfeatureExon}"
                """
}
