process salmon_quantSE {
	tag "${reads.simpleName}"
		
	publishDir "${params.out}/salmon", mode: 'copy'

	cpus "$params.mapping.cpus"
        conda 'salmon=1.3.0'

        input:
                path transcriptome_idx
                path reads
        output:
                path "$reads.simpleName"
        script:
                """
                salmon quant --libType A --validateMappings -p ${task.cpus} -i $transcriptome_idx -r ${reads} -o $reads.simpleName
                """
}
