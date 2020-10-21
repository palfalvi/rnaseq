process salmon_quant {
	tag "${sample_id}"

	publishDir "${params.out}/salmon", mode: 'copy'

	cpus "$params.mapping.cpus"
        conda 'bioconda::salmon=1.3.0'

        input:
                path transcriptome_idx
                tuple val(sample_id), file(reads)
        output:
                path sample_id
        script:
                """
                salmon quant --libType A --validateMappings -p ${task.cpus} -i $transcriptome_idx -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
                """
}
