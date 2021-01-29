process kallisto_quant {
	tag "${sample_id}"
	
	label 'small'
	publishDir "${params.out}/kallisto", mode: 'copy'

  conda "$baseDir/conda-envs/kallisto-env.yaml"

  input:
  	path transcriptome_idx
    tuple val(sample_id), file(reads)

  output:
    path sample_id

  script:
		def single = params.single ? "--single -l ${params.fragment_length} -s ${params.fragment_sd}" : ""

  	"""
    kallisto quant \
		$single \
		-t ${task.cpus} \
		-i $transcriptome_idx \
		-o $sample_id \
		${reads} \
		&> ${sample_id}.log

    mv ${sample_id}.log ${sample_id}/
    """
}
