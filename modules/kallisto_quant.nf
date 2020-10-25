process kallisto_quant {
	tag "${sample_id}"
	label 'small'
	publishDir "${params.out}/kallisto", mode: 'copy'
	cpus "$params.cpus"
  conda "$baseDir/conda-envs/kallisto-env.yaml"

  input:
  	path transcriptome_idx
    tuple val(sample_id), file(reads)
  output:
    path sample_id
  script:
  	"""
    kallisto quant \
		-t ${task.cpus} \
		-i $transcriptome_idx \
		-o $sample_id ${reads[0]} ${reads[1]} \
		&> ${sample_id}.log

    mv ${sample_id}.log ${sample_id}/
    """
}
