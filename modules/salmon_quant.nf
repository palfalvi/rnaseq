process salmon_quant {
	tag "${sample_id}"

	label 'small'
	publishDir "${params.out}/salmon", mode: 'copy'

  conda "$baseDir/conda-envs/salmon-env.yaml"
	container "quay.io/biocontainers/salmon"

  input:
  	path transcriptome_idx
    tuple val(sample_id), file(reads)

  output:
    path sample_id

  script:
	  def inputs = params.single ? "-r $reads"      : "-1 ${reads[0]} -2 ${reads[1]}"

    """
    salmon quant \
		--libType A \
		--validateMappings \
		-p ${task.cpus} \
		-i $transcriptome_idx \
		$inputs \
		-o $sample_id
    """
}
