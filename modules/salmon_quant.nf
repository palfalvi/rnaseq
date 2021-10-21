process salmon_quant {
	tag "${sample_id}"
	errorStrategy 'ignore'

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
		def inputs = params.ont ? "-a $reads"             : ( params.single ? "-r $reads --validateMappings" : "-1 ${reads[0]} -2 ${reads[1]} --validateMappings" )
		def idx    = params.ont ? "-t $transcriptome_idx" : "-i $transcriptome_idx"

    """
    salmon quant \
		--libType A \
		-p ${task.cpus} \
		$idx \
		$inputs \
		-o $sample_id
    """
}
