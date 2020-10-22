process salmon_quantSE {
	tag "${reads.simpleName}"
	publishDir "${params.out}/salmon", mode: 'copy'
	cpus "$params.cpus"
  conda "$baseDir/conda-envs/salmon-env.yaml"

  input:
  	path transcriptome_idx
    path reads
  output:
    path "$reads.simpleName"
  script:
    """
    salmon quant \
		--libType A \
		--validateMappings \
		-p ${task.cpus} \
		-i $transcriptome_idx \
		-r ${reads} \
		-o $reads.simpleName
    """
}
