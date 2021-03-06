process kallisto_quantSE {
	tag "$reads.simpleName"
	label 'small'
	publishDir "${params.out}/kallisto", mode: 'copy'
	cpus "$params.cpus"
  conda "$baseDir/conda-envs/kallisto-env.yaml"

  input:
  	path transcriptome_idx
    path reads
  output:
    path "$reads.simpleName"
  script:
    """
    kallisto quant \
		-t ${task.cpus} \
		--single \
		-l ${params.fragment_length} \
		-s ${params.fragment_sd} \
		-i $transcriptome_idx \
		-o $reads.simpleName ${reads} \
		&> ${reads.simpleName}.log

    mv ${reads.simpleName}.log ${reads.simpleName}/
    """
}
