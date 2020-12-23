process run_multiqc {
  conda "$baseDir/conda-envs/multiqc-env.yaml"
 	publishDir "${params.out}", mode: 'move'
  when:
    !params.skip_multiqc
	input:
		path('*')
		path config
  output:
    path 'multiqc*html'
  script:
    """
    export LC_ALL=en_US.utf8
    multiqc $launchDir/$config/
    """
}
