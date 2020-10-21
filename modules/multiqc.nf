process run_multiqc {
  conda "$baseDir/conda-envs/multiqc-env.yaml"
 	publishDir "${params.out}", mode: 'move'

	input:
		path('*')
		path config
  output:
    path 'multiqc*html'
  script:
                """
                multiqc .
                """
}
