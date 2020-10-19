process run_multiqc {
  conda 'multiqc=1.9'
 	publishDir "${params.out}", mode: 'move'       

	input:
		path('*')
		path config 
  output:
    path 'multiqc*html'
  script:
                """
                multiqc $config/.
                """
}
