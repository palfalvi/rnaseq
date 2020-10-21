process collect_star {
  tag "$bam"
  cpus "$params.fastqc.cpus"

	publishDir "${params.out}/featureCounts", mode: 'copy'

  conda "$baseDir/conda-envs/subread-env.yaml"

	input:
	  path bam
	  path gtf
  output:
    path "${bam}_gene.featureCounts.txt*"
  script:
    """
    featureCounts -p -T ${task.cpus} -s ${params.featureCounts_direction} -a $gtf -o ${bam}_gene.featureCounts.txt ${bam}/*.out.bam
    """

}
