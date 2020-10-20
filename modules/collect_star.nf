process collect_star {
  cpus params.fastqc.cpus
        
	publishDir "${params.out}/featureCount", mode: 'copy'

  conda 'bioconda::subread=2.0.1'
        
	input:
	  path sam
	  path gtf
  output:
    path "${bam}_gene.featureCounts.txt"
    path "${bam}_gene.summary"
  script:
    """
    featureCounts -T ${task.cpus} -s ${params.featureCounts_direction} -a $gtf -o ${bam}_gene.featureCounts.txt -p ${bam}/*.out.bam 2> ${bam}_gene.summary
    """
}