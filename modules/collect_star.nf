process collect_star {
  tag "$bam"
  cpus "$params.mapping.cpus"
  
	publishDir "${params.out}/featureCounts", mode: 'copy'

  conda 'bioconda::subread=2.0.1'
        
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