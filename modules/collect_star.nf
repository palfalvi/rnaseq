process collect_star {
  tag "$bam"
  cpus "$params.fastqc.cpus"
  
	publishDir "${params.out}/featureCounts", mode: 'copy'

  conda 'bioconda::subread=2.0.1'
        
	input:
	  path bam
	  path gtf
  output:
    path "${bam}_gene.featureCounts.txt"
    path "${bam}_gene.summary"
  script:
  if ("$params.single") 
    """
    featureCounts -T ${task.cpus} -s ${params.featureCounts_direction} -a $gtf -o ${bam}_gene.featureCounts.txt ${bam}/*.out.bam 2> ${bam}_gene.featureCounts.summary
    """
  else 
    """
    featureCounts -p -T ${task.cpus} -s ${params.featureCounts_direction} -a $gtf -o ${bam}_gene.featureCounts.txt ${bam}/*.out.bam 2> ${bam}_gene.featureCounts.summary
    """
  
}