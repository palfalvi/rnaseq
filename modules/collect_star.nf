process collect_star {
  cpus params.fastqc.cpus
        
	publishDir "${params.out}/featureCount", mode: 'copy'

  conda 'bioconda::subread=2.0.1'
        
	input:
	  path sam
	  path gtf
  output:
    path "${sam}_gene.featureCounts.txt"
    path "${sam}_gene.summary"
  script:
    """
    featureCounts -T ${task.cpus} -a $gtf -o ${sam}_gene.featureCounts.txt -p ${sam}/*.out.sam 2> ${sam}_gene.summary
    """
}