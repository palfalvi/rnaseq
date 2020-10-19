process collect_star {
  cpus params.fastqc.cpus
        
	publishDir "${params.out}/featureCount", mode: 'copy'

  conda 'bioconda::subread=2.0.1'
        
	input:
	  path bam
	  path gtf
  output:
    path "${bam}_gene.featureCounts.txt"
    path "${bam}_gene.summary"
  script:
    """
    featureCounts -T ${task.cpus} -a $gtf -o ${bam}_gene.featureCounts.txt -p -s $featureCounts_direction ${bam}/bams/*.out.bam 2> ${bam}_gene.summary
    """
}