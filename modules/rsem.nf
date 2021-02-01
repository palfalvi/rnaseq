process rsem {
  tag "$bam"
  label 'small'

	publishDir "${params.out}/rsem", mode: 'copy'

  conda "$baseDir/conda-envs/rsem-env.yaml"

	input:
	  path bam
	  path gtf
    path rsem_idx

  output:
    path "${bam}_rsem*"

  script:

    def paired = params.single ? "" : "--paired-end"
    """
    rsem-calculate-expression \
    --alignments
    --num-threads $task.cpus \
    $paired \
    --no-bam-output \
    --forward-prob 0.5 \
    --no-qualities \
    --bam \
    ${bam} \
    $rsem_idx \
    ${bam.simpleName}_rsem
    """

}
