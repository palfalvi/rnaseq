process rsem_idx {
  tag "$bam"
  label 'small'

  publishDir path: { params.save_index ? "${params.out}/rsem_index" : params.out },
              mode: 'copy', saveAs: { params.save_index ? it : null }

  conda "$baseDir/conda-envs/rsem-env.yaml"
  container "quay.io/biocontainers/rsem"

	input:
	  path genome
	  path gtf

  output:
    path "*_rsem_idx*"

  script:

    def paired = params.single ? "" : "-p"
    """
    rsem-prepare-reference \
    --num-threads $task.cpus \
    --gtf $gtf \
    $genome \
    ${genome.simpleName}_rsem_idx
    """

}
