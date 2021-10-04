process minimap2_idx {
  tag "$sample_id"
  label 'long_job'

  conda "$baseDir/conda-envs/minimap-env.yaml"

  // publishDir "${params.outdir}/bwa", mode: 'copy'

  input:
    path index
    tuple val(sample_id), file(reads)

  output:
    path "*.mmi", emit: idx

  script:
    """
    minimap2 -d ${index.simpleName}.mmi -t ${task.cpus} $index
    """
}
