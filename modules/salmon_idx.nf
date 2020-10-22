process salmon_idx {
  tag "$transcriptome"
  cpus "$params.cpus"
  publishDir "${params.out}/salmon_index",
              mode: 'copy', saveAs: { params.save_index ? it : null },
  conda "$baseDir/conda-envs/salmon-env.yaml"

  input:
    path transcriptome
  output:
    path 'index'
  script:
    """
    salmon index -p ${task.cpus} -t $transcriptome -i index
    """
}
