

process kallisto_idx {
  tag "$transcriptome"
  publishDir path: { params.save_index ? "${params.out}/kallisto_index" : params.out },
             mode: 'copy', saveAs: { params.save_index ? it : null }
  conda "$baseDir/conda-envs/kallisto-env.yaml"

  input:
    path transcriptome
  output:
    path 'index'
  script:
    """
    kallisto index -i index $transcriptome
    """
}
