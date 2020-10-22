process star_idx {
  tag "$genome"
  cpus "$params.cpus"
  publishDir path: { params.save_index ? "${params.out}/star_index" : params.out },
             mode: 'copy', saveAs: { params.save_index ? it : null }
  conda "$baseDir/conda-envs/star-env.yaml"

  input:
    path genome
    path gtf
  output:
    path "${genome.simpleName}_idx"
  script:
    """
		mkdir ${genome.simpleName}_idx
    STAR --runMode genomeGenerate \
    --runThreadN ${task.cpus} \
    --genomeDir ${genome.simpleName}_idx \
    --genomeFastaFiles $genome \
    --sjdbGTFfile $gtf \
    --sjdbGTFfeatureExon "${params.sjdbGTFfeatureExon}"
    """
}
