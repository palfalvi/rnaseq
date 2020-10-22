process star_idx {
  tag "$genome"
  cpus "$params.cpus"
  publishDir "${params.out}/star_index", mode: 'copy', enable: "${params.save_index}"
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
