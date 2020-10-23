process run_fastpSE_qc {
tag "$sample_id"
label 'small_plus'
cpus "$params.cpus"
publishDir "${params.out}/fastp_qc", mode: 'copy', pattern: '*.json'
conda "$baseDir/conda-envs/trim-env.yaml"

when:
  !skip_qc
input:
  file reads
output:
  path "*.json", emit: json
script:
  """
  fastp \
  -w ${task.cpus} \
  -i ${reads} \
  --overrepresentation_analysis \
  --json ${reads.simpleName}_fastp.json
  """
}
