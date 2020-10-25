process run_fastp_qc {
tag "$sample_id"
label 'small_plus'
cpus "$params.cpus"
publishDir "${params.out}/fastp_qc", mode: 'copy', pattern: '*.json'
conda "$baseDir/conda-envs/trim-env.yaml"

when:
  !kip_qc
input:
  tuple val(sample_id), file(reads)
output:
  path "*.json", emit: json
script:
  """
  fastp \
  -w ${task.cpus} \
  -i ${reads[0]} \
  -I ${reads[1]} \
  --overrepresentation_analysis \
  --json ${sample_id}_fastp.json
  """
}
