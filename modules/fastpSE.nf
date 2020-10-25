process run_fastpSE {
tag "$sample_id"
label 'small_plus'
cpus "$params.cpus"
publishDir "${params.out}/fastp_qc", mode: 'copy', pattern: '*.json'
conda "$baseDir/conda-envs/trim-env.yaml"

when:
  !skip_trim
input:
  file reads
output:
  path "trim_*", emit: trimmed
  path "*.json", emit: json
script:
  """
  fastp \
  -w ${task.cpus} \
  -i ${reads} \
  -o trim_${reads} \
  --overrepresentation_analysis \
  --json ${reads.simpleName}_fastp.json
  """
}
