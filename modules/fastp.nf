process run_fastp {
tag "$sample_id"
cpus "$params.cpus"
publishDir "${params.out}/fastp_qc", mode: 'copy', pattern: '*.json'
conda "$baseDir/conda-envs/trim-env.yaml"

input:
  tuple val(sample_id), file(reads)
output:
  tuple val(sample_id), file("trim_*"), emit: trimmed
  path "*.json", emit: json
script:
  """
  fastp \
  -w ${task.cpus} \
  -i ${reads[0]} \
  -I ${reads[1]} \
  -o trim_${reads[0]} \
  -O trim_${reads[1]} \
  --detect_adapter_for_pe \
  --json ${sample_id}_fastp.json
  """
}
