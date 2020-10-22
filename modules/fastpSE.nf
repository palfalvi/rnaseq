process run_fastpSE {
tag "$sample_id"
cpus "$params.fastqc.cpus"
publishDir "${params.out}/fastp_qc", mode: 'copy', pattern: '*.json'
conda "$baseDir/conda-envs/trim-env.yaml"

input:
  file reads
output:
  file "trim_*", emit: trimmed
  file "*.json", emit: json
script:
  """
  fastp \
  -w ${task.cpus} \
  -i ${reads} \
  -o trim_${reads} \
  --json ${sample_id}_fastp.json
  """
}
