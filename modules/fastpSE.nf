process run_fastpSE {
tag "$sample_id"
cpus "$params.fastqc.cpus"

conda "$baseDir/conda-envs/trim-env.yaml"

when:
  !params.skip_qc
input:
  file reads
output:
  tuple file("trim_${reads}"), emit: trimmed
  file "${sample_id}_fastp.json", emit: json
script:
  """
  fastp \
  -w ${task.cpus} \
  -i ${reads} \
  -o trim_${reads} \
  --json ${sample_id}_fastp.json
  """
}
