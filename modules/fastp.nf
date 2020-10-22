process run_fastp {
tag "$sample_id"
cpus "$params.fastqc.cpus"

conda "$baseDir/conda-envs/trim-env.yaml"

when:
  !params.skip_qc
input:
  tuple val(sample_id), file(reads)
output:
  tuple val(sample_id), file("trim_${reads[0]}", "trim_${reads[1]}"), emit: reads
  file "${sample_id}_fastp.json", emit: json
script:
  """
  fastp \
  -w ${task.cpus} \
  -i ${reads[0]} \
  -I ${reads[1]} \
  -o trim_${reads[0]} \
  -O trim_${reads[1]} \
  --detect_adapter_for_pe \
  -json ${sample_id}_fastp.json
  """
}
