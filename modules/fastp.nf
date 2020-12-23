process run_fastp {
tag "$sample_id"
label 'small_plus'
cpus "$params.cpus"
publishDir "${params.out}/fastp_qc", mode: 'copy', pattern: '*.json'
publishDir path: { params.save_trimmed ? "${params.out}/trimmed" : params.out },
            mode: 'copy', saveAs: { params.save_index ? it : null }, pattern: '*.fastq.gz'
conda "$baseDir/conda-envs/trim-env.yaml"

when:
  !skip_trim
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
  --overrepresentation_analysis \
  --json ${sample_id}_fastp.json
  """
}
