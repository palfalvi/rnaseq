process run_fastp {
tag "$sample_id"

label 'small_plus'

publishDir "${params.out}/fastp_qc", mode: 'copy', pattern: '*.json'
publishDir path: { params.save_trimmed ? "${params.out}/trimmed" : params.out },
            mode: 'copy', saveAs: { params.save_index ? it : null },
            pattern: 'trim_*'

conda "$baseDir/conda-envs/trim-env.yaml"

input:
  tuple val(sample_id), file(reads)
  
output:
  tuple val(sample_id), file("trim_*"), optional: true, emit: trimmed
  path "*.json", emit: json

script:

  def readfiles  = params.single    ? "-i $reads"      : "-i ${reads[0]} -I ${reads[1]}"
  def outfiles   = params.single    ? "-o trim_$reads" : "-o trim_${reads[0]} -O trim_${reads[1]}"
  def outff      = params.skip_trim ? ""               : "$outfiles"
  def adapter    = params.single    ? ""               : "--detect_adapter_for_pe"

  """
  fastp \
  -w ${task.cpus} \
  $readfiles \
  $outff \
  $adapter \
  --overrepresentation_analysis \
  --json ${sample_id}_fastp.json
  """
}
