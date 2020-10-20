process run_fastqc {
        tag "$sample_id"
        cpus "$params.fastqc.cpus"
        
	publishDir "${params.out}/fastqc", mode: 'copy'

        conda 'bioconda::fastqc=0.11.9 perl-app-cpanminus'
        
	when:
                !params.skip_qc
	input:
                tuple val(sample_id), file(reads)
        output:
                path 'fastqc_*'
        script:
                """
                mkdir fastqc_${sample_id}
		fastqc -t ${task.cpus} -o fastqc_${sample_id} $reads
                """
}
