process run_fastqcSE {
        cpus "$params.fastqc.cpus"
        
	publishDir "${params.out}/fastqc", mode: 'copy'

        conda 'bioconda::fastqc=0.11.9 perl-app-cpanminus'
        
        when:
		!params.skip_qc
	input:
                path reads
        output:
                path 'fastqc_*'
        script:
                """
                mkdir fastqc_${reads.simpleName}
		fastqc -t ${task.cpus} -o fastqc_${reads.simpleName} $reads
                """
}
