process run_fastqcSE {
        tag "$reads.simpleName"
        cpus "$params.fastqc.cpus"

	publishDir "${params.out}/fastqc", mode: 'copy'

        conda "$baseDir/conda-envs/fastqc-env.yaml"

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
