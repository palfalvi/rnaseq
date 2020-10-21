process kallisto_quantSE {
	tag "$reads.simpleName"

	publishDir "${params.out}/kallisto", mode: 'copy'

	cpus "$params.mapping.cpus"
        conda 'bioconda::kallisto=0.46.2'

        input:
                path transcriptome_idx
                path reads
        output:
                path "$reads.simpleName"
        script:
                """
                kallisto quant -t ${task.cpus} --single -l ${params.fragment_length} -s ${params.fragment_sd} -i $transcriptome_idx -o $reads.simpleName ${reads} &> ${reads.simpleName}/${reads.simpleName}.log
                mv ${reads.simpleName}.log ${reads.simpleName}/
                """
}
