

process kallisto_idx {
        tag "$transcriptome"
        publishDir "${params.out}/kallisto_index", mode: 'copy', enable: "${params.save_index}"
        conda 'bioconda::kallisto=0.46.2'

        input:
                path transcriptome
        output:
                path 'index'
        script:
                """
                kallisto index -i index $transcriptome
                """
}
