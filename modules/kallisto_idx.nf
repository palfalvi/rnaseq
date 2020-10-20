

process kallisto_idx {
        tag "$transcriptome"
        
        conda 'kallisto=0.46.2'

        input:
                path transcriptome
        output:
                path 'index'
        script:
                """
                kallisto index -i index $transcriptome
                """
}
