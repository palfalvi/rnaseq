

process kallisto_idx {
        tag "$transcriptome"
        publishDir "${params.out}/kallisto_index", mode: 'copy', enable: "${params.save_index}"
        conda './conda-envs/kallisto-env.yaml'

        input:
                path transcriptome
        output:
                path 'index'
        script:
                """
                kallisto index -i index $transcriptome
                """
}
