

process salmon_idx {
        cpus params.mapping.cpus
        conda 'salmon=1.3.0'

        input:
                path transcriptome
        output:
                path 'index'
        script:
                """
                salmon index -p ${task.cpus} -t $transcriptome -i index
                """
}
