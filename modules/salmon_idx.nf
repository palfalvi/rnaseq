process salmon_idx {
        tag "$transcriptome"
        cpus "$params.mapping.cpus"
        publishDir "${params.out}", mode: 'copy', enable: "${params.save_index}"
        
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
