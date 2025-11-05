process CONCAT_FILES {
    container "python:3.14.0a7-bookworm"
    tag "$meta.id"
    input:
    tuple val(meta), path(files, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}.fa*"), emit: fasta

    script:
    """
    # Check if all files have the same MIME type
    first_mime=\$(file -b --mime-type ${files[0]})
    
    for file in ${files.join(' ')}; do
        mime_type=\$(file -b --mime-type \$file)
        if [[ "\$mime_type" != "\$first_mime" ]]; then
            echo "Error: Mixed file types detected. All files must be of the same type."
            echo "First file: ${files[0]} is \$first_mime"
            echo "Current file: \$file is \$mime_type"
            exit 1
        fi
    done
    
    # Determine output extension based on MIME type
    if [[ "\$first_mime" == "application/gzip" ]]; then
        cat ${files.join(' ')} > ${meta.id}.fa.gz
    else
        cat ${files.join(' ')} > ${meta.id}.fa
    fi
    """

    stub:
    """
    # Determine output extension based on file extension
    if [[ "${files[0]}" == *.gz ]]; then
        touch ${meta.id}.fa.gz
    else
        touch ${meta.id}.fa
    fi
    """
}
