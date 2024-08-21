process QUAST {
    input:
    tuple val(meta), path(assembly)     // abacas scaffold
    tuple val(meta_2), path(bam)         // sorted bam (cleanedreads on abacas scaffold) (bai ??)
    tuple val(meta_3), path(reference_fasta)   // ref_genome (optionnel)
    // include spades ??

    output:
    ptuple val(meta), path("${meta.id}/meta.idreport.html"), emit: quast_report  // html
    tuple val(meta), path("${meta.id}")                    , emit: quast_output  // all output
    tuple val(meta), path("${meta.id}.log")                , emit: log

    script:
    def ref_flag = reference_fasta ? "--reference ${reference_fasta}" : ""

    """
    mkdir -p ${meta.id}

    # Exécuter QUAST pour évaluer la qualité de l'assemblage et du mapping
    quast.py ${assembly} \\
        -o ${meta.id} \\
        --rna-finding \\
        --glimmer \\
        --bam ${bam} \\
        ${ref_flag}

    # QUAST génère plusieurs fichiers, nous récupérons le rapport HTML comme sortie principale
    """
}
