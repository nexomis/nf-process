// assembly using 'spades' with short reads : 1 input fastq[.gz] file or 2 in case of paird-end.
// (even if compatible with the 'spades' tool, cases with several independent input files are not handled here, this inclued additional librarie or hybrid approach).

process SPADES {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/spades:4.0.0--h5fb382e_2"

  label 'cpu_high'
  label 'mem_28G'

  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("$meta.id", type: 'dir'), optional:false, emit: output_dir
  tuple val(meta), path("${meta.id}.log", type: 'file'), optional:false, emit: log
//  tuple val(meta), path('${meta.id}/${meta.id}.scaffolds.fa.gz'), optional:false, emit: scaffolds_gz
//  tuple val(meta), path('${meta.id}/${meta.id}.contigs.fa.gz'), optional:false, emit: contigs_gz
//  tuple val(meta), path('${meta.id}/${meta.id}.scaffolds.gfa'), optional:true, emit: gfa


  script:
  """
  #!/usr/bin/bash

  # run SPAdes
  spades.py --threads ${task.cpus} \\
    --memory ${task.mem} \\
    -o ./ \\
    ${ (reads.size() == 1) ? '-s ${reads}' : '-1 ${reads[0]} -2 ${reads[1]}' } \\
    ${meta.spades_args ?: ''} \\
    ${task.ext.args ?: ''} \\
    2> ${meta.id}.log

  # Move and compress files
  mkdir ${meta.id}
  mv scaffolds.fasta ${meta.id}/${meta.id}.scaffolds.fa && gzip -n ${meta.id}/${meta.id}.scaffolds.fa
  mv contigs.fasta ${meta.id}/${meta.id}.contigs.fa && gzip -n ${meta.id}/${meta.id}.contigs.fa
  mv assembly_graph_with_scaffolds.gfa ${meta.id}/${meta.id}.scaffolds.gfa

  """

  stub:
  """
  #!/usr/bin/bash

  mkdir ${meta.id}
  touch ${meta.id}.log ${meta.id}/${meta.id}.scaffolds.fa.gz ${meta.id}/${meta.id}.contigs.fa.gz ${meta.id}/${meta.id}.scaffolds.gfa
  """
}

