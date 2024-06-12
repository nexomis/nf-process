
process KALLISTO_QUANT {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/kallisto:0.50.1--h6de1650_2"

  label 'cpu_med'
  label 'mem_8G'

  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  path index

  output:
  tuple val(meta), path("$meta.id", type: 'dir'), emit: output_dir
  tuple val(meta), path("${meta.id}.log", type: 'file'), emit: log

  script:
  """
  #!/usr/bin/bash
  kallisto quant --threads ${task.cpus} \\
    --index ${index} \\
    --output-dir ${meta.id} \\
    ${meta.kallisto_args ?: ''} \\
    ${task.ext.args ?: ''} \\
    ${ (reads.size() == 1) ? '--single' : ''} ${reads} \\
    2> ${meta.id}.log

  rm ${meta.id}/abundance.tsv
    
  """

  stub:
  """
  #!/usr/bin/bash

  mkdir ${meta.id}
  touch ${meta.id}/abundance.tsv ${meta.id}/abundance.h5 ${meta.id}/run_info.json
  """
}

