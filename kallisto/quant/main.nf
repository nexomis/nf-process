
process KALLISTO_QUANT {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/kallisto:0.50.1--h6de1650_2"
  tag "$meta.id"
  label 'cpu_med'
  label 'mem_8G'

  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  tuple val(meta2), path(index, stageAs: "input_index/index")

  output:
  tuple val(meta), path("${meta.label ?: meta.id}/abundance.h5", type: 'file'), emit: h5
  tuple val(meta), path("${meta.label ?: meta.id}.log", type: 'file'), emit: log

  script:
  def name=meta.label ?: meta.id
  """
  #!/usr/bin/bash
  
  kallisto quant --threads ${task.cpus} \\
    --index ${index} \\
    --output-dir ${name} \\
    ${meta.kallisto_args ?: ''} \\
    ${task.ext.args ?: ''} \\
    ${ (reads.size() == 1) ? '--single' : ''} ${reads} \\
    2> ${name}.log

  """

  stub:
  """
  #!/usr/bin/bash

  mkdir ${meta.id}
  touch ${meta.id}/abundance.tsv ${meta.id}/abundance.h5 ${meta.id}/run_info.json
  """
}

