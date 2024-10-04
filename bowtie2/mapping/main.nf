
process BOWTIE2 {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/bowtie2:2.5.4--he20e202_3"

  label 'cpu_high'
  label 'mem_8GB' // 4GB max of foot print: https://hpc.nih.gov/apps/bowtie2.html

  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  tuple val(meta2), path(idx, arity: 1, stageAs: 'input_ref/*')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}.sam")

  script:
  """
  #!/usr/bin/bash

  idx_w_prefix=\$(ls -d ${idx}/*.rev.1.bt2 | sed "s/\\.rev\\.1\\.bt2\$//")
  if [ -z "\$idx_w_prefix" ]; then
    idx_w_prefix=\$(ls -d ${idx}/*.rev.1.bt2l | sed "s/\\.rev\\.1\\.bt2l\$//")
  fi
  if [ -z "\$idx_w_prefix" ]; then
    echo "Not found bowtie2 index files: zero math of sufix file search ('.rev.1.bt2' and '.rev.1.bt2l') on directories '${idx}'" 1>&2
    exit 1
  fi

  bowtie2 \\
    -x \$idx_w_prefix \\
    ${ (reads.size() == 1) ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}" } \\
    --threads $task.cpus \\
    -S ${meta.label ?: meta.id}.sam \\
    ${task.ext.args ?: ''}
  """

  stub:

  """
  #!/usr/bin/bash

  touch ${meta.label ?: meta.id}.sam
  """
}