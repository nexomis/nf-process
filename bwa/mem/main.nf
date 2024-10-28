process BWA_MEM {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/bwa:0.7.18--he4a0461_1"

  label 'cpu_high'
  label 'mem_med'

  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  tuple val(meta2), path(idx, arity: 1, stageAs: 'input_ref/*')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}.sam")

  script:
  """
  #!/usr/bin/bash

  idx_w_prefix=\$(ls -d ${idx}/*.bwt | sed "s/\\.bwt\$//")
  if [ -z "\$idx_w_prefix" ]; then
    echo "Not found bwa mem index files: zero math of sufix file search ('.bwt') on directories '${idx}'" 1>&2
    exit 1
  fi

  # -M (for PICARD MarkDuplicates compatibility: force single primary by reads, even if reads are splited and mapp independantly. Probably not ideal in all context)
  bwa mem \\
    -t ${task.cpus} \\
    -R "@RG\\tID:${meta.label ?: meta.id}\\tSM:${meta.label ?: meta.id}\\tPL:ILLUMINA\\tLB:${meta.label ?: meta.id}" \\
    ${task.ext.args ?: ''} \\
    \${idx_w_prefix} \\
    ${reads} \\
    > ${meta.label ?: meta.id}.sam
  """

  stub:

  """
  #!/usr/bin/bash

  touch ${meta.label ?: meta.id}.sam
  """
}
