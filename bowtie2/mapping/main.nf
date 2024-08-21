
process BOWTIE2 {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/bowtie2:2.5.4--he20e202_3"

  label 'cpu_high'
  label 'mem_high'

  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  tuple val(meta_idx), path(idx, arity: 1, stageAs: 'input_ref/*')

  output:
  tuple val(meta), path("${meta.id}.sam") , emit: sam
  tuple val(meta), path("${meta.id}.log") , emit: log

  script:

  // TODO: check and manage if necessary option to extract unmapped read on fastq file (maybe not via task.ext if depend of reads are SE or PE): '--un-conc[-gz]' '--un[-gz]'
  // in addition with extraction of unmapped reads, the default settings can be '--no-unal' (exclusion of unmapped reads from output SAM)

  """
  #!/usr/bin/bash

  # index prefix (in case of error: check if 'ls -d {idx}/{pattern}' is sufficent with symbolic links, whether add level 'ls -d {idx}/*/{pattern}' or use find use 'max-depth=2' and '-L')
  idx_w_prefix=\$(ls -d ${idx}/*.rev.1.bt2 | sed "s/\\.rev\\.1\\.bt2\$//")
  if [ -z "\$idx_w_prefix" ]; then
    idx_w_prefix=\$(ls -d ${idx}/*.rev.1.bt2l | sed "s/\\.rev\\.1\\.bt2l\$//")
  fi
  if [ -z "\$idx_w_prefix" ]; then
    echo "Not found bowtie2 index files: zero math of sufix file search ('.rev.1.bt2' and '.rev.1.bt2l') on directories '${idx}'" 1>&2
    exit 1
  fi

  bowtie2 \\
    ${task.ext.args_mapping_mode ?: '--end-to-end'} \\
    -x \$index_w_bn \\
    ${ (reads.size() == 1) ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}" } \\
    --threads $task.cpus \\
    ${task.ext.args ?: ''} \\
    -S ${meta.id}.sam \\
    2 > ${meta.id}.log
  """

  stub:

  """
  #!/usr/bin/bash

  touch ${meta.id}.sam ${meta.id}.log
  """
}