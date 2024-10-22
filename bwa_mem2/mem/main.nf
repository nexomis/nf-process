// IN PROGRESS !
process BWA_MEM2 {
  container "quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_3"  // ???

  label 'cpu_high'
  label 'mem_8GB'

  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  tuple val(meta2), path(idx, arity: 1, stageAs: 'input_ref/*')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}.sam")

  script:
  """
  #!/usr/bin/bash

  idx_w_prefix=\$(ls -d ${idx}/*.bwt.2bit.64 | sed "s/\\.bwt\\.2bit\\.64\$//")
  if [ -z "\$idx_w_prefix" ]; then
    echo "BWA-MEM2 index not found. Please ensure the index files are in '${idx}'" 1>&2
    exit 1
  fi

  bwa-mem2 mem \\ 
    -t $task.cpus \\ 
    -R "@RG\\tID:${meta.label ?: meta.id}\\tSM:${meta.label ?: meta.id}\\tPL:ILLUMINA\\tLB:${meta.label ?: meta.id}" \\ 
    \$idx_w_prefix \\ 
    ${ (reads.size() == 1) ? "${reads}" : "${reads[0]} ${reads[1]}" } \\ 
    > ${meta.label ?: meta.id}.sam
  """

  stub:
  """
  #!/usr/bin/bash

  touch ${meta.label ?: meta.id}.sam
  """
}
