// ABRA2 can also be run on a spliced-aware bam, but this may require specific parameterization which, for now, is not managed here!

process ABRA2 {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/abra2:2.24--hdcf5f25_3"

  label 'cpu_low'
  label 'mem_low'

  input:
  tuple val(meta), path(bam, stageAs: 'input_raw/*'), path(bai, stageAs: 'input_raw/*')
  tuple val(meta2), path(ref_fa, arity: 1, stageAs: 'input_ref/*')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}_abra2.bam"), path("${meta.label ?: meta.id}_abra2.bai"), emit: bam

  script:

  """
  #!/usr/bin/bash

  ## other parameters:  interested: --mapq, --single, --mad, (--undup ??) / to check: --contigs, --sa ???, --mac, --mmr, --no-edge-ci, --no-ndn, --rcf, --sua/--ssc/-sobs

  abra2 --in ${bam} \\
    --out ${meta.label ?: meta.id}_abra2.bam \\
    --ref ${ref_fa} \\
    --threads ${task.cpus} \\
    --index \\
    ${ (meta?.read_type == 'SE') ? '--single' : ''} \\
    ${task.ext.args ?: ''}
  """

  stub:

  """
  #!/usr/bin/bash

  touch ${meta.label ?: meta.id}_abra2.bam ${meta.label ?: meta.id}_abra2.bai
  """
}
