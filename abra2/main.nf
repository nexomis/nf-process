// ABRA2 can also be run on a spliced-aware bam, but this may require specific parameterization which, for now, is not managed here!
// TODO: management of the case where the bam input can be greather than 1 sample: management of the bam and bai lists, and addition of the suffix "_abra2".

process ABRA2 {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/abra2:2.24--hdcf5f25_3"

  label 'cpu_high'
  label 'mem_high'

  input:
  tuple val(meta), path(bam, stageAs: 'input_raw/*'), path(bai, stageAs: 'input_raw/*')
  tuple val(meta2), path(ref_fa, arity: 1, stageAs: 'input_ref/*')

  output:
  tuple val(meta), path("${meta.id}_abra2.bam"), path("${meta.id}_abra2.bai"), emit: bam

  script:
  //bam_out = "${bam[0].parent}/${bam[0].getBaseName()}_abra2.${bam[0].getExtension()}"


  """
  #!/usr/bin/bash

  ## parameters:
  # interested: --mapq, --single, --mad, (--undup ??)
  # to check: --contigs, --sa ???, --mac, --mmr, --no-edge-ci, --no-ndn, --rcf, --sua/--ssc/-sobs
  # note in case of splice-aware mapping, check '--help'!

  abra2 --in ${bam} \\
    --out ${meta.id}_abra2.bam \\
    --ref ${ref_fa} \\
    --threads ${task.cpus} \\
    --index \\
    ${ (meta?.read_type == 'SE') ? '--single' : ''} \\
    ${task.ext.args ?: ''}
  """

  stub:

  """
  #!/usr/bin/bash

  touch ${meta.id}_abra2.bam ${meta.id}_abra2.bam.bai
  """
}
