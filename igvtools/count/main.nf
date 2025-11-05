// need to be valided in case of multifasta reference (segmenteed virus, contigs, ...)

process IGVTOOLS_COUNT {
  container "quay.io/biocontainers/igvtools:2.3.93--0"
  tag "$meta.id"
  label 'cpu_x1'
  label 'mem_med'

  input:
  tuple val(meta), path(bam, arity: 1..2, stageAs: 'input/')   // bam[0] is bam file [recquired] and bam[1] is its index (bai) [optional but save time if its provided]
  path(genom_ref_fa, arity: 1..2, stageAs: 'input/')           // genom_ref_fa[0] is fasta file [recquired] and genom_ref_fa[1] is its index (fasta.fai) [optional but save time if its provided in case of large genome]
  val by_strand

  output:
  tuple val(meta), path("${meta.id}.wig", type: 'file') , optional:false , emit: wig

  script:

  """
  #!/usr/bin/bash

  # TODO: interest option: --includeDuplicates / --strands 'read'/'first' / --minMapQuality 0
  igvtools count \\
    -z 1 \\
    -w 1 \\
    --bases \\
    ${meta.args_igvcount ?: ''} \\
    ${task.ext.args} \\
    ${bam[0]} \\
    ${meta}.wig \\
    ${genom_ref_fa[0]}
  """

  stub:  
  """
  #!/usr/bin/bash

  touch ${meta.id}.wig
  """
}