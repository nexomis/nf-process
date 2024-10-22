process PICARD_MARK_DUPLICATES {
  container "quay.io/biocontainers/picard:3.3.0--hdfd78af_0"    // 615 MiB !!

  label 'cpu_x1'
  label 'mem_med'

  input:
  tuple val(meta), path(bam, arity: 1, stageAs: 'input_raw/*')      // Must be coordinate sorted (TODO: does not use index 'bai' for save time?) and contains minimal read group information (cf. GATK recquirement)
  tuple val(meta2), path(ref_fa, arity: 1, stageAs: 'input/*')      // Improve duplicate detection ?? make-it optional ?

  output:
  tuple val(meta), path("${meta.label ?: meta.id}_RmDup.bam"), path("${meta.label ?: meta.id}_RmDup.bai") , optional:false , emit: bam_bai
  tuple val(meta), path("${meta.label ?: meta.id}_RmDup.txt") , optional:false , emit: metrics


  script:

  """
  #!/usr/bin/bash

  picard MarkDuplicates \\
    -I ${bam} \\
    -O ${meta.label ?: meta.id}_RmDup.bam \\
    -M ${meta.label ?: meta.id}_RmDup.txt \\
    -R ${ref_fa} \\
    --CREATE_INDEX true \\
    ${task.ext.args ?: '--REMOVE_DUPLICATES true'}

  """

  stub:  
  """
  #!/usr/bin/bash

  touch ${meta.label ?: meta.id}_RmDup.bam ${meta.label ?: meta.id}_RmDup.txt
  """ 
}