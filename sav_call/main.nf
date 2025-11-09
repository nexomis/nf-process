process SAV_CALL {
  container "quay.io/nexomis/sav_call:0.2.0-py"
  tag "$meta.id"
  cpus 1
  memory 15.GB

  input:
  tuple val(meta), path(bam, arity: 1..2, stageAs: 'input/')        // bam[0] is bam file [recquired] and bam[1] is its index (bai) [by default optional, but it may be necessary (and therefore required) on other sav_call version (depending on particular pileup options: in particular region-specific options)]
  tuple val(meta2), path(ref_fa, arity: 1, stageAs: 'input/')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}_snv.both.csv", type: 'file') , optional:false , emit: snv_both
  tuple val(meta), path("${meta.label ?: meta.id}_snv.rev.csv", type: 'file') , optional:false , emit: snv_rev
  tuple val(meta), path("${meta.label ?: meta.id}_snv.fwd.csv", type: 'file') , optional:false , emit: snv_fwd
  tuple val(meta), path("${meta.label ?: meta.id}.raw_indel.csv", type: 'file') , optional:false , emit: snv_raw_indel
  tuple val(meta), path("${meta.label ?: meta.id}.base.csv", type: 'file') , optional:false , emit: base_comp

  script:
  """
  #!/usr/bin/bash
  
  sav_call \\
    --prefix-out ${meta.label ?: meta.id}_snv \\
    --bam ${bam[0]} \\
    --reference ${ref_fa} \\
    --base-csv ${meta.label ?: meta.id}.base.csv \\
    --indel-csv ${meta.label ?: meta.id}.raw_indel.csv \\
    ${task.ext.sav_call_args ?: ''}
  """

  stub:  
  """
  #!/usr/bin/bash

  touch ${meta.label ?: meta.id}_snv.both.csv ${meta.label ?: meta.id}_snv.rev.csv ${meta.label ?: meta.id}_snv.fwd.csv ${meta.label ?: meta.id}.raw_indel.csv ${meta.label ?: meta.id}.base.csv
  """ 
}
