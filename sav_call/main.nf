process SAV_CALL {
  container "quay.io/nexomis/sav_call:0.2.0"

  label 'cpu_x1'
  label 'mem_med'

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

  // TODO: manage sample specific option using '${meta.sav_call_args ?: ''}' like in kallisto process (e.g. '--R1-strand' and '--R2-strand')?

  """
  #!/usr/bin/bash
  
  sav_call \\
    --prefix-out ${meta.label ?: meta.id}_snv \\
    --bam ${bam[0]} \\
    --reference ${ref_fa} \\
    --R1-strand forward --R2-strand reverse \\
    --base-csv ${meta.label ?: meta.id}.base.csv \\
    --indel-csv ${meta.label ?: meta.id}.raw_indel.csv \\
    ${ task.ext.alt_ratio_threshold ? "--min-freq '${task.ext.alt_ratio_threshold};${task.ext.alt_ratio_threshold};${task.ext.alt_ratio_threshold}'" : '' } \\
    ${ task.ext.min_dp ? "--min-count '${task.ext.min_dp};${task.ext.min_dp};${task.ext.min_dp}'" : '' } \\
    --max-n-pileup ${task.ext.max_n_pileup ?: ''} \\
    ${task.ext.args ?: ''}

  """

  stub:  
  """
  #!/usr/bin/bash

  touch ${meta.label ?: meta.id}_snv.both.csv ${meta.label ?: meta.id}_snv.rev.csv ${meta.label ?: meta.id}_snv.fwd.csv ${meta.label ?: meta.id}.raw_indel.csv ${meta.label ?: meta.id}.base.csv
  """ 
}