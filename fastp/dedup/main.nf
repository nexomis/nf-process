process FASTP_DEDUP {
  container 'staphb/fastp:0.23.4'
  tag "$meta.id"
  cpus 8
  memory 32.GB

  input:
  tuple val(meta), path(files, arity: 1..2, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.id}*.fq.gz", arity: 1..2), emit: reads
  tuple val(meta), path("${meta.id}.json", arity: 1), emit: report

  script:
  def default_args = "-Q -L -A -G --dedup --dup_calc_accuracy 6"

  def in_args = "-i " + files[0]
  def out_args = " -o "
  if (files.size() > 1) {
    in_args += " -I " + files[1]
    in_args += " --detect_adapter_for_pe"
    out_args += meta.id + "_R1.fq.gz"
    out_args += " -O " + meta.id + "_R2.fq.gz"
  } else {
    out_args += meta.id + ".fq.gz"
  }
  
  """
  #!/usr/bin/bash

  fastp --thread $task.cpus ${task.ext.args ?: default_args} \\
  --json ${meta.id}.json \\
  $in_args $out_args

  """

  stub:
  def args_out = " "
  if (files.size() > 1) {
    args_out += meta.id + "_R1.fq.gz"
    args_out += " " + meta.id + "_R2.fq.gz"
  } else {
    args_out += meta.id + ".fq.gz"
  }
  """
  #!/usr/bin/bash

  touch $args_out
  touch ${meta.id}.json

  """

}
