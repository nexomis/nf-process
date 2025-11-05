process FASTP {
  container 'staphb/fastp:0.23.4'
  tag "$meta.id"
  label 'cpu_low'
  label 'mem_8G'

  input:
  tuple val(meta), path(files, arity: 1..2, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.id}*.fq.gz", arity: 1..2), emit: reads
  tuple val(meta), path("${meta.id}.json", arity: 1), emit: report

  script:
  def default_args = "--trim_poly_g"
  default_args += " --cut_right_window_size 4" 
  default_args += " --cut_right_mean_quality 20" 
  default_args += " --cut_right"
  default_args += " --cut_tail_window_size 4" 
  default_args += " --cut_tail_mean_quality 25" 
  default_args += " --cut_tail"
  default_args += " --average_qual 25"
  default_args += " --trim_poly_x"
  default_args += " --length_required 31"

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
