process FASTP {
  container 'staphb/fastp:0.23.4'

  label 'cpu_low'
  label 'mem_8G'

  input:
  tuple val(sample_name), path(files, arity: 1..2, stageAs: 'input_raw/*')
  tuple val(trim_poly_g), val(cut_right_window_size), val(cut_right_mean_qual), val(cut_tail_window_size), val(cut_tail_mean_qual), val(min_avg_qual), val(trim_poly_x), val(min_len)

  output:
  tuple val("${sample_name}"), path("${sample_name}*.fq.gz", arity: 1..2)

  script:
  def s_args = ""
  if (trim_poly_g) {
    s_args += " --trim_poly_g"
  }
  s_args += " --cut_right_window_size " + cut_right_window_size
  s_args += " --cut_right_mean_quality " + cut_right_mean_qual
  s_args += " --cut_right"
  s_args += " --cut_tail_window_size " + cut_tail_window_size
  s_args += " --cut_tail_mean_quality " + cut_tail_mean_qual
  s_args += " --cut_tail"
  s_args += " --average_qual " + min_avg_qual 
  if (trim_poly_x) {
    s_args += " --trim_poly_x"
  }
  def s_args_in = "-i " + files[0]
  def s_args_out = " -o "
  if (files.size() > 1) {
    s_args += " --detect_adapter_for_pe"
    s_args_in += " -I " + files[1]
    s_args_out += sample_name + "_R1.fq.gz"
    s_args_out += " -O " + sample_name + "_R2.fq.gz"
  } else {
    s_args_out += sample_name + ".fq.gz"
  }
  s_args = s_args_in + s_args_out + s_args
  
  """
  #!/usr/bin/bash

  fastp --thread $task.cpus $s_args

  """

  stub:
  def args_out = " "
  if (files.size() > 1) {
    args_out =+ sample_name + "_R1.fq.gz"
    args_out += " " + sample_name + "_R2.fq.gz"
  } else {
    args_out =+ sample_name + ".fq.gz"
  }
  """
  #!/usr/bin/bash

  touch $args_out

  """

}
