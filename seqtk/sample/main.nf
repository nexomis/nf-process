process SEQTK_SAMPLE {
  container 'staphb/seqtk:1.4'

  label 'cpu_low'
  label 'mem_8G'

  input:
  tuple val(meta), path(files, arity: 1..2, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.id}*.fq", arity: 1..2)

  script:
  def num_reads = 100000
  if (task.ext.num_reads) {
    num_reads = task.ext.num_reads
  }
  def base_cmd = "seqtk sample -s42"
  def cmd_reads1 = "${base_cmd} ${files[0]} ${num_reads}  > ${meta.id}.fq"
  def cmd_reads2 = ""
  if (files.size() > 1) {
    cmd_reads1 = "${base_cmd} ${files[0]} ${num_reads}  > ${meta.id}_R1.fq"
    cmd_reads2 = "${base_cmd} ${files[1]} ${num_reads} > ${meta.id}_R2.fq"
  }
  """
  #!/usr/bin/bash

  $cmd_reads1
  
  $cmd_reads2

  """

  stub:
  def cmd_reads1 = "touch ${meta.id}_R1.fq"
  def cmd_reads2 = ""
  if (files.size() > 1) {
    cmd_reads2 = "touch > ${meta.id}_R2.fq"
  }
  """
  #!/usr/bin/bash

  $cmd_reads1
  
  $cmd_reads2

  """

}
