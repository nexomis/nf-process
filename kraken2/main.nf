process KRAKEN2 {
  container 'ghcr.io/nexomis/kraken2:2.1.3'

  label 'cpu_med'
  label 'mem_1G_per_cpu' // mem_4G_per_cpu // TODO: params.kraken2_memory

  input:
  tuple val(meta), path(files, arity: 1..2, stageAs: 'input_raw/*')
  path db

  output:
  tuple val(meta), path("*.classified*")       , optional:true, emit: classified_reads_fastq
  tuple val(meta), path("*.unclassified*")     , optional:true, emit: unclassified_reads_fastq
  tuple val(meta), path("*classifiedreads.txt"), optional:true, emit: classified_reads_assignment
  tuple val(meta), path("*.report.txt")        , emit: report
  tuple val(meta), path("*.output.txt")        , emit: output

  // modify ext from config to enable classified and unclassified
  // ext:
  // args: --classified-out classified.#.fq --unclassified-out unclassified.#.fq

  script:
  
  def compress_args = ""
  def fq_ext = ""
  
  def r1_ext = files[0].getExtension().toLowerCase()
  if (r1_ext == 'gz' || r1_ext == 'gzip' || r1_ext == 'z') {
    compress_args = "--gzip-compressed"
    fq_ext = ".fq.gz"
  } else if (r1_ext == '.bz' || r1_ext == '.bz2' || r1_ext == '.bzip2') {
    compress_args = "--bzip2-compressed"
    fq_ext = ".fq.bz2"
  } else {
    compress_args = ""
    fq_ext = ".fq"
  }
  
  def in_args = ""

  // we need to rename reads in case in classify, kraken expect _1 or _2
  def rename_paired_inputs = ""
  
  if (files.size() > 1) {
    rename_paired_inputs = "mkdir input_renamed ; "
    rename_paired_inputs += "ln -sf \$PWD/" + files[0] + " \$PWD/input_renamed/" + meta.id + "_1" + fq_ext + "; "
    rename_paired_inputs += "ln -sf \$PWD/" + files[1] + " \$PWD/input_renamed/" + meta.id + "_2" + fq_ext
    in_args = "--paired input_renamed/" + meta.id + "_1" + fq_ext + " input_renamed/" + meta.id + "_2" + fq_ext
  } else {
    in_args += files[0]
  }
  
  def args_save_unclassified = ''
  if (task.ext.save_unclassified) {
    args_save_unclassified = (files.size() == 1) ? "--unclassified-out ${meta.id}.unclassified.fastq" : "--unclassified-out ${meta.id}.unclassified#.fastq"
  }

  """
  #!/usr/bin/bash

  $rename_paired_inputs

  kraken2 --thread $task.cpus \\
  --output ${meta.id}.output.txt --report ${meta.id}.report.txt \\
  --db $db \\
  ${args_save_unclassified} \\
  ${task.ext.args ?: ''} \\
  $compress_args $in_args
  ec=\$?

  rename 's/^(classified|unclassified)/${meta.id}.\$1/' {classified,unclassified}* 2>/dev/null && exit \$ec || exit \$ec

  """

  stub:
  def cmd_reads2 = ""
  if (files.size() > 1) {
    cmd_reads2 = "touch ${meta.id}.unclassified_R2.fq ${meta.id}.classified_R2.fq"
  }
  """
  #!/usr/bin/bash

  touch ${meta.id}.output.txt ${meta.id}.report.txt ${meta.id}.classifiedreads.txt
  touch ${meta.id}.unclassified_R1.fq ${meta.id}.classified_R1.fq
  ${cmd_reads2}
  """

}
