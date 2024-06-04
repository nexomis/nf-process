process KRAKEN2 {
  container 'ghcr.io/nexomis/kraken2:2.1.3'

  label 'cpu_medium'
  label 'mem_8G'

  input:
  tuple val(sample_name), path(files, arity: 1..2, stageAs: 'input_raw/*')
  path db

  output:
  tuple val(sample_name), path("*.classified*")       , optional:true, emit: classified_reads_fastq
  tuple val(sample_name), path("*.unclassified*")     , optional:true, emit: unclassified_reads_fastq
  tuple val(sample_name), path("*classifiedreads.txt"), optional:true, emit: classified_reads_assignment
  tuple val(sample_name), path("*.report.txt")        , emit: report
  tuple val(sample_name), path("*.output.txt")        , emit: output

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
    rename_paired_inputs += "ln -sf \$PWD/" + files[0] + " \$PWD/input_renamed/" + sample_name + "_1" + fq_ext + "; "
    rename_paired_inputs += "ln -sf \$PWD/" + files[1] + " \$PWD/input_renamed/" + sample_name + "_2" + fq_ext
    in_args = "--paired input_renamed/" + sample_name + "_1" + fq_ext + " input_renamed/" + sample_name + "_2" + fq_ext
  } else {
    in_args += files[0]
  }
  
  
  
  """
  #!/usr/bin/bash

  $rename_paired_inputs

  kraken2 --thread $task.cpus \\
  --output ${sample_name}.output.txt --report ${sample_name}.report.txt \\
  --db $db \\
  ${task.ext.args ?: ''} \\
  $compress_args $in_args
  ec=\$?

  rename 's/^(classified|unclassified)/${sample_name}.\$1/' {classified,unclassified}* 2>/dev/null && exit \$ec || exit \$ec

  """

  stub:
  """
  #!/usr/bin/bash

  touch ${sample_name}.output.txt
  touch ${sample_name}.report.txt

  """

}
