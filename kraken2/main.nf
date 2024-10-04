process KRAKEN2 {
  container 'ghcr.io/nexomis/kraken2:2.1.3'

  label 'cpu_med'
  label 'mem_4G_per_cpu' // params.kraken2_memory ?

  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  tuple val(meta2), path(db, stageAs: "inputs/index")

  output:
  tuple val(meta), path("*classified*")        , optional:true, emit: classified_reads_fastq
  tuple val(meta), path("*unclassified*")      , optional:true, emit: unclassified_reads_fastq
  tuple val(meta), path("*classifiedreads.txt"), optional:true, emit: classified_reads_assignment
  tuple val(meta), path("*.txt")               , emit: report
  tuple val(meta), path("*.gz")            , emit: output

  // modify ext from config to enable classified and unclassified
  // ext:
  // args: --classified-out classified.#.fq --unclassified-out unclassified.#.fq

  script:

  """
  #!/usr/bin/bash

  kraken2 --thread $task.cpus \\
  --report ${meta.id}.report.txt \\
  --db $db \\
  ${(reads.size() === 2)? '--paired': ''} \\
  ${task.ext.args ?: ''} \\
  $reads \\
  | gzip > ${meta.id}.gz
  
  ec=\$?

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
