process CONCAT_FQ {
  container 'ubuntu:noble-20241011'
  tag "$meta.id"
  cpus 1
  memory 2.GB

  input:
  tuple val(meta), path(fastq1, arity: 1..2, stageAs: 'input1/*')
  tuple val(meta2), path(fastq2, arity: 1..2, stageAs: 'input2/*')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}*.fq.gz", arity: 1..2)

  script:
  def cmd2 = ""
  if (fastq1.size() > 1) {
    cmd2 = "cat ${fastq1[1]} ${fastq2[1]} > ${meta.label ?: meta.id}_R2.fq.gz"
  }

  """
  cat ${fastq1[0]} ${fastq2[0]} > ${meta.label ?: meta.id}_R1.fq.gz
  $cmd2
  """ 

  stub:
  """
  touch ${meta.label ?: meta.id}_R1.fq.gz
  """

}