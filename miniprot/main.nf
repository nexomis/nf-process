process MINIPROT {
  container "quay.io/biocontainers/miniprot:0.13--he4a0461_0" // does not exists in aws
  tag "$meta.id"
  cpus 4
  memory 7.GB

  input:
  tuple val(meta), path(genome_fasta, arity: 1, stageAs: "inputs/genome.fa")
  tuple val(meta2), path(prot_fasta, arity: 1, stageAs: "inputs/proteome.fa")  

  output:
  tuple val(meta), path("${meta.label ?: meta.id}.gff")

  script:
  """
  #!/bin/bash
  miniprot ${task.ext.args ?: ''} -t $task.cpus --gff $genome_fasta $prot_fasta > ${meta.label ?: meta.id}.gff.tmp
  grep -v "##PAF" ${meta.label ?: meta.id}.gff.tmp > ${meta.label ?: meta.id}.gff
  """
}
