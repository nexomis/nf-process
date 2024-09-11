process MINIPROT {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/miniprot:0.12--he4a0461_0"

  label 'cpu_low'
  label 'mem_low'

  input:
  tuple val(meta), path(genome_fasta, arity: 1, stageAs: "inputs/genome.fa")
  tuple val(meta2), path(prot_fasta, arity: 1, stageAs: "inputs/proteome.fa")  

  output:
  tuple val("${meta.id}_${meta2.id}"), path("${meta.id}_${meta2.id}.gff")

  script:
  def name_annot = "${meta.id}_${meta2.id}"
  """
  #!/bin/bash
  miniprot ${task.ext.args ?: ''} -t $task.cpus --gff $genome_fasta $prot_fasta > ${name_annot}.gff.tmp
  grep -v "##PAF" ${name_annot}.gff.tmp > ${name_annot}.gff
  """
}