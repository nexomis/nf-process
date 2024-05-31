process MINIPROT {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/miniprot:0.12--he4a0461_0"

  label 'cpu_low'
  label 'mem_12G'

  input:
  tuple val(name_genome), path(genome_fasta, arity: 1)
  tuple val(name_proteome), path(prot_fasta, arity: 1)  

  output:
  tuple val("${name_genome}_${name_proteome}"), path("${name_genome}_${name_proteome}.gff")

  script:
  def name_annot = "${name_genome}_${name_proteome}"
  """
  #!/bin/bash
  miniprot ${task.ext.args ?: ''} -t $task.cpus --gff $genome_fasta $prot_fasta > ${name_annot}.gff.tmp
  grep -v "##PAF" ${name_annot}.gff.tmp > ${name_annot}.gff
  """
}