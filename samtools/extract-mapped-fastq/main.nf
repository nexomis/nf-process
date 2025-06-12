process EXTRACT_MAPPED_FASTQ {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/samtools:1.20--h50ea8bc_1"

  label 'cpu_high'
  label 'mem_2G_per_cpu'

  input:
  tuple val(meta), path(bam, arity: 1, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}.mapped_*.fq.gz"), optional:false , emit: mapped
  tuple val(meta), path("${meta.label ?: meta.id}.unmapped_*.fq.gz"), optional:false , emit: unmapped

  script:

  """
  #!/usr/bin/bash

  samtools sort -n $bam -o collate.bam -O bam -@ $task.cpus
  samtools view -@ $task.cpus -b -F 12 -o mapped.bam collate.bam
  samtools view -@ $task.cpus -b -f 4 -F 8 -o first.bam collate.bam
  samtools view -@ $task.cpus -b -f 8 -F 4 -o second.bam collate.bam
  samtools merge -@ $task.cpus -n -o merged.bam mapped.bam first.bam second.bam

  samtools fastq -@ $task.cpus merged.bam \\
    -s ${meta.label ?: meta.id}.mapped_R1.fq.gz
  if [ "\$(zcat ${meta.label ?: meta.id}.mapped_R1.fq.gz 2>/dev/null | wc -l)" -lt 4 ]; then
    samtools fastq -@ $task.cpus merged.bam \\
      -1 ${meta.label ?: meta.id}.mapped_R1.fq.gz \\
      -2 ${meta.label ?: meta.id}.mapped_R2.fq.gz
  fi

  samtools fastq -f 12 -@ $task.cpus merged.bam \\
    -s ${meta.label ?: meta.id}.unmapped_R1.fq.gz
  if [ "\$(zcat ${meta.label ?: meta.id}.unmapped_R1.fq.gz 2>/dev/null | wc -l)" -lt 4 ]; then
    samtools fastq -f 12 -@ $task.cpus merged.bam \\
      -1 ${meta.label ?: meta.id}.unmapped_R1.fq.gz \\
      -2 ${meta.label ?: meta.id}.unmapped_R2.fq.gz
  fi

  """

  stub:
  """
  #!/usr/bin/bash

  touch \\
    ${meta.label ?: meta.id}.mapped_R1.fq.gz \\
    ${meta.label ?: meta.id}.mapped_R2.fq.gz \\
    ${meta.label ?: meta.id}.unmapped_R1.fq.gz \\
    ${meta.label ?: meta.id}.unmapped_R2.fq.gz
  
  """
}
