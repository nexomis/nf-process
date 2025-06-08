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
    -0 /dev/null -s /dev/null \\
    -1 ${meta.label ?: meta.id}.mapped_R1.fq.gz \\
    -2 ${meta.label ?: meta.id}.mapped_R2.fq.gz

  samtools fastq -@ $task.cpus -f 12 collate.bam \\
    -0 /dev/null -s /dev/null \\
    -1 ${meta.label ?: meta.id}.unmapped_R1.fq.gz \\
    -2 ${meta.label ?: meta.id}.unmapped_R2.fq.gz

  for file in *_R2.fq.gz; do
    [ "\$(zcat \\"\$file\\" 2>/dev/null | wc -l)" -lt 4 ] && rm "\$file"
  done

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