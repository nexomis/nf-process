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

  # Sort BAM by read name for proper pairing during extraction
  samtools sort -n $bam -o collate.bam -O bam -@ $task.cpus
  
  # SAM FLAGS: 4=read unmapped, 8=mate unmapped, 12=both unmapped (4+8)
  # Extract reads where BOTH mates are mapped (-F 12: exclude flag 4+8, i.e., exclude both unmapped)
  samtools view -@ $task.cpus -b -F 12 -o mapped.bam collate.bam
  samtools view -@ $task.cpus -b -f 12 -o unmapped.bam collate.bam
  
  # Extract reads where THIS read is unmapped BUT mate is mapped 
  # (-f 4: include read unmapped, -F 8: exclude mate unmapped)
  # CONTAINS: Both R1 records (where R1 unmapped, R2 mapped) AND R2 records (where R2 unmapped, R1 mapped)
  samtools view -@ $task.cpus -b -f 4 -F 8 -o first.bam collate.bam
  
  # Extract reads where THIS read is mapped BUT mate is unmapped
  # (-f 8: include mate unmapped, -F 4: exclude read unmapped)
  # CONTAINS: Both R1 records (where R1 mapped, R2 unmapped) AND R2 records (where R2 mapped, R1 unmapped)
  samtools view -@ $task.cpus -b -f 8 -F 4 -o second.bam collate.bam
  
  # Merge all reads that have at least one mate mapped (properly paired + singletons)
  # -n flag indicates input files are already sorted by read name (from collate.bam)
  # This captures: both mapped + read unmapped/mate mapped + read mapped/mate unmapped
  samtools merge -@ $task.cpus -n -o merged.bam mapped.bam first.bam second.bam
  
  # Sort merged BAM by name for proper read pairing during FASTQ extraction
  samtools sort -n merged.bam -o merged_sorted.bam -O bam -@ $task.cpus

  # Check if data is paired-end by looking for "0 + 0 paired in sequencing" in flagstats
  # If found, use -s (single-end), otherwise use -1/-2 (paired-end)
  samtools view -h collate.bam | head -n 100000 | samtools view -h -b - | samtools flagstats - > flagstats.txt
  if cat flagstats.txt | grep -rq "^0 + 0 paired in sequencing"; then
    # Single-end data: use -s flag
    samtools fastq -@ $task.cpus merged_sorted.bam \\
      -s ${meta.label ?: meta.id}.mapped_R1.fq.gz
    samtools fastq -@ $task.cpus unmapped.bam \\
      -s ${meta.label ?: meta.id}.unmapped_R1.fq.gz
  else
    # Paired-end data: use -1/-2 flags
    samtools fastq -@ $task.cpus merged_sorted.bam \\
      -1 ${meta.label ?: meta.id}.mapped_R1.fq.gz \\
      -2 ${meta.label ?: meta.id}.mapped_R2.fq.gz
    samtools fastq -@ $task.cpus unmapped.bam \\
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
