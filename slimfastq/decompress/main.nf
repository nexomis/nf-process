process SLIMFASTQ_DECOMPRESS {
  container 'quay.io/biocontainers/slimfastq:2.04--h503566f_5'
  tag "${meta.id}"
  label 'cpu_x2'
  label 'mem_4G'

  input:
  tuple val(meta), path(files, arity: 1..2)

  output:
  tuple val(meta), path("*.fq.gz", arity: 1..2)

  script:
  def R1_sfq = files[0]
  def R2_sfq = files.size() > 1 ? files[1] : null
  def R1_out = "${meta.id}_R1.fq.gz"
  def R2_out = "${meta.id}_R2.fq.gz"
  meta.remove("read_type")
  """
  #!/usr/bin/bash

  set -e

  # Decompress and compress to gzip
  slimfastq ${R1_sfq} | gzip > ${R1_out}
  ${R2_sfq ? "slimfastq ${R2_sfq} | gzip > ${R2_out}" : ""}
  """

  stub:
  """
  #!/usr/bin/bash
  touch ${meta.id}_R1.fq.gz ${files.size() > 1 ? meta.id + '_R2.fq.gz' : ''}
  """
}
