process SLIMFASTQ_COMPRESS {
  container 'quay.io/biocontainers/slimfastq:2.04--h503566f_5'
  tag "${meta.id}"
  label 'cpu_x2'
  label 'mem_4G'

  input:
  tuple val(meta), path(files, arity: 1..2, stageAs: "inputs/*")

  output:
  tuple val(meta), path("*.sfq", arity: 1..2)

  script:
  meta.read_type = "sfq"
  def R1 = files[0]
  def R2 = files.size() > 1 ? files[1] : null
  def slimfastq = "slimfastq ${task.ext.args ?: '-l 3'}"
  """
  #!/usr/bin/bash
  set -e
  
  # Generate MD5 for input files
  if [[ "${R1}" == *.gz ]] || [[ "${R1}" == *.gzip ]]; then
    gzip -dc ${R1} | md5sum | cut -f 1 -d " " > initial.md5
    ${R2 ? "gzip -dc ${R2} | md5sum | cut -f 1 -d \" \" >> initial.md5" : ""}
    gzip -dc ${R1} | ${slimfastq} ${meta.id}_R1.sfq
    ${R2 ? "gzip -dc ${R2} | ${slimfastq} ${meta.id}_R2.sfq" : ""}
  else
    md5sum ${R1} | cut -f 1 -d " " > initial.md5
    ${R2 ? "md5sum ${R2} | cut -f 1 -d \" \" >> initial.md5" : ""}
    ${slimfastq} ${R1} ${meta.id}_R1.sfq
    ${R2 ? "${slimfastq} ${R2} ${meta.id}_R2.sfq" : ""}
  fi
  
  slimfastq ${meta.id}_R1.sfq | md5sum | cut -f 1 -d " " > decompressed.md5
  ${R2 ? "slimfastq ${meta.id}_R2.sfq | md5sum | cut -f 1 -d \" \" >> decompressed.md5" : ""}
  
  # Check if MD5s match
  if ! diff initial.md5 decompressed.md5; then
    echo "MD5s do not match"
    exit 1
  fi
  """

  stub:
  """
  #!/usr/bin/bash
  touch ${meta.id}_R1.sfq ${files.size() > 1 ? meta.id + '_R2.sfq' : ''}
  """
}
