process CALL_BATCH {
  container "quay.io/nexomis/sav_call:0.2.0-py"
  tag "$meta.id"
  label 'cpu_x1'
  label 'mem_med'

  input:
  tuple val(meta), path(base_files), path(indel_files)  // base_files and indel_files from sav_call outputs grouped by batch
  tuple val(meta_ref), path(ref_fa)
  tuple val(meta_gff), path(gff)

  output:
  tuple val(meta), path("${meta.id}.variants.csv", type: 'file') , optional:false , emit: variants
  tuple val(meta), path("${meta.id}.variants.vcf", type: 'file') , optional:false , emit: vcf
  tuple val(meta), path("${meta.id}.proteins.fa", type: 'file') , optional:false , emit: proteins

  script:
  def base_files_str = base_files instanceof List ? base_files.join(' ') : base_files
  def indel_files_str = indel_files instanceof List ? indel_files.join(' ') : indel_files
  """
  #!/usr/bin/bash

  call.py \\
    --ref ${ref_fa} \\
    --annot ${gff} \\
    --labels ${meta.labels} \\
    --out ${meta.id}.variants.csv \\
    --out_prot ${meta.id}.proteins.fa \\
    --vcf ${meta.id}.variants.vcf \\
    --base_files ${base_files_str} \\
    --indel_files ${indel_files_str} \\
    ${task.ext.call_py_args ?: ''}
  """

  stub:  
  """
  #!/usr/bin/bash

  touch ${meta.id}.variants.csv ${meta.id}.variants.vcf ${meta.id}.proteins.fa
  """ 
}
