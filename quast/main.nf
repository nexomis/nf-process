process QUAST {
  container "staphb/quast:5.2.0"

  label 'cpu_med'
  label 'mem_med'
  
  input:
  tuple val(meta), path(assembly, stageAs: "inputs/assembly??.fa")
  tuple val(meta2), path(ref_fa, stageAs: "inputs/reference.fa")
  tuple val(meta3), path(bam, stageAs: "inputs/aln??.bam"), path(bai, stageAs: "inputs/aln??.bai") 

  output:
  tuple val(meta), path("${meta.id}/report.html", type: 'file') , emit: html
  tuple val(meta), path("${meta.id}/report.tsv", type: 'file')  , emit: tsv
  tuple val(meta), path("${meta.id}", type: 'dir')              , emit: dir

  script:
  def args_ref = ref_fa.size() > 1 ? "-r $ref_fa" : ""
  def args_bam = ''
  def bamList = []
  def bam_not_empty = false
  if (bam instanceof Path){
    bamList << bam.toString()
    if (bam.size() > 1) {bam_not_empty = true}
  } else {
    bam.each{item -> 
      bamList << item.toString()
      if (item.size() > 1) {bam_not_empty = true}
    }
  }
  if (bam_not_empty) {args_bam = "--bam " + bamList.join(',')}
  
  """
  quast.py \\
    --output-dir ${meta.id} \\
    --labels ${meta.id} \\
    --threads $task.cpus \\
    $args_ref \\
    ${args_bam} \\
    ${task.ext.args ?: ''} \\
    $assembly
  """

  stub:
  """
  #!/usr/bin/bash
  mkdir ${meta.id}
  touch ${meta.id}.log ${meta.id}/report.html ${meta.id}/report.tsv ${meta.id}/report.tex ${meta.id}/icarus.html ${meta.id}/report.pdf 
  """
}
