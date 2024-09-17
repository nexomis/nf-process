process QUAST {
  container "staphb/quast:5.2.0"

  label 'cpu_med'
  label 'mem_med'
  
  input:
  tuple val(meta), path(assembly, stageAs: "inputs/assembly??.fa")
  tuple val(meta2), path(ref_fa, stageAs: "inputs/reference.fa")
  tuple val(meta3), path(bam, stageAs: "inputs/aln??.bam"), path(bai, stageAs: "inputs/aln??.bai") 

  output:
  tuple val(meta), path("${meta.sample_id ? meta.sample_id : meta.id}/report.html", type: 'file') , emit: html
  tuple val(meta), path("${meta.sample_id ? meta.sample_id : meta.id}/report.tsv", type: 'file')  , emit: tsv
  tuple val(meta), path("${meta.sample_id ? meta.sample_id : meta.id}", type: 'dir')              , emit: dir

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
  def out_dir = meta.sample_id ? meta.sample_id : meta.id
  """
  quast.py \\
    --output-dir $out_dir \\
    --labels ${meta.id} \\
    --threads $task.cpus \\
    $args_ref \\
    ${args_bam} \\
    ${task.ext.args ?: ''} \\
    $assembly
  """

  stub:
  def out_dir = meta.sample_id ? meta.sample_id : meta.id
  """
  #!/usr/bin/bash
  mkdir $out_dir
  touch $out_dir.log $out_dir/report.html $out_dir/report.tsv $out_dir/report.tex $out_dir/icarus.html ${meta.id}/report.pdf 
  """
}
