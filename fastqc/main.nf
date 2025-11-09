// inpired from https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqc/main.nf

process FASTQC {
  container 'staphb/fastqc:0.12.1'
  tag "$meta.id"
  cpus 4
  memory 12.GB

  input:
  tuple val(meta), path(files, arity: 1..2)

  output:
  tuple val(meta), path("*.html"), emit: html
  tuple val(meta), path("*.zip") , emit: zip

  script:
  def fastqc_memory = Math.min(Math.max(100, task.memory.toMega() - 512), 10000)
  
  """
  #!/usr/bin/bash

  fastqc --nogroup --threads $task.cpus --memory $fastqc_memory $files

  """

  stub:
  def zip_stub = meta.id
  def html_stub = meta.id
  if (files.size() > 1) {
    zip_stub += "_R1.zip " + meta.id + "_R2.zip"
    html_stub += "_R1.html " + meta.id + "_R2.html"
  } else {
    zip_stub += ".zip"
    html_stub += ".html"
  }
  """
  #!/usr/bin/bash

  touch $html_stub $zip_stub

  """

}
