// inpired from https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqc/main.nf

process FASTQC {
  container 'staphb/fastqc:0.12.1'

  label 'cpu_low'
  label 'mem_8G'

  input:
  tuple val(sample_name), path(files, arity: 1..2)

  output:
  tuple val(sample_name), path("*.html"), emit: html
  tuple val(sample_name), path("*.zip") , emit: zip

  script:
  def memory_in_mb = MemoryUnit.of("${task.memory}").toUnit('MB')
  def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)
  
  """
  #!/usr/bin/bash

  fastqc --nogroup --threads $task.cpus --memory $fastqc_memory $files

  """

  stub:
  def zip_stub = sample_name
  def html_stub = sample_name
  if (files.size() > 1) {
    zip_stub += "_R1.zip " + sample_name + "_R2.zip"
    html_stub += "_R1.html " + sample_name + "_R2.html"
  } else {
    zip_stub += ".zip"
    html_stub += ".html"
  }
  """
  #!/usr/bin/bash

  touch $html_stub $zip_stub

  """

}
