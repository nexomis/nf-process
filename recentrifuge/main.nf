process RECENTRIFUGE {
  container 'quay.io/biocontainers/recentrifuge:1.14.1--pyhdfd78af_0'

  label 'cpu_x1'
  label 'mem_low'

  input:
  path(files)
  path(tax_dir, stageAs: 'taxdump')

  output:
  path("report.rcf.html"), emit: html
  path("report.rcf.xlsx"), emit: xlsx

  script:
  def args_files = []
  if (files instanceof List) {
    files.each{
      args_files << it.toString()
    }
  } else {
    args_files << files.toString()
  }
  """
  #!/usr/bin/bash

  rcf -k ${args_files.join(" -k ")} -o report --sequential

  """

  stub:
  """
  #!/usr/bin/bash

  touch report.rcf.html
  touch report.rcf.xlsx

  """

}