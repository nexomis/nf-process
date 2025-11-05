process EMPTY_FILE {
  container "ubuntu:jammy"
  tag "empty_file"
  label 'cpu_low'
  label 'mem_low'

  output:
  tuple val({[id:"empty_file"]}), path("EMPTY_FILE")

  script:
  """
  touch EMPTY_FILE
  """

  stub:
  """
  touch EMPTY_FILE
  """
}