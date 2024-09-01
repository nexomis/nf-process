process EMPTY_FILE {
  container "ubuntu:jammy"

  label 'cpu_low'
  label 'mem_low'

  output:
  path "EMPTY_FILE"

  script:
  """
  touch EMPTY_FILE
  """

  stub:
  """
  touch EMPTY_FILE
  """
}