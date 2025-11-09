process EMPTY_FILE {
  container "ubuntu:jammy"
  tag "empty_file"
  cpus 4
  memory 7.GB

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