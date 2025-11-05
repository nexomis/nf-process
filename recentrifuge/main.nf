process RECENTRIFUGE {
  container 'quay.io/biocontainers/recentrifuge:1.14.1--pyhdfd78af_0'
  label 'cpu_x1'
  label 'mem_low'

  input:
  path(files)
  path(tax_dir, stageAs: 'raw_taxdump')

  output:
  path("report.rcf.html"), emit: html
  path("report.rcf.xlsx"), emit: xlsx

  script:
  def args_files = []
  if (files instanceof List) {
    files.each{ it ->
      args_files << it.toString()
    }
  } else {
    args_files << files.toString()
  }
  """
  #!/usr/bin/bash

  SOURCE_DIR=raw_taxdump
  DEST_DIR=taxdump

  mkdir -p \$DEST_DIR
  DEST_ABS_PATH=\$(realpath "\$DEST_DIR")
  for file in \$SOURCE_DIR/*; do
    if [ -f "\$file" ]; then
      filename=\$(basename "\$file")
        if echo "\$filename" | grep -q '\\.gz\$'; then
          echo "Unzipping \$filename to \$DEST_DIR"
          gunzip -c "\$file" > "\$DEST_DIR/\${filename%.gz}"
        else
          echo "Creating symlink for \$filename in \$DEST_DIR"
          ln -sf "\$(realpath "\$file")" "\$DEST_ABS_PATH/\$filename"
        fi
    fi
  done

  rcf -k ${args_files.join(" -k ")} -o report --sequential

  """

  stub:
  """
  #!/usr/bin/bash

  touch report.rcf.html
  touch report.rcf.xlsx

  """

}