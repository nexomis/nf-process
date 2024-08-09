process ABACAS {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/abacas:1.3.1--pl5321hdfd78af_3"

  label 'cpu_med'
  label 'mem_8G'

  input:
  tuple val(meta), path(scaffolds, arity: 1, stageAs: 'input_raw/*')
  path ref_fasta

  output:
  tuple val(meta), path("$meta.id", type: 'dir'), optional:false, emit: output_dir
  tuple val(meta), path("${meta.id}.log", type: 'file'), optional:false, emit: log


  script:
  """
  #!/usr/bin/bash

  # run ABACAS
  abacas.pl -r ${ref_fasta} \\
    -q ${scaffolds} \\
    -o ${meta.id} \\
    ${meta.abacas_args ?: ''} \\
    ${task.ext.args ?: ''} \\
    2> ${meta.id}.log

// intetest option : -p nucmer/promer ?  |  -N  |  -g  |  -P ?  |  -f ?  
  # Move and compress files
  # keep only interest file
  
  """

  stub:
  """
  #!/usr/bin/bash

  mkdir ${meta.id}
  touch ${meta.id}.log xxx
  # add output file
  """
}

