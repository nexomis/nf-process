process ABACAS {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/abacas:1.3.1--pl5321hdfd78af_3"

  label 'cpu_med'
  label 'mem_8G'

  input:
  tuple val(meta), path(scaffolds, arity: 1, stageAs: 'input_scaffolds.fa')
  tuple val(meta2), path (ref_genome, arity: 1, stageAs: 'input_ref.fa')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}.fasta", type: 'file')

  script:
  // useful to redefine the default mummer_program value (already defined in nextflow.config)?
  """
  #!/usr/bin/bash

  # run ABACAS (TODO, check interest of option: -N  |  -g  |  -P ?  |  -f ?  | -t (tblastx)?)

  abacas.pl -r ${ref_genome} \\
    -q ${scaffolds} \\
    -o ${meta.label ?: meta.id} \\
    ${task.ext.args ?: '-p nucmer'}

  """

  stub:
  """
  #!/usr/bin/bash
  touch ${meta.label ?: meta.id}.fasta
  """
}
