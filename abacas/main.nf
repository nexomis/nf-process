process ABACAS {
  container "${params.biocontainers_registry ?: 'quay.io/biocontainers/abacas:1.3.1--pl5321hdfd78af_3'}"

  label 'cpu_med'
  label 'mem_8G'

  input:
  tuple val(meta), path(scaffolds, arity: 1, stageAs: 'input_raw/*')
  path (ref_genome)
  val(MUMmer_program)
  val keep_on_output  //  = 'scaffold'

  output:
  tuple val(meta), path("$meta.id", type: 'dir'), optional:false, emit: output_dir
  tuple val(meta), path("${meta.id}.log", type: 'file'), optional:false, emit: log

  script:
  """
  #!/usr/bin/bash

  # run ABACAS (TODO, check interest of option: -N  |  -g  |  -P ?  |  -f ?  ) : include 'MUMmer_program' and 'keep_on_output' on task.ext.args (not on meta.abacas_args) ?
  abacas.pl -r ${ref_genome} \\
    -q ${scaffolds} \\
    -o ${meta.id} \\
    -p ${MUMmer_program}
    ${meta.abacas_args ?: ''} \\
    ${task.ext.args ?: ''} \\
    2> ${meta.id}.log

  # manage output :
  mkdir ${meta.id}/
  mv ${meta.id}.fasta ${meta.id}/${meta.id}.fasta

  case "${keep_on_output}" in
    "all")
      mv unused_contigs.out ${MUMmer_program}.* ${meta.id}.* ${meta.id}/
      ;;
    "include_extended_act_files")
      mv ${meta.id}.* ${meta.id}/
      ;;
    "include_act_files")
      mv ${meta.id}.crunch ${meta.id}.fasta ${meta.id}/
      ;;
    "scaffold")
      :
      ;;
    *)
      echo "No valid include option provided ('keep_on_output': 'scaffold', 'include_act_files', 'include_extended_act_files' or 'all'): defaulting to 'scaffold'."
      ;;
  esac
  """

  stub:
  """
  #!/usr/bin/bash

  mkdir ${meta.id}
  touch ${meta.id}.log ${meta.id}/${meta.id}.fasta 
  """
}

