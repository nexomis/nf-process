process ABACAS {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/abacas:1.3.1--pl5321hdfd78af_3"

  label 'cpu_med'
  label 'mem_8G'

  input:
  tuple val(meta), path(scaffolds, arity: 1, stageAs: 'input_raw/*')
  path (ref_genome)

  output:
  tuple val(meta), path("${meta.id}", type: 'dir')                          , optional:false, emit: output_dir
  tuple val(meta), path("${meta.id}/${meta.id}.fasta", type: 'file')        , optional:false, emit: scaffold
  tuple val(meta), path("${meta.id}/${meta.id}.*", type: 'file')            , optional:false, emit: act
  tuple val(meta), path("unused_contigs.out", type: 'file')                 , optional:false, emit: unused_contigs
  tuple val(meta), path("${meta.id}.log", type: 'file')                     , optional:false, emit: log
  // tuple val(meta), path("${task.ext.mummer_program}.*", type: 'file')    , optional:false , emit: mumer
  // tuple val(meta), path("reference.notMapped.contigs.tab", type: 'file') , optional:true  , emit: ref_tab
  // tuple val(meta), path("reference.Repeats.plot", type: 'file')          , optional:true  , emit: ref_plot


  script:
  """
  #!/usr/bin/bash

  # run ABACAS (TODO, check interest of option: -N  |  -g  |  -P ?  |  -f ?  ) : include 'MUMmer_program' and 'keep_on_output' on task.ext.args (not on meta.abacas_args) ?
  abacas.pl -r ${ref_genome} \\
    -q ${scaffolds} \\
    -o ${meta.id} \\
    -p ${task.ext.mummer_program ?: 'nucmer'} \\
    ${meta.abacas_args ?: ''} \\
    2> ${meta.id}.log
  
  mkdir ${meta.id}
  cp ${meta.id}* ${meta.id}/
  rm ${meta.id}/${meta.id}.log
  """

  stub:
  """
  #!/usr/bin/bash
  mkdir ${meta.id}
  touch ${meta.id}/${meta.id}.fasta ${meta.id}/${meta.id}.tab ${meta.id}/${meta.id}.bin ${meta.id}/${meta.id}.crunch ${meta.id}/${meta.id}.gaps
  touch ${meta.id}.log ${meta.id}.fasta ${meta.id}.tab ${meta.id}.bin ${meta.id}.crunch ${meta.id}.gaps  unused_contigs.out ${task.ext.mummer_program}.delta ${task.ext.mummer_program}.filtered.delta  ${task.ext.mummer_program}.tiling
  """
}
