process ABACAS {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/abacas:1.3.1--pl5321hdfd78af_3"

  label 'cpu_med'
  label 'mem_8G'

  input:
  tuple val(meta), path(scaffolds, arity: 1, stageAs: 'input_raw/*')
  path (ref_genome)

  output:  // in finally, no need to regroup as much because the name (emit) is not usable in publishDir ...
  tuple val(meta), path("${meta.id}/${meta.id}.fasta", type: 'file')             , optional:false, emit: scaffolds
  tuple val(meta), path("${meta.id}/${meta.id}.*", type: 'file')                 , optional:false, emit: act
  tuple val(meta), path("${meta.id}/unused_contigs.out", type: 'file')           , optional:false, emit: unused_contigs
  tuple val(meta), path("${meta.id}/${task.ext.mummer_program}.*", type: 'file') , optional:false , emit: mumer
  tuple val(meta), path("${meta.id}.log", type: 'file')                          , optional:false, emit: log
  // tuple val(meta), path("reference.notMapped.contigs.tab", type: 'file')      , optional:true  , emit: ref_tab
  // tuple val(meta), path("reference.Repeats.plot", type: 'file')               , optional:true  , emit: ref_plot
  //!!!!! tuple val(meta), path("${meta.id}", type: 'dir')                       , optional:false, emit: output_dir  // usefull ?


  script:
  // useful to redefine the default mummer_program value (already defined in nextflow.config)?
  """
  #!/usr/bin/bash

  # run ABACAS (TODO, check interest of option: -N  |  -g  |  -P ?  |  -f ?  | -t (tblastx)?)
  mkdir ${meta.id}
  cd ${meta.id}/
  abacas.pl -r ../${ref_genome} \\
    -q ../${scaffolds} \\
    -o ${meta.id} \\
    -p ${task.ext.mummer_program ?: 'nucmer'} \\
    ${task.ext.args ?: ''} \\
    2> ../${meta.id}.log
  cd ../
  """

  stub:
  """
  #!/usr/bin/bash
  mkdir ${meta.id}
  touch ${meta.id}.log ${meta.id}/${meta.id}.fasta ${meta.id}/${meta.id}.tab ${meta.id}/nucmer.tiling ${meta.id}/${meta.id}.bin ${meta.id}/${meta.id}.gaps ${meta.id}/nucmer.delta ${meta.id}/${meta.id}.crunch ${meta.id}/${meta.id}.gaps.tab ${meta.id}/nucmer.filtered.delta ${meta.id}/unused_contigs.out
  """
}
