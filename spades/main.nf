// assembly using 'spades' with short reads : 1 input fastq[.gz] file or 2 in case of paird-end.
// (even if compatible with the 'spades' tool, cases with several independent input files are not handled here, this inclued additional librarie or hybrid approach).

process SPADES {
  container "${params.biocontainers_registry ?: 'quay.io/biocontainers/spades:4.0.0--h5fb382e_2'}"

  label 'cpu_high'
  label 'mem_high'

  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  
  output:
  tuple val(meta), path("${meta.id}", type: 'dir')                                      , optional:false, emit: output_dir
  tuple val(meta), path("${meta.id}/scaffolds.fasta", type: 'file')                     , optional:true , emit: scaffolds
  tuple val(meta), path("${meta.id}/contigs.fasta", type: 'file')                       , optional:true , emit: contigs
  tuple val(meta), path("${meta.id}/transcripts.fasta", type: 'file')                   , optional:true , emit: transcripts
  tuple val(meta), path("${meta.id}/gene_clusters.fasta", type: 'file')                 , optional:true , emit: gene_clusters
  tuple val(meta), path("${meta.id}/assembly_graph_with_scaffolds.fasta", type: 'file') , optional:true , emit: gfa
  tuple val(meta), path("${meta.id}/spades.log"), type: 'file'                          , optional:false, emit: log_tools
  tuple val(meta), path("${meta.id}.log", type: 'file')                                 , optional:false, emit: log


  script:
  """
  #!/usr/bin/bash
  
  mem_int=\$(echo '${task.memory}' | sed "s/[^0-9. ]*\$//")  # task.memory must be in Gb .... | TODO: manage other unit !
  # run SPAdes : spadesMode (rnaviral) on meta.spades_args (not on task.ext.args) ?
  spades.py --threads ${task.cpus} \\
    --memory \${mem_int} \\
    -o ${meta.id} \\
    ${ (reads.size() == 1) ? "-s ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}" } \\
    ${meta.spades_args ?: ''} \\
    ${task.ext.args ?: ''} \\
    2> ${meta.id}.log

  # Select output files :
  mkdir ${meta.id}
  if [ ! -f ${meta.id}/scaffolds.fasta ]; then   # else - TODO: print big Warning Message ?!
    echo "WARNING : scaffolds file not generated,: use contigs of raw_contigs instead ?? (and contigs.gfa instead of scaffolds.gfa)
  fi
  """

  stub:
  """
  #!/usr/bin/bash

  mkdir ${meta.id}
  touch ${meta.id}.log ${meta.id}/${meta.id}.scaffolds.fa ${meta.id}/${meta.id}.contigs.fa ${meta.id}/${meta.id}.scaffolds.gfa
  """
}

