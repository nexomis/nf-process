// assembly using 'spades' with short reads : 1 input fastq[.gz] file or 2 in case of paird-end.
// (even if compatible with the 'spades' tool, cases with several independent input files are not handled here, this inclued additional librarie or hybrid approach).

process SPADES {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/spades:4.0.0--h5fb382e_2"

  label 'cpu_high'
  label 'mem_high'

  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  
  output:
  tuple val(meta), path("${meta.id}/scaffolds.fasta", type: 'file')                   , optional:true , emit: scaffolds
  tuple val(meta), path("${meta.id}/contigs.fasta", type: 'file')                     , optional:true , emit: contigs
  tuple val(meta), path("${meta.id}/assembly_graph_with_scaffolds.gfa", type: 'file') , optional:true , emit: gfa
  tuple val(meta), path("${meta.id}/spades.log", type: 'file')                        , optional:false, emit: log_tools
  tuple val(meta), path("${meta.id}.log", type: 'file')                               , optional:false, emit: log
  //!!! tuple val(meta), path("${meta.id}", type: 'dir')                              , optional:false, emit: output_dir  // usefull ?

  script:
  def memory_in_gbit = MemoryUnit.of("${task.memory}").toUnit('GB') * 8
  def memory_in_gbit_min1 = Math.max(memory_in_gbit, 1)

  """
  #!/usr/bin/bash
  
  # run SPAdes : spadesMode (rnaviral) on meta.spades_args (not on task.ext.args) ?
  spades.py --threads ${task.cpus} \\
    --memory ${memory_in_gbit_min1} \\
    -o ${meta.id} \\
    ${ (reads.size() == 1) ? "-s ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}" } \\
    ${meta.args_spades ?: ''} \\
    ${task.ext.args ?: ''} \\
    2> ${meta.id}.log

  # Select output files :
  mkdir ${meta.id}
  if [ ! -f ${meta.id}/scaffolds.fasta ]; then
    echo "WARNING : scaffolds file not generated: use contigs of raw_contigs instead ?? (and contigs.gfa instead of scaffolds.gfa)"
  fi
  """

  stub:
  def memory_in_gbit = MemoryUnit.of("${task.memory}").toUnit('GB') * 8
  def memory_in_gbit_min1 = Math.max(memory_in_gbit, 1)
  
  """
  #!/usr/bin/bash
  mkdir ${meta.id}
  touch ${meta.id}.log ${meta.id}/scaffolds.fasta ${meta.id}/contigs.fasta ${meta.id}/assembly_graph_with_scaffolds.gfa ${meta.id}/spades.log
  """
}
