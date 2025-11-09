// assembly using 'spades' with short reads : 1 input fastq[.gz] file or 2 in case of paird-end.
// (even if compatible with the 'spades' tool, cases with several independent input files are not handled here, this inclued additional librarie or hybrid approach).

process SPADES {
  container "quay.io/nexomis/spades:4.0.0-91e677"
  tag "$meta.id"
  cpus 16
  memory 30.GB

  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  
  output:
  tuple val(meta), path("${meta.label ?: meta.id}/*scaffolds.fasta", type: 'file')  , emit: scaffolds
  tuple val(meta), path("${meta.label ?: meta.id}/*contigs.fasta", type: 'file')    , emit: contigs

  script:
  def out_dir = meta.label ?: meta.id
  """
  #!/usr/bin/bash

  set -e
  
  spades.py --threads ${task.cpus} \\
    --memory ${task.memory.toGiga()} \\
    -o ${out_dir} \\
    ${ (reads.size() == 1) ? "-s ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}" } \\
    ${task.ext.args ?: ''} ${meta.args_spades ?: ''}

  if [ -f "${out_dir}/scaffolds.fasta" ] && [ -f "${out_dir}/raw_scaffolds.fasta" ]; then
    rm "${out_dir}/raw_scaffolds.fasta"
  fi

  if [ -f "${out_dir}/contigs.fasta" ] && [ -f "${out_dir}/raw_contigs.fasta" ]; then
    rm "${out_dir}/raw_contigs.fasta"
  fi

  if [ ! -f "${out_dir}/raw_scaffolds.fasta" ] && [ ! -f "${out_dir}/scaffolds.fasta" ]; then
    if [ -f "${out_dir}/contigs.fasta" ]; then
      cp "${out_dir}/contigs.fasta" "${out_dir}/scaffolds.fasta"
    elif [ -f "${out_dir}/raw_contigs.fasta" ]; then
      cp "${out_dir}/raw_contigs.fasta" "${out_dir}/scaffolds.fasta"
    else
      touch ${out_dir}/contigs.fasta
      touch ${out_dir}/scaffolds.fasta
    fi
  fi

  """

  stub:
  def out_dir = meta.label ?: meta.id
  """
  #!/usr/bin/bash
  mkdir ${out_dir}
  touch ${out_dir}.log ${out_dir}/scaffolds.fasta ${out_dir}/contigs.fasta ${out_dir}/assembly_graph_with_scaffolds.gfa
  """
}
