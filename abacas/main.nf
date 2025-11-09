process ABACAS {
  container "quay.io/biocontainers/abacas:1.3.1--pl5321hdfd78af_3"
  tag "$meta.id-$meta2.id"
  cpus 8
  memory 8.GB

  input:
  tuple val(meta), path(scaffolds, arity: 1, stageAs: 'input_scaffolds.fa')
  tuple val(meta2), path (ref_genome, arity: 1, stageAs: 'input_ref.fa')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}.fasta", type: 'file')

  script:
  // useful to redefine the default mummer_program value (already defined in nextflow.config)?
  """
  #!/usr/bin/bash

  awk '/^>/{filename="input_ref_" ++i ".fa"} {print \$0 > filename}' input_ref.fa

  for ref_file in input_ref_*.fa; do
    index=\${ref_file#input_ref_}
    index=\${index%.fa}
    output_prefix="${meta.label ?: meta.id}_\$index"

    abacas.pl -r \$ref_file \\
      -q ${scaffolds} \\
      -o "\$output_prefix" \\
      ${task.ext.args ?: '-p nucmer'}

    output_fa="\${output_prefix}.fasta"

    if [ -f "\$output_fa" ]; then
      char_count=\$(head -n 100 "\$output_fa" | grep -v '^>' | tr -d '\n' | wc -c)
      if [ "\$char_count" -le 100 ]; then
        rm "\$output_fa"
      fi
    fi
  done

  touch ${meta.label ?: meta.id}_0.fasta
  cat ${meta.label ?: meta.id}_*.fasta > ${meta.label ?: meta.id}.fasta

  """

  stub:
  """
  #!/usr/bin/bash
  touch ${meta.label ?: meta.id}.fasta
  """
}
