process IVAR_VARIANTS_ALL {
  container "quay.io/nexomis/ivar:1.4.3"

  label 'cpu_x1'
  label 'mem_med'

  input:
  tuple val(meta), path(bam, arity: 1..2, stageAs: 'input/')        // bam[0] is bam file [recquired] and bam[1] is its index (bai) [by default optional, but it may be necessary (and therefore required) depending on particular samtools mpileup options ( in particular region-specific options such as '-r')]
  tuple val(meta2), path(genome_ref, arity: 2, stageAs: 'input/')   // genome_ref: .fa and .gff/.gtf

  output:
  tuple val(meta), path("${meta.label ?: meta.id}_raw.tsv", type: 'file') , optional:false , emit: raw_iSNV_tsv
  tuple val(meta), path("${meta.label ?: meta.id}_cols_1_4.mpileup", type: 'file') , optional:false , emit: mpileup_cov

  script:
  genome_ref_fa = genome_ref[0]
  genome_ref_gff = genome_ref[1]

  """
  #!/usr/bin/bash

  samtools mpileup \\
    -aa -A -B -Q 0 \\
    --reference ${genome_ref_fa} \\
    ${task.ext.args ?: ''} \\
    ${bam[0]} |\\
    loose_ends |\\
    tee >(cut -f1-4 > ${meta.label ?: meta.id}_cols_1_4.mpileup) |\\
    ivar variants \\
      -p ${meta.label ?: meta.id}_raw \\
      -q 0 -t 0 -m 0 \\
      -r ${genome_ref_fa} \\
      -g ${genome_ref_gff} \\
       ${task.ext.args2 ?: ''}  

  """

  stub:  
  """
  #!/usr/bin/bash

  touch ${meta.label ?: meta.id}_raw.tsv ${meta.label ?: meta.id}_col1-4.mpileup
  """ 
}