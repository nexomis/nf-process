// need to be valided in case of multifasta reference (segmenteed virus, contigs, ...)

process IVAR_VARIANTS_ALL_POS {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/ivar:1.4.3--h43eeafb_0"

  label 'cpu_x1'
  label 'mem_med'

  input:
  tuple val(meta), path(bam, arity: 1..2, stageAs: 'input/')   // bam[0] is bam file [recquired] and bam[1] is its index (bai) [optional but save time if its provided]
  tuple val(meta2), path(genome_ref, arity: 2, stageAs: 'input/')
  //tuple val(meta2), path(genome_ref_fa, arity: 1..2, stageAs: 'input/'), path(genome_ref_gff, arity: 1, stageAs: 'input/')       // genom_ref_fa[0] is fasta file [recquired] and genom_ref_fa[1] is its index (fasta.fai) [optional but save time if its provided in case of large genome]

  output:
  tuple val(meta), path("${meta.id}.tsv", type: 'file') , optional:false , emit: raw_tsv
  tuple val(meta), path("${meta.id}_filter.tsv", type: 'file') , optional:false , emit: filter_iSNV_tsv
  tuple val(meta), path("${meta.id}_iSNVpos.txt", type: 'file') , optional:false , emit: filter_iSNV_pos

  script:
  genome_ref_fa = genome_ref[0]
  genome_ref_gff = genome_ref[1]

  """
  #!/usr/bin/bash

  samtools mpileup \\
    -aa -A -B -Q 0 \\
    --reference ${genome_ref_fa} \\
    ${task.ext.args} \\
    ${bam[0]} > ${meta.id}.mpileup

  cat ${meta.id}.mpileup | ivar variants \\
    -p ${meta.id}_raw \\
    -q 0 -t 0 -m 0 \\
    -r ${genome_ref_fa} \\
    -g ${genome_ref_gff} \\
    ${task.ext.args2}

  ## complete tsv result with unvariable position by adding of coverage and RefNucl

  # recup coverage and RefNucl. TODO: recup 
  nbCol = \$(awk -F'\t' 'NR==1 { print NF }' ${meta.id}_raw.tsv)
  awk -F '\t' -v nbCol="\$nbCol" '{
    printf \$1"\t"\$2"\t"\$3"\tNA\t"\$4"\tNA\tNA\tNA\tNA\tNA\t0\t"\$4; 
    for (i = 13; i <= nbCol; i++) {
      printf "\tNA";
    } 
    print "";
  }' ${meta.id}.mpileup > ${meta.id}.cov

  # add remain line
  sed '1d' ${meta.id}_raw.tsv > ${meta.id}_tmp.tsv
  grep -v "^\$(sed '1d' ${meta.id}_raw.tsv | cut -f1,2)" ${meta.id}.cov >> ${meta.id}_tmp.tsv  # string length may exeed bash limit!!!
  head -n1 ${meta.id}_raw.tsv > ${meta.id}.tsv
  sort -k1,1 -k2,2n -k8r,8n -s ${meta.id}_tmp.tsv >> ${meta.id}.tsv
  rm ${meta.id}_tmp.tsv

  ## Filter iSNVs: (Note: in case of sample anlysed individually (not in batch context), this is the real results).
  awk -F '\t' \\
    -v minDepth='${task.ext.args3 ?: '30'}' \\
    -v maxRefDepth='${task.ext.args4 ?: '0.8'}' '{
    if (\$12 >= minDepth && \$5/\$12 < maxRefDepth) 
      print \$0
    }' ${meta.id}_raw.tsv > ${meta.id}_filter.tsv

  ## Pos with iSNVs
  cut -f1,2 ${meta.id}_filter.tsv > ${meta.id}_iSNVpos.txt
  
  """

  stub:  
  """
  #!/usr/bin/bash

  touch ${meta.id}.mpileup ${meta.id}.cov ${meta.id}_raw.tsv ${meta.id}.tsv ${meta.id}_filter.tsv ${meta.id}_iSNVpos.txt
  """ 
}