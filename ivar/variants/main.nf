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
  tuple val(meta), path("${meta.id}.tsv", type: 'file') , optional:true , emit: complete_iSNV_tsv
  tuple val(meta), path("${meta.id}_raw.tsv", type: 'file') , optional:false , emit: raw_iSNV_tsv
  tuple val(meta), path("${meta.id}_filter.tsv", type: 'file') , optional:false , emit: filter_iSNV_tsv
  tuple val(meta), path("${meta.id}_iSNVpos.txt", type: 'file') , optional:false , emit: filter_iSNV_pos
  tuple val(meta), path("${meta.id}_col1-4.mpileup", type: 'file') , optional:false , emit: mpileup_cov

  script:
  genome_ref_fa = genome_ref[0]
  genome_ref_gff = genome_ref[1]

  """
  #!/usr/bin/bash

  samtools mpileup \\
    -aa -A -B -Q 0 \\
    --reference ${genome_ref_fa} \\
    ${task.ext.args ?: ''} \\
    ${bam[0]} > ${meta.id}.mpileup

  cat ${meta.id}.mpileup | ivar variants \\
    -p ${meta.id}_raw \\
    -q 0 -t 0 -m 0 \\
    -r ${genome_ref_fa} \\
    -g ${genome_ref_gff} \\
    ${task.ext.args2 ?: ''}

  ## Filter iSNVs: (Note: in case of sample anlysed individually (not in batch context), this is the real results).
  awk -F '\\t' \\
    -v minDepth='${task.ext.args3 ?: "30"}' \\
    -v maxRefDepth='${task.ext.args4 ?: "0.8"}' '
    NR == 1 { print \$0; next; } {
    if (\$12 >= minDepth && \$5/\$12 < maxRefDepth) 
      print \$0
    }' ${meta.id}_raw.tsv > ${meta.id}_filter.tsv

  ## Pos with iSNVs
  sed '1d' ${meta.id}_filter.tsv | cut -f1,2 > ${meta.id}_iSNVpos.txt
  
  ## Coverage
  cut -f1-4 ${meta.id}.mpileup > ${meta.id}_col1-4.mpileup


  #### only if requested: use boolean input? in case of batch
  ## complete tsv result with unvariable position by adding of coverage and RefNucl
  # recup coverage and Nucl. of all pos: variable and unvariable (not necessary ref.nucl ?! but in all case not impact the following because about this generate file, only unvariable region while used!) 
  nbCol=\$(awk -F'\\t' 'NR==1 { print NF }' ${meta.id}_raw.tsv)
  awk -F '\\t' -v nbCol="\$nbCol" '{
    printf \$1"\\t"\$2"\\t"\$3"\\tNA\\t"\$4"\\tNA\\tNA\\tNA\\tNA\\tNA\\t0\\t"\$4; 
    for (i = 13; i <= nbCol; i++) {
      printf "\\tNA";
    } 
    print "";
  }' ${meta.id}.mpileup > ${meta.id}.cov
  # TSV: add remain line (variable + unvaribale)
  sed '1d' ${meta.id}_raw.tsv > ${meta.id}_tmp.tsv
  sed '1d' ${meta.id}_raw.tsv | awk -F "\\t" '{print "^"\$1"\\t"\$2"\t"}' | sort -u > ${meta.id}_pattern_varPos.txt
  grep -v -f ${meta.id}_pattern_varPos.txt ${meta.id}.cov >> ${meta.id}_tmp.tsv
  head -n1 ${meta.id}_raw.tsv > ${meta.id}.tsv
  sort -k1,1 -k2,2n -k8r,8n -s ${meta.id}_tmp.tsv >> ${meta.id}.tsv
  rm ${meta.id}_tmp.tsv ${meta.id}_pattern_varPos.txt

  """

  stub:  
  """
  #!/usr/bin/bash

  touch ${meta.id}.mpileup ${meta.id}.cov ${meta.id}_raw.tsv ${meta.id}.tsv ${meta.id}_filter.tsv ${meta.id}_iSNVpos.txt
  """ 
}