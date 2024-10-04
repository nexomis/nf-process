process TRANSFERT_GFF {
  container "${params.biocontainers_registry ?: 'quay.io'}/nexomis/py-bioseq:102024"

  input:
  tuple val(meta), path(sample_fa, arity: 1, stageAs: 'input/')  // 'stageAs' and templates ?
  tuple val(meta2), path(ref, arity: 2, stageAs: 'input/')  // 'stageAs' and templates ?
  //val(save_psa)               // boolean
  //val(include_metada_in_gff)  // boolean

  output:
  tuple val(meta), path("${meta.id}_genomicCoords.csv", type: 'file') , optional:false ,  emit: genomic_coords
  tuple val(meta), path("${meta.id}_transferredAnnotation.gff", type: 'file') , optional:false ,  emit: transfered_gff
  tuple val(meta), path("${meta.id}_*_vs_*_globalAlgn.txt", type: 'file') , optional:true ,  emit: psa
  //don-t work: tuple val(meta), path([sample, "${meta.id}_transferedAnnotation.gff"], type: 'file') , optional:false ,  emit: fa_and_transfered_gff


  script:
  annot_fa = ref[0]
  annot_gff = ref[1]
  template "transfer_annot_by_global_psa.py"

  stub:
  """
  touch ${meta.id}_transferedAnnotation.gff ${meta.id}_genomicCoords.csv
  """
}
