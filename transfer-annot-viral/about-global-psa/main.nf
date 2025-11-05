process TRANSFERT_GFF {
  container "quay.io/nexomis/py-bioseq:102024"
  tag "$meta.id"
  input:
  tuple val(meta), path(sample_fa, arity: 1, stageAs: 'input/')
  tuple val(meta2), path(ref, arity: 2, stageAs: 'input/')


  output:
  tuple val(meta), path("${meta.label ?: meta.id}_genomicCoords.csv", type: 'file') , optional:false ,  emit: genomic_coords
  tuple val(meta), path("${meta.label ?: meta.id}_transferredAnnotation.gff", type: 'file') , optional:false ,  emit: transfered_gff
  tuple val(meta), path("${meta.label ?: meta.id}_*_vs_*_globalAlgn.txt", type: 'file') , optional:true ,  emit: psa


  script:
  def annot_fa = ref[0]
  def annot_gff = ref[1]
  def out_prefix = meta.label ?: meta.id
  template "transfer_annot_by_global_psa.py"

  stub:
  """
  touch ${meta.label ?: meta.id}_transferredAnnotation.gff ${meta.label ?: meta.id}_genomicCoords.csv
  """
}
