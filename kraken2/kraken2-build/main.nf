process KRAKEN2_BUILD_CUSTOM_DB {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_1"  // too big img !! Create a new one dedicated to "simple" index building from 'ghcr.io/nexomis/kraken2:2.1.3' and including missing tools (rsync, tools for remove low complexity, ... but no need to include all recquired tools (ex: peptide sequence inclusion tools))

  label 'cpu_high'
  label 'mem_4G_per_cpu'
  
  input:
  tuple val(meta), path(fa_to_add, arity: 1, stageAs: 'ref_raw/*')  // meta.id: name of db ; 'fa_to_add' can be path to file or to dir.
  val library
  // TODO (ajust this process to avoid creation of new process)
  //  - manage case where library is a list
  //  - manage case where fa_to_add is null and library is defined (sample dependant or run dependant ??) -> renome process to simply KRAKEN2_BUILD ?
  //    (use "when" before script to check at least one of library or fa_to_add are not null)

  output:
  tuple val(meta), path("${meta.id}", type: 'dir'), optional:false, emit: db
  tuple val(meta), path("${meta.id}.log", type: 'file'), optional:false, emit: log
  // output_dir usefull ?
  
  script:

  // --download-library: maybe is better with "--use-ftp" (instead of rsync)
  // --add-to-library: maybe is better if default process is '--no-masking' (e.g: for classic virus this simplify dowstream assembly whithout risk to exclude viral reads at kraken host exclussion step and probably accelerade krakn index building).
  // --build: maybe is better if default process is : '--fast-build' & '--skip-maps'
  // --clean: or "rm -rf ${meta.id}/taxonomy/ ${meta.id}/library" ?

  """ 
  #!/usr/bin/bash

  kraken2-build --db ${meta.id}/ --threads $task.cpus --download-taxonomy ${task.ext.args_dl-tax ?: ''} 2>${meta.id}.log

  ${library ? "kraken2-build --db ${meta.id}/ --threads $task.cpus --download-library $library ${task.ext.args_dl-lib ?: ''} 2>>${meta.id}.log" : "# pass '--download-library' library step because no 'library' has specified"}

  if [ -d "${fa_to_add}" ]; then
      find ${fa_to_add} -type f \\( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" -o -name "*.fa.gz" -o -name "*.fna.gz" -o -name "*.fasta.gz" \\) | while read fa_file; do
          kraken2-build --db ${meta.id}/ --add-to-library "$fa_file" ${task.ext.args_add-lib ?: ''} 2>>${meta.id}.log
      done
  else
      kraken2-build --db ${meta.id}/ --add-to-library ${fa_to_add} ${task.ext.args_add-lib ?: ''} 2>>${meta.id}.log
  fi

  kraken2-build --build --db ${meta.id}/ --threads $task.cpus ${task.ext.args_build ?: ''} 2>>${meta.id}.log

  kraken2-build --clean --db ${meta.id}/
  """

  stub:
  """
  mkdir ${meta.id}/
  touch ${meta.id}.log
  """
}

