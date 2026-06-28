process METADATA {
  tag "Generating metadata"
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq'

  //container "${ workflow.containerEngine == 'singularity' ?
  //  'docker://ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' :
  //  'ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' }"

  input:
  path sample_info_file
  path seqfiles_metadata
  path samples_to_remove

  output:
  path 'metadata.csv', emit: metadata

  script:
  """
  metadata.R --sample_info ${sample_info_file}
  """
}