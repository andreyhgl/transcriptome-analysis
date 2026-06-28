process DESEQ2_QC_PLOTS {
  tag "Generating QC plots"
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq'

  //container "${ workflow.containerEngine == 'singularity' ?
  //  'docker://ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' :
  //  'ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' }"

  input:
  path DDS

  output:
  path '*.pdf'

  """
  DESeq2_QC_plots.R
  """
}