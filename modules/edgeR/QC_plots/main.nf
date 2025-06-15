process EDGER_QC_PLOTS {
  //tag 
  //label

  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq'

  time 1.h
  memory 8.GB
  cpus 1

  input:
  path DGEList

  output:
  path '*.pdf'

  """
  edgeR_QC_plots.R
  """
}