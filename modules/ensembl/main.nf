process ENSEMBL {
  tag "Building database for ${reference_genome}"
  //label
  
  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq'
  
  time 1.h
  memory 8.GB
  cpus 1

  input:
  val reference_genome

  output:
  path 'ensembl_dataset.csv.gz', emit: ENSEMBL_DATASET

  script:
  """
  # set environment variables for biomart cache
  export BIOMART_CACHE="\$SNIC_TMP"

  ensembl.R ${reference_genome}
  """
}