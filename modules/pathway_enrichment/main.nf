process PATHWAY_ENRICHMENT {
  tag "Pathway enrichment"
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq'

  input:
  val species
  path ensembl
  path DDS
  path genexp_table
  
  output:
  path 'pathway_enrichment.Rds', emit: pathway_enrichment_table
  path 'pathway_enrichment_S2.xlsx'

  script:
  """
  # clusterProfiler needs to cache database
  export HOME=\$PWD
  export XDG_DATA_HOME=\$PWD/.local/share
  export XDG_CACHE_HOME=\$PWD/.cache
  mkdir -p \$XDG_DATA_HOME \$XDG_CACHE_HOME
  
  pathway_enrichment.R --species ${species}
  """
}