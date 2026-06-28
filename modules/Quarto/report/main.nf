process QUARTO_REPORT {
  tag "Generating report (dashboard)"

  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq'

  time 1.h
  cpus 1

  input:
  path quarto_report_file
  path quarto_assets
  path ensembl_table
  path DDS
  path genexp_table
  path pathway_enrichment

  output:
  path "*.html"

  script:
  """
  # Set environment variables needed for Quarto rendering
  export XDG_DATA_HOME="./cache/"
  export XDG_CACHE_HOME="./cache/"
  export XDG_RUNTIME_DIR="./cache/"

  # TMPDIR is used by Quarto.makeTempDirSync(), set manually to avoid HPC issues
  export TMPDIR="./cache/"
  
  #echo "\${XDG_DATA_HOME}"
  #echo "\${XDG_CACHE_HOME}"
  #echo "\${XDG_RUNTIME_DIR}"
  #echo "\${TMPDIR}"

  quarto render ${quarto_report_file}
  """
}