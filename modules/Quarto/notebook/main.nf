process QUARTO_NOTEBOOK {
  tag "Generating results report"
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq'

  input:
  path quarto_results_file
  path quarto_meta
  path quarto_css
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

  quarto render ${quarto_results_file} \
    --metadata-file ${quarto_meta} \
    --css ${quarto_css}
  """
}