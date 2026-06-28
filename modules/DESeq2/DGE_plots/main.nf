process DESEQ2_DGE_PLOTS {
  tag "Rendering DGE plots"
  label 'process_parallel'

  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq'
  //container "${ workflow.containerEngine == 'singularity' ?
  //  'docker://ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' :
  //  'ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' }"

  input:
  path DDS
  path genexp_table

  output:
  path 'DGE_plots.pdf'
  path 'list_of_genes.txt'

  script:
  def list_of_genes   = "list_of_genes.txt"
  def batchSize       = 500
  
  """
  # ---- 1. generate all plots individually -----------------------------------
  DESeq2_DGE_plots.R --cores ${task.cpus}

  
  # ---- 2. merge all plots into one ------------------------------------------
  # sort -V = sort numerically
  ls *.pdf | sort -V > ${list_of_genes}

  wc ${list_of_genes}

  mkdir batches

  # Split into batches
  split -l ${batchSize} "${list_of_genes}" "batches/batch_"

  # Merge each batch into a temp PDF
  batch_pdfs=()
  for f in batches/batch_*; do
    batch_name=\$(basename "\$f")
    out_pdf="\${f}.pdf"

    echo "Merging batch \$batch_name -> \$out_pdf"
  
    qpdf --empty --pages \$(cat "\$f") -- \$out_pdf
      
    batch_pdfs+=("\$out_pdf")
  done

  # Merge all batches into final PDF
  qpdf --empty --pages \${batch_pdfs[@]} -- DGE_plots.pdf
  """
}