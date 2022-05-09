#### Utility scripts

- `extract_txinfo.R`: extract basic informations of all the transcripts from gene annotations stored in GTF format.
  To use this R script, R packages `data.table`(cran), `GenomicFeatures`(bioconductor) and `rtracklayer` (bioconductor) should be installed in advance.
  ```bash
  Rscript --vanilla /path/to/extract_txinfo.R input.gtf output.tsv
  ```