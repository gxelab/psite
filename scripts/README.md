## Helper scripts

`extract_txinfo.R` and `extract_txinfo_ensembl.R` can be used to extract basic informations of all the transcripts from gene annotations obtained from public databases like the [Ensembl Genome Browser](https://www.ensembl.org/index.html?redirect=no) in GTF format.

```bash
Rscript --vanilla extract_txinfo_ensembl.R ensembl_gene_annotations.gtf txinfo.tsv
```

To use either script, R packages [`data.table`](https://cran.r-project.org/web/packages/data.table/index.html), [`GenomicFeatures`](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) and [`rtracklayer`](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html) should be installed in advance. There are two parameters, the first is the path of input gene annotation file in `GTF` format, and the other is the output path of transcript informations in tsv format as follows:

```
gene_id tx_name tx_id   nexon   tx_len  cds_len utr5_len        utr3_len        gene_name       gene_biotype    transcript_biotype      protein_id      chrom   strand
ENSG00000000003 ENST00000373020 234488  8       3768    738     112     2918    TSPAN6  protein_coding  protein_coding  ENSP00000362111 X       -
ENSG00000000003 ENST00000494424 234492  6       820     0       0       0       TSPAN6  protein_coding  processed_transcript            X       -
ENSG00000000003 ENST00000496771 234491  6       1025    0       0       0       TSPAN6  protein_coding  processed_transcript            X       -
ENSG00000000003 ENST00000612152 234489  7       3796    390     489     2917    TSPAN6  protein_coding  protein_coding  ENSP00000482130 X       -
ENSG00000000003 ENST00000614008 234490  7       900     411     489     0       TSPAN6  protein_coding  protein_coding  ENSP00000482894 X       -
ENSG00000000005 ENST00000373031 230753  7       1205    954     83      168     TNMD    protein_coding  protein_coding  ENSP00000362122 X       +
ENSG00000000005 ENST00000485971 230754  3       542     0       0       0       TNMD    protein_coding  processed_transcript            X       +
ENSG00000000419 ENST00000371582 220002  10      1161    864     32      265     DPM1    protein_coding  protein_coding  ENSP00000360638 20      -
ENSG00000000419 ENST00000371584 220007  10      1084    888     9       187     DPM1    protein_coding  protein_coding  ENSP00000360640 20      -
ENSG00000000419 ENST00000371588 220003  9       1054    783     9       262     DPM1    protein_coding  protein_coding  ENSP00000360644 20      -
ENSG00000000419 ENST00000413082 220013  8       672     663     9       0       DPM1    protein_coding  protein_coding  ENSP00000394921 20      -
```
