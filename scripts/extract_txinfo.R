library(data.table)
library(GenomicFeatures)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
cat(paste(c('...input gtf and output path:', args, '\n'), collapse = '\n'), file = stderr())

extract_txinfo <- function(gtf_path){
    txdb <- makeTxDbFromGFF(gtf_path)
    
    txlen <- transcriptLengths(txdb, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE)
    gtf_track <- import(gtf_path, 'gtf')
    gtf_track <- as.data.frame(gtf_track)
    scols <- c('gene_id', 'gene_name', 'gene_biotype', 'transcript_id',
               'transcript_biotype', 'protein_id', 'seqnames', 'strand')
    scols <- scols[scols %in% colnames(gtf_track)]
    txtype <- gtf_track[, scols]
    setDT(txtype)
    txtype <- unique(txtype)
    res <- merge(txlen, txtype, by.x = c('gene_id', 'tx_name'), by.y = c('gene_id', 'transcript_id'))
    res <- as.data.table(res)
    res <- res[order(gene_id, tx_name, protein_id)][!duplicated(tx_name)]
    setnames(res, 'seqnames', 'chrom')
    return(res)
}

dtt <- extract_txinfo(args[1])
fwrite(dtt, file = args[2], sep = '\t')
