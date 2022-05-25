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
    
    gtf_cmd <- paste(ifelse(endsWith(gtf_path, 'gz'), 'zcat', 'cat'), gtf_path, '| grep -v "^#"')
    gtf_dtt <- fread(cmd = gtf_cmd, header = FALSE)
    
    gtf_dtt[, tx_name := sub('.*?transcript_id "(ENS[A-Z]*T\\d+)".*', '\\1', V9)]
    gtf_dtt[!grepl("^ENS[A-Z]*T\\d+$", tx_name), tx_name := NA_character_]
    gtf_dtt[, gene_id := sub('.*?gene_id "(ENS[A-Z]*G\\d+)".*', '\\1', V9)]
    
    cds_problem <- merge(gtf_dtt[V3 == 'transcript'][grepl('cds_start_NF', V9), .(tx_name, cds_start_nf = TRUE)],
          gtf_dtt[V3 == 'transcript'][grepl('cds_end_NF', V9), .(tx_name, cds_end_nf = TRUE)],
          by = 'tx_name', all = TRUE)
    mrna_problem <- merge(gtf_dtt[V3 == 'transcript'][grepl('mRNA_start_NF', V9), .(tx_name, mrna_start_nf = TRUE)],
                         gtf_dtt[V3 == 'transcript'][grepl('mRNA_end_NF', V9), .(tx_name, mrna_end_nf = TRUE)],
                         by = 'tx_name', all = TRUE)
    txproblem <- merge(cds_problem, mrna_problem, by = 'tx_name', all = TRUE)
    res <- merge(res, txproblem, by = 'tx_name', all = TRUE)
    res[is.na(cds_start_nf), cds_start_nf := FALSE]
    res[is.na(cds_end_nf), cds_end_nf := FALSE]
    res[is.na(mrna_start_nf), mrna_start_nf := FALSE]
    res[is.na(mrna_end_nf), mrna_end_nf := FALSE]
    return(res)
}

dtt <- extract_txinfo(args[1])
fwrite(dtt, file = args[2], sep = '\t')
