# counts <- read.table("counts/counts_in_vitro_17.09.tsv", header = TRUE, row.names = 1, sep = "\t")[, c(6:11)]
# colnames(counts) <- gsub("X.media.leon.DISK2.icig.done.alignments.in_vitro.RNASEQ_", "", gsub("R.bam", "", colnames(counts)))
# coldata <- data.frame(patient = substr(colnames(counts), 1, 3),
#                       condition = substr(colnames(counts), 5, 5)) %>%
#   mutate(patient = as.factor(patient), condition = as.factor(condition))
# rownames(coldata) <- colnames(counts)
# write.table(counts, "counts/in_vitro_counts.tsv", sep = "\t")
# write.table(coldata, "counts/in_vitro_coldata.tsv", sep = "\t")

find_tf_degs <- function(counts, coldata, tfs, comparisons_df) {
  tf_counts <- counts %>%
    mutate(ensembl = rownames(.)) %>%
    inner_join(annotations_clean) %>%
    filter(symbol %in% tfs) %>% 
    distinct(symbol, .keep_all = TRUE)
  
  rownames(tf_counts) <- tf_counts$symbol
  tf_counts <- tf_counts %>%
    select(-c(ensembl, entrezid, symbol))
  dds <- DESeqDataSetFromMatrix(countData = tf_counts,
                                colData = coldata,
                                design = ~ patient) # condition + patient?
  dds <- estimateSizeFactors(dds)
  keep <- rowSums(counts(dds, normalized = TRUE) >= 5 ) >= length(unique(coldata$patient)) - 1
  dds <- dds[keep, ]
  dds <- DESeq(dds)

  cmbntns <- unique(comparisons_df$patients)
  
  tf_deseq_results <- data.frame()
  for (pair in cmbntns) {
    p1 <- gsub("(s\\w+)_s\\w+", "\\1", pair)
    p2 <- gsub("s\\w+_(s\\w+)", "\\1", pair)
    
    res_pairwise <- as.data.frame(results(dds, contrast = c("patient", p1, p2))) %>%
      mutate(symbol = rownames(.), contrast = paste0(p1, "_vs_", p2)) %>%
      filter(symbol %in% comparisons_df$tf[comparisons_df$patients == pair]) 
    
    tf_deseq_results <- bind_rows(tf_deseq_results, res_pairwise)
  }
  
  tf_deseq_results <- tf_deseq_results %>% 
    filter(padj < 0.1 & abs(log2FoldChange) > 0.5)
  return(tf_deseq_results)
}
