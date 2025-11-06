cnv <- read_csv("~/Downloads/variants_for_nstd186.csv") %>%
  filter(`Variant Region type` == "copy number variation") %>%
  filter(`Assembly Name` == "GRCh38.p12")

cnv <- cnv %>%
  mutate(start = ifelse(is.na(Start), cnv$`Inner Start`, cnv$Start),
         end = ifelse(is.na(End), cnv$`Inner End`, cnv$End)) %>%
  select(Chromosome, start, end) %>% na.omit()

View(cnv)

cnv_ranges <- GRanges(
  seqnames = paste0("chr", cnv$Chromosome),
  ranges = IRanges(start = cnv$start, end = cnv$end)  # BED is 0-based, GRanges is 1-based
)

overlaps_cs_cnv <- findOverlaps(snp_ranges, cnv_ranges)
if (length(overlaps_cs_cnv) > 0) {
  cnv_queries <- unique(queryHits(overlaps_cs_cnv))
  joined_df_no_cnv <- joined_df[-cnv_queries, , drop = FALSE]
  snp_ranges_no_cnv <- snp_ranges[-cnv_queries]
}
overlaps <- findOverlaps(snp_ranges, prmtrs)
prmrt_genes <- prmtrs$gene_id[subjectHits(overlaps)]
joined_df$ENTREZID <- NA
joined_df$ENTREZID[queryHits(overlaps)] <- prmrt_genes
joined_df$SYMBOL <- mapIds(org.Hs.eg.db,
                           keys = as.character(joined_df$ENTREZID),
                           column = "SYMBOL",
                           keytype = "ENTREZID")
joined_df <- na.omit(joined_df)
