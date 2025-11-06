set.seed(988482)

library(readr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(motifbreakR)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyr)
library(jsonlite)
library(purrr)
library(dplyr)
library(MBASED)
library(DESeq2)

add_rsids <- function(chr_vec, pos_vec, ref_vec) {
  snp_ranges <- GRanges(seqnames = as.character(chr_vec),
                        IRanges(pos_vec, width = 1), ref = ref_vec)
  matched_snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, snp_ranges)
  df_matches <- as.data.frame(matched_snps)[, c("seqnames", "pos")]
  df_matches$rsid <- mcols(matched_snps)$RefSNP_id
  df_matches <- na.omit(df_matches)
  return(df_matches)
}

########### SNP из чипсеков
dm <- read_tsv("/media/leon/Polina/diabet/merge/Damarov.vcf", show_col_types = FALSE)
an <- read_tsv("/media/leon/Polina/diabet/merge/E.Viacheslavovna.vcf", show_col_types = FALSE)
bn <- read_tsv("/media/leon/Polina/diabet/merge/N.Petrovna.vcf", show_col_types = FALSE)

dm$patient <- "D"; an$patient <- "A"; bn$patient <- "B"

dm$variant <- paste(dm$chr, dm$pos, sep = '_')
an$variant <- paste(an$chr, an$pos, sep = '_')
bn$variant <- paste(bn$chr, bn$pos, sep = '_')

all(all(colnames(dm) == colnames(an)), all(colnames(dm) == colnames(bn)))
joined_df <- rbind(an, dm, bn)
joined_df$chr <- as.character(as.numeric(joined_df$chr))
joined_df <- na.omit(joined_df)

# выбираем SNP которые есть хотя бы у 2 из 3 людей
joined_df <- joined_df %>%
  group_by(variant) %>% filter(n() > 1) %>% ungroup() %>%
  mutate(alt = ifelse(allel1 != ref, allel1, allel2)) %>%
  select(-DP, -QS1, -QS2, -allel1, -allel2)

# пересекаем с промоторами
hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
prmtrs <- promoters(hg38, columns = c("gene_id"))
exns <- exons(hg38, columns = c("exon_id", "gene_id"))

snp_ranges <- GRanges(
  seqnames = paste0("chr", joined_df$chr),
  ranges = IRanges(start = joined_df$pos, width = 1)
)
overlaps <- findOverlaps(snp_ranges, prmtrs)
prmrt_genes <- prmtrs$gene_id[subjectHits(overlaps)]
joined_df$ENTREZID <- NA
joined_df$ENTREZID[queryHits(overlaps)] <- prmrt_genes
joined_df$SYMBOL <- mapIds(org.Hs.eg.db,
                           keys = as.character(joined_df$ENTREZID),
                           column = "SYMBOL",
                           keytype = "ENTREZID")
joined_df <- na.omit(joined_df)
joined_df <- merge(joined_df, add_rsids(joined_df$chr, joined_df$pos, joined_df$ref), by.x = c("chr", "pos"), by.y = c("seqnames", "pos"))

######### in vitro RNA-seq
stats_dir <- "/media/leon/DISK2/icig/done/all_stats/"

ASE_in_vitro <- function(patient, min_DP = 100) {
  if (!is.na(as.numeric(substr(patient, 2, 3)))) {
    stats0R <- read_delim(paste0(stats_dir, patient, "_0R.stat"), show_col_types = FALSE) %>%
      filter(DP_total >= min_DP)
    stats1R <- read_delim(paste0(stats_dir, patient, "_1R.stat"), show_col_types = FALSE) %>%
      filter(DP_total >= min_DP)
  } else {
    stats0R <- read_tsv(paste0(stats_dir, patient, "_DR.stat"), show_col_types = FALSE) %>%
      filter(DP_total >= min_DP)
    stats1R <- read_tsv(paste0(stats_dir, patient, "_ER.stat"), show_col_types = FALSE) %>%
      filter(DP_total >= min_DP)
  }
  stats <- merge(stats0R, stats1R, by = c("CHR", "POS"))
  colnames(stats) <- c("chr", "pos", "id_sample_0R", "DP_0R", "ref_0R", "QS_ref_0R", "alt1_0R", "QS_alt1_0R", "alt2_0R", "QS_alt2_0R", "id_sample_1R", "DP_1R", "ref_1R", "QS_ref_1R", "alt1_1R", "QS_alt1_1R", "alt2_1R", "QS_alt2_1R")
  stats$chr <- as.character(as.numeric(stats$chr))
  
  stats <- stats[stringr::str_length(stats$ref_0R) == 1, ]
  stats <- stats[stringr::str_length(stats$alt1_0R) == 1, ]
  stats <- stats[stringr::str_length(stats$alt1_1R) == 1, ]
  
  stats[stats$alt1_0R != stats$alt1_1R, c(7, 8, 9, 10, 15, 16, 17, 18)] <- stats[stats$alt1_0R != stats$alt1_1R, c(7, 8, 9, 10, 17, 18, 15, 16)]
  
  stats$p01 <-  apply(stats[, c(4, 6, 8, 12, 14, 16)], 1, function(x) {
    data_matrix <- matrix(unlist(c(round(x[1] * x[2]), round(x[1] * x[3]), round(x[4] * x[5]), round(x[4] * x[6]))), nrow = 2, byrow = TRUE)
    v <- chisq.test(data_matrix)
    return(v$p.value)
  })
  stats <- na.omit(stats)
  
  gr_positions <- GRanges(
    seqnames = paste0("chr", stats$chr),
    ranges = IRanges(start = stats$pos, end = stats$pos)
  )
  overlaps <- findOverlaps(gr_positions, exns)
  stats <- stats[queryHits(overlaps), ]
  snp_exons <- exns[subjectHits(overlaps)]
  gene_ids <- sapply(snp_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
  stats$ENTREZID <- gene_ids
  stats$SYMBOL <- mapIds(org.Hs.eg.db,
                         keys = stats$ENTREZID,
                         keytype = "ENTREZID",
                         column = "SYMBOL")
  stats <- stats %>%
    unnest(SYMBOL) %>%
    filter(!grepl("HLA-", SYMBOL)) %>%
    unique() %>% na.omit() %>%
    mutate(patient = gsub("._.*", "", gsub("RNASEQ_s", "", id_sample_0R)),
           ref = ref_0R, alt = alt1_0R) %>%
    select(-contains("id_sample_"), -contains("alt2"), -ref_0R, -ref_1R, -alt1_0R, -alt1_1R) %>%
    mutate(log2FC_1vs0 = log2(QS_ref_1R * QS_alt1_0R / QS_alt1_1R / QS_ref_0R))
  
  if (!is.na(as.numeric(substr(patient, 2, 3)))) { stats$patient = patient }
  return(stats)
}

patients_vtr <- c("sAn", "sBn", "sDm")
setwd("/media/leon/DISK2/icig/done/all_stats/")
ase_in_vitro <- lapply(patients_vtr, function(x) {
  ASE_in_vitro(x, min_DP = 100)
})

ase <- bind_rows(ase_in_vitro)
ase$padj <- p.adjust(ase$p01, method = "BH")
ase$diffASE <- ase$padj < 0.01

# 843 генf с diffASE при FDR = 0.01 и min_DP = 100
ase %>% filter(diffASE == TRUE) %>% select(SYMBOL) %>% unique() %>% nrow()

# выбираем те SNP из чипсеков которые находятся в промоторах генов с SNP в экзонах
ase <- ase %>%
  filter(SYMBOL %in% intersect(ase$SYMBOL, joined_df$SYMBOL))

joined_df_filt <- joined_df %>%
  filter(SYMBOL %in% intersect(ase$SYMBOL, joined_df$SYMBOL))

# 3377 генов с SNP в промоторах (общих у 2+ людей) и SNP в экзонах
length(unique(joined_df_filt$SYMBOL))

joined_df_filt$variant <- paste(joined_df_filt$variant, joined_df_filt$ref, joined_df_filt$alt, sep = '_')

df_collapsed <- joined_df_filt %>%
  group_by(chr, pos, ref, alt, variant, ENTREZID, SYMBOL, rsid) %>%
  summarise(patients = paste0(sort(unique(patient)), collapse = ""), .groups = "drop")

#################### для выбранных генов с diffASE считаем ASE на уровне гена
ase$variant <- paste(ase$chr, ase$pos, sep = '_')
ase <- ase %>%
  mutate(DP_ref_0R = round(DP_0R * QS_ref_0R),
         DP_alt_0R = round(DP_0R * QS_alt1_0R),
         DP_ref_1R = round(DP_1R * QS_ref_1R),
         DP_alt_1R = round(DP_1R * QS_alt1_1R)) %>%
  filter(DP_ref_0R > 10 & DP_alt_0R > 10 & DP_ref_1R > 10 & DP_alt_1R > 10)

summarizeASEResults_2s <- function(MBASEDOutput) {
  geneOutputDF <- data.frame(
    majorAlleleFrequencyDifference = assays(MBASEDOutput)$majorAlleleFrequencyDifference[, 1],
    pValueASE = assays(MBASEDOutput)$pValueASE[, 1],
    pValueHeterogeneity = assays(MBASEDOutput)$pValueHeterogeneity[, 1]
  )
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[, 1]
  lociOutputGR$MAFDifference <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAFDifference[, 1]
  lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels = unique(lociOutputGR$aseID)))
  return(
    list(
      geneOutput = geneOutputDF,
      locusOutput = lociOutputList
    )
  )
}  

ase_sub_A <- ase %>%
  filter(patient == "A") %>%
  group_by(SYMBOL) %>%
  group_modify(~ slice_sample(.x, n = min(15, nrow(.x)))) %>%
  ungroup()
 
snps_A <- GRanges(
  seqnames = ase_sub_A$chr,
  ranges   = IRanges(ase_sub_A$pos, width = 1),
  aseID    = ase_sub_A$SYMBOL,
  allele1  = ase_sub_A$ref,
  allele2  = ase_sub_A$alt
)
names(snps_A) <- paste(ase_sub_A$SYMBOL, ase_sub_A$variant, sep = "_")

se_A_0vs1 <- SummarizedExperiment(
  assays = list(
    lociAllele1Counts = matrix(c(ase_sub_A[["DP_ref_0R"]], ase_sub_A[["DP_ref_1R"]]), ncol = 2, dimnames = list(names(snps_A), c("0", "1"))),
    lociAllele2Counts = matrix(c(ase_sub_A[["DP_alt_0R"]], ase_sub_A[["DP_alt_1R"]]), ncol = 2, dimnames = list(names(snps_A), c("0", "1")))
  ),
  rowRanges = snps_A
)
bp <- MulticoreParam(workers = 30)

start_time <- Sys.time()
results_A_0vs1 <- runMBASED(ASESummarizedExperiment = se_A_0vs1, isPhased = FALSE, numSim = 10^5, BPPARAM = bp)
end_time <- Sys.time()
print(end_time - start_time)

gene_output_A_0vs1 <- summarizeASEResults_2s(results_A_0vs1)$geneOutput
gene_output_A_0vs1$patient <- "A"
gene_output_A_0vs1$SYMBOL <- rownames(gene_output_A_0vs1)
names(gene_output_A_0vs1)[names(gene_output_A_0vs1) == "majorAlleleFrequencyDifference"] <- "genewiseMAFDiff"

locus_df_A_0vs1 <- data.frame(
  snp_id = names(unlist(summarizeASEResults_2s(results_A_0vs1)$locusOutput)),
  as.data.frame(unlist(summarizeASEResults_2s(results_A_0vs1)$locusOutput), row.names = NULL)
) %>%
  merge(gene_output_A_0vs1[, c("SYMBOL", "genewiseMAFDiff", "pValueASE", "pValueHeterogeneity")], by.x = "aseID", by.y = "SYMBOL")

full_df_A_0vs1 <- merge(locus_df_A_0vs1, ase_sub_A,
                        by.x = c("seqnames", "start", "aseID"),
                        by.y = c("chr", "pos", "SYMBOL")) %>%
  mutate(comp = "0vs1")

se_A_1vs0 <- SummarizedExperiment(
  assays = list(
    lociAllele1Counts = matrix(c(ase_sub_A[["DP_ref_1R"]], ase_sub_A[["DP_ref_0R"]]), ncol = 2, dimnames = list(names(snps_A), c("1", "0"))),
    lociAllele2Counts = matrix(c(ase_sub_A[["DP_alt_1R"]], ase_sub_A[["DP_alt_0R"]]), ncol = 2, dimnames = list(names(snps_A), c("1", "0")))
  ),
  rowRanges = snps_A
)

start_time <- Sys.time()
results_A_1vs0 <- runMBASED(ASESummarizedExperiment = se_A_1vs0, isPhased = FALSE, numSim = 10^5, BPPARAM = bp)
end_time <- Sys.time()
print(end_time - start_time)

gene_output_A_1vs0 <- summarizeASEResults_2s(results_A_1vs0)$geneOutput
gene_output_A_1vs0$patient <- "A"
gene_output_A_1vs0$SYMBOL <- rownames(gene_output_A_1vs0)
names(gene_output_A_1vs0)[names(gene_output_A_1vs0) == "majorAlleleFrequencyDifference"] <- "genewiseMAFDiff"

locus_df_A_1vs0 <- data.frame(
  snp_id = names(unlist(summarizeASEResults_2s(results_A_1vs0)$locusOutput)),
  as.data.frame(unlist(summarizeASEResults_2s(results_A_1vs0)$locusOutput), row.names = NULL)
) %>%
  merge(gene_output_A_1vs0[, c("SYMBOL", "genewiseMAFDiff", "pValueASE", "pValueHeterogeneity")], by.x = "aseID", by.y = "SYMBOL")

full_df_A_1vs0 <- merge(locus_df_A_1vs0, ase_sub_A,
                        by.x = c("seqnames", "start", "aseID"),
                        by.y = c("chr", "pos", "SYMBOL")) %>%
  mutate(comp = "1vs0")

#######################################
ase_sub_B <- ase %>%
  filter(patient == "B") %>%
  group_by(SYMBOL) %>%
  group_modify(~ slice_sample(.x, n = min(15, nrow(.x)))) %>%
  ungroup()

snps_B <- GRanges(
  seqnames = ase_sub_B$chr,
  ranges   = IRanges(ase_sub_B$pos, width = 1),
  aseID    = ase_sub_B$SYMBOL,
  allele1  = ase_sub_B$ref,
  allele2  = ase_sub_B$alt
)
names(snps_B) <- paste(ase_sub_B$SYMBOL, ase_sub_B$variant, sep = "_")

se_B_0vs1 <- SummarizedExperiment(
  assays = list(
    lociAllele1Counts = matrix(c(ase_sub_B_0vs1[["DP_ref_0R"]], ase_sub_B_0vs1[["DP_ref_1R"]]), ncol = 2, dimnames = list(names(snps_B_0vs1), c("0", "1"))),
    lociAllele2Counts = matrix(c(ase_sub_B_0vs1[["DP_alt_0R"]], ase_sub_B_0vs1[["DP_alt_1R"]]), ncol = 2, dimnames = list(names(snps_B_0vs1), c("0", "1")))
  ),
  rowRanges = snps_B_0vs1
)

start_time <- Sys.time()
results_B_0vs1 <- runMBASED(ASESummarizedExperiment = se_B_0vs1, isPhased = FALSE, numSim = 10^5, BPPARAM = bp)
end_time <- Sys.time()
print(end_time - start_time)

gene_output_B_0vs1 <- summarizeASEResults_2s(results_B_0vs1)$geneOutput
gene_output_B_0vs1$patient <- "B"
gene_output_B_0vs1$SYMBOL <- rownames(gene_output_B_0vs1)
names(gene_output_B_0vs1)[names(gene_output_B_0vs1) == "majorAlleleFrequencyDifference"] <- "genewiseMAFDiff"

locus_df_B_0vs1 <- data.frame(
  snp_id = names(unlist(summarizeASEResults_2s(results_B_0vs1)$locusOutput)),
  as.data.frame(unlist(summarizeASEResults_2s(results_B_0vs1)$locusOutput), row.names = NULL)
) %>%
  merge(gene_output_B_0vs1[, c("SYMBOL", "genewiseMAFDiff", "pValueASE", "pValueHeterogeneity")], by.x = "aseID", by.y = "SYMBOL")

full_df_B_0vs1 <- merge(locus_df_B_0vs1, ase_sub_B,
                        by.x = c("seqnames", "start", "aseID"),
                        by.y = c("chr", "pos", "SYMBOL")) %>%
  mutate(comp = "0vs1")

se_B_1vs0 <- SummarizedExperiment(
  assays = list(
    lociAllele1Counts = matrix(c(ase_sub_B[["DP_ref_1R"]], ase_sub_B[["DP_ref_0R"]]), ncol = 2, dimnames = list(names(snps_B), c("1", "0"))),
    lociAllele2Counts = matrix(c(ase_sub_B[["DP_alt_1R"]], ase_sub_B[["DP_alt_0R"]]), ncol = 2, dimnames = list(names(snps_B), c("1", "0")))
  ),
  rowRanges = snps_B
)

start_time <- Sys.time()
results_B_1vs0 <- runMBASED(ASESummarizedExperiment = se_B_1vs0, isPhased = FALSE, numSim = 10^5, BPPARAM = bp)
end_time <- Sys.time()
print(end_time - start_time)

gene_output_B_1vs0 <- summarizeASEResults_2s(results_B_1vs0)$geneOutput
gene_output_B_1vs0$patient <- "B"
gene_output_B_1vs0$SYMBOL <- rownames(gene_output_B_1vs0)
names(gene_output_B_1vs0)[names(gene_output_B_1vs0) == "majorAlleleFrequencyDifference"] <- "genewiseMAFDiff"

locus_df_B_1vs0 <- data.frame(
  snp_id = names(unlist(summarizeASEResults_2s(results_B_1vs0)$locusOutput)),
  as.data.frame(unlist(summarizeASEResults_2s(results_B_1vs0)$locusOutput), row.names = NULL)
) %>%
  merge(gene_output_B_1vs0[, c("SYMBOL", "genewiseMAFDiff", "pValueASE", "pValueHeterogeneity")], by.x = "aseID", by.y = "SYMBOL")

full_df_B_1vs0 <- merge(locus_df_B_1vs0, ase_sub_B,
                        by.x = c("seqnames", "start", "aseID"),
                        by.y = c("chr", "pos", "SYMBOL")) %>%
  mutate(comp = "1vs0")

#######################################
ase_sub_D <- ase %>%
  filter(patient == "D") %>%
  group_by(SYMBOL) %>%
  group_modify(~ slice_sample(.x, n = min(15, nrow(.x)))) %>%
  ungroup()

snps_D <- GRanges(
  seqnames = ase_sub_D$chr,
  ranges   = IRanges(ase_sub_D$pos, width = 1),
  aseID    = ase_sub_D$SYMBOL,
  allele1  = ase_sub_D$ref,
  allele2  = ase_sub_D$alt
)
names(snps_D) <- paste(ase_sub_D$SYMBOL, ase_sub_D$variant, sep = "_")

se_D_0vs1 <- SummarizedExperiment(
  assays = list(
    lociAllele1Counts = matrix(c(ase_sub_D[["DP_ref_0R"]], ase_sub_D[["DP_ref_1R"]]), ncol = 2, dimnames = list(names(snps_D), c("0", "1"))),
    lociAllele2Counts = matrix(c(ase_sub_D[["DP_alt_0R"]], ase_sub_D[["DP_alt_1R"]]), ncol = 2, dimnames = list(names(snps_D), c("0", "1")))
  ),
  rowRanges = snps_D
)

start_time <- Sys.time()
results_D_0vs1 <- runMBASED(ASESummarizedExperiment = se_D_0vs1, isPhased = FALSE, numSim = 10^5, BPPARAM = bp)
end_time <- Sys.time()
print(end_time - start_time)

gene_output_D_0vs1 <- summarizeASEResults_2s(results_D_0vs1)$geneOutput
gene_output_D_0vs1$patient <- "D"
gene_output_D_0vs1$SYMBOL <- rownames(gene_output_D_0vs1)
names(gene_output_D_0vs1)[names(gene_output_D_0vs1) == "majorAlleleFrequencyDifference"] <- "genewiseMAFDiff"

locus_df_D_0vs1 <- data.frame(
  snp_id = names(unlist(summarizeASEResults_2s(results_D_0vs1)$locusOutput)),
  as.data.frame(unlist(summarizeASEResults_2s(results_D_0vs1)$locusOutput), row.names = NULL)
) %>%
  merge(gene_output_D_0vs1[, c("SYMBOL", "genewiseMAFDiff", "pValueASE", "pValueHeterogeneity")], by.x = "aseID", by.y = "SYMBOL")

full_df_D_0vs1 <- merge(locus_df_D_0vs1, ase_sub_D,
                        by.x = c("seqnames", "start", "aseID"),
                        by.y = c("chr", "pos", "SYMBOL")) %>%
  mutate(comp = "0vs1")

se_D_1vs0 <- SummarizedExperiment(
  assays = list(
    lociAllele1Counts = matrix(c(ase_sub_D[["DP_ref_1R"]], ase_sub_D[["DP_ref_0R"]]), ncol = 2, dimnames = list(names(snps_D), c("1", "0"))),
    lociAllele2Counts = matrix(c(ase_sub_D[["DP_alt_1R"]], ase_sub_D[["DP_alt_0R"]]), ncol = 2, dimnames = list(names(snps_D), c("1", "0")))
  ),
  rowRanges = snps_D
)

start_time <- Sys.time()
results_D_1vs0 <- runMBASED(ASESummarizedExperiment = se_D_1vs0, isPhased = FALSE, numSim = 10^5, BPPARAM = bp)
end_time <- Sys.time()
print(end_time - start_time)

gene_output_D_1vs0 <- summarizeASEResults_2s(results_D_1vs0)$geneOutput
gene_output_D_1vs0$patient <- "D"
gene_output_D_1vs0$SYMBOL <- rownames(gene_output_D_1vs0)
names(gene_output_D_1vs0)[names(gene_output_D_1vs0) == "majorAlleleFrequencyDifference"] <- "genewiseMAFDiff"

locus_df_D_1vs0 <- data.frame(
  snp_id = names(unlist(summarizeASEResults_2s(results_D_1vs0)$locusOutput)),
  as.data.frame(unlist(summarizeASEResults_2s(results_D_1vs0)$locusOutput), row.names = NULL)
) %>%
  merge(gene_output_D_1vs0[, c("SYMBOL", "genewiseMAFDiff", "pValueASE", "pValueHeterogeneity")], by.x = "aseID", by.y = "SYMBOL")

full_df_D_1vs0 <- merge(locus_df_D_1vs0, ase_sub_D,
                        by.x = c("seqnames", "start", "aseID"),
                        by.y = c("chr", "pos", "SYMBOL")) %>%
  mutate(comp = "1vs0")

full_mbased_0vs1 <- bind_rows(full_df_A_0vs1, full_df_B_0vs1, full_df_D_0vs1) %>%
  select(-c(snp_id, end, width, strand, MAFDifference, QS_ref_0R, QS_alt1_0R, QS_ref_1R, QS_alt1_1R, DP_0R, DP_1R, comp)) %>%
  rename(SYMBOL = aseID)

genewise_mbased_0vs1 <- full_mbased_0vs1 %>%
  select(SYMBOL, genewiseMAFDiff,  pValueASE, pValueHeterogeneity, ENTREZID, patient) %>%
  distinct() %>%
  mutate(padj = p.adjust(pValueASE, method = "BH")) %>%
  mutate(diffASE = padj < 0.01) # ?

full_mbased_1vs0 <- bind_rows(full_df_A_1vs0, full_df_B_1vs0, full_df_D_1vs0) %>%
  select(-c(snp_id, end, width, strand, MAFDifference, QS_ref_0R, QS_alt1_0R, QS_ref_1R, QS_alt1_1R, DP_0R, DP_1R, comp)) %>%
  rename(SYMBOL = aseID)

genewise_mbased_1vs0 <- full_mbased_1vs0 %>%
  select(SYMBOL, genewiseMAFDiff,  pValueASE, pValueHeterogeneity, ENTREZID, patient) %>%
  distinct() %>%
  mutate(padj = p.adjust(pValueASE, method = "BH")) %>%
  mutate(diffASE = padj < 0.01) # ?

##################################
# для каждого гена из RNA-seq находим какие SNP из чипсеков (общие у 2+ людей) есть в промоторе этого гена
lfc_divided_0vs1 <- merge(genewise_mbased_0vs1, joined_df_filt, by = "SYMBOL", suffixes = c(".RNA", ".CS"), all.y = TRUE) %>%
  select(-c(ENTREZID.RNA, ENTREZID.CS, pValueASE, pValueHeterogeneity)) %>% 
  group_by(SYMBOL) %>%
  filter(n_distinct(diffASE) >= 2) %>% # есть один человек с diffASE и один без
  ungroup() %>%
  group_by(across(-patient.CS)) %>% 
  summarize(patient.CS = paste(sort(patient.CS), collapse = ""), .groups = "drop")%>% 
  select(SYMBOL, patient.RNA, genewiseMAFDiff) %>%
  distinct() %>% 
  group_by(SYMBOL) %>% filter(n() >= 2) %>% # не нужно?
  group_modify(~ {
    patients <- .x$patient.RNA
    # LFCs <- .x$absLFC
    MAFDiffs <- .x$genewiseMAFDiff
    combn(length(patients), 2, simplify = FALSE) %>%
      map_dfr(function(pair) {
        i <- pair[1]
        j <- pair[2]
        
        patient_1 <- patients[i]
        patient_2 <- patients[j]
        mafDiff_1 <- MAFDiffs[i]
        mafDiff_2 <- MAFDiffs[j]
        
        # считаем LFCratio только если имеется одинаковый SNP в промоторе данного гена у обоих людей
        common_snps <- df_collapsed$rsid[df_collapsed$SYMBOL == .y$SYMBOL & (df_collapsed$patients == paste0(sort(unique(patients)), collapse = "") | df_collapsed$patients == "ABD")]
        mafDiffRatio <- if(length(common_snps) > 0) max(abs(c(mafDiff_1, mafDiff_2))) / min(abs(c(mafDiff_1, mafDiff_2))) else NA
        ddMAF <- if(length(common_snps) > 0) abs(abs(mafDiff_1) - abs(mafDiff_2)) else NA
        
        tibble(
          patients = ifelse(patient_1 < patient_2, paste0(patient_1, patient_2), paste0(patient_2, patient_1)),
          mafDiffRatio = mafDiffRatio,
          ddMAF = ddMAF,
          logMafDiffRatio = if(!is.na(mafDiffRatio)) log2(mafDiffRatio) else NA,
          deltaMAF1 = mafDiff_1,
          deltaMAF2 = mafDiff_2,
          rsids = list(common_snps),
          n_common_snps = length(common_snps)
        )
      })
  }) %>% ungroup() %>% distinct() %>% filter(n_common_snps > 0)

lfc_divided_1vs0 <- merge(genewise_mbased_1vs0, joined_df_filt, by = "SYMBOL", suffixes = c(".RNA", ".CS"), all.y = TRUE) %>%
  select(-c(ENTREZID.RNA, ENTREZID.CS, pValueASE, pValueHeterogeneity)) %>% 
  group_by(SYMBOL) %>% filter(n_distinct(diffASE) >= 2) %>% ungroup() %>%
  group_by(across(-patient.CS)) %>% 
  summarize(patient.CS = paste(sort(patient.CS), collapse = ""), .groups = "drop")%>% 
  select(SYMBOL, patient.RNA, genewiseMAFDiff) %>%
  distinct() %>% 
  group_by(SYMBOL) %>% filter(n() >= 2) %>% 
  group_modify(~ {
    patients <- .x$patient.RNA
    # LFCs <- .x$absLFC
    MAFDiffs <- .x$genewiseMAFDiff
    combn(length(patients), 2, simplify = FALSE) %>%
      map_dfr(function(pair) {
        i <- pair[1]
        j <- pair[2]
        
        patient_1 <- patients[i]
        patient_2 <- patients[j]
        mafDiff_1 <- MAFDiffs[i]
        mafDiff_2 <- MAFDiffs[j]
        
        # считаем LFCratio только если имеется одинаковый SNP в промоторе данного гена у обоих людей
        common_snps <- df_collapsed$rsid[df_collapsed$SYMBOL == .y$SYMBOL & (df_collapsed$patients == paste0(sort(unique(patients)), collapse = "") | df_collapsed$patients == "ABD")]
        mafDiffRatio <- if(length(common_snps) > 0) max(abs(c(mafDiff_1, mafDiff_2))) / min(abs(c(mafDiff_1, mafDiff_2))) else NA
        ddMAF <- if(length(common_snps) > 0) abs(abs(mafDiff_1) - abs(mafDiff_2)) else NA
        
        tibble(
          patients = ifelse(patient_1 < patient_2, paste0(patient_1, patient_2), paste0(patient_2, patient_1)),
          mafDiffRatio = mafDiffRatio,
          ddMAF = ddMAF,
          logMafDiffRatio = if(!is.na(mafDiffRatio)) log2(mafDiffRatio) else NA,
          deltaMAF1 = mafDiff_1,
          deltaMAF2 = mafDiff_2,
          rsids = list(common_snps),
          n_common_snps = length(common_snps)
        )
      })
  }) %>% ungroup() %>% distinct() %>% filter(n_common_snps > 0)

# выбираем 100 генов с наибольшим log(LFCratio)
ddmaf_df <- bind_rows(lfc_divided_0vs1, lfc_divided_1vs0) %>% arrange(desc(ddMAF)) %>% distinct()
topgenes <- (ddmaf_df %>% distinct(SYMBOL) %>% pull(SYMBOL))[1:100]
ddmaf_df_filt <- ddmaf_df[1:221,] # %>% filter(SYMBOL %in% topgenes)
length(unique(ddmaf_df_filt$SYMBOL))  
snps_for_mtfbrkr <- unique(unlist(ddmaf_df_filt$rsids)) # 207
# writeLines(snps_for_mtfbrkr, "snps_for_mtfbrkr_13.09.25.txt")

########### motifbreakR
setwd("/media/leon/DISK2/icig/done/motifbreakr_12_09_25/")
result_files <- list.files(pattern = "motifbreak_.*\\.rds$")

all_results <- list()
for (file in result_files) {
  res <- readRDS(file)
  rsid <- sub("motifbreak_(.*)\\.rds", "\\1", file)
  if (!rsid %in% snps_for_mtfbrkr) { next }
  if (length(res) > 0) {
    mcols(res)$rsid <- rsid
    all_results[[rsid]] <- res
  }
}

combined_results <- do.call(c, unname(all_results))
hcmc <- combined_results[combined_results$dataSource == "HOCOMOCOv13"]

snpgenedf <-  hcmc %>% mcols() %>% as.data.frame() %>% 
  select(geneSymbol, SNP_id) %>% 
  rename(TF = geneSymbol) %>%
  distinct()

diffase_snp_tf <- ddmaf_df_filt %>%
  select(SYMBOL, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(., snpgenedf, by = "SNP_id")

########## дифэкспрессия генов ТФ
counts <- read.table("/media/leon/DISK2/icig/done/in_vitro_results/counts/counts.tsv", header = TRUE, row.names = 1, sep = "\t")
counts <- counts[, c(6, 7, 8, 9, 11, 13)]
colnames(counts) <- gsub("\\.bam", "", colnames(counts))
colnames(counts) <- gsub("RNASEQ_", "", colnames(counts))
coldata <- data.frame(patient = substr(colnames(counts), 1, 3),
                      time = substr(colnames(counts), 5, 5))
rownames(coldata) <- colnames(counts)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ time)
dds <- DESeq(dds)
res <- results(dds)

annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = rownames(res),
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "ENSEMBL")

annotations_ordered <- annotations[match(rownames(res), annotations$ENSEMBL), ]
annotations_clean <- annotations[!is.na(annotations$ENTREZID) & !is.na(annotations$SYMBOL), ]
all_symbols <- annotations_ordered$SYMBOL

TFs <- unique(snpgenedf$TF)
length(unique(snpgenedf$SNP_id)) # 198
length(TFs) # 569
length(intersect(TFs, all_symbols)) # 317

# mapping hocomoco TF names to gene symbols
json_lines <- readLines("/media/leon/DISK2/icig/done/H13CORE_annotation.jsonl")

# Parse each JSON line and extract the TF:gene_symbol mapping
tf_gene_dict <- map_dfr(json_lines, function(line) {
  data <- fromJSON(line)
  
  if (!is.null(data$masterlist_info$species$HUMAN$gene_symbol)) {
    data.frame(
      tf = data$tf,
      motif_name = data$name,
      gene_symbol = data$masterlist_info$species$HUMAN$gene_symbol,
      stringsAsFactors = FALSE
    )
  } else { NULL }
})

# Create a named vector for the dictionary (TF -> gene_symbol)
tf_to_gene <- setNames(tf_gene_dict$gene_symbol, tf_gene_dict$tf)
tf_to_gene["ZN821"] <- "ZNF821"
tf_to_gene["ZN704"] <- "ZNF704"

snpgenedf <- snpgenedf %>%
  mutate(TF = recode(TF, !!!tf_to_gene))

TFs <- unique(snpgenedf$TF)
length(TFs) # 564 (?)
length(intersect(TFs, all_symbols)) # 562

setdiff(TFs, all_symbols)

TF_counts <- counts %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  filter(SYMBOL %in% TFs)

rownames(TF_counts) <- TF_counts$ENSEMBL
TF_counts <- TF_counts[, c(2:7)]
dds <- DESeqDataSetFromMatrix(countData = TF_counts,
                              colData = coldata,
                              design = ~ patient)
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized = TRUE) >= 5 ) >= 2
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- results(dds)

###### выбираем попарно те дэги которые связаны с общими снп в промоторах генов по которым у этих двух людей разница в диффасе
resAB <- as.data.frame(results(dds, contrast = c("patient", "sAn", "sBn"))) %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "AB") %>%
  filter(SYMBOL %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "AB"])

resAD <- as.data.frame(results(dds, contrast = c("patient", "sAn", "sDm"))) %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "AD") %>%
  filter(SYMBOL %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "AD"])

resBD <- as.data.frame(results(dds, contrast = c("patient", "sBn", "sDm"))) %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "BD") %>%
  filter(SYMBOL %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "BD"])

TF_deseq_results <- bind_rows(resAB, resAD, resBD) %>%
  filter(padj < 0.1 & abs(log2FoldChange) > log2(1.5))
TF_deseq_results
all_degs <- unique(TF_deseq_results$SYMBOL)

deg_tf_rsids <- ddmaf_df_filt %>% 
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(., snpgenedf, by = "SNP_id") %>%
  filter(TF %in% all_degs)

TF_genes_df <- deg_tf_rsids %>% group_by(TF) %>% summarize(genes = list(SYMBOL))
length(TF_genes_df$TF) # 45

go_list <- list()
for (i in seq(nrow(TF_genes_df))) {
  gen <-  TF_genes_df[i, ][[1]]
  print(gen)
  try({
    tfdf <- as.data.frame(enrichGO(gene = unique(unlist(TF_genes_df[i, ]$genes)),
                                   OrgDb = org.Hs.eg.db,
                                   keyType = "SYMBOL",
                                   ont = "BP",
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   # universe = res$SYMBOL,
                                   readable = TRUE))
    if (nrow(tfdf) > 0) {
      tfdf$TF <- gen
      go_list[[gen]] <- tfdf
    }
  })
}
go_df <- bind_rows(go_list)

kegg_list <- list()
for (i in seq(nrow(TF_genes_df))) {
  gen <-  TF_genes_df[i, ][[1]]
  gens <- mapIds(org.Hs.eg.db,
                 keys = unique(unlist(TF_genes_df[i, ]$genes)),
                 keytype = "SYMBOL",
                 column = "ENTREZID")
  print(gens)
  try({
    tfdf <- as.data.frame(enrichKEGG(gene = gens,
                                     organism = "hsa",
                                     pAdjustMethod = "BH",
                                     # universe = res$ENTREZID,
                                     qvalueCutoff = 0.1))
    if (nrow(tfdf) > 0) {
      tfdf$TF <- gen
      kegg_list[[gen]] <- tfdf
    }
  })
}
kegg_df <- bind_rows(kegg_list)

causal_snps_df <- deg_tf_rsids %>% select(SNP_id) %>%
  merge(., joined_df[, c("chr", "pos", "ref", "alt", "rsid")], by.x = "SNP_id", by.y = "rsid") %>%
  distinct()

write_tsv(causal_snps_df[, c('chr', 'pos')], '/media/leon/DISK2/icig/done/79_snps.tsv', col_names = F)

check_snps <- function(path) {
  ddf <- readr::read_tsv(path, show_col_types = FALSE)
  colnames(ddf) <- c("sample", "chr", "pos", "dp", "ref", "QS_ref", "alt", "QS_alt", "alt2", "QS_alt2")
  if (max(ddf$QS_alt2) < 0.1) { ddf <- ddf %>% select(-alt2, -QS_alt2) }
  ddf <- ddf %>%
    separate(col = sample, into = c("patient", "timepoint"), sep = "_") %>%
    filter(stringr::str_length(ref) == 1) %>%
    filter(stringr::str_length(alt) == 1) %>%
    mutate(dp_ref = round(dp * QS_ref), dp_alt = round(dp * QS_alt)) %>%
    select(-QS_ref, -QS_alt, -dp)
  return(ddf)
}

chipseq_tables <- lapply(Sys.glob("/media/leon/DISK2/icig/chipseq_july/stats_for_79_snps/*.stat"), check_snps)
chipseq_table <- bind_rows(chipseq_tables) %>%
  group_by(patient, chr, pos, ref, alt) %>%
  summarize(total_dp_ref = sum(dp_ref), total_dp_alt = sum(dp_alt), .groups = 'drop') %>%
  merge(., add_rsids(chipseq_table$chr, chipseq_table$pos, chipseq_table$ref),
        by.x = c("chr", "pos"), by.y = c("seqnames", "pos"))

chipseq_table_wtfs <- deg_tf_rsids %>%
  select(SNP_id, SYMBOL, TF) %>% distinct() %>%
  group_by(SNP_id, SYMBOL) %>% summarize(TFs = list(TF), .groups = 'drop') %>%
  merge(chipseq_table, ., by.x = "rsid", by.y = 'SNP_id')



patients_vv <- c("s3", "s4", "s5", "s6", "s12")
ase_in_vivo <- lapply(patients_vv, function(x) {
  ASE_in_vitro(x, min_DP = 100)
})

ase_vv <- bind_rows(ase_in_vivo) %>%
  mutate(padj = p.adjust(p01, method = "BH")) %>%
  mutate(diffASE = padj < 0.01)

ase_vv_subset <- ase_vv %>%
  filter(SYMBOL %in% chipseq_table_wtfs$SYMBOL)

merge(chipseq_table_wtfs[, c('rsid', 'patient', 'SYMBOL', 'TFs')], ase_vv[, c('SYMBOL', 'patient', 'diffASE')], by = c('patient', 'SYMBOL')) %>% group_by(SYMBOL) %>% filter(any(diffASE == T)) %>% distinct() %>% View()
write_tsv(TF_genes_df['TF'], '~/45_tfs.tsv', col_names = F)
  