library(readr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyr)
library(dplyr)
set.seed(988482)
library(MBASED)
select <- dplyr::select

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
ASE_in_vitro <- function(patient, min_DP = 100) {
  stats0R <- read_tsv(paste0(patient, "_DR.stat"), show_col_types = FALSE) %>%
    filter(DP_total >= min_DP)
  stats1R <- read_tsv(paste0(patient, "_ER.stat"), show_col_types = FALSE) %>%
    filter(DP_total >= min_DP)
  
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

# 1565 генов с diffASE при FDR = 0.01 и min_DP = 100
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

full_mbased <- bind_rows(full_df_A_0vs1, full_df_A_1vs0, full_df_B_0vs1, full_df_B_1vs0, full_df_D_0vs1, full_df_D_1vs0) %>%
  select(-c(snp_id, end, width, strand, MAFDifference, QS_ref_0R, QS_alt1_0R, QS_ref_1R, QS_alt1_1R, DP_0R, DP_1R)) %>%
  rename(SYMBOL = aseID)

genewise_mbased <- full_mbased %>%
  select(SYMBOL, genewiseMAFDiff,  pValueASE, pValueHeterogeneity, ENTREZID, patient, comp) %>%
  distinct() %>%
  mutate(padj = p.adjust(pValueASE, method = "BH")) %>%
  mutate(diffASE = padj < 0.01) # ?

##################################
# для каждого гена из RNA-seq находим какие SNP из чипсеков (общие у 2+ людей) есть в промоторе этого гена
diffASE_w_CS <- merge(genewise_mbased, joined_df_filt, by = "SYMBOL", suffixes = c(".RNA", ".CS"), all.y = TRUE) %>%
  select(-c(ENTREZID.RNA, ENTREZID.CS, pValueASE, pValueHeterogeneity)) %>% 
  group_by(SYMBOL) %>%
  filter(n_distinct(diffASE) >= 2) %>% # есть один человек с diffASE и один без
  ungroup() # %>%
  # group_by(across(-patient.CS)) %>% # если один и тот же SNP в промоторе одного и того же гена у нескольких людей - склеиваем
  # summarize(patient.CS = paste(sort(patient.CS), collapse = ""), .groups = "drop") ####### зачем

# логарифм отношения LFC для каждой пары людей по каждому гену с ASE
lfc_divided <- diffASE_w_CS %>%
  select(SYMBOL, patient.RNA, genewiseMAFDiff) %>%
  mutate(genewiseMAFDiff = round(genewiseMAFDiff, 4)) %>%
  # mutate(genewiseMAFDiff = abs(genewiseMAFDiff)) %>%
  rename(absMAFDiff = genewiseMAFDiff) %>%
  distinct() %>% View()
  group_by(SYMBOL) %>% 
  filter(n() >= 2) %>%
  group_modify(~ {
    patients <- .x$patient
    LFCs <- .x$absLFC
    combn(length(patients), 2, simplify = FALSE) %>%
      purrr::map_dfr(function(pair) {
        i <- pair[1]
        j <- pair[2]
        
        patient_1 <- patients[i]
        patient_2 <- patients[j]
        LFC_1 <- LFCs[i]
        LFC_2 <- LFCs[j]
        
        # считаем LFCratio только если имеется одинаковый SNP в промоторе данного гена у обоих людей
        common_snps <- df_collapsed$rsid[df_collapsed$SYMBOL == .y$aseID & (df_collapsed$patients == paste0(sort(unique(patients)), collapse = "") | df_collapsed$patients == "ABD")]
        
        LFCratio <- if(length(common_snps) > 0) max(LFC_1, LFC_2)/min(LFC_1, LFC_2) else NA
        
        tibble(
          patients = ifelse(patient_1 < patient_2, paste0(patient_1, patient_2), paste0(patient_2, patient_1)),
          LFCratio = LFCratio,
          LFCdiff = abs(LFC_1 - LFC_2),
          logLFC = if(!is.na(LFCratio)) log2(LFCratio) else NA,
          LFC1 = LFC_1,
          LFC2 = LFC_2
        )
      })
  }) %>% 
  ungroup() %>% distinct()

lfcdiv <- lfc_divided %>%
  filter(aseID %in% intersect(lfc_divided$aseID, merged_result$SYMBOL)) %>%
  filter(!is.na(LFCratio))

dffff <- merge(lfcdiv[, c("aseID", "patients", "logLFC", "LFCdiff", "LFC1", "LFC2")], merged_result,
               by.x = c("aseID"), by.y = c("SYMBOL")) %>%filter(patient.CS == patients & (patient.RNA == substr(patients, 1, 1) | patient.RNA == substr(patients, 2, 2)))

# выбираем 100 генов с наибольшим log(LFCratio)
topgenes <- dffff %>% select(aseID, logLFC) %>% arrange(desc(logLFC)) %>% distinct(aseID) %>% .[1:100, "aseID"]
dffff <- dffff %>%
  filter(aseID %in% topgenes)

length(unique(result$SYMBOL)) # 806 генов
length(unique(result$rsid)) # 2958 SNP
length(unique(result[result$SYMBOL %in% topgenes,]$rsid)) # 427 SNP

write_tsv(unique(result['rsid']), "rsids_for_TFs.txt", col_names = FALSE)

########### motifbreakR
library(motifbreakR)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)

# setwd("/media/leon/DISK2/icig/done/motifbreak_results_less_pvalue//") ## threshold = 0.0005
# result_files <- list.files(pattern = "motifbreak_.*\\.rds$")
# 
# all_results <- list()
# for (file in result_files) {
#   res <- readRDS(file)
#   rsid <- sub("motifbreak_(.*)\\.rds", "\\1", file)
#   if (length(res) > 0) {
#     mcols(res)$rsid <- rsid
#     all_results[[rsid]] <- res
#   }
# }
# 
# combined_results <- do.call(c, unname(all_results))
hocomoco <- combined_results[combined_results$dataSource == "HOCOMOCOv13"]

snpgenedf <-  hocomoco %>%
  mcols() %>% 
  as.data.frame() %>% 
  select(geneSymbol, SNP_id, ALT) %>% 
  rename(TF = geneSymbol) %>%
  distinct()

result <- merge(result, full_mbased[, c("aseID", "patient", "absLFC")], by.x = c("SYMBOL", "patient.RNA"), by.y = c("aseID", "patient"))
merged_result <- merge(snpgenedf, result, by.x = c("SNP_id", "ALT"), by.y = c("rsid", "alt")) %>% distinct()

expanded_rows <- merged_result %>%
  filter(patient.CS == "ABD") %>% 
  dplyr::slice(rep(1:n(), each = 3)) %>%
  mutate(patient.CS = rep(c("AB", "BD", "AD"), times = n()/3))

merged_result <- merged_result %>%
  filter(patient.CS != "ABD") %>% 
  bind_rows(expanded_rows) %>%
  arrange(SNP_id, TF, patient.RNA)

rownames(merged_result) <- NULL





length(unique(dffff$SNP_id)) # 251
length(unique(dffff$TF))     # 887

TFdct <- c("Z780A" = "ZNF780A", "THB" = "THRB", "ZSC22" = "ZSCAN22", "ZKSC8" = "ZKSCAN8", "Z585B" = "ZNF585B", "UBIP1" = "UBP1", "ZBT14" = "ZBTB14", "ANDR" = "AR",
           "ZSC31" = "ZSCAN31", "NKX25" = "NKX2-5", "TF2LY" = "TGIF2LY", "PKNX2" =  "PKNOX2", "UBIP1" = "UBP1", "MGAP" = "MGA", "HEN2" = "NHLH2", "HTF4" = "TCF12", "ITF2" = "TCF4")
dffff <- dffff %>%
  mutate(TF = ifelse(grepl("ZNF", TF), TF, gsub("^ZN", "ZNF", TF))) %>%
  mutate(TF = recode(TF, !!!TFdct))

# tf_snp_count <- dffff %>%
#   select(TF, SNP_id) %>%       
#   group_by(TF) %>%             
#   summarize(SNP_count = n_distinct(SNP_id)) %>%
#   arrange(desc(SNP_count))
# 
# tf_gene_count <- dffff %>%
#   select(TF, aseID) %>%       
#   group_by(TF) %>%             
#   summarize(gene_count = n_distinct(aseID)) %>%
#   arrange(desc(gene_count))

# TFs <- tf_snp_count[tf_snp_count$SNP_count > 9, ][["TF"]] %>% unique()
TFs <- dffff$TF
length(intersect(TFs, all_symbols)) # 627 из 887
setdiff(TFs, all_symbols)

########## дифэкспрессия генов ТФ
TF_counts <- counts %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  filter(SYMBOL %in% TFs)

rownames(TF_counts) <- TF_counts$ENSEMBL
TF_counts <- TF_counts[, c(7, 8, 9, 10, 12, 14)]
colnames(TF_counts) <- gsub('R', '', colnames(TF_counts))
TF_counts <- TF_counts[, match(rownames(coldata), colnames(TF_counts))]
all(rownames(coldata) == colnames(TF_counts))
dds <- DESeqDataSetFromMatrix(countData = TF_counts,
                              colData = coldata,
                              design = ~ patient)
dds <- DESeq(dds)

resAB <- as.data.frame(results(dds, contrast = c("patient", "sAn", "sBn")))
resAB <- resAB %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "AB")

resAD <- as.data.frame(results(dds, contrast = c("patient", "sAn", "sDm")))
resAD <- resAD %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "AD")

resBD <- as.data.frame(results(dds, contrast = c("patient", "sBn", "sDm")))
resBD <- resBD %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "BD")

TF_deseq_results <- bind_rows(resAB, resAD, resBD) %>%
  filter(padj < 0.1 & abs(log2FoldChange) > log2(1.5))
TF_deseq_results
all_degs <- unique(TF_deseq_results$SYMBOL)

res_table <- dffff %>%
  filter(TF %in% TFs) %>%
  group_by(across(-TF)) %>% summarize(TF = list(TF), .groups = 'drop') %>%
  select(aseID, patients, logLFC, LFCdiff, LFC1, LFC2, SNP_id, TF) %>%
  mutate(TF = lapply(TF, sort),
         TF_number = lapply(TF, length),
         logLFC = round(logLFC, 2), LFC2 = round(LFC2, 2), LFCdiff = round(LFCdiff, 2), LFC1 = round(LFC1, 2)) %>%
  distinct()

deg_res_table <- dffff %>%
  filter(TF %in% all_degs) %>%
  group_by(across(-TF)) %>% summarize(TF = list(TF), .groups = 'drop') %>%
  select(aseID, patients, logLFC, LFCdiff, LFC1, LFC2, SNP_id, TF) %>%
  mutate(TF = lapply(TF, sort),
         TF_number = lapply(TF, length),
         logLFC = round(logLFC, 2), LFC2 = round(LFC2, 2), LFCdiff = round(LFCdiff, 2), LFC1 = round(LFC1, 2)) %>%
  distinct() %>% mutate(TF = sapply(TF, function(x) paste(x, collapse = ","))) %>% mutate(TF_number = sapply(TF_number, function(x) paste(x, collapse = ","))) 

TF_genes_df <- dffff %>% filter(TF %in% all_degs) %>% group_by(TF) %>% summarize(genes = list(aseID))
length(TF_genes_df$TF) # 140


library(clusterProfiler)
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
                                   universe = res$SYMBOL,
                                   readable = TRUE))
    if (nrow(tfdf) > 0) {
      tfdf$TF <- gen
      go_list[[gen]] <- tfdf
    }
  })
}
bind_rows(go_list) %>% View()

kegg_list <- list()
for (i in seq(nrow(TF_genes_df))) {
  gen <-  TF_genes_df[i, ][[1]]
  print(gen)
  gens <- mapIds(org.Hs.eg.db,
                 keys = unique(unlist(TF_genes_df[i, ]$genes)),
                 keytype = "SYMBOL",
                 column = "ENTREZID")
  try({
    tfdf <- as.data.frame(enrichKEGG(gene = gens,
                                     organism = "hsa",
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.1,
                                     universe = res$ENTREZID))
    if (nrow(tfdf) > 0) {
      tfdf$TF <- gen
      kegg_list[[gen]] <- tfdf
    }
  })
}
bind_rows(kegg_list) %>% View()


