library(readr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyr)

select <- dplyr::select

########### SNP из чипсеков
dm <- read_tsv("/media/leon/Polina/diabet/merge/Damarov.vcf")
an <- read_tsv("/media/leon/Polina/diabet/merge/E.Viacheslavovna.vcf")
bn <- read_tsv("/media/leon/Polina/diabet/merge/N.Petrovna.vcf")

dm$patient <- "D"; an$patient <- "A"; bn$patient <- "B"

dm$variant <- paste(dm$chr, dm$pos, sep = '_')
an$variant <- paste(an$chr, an$pos, sep = '_')
bn$variant <- paste(bn$chr, bn$pos, sep = '_')

all(all(colnames(dm) == colnames(an)), all(colnames(dm) == colnames(bn)), all(colnames(bn) == colnames(an)))
joined_df <- rbind(an, dm, bn)
joined_df$chr <- as.character(as.numeric(joined_df$chr))
joined_df <- na.omit(joined_df)

# выбираем SNP которые есть хотя бы у 2 из 3 людей
common_vars <- joined_df %>% dplyr::count(variant) %>% dplyr::filter(n > 1) %>% select(variant) %>% .[[1]] %>% unique()
joined_df <- joined_df %>% filter(variant %in% common_vars)

# пересекаем с промоторами
hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
prmtrs <- promoters(hg38, columns = c("gene_id"))
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

add_rsids <- function(chr_vec, pos_vec, ref_vec) {
  snp_ranges <- GRanges(seqnames = as.character(chr_vec),
                        IRanges(pos_vec, width = 1), ref = ref_vec)
  matched_snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, snp_ranges)
  df_matches <- as.data.frame(matched_snps)[, c("seqnames", "pos")]
  df_matches$rsid <- mcols(matched_snps)$RefSNP_id
  df_matches <- na.omit(df_matches)
}
joined_df <- merge(joined_df, add_rsids(joined_df$chr, joined_df$pos, joined_df$ref), by.x = c("chr", "pos"), by.y = c("seqnames", "pos"))


######### in vitro RNA-seq
ASE_in_vitro <- function(patient, min_DP = 100) {
  stats0R <- read_tsv(paste0("/media/leon/DISK2/icig/done/all_stats/", patient, "_DR.stat"), show_col_types = FALSE) %>%
    filter(DP_total >= min_DP)
  stats1R <- read_tsv(paste0("/media/leon/DISK2/icig/done/all_stats/", patient, "_ER.stat"), show_col_types = FALSE) %>%
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
    select(-contains("id_sample_"), -contains("alt2"), -ref_0R, -ref_1R) %>%
    mutate(log2FC_1vs0 = log2(QS_ref_1R * QS_alt1_0R / QS_alt1_1R / QS_ref_0R))
  return(stats)
}

patients_vtr <- c("sAn", "sBn", "sDm")
ase_in_vitro <- lapply(patients_vtr, function(x) {
  ASE_in_vitro(x, min_DP = 100)
})
ase <- bind_rows(ase_in_vitro)

ase$padj <- p.adjust(ase$p01, method = "BH")
ase$diffASE <- ase$padj < 0.01

# write_tsv(ase, "/media/leon/DISK2/icig/done/three_ppl_invtr_ase.tsv")
# ase <- read_tsv("/media/leon/DISK2/icig/done/three_ppl_invtr_ase.tsv")

# 1565 генов с diffASE при FDR = 0.01 и min_DP = 100
ase %>% filter(diffASE == TRUE) %>% select(SYMBOL) %>% unique() %>% nrow()

# всего 7253 гена с достаточным покрытием
ase %>% select(SYMBOL) %>% unique() %>% nrow()

# из них 843 с diffASE и SNP в промоторах
ase <- ase %>%
  filter(SYMBOL %in% joined_df$SYMBOL)
ase %>% filter(diffASE == TRUE) %>% select(SYMBOL) %>% unique() %>% nrow()

# выбираем те SNP из чипсеков которые находятся в промоторах генов с SNP в экзонах
joined_df_filt <- joined_df %>%
  filter(SYMBOL %in% ase$SYMBOL)

# 3377 генов с SNP в промоторах (общих у 2+ людей) и SNP в экзонах
length(unique(joined_df_filt$SYMBOL))

# из всех SNP в экзонах одного гена выбираем один с наибольшим покрытием (?)
ase_filtered <- ase %>%
  mutate(total_DP = DP_0R + DP_1R) %>%
  group_by(SYMBOL, patient) %>%
  filter(total_DP == max(total_DP)) %>%
  ungroup()

# для каждого гена из RNA-seq находим какие SNP из чипсеков (общие у 2+ людей) есть в промоторе этого гена
df <- merge(ase_filtered, joined_df_filt, by = "SYMBOL", suffixes = c(".RNA", ".CS"), all.y = TRUE) %>%
  select(SYMBOL, patient.RNA, diffASE, log2FC_1vs0, padj, chr.RNA, pos.RNA,
                ref.RNA, alt, DP_0R, QS_ref_0R, DP_1R, QS_ref_1R, rsid, patient.CS) %>%
  # оставляем гены для которых хотя бы у одного человека есть diffASE
  group_by(SYMBOL) %>% filter(any(diffASE == TRUE)) %>% ungroup()

result <- df %>%
  group_by(chr.RNA, pos.RNA, alt) %>%
  mutate(diffdiffASE = all(diffASE == TRUE) & any(log2FC_1vs0 > 0) & any(log2FC_1vs0 < 0)) %>%
  ungroup() %>%
  group_by(SYMBOL) %>%
  filter(
    n_distinct(diffASE) >= 2 | # есть один человек с diffASE и один без
    diffdiffASE                # либо diffASE у обоих но в разную сторону
  ) %>% ungroup() %>%
  group_by(across(-patient.CS)) %>%
  summarize(patient.CS = paste(sort(patient.CS), collapse = ""), .groups = "drop")

View(result)
length(unique(result$SYMBOL)) # 356 генов
length(unique(result$rsid)) # 1237 SNP

write_tsv(unique(result['rsid']), "rsids_for_TFs.txt", col_names = FALSE)


########### motifbreakR
library(motifbreakR)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("/media/leon/DISK2/icig/done/motifbreak_results_3000rsids_pvalue_10e4/") ## threshold = 0.0001
result_files <- list.files(pattern = "motifbreak_.*\\.rds$")
i <- 0
all_results <- list()
for (file in result_files) {
  i <- i + 1
  print(i/length(result_files)*100)
  res <- readRDS(file)
  rsid <- sub("motifbreak_(.*)\\.rds", "\\1", file)
  if (length(res) > 0) {
    mcols(res)$rsid <- rsid
    all_results[[rsid]] <- res
  }
}

combined_results <- do.call(c, unname(all_results))
hocomoco <- combined_results[combined_results$dataSource == "HOCOMOCOv13"]

snpgenecount <- hocomoco %>%
  mcols() %>% 
  as.data.frame() %>% 
  select(geneSymbol, SNP_id) %>% 
  distinct() %>% 
  group_by(geneSymbol) %>% 
  summarize(SNP_count = n_distinct(SNP_id)) %>% 
  arrange(desc(SNP_count))

tfs_w_10_snps <- snpgenecount[snpgenecount$SNP_count > 10, ][["geneSymbol"]]
mean(snpgenecount[snpgenecount$SNP_count > 10, ]$SNP_count) # 39 SNP на ТФ

snpgenedf <-  hocomoco %>%
  mcols() %>% 
  as.data.frame() %>% 
  select(geneSymbol, SNP_id, ALT) %>% 
  filter(geneSymbol %in% tfs_w_10_snps) %>%
  rename(TF = geneSymbol) %>%
  distinct()

snpgenedf %>% group_by(TF) %>%             
  summarize(SNP_count = n_distinct(SNP_id)) %>% View()

####### для выбранных генов с diffASE считаем ASE на уровне гена
set.seed(988482)
library(MBASED)

ase$variant <- paste(ase$chr, ase$pos, sep = '_')
ase <- ase %>%
  mutate(DP_ref_0R = round(DP_0R * QS_ref_0R),
         DP_alt_0R = round(DP_0R * QS_alt1_0R),
         DP_ref_1R = round(DP_1R * QS_ref_1R),
         DP_alt_1R = round(DP_1R * QS_alt1_1R))

aseDm <- ase %>% filter(patient == "D")
snpsDm <- GRanges(
  seqnames = aseDm$chr, ranges = IRanges(aseDm$pos, width = 1), aseID = aseDm$SYMBOL, allele1 = aseDm$ref,allele2 = aseDm$alt
)
names(snpsDm) <- paste(aseDm$SYMBOL, aseDm$variant, sep = "_")

dm0R <- SummarizedExperiment(assays = list(lociAllele1Counts = matrix(aseDm$DP_ref_0R, ncol = 1,dimnames = list(names(snpsDm), "dm0R")),
    lociAllele2Counts = matrix(aseDm$DP_alt_0R, ncol = 1, dimnames = list(names(snpsDm), "dm0R"))), rowRanges = snpsDm)

dm0R_mbased_results <- runMBASED(ASESummarizedExperiment = dm0R, isPhased = FALSE, numSim = 10^6, BPPARAM = SerialParam())

summarizeASEResults <- function(MBASEDOutput) {
  geneOutputDF <- data.frame(
    majorAlleleFrequency = assays(MBASEDOutput)$majorAlleleFrequency[, 1],
    pValueASE = assays(MBASEDOutput)$pValueASE[, 1],
    pValueHeterogeneity = assays(MBASEDOutput)$pValueHeterogeneity[, 1]
  )
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[, 1]
  lociOutputGR$MAF <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAF[, 1]
  lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels = unique(lociOutputGR$aseID)))
  return(
    list(
      geneOutput = geneOutputDF,
      locusOutput = lociOutputList
    )
  )
}
dm0R_mbased <- summarizeASEResults(dm0R_mbased_results)$geneOutput

dm1R <- SummarizedExperiment(assays = list(lociAllele1Counts = matrix(aseDm$DP_ref_1R, ncol = 1, dimnames = list(names(snpsDm), "dm1R")),
    lociAllele2Counts = matrix(aseDm$DP_alt_1R, ncol = 1, dimnames = list(names(snpsDm), "dm1R"))), rowRanges = snpsDm)
dm1R_mbased_results <- runMBASED(ASESummarizedExperiment = dm1R, isPhased = FALSE, numSim = 10^6, BPPARAM = SerialParam())
dm1R_mbased <- summarizeASEResults(dm1R_mbased_results)$geneOutput

aseAn <- ase %>% filter(patient == "A")
snpsAn <- GRanges(seqnames = aseAn$chr, ranges = IRanges(aseAn$pos, width = 1), aseID = aseAn$SYMBOL, allele1 = aseAn$ref, allele2 = aseAn$alt)
names(snpsAn) <- paste(aseAn$SYMBOL, aseAn$variant, sep = "_")

an0R <- SummarizedExperiment(assays = list(lociAllele1Counts = matrix(aseAn$DP_ref_0R, ncol = 1, dimnames = list(names(snpsAn), "an0R")),
    lociAllele2Counts = matrix(aseAn$DP_alt_0R, ncol = 1, dimnames = list(names(snpsAn), "an0R"))), rowRanges = snpsAn)
an0R_mbased_results <- runMBASED(ASESummarizedExperiment = an0R, isPhased = FALSE, numSim = 10^6, BPPARAM = SerialParam())
an0R_mbased <- summarizeASEResults(an0R_mbased_results)$geneOutput

an1R <- SummarizedExperiment(assays = list(lociAllele1Counts = matrix(aseAn$DP_ref_1R, ncol = 1, dimnames = list(names(snpsAn), "an1R")),
    lociAllele2Counts = matrix(aseAn$DP_alt_1R, ncol = 1, dimnames = list(names(snpsAn), "an1R"))), rowRanges = snpsAn)
an1R_mbased_results <- runMBASED(ASESummarizedExperiment = an1R, isPhased = FALSE, numSim = 10^6, BPPARAM = SerialParam())
an1R_mbased <- summarizeASEResults(an1R_mbased_results)$geneOutput

aseBn <- ase %>% filter(patient == "B")
snpsBn <- GRanges(seqnames = aseBn$chr, ranges = IRanges(aseBn$pos, width = 1), aseID = aseBn$SYMBOL, allele1 = aseBn$ref, allele2 = aseBn$alt)
names(snpsBn) <- paste(aseBn$SYMBOL, aseBn$variant, sep = "_")

bn0R <- SummarizedExperiment(assays = list(lociAllele1Counts = matrix(aseBn$DP_ref_0R, ncol = 1, dimnames = list(names(snpsBn), "bn0R")),
    lociAllele2Counts = matrix(aseBn$DP_alt_0R, ncol = 1, dimnames = list(names(snpsBn), "bn0R"))), rowRanges = snpsBn)
bn0R_mbased_results <- runMBASED(ASESummarizedExperiment = bn0R, isPhased = FALSE, numSim = 10^6, BPPARAM = SerialParam())
bn0R_mbased <- summarizeASEResults(bn0R_mbased_results)$geneOutput

bn1R <- SummarizedExperiment(assays = list(lociAllele1Counts = matrix(aseBn$DP_ref_1R, ncol = 1, dimnames = list(names(snpsBn), "bn1R")),
    lociAllele2Counts = matrix(aseBn$DP_alt_1R, ncol = 1, dimnames = list(names(snpsBn), "bn1R"))), rowRanges = snpsBn)
bn1R_mbased_results <- runMBASED(ASESummarizedExperiment = bn1R, isPhased = FALSE, numSim = 10^6, BPPARAM = SerialParam())
bn1R_mbased <- summarizeASEResults(bn1R_mbased_results)$geneOutput

dm0R_mbased$patient <- "D"; dm0R_mbased$timepoint <- "0"; dm0R_mbased$SYMBOL <- rownames(dm0R_mbased)
dm1R_mbased$patient <- "D"; dm1R_mbased$timepoint <- "1"; dm1R_mbased$SYMBOL <- rownames(dm1R_mbased)
an0R_mbased$patient <- "A"; an0R_mbased$timepoint <- "0"; an0R_mbased$SYMBOL <- rownames(an0R_mbased)
an1R_mbased$patient <- "A"; an1R_mbased$timepoint <- "1"; an1R_mbased$SYMBOL <- rownames(an1R_mbased)
bn0R_mbased$patient <- "B"; bn0R_mbased$timepoint <- "0"; bn0R_mbased$SYMBOL <- rownames(bn0R_mbased)
bn1R_mbased$patient <- "B"; bn1R_mbased$timepoint <- "1"; bn1R_mbased$SYMBOL <- rownames(bn1R_mbased)

dm0gr <- unlist(summarizeASEResults(dm0R_mbased_results)$locusOutput)
dm0df <- data.frame(snp_id = names(dm0gr), as.data.frame(dm0gr, row.names = NULL))
dm0df <- merge(dm0df, dm0R_mbased[, c("SYMBOL", "majorAlleleFrequency")], by.x = "aseID", by.y = "SYMBOL")
dm1gr <- unlist(summarizeASEResults(dm1R_mbased_results)$locusOutput)
dm1df <- data.frame(snp_id = names(dm1gr), as.data.frame(dm1gr, row.names = NULL))
dm1df <- merge(dm1df, dm1R_mbased[, c("SYMBOL", "majorAlleleFrequency")], by.x = "aseID", by.y = "SYMBOL")
dm_full <- merge(dm0df, aseDm, by.x = c("seqnames", "start", "aseID"), by.y = c("chr", "pos", "SYMBOL"))
dm_full <- merge(dm_full, dm1df, by = c("seqnames", "start", "aseID"))

# отношение аллеля, преобладающего в точке 1, к альтернативному
dm_full <- dm_full %>%
  mutate(mbased_lfc = ifelse(allele1IsMajor.x == allele1IsMajor.y,
                             log2((majorAlleleFrequency.y / (1 - majorAlleleFrequency.y)) / (majorAlleleFrequency.x / (1 - majorAlleleFrequency.x))),
                             log2((majorAlleleFrequency.y / (1 - majorAlleleFrequency.y)) / ((1 - majorAlleleFrequency.x) / majorAlleleFrequency.x))))

bn0gr <- unlist(summarizeASEResults(bn0R_mbased_results)$locusOutput)
bn0df <- data.frame(snp_id = names(bn0gr), as.data.frame(bn0gr, row.names = NULL))
bn0df <- merge(bn0df, bn0R_mbased[, c("SYMBOL", "majorAlleleFrequency")], by.x = "aseID", by.y = "SYMBOL")
bn1gr <- unlist(summarizeASEResults(bn1R_mbased_results)$locusOutput)
bn1df <- data.frame(snp_id = names(bn1gr), as.data.frame(bn1gr, row.names = NULL))
bn1df <- merge(bn1df, bn1R_mbased[, c("SYMBOL", "majorAlleleFrequency")], by.x = "aseID", by.y = "SYMBOL")
bn_full <- merge(bn0df, aseBn, by.x = c("seqnames", "start", "aseID"), by.y = c("chr", "pos", "SYMBOL"))
bn_full <- merge(bn_full, bn1df, by = c("seqnames", "start", "aseID"))
bn_full <- bn_full %>%
  mutate(mbased_lfc = ifelse(allele1IsMajor.x == allele1IsMajor.y,
                             log2((majorAlleleFrequency.y / (1 - majorAlleleFrequency.y)) / (majorAlleleFrequency.x / (1 - majorAlleleFrequency.x))),
                             log2((majorAlleleFrequency.y / (1 - majorAlleleFrequency.y)) / ((1 - majorAlleleFrequency.x) / majorAlleleFrequency.x))))

an0gr <- unlist(summarizeASEResults(an0R_mbased_results)$locusOutput)
an0df <- data.frame(snp_id = names(an0gr), as.data.frame(an0gr, row.names = NULL))
an0df <- merge(an0df, an0R_mbased[, c("SYMBOL", "majorAlleleFrequency")], by.x = "aseID", by.y = "SYMBOL")
an1gr <- unlist(summarizeASEResults(an1R_mbased_results)$locusOutput)
an1df <- data.frame(snp_id = names(an1gr), as.data.frame(an1gr, row.names = NULL))
an1df <- merge(an1df, an1R_mbased[, c("SYMBOL", "majorAlleleFrequency")], by.x = "aseID", by.y = "SYMBOL")
an_full <- merge(an0df, aseAn, by.x = c("seqnames", "start", "aseID"), by.y = c("chr", "pos", "SYMBOL"))
an_full <- merge(an_full, an1df, by = c("seqnames", "start", "aseID"))
an_full <- an_full %>%
  mutate(mbased_lfc = ifelse(allele1IsMajor.x == allele1IsMajor.y,
                             log2((majorAlleleFrequency.y / (1 - majorAlleleFrequency.y)) / (majorAlleleFrequency.x / (1 - majorAlleleFrequency.x))),
                             log2((majorAlleleFrequency.y / (1 - majorAlleleFrequency.y)) / ((1 - majorAlleleFrequency.x) / majorAlleleFrequency.x))))

full_mbased <- bind_rows(dm_full, an_full, bn_full) %>%
  select(-c(snp_id.x, end.x, width.x, strand.x, snp_id.y, end.y, width.y, strand.y)) %>%
  filter(DP_ref_0R > 10 & DP_alt_0R > 10 & DP_ref_1R > 10 & DP_alt_1R > 10) %>% 
  mutate(absLFC = abs(mbased_lfc))
  
result <- merge(result, full_mbased[, c("aseID", "patient", "absLFC")], by.x = c("SYMBOL", "patient.RNA"), by.y = c("aseID", "patient"))
merged_result <- merge(snpgenedf, result, by.x = c("SNP_id", "ALT"), by.y = c("rsid", "alt")) %>% distinct()
View(merged_result)

expanded_rows <- merged_result %>%
  filter(patient.CS == "ABD") %>% 
  slice(rep(1:n(), each = 3)) %>%
  mutate(patient.CS = rep(c("AB", "BD", "AD"), times = n()/3))

merged_result <- merged_result %>%
  filter(patient.CS != "ABD") %>% 
  bind_rows(expanded_rows) %>%
  arrange(SNP_id, TF, patient.RNA)

rownames(merged_result) <- NULL

lfc_divided <- full_mbased %>%
  select(aseID, patient, absLFC) %>% 
  distinct() %>%
  group_by(aseID) %>% 
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
          tibble(
            patients = ifelse(patient_1 < patient_2, paste0(patient_1, patient_2), paste0(patient_2, patient_1)),
            LFCratio = max(LFC_1, LFC_2) / min(LFC_1, LFC_2),
            LFCdiff = abs(LFC_1 - LFC_2),
            LFC_patient1 = LFC_1,
            LFC_patient2 = LFC_2
          )
        })
    }) %>%
  ungroup() %>%
  distinct()
  
lfcdiv <- lfc_divided %>%
  filter(aseID %in% intersect(lfc_divided$aseID, merged_result$SYMBOL)) %>%
  arrange(desc(LFCdiff))

dff <- merge(lfcdiv[, c("aseID", "patients", "LFCdiff")], merged_result, by.x = c("aseID"), by.y = c("SYMBOL"))

topgenes <- lfcdiv[1:100, ][["aseID"]]
dfff <- dff %>%
  filter(aseID %in% topgenes)

length(unique(dfff$SNP_id))
length(unique(dfff$TF))

tf_snp_count <- dfff %>%
  select(TF, SNP_id) %>%       
  group_by(TF) %>%             
  summarize(SNP_count = n_distinct(SNP_id)) %>%
  arrange(desc(SNP_count))

tf_snp_count %>% View()
plot_df <- tf_snp_count[tf_snp_count$SNP_count > 9, ]
# plot(plot_df$SNP_count)
# axis(1, at = 1:nrow(plot_df), labels = plot_df$TF)

TFs <- tf_snp_count[tf_snp_count$SNP_count > 9, ][["TF"]]
TFs <- ifelse(grepl("ZNF", TFs), TFs, gsub("ZN", "ZNF", TFs))
TFs[match(c("Z780A", "Z585B", "ZSC22", "ZSC31", "NKX25", "ZBT14", "THB"), TFs)] <- c("ZNF780A", "ZNF585B", "ZSCAN22", "ZSCAN31", "NKX2-5", "ZBTB14", "THRB")

length(intersect(TFs, counts$SYMBOL))
setdiff(TFs, counts$SYMBOL)

TF_counts <- counts %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  filter(SYMBOL %in% TFs)

rownames(TF_counts) <- TF_counts$ENSEMBL
TF_counts <- TF_counts[, 2:7]

dds <- DESeqDataSetFromMatrix(countData = TF_counts,
                              colData = coldata,
                              design = ~ patient)
dds <- DESeq(dds)

resAB <- as.data.frame(results(dds, contrast = c("patient", "sAn", "sBn")))
resAB <- resAB %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL")
resAB <- resAB %>%
  filter(padj < 0.1 & abs(log2FoldChange) > log2(1.5)) %>%
  mutate(contrast = "AB")

resAD <- as.data.frame(results(dds, contrast = c("patient", "sAn", "sDm")))
resAD <- resAD %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL")
resAD <- resAD %>%
  filter(padj < 0.1 & abs(log2FoldChange) > log2(1.5)) %>%
  mutate(contrast = "AD")

resBD <- as.data.frame(results(dds, contrast = c("patient", "sBn", "sDm")))
resBD <- resBD %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL")
resBD <- resBD %>%
  filter(padj < 0.1 & abs(log2FoldChange) > log2(1.5)) %>%
  mutate(contrast = "BD")

TF_deseq_results <- bind_rows(resAB, resAD, resBD)
TF_deseq_results
all_degs <- TF_deseq_results$SYMBOL

ABupreg <- (resAB %>% filter(log2FoldChange > 0))$SYMBOL
ABdownreg <- (resAB %>% filter(log2FoldChange < 0))$SYMBOL
ADupreg <- (resAD %>% filter(log2FoldChange > 0))$SYMBOL
ADdownreg <- (resAD %>% filter(log2FoldChange < 0))$SYMBOL
BDupreg <- (resBD %>% filter(log2FoldChange > 0))$SYMBOL
BDdownreg <- (resBD %>% filter(log2FoldChange < 0))$SYMBOL
                                                                        



