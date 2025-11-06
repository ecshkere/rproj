library(readr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library()
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
common_vars <- joined_df %>% dplyr::count(variant) %>% dplyr::filter(n > 1) %>% dplyr::select(variant) %>% .[[1]] %>% unique()
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
    dplyr::select(-contains("id_sample_"), -contains("alt2"), -ref_0R, -ref_1R) %>%
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
ase %>% filter(diffASE == TRUE) %>% dplyr::select(SYMBOL) %>% unique() %>% nrow()

# всего 7253 гена с достаточным покрытием
ase %>% dplyr::select(SYMBOL) %>% unique() %>% nrow()

# из них 843 с diffASE и SNP в промоторах
ase <- ase %>%
  filter(SYMBOL %in% joined_df$SYMBOL)
ase %>% filter(diffASE == TRUE) %>% dplyr::select(SYMBOL) %>% unique() %>% nrow()

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
  dplyr::select(SYMBOL, patient.RNA, diffASE, log2FC_1vs0, padj, chr.RNA, pos.RNA,
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

## threshold = 0.0005
# setwd("/media/leon/DISK2/icig/done/motifbreak_results_less_pvalue//")
# result_files <- list.files(pattern = "motifbreak_.*\\.rds$")
# 
# all_results <- list()
# 
# for (file in result_files) {
#   res <- readRDS(file)
#   rsid <- sub("motifbreak_(.*)\\.rds", "\\1", file)
#   
#   if (length(res) > 0) {
#     mcols(res)$rsid <- rsid
#     all_results[[rsid]] <- res
#   }
# }
# 
# combined_results <- do.call(c, unname(all_results))
hocomoco <- combined_results[combined_results$dataSource == "HOCOMOCOv13"]

snpgenecount <- hocomoco %>%
  mcols() %>% 
  as.data.frame() %>% 
  dplyr::select(geneSymbol, SNP_id) %>% 
  distinct() %>% 
  group_by(geneSymbol) %>% 
  summarise(SNP_count = n_distinct(SNP_id)) %>% 
  arrange(desc(SNP_count))

tfs_w_10_snps <- snpgenecount[snpgenecount$SNP_count > 10,][["geneSymbol"]]
mean(snpgenecount[snpgenecount$SNP_count > 10,]$SNP_count) # 39 SNP на ТФ

snpgenedf <-  hocomoco %>%
  mcols() %>% 
  as.data.frame() %>% 
  dplyr::select(geneSymbol, SNP_id, ALT) %>% 
  filter(geneSymbol %in% tfs_w_10_snps) %>%
  rename(TF = geneSymbol) %>%
  distinct()

snpgenedf %>% group_by(TF) %>%             
  summarise(SNP_count = n_distinct(SNP_id)) %>% View()

merged_result <- merge(snpgenedf, result, by.x = c("SNP_id", "ALT"), by.y = c("rsid", "alt")) %>% unique()
View(merged_result)

lfc_divided <-
  merged_result %>%
  group_by(chr.RNA, pos.RNA, ref.RNA, ALT, SNP_id, TF, SYMBOL, patient.CS) %>%
  filter(n() >= 2) %>%
    group_modify(~ {
      patients <- .x$patient.RNA
      log2FCs <- .x$log2FC_1vs0
      combn(length(patients), 2, simplify = FALSE) %>%
        purrr::map_dfr(function(pair) {
          i <- pair[1]
          j <- pair[2]
          if (log2FCs[i] >= log2FCs[j]) {
            patient_high <- patients[i]
            patient_low <- patients[j]
            log2FC_high <- log2FCs[i]
            log2FC_low <- log2FCs[j]
          } else {
            patient_high <- patients[j]
            patient_low <- patients[i]
            log2FC_high <- log2FCs[j]
            log2FC_low <- log2FCs[i]
          }
          tibble(
            patient1 = patient_high,
            patient2 = patient_low,
            LFCratio = abs(log2FC_high) / abs(log2FC_low),
            LFCdiff = abs(log2FC_high - log2FC_low),
            log2FC_patient1 = log2FC_high,
            log2FC_patient2 = log2FC_low
          )
        })
    }) %>%
  ungroup() %>%
  unique()
  
dffff <- lfc_divided %>%
  group_by(TF) %>%
  mutate(TF_CV = mean(LFCdiff) / sd(LFCdiff) * 100) %>%
  select(TF, TF_CV) %>% unique() %>%
  arrange(desc(TF_CV))


toptfs <- unique(dffff[1:100, ]$TF)
merged_result %>%
  filter(TF %in% toptfs) %>%  
  select(TF, SNP_id) %>%       
  group_by(TF) %>%             
  summarise(SNP_count = n_distinct(SNP_id)) %>%
  arrange(desc(SNP_count)) %>% View()

# snpMart <- useEnsembl(biomart = "snps", dataset = "hsapiens_snp")
# SNP_MAF <- getBM(attributes = c('refsnp_id', 'minor_allele', 'minor_allele_freq'),
#                  filters = 'snp_filter', values = unique(merged_result$SNP_id), mart = snpMart)
# 
# stat_shot <- merge(stat_shot, SNP_MAF, by.x = 'rsid', by.y = 'refsnp_id')


cvs <- merged_result %>%
  group_by(geneSymbol, patient.CS) %>%
  summarise(cv = sd(abs(log2FC_1vs0)) / mean(abs(log2FC_1vs0)) * 100)

merged_result <- merge(merged_result, cvs, by = "geneSymbol")


strong_results <- combined_results[combined_results$effect == "strong"]
cleaned_rownames <- gsub("\\..*", "", names(strong_results))

AD_snps <- result %>% filter(patient.CS == "AD") %>% dplyr::select(rsid) %>% unique()
AD_results <- strong_results[cleaned_rownames %in% AD_snps[[1]]]
AD_genes <- unique(AD_results[AD_results$dataSource == "HOCOMOCOv13"]$geneSymbol)
length(AD_genes)

BD_snps <- result %>% filter(patient.CS == "BD") %>% dplyr::select(rsid) %>% unique()
BD_results <- strong_results[cleaned_rownames %in% BD_snps[[1]]]
BD_genes <- unique(BD_results[BD_results$dataSource == "HOCOMOCOv13"]$geneSymbol)
length(BD_genes)

AB_snps <- result %>% filter(patient.CS == "AB") %>% dplyr::select(rsid) %>% unique()
AB_results <- strong_results[cleaned_rownames %in% AB_snps[[1]]]
AB_genes <- unique(AB_results[AB_results$dataSource == "HOCOMOCOv13"]$geneSymbol)
length(AB_genes)

ABD_snps <- result %>% filter(patient.CS == "ABD") %>% dplyr::select(rsid) %>% unique()
ABD_results <- strong_results[cleaned_rownames %in% ABD_snps[[1]]]
ABD_genes <- unique(ABD_results[ABD_results$dataSource == "HOCOMOCOv13"]$geneSymbol)
length(ABD_genes)

# sig_results <- combined_results[abs(combined_results$alleleEffectSize) > 0.25]
# tf_genes <- sig_results$geneSymbol %>% unique()
# length(tf_genes)
# sigs <- na.omit(combined_results)[elementMetadata(na.omit(combined_results))$Refpvalue < 0.05]

disruption_counts <- as.data.frame(table(strong_results$rsid))
colnames(disruption_counts) <- c("rsid", "disruption_count")
disruption_counts <- disruption_counts[order(-disruption_counts$disruption_count), ]

tf_counts <- as.data.frame(table(strong_results$geneSymbol))
tf_counts <- tf_counts[order(-tf_counts$Freq), ]
head(tf_counts, 20)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Make row names unique
names(strong_results) <- make.names(names(strong_results), unique = TRUE)

# Now run annotation
annotated <- annotatePeak(
  strong_results,
  tssRegion = c(-2000, 2000),
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  annoDb = "org.Hs.eg.db"
)

# Convert to data frame
annotated_df <- as.data.frame(annotated)

# Filter for known disease-associated genes
library(gwasrapidd)
disease_genes <- get_traits()#$trait  # Get disease-associated genes
prioritized <- annotated_df[annotated_df$SYMBOL %in% disease_genes, ]

# TF disruption landscape
library(ggplot2)
ggplot(as.data.frame(tf_counts[1:20, ]), 
       aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Transcription Factor", y = "Disruption Count")

snps <- snps.from.rsid(rsid = "rs4759021",
                       dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                       search.genome = BSgenome.Hsapiens.UCSC.hg38)

results <- motifbreakR(snpList = snps,
                       filterp = TRUE,              # Filter by p-value
                       pwmList = MotifDb,           # Default motif database
                       threshold = 0.01,            # P-value cutoff
                       method = "ic",               # Use "information content"
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25), # Uniform background
                       # BPPARAM = BiocParallel::bpparam()
) # Parallel processing

# Filter for significant disruptions (p < 1e-4)
significant_results <- results[elementMetadata(results)$pvalue < 0.01]
# 
# # Save results to CSV
# # write.csv(as.data.frame(significant_results), "motifbreakR_results.csv")
# 
# # Plot top 10 most disruptive SNPs
plotMB(results = significant_results, rsid = names(significant_results)[1:10], effect = "strong")