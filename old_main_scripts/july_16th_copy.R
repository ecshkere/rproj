library(readr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyr)

select <- dplyr::select

########### SNP из чипсеков
read_stats <- function(path) {
  # path <- "s15_0C.stat"; patient <- "s15"
  df <- readr::read_tsv(path, show_col_types = FALSE)
  colnames(df) <- c("sample", "chr", "pos", "dp", "ref", "QS_ref", "alt", "QS_alt", "alt2", "QS_alt2")
  df <- df %>% separate(col = sample, into = c("patient", "timepoint"), sep = "_")
  df$chr <- as.character(as.numeric(df$chr))
  df <- na.omit(df) 
  df <- df[stringr::str_length(df$ref) == 1, ]
  df <- df[stringr::str_length(df$alt) == 1, ]
  df$variant <- paste(df$chr, df$pos, sep = '_')
  return(df)
}

chipseq_tables <- lapply(Sys.glob("*C.stat"), read_stats)
joined_df <- bind_rows(chipseq_tables)

# выбираем SNP которые есть хотя бы у 2 из 3 людей
common_vars <- joined_df %>% dplyr::count(variant) %>% dplyr::filter(n > 1) %>% pull(variant) %>% unique()
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

ptnts <- unique(joined_df$patient)
ase <- ase_all %>%
  filter(patient %in% ptnts) %>%
  select(-c("DP_30R", "QS_ref_30R", "QS_alt1_30R", "DP_90R", "QS_ref_90R", "QS_alt1_90R", "p30", "p90", "log2FC_30vs0", "log2FC_90vs0", "padj30vs0", "padj90vs0", "ptnt_smbl")) %>%
  mutate(diffASE = padj1vs0 < 0.01)
         
    
# 5201 ген с diffASE при FDR = 0.01 и min_DP = 100
ase %>% filter(diffASE == TRUE) %>% select(SYMBOL) %>% unique() %>% nrow()

# всего 9503 гена с достаточным покрытием
ase %>% select(SYMBOL) %>% unique() %>% nrow()

# из них 3849 с diffASE и SNP в промоторах
ase <- ase %>%
  filter(SYMBOL %in% joined_df$SYMBOL)
ase %>% filter(diffASE == TRUE) %>% select(SYMBOL) %>% unique() %>% nrow()

# выбираем те SNP из чипсеков которые находятся в промоторах генов с SNP в экзонах
joined_df_filt <- joined_df %>%
  filter(SYMBOL %in% ase$SYMBOL)

# 6444 генов с SNP в промоторах (общих у 2+ людей) и SNP в экзонах
length(unique(joined_df_filt$SYMBOL))


#################### для выбранных генов с diffASE считаем ASE на уровне гена
set.seed(988482)
library(MBASED)

ase$variant <- paste(ase$chr, ase$pos, sep = '_')
ase <- ase %>%
  mutate(DP_ref_0R = round(DP_0R * QS_ref_0R),
         DP_alt_0R = round(DP_0R * QS_alt1_0R),
         DP_ref_1R = round(DP_1R * QS_ref_1R),
         DP_alt_1R = round(DP_1R * QS_alt1_1R))

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

mbased_processing <- function(patient) {
  df <- ase %>% filter(patient == patient)
  snps <- GRanges(seqnames = df$chr, ranges = IRanges(df$pos, width = 1), aseID = df$SYMBOL, allele1 = df$ref,allele2 = df$alt)
  names(snps) <- paste(df$SYMBOL, df$variant, sep = "_")
  df0R <- SummarizedExperiment(
    assays = list(lociAllele1Counts = matrix(df$DP_ref_0R, ncol = 1, dimnames = list(names(snps), paste0(patient, "0R"))),
                  lociAllele2Counts = matrix(df$DP_alt_0R, ncol = 1, dimnames = list(names(snps), paste0(patient, "0R")))),
    rowRanges = snps)
  df0R_mbased_results <- runMBASED(ASESummarizedExperiment = df0R, isPhased = FALSE, numSim = 10^3, BPPARAM = SerialParam())
  df0R_mbased <- summarizeASEResults(df0R_mbased_results)$geneOutput
  df1R <- SummarizedExperiment(
    assays = list(lociAllele1Counts = matrix(df$DP_ref_1R, ncol = 1, dimnames = list(names(snps), paste0(patient, "1R"))),
                  lociAllele2Counts = matrix(df$DP_alt_1R, ncol = 1, dimnames = list(names(snps), paste0(patient, "1R")))),
    rowRanges = snps)
  df1R_mbased_results <- runMBASED(ASESummarizedExperiment = df1R, isPhased = FALSE, numSim = 10^3, BPPARAM = SerialParam())
  df1R_mbased <- summarizeASEResults(df1R_mbased_results)$geneOutput
  
  df0R_mbased$patient <- patient; df0R_mbased$timepoint <- "0"; df0R_mbased$SYMBOL <- rownames(df0R_mbased)
  df1R_mbased$patient <- patient; df1R_mbased$timepoint <- "1"; df1R_mbased$SYMBOL <- rownames(df1R_mbased)

  df0gr <- unlist(summarizeASEResults(df0R_mbased_results)$locusOutput)
  df0df <- data.frame(snp_id = names(df0gr), as.data.frame(df0gr, row.names = NULL))
  df0df <- merge(df0df, df0R_mbased[, c("SYMBOL", "majorAlleleFrequency")], by.x = "aseID", by.y = "SYMBOL")
  df1gr <- unlist(summarizeASEResults(df1R_mbased_results)$locusOutput)
  df1df <- data.frame(snp_id = names(df1gr), as.data.frame(df1gr, row.names = NULL))
  df1df <- merge(df1df, df1R_mbased[, c("SYMBOL", "majorAlleleFrequency")], by.x = "aseID", by.y = "SYMBOL")
  df_full <- merge(df0df, df, by.x = c("seqnames", "start", "aseID"), by.y = c("chr", "pos", "SYMBOL"))
  df_full <- merge(df_full, df1df, by = c("seqnames", "start", "aseID"))
  df_full <- df_full %>% # отношение аллеля, преобладающего в точке 1, к альтернативному
    mutate(mbased_lfc = ifelse(allele1.x == allele1.y, 
                               log2((majorAlleleFrequency.y / (1 - majorAlleleFrequency.y)) / (majorAlleleFrequency.x / (1 - majorAlleleFrequency.x))),
                               log2((majorAlleleFrequency.y / (1 - majorAlleleFrequency.y)) / (1 - majorAlleleFrequency.x) / (majorAlleleFrequency.x))))
  print(paste(patient, "done"))
  return(df_full)
}

mbased_per_patient_list <- lapply(ptnts, mbased_processing)
full_mbased <- bind_rows(mbased_per_patient_list) %>%
  select(-c(snp_id.x, end.x, width.x, strand.x, snp_id.y, end.y, width.y, strand.y)) %>%
  filter(DP_ref_0R > 10 & DP_alt_0R > 10 & DP_ref_1R > 10 & DP_alt_0R > 10) %>% 
  mutate(absLFC = abs(mbased_lfc))

#####################################
# из всех SNP в экзонах одного гена выбираем один с наибольшим покрытием (?)
ase_filtered <- ase %>%
  mutate(total_DP = DP_0R + DP_1R) %>%
  group_by(SYMBOL, patient) %>%
  filter(total_DP == max(total_DP)) %>%
  ungroup()

# для каждого гена из RNA-seq находим какие SNP из чипсеков (общие у 2+ людей) есть в промоторе этого гена
df <- merge(ase_filtered, joined_df_filt, by = "SYMBOL", suffixes = c(".RNA", ".CS"), all.y = TRUE) %>%
  select(SYMBOL, patient.RNA, diffASE, log2FC_1vs0, padj1vs0, chr.RNA, pos.RNA,
         ref.RNA, alt.RNA, DP_0R, QS_ref_0R, DP_1R, QS_ref_1R, rsid, patient.CS) %>%
  # оставляем гены для которых хотя бы у одного человека есть diffASE
  group_by(SYMBOL) %>% filter(any(diffASE == TRUE)) %>% ungroup()

result <- df %>%
  # group_by(chr.RNA, pos.RNA, alt) %>%
  # mutate(diffdiffASE = all(diffASE == TRUE) & any(log2FC_1vs0 > 0) & any(log2FC_1vs0 < 0)) %>%
  # ungroup() %>%
  group_by(SYMBOL) %>%
  filter(n_distinct(diffASE) >= 2) %>% # есть один человек с diffASE и один без
  ungroup() %>%
  group_by(across(-patient.CS)) %>%
  summarize(patient.CS = paste(sort(unique(patient.CS)), collapse = ""), .groups = "drop")

View(result)
length(unique(result$SYMBOL)) # 356 генов
length(unique(result$rsid)) # 1237 SNP

write_tsv(unique(result['rsid']), "25.07.16.in_vivo_rsids.txt", col_names = FALSE)

########### motifbreakR
library(motifbreakR)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("/media/leon/DISK2/icig/done/motifbreakr_for_in_vivo/") ## threshold = 0.0005
result_files <- list.files(pattern = "motifbreak_.*\\.rds$")

all_results <- list()
for (file in result_files) {
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

tfs_w_10_snps <- snpgenecount[snpgenecount$SNP_count > 9, ][["geneSymbol"]]
mean(snpgenecount[snpgenecount$SNP_count > 9, ]$SNP_count) # в среднем 37 SNP на ТФ

snpgenedf <-  hocomoco %>%
  mcols() %>% 
  as.data.frame() %>% 
  select(geneSymbol, SNP_id, ALT) %>% 
  filter(geneSymbol %in% tfs_w_10_snps) %>%
  rename(TF = geneSymbol) %>%
  distinct()

snpgenedf %>% group_by(TF) %>%             
  summarize(SNP_count = n_distinct(SNP_id)) %>% View()

snpgenedf <- readRDS("/media/leon/DISK2/icig/done/motifbreakr_for_in_vivo/invivomtfbrkr_snpgenedf.rds")

result <- merge(result, full_mbased[, c("aseID", "patient", "absLFC")], by.x = c("SYMBOL", "patient.RNA"), by.y = c("aseID", "patient"))
merged_result <- merge(snpgenedf, result, by.x = c("SNP_id", "ALT"), by.y = c("rsid", "alt.RNA")) %>% distinct()
View(merged_result)

expanded_rows <- merged_result %>%
  filter(patient.CS == "ABD") %>% 
  dplyr::slice(rep(1:n(), each = 3)) %>%
  mutate(patient.CS = rep(c("AB", "BD", "AD"), times = n()/3))

merged_result <- merged_result %>%
  filter(patient.CS != "ABD") %>% 
  bind_rows(expanded_rows) %>%
  arrange(SNP_id, TF, patient.RNA)

rownames(merged_result) <- NULL

# логарифм отношения LFC для каждой пары людей по каждому гену с ASE
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
          logLFC = log2(LFCratio),
          LFC1 = LFC_1,
          LFC2 = LFC_2
        )
      })
  }) %>%
  ungroup() %>%  distinct()

lfcdiv <- lfc_divided %>%
  filter(aseID %in% intersect(lfc_divided$aseID, merged_result$SYMBOL))

dff <- merge(lfcdiv[, c("aseID", "patients", "logLFC", "LFCdiff", "LFC1", "LFC2")], merged_result, by.x = c("aseID"), by.y = c("SYMBOL"))
dffff <- dff %>%
  filter(patient.CS == patients & (patient.RNA == substr(patients, 1, 1) | patient.RNA == substr(patients, 2, 2)))

# выбираем 100 генов с наибольшим log(LFCratio)
topgenes <- dffff %>% select(aseID, logLFC) %>% arrange(desc(logLFC)) %>% distinct(aseID) %>% .[1:100, "aseID"]
dffff <- dffff %>%
  filter(aseID %in% topgenes)

length(unique(dffff$SNP_id)) # 259
length(unique(dffff$TF))     # 805

TFdct <- c("Z780A" = "ZNF780A", "THB" = "THRB", "ZSC22" = "ZSCAN22", "ZKSC8" = "ZKSCAN8", "Z585B" = "ZNF585B", "UBIP1" = "UBP1", "ZBT14" = "ZBTB14", "ANDR" = "AR",
           "ZSC31" = "ZSCAN31", "NKX25" = "NKX2-5", "TF2LY" = "TGIF2LY", "PKNX2" =  "PKNOX2", "UBIP1" = "UBP1", "MGAP" = "MGA", "HEN2" = "NHLH2", "HTF4" = "TCF12", "ITF2" = "TCF4")
dffff <- dffff %>%
  mutate(TF = ifelse(grepl("ZNF", TF), TF, gsub("^ZN", "ZNF", TF))) %>%
  mutate(TF = recode(TF, !!!TFdct))

tf_snp_count <- dffff %>%
  select(TF, SNP_id) %>%       
  group_by(TF) %>%             
  summarize(SNP_count = n_distinct(SNP_id)) %>%
  arrange(desc(SNP_count))
# tf_snp_count %>% View()

tf_gene_count <- dffff %>%
  select(TF, aseID) %>%       
  group_by(TF) %>%             
  summarize(gene_count = n_distinct(aseID)) %>%
  arrange(desc(gene_count))

# plot_df <- tf_snp_count[tf_snp_count$SNP_count > 1, ]
# plot(plot_df$SNP_count)
# axis(1, at = 1:nrow(plot_df), labels = plot_df$TF)

TFs <- tf_snp_count[tf_snp_count$SNP_count > 9, ][["TF"]] %>% unique()
length(TFs) # 148
tf_gene_count %>% filter(TF %in% TFs) %>% .[["gene_count"]] %>% mean() # 12.5
length(intersect(TFs, all_symbols)) # 136
setdiff(fixedTFs, all_symbols)

########## дифэкспрессия генов ТФ
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
  distinct()

deg_res_table %>% View()

TF_genes_df <- dffff %>% filter(TF %in% all_degs) %>% group_by(TF) %>% summarize(genes = list(aseID))

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