library(readr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

setwd("/media/leon/Polina/diabet/merge/")

# SNP из чипсеков
dm <- read_tsv("/media/leon/Polina/diabet/merge/Damarov.vcf")
an <- read_tsv("/media/leon/Polina/diabet/merge/E.Viacheslavovna.vcf")
bn <- read_tsv("/media/leon/Polina/diabet/merge/N.Petrovna.vcf")

dm$patient <- "sDm"; an$patient <- "sAn"; bn$patient <- "sBn"

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

add_rsids <- function(df) {
  snp_ranges <- GRanges(seqnames = as.character(df$chr),
                        IRanges(df$pos, width = 1), ref = df$ref)
  matched_snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, snp_ranges)
  df_matches <- as.data.frame(matched_snps)[, c("seqnames", "pos")]
  df_matches$rsid <- mcols(matched_snps)$RefSNP_id
  df_matches <- na.omit(df_matches)
  df <- merge(df, df_matches, by.x = c("chr", "pos"), by.y = c("seqnames", "pos"))
}
joined_df <- add_rsids(joined_df)

###### in vitro RNA-seq

ASE_in_vitro <- function(patient, FDR = 0.01, min_DP = 50) {
  # patient <- "sAn"
  stats0R <- read_delim(paste0("/media/leon/DISK2/icig/done/all_stats/", patient, "_DR.stat"), delim = "\t", trim_ws = TRUE) %>%
    filter(DP_total >= min_DP)
  stats1R <- read_delim(paste0("/media/leon/DISK2/icig/done/all_stats/", patient, "_ER.stat"), delim = "\t", trim_ws = TRUE) %>%
    filter(DP_total >= min_DP)
  
  stats <- merge(stats0R, stats1R, by = c("CHR", "POS"))
  colnames(stats) <- c("chr", "pos", "id_sample_0R", "DP_0R", "ref_0R", "QS_ref_0R", "alt1_0R", "QS_alt1_0R", "alt2_0R", "QS_alt2_0R", "id_sample_1R", "DP_1R", "ref_1R", "QS_ref_1R", "alt1_1R", "QS_alt1_1R", "alt2_1R", "QS_alt2_1R")
  stats$chr <- as.character(as.numeric(stats$chr))
  
  stats <- stats[str_length(stats$ref_0R) == 1,]
  stats <- stats[str_length(stats$alt1_0R) == 1,]
  stats <- stats[str_length(stats$alt1_1R) == 1,]
  
  stats[stats$alt1_0R != stats$alt1_1R, c(7, 8, 9, 10, 15, 16, 17, 18)] <- stats[stats$alt1_0R != stats$alt1_1R, c(7, 8, 9, 10, 17, 18, 15, 16)]
  
  stats$p01 <-  apply(stats[, c(4, 6, 8, 12, 14, 16)], 1, function(x) {
    data_matrix <- matrix(unlist(c(round(x[1] * x[2]), round(x[1] * x[3]), round(x[4] * x[5]), round(x[4] * x[6]))), nrow = 2, byrow = TRUE)
    v <- chisq.test(data_matrix)
    return(v$p.value)
  })
  
  stats$p01adj <- p.adjust(stats$p01, method = "BH")
  # stats <- stats[stats$p01adj < FDR, ]
  stats$diffASE <- stats$p01adj < FDR
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
  stats <- unique(stats)
  stats <- na.omit(stats)
  stats$SYMBOL <- mapIds(org.Hs.eg.db,
                         keys = stats$ENTREZID,
                         keytype = "ENTREZID",
                         column = "SYMBOL")
  stats <- stats %>% filter(!grepl("HLA-", SYMBOL))
  
  gene_names <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = as.character(stats$ENTREZID),
                                      columns = "GENENAME",
                                      keytype = "ENTREZID")
  
  stats <- merge(stats, gene_names, by = "ENTREZID")
  stats <- unique(stats)
  stats <- stats %>%
    mutate(DP_ref_0R = round(DP_0R * QS_ref_0R),
           DP_alt_0R = round(DP_0R * QS_alt1_0R),
           DP_ref_1R = round(DP_1R * QS_ref_1R),
           DP_alt_1R = round(DP_1R * QS_alt1_1R),
           patient = gsub("_.*", "", gsub("RNASEQ_", "", id_sample_0R)),
           ref = ref_0R, alt = alt1_0R) %>%
    dplyr::select(-contains("id_sample_"), -contains("alt2"), -contains("QS"), -contains("alt1"), -ref_0R, -ref_1R)
  
  stats$ase0pval <- apply(stats[, c("DP_ref_0R", "DP_alt_0R")], 1, function(x) {
    binom.test(c(x[1], x[2]), 0.5)$p.value
  })
  stats$ase0padj <- p.adjust(stats$ase0pval, method = "BH")
  stats$ase0 <- stats$ase0padj < FDR
  
  stats$ase1pval <- apply(stats[, c("DP_ref_1R", "DP_alt_1R")], 1, function(x) {
    binom.test(c(x[1], x[2]), 0.5)$p.value
  })
  stats$ase1padj <- p.adjust(stats$ase1pval, method = "BH")
  stats$ase1 <- stats$ase1padj < FDR
  
  stats <- stats %>% mutate(
    log2FC_0 = log2(DP_ref_0R / DP_alt_0R),
    log2FC_1 = log2(DP_ref_1R / DP_alt_1R),
    log2FC_1vs0 = log2(DP_ref_1R * DP_alt_0R / DP_alt_1R / DP_ref_0R),
  )
  return(stats)
}

patients_vtr <- c("sAn", "sBn", "sDm")
ase_in_vitro <- lapply(patients_vtr, function(x) {
  ASE_in_vitro(x, FDR = 0.01)
})
ase <- bind_rows(ase_in_vitro)
# write_csv(ase, "/media/leon/DISK2/icig/done/three_people_vtr_ase.csv")
# ase <- read_csv("/media/leon/DISK2/icig/done/three_people_vtr_ase.csv")

# 2501 ген с diffASE
ase %>% filter(diffASE == TRUE) %>% select(SYMBOL) %>% nrow()
ase <- ase %>%
  filter(SYMBOL %in% joined_df$SYMBOL) # 1394 с SNP в промоторах

# выбираем те SNP из чипсеков которые находятся в промоторах генов с diffASE 
joined_df_filt <- joined_df %>%
  filter(SYMBOL %in% ase$SYMBOL)

# 888 генов с diffASE и SNP в промоторах общих у 2+ людей
length(unique(joined_df_filt$SYMBOL))

# # ase <- ase %>%
# #   mutate(DP_ref_0R = round(DP_0R * QS_ref_0R),
# #          DP_alt_0R = round(DP_0R * QS_alt1_0R),
# #          DP_ref_1R = round(DP_1R * QS_ref_1R),
# #          DP_alt_1R = round(DP_1R * QS_alt1_1R),
# #          ref = ref_0R, alt = alt1_0R) %>%
# #   rename(ENTREZID = gene_id) %>%
# #   dplyr::select(-contains("id_sample_"), -contains("alt2"), -contains("QS"), -contains("alt1"), -ref_0R, -ref_1R)
# # 
# # ase <- ase[, c("chr", "pos", "ref", "alt", "DP_0R", "DP_ref_0R", "DP_alt_0R", "DP_1R", "DP_ref_1R", "DP_alt_1R", "ENTREZID", "SYMBOL", "GENENAME", "patient", "p01", "p01adj", "log2FC1")]

# из всех SNP в экзонах одного гена выбираем один с наибольшим покрытием (?)
ase_filtered <- ase %>%
  mutate(total_DP = DP_0R + DP_1R) %>%
  group_by(SYMBOL) %>%
  filter(total_DP == max(total_DP)) %>%
  ungroup()

# для каждого гена с diffASE находим какие SNP из чипсеков (общие у 2+ людей) есть в промоторе этого гена
df <- merge(ase_filtered, joined_df_filt, by = c("SYMBOL", "patient"), suffixes = c(".RNA", ".CS"), all.y = TRUE)
length(unique(df[df$diffASE == TRUE, ]$rsid)) # 600 SNP в чипсеках общих у 2+ людей в промоторах генов с diffASE

dff <- df %>% dplyr::select(SYMBOL, patient, diffASE, log2FC_1vs0, chr.RNA, pos.RNA, ref.RNA, alt, DP_0R, DP_1R, rsid)
dfff <- dff %>% group_by(SYMBOL) %>% filter(any(diffASE == TRUE)) %>% na.omit()
length(unique(dfff$rsid))

# 1 ген такой чтобы у двух людей было diffASE и одинаковые SNP в промоторах
rsids <- dfff %>%  dplyr::count(rsid) %>% dplyr::filter(n > 1) %>% dplyr::select(rsid) %>% .[[1]] %>% unique()
length(rsids)
fd <- dfff %>% filter(rsid %in% rsids) %>% na.omit() 




set.seed(988482)
library(MBASED)

aseDm <- ase %>% filter(patient == "sDm")
snpsDm <- GRanges(
  seqnames = aseDm$chr,
  ranges = IRanges(aseDm$pos, width = 1),
  aseID = aseDm$SYMBOL,
  allele1 = aseDm$ref,
  allele2 = aseDm$alt
)
names(snpsDm) <- paste(aseDm$SYMBOL, aseDm$rsid, sep = "_")

dm0R <- SummarizedExperiment(
  assays = list(
    lociAllele1Counts = matrix(
      aseDm$DP_ref_0R,
      ncol = 1,
      dimnames = list(names(snpsDm), "dm0R")
    ),
    lociAllele2Counts = matrix(
      aseDm$DP_alt_0R,
      ncol = 1,
      dimnames = list(names(snpsDm), "dm0R")
    )
  ),
  rowRanges = snpsDm
)

dm0R_mbased_results <- runMBASED(
  ASESummarizedExperiment = dm0R,
  isPhased = FALSE,
  numSim = 10^6,
  BPPARAM = SerialParam()
)

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
dm0R_mbased <- summarizeASEResults(dm0R_mbased_results)

dm1R <- SummarizedExperiment(
  assays = list(
    lociAllele1Counts = matrix(
      aseDm$DP_ref_1R,
      ncol = 1,
      dimnames = list(names(snpsDm), "dm1R")
    ),
    lociAllele2Counts = matrix(
      aseDm$DP_alt_1R,
      ncol = 1,
      dimnames = list(names(snpsDm), "dm1R")
    )
  ),
  rowRanges = snpsDm
)

dm1R_mbased_results <- runMBASED(
  ASESummarizedExperiment = dm1R,
  isPhased = FALSE,
  numSim = 10^6,
  BPPARAM = SerialParam()
)
dm1R_mbased <- summarizeASEResults(dm1R_mbased_results)

dm0R_mbased$geneOutput

counts <- read.table("/media/leon/DISK2/icig/done/in_vitro_results/counts/counts.tsv", header = TRUE, row.names = 1, sep = "\t")
counts <- counts[, 6:length(counts)]
colnames(counts) <- gsub("\\.bam", "", colnames(counts))
colnames(counts) <- gsub("RNASEQ_", "", colnames(counts))
colnames(counts) <- gsub("DR", "0R", colnames(counts))
colnames(counts) <- gsub("ER", "1R", colnames(counts))
counts <- counts %>%
  dplyr::select(-contains("Kh"))
counts$ENSEMBL <- rownames(counts)

annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = rownames(counts),
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "ENSEMBL")

annotations_ordered <- annotations[match(rownames(counts), annotations$ENSEMBL), ]
annotations_clean <- annotations[!is.na(annotations$ENTREZID) & !is.na(annotations$SYMBOL), ]
counts <- merge(counts, annotations_ordered, by = "ENSEMBL")
countsDm <- counts %>% filter(SYMBOL %in% aseDm$SYMBOL)

gfhfsf <- aseDm %>% dplyr::count(SYMBOL) %>% dplyr::filter(n > 1) %>% dplyr::select(SYMBOL) %>% .[[1]]
