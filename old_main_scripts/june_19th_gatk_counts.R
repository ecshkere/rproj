library(readr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)

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

add_rsids <- function(chr_vec, pos_vec, ref_vec) {
  snp_ranges <- GRanges(seqnames = as.character(chr_vec),
                        IRanges(pos_vec, width = 1), ref = ref_vec)
  matched_snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, snp_ranges)
  df_matches <- as.data.frame(matched_snps)[, c("seqnames", "pos")]
  df_matches$rsid <- mcols(matched_snps)$RefSNP_id
  df_matches <- na.omit(df_matches)
}
joined_df <- merge(joined_df, add_rsids(joined_df$chr, joined_df$pos, joined_df$ref), by.x = c("chr", "pos"), by.y = c("seqnames", "pos"))


###### in vitro RNA-seq
ASE_in_vitro <- function(patient, FDR = 0.01, min_DP = 50) {
  stats0R <- read_tsv(paste0("/media/leon/DISK2/icig/done/gatk_counts_filt/", patient, "_DR.stat"), show_col_types = FALSE) %>%
    dplyr::select(contig, position, refAllele, altAllele, refCount, altCount, totalCount)
  stats1R <- read_tsv(paste0("/media/leon/DISK2/icig/done/gatk_counts_filt/", patient, "_ER.stat"), show_col_types = FALSE) %>%
    dplyr::select(contig, position, refAllele, altAllele, refCount, altCount, totalCount)

  stats <- merge(stats0R, stats1R, by = c("contig", "position", "refAllele", "altAllele"), suffixes = c(".0R", ".1R"))
  stats$patient <- patient
  stats <- stats %>% relocate(patient, .after = position )
  colnames(stats) <- c("chr", "pos", "patient", "ref", "alt", "DP_ref_0R", "DP_alt_0R", "DP_0R", "DP_ref_1R", "DP_alt_1R", "DP_1R")
  stats <- stats %>% filter(DP_0R >= min_DP & DP_1R >= min_DP &
                              (DP_ref_0R / DP_0R > 0.1 & DP_ref_0R / DP_0R < 0.9) &
                              (DP_ref_1R / DP_1R > 0.1 & DP_ref_1R / DP_1R < 0.9))
  
  stats$chr <- as.character(as.numeric(stats$chr))
  
  stats$p01 <- apply(stats[, c("DP_ref_0R", "DP_alt_0R", "DP_ref_1R", "DP_alt_1R")], 1, function(x) {
    data_matrix <- matrix(unlist(c(x[1], x[2], x[3], x[4])), nrow = 2, byrow = TRUE)
    v <- chisq.test(data_matrix)
    return(v$p.value)
  })
  
  stats$p01adj <- p.adjust(stats$p01, method = "BH")
  stats <- stats[stats$p01adj < FDR, ]
  stats <- na.omit(stats)
  
  gr_positions <- GRanges(
    seqnames = paste0("chr", stats$chr),
    ranges = IRanges(start = stats$pos, end = stats$pos)
  )
  
  overlaps <- findOverlaps(gr_positions, exns)
  stats <- stats[queryHits(overlaps), ]
  snp_exons <- exns[subjectHits(overlaps)]
  gene_ids <- sapply(snp_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
  stats$gene_id <- gene_ids
  stats$SYMBOL <- mapIds(org.Hs.eg.db,
                         keys = stats$gene_id,
                         keytype = "ENTREZID",
                         column = "SYMBOL")
  stats <- stats %>% filter(!grepl("HLA-", SYMBOL))
  stats <- unique(stats)
  stats <- na.omit(stats)
  
  gene_names <- AnnotationDbi::select(org.Hs.eg.db, keys = as.character(stats$gene_id), columns = "GENENAME", keytype = "ENTREZID")
  
  stats <- merge(stats, gene_names, by.x = "gene_id", by.y = "ENTREZID")
  stats <- unique(stats)
  stats <- stats %>% mutate(log2FC1 = log2(DP_ref_1R * DP_alt_0R / DP_alt_1R / DP_ref_0R))  
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
    dplyr::select(-contains("id_sample_"), -contains("alt2"), -ref_0R, -ref_1R)
  
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
  ASE_in_vitro(x, FDR = 0.01, min_DP = 100)
})
ase <- bind_rows(ase_in_vitro)
# write_tsv(ase, "/media/leon/DISK2/icig/done/three_people_vtr_ase.tsv")
# ase <- read_tsv("/media/leon/DISK2/icig/done/three_people_vtr_ase.tsv")
ase$padj_all <- p.adjust(ase$p01, method = "BH")
ase$diffASEall <- ase$padj_all < FDR

# 2783 генов с diffASE при FDR = 0.01 и min_DP = 100
ase %>% filter(diffASEall == TRUE) %>% dplyr::select(SYMBOL) %>% unique() %>% nrow()

# всего 7254 гена с достаточным покрытием
ase %>% dplyr::select(SYMBOL) %>% unique() %>% nrow()

# из них 1468 с SNP в промоторах
ase <- ase %>%
  filter(SYMBOL %in% joined_df$SYMBOL)

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
  dplyr::select(SYMBOL, ENTREZID.RNA, patient.RNA, diffASEall, log2FC_1vs0, p01adj, chr.RNA, pos.RNA,
                ref.RNA, alt, DP_0R, QS_ref_0R, DP_1R, QS_ref_1R, rsid, patient.CS) %>%
  mutate(p01adj = round(p01adj, 2), log2FC_1vs0 = round(log2FC_1vs0, 2)) %>%
  # оставляем гены для которых хотя бы у одного человека есть diffASE
  group_by(SYMBOL) %>% filter(any(diffASEall == TRUE)) %>% ungroup()

# 865 генов таких чтобы было diffASE хотя бы у одного человека и совпадающие SNP в промоторах хотя бы у двух
length(unique(df$SYMBOL))
length(unique(df$rsid)) # 2855 SNP

result <- df %>%
  group_by(SYMBOL) %>% filter(n_distinct(diffASEall) >= 2) %>% ungroup() %>%
  group_by(across(-patient.CS)) %>%
  summarize(patient.CS = paste(patient.CS, collapse = ", "), .groups = "drop")

View(result)
length(unique(result$SYMBOL)) # 716 генов
length(unique(result$rsid)) # 2459 SNP











exns_df <- as.data.frame(exns)
exons_coords <- exns_df %>%
  filter(gene_id %in% df$ENTREZID.RNA) %>%
  dplyr::select(seqnames, start, end) %>%
  mutate(seqnames = substr(seqnames, 4, length(seqnames)))
write_tsv(exons_coords, "/media/leon/DISK2/icig/done/exons.tsv", col_names = FALSE)

write_tsv(df %>% group_by(rsid) %>% filter(!any(patient.RNA == "sBn")) %>% dplyr::select(chr.RNA, pos.RNA) %>% unique(),
          "/media/leon/DISK2/icig/done/rna_snps_BN.tsv",
          col_names = FALSE)

write_tsv(df %>% group_by(rsid) %>% filter(!any(patient.RNA == "sAn")) %>% dplyr::select(chr.RNA, pos.RNA) %>% unique(),
          "/media/leon/DISK2/icig/done/rna_snps_AN.tsv",
          col_names = FALSE)

write_tsv(df %>% group_by(rsid) %>% filter(!any(patient.RNA == "sDm")) %>% dplyr::select(chr.RNA, pos.RNA) %>% unique(),
          "/media/leon/DISK2/icig/done/rna_snps_DM.tsv",
          col_names = FALSE)


library(ASEP)

dmd <- read_tsv("/media/leon/DISK2/icig/done/alignments/in_vitro/sDm_DR_exons_ase_counts.tsv")
dmd$id <- "sDm"
dmd$group <- "D"
dmd <- merge(dmd, add_rsids(dmd$contig, dmd$position, dmd$refAllele), by.x = c("contig", "position"), by.y = c("seqnames", "pos"))
dmd <- dmd %>%
  rename(snp = rsid, total = totalCount, ref = refCount,
         chr = contig, pos = position) %>%
  filter(total > 50)

dmd <- dmd %>%
  dplyr::select(chr, pos, snp, refAllele, altAllele, ref, total, id, group)

gr_positions <- GRanges(
  seqnames = paste0("chr", dmd$chr),
  ranges = IRanges(start = dmd$pos, end = dmd$pos)
)
overlaps <- findOverlaps(gr_positions, exns)
dmd <- dmd[queryHits(overlaps), ]
snp_exons <- exns[subjectHits(overlaps)]
gene_ids <- sapply(snp_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
dmd$ENTREZID <- gene_ids
dmd <- na.omit(dmd)
dmd$gene <- mapIds(org.Hs.eg.db,
                       keys = dmd$ENTREZID,
                       keytype = "ENTREZID",
                       column = "SYMBOL")

dme <- read_tsv("/media/leon/DISK2/icig/done/alignments/in_vitro/sDm_ER_exons_ase_counts.tsv")
dme$id <- "sDm"
dme$group <- "E"
dme <- merge(dme, add_rsids(dme$contig, dme$position, dme$refAllele), by.x = c("contig", "position"), by.y = c("seqnames", "pos"))
dme <- dme %>%
  rename(snp = rsid, total = totalCount, ref = refCount,
         chr = contig, pos = position) %>%
  filter(total > 50)

dme <- dme %>%
  dplyr::select(chr, pos, snp, refAllele, altAllele, ref, total, id, group)

gr_positions <- GRanges(
  seqnames = paste0("chr", dme$chr),
  ranges = IRanges(start = dme$pos, end = dme$pos)
)
overlaps <- findOverlaps(gr_positions, exns)
dme <- dme[queryHits(overlaps), ]
snp_exons <- exns[subjectHits(overlaps)]
gene_ids <- sapply(snp_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
dme$ENTREZID <- gene_ids
dme <- na.omit(dme)
dme$gene <- mapIds(org.Hs.eg.db,
                     keys = dme$ENTREZID,
                     keytype = "ENTREZID",
                     column = "SYMBOL")


and <- read_tsv("/media/leon/DISK2/icig/done/alignments/in_vitro/sAn_DR_exons_ase_counts.tsv")
and$id <- "sDm"
and$group <- "D"
and <- merge(and, add_rsids(and$contig, and$position, and$refAllele), by.x = c("contig", "position"), by.y = c("seqnames", "pos"))
and <- and %>%
  rename(snp = rsid, total = totalCount, ref = refCount,
         chr = contig, pos = position) %>%
  filter(total > 50)

and <- and %>%
  dplyr::select(chr, pos, snp, refAllele, altAllele, ref, total, id, group)

gr_positions <- GRanges(
  seqnames = paste0("chr", and$chr),
  ranges = IRanges(start = and$pos, end = and$pos)
)
overlaps <- findOverlaps(gr_positions, exns)
and <- and[queryHits(overlaps), ]
snp_exons <- exns[subjectHits(overlaps)]
gene_ids <- sapply(snp_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
and$ENTREZID <- gene_ids
and <- na.omit(and)
and$gene <- mapIds(org.Hs.eg.db,
                   keys = and$ENTREZID,
                   keytype = "ENTREZID",
                   column = "SYMBOL")

ane <- read_tsv("/media/leon/DISK2/icig/done/alignments/in_vitro/sAn_ER_exons_ase_counts.tsv")
ane$id <- "sDm"
ane$group <- "E"
ane <- merge(ane, add_rsids(ane$contig, ane$position, ane$refAllele), by.x = c("contig", "position"), by.y = c("seqnames", "pos"))
ane <- ane %>%
  rename(snp = rsid, total = totalCount, ref = refCount,
         chr = contig, pos = position) %>%
  filter(total > 50)

ane <- ane %>%
  dplyr::select(chr, pos, snp, refAllele, altAllele, ref, total, id, group)

gr_positions <- GRanges(
  seqnames = paste0("chr", ane$chr),
  ranges = IRanges(start = ane$pos, end = ane$pos)
)
overlaps <- findOverlaps(gr_positions, exns)
ane <- ane[queryHits(overlaps), ]
snp_exons <- exns[subjectHits(overlaps)]
gene_ids <- sapply(snp_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
ane$ENTREZID <- gene_ids
ane <- na.omit(ane)
ane$gene <- mapIds(org.Hs.eg.db,
                   keys = ane$ENTREZID,
                   keytype = "ENTREZID",
                   column = "SYMBOL")

bnd <- read_tsv("/media/leon/DISK2/icig/done/alignments/in_vitro/sBn_DR_exons_ase_counts.tsv")
bnd$id <- "sDm"
bnd$group <- "D"
bnd <- merge(bnd, add_rsids(bnd$contig, bnd$position, bnd$refAllele), by.x = c("contig", "position"), by.y = c("seqnames", "pos"))
bnd <- bnd %>%
  rename(snp = rsid, total = totalCount, ref = refCount,
         chr = contig, pos = position) %>%
  filter(total > 50)

bnd <- bnd %>%
  dplyr::select(chr, pos, snp, refAllele, altAllele, ref, total, id, group)

gr_positions <- GRanges(
  seqnames = paste0("chr", bnd$chr),
  ranges = IRanges(start = bnd$pos, end = bnd$pos)
)
overlaps <- findOverlaps(gr_positions, exns)
bnd <- bnd[queryHits(overlaps), ]
snp_exons <- exns[subjectHits(overlaps)]
gene_ids <- sapply(snp_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
bnd$ENTREZID <- gene_ids
bnd <- na.omit(bnd)
bnd$gene <- mapIds(org.Hs.eg.db,
                   keys = bnd$ENTREZID,
                   keytype = "ENTREZID",
                   column = "SYMBOL")

bne <- read_tsv("/media/leon/DISK2/icig/done/alignments/in_vitro/sBn_ER_exons_ase_counts.tsv")
bne$id <- "sDm"
bne$group <- "E"
bne <- merge(bne, add_rsids(bne$contig, bne$position, bne$refAllele), by.x = c("contig", "position"), by.y = c("seqnames", "pos"))
bne <- bne %>%
  rename(snp = rsid, total = totalCount, ref = refCount,
         chr = contig, pos = position) %>%
  filter(total > 50)

bne <- bne %>%
  dplyr::select(chr, pos, snp, refAllele, altAllele, ref, total, id, group)

gr_positions <- GRanges(
  seqnames = paste0("chr", bne$chr),
  ranges = IRanges(start = bne$pos, end = bne$pos)
)
overlaps <- findOverlaps(gr_positions, exns)
bne <- bne[queryHits(overlaps), ]
snp_exons <- exns[subjectHits(overlaps)]
gene_ids <- sapply(snp_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
bne$ENTREZID <- gene_ids
bne <- na.omit(bne)
bne$gene <- mapIds(org.Hs.eg.db,
                   keys = bne$ENTREZID,
                   keytype = "ENTREZID",
                   column = "SYMBOL")


df_asep <- rbind(dmd, dme, and, ane, bnd, bne)
df_asep$ref_condition <- "D"
df_asep <- df_asep %>%
  # filter(gene %in% intersect(dmd$gene, dme$gene)) %>%
  filter(ref/total > 0.1 & ref/total < 0.9) %>%
  na.omit() %>% unique()
df_asep <- df_asep %>% group_by(snp) %>% filter(n_distinct(group) >= 2) %>% ungroup()
differential_ASE_detection(df_asep, adaptive=FALSE, n_resample=10, parallel=TRUE, n_core = 10, save_out=FALSE)

ASE_in_vitro_GATK <- function(patient, FDR = 0.01, min_DP = 100) {
  df0R <- read_tsv(paste0("/media/leon/DISK2/icig/done/gatk_counts_filt/", patient, "_DR.stat")) %>%
    dplyr::select(contig, position, refAllele, altAllele, refCount, altCount, totalCount)
  df1R <- read_tsv(paste0("/media/leon/DISK2/icig/done/gatk_counts_filt/", patient, "_ER.stat")) %>%
    dplyr::select(contig, position, refAllele, altAllele, refCount, altCount, totalCount)
  
  df <- merge(df0R, df1R, by = c("contig", "position", "refAllele", "altAllele"), suffixes = c(".0R", ".1R"))
  df$patient <- patient
  df <- df %>% relocate(patient, .after = position )
  colnames(df) <- c("chr", "pos", "patient", "ref", "alt", "DP_ref_0R", "DP_alt_0R", "DP_0R", "DP_ref_1R", "DP_alt_1R", "DP_1R")
  df <- df %>% filter(DP_0R >= min_DP & DP_1R >= min_DP &
                      (DP_ref_0R / DP_0R > 0.1 & DP_ref_0R / DP_0R < 0.9) &
                      (DP_ref_1R / DP_1R > 0.1 & DP_ref_1R / DP_1R < 0.9))
  
  df$chr <- as.character(as.numeric(df$chr))
  df <- na.omit(df)
  
  df$p01 <- apply(df[, c("DP_ref_0R", "DP_alt_0R", "DP_ref_1R", "DP_alt_1R")], 1, function(x) {
    data_matrix <- matrix(unlist(c(x[1], x[2], x[3], x[4])), nrow = 2, byrow = TRUE)
    v <- chisq.test(data_matrix)
    return(v$p.value)
  })

  df$p01adj <- p.adjust(df$p01, method = "BH")
  df$diffASE <- df$p01adj < FDR
  
  df$p01[is.na(df$p01)] <- "2"
  df$p01adj[is.na(df$p01adj)] <- "2"
  
  gr_positions <- GRanges(
    seqnames = paste0("chr", df$chr),
    ranges = IRanges(start = df$pos, end = df$pos)
  )
  
  overlaps <- findOverlaps(gr_positions, exns)
  df <- df[queryHits(overlaps), ]
  snp_exons <- exns[subjectHits(overlaps)]
  gene_ids <- sapply(snp_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
  df$ENTREZID <- gene_ids
  df$SYMBOL <- mapIds(org.Hs.eg.db,
                         keys = df$ENTREZID,
                         keytype = "ENTREZID",
                         column = "SYMBOL")
  df <- df %>% filter(!grepl("HLA-", SYMBOL))
  df <- unique(df)
  df <- na.omit(df)
  
  gene_names <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = as.character(df$ENTREZID),
                                      columns = "GENENAME",
                                      keytype = "ENTREZID")
  
  df <- merge(df, gene_names, by = "ENTREZID")
  df <- unique(df)
  df <- df %>% mutate(log2FC = log2(DP_ref_1R * DP_alt_0R / DP_alt_1R / DP_ref_0R))
  return(df)
}
dm_exons_ase <- ASE_in_vitro_GATK("sDm", FDR = 0.05, min_DP = 100)
bn_exons_ase <- ASE_in_vitro_GATK("sBn", FDR = 0.05, min_DP = 100)
an_exons_ase <- ASE_in_vitro_GATK("sAn", FDR = 0.05, min_DP = 100)
mean(dm_exons_ase$diffASE)
mean(bn_exons_ase$diffASE)
mean(an_exons_ase$diffASE)
sum(dm_exons_ase$diffASE)
sum(an_exons_ase$diffASE)
sum(bn_exons_ase$diffASE)

all_ase <- rbind(dm_exons_ase, bn_exons_ase, an_exons_ase)
all_ase$padj <- p.adjust(all_ase$p01, method = "BH")
all_ase$diffASE_all <- all_ase$padj < 0.05
sum(all_ase$diffASE_all)
all_ase %>% filter(patient == "sDm") %>% dplyr::select(diffASE_all) %>% sum()
all_ase %>% filter(patient == "sBn") %>% dplyr::select(diffASE_all) %>% sum()
all_ase %>% filter(patient == "sAn") %>% dplyr::select(diffASE_all) %>% sum()
