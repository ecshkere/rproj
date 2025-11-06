library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(org.Hs.eg.db)
# library(AnnotationDbi)
library(clusterProfiler)
library(motifbreakR)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(jsonlite)
library(DESeq2)
library(biomaRt)
library(GenomicRanges)
library(rtracklayer)
library(vcfR)
library(tidyverse)

source("scripts/functions.R")

## in vitro chip-seq snps
dm <- read_tsv("data/Damarov.vcf")
an <- read_tsv("data/E.Viacheslavovna.vcf")
bn <- read_tsv("data/N.Petrovna.vcf")

dm$patient <- "D"; an$patient <- "A"; bn$patient <- "B"

all(all(colnames(dm) == colnames(an)), all(colnames(dm) == colnames(bn)))

joined_df <- bind_rows(an, dm, bn) %>% 
  mutate(chr = as.character(as.numeric(chr))) %>% # removing X Y MT and alternative contigs
  drop_na() %>%
  mutate(alt = ifelse(allel1 != ref, allel1, allel2)) %>%
  select(-c(DP, QS1, QS2, allel1, allel2)) %>% 
  mutate(variant = paste(chr, pos, ref, alt, sep = '_')) %>% 
  group_by(variant) %>% filter(n_distinct(patient) > 1) %>% ungroup() # choosing snps present in at least 2 out of 3 people

# choosing snps in promoters
joined_df <- assign_genes(joined_df, "chr", "pos", promoters = TRUE)
joined_df <- add_rsids(joined_df, "chr", "pos", "ref") 

stats_dir <- "data/"

#### in vitro RNA-seq ------------------------------------------------------------------------
read_stats <- function(patient, min_DP = 100, gatk = FALSE) {
  if (gatk) {
    stats_dir <- "/media/leon/DISK2/icig/done/gatk_counts_filt/"
    if (!is.na(as.numeric(substr(patient, 2, 3)))) {
      stats0R <- read_tsv(paste0(stats_dir, patient, "_0R.stat")) 
      stats1R <- read_tsv(paste0(stats_dir, patient, "_1R.stat"))
    } else {
      stats0R <- read_tsv(paste0(stats_dir, patient, "_DR.stat"))
      stats1R <- read_tsv(paste0(stats_dir, patient, "_ER.stat"))
    }
    stats0R <- stats0R %>% 
      select(contig, position, refAllele, altAllele, refCount, altCount, totalCount) %>% 
      rename(chr = contig, pos = position, ref = refAllele, alt = altAllele, DP_ref = refCount, DP_alt = altCount, DP = totalCount)
    
    stats1R <- stats1R %>% 
      select(contig, position, refAllele, altAllele, refCount, altCount, totalCount) %>% 
      rename(chr = contig, pos = position, ref = refAllele, alt = altAllele, DP_ref = refCount, DP_alt = altCount, DP = totalCount)
    
    stats <- inner_join(stats0R, stats1R, by = c("chr", "pos", "ref", "alt"), suffix = c("_0R", "_1R")) %>% 
      filter(DP_ref_0R > 5, DP_alt_0R > 5, DP_ref_1R > 5, DP_alt_1R > 5) %>% 
      filter(chr %in% as.character(seq(1, 22)))
    
    stats <- assign_genes(stats, "chr", "pos", exons_only = TRUE) %>%
      filter(!grepl("HLA-", symbol)) %>%
      mutate(patient = patient) %>%
      filter(DP_ref_0R / DP_0R > 0.1 & DP_alt_0R / DP_0R > 0.1 & DP_ref_1R / DP_1R > 0.1 & DP_alt_1R / DP_1R > 0.1) %>%
      mutate(log2FC_1vs0 = log2(DP_ref_1R * DP_alt_0R / DP_alt_1R / DP_ref_0R),
             ref_frac_diff = DP_ref_1R / DP_1R - DP_ref_0R / DP_0R)
    return(stats)
    }
  if (!is.na(as.numeric(substr(patient, 2, 3)))) {
    stats0R <- read_tsv(paste0(stats_dir, patient, "_0R.stat"))
    stats1R <- read_tsv(paste0(stats_dir, patient, "_1R.stat"))
  } else {
    stats0R <- read_tsv(paste0(stats_dir, patient, "_DR.stat"))
    stats1R <- read_tsv(paste0(stats_dir, patient, "_ER.stat"))
  }
  stats <- inner_join(stats0R, stats1R, by = c("chr", "pos"), suffix = c("_0R", "_1R")) %>% 
    filter(DP_0R >= min_DP, DP_1R >= min_DP) %>%
    # filter(DP_ref_0R > 5, DP_alt_0R > 5, DP_ref_1R > 5, DP_alt_1R > 5) %>% 
    filter(chr %in% as.character(seq(1, 22))) %>% 
    filter(str_length(ref_0R) == 1, str_length(alt1_0R) == 1, str_length(alt1_1R) == 1) 
   
  stats[stats$alt1_0R == stats$alt2_1R, c("alt1_1R", "QS_alt1_1R", "alt2_1R", "QS_alt2_1R")] <- stats[stats$alt1_0R == stats$alt2_1R, c("alt2_1R", "QS_alt2_1R", "alt1_1R", "QS_alt1_1R")]  
  
  stats <- assign_genes(stats, "chr", "pos", exons_only = TRUE) %>%
    filter(!grepl("HLA-", symbol)) %>%
    mutate(patient = gsub("._.*", "", gsub("RNASEQ_s", "", id_sample_0R)),
           ref = ref_0R, alt = alt1_0R,
           DP_ref_0R = round(DP_0R * QS_ref_0R),
           DP_alt_0R = round(DP_0R * QS_alt1_0R),
           DP_ref_1R = round(DP_1R * QS_ref_1R),
           DP_alt_1R = round(DP_1R * QS_alt1_1R)) %>%
    filter(DP_ref_0R / DP_0R > 0.1 & DP_alt_0R / DP_0R > 0.1 & DP_ref_1R / DP_1R > 0.1 & DP_alt_1R / DP_1R > 0.1) %>%
    mutate(log2FC_1vs0 = log2(QS_ref_1R * QS_alt1_0R / QS_alt1_1R / QS_ref_0R),
           ref_frac_diff = QS_ref_1R - QS_ref_0R) %>% 
    select(-c(contains("id_sample_"), contains("alt2"), ref_0R, ref_1R, alt1_0R, alt1_1R, contains("QS")))
  
  if (!is.na(as.numeric(substr(patient, 2, 3)))) {
    stats$patient <- patient
  }
  return(stats)
}

patients_vtr <- c("sAn", "sBn", "sDm")
ase_in_vitro <- lapply(patients_vtr, function(x) {
  read_stats(x, min_DP = 50)
})

gatk_ase <- lapply(patients_vtr, function(x) {
  read_stats(x, min_DP = 100, gatk = TRUE)
})

pthrshld <- 0.05
gatkk <- bind_rows(gatk_ase) %>% 
  mutate(patient = substr(patient, 2, 2)) %>%
  filter(symbol %in% joined_df$symbol) %>%
  mutate(total_DP = DP_0R + DP_1R) %>% 
  group_by(symbol, chr, pos) %>%  # choosing the one most common snp for each gene across patients
  mutate(n = n(), total_across_dp = sum(total_DP), patients = paste(sort(unique(patient)), collapse = "")) %>% 
  filter(n > 1) %>% ungroup() %>% 
  rowwise() %>% # calculating chi-square and adjusting p-value
  mutate(p01 = chisq.test(matrix(c(DP_ref_0R, DP_alt_0R, DP_ref_1R, DP_alt_1R), nrow = 2, byrow = TRUE))$p.value) %>%
  ungroup() %>% drop_na() %>% 
  mutate(padj = p.adjust(p01, method = "BH")) %>% 
  mutate(diffASE = padj < pthrshld) %>%
  mutate(patients = if_else(patients == "ABD", "AB,BD,AD", patients)) %>% # про каждый снп записываем у каких попарных комбинаций людей он есть
  separate_rows(patients, sep = ",") %>% 
  ungroup() %>% distinct()

ase <- bind_rows(ase_in_vitro) %>% 
  filter(symbol %in% joined_df$symbol) %>%
  mutate(total_DP = DP_0R + DP_1R) %>% 
  group_by(symbol, chr, pos) %>%  # choosing the one most common snp for each gene across patients
  mutate(n = n(), total_across_dp = sum(total_DP), patients = paste(sort(unique(patient)), collapse = "")) %>% 
  filter(n > 1) %>% ungroup() %>% 
  rowwise() %>% # calculating chi-square and adjusting p-value
  mutate(p01 = chisq.test(matrix(c(DP_ref_0R, DP_alt_0R, DP_ref_1R, DP_alt_1R), nrow = 2, byrow = TRUE))$p.value) %>%
  ungroup() %>% drop_na() %>% 
  mutate(padj = p.adjust(p01, method = "BH")) %>% 
  mutate(diffASE = padj < pthrshld) %>%
  mutate(patients = if_else(patients == "ABD", "AB,BD,AD", patients)) %>% # про каждый снп записываем у каких попарных комбинаций людей он есть
  separate_rows(patients, sep = ",") %>% 
  ungroup() %>% distinct()

# ase <- gatkk
# для каждого гена из RNA-seq находим какие SNP из чипсеков у каких попарных комбинаций людей есть в промоторе этого гена
joined_df_filt <- joined_df %>% filter(symbol %in% ase$symbol)

df_collapsed <- joined_df_filt %>%
  group_by(across(-patient)) %>%
  summarise(patients = paste(sort(unique(patient)), collapse = ""), .groups = "drop") %>% 
  mutate(patients = if_else(patients == "ABD", "AB,BD,AD", patients)) %>%
  separate_rows(patients, sep = ",") %>% distinct()

# остались только те строки где у двух людей есть одинаковые снп и в экзоне и в промоторе
ddmaf_df <- inner_join(ase, distinct(df_collapsed[, c("symbol", "patients")])) %>%
  group_by(symbol) %>%
  filter(any(diffASE == TRUE) & n_distinct(patient) > 1) %>%
  filter(n == max(n)) %>%
  filter(total_across_dp == max(total_across_dp)) %>%
  select(symbol, chr, pos, patient, ref_frac_diff, padj, n, patients) %>%
  distinct() %>%
  group_by(symbol) %>%
  group_modify(~ {
    patients <- .x$patient
    comparisons <- .x$patients
    ref_frac_diffs <- .x$ref_frac_diff
    padjes <- .x$padj
    combn(length(patients), 2, simplify = FALSE) %>%
      map_dfr(function(pair) {
        i <- pair[1]; j <- pair[2]
        p1 <- patients[i]; p2 <- patients[j]
        drf1 <- ref_frac_diffs[i]; drf2 <- ref_frac_diffs[j]
        pvl1 <- padjes[i]; pvl2 <- padjes[j]
        ptnts <- ifelse(p1 < p2, paste0(p1, p2), paste0(p2, p1))

        # есть хотя бы один общий снп + хотя бы у одного достоверно диффасе => считаем модуль разницы изменения доли референсного аллеля
        dd_ref_frac <- if (ptnts %in% comparisons & min(pvl1, pvl2) < pthrshld) abs(drf1 - drf2) else NA
    
        tibble(
          patients = ptnts,
          dd_ref_frac = dd_ref_frac,
          drf1 = drf1, drf2 = drf2,
          padj1 = pvl1, padj2 = pvl2
        )})}) %>%
  ungroup() %>%
  distinct() %>%
  filter(!is.na(dd_ref_frac)) %>%
  arrange(desc(dd_ref_frac))

# выбираем 100 генов с наибольшим модулем разницы изменения доли референсного аллеля
topgenes <- unique(ddmaf_df %>% pull(symbol))[1:100]
length(unique(ddmaf_df$symbol)) # 101
summary(ddmaf_df$dd_ref_frac)
length(unique(ddmaf_df[1:110, ]$symbol)) # 101
ddmaf_df[[109, "dd_ref_frac"]] # 0.201  # 0.143773
ddmaf_df_filt <- ddmaf_df[1:109, ] %>% 
  left_join(df_collapsed %>% select(symbol, patients, rsid))

snps_for_mtfbrkr <- unique(ddmaf_df_filt$rsid)
length(snps_for_mtfbrkr) # 272

######## motifbreakR --------------------------------------------
mtfbrkr_dir <- "motifbreakr"
result_files <- list.files(mtfbrkr_dir, pattern = ".*\\.rds$") 
rsids_in_folder <- gsub("motifbreak_", "", gsub("\\.rds", "", result_files))
rsids_to_process <- setdiff(snps_for_mtfbrkr, rsids_in_folder)

## writeLines(rsids_to_process, "/media/leon/DISK2/icig/done/snps_for_mtfbrkr_to_add.txt")
if (length(rsids_to_process)) {
  batch_size <- 32
  batches <- split(rsids_to_process, ceiling(seq_along(rsids_to_process) / batch_size))
  
  counter <- 0
  for (batch in batches) {
    snp_batch <- snps.from.rsid(rsid = batch,
                                dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                search.genome = BSgenome.Hsapiens.UCSC.hg38)
    
    res_batch <- motifbreakR(snpList = snp_batch,
                             pwmList = MotifDb,
                             filterp = TRUE,
                             threshold = 0.0001,
                             method = "ic")
    
    for (rsid in batch) {
      output_file <- file.path(mtfbrkr_dir, paste0(rsid, ".rds"))
      snp_res <- res_batch[res_batch$SNP_id == rsid]
      if (length(snp_res) > 0) {
        saveRDS(snp_res, output_file)
      }
    }
    counter <- counter + 1
    cat("batch", counter, "out of", length(batches), "done\n")
  }
}
result_files <- list.files(mtfbrkr_dir, pattern = ".*\\.rds$") 
rsids_in_folder <- gsub(".rds", "", result_files)
result_files <- result_files[rsids_in_folder %in% snps_for_mtfbrkr]

all_results <- list()
for (file in result_files) {
  res <- readRDS(paste(mtfbrkr_dir, file, sep = '/'))
  rsid <- sub("(.*)\\.rds", "\\1", file)
  if (!rsid %in% snps_for_mtfbrkr) { next }
  if (length(res) > 0) {
    mcols(res)$rsid <- rsid
    all_results[[rsid]] <- res
  }
}
combined_results <- do.call(c, unname(all_results))
hcmc <- combined_results[combined_results$dataSource == "HOCOMOCOv13"]

snpgenedf <- hcmc %>%
  mcols() %>% as.data.frame() %>% 
  select(geneSymbol, SNP_id) %>% 
  rename(TF = geneSymbol) %>%
  distinct()

##### дифэкспрессия ТФ in vitro -------------------------------------------------------------------
counts <- read.table("counts_in_vitro_17.09.tsv", header = TRUE, row.names = 1, sep = "\t")
counts <- counts[, c(6:11)]
colnames(counts) <- c("D0", "D1", "B0", "B1", "A0", "A1")
coldata <- data.frame(patient = substr(colnames(counts), 1, 1),
                      time = substr(colnames(counts), 2, 2)) %>%
  mutate(patient = as.factor(patient), time = as.factor(time))
rownames(coldata) <- colnames(counts)

# список всех генов которые экспрессируются вообще (background для GO)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized = TRUE) >= 5 ) >= 2
dds <- dds[keep, ]
dds <- DESeq(dds)
res_all <- results(dds)

annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = rownames(res_all),
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "ENSEMBL") %>% 
  rename(entrezid = ENTREZID, symbol = SYMBOL, ensembl = ENSEMBL)

annotations_ordered <- annotations[match(rownames(res_all), annotations$ensembl), ]
annotations_clean <- annotations_ordered %>%
  filter(!is.na(entrezid) & !is.na(symbol))
all_symbols <- unique(annotations_clean$symbol)

TFs <- unique(snpgenedf$TF)
length(unique(snpgenedf$SNP_id)) # 263
length(TFs) # 658
length(intersect(TFs, all_symbols)) # 318

snpgenedf <- snpgenedf %>%
  mutate(TF = recode(TF, !!!tf_to_gene))

TFs <- unique(snpgenedf$TF)
length(TFs) # 652
length(intersect(TFs, all_symbols)) # 545

TF_counts <- counts %>%
  mutate(ensembl = rownames(.)) %>%
  inner_join(annotations_clean) %>%
  filter(symbol %in% TFs)

rownames(TF_counts) <- TF_counts$ensembl
TF_counts <- TF_counts %>%
  select(-ensembl, -entrezid, -symbol)
dds <- DESeqDataSetFromMatrix(countData = TF_counts,
                              colData = coldata,
                              design = ~ patient) # time + patient?
dds <- DESeq(dds)

### выбираем попарно те дэги которые связаны с общими снп в промоторах генов по которым у этих двух людей разница в диффасе
diffase_snp_tf <- ddmaf_df_filt %>%
  select(symbol, patients, rsid) %>%
  rename(SNP_id = rsid) %>%
  inner_join(snpgenedf)

resAB <- as.data.frame(results(dds, contrast = c("patient", "A", "B"))) %>%
  mutate(ensembl = rownames(.)) %>%
  inner_join(annotations_clean) %>%
  mutate(contrast = "AB") %>%
  filter(symbol %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "AB"])

resAD <- as.data.frame(results(dds, contrast = c("patient", "A", "D"))) %>%
  mutate(ensembl = rownames(.)) %>%
  inner_join(annotations_clean) %>%
  mutate(contrast = "AD") %>%
  filter(symbol %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "AD"])

resBD <- as.data.frame(results(dds, contrast = c("patient", "B", "D"))) %>%
  mutate(ensembl = rownames(.)) %>%
  inner_join(annotations_clean) %>%
  mutate(contrast = "BD") %>%
  filter(symbol %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "BD"])

TF_deseq_results <- bind_rows(resAB, resAD, resBD) %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 0.5) 
all_degs <- sort(unique(TF_deseq_results$symbol))
length(all_degs) # 118
length(intersect(read_tsv("output/newer_list_of_degtfs.txt", col_names = F)$X1, all_degs))
setdiff(read_tsv("output/newer_list_of_degtfs.txt", col_names = F)$X1, all_degs)
setdiff(all_degs, read_tsv("output/newer_list_of_degtfs.txt", col_names = F)$X1)

# writeLines(all_degs, "output/newerer_list_of_degtfs.txt")

deg_tf_rsids <- ddmaf_df_filt %>% 
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id") %>%
  filter(TF %in% all_degs)

TF_genes_df <- deg_tf_rsids %>% group_by(TF) %>% summarize(genes = list(symbol))


#################### MISSENSES
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
attrs <- c("external_gene_name", "chromosome_name", "exon_chrom_start", "exon_chrom_end")

exons_df <- getBM(attributes = attrs,
                  filters = "external_gene_name",
                  values = TFs,
                  mart = ensembl) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name),
         start = exon_chrom_start,
         end = exon_chrom_end) %>%
 select(chromosome_name, exon_chrom_start, exon_chrom_end)

write.table(exons_df, "~/in_vivo/exons_for_all_tfs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## bedtools intersect -a /media/leon/Polina/atac_rna/dbSnp153.bed -b ~/exons_for_all_tfs.bed | awk '{print $1, $2, $3, $4}' - > ~/tf_exon_rsids_153.bed

snps_gr <- read_table("~/all_tfs_exon_rsids_153.bed", col_names = F)
the_rsids <- unique(snps_gr$X4)
length(the_rsids) # 1245398

snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
snps_cons <- readRDS("/media/leon/DISK2/icig/biomart_third.RDS")
the_rsids <- setdiff(the_rsids, snps_cons$refsnp_id)
for (i in seq(1, length(the_rsids), by = 10000)) {
  snps_cons <- rbind(snps_cons,
                     getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "consequence_type_tv"),
                           filters = 'snp_filter',
                           values = the_rsids[i:(i+10000)],
                           mart = snp_mart)
  )
  print(i)
}
# saveRDS(snps_cons, "/media/leon/DISK2/icig/biomart4.RDS")

snps_missense <- snps_cons %>% 
  filter(refsnp_id %in% the_rsids) %>% 
  filter(consequence_type_tv == "missense_variant")

missense_ranges <- GRanges(
  seqnames = paste0("chr", snps_missense$chr_name),
  ranges = IRanges(start = snps_missense$chrom_start,
                   end = snps_missense$chrom_start),
  rsid = snps_missense$refsnp_id,
  consequence = snps_missense$consequence_type_tv
)

mssns_ovrlps <- findOverlaps(missense_ranges, exns)
mssns_snps <- snps_missense[queryHits(mssns_ovrlps), ]
mssns_exons <- exns[subjectHits(mssns_ovrlps)]
mssns_genes <- sapply(mssns_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
mssns_snps$entrezid <- mssns_genes
mssns_snps$symbol <- mapIds(org.Hs.eg.db,
                            keys = mssns_snps$entrezid,
                            column = "SYMBOL",
                            keytype = "ENTREZID")

mssns_snps <- mssns_snps %>% 
  filter(symbol %in% TFs) %>%  select(-entrezid) %>%  distinct() %>%
  rename(TF = symbol)

mssns_snps$TF <- unlist(mssns_snps$TF)

mssns_snps_bed <- mssns_snps %>%
  mutate(end = chrom_start) %>%
  mutate(chrom_start = chrom_start - 1) %>%
  select(chr_name, chrom_start, end) %>% distinct() 

write.table(mssns_snps_bed, "~/13.10.2025.5.23.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# awk 'NF==3 && $1!="NA"' ~/in_vitro_missense_snps.bed > ~/in_vitro_missense_snps.clean.bed
# cd /media/leon/DISK2/icig/done/alignments/in_vitro/
# bcftools mpileup -R  ~/in_vitro_missense_snps.clean.bed -f /media/leon/Polina/Genomes/hg38.chromFa/hg38_FOR_HISAT.fa RNASEQ_sBn.merged.bam RNASEQ_sAn.merged.bam RNASEQ_sDm.merged.bam -d 1000 | bcftools call -mv -Oz -o final_in_vitro_missense_genotypes.vcf
# bcftools filter -i 'MIN(DP)>=20' final_in_vitro_missense_genotypes.vcf -Oz -o final_in_vitro_missense_genotypes.dp20.vcf; mv final_in_vitro_missense_genotypes.dp20.vcf final_in_vitro_missense_genotypes.vcf

translate_gt <- function(gt, ref, alt) {
  case_when(
    gt == "0/0" ~ paste0(ref, ref),
    gt == "0/1" ~ paste0(ref, alt),
    gt == "1/1" ~ paste0(alt, alt),
    gt == "./." ~ NA_character_,
    TRUE ~ NA_character_
  )
}

add_comparisons <- function(patients_vec) {
  comparison_map <- combn(patients_vec, 2, simplify = FALSE) %>%
    map_df(~ tibble(
      patients = paste0(sort(c(.x[1], .x[2])), collapse = ''),
      patient1 = .x[1],
      patient2 = .x[2]
    ))
  return(comparison_map)
}

vcf <- read.vcfR("/media/leon/DISK2/icig/done/alignments/s4_bn_an_dm_filtered.vcf.gz")

refalt <- as.data.frame(vcf@fix) %>%
  select(CHROM, POS, REF, ALT) %>%
  rename(chr = CHROM, pos = POS) %>%
  separate_rows(ALT, sep = ",")

gt <- as.data.frame(extract.gt(vcf)) %>%
  rownames_to_column("variant") %>%
  separate(variant, into = c("chr", "pos"), sep = "_", remove = TRUE) %>%
  distinct() %>%
  merge(mssns_snps, by.x = c("chr", "pos"), by.y = c("chr_name", "chrom_start")) %>%
  left_join(refalt, by = c("chr", "pos"), relationship = "many-to-many") %>%
  filter(str_length(REF) == 1) %>%
  filter(str_length(ALT) == 1) %>%
  rename(B = sBn, A = sAn, D = sDm) %>% 
  select(-RNASEQ_s4_0R)

diffase_snp_tf <- ddmaf_df_filt %>%
  select(symbol, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id")

ptnts <- c("A", "B", "D")

result_missenses <- diffase_snp_tf %>%
  select(-SNP_id, -symbol) %>% distinct() %>% 
  inner_join(add_comparisons(ptnts), by = "patients") %>%
  left_join(gt, by = "TF", relationship = "many-to-many") %>%
  rowwise() %>%
  mutate(allele1 = get(patient1), allele2 = get(patient2)) %>%
  filter(!is.na(allele1), !is.na(allele2)) %>%
  filter(allele1 != allele2) %>%
  # filter(allele1 != "0/1" & allele2 != "0/1") %>%
  ungroup() %>%
  distinct(patients, TF, refsnp_id, .keep_all = TRUE) %>%
  select(TF, refsnp_id, chr, pos, REF, ALT, patient1, allele1, patient2, allele2, consequence_type_tv) %>%
  mutate(allele1 = translate_gt(allele1, REF, ALT),
         allele2 = translate_gt(allele2, REF, ALT))  

missense_tfs <- sort(unique(result_missenses$TF))
length(missense_tfs) # 194
length(unique(na.omit(gt_vtr)$TF)) # 212
length(unique(missense_tfs)) / length(unique(gt$TF)) # 0.85
writeLines(missense_tfs, "~/less_in_vitro_missenses.txt")


## choosing snps in promoters of differentially expressed or missense tfs
causal_tf_rsids <- ddmaf_df_filt %>% 
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id") %>%
  filter(TF %in% union(all_degs, missense_tfs))

causal_snps_df <- causal_tf_rsids %>% select(SNP_id) %>%
  merge(joined_df[, c("chr", "pos", "ref", "alt", "rsid")], by.x = "SNP_id", by.y = "rsid") %>%
  distinct()

write_tsv(causal_snps_df[, c('chr', 'pos')], '/media/leon/DISK2/icig/done/85_snps.tsv', col_names = F) # 260

# выбираем пациентов ин виво гетерозигот по снп связанным с дэгами или миссенс тфами
# cd /media/leon/DISK2/icig/chipseq_july/stats_for_promoter_snps/; for i in 3 4 5 6 8 9 12 15; do perl /media/leon/DISK2/icig/done/nextflow/bin/stats.pl s${i} ../alignments/CHIPSEQ_s${i}.bam; done

check_snps <- function(path) {
  ddf <- read_tsv(path, show_col_types = FALSE)
  colnames(ddf) <- c("patient", "chr", "pos", "dp", "ref", "QS_ref", "alt", "QS_alt", "alt2", "QS_alt2")
  ddf <- ddf %>%
    select(-alt2, -QS_alt2) %>%
    filter(str_length(ref) == 1) %>%
    filter(str_length(alt) == 1) %>%
    mutate(dp_ref = round(dp * QS_ref), dp_alt = round(dp * QS_alt)) %>%
    select(-QS_ref, -QS_alt, -dp)
  return(ddf)
}

# joining dfs with the coverage of heterozygous snps motifbreaking DE or missense tfs in in vivo chipseq
chipseq_tables <- lapply(Sys.glob("/media/leon/DISK2/icig/chipseq_july/stats_for_promoter_snps/*.stat"), check_snps)
chipseq_table <- bind_rows(chipseq_tables) %>%
  mutate(chr = factor(chr)) %>% 
  left_join(joined_df[, c('chr', 'pos', 'symbol', 'rsid')], by = c("chr", "pos")) %>% 
  distinct()

chipseq_table_tfs <- causal_tf_rsids %>%
  select(SNP_id, symbol, TF) %>% distinct() %>%
  group_by(SNP_id, symbol) %>% summarize(TFs = list(TF), .groups = 'drop') %>% # for each snp and corresponding diffase gene collecting tfs into a list
  rename(rsid = SNP_id) %>%
  merge(chipseq_table, by = c('rsid', 'symbol'))

# считаем ddMAF для in vivo RNA-seq
patients_vv <- c("s3", "s4", "s5", "s6", "s8", "s9", "s12", "s15")
ase_in_vivo <- lapply(patients_vv, function(x) {
  read_stats(x, min_DP = 100)
})

# выбираем только те гены в промоторах которых и находятся итоговые SNP (== гены по которым была наибольшая ddMAF)
ase_vv <- bind_rows(ase_in_vivo) %>% 
  filter(symbol %in% chipseq_table_tfs$symbol) %>%
  mutate(padj = p.adjust(p01, method = "BH")) %>% 
  mutate(diffASE  = padj < 0.01) %>%
  mutate(total_DP = DP_0R + DP_1R) %>%
  group_by(symbol, chr, pos) %>%
  mutate(n = n(), total_across_dp = sum(total_DP),
         patients_vec = list(unique(patient))) %>% 
  filter(n > 1) %>%
  mutate(
    patients = map(patients_vec, ~ {
      combn(sort(.x), 2, FUN = function(p) paste(p, collapse = ""), simplify = TRUE)
    })
  ) %>%
  unnest(patients) %>%
  distinct()

df_collapsed_vv <- chipseq_table_tfs %>%
  group_by(chr, pos, ref, alt, symbol, rsid) %>%
  summarise(patients_vec = list(sort(unique(patient))), .groups = "drop") %>%
  mutate(
    patients = map(patients_vec, ~ {
      if (length(.x) > 1) {
        combn(.x, 2, FUN = function(p) paste(p, collapse = ""), simplify = TRUE)
      } else {
        .x  # keep single patient as-is if only one
      }
    })
  ) %>%
  select(-patients_vec) %>%
  unnest(patients) %>%
  distinct()

ddmaf_df_vv <- merge(ase_vv, df_collapsed_vv[, c("symbol", "rsid", "patients")], by = c("symbol", "patients"), suffixes = c(".RNA", ".CS")) %>% 
  group_by(symbol) %>% filter(any(diffASE == TRUE) & n() > 1) %>% 
  filter(n == max(n)) %>% filter(total_across_dp == max(total_across_dp)) %>% 
  group_by(across(-rsid)) %>% summarize(rsids = list(unique(rsid)), .groups = 'drop') %>%
  group_by(across(-patients)) %>% summarize(patients = list(unique(patients)), .groups = 'drop') %>%
  select(symbol, chr, pos, patient, ref_frac_diff, padj, rsids, n, patients) %>% distinct() %>%
  group_by(symbol) %>% 
  group_modify(~ { 
      patients <- .x$patient 
      comparisons <- .x$patients
      ref_frac_diffs <- .x$ref_frac_diff 
      pValues <- .x$padj
      common_snps <- .x$rsids
      combn(length(patients), 2, simplify = FALSE) %>%
        map_dfr(function(pair) {
          i <- pair[1]
          j <- pair[2]
          p1 <- patients[i]
          p2 <- patients[j]
          ref_frac_diff_1 <- ref_frac_diffs[i]
          ref_frac_diff_2 <- ref_frac_diffs[j]
          pValue_1 <- pValues[i]
          pValue_2 <- pValues[j]
          ptnts <- ifelse(p1 < p2, paste0(p1, p2), paste0(p2, p1))
          
          # есть хотя бы один общий снп + хотя бы у одного достоверно диффасе => считаем модуль разницы изменения доли референсного аллеля
          dd_ref_frac <- if(ptnts %in% comparisons & (min(c(pValue_1, pValue_2)) < pthrshld)) abs(ref_frac_diff_1 - ref_frac_diff_2) else NA
          tibble(
            patients = ptnts, 
            dd_ref_frac = dd_ref_frac,
            ref_frac_diff_1 = ref_frac_diff_1,
            ref_frac_diff_2 = ref_frac_diff_2,
            padj1 = pValues[i],
            padj2 = pValues[j] 
          )})}) %>%
  ungroup() %>% distinct() %>%
  filter(!is.na(dd_ref_frac)) %>%
  arrange(desc(dd_ref_frac))
    
ddmaf_df_filt_vv <- ddmaf_df_vv %>% filter(abs(dd_ref_frac) > 0.136)
length(unique(ddmaf_df_filt_vv$symbol)) # 22 из 100
length(unique(unlist(ddmaf_df_filt_vv$rsids))) # 37

TFs_to_check_invivo <- causal_tf_rsids %>%
  filter(symbol %in% ddmaf_df_filt_vv$symbol) %>% pull(TF) %>% unique()

length(TFs_to_check_invivo) # 137


####################### DESEQ IN VIVO
counts_vv <- read.table("~/counts_in_vivo_17.09.tsv", header = TRUE, row.names = 1, sep = "\t")
colnames(counts_vv) <- gsub("X.media.leon.DISK2.icig.done.alignments.RNASEQ_", "", gsub("R.bam", "", colnames(counts_vv)))

coldata_vv <- read.csv("/media/leon/DISK2/icig/done/coldata.csv", header = TRUE, row.names = 1, sep = ';')

rownames(coldata_vv) <- gsub('R.bam', '', rownames(coldata_vv))
rownames(coldata_vv) <- gsub('RNASEQ_', '', rownames(coldata_vv))

rename_map <- c("s6_0" = "s7_0", "s7_0" = "s6_0", "s6_1" = "s7_1", "s7_1" = "s6_1", "s6_30" = "s7_30", "s7_30" = "s6_30", "s6_90" = "s7_90", "s7_90" = "s6_90")

rownames(coldata_vv) <- ifelse(
  rownames(coldata_vv) %in% names(rename_map),
  rename_map[rownames(coldata_vv)],
  rownames(coldata_vv)
)

coldata_vv <- coldata_vv %>% 
  filter(patient %in% patients_vv) %>% 
  select(time, patient)

counts_vv <- counts_vv[, match(rownames(coldata_vv), colnames(counts_vv))]
all(rownames(coldata_vv) == colnames(counts_vv))

TF_counts_vv <- counts_vv %>%
  mutate(ensembl = rownames(.)) %>%
  merge(annotations_ordered, by = "ensembl") %>%
  filter(symbol %in% TFs_to_check_invivo)

rownames(TF_counts_vv) <- TF_counts_vv$symbol

TF_counts_vv <- TF_counts_vv %>% select(-entrezid, -symbol, -ensembl) 

ddsvv <- DESeqDataSetFromMatrix(countData = TF_counts_vv,
                                colData = coldata_vv,
                                design = ~ patient)
ddsvv <- estimateSizeFactors(ddsvv)
keepvv <- rowSums(counts(ddsvv, normalized = TRUE) >= 5 ) >= 2
ddsvv <- ddsvv[keepvv, ]
ddsvv <- DESeq(ddsvv)
resvv <- results(ddsvv)

rs_gene_ptnt_tf_vv <- ddmaf_df_filt_vv %>%   
  select(symbol, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id") 

cmbntns <- unique(rs_gene_ptnt_tf_vv$patients)

TF_deseq_results_vv <- data.frame()
for (pair in cmbntns) {
  p1 <- gsub("(s\\d+)(s\\d+)", "\\1", pair)
  p2 <- gsub("(s\\d+)(s\\d+)", "\\2", pair)
  
  res_pairwise <- as.data.frame(results(ddsvv, contrast = c("patient", p1, p2))) %>%
    mutate(symbol = rownames(.)) %>%
    mutate(contrast = paste0(p1, "_vs_", p2)) %>%
    filter(symbol %in% rs_gene_ptnt_tf_vv$TF[rs_gene_ptnt_tf_vv$patients == pair]) 
  TF_deseq_results_vv <- bind_rows(TF_deseq_results_vv, res_pairwise)
}

TF_deseq_results_vv <- TF_deseq_results_vv %>% 
  filter(padj < 0.1 & abs(log2FoldChange) > 0.5) 
all_degs_vv <- unique(TF_deseq_results_vv$symbol)
length(all_degs_vv) # 21
sort(all_degs_vv)

length(intersect(all_degs, all_degs_vv)) # 16

########## MISSENSES IN VIVO

in_vivo_missense_snps <- mssns_snps %>%
  filter(TF %in% TFs_to_check_invivo)

in_vivo_missense_snps %>%
  mutate(end = chrom_start) %>%
  mutate(chrom_start = chrom_start - 1) %>%
  select(chr_name, chrom_start, end) %>% distinct() %>%
  write.table("~/new_in_vivo_missense_snps.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# awk 'NF==3 && $1!="NA"' ~/in_vivo_missense_snps.bed > ~/in_vivo_missense_snps.clean.bed
# bcftools mpileup -R ~/in_vivo_missense_snps.clean.bed -f /media/leon/Polina/Genomes/hg38.chromFa/hg38_FOR_HISAT.fa RNASEQ_s*.merged.bam -d 1000 | bcftools call -mv -Oz -o in_vivo_missenses.more.vcf
# bcftools filter -i 'MIN(DP)>=20' in_vivo_missenses.more.vcf -Oz -o in_vivo_missenses.more.dp20.vcf; mv in_vivo_missenses.more.dp20.vcf in_vivo_missenses.more.vcf

vcf_vv <- read.vcfR("/media/leon/DISK2/icig/done/alignments/ALL_DIABETICS_VARIANTS.FILTERED.vcf.gz")

refalt_vv <- as.data.frame(vcf_vv@fix) %>%
  select(CHROM, POS, REF, ALT) %>%
  rename(chr = CHROM, pos = POS) %>%
  separate_rows(ALT, sep = ",")

gt_vv <- as.data.frame(extract.gt(vcf_vv)) %>%
  rownames_to_column("variant") %>%
  separate(variant, into = c("chr", "pos"), sep = "_", remove = TRUE) %>%
  distinct() %>%
  merge(in_vivo_missense_snps, by.x = c("chr", "pos"), by.y = c("chr_name", "chrom_start")) %>%
  left_join(refalt_vv, by = c("chr", "pos"), relationship = "many-to-many") %>%
  filter(str_length(REF) == 1) %>%
  filter(str_length(ALT) == 1) %>%
  rename(s3 = RNASEQ_s3_0R, s4 = RNASEQ_s4_0R, s5 = RNASEQ_s5_0R, s6 = RNASEQ_s6_0R, s8 = RNASEQ_s8_0R, s9 = RNASEQ_s9_0R, s12 = RNASEQ_s12_0R, s15 = RNASEQ_s15_0R)

diffase_snp_tf_vv <- ddmaf_df_filt_vv %>%
  select(symbol, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id")

result_missenses_vv <- diffase_snp_tf_vv %>%
  select(-SNP_id, -symbol) %>% distinct() %>% 
  inner_join(add_comparisons(patients_vv), by = "patients") %>%
  left_join(gt_vv, by = "TF", relationship = "many-to-many") %>%
  rowwise() %>%
  mutate(allele1 = get(patient1), allele2 = get(patient2)) %>%
  filter(!is.na(allele1), !is.na(allele2)) %>%
  filter(allele1 != allele2) %>%
  # filter(allele1 != "0/1" & allele2 != "0/1") %>%
  ungroup() %>%
  distinct(patients, TF, refsnp_id, .keep_all = TRUE) %>%
  select(TF, refsnp_id, chr, pos, REF, ALT, patient1, allele1, patient2, allele2, consequence_type_tv) %>%
  mutate(allele1 = translate_gt(allele1, REF, ALT),
         allele2 = translate_gt(allele2, REF, ALT))  

vv_missenses <- sort(unique(result_missenses_vv$TF))
length(vv_missenses)
both_missenses <- union(vv_missenses, missense_tfs)

writeLines(result_missenses$refsnp_id, "~/in_vivo/missense_snps.txt")

########## annotating missense snps with uniprot
### https://www.uniprot.org/id-mapping/
### Gene_Name -> UniProtKB, 9606
uniprot_ids <- read_tsv("~/in_vivo/idmapping_2025_10_14.tsv", show_col_types = FALSE) %>% 
  filter(Reviewed == "reviewed") %>% 
  select(From, Entry)

result_missenses %>% 
  select(TF, refsnp_id) %>%
  distinct() %>% 
  left_join(uniprot_ids, by = c("TF" = "From")) %>% 
  select(Entry, refsnp_id) %>% 
  write_tsv("~/in_vivo/tf_snp_df.tsv")

## /media/leon/DISK2/icig/done/uniprot/uniprot.ipynb

uniprot <- read_tsv("~/in_vivo/uniprot_isoforms_with_domains_2.tsv") %>% 
  left_join(uniprot_ids, by = c("uniprot_id" = "Entry")) %>% 
  rename(symbol = From)

length(unique(both_missenses)) # 209
length(intersect(both_missenses, uniprot$symbol)) # 204

uniprot_with_domains <- uniprot %>% 
  filter(!is.na(Domain))

length(unique(uniprot_with_domains$symbol)) # 154

uniprot_ordered <- uniprot_with_domains %>% 
  filter(Domain != "REGION|Disordered")

length(unique(uniprot_ordered$symbol)) # 123

####### FUNCTIONAL ANNOTATION
all_symbols <- unique(all_symbols[!is.na(all_symbols)])
all_entrezids <- mapIds(org.Hs.eg.db,
                        keys = all_symbols,
                        keytype = "SYMBOL",
                        column = "ENTREZID")

go_list <- list()
for (i in seq(nrow(TF_genes_df))) {
  gen <-  TF_genes_df[i, ][[1]]
  print(gen)
  try({
    tfdf <- as.data.frame(enrichGO(gene = union(all_degs, missense_tfs),
                                   OrgDb = org.Hs.eg.db,
                                   keyType = "SYMBOL",
                                   ont = "BP",
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   universe = all_symbols,
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
    tfdf <- as.data.frame(enrichKEGG(gene = all_entrezids[names(all_entrezids) %in% union(all_degs, missense_tfs)],
                                     organism = "hsa",
                                     pAdjustMethod = "BH",
                                     universe = all_entrezids,
                                     qvalueCutoff = 0.1))
    if (nrow(tfdf) > 0) {
      tfdf$TF <- gen
      kegg_list[[gen]] <- tfdf
    }
  })
}
kegg_df <- bind_rows(kegg_list)

########## SUPPLEMENTARY TABLE
# rsid  diffase_gene  patients  ddmaf_in_vivo tfs n_tfs
table1 <- ddmaf_df_filt %>%
  unnest(rsids) %>%
  select(rsids, symbol, patients, dd_ref_frac) %>%
  mutate(dd_ref_frac = round(dd_ref_frac, 3)) %>%
  group_by(rsids, symbol, patients) %>%
  slice_max(ddMAF, n = 1, with_ties = TRUE) %>%
  ungroup() %>% select(-comp) %>% distinct() %>% 
  rename(rsid = rsids, diffase_gene = symbol, ddref_in_vitro = dd_ref_frac) %>% 
  inner_join(snpgenedf, by = c("rsid" = "SNP_id")) %>%
  group_by(rsid, diffase_gene, patients, ddref_in_vitro) %>% 
  reframe(tfs = list(unique(TF)), n_tfs = n_distinct(TF))

table1_vv <- ddmaf_vv_df %>% drop_na() %>% 
  unnest(rsids) %>%
  select(rsids, symbol, patients, ddMAF, comp) %>%
  mutate(ddMAF = round(ddMAF, 3)) %>%
  group_by(rsids, symbol, patients) %>%
  slice_max(ddMAF, n = 1, with_ties = TRUE) %>%
  ungroup() %>% select(-comp) %>% distinct() %>% 
  rename(rsid = rsids, diffase_gene = symbol, ddmaf_in_vivo = ddMAF) %>% 
  inner_join(snpgenedf, by = c("rsid" = "SNP_id")) %>%
  group_by(rsid, diffase_gene, patients, ddmaf_in_vivo) %>% 
  reframe(tfs = list(unique(TF)), n_tfs = n_distinct(TF))

snp_patients <- joined_df_filt %>%
  group_by(rsid) %>% summarize(present_at = paste0(sort(c(patient)), collapse = "")) %>%
  distinct()

# rsid  present_at
snp_patients_vv <- chipseq_table %>%
  group_by(rsid) %>%
  summarize(n_patients = n_distinct(patient),
            present_at = paste0(sort(c(patient)), collapse = "")) %>%
  distinct()

table2 <- table1 %>% 
  left_join(snp_patients, by = "rsid") %>%
  group_by(rsid, diffase_gene) %>%
  summarize(
    patients_ddmaf = paste0(patients, ":", ddmaf_in_vitro, collapse = ", "),
    tfs = list(sort(unique(unlist(tfs)))),
    n_tfs = length(unique(unlist(tfs))),
    present_at = first(present_at),
    .groups = "drop"
  ) %>%
  mutate(
    n_degs = vapply(tfs, function(x) length(intersect(x, both_degs)), integer(1)),
    n_mssnss = vapply(tfs, function(x) length(intersect(x, both_missenses)), integer(1)),
    n_fnctnl_mssnss = vapply(tfs, function(x) length(intersect(x, uniprot_ordered$symbol)), integer(1))
  ) %>% 
  # mutate(str_tfs = sapply(tfs, function(x) paste(sort(x), collapse = ", "))) %>% 
  left_join(snp_patients_vv[, c('rsid', 'n_patients')], by = 'rsid') %>% 
  rename(n_ptnts_in_vivo = n_patients) %>% 
  mutate(n_ptnts_in_vivo = ifelse(is.na(table2$n_ptnts_in_vivo), 0, table2$n_ptnts_in_vivo)) %>% 
  mutate(n_degs_in_vivo = vapply(tfs, function(x) length(intersect(x, all_degs_vv)), integer(1)),
         n_mssnss_in_vivo = vapply(tfs, function(x) length(intersect(x, vv_missenses)), integer(1)))

table2_vv <- table1_vv %>% 
  left_join(snp_patients_vv, by = "rsid") %>%
  group_by(rsid, diffase_gene) %>%
  summarize(
    patients_ddmaf = paste0(patients, ":", ddmaf_in_vivo, collapse = ", "),
    tfs = list(sort(unique(unlist(tfs)))),
    n_tfs = length(unique(unlist(tfs))),
    present_at = first(present_at),
    ddmaf_signif = sum(ddmaf_in_vivo > 0.121),
    max_ddmaf = max(ddmaf_in_vivo),
    mean_ddmaf = round(mean(ddmaf_in_vivo), 3),
    .groups = "drop"
  ) %>%
  mutate(
    n_degs = vapply(tfs, function(x) length(intersect(x, both_degs)), integer(1)),
    n_missenses = vapply(tfs, function(x) length(intersect(x, both_missenses)), integer(1))
  ) %>% 
  # mutate(str_tfs = sapply(tfs, function(x) paste(sort(x), collapse = ", "))) %>% 
  arrange(desc(n_tfs))

table2 <- table2 %>% 
  left_join(table2_vv[, c('rsid', 'ddmaf_signif', 'max_ddmaf', 'mean_ddmaf')], by = 'rsid') 

table2 <- table2 %>% 
  mutate(ddmaf_signif = ifelse(is.na(table2$ddmaf_signif), 0, table2$ddmaf_signif)) %>% 
  arrange(desc(max_ddmaf)) %>% 
  mutate(tfs = sapply(tfs, function(x) paste(sort(x), collapse = ", "))) %>% 
  select(rsid, diffase_gene, present_at, patients_ddmaf, n_tfs, n_degs, n_mssnss, n_fnctnl_mssnss, n_ptnts_in_vivo, ddmaf_signif, max_ddmaf, mean_ddmaf, n_degs_in_vivo, n_mssnss_in_vivo, tfs)

# write_tsv(table2, "~/supplementary_table_for_snps.tsv")
# write_tsv(table2[, "rsid"], "~/snps_to_check_in_vep.tsv", col_names = F)

vep <- read_tsv("~/Downloads/promoter_snps_vep_results.txt") %>% 
  rename(rsid = `#Uploaded_variation`) %>% 
  inner_join(table2[, c('rsid', 'diffase_gene')], by = c('rsid' = 'rsid', 'symbol' = 'diffase_gene')) %>% 
  select(rsid, symbol, Consequence, CLIN_SIG, PUBMED) %>% 
  distinct()
