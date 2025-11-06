# packages <- c("jsonlite", "tidyverse", "BiocManager")
# install.packages(setdiff(packages, rownames(installed.packages())))
# 
# bioc <- c(
#   "DESeq2", "org.Hs.eg.db", "biomaRt", "AnnotationDbi", "motifbreakR",
#   "MotifDb", "clusterProfiler", "TxDb.Hsapiens.UCSC.hg38.knownGene",
#   "BSgenome.Hsapiens.UCSC.hg38", "SNPlocs.Hsapiens.dbSNP155.GRCh38"
# )
# 
# for (pkg in bioc) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     BiocManager::install(pkg) # версии
#   }
# }

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(biomaRt)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(DESeq2)
library(clusterProfiler)
library(tidyverse)

source("scripts/1_functions.R")
source("scripts/2_snp_stats_31.10.R")
source("scripts/4_deseq.R")
source("scripts/5_missenses.R")

## in vitro chip-seq snps
if (!exists("IN_VIVO")) {
  patients <- c("sAn", "sBn", "sDm")
  counts_path <- "counts/in_vitro_counts.tsv"
  coldata_path <- "counts/in_vitro_coldata.tsv"
} else {
  patients <- c("s1", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s12", "s15")
  counts_path <- "counts/in_vivo_counts.tsv"
  coldata_path <- "counts/in_vivo_coldata.tsv"
}

# joined_df_vv <- joined_df; ase_vv <- ase; chipseq_df_vv <- chipseq_df; ddmaf_df_vv <- ddmaf_df; ddmaf_df_filt_vv <- ddmaf_df_filt; snpgenedf_vv <- snpgenedf; counts_vv <- counts; coldata_vv <- coldata; patients_vv <- patients

joined_df <- lapply(patients, read_chipseq) %>% bind_rows() %>%
  mutate(variant = paste(chr, pos, ref, alt, sep = '_')) %>%
  group_by(variant) %>% filter(n_distinct(patient) > 1) %>% ungroup() # filtering snps found in at least 2 patients

# filtering snps in promoters
joined_df <- assign_genes(joined_df, "chr", "pos", promoters = TRUE)
joined_df <- add_rsids(joined_df, "chr", "pos", "ref")

#### in vitro RNA-seq
stats_dir <- "data/"
ase_dfs <- lapply(patients, function(x) {
  read_rnaseq(x, min_DP = 50)
})

ase <- prepare_rna(ase_dfs, joined_df)
chipseq_df <- prepare_cs(joined_df, ase)
input_for_glm <- prepare_for_py(ase, chipseq_df)

subset_for_glm <- input_for_glm %>% select(symbol, patients, contains("DP"), frac_common) %>% 
  drop_na() # откуда

res_df1 <- run_glm(subset_for_glm)

if ("n_common" %in% names(subset_for_glm)) {
  swapped <- subset_for_glm %>% filter(n_common == 0)
} else if ("frac_common" %in% names(subset_for_glm)) {
  swapped <- subset_for_glm %>% filter(frac_common < 0.5)
}

# Swap patient 2 columns for genes with not enough common SNPs
for (pref in c("DP_maj_gene", "DP_min_gene")) {
  tmp <- swapped[[paste0(pref, "_0_p2")]]
  swapped[[paste0(pref, "_0_p2")]] <- swapped[[paste0(pref, "_1_p2")]]
  swapped[[paste0(pref, "_1_p2")]] <- tmp
}

res_df2 <- run_glm(swapped)

merged <- bind_rows(res_df1, res_df2) %>%
  arrange(desc(p_interaction)) %>%
  distinct(symbol, patients, .keep_all = TRUE) %>% 
  mutate(padj = p.adjust(p_interaction, method = "BH"))

ddmaf_df_filt <- merged %>% 
  inner_join(chipseq_df %>% select(symbol, patients, rsid)) %>% 
  filter(padj < 0.1)

length(unique(ddmaf_df_filt$symbol)) # 116

snps_for_mtfbrkr <- unique(ddmaf_df_filt$rsid)
length(snps_for_mtfbrkr) # 240

#### motifbreakR 
mtfbrkr_dir <- "motifbreakr"
result_files <- list.files(mtfbrkr_dir, pattern = ".*\\.rds$") 
rsids_in_folder <- gsub("\\.rds", "", result_files)
rsids_to_process <- setdiff(snps_for_mtfbrkr, rsids_in_folder)

# writeLines(rsids_to_process, "snps_for_mtfbrkr_to_add.txt")
## Rscript scripts/motifbreakr.R snps_for_mtfbrkr_to_add.txt motifbreakr/ 30 

if (length(rsids_to_process)) {
  run_motifbreakr(rsids_to_process, mtfbrkr_dir, batch_size = 10) 
}
result_files <- list.files(mtfbrkr_dir, pattern = ".*\\.rds$") 
rsids_in_folder <- gsub(".rds", "", result_files)
result_files <- result_files[rsids_in_folder %in% snps_for_mtfbrkr]

all_results <- list()
for (file in result_files) {
  res <- readRDS(paste(mtfbrkr_dir, file, sep = '/'))
  rsid <- sub(".rds", "", file)
  if (!rsid %in% snps_for_mtfbrkr) { next }
  if (length(res) > 0) {
    mcols(res)$rsid <- rsid
    all_results[[rsid]] <- res
  }
}
combined_results <- do.call(c, unname(all_results))
hcmc <- combined_results[combined_results$dataSource == "HOCOMOCOv13"]

# TFs motifbroken by SNPs in promoters of genes with highest ddMAF
snpgenedf <- hcmc %>% mcols() %>% as.data.frame() %>% 
  select(geneSymbol, SNP_id) %>% rename(tf = geneSymbol) %>% distinct() %>% 
  mutate(tf = recode(tf, !!!tf_to_gene))

tfs <- unique(snpgenedf$tf)
length(tfs) # 640
length(unique(snpgenedf$SNP_id)) # 234

counts <- read.table(counts_path, header = TRUE, sep = "\t")
coldata <- read.table(coldata_path, header = TRUE, sep = "\t", stringsAsFactors = TRUE)

comparisons_df <- ddmaf_df_filt %>% select(symbol, patients, rsid) %>% # gene with high ddMAF => patients having high ddMAF at this gene => promoter SNP => corresponding TF
  rename(SNP_id = rsid) %>% inner_join(snpgenedf, relationship = "many-to-many") %>% distinct(patients, tf)

deseq_df <- find_tf_degs(counts, coldata, tfs, comparisons_df)
all_degs <- sort(unique(deseq_df$symbol))

length(all_degs) # 102

check <- readLines("output/2025.10.28.in_vitro_degs.dp100_padj01_logROR01.txt")
length(intersect(check, all_degs))
setdiff(check, all_degs)
setdiff(all_degs, check)

result_missenses <- find_missense_tfs(patients, tfs, comparisons_df)
homo_missenses <- result_missenses %>%
  filter(substr(allele1, 1, 1) == substr(allele1, 2, 2) & substr(allele2, 1, 1) == substr(allele2, 2, 2)) %>% 
  pull(tf) %>% unique() %>% sort()

hetero_missenses <- result_missenses %>%
  filter(substr(allele1, 1, 1) != substr(allele1, 2, 2) | substr(allele2, 1, 1) != substr(allele2, 2, 2)) %>% 
  pull(tf) %>% unique() %>% sort()

missense_tfs <- sort(unique(result_missenses$tf))
length(missense_tfs)                      
length(intersect(missense_tfs, all_degs))

# writeLines(all_degs, "output/2025.10.31.in_vitro_degs.dp50_padj01_half_common_snps.txt")
# writeLines(missense_tfs, "output/2025.10.31.in_vitro_missenses.dp50_padj01_half_common_snps.txt")

## filtering snps found at IN VIVO DIABETICS in promoters of TFs that were differentially expressed or carried missense variants in IN VITRO dataset
if (!exists("IN_VIVO")) {
  causal_snps <- snpgenedf %>% filter(tf %in% union(all_degs, missense_tfs))
  degs_in_vitro <- all_degs; missenses_in_vitro <- missense_tfs
  missense_df_vtr <- result_missenses
  IN_VIVO <- TRUE # run again up to all_degs ## сделать нормально
}

chipseq_subset <- chipseq_df %>% filter(rsid %in% causal_snps$SNP_id)
ase_subset <- ase %>% filter(symbol %in% chipseq_subset$symbol)

input_for_glm_ss <- prepare_for_py(ase_subset, chipseq_subset)

subset_for_glm_ss <- input_for_glm_ss %>% select(symbol, patients, contains("DP"), frac_common) %>% 
  drop_na()

res_df1_ss <- run_glm(subset_for_glm_ss)

if ("n_common" %in% names(subset_for_glm_ss)) {
  swapped_ss <- subset_for_glm_ss %>% filter(n_common == 0)
} else if ("frac_common" %in% names(subset_for_glm_ss)) {
  swapped_ss <- subset_for_glm_ss %>% filter(frac_common < 0.5)
}

# Swap patient 2 columns for genes with not enough common SNPs
for (pref in c("DP_maj_gene", "DP_min_gene")) {
  tmp <- swapped_ss[[paste0(pref, "_0_p2")]]
  swapped_ss[[paste0(pref, "_0_p2")]] <- swapped_ss[[paste0(pref, "_1_p2")]]
  swapped_ss[[paste0(pref, "_1_p2")]] <- tmp
}

res_df2_ss <- run_glm(swapped_ss)

merged_ss <- bind_rows(res_df1_ss, res_df2_ss) %>%
  arrange(desc(p_interaction)) %>%
  distinct(symbol, patients, .keep_all = TRUE) %>% 
  mutate(padj = p.adjust(p_interaction, method = "BH"))

ddmaf_df_filt_ss <- merged_ss %>% 
  inner_join(chipseq_subset %>% select(symbol, patients, rsid)) %>% 
  filter(padj < 0.1) %>% 
  inner_join(causal_snps, by = c("rsid" = "SNP_id"), relationship = "many-to-many") %>% distinct()

length(unique(ddmaf_df_filt_ss$symbol)) # 40

tfs_to_check <- unique(ddmaf_df_filt_ss$tf)
comparisons_ss <- comparisons_df %>% filter(tf %in% tfs_to_check)
deseq_ss <- find_tf_degs(counts, coldata, tfs_to_check, comparisons_ss)
degs_ss <- sort(unique(deseq_ss$symbol))
degs_conf <- intersect(degs_ss, degs_in_vitro)
length(degs_conf) # 50

missense_ss_df <- find_missense_tfs(patients, tfs_to_check, comparisons_ss)
missense_ss <- sort(unique(missense_ss_df$tf)) # 63
missense_conf <- intersect(missense_ss, missenses_in_vitro)
length(missense_conf) # 45

missense_conf_homo <- missense_ss_df %>%
  filter(substr(allele1, 1, 1) == substr(allele1, 2, 2) & substr(allele2, 1, 1) == substr(allele2, 2, 2)) %>% 
  pull(tf) %>% unique() %>% sort() %>%
  intersect(missense_conf)

# writeLines(degs_conf, "output/2025.10.31.degs_confirmed_in_vivo.dp50_padj01_half_common_snps.txt")
# writeLines(missense_conf, "output/2025.10.31.missenses_confirmed_in_vivo.dp50_padj01_half_common_snps.txt")
