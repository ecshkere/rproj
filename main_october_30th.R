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
source("scripts/2_snp_stats_27.10.R")
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

joined_df <- lapply(patients, read_chipseq) %>% bind_rows() %>%
  mutate(variant = paste(chr, pos, ref, alt, sep = '_')) %>%
  group_by(variant) %>% filter(n_distinct(patient) > 1) %>% ungroup() # filtering snps found in at least 2 patients

# filtering snps in promoters
joined_df <- assign_genes(joined_df, "chr", "pos", promoters = TRUE)
joined_df <- add_rsids(joined_df, "chr", "pos", "ref")

#### in vitro RNA-seq
stats_dir <- "data/"
ase_dfs <- lapply(patients, function(x) {
  read_rnaseq(x, min_DP = 100)
})

ase <- prepare_rna(ase_dfs, joined_df)
chipseq_df <- prepare_cs(joined_df, ase)
ddmaf_df <- prepare_for_glm(ase, chipseq_df)

glm_res <- run_glm(select(ddmaf_df, symbol, patients, starts_with("DP"))) %>%
  mutate(padj = p.adjust(p_interaction, method = "BH")) %>%
  select(symbol, patients, padj, coef_interaction)

length(unique(glm_res$symbol[glm_res$padj < 0.1])) # 111

ddmaf_df_filt <- ddmaf_df %>%
  inner_join(glm_res) %>% 
  filter(padj < 0.1, log2ROR > 0.1) %>%
  inner_join(chipseq_df %>% select(symbol, patients, rsid))

length(unique(ddmaf_df_filt$symbol)) # 99

snps_for_mtfbrkr <- unique(ddmaf_df_filt$rsid)
length(snps_for_mtfbrkr) # 242

#### motifbreakR 
mtfbrkr_dir <- "motifbreakr"
result_files <- list.files(mtfbrkr_dir, pattern = ".*\\.rds$") 
rsids_in_folder <- gsub("\\.rds", "", result_files)
rsids_to_process <- setdiff(snps_for_mtfbrkr, rsids_in_folder)

# writeLines(rsids_to_process, "snps_for_mtfbrkr_to_add.txt")
## Rscript scripts/motifbreakr.R snps_for_mtfbrkr_to_add.txt motifbreakr/ 30 

if (length(rsids_to_process)) {
  run_motifbreakr(rsids_to_process, mtfbrkr_dir, batch_size = 26) 
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
length(tfs) # 651
length(unique(snpgenedf$SNP_id)) # 232

counts <- read.table(counts_path, header = TRUE, sep = "\t")
coldata <- read.table(coldata_path, header = TRUE, sep = "\t", stringsAsFactors = TRUE)
  
comparisons_df <- ddmaf_df_filt %>% select(symbol, patients, rsid) %>% # gene with high ddMAF => patients having high ddMAF at this gene => promoter SNP => corresponding TF
  rename(SNP_id = rsid) %>% inner_join(snpgenedf, relationship = "many-to-many") %>% distinct(patients, tf)
  
deseq_df <- find_tf_degs(counts, coldata, tfs, comparisons_df)
all_degs <- sort(unique(deseq_df$symbol))
length(all_degs) # 107

result_missenses <- find_missense_tfs(patients, tfs, comparisons_df)
homo_missenses <- result_missenses %>%
  filter(qual == "homo") %>% pull(tf) %>% unique() %>% sort()

hetero_missenses <- result_missenses %>%
  filter(qual == "hetero") %>% pull(tf) %>% unique() %>% sort()

missense_tfs <- sort(unique(result_missenses$tf))
length(missense_tfs) # 82
length(intersect(missense_tfs, all_degs)) # 8

# writeLines(all_degs, "output/2025.10.28.in_vitro_degs.dp100_padj01_logROR01.txt")
# writeLines(missense_tfs, "output/2025.10.28.in_vitro_missenses.dp100_padj01_logROR01.txt")

## filtering snps found at IN VIVO DIABETICS in promoters of TFs that were differentially expressed or carried missense variants in IN VITRO dataset
if (!exists("IN_VIVO")) {
  causal_snps <- snpgenedf %>% filter(tf %in% union(all_degs, missense_tfs))
  degs_in_vitro <- all_degs; missenses_in_vitro <- missense_tfs; result_missenses_vtr <- result_missenses; deseq_df_vtr <- deseq_df
  joined_df_vtr <- joined_df; ase_vtr <- ase; chipseq_df_vtr <- chipseq_df; ddmaf_df_vtr <- ddmaf_df; ddmaf_df_filt_vtr <- ddmaf_df_filt
  snpgenedf_vtr <- snpgenedf; counts_vtr <- counts; coldata_vtr <- coldata; patients_vtr <- patients; comparisons_df_vtr <- comparisons_df
  IN_VIVO <- TRUE # run again up to all_degs ## сделать нормально
}

chipseq_subset <- chipseq_df %>% filter(rsid %in% causal_snps$SNP_id)
ase_subset <- ase %>% filter(symbol %in% chipseq_subset$symbol)

ddmaf_df_ss <- prepare_for_glm(ase_subset, chipseq_subset)

glm_res_ss <- run_glm(select(ddmaf_df_ss, symbol, patients, starts_with("DP"))) %>%
  mutate(padj = p.adjust(p_interaction, method = "BH")) %>%
  select(symbol, patients, padj, coef_interaction)

length(unique(glm_res_ss$symbol[glm_res_ss$padj < 0.1])) # 25

ddmaf_df_filt_ss <- ddmaf_df_ss %>%
  inner_join(glm_res_ss) %>% 
  filter(padj < 0.1, log2ROR > 0.1) %>%
  inner_join(chipseq_subset %>% select(symbol, patients, rsid)) %>% 
  inner_join(causal_snps, by = c("rsid" = "SNP_id"), relationship = "many-to-many") %>% distinct()

length(unique(ddmaf_df_filt_ss$symbol)) # 23 out of 100

tfs_to_check <- unique(ddmaf_df_filt_ss$tf)
comparisons_ss <- comparisons_df %>% filter(tf %in% tfs_to_check)
deseq_ss <- find_tf_degs(counts, coldata, tfs_to_check, comparisons_ss)
degs_ss <- sort(unique(deseq_ss$symbol))
degs_conf <- intersect(degs_ss, degs_in_vitro)
length(degs_conf) # 51

missense_ss_df <- find_missense_tfs(patients, tfs_to_check, comparisons_ss) %>% 
  filter(tf %in% missenses_in_vitro)
missense_conf <- sort(unique(missense_ss_df$tf))
length(missense_conf) # 41

missense_conf_homo <- missense_ss_df %>%
  filter(qual == "homo") %>% pull(tf) %>% unique() %>% sort()

missense_conf_hetero <- missense_ss_df %>%
  filter(qual == "hetero") %>% pull(tf) %>% unique() %>% sort()

# writeLines(degs_conf, "output/2025.10.28.degs_confirmed_in_vivo.p01logror01.txt")
# writeLines(missense_conf, "output/2025.11.06.missenses_confirmed_in_vivo.p01logror01.txt")

########## annotating missense snps with uniprot
uniprot_ids <- read_tsv("idmapping_2025_10_28.tsv") %>% 
  filter(Reviewed == "reviewed") %>% 
  select(From, Entry)

missense_df_vtr %>% select(tf, rsid) %>% distinct() %>% 
  left_join(uniprot_ids, by = c("tf" = "From"), relationship = "many-to-many") %>% select(Entry, rsid) %>% 
  write_tsv("for_uniprot_py_script.tsv")

## python3 ./scripts/uniprot_parser.py for_uniprot_py_script.tsv output/2025.10.28.uniprot_output_in_vitro.tsv

uniprot <- read_tsv("output/2025.10.28.uniprot_output_in_vitro.tsv") %>% 
  left_join(uniprot_ids, by = c("uniprot_id" = "Entry")) %>% 
  rename(symbol = From)

uniprot_with_domains <- uniprot %>% filter(!is.na(Domain)) %>% arrange(symbol)

write_tsv(uniprot_with_domains %>% select(symbol, rsid, Domain), "output/missense_snps_uniprot_annotation.tsv")

uniprot_ordered <- uniprot_with_domains %>% filter(Domain != "REGION|Disordered")
length(unique(uniprot_ordered$symbol)) 

####### FUNCTIONAL ANNOTATION
## background 
counts <- read.table("counts/in_vivo_counts.tsv", header = TRUE, sep = "\t")
coldata <- read.table("counts/in_vivo_coldata.tsv", header = TRUE, sep = "\t", stringsAsFactors = TRUE)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ patient)
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized = TRUE) >= 5 ) >= 15
background_ens <- names(keep[keep])
background_symb <- annotations_clean$symbol[annotations_clean$ensembl %in% background_ens]
background_entr <- annotations_clean$entrezid[annotations_clean$ensembl %in% background_ens]

go_degs <- go(all_degs)
kegg_degs <- kegg(all_degs)

go_missenses <- go(missense_tfs)
kegg_missenses <- kegg(missense_tfs)

go_ddase <- ddmaf_df_filt %>%
  # filter(rsid %in% (snpgenedf %>% filter(tf %in% union(all_degs, missense_tfs)) %>% pull(SNP_id))) %>%
  pull(symbol) %>% unique() %>% go()

kegg_ddase <- ddmaf_df_filt %>% 
  # filter(rsid %in% (snpgenedf %>% filter(tf %in% union(all_degs, missense_tfs)) %>% pull(SNP_id))) %>%
  pull(symbol) %>% unique() %>% kegg()

write_tsv(go_degs, "output/2025.10.28.GO_degs.tsv")
write_tsv(kegg_degs, "output/2025.10.28.KEGG_degs.tsv")
# write_tsv(go_missenses, "output/2025.10.28.GO_missenses.tsv")
write_tsv(kegg_missenses, "output/2025.10.28.KEGG_missenses.tsv")
write_tsv(go_ddase, "output/2025.10.28.GO_ddASE_genes.tsv")
# write_tsv(kegg_ddase, "output/2025.10.28.KEGG_ddASE_genes.tsv")










########## SUPPLEMENTARY TABLE
# rsid  diffase_gene  patients  ddmaf_in_vivo tfs n_tfs
table1 <- ddmaf_df_filt_vtr %>%
  select(rsid, symbol, patients, padj, log2ROR) %>%
  mutate(padj = round(padj, 2), log2ROR = round(log2ROR, 2)) %>% distinct() %>% 
  rename(diffase_gene = symbol, padj_in_vitro = padj, logROR_in_vitro = log2ROR, patients_in_vitro = patients) %>% 
  inner_join(snpgenedf_vtr, by = c("rsid" = "SNP_id"), relationship = "many-to-many") %>% 
  group_by(across(-tf)) %>% 
  reframe(tfs = list(unique(tf)), n_tfs = n_distinct(tf)) %>% distinct()

table1_vv <- ddmaf_df_filt_ss %>%
  select(rsid, symbol, patients, padj, log2ROR) %>%
  mutate(padj = round(padj, 2), log2ROR = round(log2ROR, 2)) %>% distinct() %>% 
  rename(diffase_gene = symbol, padj_in_vivo = padj, logROR_in_vivo = log2ROR, patients_in_vivo = patients) %>% 
  inner_join(snpgenedf, by = c("rsid" = "SNP_id"), relationship = "many-to-many") %>% 
  group_by(across(-tf)) %>% 
  reframe(tfs = list(unique(tf)), n_tfs = n_distinct(tf)) %>% distinct()

snp_patients <- joined_df_vtr %>% filter(symbol %in% ddmaf_df_filt_vtr$symbol) %>%
  group_by(rsid) %>% summarize(present_at = paste0(sort(c(patient)), collapse = "_")) %>%
  distinct()

snp_patients_vv <- joined_df %>% filter(symbol %in% ddmaf_df_filt_ss$symbol) %>%
  group_by(rsid) %>% summarize(present_at = paste0(sort(c(patient)), collapse = "_")) %>%
  distinct()


degs_freq_vtr <- enframe(sort(table(deseq_df_vtr[deseq_df_vtr$symbol %in% degs_in_vitro,]$symbol), decreasing = T), name = "tf", value = "times")
mssns_freq_vtr <- enframe(sort(table(result_missenses_vtr %>% group_by(tf, patient1, patient2) %>% summarize(.groups = "drop") %>%
                                       filter(tf %in% missenses_in_vitro) %>% pull(tf)), decreasing = T), name = "tf", value = "times")
mssns_freq_vtr <- enframe(sort(table(result_missenses_vtr %>% group_by(tf, patient1, patient2) %>% summarize(.groups = "drop") %>%
                                       filter(tf %in% missenses_in_vitro) %>% pull(tf)), decreasing = T), name = "tf", value = "times")

degs_freq_vv <- enframe(sort(table(deseq_ss[deseq_ss$symbol %in% degs_conf,]$symbol), decreasing = T), name = "tf", value = "times")
mssns_freq_vv <- enframe(sort(table(missense_ss_df %>% group_by(tf, patient1, patient2) %>% summarize(.groups = "drop") %>%
             filter(tf %in% missense_conf) %>% pull(tf)), decreasing = T), name = "tf", value = "times")

degs_freq <- full_join(degs_freq_vtr, degs_freq_vv, by = "tf", suffix = c("_vtr", "_vv"))
mssns_freq <- full_join(mssns_freq_vtr, mssns_freq_vv, by = "tf", suffix = c("_vtr", "_vv"))

ress <- full_join(degs_freq, mssns_freq, by = "tf", suffix = c("_deg", "_mssns")) %>% 
  rowwise() %>%
  mutate(sum = sum(c_across(2:5), na.rm = T)) %>%
  ungroup() %>%
  arrange(desc(sum))

table2 <- table1 %>% 
  left_join(snp_patients, by = "rsid") %>%
  group_by(rsid, diffase_gene) %>%
  summarize(
    ptnts_logROR_padj = paste0(patients_in_vitro, " : ", logROR_in_vitro, ", ", padj_in_vitro, collapse = "; "),
    tfs = list(sort(unique(unlist(tfs)))),
    n_tfs = length(unique(unlist(tfs))),
    present_at = first(present_at),
    .groups = "drop"
  ) %>%
  mutate(
    # n_degs = vapply(tfs, function(x) length(intersect(x, all_degs)), integer(1)),
    # n_mssnss = vapply(tfs, function(x) length(intersect(x, missense_tfs)), integer(1)),
    # n_fnctnl_mssnss = vapply(tfs, function(x) length(intersect(x, uniprot_ordered$symbol)), integer(1))
    degs = sapply(tfs, function(x) paste(sort(intersect(x, all_degs)), collapse = ", ")),
    mssnss = sapply(tfs, function(x) paste(sort(intersect(x, missense_tfs)), collapse = ", ")),
    fnctnl_mssnss = sapply(tfs, function(x) paste(sort(intersect(x, uniprot_ordered$symbol)), collapse = ", "))
  ) %>% 
  # mutate(tfs = sapply(tfs, function(x) paste(sort(x), collapse = ", "))) %>%
  # left_join(snp_patients_vv[, c('rsid', 'n_patients')], by = 'rsid') %>% 
  # rename(n_ptnts_in_vivo = n_patients) %>% 
  # mutate(n_ptnts_in_vivo = ifelse(is.na(table2$n_ptnts_in_vivo), 0, table2$n_ptnts_in_vivo)) %>% 
  # mutate(n_degs_in_vivo = vapply(tfs, function(x) length(intersect(x, all_degs_vv)), integer(1)),
  #        n_mssnss_in_vivo = vapply(tfs, function(x) length(intersect(x, vv_missenses)), integer(1)))
  select(-tfs)

table2_vv <- table1_vv %>% 
  left_join(snp_patients_vv, by = "rsid") %>%
  group_by(rsid, diffase_gene) #%>%
  summarize(
    patients_log2ROR = paste0(patients, ":", ddmaf_in_vivo, collapse = ", "),
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


