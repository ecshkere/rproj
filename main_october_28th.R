# packages <- c("jsonlite", "tidyverse", "BiocManager")
# install.packages(setdiff(packages, rownames(installed.packages())))
# 
# bioc <- c("DESeq2", "org.Hs.eg.db", "biomaRt", "AnnotationDbi", "motifbreakR",
#           "MotifDb", "clusterProfiler", "TxDb.Hsapiens.UCSC.hg38.knownGene",
#           "BSgenome.Hsapiens.UCSC.hg38", "SNPlocs.Hsapiens.dbSNP155.GRCh38")
# 
# for (pkg in bioc) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     BiocManager::install(pkg)
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
source("scripts/3_deseq.R")
source("scripts/4_missenses.R")

if (!exists("IN_VIVO")) {
  patients <- c("sAn", "sBn", "sDm")
  counts_path <- "counts/in_vitro_counts.tsv"
  coldata_path <- "counts/in_vitro_coldata.tsv"
} else {
  patients <- c("s1", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s12", "s15")
  counts_path <- "counts/in_vivo_counts.tsv"
  coldata_path <- "counts/in_vivo_coldata.tsv"
}

# stats_dir <- "data/"
stats_dir <- "data/unfiltered_stats/"

# in vitro chip-seq 
joined_df <- map_dfr(patients, read_chipseq) %>%
  mutate(variant = paste(chr, pos, ref, alt, sep = '_')) %>%
  filter(n_distinct(patient) > 1, .by = variant) # filtering snps common for 2+ patients

# filtering snps in promoters
joined_df <- assign_genes(joined_df, "chr", "pos", promoters = TRUE)
joined_df <- add_rsids(joined_df, "chr", "pos", "ref")

# in vitro RNA-seq
ase_dfs <- map_dfr(patients, read_rnaseq, min_DP = 100)
ase <- prepare_rna(ase_dfs, joined_df)
chipseq_df <- prepare_cs(joined_df, ase)
ddmaf_df <- prepare_for_glm(ase, chipseq_df)

# comparing diffASE
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

# writeLines(rsids_to_process, "output/snps_for_mtfbrkr_to_add.txt")
## Rscript scripts/motifbreakr.R output/snps_for_mtfbrkr_to_add.txt motifbreakr/ 30

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

# writeLines(all_degs, "results/2025.11.10.in_vitro_degs.dp100_padj01_logROR01.txt")
# writeLines(missense_tfs, "results/2025.11.10.in_vitro_missenses.dp100_padj01_logROR01.txt")

## filtering snps found among in vivo diabetics in promoters of TFs that were differentially expressed or carried missense variants in in vitro dataset
if (!exists("IN_VIVO")) {
  causal_snps <- snpgenedf %>% filter(tf %in% union(all_degs, missense_tfs))
  degs_in_vitro <- all_degs; missenses_in_vitro <- missense_tfs; result_missenses_vtr <- result_missenses; deseq_df_vtr <- deseq_df
  joined_df_vtr <- joined_df; ase_vtr <- ase; chipseq_df_vtr <- chipseq_df; ddmaf_df_vtr <- ddmaf_df; ddmaf_df_filt_vtr <- ddmaf_df_filt
  snpgenedf_vtr <- snpgenedf; counts_vtr <- counts; coldata_vtr <- coldata; patients_vtr <- patients; comparisons_df_vtr <- comparisons_df
  IN_VIVO <- TRUE # run again ## сделать нормально
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

# writeLines(degs_conf, "results/2025.11.10.degs_confirmed_in_vivo.p01logror01.txt")
# writeLines(missense_conf, "results/2025.11.10.missenses_confirmed_in_vivo.p01logror01.txt")

## annotating missense snps with uniprot
uniprot_ids <- read_tsv("output/idmapping_2025_11_10.tsv") %>% 
  filter(Reviewed == "reviewed") %>% 
  select(From, Entry)

result_missenses_vtr %>% select(tf, rsid) %>% distinct() %>% 
  left_join(uniprot_ids, by = c("tf" = "From"), relationship = "many-to-many") %>% select(Entry, rsid) #%>% 
# write_tsv("output/2025.11.10.for_uniprot_py_script.tsv")

## python3 ./scripts/uniprot_parser.py for_uniprot_py_script.tsv output/2025.10.28.uniprot_output_in_vitro.tsv

uniprot <- read_tsv("output/2025.11.10.uniprot_output_in_vitro.tsv") %>% 
  left_join(uniprot_ids, by = c("uniprot_id" = "Entry")) %>% 
  rename(symbol = From) %>% filter(symbol %in% missenses_in_vitro)

uniprot_with_domains <- uniprot %>% filter(!is.na(Domain)) %>% arrange(symbol)

# write_tsv(uniprot_with_domains %>% select(symbol, rsid, Domain), "results/2025.11.10.missense_snps_uniprot_annotation.tsv")

uniprot_ordered <- uniprot_with_domains %>% filter(Domain != "REGION|Disordered")
length(unique(uniprot_ordered$symbol)) 

####### FUNCTIONAL ANNOTATION
# background 
go_degs <- go(degs_in_vitro)
kegg_degs <- kegg(degs_in_vitro)

go_missenses <- go(missenses_in_vitro)
kegg_missenses <- kegg(missenses_in_vitro)

# genes having high ddASE and associated through motifbreaking SNPs in their promoters with DE or missense TFs
ddase_genes <- ddmaf_df_filt_vtr %>%
  filter(rsid %in% (snpgenedf_vtr %>% filter(tf %in% union(degs_in_vitro, missenses_in_vitro)) %>% pull(SNP_id))) %>%
  pull(symbol) %>% unique() 

ddase_genes_ss <- ddmaf_df_filt_ss %>%
  filter(rsid %in% (snpgenedf %>% filter(tf %in% union(degs_conf, missense_conf)) %>% pull(SNP_id))) %>%
  pull(symbol) %>% unique() 

go_ddase <- go(ddase_genes)
kegg_ddase <- kegg(ddase_genes)

# write_tsv(go_degs, "results/2025.10.28.GO_degs.tsv")
# write_tsv(kegg_degs, "results/2025.10.28.KEGG_degs.tsv")
# # write_tsv(go_missenses, "results/2025.10.28.GO_missenses.tsv")
# write_tsv(kegg_missenses, "results/2025.10.28.KEGG_missenses.tsv")
# write_tsv(go_ddase, "results/2025.10.28.GO_ddASE_genes.tsv")
# # write_tsv(kegg_ddase, "results/2025.10.28.KEGG_ddASE_genes.tsv")


########## SUPPLEMENTARY TABLE
# choose significant ddase genes having promoter snps motifbreaking differentially expressed or missense tfs
tableee <- ddmaf_df_filt_vtr %>% select(rsid, symbol, patients, padj, log2ROR) %>%
  inner_join(snpgenedf_vtr, by = c("rsid" = "SNP_id"), relationship = "many-to-many") %>%
  filter(tf %in% union(degs_in_vitro, missenses_in_vitro)) %>%
  # collect all tfs motifbroken by any snp in the promoter of the gene
  mutate(tfs = list(unique(tf)), .by = symbol) %>% select(-rsid, -tf) %>% distinct() %>%
  # calculate mean p-value and log2ROR across all comparisons with padj < 0.1
  summarize(padj = mean(padj), log2ROR = mean(log2ROR), .by = c(symbol, tfs)) %>%
  rowwise() %>% mutate(n_tfs = length(tfs), n_degs = length(intersect(tfs, degs_in_vitro)), n_missenses = length(intersect(tfs, missenses_in_vitro)))

# same for in vivo subset
tableee_ss <- ddmaf_df_filt_ss %>% select(rsid, symbol, patients, padj, log2ROR) %>%
  inner_join(snpgenedf, by = c("rsid" = "SNP_id"), relationship = "many-to-many") %>%
  filter(tf %in% union(degs_conf, missense_conf)) %>%
  mutate(tfs = list(unique(tf)), .by = symbol) %>% select(-rsid, -tf) %>% distinct() %>%
  summarize(padj = mean(padj), log2ROR = mean(log2ROR), .by = c(symbol, tfs)) %>%
  rowwise() %>% mutate(n_tfs = length(tfs), n_degs = length(intersect(tfs, degs_conf)), n_missenses = length(intersect(tfs, missense_conf)))

# how many pairs of people have this TF as differentially expressed
# in vitro
degs_freq_vtr <- enframe(table(deseq_df_vtr$symbol), name = "tf", value = "times_deg") %>% 
  # divide by how many could have had it diffexpressed (i.e. had a common snp motifbreaking this tf)
  inner_join(enframe(table(comparisons_df_vtr$tf), name = "tf", value = "times_all")) %>% 
  filter(tf %in% degs_in_vitro) %>% mutate(freq_in_vitro = times_deg / times_all) 

# in vivo
degs_freq_vv <- enframe(table(deseq_ss$symbol), name = "tf", value = "times_deg") %>%
  inner_join(enframe(table(comparisons_ss$tf), name = "tf", value = "times_all")) %>% 
  filter(tf %in% degs_conf) %>% mutate(freq_in_vivo = times_deg / times_all)

deg_tfs_freq_table <- full_join(degs_freq_vtr, degs_freq_vv, by = "tf", suffix = c("_vtr", "_vv")) %>% 
  arrange(desc(freq_in_vitro))

# missenses
mssns_freq_vtr <- result_missenses_vtr %>%
  left_join(distinct(uniprot[, c("rsid", "Domain")]), relationship = "many-to-many") %>%
  # divide by 4 if the missense is heterozygous
  mutate(value = ifelse(qual == "hetero", 0.25, 1)) %>% 
  # divide again is the mutation is non-functional (?)
  mutate(value = ifelse(is.na(Domain) | grepl("REGION|Disordered", Domain), value / 4, value)) %>% 
  # sum the values for each tf (all rsids and all comparisons (may be greater than 3))
  summarize(sum_value = sum(value), .by = tf) %>% 
  # divide by how many comparisons have been made (gene-wise)
  inner_join(enframe(table(comparisons_df_vtr$tf), name = "tf", value = "times_all")) %>% 
  mutate(freq_in_vitro = sum_value / times_all) %>% arrange(desc(freq_in_vitro))

mssns_freq_vv  <- missense_ss_df %>%
  filter(tf %in% missense_conf) %>% 
  left_join(distinct(uniprot[, c("rsid", "Domain")]), relationship = "many-to-many") %>%
  mutate(value = ifelse(qual == "hetero", 0.25, 1)) %>% 
  mutate(value = ifelse(is.na(Domain) | grepl("REGION|Disordered", Domain), value / 4, value)) %>% 
  summarize(sum_value = sum(value), .by = tf) %>% 
  inner_join(enframe(table(comparisons_ss$tf), name = "tf", value = "times_all")) %>% 
  mutate(freq_in_vivo = sum_value / times_all) %>% arrange(desc(freq_in_vivo))

mssns_tfs_freq_table <- full_join(mssns_freq_vtr, mssns_freq_vv, by = "tf", suffix = c("_vtr", "_vv"))

# degs and missenses in vitro and in vivo
ress <- full_join(deg_tfs_freq_table[, c("tf", "freq_in_vitro", "freq_in_vivo")],
                  mssns_tfs_freq_table[, c("tf", "freq_in_vitro", "freq_in_vivo")],
                  by = "tf", suffix = c("_deg", "_mssns"))

tablee <- tableee %>% unnest(tfs) %>% rename(tf = tfs) %>%
  left_join(ress[, c("tf", "freq_in_vitro_deg", "freq_in_vitro_mssns")]) %>%
  summarize(
    tfs = list(tf),
    freq_in_vitro_deg = sum(freq_in_vitro_deg, na.rm = TRUE),
    freq_in_vitro_mssns = sum(freq_in_vitro_mssns, na.rm = TRUE),
    .by = c(symbol, padj, log2ROR, n_tfs, n_degs, n_missenses)) %>%
  # divide the sum of degs' values by total number of degs (meaning if all of them were differentially expressed in all possible comparisons)
  mutate(deg_freq_vtr_norm = freq_in_vitro_deg / n_degs, mssns_freq_vtr_norm = freq_in_vitro_mssns / n_missenses)

tablee_ss <- tableee_ss %>% unnest(tfs) %>% rename(tf = tfs) %>%
  left_join(ress[, c("tf", "freq_in_vivo_deg", "freq_in_vivo_mssns")]) %>%
  summarize(tfs = list(tf), freq_in_vivo_deg = sum(freq_in_vivo_deg, na.rm = TRUE), freq_in_vivo_mssns = sum(freq_in_vivo_mssns, na.rm = TRUE),
            .by = c(symbol, padj, log2ROR, n_tfs, n_degs, n_missenses)) %>% 
  mutate(deg_freq_vv_norm = freq_in_vivo_deg / n_degs, mssns_freq_vv_norm = freq_in_vivo_mssns / n_missenses)

# calculate the total rank based on proportion of degs, missenses, mean p-value and odds-ratio difference
the_table <- tablee %>% 
  mutate(
    rank_padj = dense_rank(padj),
    rank_ror = dense_rank(-log2ROR),
    rank_degs = dense_rank(-deg_freq_vtr_norm),
    rank_mssns = dense_rank(-mssns_freq_vtr_norm)
  ) %>% rowwise() %>%
  mutate(
    total_rank = sum(rank_padj, rank_ror, rank_degs, rank_mssns, na.rm = TRUE),
    degs = paste(intersect(tfs, degs_in_vitro), collapse = ", "),
    missenses = paste(intersect(tfs, missenses_in_vitro), collapse = ", ")
  ) %>% ungroup() %>% arrange(total_rank) %>% select(-tfs)

the_table_ss <- tablee_ss %>% 
  mutate(
    rank_padj = dense_rank(padj),
    rank_ror = dense_rank(-log2ROR),
    rank_degs = dense_rank(-deg_freq_vv_norm),
    rank_mssns = dense_rank(-mssns_freq_vv_norm)
  ) %>% rowwise() %>%
  mutate(
    total_rank = sum(rank_padj, rank_ror, rank_degs, rank_mssns, na.rm = TRUE),
    degs = paste(intersect(tfs, degs_conf), collapse = ", "),
    missenses = paste(intersect(tfs, missense_conf), collapse = ", ")
  ) %>% ungroup() %>% arrange(total_rank) %>% select(-tfs)

final_in_vitro_table <- the_table %>%
  mutate(across(c(padj, log2ROR, freq_in_vitro_deg, freq_in_vitro_mssns, deg_freq_vtr_norm, mssns_freq_vtr_norm), ~ round(., 2))) %>%
  select(-freq_in_vitro_deg, -freq_in_vitro_mssns, -rank_padj, -rank_ror, -rank_degs, -rank_mssns)

final_in_vivo_table <- the_table_ss %>%
  mutate(across(c(padj, log2ROR, freq_in_vivo_deg, freq_in_vivo_mssns, deg_freq_vv_norm, mssns_freq_vv_norm), ~ round(., 2))) %>%
  select(-freq_in_vivo_deg, -freq_in_vivo_mssns, -rank_padj, -rank_ror, -rank_degs, -rank_mssns)