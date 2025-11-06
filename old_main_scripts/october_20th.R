source("scripts/1_functions.R")
source("scripts/2_snp_stats_21.10.R")
source("scripts/3_motifbreakr.R")
# source("scripts/4_deseq.R")
# source("scripts/5_missenses.R")

## in vitro chip-seq snps
patients_vtr <- c("sAn", "sBn", "sDm")
joined_df <- lapply(patients_vtr, read_chipseq) %>% bind_rows() %>% 
  mutate(variant = paste(chr, pos, ref, alt, sep = '_')) %>% 
  group_by(variant) %>% filter(n_distinct(patient) > 1) %>% ungroup() # selecting snps present in at least 2 out of 3 people

# selecting snps in promoters
joined_df <- assign_genes(joined_df, "chr", "pos", promoters = TRUE)
joined_df <- add_rsids(joined_df, "chr", "pos", "ref") 

#### in vitro RNA-seq ------------------------------------------------------------------------
stats_dir <- "data/"
p <- 0.1
ase_in_vitro <- lapply(patients_vtr, function(x) {
  read_rnaseq(x, min_DP = 100)
})

ase <- prepare_rna(ase_in_vitro, joined_df)
df_collapsed <- prepare_cs(joined_df, ase)
ddmaf_df <- compute_ddmaf(ase, df_collapsed)

# выбираем 100 генов с наибольшим модулем разницы изменения доли референсного аллеля
length(unique(ddmaf_df$symbol)) # 257
length(unique(ddmaf_df[1:128, ]$symbol)) # 101
ddmaf_df[[127, "ddMAF"]] # 0.145
ddmaf_df_filt <- ddmaf_df[1:127, ] %>% 
  inner_join(df_collapsed %>% select(symbol, patients, rsid))

snps_for_mtfbrkr <- unique(ddmaf_df_filt$rsid)
length(snps_for_mtfbrkr) # 228

#### motifbreakR 
mtfbrkr_dir <- "motifbreakr"
result_files <- list.files(mtfbrkr_dir, pattern = ".*\\.rds$") 
rsids_in_folder <- gsub("\\.rds", "", result_files)
rsids_to_process <- setdiff(snps_for_mtfbrkr, rsids_in_folder)

## writeLines(rsids_to_process, "snps_for_mtfbrkr_to_add.txt")
## Rscript motifbreakr.R snps_for_mtfbrkr_to_add.txt 

if (length(rsids_to_process)) {
  # creates an RDS object for each SNP in the "./motifbreakr" directory
  run_motifbreakr(rsids_to_process, batch_size = 31) 
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

snpgenedf <- hcmc %>% mcols() %>% as.data.frame() %>% 
  select(geneSymbol, SNP_id) %>% rename(TF = geneSymbol) %>% distinct() %>% 
  mutate(TF = recode(TF, !!!tf_to_gene))

TFs <- unique(snpgenedf$TF)
length(TFs) # 615
length(unique(snpgenedf$SNP_id)) # 220

source("scripts/4_deseq.R")

deg_tf_rsids <- ddmaf_df_filt_vv %>% 
  inner_join(snpgenedf, by = c("rsid" = "SNP_id")) %>%
  filter(TF %in% all_degs)

TF_genes_df <- deg_tf_rsids %>%
  group_by(TF) %>% summarize(genes = list(symbol))

source("scripts/5_missenses.R")

## choosing snps in promoters of differentially expressed or missense tfs

# causal_tf_rsids <- ddmaf_df_filt %>% 
#   inner_join(snpgenedf, by = c("rsid" = "SNP_id")) %>%
#   filter(TF %in% union(all_degs, missense_tfs))
# 
# causal_snps_df <- causal_tf_rsids %>% select(rsid) %>%
#   inner_join(joined_df[, c("chr", "pos", "ref", "alt", "rsid")]) %>%
#   distinct()

# write_tsv(causal_snps_df[, c('chr', 'pos')], 'output/snps_to_.tsv', col_names = F) # 260

causal_snps <- snpgenedf %>% filter(TF %in% union(all_degs, missense_tfs)) %>% distinct(SNP_id) %>% pull(SNP_id)
writeLines(causal_snps, "causal_snps.tsv")




# ################
# joined_df_ss <- chipseq_table %>% filter(rsid %in% causal_snps)
################













# выбираем пациентов ин виво гетерозигот по снп связанным с дэгами или миссенс тфами
# cd /media/leon/DISK2/icig/chipseq_july/stats_for_promoter_snps/; for i in 3 4 5 6 8 9 12 15; do perl /media/leon/DISK2/icig/done/nextflow/bin/stats.pl s${i} ../alignments/CHIPSEQ_s${i}.bam; done

check_snps <- function(path) {
  ddf <- read_tsv(path)
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
uniprot_ids <- read_tsv("~/in_vivo/idmapping_2025_10_14.tsv") %>% 
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

snp_patients <- joined_df %>% filter(symbol %in% ase$symbol) %>%
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
