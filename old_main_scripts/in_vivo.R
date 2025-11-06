# all in vivo promoter chip-seq snps
patients_vv <- c("s1", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s12", "s15")
chipseq_tables <- lapply(patients_vv, read_chipseq)

chipseq_table <- bind_rows(chipseq_tables) %>% 
  mutate(variant = paste(chr, pos, ref, alt, sep = '_')) %>% 
  group_by(variant) %>% filter(n_distinct(patient) > 1) %>% ungroup() 

chipseq_table <- assign_genes(chipseq_table, "chr", "pos", promoters = TRUE)
chipseq_table <- add_rsids(chipseq_table, "chr", "pos", "ref")

# in vivo RNA-seq
ase_in_vivo <- lapply(patients_vv, function(x) {
  read_rnaseq(x, min_DP = 100)
})

ase_vv <- prepare_rna(ase_in_vivo)
df_collapsed_vv <- prepare_cs(chipseq_table ) 
ddmaf_df_vv <- compute_ddmaf(ase_vv, df_collapsed_vv)

topgenes_vv <- (ddmaf_df_vv %>% distinct(SYMBOL) %>% pull(SYMBOL))[1:100]
length(unique(ddmaf_df_vv[1:658, ]$SYMBOL)) # 101
ddmaf_df_vv[657, "ddMAF"] # 0.302
ddmaf_df_filt_vv <- ddmaf_df_vv[1:657, ]
# ddmaf_df_filt <- ddmaf_df %>% filter(dd_ref_frac >= 0.136)

snps_for_mtfbrkr <- unique(unlist(ddmaf_df_filt$rsids))
length(snps_for_mtfbrkr) # 408

mtfbrkr_dir <- "/media/leon/DISK2/icig/done/motifbreakr_16_09_25"

result_files <- list.files(mtfbrkr_dir, pattern = ".*\\.rds$") 
rsids_in_folder <- gsub("motifbreak_", "", gsub("\\.rds", "", result_files))
rsids_to_process <- setdiff(snps_for_mtfbrkr, rsids_in_folder)
# writeLines(rsids_to_process, "/media/leon/DISK2/icig/done/snps_for_mtfbrkr_to_add.txt")

combined_results <- do.call(c, unname(all_results))
hcmc <- combined_results[combined_results$dataSource == "HOCOMOCOv13"]

snpgenedf <-  hcmc %>% mcols() %>% as.data.frame() %>% 
  select(geneSymbol, SNP_id) %>% 
  rename(TF = geneSymbol) %>%
  distinct()


##################### дифэкспрессия
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

# список всех генов которые экспрессируются вообще для функциональной аннотации и тп
dds <- DESeqDataSetFromMatrix(countData = counts_vv,
                              colData = coldata_vv,
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized = TRUE) >= 5 ) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- results(dds)

annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = rownames(res),
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "ENSEMBL")

annotations_ordered <- annotations[match(rownames(res), annotations$ENSEMBL), ]
annotations_clean <- annotations[!is.na(annotations$ENTREZID) & !is.na(annotations$SYMBOL), ]
all_symbols <- annotations_ordered$SYMBOL

TFs <- unique(snpgenedf$TF)
length(unique(snpgenedf$SNP_id)) # 401
length(TFs) # 764
length(intersect(TFs, all_symbols)) 

# mapping hocomoco TF names to gene symbols
json_lines <- readLines("/media/leon/DISK2/icig/done/H13CORE_annotation.jsonl")

tf_gene_dict <- map_dfr(json_lines, function(line) {
  data <- fromJSON(line)
  
  if (!is.null(data$masterlist_info$species$HUMAN$gene_symbol)) {
    data.frame(
      tf = data$tf,
      motif_name = data$name,
      gene_symbol = data$masterlist_info$species$HUMAN$gene_symbol,
      stringsAsFactors = FALSE
    )
  } else { NULL }
})

tf_to_gene <- setNames(tf_gene_dict$gene_symbol, tf_gene_dict$tf)
tf_to_gene["ZN821"] <- "ZNF821"; tf_to_gene["ZN704"] <- "ZNF704"

snpgenedf <- snpgenedf %>%
  mutate(TF = recode(TF, !!!tf_to_gene))

TFs <- unique(snpgenedf$TF)
length(TFs) # 764
length(intersect(TFs, all_symbols)) # 648

TF_counts_vv <- counts_vv %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  filter(SYMBOL %in% TFs)

rownames(TF_counts_vv) <- TF_counts_vv$SYMBOL
TF_counts_vv <- TF_counts_vv %>% select(-ENSEMBL, -ENTREZID, -SYMBOL)
dds <- DESeqDataSetFromMatrix(countData = TF_counts_vv,
                              colData = coldata_vv,
                              design = ~ patient) # time + patient?
dds <- DESeq(dds)
res <- results(dds)

diffase_snp_tf <- ddmaf_df_filt %>%
  select(SYMBOL, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id") 

rs_gene_ptnt_tf <- ddmaf_df_filt %>%   
  select(SYMBOL, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id") 

cmbntns <- unique(rs_gene_ptnt_tf$patients)

TF_deseq_results <- data.frame()
for (pair in cmbntns) {
  p1 <- gsub("(s\\d+)s\\d+", "\\1", pair)
  p2 <- gsub("s\\d+(s\\d+)", "\\1", pair)
  
  res_pairwise <- as.data.frame(results(dds, contrast = c("patient", p1, p2))) %>%
    mutate(SYMBOL = rownames(.)) %>%
    mutate(contrast = paste0(p1, "_vs_", p2)) %>%
    filter(SYMBOL %in% rs_gene_ptnt_tf$TF[rs_gene_ptnt_tf$patients == pair]) 

  TF_deseq_results <- bind_rows(TF_deseq_results, res_pairwise)
}

TF_deseq_results <- TF_deseq_results %>% 
  filter(padj < 0.1 & abs(log2FoldChange) > 0.5) 
all_degs <- sort(unique(TF_deseq_results$SYMBOL))
length(all_degs) # 150
writeLines(all_degs, "~/in_vivo_initial_degs.txt")
length(unique(read_tsv("~/full_degs.txt", col_names = F)$X1)) # 121
length(intersect(all_degs, read_tsv("~/full_degs.txt", col_names = F)$X1)) # 48


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

## bedtools intersect -a /media/leon/Polina/atac_rna/dbSnp153.bed -b ~/in_vivo/exons_for_all_tfs.bed | awk '{print $1, $2, $3, $4}' - > ~/in_vivo/tf_exon_rsids_153.bed

snps_gr <- read_table("~/in_vivo/tf_exon_rsids_153.bed", col_names = F)
the_rsids <- unique(snps_gr$X4)
length(the_rsids) # 1245398

snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
snps_cons <- readRDS("/media/leon/DISK2/icig/biomart4.RDS")
rsids_left <- setdiff(the_rsids, snps_cons$refsnp_id)
for (i in seq(1, length(rsids_left), by = 10000)) {
  snps_cons <- rbind(snps_cons,
                     getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "consequence_type_tv"),
                           filters = 'snp_filter',
                           values = the_rsids[i:(i+10000)],
                           mart = snp_mart)
  )
  print(i)
}
# saveRDS(snps_cons, "/media/leon/DISK2/icig/biomart5.RDS")

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
mssns_snps$ENTREZID <- mssns_genes
mssns_snps$SYMBOL <- mapIds(org.Hs.eg.db,
                            keys = mssns_snps$ENTREZID,
                            column = "SYMBOL",
                            keytype = "ENTREZID")

mssns_snps <- mssns_snps %>% 
  filter(SYMBOL %in% TFs) %>%  select(-ENTREZID) %>%  distinct() %>%
  rename(TF = SYMBOL)

mssns_snps$TF <- unlist(mssns_snps$TF)

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

vcf <- read.vcfR("/media/leon/DISK2/icig/done/alignments/eight_people_in_vivo.filtered.vcf.gz")

refalt <- as.data.frame(vcf@fix) %>%
  select(CHROM, POS, REF, ALT) %>%
  rename(chr = CHROM, pos = POS) %>%
  separate_rows(ALT, sep = ",")

gt <- as.data.frame(extract.gt(vcf)) %>%
  rownames_to_column("variant") %>%
  separate(variant, into = c("chr", "pos", "n"), sep = "_", remove = TRUE) %>%
  select(-n) %>% distinct() %>% 
  mutate(pos = as.integer(pos)) %>% 
  inner_join(mssns_snps, by = c("chr" = "chr_name", "pos" = "chrom_start")) %>%
  mutate(pos = as.character(pos)) %>% 
  left_join(refalt, relationship = "many-to-many") %>%
  filter(str_length(REF) == 1, str_length(ALT) == 1) %>%
  rename(s3 = RNASEQ_s3_0R, s4 = RNASEQ_s4_0R, s5 = RNASEQ_s5_0R, s6 = RNASEQ_s6_0R, s8 = RNASEQ_s8_0R, s9 = RNASEQ_s9_0R, s12 = RNASEQ_s12_0R, s15 = RNASEQ_s15_0R) %>% 
  unnest_longer(TF)

diffase_snp_tf <- ddmaf_df_filt %>%
  select(SYMBOL, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id")

result_missenses <- diffase_snp_tf %>%
  select(-SNP_id, -SYMBOL) %>% distinct() %>% 
  inner_join(add_comparisons(patients_vv), by = "patients") %>%
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
length(missense_tfs) # 217
length(unique(na.omit(gt)$TF)) # 237
length(unique(missense_tfs)) / length(unique(gt$TF)) # 0.88
writeLines(missense_tfs, "~/in_vivo/missense_tf.txt")

# writeLines(union(result_missenses$refsnp_id, result_missenses_vv$refsnp_id), "~/missense_snps.txt")
