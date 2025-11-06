library(biomaRt)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(vcfR)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
attrs <- c("external_gene_name", "chromosome_name", "exon_chrom_start", "exon_chrom_end")

exons_df <- getBM(attributes = attrs,
                  filters = "external_gene_name",
                  values = setdiff(TFs, read_tsv("~/tfs_for_missense_prev.tsv", col_names = F)$X1),
                  mart = ensembl) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name),
         start = exon_chrom_start,
         end = exon_chrom_end) %>%
  select(chromosome_name, exon_chrom_start, exon_chrom_end)

write.table(exons_df, "~/exons_for_all_tfs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## bedtools intersect -a /media/leon/Polina/atac_rna/dbSnp153.bed -b ~/exons_for_all_tfs.bed > ~/tf_exon_rsids_153.bed
## awk '{print $1, $2, $3, $4}' ~/tf_exon_rsids_153.bed > ~/tf_exon_rsids_153clean.bed

snps_gr <- read_table("~/tf_exon_rsids_153clean.bed", col_names = F)
the_rsids <- unique(snps_gr$X4)
length(the_rsids) # 1025126

snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
snps_cons <- readRDS("/media/leon/DISK2/icig/biomart_third.RDS")
for (i in seq(1, length(the_rsids), by = 10000)) {
  snps_cons <- rbind(snps_cons,
  getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "consequence_type_tv"),
        filters = 'snp_filter',
        values = the_rsids[i:(i+10000)],
        mart = snp_mart)
  )
  print(i)
}

# snps_cons <- readRDS("/media/leon/DISK2/icig/biomart.RDS")
# saveRDS(snps_cons, "/media/leon/DISK2/icig/biomart645.RDS")
# saveRDS(snps_cons, "/media/leon/DISK2/icig/biomart_third.RDS")

snps_missense <- snps_cons %>% 
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

mssns_snps_bed <- mssns_snps %>%
  mutate(end = chrom_start) %>%
  mutate(chrom_start = chrom_start - 1) %>%
  select(chr_name, chrom_start, end) %>% distinct() 

old_missense_bed <- read_tsv("~/in_vitro_missense_snps.bed", col_names = F)
mssns_snps_bed <- mssns_snps_bed %>% 
  filter(!paste(chr_name, chrom_start, end) %in% paste(old_missense_bed$X1, old_missense_bed$X2, old_missense_bed$X3))

write.table(mssns_snps_bed, "~/in_vitro_missense_snps.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# awk 'NF==3 && $1!="NA"' ~/in_vitro_missense_snps.bed > ~/in_vitro_missense_snps.clean.bed
# bcftools mpileup -R  ~/in_vitro_missense_snps.clean.bed -f /media/leon/Polina/Genomes/hg38.chromFa/hg38_FOR_HISAT.fa RNASEQ_sBn.merged.bam RNASEQ_sAn.merged.bam  RNASEQ_sDm.merged.bam -d 10000 | bcftools call -mv -Oz -o final_in_vitro_missense_genotypes.vcf
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

vcf <- read.vcfR("/media/leon/DISK2/icig/done/alignments/in_vitro/in_vitro_merged.vcf.gz")

refalt <- as.data.frame(vcf@fix) %>%
  select(CHROM, POS, REF, ALT) %>%
  rename(chr = CHROM, pos = POS) %>%
  separate_rows(ALT, sep = ",")

gt <- as.data.frame(extract.gt(vcf)) %>%
  tibble::rownames_to_column("variant") %>%
  separate(variant, into = c("chr", "pos"), sep = "_", remove = FALSE) %>%
  select(-variant) %>% distinct() %>%
  merge(mssns_snps, by.x = c("chr", "pos"), by.y = c("chr_name", "chrom_start")) %>%
  left_join(refalt, by = c("chr", "pos"), relationship = "many-to-many") %>%
  filter(stringr::str_length(REF) == 1) %>%
  filter(stringr::str_length(ALT) == 1) %>%
  rename(B = sBn, A = sAn, D = sDm)

diffase_snp_tf <- ddmaf_df_filt %>%
  select(SYMBOL, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id")

ptnts <- c("A", "B", "D")

result_missenses <- diffase_snp_tf %>%
  select(-SNP_id, -SYMBOL) %>% distinct() %>% # 1248
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

missense_tfs <- unique(result_missenses$TF)
length(missense_tfs) # 200
length(unique(gt$TF)) # 236
length(unique(missense_tfs)) / length(unique(gt$TF)) # 0.85
writeLines(missense_tfs, "~/new_in_vitro_missenses.txt")

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

vcf_vv <- read.vcfR("/media/leon/DISK2/icig/done/alignments/in_vivo_missenses.more.vcf")
vcf_vv <- read.vcfR("/media/leon/DISK2/icig/done/vcfs/unfiltered/merged_samples.dp20.vcf.gz")

refalt_vv <- as.data.frame(vcf_vv@fix) %>%
  select(CHROM, POS, REF, ALT) %>%
  rename(chr = CHROM, pos = POS) %>%
  separate_rows(ALT, sep = ",")
# 
# ########
# gt_vv <- as.data.frame(extract.gt(vcf_vv)) %>%
#   tibble::rownames_to_column("variant") %>%
#   separate(variant, into = c("chr", "pos"), sep = "_", remove = FALSE) %>%
#   select(-variant) %>% distinct() %>%
#   mutate(pos = as.integer(pos)) %>% 
#   full_join(in_vivo_missense_snps, by = c("chr" = "chr_name", "pos" = "chrom_start")) %>%
#   left_join(refalt_vv, by = c("chr", "pos"), relationship = "many-to-many") %>%
#   filter(stringr::str_length(REF) == 1) %>%
#   filter(stringr::str_length(ALT) == 1) %>%
#   rename(s3 = RNASEQ_s3_0R, s4 = RNASEQ_s4_0R, s5 = RNASEQ_s5_0R, s6 = RNASEQ_s6_0R, s8 = RNASEQ_s8_0R, s9 = RNASEQ_s9_0R, s12 = RNASEQ_s12_0R, s15 = RNASEQ_s15_0R)
# 
# ##########



gt_vv <- as.data.frame(extract.gt(vcf_vv)) %>%
  tibble::rownames_to_column("variant") %>%
  separate(variant, into = c("chr", "pos"), sep = "_", remove = FALSE) %>%
  select(-variant) %>% distinct() %>%
  merge(in_vivo_missense_snps, by.x = c("chr", "pos"), by.y = c("chr_name", "chrom_start")) %>%
  left_join(refalt_vv, by = c("chr", "pos"), relationship = "many-to-many") %>%
  filter(stringr::str_length(REF) == 1) %>%
  filter(stringr::str_length(ALT) == 1) %>%
  rename(s3 = RNASEQ_s3_0R, s4 = RNASEQ_s4_0R, s5 = RNASEQ_s5_0R, s6 = RNASEQ_s6_0R, s8 = RNASEQ_s8_0R, s9 = RNASEQ_s9_0R, s12 = RNASEQ_s12_0R, s15 = RNASEQ_s15_0R)

diffase_snp_tf_vv <- ddmaf_df_filt_vv %>%
  select(SYMBOL, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id")

result_missenses_vv <- diffase_snp_tf_vv %>%
  select(-SNP_id, -SYMBOL) %>% distinct() %>% # 365
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

vv_missenses <- unique(result_missenses_vv$TF)
length(vv_missenses)
both_missenses <- union(vv_missenses, missense_tfs)

writeLines(union(result_missenses$refsnp_id, result_missenses_vv$refsnp_id), "~/missense_snps.txt")

### https://www.uniprot.org/id-mapping/
uniprot_ids <- read_tsv("/media/leon/DISK2/icig/done/uniprot/idmapping_2025_10_03.tsv") %>% 
  filter(Reviewed == "reviewed") %>% 
  select(From, Entry)

bind_rows(result_missenses[, c("TF", "refsnp_id")],
          result_missenses_vv[, c("TF", "refsnp_id")]) %>%
  distinct() %>% 
  left_join(uniprot_ids, by = c("TF" = "From")) %>% 
  select(Entry, refsnp_id) %>% 
  write_tsv("/media/leon/DISK2/icig/done/uniprot/tf_snp_df.tsv")

# vep_table <- read_tsv("~/Downloads/vep_423.txt") %>% 
#   rename(rsid = `#Uploaded_variation`)
# 
# vep <- vep_table %>% 
#   select(rsid, SYMBOL, Protein_position, UNIPROT_ISOFORM) %>% 
#   filter(UNIPROT_ISOFORM != "-") %>% distinct() %>%
#   mutate(isoform_number = gsub(".*-(.*)", "\\1", UNIPROT_ISOFORM)) %>%
#   filter(isoform_number == "1") %>% 
#   mutate(UNIPROT_ISOFORM = gsub("(.*)-.*", "\\1", UNIPROT_ISOFORM)) %>% 
#   filter(Protein_position != '-') %>% 
#   select(-isoform_number)
# 
# write_tsv(vep, "/media/leon/DISK2/icig/done/uniprot/uniprot_isoforms_to_parse.tsv")
# 
# length(unique(both_missenses)) # 209
# length(intersect(both_missenses, uniprot$SYMBOL)) # 118
# length(unique(uniprot %>% filter(!is.na(Domain)) %>% pull(SYMBOL))) # 88
# got_uniprot_id <- unique((vep_table %>% filter(Consequence == "missense_variant") %>%  filter(UNIPROT_ISOFORM != "-"))$SYMBOL) # 121

uniprot <- read_tsv("/media/leon/DISK2/icig/done/uniprot/uniprot_isoforms_with_domains_2.tsv") %>% 
  left_join(uniprot_ids, by = c("uniprot_id" = "Entry")) %>% 
  rename(SYMBOL = From)

length(unique(both_missenses)) # 209
length(intersect(both_missenses, uniprot$SYMBOL)) # 204

uniprot_with_domains <- uniprot %>% 
  filter(!is.na(Domain))

length(unique(uniprot_with_domains$SYMBOL)) # 154

uniprot_ordered <- uniprot_with_domains %>% 
  filter(Domain != "REGION|Disordered")

length(unique(uniprot_ordered$SYMBOL)) # 123
