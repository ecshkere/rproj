library(tidyverse)

read_chipseq <- function(patient) {
  df <- read_tsv(paste0(stats_dir, patient, ".stat")) %>%
    rename(alt = alt1, patient = id_sample) %>% 
    filter(str_length(ref) == 1, str_length(alt) == 1,
           chr %in% as.character(seq(1, 22))) %>%
    mutate(DP_ref = round(DP * QS_ref), DP_alt = round(DP * QS_alt1)) %>%
    filter(DP >= 20, DP_ref / DP >= 0.1, DP_ref / DP <= 0.9) %>% 
    select(-c(QS_ref, QS_alt1, DP, alt2, QS_alt2))
  return(df)
}

read_rnaseq <- function(patient, min_DP = 100) {
  if (!is.na(as.numeric(substr(patient, 2, 3)))) {
    stats0R <- read_tsv(paste0(stats_dir, patient, "_0R.stat"))
    stats1R <- read_tsv(paste0(stats_dir, patient, "_1R.stat"))
  } else {
    stats0R <- read_tsv(paste0(stats_dir, patient, "_DR.stat"))
    stats1R <- read_tsv(paste0(stats_dir, patient, "_ER.stat"))
  }
  stats <- inner_join(stats0R, stats1R, by = c("chr", "pos"), suffix = c("_0R", "_1R"), relationship = "many-to-many") %>% 
    filter(DP_0R >= min_DP, DP_1R >= min_DP) %>%
    filter(chr %in% as.character(seq(1, 22))) %>% 
    filter(str_length(ref_0R) == 1, str_length(alt1_0R) == 1, str_length(alt1_1R) == 1) 
  
  stats[stats$alt1_0R == stats$alt2_1R, c("alt1_1R", "QS_alt1_1R", "alt2_1R", "QS_alt2_1R")] <- stats[stats$alt1_0R == stats$alt2_1R, c("alt2_1R", "QS_alt2_1R", "alt1_1R", "QS_alt1_1R")]  
  
  stats <- assign_genes(stats, "chr", "pos", exons_only = TRUE) %>%
    filter(!grepl("HLA-", symbol)) %>%
    mutate(patient = patient,
           ref = ref_0R, alt = alt1_0R,
           DP_ref_0R = round(DP_0R * QS_ref_0R),
           DP_alt_0R = round(DP_0R * QS_alt1_0R),
           DP_ref_1R = round(DP_1R * QS_ref_1R),
           DP_alt_1R = round(DP_1R * QS_alt1_1R)) %>%
    filter(DP_ref_0R / DP_0R > 0.1 & DP_alt_0R / DP_0R > 0.1 & DP_ref_1R / DP_1R > 0.1 & DP_alt_1R / DP_1R > 0.1) %>%
    mutate(log2FC_1vs0 = log2(QS_ref_1R * QS_alt1_0R / QS_alt1_1R / QS_ref_0R),
           ref_frac_diff = QS_ref_1R - QS_ref_0R) %>% 
    select(-c(contains("id_sample_"), contains("alt2"), ref_0R, ref_1R, alt1_0R, alt1_1R, contains("QS")))
  return(stats)
}

prepare_rna <- function(ase_dfs, chipseq_df) {
  ase <- ase_dfs %>% 
    filter(symbol %in% chipseq_df$symbol) %>%
    filter(n_distinct(patient) > 1, .by = symbol) %>% 
    mutate(
      snp = paste(chr, pos, sep = "_"),
      ref_major = DP_ref_0R > DP_alt_0R,
      major = ifelse(ref_major, ref, alt),
      minor = ifelse(ref_major, alt, ref),
      DP_maj_0R = ifelse(ref_major, DP_ref_0R, DP_alt_0R),
      DP_min_0R = ifelse(ref_major, DP_alt_0R, DP_ref_0R),
      DP_maj_1R = ifelse(ref_major, DP_ref_1R, DP_alt_1R),
      DP_min_1R = ifelse(ref_major, DP_alt_1R, DP_ref_1R)
    ) %>% mutate(
      n_snps = n_distinct(snp),
      snps = list(sort(unique(snp))),
      .by = c(symbol, patient)
    ) %>% select(-c(entrezid, DP_0R, DP_1R, log2FC_1vs0, ref_frac_diff))
  return(ase)
}

# для каждого гена из RNA-seq находим какие SNP из чипсеков у каких попарных комбинаций людей есть в промоторе этого гена
prepare_cs <- function(chipseq_dfs, ase) {
  chipseq_df <- chipseq_dfs %>%
    filter(symbol %in% ase$symbol) %>% select(-DP_ref, -DP_alt) %>% 
    mutate(variant = paste(chr, pos, ref, alt, sep = '_')) %>% 
    group_by(across(-patient)) %>% filter(n_distinct(patient) > 1) %>%
    mutate(patients = {
      combs <- combn(sort(unique(patient)), 2, simplify = FALSE)
      paste(sapply(combs, function(x) paste0(x, collapse = "_")), collapse = ",")
    }) %>% ungroup() %>% 
    separate_rows(patients, sep = ",") %>% distinct()
  return(chipseq_df)
}

#### at least one common snp (commented out)
## half the snps are common
prepare_for_glm <- function(ase, chipseq_df) {
  df_pairs <- chipseq_df %>%
    distinct(symbol, patients)

  valid_symbols <- ase %>% 
    filter(n_distinct(patient) > 1, .by = symbol) %>% pull(symbol)
  
  ase_sub <- ase %>% filter(symbol %in% valid_symbols)

  res <- ase_sub %>%
    inner_join(ase_sub, by = "symbol", suffix = c("_p1", "_p2"), relationship = "many-to-many") %>%
    rename(p1 = patient_p1, p2 = patient_p2) %>%
    filter(p1 < p2) %>%
    mutate(patients = paste(p1, p2, sep = "_")) %>%
    rowwise %>%
    mutate(n_common = length(intersect(snps_p1, snps_p2)),
           frac_common = round(n_common / length(union(snps_p1, snps_p2)), 2)) %>%
    ungroup() %>%
    inner_join(df_pairs) %>% distinct() %>%
    filter(case_when(
      # n_common > 0 ~ snp_p1 == snp_p2,
      # n_common == 0 ~ snp_p1 != snp_p2)) # drop unique snps for genes where at least one common snp is present
      
      frac_common >= 0.5 ~ snp_p1 == snp_p2, # drop unique snps for genes where at least half the snps are common 
      frac_common < 0.5 ~ ((snp_p1 == snp_p2) | (snp_p1 != snp_p2))), # otherwise keep them all but without being able to compare the direction of asymmetry shift later
      .by = c(symbol, patients)
    ) %>% 
    mutate(same_major = case_when( # swap DP_maj and DP_min for common snps if major_p2 differs from major_p1
      # snp_p1 == snp_p2 & major_p1 == major_p2 ~ "TRUE",
      # snp_p1 == snp_p2 & major_p1 == minor_p2 ~ "FALSE",
      # snp_p1 != snp_p2 ~ "NA")
      
      frac_common >= 0.5 & snp_p1 == snp_p2 & major_p1 == major_p2 ~ "T",
      frac_common >= 0.5 & snp_p1 == snp_p2 & major_p1 == minor_p2 ~ "F",
      frac_common < 0.5 ~ "NA") # will be calculating mean coverage for each patient individually
    ) %>%
    mutate(
      major_p2_tmp = ifelse(same_major == "F", minor_p2, major_p2),
      minor_p2_tmp = ifelse(same_major == "F", major_p2, minor_p2),
      DP_maj_0R_p2_tmp = ifelse(same_major == "F", DP_min_0R_p2, DP_maj_0R_p2),
      DP_min_0R_p2_tmp = ifelse(same_major == "F", DP_maj_0R_p2, DP_min_0R_p2),
      DP_maj_1R_p2_tmp = ifelse(same_major == "F", DP_min_1R_p2, DP_maj_1R_p2),
      DP_min_1R_p2_tmp = ifelse(same_major == "F", DP_maj_1R_p2, DP_min_1R_p2),
    ) %>%
    mutate(major_p2 = major_p2_tmp, minor_p2 = minor_p2_tmp, DP_maj_0R_p2 = DP_maj_0R_p2_tmp, DP_min_0R_p2 = DP_min_0R_p2_tmp, DP_maj_1R_p2 = DP_maj_1R_p2_tmp, DP_min_1R_p2 = DP_min_1R_p2_tmp) %>%
    select(-contains("_tmp")) %>%
    mutate(
      DP_maj_gene_0_p1 = round(mean(DP_maj_0R_p1)),
      DP_min_gene_0_p1 = round(mean(DP_min_0R_p1)),
      DP_maj_gene_1_p1 = round(mean(DP_maj_1R_p1)),
      DP_min_gene_1_p1 = round(mean(DP_min_1R_p1)),
      DP_maj_gene_0_p2 = round(mean(DP_maj_0R_p2)),
      DP_min_gene_0_p2 = round(mean(DP_min_0R_p2)),
      DP_maj_gene_1_p2 = round(mean(DP_maj_1R_p2)),
      DP_min_gene_1_p2 = round(mean(DP_min_1R_p2)),
      OR_p1 = (DP_maj_gene_1_p1 / DP_min_gene_1_p1) / (DP_maj_gene_0_p1 / DP_min_gene_0_p1),
      OR_p2 = (DP_maj_gene_1_p2 / DP_min_gene_1_p2) / (DP_maj_gene_0_p2 / DP_min_gene_0_p2),
      log2ROR = ifelse(
        # n_common == 0,
        frac_common < 0.5,
        abs(abs(log2(OR_p1)) - abs(log2(OR_p2))),
        abs(log2(OR_p1) - log2(OR_p2))), # ignore direction if snps were different
        .by = c(symbol, patients)
    ) %>% select(symbol, patients, n_common, frac_common, n_snps_p1, n_snps_p2, p1, OR_p1, OR_p2, log2ROR, DP_maj_gene_0_p1, DP_min_gene_0_p1, DP_maj_gene_1_p1, DP_min_gene_1_p1, p2, DP_maj_gene_0_p2, DP_min_gene_0_p2, DP_maj_gene_1_p2, DP_min_gene_1_p2) %>%
    distinct() %>% arrange(frac_common)

  return(res)
}

## logistic regression
run_glm <- function(df) {
  results <- lapply(seq_len(nrow(df)), function(i) {
    r <- df[i, ]
    patients <- unlist(strsplit(r$patients, "_"))
    p1 <- patients[1]; p2 <- patients[2]
    
    dat <- data.frame(
      patient = c(p1, p2, p1, p2),
      condition = factor(c(0, 0, 1, 1)),
      maj = c(r$DP_maj_gene_0_p1, r$DP_maj_gene_0_p2,
              r$DP_maj_gene_1_p1, r$DP_maj_gene_1_p2),
      min = c(r$DP_min_gene_0_p1, r$DP_min_gene_0_p2,
              r$DP_min_gene_1_p1, r$DP_min_gene_1_p2)
    )
    dat <- dat %>%
      mutate(total = maj + min) %>%
      pivot_longer(cols = c(maj, min), names_to = "allele_type", values_to = "count") %>%
      mutate(is_maj = ifelse(allele_type == "maj", 1, 0))
    
    fit <- stats::glm(is_maj ~ condition * patient,
                      family = binomial(),
                      weights = count,
                      data = dat)
    
    inter <- grep(":", names(coef(fit)), value = TRUE)
    if (length(inter) > 0) {
      tibble(
        symbol = r$symbol,
        patients = r$patients,
        p_interaction = summary(fit)$coefficients[inter[length(inter)], 4],
        coef_interaction = coef(fit)[inter[length(inter)]]
      )
    } else { NULL }})
  
  bind_rows(results) %>% drop_na(p_interaction)
}
