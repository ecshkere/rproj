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

prepare_rna <- function(ase_dfs, chipseq_df, p = 0.1) {
  ase <- ase_dfs %>% 
    filter(symbol %in% chipseq_df$symbol) %>%
    filter(n_distinct(patient) > 1, .by = symbol) %>% 
    mutate(ref_major = DP_ref_0R > DP_alt_0R,
           major = ifelse(ref_major, ref, alt),
           minor = ifelse(ref_major, alt, ref),
           DP_maj_0R = ifelse(ref_major, DP_ref_0R, DP_alt_0R),
           DP_min_0R = ifelse(ref_major, DP_alt_0R, DP_ref_0R),
           DP_maj_1R = ifelse(ref_major, DP_ref_1R, DP_alt_1R),
           DP_min_1R = ifelse(ref_major, DP_alt_1R, DP_ref_1R)) %>% 
    summarize(
      DP_maj_gene_0 = round(mean(DP_maj_0R)),
      DP_min_gene_0 = round(mean(DP_min_0R)),
      DP_maj_gene_1 = round(mean(DP_maj_1R)),
      DP_min_gene_1 = round(mean(DP_min_1R)),
      dMAF = abs(DP_maj_gene_1 / (DP_maj_gene_1 + DP_min_gene_1) - DP_maj_gene_0 / (DP_maj_gene_0 + DP_min_gene_0)), 
      OR = (DP_maj_gene_1 / DP_min_gene_1) / (DP_maj_gene_0 / DP_min_gene_0),
      .by = c(symbol, patient))
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

prepare_for_glm <- function(ase, chipseq_df) {
  df_pairs <- chipseq_df %>%
    distinct(symbol, patients)
  
  valid_symbols <- ase %>%
    filter(n_distinct(patient) > 1, .by = symbol) %>% pull(symbol)
  
  ase_sub <- ase %>% filter(symbol %in% valid_symbols)
  
  pairs <- ase_sub %>% 
    inner_join(ase_sub, by = "symbol", suffix = c("_p1", "_p2"), relationship = "many-to-many") %>%
    rename(p1 = patient_p1, p2 = patient_p2) %>% 
    filter(p1 < p2) %>% 
    mutate(patients = paste(p1, p2, sep = "_"),
           log2OR1 = log2(OR_p1), log2OR2 = log2(OR_p2),
           log2ROR = abs(abs(log2(OR_p1)) - abs(log2(OR_p2))),
           ddMAF = abs(dMAF_p1 - dMAF_p2))
  
  res <- pairs %>%
    inner_join(df_pairs, by = c("symbol", "patients"), relationship = "many-to-many") %>%
    distinct() %>% arrange(desc(ddMAF))
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
