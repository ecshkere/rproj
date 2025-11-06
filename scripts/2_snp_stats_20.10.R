library(tidyverse)

## FOR IN VITRO CHIP-SEQ
read_old_chipseq <- function(patient) {
  df <- read_tsv(paste0("data/", patient, ".stat"), show_col_types = FALSE) %>%
    mutate(patient = patient,
           alt = ifelse(allel1 == ref, allel2, allel1), 
           dp_ref = ifelse(allel1 == ref, round(DP * QS1), round(DP * QS2)), 
           dp_alt = ifelse(allel1 == ref, round(DP * QS2), round(DP * QS1))) %>% 
    filter(str_length(ref) == 1, str_length(alt) == 1,
           chr %in% as.character(seq(1, 22))) %>%
    select(-c(QS1, QS2, DP, allel1, allel2)) %>% 
    filter(dp_ref + dp_alt >= 20, dp_ref / (dp_ref + dp_alt) >= 0.1, dp_ref / (dp_ref + dp_alt) <= 0.9)
  return(df)
}

## FOR IN VIVO CHIP-SEQ
read_chipseq <- function(sample) {
  df <- read_tsv(paste0("data/", sample, "C.stat"), show_col_types = FALSE) %>%
    rename(alt = alt1, QS_alt = QS_alt1) %>% 
    filter(str_length(ref) == 1, str_length(alt) == 1,
           chr %in% as.character(seq(1, 22))) %>%
    mutate(dp_ref = round(DP * QS_ref), dp_alt = round(DP * QS_alt), 
           condition = gsub("s\\d+_", "", id_sample),
           patient = gsub("_\\d+C", "", id_sample)) %>%
    select(-c(QS_ref, QS_alt, DP, alt2, QS_alt2)) %>% 
    group_by(patient, chr, pos, ref, alt) %>%
    summarize(dp_ref = sum(dp_ref), dp_alt = sum(dp_alt), .groups = 'drop') %>% 
    filter(dp_ref + dp_alt >= 20, dp_ref / (dp_ref + dp_alt) >= 0.1, dp_ref / (dp_ref + dp_alt) <= 0.9)
  return(df)
}

## FOR RNA-SEQ
read_stats <- function(patient, min_DP = 100) {
  if (!is.na(as.numeric(substr(patient, 2, 3)))) {
    stats0R <- read_tsv(paste0(stats_dir, patient, "_0R.stat"), show_col_types = FALSE)
    stats1R <- read_tsv(paste0(stats_dir, patient, "_1R.stat"), show_col_types = FALSE)
  } else {
    stats0R <- read_tsv(paste0(stats_dir, patient, "_DR.stat"), show_col_types = FALSE)
    stats1R <- read_tsv(paste0(stats_dir, patient, "_ER.stat"), show_col_types = FALSE)
  }
  stats <- inner_join(stats0R, stats1R, by = c("chr", "pos"), suffix = c("_0R", "_1R"), relationship = "many-to-many") %>% 
    filter(DP_0R >= min_DP, DP_1R >= min_DP) %>%
    # filter(DP_ref_0R > 5, DP_alt_0R > 5, DP_ref_1R > 5, DP_alt_1R > 5) %>%
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

process_rna_stats <- function(ase_dfs, chipseq_df, p = 0.01) {
  # ase_dfs <- ase_in_vitro; chipseq_df <- joined_df; p = 0.01
  ase_dfs <- ase_in_vivo; chipseq_df <- chipseq_table; p = 0.01
  ase <- bind_rows(ase_dfs) %>% 
    filter(symbol %in% chipseq_df$symbol) %>%
    mutate(total_DP = DP_0R + DP_1R) %>% 
    group_by(symbol, chr, pos) %>%  # choosing the one most common snp for each gene across patients
    mutate(n = n(), total_across_dp = sum(total_DP)) %>% 
    filter(n > 1, n_distinct(patient) > 1) %>% ungroup() %>% 
    rowwise() %>% # calculating chi-square and adjusting p-value
    mutate(p01 = chisq.test(matrix(c(DP_ref_0R, DP_alt_0R, DP_ref_1R, DP_alt_1R), nrow = 2, byrow = TRUE))$p.value) %>%
    ungroup() %>% drop_na() %>% 
    mutate(padj = p.adjust(p01, method = "BH")) %>% 
    mutate(diffASE = padj < p) %>%
    group_by(symbol, chr, pos) %>% 
    mutate(patients = {
      combs <- combn(sort(unique(patient)), 2, simplify = FALSE)
      paste(sapply(combs, function(x) paste0(x, collapse = "")), collapse = ",")
    }) %>% ungroup() %>% 
    separate_rows(patients, sep = ",") %>%
    distinct()
  
  # для каждого гена из RNA-seq находим какие SNP из чипсеков у каких попарных комбинаций людей есть в промоторе этого гена
  df_collapsed <- chipseq_df %>%
    filter(symbol %in% ase$symbol) %>% 
    mutate(variant = paste(chr, pos, ref, alt, sep = '_')) %>% 
    group_by(variant) %>% filter(n_distinct(patient) > 1) %>% ungroup() %>%
    mutate(patients = {
      combs <- combn(sort(unique(patient)), 2, simplify = FALSE)
      paste(sapply(combs, function(x) paste0(x, collapse = "")), collapse = ",")
    }) %>% ungroup() %>% 
    separate_rows(patients, sep = ",") %>%
    distinct()
  
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
          
          tibble(patients = ptnts,
                 dd_ref_frac = dd_ref_frac,
                 drf1 = drf1, drf2 = drf2,
                 padj1 = pvl1, padj2 = pvl2)
        })}) %>%
    ungroup() %>% distinct() %>%
    filter(!is.na(dd_ref_frac)) %>%
    arrange(desc(dd_ref_frac))
  
  return(list(ase, df_collapsed, ddmaf_df))
}