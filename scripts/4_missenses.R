read_missense <- function(patient) {
  df <- read_tsv(paste0("data/", patient, ".mssns.tsv")) %>%
    filter(str_length(ref) == 1, str_length(alt1) == 1,
           as.character(chr) %in% as.character(seq(1, 22))) %>% 
    mutate(gt = case_when(QS_ref > 0.9 ~ "0/0",
                          QS_alt1 > 0.9 ~ "1/1",
                          QS_ref > 0.1 & QS_ref < 0.9 & QS_alt1 > 0.1 & QS_alt1 < 0.9 ~ "0/1",
                          .default = NA_character_),
           chr = as.character(chr)) %>% 
    rename(alt = alt1, patient = id_sample) %>% 
    select(patient, chr, pos, ref, alt, gt) %>% 
    mutate(chr = as.character(chr))
  return(df)
}

join_missense <- function(patients) {
  lapply(patients, read_missense) %>% 
    reduce(full_join, by = c("chr", "pos", "ref", "alt")) %>% 
    select(-contains("patient")) %>% 
    set_names(c("chr", "pos", "ref", "alt", patients)) %>%
    filter(rowSums(is.na(.)) <= length(patients) - 2)
}

add_comparisons <- function(patients_vec) {
  comparison_map <- combn(patients_vec, 2, simplify = FALSE) %>%
    map_df(~ tibble(
      patients = paste0(sort(c(.x[1], .x[2])), collapse = '_'),
      patient1 = .x[1],
      patient2 = .x[2]
    ))
  return(comparison_map)
}

translate_gt <- function(gt, ref, alt) {
  case_when(
    gt == "0/0" ~ paste0(ref, ref),
    gt == "0/1" ~ paste0(ref, alt),
    gt == "1/1" ~ paste0(alt, alt),
    gt == "./." ~ NA_character_,
    TRUE ~ NA_character_
  )
}

find_missense_tfs <- function(patients, tfs, patients_tf) {
  missense_gts <- join_missense(patients)

  missense_gts <- assign_genes(missense_gts, "chr", "pos", exons_only = TRUE) %>% 
   filter(symbol %in% tfs) %>% rename(tf = symbol)

  missense_gts <- add_rsids(missense_gts, "chr", "pos", "ref")

  result_missenses <- patients_tf %>% 
    inner_join(add_comparisons(patients)) %>%
    left_join(missense_gts, relationship = "many-to-many") %>%
    filter(!is.na(rsid)) %>% 
    rowwise() %>%
    mutate(allele1 = get(patient1), allele2 = get(patient2)) %>%
    filter(!is.na(allele1), !is.na(allele2)) %>%
    filter(allele1 != allele2) %>%
    # filter((allele1 == "0/0" & allele2 == "1/1") | (allele2 == "0/0" & allele1 == "1/1")) %>%
    ungroup() %>%
    select(tf, rsid, chr, pos, ref, alt, patient1, allele1, patient2, allele2) %>%
    mutate(allele1 = translate_gt(allele1, ref, alt),
           allele2 = translate_gt(allele2, ref, alt)) %>% 
    mutate(qual = ifelse(substr(allele1, 1, 1) == substr(allele1, 2, 2) & substr(allele2, 1, 1) == substr(allele2, 2, 2), "homo", "hetero"))
  return(result_missenses)
}
