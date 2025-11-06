set.seed(988482)

library(readr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(motifbreakR)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyr)
library(jsonlite)
library(purrr)
library(MBASED)
library(DESeq2)
library(biomaRt)
library(GenomicRanges)
library(rtracklayer)
library(vcfR)
library(dplyr)

########### SNP из чипсеков
dm <- read_tsv("/media/leon/Polina/diabet/merge/Damarov.vcf", show_col_types = FALSE)
an <- read_tsv("/media/leon/Polina/diabet/merge/E.Viacheslavovna.vcf", show_col_types = FALSE)
bn <- read_tsv("/media/leon/Polina/diabet/merge/N.Petrovna.vcf", show_col_types = FALSE)

dm$patient <- "D"; an$patient <- "A"; bn$patient <- "B"

dm$variant <- paste(dm$chr, dm$pos, sep = '_')
an$variant <- paste(an$chr, an$pos, sep = '_')
bn$variant <- paste(bn$chr, bn$pos, sep = '_')

all(all(colnames(dm) == colnames(an)), all(colnames(dm) == colnames(bn)))
joined_df <- rbind(an, dm, bn)
joined_df$chr <- as.character(as.numeric(joined_df$chr)) # removing X Y MT and alternative contigs
joined_df <- na.omit(joined_df)

# выбираем SNP которые есть хотя бы у 2 из 3 людей
joined_df <- joined_df %>%
  group_by(variant) %>% filter(n() > 1) %>% ungroup() %>%
  mutate(alt = ifelse(allel1 != ref, allel1, allel2)) %>%
  dplyr::select(-DP, -QS1, -QS2, -allel1, -allel2)

# пересекаем с промоторами
hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
prmtrs <- promoters(hg38, columns = c("gene_id"))
exns <- exons(hg38, columns = c("exon_id", "gene_id"))
# blacklist <- read_tsv("/media/leon/DISK2/icig/ENCFF356LFX.bed", col_names = c("chr", "start", "end"))
# blacklist <- GRanges(
#   seqnames = blacklist$chr,
#   ranges = IRanges(start = blacklist$start + 1, end = blacklist$end)  # BED is 0-based, GRanges is 1-based
# )
# snp_ranges <- GRanges(
#   seqnames = paste0("chr", joined_df$chr),
#   ranges = IRanges(start = joined_df$pos, width = 1)
# )
# 
# overlaps_cs_blacklist <- findOverlaps(snp_ranges, blacklist) ## excluding ENCODE blacklisted regions (?)
# if (length(overlaps_cs_blacklist) > 0) {
#   blacklist_queries <- unique(queryHits(overlaps_cs_blacklist))
#   joined_df <- joined_df[-blacklist_queries, , drop = FALSE]
#   snp_ranges <- snp_ranges[-blacklist_queries]
# }

overlaps <- findOverlaps(snp_ranges, prmtrs)
prmrt_genes <- prmtrs$gene_id[subjectHits(overlaps)]
joined_df$ENTREZID <- NA
joined_df$ENTREZID[queryHits(overlaps)] <- prmrt_genes
joined_df$SYMBOL <- mapIds(org.Hs.eg.db,
                           keys = as.character(joined_df$ENTREZID),
                           column = "SYMBOL",
                           keytype = "ENTREZID")
joined_df <- na.omit(joined_df)

add_rsids <- function(chr_vec, pos_vec, ref_vec) {
  snp_ranges <- GRanges(seqnames = as.character(chr_vec),
                        IRanges(pos_vec, width = 1), ref = ref_vec)
  matched_snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, snp_ranges)
  df_matches <- as.data.frame(matched_snps)[, c("seqnames", "pos")]
  df_matches$rsid <- mcols(matched_snps)$RefSNP_id
  df_matches <- na.omit(df_matches)
  return(df_matches)
}

joined_df <- merge(joined_df, add_rsids(joined_df$chr, joined_df$pos, joined_df$ref), by.x = c("chr", "pos"), by.y = c("seqnames", "pos"))

######### in vitro RNA-seq - выбираем гетерозиготные SNP с покрытием > 100 и привязываем к генам
stats_dir <- "/media/leon/DISK2/icig/done/all_stats/"

ASE_in_vitro <- function(patient, min_DP = 100) {
  if (!is.na(as.numeric(substr(patient, 2, 3)))) {
    stats0R <- read_delim(paste0(stats_dir, patient, "_0R.stat"), show_col_types = FALSE) %>%
      filter(DP_total >= min_DP)
    stats1R <- read_delim(paste0(stats_dir, patient, "_1R.stat"), show_col_types = FALSE) %>%
      filter(DP_total >= min_DP)
  } else {
    stats0R <- read_tsv(paste0(stats_dir, patient, "_DR.stat"), show_col_types = FALSE) %>%
      filter(DP_total >= min_DP)
    stats1R <- read_tsv(paste0(stats_dir, patient, "_ER.stat"), show_col_types = FALSE) %>%
      filter(DP_total >= min_DP)
  }
  stats <- merge(stats0R, stats1R, by = c("CHR", "POS"))
  colnames(stats) <- c("chr", "pos", "id_sample_0R", "DP_0R", "ref_0R", "QS_ref_0R", "alt1_0R", "QS_alt1_0R", "alt2_0R", "QS_alt2_0R", "id_sample_1R", "DP_1R", "ref_1R", "QS_ref_1R", "alt1_1R", "QS_alt1_1R", "alt2_1R", "QS_alt2_1R")
  stats$chr <- as.character(as.numeric(stats$chr)) # removing X Y MT and alternative contigs
  
  stats <- stats[stringr::str_length(stats$ref_0R) == 1, ]
  stats <- stats[stringr::str_length(stats$alt1_0R) == 1, ]
  stats <- stats[stringr::str_length(stats$alt1_1R) == 1, ]
  
  stats[stats$alt1_0R != stats$alt1_1R, c(7, 8, 9, 10, 15, 16, 17, 18)] <- stats[stats$alt1_0R != stats$alt1_1R, c(7, 8, 9, 10, 17, 18, 15, 16)]
  
  gr_positions <- GRanges(
    seqnames = paste0("chr", stats$chr),
    ranges = IRanges(start = stats$pos, end = stats$pos)
  )
  # overlaps_cs_cnv <- findOverlaps(gr_positions, cnv_ranges) ## excluding ENCODE blacklisted regions (?)
  # 
  # if (length(overlaps_cs_cnv) > 0) {
  #   cnv_queries <- unique(queryHits(overlaps_cs_cnv))
  #   stats <- stats[-cnv_queries, , drop = FALSE]
  #   gr_positions <- gr_positions[-cnv_queries]
  # }
  overlaps <- findOverlaps(gr_positions, exns)
  stats <- stats[queryHits(overlaps), ]
  snp_exons <- exns[subjectHits(overlaps)]
  gene_ids <- sapply(snp_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
  stats$ENTREZID <- gene_ids
  stats$SYMBOL <- mapIds(org.Hs.eg.db,
                         keys = stats$ENTREZID,
                         keytype = "ENTREZID",
                         column = "SYMBOL")
  stats <- stats %>%
    unnest(SYMBOL) %>%
    filter(!grepl("HLA-", SYMBOL)) %>%
    unique() %>% na.omit() %>%
    mutate(patient = gsub("._.*", "", gsub("RNASEQ_s", "", id_sample_0R)),
           ref = ref_0R, alt = alt1_0R) %>%
    dplyr::select(-contains("id_sample_"), -contains("alt2"), -ref_0R, -ref_1R, -alt1_0R, -alt1_1R) %>%
    mutate(log2FC_1vs0 = log2(QS_ref_1R * QS_alt1_0R / QS_alt1_1R / QS_ref_0R))
  
  if (!is.na(as.numeric(substr(patient, 2, 3)))) { stats$patient = patient }
  return(stats)
}

patients_vtr <- c("sAn", "sBn", "sDm")
ase_in_vitro <- lapply(patients_vtr, function(x) {
  ASE_in_vitro(x, min_DP = 100)
})

ase <- bind_rows(ase_in_vitro)

# выбираем такие гены у которых есть и SNP в экзонах, и SNP в промоторах у 2+ человек
ase <- ase %>% filter(SYMBOL %in% intersect(ase$SYMBOL, joined_df$SYMBOL))
joined_df_filt <- joined_df %>% filter(SYMBOL %in% intersect(ase$SYMBOL, joined_df$SYMBOL))
length(unique(joined_df_filt$SYMBOL)) # 3377 

joined_df_filt$variant <- paste(joined_df_filt$variant, joined_df_filt$ref, joined_df_filt$alt, sep = '_')
joined_df_filt <- joined_df_filt %>% dplyr::select(-ENTREZID)

#################### для всех генов с SNP в экзонах с достаточным покрытием считаем diffASE на уровне гена
ase$variant <- paste(ase$chr, ase$pos, sep = '_')
ase <- ase %>%
  mutate(DP_ref_0R = round(DP_0R * QS_ref_0R),
         DP_alt_0R = round(DP_0R * QS_alt1_0R),
         DP_ref_1R = round(DP_1R * QS_ref_1R),
         DP_alt_1R = round(DP_1R * QS_alt1_1R)) %>%
  filter(DP_ref_0R > 10 & DP_alt_0R > 10 & DP_ref_1R > 10 & DP_alt_1R > 10)

summarizeASEResults_2s <- function(MBASEDOutput) {
  geneOutputDF <- data.frame(
    majorAlleleFrequencyDifference = assays(MBASEDOutput)$majorAlleleFrequencyDifference[, 1],
    pValueASE = assays(MBASEDOutput)$pValueASE[, 1],
    pValueHeterogeneity = assays(MBASEDOutput)$pValueHeterogeneity[, 1]
  )
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[, 1]
  lociOutputGR$MAFDifference <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAFDifference[, 1]
  lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels = unique(lociOutputGR$aseID)))
  return(
    list(
      geneOutput = geneOutputDF,
      locusOutput = lociOutputList
    )
  )
}  
run_mbased_analysis <- function(ase_data, patient_id, num_sim = 10^6, workers = 20) {
  
  ase_sub <- ase_data %>%
    filter(patient == patient_id) %>%
    group_by(SYMBOL) %>%
    group_modify(~ slice_sample(.x, n = min(15, nrow(.x)))) %>%
    ungroup()
  
  snps <- GRanges(
    seqnames = ase_sub$chr,
    ranges = IRanges(ase_sub$pos, width = 1),
    aseID = ase_sub$SYMBOL,
    allele1 = ase_sub$ref,
    allele2 = ase_sub$alt
  )
  names(snps) <- paste(ase_sub$SYMBOL, ase_sub$variant, sep = "_")  
  
  run_comparison <- function(comparison_type) {
    if (comparison_type == "0vs1") {
      
      se <- SummarizedExperiment(
        assays = list(
          lociAllele1Counts = matrix(
            c(
              ase_sub[["DP_ref_0R"]],
              ase_sub[["DP_ref_1R"]]
            ),
            ncol = 2,
            dimnames = list(
              names(snps),
              c("0", "1")
            )
          ),
          lociAllele2Counts = matrix(
            c(
              ase_sub[["DP_alt_0R"]],
              ase_sub[["DP_alt_1R"]]
            ),
            ncol = 2,
            dimnames = list(
              names(snps),
              c("0", "1")
            )
          )
        ),
        rowRanges = snps
      )
    } else {
      
      se <- SummarizedExperiment(
        assays = list(
          lociAllele1Counts = matrix(
            c(
              ase_sub[["DP_ref_1R"]],
              ase_sub[["DP_ref_0R"]]
            ),
            ncol = 2,
            dimnames = list(
              names(snps),
              c("1", "0")
            )
          ),
          lociAllele2Counts = matrix(
            c(
              ase_sub[["DP_alt_1R"]],
              ase_sub[["DP_alt_0R"]]
            ),
            ncol = 2,
            dimnames = list(
              names(snps),
              c("1", "0")
            )
          )
        ),
        rowRanges = snps
      )
    }
    
    bp <- MulticoreParam(workers = workers)
    start_time <- Sys.time()
    results <- runMBASED(
      ASESummarizedExperiment = se,
      isPhased = FALSE,
      numSim = num_sim,
      BPPARAM = bp
    )
    end_time <- Sys.time()
    print(paste("Time for", patient_id, comparison_type, ":", round(end_time - start_time, 1)))
    
    # Process results
    gene_output <- summarizeASEResults_2s(results)$geneOutput
    gene_output$patient <- patient_id
    gene_output$SYMBOL <- rownames(gene_output)
    names(gene_output)[names(gene_output) == "majorAlleleFrequencyDifference"] <- "genewiseDeltaMAF"
    
    locus_df <- data.frame(
      snp_id = names(unlist(summarizeASEResults_2s(results)$locusOutput)),
      as.data.frame(unlist(summarizeASEResults_2s(results)$locusOutput), row.names = NULL)) %>%
      merge(gene_output[, c("SYMBOL", "genewiseDeltaMAF", "pValueASE", "pValueHeterogeneity")], 
            by.x = "aseID", by.y = "SYMBOL")
    
    # Merge with original data
    full_df <- merge(locus_df, ase_sub,
                     by.x = c("seqnames", "start", "aseID"),
                     by.y = c("chr", "pos", "SYMBOL")) %>%
      mutate(comp = comparison_type)
    
    return(full_df)
  }
  
  # Run both comparisons
  df_0vs1 <- run_comparison("0vs1")
  df_1vs0 <- run_comparison("1vs0")
  
  # Post-processing for each comparison type
  process_comparison <- function(comp_df) {
    
    full_results <- comp_df %>%
      select(-c(snp_id, end, width, strand, QS_ref_0R, QS_alt1_0R, QS_ref_1R, QS_alt1_1R, DP_0R, DP_1R)) %>%
      rename(SYMBOL = aseID)
    
    genewise_results <- full_results %>%
      select(SYMBOL, genewiseDeltaMAF, pValueASE, pValueHeterogeneity, patient) %>%
      distinct()
    
    return(list(full = full_results, genewise = genewise_results))
  }
  
  # Process both comparisons
  results_0vs1 <- process_comparison(df_0vs1)
  results_1vs0 <- process_comparison(df_1vs0)
  
  # Return results as a list
  return(list(
    full_0vs1 = results_0vs1$full,
    genewise_0vs1 = results_0vs1$genewise,
    full_1vs0 = results_1vs0$full,
    genewise_1vs0 = results_1vs0$genewise
  ))
}

# Run for all patients
results_A <- run_mbased_analysis(ase, "A")
results_B <- run_mbased_analysis(ase, "B")
results_D <- run_mbased_analysis(ase, "D")

# # Combine all results for each comparison type
full_mbased_0vs1 <- bind_rows(results_A$full_0vs1, results_B$full_0vs1, results_D$full_0vs1) ## не нужны
full_mbased_1vs0 <- bind_rows(results_A$full_1vs0, results_B$full_1vs0, results_D$full_1vs0)

genewise_mbased_0vs1 <- bind_rows(results_A$genewise_0vs1, results_B$genewise_0vs1, results_D$genewise_0vs1) %>%
  mutate(padj = p.adjust(pValueASE, method = "BH")) %>%
  mutate(diffASE = padj < 0.01)

genewise_mbased_1vs0 <- bind_rows(results_A$genewise_1vs0, results_B$genewise_1vs0, results_D$genewise_1vs0) %>%
  mutate(padj = p.adjust(pValueASE, method = "BH")) %>%
  mutate(diffASE = padj < 0.01)

##################################
# для каждого гена из RNA-seq находим какие SNP из чипсеков (общие у 2+ людей) есть в промоторе этого гена
df_collapsed <- joined_df_filt %>%
  group_by(chr, pos, ref, alt, variant, SYMBOL, rsid) %>%
  summarise(patients = list(unique(patient)), .groups = "drop")

calculate_ddMAF <- function(genewise_mbased_df, chipseq_df, collapsed_cs_df) {
  ddmaf_df <- merge(genewise_mbased_df, chipseq_df, by = "SYMBOL", suffixes = c(".RNA", ".CS")) %>%
    select(-pValueHeterogeneity) %>% group_by(SYMBOL) %>% filter(any(diffASE)) %>% ungroup() %>% # есть хотя бы один человек с diffASE 
    group_by(across(-patient.CS)) %>% summarize(patient.CS = list(unique(patient.CS)), .groups = "drop") %>% distinct() %>% 
    select(SYMBOL, patient.RNA, genewiseDeltaMAF, padj) %>% group_by(SYMBOL) %>% filter(n() >= 2) %>% 
    group_modify(~ { 
      patients <- .x$patient.RNA 
      deltaMAFs <- .x$genewiseDeltaMAF 
      pValues <- .x$padj
      combn(length(patients), 2, simplify = FALSE) %>%
        map_dfr(function(pair) {
          i <- pair[1]
          j <- pair[2]
          patient_1 <- patients[i]
          patient_2 <- patients[j]
          deltaMAF_1 <- deltaMAFs[i]
          deltaMAF_2 <- deltaMAFs[j]
          pValue_1 <- pValues[i]
          pValue_2 <- pValues[j]
          
          common_snps <- collapsed_cs_df %>%
            filter(SYMBOL == .y$SYMBOL) %>%
            filter(map_lgl(patients, ~ all(c(patient_1, patient_2) %in% .x))) %>%
            pull(rsid)
          
          ddMAF <- if((length(common_snps) > 0) & (min(c(pValue_1, pValue_2)) < 0.01) & (patient_1 != patient_2)) abs(abs(deltaMAF_1) - abs(deltaMAF_2)) else NA
          tibble(
            patients = ifelse(patient_1 < patient_2, paste0(patient_1, patient_2), paste0(patient_2, patient_1)), 
            ddMAF = ddMAF,
            deltaMAF1 = deltaMAF_1,
            deltaMAF2 = deltaMAF_2,
            rsids = list(common_snps),
            n_common_snps = length(common_snps),
            padj1 = pValues[i],
            padj2 = pValues[j] 
            )
          })
      }) %>% ungroup() %>% distinct() %>% filter(n_common_snps > 0) #%>% na.omit()
}
ddmaf_0vs1 <- calculate_ddMAF(genewise_mbased_0vs1, joined_df_filt, df_collapsed)
ddmaf_1vs0 <- calculate_ddMAF(genewise_mbased_1vs0, joined_df_filt, df_collapsed)

# выбираем 100 генов с наибольшим ddMAF
ddmaf_df <- bind_rows(ddmaf_0vs1 %>% mutate(comp = "0vs1"), ddmaf_1vs0 %>% mutate(comp = "1vs0")) %>%
  arrange(desc(ddMAF)) %>% distinct() 

topgenes <- (ddmaf_df %>% distinct(SYMBOL) %>% pull(SYMBOL))[1:100]
length(unique(ddmaf_df[1:202,]$SYMBOL)) # 101
ddmaf_df[201, "ddMAF"] # 0.121
ddmaf_df_filt <- ddmaf_df[1:201,]

snps_for_mtfbrkr <- unique(unlist(ddmaf_df_filt$rsids))
length(snps_for_mtfbrkr) # 271

setwd("/media/leon/DISK2/icig/done/motifbreakr_16_09_25/") ## сделать нормально
result_files <- list.files(pattern = "motifbreak_.*\\.rds$") 
rsids_in_folder <- gsub("motifbreak_", "", gsub(".rds", "", result_files))

writeLines(setdiff(snps_for_mtfbrkr, rsids_in_folder), "/media/leon/DISK2/icig/done/snps_for_mtfbrkr_to_add.txt")

########### motifbreakR

result_files <- result_files[rsids_in_folder %in% snps_for_mtfbrkr]

all_results <- list()
for (file in result_files) {
  res <- readRDS(file)
  rsid <- sub("motifbreak_(.*)\\.rds", "\\1", file)
  if (!rsid %in% snps_for_mtfbrkr) { next }
  if (length(res) > 0) {
    mcols(res)$rsid <- rsid
    all_results[[rsid]] <- res
  }
}

combined_results <- do.call(c, unname(all_results))
hcmc <- combined_results[combined_results$dataSource == "HOCOMOCOv13"]

snpgenedf <-  hcmc %>% mcols() %>% as.data.frame() %>% 
  select(geneSymbol, SNP_id) %>% 
  rename(TF = geneSymbol) %>%
  distinct()

########## дифэкспрессия генов ТФ
counts <- read.table("~/counts_in_vitro_17.09.tsv", header = TRUE, row.names = 1, sep = "\t")
counts <- counts[, c(6:11)]
colnames(counts) <- c("D0", "D1", "B0", "B1", "A0", "A1")
coldata <- data.frame(patient = substr(colnames(counts), 1, 1),
                      time = substr(colnames(counts), 2, 2)) %>%
  mutate(patient = as.factor(patient), time = as.factor(time))
rownames(coldata) <- colnames(counts)

# список всех генов которые экспрессируются вообще для функциональной аннотации и тп
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized = TRUE) >= 5 ) >= 2
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
length(unique(snpgenedf$SNP_id)) # 261
length(TFs) # 661
length(intersect(TFs, all_symbols)) # 317

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
tf_to_gene["ZN821"] <- "ZNF821"
tf_to_gene["ZN704"] <- "ZNF704"

snpgenedf <- snpgenedf %>%
  mutate(TF = recode(TF, !!!tf_to_gene))

TFs <- unique(snpgenedf$TF)
length(TFs) # 654
length(intersect(TFs, all_symbols)) # 545
setdiff(TFs, all_symbols)

TF_counts <- counts %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  filter(SYMBOL %in% TFs)

rownames(TF_counts) <- TF_counts$ENSEMBL
TF_counts <- TF_counts %>% select(-ENSEMBL, -ENTREZID, -SYMBOL)
dds <- DESeqDataSetFromMatrix(countData = TF_counts,
                              colData = coldata,
                              design = ~ patient) # добавить time чтобы учитывать только межиндивидуальную вариабельность а не эмпаглифлозин (??) 
dds <- DESeq(dds)
res <- results(dds)

###### выбираем попарно те дэги которые связаны с общими снп в промоторах генов по которым у этих двух людей разница в диффасе
diffase_snp_tf <- ddmaf_df_filt %>%
  select(SYMBOL, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id") 

res <- as.data.frame(results(dds)) %>% mutate(ENSEMBL = rownames(.)) %>% merge(annotations_ordered, by = "ENSEMBL")
length(unique(res$SYMBOL))

resAB <- as.data.frame(results(dds, contrast = c("patient", "A", "B"))) %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "AB") %>%
  filter(SYMBOL %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "AB"])

resAD <- as.data.frame(results(dds, contrast = c("patient", "A", "D"))) %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "AD") %>%
  filter(SYMBOL %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "AD"])

resBD <- as.data.frame(results(dds, contrast = c("patient", "B", "D"))) %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "BD") %>%
  filter(SYMBOL %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "BD"])

TF_deseq_results <- bind_rows(resAB, resAD, resBD) %>% filter(padj < 0.1 & abs(log2FoldChange) > 0.5) 
all_degs <- unique(TF_deseq_results$SYMBOL)
length(all_degs) # 117 # 108 если добавить time в дизайн
writeLines(all_degs, "~/final_list_of_degtfs.txt")

deg_tf_rsids <- ddmaf_df_filt %>% 
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id") %>%
  filter(TF %in% all_degs)

TF_genes_df <- deg_tf_rsids %>% group_by(TF) %>% summarize(genes = list(SYMBOL))

causal_tf_rsids <- ddmaf_df_filt %>% 
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id") %>%
  filter(TF %in% union(all_degs, missense_tfs))

# choosing snps in promoters of DE or missense tfs
causal_snps_df <- causal_tf_rsids %>% select(SNP_id) %>%
  merge(joined_df[, c("chr", "pos", "ref", "alt", "rsid")], by.x = "SNP_id", by.y = "rsid") %>%
  distinct()

write_tsv(causal_snps_df[, c('chr', 'pos')], '/media/leon/DISK2/icig/done/85_snps.tsv', col_names = F) # 242

# выбираем пациентов ин виво гетерозигот по снп связанным с дэгами или миссенс тфами
# cd /media/leon/DISK2/icig/chipseq_july/stats_for_promoter_snps/; for i in 3 4 5 6 8 9 12 15; do perl /media/leon/DISK2/icig/done/nextflow/bin/stats.pl s${i} ../alignments/CHIPSEQ_s${i}.bam; done

check_snps <- function(path) {
  ddf <- read_tsv(path, show_col_types = FALSE)
  colnames(ddf) <- c("patient", "chr", "pos", "dp", "ref", "QS_ref", "alt", "QS_alt", "alt2", "QS_alt2")
  if (max(ddf$QS_alt2) < 0.1) { ddf <- ddf %>% select(-alt2, -QS_alt2) } else { print("alarm") }
  ddf <- ddf %>%
    filter(stringr::str_length(ref) == 1) %>%
    filter(stringr::str_length(alt) == 1) %>%
    mutate(dp_ref = round(dp * QS_ref), dp_alt = round(dp * QS_alt)) %>%
    select(-QS_ref, -QS_alt, -dp)
  return(ddf)
}

# joining dfs with coverage of heterozygous snps motifbreaking DE or missense tfs in in vivo chipseq
chipseq_tables <- lapply(Sys.glob("/media/leon/DISK2/icig/chipseq_july/stats_for_promoter_snps/*.stat"), check_snps)
chipseq_table <- bind_rows(chipseq_tables) %>%
  mutate(chr = factor(chr)) %>% 
  left_join(joined_df[, c('chr', 'pos', 'SYMBOL', 'rsid')], by = c("chr", "pos")) %>% 
  distinct()

chipseq_table_tfs <- causal_tf_rsids %>%
  select(SNP_id, SYMBOL, TF) %>% distinct() %>%
  group_by(SNP_id, SYMBOL) %>% summarize(TFs = list(TF), .groups = 'drop') %>% # for each snp and corresponding diffase gene collecting tfs into a list
  rename(rsid = SNP_id) %>%
  merge(chipseq_table, by = 'rsid')

# считаем ddMAF для in vivo RNA-seq
patients_vv <- c("s3", "s4", "s5", "s6", "s8", "s9", "s12", "s15")
ase_in_vivo <- lapply(patients_vv, function(x) {
  ASE_in_vitro(x, min_DP = 100)
})

ase_vv <- bind_rows(ase_in_vivo)

ase_vv <- ase_vv %>%
  mutate(DP_ref_0R = round(DP_0R * QS_ref_0R),
         DP_alt_0R = round(DP_0R * QS_alt1_0R),
         DP_ref_1R = round(DP_1R * QS_ref_1R),
         DP_alt_1R = round(DP_1R * QS_alt1_1R)) %>%
  filter(DP_ref_0R > 10 & DP_alt_0R > 10 & DP_ref_1R > 10 & DP_alt_1R > 10)

ase_vv$variant <- paste(ase_vv$chr, ase_vv$pos, sep = '_')

# выбираем только те гены в промоторах которых и находятся итоговые SNP (== гены по которым была наибольшая ddMAF)
ase_vv_subset <- ase_vv %>%
  filter(SYMBOL %in% chipseq_table_tfs$SYMBOL)

results_s3 <- run_mbased_analysis(ase_vv_subset, "s3")
results_s4 <- run_mbased_analysis(ase_vv_subset, "s4")
results_s5 <- run_mbased_analysis(ase_vv_subset, "s5")
results_s6 <- run_mbased_analysis(ase_vv_subset, "s6")
results_s8 <- run_mbased_analysis(ase_vv_subset, "s8")
results_s9 <- run_mbased_analysis(ase_vv_subset, "s9")
results_s12 <- run_mbased_analysis(ase_vv_subset, "s12")
results_s15 <- run_mbased_analysis(ase_vv_subset, "s15")

full_mbased_vv_0vs1 <- bind_rows(results_s3$full_0vs1, results_s4$full_0vs1, results_s5$full_0vs1, results_s6$full_0vs1, results_s8$full_0vs1, results_s9$full_0vs1, results_s12$full_0vs1, results_s15$full_0vs1)
full_mbased_vv_1vs0 <- bind_rows(results_s3$full_1vs0, results_s4$full_1vs0, results_s5$full_1vs0, results_s6$full_1vs0, results_s8$full_1vs0, results_s9$full_1vs0, results_s12$full_1vs0, results_s15$full_1vs0)

genewise_mbased_vv_0vs1 <- bind_rows(results_s3$genewise_0vs1, results_s4$genewise_0vs1, results_s5$genewise_0vs1, results_s6$genewise_0vs1, results_s8$genewise_0vs1, results_s9$genewise_0vs1, results_s12$genewise_0vs1, results_s15$genewise_0vs1) %>%
  mutate(padj = p.adjust(pValueASE, method = "BH")) %>%
  mutate(diffASE = padj < 0.01)

genewise_mbased_vv_1vs0 <- bind_rows(results_s3$genewise_1vs0, results_s4$genewise_1vs0, results_s5$genewise_1vs0, results_s6$genewise_1vs0, results_s8$genewise_1vs0, results_s9$genewise_1vs0, results_s12$genewise_1vs0, results_s15$genewise_1vs0) %>%
  mutate(padj = p.adjust(pValueASE, method = "BH")) %>%
  mutate(diffASE = padj < 0.01)

df_collapsed_vv <- chipseq_table_tfs %>%
  group_by(chr, pos, ref, alt, SYMBOL, rsid) %>%
  summarise(patients = list(unique(patient)), .groups = "drop")

ddmaf_vv_0vs1 <- calculate_ddMAF(genewise_mbased_vv_0vs1, chipseq_table_tfs[, 1:8], df_collapsed_vv)
ddmaf_vv_1vs0 <- calculate_ddMAF(genewise_mbased_vv_1vs0, chipseq_table_tfs[, 1:8], df_collapsed_vv)

ddmaf_vv_df <- bind_rows(ddmaf_vv_0vs1 %>% mutate(comp = "0vs1"), ddmaf_vv_1vs0 %>% mutate(comp = "1vs0")) %>%
  arrange(desc(ddMAF)) %>% distinct()

ddmaf_df_filt_vv <- ddmaf_vv_df %>% filter(abs(ddMAF) > 0.121)
length(unique(ddmaf_df_filt_vv$SYMBOL)) # 15 из 100
length(unique(unlist(ddmaf_df_filt_vv$rsids))) # 23

TFs_to_check_invivo <- causal_tf_rsids %>%
  filter(SYMBOL %in% ddmaf_df_filt_vv$SYMBOL) %>% pull(TF) %>% unique()

length(TFs_to_check_invivo) # 108

########## DESEQ IN VIVO
counts_vv <- read.table("~/counts_in_vivo_17.09.tsv", header = TRUE, row.names = 1, sep = "\t")
colnames(counts_vv) <- gsub("X.media.leon.DISK2.icig.done.alignments.RNASEQ_", "", gsub("R.bam", "", colnames(counts_vv)))
counts_vv

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
  mutate(ENSEMBL = rownames(.)) %>%
  merge(annotations_ordered, by = "ENSEMBL") %>%
  filter(SYMBOL %in% TFs_to_check_invivo)

rownames(TF_counts_vv) <- TF_counts_vv$SYMBOL

TF_counts_vv <- TF_counts_vv %>% select(-ENTREZID, -SYMBOL, -ENSEMBL) 

ddsvv <- DESeqDataSetFromMatrix(countData = TF_counts_vv,
                                colData = coldata_vv,
                                design = ~ patient)
ddsvv <- estimateSizeFactors(ddsvv)
keepvv <- rowSums(counts(ddsvv, normalized = TRUE) >= 5 ) >= 2
ddsvv <- ddsvv[keepvv, ]
ddsvv <- DESeq(ddsvv)
resvv <- results(ddsvv)

rs_gene_ptnt_tf_vv <- ddmaf_df_filt_vv %>%   
  select(SYMBOL, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id") 

cmbntns <- unique(rs_gene_pat_tf_vv$patients)
TF_deseq_results_vv <- data.frame()
for (pair in cmbntns) {
  p1 <- gsub("(s\\d+)(s\\d+)", "\\1", pair)
  p2 <- gsub("(s\\d+)(s\\d+)", "\\2", pair)

 res_pairwise <- as.data.frame(results(ddsvv, contrast = c("patient", p1, p2))) %>%
  mutate(SYMBOL = rownames(.)) %>%
  mutate(contrast = paste0(p1, "_vs_", p2)) %>%
  filter(SYMBOL %in% rs_gene_pat_tf_vv$TF[rs_gene_pat_tf_vv$patients == pair]) 
 TF_deseq_results_vv <- bind_rows(TF_deseq_results_vv, res_pairwise)
}

TF_deseq_results_vv <- TF_deseq_results_vv %>% 
  filter(padj < 0.1 & abs(log2FoldChange) > 0.5) 
all_degs_vv <- unique(TF_deseq_results_vv$SYMBOL)
length(all_degs_vv) # 11

length(intersect(all_degs, all_degs_vv)) # 9
both_degs <- union(all_degs, all_degs_vv) # 119

deg_tf_rsids_vv <- ddmaf_df_filt_vv %>% # ?
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(snpgenedf, by = "SNP_id") %>%
  filter(TF %in% TF_deseq_results_vv$SYMBOL)

### write_tsv(TF_genes_df['TF'], '~/45_tfs.tsv', col_names = F)


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
  select(rsids, SYMBOL, patients, ddMAF, comp) %>%
  mutate(ddMAF = round(ddMAF, 3)) %>%
  group_by(rsids, SYMBOL, patients) %>%
  slice_max(ddMAF, n = 1, with_ties = TRUE) %>%
  ungroup() %>% select(-comp) %>% distinct() %>% 
  rename(rsid = rsids, diffase_gene = SYMBOL, ddmaf_in_vitro = ddMAF) %>% 
  inner_join(snpgenedf, by = c("rsid" = "SNP_id")) %>%
  group_by(rsid, diffase_gene, patients, ddmaf_in_vitro) %>% 
  reframe(tfs = list(unique(TF)), n_tfs = n_distinct(TF))

table1_vv <- ddmaf_vv_df %>% na.omit() %>% 
  unnest(rsids) %>%
  select(rsids, SYMBOL, patients, ddMAF, comp) %>%
  mutate(ddMAF = round(ddMAF, 3)) %>%
  group_by(rsids, SYMBOL, patients) %>%
  slice_max(ddMAF, n = 1, with_ties = TRUE) %>%
  ungroup() %>% select(-comp) %>% distinct() %>% 
  rename(rsid = rsids, diffase_gene = SYMBOL, ddmaf_in_vivo = ddMAF) %>% 
  inner_join(snpgenedf, by = c("rsid" = "SNP_id")) %>%
  group_by(rsid, diffase_gene, patients, ddmaf_in_vivo) %>% 
  reframe(tfs = list(unique(TF)), n_tfs = n_distinct(TF))

snp_patients <- joined_df_filt %>%
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
    n_fnctnl_mssnss = vapply(tfs, function(x) length(intersect(x, uniprot_ordered$SYMBOL)), integer(1))
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
  inner_join(table2[, c('rsid', 'diffase_gene')], by = c('rsid' = 'rsid', 'SYMBOL' = 'diffase_gene')) %>% 
  select(rsid, SYMBOL, Consequence, CLIN_SIG, PUBMED) %>% 
  distinct()
  