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
library(dplyr)
library(MBASED)
library(DESeq2)

add_rsids <- function(chr_vec, pos_vec, ref_vec) {
  snp_ranges <- GRanges(seqnames = as.character(chr_vec),
                        IRanges(pos_vec, width = 1), ref = ref_vec)
  matched_snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, snp_ranges)
  df_matches <- as.data.frame(matched_snps)[, c("seqnames", "pos")]
  df_matches$rsid <- mcols(matched_snps)$RefSNP_id
  df_matches <- na.omit(df_matches)
  return(df_matches)
}

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
joined_df$chr <- as.character(as.numeric(joined_df$chr))
joined_df <- na.omit(joined_df)

# выбираем SNP которые есть хотя бы у 2 из 3 людей
joined_df <- joined_df %>%
  group_by(variant) %>% filter(n() > 1) %>% ungroup() %>%
  mutate(alt = ifelse(allel1 != ref, allel1, allel2)) %>%
  select(-DP, -QS1, -QS2, -allel1, -allel2)

# пересекаем с промоторами
hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
prmtrs <- promoters(hg38, columns = c("gene_id"))
exns <- exons(hg38, columns = c("exon_id", "gene_id"))

snp_ranges <- GRanges(
  seqnames = paste0("chr", joined_df$chr),
  ranges = IRanges(start = joined_df$pos, width = 1)
)
overlaps <- findOverlaps(snp_ranges, prmtrs)
prmrt_genes <- prmtrs$gene_id[subjectHits(overlaps)]
joined_df$ENTREZID <- NA
joined_df$ENTREZID[queryHits(overlaps)] <- prmrt_genes
joined_df$SYMBOL <- mapIds(org.Hs.eg.db,
                           keys = as.character(joined_df$ENTREZID),
                           column = "SYMBOL",
                           keytype = "ENTREZID")
joined_df <- na.omit(joined_df)
joined_df <- merge(joined_df, add_rsids(joined_df$chr, joined_df$pos, joined_df$ref), by.x = c("chr", "pos"), by.y = c("seqnames", "pos"))

######### in vitro RNA-seq
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
  stats$chr <- as.character(as.numeric(stats$chr))
  
  stats <- stats[stringr::str_length(stats$ref_0R) == 1, ]
  stats <- stats[stringr::str_length(stats$alt1_0R) == 1, ]
  stats <- stats[stringr::str_length(stats$alt1_1R) == 1, ]
  
  stats[stats$alt1_0R != stats$alt1_1R, c(7, 8, 9, 10, 15, 16, 17, 18)] <- stats[stats$alt1_0R != stats$alt1_1R, c(7, 8, 9, 10, 17, 18, 15, 16)]
  
  stats$p01 <-  apply(stats[, c(4, 6, 8, 12, 14, 16)], 1, function(x) {
    data_matrix <- matrix(unlist(c(round(x[1] * x[2]), round(x[1] * x[3]), round(x[4] * x[5]), round(x[4] * x[6]))), nrow = 2, byrow = TRUE)
    v <- chisq.test(data_matrix)
    return(v$p.value)
  })
  stats <- na.omit(stats)
  
  gr_positions <- GRanges(
    seqnames = paste0("chr", stats$chr),
    ranges = IRanges(start = stats$pos, end = stats$pos)
  )
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
    select(-contains("id_sample_"), -contains("alt2"), -ref_0R, -ref_1R, -alt1_0R, -alt1_1R) %>%
    mutate(log2FC_1vs0 = log2(QS_ref_1R * QS_alt1_0R / QS_alt1_1R / QS_ref_0R))
  
  if (!is.na(as.numeric(substr(patient, 2, 3)))) { stats$patient = patient }
  return(stats)
}

patients_vtr <- c("sAn", "sBn", "sDm")
setwd("/media/leon/DISK2/icig/done/all_stats/")
ase_in_vitro <- lapply(patients_vtr, function(x) {
  ASE_in_vitro(x, min_DP = 100)
})

ase <- bind_rows(ase_in_vitro)

# выбираем такие гены у которых есть и SNP в экзонах, и SNP в промоторах у 2+ человек
ase <- ase %>% filter(SYMBOL %in% intersect(ase$SYMBOL, joined_df$SYMBOL))
joined_df_filt <- joined_df %>% filter(SYMBOL %in% intersect(ase$SYMBOL, joined_df$SYMBOL))
length(unique(joined_df_filt$SYMBOL)) # 3377 

joined_df_filt$variant <- paste(joined_df_filt$variant, joined_df_filt$ref, joined_df_filt$alt, sep = '_')

df_collapsed <- joined_df_filt %>%
  group_by(chr, pos, ref, alt, variant, ENTREZID, SYMBOL, rsid) %>%
  summarise(patients = paste0(sort(unique(patient)), collapse = ""), .groups = "drop")

#################### для выбранных генов с diffASE считаем ASE на уровне гена
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
  # Filter data for patient and sample 15 SNPs per gene
  ase_sub <- ase_data %>%
    filter(patient == patient_id) %>%
    group_by(SYMBOL) %>%
    group_modify(~ slice_sample(.x, n = min(15, nrow(.x)))) %>%
    ungroup()
  
  # Create GRanges object for SNPs
  snps <- GRanges(
    seqnames = ase_sub$chr,
    ranges = IRanges(ase_sub$pos, width = 1),
    aseID = ase_sub$SYMBOL,
    allele1 = ase_sub$ref,
    allele2 = ase_sub$alt
  )
  names(snps) <- paste(ase_sub$SYMBOL, ase_sub$variant, sep = "_")
  
  # Define comparison function
  run_comparison <- function(comparison_type) {
    # Create assay matrices based on comparison type
    if (comparison_type == "0vs1") {
      assay1 <- matrix(c(ase_sub[["DP_ref_0R"]], ase_sub[["DP_ref_1R"]]), ncol = 2)
      assay2 <- matrix(c(ase_sub[["DP_alt_0R"]], ase_sub[["DP_alt_1R"]]), ncol = 2)
      colnames(assay1) <- colnames(assay2) <- c("0", "1")
    } else {
      assay1 <- matrix(c(ase_sub[["DP_ref_1R"]], ase_sub[["DP_ref_0R"]]), ncol = 2)
      assay2 <- matrix(c(ase_sub[["DP_alt_1R"]], ase_sub[["DP_alt_0R"]]), ncol = 2)
      colnames(assay1) <- colnames(assay2) <- c("1", "0")
    }
    
    # Create SummarizedExperiment
    se <- SummarizedExperiment(
      assays = list(
        lociAllele1Counts = assay1,
        lociAllele2Counts = assay2
      ),
      rowRanges = snps
    )
    
    # Run MBASED
    bp <- MulticoreParam(workers = workers)
    start_time <- Sys.time()
    results <- runMBASED(
      ASESummarizedExperiment = se,
      isPhased = FALSE,
      numSim = num_sim,
      BPPARAM = bp
    )
    end_time <- Sys.time()
    print(paste("Time for", patient_id, comparison_type, ":", end_time, "\n", start_time))
    
    # Process results
    gene_output <- summarizeASEResults_2s(results)$geneOutput
    gene_output$patient <- patient_id
    gene_output$SYMBOL <- rownames(gene_output)
    names(gene_output)[names(gene_output) == "majorAlleleFrequencyDifference"] <- "genewiseMAFDiff"
    
    locus_df <- data.frame(
      snp_id = names(unlist(summarizeASEResults_2s(results)$locusOutput)),
      as.data.frame(unlist(summarizeASEResults_2s(results)$locusOutput), row.names = NULL)) %>%
      merge(gene_output[, c("SYMBOL", "genewiseMAFDiff", "pValueASE", "pValueHeterogeneity")], 
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
  process_comparison <- function(df, comp_type) {
    # Filter for the specific comparison type
    comp_df <- df %>% filter(comp == comp_type)
    
    # Create full results
    full_results <- comp_df %>%
      select(-c(snp_id, end, width, strand, MAFDifference, QS_ref_0R, QS_alt1_0R, QS_ref_1R, QS_alt1_1R, DP_0R, DP_1R, comp)) %>%
      rename(SYMBOL = aseID)
    
    # Create genewise results
    genewise_results <- full_results %>%
      select(SYMBOL, genewiseMAFDiff, pValueASE, pValueHeterogeneity, ENTREZID, patient) %>%
      distinct() # %>%
      # mutate(padj = p.adjust(pValueASE, method = "BH")) %>%
      # mutate(diffASE = padj < 0.01)
    
    return(list(full = full_results, genewise = genewise_results))
  }
  
  # Process both comparisons
  results_0vs1 <- process_comparison(bind_rows(df_0vs1, df_1vs0), "0vs1")
  results_1vs0 <- process_comparison(bind_rows(df_0vs1, df_1vs0), "1vs0")
  
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

# Combine all results for each comparison type
full_mbased_0vs1 <- bind_rows(results_A$full_0vs1, results_B$full_0vs1, results_D$full_0vs1)
genewise_mbased_0vs1 <- bind_rows(results_A$genewise_0vs1, results_B$genewise_0vs1, results_D$genewise_0vs1)

full_mbased_1vs0 <- bind_rows(results_A$full_1vs0, results_B$full_1vs0, results_D$full_1vs0)
genewise_mbased_1vs0 <- bind_rows(results_A$genewise_1vs0, results_B$genewise_1vs0, results_D$genewise_1vs0)

genewise_mbased_0vs1 <- genewise_mbased_0vs1 %>%
  mutate(padj = p.adjust(pValueASE, method = "BH")) %>%
  mutate(diffASE = padj < 0.01)

genewise_mbased_1vs0 <- genewise_mbased_1vs0 %>%
  mutate(padj = p.adjust(pValueASE, method = "BH")) %>%
  mutate(diffASE = padj < 0.01)

##################################
# для каждого гена из RNA-seq находим какие SNP из чипсеков (общие у 2+ людей) есть в промоторе этого гена
lfc_divided_0vs1 <- merge(genewise_mbased_0vs1, joined_df_filt, by = "SYMBOL", suffixes = c(".RNA", ".CS"), all.y = TRUE) %>%
  select(-c(ENTREZID.RNA, ENTREZID.CS, pValueHeterogeneity)) %>% 
  group_by(SYMBOL) %>% filter(any(diffASE)) %>% ungroup() %>% # есть хотя бы один человек с diffASE ######## НАДО ПРОВЕРЯТЬ КОГДА ПОПАРНО СЧИТАЮ ОТНОШЕНИЕ ЧТО ХОТЯ БЫ У ОДНОГО ЧЕЛОВЕКА ДОСТОВЕРНОЕ ДИФФАСЕ
  group_by(across(-patient.CS)) %>% 
  summarize(patient.CS = paste(sort(patient.CS), collapse = ""), .groups = "drop")%>% 
  select(SYMBOL, patient.RNA, genewiseMAFDiff, padj) %>% distinct() %>% 
  group_by(SYMBOL) %>% filter(n() >= 2) %>% 
  group_modify(~ {
    patients <- .x$patient.RNA
    MAFDiffs <- .x$genewiseMAFDiff
    pValues <- .x$padj
    combn(length(patients), 2, simplify = FALSE) %>%
      map_dfr(function(pair) {
        i <- pair[1]
        j <- pair[2]
        
        patient_1 <- patients[i]
        patient_2 <- patients[j]
        mafDiff_1 <- MAFDiffs[i]
        mafDiff_2 <- MAFDiffs[j]
        
        # считаем ddMAF только если имеется одинаковый SNP в промоторе данного гена у обоих людей И ХОТЯ БЫ У ОДНОГО ЧЕЛОВЕКА ДИФФАСЕ
        # у второго либо нет вообще либо есть но величина дельта-маф сильно меньше
        common_snps <- df_collapsed$rsid[df_collapsed$SYMBOL == .y$SYMBOL & (df_collapsed$patients == paste0(sort(unique(patients)), collapse = "") | df_collapsed$patients == "ABD")]
        mafDiffRatio <- if(length(common_snps) > 0) max(abs(c(mafDiff_1, mafDiff_2))) / min(abs(c(mafDiff_1, mafDiff_2))) else NA
        ddMAF <- if((length(common_snps) > 0) & (min(pValues) < 0.01)) abs(abs(mafDiff_1) - abs(mafDiff_2)) else NA
        
        tibble(
          patients = ifelse(patient_1 < patient_2, paste0(patient_1, patient_2), paste0(patient_2, patient_1)),
          mafDiffRatio = mafDiffRatio,
          ddMAF = ddMAF,
          logMafDiffRatio = if(!is.na(mafDiffRatio)) log2(mafDiffRatio) else NA,
          deltaMAF1 = mafDiff_1,
          deltaMAF2 = mafDiff_2,
          padj_1 = pValues[i],  padj_2 = pValues[j],
          rsids = list(common_snps),
          n_common_snps = length(common_snps)
        )
      })
  }) %>% ungroup() %>% distinct() %>% filter(n_common_snps > 0)

lfc_divided_1vs0 <- merge(genewise_mbased_1vs0, joined_df_filt, by = "SYMBOL", suffixes = c(".RNA", ".CS"), all.y = TRUE) %>%
  select(-c(ENTREZID.RNA, ENTREZID.CS, pValueHeterogeneity)) %>% 
  group_by(SYMBOL) %>% filter(any(diffASE)) %>% ungroup() %>% 
  group_by(across(-patient.CS)) %>% 
  summarize(patient.CS = paste(sort(patient.CS), collapse = ""), .groups = "drop")%>% 
  select(SYMBOL, patient.RNA, genewiseMAFDiff, padj) %>% distinct() %>% 
  group_by(SYMBOL) %>% filter(n() >= 2) %>% 
  group_modify(~ {
    patients <- .x$patient.RNA
    MAFDiffs <- .x$genewiseMAFDiff
    pValues <- .x$padj
    combn(length(patients), 2, simplify = FALSE) %>%
      map_dfr(function(pair) {
        i <- pair[1]
        j <- pair[2]
        
        patient_1 <- patients[i]
        patient_2 <- patients[j]
        mafDiff_1 <- MAFDiffs[i]
        mafDiff_2 <- MAFDiffs[j]
        
        # считаем LFCratio только если имеется одинаковый SNP в промоторе данного гена у обоих людей
        common_snps <- df_collapsed$rsid[df_collapsed$SYMBOL == .y$SYMBOL & (df_collapsed$patients == paste0(sort(unique(patients)), collapse = "") | df_collapsed$patients == "ABD")]
        mafDiffRatio <- if(length(common_snps) > 0) max(abs(c(mafDiff_1, mafDiff_2))) / min(abs(c(mafDiff_1, mafDiff_2))) else NA
        ddMAF <- if((length(common_snps) > 0) & (min(pValues) < 0.01)) abs(abs(mafDiff_1) - abs(mafDiff_2)) else NA
        
        tibble(
          patients = ifelse(patient_1 < patient_2, paste0(patient_1, patient_2), paste0(patient_2, patient_1)),
          mafDiffRatio = mafDiffRatio,
          ddMAF = ddMAF,
          logMafDiffRatio = if(!is.na(mafDiffRatio)) log2(mafDiffRatio) else NA,
          deltaMAF1 = mafDiff_1,
          deltaMAF2 = mafDiff_2,
          padj_1 = pValues[i],  padj_2 = pValues[j],
          rsids = list(common_snps),
          n_common_snps = length(common_snps)
        )
      })
  }) %>% ungroup() %>% distinct() %>% filter(n_common_snps > 0)

# выбираем 100 генов с наибольшим log(LFCratio)
ddmaf_df <- bind_rows(lfc_divided_0vs1, lfc_divided_1vs0) %>% arrange(desc(ddMAF)) %>% distinct()
topgenes <- (ddmaf_df %>% distinct(SYMBOL) %>% pull(SYMBOL))[1:100]
ddmaf_df_filt <- ddmaf_df[1:222, ] ## %>% filter(SYMBOL %in% topgenes) ## не так
length(unique(ddmaf_df_filt$SYMBOL))  
snps_for_mtfbrkr <- unique(unlist(ddmaf_df_filt$rsids))
length(snps_for_mtfbrkr) # 208
writeLines(snps_for_mtfbrkr, "/media/leon/DISK2/icig/done/snps_for_mtfbrkr_16.09.25.txt")

########### motifbreakR
setwd("/media/leon/DISK2/icig/done/motifbreakr_12_09_25/")
result_files <- list.files(pattern = "motifbreak_.*\\.rds$")

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

diffase_snp_tf <- ddmaf_df_filt %>%
  select(SYMBOL, patients, rsids) %>%
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(., snpgenedf, by = "SNP_id")

########## дифэкспрессия генов ТФ
counts <- read.table("/media/leon/DISK2/icig/done/in_vitro_results/counts/counts.tsv", header = TRUE, row.names = 1, sep = "\t")
counts <- counts[, c(6, 7, 8, 9, 11, 13)]
colnames(counts) <- gsub("\\.bam", "", colnames(counts))
colnames(counts) <- gsub("RNASEQ_", "", colnames(counts))
coldata <- data.frame(patient = substr(colnames(counts), 1, 3),
                      time = substr(colnames(counts), 5, 5))
rownames(coldata) <- colnames(counts)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ time)
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
length(unique(snpgenedf$SNP_id)) # 198
length(TFs) # 569
length(intersect(TFs, all_symbols)) # 317

# mapping hocomoco TF names to gene symbols
json_lines <- readLines("/media/leon/DISK2/icig/done/H13CORE_annotation.jsonl")

# Parse each JSON line and extract the TF:gene_symbol mapping
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

# Create a named vector for the dictionary (TF -> gene_symbol)
tf_to_gene <- setNames(tf_gene_dict$gene_symbol, tf_gene_dict$tf)
tf_to_gene["ZN821"] <- "ZNF821"
tf_to_gene["ZN704"] <- "ZNF704"

snpgenedf <- snpgenedf %>%
  mutate(TF = recode(TF, !!!tf_to_gene))

TFs <- unique(snpgenedf$TF)
length(TFs) # 564 (?)
length(intersect(TFs, all_symbols)) # 562

setdiff(TFs, all_symbols)

TF_counts <- counts %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  filter(SYMBOL %in% TFs)

rownames(TF_counts) <- TF_counts$ENSEMBL
TF_counts <- TF_counts[, c(2:7)]
dds <- DESeqDataSetFromMatrix(countData = TF_counts,
                              colData = coldata,
                              design = ~ patient)
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized = TRUE) >= 5 ) >= 2
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- results(dds)

###### выбираем попарно те дэги которые связаны с общими снп в промоторах генов по которым у этих двух людей разница в диффасе
resAB <- as.data.frame(results(dds, contrast = c("patient", "sAn", "sBn"))) %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "AB") %>%
  filter(SYMBOL %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "AB"])

resAD <- as.data.frame(results(dds, contrast = c("patient", "sAn", "sDm"))) %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "AD") %>%
  filter(SYMBOL %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "AD"])

resBD <- as.data.frame(results(dds, contrast = c("patient", "sBn", "sDm"))) %>%
  mutate(ENSEMBL = rownames(.)) %>%
  merge(., annotations_ordered, by = "ENSEMBL") %>%
  mutate(contrast = "BD") %>%
  filter(SYMBOL %in% diffase_snp_tf$TF[diffase_snp_tf$patients == "BD"])

TF_deseq_results <- bind_rows(resAB, resAD, resBD) %>%
  filter(padj < 0.1 & abs(log2FoldChange) > log2(1.5))
TF_deseq_results
all_degs <- unique(TF_deseq_results$SYMBOL)

deg_tf_rsids <- ddmaf_df_filt %>% 
  unnest(rsids) %>%
  rename(SNP_id = rsids) %>%
  merge(., snpgenedf, by = "SNP_id") %>%
  filter(TF %in% all_degs)

TF_genes_df <- deg_tf_rsids %>% group_by(TF) %>% summarize(genes = list(SYMBOL))
length(TF_genes_df$TF) # 45

go_list <- list()
for (i in seq(nrow(TF_genes_df))) {
  gen <-  TF_genes_df[i, ][[1]]
  print(gen)
  try({
    tfdf <- as.data.frame(enrichGO(gene = unique(unlist(TF_genes_df[i, ]$genes)),
                                   OrgDb = org.Hs.eg.db,
                                   keyType = "SYMBOL",
                                   ont = "BP",
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   # universe = res$SYMBOL,
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
    tfdf <- as.data.frame(enrichKEGG(gene = gens,
                                     organism = "hsa",
                                     pAdjustMethod = "BH",
                                     # universe = res$ENTREZID,
                                     qvalueCutoff = 0.1))
    if (nrow(tfdf) > 0) {
      tfdf$TF <- gen
      kegg_list[[gen]] <- tfdf
    }
  })
}
kegg_df <- bind_rows(kegg_list)

causal_snps_df <- deg_tf_rsids %>% select(SNP_id) %>%
  merge(., joined_df[, c("chr", "pos", "ref", "alt", "rsid")], by.x = "SNP_id", by.y = "rsid") %>%
  distinct()

# write_tsv(causal_snps_df[, c('chr', 'pos')], '/media/leon/DISK2/icig/done/79_snps.tsv', col_names = F)

check_snps <- function(path) {
  ddf <- readr::read_tsv(path, show_col_types = FALSE)
  colnames(ddf) <- c("sample", "chr", "pos", "dp", "ref", "QS_ref", "alt", "QS_alt", "alt2", "QS_alt2")
  if (max(ddf$QS_alt2) < 0.1) { ddf <- ddf %>% select(-alt2, -QS_alt2) }
  ddf <- ddf %>%
    separate(col = sample, into = c("patient", "timepoint"), sep = "_") %>%
    filter(stringr::str_length(ref) == 1) %>%
    filter(stringr::str_length(alt) == 1) %>%
    mutate(dp_ref = round(dp * QS_ref), dp_alt = round(dp * QS_alt)) %>%
    select(-QS_ref, -QS_alt, -dp)
  return(ddf)
}

chipseq_tables <- lapply(Sys.glob("/media/leon/DISK2/icig/chipseq_july/stats_for_79_snps/*.stat"), check_snps)
chipseq_table <- bind_rows(chipseq_tables)
chipseq_table <- merge(chipseq_table,
                       add_rsids(chipseq_table$chr, chipseq_table$pos, chipseq_table$ref),
                       by.x = c("chr", "pos"), by.y = c("seqnames", "pos"))
chipseq_table <- chipseq_table %>%
  filter(rsid %in% causal_snps_df$SNP_id) %>%
  group_by(patient, chr, pos, ref, alt, rsid) %>%
  summarize(total_dp_ref = sum(dp_ref), total_dp_alt = sum(dp_alt), .groups = 'drop')

# chipseq_table <- merge(chipseq_table,
#                        add_rsids(chipseq_table$chr, chipseq_table$pos, chipseq_table$ref),
#                        by.x = c("chr", "pos"), by.y = c("seqnames", "pos"))

chipseq_table_tfs <- deg_tf_rsids %>%
  select(SNP_id, SYMBOL, TF) %>% distinct() %>%
  group_by(SNP_id, SYMBOL) %>% summarize(TFs = list(TF), .groups = 'drop') %>%
  merge(chipseq_table, ., by.x = 'rsid', by.y = 'SNP_id')

patients_vv <- c("s3", "s4", "s5", "s6", "s8", "s9", "s12", "s15")
ase_in_vivo <- lapply(patients_vv, function(x) {
  ASE_in_vitro(x, min_DP = 100)
})

ase_vv <- bind_rows(ase_in_vivo) %>%
  mutate(padj = p.adjust(p01, method = "BH")) %>%
  mutate(diffASE = padj < 0.01)

ase_vv <- ase_vv %>%
  mutate(DP_ref_0R = round(DP_0R * QS_ref_0R),
         DP_alt_0R = round(DP_0R * QS_alt1_0R),
         DP_ref_1R = round(DP_1R * QS_ref_1R),
         DP_alt_1R = round(DP_1R * QS_alt1_1R)) %>%
  filter(DP_ref_0R > 10 & DP_alt_0R > 10 & DP_ref_1R > 10 & DP_alt_1R > 10)

ase_vv_subset <- ase_vv %>%
  filter(SYMBOL %in% chipseq_table_tfs$SYMBOL)

ase_vv$variant <- paste(ase_vv$chr, ase_vv$pos, sep = '_')
results_s3 <- run_mbased_analysis(ase_vv_subset, "s3")
results_s4 <- run_mbased_analysis(ase_vv_subset, "s4")
results_s5 <- run_mbased_analysis(ase_vv_subset, "s5")
results_s6 <- run_mbased_analysis(ase_vv_subset, "s6")
results_s8 <- run_mbased_analysis(ase_vv_subset, "s8")
results_s9 <- run_mbased_analysis(ase_vv_subset, "s9")
results_s12 <- run_mbased_analysis(ase_vv_subset, "s12")
results_s15 <- run_mbased_analysis(ase_vv_subset, "s15")

full_mbased_vv_0vs1 <- bind_rows(
  results_s3$full_0vs1,
  results_s4$full_0vs1,
  results_s5$full_0vs1,
  results_s6$full_0vs1,
  results_s8$full_0vs1,
  results_s9$full_0vs1,
  results_s12$full_0vs1,
  results_s15$full_0vs1
)

genewise_mbased_vv_0vs1 <- bind_rows(
  results_s3$genewise_0vs1,
  results_s4$genewise_0vs1,
  results_s5$genewise_0vs1,
  results_s6$genewise_0vs1,
  results_s8$genewise_0vs1,
  results_s9$genewise_0vs1,
  results_s12$genewise_0vs1,
  results_s15$genewise_0vs1
)

full_mbased_vv_1vs0 <- bind_rows(
  results_s3$full_1vs0,
  results_s4$full_1vs0,
  results_s5$full_1vs0,
  results_s6$full_1vs0,
  results_s8$full_1vs0,
  results_s9$full_1vs0,
  results_s12$full_1vs0,
  results_s15$full_1vs0
)

genewise_mbased_vv_1vs0 <- bind_rows(
  results_s3$genewise_1vs0,
  results_s4$genewise_1vs0,
  results_s5$genewise_1vs0,
  results_s6$genewise_1vs0,
  results_s8$genewise_1vs0,
  results_s9$genewise_1vs0,
  results_s12$genewise_1vs0,
  results_s15$genewise_1vs0
)

genewise_mbased_vv_0vs1$comp <- '0vs1'
genewise_mbased_vv_1vs0$comp <- '1vs0'


merge(chipseq_table_tfs[, c('rsid', 'patient', 'SYMBOL', 'TFs')],
      ase_vv[, c('SYMBOL', 'patient', 'diffASE')],
      by = c('patient', 'SYMBOL')) %>%
  group_by(SYMBOL) %>% filter(any(diffASE == T)) %>% distinct() %>% View()

### write_tsv(TF_genes_df['TF'], '~/45_tfs.tsv', col_names = F)
  
merge(chipseq_table_tfs[, c('rsid', 'patient', 'SYMBOL', 'TFs')],
      bind_rows(genewise_mbased_vv_0vs1, genewise_mbased_vv_1vs0)[, c('SYMBOL', 'patient', 'comp', 'padj', 'genewiseMAFDiff', 'diffASE')],
      by = c('patient', 'SYMBOL')) %>%
  distinct() %>% group_by(SYMBOL) %>% filter(any(diffASE == T)) %>% ungroup() %>%
  # group_by(SYMBOL, patient) %>% filter(diffASE == T) %>% 
  View()
      



dmm <- read_tsv("/media/leon/DISK2/icig/done/for_missense/sDm.stat", show_col_types = FALSE)
anm <- read_tsv("/media/leon/DISK2/icig/done/for_missense/sAn.stat", show_col_types = FALSE)
bnm <- read_tsv("/media/leon/DISK2/icig/done/for_missense/sBn.stat", show_col_types = FALSE)

all(all(colnames(dmm) == colnames(anm)), all(colnames(dmm) == colnames(bnm)))
joined_dfm <- rbind(anm, dmm, bnm)
joined_dfm$chr <- as.character(as.numeric(joined_dfm$CHR))
joined_dfm <- na.omit(joined_dfm)
joined_dfm <- joined_dfm[stringr::str_length(joined_dfm$REF) == 1, ]
joined_dfm <- joined_dfm[stringr::str_length(joined_dfm$ALT1) == 1, ]

# пересекаем с экзонами
hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene

snp_ranges <- GRanges(
  seqnames = paste0("chr", joined_dfm$CHR),
  ranges = IRanges(start = joined_dfm$POS, width = 1)
)
overlaps <- findOverlaps(snp_ranges, exns)
joined_dfm <- joined_dfm[queryHits(overlaps), ]

snp_exons <- exns[subjectHits(overlaps)]
gene_ids <- sapply(snp_exons$gene_id, function(x) ifelse(length(x) == 0, NA, unlist(x)))
joined_dfm$ENTREZID <- gene_ids
joined_dfm$SYMBOL <- mapIds(org.Hs.eg.db,
                       keys = joined_dfm$ENTREZID,
                       keytype = "ENTREZID",
                       column = "SYMBOL")

joined_dfm <- na.omit(joined_dfm)
joined_dfm <- merge(joined_dfm, add_rsids(joined_dfm$CHR, joined_dfm$POS, joined_dfm$REF), by.x = c("CHR", "POS"), by.y = c("seqnames", "pos"))
# writeLines(unique(joined_dfm$rsid), "~/potentially_missense_all_tfs.tsv")

tf_patient_df <- ddmaf_df %>% dplyr::select(rsids, patients) %>% unnest(rsids) %>%
  merge(., snpgenedf, by.x = 'rsids', by.y = 'SNP_id')

all_missenses <- read_csv("~/rsid_gene_missense_all_tfs.csv") %>% 
  filter(is_missense == TRUE) %>%
  tidyr::separate_rows(gene, sep = ";") %>%
  merge(., joined_dfm[, c('SAMPLE', 'rsid')], by = "rsid") %>%
  distinct() %>%
  rename(patient = SAMPLE) %>%
  filter((patient == 'sDm' & gene %in% tf_patient_df$TF[tf_patient_df$patients == 'AD' | tf_patient_df$patients == 'BD']) | (patient == 'sAn' & gene %in% tf_patient_df$TF[tf_patient_df$patients == 'AD' | tf_patient_df$patients == 'AB']) | (patient == 'sBn' & gene %in% tf_patient_df$TF[tf_patient_df$patients == 'BD' | tf_patient_df$patients == 'AB']))

length(unique(all_missenses$gene))

write_tsv(all_missenses, '~/all_missenses_filtered_to_check.tsv')
