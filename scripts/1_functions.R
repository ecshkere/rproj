options(readr.show_col_types = FALSE)

## genome annotation
hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
prmtrs <- promoters(hg38, columns = "gene_id")
exns <- exons(hg38, columns = c("exon_id", "gene_id"))
gens <- genes(hg38)

# returns the input dataframe with "gene_id" and "symbol" columns added
assign_genes <- function(df, chr_col, start_col, end_col = start_col,
                         promoters = FALSE, exons_only = FALSE, omit_na = TRUE) {
  seqnames_col <- gsub("^([0-9XY])", "chr\\1", df[[chr_col]])
  ranges <- GRanges(seqnames = seqnames_col,
                    ranges = IRanges(start = df[[start_col]],
                                     end = df[[end_col]]))
  df <- df %>% mutate(entrezid = NA_character_)
  
  if (promoters) {
    overlaps <- findOverlaps(ranges, prmtrs)
    df$entrezid[queryHits(overlaps)] <- prmtrs$gene_id[subjectHits(overlaps)]
    df <- df %>% unnest_longer(entrezid)
  } else if (exons_only) {
    overlaps <- findOverlaps(ranges, exns)
    df$entrezid[queryHits(overlaps)] <- exns$gene_id[subjectHits(overlaps)]
    df <- df %>% unnest_longer(entrezid) %>% 
      filter(!(is.na(entrezid) | sapply(entrezid, length) == 0))
  } else {
    overlaps <- findOverlaps(ranges, gens)
    df$entrezid[queryHits(overlaps)] <- gens$gene_id[subjectHits(overlaps)]
  }
  df <- df %>%
    mutate(symbol = mapIds(org.Hs.eg.db,
                           keys = entrezid,
                           column = "SYMBOL",
                           keytype = "ENTREZID")) %>%
    unnest_longer(symbol)
  
  if (omit_na) { df <- drop_na(df, symbol, entrezid) }
  if ("symbol_id" %in% colnames(df)) { df <- df %>% select(-symbol_id) } # ?
  return(df)
}

# returns the input dataframe with "rsid" column added  
# omits rows without rsid 
# todo gsub "chr" with nothing
add_rsids <- function(df, chr_col, pos_col, ref_col, omit_na = TRUE) {
  ranges <- GRanges(seqnames = as.character(df[[chr_col]]),
                    ranges = IRanges(df[[pos_col]], width = 1),
                    ref = df[[ref_col]])
  snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, ranges)
  
  snp_df <- as.data.frame(snps)[, c("seqnames", "pos")] %>% 
    mutate(rsid = mcols(snps)$RefSNP_id) %>% drop_na()
  
  df <- left_join(df, snp_df,
                  by = set_names(c("seqnames", "pos"), c(chr_col, pos_col)),
                  relationship = "many-to-many")
  
  if (omit_na) { df <- drop_na(df, rsid) }
  return(df)
}


# takes a vector of rsids as input
# splits it in batches
# creates an RDS object for each SNP in the $mtfbrkr_dir directory
run_motifbreakr <- function(rsids_to_process, mtfbrkr_dir, batch_size = 25) {
  batches <- split(rsids_to_process, ceiling(seq_along(rsids_to_process) / batch_size))
  cat("Starting motifbreakR analysis for", length(rsids_to_process), "rsIDs in", length(batches), "batches\n")
  
  counter <- 0
  for (batch in batches) {
    snp_batch <- snps.from.rsid(rsid = batch,
                                dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                search.genome = BSgenome.Hsapiens.UCSC.hg38)
    
    res_batch <- motifbreakR(snpList = snp_batch,
                             pwmList = MotifDb,
                             filterp = TRUE,
                             threshold = 0.0001,
                             method = "ic")
    
    for (rsid in batch) {
      output_file <- file.path(mtfbrkr_dir, paste0(rsid, ".rds"))
      snp_res <- res_batch[res_batch$SNP_id == rsid]
      if (length(snp_res) > 0) {
        saveRDS(snp_res, output_file)
      }
    }
    counter <- counter + 1
    cat("Batch", counter, "out of", length(batches), "done\n")
  }
}

# mapping hocomoco TF names to gene symbols
json_lines <- readLines("data/H13CORE_annotation.jsonl")

tf_gene_dict <- map_dfr(json_lines, function(line) {
  data <- jsonlite::fromJSON(line)
  if (!is.null(data$masterlist_info$species$HUMAN$gene_symbol)) {
    data.frame(tf = data$tf,
               motif_name = data$name,
               gene_symbol = data$masterlist_info$species$HUMAN$gene_symbol,
               stringsAsFactors = FALSE)
  } else { NULL }
})
tf_to_gene <- setNames(tf_gene_dict$gene_symbol, tf_gene_dict$tf)
tf_to_gene["ZN821"] <- "ZNF821"
tf_to_gene["ZN704"] <- "ZNF704"

## VECTOR OF ALL GENE SYMBOLS
all_genes <- read_tsv("counts/in_vivo_counts.tsv", skip = 1, col_names = F, col_select = 1)[[1]]
annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = all_genes,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "ENSEMBL") %>% 
  rename(entrezid = ENTREZID, symbol = SYMBOL, ensembl = ENSEMBL)

annotations_ordered <- annotations[match(all_genes, annotations$ensembl), ]
annotations_clean <- annotations_ordered %>%
  filter(!is.na(entrezid) & !is.na(symbol))
all_symbols <- unique(annotations_clean$symbol)
all_entrezids <- unique(annotations_clean$entrezid)

all_counts <- read.table("counts/in_vivo_counts.tsv", header = TRUE, sep = "\t")
all_coldata <- read.table("counts/in_vivo_coldata.tsv", header = TRUE, sep = "\t", stringsAsFactors = TRUE)

dds <- DESeqDataSetFromMatrix(countData = all_counts,
                              colData = all_coldata,
                              design = ~ patient)
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized = TRUE) >= 5 ) >= 15
background_ens <- names(keep[keep])
background_symb <- annotations_clean$symbol[annotations_clean$ensembl %in% background_ens]
background_entr <- annotations_clean$entrezid[annotations_clean$ensembl %in% background_ens]

symbol_to_entrez <- function(x) {
  if (is.na(as.numeric(x[1]))) {
    mapIds(org.Hs.eg.db,
           keys = x,
           keytype = "SYMBOL",
           column = "ENTREZID")
  }
}

go <- function(x, q_v = 0.1) {
  enrichGO(gene = x,
           OrgDb = org.Hs.eg.db,
           keyType = "SYMBOL",
           universe = background_symb,
           ont = "BP",
           pAdjustMethod = "BH",
           qvalueCutoff = q_v,
           readable = TRUE) %>%
    as.data.frame()
}

kegg <- function(x, q_v = 0.1) {
  if (is.na(as.numeric(x[1]))) { x <- symbol_to_entrez(x) }
  k <- enrichKEGG(gene = x,
                  organism = "hsa",
                  pAdjustMethod = "BH",
                  universe = background_entr,
                  qvalueCutoff = q_v) %>% as.data.frame()
  
  if (nrow(k) > 0) {
    entrezes <- strsplit(k$geneID, "/")
    
    entrez_to_symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                              keys = unique(unlist(entrezes)),
                                              columns = "SYMBOL",
                                              keytype = "ENTREZID")
    
    symbols <- lapply(entrezes, function(ids) {
      res <- entrez_to_symbol$SYMBOL[match(ids, entrez_to_symbol$ENTREZID)]
      paste(res, collapse = "/")})
    
    k$geneSYMBOL <- unlist(symbols)
  }
  return(k)
}

## ALL MISSENSE SNPS IN TRANSCRIPTION FACTORS
if (FALSE) {
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  attrs <- c("external_gene_name", "chromosome_name", "exon_chrom_start", "exon_chrom_end")
  
  exons_df <- getBM(attributes = attrs,
                    filters = "external_gene_name",
                    values = tf_to_gene,
                    mart = ensembl) %>%
    mutate(chromosome_name = paste0("chr", chromosome_name),
           start = exon_chrom_start,
           end = exon_chrom_end) %>%
    select(chromosome_name, exon_chrom_start, exon_chrom_end)
  
  write_tsv(exons_df, "output/exons_for_all_tfs.bed", col_names = FALSE)
  
  ## bedtools intersect -a /media/leon/Polina/atac_rna/dbSnp153.bed -b output/exons_for_all_tfs.bed | awk '{print $1, $2, $3, $4}' - > output/all_tfs_exon_rsids.bed
  
  snps_gr <- read_table("output/all_tfs_exon_rsids.bed", col_names = F) %>% 
    filter(X3 - X2 == 1)
  the_rsids <- unique(snps_gr$X4)

  snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
  snps_cons <- data.frame()
  for (i in seq(1, length(the_rsids), by = 10000)) {
    snps_cons <- rbind(snps_cons,
                       getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "consequence_type_tv"),
                             filters = "snp_filter",
                             values = the_rsids[i:(i+10000)],
                             mart = snp_mart))
    print(i)
  }
  # saveRDS(snps_cons, "output/biomart.RDS")
  # snps_cons <- readRDS("output/biomart.RDS")
  
  snps_missense <- snps_cons %>% filter(consequence_type_tv == "missense_variant")
  
  snps_missense <- assign_genes(snps_missense, "chr_name", "chrom_start", exons_only = TRUE) %>% 
    filter(symbol %in% tf_to_gene, !chr_name %in% c("X", "Y")) %>%
    select(-entrezid) %>% distinct() %>% rename(tf = symbol)
  
  snps_missense %>% select(chr_name, chrom_start) %>% distinct() %>% write_tsv("output/2025.10.27.all_missense_positions.tsv", col_names = F) # for missenses.pl
}
