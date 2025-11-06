#!/usr/bin/env Rscript

# motifbreakr.R
# Usage: Rscript motifbreakr.R <input_file.txt> <output_dir> <batch_size>

suppressPackageStartupMessages({
  library(motifbreakR)
  library(MotifDb)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
})

args <- commandArgs(trailingOnly = TRUE)

# Check input arguments
if (length(args) != 3) {
  stop("Usage: Rscript motifbreakr.R <input_file.txt> <output_dir> <batch_size>", call. = FALSE)
}

input_file <- args[1]
output_dir <- args[2]
batch_size <- as.integer(args[3])
dir.create(output_dir, showWarnings = FALSE)

# Read rsIDs from input file
rsids <- unique(readLines(input_file))
rsids <- rsids[rsids != ""]  # Remove empty lines
rsids <- grep("^rs\\d+$", rsids, value = TRUE)  # Keep valid rsIDs

if (length(rsids) == 0) {
  stop("No valid rsIDs found in input file", call. = FALSE)
}

# Batch processing parameters
total <- length(rsids)
batches <- split(rsids, ceiling(seq_along(rsids)/batch_size))
num_batches <- length(batches)

cat("Starting motifbreakR analysis for", total, "rsIDs in", num_batches, "batches\n")
cat("============================================\n\n")

# Initialize counters and error log
processed <- 0
batch_counter <- 0
errors <- character(0)

for (batch in batches) {
  batch_counter <- batch_counter + 1
  batch_rsids <- unlist(batch)
  current_batch_size <- length(batch_rsids)
  processed <- processed + current_batch_size
  
  cat("\n[ BATCH", batch_counter, "/", num_batches, "]\n")
  cat("Processing", current_batch_size, "rsIDs:", paste(batch_rsids, collapse = ", "), "\n")
  
  tryCatch({
    # Retrieve SNP data in batch
    snp_batch <- snps.from.rsid(
      rsid = batch_rsids,
      dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
      search.genome = BSgenome.Hsapiens.UCSC.hg38
    )
    
    # Run motifbreakR analysis on batch
    res_batch <- motifbreakR(
      snpList = snp_batch,
      filterp = TRUE,
      pwmList = MotifDb,
      threshold = 0.0001,
      method = "ic",
      bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
      verbose = FALSE
    )
    
    # Save individual results
    for (rsid in batch_rsids) {
      output_file <- file.path(output_dir, paste0(rsid, ".rds"))
      res_single <- res_batch[res_batch$SNP_id == rsid]
      
      if (length(res_single) > 0) {
        saveRDS(res_single, file = output_file)
        cat("  ✓ Saved:", rsid, "\n")
      } else {
        msg <- paste("No results for", rsid)
        errors <<- c(errors, msg)
        cat("  !! WARNING:", msg, "\n")
      }
    }
    
    # Explicit memory cleanup
    rm(snp_batch, res_batch)
    gc(full = TRUE)
    
  }, error = function(e) {
    msg <- paste("Batch", batch_counter, "failed:", conditionMessage(e))
    errors <<- c(errors, msg)
    cat("  !! BATCH ERROR:", msg, "\n")
    
    # Try to save partial results if available
    if (exists("res_batch")) {
      for (rsid in batch_rsids) {
        output_file <- file.path(output_dir, paste0(rsid, ".rds"))
        if (rsid %in% res_batch$SNP_id) {
          res_single <- res_batch[res_batch$SNP_id == rsid]
          saveRDS(res_single, file = output_file)
          cat("  ✓ Partial save:", rsid, "\n")
        }
      }
    }
  })
  

  cat("Completed:", processed, "/", total, "rsIDs\n")
}

# Final summary
cat("\n============================================\n")
cat("Processing complete\n")
cat("Total rsIDs processed:", processed, "/", total, "\n")
cat("Successful batches:", num_batches - length(errors), "/", num_batches, "\n")
cat("Failures:", length(errors), "\n")

if (length(errors) > 0) {
  cat("\nError details:\n")
  writeLines(errors)
  writeLines(errors, "motifbreak_errors.log")
}
