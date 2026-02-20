#!/usr/bin/env Rscript

library(dada2)
library(here)

all_track <- list()

# -------------------------------------------------------------------
# 1. Discover samples
# -------------------------------------------------------------------

# Find all forward reads in the length_filtered directory
filtFs <- sort(list.files(
  here("trim_clean_qc", "length_filtered"),
  pattern = "_trnL_R1.lenfilt.fastq.gz$",
  full.names = TRUE
))

# Derive reverse reads by replacing R1 → R2
filtRs <- gsub("_R1.lenfilt.fastq.gz$", "_R2.lenfilt.fastq.gz", filtFs)

# Extract sample names (everything before _trnL_R1...)
sample_names <- sub("_trnL_R1.lenfilt.fastq.gz$", "", basename(filtFs))

cat("Discovered samples:\n")
print(sample_names)

# Create top-level output directory
top_outdir <- "dada2"
if (!dir.exists(top_outdir)) {
  dir.create(top_outdir, recursive = TRUE)
}

# -------------------------------------------------------------------
# 2. Learn pooled error models (across all samples)
# -------------------------------------------------------------------

cat("\n=========================================\n")
cat("Learning pooled error models across all samples...\n")
cat("=========================================\n\n")

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Save pooled error model plots once at the top level
pdf(file.path(top_outdir, "trnL_pooled_error_model_forward.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf(file.path(top_outdir, "trnL_pooled_error_model_reverse.pdf"))
plotErrors(errR, nominalQ = TRUE)
dev.off()

# -------------------------------------------------------------------
# 3. Loop over samples
# -------------------------------------------------------------------

for (i in seq_along(sample_names)) {
  
  sample <- sample_names[i]
  filtF <- filtFs[i]
  filtR <- filtRs[i]
  
  cat("\n=========================================\n")
  cat("Processing sample:", sample, "\n")
  cat("Forward:", filtF, "\n")
  cat("Reverse:", filtR, "\n")
  cat("=========================================\n\n")
  
  # Create per-sample output directory
  outdir <- file.path(top_outdir, sample)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # -------------------------------------------------------------------
  # Dereplication
  # -------------------------------------------------------------------
  cat("Dereplicating...\n")
  derepF <- derepFastq(filtF, verbose = TRUE)
  derepR <- derepFastq(filtR, verbose = TRUE)
  
  cat("Forward derep unique sequences:", length(derepF$uniques), "\n")
  cat("Reverse derep unique sequences:", length(derepR$uniques), "\n")
  
  # -------------------------------------------------------------------
  # Denoising using pooled error models
  # -------------------------------------------------------------------
  cat("Running dada() with pooled error models...\n")
  dadaF <- dada(derepF, err = errF, multithread = TRUE)
  dadaR <- dada(derepR, err = errR, multithread = TRUE)
  
  # -------------------------------------------------------------------
  # Merge paired reads
  # -------------------------------------------------------------------
  cat("Merging pairs...\n")
  mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose = TRUE)
  
  cat("Total input pairs:", sum(derepF$uniques), "\n")
  cat("Merged pairs:", nrow(mergers), "\n")
  
  # -------------------------------------------------------------------
  # Build sequence table
  # -------------------------------------------------------------------
  seqtab <- makeSequenceTable(mergers)
  
  cat("Sequence table dimensions (samples x ASVs):\n")
  print(dim(seqtab))
  
  # -------------------------------------------------------------------
  # ASV length distribution
  # -------------------------------------------------------------------
  asv_lengths <- nchar(getSequences(seqtab))
  cat("ASV length distribution:\n")
  print(table(asv_lengths))
  
  pdf(file.path(outdir, "trnL_asv_length_histogram.pdf"))
  hist(asv_lengths, breaks = 30, main = "trnL ASV Length Distribution", xlab = "bp")
  dev.off()
  
  # -------------------------------------------------------------------
  # Chimera removal
  # -------------------------------------------------------------------
  seqtab_nochim <- removeBimeraDenovo(
    seqtab,
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE
  )
  
  cat("Reads retained after chimera removal:",
      sum(seqtab_nochim), "of", sum(seqtab), "\n")
  
  # -------------------------------------------------------------------
  # Read retention tracking
  # -------------------------------------------------------------------
  
  getN <- function(x) sum(getUniques(x))
  
  track <- data.frame(
    input      = sum(derepF$uniques),   # total input pairs
    denoisedF  = getN(dadaF),
    denoisedR  = getN(dadaR),
    merged     = getN(mergers),
    nonchim    = sum(seqtab_nochim)
  )
  
  rownames(track) <- sample
  
  all_track[[sample]] <- track
  
  # -------------------------------------------------------------------
  # Save outputs
  # -------------------------------------------------------------------
  saveRDS(dadaF, file.path(outdir, "trnL_dadaF.rds"))
  saveRDS(dadaR, file.path(outdir, "trnL_dadaR.rds"))
  saveRDS(mergers, file.path(outdir, "trnL_merged_pairs.rds"))
  saveRDS(seqtab_nochim, file.path(outdir, "trnL_seqtab_nochim.rds"))
  
  cat("Finished sample:", sample, "\n")
}

combined <- do.call(rbind, all_track)

write.table(
  combined,
  file = file.path(top_outdir, "trnL_read_retention_all_samples.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cat("\nAll samples processed.\n")

# After the loop, combine all seqtab_nochim objects
all_seqtabs <- list.files(
  path = "dada2",
  pattern = "trnL_seqtab_nochim.rds$",
  full.names = TRUE,
  recursive = TRUE
)

seqtabs <- lapply(all_seqtabs, readRDS)

# Extract sample names from directory names
sample_names <- basename(dirname(all_seqtabs))

# Assign rownames to each table
for (i in seq_along(seqtabs)) {
  rownames(seqtabs[[i]]) <- sample_names[i]
}

# Merge safely
if (length(seqtabs) == 1) {
  seqtab_all <- seqtabs[[1]]
} else {
  seqtab_all <- do.call(mergeSequenceTables, seqtabs)
}

# Extract ASV sequences
asv_seqs <- colnames(seqtab_all)

asv_ids <- paste0("ASV_", seq_along(asv_seqs))
colnames(seqtab_all) <- asv_ids

# Write FASTA
asv_fasta <- "dada2/trnL_ASVs.fasta"
sink(asv_fasta)
for (i in seq_along(asv_seqs)) {
  cat(paste0(">ASV_", i, "\n"))
  cat(asv_seqs[i], "\n")
}
sink()

# Export AVS table for Qiime2
write.table(
  seqtab_all,
  file = "dada2/trnL_seqtab_all.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


