#!/usr/bin/env Rscript

library(dada2)
library(here)

# 1. Set file paths

filtF <- here("trim_clean_qc", "length_filtered", "23412_S1_trnL_R1.lenfilt.fastq.gz") 
filtR <- here("trim_clean_qc", "length_filtered", "23412_S1_trnL_R2.lenfilt.fastq.gz")

sample_name <- "23412_S1"

# 2. Dereplication

cat("Dereplicating...\n")
derepF <- derepFastq(filtF, verbose = TRUE) 
derepR <- derepFastq(filtR, verbose = TRUE) 

# Inspect derep sizes 
cat("Forward derep unique sequences:", length(derepF$uniques), "\n") 
cat("Reverse derep unique sequences:", length(derepR$uniques), "\n")

# 3. Learn error models

cat("Learning error models...\n")
errF <- learnErrors(filtF, multithread = TRUE) 
errR <- learnErrors(filtR, multithread = TRUE)

# Plot error models 
pdf("dada2_error_model_forward.pdf") 
plotErrors(errF, nominalQ = TRUE) 
dev.off() 
pdf("dada2_error_model_reverse.pdf") 
plotErrors(errR, nominalQ = TRUE) 
dev.off()

# 4. Denoising

cat("Running dada()...\n")
dadaF <- dada(derepF, err = errF, multithread = TRUE) 
dadaR <- dada(derepR, err = errR, multithread = TRUE)

# 5. Merge paired reads

cat("Merging pairs...\n")
mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose = TRUE)

# Inspect merge rates
cat("Total input pairs:", sum(derepF$uniques), "\n")
cat("Merged pairs:", nrow(mergers), "\n")

# 6. Build sequence table

seqtab <- makeSequenceTable(mergers)

cat("Sequence table dimensions (samples x ASVs):\n")
print(dim(seqtab))

# 7. Inspect ASV length distribution

asv_lengths <- nchar(getSequences(seqtab))
cat("ASV length distribution:\n")
print(table(asv_lengths))

pdf("dada2_asv_length_histogram.pdf")
hist(asv_lengths, breaks = 30, main = "ASV Length Distribution", xlab = "bp")
dev.off()

# 8. Chimera removal

seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)

cat("Reads retained after chimera removal:",
    sum(seqtab_nochim), "of", sum(seqtab), "\n")


# 9. Save outputs

saveRDS(seqtab_nochim, file = "dada2_trnL_seqtab_nochim.rds")
saveRDS(errF, file = "dada2_trnL_errF.rds")
saveRDS(errR, file = "dada2_trnL_errR.rds")

cat("Done.\n")

