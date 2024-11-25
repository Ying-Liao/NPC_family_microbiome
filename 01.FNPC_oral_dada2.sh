# 1. library packages
library(dada2);paste("version of dada2 ", packageVersion("dada2"))
library(ggplot2);paste("version of ggplot2 ", packageVersion("ggplot2"))
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(reshape2);paste("version of reshape2 ", packageVersion("reshape2"))
library(RColorBrewer);paste("version of RColorBrewer ", packageVersion("RColorBrewer"))
library(optparse);paste("version of optparse ", packageVersion("optparse"))
library(gridExtra); packageVersion("gridExtra")

# 2. read in the names of the fastq files
# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(paste(path, "/rawdata", sep=""), pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(paste(path, "/rawdata", sep=""), pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_R#.fastq
sample.names <- sub('.........$','',basename(fnFs))
# write out sample.names
write.table(sample.names, paste(path.file,"/sample.names.txt", sep=""), sep="\t")

# 3. Inspect read quality profiles
#   3.1 visualizing the quality profiles of the forward reads
QC_F <- plotQualityProfile(fnFs[7:12])
ggsave(QC_F, file=paste(path.fig,'QC_Forward.pdf',sep=""))
#   3.2 visualizing the quality profiles of the reverse reads
QC_R <- plotQualityProfile(fnRs[7:12])
ggsave(QC_R, file=paste(path.fig,'QC_Reverse.pdf',sep=""))

# 4. Filter and trim
#   4.1 Assign the filenames for the filtered fastq.gz files
#       Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered_data", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered_data", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#   4.2 Filter and trim (using standard filtering parameters)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,200),trimLeft=c(31,32),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)

# 5. Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plot_er <- plotErrors(errF, nominalQ=TRUE)
ggsave(plot_er, file=paste(path.fig,'plot_error_rates.pdf',sep=""))


# 6  Sample Inference
filt0 <- rownames(out[which(out[,2]==0) ,])
filt0 <- sub('.........$','',filt0)

filtFs_f <-filtFs[setdiff(names(filtFs), filt0)]
filtRs_f <-filtRs[setdiff(names(filtRs), filt0)]
dadaFs <- dada(filtFs_f, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs_f, err=errR, multithread=TRUE)

saveRDS(dadaFs, paste(path.rds, "/dadaFs.rds", sep=""))
saveRDS(dadaRs, paste(path.rds, "/dadaRs.rds", sep=""))
dadaFs[[1]]
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# 7. Merge paired reads to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, filtFs_f, dadaRs, filtRs_f, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
saveRDS(mergers, paste(path.rds, "/mergers.rds", sep=""))

# 8. Construct sequence table -- amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#   Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
write.table(seqtab, file=paste(path.file,'/ASVtable_raw.txt',sep=""), sep="\t")

# 9. Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

write.table(seqtab.nochim, file=paste(path.file,'/ASVtable.txt',sep=""))
saveRDS(seqtab.nochim, paste(path.rds, "/seqtab.nochim.rds", sep=""))

# 10. Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, file=paste(path.file,'/dada2_track_stats.txt',sep=""))