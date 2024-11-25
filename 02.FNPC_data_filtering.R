library(ggplot2);paste("version of ggplot2 ", packageVersion("ggplot2"))
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")

# Import data -----
seqtab <- readRDS("seqtab.nochim.rds")
taxa <- readRDS("taxa.rds")
samdf <- read.table("./files/samplelist.txt", sep = "\t", row.names = 1, header = TRUE)

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa))

ps

dna <- Biostrings::DNAStringSet(taxa_names(ps)) 
colnames(otu_table(seqtab, taxa_are_rows=FALSE)) == taxa_names(ps) 
row.names(taxa) == taxa_names(ps)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
names(dna) <- taxa_names(ps)
colnames(seqtab) <- taxa_names(ps)
row.names(taxa) <- taxa_names(ps)
head(seqtab)
head(dna)
head(taxa)

ps_names <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
                     sample_data(samdf),
                     tax_table(taxa)) 
ps_names <- merge_phyloseq(ps_names, dna)
saveRDS(ps_names, paste(path.rds, "/ps_names.rds", sep=""))

ps.filter1 <- prune_taxa(taxa_sums(ps_noqc) > 20, ps_names)
ps.filter1

prevelancedf = apply(X = otu_table(ps_noqc, taxa_are_rows=FALSE),
                     MARGIN = 2,
                     FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevelancedf = data.frame(Prevalence = prevelancedf,
                          preRate = prevelancedf/nsamples(ps_noqc), 
                          TotalAbundance = taxa_sums(otu_table(ps_noqc, taxa_are_rows=FALSE)),
                          tax_table(ps_noqc))

write.table(prevelancedf, file = paste(path.file,'/taxa_preRate.txt',sep=""), sep="\t")

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevelancedf)[(prevelancedf$preRate >= 0.01)]
length(keepTaxa)

ps.filter2 = prune_taxa(keepTaxa, ps.filter1)

saveRDS(ps.filter2, paste(path.rds, "/ps.filter2.rds", sep=""))

seqdepth = apply(X = otu_table(ps.filter2, taxa_are_rows=FALSE),
                 MARGIN = 1,
                 FUN = function(x){sum(x)})
summary(seqdepth)

seqdepth = data.frame(seqDepth = seqdepth,
                      sample_data(ps.filter2))

write.table(seqdepth, paste(path.file, "/sample_seq_depth.txt", sep = ""), sep = "\t")


keepsample = rownames(seqdepth)[(seqdepth$seqDepth >= 5000)]
ps.filter3 <- prune_samples(keepsample, ps.filter2)
ps.filter3

saveRDS(ps.filter3, paste(path.rds, "/ps.filter3.rds", sep=""))

# Normalize number of reads in each sample.
sample_sums(ps.filter3)
set.seed(13)
ps.filter4 = rarefy_even_depth(ps.filter3,sample.size = min(sample_sums(ps.filter3)),
                               rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

sample_sums(ps.filter4) 
saveRDS(ps.filter4, paste(path.rds, "/ps.filter4.rds", sep=""))


