library("DESeq2")

# import data
meta <- read.csv("./case_control.csv")
meta$group <-as.factor(meta$group)
meta$sex <-as.factor(meta$sex)
meta$smoke <-as.factor(meta$smoke)
asv <- read.table("./coreASVs_table.txt", row.names = 1, header = TRUE, sep = "\t")

# Perform DESeq analysis
dds <- DESeqDataSetFromMatrix(countData = asv, colData = meta, 
                              design = ~age+sex+smoke+group)
dds$group <- relevel( dds$group, "0" )
dds <- DESeq(dds)
suppressMessages(dds)
cooks <- assays(dds)[["cooks"]]
maxCooks <- as.data.frame(apply(cooks,1,max))
rep <- assays(dds)[["replaceCounts"]]
res <- results(dds, contrast=c('group', '1', '0'))
deseq_res <- as.data.frame(res[order(res$padj), ])
deseq_res$otu_id <- rownames(deseq_res)
deseq_res$baseMean <- deseq_res$baseMean / sum(deseq_res$baseMean)
write.table(deseq_res, './result_DESeq2.txt', row.names = FALSE, sep = '\t', quote = FALSE)