library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ape); packageVersion("ape")

phyloseq <- ps.filter4

# FNPC beta diversity------------
# unifrac, wunifrac, euclidean,bray,jaccard
dist_list <- c("unifrac", "wunifrac" ,"euclidean", "bray", "jaccard")
dist_name <- c("unifrac_dist", "wunifrac_dist" ,"euclidean_dist", "bray_dist", "jaccard_dist")

# phy_tree-----
random_tree <- rtree(ntaxa(phyloseq), rooted=TRUE, tip.label = taxa_names(phyloseq))
tree_plot <- plot(random_tree)

phyloseq <- merge_phyloseq(phyloseq, random_tree)


# calculate distance matrix--------
unifrac_dist <- phyloseq::distance(phyloseq, method="unifrac")
wunifrac_dist <- phyloseq::distance(phyloseq, method="wunifrac")
euclidean_dist <- phyloseq::distance(phyloseq, method="euclidean")
bray_dist <- phyloseq::distance(phyloseq, method="bray")
jaccard_dist <- phyloseq::distance(phyloseq, method="jaccard")

saveRDS(unifrac_dist, paste(path.rds, "unifrac_dist.rds", sep=""))
saveRDS(wunifrac_dist, paste(path.rds, "wunifrac_dist.rds", sep=""))
saveRDS(euclidean_dist, paste(path.rds, "euclidean_dist", sep=""))
saveRDS(bray_dist, paste(path.rds, "bray_dist.rds", sep=""))
saveRDS(jaccard_dist, paste(path.rds, "jaccard_dist.rds", sep=""))

m.unifrac_dist <- as.matrix(unifrac_dist)
write.table(m.unifrac_dist,paste(path.file, "unifrac_dist.txt",sep = ""),sep = "\t")
m.wunifrac_dist <- as.matrix(wunifrac_dist)
write.table(m.wunifrac_dist,paste(path.file, "wunifrac_dist.txt",sep = ""),sep = "\t")
m.euclidean_dist <- as.matrix(euclidean_dist)
write.table(m.euclidean_dist,paste(path.file, "euclidean_dist.txt",sep = ""),sep = "\t")
m.bray_dist <- as.matrix(bray_dist)
write.table(m.bray_dist,paste(path.file, "bray_dist.txt",sep = ""),sep = "\t")
m.jaccard_dist <- as.matrix(jaccard_dist)
write.table(m.jaccard_dist,paste(path.file, "jaccard_dist.txt",sep = ""),sep = "\t")

# calculate alpha diversoty indice--------
alpha_div <- estimate_richness(phyloseq, split = TRUE, measures = NULL)
write.table(alpha_div, paste(path.file, "alpha_div.txt", sep = ""), sep = "\t")