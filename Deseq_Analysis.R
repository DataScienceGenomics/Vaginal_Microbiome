#DESe2 Analysis

library("DESeq2")
otu_table(carbomVtaxa)<-otu_table(carbomVtaxa)+1
diagdds = phyloseq_to_deseq2(carbomVtaxa, ~ Labor_Fever)##FebrileVSAfebrile
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(carbomVtaxa)[rownames(sigtab), ]
                                            , "matrix"))
head(sigtab)
View(sigtab)

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Species
x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
x = sort(x, TRUE)
sigtab$Species = factor(as.character(sigtab$Species), levels=names(x))

#Adjusting for Site and community
TreeOfCol <- factor(cutree(p$tree_col, k = 5))

sample_data(carbomVtaxa)$TreeOfCol <- TreeOfCol#tree of col is a coloumn generated from assignment of samples from hierarchical clustering
#WithCommunity
library("DESeq2")
otu_table(carbomVtaxa)<-otu_table(carbomVtaxa)+1
diagdds = phyloseq_to_deseq2(carbomVtaxa, ~ Labor_Fever+TreeOfCol+SSITE)##FebrileVSAfebrile
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#resultsNames(diagdds)
res= results(diagdds, contrast=c("Labor_Fever","Febrile","Afebrile"), cooksCutoff = FALSE)

alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(carbomVtaxa)[rownames(sigtab), ]
                                            , "matrix"))
head(sigtab)
View(sigtab)
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "phylum_colors", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#SpeciesLevel
x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
x = sort(x, TRUE)
sigtab$Species = factor(as.character(sigtab$Species), levels=names(x))


