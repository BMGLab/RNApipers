#Downstream analysis for Sven's Head and Neck cancer RNAseq
library(tidyverse)

setwd("/Users/yasinkaymaz/Documents/data/Sven/results.1")
pcg <- read.delim("protein_coding_genes.Ensembl.ids.txt", header = F)
head(pcg)

tpm <- read.delim("salmon.merged.gene_tpm.tsv", header = T, row.names = 1)
head(tpm)
tpm %>% select(-gene_name) -> tpm

table(rownames(tpm) %in% pcg$V1)

tpm[which(rownames(tpm) %in% pcg$V1),] -> tpm

tpm[which(rowMeans(tpm) > 1.0),] -> tpm

tpm[order(apply(tpm, 1, sd, na.rm=TRUE) / apply(tpm, 1, mean, na.rm=TRUE), decreasing = T),] %>% head(2000) -> tpm.top.cv

pheatmap(cor(tpm.top.cv))
head(tpm)
tpm.top.cv %>% reshape2::melt() %>% ggplot(aes(variable, log(value+.01) )) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Counts
sscounts <- read.table("salmon.merged.gene_counts.tsv", row.names = 1, header = T)
sscounts %>% select(gene_name) -> geneSymbol.dic

sscounts %>% select(-gene_name) -> cts
names(cts)
#Keep only Protein coding genes:
cts[which(rownames(cts) %in% pcg$V1),] -> cts
dim(cts)
#Convert to integers
cts <- round(cts, digits = 0)

cts[which(rowMeans(cts) > 10.0),] -> cts

coldata <- read.delim("meta.data.table.txt", row.names = 1)
coldata
coldata$Condition <- factor(coldata$Condition)
coldata$PatientID <- factor(coldata$PatientID)


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Condition)
dds
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 100) >= smallestGroupSize
dds <- dds[keep,]

dds$Condition <- factor(dds$Condition, levels = c("Post-treatment", "Pre-treatment"))
dds$Condition <- relevel(dds$Condition, "Pre-treatment")

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
summary(res)
#res[order(res$padj, decreasing = F),] %>% head(350)
resultsNames(dds)

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Condition", 
                returnData=TRUE)
#library("ggplot2")
ggplot(d, aes(x=Condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) + ggtitle(geneSymbol.dic[rownames(dds[which.min(res$padj),]),])

#res <- results(dds, name="Condition_Pre.treatment_vs_Post.treatment")
#res <- results(dds, contrast=c("condition","treated","untreated"))
resLFC <- lfcShrink(dds, coef="Condition_Pre.treatment_vs_Post.treatment", type="apeglm")
resLFC
resOrdered <- res[order(res$pvalue),]
summary(res)


library("RColorBrewer")
library(pheatmap)
vsd <- vst(dds, blind=FALSE)

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$PatientID, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("Condition"))


sgenes <- resOrdered[order(abs(resOrdered$log2FoldChange),decreasing = T),] %>% head(500) %>% rownames()


rownames(dds) <- geneSymbol.dic[rownames(dds),"gene_name"]

rld <- vst(dds, blind=FALSE)

de_mat <- assay(rld)[sgenes,]
pheatmap(t(scale(t(de_mat))),show_rownames = F,show_colnames = T,annotation_col =coldata)


