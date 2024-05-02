library(dplyr)
library(edgeR)
library(tidyverse)
library(pheatmap)
#install.packages("gplots")
library(gplots)
#install.packages("RColorBrewer")
library(RColorBrewer)
library(limma)

setwd("/Users/edyilmaz/dn")

#configFiles
pcg_file <- "/Users/edyilmaz/dn/protein_coding_genes.Ensembl.ids.txt"
tpm_file <- "/Users/edyilmaz/dn/salmon.merged.gene_tpm.tsv"
sscounts_file <- "/Users/edyilmaz/dn/salmon.merged.gene_counts.tsv"
coldata_file <- "/Users/edyilmaz/dn/meta.data.table.txt"
#list: (gene_id	gene_name	X24.004.1_S56	X24.004.10_S20	X24.004.11_S18	X24.004.2_S52	X24.004.3_S41	X24.004.4_S27	X24.004.5_S19	X24.004.6_S17	X24.004.7_S14	X24.004.8_S53	X24.004.9_S42	X24004.12_S15)
pcg <- read.delim(pcg_file, header = F)
tpm <- read.delim(tpm_file, header = T, row.names = 1)
tpm %>% select(-gene_name) -> tpm

table(rownames(tpm) %in% pcg$V1)

tpm[which(rownames(tpm) %in% pcg$V1),] -> tpm
tpm[which(rowMeans(tpm) > 1.0),] -> tpm
tpm[order(apply(tpm, 1, sd, na.rm=TRUE) / apply(tpm, 1, mean, na.rm=TRUE), decreasing = T),] %>% head(2000) -> tpm.top.cv


pdf(file="Rplot_Figure1_heatmap.pdf")
pheatmap(cor(tpm.top.cv))
dev.off()

pdf(file="Rplot_Figure2_boxplot.pdf")
tpm.top.cv %>% reshape2::melt() %>% ggplot(aes(variable, log(value+.01) )) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#Counts
sscounts <- read.table(sscounts_file, row.names = 1, header = T)
sscounts %>% select(gene_name) -> geneSymbol.dic
sscounts %>% select(-gene_name) -> cts
names(cts)

#Keep only Protein coding genes:
cts[which(rownames(cts) %in% pcg$V1),] -> cts
dim(cts)

#Convert to integers
cts <- round(cts, digits = 0)
cts[which(rowMeans(cts) > 10.0),] -> cts
coldata <- read.delim(coldata_file, row.names = 1)

coldata$Condition <- factor(coldata$Condition)
coldata$PatientID <- factor(coldata$RawDataID)
cts <- cts[, coldata$RawDataID]
cts


#EdgeR using:
d <- DGEList(counts = cts, group = (coldata$Condition))
d

dim(d)

d.full <- d # keep the old one in case we mess up
#head(d$counts)
#head(cpm(d))
apply(d$counts, 2, sum)

keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)

d$samples$lib.size <- colSums(d$counts)
d$samples

d <- calcNormFactors(d)
d

#needed regulation:
pdf("Rplot_Figure3_location.pdf")
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:4, pch=20)
dev.off()

#estimating the dispersion
pdf(file="Rplot_Figure4_estimatingDispersion.pdf")
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
plotBCV(d1)
dev.off()
#glm esitmates of dispersion

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
pdf(file="Rplot_Figure5_gmlEstimateBCV.pdf")
plotBCV(d2)
dev.off()

################## Differantial Expression Analysis

et24 <- exactTest(d1, pair=c(2,4)) # compare groups 2 and 4
et24t = topTags(et24,n=50)
write.csv(et24, "toptags_output_et24.csv", row.names=TRUE)
write.csv(et24t, "toptags_output_et24_wFDR.csv", row.names=TRUE)

et21 <- exactTest(d1, pair=c(2,1)) # compare groups 2 and 1
et21t = topTags(et21,n=50)
write.csv(et21, "toptags_output_et21.csv", row.names=TRUE)
write.csv(et21t, "toptags_output_et21_wFDR.csv", row.names=TRUE)

et43 <- exactTest(d1, pair=c(4,3)) # compare groups 4 and 3
et43t = topTags(et43,n=50)
write.csv(et43, "toptags_output_et43.csv", row.names=TRUE)
write.csv(et43t, "toptags_output_et43_wFDR.csv", row.names=TRUE)

et23 <- exactTest(d1, pair=c(2,3)) # compare groups 2 and 3
et23t = topTags(et23,n=50)
write.csv(et23, "toptags_output_et23.csv", row.names=TRUE)
write.csv(et23t, "toptags_output_et23_wFDR.csv", row.names=TRUE)

de1 <- decideTestsDGE(et24, adjust.method="BH", p.value=0.05)
summary(de1)
de1tags24 <- rownames(d1)[as.logical(de1)] 
plotSmear(et24, de.tags=de1tags24)
abline(h = c(-2, 2), col = "blue")

de2 <- decideTestsDGE(et21, adjust.method="BH", p.value=0.05)
summary(de2)
de2tags21 <- rownames(d2)[as.logical(de2)] 
plotSmear(et21, de.tags=de2tags21)
abline(h = c(-2, 2), col = "blue")

de3 <- decideTestsDGE(et23, adjust.method="BH", p.value=0.05)
summary(de3)
de3tags23 <- rownames(d3)[as.logical(de3)] 
plotSmear(et23, de.tags=de3tags23)
abline(h = c(-2, 2), col = "blue")

de4 <- decideTestsDGE(et43, adjust.method="BH", p.value=0.05)
summary(de4)
de4tags43 <- rownames(d4)[as.logical(de4)] 
plotSmear(et43, de.tags=de4tags43)
abline(h = c(-2, 2), col = "blue")


#glm testing for differential expression
dds <- d

library("RColorBrewer")
library(pheatmap)
voom_transformed_data <- voom(dds$counts, design.mat)
vsd <- voom_transformed_data

voom_exprs <- voom_transformed_data$E
sampleDists <- dist(t(voom_exprs))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(voom_exprs), sep="-")
colnames(sampleDistMatrix) <- NULL

pdf(file="Rplot_Figure6_MDSplot.pdf")
my_colors <- c("red", "blue", "darkgreen", "orange", "purple", "darkblue", "brown", "darkmagenta", "black", "darkred")
plotMDS(voom_exprs, labels = colnames(voom_exprs), col = my_colors, cex = 0.5)
par(cex.lab = 1)
title(main = "MDS Plot", cex.main = 1, font.main=2)
dev.off()

pdf(file="Rplot_Figure7_pheatmap.pdf")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

pdf(file="Rplot_Figure8_pcaPlot.pdf")
pca_result <- prcomp(t(voom_exprs))
plot(pca_result$x[,1], pca_result$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA Plot")
text(pca_result$x[,1], pca_result$x[,2], labels = colnames(voom_transformed_data$E), pos = 1, cex=0.6)
dev.off()

logcounts <- cpm(d,log=TRUE)
var_genes <- apply(logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:30]
head(select_var)

highly_variable_lcpm <- logcounts[select_var, ]
dim(highly_variable_lcpm)

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("red", "blue", "darkgreen", "orange")[d$samples$group]
# Plot the heatmap
pdf(file="Rplot_Figure9_lastHeatmap.pdf")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top variable genes across conditions Samples",ColSideColors=col.cell,scale="row",margins=c(9,9))
dev.off()

################################Voom Trans.

voom_transformed_data <- voom(counts = dds$counts, design = design.mat)
fit <- lmFit(voom_transformed_data, design=design.mat)
fit <- eBayes(fit)
top_tags <- topTable(fit, coef=1, n=500, sort.by="logFC")
sgenes <- rownames(top_tags)

highly_variable_lcpm2 <- logcounts[sgenes, ]
dim(highly_variable_lcpm2)

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("red", "blue", "darkgreen", "orange")[d$samples$group]
# Plot the heatmap
pdf(file="Rplot_Figure10_DE_heatmap.pdf")
heatmap.2(highly_variable_lcpm2,col=rev(morecols(50)),trace="none", main="Top variable genes across conditions Samples with Voom",ColSideColors=col.cell,scale="row",margins=c(9,9))
dev.off()