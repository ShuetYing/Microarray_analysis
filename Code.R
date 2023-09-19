# load the libraries
library(affy)
library(limma)

# read annotation data 
annotation <-read.table("Annotation.txt",header=T,as.is=T)
row.names(annotation)<-annotation$Name

# read CEL files 
data <- ReadAffy(filenames=annotation$Filename)
sampleNames(data) <- annotation$Name

# check quality of the data 
# log intensity plot
png("log_intensity_plot.png")
hist(data, main="Log intensity plot", ylab="Density", xlab="Log intensity")
dev.off()

# box plot 
png("boxplot_before_normalization.png")
par(mar=c(5, 5, 4, 2) + 0.7)
boxplot(data, col=annotation$Color, names=annotation$Name, las=2, main="Boxplot of signal intensity before normalization", ylab="Log intensity", xlab="")
mtext("Samples", side=1, line=4, padj=1)
dev.off()

#scatter plot
png("Mn1wt_vs_Mn1null_before_normalization.png")
Mn1wt <- exprs(data[,c("Mn1wt.1", "Mn1wt.2", "Mn1wt.3")])
Mn1null <- exprs(data[,c("Mn1null.1", "Mn1null.2", "Mn1null.3")])
plot(Mn1wt, Mn1null, pch=1, col=c("red","blue"), xlab="Mn1wt cell", ylab="Mn1null cell", main="Mn1wt cell vs Mn1null cell")
legend("topright", legend=c("Mn1wt", "Mn1null"), pch=1, col=c("red","blue"))
dev.off()
Mn1wt <- log2(Mn1wt)
Mn1null <- log2(Mn1null)
png("log2_Mn1wt_vs_Mn1null_before_normalization.png")
plot(Mn1wt, Mn1null, pch=1, col=c("red","blue"), xlab="log2 Mn1wt cell", ylab=" log2 Mn1null cell", main="Mn1wt cell vs Mn1null cell")
abline(-1, 1, col="black", lty="dashed")
abline(1, 1, col="black", lty="dashed")
dev.off()

# normalize data using Robust Multichip Average (RMA)
eset <- rma(data)

# annotate the data 
library(mouse4302.db)
sampleNames(eset)<-annotation$Name
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "mouse4302.db")
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, stringsAsFactors=F)
tmp[tmp=="NA"] <- NA
fData(eset) <- tmp

# extract normalized expression values
values <- exprs(eset) 

# box plot
png("boxplot_after_normalization.png")
par(mar=c(5, 5, 4, 2) + 0.7)
boxplot(values, col=annotation$Color, names=annotation$Name, las=2, main="Boxplot of signal intensity after normalization", ylab="Log intensity", xlab="")
mtext("Samples", side=1, line=4, padj=1)
dev.off()

#scatter plot
png("Mn1wt_vs_Mn1null_after_normalization.png")
Mn1wt <- exprs(eset[,c("Mn1wt.1", "Mn1wt.2", "Mn1wt.3")])
Mn1null <- exprs(eset[,c("Mn1null.1", "Mn1null.2", "Mn1null.3")])
plot(Mn1wt, Mn1null, pch=1, col=c("red","blue"), xlab="log2 Mn1wt cell", ylab="log2 Mn1null cell", main="Mn1wt cell vs Mn1null cell")
abline(-1, 1, col="black", lty="dashed")
abline(1, 1, col="black", lty="dashed")
legend("topright", legend=c("Mn1wt", "Mn1null"), pch=1, col=c("red","blue"))
dev.off()

# quality control for normalized data 
png("MvA_plot.png")
mva.pairs(values, cex=0.8, main="MvA plot of gene expression data after normalization")
dev.off()

# assess variation between replicates in Mn1wt using correlation
rep1 <- exprs(data[,"Mn1wt.1"])
rep2 <- exprs(data[,"Mn1wt.2"])
rep3 <- exprs(data[,"Mn1wt.3"])
correlation <- cor(cbind(rep1, rep2, rep3))
write.csv(correlation,file="correlation.csv")

# scatter plot btw replicates of Mn1wt
png("replicate_1_&_2.png")
plot(rep1, rep2, pch=1, xlab="Mn1wt.1", col="black", ylab="Mn1wt.2", main="Mn1wt.1 vs Mn1wt.2")
abline(0, 1, col="red")
dev.off()
png("replicate_1_&_3.png")
plot(rep1, rep3, pch=1, xlab="Mn1wt.1", col="black", ylab="Mn1wt.3", main="Mn1wt.1 vs Mn1wt.3")
abline(0, 1, col="red")
dev.off()
png("replicate_2_&_3.png")
plot(rep2, rep3, pch=1, xlab="Mn1wt.2", col="black", ylab="Mn1wt.3", main="Mn1wt.2 vs Mn1wt.3")
abline(0, 1, col="red")
dev.off()

# perform hierarchical clustering 
colnames(values) <- annotation$Name # assign sample name as col name 
# convert corr matrix into distance matrix and calculate the distances btw clusters 
hc<-hclust(as.dist(1-cor(values, method="pearson")), method="average")
# output=dendrogram, samples with similar expression patterns will be grouped together into clusters
png("dendogram.png")
plot(hc)
dev.off()

# perform Principle Components Analysis (PCA)
library(scatterplot3d)
pca <- prcomp(t(values), scale=T) # work on transpose of the matrix 
# coordinates of the plotted points 
png("PCA.png")
s3d<-scatterplot3d(pca$x[,1:3], pch=19, color=annotation$Color) # the first 3 PCs usually capture the largest amount of variation 
# convert into screen coordinates
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels=colnames(values), pos=3, offset=0.5)
dev.off()

# convert expression values into linear scale 
values <- exprs(eset)
exprsvals10 <- 2^values

# calculates the ratio of the mean expression levels between samples 
wt.mean <- apply(exprsvals10[,c("Mn1wt.1", "Mn1wt.2", "Mn1wt.3")], 1, mean)
null.mean <- apply(exprsvals10[,c("Mn1null.1", "Mn1null.2", "Mn1null.3")], 1, mean)
fold_change <- wt.mean/null.mean
write.csv(fold_change,file="group_means.csv")

# build design matrix 
design <- model.matrix(~-1+factor(c(1,1,1,2,2,2)))
colnames(design) <- c("Mn1wt", "Mn1null")

# build contrast matrix 
contrastmatrix <- makeContrasts(Mn1wt-Mn1null, levels=design)

# fit the model 
fit_model <- lmFit(eset, design) 
# calculate the log-fold changes and standard errors for each contrast specified in the contrast matrix
fit_contrast <- contrasts.fit(fit_model, contrastmatrix)
# perform empirical Bayes moderation of the standard errors
fit_contrast <- eBayes(fit_contrast)
# look at the summary at log fold change of 1 (2 fold change in expression level in either direction)
summary(decideTests(fit_contrast, lfc=1))

# obtain statistical summary of all the differentially expressed genes
results <- topTable(fit_contrast, adjust="fdr", number=nrow(eset))
write.csv(results,"statistical_results.csv")

# obtain top 10 differentially expressed genes 
# genes by log odds
top_genes_limma <- topTable(fit_contrast, n=10)
write.csv(top_genes_limma,"top_10_genes_by_limma.csv")
# genes by fold change
top_genes_fc <- topTable(fit_contrast, sort.by="logFC")
top_genes_fc <- top_genes_fc[order(-top_genes_fc$logFC),]
top_genes_fc <- head(top_genes_fc, 10)
write.csv(top_genes_fc,"top_10_genes_by_fc.csv")

# determine the enriched pathway for the top 10 differentially expressed genes identified by limma
library(clusterProfiler)
library(org.Mm.eg.db)
limma_genes <- top_genes_limma$SYMBOL
limma_result <- enrichGO(limma_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05, universe = NULL)
png("limma_function.png")
dotplot(limma_result)
dev.off

# determine the enriched pathway for the top 10 differentially expressed genes identified by logFC
logFC_genes <- top_genes_fc$SYMBOL
logFC_result <- enrichGO(logFC_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05, universe = NULL)
png("logFC_function.png")
dotplot(logFC_result)
dev.off()

# obtain annotation for the probes 
res <- select(mouse4302.db, keys = rownames(eset), columns = c("ENTREZID", "ENSEMBL","SYMBOL"), keytype="PROBEID")
idx <- match(rownames(eset), res$PROBEID)
fData(eset) <- res[idx, ]
head(fData(eset), 10)

# remove probe w/o ENTREZID
eset_removed <- eset[is.na(fData(eset)$ENTREZID)==0,]

# convert ENTREZID to indices in Mm.h
Mm.H <- readRDS("Mm.h.all.v7.1.entrez.rds")
H.indices <- ids2indices(Mm.H,fData(eset_removed)$ENTREZID)
# use Camera to run functional enrichment analysis 
results <- camera(eset_removed, index=H.indices, design=design, contrast=contrastmatrix)
write.csv(results,"functional_enrichment.csv")
# find the row with HALLMARK_APOPTOSIS
which(rownames(results) == "HALLMARK_APOPTOSIS")
apoptosis <- results[34,]
write.csv(apoptosis,"functional_enrichment_apoptosis.csv")

# obtain differentially expressed genes from top table & top treat using default method: BH
top_table <- topTable(fit_contrast, number=nrow(eset))
treat_fit_model <- lmFit(eset, design)
treat_fit_contrast <-contrasts.fit(treat_fit_model, contrastmatrix)
treat_fit_contrast <- treat(treat_fit_contrast)
top_treat <- topTreat(treat_fit_contrast, number=nrow(eset))
# filter by fold-change and FDR
filtered_top_table <- top_table[abs(top_table$logFC) > 1 & top_table$adj.P.Val < 0.01, ]
write.csv(filtered_top_table, "top_table.csv")
filtered_top_treat <- top_treat[abs(top_treat$logFC) > 1 & top_treat$adj.P.Val < 0.01, ]
write.csv(filtered_top_treat, "top_treat.csv")
# find the difference between the two methods 
library(dplyr)
diff <- anti_join(filtered_top_table, filtered_top_treat, by = "SYMBOL")
write.csv(diff, "diff_btw_table_treat.csv") # return genes that are not found in the top treat table

# volcano plot
# code modified from GTK-teaching on GitHub: https://gtk-teaching.github.io/Microarrays-R/09-downstreamAnalysis/index.html
png("volcano_plot.png")
volcanoplot(fit_contrast, style = "p-value", highlight = 0, names = fit_contrast$genes$SYMBOL, xlab = "Log2 Fold Change", ylab = "-log10(p-value)", main = "Volcano plot", pch=16, cex=0.45)
dev.off()
png("volcano_plot_with_DEGs.png")
volcanoplot(fit_contrast, style = "p-value", highlight = 0, names = fit_contrast$genes$SYMBOL, xlab = "Log2 Fold Change", ylab = "-log10(p-value)", main = "Volcano plot", pch=16, cex=0.45)
DEGs <- topTable(fit_contrast, p.value = 0.05, lfc = 4, number = Inf)
points(DEGs$logFC,-log10(DEGs$P.Value),col='red')
dev.off()
