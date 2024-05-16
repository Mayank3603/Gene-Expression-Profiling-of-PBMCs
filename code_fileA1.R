library(GEOquery)
library(ggplot2)
library(fgsea)
library(limma)
library(umap)


# download the microarray data
gset <- getGEO("GSE48556", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6947", attr(gset, "names")) else idx <- 1
gset<- (gset[[idx]])


# normailze the dataset
exprs(gset) <- normalizeBetweenArrays(exprs(gset)) 


# pdata and fdata and attributes.
pdata <- pData(gset)
fdata <- fData(gset)
sampleNames(gset)

print(pdata)
print(fdata)



# EDA and pre processing
# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "000000111111111111111111111111111111111")
sml <- strsplit(gsms, split="")[[1]]

ex <- exprs(gset)
ex[which(ex <= 0)] <- NaN
# log2 transformation
exprs(gset) <- log2(ex) 

exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data

# assign samples to groups and set up design matrix
gs <- factor(sml)

groups <- make.names(c("osteoarthritis","control"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

nall <- nrow(gset)
gset <- gset[complete.cases(exprs(gset)), ]



# t_test , holm correction and holm test
t_test <- apply(gset, 1, function(x) t.test(x[gs=="osteoarthritis"], x[gs=="control"]))


p_value <- sapply(t_test, function(x) x$p.value)

log_FoldChange<- apply(gset, 1, function(x) log2(mean(x[gs=="osteoarthritis"])/mean(x[gs=="control"])))

# using Holm correction we correct the p-value
holm_correction <- p.adjust(p_value, method = "holm")

# Create a volcano plot to visualize the results
plot(log_FoldChange, -log10(holm_correction), pch=20, main="volcano_plot", xlab="log2 fold change", ylab="-log10(pvalue)",
     xlim=c(-6,6), ylim=c(0,10), col=ifelse(abs(log_FoldChange) > 2 & holm_correction < 0.05, "red", "black"))
abline(h=-log10(0.05), col="blue", lty=2)
abline(v=c(-2,2), col="blue", lty=2)



# Volcano plot using limma package

fit  <- lmFit(gset,design)

#recalculate model coefficients
vector1 <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=vector1, levels=design)
new_cordinate <- contrasts.fit(fit, cont.matrix)

# compuation of top significant values
new_cordinate <- eBayes(new_cordinate, 0.01)
top_table <- topTable(new_cordinate, adjust="holm", sort.by="B", number=250)
dT <- decideTests(new_cordinate, adjust.method="holm", p.value=0.05)
colnames(new_cordinate) # list contrast names

volcanoplot(new_cordinate, coef=1, main=colnames(new_cordinate)[1], pch=20, highlight=length(which(dT[,1]!=0)), names=rep('+', nrow(new_cordinate)))      

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

expressed_genes <- c(top_table$"Gene.symbol")
genes_mapped <- bitr(expressed_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)


# perform gene set enrichment analysis using enrichGO
ans <- enrichGO(
  gene=genes_mapped$ENTREZID,
  keyType = "ENTREZID",
  OrgDb=org.Hs.eg.db,
  ont="ALL",
  pvalueCutoff=0.05,
  qvalueCutoff=0.05,
  minGSSize=5, 
  maxGSSize=500, 
  
)


#Bar Graph
most_significant_10 <- ans[1:10]$pvalue
barplot(most_significant_10,names.arg =ans[1:10]$ID, las = 2, col = "blue", ylim = c(0, max(most_significant_10)*1.2))

# dot-plot
dotplot(ans, showCategory = 25)

# enrichment map plot
pairwise_termsim_res <- pairwise_termsim(ans)
emapplot(pairwise_termsim_res, enrich=res, showCategory=25)

# category netplot
cnetplot(ans, foldChange = de_genes_ranked, showCategory = 25)

# To view different pathways in gene enrichment analysis
View(ans@result)


