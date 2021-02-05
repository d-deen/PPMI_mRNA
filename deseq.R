################################
########## PACKAGES ###########
###############################

library(tximport)
library(DESeq2)
library(AnnotationHub)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ensembldb)
library(RColorBrewer)
library(PCAtools)
library(pheatmap)
library(ggrepel)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(DEGreport)
library(Cairo)
library(ggpubr)
library(EnhancedVolcano)
library(tidyverse)

###### graph theme
graph_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20,face="bold")
  )

# import sample data sheet
samples <- read_csv("~/OneDrive - Newcastle University/PhD/PPMI/merged_BL.csv")

#### filter samples
# import flagged samples
flagged <- read_csv("~/OneDrive - Newcastle University/PhD/PPMI/flagged_samples.csv", col_names = FALSE)

# filter samples removing those where the clin diagnosis doesn't match imaging, those that have withdrawn from study
samples <- samples %>% 
  filter(clin_imag_agree == "Y" & ppmi_enroll != "Withdrew")

#vector of all filenames including path
files <- file.path("~/BL",samples$filename_long)
# rename files with patno
names(files) <- pull(samples, PATNO)

# creating tx2gene
txdb <- makeTxDbFromGFF(file="~/OneDrive - Newcastle University/PhD/PPMI/gencode.v19.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME") 

# tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)

# create metadata file and check both files match
metadata <- data.frame(samples, row.names = colnames(txi$counts))
all(colnames(txi$counts) %in% rownames(metadata))
all(colnames(txi$counts) == rownames(metadata))

# convert sex to character
metadata$sex <- as.character(as.numeric(metadata$sex))

# create DESeqDataSet
dds <- DESeqDataSetFromTximport(txi,
                                colData = metadata,
                                design = ~ status_clin)

#########################################
#### EXPLORATORY DATA ANALYSIS ##########
#########################################

# transform counts via vst and run pca
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

##### principal component analysis

# run pca
pca <- prcomp(t(vsd_mat))

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
percentVar <- percentVar * 100
percentVar <- round(percentVar, digits = 2)

# combine with metadata
p <- pca(vsd_mat, metadata = metadata)

# scree plot
screeplot(p, axisLabSize = 18, titleLabSize = 22)

### PCA biplots
df <- cbind(metadata, pca$x)
PCA_plot <- ggplot(df) + geom_point(aes(x = PC1, y = PC2, color = sex), size = 5, alpha = 0.7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  graph_theme
PCA_plot


##### HIERARCHICAL CLUSTERING
# compute pairwise correlation
vsd_cor <- cor(vsd_mat)
head(vsd_cor)


# plot as heatmap
df <- as.data.frame(colData(dds)[,"status_clin"])
rownames(df) <- colnames(vsd_cor)
heat.colors <- brewer.pal(6, "Blues")
pheatmap(vsd_cor, 
         color = heat.colors, 
         annotation_col = df,
         annotation_names_col = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)

## DIFFERENTIAL GENE EXPRESSION ANALYSIS
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
results <- results(dds)

# Get normalised counts and save to file
normalized_counts <- counts(dds, normalized = TRUE)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)
