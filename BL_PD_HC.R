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
library(viridis) 
library(tidyverse)

###### graph theme
graph_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15)
  )

# import sample data sheet
samples <- read.csv("~/OneDrive - Newcastle University/PhD/PPMI/merged_BL.csv", header = T)

# filter samples removing those where the clin diagnosis doesn't match imaging, those that have withdrawn from study
samples <- samples %>% 
  filter(status_clin == "PD" | status_clin == "HC") %>% 
  filter(EVENT_ID == "BL" & clin_imag_agree == "Y" & ppmi_enroll != "Withdrew")
  

#vector of all filenames including path
files <- file.path("~/BL/",samples$filename_long)
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

# convert sex and site to character
metadata$sex <- as.character(as.numeric(metadata$sex))
metadata$SITE <- sapply(metadata$SITE, as.factor)

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
screeplot(p, axisLabSize = 8, titleLabSize = 22) +
  graph_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
    )


  

### PCA biplots
# condition
df <- cbind(metadata, pca$x)
PCA_condition_1_2 <- ggplot(df) + geom_point(aes(x = PC1, y = PC2, color = status_clin), size = 5, alpha = 0.75) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_brewer(palette = "Set1") +
  graph_theme
PCA_condition_1_2

df <- cbind(metadata, pca$x)
PCA_condition_3_4 <- ggplot(df) + geom_point(aes(x = PC3, y = PC4, color = status_clin), size = 5, alpha = 0.75) +
  xlab(paste0("PC3: ", percentVar[3], "% variance")) +
  ylab(paste0("PC4: ", percentVar[4], "% variance")) +
  scale_color_brewer(palette = "Set1") +
  graph_theme
PCA_condition_3_4

ggarrange(PCA_condition_1_2, PCA_condition_3_4,
          common.legend = T
          )

# sex
df <- cbind(metadata, pca$x)
PCA_sex_1_2 <- ggplot(df) + geom_point(aes(x = PC1, y = PC2, color = sex), size = 5, alpha = 0.75) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_brewer(palette = "Set1") +
  graph_theme
PCA_sex_1_2

df <- cbind(metadata, pca$x)
PCA_sex_3_4 <- ggplot(df) + geom_point(aes(x = PC3, y = PC4, color = sex), size = 5, alpha = 0.75) +
  xlab(paste0("PC3: ", percentVar[3], "% variance")) +
  ylab(paste0("PC4: ", percentVar[4], "% variance")) +
  scale_color_brewer(palette = "Set1") +
  graph_theme
PCA_sex_3_4

ggarrange(PCA_sex_1_2, PCA_sex_3_4,
          common.legend = T
)

# age 
df <- cbind(metadata, pca$x)
PCA_age_1_2 <- ggplot(df) + geom_point(aes(x = PC1, y = PC2, color = age_bl), size = 5, alpha = 0.75) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_colour_viridis() +
  graph_theme
PCA_age_1_2

df <- cbind(metadata, pca$x)
PCA_age_3_4 <- ggplot(df) + geom_point(aes(x = PC3, y = PC4, color = age_bl), size = 5, alpha = 0.75) +
  xlab(paste0("PC3: ", percentVar[3], "% variance")) +
  ylab(paste0("PC4: ", percentVar[4], "% variance")) +
  scale_colour_viridis() +
  graph_theme
PCA_age_3_4

ggarrange(PCA_age_1_2, PCA_age_3_4,
          common.legend = T
)


# site 
df <- cbind(metadata, pca$x)
PCA_site_1_2 <- ggplot(df) + geom_point(aes(x = PC1, y = PC2, color = SITE), size = 5, alpha = 0.75) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_colour_viridis(discrete = T) +
  graph_theme
PCA_site_1_2

df <- cbind(metadata, pca$x)
PCA_site_3_4 <- ggplot(df) + geom_point(aes(x = PC3, y = PC4, color = SITE), size = 5, alpha = 0.75) +
  xlab(paste0("PC3: ", percentVar[3], "% variance")) +
  ylab(paste0("PC4: ", percentVar[4], "% variance")) +
  scale_colour_viridis(discrete= T) +
  graph_theme
PCA_site_3_4

ggarrange(PCA_site_1_2, PCA_site_3_4,
          common.legend = T
)


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

# plot disperal estimates
plotDispEsts(dds)

## building results table
contrast <- c("status_clin", "PD", "HC")
results_sig <- results(dds,
                   contrast = contrast,
                   alpha = 0.05)
                   #lfcThreshold = 0.58)

results_sig %>% data.frame() %>% View()

results_shrink <- lfcShrink(dds,
                     contrast = contrast,
                     res = results_sig,
                     type = "ashr")

plotMA(results_shrink)