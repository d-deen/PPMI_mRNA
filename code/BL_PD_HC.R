

# PACKAGES ----------------------------------------------------------------

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
library(cowplot)
library(tidyverse)



# GRAPH THEME -------------------------------------------------------------

graph_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_text(size = 13, colour = "black"),
    axis.title = element_text(size = 13),
    aspect.ratio = 1
  )



# SETTING UP SAMPLE INFO AND METADATA -------------------------------------

samples <- read.csv("~/OneDrive - Newcastle University/PhD/PPMI/metadata.csv", header = T)

# filter samples removing those where the clin diagnosis doesn't match imaging, those that have withdrawn from study
samples_filtered <- samples %>% 
  filter(status_clin == "PD" | status_clin == "HC") %>% 
  filter(EVENT_ID == "BL" & clin_imag_agree == "Y" & ppmi_enroll != "Withdrew")

nrow(samples_filtered %>% 
       filter(status_clin == "PD"))       # number of PD samples

nrow(samples_filtered %>% 
       filter(status_clin == "HC"))       # number of HC samples
  

#vector of all filenames including path
files <- file.path("~/BL/",samples_filtered$file_path)
# rename files with patno
names(files) <- pull(samples_filtered, PATNO)

# creating tx2gene
txdb <- makeTxDbFromGFF(file="~/OneDrive - Newcastle University/PhD/PPMI/gencode.v19.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME") 
tx2gene$TXNAME <- sub("\\.\\d+", "", tx2gene$TXNAME)
tx2gene$GENEID <- sub("\\.\\d+", "", tx2gene$GENEID)


# tximport
txi <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM", 
                ignoreAfterBar = TRUE,
                ignoreTxVersion = TRUE)

# create metadata file and check both files match
metadata <- data.frame(samples_filtered, row.names = colnames(txi$counts))
all(colnames(txi$counts) %in% rownames(metadata))
all(colnames(txi$counts) == rownames(metadata))

# convert sex and site to character
metadata$sex <- as.character(as.numeric(metadata$sex))
metadata$SITE <- sapply(metadata$SITE, as.factor)

# create DESeqDataSet
dds <- DESeqDataSetFromTximport(txi,
                                colData = metadata,
                                design = ~ status_clin)



#### EXPLORATORY DATA ANALYSIS ------------------------------------------------------

# transform counts via vst and run pca
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)


## principal component analysis
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

## PCA biplots
# condition
df <- cbind(metadata, pca$x)
PCA_condition_1_2 <- ggplot(df) + geom_point(aes(x = PC1, y = PC2, color = status_clin), size = 1, alpha = 0.75) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_color_brewer(palette = "Set2") +
  graph_theme
PCA_condition_1_2

df <- cbind(metadata, pca$x)
PCA_condition_3_4 <- ggplot(df) + geom_point(aes(x = PC3, y = PC4, color = status_clin), size = 1, alpha = 0.75) +
  xlab(paste0("PC3: ", percentVar[3], "% variance")) +
  ylab(paste0("PC4: ", percentVar[4], "% variance")) +
  coord_fixed() +
  scale_color_brewer(palette = "Set2") +
  graph_theme
PCA_condition_3_4

PCA_condition <- ggarrange(PCA_condition_1_2, PCA_condition_3_4,
          common.legend = T, 
          legend = "top",
          align = "hv"
          )

ggsave("figures/PCA_condition.png")

# sex
df <- cbind(metadata, pca$x)
PCA_sex_1_2 <- ggplot(df) + geom_point(aes(x = PC1, y = PC2, color = sex), size = 1, alpha = 0.75) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_brewer(palette = "Set1") +
  coord_fixed() +
  graph_theme
PCA_sex_1_2

df <- cbind(metadata, pca$x)
PCA_sex_3_4 <- ggplot(df) + geom_point(aes(x = PC3, y = PC4, color = sex), size = 1, alpha = 0.75) +
  xlab(paste0("PC3: ", percentVar[3], "% variance")) +
  ylab(paste0("PC4: ", percentVar[4], "% variance")) +
  scale_color_brewer(palette = "Set1") +
  coord_fixed() +
  graph_theme
PCA_sex_3_4

PCA_sex <- ggarrange(PCA_sex_1_2, PCA_sex_3_4,
          common.legend = T,
          legend = "bottom",
          align = "hv"
)

ggsave("figures/PCA_sex.png")

# age 
df <- cbind(metadata, pca$x)
PCA_age_1_2 <- ggplot(df) + geom_point(aes(x = PC1, y = PC2, color = age_bl), size = 1, alpha = 0.75) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_colour_viridis() +
  coord_fixed() +
  graph_theme
PCA_age_1_2

df <- cbind(metadata, pca$x)
PCA_age_3_4 <- ggplot(df) + geom_point(aes(x = PC3, y = PC4, color = age_bl), size = 1, alpha = 0.75) +
  xlab(paste0("PC3: ", percentVar[3], "% variance")) +
  ylab(paste0("PC4: ", percentVar[4], "% variance")) +
  scale_colour_viridis() +
  coord_fixed() +
  graph_theme
PCA_age_3_4

PCA_age <- ggarrange(PCA_age_1_2, PCA_age_3_4,
          common.legend = T,
          legend = "right",
          align = "hv"
)

ggsave("figures/PCA_age.png")

# site 
df <- cbind(metadata, pca$x)
PCA_site_1_2 <- ggplot(df) + geom_point(aes(x = PC1, y = PC2, color = SITE), size = 1, alpha = 0.75) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_colour_viridis(discrete = T) +
  coord_fixed() +
  graph_theme
PCA_site_1_2

df <- cbind(metadata, pca$x)
PCA_site_3_4 <- ggplot(df) + geom_point(aes(x = PC3, y = PC4, color = SITE), size = 1, alpha = 0.75) +
  xlab(paste0("PC3: ", percentVar[3], "% variance")) +
  ylab(paste0("PC4: ", percentVar[4], "% variance")) +
  scale_colour_viridis(discrete= T) +
  coord_fixed() +
  graph_theme
PCA_site_3_4

PCA_site <- ggarrange(PCA_site_1_2, PCA_site_3_4,
          common.legend = T,
          legend = "right",
          align = "hv"
)

ggsave("figures/PCA_site.png")

# combine condition and sex PC plots
PCA_combined_condition_sex <- ggarrange(PCA_condition, PCA_sex, 
                                        nrow = 2,
                                        align = "v")
PCA_combined_condition_sex
ggsave("figures/PCA_combined_condition_sex.png")

## HIERARCHICAL CLUSTERING
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


# DIFFERENTIAL EXPRESSION ANALYSIS ----------------------------------------
# create DESeqDataSet
dds <- DESeqDataSetFromTximport(txi,
                                colData = metadata,
                                design = ~ status_clin)

dds <- DESeq(dds)

# library sizes
library_sizes <- colSums(counts(dds))

library_sizes <- as.data.frame(library_sizes)
library_sizes <- library_sizes %>% 
  rownames_to_column(var = "patno")
  
lib_bar <- ggplot(library_sizes, aes(x = patno, y = library_sizes)) +
  geom_col(colour = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Sample") +
  ylab("Total number of counts") +
  graph_theme +
  theme(axis.text.x = element_blank(),
        aspect.ratio = 1)
lib_bar

ggsave("figures/library_size.png")

# write unfiltered normalized counts to file
normalized_counts_unfiltered <- counts(dds, normalized = TRUE)
write.csv(normalized_counts_unfiltered, "output/normalized_counts_unfiltered.csv")
view(normalized_counts_unfiltered)

# average gene count
ncu_mean_gene <- as.data.frame(colMeans(normalized_counts_unfiltered))
ncu_mean_gene <- ncu_mean_gene %>% 
  rownames_to_column (var = "patno") %>% 
  rename(means = `colMeans(normalized_counts_unfiltered)`)

# number of genes
ncu_observations <- as.data.frame(colSums(normalized_counts_unfiltered != 0))
ncu_observations <- ncu_observations %>% 
  rownames_to_column (var = "patno") %>% 
  rename(observations = `colSums(normalized_counts_unfiltered != 0)`)

# library size
ncu_sum_gene <- as.data.frame(colSums(normalized_counts_unfiltered))
ncu_sum_gene <- ncu_sum_gene %>% 
  rownames_to_column (var = "patno") %>% 
  rename(lib_size = `colSums(normalized_counts_unfiltered)`)

# merge datasets togther
df <- merge(ncu_mean_gene, ncu_observations, by = "patno")
df <- merge(df, ncu_sum_gene, by = "patno")

patno.status_clin <- samples_filtered %>% 
  select(PATNO, status_clin) %>% 
  rename(patno = PATNO) 
patno.status_clin[] <- lapply(patno.status_clin, as.character)

gene_counts <- left_join(df, patno.status_clin, by = "patno")

# plot no of genes vs mean gene count
a <- ggplot(gene_counts, aes(x = observations, y = means)) +
  geom_point(shape = 1) +
  geom_smooth(method=lm) +
  xlab("Number of genes detected") +
  ylab("Average gene count") +
  graph_theme
a

# plot library size vs no of genes
b <- ggplot(gene_counts, aes(x = observations, y = lib_size)) +
  geom_point(aes(colour = status_clin)) +
  geom_smooth(method=lm) +
  stat_cor(p.accuracy = 0.001) +
  xlab("Number of genes detected") +
  ylab("Total number of counts") +
  scale_color_brewer(palette = "Set2") +
  graph_theme +
  theme(legend.position = "bottom")
b

# plot lib size vs mean gene count
c <- ggplot(gene_counts, aes(x = means, y = lib_size)) +
  geom_point(shape = 1) +
  xlab("Average gene count") +
  ylab("Total number of counts") +
  stat_cor(p.accuracy = 0.001) +
  graph_theme +
  theme(aspect.ratio = 1)
c

lib_scatter <- plot_grid(c, b, nrow = 2, align = "h", axis = "b")

plot_grid(lib_bar, lib_scatter, ncol = 2, align = "hv")
ggsave("figures/library size plots.pdf")

# pre filtering
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

dds <- DESeq(dds)
results <- results(dds)

# Get normalised counts and save to file
normalized_counts <- counts(dds, normalized = TRUE)
#write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# plot disperal estimates
plotDispEsts(dds)


## building results table
contrast <- c("status_clin", "PD", "HC")
results <- results(dds,
                   contrast = contrast,
                   alpha = 0.1)
                   
                   #lfcThreshold = 0.58)

#results %>% data.frame() %>% View()
#results_csv <- results %>% data.frame
#write.csv(results_csv, "results.csv")

results_shrink <- lfcShrink(dds,
                     contrast = contrast,
                     res = results,
                     type = "ashr")

plotMA(results_shrink)

# histogram of p values
use <- results$baseMean > metadata(results)$filterThreshold
h1 <- hist(results$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(results$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

## investigating independent filtering
metadata(results)$filterThreshold



### plot counts of ICICLE genes

# HLA-E
HLA_E.ppmi <- plotCounts(dds, gene = "ENSG00000204592", intgroup = "status_clin", returnData = T)
HLA_E.ppmi <- HLA_E.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(trans = "log2") +
    ggtitle("HLA-E") +
    ylab("Log2 Normalised Counts") +
    xlab("Disease State") +
    graph_theme +
    theme(legend.position = "none",
          aspect.ratio = 1, 
          axis.title.x = element_blank(),
          plot.title = element_text(size = 8),
          axis.title = element_text(size = 5),
          axis.text = element_text(size = 5))

# CLK2
CLK2.ppmi <- plotCounts(dds, gene = "ENSG00000176444", intgroup = "status_clin", returnData = T)
CLK2.ppmi <- CLK2.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(trans = "log2") +
  ggtitle("CLK2") +
  ylab("Log2 Normalised Counts") +
  xlab("Disease State") +
  graph_theme +
  theme(legend.position = "none",
        aspect.ratio = 1, 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5))


# AGER
AGER.ppmi <- plotCounts(dds, gene = "ENSG00000204305", intgroup = "status_clin", returnData = T)
AGER.ppmi <- AGER.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(trans = "log2") +
  ggtitle("AGER") +
  ylab("Log2 Normalised Counts") +
  xlab("Disease State") +
  graph_theme +
  theme(legend.position = "none",
        aspect.ratio = 1, 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5))

#VARS1 - not found

#HLA-DMB
HLA_DMB.ppmi <- plotCounts(dds, gene = "ENSG00000242574", intgroup = "status_clin", returnData = T)
HLA_DMB.ppmi <- HLA_DMB.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(trans = "log2") +
  ggtitle("HLA-DMB") +
  ylab("Log2 Normalised Counts") +
  xlab("Disease State") +
  graph_theme +
  theme(legend.position = "none",
        aspect.ratio = 1, 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5))

#DEFA1B
DEFA1B.ppmi <- plotCounts(dds, gene = "ENSG00000240247", intgroup = "status_clin", returnData = T)
DEFA1B.ppmi <- DEFA1B.ppmi <- DEFA1B.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(trans = "log2") +
  ggtitle("DEFA1B") +
  ylab("Log2 Normalised Counts") +
  xlab("Disease State") +
  graph_theme +
  theme(legend.position = "none",
        aspect.ratio = 1, 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5))

#MUC5B
MUC5B.ppmi <- plotCounts(dds, gene = "ENSG00000117983", intgroup = "status_clin", returnData = T)
MUC5B.ppmi <- MUC5B.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(trans = "log2") +
  ggtitle("MUC5B") +
  ylab("Log2 Normalised Counts") +
  xlab("Disease State") +
  graph_theme +
  theme(legend.position = "none",
        aspect.ratio = 1, 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5))

#C6orf47
C6orf47.ppmi <- plotCounts(dds, gene = "ENSG00000204439", intgroup = "status_clin", returnData = T)
C6orf47.ppmi <- C6orf47.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(trans = "log2") +
  ggtitle("C6orf47") +
  ylab("Log2 Normalised Counts") +
  xlab("Disease State") +
  graph_theme +
  theme(legend.position = "none",
        aspect.ratio = 1, 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5))

#GATAD2B
GATAD2B.ppmi <- plotCounts(dds, gene = "ENSG00000143614", intgroup = "status_clin", returnData = T)
GATAD2B.ppmi <- GATAD2B.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(trans = "log2") +
  ggtitle("GATAD2B") +
  ylab("Log2 Normalised Counts") +
  xlab("Disease State") +
  graph_theme +
  theme(legend.position = "none",
        aspect.ratio = 1, 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5))

# DHX16
DHX16.ppmi <- plotCounts(dds, gene = "ENSG00000204560", intgroup = "status_clin", returnData = T)
DHX16.ppmi <- DHX16.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(trans = "log2") +
  ggtitle("DHX16") +
  ylab("Log2 Normalised Counts") +
  xlab("Disease State") +
  graph_theme +
  theme(legend.position = "none",
        aspect.ratio = 1, 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5))

#ABHD16A
ABHD16A.ppmi <- plotCounts(dds, gene = "ENSG00000204427", intgroup = "status_clin", returnData = T)
ABHD16A.ppmi <- ABHD16A.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(trans = "log2") +
  ggtitle("ABHD16A") +
  ylab("Log2 Normalised Counts") +
  xlab("Disease State") +
  graph_theme +
  theme(legend.position = "none",
        aspect.ratio = 1, 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5))

# APOBEC3A
APOBEC3A.ppmi <- plotCounts(dds, gene = "ENSG00000128383", intgroup = "status_clin", returnData = T)
APOBEC3A.ppmi <- APOBEC3A.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(trans = "log2") +
  ggtitle("APOBEC3A") +
  ylab("Log2 Normalised Counts") +
  xlab("Disease State") +
  graph_theme +
  theme(legend.position = "none",
        aspect.ratio = 1, 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5))

#PPP1R3B
PPP1R3B.ppmi <- plotCounts(dds, gene = "ENSG00000173281", intgroup = "status_clin", returnData = T)
PPP1R3B.ppmi <- PPP1R3B.ppmi %>% 
  mutate(status_clin = fct_relevel(status_clin, "PD", "HC")) %>% 
  ggplot(aes(x = status_clin, y = count, fill = status_clin)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(trans = "log2") +
  ggtitle("PPP1R3B") +
  ylab("Log2 Normalised Counts") +
  graph_theme +
  theme(legend.position = "none", 
        aspect.ratio = 1, 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5)
        )

cowplot_icicle_in_ppmi <- plot_grid(ABHD16A.ppmi,
          AGER.ppmi,
          APOBEC3A.ppmi,
          C6orf47.ppmi,
          CLK2.ppmi,
          DEFA1B.ppmi,
          DHX16.ppmi,
          GATAD2B.ppmi,
          HLA_DMB.ppmi,
          HLA_E.ppmi,
          MUC5B.ppmi,
          PPP1R3B.ppmi,
          nrow = 4,
          ncol = 3,
          align = "hv")

cowplot_icicle_in_ppmi
ggsave("figures/cowplot of icicle genes in ppmi.pdf")