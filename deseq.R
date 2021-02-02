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

# filter samples
samples <- samples %>% 
  filter(clin_imag_agree == "Y" & ppmi_enroll != "Withdrew")

#vector of all filenames including path
files <- file.path("~/BL",samples$filename_long)
# rename files with patno
names(files) <- pull(samples, PATNO)

# creating tx2gene
txdb <- makeTxDbFromGFF(file="~/OneDrive - Newcastle University/PhD/PPMI/gencode.v19.annotation.gtf")
saveDb(x=txdb, file = "gencode.v19.annotation.TxDb")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME") 

# tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)

# create metadata file and check both files match
metadata <- data.frame(samples, row.names = colnames(txi$counts))
all(colnames(txi$counts) %in% rownames(metadata))
all(colnames(txi$counts) == rownames(metadata))

# deseq
dds <- DESeqDataSetFromTximport(txi,
                                colData = metadata,
                                design = ~ status_clin)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
results <- results(dds)
