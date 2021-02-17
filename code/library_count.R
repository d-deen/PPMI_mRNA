dds <- DESeqDataSetFromTximport(txi,
                                colData = metadata,
                                design = ~ status_clin)

dds <- DESeq(dds)

# library sizes
library_sizes <- colSums(counts(dds))

library_sizes <- as.data.frame(library_sizes)
library_sizes <- library_sizes %>% 
  rownames_to_column(var = "patno")

ggplot(library_sizes, aes(x = patno, y = library_sizes)) +
  geom_col(colour = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Sample") +
  ylab("Library Size") +
  coord_flip() +
  graph_theme +
  theme(axis.text.x = element_blank())


normalized_counts_unfiltered <- counts(dds, normalized = TRUE)
write.csv(normalized_counts_unfiltered, "output/normalized_counts_unfiltered.csv")
view(normalized_counts_unfiltered)

# avg lib size vs gene counts 
ncu_mean_gene <- as.data.frame(colMeans(normalized_counts_unfiltered))
ncu_mean_gene <- ncu_mean_gene %>% 
  rownames_to_column (var = "patno") %>% 
  rename(means = `colMeans(normalized_counts_unfiltered)`)

ncu_observations <- as.data.frame(colSums(normalized_counts_unfiltered != 0))
ncu_observations <- ncu_observations %>% 
  rownames_to_column (var = "patno") %>% 
  rename(observations = `colSums(normalized_counts_unfiltered != 0)`)

ncu_sum_gene <- as.data.frame(colSums(normalized_counts_unfiltered))
ncu_sum_gene <- ncu_sum_gene%>% 
  rownames_to_column (var = "patno") %>% 
  rename(lib_size = `colSums(normalized_counts_unfiltered)`)



df <- merge(ncu_mean_gene, ncu_observations, by = "patno")
df <- merge(df, ncu_sum_gene, by = "patno")
a <- ggplot(df, aes(x = observations, y = means)) +
  geom_point(shape = 1) +
  geom_smooth(method=lm) +
  xlab("Number of genes detected") +
  ylab("Average gene count") +
  graph_theme

b <- ggplot(df, aes(x = observations, y = lib_size)) +
  geom_point(shape = 1) +
  geom_smooth(method=lm) +
  stat_cor(p.accuracy = 0.001) +
  xlab("Number of genes detected") +
  ylab("Library size") +
  graph_theme

c <- ggplot(df, aes(x = means, y = lib_size)) +
  geom_point(shape = 1) +
  xlab("Average gene count") +
  ylab("Library size") +
  stat_cor(p.accuracy = 0.001) +
  graph_theme
c

p <- ggarrange(c, b)
annotate_figure(p, top = "Normalised")