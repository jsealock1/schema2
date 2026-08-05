

# a - GTEx by tissue
# b - psychencode single cell clusters
# c - brainspan heatmap


# A Gtex tissue expression
dat = read.delim(FUMA)
dat = subset(dat, Category=='DEG.up')
dat$color = ifelse(dat$adjP < 0.05, 'darkred','cornflowerblue')

dat$GeneSet = gsub('_',' ',dat$GeneSet)

plot_a = ggplot(dat, aes(x = reorder(GeneSet, p), y = -log10(p), fill=color)) +
    geom_col() +
    theme_classic() +
    labs(x = NULL, y = '-log10(p)') +  
    scale_fill_identity() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size=14), 
          axis.text.y = element_text(size=16),
          axis.title.y = element_text(size=16), 
        plot.margin = unit(c(1, 1, 1, 4), "cm"))

# B single cell clusters

library(ggplot2)
library(ggpubr)
library(FSA)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)


df <- do.call(rbind, lapply(files, read.table, sep = "\t", header = TRUE))

df$Cell.type <- factor(df$Cell.type, levels = unique(df$Cell.type))
df$Gene <- factor(df$Gene, levels = unique(df$Gene))

# Step 1: Average expression across datasets
df_avg <- df %>%
  group_by(Gene, Cell.type) %>%
  summarise(AvgExpr = mean(Average.Expression, na.rm = TRUE), .groups = "drop")

# Step 2: Reshape to matrix form
expr_matrix_raw <- df_avg %>%
  pivot_wider(names_from = Cell.type, values_from = AvgExpr) %>%
  column_to_rownames("Gene")

# Step 3: Normalize expression per gene (z-score)
expr_matrix_z <- t(scale(t(expr_matrix_raw)))  # row-wise z-score normalization

# Optional: Check for any NA rows from all-zero genes and remove
expr_matrix_z <- expr_matrix_z[complete.cases(expr_matrix_z), ]

# Step 4: Elbow method to choose optimal k
set.seed(1)
wss <- sapply(1:10, function(k) {
  kmeans(expr_matrix_z, centers = k, nstart = 25)$tot.withinss
})


plot(1:10, wss, type = "b", pch = 19,
     xlab = "Number of Clusters (k)",
     ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method for K-means")

# Step 5: Final clustering with chosen k
set.seed(42)
k_clust <- kmeans(expr_matrix_z, centers = 4)  # adjust k as needed

# Step 6: Add cluster labels
gene_clusters <- data.frame(Gene = rownames(expr_matrix_z), Cluster = k_clust$cluster)

# Optional: Join back with original data
df_clustered <- df %>%
  left_join(gene_clusters, by = "Gene")


# Convert z-score matrix to long format and join clusters
heatmap_df <- expr_matrix_z %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Cell.type", values_to = "Norm.Expression") %>%
  left_join(gene_clusters, by = "Gene")

# Order genes within cluster for prettier plotting
heatmap_df <- heatmap_df %>%
  group_by(Cluster) %>%
  arrange(Cluster, Gene) %>%
  mutate(Gene = factor(Gene, levels = unique(Gene))) %>%
  ungroup()

pvalues = read.csv('case_control_minp_ptv_mismean_rank93_2026-01-07.csv')
colnames(pvalues)[2] = 'pvalue_out'
pvalues = pvalues[,c('gene_symbol','pvalue_out')]
bonf = 1.31e-6
fdr = 6.84e-05

heatmap_df = merge(heatmap_df, pvalues, by.x='Gene', by.y='gene_symbol')
heatmap_df = heatmap_df[order(heatmap_df$pvalue_out),]
heatmap_df$Gene <- factor(heatmap_df$Gene, levels = unique(heatmap_df$Gene))

heatmap_df$Cluster = paste0('Cluster ', heatmap_df$Cluster)


library(ggplot2)
library(patchwork)


cluster_sizes <- heatmap_df %>%
  distinct(Cluster, Gene) %>%
  count(Cluster, name = "n_genes") %>%
  arrange(desc(n_genes)) %>%       # largest → smallest
  mutate(new_label = paste0("Cluster ", row_number()))

# Mapping: old → new labels
cluster_map <- setNames(cluster_sizes$new_label, cluster_sizes$Cluster)

# --- 2. Apply new cluster labels to the dataframe ---
heatmap_df$Cluster <- cluster_map[heatmap_df$Cluster]

# --- 3. Make clusters an ordered factor: Cluster 1 → Cluster 4 ---
heatmap_df$Cluster <- factor(
  heatmap_df$Cluster,
  levels = paste0("Cluster ", seq_len(nrow(cluster_sizes)))
)

clusters <- levels(heatmap_df$Cluster)



plots <- lapply(seq_along(clusters), function(i) {
  clust <- clusters[i]
  df_sub <- subset(heatmap_df, Cluster == clust)

  p <- ggplot(df_sub, aes(x = Gene, y = Cell.type, fill = Norm.Expression)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      colors = rev(c("#8B0000", "#FF6347", "#FFD700", "#ADD8E6", "#00008B")),
      name = "Z-score",
      limits = c(-3, 3)
    ) +
    coord_fixed() +
    ggtitle(clust) +
    theme_bw() +
    labs(x = "", y = "") +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      legend.position = "none"      # ← REMOVE LEGEND FROM ALL PLOTS
    )

  # Remove y-axis & legend from all but first plot
  if (i != 1) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
      )
  }

  return(p)
})


# Stitch plots and collect shared legend
plot_b = wrap_plots(plots, nrow = 1) + plot_layout(guides = "collect")




## C - brainspan heatmap 

#### heatmap 
dat2 = read.delim(BRAINSPAN, header=T, sep="\t")

# divide into larger age groups 
first_trimester = c('8 pcw', '9 pcw', '12 pcw')
second_trimester = c('13 pcw', '16 pcw', '17 pcw', '19 pcw','21 pcw', '24 pcw', '25 pcw', '26 pcw')
third_trimester = c('35 pcw', '37 pcw')

infancy = c('4 mos', '10 mos')
childhood = c('1 yrs','2 yrs', '3 yrs', '4 yrs', '8 yrs', '11 yrs')
adolescence = c('13 yrs', '15 yrs', '18 yrs', '19 yrs')
adulthood = c('21 yrs', '23 yrs', '30 yrs', '36 yrs', '37 yrs', '40 yrs')

dat2$age_group = ifelse(dat2$age %in% first_trimester, 'first trimester',
                ifelse(dat2$age %in% second_trimester, 'second trimester',
                ifelse(dat2$age %in% third_trimester, 'third trimester',
                  ifelse(dat2$age %in% infancy, 'infancy', 
                    ifelse(dat2$age %in% childhood, 'childhood', 
                      ifelse(dat2$age %in% adolescence, 'adolescence', 
                        ifelse(dat2$age %in% adulthood, 'adulthood', 'NA')))))))

dat2$age_group = factor(dat2$age_group, levels=c('first trimester', 'second trimester','third trimester','infancy','childhood','adolescence','adulthood'))


pvalues = read.csv(OUT)
colnames(pvalues)[2] = 'pvalue_out'
fdr_genes = subset(pvalues, pvalue_out <= fdr)

genes = c("SETD1A",  "ZMYM2",   "HERC1",   "RB1CC1",  "SP4",     "SCAF1" ,  "XPO7",    "FYN",     "PPP3CA",  "CUL1",
        "HDAC9",   "JARID2",  "ATP9A",   "PTK2",    "STAG1",   "SCN2A",   "CHRM4",   "PHIP",    "PPFIA4",  "TRIO", 
        "ADGRL1",  "ATP1A1",  "ADAM23",  "CRAMP1",  "ZMYND11", "DLGAP3",  "ZC3H12B", "TNPO3",   "AKAP11",  "RFX3",   
        "DAGLA",   "UBE2E3",  "LRP1",    "NRXN1",   "CACNA1B", "STK38L",  "WDR78",   "KDM5B",   "IQSEC1",  "LRRC4"
)

fdr_dat = dat2[(dat2$gene_symbol %in% genes),]
## ADGRL1 and CRAMP1 missing

## z-score scale within gene
library(dplyr)

fdr_dat <- fdr_dat %>%
  group_by(gene_symbol) %>%
  mutate(TPM_zscore = scale(tpm)) %>%
  ungroup()

fdr_dat$age = factor(fdr_dat$age, levels=c('8 pcw', '9 pcw', '12 pcw','13 pcw', '16 pcw', '17 pcw', '19 pcw','21 pcw', '24 pcw', '25 pcw', '26 pcw',
                                          '35 pcw', '37 pcw', '4 mos', '10 mos', '1 yrs','2 yrs', '3 yrs', '4 yrs', '8 yrs', '11 yrs', 
                                          '13 yrs', '15 yrs', '18 yrs', '19 yrs', '21 yrs', '23 yrs', '30 yrs', '36 yrs', '37 yrs', '40 yrs'))
                                

fdr_dat$age_group = factor(fdr_dat$age_group, levels=rev(c('first trimester', 'second trimester','third trimester',
                                                          'infancy','childhood','adolescence','adulthood')))


heatmap_data <- fdr_dat %>%
  group_by(gene_symbol, age_group) %>%
  summarize(mean_TPM_z = mean(TPM_zscore, na.rm = TRUE), .groups = "drop")

heatmap_matrix <- heatmap_data %>%
  pivot_wider(names_from = age_group, values_from = mean_TPM_z) %>%
  column_to_rownames("gene_symbol")

row_order <- hclust(dist(heatmap_matrix))$order
gene_order <- rownames(heatmap_matrix)[row_order]

heatmap_data$gene_symbol <- factor(heatmap_data$gene_symbol, levels = gene_order)



early_ordered <- c("SETD1A", "ZMYM2", "HERC1", "RB1CC1",  "SP4", 
                  "XPO7", "FYN", "CUL1", "HDAC9","JARID2", "ATP9A",
                  "STAG1","PHIP","TRIO","ZMYND11", "TNPO3","KDM5B",
                    "RFX3","UBE2E3")

mid_ordered <- c(
  "SCAF1", "WDR78","LRP1","CACNA1B","NRXN1","STK38L")

late_ordered <- c("DLGAP3","PPP3CA","PTK2","SCN2A","CHRM4","PPFIA4",
                "ATP1A1","ADAM23","AKAP11","ZC3H12B","DAGLA",
                "IQSEC1", 'LRRC4')

all = c(early_ordered, mid_ordered, late_ordered)



heatmap_data$gene_symbol = factor(heatmap_data$gene_symbol, levels = rev(all))
heatmap_data$age_group <- factor(heatmap_data$age_group, levels = c(
  'first trimester', 'second trimester', 
  'third trimester', 'infancy', 'childhood', 'adolescence', 'adulthood'
))



plot_c = ggplot(heatmap_data, aes(x = age_group, y = gene_symbol, fill = mean_TPM_z)) +
  geom_tile() +
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      scale_fill_gradientn(
      colors = rev(c("#8B0000", "#FF6347", "#FFD700", "#ADD8E6", "#00008B")),
      name = "Z-score",
      limits = c(-1.05, 1.5)  # adjust as needed
    ) +
  labs(x = "Developmental Stage", y = " ", fill = "Mean TPM Z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=14),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16)) +
  coord_fixed()


png('figure_3c_2.png', width=20, height=8, unit='in', res=1000)
print(plot_c)
dev.off()
