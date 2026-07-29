setwd('/Users/juliasealock/Desktop/schema/add_dragen_data_nov25/manuscript/figures/figure2')

## A. enrichment in asd, dd, bip

plot_dat = read.csv('enrichment_by_pheno_group_and_anno_meta_updated.csv')
plot_dat = subset(plot_dat, pop == 'Meta-analysis')
plot_dat$Phenotype = ifelse(plot_dat$Phenotype == 'DD/ID', 'DD/ID\n(261 genes)', 'ASD\n(257 genes)')
plot_dat$Phenotype = factor(plot_dat$Phenotype, levels=c('DD/ID', 'DD/ID\n(261 genes)','ASD', 'ASD\n(257 genes)'))#,'BIP\n(27 genes)'))
plot_dat$x = ' '


plot_dat$anno = gsub('pLoF','PTV',plot_dat$anno)
plot_dat$anno = factor(plot_dat$anno, levels=c('PTV', 'PTV + Missense Rank 93%', 'Missense Rank 93%', 'Synonymous'))

## All mac plot

cols = c('PTV' = '#d90429', 'PTV + Missense Rank 93%' = "#fc9f5b", 'Missense Rank 93%' = "#129490", 'Synonymous' = '#9fb1bc') # #629677


plot_a = ggplot(data=plot_dat, aes(x=x, y=OR, ymin=Lower.CI, ymax=Upper.CI, fill=anno, colour=anno)) +
    geom_pointrange(position=position_dodge(width=2), size=1, linewidth=1.1) + 
    geom_hline(yintercept=1, lty=2) + 
    ggtitle(" ") +
    theme_bw() +
    xlab("") +
    theme(axis.text.x=element_text(size=12), axis.title.x=element_text(size=12)) +
    theme(axis.text.y=element_text(size=12), axis.title.y=element_text(size=12)) +
    theme(legend.title = element_text(size = 12), legend.text = element_text(size=12)) +
    theme(strip.text.x = element_text(size = 12, face='bold')) + 
    theme(axis.ticks.x = element_blank()) +
    ylab("Odds Ratio") +
    scale_colour_manual(values = cols) +
    scale_y_continuous(limits = c(1, 1.85), n.breaks=10) +
    facet_wrap(~Phenotype, nrow=1) +
    guides(colour=guide_legend(title="Annotation"), fill='none') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    	    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) 
  

## B - upset plot

df = read.csv('upset_plot_data.csv')
df = filtered_data2


sczp = read.csv('case_control_minp_ptv_mismean_rank93_2026-01-07.csv')[c(1,2)]
sczp = subset(sczp, cc_p_value <= fdr)
colnames(sczp)[2] = 'pvalue_out'
sczp = sczp[,c('gene_symbol','pvalue_out')]
df = merge(df, sczp, by.x='Gene', by.y='gene_symbol')
df = df[order(df$pvalue_out),]

df <- df %>%
  group_by(Gene) %>%
  mutate(
    # collapse the diseases with count > 0 into a single string per gene
    overlap_group = paste(sort(Disease[Presence > 0]), collapse = "&")
  ) %>%
  ungroup()

overlap_levels <- c(
  "ASD&BIP&DD&SCZ",
  "ASD&DD&SCZ",
  "ASD&SCZ",
  "BIP&SCZ",
  "DD&SCZ",
  "SCZ"              # SCZ only at the bottom
)


df <- df %>%
  mutate(
    overlap_group = factor(overlap_group, levels = overlap_levels)
  )

df <- df %>%
  arrange(overlap_group, pvalue_out) %>%
  mutate(
    Gene    = factor(Gene, levels = unique(Gene)),
    Disease = factor(Disease, levels = c("SCZ", "DD", "ASD", "BIP"))
  )

plot_b = ggplot(df, aes(y = Disease, x = Gene, group = Gene)) +
  geom_point(size = 4.5, aes(color = color), alpha = 1) +  # Move alpha outside of aes
  geom_line(linewidth = 0.75, aes(color = color)) +
  scale_color_identity() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.6, hjust = 1.25, size = 12),
        axis.text.y = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  guides(color = 'none') +
  xlab('') + 
  ylab('')



### C - heatmap of ORs between scz/bip/dd/asd
bonf = 1.31e-6
fdr = 6.84e-5

colnames(scz)[2] = 'pvalue_out'


scz = scz[,c('gene_symbol','pvalue_out','OR')]
bip = bip[,c('gene_symbol','OR','RES.p_value')]
asd = asd[,c('gene','or.odds.ratio','p')]

scz = subset(scz, pvalue_out <= fdr)

colnames(scz)[c(2:3)] = c('scz_pvalue','scz_or')
colnames(bip)[c(2:3)] = c('bip_or','bip_pvalue')
colnames(asd)[c(2:3)] = c('asd_or','asd_pvalue')

scz_bip = merge(scz, bip, by='gene_symbol')
dat = merge(scz_bip, asd, by.x='gene_symbol', by.y='gene')

dat = dat[order(dat$scz_pvalue),]



# Reshape to long format
library(tidyr)
library(dplyr)

df_long <- dat %>%
  pivot_longer(
    cols = ends_with("_or"),
    names_to = "condition",
    values_to = "OR"
  ) %>%
  mutate(condition = gsub("_or", "", condition))
df_long <- df_long %>%
  mutate(OR_plot = ifelse(is.infinite(OR),
                          max(OR[is.finite(OR)], na.rm = TRUE) * 1.1,
                          OR))
df_long$condition = toupper(df_long$condition)
df_long$condition = factor(df_long$condition, levels=c('SCZ','ASD','BIP'))

df_long$gene_symbol2 = factor(df_long$gene_symbol, levels=c("KDM5B", "SETD1A", "SCN2A", "TRIO", "ZMYM2", "SCAF1", "PPP3CA", "JARID2",
                                                            "PTK2", "RFX3", "NRXN1", "LRRC4", "HERC1", "RB1CC1", "SP4", "CUL1",
                                                            "ATP9A", "AKAP11", "PHIP", "ZMYND11", "XPO7", "FYN", "HDAC9", "STAG1",
                                                            "CHRM4", "PPFIA4", "ADGRL1", "ATP1A1", "ADAM23", "CRAMP1", "DLGAP3", "ZC3H12B",
                                                            "TNPO3", "DAGLA", "UBE2E3", "LRP1", "CACNA1B", "STK38L", "WDR78", "IQSEC1"))


plot_c = ggplot(df_long, aes(y = condition, x = gene_symbol2, fill = OR_plot)) +
  geom_tile(color = "white") +
    scale_fill_viridis_c(
    option = "D",
    trans = "log10",
    breaks = c(0.8, 1, 1.5, 2, 3, 5, 10, 20),
    labels = c("0.8", "1", '', "2", "3", "5", "10", "20")
  ) +
  theme_minimal() +
  labs(
       x = " ",
       y = "",
       fill = "Odds Ratio") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.25, size = 12),
        axis.text.y = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14)) + 
    # coord_equal() +
    guides(
    fill = guide_colorbar(
      barwidth = 12,
      barheight = 0.6,
      title.position = "top",
      title.hjust = 0.5
    )
  )


png("figure_2a.png", width = 10, height = 4, units = 'in', res = 1000)
print(plot_a)
dev.off()

png("figure_2b.png", width = 10, height = 4, units = 'in', res = 1000)
print(plot_b)
dev.off()

png("figure_2c.png", width = 10, height = 4, units = 'in', res = 1000)
print(plot_c)
dev.off()


