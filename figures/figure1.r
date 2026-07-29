## figure 1
# 1a - Manhattan plot
# 1b - QQ plot
# 1c - sign test between schema 1 and 2


library(ggplot2)
library(ggrepel)
library(qqman)
library(tidyverse)

# A
bonf = 1.31e-6
fdr = 6.84e-05


pli = read.delim("gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
pc = subset(pli, gene_type=="protein_coding")
loc = pc[c(1,75:77)]
loc$chromosome = ifelse(loc$chromosome=="X",23, ifelse(loc$chromosome=="Y",24, loc$chromosome))
loc$chromosome = as.numeric(loc$chromosome)
loc$start_position = as.numeric(loc$start_position)
colnames(loc)[2:3] = c('CHR','BP')

gene = c('ADGRL1','PCNX3','CRAMP1', 'LINC01999', 'DOP1B',"IQSEC1")
CHR = c(19,11,16,17,21,3)
BP = c(14147743,65615776,1612337,60564547,36156782,12897043)
end_position = c(14206187,65637439,1677908,60586623,36294274,13283281)
test = data.frame(gene, CHR, BP, end_position)

loc = rbind(loc, test)



man_plot_dat = function(out=out){
  colnames(out)[1] = 'gene'
  out2 = merge(out, loc, by="gene")
  out2 = out2[!duplicated(out2$gene),]
  out2$fdr_sig = ifelse(out2$pvalue_out <= fdr, 'true','false')
  fdr_sig_genes = subset(out2, fdr_sig=='true')$gene
  n_fdr = length(fdr_sig_genes)
  print(paste0('fdr n:', n_fdr))
  n_bonf = nrow(subset(out2, pvalue_out < bonf))
  print(paste0('bonf n:', n_bonf))
  gwasResults = out2
  don <- gwasResults %>% 
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    mutate( is_annotate=ifelse(gene %in% fdr_sig_genes, "yes", "no")) #%>%
  return(don)
}

out = read.csv('case_control_minp_ptv_mismean_rank93_2026-01-07.csv')
colnames(out)[2] = 'pvalue_out'
don = man_plot_dat(out = out)

don$CHR = ifelse(don$CHR==23,'X',don$CHR)
don = subset(don, CHR != '24')

chr1 = c(1,5,9,13,17,21)
chr2 = c(2,6,10,14,18,22)
chr3 = c(3,7,11,15,19,'X')
chr4 = c(4,8,12,16,20)


don$colors = ifelse(don$pvalue_out <= bonf, "darkgoldenrod2",
                  ifelse(don$pvalue_out <= fdr, "darkgoldenrod2", 
                    ifelse(don$pvalue_out >= fdr & don$CHR %in% chr1, '#2596be', 
                      ifelse(don$pvalue_out >= fdr & don$CHR %in% chr2, '#3ba1c5',
                        ifelse(don$pvalue_out >= fdr & don$CHR %in% chr3, '#51abcb',
                          ifelse(don$pvalue_out >= fdr & don$CHR %in% chr4, '#66b6d2',
                        'black'))))))



axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Adjust point size based on pvalue_out relative to bonf
don$point_size = ifelse(don$pvalue_out <= bonf, 8, 
                  ifelse(don$pvalue_out <= fdr, 7, 1.3))
don$alpha_var = ifelse(don$pvalue_out <= fdr, 1, 0.9)

# Adjust label boldness based on pvalue_out relative to bonf
don$label_fontface = ifelse(don$pvalue_out < bonf, "bold", "plain")

nudge_df <- data.frame(
  gene = c("SCN2A", "UBE2E3", "STAG1", "KDM5B", "WDR78", "DLGAP3", "ATP1A1", "PPFIA4", "NRXN1", "TRIO", "PHIP"),
  nudge_x = c(-0.5e8,  0.3e8,   0.5e8,   -0.3e8,  -0.5e8,  -0.2e8,   -0.4e8,   -0.6e8,   0.4e8, 0.2e8, 0.2e8),
  nudge_y = c(0.3,    -0.8,     0.3,     -0.5,    -0.8,    -0.6,      0.3,     0.5,      0.5, -0.5, -0.5)
)

don <- merge(don, nudge_df, by = "gene", all.x = TRUE)
don$nudge_x[is.na(don$nudge_x)] <- 0
don$nudge_y[is.na(don$nudge_y)] <- 0

# subsets used in both geom_label_repel calls
bonf_dat <- subset(don, is_annotate=="yes" & pvalue_out <= bonf)
fdr_dat  <- subset(don, is_annotate=="yes" & pvalue_out > bonf & pvalue_out <= fdr)


plot_a = ggplot(don, aes(x=BPcum, y=-log10(pvalue_out))) + 
    geom_point( aes(color=colors, size=point_size)) +
    scale_color_identity() +
      geom_hline(yintercept = -log10(bonf), linetype=2) + 
      geom_hline(yintercept = -log10(fdr), linetype=2, color="darkgray") + 
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center, expand=c(0,0)) +
      scale_y_continuous(expand = c(0, 1) ) +     # remove space between plot area and x axis
geom_text_repel(data = bonf_dat, aes(label = gene), nudge_x = bonf_dat$nudge_x, nudge_y = bonf_dat$nudge_y, fontface = "bold",
    size = 6.5, box.padding = 0.5, point.padding = 0.3, force = 2, force_pull = 0.5, max.overlaps = Inf, 
    segment.color = "grey50", segment.size = 0.3, min.segment.length = 0.8) +  
# FDR only hits - plain
geom_text_repel(data = fdr_dat, aes(label = gene), nudge_x = fdr_dat$nudge_x, nudge_y = fdr_dat$nudge_y, fontface = "plain",
    size = 6, box.padding = 0.5, point.padding = 0.3, force = 2, force_pull = 0.5, max.overlaps = Inf,
    segment.color = "grey50", segment.size = 0.3, min.segment.length = 0.8) +
      theme_bw() +
      theme( legend.position="top", legend.text = element_text(size = 16), legend.title = element_text(size = 16),  
        panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text = element_text(size=14), axis.title = element_text(size=14) ) +
      guides(color=guide_legend(override.aes = list(size=5)), size='none') +
      xlab("") + 
      ylab("-log10(P)") 

png('figure_1a_2.png', width=20, height=10, unit='in', res=1000)
print(plot_a)
dev.off()


## B - qqplot
out = read.csv('case_control_minp_ptv_mismean_rank93_2026-01-07.csv')
get_qqp_data = function(input=input, group=group){
    colnames(input)[2] = 'pvalue_out'
    input = input[!duplicated(input$gene_symbol),]
    input$pvector = input$pvalue_out
    input <- subset(input, !is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector < 1 & pvector > 0)
    # Observed and expected
    input$o = -log10(input$pvector)
    input = input[order(-input$o),]
    input$e = -log10(ppoints(length(input$pvector) ))
    # input$group = paste0(group)
    input$group = paste0(group)
    return(input)
}

qqp = get_qqp_data(input=out, group='SCHEMA2')

qqp$is_annotate = ifelse(qqp$pvalue_out <= bonf, 'yes','no')

don_colors = don[,c('gene','colors')]
qqp = merge(qqp, don_colors, by.x='gene_symbol', by.y='gene')

qqp$point_size = ifelse(qqp$is_annotate == 'yes', 'big', 'small')

    
plot_b = ggplot(qqp, aes(x = e, y = o, color=colors)) +
    geom_hline(yintercept=-log10(bonf), linetype="dashed", color = "darkgray") +
    geom_hline(yintercept=-log10(fdr), linetype="dashed", color = "gray") +
    geom_point(aes(size = point_size)) +
    scale_size_manual(values = c('small' = 3, 'big' = 6)) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
       x = "Expected -log10(p-value)",
       y = "Observed -log10(p-value)") +
    scale_color_identity() +
    theme_classic() + 
    theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
    guides(color=guide_legend(title="Annotation", override.aes=list(size=4))) +
    theme(axis.text = element_text(size=14), axis.title=element_text(size=14)) +
    guides(color='none', size='none') 


### OR x freq plot
ptv = read.delim('2026_01_07_ALL_SCHEMA_PTV_MAC15_GQ25_ptv_AC15_CMH_pc_matched.tsv')[c(1,18,19)]
out1 = subset(out, cc_p_value <= fdr)
ptv = ptv[(ptv$gene_symbol %in% out1$gene_symbol),]
write.csv(ptv, 'schema2_fdr_genes_OR_gnomad_freq.csv', row.names=F, quote=F)

## plot
dat = read.csv('schema2_fdr_genes_OR_gnomad_freq.csv')

dat$gnomad_freq = dat$gnomad_obs/(807162*2)

dat$OR_plot <- ifelse(is.infinite(dat$OR), 100, dat$OR)

dat$gnomad_freq = as.numeric(dat$gnomad_freq)

## add gwas loci
gwas = read.delim('PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv')
gwas = subset(gwas, PVAL < 5e-8)
gwas$OR = exp(gwas$BETA)
gwas$Cite = 'Trubetskoy et al'
gwas$old = 0
gwas$gnomad_obs = 0
gwas$gnomad_freq = gwas$FCON 
gwas$Class = 'GWAS'
gwas$Lower_OR = gwas$OR
gwas$Upper_OR = gwas$OR
gwas$Lower_OR <- exp(gwas$BETA - 1.96 * gwas$SE)
gwas$Upper_OR <- exp(gwas$BETA + 1.96 * gwas$SE)

gwas = gwas[,c('ID','PVAL', 'OR', 'gnomad_obs', 'gnomad_freq')]

colnames(gwas)[1] = 'gene_symbol'
colnames(gwas)[2] = 'SCHEMA2_pvalue'
gwas$OR_plot = gwas$OR
gwas = subset(gwas, gnomad_freq <= 0.5)
gwas$Class = 'GWAS'
dat$Class = 'PTV'

colnames(dat)[2] = 'SCHEMA2_pvalue'

dat = rbind(dat, gwas)
dat$label = ifelse(dat$Class=='GWAS','false','true')

dat$Class = factor(dat$Class, levels=c('PTV','GWAS'))
colors = c('PTV' = '#d90429','GWAS' = '#4e6766')


custom_labels <- c(
    "0.5", "0.1", "0.05", "0.01", "0.001", "0.0001",
    expression("1" %*% 10^-5),
    expression("1" %*% 10^-6)
)

plot_c = ggplot(dat, aes(x = gnomad_freq, y = OR_plot, color = Class, label = gene_symbol)) +
  geom_jitter(aes(size = OR_plot)) +
  scale_x_log10(
    breaks = c(0.5, 0.1, 0.05, 0.01, 0.001, 0.0001, 1e-5, 1e-6),
    labels = custom_labels
  ) +
  scale_y_log10() +
  labs(
    x = "Minor Allele Frequency in General Population",
    y = "Odds ratio",
    size = "Odds Ratio"
  ) +
  theme_minimal() +
  geom_text_repel(
    data = subset(dat, label == 'true'),
    aes(label = gene_symbol),
    size = 5,
    max.overlaps = Inf,
    color = 'black'
  ) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  ) +
  scale_color_manual(
    values = colors,
    name = 'Variant Type',
    guide = guide_legend(override.aes = list(shape = 16, size = 4))
  ) +
  scale_size_continuous(range = c(3, 8)) +
  guides(size = 'none')



# figure D - forest plot
results = read.delim('schema2_case_control_results_05-19-2026.tsv')
results$minp_group <- ifelse(
  results$PTV.Pvalue == results$PTV...Missense.Pvalue, 'Both',
  ifelse(results$PTV.Pvalue < results$PTV...Missense.Pvalue, 'PTV', 'PTV + Missense Rank 93%')
)
results = results[,c('Gene','minp_group','PTV...Missense.Pvalue','PTV.Pvalue')]
results = results[(results$Gene %in% s2_genes$gene),]

ptv_label_genes = subset(results, minp_group == 'PTV' | minp_group == 'Both')$Gene
ptv_mis_label_genes = subset(results, minp_group == 'PTV + Missense Rank 93%' | minp_group == 'Both')$Gene

s2_ptv = read.delim('2026_06_12_ALL_SCHEMA_PTV_MAC15_GQ25_ptv_AC15_CMH_pc_matched.tsv')
s2_ptv = s2_ptv[,c('gene_symbol','OR','CI_lower','CI_upper')]
colnames(s2_ptv) = c('gene', 'OR_PTV','PTV_lower','PTV_upper')
s2_ptv$group = 'PTV'
s2_ptv$label = ifelse(s2_ptv$gene %in% ptv_label_genes, '*', ' ')

s2_ptv_mis = read.delim('2026_06_12_ALL_SCHEMA_PTV_MIS_MAC15_GQ25_ptv-mis_mean_rank93_AC15_CMH_pc_matched.tsv')
s2_ptv_mis = s2_ptv_mis[,c('gene_symbol','OR','CI_lower','CI_upper')]
colnames(s2_ptv_mis) = c('gene', 'OR_PTV','PTV_lower','PTV_upper')
s2_ptv_mis$group = 'PTV + Missense Rank 93%'
s2_ptv_mis$label = ifelse(s2_ptv_mis$gene %in% ptv_mis_label_genes, '*', ' ')


schema2_genes = c('SETD1A', 'ZMYM2', 'HERC1', 'RB1CC1', 'SCAF1', 'XPO7', 'SP4', 'FYN', 'PPP3CA', 'CUL1', 'HDAC9', 'JARID2', 'ATP9A', 'PTK2', 'STAG1', 'SCN2A')

dat_ptv = rbind(s2_ptv, s2_ptv_mis)

s2_genes = dat_ptv[(dat_ptv$gene %in% schema2_genes),]

library(dplyr)
library(ggplot2)

# Choose a display cap — log(100) is reasonable given your data
cap <- 100

s2_genes <- s2_genes %>%
  mutate(
    OR_plot    = ifelse(is.infinite(OR_PTV),  cap, pmin(OR_PTV, cap)),
    lower_plot = ifelse(is.infinite(PTV_lower), NA,  PTV_lower),
    upper_plot = ifelse(is.infinite(PTV_upper), cap, pmin(PTV_upper, cap)),
    capped     = OR_PTV > cap | PTV_upper > cap | is.infinite(OR_PTV)  # flag for annotation
  )


  # Manually assign dodge offsets
dodge_offset <- 0.7 / 2 * 0.5  # half the dodge width

capped_data <- filter(s2_genes, capped) %>%
  mutate(x_dodge = as.numeric(gene) + ifelse(group == "SCHEMA 1", dodge_offset, -dodge_offset))



s2_genes$gene = factor(s2_genes$gene, levels=rev(c('SETD1A', 'ZMYM2', 'HERC1', 'RB1CC1', 'SCAF1', 'XPO7', 'SP4', 'FYN', 'PPP3CA', 'CUL1', 'HDAC9', 'JARID2', 'ATP9A', 'PTK2', 'STAG1', 'SCN2A')))
s2_genes$group = factor(s2_genes$group, levels=c("PTV + Missense Rank 93%", "PTV"))

cols = c('PTV' = '#d90429', 'PTV + Missense Rank 93%' = "#fc9f5b") 

plot_d = ggplot(s2_genes, aes(x = gene, y = log(OR_plot), ymin = log(lower_plot), ymax = log(upper_plot),
         color = group, fill = group)) +
geom_pointrange(size = 0.8, linewidth = 1, position = position_dodge(width = 0.7)) +
geom_segment(
    data = capped_data,
    aes(x = x_dodge, xend = x_dodge, y = log(OR_plot), yend = log(cap) + 0.2,
        color = group),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth= 1.2,
    show.legend = FALSE
  ) +
  #  scale_color_npg() + 
  scale_colour_manual(values = cols) +
  #scale_shape_manual(values = c("PTV" = 16, "PTV + Missense" = 17)) +
  guides(
    shape = guide_legend(reverse = TRUE),  # flips legend back to SCHEMA 1 first
    color = guide_legend(reverse = TRUE),
    fill  = "none"
  ) +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = NULL, y = "PTV OR", color = NULL, shape = NULL) +
  scale_y_continuous(
    limits = c(log(0.1), log(cap) + 0.3),
    breaks = log(c(0.1, 0.5, 1, 2, 5, 10, 25, 50, 100)),
    labels =     c("0.1", "0.5", "1", "2", "5", "10", "25", "50", "100")
  ) +
  theme_bw() +  
theme(
    legend.position   = "top", 
    legend.box        = "vertical", 
    legend.spacing.y  = unit(0, "pt"),
    legend.key.width  = unit(1.5, "cm"),
    axis.text         = element_text(size = 14), 
    axis.title        = element_text(size = 14), 
    legend.title      = element_text(size = 14, face = "bold"), 
    legend.text       = element_text(size = 14), 
    strip.background  = element_rect(fill = "white"), 
    strip.text        = element_text(face = "bold", size = 14, hjust = 0),
    plot.title = element_text(face = "bold")
  ) +
  ylab('Odds Ratio')





png('figure_1a_2.png', width=20, height=10, unit='in', res=1000)
print(plot_a)
dev.off()

png('figure_1b.png', width=10, height=10, unit='in', res=1000)
print(plot_b)
dev.off()

png('figure_1c.png', width=20, height=10, unit='in', res=1000)
print(plot_c)
dev.off()


png('figure_1d_ptv_mis_forest.png', width=10, height=10, unit='in', res=1000)
print(plot_d)
dev.off()

