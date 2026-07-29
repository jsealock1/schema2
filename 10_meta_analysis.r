
setwd('/Users/juliasealock/Desktop/schema/add_dragen_data_nov25/cmh/acat')

library(dplyr)
library(data.table)
date = Sys.Date()

library(data.table)

ptv = read.delim('2026_01_07_ALL_SCHEMA_PTV_MAC15_GQ25_ptv_AC15_CMH_pc_matched.tsv')[c(1,18,19)]
ptv_mis = read.delim('2026_01_07_ALL_SCHEMA_PTV_MIS_MAC15_GQ25_ptv-mis_mean_rank93_AC15_CMH_pc_matched.tsv')[c(1,18,19)]


ptv = subset(ptv, RES.p_value > 0)
ptv_mis = subset(ptv_mis, RES.p_value > 0)

## find minp of case-controls
ptv_mis$annotation = 'ptv_mis'
ptv$annotation = 'ptv'

cc = rbind(ptv_mis, ptv)
colnames(cc)[2] = 'p_value'
cc = subset(cc, p_value > 0)
cc_minp = cc %>% group_by(gene_symbol) %>% slice(which.min(p_value))
table(cc_minp$annotation)

colnames(cc_minp)[2] = 'cc_p_value'
colnames(cc_minp)[4] = 'cc_annotation'

write.csv(cc_minp, paste0('case_control_minp_ptv_mismean_rank93_',date,'.csv'), row.names=F, quote=F)


# find minp of de novo 
de_novo_v2_results = '/Users/juliasealock/Desktop/schema/schema1_de_novo/july_2024/schema1_de_novo_poisson_results_with_gnomad_v2_rates_with_misrank.tsv.bgz'

# find minp of de novo 
de_novo = read.delim(de_novo_v2_results, header=T)
de_novo1 = subset(de_novo, annotation=='pLoF' | annotation=='mis' | annotation=='pLoF_or_mis')
de_novo1 = de_novo1 %>% group_by(gene_symbol) %>% slice(which.min(p_value))

dat_ptv_mis = merge(de_novo1, cc_minp, by="gene_symbol")


## case-control + de novo stouffer
run_stouffer = function(dat = dat, de_novo = de_novo, case_control=case_control){
    p1 = setDT(dat)[,c('gene_symbol','cc_p_value','p_value')]
    colnames(p1)[3] = 'de_novo_pvalue'
    p1[p1 == 1] <- 0.99999

    w1 = 10
    w2 = 1
    w = c(w1, w2)
    p1$cc_p_value1 = qnorm(p1$cc_p_value/2) ## two sided with p/2
    p1$de_novo_pvalue1 = qnorm(p1$de_novo_pvalue/2)
    p1$stoufferP <- 2*pnorm(-abs(((p1$cc_p_value1*w1) + (p1$de_novo_pvalue1*w2))/sqrt(sum(w^2))))
  
    # add back in genes not tested in the meta 
    p2 = p1[,c('gene_symbol','stoufferP')]
    p2 = merge(p2, case_control, by='gene_symbol', all=TRUE)
    p3 = merge(p2, de_novo, by="gene_symbol", all=TRUE)

    p3$pvalue_out = ifelse(is.na(p3$stoufferP), p3$cc_p_value, p3$stoufferP)
    p3 = p3[!is.na(p3$gene_symbol),]
    p3 = p3[order(p3$pvalue_out),]
    p3 = subset(p3, pvalue_out>0)
    # colnames(p3)[4:9] = c('cc_OR','n_de_novo','de_novo_mu','de_novo_expected','de_novo_pvalue','de_novo_annotation')
    return(p3)
}

out_ptv_mis_dn = run_stouffer(dat=dat_ptv_mis, case_control=cc_minp, de_novo=de_novo1)
write.csv(out_ptv_mis_dn, paste0("minp_ptv_mis_misrank93_stouffer_meta_de_novo_v2_rates_mac15_results_cc_weight10_",date,".csv"), row.names=F, quote=F)






# setDT(ptv)
# setDT(ptv_mis)

# # Unweighted ACAT
# acat <- function(p, eps = 1e-15) {
#   p <- pmax(pmin(p, 1 - eps), eps)
#   if (!length(p)) return(NA_real_)
#   T <- mean(tan((0.5 - p) * pi))
#   0.5 - atan(T) / pi
# }

# res_acat <- merge(
#   ptv[, .(gene_symbol, p_ptv = RES.p_value)],
#   ptv_mis[, .(gene_symbol, p_mis = RES.p_value)],
#   by = "gene_symbol",
#   all = TRUE
# )[
#   , .(
#       p_ptv,
#       p_mis,
#       p_acat = {
#         p_vec <- c(p_ptv, p_mis)
#         p_vec <- p_vec[is.finite(p_vec) & !is.na(p_vec)]
#         acat(p_vec)
#       }
#     ),
#   by = gene_symbol
# ]
# head(res_acat[order(res_acat$p_acat),])


library(ggplot2)
library(ggrepel)
library(qqman)
# detach("package:plyr", unload=TRUE)



bonf = 1.31e-6
fdr = 6.84e-05








setwd('/Users/juliasealock/Desktop/schema/add_dragen_data_nov25/cmh/acat')

library(dplyr)
library(data.table)
library(ggplot2)
library(ggsci)

date = Sys.Date()

ptv = read.delim('2026_01_07_ALL_SCHEMA_PTV_MAC15_GQ25_ptv_AC15_CMH_pc_matched.tsv')[c(1,18)]
ptv_mis = read.delim('2026_01_07_ALL_SCHEMA_PTV_MIS_MAC15_GQ25_ptv-mis_mean_rank93_AC15_CMH_pc_matched.tsv')[c(1,18)]
minp = read.csv('case_control_minp_ptv_mismean_rank93_2026-01-07.csv')[c(1,2)]
minp_dn = read.csv('minp_ptv_mis_misrank93_stouffer_meta_de_novo_v2_rates_mac15_results_cc_weight10_2026-01-07.csv')[c(1,11)]

ptv = subset(ptv, RES.p_value > 0)
ptv_mis = subset(ptv_mis, RES.p_value > 0)

get_qqp_data = function(input=input, group=group){
    colnames(input)[2] = 'pvalue_out'
    input = input[!duplicated(input$gene_symbol),]
    input$pvector = input$pvalue_out
    input <- subset(input, !is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector < 1 & pvector > 0)
    input$o = -log10(input$pvector)
    input = input[order(-input$o),]
    input$e = -log10(ppoints(length(input$pvector) ))
    input$group = paste0(group)
    return(input)
}

ptv_qqp = get_qqp_data(input=ptv, group='PTV')
ptv_mis_qqp = get_qqp_data(input=ptv_mis, group='PTV + Missense')
minp_qqp = get_qqp_data(input=minp, group='Lowest Pvalue')
minp_dn_qqp = get_qqp_data(input=minp_dn, group='De Novo Meta-Analysis')


# qq_dat = rbind(ptv_qqp, ptv_mis_qqp, acat_qqp, minp_qqp)
qq_dat = rbind(ptv_qqp, ptv_mis_qqp, minp_qqp, minp_dn_qqp)


bonf = 1.31e-6
fdr = 6.84e-05


qq_dat$is_annotate = ifelse(qq_dat$pvalue_out <= bonf, 'yes', 'no')

qq_dat$group = factor(qq_dat$group, levels=c('PTV', 'PTV + Missense', 'Lowest Pvalue', 'De Novo Meta-Analysis'))

pdf('schema2_qqp_by_anno_wide_2026-01-07.pdf', width=18, height=6)
ggplot(qq_dat, aes(x = e, y = o, color=group)) +
    geom_hline(yintercept=-log10(bonf), linetype="dashed", color = "darkgray") +
    geom_hline(yintercept=-log10(fdr), linetype="dashed", color = "gray") +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
       x = "Expected -log10(p-value)",
       y = "Observed -log10(p-value)") +
    scale_colour_d3() +
    theme_bw() + 
    facet_wrap(~group, nrow=1) +
    theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
    guides(color=guide_legend(title="Annotation", override.aes=list(size=4))) + 
    geom_text_repel(data=subset(qq_dat, is_annotate=="yes"), aes(label=gene_symbol), position = position_nudge_repel(y=0.3, x=0.1), size=4, max.overlaps = Inf) +
    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA),
            strip.text = element_text(size=14, face='bold')) +
    guides(color='none')
dev.off()

png('schema2_qqp_by_anno2_2026-01-07.png', width=8, height=8, unit='in', res=300)
ggplot(qq_dat, aes(x = e, y = o, color=group)) +
    geom_hline(yintercept=-log10(bonf), linetype="dashed", color = "darkgray") +
    geom_hline(yintercept=-log10(fdr), linetype="dashed", color = "gray") +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
       x = "Expected -log10(p-value)",
       y = "Observed -log10(p-value)") +
    scale_colour_d3() +
    theme_bw() + 
    facet_wrap(~group, nrow=2) +
    theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
    guides(color=guide_legend(title="Annotation", override.aes=list(size=4))) + 
    geom_text_repel(data=subset(qq_dat, is_annotate=="yes"), aes(label=gene_symbol), position = position_nudge_repel(y=0.3, x=0.1), size=4, max.overlaps = Inf) +
    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA),
            strip.text = element_text(size=14, face='bold')) +
    guides(color='none')
dev.off()








## AF filtering

## gnomad
ptv = read.delim('2026_01_07_ALL_SCHEMA_PTV_MAC15_gnomad_0.001_popmax_threshold_CMH_out.tsv')[c(1,18)]
ptv_mis = read.delim('2026_01_07_ALL_SCHEMA_PTV_MIS_MAC15_gnomad_0.001_popmax_threshold_CMH_out.tsv')[c(1,18)]

ptv = subset(ptv, RES.p_value > 0)
ptv_mis = subset(ptv_mis, RES.p_value > 0)

cc = rbind(ptv_mis, ptv)
colnames(cc)[2] = 'p_value'
cc = subset(cc, p_value > 0)
cc_minp = cc %>% group_by(gene_symbol) %>% slice(which.min(p_value))

colnames(cc_minp)[2] = 'cc_p_value'

gnomad_filtered = cc_minp

## rgc
ptv = read.delim('2026_01_07_ALL_SCHEMA_PTV_MAC15_rgc_0.001_popmax_threshold_CMH_out.tsv')[c(1,18)]
ptv_mis = read.delim('2026_01_07_ALL_SCHEMA_PTV_MIS_MAC15_rgc_0.001_popmax_threshold_CMH_out.tsv')[c(1,18)]

ptv = subset(ptv, RES.p_value > 0)
ptv_mis = subset(ptv_mis, RES.p_value > 0)

cc = rbind(ptv_mis, ptv)
colnames(cc)[2] = 'p_value'
cc = subset(cc, p_value > 0)
cc_minp = cc %>% group_by(gene_symbol) %>% slice(which.min(p_value))

colnames(cc_minp)[2] = 'cc_p_value'

rgc_filtered = cc_minp


# unfiltered
unfiltered = read.csv('case_control_minp_ptv_mismean_rank93_2026-01-07.csv')[c(1,2)]

get_qqp_data = function(input=input, group=group){
    colnames(input)[2] = 'pvalue_out'
    input = input[!duplicated(input$gene_symbol),]
    input$pvector = input$pvalue_out
    input <- subset(input, !is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector < 1 & pvector > 0)
    input$o = -log10(input$pvector)
    input = input[order(-input$o),]
    input$e = -log10(ppoints(length(input$pvector) ))
    input$group = paste0(group)
    return(input)
}

gnomad_qqp = get_qqp_data(input=gnomad_filtered, group='gnomAD popmax < 0.1%')
rgc_qqp = get_qqp_data(input=rgc_filtered, group='RGC popmax < 0.1%')
unfiltered = get_qqp_data(input=unfiltered, group='Unfiltered')

qq_dat = rbind(gnomad_qqp, rgc_qqp, unfiltered)

qq_dat$is_annotate = ifelse(qq_dat$pvalue_out <= bonf, 'yes', 'no')

qq_dat$group = factor(qq_dat$group, levels=c('Unfiltered', 'gnomAD popmax < 0.1%', 'RGC popmax < 0.1%'))

pdf('schema2_minp_comparison_with_popmax_qqp.pdf', width=15, height=8)
ggplot(qq_dat, aes(x = e, y = o, color=group)) +
    geom_hline(yintercept=-log10(bonf), linetype="dashed", color = "darkgray") +
    geom_hline(yintercept=-log10(fdr), linetype="dashed", color = "gray") +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
       x = "Expected -log10(p-value)",
       y = "Observed -log10(p-value)") +
    scale_colour_d3() +
    theme_bw() + 
    facet_wrap(~group, nrow=1) +
    theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
    guides(color=guide_legend(title="Annotation", override.aes=list(size=4))) + 
    geom_text_repel(data=subset(qq_dat, is_annotate=="yes"), aes(label=gene_symbol), position = position_nudge_repel(y=0.3, x=0.1), size=4, max.overlaps = Inf) +
    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA),
            strip.text = element_text(size=14, face='bold')) +
    guides(color='none')
dev.off()










## MAC 1

setwd('/Users/juliasealock/Desktop/schema/add_dragen_data_nov25/cmh/acat/mac1')

library(dplyr)
library(data.table)
date = Sys.Date()

library(data.table)

ptv = read.delim('2026_01_07_ALL_SCHEMA_PTV_MAC15_GQ25_ptv_AC1_CMH_pc_matched.tsv')[c(1,18,19)]
ptv_mis = read.delim('2026_01_07_ALL_SCHEMA_PTV_MIS_MAC15_GQ25_ptv-mis_mean_rank93_AC1_CMH_pc_matched.tsv')[c(1,18,19)]


ptv = subset(ptv, RES.p_value > 0)
ptv_mis = subset(ptv_mis, RES.p_value > 0)

## find minp of case-controls
ptv_mis$annotation = 'ptv_mis'
ptv$annotation = 'ptv'

cc = rbind(ptv_mis, ptv)
colnames(cc)[2] = 'p_value'
cc = subset(cc, p_value > 0)
cc_minp = cc %>% group_by(gene_symbol) %>% slice(which.min(p_value))
table(cc_minp$annotation)

colnames(cc_minp)[2] = 'cc_p_value'
colnames(cc_minp)[4] = 'cc_annotation'

write.csv(cc_minp, paste0('case_control_minp_ptv_mismean_rank93_mac1_',date,'.csv'), row.names=F, quote=F)

mac1 = read.csv('case_control_minp_ptv_mismean_rank93_mac1_2026-01-12.csv')[c(1,2,3)]
mac15 = read.csv('../case_control_minp_ptv_mismean_rank93_2026-01-07.csv')[c(1,2,3)]

colnames(mac1)[2] = 'mac1_pvalue'
colnames(mac15)[2] = 'mac15_pvalue'

colnames(mac1)[3] = 'mac1_or'
colnames(mac15)[3] = 'mac15_or'

dat = merge(mac1, mac15, by='gene_symbol')

pvalue_cor = cor.test(dat$mac1_pvalue, dat$mac15_pvalue, method='spearman')
# rho = 0.387; pvalue < 2.2e-16

bonf = 1.31e-6
fdr = 6.84e-05



a = ggplot(dat, aes(x = -log10(mac15_pvalue), y = -log10(mac1_pvalue))) + 
    geom_point() + 
    theme_bw() + 
    geom_abline(intercept = 0, slope = 1, color='red') + 
    geom_text_repel(data=subset(dat, mac15_pvalue <= fdr), aes(label=gene_symbol), position = position_nudge_repel(y=0.3, x=0.1), size=4, max.overlaps = Inf) + 
    geom_vline(xintercept = -log10(fdr), linetype='dashed') + 
    geom_hline(yintercept = -log10(fdr), linetype='dashed') +
    xlab('-log10(MAC <= 15 pvalue)') + ylab('-log10(Singleton pvalue)') +
    theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12))

get_qqp_data = function(input=input, group=group){
    colnames(input)[2] = 'pvalue_out'
    input = input[!duplicated(input$gene_symbol),]
    input$pvector = input$pvalue_out
    input <- subset(input, !is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector < 1 & pvector > 0)
    input$o = -log10(input$pvector)
    input = input[order(-input$o),]
    input$e = -log10(ppoints(length(input$pvector) ))
    input$group = paste0(group)
    return(input)
}

mac1 = read.csv('case_control_minp_ptv_mismean_rank93_mac1_2026-01-12.csv')[c(1,2)]
mac15 = read.csv('../case_control_minp_ptv_mismean_rank93_2026-01-07.csv')[c(1,2)]

mac1_qqp = get_qqp_data(input=mac1, group='Singletons')
mac15_qqp = get_qqp_data(input=mac15, group='MAC <= 15')

qq_dat = rbind(mac1_qqp, mac15_qqp)

qq_dat$is_annotate = ifelse(qq_dat$pvalue_out <= bonf, 'yes', 'no')


b = ggplot(qq_dat, aes(x = e, y = o, color=group)) +
    geom_hline(yintercept=-log10(bonf), linetype="dashed", color = "darkgray") +
    geom_hline(yintercept=-log10(fdr), linetype="dashed", color = "gray") +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
       x = "Expected -log10(p-value)",
       y = "Observed -log10(p-value)") +
    scale_colour_d3() +
    theme_bw() + 
    facet_wrap(~group, nrow=1) +
    theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
    guides(color=guide_legend(title="Annotation", override.aes=list(size=4))) + 
    geom_text_repel(data=subset(qq_dat, is_annotate=="yes"), aes(label=gene_symbol), position = position_nudge_repel(y=0.3, x=0.1), size=4, max.overlaps = Inf) +
    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA),
            strip.text = element_text(size=14, face='bold')) +
    guides(color='none')


png('mac1_vs_mac15.png', width=8, height=12, unit='in', res=300)
ggarrange(a, b, labels=c('A.','B.'), nrow=1)
dev.off()



