
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
de_novo_v2_results = 'schema1_de_novo_poisson_results_with_gnomad_v2_rates_with_misrank.tsv.bgz'

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

qq_dat = rbind(ptv_qqp, ptv_mis_qqp, minp_qqp, minp_dn_qqp)


bonf = 1.31e-6
fdr = 6.84e-05


qq_dat$is_annotate = ifelse(qq_dat$pvalue_out <= bonf, 'yes', 'no')

qq_dat$group = factor(qq_dat$group, levels=c('PTV', 'PTV + Missense', 'Lowest Pvalue', 'De Novo Meta-Analysis'))

pdf('schema2_qqp.pdf', width=18, height=6)
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





