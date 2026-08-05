

## run capture x pop and meta

pcs$scores = gsub("\\[|\\]", "", pcs$scores)
pcs = separate(data = pcs, col = scores, into = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), sep = ",")
pcs[,2:11] <- sapply(pcs[,2:11],as.numeric)
meta = merge(meta, pcs, by='s')

mac1$ptv_misrank93 = mac1$total_plof_os + mac1$total_misrank93
mac2_5$ptv_misrank93 = mac2_5$total_plof_os + mac2_5$total_misrank93
mac6_10$ptv_misrank93 = mac6_10$total_plof_os + mac6_10$total_misrank93
mac11_15$ptv_misrank93 = mac11_15$total_plof_os + mac11_15$total_misrank93

meta$pop = ifelse(meta$pop == 'SAS', 'CSA', meta$pop)
meta$group2 = paste0(meta$pop,'_',meta$capture)

## functions

run_glms2 = function(dat=dat, strata=strata){
    dat1 = merge(dat, meta, by='s')
    dat1 = subset(dat1, group2==strata)
    dat1$scz_status = ifelse(dat1$status=="schizophrenia",1,0)
    dat1$scz_status = as.factor(dat1$scz_status)

    plof_os_glm = summary(glm(scz_status ~ total_plof_os + n_total_singleton + PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family=binomial))
    plof_misrank93_glm = summary(glm(scz_status ~ ptv_misrank93 + n_total_singleton + PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family=binomial))
    misrank93_glm = summary(glm(scz_status ~ total_misrank93 +  n_total_singleton +  PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family=binomial))
    syn_glm = summary(glm(scz_status ~ total_syn + n_total_singleton + PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family=binomial))

    plof_os = coef(plof_os_glm)[2,]
    plof_misrank93 = coef(plof_misrank93_glm)[2,]
    misrank93 = coef(misrank93_glm)[2,]
    syn = coef(syn_glm)[2,]

    out = rbind(plof_misrank93, plof_os, misrank93, syn)
    out = as.data.frame(out)
    colnames(out) = c("Beta","SE","Z","Pvalue")
    out$OR = exp(out$Beta)
    out$Lower.CI = exp(out$Beta - (1.96*out$SE))
    out$Upper.CI = exp(out$Beta + (1.96*out$SE))
    out$annotation = rownames(out)
    out$strata = paste0(strata)
    return(out)
}

### mac1
results_list <- list()
for(strata in unique(meta$group2)){
    message("Running: ", strata)
    results_list[[strata]] <- run_glms2(
        dat = mac1,
        strata = strata
    )
}
results_mac1 <- do.call(rbind, results_list)


## mac2-5
results_list <- list()
for(strata in unique(meta$group2)){
    message("Running: ", strata)
    results_list[[strata]] <- run_glms2(
        dat = mac2_5,
        strata = strata
    )
}
results_mac2_5 <- do.call(rbind, results_list)


## mac6-10
results_list <- list()
for(strata in unique(meta$group2)){
    message("Running: ", strata)
    results_list[[strata]] <- run_glms2(
        dat = mac6_10,
        strata = strata
    )
}
results_mac6_10 <- do.call(rbind, results_list)

## mac11-15
results_list <- list()
for(strata in unique(meta$group2)){
    message("Running: ", strata)
    results_list[[strata]] <- run_glms2(
        dat = mac11_15,
        strata = strata
    )
}
results_mac11_15 <- do.call(rbind, results_list)


## make MAC 1 - 15 file
colnames(mac1)[2:5] = paste0('mac1_',colnames(mac1)[2:5])
colnames(mac2_5)[2:5] = paste0('mac2_5_',colnames(mac2_5)[2:5])
colnames(mac6_10)[2:5] = paste0('mac6_10_',colnames(mac6_10)[2:5])
colnames(mac11_15)[2:5] = paste0('mac11_15_',colnames(mac11_15)[2:5])

mac1_5 = merge(mac1, mac2_5, by='s')
mac1_10 = merge(mac1_5, mac6_10, by='s')
mac1_15 = merge(mac1_10, mac11_15, by='s')

mac1_15$total_plof_os = mac1_15$mac1_total_plof_os + mac1_15$mac2_5_total_plof_os + mac1_15$mac6_10_total_plof_os + mac1_15$mac11_15_total_plof_os
mac1_15$total_misrank93 = mac1_15$mac1_total_misrank93 + mac1_15$mac2_5_total_misrank93 + mac1_15$mac6_10_total_misrank93 + mac1_15$mac11_15_total_misrank93
mac1_15$ptv_misrank93 = mac1_15$mac1_ptv_misrank93+ mac1_15$mac2_5_ptv_misrank93 + mac1_15$mac6_10_ptv_misrank93+ mac1_15$mac11_15_ptv_misrank93
mac1_15$total_syn = mac1_15$mac1_total_syn + mac1_15$mac2_5_total_syn + mac1_15$mac6_10_total_syn + mac1_15$mac11_15_total_syn

mac1_15 = mac1_15[,c('s', 'total_plof_os', 'total_misrank93', 'ptv_misrank93', 'total_syn')]

## mac1-15
results_list <- list()
for(strata in unique(meta$group2)){
    message("Running: ", strata)
    results_list[[strata]] <- run_glms2(
        dat = mac1_15,
        strata = strata
    )
}
results_mac1_15 <- do.call(rbind, results_list)

## write
results_mac1$MAC = 'MAC1'
results_mac2_5$MAC = 'MAC2_5'
results_mac6_10$MAC = 'MAC6_10'
results_mac11_15$MAC = 'MAC11_15'
results_mac1_15$MAC = 'MAC1_15'

all_results = rbind(results_mac1, results_mac2_5, results_mac6_10, results_mac11_15, results_mac1_15)
write.csv(all_results, 'enrichment_results_by_mac_x_anno_x_capture_ancestry.csv', row.names=F, quote=F, col.names=F)


library(meta)
run_meta = function(data, cat, anc){
    # Subset to one ancestry (pop) and one annotation
    studies = subset(data, annotation == cat & pop == anc)
    if (nrow(studies) == 0) return(NULL)
    meta_analysis <- metagen(
        TE    = studies$Beta,     # log(OR) or Beta
        seTE  = studies$SE,       # SE of Beta
        studlab = studies$capture,  # Label by capture platform
        sm    = "SMD",            # keep as in your code (or change if Beta is log(OR))
        method.tau = "DL",
        fixed = TRUE,
        random = FALSE
    )
    Pvalue   <- meta_analysis$pval.fixed
    Beta     <- meta_analysis$TE.fixed
    SE       <- meta_analysis$seTE.fixed
    Lower.CI <- meta_analysis$lower.fixed
    Upper.CI <- meta_analysis$upper.fixed
    out <- data.frame(
        annotation = cat,
        pop        = anc,
        Pvalue     = Pvalue,
        Beta       = Beta,
        SE         = SE,
        Lower.CI   = Lower.CI,
        Upper.CI   = Upper.CI,
        stringsAsFactors = FALSE
    )
    # Convert Beta to OR
    out$OR        <- exp(out$Beta)
    out$Lower.CI  <- exp(out$Lower.CI)
    out$Upper.CI  <- exp(out$Upper.CI)
    return(out)
}

combine_meta = function(input2){
    # annotations you care about
    annos <- c("plof_os", "misrank93", "syn", "plof_misrank93")
    # ancestries / pops present
    pops  <- unique(input2$pop)
    out_list <- list()
    k <- 1

    for (anc in pops){
        for (cat in annos){
            tmp <- run_meta(data = input2, cat = cat, anc = anc)
            if (!is.null(tmp)){
                out_list[[k]] <- tmp
                k <- k + 1
            }
        }
    }
    meta_out <- do.call(rbind, out_list)
    # Nice factor ordering for ancestry if you want
    meta_out$pop <- factor(
        meta_out$pop,
        levels = c("NFE","FIN","EAS","AFR","AMR","CSA","ASJ")
    )
    # Add any convenience columns you were using before
    meta_out$x <- "pLI > 0.9"
    meta_out$anno <- ifelse(meta_out$annotation == "plof_os",        "PTV",
                     ifelse(meta_out$annotation == "misrank93",      "Missense Rank 93%",
                     ifelse(meta_out$annotation == "plof_misrank93", "PTV + Missense Rank 93%",
                                                                  "Synonymous")))
    meta_out$anno <- factor(
        meta_out$anno,
        levels = c("PTV + Missense Rank 93%", "PTV", "Missense Rank 93%", "Synonymous")
    )
    return(meta_out)
}



results_mac1 = separate(data = results_mac1, col = strata, into = c("pop","capture"), sep = "_")
meta_results_mac1 <- combine_meta(results_mac1)
meta_results_mac1$MAC = 'MAC = 1'

results_mac2_5 = separate(data = results_mac2_5, col = strata, into = c("pop","capture"), sep = "_")
meta_results_mac2_5 <- combine_meta(results_mac2_5)
meta_results_mac2_5$MAC = '2 <= MAC <= 5'

results_mac6_10 = separate(data = results_mac6_10, col = strata, into = c("pop","capture"), sep = "_")
meta_results_mac6_10 <- combine_meta(results_mac6_10)
meta_results_mac6_10$MAC = '6 <= MAC <= 10'

results_mac11_15 = separate(data = results_mac11_15, col = strata, into = c("pop","capture"), sep = "_")
meta_results_mac11_15 <- combine_meta(results_mac11_15)
meta_results_mac11_15$MAC = '11 <= MAC <= 15'


results_mac1_15 = separate(data = results_mac1_15, col = strata, into = c("pop","capture"), sep = "_")
meta_results_mac1_15 <- combine_meta(results_mac1_15)
meta_results_mac1_15$MAC = 'MAC <= 15'


meta_ancestry = rbind(meta_results_mac1, meta_results_mac2_5, meta_results_mac6_10, meta_results_mac11_15, meta_results_mac1_15)
write.csv(meta_ancestry, ENRICHMENT_OUT, row.names=F, quote=F)




meta_across_ancestries <- function(data, mac = NULL){
    annos <- unique(data$annotation)
    out_list <- list()
    k <- 1
    for (cat in annos) {
        studies <- subset(data, annotation == cat)
        if (nrow(studies) == 0) next
        meta_analysis <- metagen(
            TE      = studies$Beta,
            seTE    = studies$SE,
            studlab = studies$pop,   # each ancestry is a study
            sm      = "SMD",
            method.tau = "DL",
            fixed   = TRUE,
            random  = FALSE
        )
        Pvalue   <- meta_analysis$pval.fixed
        Beta     <- meta_analysis$TE.fixed
        SE       <- meta_analysis$seTE.fixed
        Lower.CI <- meta_analysis$lower.fixed
        Upper.CI <- meta_analysis$upper.fixed
        out <- data.frame(
            annotation = cat,
            Pvalue     = Pvalue,
            Beta       = Beta,
            SE         = SE,
            Lower.CI   = Lower.CI,
            Upper.CI   = Upper.CI,
            stringsAsFactors = FALSE
        )
        if (!is.null(mac)) out$MAC <- mac
        # OR scale
        out$OR        <- exp(out$Beta)
        out$Lower.CI  <- exp(out$Lower.CI)
        out$Upper.CI  <- exp(out$Upper.CI)
        out_list[[k]] <- out
        k <- k + 1
    }
    if (length(out_list) == 0) return(NULL)
    meta_out <- do.call(rbind, out_list)
    return(meta_out)
}

meta_results_mac1_global <- meta_across_ancestries(meta_results_mac1)
meta_results_mac2_5_global <- meta_across_ancestries(meta_results_mac2_5)
meta_results_mac6_10_global <- meta_across_ancestries(meta_results_mac6_10)
meta_results_mac11_15_global <- meta_across_ancestries(meta_results_mac11_15)
meta_results_mac1_15_global <- meta_across_ancestries(meta_results_mac1_15)

meta_results_mac1_global$MAC = 'MAC = 1'
meta_results_mac2_5_global$MAC = '2 <= MAC <= 5'
meta_results_mac6_10_global$MAC = '6 <= MAC <= 10'
meta_results_mac11_15_global$MAC = '11 <= MAC <= 15'
meta_results_mac1_15_global$MAC = 'MAC <= 15'

meta_results = rbind(meta_results_mac1_global, meta_results_mac2_5_global, meta_results_mac6_10_global, 
                        meta_results_mac11_15_global, meta_results_mac1_15_global)
meta_results$anno <- ifelse(meta_results$annotation == "plof_os",        "PTV",
                    ifelse(meta_results$annotation == "misrank93",      "Missense Rank 93%",
                    ifelse(meta_results$annotation == "plof_misrank93", "PTV + Missense Rank 93%",
                                                                "Synonymous")))

write.csv(meta_results, ENRICHMENT_OUT2, row.names=F, quote=F)



cols1 = c('PTV' = '#d90429', 'PTV + Missense Rank 93%' = "#fc9f5b", 'Missense Rank 93%' = "#129490", 'Synonymous' = '#9fb1bc') 
meta_results$anno = factor(meta_results$anno, levels=c('PTV','PTV + Missense Rank 93%', 'Missense Rank 93%','Synonymous'))
meta_results$MAC = factor(meta_results$MAC, levels=c('MAC <= 15', 'MAC = 1', '2 <= MAC <= 5', '6 <= MAC <= 10', '11 <= MAC <= 15'))
meta_results$x = 'pLI > 0.9'

png('enrichment_analysis_meta_across_ancestry.png', width=14, height=6, unit='in', res=300)
ggplot(data=meta_results, aes(x=x, y=OR, ymin=Lower.CI, ymax=Upper.CI, fill=anno, colour=anno)) +
    geom_pointrange(position=position_dodge(width=2), size=1, linewidth=1.1) + 
    geom_hline(yintercept=1, lty=2) + 
    ggtitle(" ") +
    theme_bw() +
    xlab("") +
    theme(axis.text.x=element_text(size=12), axis.title.x=element_text(size=12)) +
    theme(axis.text.y=element_text(size=12), axis.title.y=element_text(size=12)) +
    theme(legend.title = element_text(size = 12), legend.text = element_text(size=12)) +
    theme(strip.text.x = element_text(size = 12, face='bold')) + 
    ylab("Odds Ratio") +
    # scale_color_npg() +
    scale_colour_manual(values = cols1) +
    scale_y_continuous(n.breaks=10) +
    facet_wrap(~MAC, nrow=1) +
    guides(colour=guide_legend(title="Annotation"), fill='none') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    	    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) 
dev.off()



## supplement: 

library(tidyverse)
library(ggplot2)
library(ggpubr)

# A. by MAC

cols1 = c('PTV' = '#d90429', 'PTV + Missense Rank 93%' = "#fc9f5b", 'Missense Rank 93%' = "#129490", 'Synonymous' = '#9fb1bc') 
meta_results$anno = factor(meta_results$anno, levels=c('PTV','PTV + Missense Rank 93%', 'Missense Rank 93%','Synonymous'))
meta_results$MAC = factor(meta_results$MAC, levels=c('MAC <= 15', 'MAC = 1', '2 <= MAC <= 5', '6 <= MAC <= 10', '11 <= MAC <= 15'))
meta_results$x = 'pLI > 0.9'

a = ggplot(data=meta_results, aes(x=x, y=OR, ymin=Lower.CI, ymax=Upper.CI, fill=anno, colour=anno)) +
    geom_pointrange(position=position_dodge(width=2), size=1, linewidth=1.1) + 
    geom_hline(yintercept=1, lty=2) + 
    ggtitle(" ") +
    theme_bw() +
    xlab("") +
    theme(axis.text.x=element_text(size=12), axis.title.x=element_text(size=12)) +
    theme(axis.text.y=element_text(size=12), axis.title.y=element_text(size=12)) +
    theme(legend.title = element_text(size = 12), legend.text = element_text(size=12)) +
    theme(strip.text.x = element_text(size = 10, face='bold')) + 
    ylab("Odds Ratio") +
    # scale_color_npg() +
    scale_colour_manual(values = cols1) +
    scale_y_continuous(n.breaks=10) +
    facet_wrap(~MAC, nrow=1) +
    guides(colour=guide_legend(title="Annotation", position = 'top'), fill='none') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    	    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) 



# B. by ancestry 
cols <- c("Meta-analysis" = "black", "NFE" = "#1F77B4FF", "FIN" = "#FF7F0EFF", "EAS" = "#2CA02CFF", 'AFR' = '#D62728FF', 'AMR' = '#9467BDFF', 'CSA' = '#8C564BFF', 'ASJ' = '#E377C2FF', 'MID' = '#BCBD22FF')

mac1_pop_meta = subset(meta_ancestry, MAC == 'MAC = 1')
mac1_pop_meta$anno = factor(mac1_pop_meta$anno, levels=c('PTV','PTV + Missense Rank 93%', 'Missense Rank 93%','Synonymous'))
mac1_pop_meta$pop = factor(mac1_pop_meta$pop, levels=c('NFE','FIN','EAS','AFR','AMR','CSA','ASJ'))

b = ggplot(data=mac1_pop_meta, aes(x=x, y=OR, ymin=Lower.CI, ymax=Upper.CI, fill=pop, colour=pop)) +
    geom_pointrange(position=position_dodge(width=2), size=1, linewidth=1.1) + 
    geom_hline(yintercept=1, lty=2) + 
    ggtitle(" ") +
    theme_bw() +
    xlab("") +
    theme(axis.text.x=element_text(size=12), axis.title.x=element_text(size=12)) +
    theme(axis.text.y=element_text(size=12), axis.title.y=element_text(size=12)) +
    theme(legend.title = element_text(size = 12), legend.text = element_text(size=12)) +
    theme(strip.text.x = element_text(size = 10, face='bold')) + 
    ylab("Odds Ratio") +
    scale_colour_manual(values = cols) +
    scale_y_continuous(n.breaks=10) +
    facet_wrap(~anno, nrow=1) +
    guides(colour=guide_legend(title="Ancestry", position='top'), fill='none') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    	    theme(
	    	panel.grid.major = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) 


# C. by capture 
mac1_ptv = subset(all_results, MAC == 'MAC1' & annotation =='plof_os')
cols <- c("NFE" = "#1F77B4FF", "FIN" = "#FF7F0EFF", "EAS" = "#2CA02CFF", 'AFR' = '#D62728FF', 'AMR' = '#9467BDFF', 'CSA' = '#8C564BFF', 'ASJ' = '#E377C2FF', 'MID' = '#BCBD22FF')

mac1_ptv = separate(data = mac1_ptv, col = strata, into = c("pop", 'capture'), sep = "_")
mac1_ptv$capture = ifelse(mac1_ptv$capture=='BGE', 'BGE-GATK', ifelse(mac1_ptv$capture == 'DRAGEN', 'BGE-DRAGEN', mac1_ptv$capture))
mac1_ptv$capture = factor(mac1_ptv$capture, levels=c('Twist', 'Nextera', 'Agilent', 'WGS', 'BGE-GATK', 'BGE-DRAGEN'))
mac1_ptv$x = 'pLI > 0.9'
c = ggplot(data=mac1_ptv, aes(x=x, y=OR, ymin=Lower.CI, ymax=Upper.CI, fill=pop, colour=pop)) +
    geom_pointrange(position=position_dodge(width=2), size=1, linewidth=1.1) + 
    geom_hline(yintercept=1, lty=2) + 
    ggtitle(" ") +
    theme_bw() +
    xlab("") +
    theme(axis.text.x=element_text(size=12), axis.title.x=element_text(size=12)) +
    theme(axis.text.y=element_text(size=12), axis.title.y=element_text(size=12)) +
    theme(legend.title = element_text(size = 12), legend.text = element_text(size=12)) +
    theme(strip.text.x = element_text(size = 12, face='bold')) + 
    ylab("Odds Ratio") +
    scale_colour_manual(values = cols) +
    scale_y_continuous(n.breaks=10) +
    # facet_wrap(~anno, nrow=2) +
    facet_wrap(~capture) +
    guides(colour=guide_legend(title="Ancestry"), fill='none') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    	    theme(
	    	panel.grid.major = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) 

