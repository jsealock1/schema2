

## run capture x pop and meta
# setwd('')
library(tidyverse)
library(ggplot2)


pcs$scores = gsub("\\[|\\]", "", pcs$scores)
pcs = separate(data = pcs, col = scores, into = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), sep = ",")
pcs[,2:11] <- sapply(pcs[,2:11],as.numeric)
meta = merge(meta, pcs, by='s')

asd$ptv_misrank93 = asd$total_plof_os + asd$total_misrank93
dd$ptv_misrank93 = dd$total_plof_os + dd$total_misrank93

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

### asd
results_list <- list()
for(strata in unique(meta$group2)){
    message("Running: ", strata)
    results_list[[strata]] <- run_glms2(
        dat = asd,
        strata = strata
    )
}
results_asd <- do.call(rbind, results_list)


## dd
results_list <- list()
for(strata in unique(meta$group2)){
    message("Running: ", strata)
    results_list[[strata]] <- run_glms2(
        dat = dd,
        strata = strata
    )
}
results_dd <- do.call(rbind, results_list)






## write
results_asd$phenotype = 'ASD'
results_dd$phenotype = 'DD/ID'

all_results = rbind(results_asd, results_dd)
write.csv(all_results, DDID_OUT, row.names=F, quote=F, col.names=F)


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



results_asd = separate(data = results_asd, col = strata, into = c("pop","capture"), sep = "_")
meta_results_asd <- combine_meta(results_asd)
meta_results_asd$phenotype = 'ASD'

results_dd = separate(data = results_dd, col = strata, into = c("pop","capture"), sep = "_")
meta_results_dd <- combine_meta(results_dd)
meta_results_dd$phenotype= 'DD/ID'


meta_ancestry = rbind(meta_results_asd, meta_results_dd)
write.csv(meta_ancestry, META_OUT, row.names=F, quote=F)




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

meta_results_asd_global <- meta_across_ancestries(meta_results_asd)
meta_results_dd_global <- meta_across_ancestries(meta_results_dd)

meta_results_asd_global$phenotype = 'ASD'
meta_results_dd_global$phenotype = 'DD/ID'

meta_results = rbind(meta_results_asd_global, meta_results_dd_global)

meta_results$anno <- ifelse(meta_results$annotation == "plof_os",        "PTV",
                    ifelse(meta_results$annotation == "misrank93",      "Missense Rank 93%",
                    ifelse(meta_results$annotation == "plof_misrank93", "PTV + Missense Rank 93%",
                                                                "Synonymous")))

write.csv(meta_results, META_OUT2, row.names=F, quote=F)



cols1 = c('PTV' = '#d90429', 'PTV + Missense Rank 93%' = "#fc9f5b", 'Missense Rank 93%' = "#129490", 'Synonymous' = '#9fb1bc') 
meta_results$anno = factor(meta_results$anno, levels=c('PTV','PTV + Missense Rank 93%', 'Missense Rank 93%','Synonymous'))
meta_results$phenotype = factor(meta_results$phenotype, levels=c('DD/ID', 'ASD'))
meta_results$x = ''

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
    facet_wrap(~phenotype, nrow=1) +
    guides(colour=guide_legend(title="Annotation"), fill='none') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    	    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) 
