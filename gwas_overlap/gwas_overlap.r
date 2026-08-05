bonf = 1.31e-6
fdr = 6.84e-05


all_meta = read.csv(META)

fdr_genes = subset(all_meta, cc_p_value <= fdr)

gwas_index = read.delim(GWAS)
gwas_index2 = gwas_index[,c('locus','alleles','top.index')]

gwas_index2$locus = gsub('chr', '', gwas_index2$locus)
gwas_index2 = separate(data = gwas_index2, col = locus, into = c("chr",'bp'), sep = ":")
colnames(gwas_index2)[1] = 'chromosome'
gwas_index2$chromosome = as.integer(gwas_index2$chromosome)
gwas_index2$bp = as.numeric(gwas_index2$bp)

pli = read.delim("gnomad.v4.0.constraint_metrics.tsv", header=T)
pli = pli[,c('gene','transcript','mane_select')]
pli = subset(pli, mane_select=='true')
pli$transcript = gsub("\\..*", "", pli$transcript)


loc = read.csv("refseq_grch38_gene_locations", header=T)
colnames(loc)[1] = 'transcript'
colnames(loc)[7] = 'gene'

loc = merge(loc, pli, by=c('transcript','gene'))
loc = loc[,c('gene','chrom','txStart','txEnd')]
loc$chrom = gsub("\\_.*", "", loc$chrom)

loc$chromosome = gsub('chr','',loc$chr)
loc$chromosome = as.integer(loc$chromosome)
loc$txStart = as.numeric(loc$txStart)
loc$txEnd = as.numeric(loc$txEnd)



result <- gwas_index2 %>%
  inner_join(loc, by = "chromosome") %>%
  filter((bp >= txStart & bp <= txEnd) |  # bp is within gene range
         (abs(bp - txStart) <= 5e5 |          # bp is within .5 megabase of start
          abs(bp - txEnd) <= 5e5))  

result_fdr = result[(result$gene %in% fdr_genes$gene_symbol),]

gwas = read.delim('supp_table3_ld_indep_loci.txt', header=T) # index snps from extended gwas 
gwas = gwas[,c('top.index','top.P','top.alleles','top.OR','top.SE')]


## permute 
set.seed(42)

eligible_genes <- loc

# Sample from unique gene names to avoid over-weighting genes with multiple rows/transcripts
gene_universe <- unique(eligible_genes$gene)

# Function to count overlaps for a random set of genes
count_random_overlaps2 <- function(n_genes = 40, window_bp = 5e5) {

  sampled_genes <- sample(gene_universe, n_genes, replace = FALSE)
  sampled_locs  <- eligible_genes[eligible_genes$gene %in% sampled_genes, ]

  # Join GWAS index SNPs to sampled gene coordinates (same chromosome)
  joined <- gwas_index2 %>%
    inner_join(sampled_locs, by = "chromosome") %>%
    filter(
      (bp >= txStart & bp <= txEnd) |
      (abs(bp - txStart) <= window_bp) |
      (abs(bp - txEnd)   <= window_bp)
    )

  # Count unique overlapping genes
  length(unique(joined$gene))
}


# Run simulation
n_sim <- 10000
overlap_counts <- replicate(n_sim, count_random_overlaps2())

# Compare to observed
# observed_count <- length(unique(result_fdr$gene))
observed_count = 6

p_value <- (sum(overlap_counts >= observed_count) + 1) / (length(overlap_counts) + 1)


hist(overlap_counts, main = "Expected Overlaps from 10k Overlap Simulations",
     xlab = "Number of Overlapping Genes", breaks = 50)
abline(v = observed_count, col = "red", lwd = 2)


