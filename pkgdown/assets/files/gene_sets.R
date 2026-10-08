# Small controlled prioritization experiment; no annotation downloads.
library(gsim)
simulation <- gsim(n = 64, m = 24, n_causal = 4, seed = 812)
# Explicitly synthetic, one gene per SNP: four causal and twenty noncausal genes.
# Substitute a real data.frame(snp = ..., gene = ...) with matching SNP IDs.
# One SNP may occur in several mapping rows for different genes.
snp_gene_map <- data.frame(snp = rownames(simulation$B),
                            gene = paste0("synthetic_gene_", seq_len(nrow(simulation$B))))
gene_sets <- gsim_gene_sets(snp_gene_map, simulation$causal_rsids,
                            set_size = 12, n_causal = c(0, 2, 4),
                            n_sets = 3, seed = 813)
# Sets are deliberately conditional on causal truth, not significance labels.
print(gene_sets$truth)
# Two size-12 null sets from twenty noncausal genes must overlap.
print(intersect(gene_sets$sets[[1]], gene_sets$sets[[2]]))
# simulation$B retains all marker effects, unchanged by set construction.
