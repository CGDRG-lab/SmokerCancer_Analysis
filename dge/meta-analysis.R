library(dplyr)
library(readr)



setwd("~/Reuben/Biomarker Finding/for_cgdrg/dge")


m_result <- read_tsv('../mutual-data/emtab_top_table.tsv')
g_result <- read_tsv('../mutual-data/gse_top_table.tsv')


merged_df <- inner_join(m_result, g_result, by = "Gene", suffix = c("_m", "_g"))

## ------------------------ Stouffer’s Z-Score Method
weight_m <- 1 #3.26,
weight_g <- 1 #2.97,

merged_df <- merged_df %>%
	mutate(
		logFC_mean = (logFC_m + logFC_g) / 2,
		z_dir_m = sign(logFC_m) * qnorm(1 - P.Value_m/2),
		z_dir_g = sign(logFC_g) * qnorm(1 - P.Value_g/2),
		stouffer_z = (weight_m * z_dir_m + weight_g * z_dir_g) / sqrt(weight_m^2 + weight_g^2),
		stouffer_pval = 2 * pnorm(-abs(stouffer_z)),
		stouffer_fdr = p.adjust(stouffer_pval, method = "fdr")
	)
## -------------------------- Just plotting the pvalues
par(mfrow=c(1,1))
hist(merged_df$stouffer_pval, breaks=100)
hist(merged_df$P.Value_m, breaks = 20)
hist(merged_df$P.Value_g, breaks = 20)
hist(merged_df$stouffer_fdr, breaks = 20)

## -------------------------- Filter
sig_genes_meta <- merged_df %>%
	# filter(stouffer_fdr < 0.1) %>%
	filter(stouffer_pval < 0.05) %>%
	pull(Gene)

sort(sig_genes_meta)



## ------------------write the stouffer sig genes
writeLines(sort(sig_genes_meta), 'sig_genes_meta.txt')

## ------------------write the stouffer results
merged_df %>%
	select(
		Gene, logFC_m, logFC_g, P.Value_m, P.Value_g, logFC_mean,
		stouffer_z, stouffer_pval, stouffer_fdr
		#z_dir_m, z_dir_g
	) %>%
	arrange(stouffer_pval) %>%
	mutate(
		across(where(is.numeric), function(x){round(x, 4)})
	) %>%
	writexl::write_xlsx('stouffer_result.xlsx')














