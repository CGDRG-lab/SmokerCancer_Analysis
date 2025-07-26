library(dendextend)
library(readr)
library(tibble)

setwd("~/Reuben/Biomarker Finding/for_cgdrg/preprocessing")

source('../algorithms/analysis2.R')

## =============================================================================
## ---------------------------- DATA LOADING -------------------------------
## =============================================================================
g_edata <- read_edata('../mutual-data/gse_raw_gene_collapsed.tsv')
g_mdata <- read_mdata('../mutual-data/gse_raw_mdata.tsv')

m_edata <- read_edata('../mutual-data/emtab_raw_gene_collapsed.tsv')
m_mdata <- read_mdata('../mutual-data/emtab_raw_mdata.tsv')



## ---------- OUTLIEAR FINDING EMTAB

# plot data before filtering
par(mfrow=c(2, 2))
clst <- hclust(dist(t(scale(m_edata))), method='average')
dnd <- as.dendrogram(clst, hang = 0.1)
dnd <- color_labels(dnd, col = ifelse(m_mdata$condition[clst$order]=='SC', 1, 2))
plot(dnd, cex = 0.1, ylab='Height', main='Cluster Dendogram E-MTAB-1690 (Before)')
legend('topleft', c('SC', 'SNC'), col=c(1, 2), pch=19, cex = 0.8)

plot_pca(
	m_edata, m_mdata$condition, cex=0.8, lpos = 'bottomleft',
	main='PCA plot (Before)'
)

# find outlier samples
m_outlier <- find_outliar_samples(m_edata, m_mdata$condition, f=2)
m_outlier
m_mdata_filt <- m_mdata[! rownames(m_mdata) %in% m_outlier, ]
table(m_mdata_filt$condition)

m_edata_filt <- m_edata[, rownames(m_mdata_filt)]

# plot data after filtering
clst <- hclust(dist(t(scale(m_edata_filt))), method='average')
dnd <- as.dendrogram(clst, hang = 0.1)
dnd <- color_labels(dnd, col = ifelse(m_mdata_filt$condition[clst$order]=='SC', 1, 2))
plot(dnd, cex = 0.1, ylab='Height', main='Cluster Dendogram E-MTAB-1690 (After)')
legend('topright', c('SC', 'SNC'), col=c(1, 2), pch=19, cex = 0.8)

plot_pca(
	m_edata_filt, m_mdata_filt$condition, cex=0.8, lpos = 'bottomright',
	main='PCA plot (After)'
)




## ---------- OUTLIEAR FINDING GSE

# plot data before filtering
par(mfrow=c(2, 2))
clst <- hclust(dist(t(scale(g_edata))), method='average')
dnd <- as.dendrogram(clst, hang = 0.1)
dnd <- color_labels(dnd, col = ifelse(g_mdata$condition[clst$order]=='SC', 1, 2))
plot(dnd, cex = 0.1, ylab='Height', main='Cluster Dendogram GSE19027 (Before)')
legend('topleft', c('SC', 'SNC'), col=c(1, 2), pch=19, cex = 0.8)

plot_pca(
	g_edata, g_mdata$condition, cex=0.8, lpos = 'topleft',
	main='PCA plot (Before)'
)

# find outlier samples
g_outlier <- find_outliar_samples(g_edata, g_mdata$condition, f=2)
g_outlier
g_outlier <- c(g_outlier, "GSM470865", "GSM470866")

g_mdata_filt <- g_mdata[! rownames(g_mdata) %in% g_outlier, ]
table(g_mdata_filt$condition)

g_edata_filt <- g_edata[, rownames(g_mdata_filt)]

# plot data after filtering
clst <- hclust(dist(t(scale(g_edata_filt))), method='average')
dnd <- as.dendrogram(clst, hang = 0.1)
dnd <- color_labels(dnd, col = ifelse(g_mdata_filt$condition[clst$order]=='SC', 1, 2))
plot(dnd, cex = 0.1, ylab='Height', main='Cluster Dendogram GSE19027 (After)')
legend('topright', c('SC', 'SNC'), col=c(1, 2), pch=19, cex = 0.8)



plot_pca(
	g_edata_filt, g_mdata_filt$condition, cex=0.8, lpos = 'bottomleft', 
	main='PCA plot (After)'
)


## ---------- WRITE DATA
g_edata_filt %>%
	as.data.frame() %>%
	rownames_to_column('Gene') %>%
	write_tsv('../mutual-data/gse_prepr_edata.tsv')

m_edata_filt %>%
	as.data.frame() %>%
	rownames_to_column('Gene') %>%
	write_tsv('../mutual-data/emtab_prepr_edata.tsv')

g_mdata_filt %>%
	rownames_to_column('Sample') %>%
	write_tsv('../mutual-data/gse_prepr_mdata.tsv')

m_mdata_filt %>%
	rownames_to_column('Sample') %>%
	write_tsv('../mutual-data/emtab_prepr_mdata.tsv')





















