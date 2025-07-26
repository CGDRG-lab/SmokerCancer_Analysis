library(sva)
library(dplyr)
library(readr)
library(tibble)
library(limma)

## set working directory
setwd("~/Reuben/Biomarker Finding/for_cgdrg/preprocessing")

source('../algorithms/analysis2.R')
## =============================================================================
## ---------------------------- DATA LOADING GSE -------------------------------
## =============================================================================
g_edata <- read_edata('../mutual-data/gse_prepr_edata.tsv')
g_mdata <- read_mdata('../mutual-data/gse_prepr_mdata.tsv')

m_edata <- read_edata('../mutual-data/emtab_prepr_edata.tsv')
m_mdata <- read_mdata('../mutual-data/emtab_prepr_mdata.tsv')


all(rownames(m_edata) == rownames(g_edata))

table(m_mdata$condition)
table(g_mdata$condition)

# --------- CHECK FOR OUTLIEAR SAMPLES
find_outliar_samples(g_edata, g_mdata$condition)
find_outliar_samples(m_edata, m_mdata$condition)

## ------------- checking the data distributions
par(mfrow=c(1,1))
boxplot(g_edata)
plotDensities(g_edata, group = g_mdata$condition)
plotDensities(m_edata, group = m_mdata$condition)

# g_edata <- normalizeBetweenArrays(g_edata, method = 'quantile')
# m_edata <- normalizeBetweenArrays(m_edata, method = 'quantile')

## ------------- Adding batch info
m_mdata$batch <- 1
g_mdata$batch <- 2


## ------------- merge the data
edata <- cbind(m_edata, g_edata)
mdata <- rbind(
	m_mdata[,c('age', 'sex', 'condition', 'batch')],
	g_mdata[,c('age', 'sex', 'condition', 'batch')]
)
table(mdata$condition)

all(colnames(edata) == rownames(mdata))

## ------------- checking the data distributions of merged data
## before batch correction
# par(mar=c(6, 4, 3, 2))
par(mfrow=c(2, 2), mar=c(6, 5, 2, 1))
boxplot(edata, col = mdata$batch+1, las=2, main = "Boxplot of Merged Data (Before)")
plotDensities(edata, group = mdata$batch, main = "Distribution of Merged Data (Before)")


plot_pca(edata, mdata$condition, xlim=c(-65, 110), lpos = 'topleft', cex=0.8,
		 main = "PCA of Merged Data (Before)")

clst <- hclust( dist(t(edata)), method = 'average' )
plot(clst, main = "Heirarchical clustering of Merged Data (Before)", cex=0.8)

## ------------- batch effect correction
## design matrix
design_mat <- model.matrix(~condition, data = mdata)

edata_combat <- ComBat(
	dat = edata,
	batch = mdata$batch, 
	mod = design_mat, 
	par.prior = TRUE,
	prior.plots = FALSE
)


## ------------- checking outliear samples in merged data
find_outliar_samples(edata_combat, mdata$condition)

## ------------- checking the data distributions of merged data
## after batch correction
# par(mar=c(6, 4, 3, 2))
par(mfrow=c(2, 2), mar=c(6, 5, 2, 1))
boxplot(edata_combat, col = mdata$batch+1, las=2, main = "Boxplot of Merged Data (After)")
plotDensities(edata_combat, group = mdata$batch, main = "Distribution of Merged Data (After)")

plot_pca(edata_combat, mdata$condition, xlim=c(-100, 100), lpos = 'topleft',
		 main = "PCA of Merged Data (After)", cex=0.8)


clst <- hclust( dist(t(edata_combat)), method = 'average' )
plot(clst, main = "Heirarchical clustering of Merged Data (After)", cex=0.8)



## ------------- just checking expression of some gene before and after correction
h_genes <- c('GAPDH', 'ACTB', 'B2M', 'TBP')
par(mfrow=c(2, 4), mar=c(3, 5, 4, 3))
for(g in h_genes){
	boxplot(edata[g, ] ~ mdata$batch, main = paste(g, " - before"))
	boxplot(edata_combat[g, ] ~ mdata$batch, main = paste(g, " - after"))
}


## ------------- differential gene expression
w <- arrayWeights(edata_combat, design = design_mat)
fit <- lmFit(edata_combat, design = design_mat)
fit2 <- eBayes(fit)



result <- topTable(fit2, coef = 2, number = Inf) %>% 
	rownames_to_column('Gene')


plot_volcano(result, lab=T, minFC = -0.8, maxFC = 0.8, by='P.Value')

par(mfrow=c(1, 1))
plotMD(fit2)
abline(h=0, lwd=2, col='red')

addmargins( table(mdata$batch ,mdata$condition) )

## --------------------------------- Write The Results
edata_combat %>% 
	as.data.frame() %>% 
	rownames_to_column('Gene') %>% 
	write_tsv('../mutual-data/merged_edata.tsv')

mdata %>% 
	as.data.frame() %>% 
	rownames_to_column('Sample') %>% 
	write_tsv('../mutual-data/merged_mdata.tsv')
















