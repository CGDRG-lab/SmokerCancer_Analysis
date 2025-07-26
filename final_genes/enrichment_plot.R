library(dplyr)
library(readr)
library(readxl)

setwd("~/Reuben/Biomarker Finding/for_cgdrg/final_genes")

source('../algorithms/analysis2.R')

## ====================== Plot for darkturquoise module genes
res2 <- read_tsv("~/Reuben/Biomarker Finding/trial4/wgcna/darkt.tsv")

res2 %>% 
	# arrange(FDR) %>%
	filter(FDR <= 0.05) %>% 
	arrange(desc(Count)) %>%
	slice_head(n=12) %>% 
	lolipop_plot('Darkturquoise Module',
				 textwrap = 40,
				 generatio = F,
				 xmin = 0,
				 foldEn = "Fold Enrichment",
				 nGene = "Count",
				 pval = "PValue",
				 term = "Term")


## ====================== Enrichment Plot for Final 6 genes
# res <- read_tsv('./david_enrichment.tsv')
res <- read_excel('./ClueGOResultTable-mod.xls')

colnames(res)

res <- res %>% rename(Count = `Nr. Genes`,
					  FDR = `Term PValue Corrected with Benjamini-Hochberg`) %>% 
	mutate(Term = paste(ID, Term, sep=' ~ '))

colnames(res)

lolipop_plot(
	res,
	title = 'Pathway Enrichment',
	textwrap = 40,
	generatio = F,
	xmin = 0,
	foldEn = "% Associated Genes",
	nGene = "Count",
	pval = "FDR" ,
	term = "Term"
)

## ====================== Gene expression plot for Final 6 genes
edata <- read_edata('../mutual-data/merged_edata.tsv')
mdata <- read_mdata('../mutual-data/merged_mdata.tsv')



strip_plot_genes(edata, mdata$condition, ncol=9,gene_names = c(
	'HYOU1', 'NDUFV1', 'POLR2E', 'ECH1', 'BAG6', 
	'GCLM', 'ALDH3A1', 'ETFB', 'TXNRD1'
)) +
	theme(
		legend.position = "top"
	)












