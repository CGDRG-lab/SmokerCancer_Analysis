# ---------------------------------------------------------------------------------
# Author : Ridwan
# Email  : reubenridwan@gmail.com
# License: MIT License (see LICENSE file)
# Date   : 2024
# ---------------------------------------------------------------------------------

library(ggrepel)
library(ggplot2)
library(stringr)
library(dplyr)


## --------------------- read data in from tsv
read_edata <- function(file, index_col='Gene'){
	expr_data <- readr::read_tsv(file = file) %>% 
		tibble::column_to_rownames(index_col) %>% 
		as.matrix()
	return(expr_data)
}

read_mdata <- function(file, index_col='Sample'){
	expr_data <- readr::read_tsv(file = file) %>% 
		tibble::column_to_rownames(index_col)
	return(expr_data)
}

## --------------------------------- find outlier samples
find_outliar_samples <- function(t_edata, groups, f=2){
	outliear_list <- character(0)
	unique_grp <- unique(groups)
	
	for(grp in unique_grp){
		tmp_df <- t_edata[, groups==grp]
		while(TRUE){
			iac_matrix <- cor(tmp_df, method = "pearson")
			avg_iac <- rowMeans(iac_matrix, na.rm = TRUE)
			mean_iac <- mean(avg_iac, na.rm = TRUE)
			sd_iac <- sd(avg_iac, na.rm = TRUE)
			outliers <- names(avg_iac)[avg_iac < (mean_iac - f * sd_iac)]
			if(length(outliers) == 0){ break }
			outliear_list <- c(outliear_list, outliers)
			tmp_df <- tmp_df[, ! colnames(tmp_df) %in% outliers]
		}
	}
	
	return(outliear_list)
}

## ================================================================================
## --------------------------------------- PCA Genes ------------------------------
## ================================================================================
plot_pca <- function(edata, groups, main="PCA of Samples", cex=1, lpos=NA, ...){
	pca <- prcomp(scale(t(edata)))
	e_var <- summary(pca)$importance[2, ]
	
	
	plot(pca$x[, 1], pca$x[, 2], 
		 # col = as.numeric(as.factor(m_mdata_filt$condition)), 
		 pch = "",
		 xlab = paste("PC1 (", round(e_var[1] * 100, 2), "%)", sep = ""),
		 ylab = paste("PC2 (", round(e_var[2] * 100, 2), "%)", sep = ""), 
		 main=main,
		 ...
	)
	
	# Adding the sample names instead of points
	text(pca$x[, 1], pca$x[, 2], 
		 labels = colnames(edata), 
		 cex = cex,  # Adjust text size as needed
		 col = as.numeric(as.factor(groups)))
	
	if(!is.na(lpos)){
		legend(
			lpos, levels(as.factor(groups)), col=1:length(unique(groups)), pch=19,
			cex=cex
		)
	}
}


# ecnrich-plot clusterprofile
lolipop_plot <- function(
	result, title="Dot Plot", xmin=1, textwrap=30, ysize=13, generatio=FALSE,
	foldEn = 'FoldEnrichment', nGene = 'nGene', term='Description', pval='p.adjust'
){
	if(nrow(result) == 0) {
		print("No Row")
		return()
	}
	x_val <- foldEn
	size_val <- nGene
	if(generatio){ x_val <- nGene; size_val <- foldEn }
	
	xmax <- max(result[, x_val]) * 1.1
	
	fig <- result %>% 
		ggplot(aes(x = !!sym(x_val), y = reorder(!!sym(term), !!sym(x_val)), color = !!sym(pval))) +
		geom_segment(aes(x=xmin, xend = !!sym(x_val)), linewidth = 1) +
		geom_point(aes(size = !!sym(size_val)), shape=19, show.legend = FALSE) +
		geom_point(aes(size = !!sym(size_val)), shape=1, color='black') +
		scale_size_continuous(range = c(4, 9)) +
		scale_color_gradientn(colors = c( "tomato", "dodgerblue")) +
		scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax)) +
		scale_y_discrete(labels = function(x) str_wrap(x, width = textwrap)) +
		ggtitle(title) +
		theme_minimal() +
		theme(
			plot.title = element_text(size = 20),
			panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
			axis.title.y = element_blank(),
			axis.text.y = element_text(size = ysize),
			
		)
	
	return(fig)
}


# ----------------------- Genes
get_summary <- function(enrichr, n=10, generatio=FALSE){
	x_val <- 'FoldEnrichment'
	if(generatio){ x_val <- 'nGene'}
	
	enrich_result <- as.data.frame(enrichr) %>% 
		dplyr::mutate(
			nGene = as.numeric(str_replace(GeneRatio, "/.*$", "")),
			pathway_genes = as.numeric(str_replace(BgRatio, "/.*$", "")),
			total_gene = as.numeric(str_replace(GeneRatio, "^.*/", "")),
			background_genes = as.numeric(str_replace(BgRatio, "^.*/", ""))
		) %>%
		dplyr::arrange(desc( !!sym(x_val) ), p.adjust) %>% 
		slice_head(n=n) %>% 
		dplyr::select(Description, p.adjust, FoldEnrichment, nGene, pathway_genes, total_gene, background_genes) %>% 
		tibble::rownames_to_column('ID')
	
	return(enrich_result)
}

## ================================================================================
## -------------------------------- Strip Plot Genes ------------------------------
## ================================================================================
strip_plot_genes <- function(
	edata, condition, gene_names, title="", ncol=5,
	cmap = NULL
){
	sample_mapper <- data.frame(
		Sample_ID = colnames(edata),
		condition = condition
	)
	
	fig <- edata %>%
		as.data.frame() %>% 
		tibble::rownames_to_column("Gene") %>% 
		dplyr::filter(Gene %in% gene_names) %>% 
		tidyr::pivot_longer(
			cols = -'Gene',
			names_to = 'Sample_ID',
			values_to = 'Expression'
		) %>% 
		dplyr::left_join(
			sample_mapper,
			by='Sample_ID'
		) %>% 
		ggplot(
			aes(x = condition, y=Expression, color=condition)
		) +
		geom_boxplot() +
		geom_jitter(width=0.2, alpha=0.6) +
		facet_wrap(
			"Gene", 
			scales = "free_y",
			ncol = ncol
		) +
		ggtitle(title) +
		theme_bw() +
		theme(
			plot.title = element_text(size=16, hjust=0.5, margin = margin(b = 10)),
			axis.title.x = element_blank(),
			axis.text.x = element_text(angle = 45, vjust = 0.5)
		)
	
	if(!is.null(cmap)){
		fig <- fig + 
			scale_color_manual(
				values = cmap
			)
	}
	
	return(fig)
}



## ================================================================================
## ---------------------------------- plot valocano ------------------------------
## ================================================================================
plot_volcano <- function(
		results, minFC=-1, maxFC=1, alpha=0.05, by="P.Value", lab=FALSE,
		highlight = NA, show_line = FALSE, byBvalue = FALSE, a=0.8, 
		fsize = 3, gene_list = NA, show_grid=TRUE, msize=2.5,
		mcolor='black', legend=FALSE, n_overlap = 10
){
	# by = "P.Value"
	# alpha = 0.01
	
	min_fc_value <- min(results$logFC)
	max_fc_value <- max(results$logFC)
	max_log_pval <- max(-log10(results[[by]]))
	
	plot_data <- results %>%
		dplyr::select(Gene, logFC, !!sym(by), B) %>%
		dplyr::mutate(logFC := as.numeric(logFC), !!sym(by) := as.numeric(!!sym(by))) %>% 
		dplyr::mutate(
			state = case_when(
				logFC >= maxFC & !!sym(by) <= alpha ~ 'up',
				logFC <= minFC & !!sym(by) <= alpha ~ 'down',
				TRUE ~ 'same'
			)
		)
	
	figure <- ggplot(plot_data, aes(
		x = logFC, 
		y = -log10( !!sym(by) ), 
		color = state)) +
		geom_hline(yintercept = 0, color = "black", linewidth=0.5)
	
	if(show_line){
		figure <- figure +
			geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "red") +
			geom_vline(xintercept = c(minFC, maxFC), linetype = "dashed", color = "red")
	}
	
	figure <- figure +
		geom_point(show.legend = legend, alpha=a) +
		scale_color_manual(values = c(
			'up'= 'blue',
			'down' = 'red',
			'same' = 'gray'),
			breaks = c("down", "same", "up"),
			labels = c("Down-regulated", "Not-Significant", "Up-regulated")
		)
	
	
	if(!is.na(highlight)){
		if(byBvalue) { plot_data <- arrange(plot_data, desc(B)) }
		else { plot_data <- arrange( plot_data, desc(abs(logFC)), !!sym(by) ) }
		
		figure <- figure + 
			# geom_point(
			# 	data = plot_data[1:highlight, ],
			# 	# color = mcolor,
			# 	# size = msize
			# ) +
			geom_text_repel(
				data = plot_data[1:highlight, ],
				aes(label = Gene),
				size = fsize,
				fontface = 'bold',
				max.overlaps = n_overlap,
				colour = 'black',
				show.legend = NA
			)
	}
	
	if (!all(is.na(gene_list))){
		mark_gene_df <- dplyr::filter(plot_data, Gene %in% gene_list)
		
		figure <- figure + 
			geom_point(
				data = mark_gene_df,
				color = mcolor,
				size = msize
			) +
			geom_text_repel(
				data = mark_gene_df,
				aes(label = Gene),
				size = fsize,
				fontface = 'bold',
				max.overlaps = n_overlap,
				colour = 'black',
				show.legend = FALSE
			)
	}
	
	if (lab){
		figure <- figure + geom_text_repel(
			data = . %>% dplyr::filter(state !=  'same'),
			aes(label = Gene),
			size = fsize,
			fontface = 'bold',
			max.overlaps = n_overlap,
			colour = 'black',
			show.legend = FALSE
		)
	}
	
	
	figure <- figure +
		labs(x = "logFC", y = paste0("-log10(", by, ")", collapse = '')) +
		scale_y_continuous(limits = c(0, NA), breaks = seq(0, ceiling(max_log_pval), 1)) +  # Step of 1 for y-axis
		scale_x_continuous(breaks = seq(floor(min_fc_value), ceiling(max_fc_value), 1)) +  # Step of 1 for x-axis
		theme_bw() +
		theme(
			legend.position = "top",
			legend.title = element_blank(),
			legend.text = element_text(size = 12),
			panel.border = element_rect(color = 'black', fill = NA, linewidth = 1.5)
		)
	
	if(!show_grid){
		figure <- figure +
		theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
		)
	}
	
	return(figure)
}

