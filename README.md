# SmokerCancer_Analysis

Comparative transcriptomic analysis of smokers with and without lung cancer using GEO datasets, WGCNA, and machine learning.

**Objective**:  
To investigate potential protective molecular mechanisms in long-term smokers who do not develop lung cancer, despite similar environmental exposure as cancer patients.

**Key Approaches**:
- Microarray preprocessing and merging from NCBI GEO
- Differential gene expression analysis with meta-analysis (Stouffer’s method)
- Weighted gene co-expression network analysis (WGCNA)
- Protein-protein interaction network construction
- Gene enrichment analysis
- Machine learning-based feature selection (PLS-DA, Genetic Algorithm, Logistic Regression)


## 🔬 Methods Summary

Publicly available microarray datasets (GSE and E-MTAB series) were downloaded and independently preprocessed. After log-transformation and probe-to-gene collapsing, a meta-analysis was conducted using Stouffer's method to identify consensus differentially expressed genes (DEGs). WGCNA was applied to identify co-expression modules, followed by STRING-based PPI network construction. Enrichment was performed via DAVID. Finally, multiple machine learning methods were used to prioritize informative genes.


## 📈 Key Findings

Smokers without lung cancer exhibited higher and more coordinated expression of NRF2-related detoxification and redox regulatory genes (e.g., **CYP1A1**, **CYP1B1**, **ALDH3A1**, **TXNRD1**). This suggests a potential role for oxidative stress buffering in protecting against tumorigenesis in some smokers.


## 📊 Data

- Fully preprocessed expression and metadata files, along with differential gene expression results, are available in the `mutual-data/` directory.
- The meta-analysis Stouffer’s combined results can be found in the `dge/` directory.
- Weighted Gene Co-Expression Network Analysis (WGCNA) outputs are stored in the `wgcna/` directory.
- Protein-Protein Interaction (PPI) network files and related data can be found in the `network/` directory.

## 🧪 Tools & Libraries Used

- R (limma, metaMA, WGCNA, clusterProfiler)
- Python (scikit-learn, pandas, numpy)
- Cytoscape & STRINGdb
- DAVID Functional Annotation Tool

## 📜 License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.

<!--
## 📣 Citation

If you use this code or analysis in your research, please cite:
-->

## 📬 Contact

For questions, please reach out via [reubenridwan@gmail.com](mailto:reubenridwan@gmail.com)


## 🙏 Acknowledgements

Special thanks to my thesis supervisor Dr. Mustak Ibn Ayub (DPhil, University of Oxford), and to the GEB Department at University of Dhaka.
