# Benchmarking integration of single-cell differential expression
-------------------------------------------------------------------------------------------------------------------------

### General performance of methods
* This analysis demonstrates the high resolution and efficacy of integrative DE analysis for specific cell type as compared to the analysis of bulk sequencing data. Overall, this project covers over 46 integrative methods for scRNA-seq DE analysis and several different levels of batch effects, thus provides a guideline to integrating DE analysis of scRNA-seq data. 
<img src="data/sum_features.png" width="1350"> 
<br>

### Performance of selected methods on specific datasets
* The general comparison of interested methods on 3 datasets using (A) Splatter  and model-free simulations: (B) MCA B-cell and (C) Pancreas data. (D-F) The correspoding evaluation of batch effect on each dataset using principal variance component analysis (PVCA). (G-I) The performance of F-beta score with β=0.5 to emphasize the role of precision. (J-L) The corresponding precision-recall curves of considering methods with pAUPR ranking.

| **Splatter simulation** | **MCA (T cells)** | **Pancreas (Alpha cells)** |
| --- | --- | --- |
|(A) <img src="data/2b_splatter_tsne.png" width="400"> |(B) <img src="data/tcell_tsne.png" width="400"> |(C) <img src="data/pan_tsne.png" width="400"> |
|(D) <img src="data/splatter_80375_pvca.png" width="400"> |(E) <img src="data/tcell_pvca.png" width="400"> |(F) <img src="data/pan_alpha_pvca.png" width="400"> |
|(G) <img src="data/2b_fbeta.png" width="400"> |(H) <img src="data/tcell_fbeta.png" width="400"> |(I) <img src="data/pan_fbeta.png" width="400"> |
|(J) <img src="data/2b_pr.png" width="400"> |(K) <img src="data/tcell_pr.png" width="400"> |(L) <img src="data/pan_pr.png" width="400"> |


### System requirements

All experiments were tested with the following softwares and packages:
- R (>=4.1.2)
- Splatter (1.18.2)
- sva (3.38.0)
- batchelor (1.6.3)
- scMerge (1.6.0)
- limma (3.46.0)
- Seurat (4.0.2)
- Seurat Data (3.0.2) 
- MAST (1.16.0)
- DESeq2 (1.30.1)
- edgeR (3.32.1)
- ZINB-WaVE (1.12.0) 
- ... (dependencies)



All requirement libraries used for testing Python code are listed in the 'requirements.txt' including:
- anndata==0.8.0
- helpers==0.2.0
- matplotlib==3.5.3
- numpy==1.23.1
- pandas==1.4.4
- scanorama==1.7.2
- scanpy==1.9.1
- scgen==2.1.0
- scipy==1.9.1
- scvi==0.6.8
- seaborn==0.12.1
- torch==1.12.1
