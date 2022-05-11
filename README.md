# Benchmarking integration of single-cell RNA-seq differential analysis
-------------------------------------------------------------------------------------------------------------------------

### General performance of methods
* This analysis demonstrates the high resolution and efficacy of integrative DE analysis for specific cell type as compared to the analysis of bulk sequencing data. Overall, this project covers over 40 integrative methods for scRNA-seq DE analysis and several different levels of batch effects, thus provides a guideline to integrating DE analysis of scRNA-seq data. 
<img src="data/summary_220407.png" width="1350"> 
<br>
### Performance of selected methods on specific datasets
* The general comparison of interested methods on 3 datasets using (A) Splatter  and model-free simulations: (B) MCA B-cell and (C) Pancreas data. (D-F) The correspoding evaluation of batch effect on each dataset using principal variance component analysis (PVCA). (G-I) The performance of F-beta score with β=0.5 to emphasize the role of precision. (J-L) The corresponding precision-recall curves of considering methods with pAUPR ranking.

| **Splatter simulation** | **MCA (T cells)** | **Pancreas (Alpha cells)** |
| --- | --- | --- |
|(A) <img src="data/2b_splatter_tsne.png" width="400"> |(B) <img src="data/tcell_tsne.png" width="400"> |(C) <img src="data/pan_tsne.png" width="400"> |
|(D) <img src="data/splatter_80375_pvca.png" width="400"> |(E) <img src="data/tcell_pvca.png" width="400"> |(F) <img src="data/pan_alpha_pvca.png" width="400"> |
|(G) <img src="data/sp80_fbeta.png" width="400"> |(H) <img src="data/tcell_fbeta.png" width="400"> |(I) <img src="data/pan_fbeta.png" width="400"> |
|(J) <img src="data/sp80_PR.png" width="400"> |(K) <img src="data/tcell_PR.png" width="400"> |(L) <img src="data/pan_PR.png" width="400"> |

*To whom correspondence should be addressed. Tel: +82-52-217-2525; Fax: +82-52-217-2639; Email: dougnam@unist.ac.kr
|-------------------------------------------------------------------------------------------------------------------------|
