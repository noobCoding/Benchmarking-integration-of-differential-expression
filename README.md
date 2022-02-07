# Benchmarking integration of single-cell RNA-seq differential analysis
#### Hai C. T. Nguyen†, Bukyung Baik†, Sora Yoon, Hae-Ock Lee, Taesung Park, Dougu Nam*
-------------------------------------------------------------------------------------------------------------------------

### General performance of methods
* This analysis demonstrates the high resolution and efficacy of integrative DE analysis for specific cell type as compared to the analysis of bulk sequencing data. Overall, this project covers over 40 integrative methods for scRNA-seq DE analysis and several different levels of batch effects, thus provides a guideline to integrating DE analysis of scRNA-seq data. 
<img src="data/summary_2.png" width="1350"> 

### Performance of selected methods on specific datasets
* The general comparison of interested methods on 3 datasets using (A) Splatter  and model-free simulations: (B) MCA B-cell and (C) Pancreas data. (D-F) The correspoding evaluation of batch effect on each dataset using principal variance component analysis (PVCA). (G-I) The performance of F-beta score with β=0.5 to emphasize the role of precision. (J-L) The corresponding precision-recall curves of considering methods with pAUPR ranking.

| **Splatter simulation** | **MCA (T cells)** | **Pancreas (Alpha cells)** |
| --- | --- | --- |
|(A) <img src="data/splatter_dist.png" width="400"> |(B) <img src="data/tcell_dist.png" width="400"> |(C) <img src="data/pan_dist.png" width="400"> |
|(D) <img src="data/pvca_splatter.png" width="400"> |(E) <img src="data/pvca_tcell.png" width="400"> |(F) <img src="data/pvca_pan.png" width="400"> |
|(G) <img src="data/8095_37_fscore.png" width="400"> |(H) <img src="data/tcell_95_fscore_zbw_psc.png" width="400"> |(I) <img src="data/pan95_fscore.png" width="400"> |
|(J) <img src="data/80375_RCC_50.png" width="400"> |(K) <img src="data/tcell_RCC_50.png" width="400"> |(L) <img src="data/pan_alpha_RCC_50.png" width="400"> |

*To whom correspondence should be addressed. Tel: +82-52-217-2525; Fax: +82-52-217-2639; Email: dougnam@unist.ac.kr
|-------------------------------------------------------------------------------------------------------------------------|
