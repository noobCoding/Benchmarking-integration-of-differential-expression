methodsuse=('combat++findmarkers' 'limma_bec++findmarkers' 'mnn_opt++findmarkers' 'scMerge++findmarkers' "Seurat_++findmarkers"  'zinbwave++findmarkers'  'raw++findmarkers' 'raw++MAST' 'raw++MAST_Cov' 'raw++DESeq2' 'raw++DESeq2_Cov' 'zinbwave++DESeq2_pseudo' 'zinbwave++DESeq2_pseudo_Cov' "raw+ind.DESeq2.pval+Fisher" 'raw+ind.DESeq2.pval+wFisher'  'raw++edgeR_DetRate' 'raw++edgeR_DetRate_Cov' 'raw++edgeR' 'raw++edgeR_Cov' 'zinbwave++edgeR' 'zinbwave++edgeR_Cov' "raw+ind.edgeR.pval+Fisher" 'raw+ind.edgeR.pval+wFisher' "raw++limma" "raw++limma_Cov" 'raw++limma_trend' 'raw++limma_trend_Cov' 'combat++limma_trend' 'mnn_opt++limma_trend' 'scMerge++limma_trend' 'LogNormalize+ind.limma_trend+Fisher' 'LogNormalize+ind.limma_trend+wFisher' 'raw+ind.DESeq2.ES+FEM' "voom+ind.ES+FEM" "LogNormalize+ind.ES+FEM" "raw+ind.DESeq2.ES+REM"  "voom+ind.ES+REM" "LogNormalize+ind.ES+REM" 'voom+ind.modt+Fisher'  'voom+ind.modt+wFisher')

seeds=(1234)
dowhich=$1
FP=$2
ct=$3
subct=$4
for i in "${methodsuse[@]}"
do
   for seed in ${seeds[@]}
   do
       cd "~/script_bk"
       id=${i}
       echo ${id}
       echo ${seed}
       Rscript run_integrated_analysis_terminal.R ${id} ${dowhich} ${FP} ${seed} ${ct} ${subct}
       rm -rf ~/.local/share/Trash/*
    done
done
