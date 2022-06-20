methodsuse=('combat++findmarkers' 'limma_bec++findmarkers' 'mnn_opt++findmarkers' 'scMerge++findmarkers' "Seurat_++findmarkers" 'zinbwave++findmarkers' 'raw++findmarkers' "raw+ind.DESeq2.pval+Fisher" "raw+ind.edgeR.pval+Fisher" 'LogNormalize+ind.limma_trend+Fisher' 'raw+ind.DESeq2.ES+REM' "voom+ind.ES+REM" "LogNormalize+ind.ES+REM" 'voom+ind.modt+Fisher')

seeds=(1234)
dowhich=$1
FP=$2
ct=$3
subct=$4
#seed=$3
#seed=1234
#FP='none'
for i in "${methodsuse[@]}"
do
   for seed in ${seeds[@]}
   do
       cd "/hdd2/SC lung/script"
       id=${i}
       echo ${id}
       echo ${seed}
       Rscript run_integrated_analysis_terminal.R ${id} ${dowhich} ${FP} ${seed} ${ct} ${subct}
       rm -rf ~/.local/share/Trash/*
    done
done
