
methodsuse=('combat++findmarkers' 'limma_bec++findmarkers' 'mnn_opt++findmarkers' 'scMerge++findmarkers' "Seurat_++findmarkers" 'zinbwave++findmarkers' 'raw++findmarkers' "raw+ind.DESeq2.pval+Fisher" 'LogNormalize+ind.limma_trend+Fisher' "voom+ind.ES+REM")

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
