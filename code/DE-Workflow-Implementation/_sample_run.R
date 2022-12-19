dir_refscript='ref_script'
files.sources = list.files(dir_refscript,full.names = T)

sapply(files.sources, source)

# example-1a on processing a BEC method
run_combat(count,cellinfo)
load('combat.rda')
run_wilcox(processed, cellinfo=cellinfo,is.log=T,former.meth = 'combat')
load('combat+wilcox.rda')

# example-1b on processing a BEC method using Python including 'scvi', 'scgen', 'scanorama'
run_format_python(cellinfo, meth='scvi')
run_wilcox(processed, cellinfo=cellinfo,is.log=T,former.meth = 'scvi')
load('scvi+wilcox.rda')

#example-2 on processing a COV method 
run_limmavoom(count,cellinfo,cov=T)
load('limmavoom.rda')

#example-3 on processing a META method
run_LogNormalize(count,cellinfo,separate = T,former.meth = '')
load('LogNormalize_sep.rda')
run_limmatrend_sep(processed=processed,cellinfo=cellinfo,former.meth = 'LogNormalize')
load('LogNormalize_sep+limmatrend_sep.rda')
run_wFisher(res,processed=processed,cellinfo=cellinfo,former.meth='LogNormalize_sep+limmatrend')
load('LogNormalize_sep+limmatrend_sep+wfisher.rda')
