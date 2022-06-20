args = commandArgs(trailingOnly = T)
if (length(args) < 5)
  stop("Must set methods and do.which", call. = F)
methods.use = args[1]
do.which = args[2]
FP = args[3]
seed = as.numeric(args[4])
ct=args[5]
sub.ct=args[6]
base_dir='~/'
source(paste0(base_dir,'script_bk/do.analysis.R'))
setwd(base_dir)
ct%<>%gsub(pattern='_',replacement=' ')
Cell.types<-ct
if(sub.ct=='none'){
  sub.ct=NULL
}else{
  sub.ct%<>%gsub(pattern='_',replacement=' ')
}
for(filter.rate in c(0.01,0.05)){
  for(ct in Cell.types){
    for(meth.comb in methods.use){
      # sub.ct=NULL
      # sub.ct='CD4+ Th'
      # sub.ct='tLung Alveolar Mac vs mo-Mac'
      if(FP=='tumor'){
        filt_dir_origin=paste0('/fp_filtered/')
      }else if(FP=='normal'){
        filt_dir_origin=paste0('/fp_filtered2/')
      }else if(FP=='none'){
        filt_dir_origin='/filtered/'
      }
      if(!is.null(sub.ct)){
        filt_dir=ifelse(is.null(filter.rate),'',paste0(filt_dir_origin,filter.rate,'_',sub.ct))
        count.base=paste0('tLung_state1-nLung.',sub.ct,'_')
      }else{
        filt_dir=ifelse(is.null(filter.rate),'',paste0(filt_dir_origin,filter.rate))
        count.base=paste0('tLung_state1-nLung.')
      }
      if(FP %ni% c('none')){
        filt_dir<-paste0(filt_dir,'/',seed)
      }
      
      if(do.which=='filter'){
        if(FP %in% c('normal','tumor')){
          do.filter.fp(ct,base_dir,filter.rate = filter.rate, filt_dir=filt_dir, count.base=count.base,FP=FP,seed=seed)
        }else if(FP=='none'){
          do.filter(ct,base_dir,filter.rate = filter.rate, filt_dir=filt_dir, count.base=count.base)
        }
      }else if(do.which=='first'){
        do.first(ct,meth.comb,base_dir,filter.rate = filter.rate, filt_dir=filt_dir, count.base=count.base,FP=FP,seed=seed, vs="original")
      }else if(do.which=='second'){
        do.second(ct,meth.comb,base_dir,filter.rate=filter.rate, filt_dir=filt_dir, count.base=count.base,FP=FP,seed=seed)
      }else if(do.which=='third'){
        do.third(ct,meth.comb,base_dir,filter.rate=filter.rate, filt_dir=filt_dir, count.base=count.base,FP=FP,seed=seed)
      }
      
      
    }
  }
}