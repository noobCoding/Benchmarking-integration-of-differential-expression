#raw count as input count=rawcount
#Python methods include 'scvi', 'scgen', 'scanorama'

run_format_python<-function(count=NULL, cellinfo=NULL, meth='scvi', former.meth=''){
  
  
  # System run Python method
  python_corrected_output<- read.csv(file = paste0( meth , '_corrected_data.csv'),header=T,row.names=1,check.names = F)
  
  res<-python_corrected_output
  processed<-python_corrected_output
  save(res,processed,cellinfo,file=paste0('./',ifelse(former.meth=='','',paste0(former.meth,'+')), meth, '.rda'))
}