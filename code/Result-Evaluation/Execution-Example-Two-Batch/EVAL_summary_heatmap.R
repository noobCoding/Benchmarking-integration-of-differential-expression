library(magrittr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stringr)
select_color = function(method_color)
{
  answer = switch(method_color,
                  # "none"= '#FFFFFF',
                  # "high"= '#FFC000',
                  # 'medium'='#37CBFF',
                  # 'low'= '#0070C0'
                  "0"= '#FFFFFF',
                  "Top good"= '#FFC000',
                  "Good"='#37CBFF',
                  "Poor"= '#0070C0'
  )
  return(answer)
}

heatmap <- readRDS("heatmapdata.RDS")
heatmap[heatmap=='1'] <- 'Top good'
heatmap[heatmap=='2'] <- 'Good'
heatmap[heatmap=='3'] <- 'Poor'

Methods=colnames(heatmap)
Methods%>%length()
y_order<-colnames(heatmap)
x_order<-rownames(heatmap)
dff <- as.data.frame(t(heatmap))

dff$Methods=Methods

dff<-melt(dff,id.vars =c("Methods"),  measure.vars= setdiff(colnames(dff),c('Methods')))
dff$Methods%<>%factor(levels=y_order)
dff$variable%<>%factor(levels=rev(x_order))
dff$value <- as.character(dff$value)
colnames(dff) <- c("Methods", "Features", "value")

p<-ggplot(dff, mapping=aes(x=Methods, y=Features))+
  geom_tile(aes(fill=value),colour='white', lwd=0.25, linetype=1)+
  scale_fill_manual(name='Performance', values=sapply(X=as.character(c('0','Top good','Good','Poor')),FUN=select_color))+
  theme(axis.text.x = element_text(angle = 45, hjust=0, size=13),
        axis.text.y = element_text(size=13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
        )+ scale_x_discrete(position = "top")
p  

pdf("summary_heatmap.pdf", width=20, height=4)
p
dev.off()
