library(dplyr)
library(tidyr)
library(ggplot2)
library(bracer)
library(magrittr)
library(patchwork)


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <filePath>  <outputPlot>', call.=FALSE)
}

# concat all ld files and add chrom (chromosome) and type (p-arm, q-arm, recombinant) columns 
new = data.frame()
fileList=glob(args[1], engine = "r")
for (f in fileList){
  myd<-read.table(f, header = T, comment = '')
  chromosome = strsplit(basename(f), ".", fixed = T) %>% sapply(extract2, 1)
  length = strsplit(basename(f), ".", fixed = T) %>% sapply(extract2, 2)
  all<-myd %>% mutate(chrom = paste(chromosome), type = paste(length))
  new <-rbind(new,all)
}

#new %>% write.table('ciro.txt')
colpalette<-c('#C499BA', '#FF7396', '#0096FF' ) #'#EF5B0C') # '#FFB562')

p1 <- new %>% group_by(type, CHR) %>%  ggplot(aes(F_MISS, fill=type) )+ geom_histogram(aes(y=..density..), alpha=0.5, binwidth=0.01) +theme_light() + xlab('Fraction of missing genotypes') + ggtitle('Missing genotypes per site')+facet_wrap (CHR ~ .,nrow=5,  scales='free_y') + scale_fill_manual(values =colpalette )+ theme(legend.position='top', legend.title=element_blank())  


#allele frequencies file 
myfreq=read.table(args[2], header=T , comment='' )
#myphr %>% write.table('luigi.txt')

#bedfile of putative homologous regions 
myphr=read.table(args[3], header=F, comment='' )

#  chromosome position 
mypos=read.table(args[4], header=T, comment='' ) 



### PLOTS 

## AFS 
p2<- myfreq %>% filter(MAF>0)%>% ggplot(aes(MAF))+geom_histogram(  position='dodge',  binwidth=0.01)  +facet_wrap(CHR ~ . , nrow=5 )+theme_light()+ggtitle('Allele frequency spectrum')

# Size of the Pseudo Homologous Regions
p3<- myphr %>% ggplot (aes(V1, (V3-V2)/1000) )+ geom_violin()  + ylab('Segment size (kb)' ) +xlab('')+ ggtitle ('Size of the Pseudo Homologous Regions' )+ theme_light() + coord_flip() 

#Distance between consecutive markers
p4<- mypos %>% group_by(X.CHROM) %>% mutate (Distance = POS - lag(POS)) %>% ggplot(aes(X.CHROM, Distance))+ geom_violin() + scale_y_log10() + coord_flip() +theme_light() +ggtitle('Distance between consecutive markers') +ylab('Distance (bp)')+ xlab('')
#
p1 |  p2 | p3/p4 + plot_annotation(tag_levels ='A') 
 
ggsave(args[5], width = 20, heigh = 10)



