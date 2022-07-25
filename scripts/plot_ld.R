library(tidyverse)
library(bracer)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <filePath>  <outputPlot>', call.=FALSE)
}

new = data.frame()
fileList=glob(args[1], engine = "r")
for (f in fileList){
myd<-read.table(f, header = T, comment = '')
chromosome = strsplit(basename(f), ".", fixed = T) %>% sapply(extract2, 1)
length = strsplit(basename(f), ".", fixed = T) %>% sapply(extract2, 2)
all<-myd %>% mutate(chrom = paste(chromosome), type = paste(length))
new <-rbind(new,all)
}

new %>% mutate (ranges=cut(BP_B-BP_A, seq(min(BP_B-BP_A), max(BP_B-BP_A), 20 ) ) ) %>% group_by(ranges, type, chrom ) %>% summarize( d_count= length (R2), d_stat=mean(R2), d_sd=sd(R2) , ci= 1.96*(d_sd/sqrt(d_count)), upper= d_stat+ci ,  lower= d_stat-ci ) %>% filter(!is.na(ranges)) %>% separate(ranges, into = c('start', 'end'), sep = "," ,remove = FALSE, convert = TRUE) %>% mutate(realStart = gsub("\\(", "", start)) %>% mutate(realEnd = gsub("\\]", "", end)) %>% select(-start, -end) %>% mutate(start = as.numeric(as.character(realStart)), end = as.numeric(as.character(realEnd))) %>% select(-realStart, -realEnd) %>% mutate(midpoint = (start+end)/2) %>% filter(midpoint < 1000) %>% mutate(end = as.factor(end)) %>% ggplot(aes(end, d_stat, color = type, size = d_count)) + geom_point() + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + scale_x_discrete(guide = guide_axis(angle = 90)) + xlab("ranges") + facet_wrap( .~ chrom, nrow = 5) + theme_bw() + theme(legend.text=element_text(size=14), axis.text.x=element_text(size=14), axis.text.y = element_text(size = 14) ,axis.title=element_text(size = 14), strip.text = element_text(size = 14)) + scale_y_continuous(0,1.0)
ggsave(args[2], width = 20, heigh = 20)
