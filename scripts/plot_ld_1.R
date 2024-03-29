library(dplyr)
library(tidyr)
library(ggplot2)
library(bracer)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify four arguments: <filePath> <size of binning in bp> <max size to plot in bp> <outputPlot>', call.=FALSE)
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

# calculate r2 statistics in binsizebp windows and plot
# x = ranges = distance
# y =  R2 statistics 
# d_count = number of observationin a bin 
# type = Type
# recombinant = PHR

binsize=as.numeric(args[2]) #in base pairs 
endSize=as.numeric(args[3]) #in base pairs 

colpalette<-c('#C499BA', '#0096FF', '#FF7396') #'#EF5B0C') # '#FFB562')

new %>%
  filter(R2>0) %>%  
  mutate(kb_dist=BP_B-BP_A) %>% 
  mutate(ranges=cut(kb_dist, seq(min(kb_dist), max(kb_dist), binsize ) ) ) %>%
  group_by(ranges, type, chrom ) %>%
  summarize( d_count= length (R2), d_stat=mean(R2), d_sd=sd(R2) , ci= 1.96*(d_sd/sqrt(d_count)), upper= d_stat+ci ,  lower= d_stat-ci ) %>% 
  filter(!is.na(ranges)) %>% 
  separate(ranges, into = c('start', 'end'), sep = "," ,remove = FALSE, convert = TRUE) %>%
  mutate(realStart = as.numeric(gsub('\\(', '', start))) %>% 
  mutate(realEnd = as.numeric(gsub('\\]', '', end))) %>% 
  select(-start, -end) %>% 
  filter(realEnd <= endSize) %>% 
  mutate(tick=paste(round(realStart/1000,2), '-', round(realEnd/1000,2), sep='')) %>% 
  #mutate(midpoint = (realStart+realEnd)/2) %>% 
  #filter(midpoint  < midpointsize) %>% 
  #mutate(f_mid=midpoint/1000) %>% 
  #ggplot(aes(as.factor(f_mid), d_count, color = type, size = d_count)) +
  ggplot(aes(as.factor(tick), d_stat, color = type, size = d_count)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + 
    xlab('Distance between markers (kb) ') + ylab('Average r^2') +
    scale_colour_manual(values =colpalette ) +
    facet_wrap( .~ chrom, nrow = 5) + 
    theme_light() + 
    theme(
      legend.text=element_text(size=16),
      axis.text.x=element_text(size=12, angle = 90),
      axis.text.y = element_text(size = 14),
      axis.title=element_text(size = 14), strip.text = element_text(size = 14), 
      legend.title=element_blank()) +
    ylim(0,1) 
ggsave(args[4], width = 12, heigh = 14)

