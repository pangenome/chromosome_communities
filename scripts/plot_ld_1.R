library(dplyr)
library(tidyr)
library(ggplot2)
library(bracer)
library(magrittr)

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

# calculate r2 statistics in binsizebp windows and plot
# x = ranges = distance
# y =  R2 statistics 
# d_count = number of observationin a bin 
# type = Type
# recombinant = PHR

binsize=500 #in base pairs 
midpointsize=5000  # in base pairs 

colpalette<-c('#C499BA', '#FF7396', '#0096FF' ) #'#EF5B0C') # '#FFB562')

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
  mutate(midpoint = (realStart+realEnd)/2) %>% 
  filter(midpoint  < midpointsize) %>% 
  ggplot(aes(as.factor(midpoint), d_count, color = type, size = d_count)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + 
    xlab('Distances (kb) ') + ylab('Average r2') +
    scale_colour_manual(values =colpalette ) +
    facet_wrap( .~ chrom, nrow = 5) + 
    theme_light() + 
    theme(
      legend.text=element_text(size=16),
      axis.text.x=element_text(size=14),
      axis.text.y = element_text(size = 14),
      axis.title=element_text(size = 14), strip.text = element_text(size = 14), 
      legend.title=element_blank()) +
    ylim(0,1) 
ggsave(args[2], width = 16, heigh = 20)

