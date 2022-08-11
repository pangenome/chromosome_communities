
require(tidyverse)

x <- read.delim('bisnps.tsv.gz')
x$graph <- as.factor(x$graph)
x$chrom <- as.factor(x$chrom)
summary(x$graph)
for (s in levels(x$chrom)) { print(s); ggplot(subset(x, chrom==s), aes(x=pos, fill=as.factor(lv))) + geom_histogram(binwidth=100000) + scale_fill_manual("LV", values=c('#888888','#d53e4f','#fc8d59','#fee08b','#99d594','#3288bd','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8','#3713d8')) + ggtitle(paste("biallelic SNPs",s)) + facet_grid(graph ~ .) ; ggsave(paste("bisnps.",s,".dens.lv_by.pdf",sep=""), height=5, width=15); }