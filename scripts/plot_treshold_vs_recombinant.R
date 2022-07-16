library(ggplot2)


x <- read.delim('/home/guarracino/git/chromosome_communities/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.recombinant_regions.table.tsv', header = F)
x$V3 <- x$V3/1000/1000 # Convert in Mbps

colnames(x) <- c('SelfCovThreshold', 'Threshold', 'Mbps')

x$SelfCovThreshold <- as.factor(x$SelfCovThreshold)
x$SelfCovThreshold <- factor(x$SelfCovThreshold, levels = c("0", "1.5", "1"))
options(scipen=10000) # Disable scientific notation on axes

ggplot(x, aes(x=Threshold, y=Mbps)) + 
  geom_line() + 
  theme_bw() + 
  facet_grid(~SelfCovThreshold) +
  scale_x_continuous(breaks = round(seq(min(x$Threshold), max(x$Threshold), by = 0.005), 3)) +
  scale_y_continuous(breaks = round(seq(0, max(x$Mbps), by = 0.5), 3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  