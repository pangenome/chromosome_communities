args <- commandArgs()
path_community_table <- '/home/guarracino/Downloads/Pangenomics/chromosome_communities/paper/HPRCy1v2genbank.self.s20k.l100k.p98.n93.h0001.l1000000.paf.edges.weights.txt.community.leiden.tsv'

library(ggplot2)

x <- read.table(path_community_table, header = T)

mylevels <- x$community.of
x$community.of <- factor(x$community.of, levels=mylevels)

#x$community.of <- factor(x$community.of, levels=unique(x[order(as.integer(gsub("[^0-9]", "", x$community.of))),'community.of']))

library(reshape2)
co <- melt(x, id.vars = 'community.of')

num_contigs <- sum(co$value)


p <- ggplot(
  co,
  aes(
    x = community.of,
    y = variable,
    fill = value / num_contigs
  )
) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1),
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "right",
    
    #axis.title.y=element_blank()
  ) +
  labs(x = "Community", y = 'Chromosome', fill = '% of contigs') + theme(aspect.ratio=1)

p
