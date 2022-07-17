args <- commandArgs()
path_community_table <- args[6]
path_community_2_size <- args[7]
total_sequence_content_bp <- as.numeric(args[8])
path_output <- args[9]

library(ggplot2)
library(tidyverse)

x <- read.table(path_community_table, header = T)

# Respect and reverted the order in the file
mylevels <- x$community.of
x$community.of <- factor(x$community.of, levels=rev(mylevels))

y <- read.table(path_community_2_size, header = F)
colnames(y) <- c('num.community', 'variable', 'sequence.content.bp')
y$variable <- as.character(y$variable)
y[y$variable == 'unmapped',]$variable <- "Not partitioned"

library(reshape2)
xy <- merge(
  melt(x, id=c("num.community" ,"community.of")),
  y,
  by = c('num.community', 'variable')
) 

# Info grouped by community
#xy %>% group_by(num.community) %>% summarise(num.contigs.comm = sum(value)) %>% View()
#xy %>% group_by(num.community) %>% summarise(sequence.content.comm = sum(sequence.content.bp)) %>% View()

#num_contigs <- sum(xy$value)
#xy$ratio <- xy$sequence.content.bp / total_sequence_content_bp * 100

p <- ggplot(
  xy,
  aes(
    x = variable,
    y = community.of,
    fill = sequence.content.bp / total_sequence_content_bp * 100.0
  )
) +
  geom_tile() + 
  scale_fill_gradient(low = "#FAFAFA", high = "red") +
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
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  ) +
  labs(x = "Chromosome", y = 'Community', fill = '% of sequence') + 
  theme(aspect.ratio=1)
#ggsave(
#  plot = p,
#  path_output,
#  width = 25, height = 25, units = "cm", dpi = 300, bg = "transparent", limitsize = F)

p <- p +
  #geom_text(aes(label=round(sequence.content.bp / total_sequence_content_bp, digits = 10)), size=2.5, color="red") + 
  geom_text(aes(label=value), size=2.6, color="black")
ggsave(
  plot = p,
  paste0(path_output),
  width = 25, height = 25, units = "cm", dpi = 300, bg = "transparent", limitsize = F)

# OLD CODE
#x$community.of <- factor(x$community.of, levels=unique(x[order(as.integer(gsub("[^0-9]", "", x$community.of))),'community.of']))
#co <- melt(x%>% select(-num.community), id.vars = 'community.of')
#num_contigs <- sum(co$value)
