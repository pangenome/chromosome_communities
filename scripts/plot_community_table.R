library(ggplot2)
library(tidyverse)

# install ggrepel from git
#install.packages("devtools")
#devtools::install_github("slowkow/ggrepel")


#path_community_table <- 'HPRCy1v2genbank.self.s50k.l250k.p95.n93.h0001.l1000000.paf.community.leiden.tsv'
#path_community_2_size <- 'HPRCy1v2genbank.self.s50k.l250k.p95.n93.h0001.l1000000.paf.community2size.tsv'
#total_sequence_content_bp <- 2.83434e+11
#path_output <- 'x.pdf'

# read in data
x <- read.table(path_community_table, header = T) %>%
  rename(`Unassigned` = not.partitioned)

x$community.of <- gsub("not.partitioned_", "Unassigned - ", x$community.of)
x$community.of <- gsub("_", " - ", x$community.of)

# Respect and reverted the order in the file
mylevels <- x$community.of
x$community.of <- factor(x$community.of, levels=rev(mylevels))

y <- read.table(path_community_2_size, header = F)
colnames(y) <- c('num.community', 'variable', 'sequence.content.bp')
y$variable <- as.character(y$variable)
y[y$variable == 'unmapped',]$variable <- "Unassigned"
y$variable <- as.factor(y$variable)

library(reshape2)
xy <- merge(
  melt(x, id=c("num.community" ,"community.of")),
  y,
  by = c('num.community', 'variable')
)

summary(xy)

#xy$value[xy$value > 100] <- ''

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
    #plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 21),
    axis.text.x = element_text(size = 19, angle = 90, vjust = 0.5, hjust=1),
    axis.text.y = element_text(size = 19),
    
    legend.title = element_text(size = 19),
    legend.text = element_text(size = 19),
    legend.position = "right",
    
    #axis.title.y=element_blank()
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(margin = margin(t = -20)),
    axis.title.y = element_text(margin = margin(r = -20))
  ) +
  labs(x = "Chromosome", y = 'Community', fill = '% seq') + 
  #theme(aspect.ratio=1) + 
  theme(
    plot.margin = unit(c(0,0,0,0), "cm")
  )
p <- p +
  #geom_text(aes(label=round(sequence.content.bp / total_sequence_content_bp, digits = 10)), size=2.5, color="red") +
  # angled text labels
  with(xy, geom_text(aes(label=value),
                     size=6.5, color="black",
                     # the angle should depend on if we're in the last column to the right
                     angle = ifelse(variable == "Unassigned", 0, 45),
                     hjust=0.5, vjust=0.5))
ggsave(
  plot = p,
  path_output,
  width = 34, height = 30, units = "cm", dpi = 300, bg = "transparent", limitsize = F)

# OLD CODE
#x$community.of <- factor(x$community.of, levels=unique(x[order(as.integer(gsub("[^0-9]", "", x$community.of))),'community.of']))
#co <- melt(x%>% select(-num.community), id.vars = 'community.of')
#num_contigs <- sum(co$value)

