library(ggplot2)
require(tidyverse)

fst <- read.table("/home/guarracino/git/chromosome_communities/fst/wc.2000.fst.txt",sep="\t",header=T,comment.char = "?")

# Find the chromosome names and order them: first numerical order, then any non-numerical chromosomes
#   e.g., chr1, chr2, chr22, chrX
chroms <- unique(fst$chromosome)
chrOrder <- sort(chroms)
fst$chrOrder <- factor(fst$chromosome,levels=chrOrder)

fst$pops <- paste(fst$pop1, "vs", fst$pop2)

w <- 2000
stat_type <- 'Weir and Cockerham 1984'
#stat_type<- 'Hudson 1992, Bhatia et al. 2013'

p <- fst %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1000000) %>%
  ggplot(aes(x = chr_position, y = avg_wc_fst, color = chromosome))+
  #facet_wrap(. ~ pops)+
  facet_grid(pops ~ chromosome) +
  geom_point(size = 0.1)+
  xlab("Position on Chromosome (Mb)")+
  ylab("Fst")+
  theme_bw()+
  labs(title=paste0(stat_type, " (window of ", w, " bps)"))+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(0.1, "cm"),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0),
    strip.placement = "outside"
    #legend.position = "none"
  )+
  #scale_x_continuous(expand = c(0, 0))+
  #scale_y_continuous(expand = c(0, 0))+
  xlim(0,20) +
  ylim(-0.2, 1) +
  scale_color_brewer(palette = "Set1")

# Add rDNA annotation
p <- p +
  geom_rect(
    data = data.frame(chromosome = "chm13#chr13"),
    aes(xmin = 5770548/1000000, xmax = 9348041/1000000, ymin = -0.15, ymax = 0.15),
    alpha = 0.5, fill="pink", inherit.aes = FALSE) +
  geom_rect(
    data = data.frame(chromosome = "chm13#chr14"),
    aes(xmin = 2099537/1000000, xmax = 2817811/1000000, ymin = -0.15, ymax = 0.15),
    alpha = 0.5, fill="pink", inherit.aes = FALSE) +
  geom_rect(
    data = data.frame(chromosome = "chm13#chr15"),
    aes(xmin = 2506442/1000000, xmax = 4707485/1000000, ymin = -0.15, ymax = 0.15),
    alpha = 0.5, fill="pink", inherit.aes = FALSE) +
  geom_rect(
    data = data.frame(chromosome = "chm13#chr21"),
    aes(xmin = 3108298/1000000, xmax = 5612715/1000000, ymin = -0.15, ymax = 0.15),
    alpha = 0.5, fill="pink", inherit.aes = FALSE) +
  geom_rect(
    data = data.frame(chromosome = "chm13#chr22"),
    aes(xmin = 4793794/1000000, xmax = 5720650/1000000, ymin = -0.15, ymax = 0.15),
    alpha = 0.5, fill="pink", inherit.aes = FALSE)

# Add centromere annotation
p <- p +
  geom_rect(
    data = data.frame(chromosome = "chm13#chr13"),
    aes(xmin = 13941594/1000000, xmax = 17573031/1000000, ymin = -0.15, ymax = 0.15),
    alpha = 0.5, fill="black", inherit.aes = FALSE) +
  geom_rect(
    data = data.frame(chromosome = "chm13#chr14"),
    aes(xmin = 10067897/1000000, xmax = 12776582/1000000, ymin = -0.15, ymax = 0.15),
    alpha = 0.5, fill="black", inherit.aes = FALSE) +
  geom_rect(
    data = data.frame(chromosome = "chm13#chr15"),
    aes(xmin = 15412039/1000000, xmax = 17709803/1000000, ymin = -0.15, ymax = 0.15),
    alpha = 0.5, fill="black", inherit.aes = FALSE) +
  geom_rect(
    data = data.frame(chromosome = "chm13#chr21"),
    aes(xmin = 10816799/1000000, xmax = 11340698/1000000, ymin = -0.15, ymax = 0.15),
    alpha = 0.5, fill="black", inherit.aes = FALSE) +
  geom_rect(
    data = data.frame(chromosome = "chm13#chr22"),
    aes(xmin = 12784333/1000000, xmax = 16188127/1000000, ymin = -0.15, ymax = 0.15),
    alpha = 0.5, fill="black", inherit.aes = FALSE)

p
