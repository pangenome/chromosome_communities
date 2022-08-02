args <- commandArgs()
path_support_dedup_tsv <- args[6]
path_entropy_tsv <- args[7]
dir_annotation <- args[8]
path_output <- args[9]

x_min <- 0
x_max <- 25000000
width <- 90
height <- 20


# library
library(ggridges)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(scales) # for pretty_breaks()

library(png)
library(grid)
library(ggpubr)


# Collapsed untangled output
s <- read.delim(path_support_dedup_tsv)
colnames(s) <- c("ground.target", "start", "end", "chr13", "chr14", "chr15", "chr21", "chr22")


# Entropy
e <- read.delim(path_entropy_tsv)
e[e$shannon_div_index == -1, ]$shannon_div_index <- NA

# Compute average SDI by window by ground.target
e_average <- e %>% 
  dplyr::group_by(ground.target, start.pos, end.pos) %>% 
  dplyr::summarize(average_sdi = mean(shannon_div_index, na.rm = TRUE))


p_panels <- c()

for (num_chr in c(13, 14, 15, 21, 22)) {
  print(num_chr)

  if (num_chr == 13) {
    colors <- c("#F8766D")
  } else if (num_chr == 14) {
    colors <- c("#A3A500")
  } else if (num_chr == 15) {
    colors <- c("#00BF7D")
  } else if (num_chr == 21) {
    colors <- c("#00B0F6")
  } else if (num_chr == 22) {
    colors <- c("#E76BF4")
  }else {
    colors <- c("#000000")
  }
  
  # Annotation
  path_annotation <- file.path(dir_annotation, paste0('hgt_genome_euro_chr', num_chr, '_0_25Mbp.png'))
  img <- readPNG(path_annotation)
  p_annotation <- ggplot() +
    annotation_custom(
      rasterGrob(img, width = 1, height = 1),
      xmin = - Inf, xmax = Inf,
      ymin = - Inf, ymax = Inf
    ) + theme(
      plot.margin = unit(c(0,1,0.5,0), "cm")
    )
  

  chr <- paste0('chm13#chr', num_chr)
  
  # Apply filters
  s_chr <- s[s$ground.target == chr, ]
  s_chr_long <- pivot_longer(s_chr, chr13:chr22, "chromosome")
  
  p_collapsed_untangle <- ggplot(s_chr_long, aes(x=start, y=value, color=chromosome)) +
    geom_step(alpha=0.5) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0, 0, 0)) +
    scale_y_continuous(limits = c(0, max(10, max(s_chr_long$value))), breaks=pretty_breaks(n=6)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      
      text = element_text(size = 32),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "top",
      
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      plot.margin = unit(c(0,1.03,0,6.74), "cm")
    ) + labs(
      x = paste(''),
      y = paste('# contigs\n')
    ) +
    guides(
      colour = guide_legend(title="Chromosome", override.aes = list(size=10))
    )
  # Gray-out regions with missing aligned contigs
  p_collapsed_untangle <- p_collapsed_untangle +
    annotate("rect",
             xmin = 0, xmax = s_chr[1,]$end,
             ymin = 0, ymax = max(10, max(s_chr_long$value)),
             fill = "#444444", alpha = .1, color = "#ffffff", size = 0.1)
  
  # Apply filters
  e_chr_average <- e_average[e_average$ground.target == chr,]
  
  p_entropy_average <- ggplot(e_chr_average, aes(x=start.pos, y=average_sdi, color=ground.target)) +
    geom_step(alpha=0.8) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0, 0, 0)) +
    scale_y_continuous(limits = c(0, min(2, max(e_chr_average$average_sdi, na.rm = T))), breaks=pretty_breaks(n=6)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      
      text = element_text(size = 32),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "top",
      
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      plot.margin = unit(c(0,1.03,0,6.7), "cm")
    ) + labs(
      x = paste('Position'),
      y = paste('Entropy\n')
    ) +
    scale_color_manual(values=colors) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  # Gray-out regions with missing aligned contigs
  p_entropy_average <- p_entropy_average +
    annotate("rect",
             xmin = 0, xmax = e_chr_average[!is.na(e_chr_average$average_sdi),][1,]$start.pos,
             ymin = 0, ymax = min(2, max(e_chr_average$average_sdi, na.rm = T)),
             fill = "#444444", alpha = .1, color = "#ffffff", size = 0.1)
  
  # Final panel
  p_panel <- ggpubr::ggarrange(
    p_annotation, p_collapsed_untangle, p_entropy_average,
    labels=c('', '', ''), font.label = list(size = 40, color = "black", face = "bold", family = NULL),
    heights = c(1.2, 1, 1),
    legend = "right", # legend position,
    common.legend = T,
    nrow = 3
  )

  #ggsave(
  #  plot = p_panel,
  #  file.path(path_output, paste0(num_chr, '.png')),
  #  width = width, height = height,
  #  units = "cm",
  #  dpi = 100, bg = "white",
  #  limitsize = FALSE
  #)
  
  p_panels[[length(p_panels)+1]] <- p_panel
}


# Put panels together
p_figure <- ggpubr::ggarrange(
  p_panels[[1]], NULL, p_panels[[2]], NULL, p_panels[[3]], NULL, p_panels[[4]], NULL, p_panels[[5]],
  #labels=c('A', 'B', 'C', 'D', 'E'), font.label = list(size = 40, color = "black", face = "bold", family = NULL),
  heights = c(1,0.05,1,0.05,1,0.05,1,0.05,1),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 9, ncol = 1
)

ggsave(
  plot = p_figure,
  path_output,
  width = width, height = height*5 + 4*height/20, # There are 4 NULL plots that have height == height/20
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)

