
args <- commandArgs()
path_msmc2_final_txt <- args[6]
population_name <- args[7]
title <- args[8]

library(ggplot2)
library(ggforce)

mu <- 1.25e-8
gen <- 30
x <- read.table(path_msmc2_final_txt, sep = '\t', header = T)


x$pop <- population_name

x$x <- x$left_time_boundary/mu*gen
x$y <- (1/x$lambda)/(2*mu)


p <- ggplot(x, aes(x = x, y = y, color = pop)) + 
  geom_step() + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 20),#, family="NimbusSan"),
    axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),#, family="NimbusSan"),
    axis.text.y = element_text( size = 16),#, family="NimbusSan"),
    
    legend.title = element_text(size=20),#, family="NimbusSan"),
    legend.text= element_text(size = 16),#, family="NimbusSan"),
    legend.position = "top",
    legend.box="vertical",
    legend.margin=margin(14.2,0,14.2,0),
    
    panel.grid.minor = element_line(size = 0.125),
    panel.grid.major = element_line(size = 0.25)
  ) +
  #xlim(c(0, 300)) +
  #ylim(c(0, 500)) +
  scale_x_continuous(trans='log10') +
  xlab("Years ago") + 
  ylab("effective population size") +
  ggtitle(title)

ggsave(plot = p, paste0(path_msmc2_final_txt, '.pdf'), width = 60, height = 50, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
