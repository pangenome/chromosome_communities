library(tidyverse)

mu <- 1.25e-8
gen <- 30
x1 <- read.table('~/git/chromosome_communities/msmc/15.pq_arms.phased.chr13.hetsep.p_arms.final.txt', sep = '\t', header = T)
x2 <- read.table('~/git/chromosome_communities/msmc/15.pq_arms.phased.chr13.hetsep.q_arms.final.txt', sep = '\t', header = T)

x1$pop <- 'p-arms'
x2$pop <- 'q-arms'

x <- rbind(x1, x2)

x$x <- x$left_time_boundary/mu*gen
x$y <- (1/x$lambda)/(2*mu)

x %>%
  ggplot(aes(x = x, y = y, color = pop)) + 
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
  scale_x_continuous(trans='log2') +
  xlab("Years ago") + 
  ylab("effective population size") +
  ggtitle('75 contigs (15 for each acro)')



mu <- 1.25e-8
gen <- 30
crossPopDat<-read.table("~/git/chromosome_communities/msmc/15.pq_arms.phased.chr13.hetsep.combined.final.txt", header=TRUE)

crossPopDat$x <- crossPopDat$left_time_boundary/mu*gen
crossPopDat$y <- 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11)

crossPopDat %>%
  ggplot(aes(x = x, y = y)) +
  geom_step()  +
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
  scale_x_continuous(trans='log2') +
  xlab("Years ago") + 
  ylab("relative cross-coalescence rate") +
  ggtitle('75 contigs (15 for each acro)')
