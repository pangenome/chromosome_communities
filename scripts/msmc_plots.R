library(tidyverse)

mu <- 1.25e-8
gen <- 30
x1 <- read.table('~/git/chromosome_communities/msmc/HG002.pq_arms.phased.chr13.hetsep.p_arms.final.txt', sep = '\t', header = T)
x2 <- read.table('~/git/chromosome_communities/msmc/HG002.pq_arms.phased.chr13.hetsep.q_arms.final.txt', sep = '\t', header = T)
plot(
  x1$left_time_boundary/mu*gen,
  (1/x1$lambda)/(2*mu),
  log="x",
 # ylim=c(0,100000),
  type="n",
  xlab="Years ago", ylab="effective population size",
  main="HG002 - All acros"
)
lines(x1$left_time_boundary/mu*gen, (1/x1$lambda)/(2*mu), type="s", col="red")
lines(x2$left_time_boundary/mu*gen, (1/x2$lambda)/(2*mu), type="s", col="blue")
legend("topright",legend=c("p-arm", "q-arm"), col=c("red", "blue"), lty=c(1,1), pch=c(1,3), cex=0.5)





mu <- 1.25e-8
gen <- 30
crossPopDat<-read.table("~/git/chromosome_communities/msmc/HG002.pq_arms.phased.chr13.hetsep.combined.final.txt", header=TRUE)
plot(
  crossPopDat$left_time_boundary/mu*gen,
  2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11),
  #xlim=c(1000,700000),
  #ylim=c(0,1),
  type="n",
  xlab="Years ago", ylab="relative cross-coalescence rate",
  main="HG002 - All acros\np-arms vs p-arms"
  )
lines(crossPopDat$left_time_boundary/mu*gen, 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11), type="s")

