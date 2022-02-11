# Usage: Rscript get_recombination_rate_plot.R S1-indexes/ $SEG_LENGTH $N $LENGTH_OF_SEQ out ./

args <- commandArgs()
dir_with_indexes_and_sumsldhat <- args[6]
segLength <- as.numeric(args[7])
nn <- as.numeric(args[8])
ll <- as.numeric(args[9])
out <- args[10]

segs = ll/segLength

impute <- function(data, index, two, segs)
{
  while (!(length(index) == 0) || (length(index) == segs)) {
    data.to.impute = sapply(index, get_impute_data, data,
                            two = two, segs = segs)
    help.ind = which(colSums(apply(data.to.impute, 2, is.na)) ==
                       0)
    if (length(help.ind) > 0) {
      data[index[help.ind]] = apply(matrix(data.to.impute[,
                                             help.ind], nrow = 2), 2, mean)
      index = as.numeric(which(is.na(data)))
    }
    else {
      data.to.impute = sapply(index, get_impute_data,
                              data, two = F, segs = segs)
      ct = 1
      help.ind = c()
      while (length(help.ind) == 0) {
        help.ind = which(colSums(apply(data.to.impute,
                                       2, is.na)) == ct)
        ct = ct + 1
      }
      if (length(help.ind) == 1) {
        imp.ind = help.ind
      }
      else {
        imp.ind = help.ind[colSums(apply(data.to.impute[1:2,
                                                        help.ind], 2, is.na)) < 2][1]
        if (is.na(imp.ind)) {
          imp.ind = help.ind[1]
        }
      }
      data[index[imp.ind]] = weighted.mean(data.to.impute[,
                                             imp.ind], c(1 / 3, 1 / 3, 1 / 6, 1 / 6), na.rm = T)
      index = as.numeric(which(is.na(data)))
    }
  }
  return(data)
}

get_smuce <- function(help, segs, alpha, ll, quant = 0.35, rescale, constant, demography, regMod)
{
  gam <- 0.5
  eps <- 0
  help$MaxChi[is.infinite(help$MaxChi)] <- NA

  if (length(regMod) == 1 && !demography) {
    pr1 = predict(object = LDJump::mod.full, newdata = help)
  }
  # ToDo: I can''t install LDJump on octopus, but for now I am using regMod = '' and demography = F
  # if (length(regMod) == 1 && demography) {
  #   pr1 = predict(object = LDJump::mod.full.demo, newdata = help)
  # }
  # if (length(regMod) > 1) {
  #   pr1 = predict(object = regMod, newdata = help)
  # }
  pr1[is.na(pr1)] = -1 / gam
  ind.na = is.na(pr1)
  ind.q = which(round(seq(0.1, 0.5, by = 0.05), 2) == quant)
  pr.cor = predict(object = LDJump::list.quantile.regs[[ind.q]], newdata = data.frame(x = pr1))
  pr.cor = ifelse(pr.cor < -1 / gam, -1 / gam, pr.cor)
  pr.cor.nat = (pr.cor * gam + 1)^(1 / gam) - eps
  help = help[, 1:8]
  if (length(regMod) == 0) {
    ind = as.numeric(which(is.na(rowSums(help))))
  }
  else {
    ind = as.numeric(which(is.na(rowSums(help[, -c(ncol(help) -
                                                     1, ncol(help))]))))
  }
  ind = c(ind, as.numeric(which(is.infinite(rowSums(help)))),
          as.numeric(which(rowSums(help) == 0)))
  ind = sort(unique(ind))
  pr.cor.nat[ind] = NA
  pr.cor.nat = impute(pr.cor.nat, ind, two = T, segs = segs)
  if (constant) {
    return(list(pr.cor.nat))
  }
  else {
    idx = c(1, 2, 4:9, 11, 13)
    if (length(alpha) > 1) {
      temp.cor.back = list()
      for (i in 1:length(alpha)) {
        temp = stepR::stepFit(pr.cor.nat, alpha = alpha[i],
                              family = "gauss", confband = T)
        if (rescale) {
          for (idxs in idx) {
            temp[[idxs]] = temp[[idxs]] / segs * ll
          }
        }
        temp.cor.back[[i]] = temp
      }
    }
    else {
      temp.cor.back = stepR::stepFit(pr.cor.nat, alpha = alpha,
                                     family = "gauss", confband = T)
      if (rescale) {
        for (idxs in idx) {
          temp.cor.back[[idxs]] = temp.cor.back[[idxs]] / segs *
            ll
        }
      }
    }
    return(list(temp.cor.back, pr.cor.nat, ind))
  }
}

helper_new <- c()
demography <- F
regMod <- ""
alpha <- 0.05
quant <- 0.35
rescale <- F
constant <- F
status <- T
polyThres <- 0

setwd(dir_with_indexes_and_sumsldhat)

helper <- read.table('indexes.tsv', sep = '\t')

hahe <- helper[, 1]
tajd <- helper[, 2]
haps <- helper[, 3] / segLength / nn

l.files <- list.files(path=dir_with_indexes_and_sumsldhat, pattern = paste0("Sums_part_main_", out))
l.files <- l.files[
  order(as.numeric(regmatches(l.files, regexpr("[0-9]+", l.files))))
]
named.list <- lapply(l.files, read.table, sep = "=")
sums <- data.table::rbindlist(named.list)

indices <- seq(1, nrow(sums), by = 6)
apwd <- sums[indices + 1, 2] / segLength
vapw <- sums[indices + 5, 2] / segLength
colnames(apwd) <- "apwd"
colnames(vapw) <- "vapw"

wath <- helper[, 6] / segLength
MaxChi <- helper[, 7]
NSS <- helper[, 8]
phi.mean <- helper[, 9]
phi.var <- helper[, 10]
helper <- data.frame(
  cbind(hahe, tajd, haps, apwd, vapw, wath, MaxChi, NSS, phi.mean, phi.var),
  row.names = 1:nrow(helper)
)
helper_new <- rbind(helper_new, helper)
full.list <- get_smuce(
  helper_new, segs, alpha, quant = quant,
  ll, rescale = rescale,
  constant = constant, demography = demography,
  regMod = regMod
)
if (!constant) {
  seq.full.cor <- full.list[[1]]
  pr.full.cor <- full.list[[2]]
  ind <- full.list[[3]]
  if (length(ind) > 0) {
    warning(paste(
      "Recombination rates were imputed for the following segments:",
      toString(ind), sep = "")
    )
  }
  results <- list(`Estimated recombination map` = seq.full.cor,
                  `Constant estimates:` = pr.full.cor, `Summary Statistics` = helper_new,
                  alpha = alpha, `Sample Size` = nn, `Sequence Length` = ll,
                  `Segment Length` = segLength, `Imputed Segments` = ind)
} else {
  pr.full.cor <- full.list[[1]]
  results <- list(`Constant estimates` = pr.full.cor, `Summary Statistics` = helper_new,
                  alpha = alpha, `Sample Size` = nn, `Sequence Length` = ll,
                  `Segment Length` = segLength)
}

postscript("Results.eps", horiz = F)
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()


#system(paste0("for f in Sums_part_main_", out, "*; do rm -f $f; done"))
