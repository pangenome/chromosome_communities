args <- commandArgs()
x <- args[6]
seqNamePart <- args[7]
pathPhi <- args[8]
prefix <- args[9]
output_dir <- args[10]

if (file.exists(seqNamePart)) {
  getPhi <- function (seqName, pathPhi, prefix)
  {
    phiName <- paste0(output_dir, "/", prefix, ".phi.log")
    system(paste0(pathPhi, " -f ", seqName, " -o -v >", phiName))

    results <- system(paste0("grep 'Value of maximum breakpoint is:' ", phiName, ' -c'), intern = T)[[1]]
    if (results != "0") {
      MaxChi <- as.numeric(
        strsplit(system(paste0("grep 'Value of maximum breakpoint is:' ", phiName), intern = T), ": ")[[1]][2]
      )
      if (is.infinite(MaxChi)) {
        MaxChi <- NA
      }
      NSS <- as.numeric(strsplit(system(paste("grep 'The Neighbour Similarity score is ' ", phiName), intern = T), " ")[[1]][6])
      phi.mean <- unlist(strsplit(system(paste0("grep 'Mean: ' ", phiName), intern = T), " "))
      phi.mean <- as.numeric(phi.mean[phi.mean != ""][2])
      phi.var <- unlist(strsplit(system(paste0("grep 'Variance: ' ", phiName), intern = T), " "))
      phi.var <- as.numeric(phi.var[phi.var != ""][2])
    } else {
      MaxChi <- NSS <- phi.mean <- phi.var <- 0
    }

    return(c(MaxChi, NSS, phi.mean, phi.var))
  }

  apwd <- vapw <- NA

  samp <- ape::read.dna(seqNamePart, as.matrix = T, format = "fasta")
  tajd <- pegas::tajima.test(samp)$D
  temp <- try(adegenet::DNAbin2genind(samp, polyThres = 0))
  hahe <- adegenet::Hs(temp)

  pathTmpFile <- paste0(output_dir, "/", prefix, ".tmpfile.", x, ".fa")
  system(paste0("sed '1d' ", seqNamePart, " > ", pathTmpFile))
  phis <- getPhi(
    seqName = pathTmpFile,
    pathPhi = pathPhi,
    prefix = paste0(prefix, ".", x)
  )
  system((paste0("rm ", pathTmpFile)))

  haps <- length(print(pegas::haplotype(samp)))

  wath <- pegas::theta.s(samp)
} else {
  haps <- hahe <- tajd <- apwd <- vapw <- wath <- 0
  phis <- rep(0, 4)
}

cat(paste(c(hahe, tajd, haps, apwd, vapw, wath, phis, "\n"), collapse = '\t'))
