options(stringsAsFactors = F)
library(extraDistr)

overlap <- function(pthres, z1, z2){
  zthres <- abs(qnorm(p = pthres))
  #hi1 <- z1 > zthres | z1 < -zthres
  #hi2 <- z2 > zthres | z2 < -zthres
  #overlap <- sum(hi1 & hi2 & (sign(z1) == sign(z2)))
  # is equivalent to:
  red_urn <- z1 > zthres
  red_draw <- red_urn & (z2 > zthres)
  blue_urn <- z1 < -zthres
  blue_draw <- blue_urn & (z2 < -zthres)
  overlap <- sum(red_draw | blue_draw)
  sum_black_urn <- 10437951 - sum(red_urn) - sum(blue_urn)
  nsign_draw <- sum(z2 > zthres | z2 < -zthres)
  sum_black_draw <- nsign_draw - overlap
  x <- c(sum(red_draw), sum(blue_draw), sum_black_draw)
  n <- c(sum(red_urn), sum(blue_urn), sum_black_urn)
  k <- nsign_draw
  pval <- dmvhyper(x, n, k, log = FALSE)
  names(overlap) <- pval
}

# overlap <- function(pthres, z1, z2){
#   p1 <- 2 * pnorm(-abs(z1))
#   p2 <- 2 * pnorm(-abs(z2))
#   sign1 <- p1 < (1*10^(-p))
#   sign2 <- p2 < (1*10^(-p))
# }
# 
# overlap <- function(p, i, j){
#   sign_i <- (mat[,i] < (1*10^(-p)))
#   sign_j <- (mat[,j] < (1*10^(-p)))
#   #n <- sum(sign_i)
#   overlap <- sum(sign_i & sign_j)
#   m <- sum(sign_i)
#   n <- 10437951-m
#   k <- sum(sign_j)
#   #pval <- binom.test(overlap, sum(sign_j), p)$p.value
#   pval <- phyper(q = overlap, m = m, n = n, k = k, lower.tail = F, log.p = T)
#   names(overlap) <- pval
#   return(overlap)
#  # return(list(overlap = overlap, log.pval = pval))
# }

MakeMatAndCalcOverlap <- function(metric, subtype, phenos, minp = 20){
  #phenos <- c("CCSc", "CCSp", "TOAST", "union_original", "intersect")
  m <- metric
  s <- subtype
  n <- length(phenos)
  for (p in phenos){
    df <- scan(sprintf("%s_%s_%s", m, p, s), what = "character")
    df <- df[df != df[1]]
    print(length(df))
    assign(sprintf("%s_%s", p, m), df)
  }
  mat <- matrix(as.numeric(unlist(mget(paste(phenos, m, sep = "_")))), ncol = n, byrow = F)
  for(i in 1:n){
    for(j in i:n){
      assign(
        sprintf("%s%sX%s%s", phenos[i], s, phenos[j], s), 
        lapply(as.list(1:minp), overlap, i, j),
        # lapply(as.list(1:minp), function(x){return(sum((mat[,i] < (1*10^(-x))) & (mat[,j] < (1*10^(-x)))))}),
        envir = .GlobalEnv)
    }
  }
}

MakeMatAndCalcOverlapClumps <- function(subtype, phenos, minp = 20, pval = "P"){
  s <- subtype
  for(i in 1:5){
    for(j in i:5){
      df_i <- get(sprintf("%s_%s_clumps", phenos[i], s))
      df_j <- get(sprintf("%s_%s_clumps", phenos[j], s))
      overlap <- lapply(as.list(1:minp), function(x){return(length(intersect(df_i$SNP[df_i$P < (1*10^(-x))], df_j$SNP[df_j$P < (1*10^(-x))])))})
      assign(
        sprintf("%s%sX%s%s", phenos[i], s, phenos[j], s),
        overlap,
        envir = .GlobalEnv)
    }
  }
}

MakeOverlapPlot <- function(ind, subtype, n, phenos, save = F, fn = NULL, ff = png, log = T, signline = T){
  nphenos <- length(phenos)
  cols <- rainbow(nphenos)
  if(save){ff(fn)}
  s <- subtype
  phe <- phenos[ind]
  #yl <- ifelse(log, as.vector(1, 1e6), as.vector(1,70))
  if(log){ylim <- c(1, 1e6)}else{ylim <- c(1, 70)}
  plot(
    unlist(get(sprintf("%s%sX%s%s", phe, s, phe, s))), 
    log = ifelse(log, "y", ""), 
    pch = 16, 
    cex = 0.8, 
    col = cols[ind*ind], 
    ylim = ylim, 
    type = "n",
    main = sprintf("Overlap in significant hits between %s and all phenotypes \n for %s", phe, s),
    xlab = "P-value significance threshold (1e-x)",
    ylab = sprintf("Overlap in significant hits between %s and all phenotypes", phe, s))
  if(signline){
    abline(v = 8)
  }
  logpvals <- c()
  for(j in 1:nphenos){
    objname <- sprintf("%s%sX%s%s", phenos[ind], s, phenos[j], s)
    if(exists(objname)){obj <- get(objname)}else{objname <- sprintf("%s%sX%s%s", phenos[j], s, phenos[ind], s); obj <- get(objname)}
    logpvals <- c(logpvals, names(obj[[8]]))
    lines(unlist(obj), type = "l", col = cols[j], lwd = n[j]/1000)
  }
  logpvals <- as.integer(logpvals)
  legend("topright", col = cols[1:nphenos], legend = paste(phenos, "(", formatC(logpvals, digits = 0, ), ")"), fill = cols[1:nphenos], border = cols[1:nphenos])
  if(save){dev.off()}
}
