options(stringsAsFactors = F)

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
        lapply(as.list(1:minp), function(x){return(sum((mat[,i] < (1*10^(-x))) & (mat[,j] < (1*10^(-x)))))}),
        envir = .GlobalEnv)
    }
  }
}

MakeMatAndCalcOverlapClumps <- function(subtype, phenos, minp = 20, pval = "P"){
  s <- subtype
  
  # for (p in phenos){
  #   df <- get(sprintf("%s_%s_clumps", p, s))
  #   df <- df[,pval]
  #   assign(sprintf("%s_%s", p, "pvals"), df)
  # }
  # mat <- matrix(as.numeric(unlist(mget(paste(phenos, "pvals", sep = "_")))), ncol = 5, byrow = F)
  for(i in 1:5){
    for(j in i:5){
      df_i <- get(sprintf("%s_%s_clumps", phenos[i], s))
      # hits_i <- df_i$SNP[df_i$P < (1*10^(-x))]
      df_j <- get(sprintf("%s_%s_clumps", phenos[j], s))
      # hits_j <- df_j$SNP[df_j$P < (1*10^(-x))]
      # overlap <- intersect(hits_i, hits_j)
      overlap <- lapply(as.list(1:minp), function(x){return(length(intersect(df_i$SNP[df_i$P < (1*10^(-x))], df_j$SNP[df_j$P < (1*10^(-x))])))})
      assign(
        sprintf("%s%sX%s%s", phenos[i], s, phenos[j], s),
        overlap,
        #lapply(as.list(1:minp), function(x){return(sum((mat[,i] < (1*10^(-x))) & (mat[,j] < (1*10^(-x)))))}),
        envir = .GlobalEnv)
    }
  }
}
# s <- "sao"
# df <- get(sprintf("%s_%s_clumps", phenos[1], s))
# overlap <- df$SNP#[df$P < 10e-8]
# for(i in 2:5){
#   df <- get(sprintf("%s_%s_clumps", phenos[i], s))
#   overlap <- intersect(overlap, df$SNP)#[df$P < 10e-8])
# }
# 
# # sao...
# s <- "sao"
# df <- get(sprintf("%s_%s_clumps", phenos[1], s))
# overlap <- df$CHR[df$P < 10e-8]
# for(i in 2:5){
#   df <- get(sprintf("%s_%s_clumps", phenos[i], s))
#   overlap <- intersect(overlap, df$CHR)[df$P < 10e-8]
# }

MakeOverlapPlot <- function(ind, subtype, n, phenos, save = F, fn = NULL, ff = png, log = T){
  nphenos <- length(phenos)
  cols <- rainbow(nphenos)
  if(save){ff(fn)}
  s <- subtype
  phe <- phenos[ind]
  plot(
    unlist(get(sprintf("%s%sX%s%s", phe, s, phe, s))), 
    log = ifelse(log, "y", ""), 
    pch = 16, 
    cex = 0.8, 
    col = cols[ind*ind], 
    ylim = c(1, 1e6), 
    type = "n",
    main = sprintf("Overlap in significant hits between %s and all phenotypes \n for %s", phe, s),
    xlab = "P-value significance threshold (1e-x)",
    ylab = sprintf("Overlap in significant hits between %s and all phenotypes", phe, s))
  for(j in 1:nphenos){
    #if(i == j){pch = 18; col = "black"}else{pch = 16; col = cols[i*j]}
    objname <- sprintf("%s%sX%s%s", phenos[ind], s, phenos[j], s)
    if(exists(objname)){obj <- get(objname)}else{objname <- sprintf("%s%sX%s%s", phenos[j], s, phenos[ind], s); obj <- get(objname)}
    lines(unlist(obj), type = "l", col = cols[j], lwd = n[j]/1000)
  }
  legend("topright", col = cols[1:nphenos], legend = phenos, fill = cols[1:nphenos], border = cols[1:nphenos])
  if(save){dev.off()}
}
