library(qqman)

endsWith <- function(x, suffix){
  n <- nchar(x)
  return(substr(x, n - nchar(suffix) + 1L, n) == suffix)
}

LoadData <- function(nchroms=22, df_name_base="imputed_chr_", path="bolt.stats.group4.EUR/bolt.chr", file_suffix=".stats"){
  print(file_suffix)
  if (endsWith(file_suffix, ".gz")) {func <- gzfile; print("this is a .gz file, will use gzfile() function")} else {func <- file}
  print(func)
  chrs <- 1:nchroms
  for(c in chrs){
    print(c)
    df_name <- sprintf("%s%i", df_name_base, c)
    assign(
      df_name, 
      read.table(func(description=sprintf("%s%i%s", path, c, file_suffix), open = "r"), header = T),
      envir = .GlobalEnv
    )
  }
}

MakeManhattan <- function(concat, save = FALSE, fn = NA, pcol = "P_BOLT_LMM_INF", rmsnps = NULL, snpcolnr = 1){
  if(save){png(fn)}
  if(!is.null(rmsnps)){concat <- concat[!(concat[,snpcolnr] %in% rmsnps),]}
  manhattan(concat, p = pcol, ylim = c(0,30), col = c(rgb(154, 137, 168, maxColorValue = 255), "lightgrey"))
  if(save){dev.off()}
}

## QQ Plot function ##
plotQQ <- function(z,color,cex){
  p <- 2*pnorm(-abs(z))
  p <- sort(p)
  expected <- c(1:length(p))
  lobs <- -(log10(p))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  ## plots all points with p < 1
  p_sig = subset(p,p<1)
  points(lexp[1:length(p_sig)], lobs[1:length(p_sig)], pch=23, cex=.3, col=color, bg=color)
  
  ## samples 5000 points from p > 0.01
  n=5001
  i<- c(length(p)- c(0,round(log(2:(n-1))/log(n)*length(p))),1)
  lobs_bottom=subset(lobs[i],lobs[i] <= 2)
  lexp_bottom=lexp[i[1:length(lobs_bottom)]]
  points(lexp_bottom, lobs_bottom, pch=23, cex=cex, col=color, bg=color)
}


MakeMAFQQ <- function(concat, save = FALSE, fn = NA, mx = 25, freqcol = "A1FREQ", pcol = "P_BOLT_LMM_INF"){
  S <- concat
  S$FRQ <- S[,freqcol]
  S$P <- S[,pcol]
  
  pvals_lo1=subset(S,(S$FRQ > 0.20 & S$FRQ < 0.8))
  pvals_lo2=subset(S,((S$FRQ < 0.20 & S$FRQ > 0.05) | (S$FRQ > 0.8 & S$FRQ < 0.95)))
  pvals_lo3=subset(S,((S$FRQ < 0.05 & S$FRQ > 0.01) | (S$FRQ > 0.95 & S$FRQ < 0.99)))
  pvals_lo4=subset(S,((S$FRQ < 0.01 & S$FRQ > 0.001) | (S$FRQ > 0.99 & S$FRQ < 0.999)))
  pvals_lo5=subset(S,(S$FRQ < 0.001 | S$FRQ > 0.999))
  
  z=qnorm(S$P/2)
  z_lo1=qnorm(pvals_lo1$P/2)
  z_lo2=qnorm(pvals_lo2$P/2)
  z_lo3=qnorm(pvals_lo3$P/2)
  z_lo4=qnorm(pvals_lo4$P/2)
  z_lo5=qnorm(pvals_lo5$P/2)
  
  lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
  l1 = round(median(z_lo1^2,na.rm=T)/qchisq(0.5,df=1),3)
  l2 = round(median(z_lo2^2,na.rm=T)/qchisq(0.5,df=1),3)
  l3 = round(median(z_lo3^2,na.rm=T)/qchisq(0.5,df=1),3)
  l4 = round(median(z_lo4^2,na.rm=T)/qchisq(0.5,df=1),3)
  l5 = round(median(z_lo5^2,na.rm=T)/qchisq(0.5,df=1),3)
  
  if(save){png(fn)}
  
  plot(c(0,mx), c(0,mx), col="gray25", lwd=3, type="l", xlab="Expected Distribution (-log10 of P value)", ylab="Observed Distribution (-log10 of P value)", xlim=c(0,mx), ylim=c(0,mx), las=1, xaxs="i", yaxs="i", bty="l",main=c(substitute(paste("QQ plot: ",lambda," = ", lam),list(lam = lambda)),expression()))
  
  plotQQ(z,rgb(255,79,0,maxColorValue=255),0.4)
  plotQQ(z_lo5,"lightpink2",0.3)
  plotQQ(z_lo4,"purple",0.3)
  plotQQ(z_lo3,"deepskyblue1",0.3)
  plotQQ(z_lo2,"slateblue3",0.3)
  plotQQ(z_lo1,"olivedrab2",0.3)
  
  ## provides legend
  legend("bottomright",legend=c("Expected (null)","Observed",
                                substitute(paste("MAF > 0.20 [", lambda," = ", lam, "]"),list(lam = l1)),expression(),
                                substitute(paste("0.05 < MAF < 0.20 [", lambda," = ", lam, "]"),list(lam = l2)),expression(),
                                substitute(paste("0.01 < MAF < 0.05 [", lambda," = ", lam, "]"),list(lam = l3)),expression(),
                                substitute(paste("0.001 < MAF < 0.01 [", lambda," = ", lam, "]"),list(lam = l4)),expression(),
                                substitute(paste("MAF < 0.001 [", lambda," = ", lam, "]"),list(lam = l5)),expression()),
         pch=c((vector("numeric",6)+1)*23), cex=1.1, pt.cex=1.5, pt.bg=c("grey25",rgb(255,79,0,maxColorValue=255),"olivedrab2","slateblue3","deepskyblue1","purple","lightpink2"))
  
  if(save){dev.off()}
}

MakeInfoQQ <- function(concat, save = FALSE, fn = NA, mx = 25, infocol = "INFO", pcol = "P_BOLT_LMM_INF"){
  S <- concat
  S$INFO <- S[,infocol]
  S$P <- S[,pcol]
  
  pvals_lo0=subset(S,(S$INFO > 1.1))
  pvals_lo1=subset(S,(S$INFO > 0.9 & S$INFO < 1.1))
  pvals_lo2=subset(S,((S$INFO <= 0.9 & S$INFO > 0.8)))
  pvals_lo3=subset(S,((S$INFO <= 0.8 & S$INFO > 0.7)))
  pvals_lo4=subset(S,((S$INFO <= 0.7 & S$INFO > 0.6)))
  pvals_lo5=subset(S,(S$INFO <= 0.6))
  
  z=qnorm(S$P/2)
  z_lo0=qnorm(pvals_lo0$p/2)
  z_lo1=qnorm(pvals_lo1$P/2)
  z_lo2=qnorm(pvals_lo2$P/2)
  z_lo3=qnorm(pvals_lo3$P/2)
  z_lo4=qnorm(pvals_lo4$P/2)
  z_lo5=qnorm(pvals_lo5$P/2)
  
  
  ## calculates lambda
  lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
  l0 = round(median(z_lo0^2,na.rm=T)/qchisq(0.5,df=1),3)
  l1 = round(median(z_lo1^2,na.rm=T)/qchisq(0.5,df=1),3)
  l2 = round(median(z_lo2^2,na.rm=T)/qchisq(0.5,df=1),3)
  l3 = round(median(z_lo3^2,na.rm=T)/qchisq(0.5,df=1),3)
  l4 = round(median(z_lo4^2,na.rm=T)/qchisq(0.5,df=1),3)
  l5 = round(median(z_lo5^2,na.rm=T)/qchisq(0.5,df=1),3)
  
  if(save){png(fn)}
  
  plot(c(0,mx), c(0,mx), col="black", lwd=3, type="l", xlab="Expected Distribution (-log10 of P value)", ylab="Observed Distribution (-log10 of P value)", xlim=c(0,mx), ylim=c(0,mx), las=1, xaxs="i", yaxs="i", bty="l",main=c(substitute(paste("QQ plot: ",lambda," = ", lam),list(lam = lambda)),expression()))
   
  #lapply(list("lightpink2", "purple", "deepskyblue1", "slateblue3", "olivedrab2"), plotQQ)
  
  plotQQ(z,rgb(255,79,0,maxColorValue=255),0.4)
  plotQQ(z_lo5,"lightpink2",0.3)
  plotQQ(z_lo4,"purple",0.3)
  plotQQ(z_lo3,"deepskyblue1",0.3)
  plotQQ(z_lo2,"slateblue3",0.3)
  plotQQ(z_lo1,"olivedrab2",0.3)
  
  ## provides legend
  legend("bottomright",legend=c("Expected (null)","Observed",
                        substitute(paste("Imp Qual > 1.1 [", lambda," = ", lam, "]"),list(lam = l0)),expression(),
                        substitute(paste("0.9 < Imp Qual < 1.1 [", lambda," = ", lam, "]"),list(lam = l1)),expression(),
                        substitute(paste("0.8 < Imp Qual < 0.9 [", lambda," = ", lam, "]"),list(lam = l2)),expression(),
                        substitute(paste("0.7 < Imp qual < 0.8 [", lambda," = ", lam, "]"),list(lam = l3)),expression(),
                        substitute(paste("0.6 < Imp qual < 0.7 [", lambda," = ", lam, "]"),list(lam = l4)),expression(),
                        substitute(paste("Imp qual < 0.6 [", lambda," = ", lam, "]"),list(lam = l5)),expression()),
         pch=c((vector("numeric",6)+1)*23), cex=1.1, pt.cex=1.5, pt.bg=c("black",rgb(255,79,0,maxColorValue=255),"palevioletred3","olivedrab2","slateblue3","deepskyblue1","purple","lightpink2"))
  
  if(save){dev.off()}
  
}

h2_and_snpcount <- function(corr_factor=5e6, snp_counts=NULL, ylim=c(0,20)){
  if(is.null(snp_counts)){
    snp_counts <- scan("/Users/jvonberg/git/LMM/snp_counts.txt")
  }
  par(mar = c(5,5,5,5))
  h2_bar <- barplot(h2s, ylim = ylim, names.arg = 1:22, xlab = "chromosome", ylab = "heritability (h2)", main = "Heritability for a certain left-out chromosome")
  lines(x = h2_bar, y = snp_counts/corr_factor, col = "purple", lwd = 2)
  points(x = h2_bar, y = snp_counts/corr_factor, col = "purple", cex = 0.6)
  axis(4, at = seq(0, 0.25, by = 0.05), labels = seq(0,0.25*corr_factor,by = 0.05*corr_factor))
  mtext(text = "SNP-count", side = 4, line = 2)
}

NrOfSignificantHits <- function(dfname, threshold = 1e-8, pcol = "P_BOLT_LMM_INF"){
  return(sum(get(dfname)[,pcol] < threshold))
}

NrOfSignificantHits_vectorized <- Vectorize(NrOfSignificantHits, vectorize.args = "dfname")

PvalAtKnownHits <- function(dfname, hits = 1:6, pcol = "P_BOLT_LMM_INF", minuslog = F){
  if(minuslog){func <- function(x){return(-log10(x))}}
  else{func <- function(x){return(x)}}
  return(func(get(dfname)[hits,pcol]))
}

PvalAtKnownHits_vectorized <- Vectorize(PvalAtKnownHits, vectorize.args = "dfname")