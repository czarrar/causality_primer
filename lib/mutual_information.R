# Slow approach to getting mutual information
# but maybe more instructive?
# borrowed from http://stackoverflow.com/questions/20491028/optimal-way-for-calculating-columnwise-mutual-information-using-numpy
library(entropy)

X <- xy[,1]; Y <- xy[,2]
calc_mi <- function(X, Y, ...) {
  ###
  # slow way to get mutual information
  # but maybe more learning friendly
  ###
  
  # again this is more for learning
  shannon.entropy <- function(counts) {
    freqs <- counts/sum(counts)
    freqs <- freqs[freqs>0] # want only positive freqs/counts
    H     <- -sum(freqs * log2(freqs))
    return(H)
  }
  
  # also have hist2d and freq2d for 2d histogram
  c_XY  <- discretize2d(X, Y, numBins1=10, numBins2=10)
  # note theoretical maximum log(numBins1*numBins2)
  c_X   <- rowSums(c_XY)
  c_Y   <- colSums(c_XY)
  
  H_X   <- shannon.entropy(c_X)
  H_Y   <- shannon.entropy(c_Y)
  H_XY  <- shannon.entropy(c_XY)
  
  MI    <- H_X + H_Y - H_XY
  
  H_X.Y <- H_X - MI
  H_Y.X <- H_Y - MI
  
  
  entropy(c_XY)
  
  p_X   <- freqs.empirical(c_X)
  p_Y   <- freqs.empirical(c_Y)
  p_XY  <- freqs.empirical(c_XY)
  
  #tmp   <- matrix(rep(p_X, 10), 10, 10)
  #p_X.Y <- p_XY/tmp
  #p_X.Y <- rowSums(p_X.Y)/rowSums(p_X.Y!=0)
  #p_X.Y[is.na(p_X.Y)] <- 0
  #p_X.Y
  
  # p(X|Y)
  p_X.Y <- rowSums(p_XY)/rowSums(p_XY!=0)
  p_X.Y <- p_X.Y/p_X
  p_X.Y[is.na(p_X.Y)] <- 0
  
  # p(Y|X)
  p_Y.X <- colSums(p_XY)/colSums(p_XY!=0)
  p_Y.X <- p_Y.X/p_Y
  p_Y.X[is.na(p_Y.X)] <- 0
  
  # Plot
  library(ggplot2)
  nbins <- 10
  df <- data.frame(
    x = rep(1:10, 4), 
    id = rep(c("p(X)", "p(Y)", "p(X|Y)", "p(Y|X)"), each=nbins), 
    node = rep(c("X", "Y", "X", "Y"), each=nbins), 
    conditional = rep(c("0", "0", "1", "1"), each=nbins), 
    probs = c(p_X/sum(p_X), p_Y/sum(p_Y), p_X.Y/sum(p_X.Y), p_Y.X/sum(p_Y.X))
  )
  
  # Distribution Measures
  ddply(df, .(id), function(sdf) {
    x   <- sdf$probs
    res <- c(var=var(x), skewness=skewness(x), kurtosis=kurtosis(x), 
             entropy=entropy.empirical(x, "log2"))
    round(res, 3)
  })
  
  qplot(x=1:10, y=p_X.Y-p_Y.X, xlab="Bins", ylab="p(X|Y) - p(Y|X)") + 
    geom_hline(yintercept=0, linetype=2) + 
    geom_line()
  
  # a faster way to get the mutual information
  # mi.empirical(c_XY)
}
