# MENTAPPR: Compute MaxEnt approximation of negentropy and differential entropy
# Aapo Hyvarinen, May 2012
# and adapted to R by Zarrar Shehzad, July 2014
# Based on NIPS*97 paper, www.cs.helsinki.fi/u/ahyvarin/papers/NIPS97.pdf
#  but using a new nonlinearity
# Input: sample of continous-valued random variable as a vector. 
#         For matrices, entropy computed for each column
# Output: (differential) entropy with negentropy as an attribute

mentappr <- function(x, standardize=FALSE) {
  if (standardize == T) {
    x <- scale(x)
    xstd <- attr(x, "scaled:scale")
  } else {
    xstd <- sd(x)
  }
  
  # Constants we need (TODO: precompute)
  k1 <- 36/(8*sqrt(3)-9)
  gamma <- 0.37457; 
  k2 <- 79.047;
  gaussianEntropy <- log(2*pi)/2+1/2
  
  # This is negentropy
  negentropy <- k2*(mean(log(cosh(x)))-gamma)^2+k1*mean(x*exp(-(x^2)/2))^2
  
  # This is entropy
  entropy <- gaussianEntropy - negentropy + log(xstd)
  attr(entropy, "negentropy") <- negentropy
  
  return(entropy)
}

