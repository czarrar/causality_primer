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


# PWLING: pairwise causality measures in linear non-Gaussian model
# Version 1.2, Aapo Hyvarinen, Feb 2013
# Input: Data matrix with variables as rows, and index of method [1...5]
# Output: Matrix LR with likelihood ratios
#         If entry (i,j) in that matrix is positive, 
#         estimate of causal direction is i -> j
# Methods 1...5 are as follows:
#   1: General entropy-based method, for variables of any distribution
#   2: First-order approximation of LR by tanh, for sparse variables
#   3: Basic skewness measure, for skewed variables
#   4: New skewness-based measure, robust to outliers
#   5: Dodge-Rousson measure, for skewed variables
#   If you want to use method 3 or 4 without skewness correction, 
#      input -3 or -4 as the method.
# See http://www.cs.helsinki.fi/u/ahyvarin/code/pwcausal/ for more information

pwling <- function(X, method, C=NULL, verbose=TRUE) {
  library(moments) # for measure of skewness
  
  if (verbose) {
    suppressMessages(library(niftir)) # for progressbar
    vcat <- function(msg, ...) cat(sprintf(msg, ...), "\n")
  } else {
    vcat <- function(x) invisible(NULL)
  }
  
  # Get size parameters
  nvars  <- nrow(X)
  nsamps <- ncol(X)
  
  # Standardize each variable
  vcat("standardizing")
  X <- t(scale(t(X)))
  
  # If using skewness measures with skewness correction, make skewnesses positive
  if (method == 3 || method == 4) {
    vcat("making skewness positive")
    for (i in 1:nvars) X[i,] <- X[i,] * sign(skewness(X[i,]))
  }
  
  # Compute covariance matrix
  vcat("computing covariance")
  if (is.null(C)) C <- cov(t(X))
  
  # Compute causality measures
  ############################
  
  vcat("computing causality")
  
  # General entropy-based method, for variables of any distribution
  if (method == 1) {
    vcat("...using EB method")
    
    # Initialize output matrix
    LR <- matrix(0, nvars, nvars)
    
    # Loop throgh pairs
    if (verbose) pb <- progressbar(nvars)
    for (i in 1:nvars) {
      for (j in 1:nvars) {
        if (i == j) next
        res1 <- X[j,] - C[j,i]*X[i,]
        res2 <- X[i,] - C[i,j]*X[j,]
        LR[i,j] <- mentappr(X[j,]) - mentappr(X[i,]) - mentappr(res1) + mentappr(res2)
      }
      if (verbose) update(pb, i)
    }
    if (verbose) end(pb)
  } 
  # First-order approximation of LR by tanh, for sparse variables
  else if (method == 2) {
    vcat("...using Tanh method")
    LR <- C * (X %*% tanh(t(X)) - tanh(X) %*% t(X))/nsamps
  }
  # Basic skewness measure, for skewed variables
  else if (abs(method) == 3) {
    if (sign(method) == 1) vcat("...using Skew method")
    else vcat("...using SkewE method")
    LR <- C * (-X %*% t(X)^2 + X^2 %*% t(X))/nsamps
  }
  # New skewed measure, robust to outliers
  else if (abs(method) == 4) {
    if (sign(method) == 1) vcat("...using RSkew method")
    else vcat("...using RSkewE method")
    Xthr <- X
    Xthr[X<0] <- 0
    gX <- log(cosh(Xthr))
    LR <- C * (-((X %*% t(gX))/nsamps) + ((gX %*% t(X))/nsamps))
  }
  # Dodge-Rousson measure, for skewed variables
  else if (method == 5) {
    vcat("...using DRSkew measure")
    LR <- (-(X %*% t(X)^2)^2 + (X^2 %*% t(X))^2)/nsamps
  }
  
  return(LR)
}
