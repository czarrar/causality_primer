# This is a streamlined adaptation of AFNI's 1dSVAR.R
# Here I do only individual-level analyses

## SETUP ######################################################################

suppressMessages(require("gsl"))      # Legendre polynomials
suppressMessages(require("vars"))     # VAR modeling 
vcat <- function(msg, ...) cat(sprintf(msg, ...), "\n")

fn_dat <- "data/sim01/sub01_ts.txt"
fn_cov <- NULL
fn_cnk <- NULL # block/trial chunk file
fn_inm <- "data/sim01/sub01_net.txt"

myData  <- read.table(fn_dat)
instMat <- read.table(fn_inm)

#diag(instMat) <- 0
#instMat <- instMat + t(instMat)
#instMat[instMat!=0] <- NA
#

#setwd("~/Dropbox/Research/yale/causality_primer/netsim")
#net   <- as.matrix(read.table("data/sim01/sub01_net.txt")) # reference
#ts    <- as.matrix(read.table("data/sim01/sub01_ts.txt"))


# myData are your ROIs
nTotal  <- nrow(myData)
nROIs   <- ncol(myData)

# Number of breaks/chunks in time-series
if (!is.null(fn_cnk)) {
  nPts <- as.vector(as.matrix(read.table(fn_cnk)))
} else {
  nPts <- nTotal
}
nChunks <- length(nPts)

# Prior detrending is NOT recommended due to potential complications
# Trend can be modeled through specifying the order of a polynomial here,
# or can be included as part of covariates later on. If you plan to model the
# trend with your own regressors or don't need to model it, choose -1 here.
# If trend has already been removed (not recommended), choose 0 here.
nPoly <- 2

maxLags <- 10

scaleTS <- TRUE

qualityCheck <- TRUE


## Covariates ##################################################################

# Any additional covariates
if (is.null(fn_cov)) {
  exData <- NULL
  nCOVs  <- 0 # TODO: is this updated later?
} else {
  exData <- read.table(fn_cov, header=TRUE)
  nCOVs  <- ncol(exData)
}

# Create exogenous variables with Legendre polynomials from gsl
# e.g. mean, linear, quadratic for nPoly = 2
gen_trends <- function(nTotal, nPoly, nPts) {
  nChunks   <- length(nPts)
  trendMat  <- as.data.frame(array(0, dim = c(nTotal, (nPoly+1)*nChunks)))
  jumpPts   <- 0
  
  for (ii in 1:nChunks) {
    chunk.inds <- (jumpPts+1):(jumpPts+nPts[ii])
    poly.inds  <- (1+(nPoly+1)*(ii-1)):((nPoly+1)*ii)
    
    trendMat[chunk.inds,poly.inds] <- 
      t(legendre_Pl_array(nPoly, seq(-1,1,len=nPts[ii])))
    
    names(trendMat)[poly.inds] <- sprintf("Run%iTrend%i", ii, seq(nPoly+1)-1)
    
    if (ii < nChunks) jumpPts <- jumpPts+nPts[ii]
  }
  
  return(trendMat)
}

if (nPoly > -1) {
  trendMat <- gen_trends(nTotal, nPoly, nPts)
  if (is.null(exData)) {
    exMat <- trendMat
  } else {
    exMat <- cbind(trendMat, exData)
  }
} else {
  exMat <- exData # if no baseline and trend, do nothing
}


## Normalization ##############################################################

# Divide each chunk (block/trial) by it's mean
normalize_chunks <- function(myData, nPts) {
  nChunks <- length(nPts)
  nROIs   <- ncol(myData)
  
  ## we don't want any negative means
  ## so set shift everything by +1000
  ## or +1000 + minimum value (if min val < -1000)
  gshift  <- 1000 + min(0, 1000 + min(myData) - 1e-6)
  newData <- myData + gshift
  
  jumpPts <- 0
  for (ii in 1:nChunks) {
    for (jj in 1:nROIs) {
      chunk.inds <- (jumpPts+1):(jumpPts+nPts[ii])
      newData[chunk.inds,jj] <- newData[chunk.inds,jj]/mean(newData[chunk.inds,jj])
    }
    if (ii < nChunks) jumpPts <- jumpPts+nPts[ii]
  }
  
  return(newData)
}

# Perform Normalization (Divide each time chunk by the mean)
if (scaleTS) {
  newData <- normalize_chunks(myData, nPts)
} else {
  newData <- myData
}


## AR Model ###################################################################

# Here I will have a max lags variable and then actually default to the maximum number of recommended lags if that variable was set
#maxLags <- 3 # SET AT TOP

critSel <- VARselect(newData, lag.max = maxLags, type = "none", exogen=exMat)

# This will give suggested orders for VAR
# including: AIC, HQ, SC, and FPE
# 
# Usually consistency exists between AIC and FPE, 
# and between HQ and QC
# 
# AIC and FPE tend to overestimate the "true order"
# Since there is no universally best criterion it might be good
# to try different analyses with various orders within the range
# covered by the 4 criteria

# Here, we will always choose the maximum one
print(critSel$selection)
nLags <- max(critSel$selection)
# TODO: also check it apply(critSel$criteria, 2, function(x) any(is.infinite(x)))

# Now we try to add the dummy variables to account for any cross run/block breaks
# this would be a value of 1 between each run/block indicating the offset there. 
if (nChunks > 1) {
  cat("iterating through number of lags taking into account blocks/runs\n")
  
  exMatMod <- ifelse(is.null(exMat))
  redo_lags <- TRUE
  max_iters <- 10
  niters    <- 0
  
  while (redo_lags || niters < max_iters) {    
    breakMat <- as.data.frame(array(0, dim = c(nTotal, (nChunks-1)*nLags)))
    jumpPts <- 0
    
    for (ii in 1:(nChunks-1)) {
      jumpPts <- jumpPts+nPts[ii]
      for (jj in 1:nLags) {
        rind <- jumpPts+jj
        cind <- (ii-1)*nLags+jj
        breakMat[rind,cind] <- 1
        #breakMat[,cind] <- c(rep(0, jumpPts+jj-1), 1, rep(0, nTotal-jumpPts-jj))
        names(breakMat)[cind] <- sprintf("Run%iLag%i", ii, jj)
      }
    }
    
    if (is.null(exMat)) {
      exMatMod <- breakMat
    } else {
      exMatMod <- cbind(breakMat, exMat)
    }
    
    critSel <- VARselect(newData, lag.max=maxLags, type="none", exogen=exMatMod)
    prev_nLags <- nLags
    nLags <- max(critSel$selection)
    if (prev_nLags == nLags) {
      vcat("choosing %i lags\n", nLags)
      redo_lags <- FALSE
    } else {
      vcat("retrying with %i lags\n", nLags)
    }
    
    niters <- niters + 1
  }
} else {
  exMatMod <- exMat
  vcat("choosing %i lags\n", nLags)
}

# Our AR model!
fm <- VAR(newData, p=nLags, type="none", exogen=exMatMod)
# what if normality tests and the different correlations do poorly
# should set to nLags - 1?
# normality.test(fm)

if (qualityCheck) {
  # the modulus of the eigenvalues (presumably less than 1 as stable condition) in the reverse characteristic polynomial; stable process is stationary, but the converse is not true
  #print("Quality check of the model:")
  if (prod(roots(fm)<1)) print("Eigenvalues of the companion coefficient matrix indciate that the VAR(p) process is stable and thus stationary") else print("The VAR(p) process seems unstable and thus is not stationary")
  print("-----------------")
  print("Normality testing of the residuals")
  print(normality.test(fm))
  print("-----------------")
  print("Serial correlation test:")
  print(serial.test(fm, lags.pt=16, lags.bg=5, type=c("PT.asymptotic")))
  print(serial.test(fm, lags.pt=16, lags.bg=5, type=c("PT.adjusted")))
  print(serial.test(fm, lags.pt=16, lags.bg=5, type=c("BG")))
  print(serial.test(fm, lags.pt=16, lags.bg=5, type=c("ES")))
  print("-----------------")
  print("Autoregressive conditional heteroskedasticity (ARCH) test")
  print(arch.test(fm))
}

# spill out the original path matrix with direction going from rows to columns
netMatR <- array(data=NA, dim=c(nLags, nROIs, nROIs))   # original path coefficient matrix
netMatT <- array(data=NA, dim=c(nLags, nROIs, nROIs))   # t values matrix
for (ii in 1:nROIs) for (jj in 1:nROIs) for (kk in 1:nLags)  { # ii: target, jj: source, kk: lag
  netMatR[kk,jj,ii] <- coef(fm)[[ii]][jj+nROIs*(kk-1), 1]  # path coefficients
  netMatT[kk,jj,ii] <- coef(fm)[[ii]][jj+nROIs*(kk-1), 3]  # t values
}

aa <- diag(nROIs)-t(as.matrix(instMat)) # A0
bb <- diag(nROIs); diag(bb) <- NA  # make B as a diagonal matrix with unknown elments

# Will first try 'iterative' (0) as per AFNI suggestion
# if that fails will try 'direct' (1)
estMeth <- 0

# TODO: should have try loop around the scoring SVAR
# if it fails, try the direct method
if (estMeth==0) {
  fm2 <- SVAR(x = fm, estmethod = "scoring", Amat = aa, Bmat = bb, 
              max.iter = 100, maxls = 1000, conv.crit = 1.0e-6)
} else if (estMeth==1) {
  fm2 <- SVAR(x = fm, estmethod = "direct", Amat = aa, Bmat = bb, 
              max.iter = 100, maxls = 1000, conv.crit = 1.0e-6)
}

instMatP <- diag(nROIs) - fm2$A   # instantaneous effects A0
nZ <- function(x) ifelse(abs(x)>1e-7, x, 1e-7)  # avoid 0 standard error division
instMatPt <- instMatP/apply(fm2$Ase, c(1,2), nZ)   # t-value for instantaneous effects A0

cat("\nInstantaneous effect matrix:\n")
print(t(instMatP))
print(fm2$LR)
if (saveInstMatP) {
  instMatPName <- as.character(readline("File name prefix for instantaneous matrix? "))
  write.table(t(instMatP), file=sprintf("%s.1D", instMatPName), append=FALSE, row.names=names(myData), col.names=names(myData))
}
cat("\nT-values for instantaneous effect matrix:\n")
print(t(instMatPt))
saveInstMatPt <- as.integer(readline("\nSave t-values for instantaneous effect matrix (0: no; 1: yes)? "))
if (saveInstMatPt) {   
  instMatPtName <- as.character(readline("File name prefix for t-values of instantaneous matrix? "))
  write.table(t(instMatPt), file=sprintf("%sSEMt.1D", instMatPtName), append=FALSE, row.names=names(myData), col.names=names(myData))
}

x.y.twice <- c(lm(x[1:50]~y[1:50])$residuals + mean(x[1:50]), lm(x[51:100]~y[51:100])$residuals + mean(x[51:100]), lm(x[101:150]~y[101:150])$residuals + mean(x[101:150]), lm(x[151:200]~y[151:200])$residuals + mean(x[151:200]))
var(x.y.twice)

y.x.twice <- c(lm(y[1:50]~x[1:50])$residuals + mean(y[1:50]), lm(y[51:100]~x[51:100])$residuals + mean(y[51:100]), lm(y[101:150]~x[101:150])$residuals + mean(y[101:150]), lm(y[151:200]~x[151:200])$residuals + mean(y[151:200]))
var(y.x.twice)


for (ii in 1:nLags) {
  cat("\n"); print(sprintf("Path coefficient matrix with a lag of %i (direction goes from row to column):", ii))
  lagMat <- netMatR[ii,,]%*%t(fm2$A); rownames(lagMat) <- colnames(instMatP)
  print(lagMat)
  cat("\n"); saveLagMat <- as.integer(readline("Save above path matrix (0: no; 1: yes)? "))
  if (saveLagMat) {
    cat("\n"); lagMatName <- as.character(readline("Filename prefix (e.g., PathLag1Subj1)? "))
    write.table(lagMat, file=sprintf("%s.1D", lagMatName), append=FALSE, row.names=names(myData), col.names=names(myData))
  }     
  print("-----------------") # Need to work out the t-values for the lagged effects
  #print(sprintf("Matrix of t values with a lag of %i (direction goes from row to column):", ii))
  #print(matrix(netMatT[ii,,], nrow = nROIs, ncol = nROIs, dimnames = list(names(myData), names(myData))))
  #print(sprintf("DFs = %i for null hypothesis H_0: a path coefficient = 0.", summary(fm)$varresult[[1]]$df[2]))
  #saveMatT <- as.integer(readline("Save matrix of t values for group analysis (0: no; 1: yes)? "))
  #if (saveMatT) {
  #   matName <- as.character(readline("File name prefix (e.g., TLag1Subj1)? "))
  #   write.table(netMatT[ii,,], file=sprintf("%s.1D", matName), append=FALSE, row.names=names(myData), col.names=names(myData))
  #}   	
  #print("-----------------")
}

