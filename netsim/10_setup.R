#!/usr/bin/env Rscript

# This script will first download the netsim data (if it isn't already there)
# and then convert the time-series and ground-truth matrices from matlab format to text files
# * ts: Ntimepoints x Nnodes (contains all subjects' timeseries concatenated)
# * net: Nsubjects x Nnodes x Nnodes (contains the ground truth networks)

library(R.matlab)

# Create output data directory
datadir <- "data"
if (!file.exists(datadir)) dir.create(datadir)

# Download the data and extract
download.file("http://fsl.fmrib.ox.ac.uk/analysis/netsim/sims.tar.gz", "data/sims.tar.gz", method="auto")
system(sprintf("cd %s; tar -xzvf sims.tar.gz", datadir))
file.remove(file.path(datadir, "sims.tar.gz"))

# First we move any 1-9 files
for (i in 1:9) {
    infile  <- file.path(datadir, sprintf("sim%i.mat", i))
    outfile <- file.path(datadir, sprintf("sim%02i.mat", i))
    if (file.exists(infile)) {
        cat(sprintf("%s -> %s", infile, outfile), "\n")
        file.rename(infile, outfile)
    }
}

# Second we read in each of the simulations and save them in text file format
for (i in 2:28) {
    cat("simulation", i, "\n")
    
    ## Paths
    cat("...paths\n")
    infile  <- file.path(datadir, sprintf("sim%02i.mat", i))
    outfile <- file.path(datadir, sprintf("sim%02i.rda", i))
    outdir  <- file.path(datadir, sprintf("sim%02i", i))
    if (!file.exists(outdir)) dir.create(outdir)
    
    ## Read in data
    cat("...read data\n")
    mlab    <- readMat(infile)
    nsubs   <- as.integer(mlab$Nsubjects)
    nnodes  <- as.integer(mlab$Nnodes)
    ntpts   <- as.integer(mlab$Ntimepoints)
    ts      <- mlab$ts
    net     <- mlab$net
    
    ## Reshape
    cat("...reshape\n")
    dim(ts) <- c(nsubs, ntpts, nnodes)
    
    ## Save
    cat("...save\n")
    save(nsubs, nnodes, ntpts, ts, net, file=outfile)
    for (j in 1:nsubs) {
        sub.ts  <- ts[j,,]
        tsfile  <- file.path(outdir, sprintf("sub%02i_ts.txt", j))
        write.table(sub.ts, file=tsfile, row.names=F, col.names=F, quote=F)
        
        sub.net <- net[j,,]
        netfile <- file.path(outdir, sprintf("sub%02i_net.txt", j))
        write.table(sub.net, file=netfile, row.names=F, col.names=F, quote=F)
    }
}
