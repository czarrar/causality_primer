#' This script will attempt to run 1dGC.R although not this requires user intervention.
#+ run
datadir <- "~/data/fsl_sims/sim01"
setwd(datadir)

program <- system("which 1dGC.R", intern=T)
source(program)
