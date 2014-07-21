#!/usr/bin/env Rscript

# This simulation will try to mimic the Biographical Face Memory task
suppressMessages(library(niftir))
suppressMessages(library(plyr))
suppressMessages(library(neuRosim))
basedir     <- "/data/psych/faceMemoryMRI"

# Regions of Interest (should be in same space as MNI brain)
# Base things (and dimensions off of this)
roi.file1   <- file.path(basedir, "scripts/rois/facescene+facehouse_stat3_peaks_sep12.nii.gz")
roi.hdr     <- read.nifti.header(roi.file1)
roi.img     <- read.nifti.image(roi.file1)
names(roi.txt) <- c("x", "y", "z", "ind")
roi.coords  <- coords(roi.hdr$dim)[roi.img!=0,]
## below are an alternative way to get the coordinates based on the text file
## output by 3dExtrema
roi.file2   <- file.path(basedir, "scripts/rois/face_gt_house+scene.txt")
roi.txt     <- read.table(roi.file2)
alt.coords  <- roi.hdr$qto.ijk %*% t(cbind(roi.txt[,-4], rep(1, nrow(roi.txt))))
alt.coords  <- t(alt.coords[-4,]) + 1
# apply(alt.coords, 1, function(x) roi.img[x[1],x[2],x[3]])

# Sample participant (to get ntpts and task timing)
sub.file1   <- file.path(basedir, "data/nifti/tb9226/tb9226_FaceMemory01_withQ_run01.nii.gz")
sub.file2   <- file.path(basedir, "command/timing/faceMemory01_tb9226_Questions_run01_bio")
sub.file3   <- file.path(basedir, "command/timing/faceMemory01_tb9226_Questions_run01_phys")
sub.hdr     <- read.nifti.header(sub.file1)
sub.taskA   <- as.vector(as.matrix(read.table(sub.file2)))
sub.taskB   <- as.vector(as.matrix(read.table(sub.file3)))

# The above should actually just be in standard space
# but since it isn't this is a work around to get the baseline image
base.file   <- file.path(basedir, "analysis/fsl/Questions/tb9226/run01.feat/reg_standard/example_func.nii.gz")
mask.file   <- file.path(basedir, "analysis/fsl/Questions/tb9226/run01.feat/reg_standard/mask.nii.gz")
base.img    <- read.nifti.image(base.file)
mask.img    <- 1*(read.nifti.image(mask.file)!=0)


###
# PARAMETERS
###

effects     <- 1
hrf         <- "Ballon"
snr         <- 1
rho.temp    <- c(0.142, 0.108, 0.084) # AR coefficients - lag 3
rho.spat    <- 0.4 # spatial neighoor correlation
noise.wts   <- c(0.05, 0.1, 0.01, 0.09, 0.05, 0.7) # system, temporal, drift, physio, task-related, ?

dims        <- roi.hdr$dim
nscans      <- sub.hdr$dim[4]
tr          <- sub.hdr$pixdim[4]/1000
total.time  <- nscans * tr
nregions    <- nrow(roi.txt)

baseline    <- as.array(base.img)
baseline.bin<- as.array(mask.img)

## a
nums        <- rle(sub.taskA)
inds        <- which(nums$values==1)
onsetsA     <- sapply(inds, function(oi) sum(nums$lengths[1:oi]))
durationsA  <- nums$lengths[inds]
## b
nums        <- rle(sub.taskB)
inds        <- which(nums$values==1)
onsetsB     <- sapply(inds, function(oi) sum(nums$lengths[1:oi]))
durationsB  <- nums$lengths[inds]
## combined
onsets      <- list(onsetsA*tr, onsetsB*tr)
durations   <- list(durationsA*tr, durationsB*tr)

# coordinates of regions
rcoords     <- alply(roi.coords, 1, as.numeric)
radius      <- 6 # 6 voxels or 12mm

# duplicate for regions
onsets.regions <- llply(1:nregions, function(i) onsets) # note: region can have diff onsets
durations.regions <- llply(1:nregions, function(i) durations)
effects.regions <- llply(1:nregions, function(i) effects)


###
# SIM
###

# Temporal Properties
design  <- simprepTemporal(regions=nregions, onsets=onsets.regions, 
                           durations=durations.regions, hrf=hrf, TR=tr, 
                           totaltime=total.time, effectsize = effects.regions)

# Spatial Properties
spatial <- simprepSpatial(regions=nregions, coord=rcoords, radius=radius, fading=0.2, 
                          form="sphere")

# Put it together
sim.dat <- simVOLfmri(design=design, image=spatial, dim=dims, SNR=snr, 
                      noise="mixture", type="rician", rho.temp=rho.temp, 
                      rho.spat=rho.spat, w=noise.wts, nscan=nscans, vee=0.1, 
                      template=baseline.bin, base=1000, spat="gaussRF")
# not sure why getting this error:
# base should be a single number or an array with dimensions corresponding to dim


onsets.regions <- list(onsets, lapply(onsets, function(x) x+1))
durations.regions <- durations.regions[1:2]
effects.regions <- effects.regions[1:2]
design  <- simprepTemporal(regions=2, onsets=onsets.regions, 
                           durations=durations.regions, hrf=hrf, TR=tr, 
                           totaltime=total.time, effectsize = effects.regions)
ts1 <- simTSfmri(design[1], base=100, nscan=nscans, TR=tr, SNR=snr, noise="mixture", 
                type="rician", rho=rho.temp, weights=noise.wts[1:5]/sum(noise.wts[1:5]), 
                vee=0.1)
ts2 <- simTSfmri(design[2], base=100, nscan=nscans, TR=tr, SNR=snr, noise="mixture", 
                type="rician", rho=rho.temp, weights=noise.wts[1:5]/sum(noise.wts[1:5]), 
                vee=0.1)