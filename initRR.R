library(TMB)
library(RColorBrewer)
library(colorRamps)
library(dplyr)
library(reshape2)
library(parallel)
library(here)
library(mvtnorm)

source("simRR.R")   # Simulate run reconstruction
source("fitRR.R")   # Fit run reconstruction
source("fitSim.R")  # Fit run reconstruction to simmed data
source("calcRunTiming.R")
source("tools.R")   # Background functions
source("plot.R")    # Plotting functions

# Compile TMB objective function
compile("yukonChinookRunRecon.cpp")