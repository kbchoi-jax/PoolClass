#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Four arguments must be supplied (<cntfile> <outfile> <cores> <seed>).n", call.=FALSE)
}

library(scRATE)

cntfile <- args[1]
outfile <- args[2]
nCores <- as.integer(args[3])
seed <- as.integer(args[4])

scRATE::get_bestmodel_fit(cntfile, nCores, seed, outfile)
