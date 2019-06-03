#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("Five arguments must be supplied (<loomfile> <num_chunks> <outdir> <unique?> <dryrun?>)", call.=FALSE)
}

library(loomR)
nCores <- parallel::detectCores()

loomfile <- args[1]
num_chunks <- as.integer(args[2])
outdir <-  args[3]
unique <- as.integer(args[4])
dryrun <- as.integer(args[5])

ds <- connect(loomfile, mode = 'r+')
if(unique) {
  dmat <- ds$layers$unique[,]
  cat('[submit_jobs] Unique read counts will be used.\n')
} else {
  dmat <- ds$matrix[,]
  cat('[submit_jobs] EMASE counts will be used.\n')
}
num_cells <- dim(dmat)[1]
num_genes <- dim(dmat)[2]
gname <- ds$row.attrs$GeneID[]
cname <- ds$col.attrs$CellID[]
selected <- ds$row.attrs$`Selected:EM`[]
if(length(selected) == 0) {
  selected <- ds$row.attrs$`Selected`[]
}
cat(sprintf('[submit_jobs] %d genes will be processed.\n', sum(selected)))
ds$close_all()

idx_gsurv <- which(selected > 0)
num_gsurv <- length(idx_gsurv)
chunk_sz <- ceiling(num_gsurv / num_chunks)
gsets <- split(idx_gsurv, ceiling(seq_along(idx_gsurv) / chunk_sz))

g_ends <- c()
for (k in 1:length(gsets)) {
  g_ends <- c(g_ends, tail(gsets[[k]], 1))
}
g_starts <- g_ends + 1
g_starts <- c(1, g_starts)
g_starts <- g_starts[-length(g_starts)]
g_ends[length(g_ends)] <- num_genes

dmat <- as.data.frame(t(dmat))
rownames(dmat) <- gname
colnames(dmat) <- cname
csize <- as.vector(colSums(dmat))

scriptfile <- '/home/kbchoi/src/utils/submit_model_comparison_on_cluster.sh'
rfile <- 'run_model_comparison.R'

for (k in 1:length(gsets)) {
  s <- g_starts[k]
  e <- g_ends[k]
  cntmat <- dmat[s:e,]
  gsurv  <- selected[s:e]
  ifile <- sprintf('%s/_chunk.%05d-%05d.rds', outdir, s, e)
  ofile <- sprintf('%s/_poolclass_compare_models.%05d-%05d.rds', outdir, s, e)
  cmdstr <- sprintf('qsub -o %s -e %s -v RFILE=%s,INFILE=%s,OUTFILE=%s %s', outdir, outdir, rfile, ifile, ofile, scriptfile)
  if(!dryrun) {
    save(cntmat, gsurv, csize, file = ifile)
    cat(cmdstr, '\n')
    system(cmdstr)
    Sys.sleep(1)
  } else {
    cat(cmdstr, '\n')
    Sys.sleep(1)
  }
}
