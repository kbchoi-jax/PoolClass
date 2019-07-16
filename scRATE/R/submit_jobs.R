submit_jobs <- function(loomfile, num_chunks, outdir, dryrun, scriptfile, rfile, layer=NULL, g_start=NULL, g_end=NULL, chunk_start=NULL, chunk_end=NULL) {
  nCores <- min(4, parallel::detectCores())
  ds <- connect(loomfile, mode = 'r+')
  if(is.null(layer)) {
    dmat <- ds$matrix[,]
    cat('[submit_jobs] Counts from main layer will be loaded.\n')
  } else {
    dmat <- ds$layers[[layer]][,]
    cat(sprintf('[submit_jobs] Counts from %s layer will be loaded.\n', layer))
  }
  num_cells <- dim(dmat)[1]
  num_genes <- dim(dmat)[2]
  gname <- ds$row.attrs$GeneID[]
  cname <- ds$col.attrs$CellID[]
  selected <- ds$row.attrs$`Selected:EM`[]
  if(length(selected) == 0) {
    selected <- ds$row.attrs$`Selected`[]
  }
  ds$close_all()

  if(is.null(g_start)) {
    gidx1 <- 1
  } else {
    gidx1 <- g_start
  }
  if(is.null(g_end)) {
    gidx2 <- num_genes
  } else {
    gidx2 <- g_end
  }

  idx_gsurv <- which(selected > 0)
  idx_gsurv <- idx_gsurv[idx_gsurv >= gidx1 & idx_gsurv <= gidx2]
  num_gsurv <- length(idx_gsurv)
  chunk_sz <- ceiling(num_gsurv / num_chunks)
  gsets <- split(idx_gsurv, ceiling(seq_along(idx_gsurv) / chunk_sz))

  cat(sprintf('[submit_jobs] %d genes (between Gene %d and %d) will be processed.\n', num_gsurv, gidx1, gidx2))

  g_ends <- c()
  for (k in 1:length(gsets)) {
    g_ends <- c(g_ends, tail(gsets[[k]], 1))
  }
  g_starts <- g_ends + 1
  g_starts <- c(gidx1, g_starts)
  g_starts <- g_starts[-length(g_starts)]
  g_ends[length(g_ends)] <- gidx2

  dmat <- as.data.frame(t(dmat))
  rownames(dmat) <- gname
  colnames(dmat) <- cname
  csize <- as.vector(colSums(dmat))

  #scriptfile <- '/home/kbchoi/src/utils/submit_model_comparison_on_cluster.sh'
  #rfile <- 'run_model_comparison.R'

  if (is.null(chunk_start)) {
    cidx1 <- 1
  } else {
    cidx1 <- chunk_start
  }
  if (is.null(chunk_end)) {
    cidx2 <- length(gsets)
    if (num_chunks != cidx2) {
      stop(sprintf("The computed number of chunks %d is not equal to the requested number %d.", num_chunks, cidx2), call. = FALSE)
    }
  } else {
    if (chunk_end > num_chunks) {
      stop(sprintf("The computed number of chunks %d is not equal to the requested number %d.", num_chunks, chunk_end), call. = FALSE)
    }
    cidx2 <- chunk_end
  }

  cat(sprintf('[submit_jobs] Chunk %d to %d (out of %d) will be processed.\n', cidx1, cidx2, num_chunks))

  for (k in cidx1:cidx2) {
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
    }
  }
}
