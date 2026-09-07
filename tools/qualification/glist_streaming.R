# Bounded synthetic resource probe; run from the repository root with the
# isolated installed gsim loaded. No external data, qgg, or dense panel.
# source('tools/qualification/glist_streaming.R'); result <- glist_streaming()
glist_streaming <- function(output = NULL) {
  stopifnot(requireNamespace('gsim', quietly = TRUE))
  root <- tempfile('gsim-resource-'); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  n <- 1024L; m <- 4096L; nt <- 3L; block <- 32L
  ids <- paste0('id', seq_len(n)); markers <- paste0('m', seq_len(m))
  prefix <- file.path(root, 'synthetic')
  con <- file(paste0(prefix, '.bed'), 'wb')
  tryCatch({
    writeBin(as.raw(c(108,27,1)), con)
    # Independently encoded one marker at a time. Synthetic dosage patterns
    # are reproducible, deliberately not biological LD or a HAPNEST fixture.
    set.seed(9107)
    for (j in seq_len(m)) {
      x <- rbinom(n, 2, .1 + .7 * (j %% 101) / 100)
      code <- c(3L,2L,0L)[x+1L]
      writeBin(as.raw(colSums(matrix(code * rep(c(1,4,16,64),n/4),4))), con)
    }
  }, finally = close(con))
  write.table(data.frame(1,markers,0,seq_len(m),'G','A'),paste0(prefix,'.bim'),
              row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(data.frame('fam',ids,0,0,0,-9),paste0(prefix,'.fam'),
              row.names=FALSE,col.names=FALSE,quote=FALSE)
  gl <- list(ids=ids,n=n,rsids=list(markers),mchr=m,bedfiles=paste0(prefix,'.bed'),
             bimfiles=paste0(prefix,'.bim'),famfiles=paste0(prefix,'.fam'))
  rows <- lapply(c(64L,512L,4096L), function(causal) {
    q <- setNames(as.numeric(seq_len(m) <= causal),markers)
    gc()
    log <- file.path(root,'allocations.log')
    Rprofmem(log)
    elapsed <- tryCatch(system.time({
      sim <- gsim::gsim(Glist=gl,causal_probability=q,nt=nt,seed=552,
                         chunk_size=block,compute_sumstats=FALSE)
    })[['elapsed']], finally=Rprofmem(NULL))
    allocations <- suppressWarnings(as.numeric(sub(' .*','',readLines(log))))
    largest <- max(allocations,na.rm=TRUE)
    stopifnot(sim$settings$genotype_stream$backend == 'native_bed_scalar',
              sim$settings$n_causal == causal,is.null(sim$W_causal))
    # At the larger causal counts a single full causal matrix would exceed
    # every observed R allocation. Native code retains one packed record.
    if (causal >= 512L) stopifnot(largest < 8*n*causal)
    data.frame(n=n,m=m,causal=causal,traits=nt,block=block,elapsed_seconds=elapsed,
      decode_capacity_bytes=8*n*block,packed_record_bytes=ceiling(n/4),
      full_causal_payload_avoided=8*n*causal,
      full_B_bytes=as.numeric(object.size(sim$B)),
      causal_B_bytes=as.numeric(object.size(sim$B_causal)),
      G_E_Y_bytes=sum(vapply(sim[c('G','E','Y')],object.size,numeric(1))),
      input_metadata_bytes=as.numeric(object.size(gl)),
      result_bytes=as.numeric(object.size(sim)),largest_R_allocation_bytes=largest)
  })
  result <- do.call(rbind,rows)
  if (!is.null(output)) write.csv(result,output,row.names=FALSE)
  print(result,row.names=FALSE)
  invisible(result)
}
