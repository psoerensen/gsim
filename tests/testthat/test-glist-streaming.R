testthat::test_that("native BED accumulation matches independent dense phenotypes", {
  root <- tempfile("gsim-stream-parity-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  f <- .stream_fixture(root)
  markers <- unlist(f$Glist$rsidsLD, use.names = FALSE)
  W <- f$W[f$ids, markers, drop = FALSE]
  configs <- list(
    list(architecture = "bayesr", n_causal = 6L, nt = 1L),
    list(architecture = "bayesr", causal_probability = 1, nt = 2L,
         rg = matrix(c(1,.3,.3,1),2), re = matrix(c(1,-.2,-.2,1),2)),
    list(architecture = "maf_dependent", n_causal = 7L, nt = 2L, scale_effects = FALSE),
    list(architecture = "fixed", beta = setNames(c(.8,-.2), c("m3","m20")),
         standardize_W = FALSE),
    list(architecture = "bayesr", causal_probability = 1, a = -.4, b = -.2,
         ld_score = setNames(seq(1,2,length.out=length(markers)),markers),
         maf = setNames(rep(.31,length(markers)),markers)),
    list(architecture = "bayesc", n_causal = 6L, nt = 2L,
         marker_multipliers = setNames(seq(.5,2,length.out=length(markers)),markers),
         scale_effects = FALSE)
  )
  for (config in configs) {
    dense <- do.call(gsim, c(list(W = W, seed = 444, compute_sumstats = TRUE), config))
    dense_rng <- .Random.seed
    raw_config <- config; raw_config$scale_effects <- FALSE
    raw_dense <- do.call(gsim,c(list(W=W,seed=444),raw_config))
    raw_native <- do.call(gsim,c(list(Glist=f$Glist,ids=f$ids,seed=444),raw_config))
    testthat::expect_identical(raw_native$B,raw_dense$B)
    for (block in c(1L, 3L, 64L)) {
      stream <- do.call(gsim, c(list(Glist = f$Glist, ids = f$ids,
                                    seed = 444, chunk_size = block, compute_sumstats = TRUE), config))
      testthat::expect_identical(.Random.seed, dense_rng)
      for (field in c("component", "causal_rsids", "marker_probabilities", "causal_probability", "marker_multipliers"))
        testthat::expect_identical(stream[[field]], dense[[field]])
      # Scalar marker-order sums versus BLAS; 1e-12 relative tolerance is
      # conservative for these <= 8-term dot products and their calibration.
      for (field in c("B", "B_causal", "G", "E", "Y", "vg_observed", "Sigma_e", "rg_observed", "re_observed", "h2_observed"))
        testthat::expect_equal(stream[[field]], dense[[field]], tolerance = 1e-12)
      if (isFALSE(config$scale_effects)) testthat::expect_identical(stream$B, dense$B)
      align_stats <- function(x) x[order(x$trait, x$rsid), c("rsid","trait","beta","se","z","n")]
      a <- align_stats(stream$sumstats); b <- align_stats(dense$sumstats)
      rownames(a) <- rownames(b) <- NULL
      testthat::expect_equal(a, b, tolerance = 1e-12)
      testthat::expect_identical(stream$settings$genotype_stream$backend, "native_bed_scalar")
      testthat::expect_false(stream$settings$genotype_stream$full_causal_matrix)
      testthat::expect_null(stream$W_causal)
    }
  }
})

testthat::test_that("BED physical offsets survive filtered and reordered metadata", {
  root <- tempfile("gsim-stream-offsets-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  f <- .stream_fixture(root)
  g <- f$Glist
  g$rsids <- list(g$rsids[[1]][c(8,3,1)],g$rsids[[2]][c(9,5,2)])
  g$rsidsLD <- NULL; g$mchr <- lengths(g$rsids)
  # A second file may have its own FAM row order. IDs, not Glist row positions,
  # define the physical rows. The old qgg decoder did not handle this mismatch.
  .stream_write_bed(f$W[32:1,13:24], sub("[.]bed$","",g$bedfiles[2]), "2")
  selected <- c("m21","m1","m3","m14","m8","m17")
  native <- gsim(Glist=g, ids=f$ids, rsids=selected, causal_probability=1,
                 seed=445, scale_effects=FALSE, nt=2)
  dense <- gsim(W=f$W[f$ids,selected], causal_probability=1,
                seed=445, scale_effects=FALSE, nt=2)
  testthat::expect_identical(native$B,dense$B)
  testthat::expect_equal(native$G,dense$G,tolerance=1e-12)
  testthat::expect_identical(rownames(native$G),f$ids)
  testthat::expect_identical(native$causal_rsids,selected)
})

testthat::test_that("missing, invariant and malformed BED inputs fail explicitly", {
  root <- tempfile("gsim-stream-errors-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  f <- .stream_fixture(root)
  run <- function(...) gsim(Glist=f$Glist, ids=f$ids, causal_probability=1, seed=4, ...)
  testthat::expect_error(run(rsids="m24"),"only missing")
  testthat::expect_error(run(rsids=c("m3","m12")),"zero/non-finite")
  x <- run(rsids=c("m3","m12"),standardize_W=FALSE)
  testthat::expect_true(all(is.finite(x$G)))
  testthat::expect_error(run(return_genotypes=TRUE),"return_genotypes")
  bad <- f$Glist; bad$rsids[[1]][1] <- "absent"; bad$rsidsLD <- NULL
  testthat::expect_error(gsim(Glist=bad, n_causal=2, seed=4),"absent from BIM")
  bad <- f$Glist; bad$ids[2] <- bad$ids[1]
  testthat::expect_error(gsim(Glist=bad, n_causal=2, seed=4),"unique")
  bad <- f$Glist; bad$a1 <- list(rep("A",12),rep("A",12))
  testthat::expect_error(gsim(Glist=bad,n_causal=2,seed=4),"orientation")
  bad <- f$Glist; bad$bedfiles <- NULL
  testthat::expect_error(gsim(Glist=bad,n_causal=2,seed=4),"requires SNP-major")
  con <- file(f$Glist$bedfiles[1],"r+b"); writeBin(as.raw(c(108,27,0)),con); close(con)
  testthat::expect_error(run(),"header|mode")
  con <- file(f$Glist$bedfiles[1],"wb"); writeBin(as.raw(c(108,27,1,0)),con); close(con)
  testthat::expect_error(run(),"file size mismatch")
})

testthat::test_that("explicit custom readers have capped requests and no dense return", {
  root <- tempfile("gsim-stream-custom-"); dir.create(root)
  on.exit(unlink(root, recursive=TRUE),add=TRUE)
  f <- .stream_fixture(root); calls <- integer()
  reader <- function(Glist,rsids,ids,chr=NULL,impute=TRUE,scale=FALSE) {
    calls <<- c(calls,length(rsids)); f$W[ids,rsids,drop=FALSE]
  }
  custom <- gsim(Glist=f$Glist,ids=f$ids,seed=44,causal_probability=1,
                 getG_fun=reader,chunk_size=3,compute_sumstats=TRUE)
  native <- gsim(Glist=f$Glist,ids=f$ids,seed=44,causal_probability=1,chunk_size=3)
  testthat::expect_lte(max(calls),3L)
  testthat::expect_equal(custom$G,native$G,tolerance=1e-12)
  testthat::expect_identical(custom$settings$genotype_stream$backend,"bounded_custom_reader")
})


testthat::test_that("native boundary checks and legacy summary ordering survive the cap", {
  root <- tempfile("gsim-stream-cap-"); dir.create(root)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  set.seed(74)
  W <- matrix(rbinom(32*130,2,.4),32,dimnames=list(paste0("s",1:32),paste0("v",1:130)))
  prefix <- file.path(root,"chr1"); .stream_write_bed(W,prefix)
  gl <- list(ids=rownames(W),n=32L,rsids=list(colnames(W)),mchr=130L,bedfiles=paste0(prefix,".bed"))
  symbol <- getFromNamespace("C_gsim_bed_read_selected","gsim")
  testthat::expect_error(.Call(symbol,gl$bedfiles,32L,130L,0L,1L),"sample index")
  testthat::expect_error(.Call(symbol,gl$bedfiles,32L,130L,1L,131L),"marker index")
  testthat::expect_error(.Call(symbol,gl$bedfiles,32L,130L,1L,rep(1L,65)),"64 columns")
  sim <- gsim(Glist=gl,n_causal=90,nt=2,seed=78,compute_sumstats=TRUE)
  dense <- gsim(W=W,n_causal=90,nt=2,seed=78,compute_sumstats=TRUE)
  # One legacy default chunk spans >64 markers; retain exact IDs and row names.
  testthat::expect_identical(sim$sumstats[c("rsid","trait","n")],dense$sumstats[c("rsid","trait","n")])
  testthat::expect_equal(sim$sumstats,dense$sumstats,tolerance=1e-12)
})
