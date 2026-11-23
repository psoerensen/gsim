testthat::test_that("shared simulation workflow validates recipes without changing caller RNG", {
  set.seed(338);before<-.Random.seed
  bad<-list(schema="gsim.simulation/v1",inputs=list(),traits=c("t1","t2"))
  testthat::expect_error(gsim_simulate(bad,tempfile()),"inputs")
  bad$unexpected<-TRUE
  testthat::expect_error(gsim_simulate(bad,tempfile()),"schema|unknown")
  testthat::expect_true(is.function(gsim_simulate))
  testthat::expect_identical(before,.Random.seed)
})


test_that("prepared Gaussian simulation shares deterministic draws and checks", {
  root <- tempfile("gsim-gaussian-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  binary <- function(name, x) {
    f <- file(file.path(root, name), "wb"); on.exit(close(f))
    writeBin(charToRaw("GSVEC001"), f)
    writeBin(as.integer(length(x)), f, size = 4L, endian = "little")
    writeBin(as.raw(rep(0, 4)), f)
    writeBin(as.double(x), f, size = 8L, endian = "little")
  }
  binary("K.bin", 2); binary("T.bin", c(0, 1))
  binary("mean.bin", c(1, 0)); binary("D.bin", c(1, 1))
  binary("R.bin", .7)
  writeLines(c("GSIM_GAUSSIAN_LABELS 1", '"r1" "t1" "id1"'), file.path(root, "labels.txt"))
  writeLines("{}", file.path(root, "provenance.json"))
  writeLines(c("GSIM_GAUSSIAN_INPUTS 1", 'latent_dimension "1"',
    'coefficient_dimension "2"', 'primary_dimension "2"', 'record_count "1"',
    'prior_covariance "K.bin"', 'coefficient_projection "T.bin"',
    'coefficient_mean "mean.bin"', 'record_design "D.bin"',
    'residual_covariance "R.bin"', 'record_labels "labels.txt"',
    'provenance "provenance.json"'), file.path(root, "spec.txt"))
  recipe <- list(schema = "gsim.gaussian_model/v1",
    inputs = list(specification = "spec.txt"), seed = 4401)
  set.seed(81); before <- .Random.seed
  gsim_simulate(recipe, file.path(root, "out"), base_directory = root)
  expect_identical(.Random.seed, before)
  b <- read.csv(file.path(root, "out", "coefficients.csv"))$value
  z <- read.csv(file.path(root, "out", "latent_innovations.csv"))$standard_normal
  e <- read.csv(file.path(root, "out", "residual_innovations.csv"))$standard_normal
  y <- read.csv(file.path(root, "out", "record_truth.csv"))$phenotype
  expect_equal(b, c(1, sqrt(2) * z), tolerance = 1e-14)
  expect_equal(y, sum(b) + sqrt(.7) * e, tolerance = 1e-14)
  expect_error(gsim_simulate(recipe, file.path(root, "out"), root), "existing")
  recipe$maximum_memory_bytes <- 1
  expect_error(gsim_simulate(recipe, file.path(root, "bad"), root), "memory")
  expect_false(dir.exists(file.path(root, "bad")))
  recipe$maximum_memory_bytes <- NULL; binary("K.bin", NaN)
  expect_error(gsim_simulate(recipe, file.path(root, "nonfinite"), root), "nonfinite")
})
