test_that("site scanning reads fixed fields and subsets gzip VCF by callback", {
  vcf <- tempfile(fileext = ".vcf")
  subset <- tempfile(fileext = ".vcf.gz")
  header <- c("##fileformat=VCFv4.2", "##reference=GRCh38",
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1")
  records <- c("1\t10\trs10\tA\tG\t30\tPASS\tAC=1;AN=2;FLAG\tGT\t0|1",
               "1\t11\t.\tC\tA,G\t.\t.\tAC=0,1;AN=2\tGT\t0|2",
               "2\t10\trs20\tC\tG\t.\tPASS\t.\tGT\t0|0")
  writeLines(c(header, records), vcf)
  batches <- list()
  observed_header <- gsim_scan_vcf(vcf, function(batch) {
    batches[[length(batches) + 1L]] <<- batch
    batch$filter == "PASS"
  }, info = c("AC", "AN", "FLAG"), chromosome = "1", region = c(10, 11),
  batch_size = 1L, output = subset)

  expect_identical(observed_header, header)
  expect_length(batches, 2L)
  expect_identical(batches[[1L]]$info_AC, "1")
  expect_identical(batches[[1L]]$info_FLAG, "TRUE")
  expect_identical(batches[[1L]]$qual, "30")
  expect_identical(batches[[2L]]$info_AC, "0,1")
  expect_identical(batches[[2L]]$info_FLAG, NA_character_)
  connection <- gzfile(subset, "rt")
  expect_identical(readLines(connection), c(header, records[1L]))
  close(connection)
  reread <- NULL
  gsim_scan_vcf(subset, function(batch) reread <<- batch, info = "AN")
  expect_identical(reread$pos, 10)
  expect_identical(reread$info_AN, "2")

  failed <- tempfile(fileext = ".gz")
  expect_error(gsim_scan_vcf(vcf, function(batch) stop("callback failed"),
                             output = failed), "callback failed")
  expect_false(file.exists(failed))
})

test_that("site scanning supports VCFs with no sample columns", {
  vcf <- tempfile()
  writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
               "1\t1\t.\tA\tC\t.\tPASS\tAC=2"), vcf)
  site <- NULL
  gsim_scan_vcf(vcf, function(batch) site <<- batch, info = "AC", batch_size = 1L)
  expect_identical(site$info_AC, "2")
  cat("1\t2\t.\n", file = vcf, append = TRUE)
  expect_error(gsim_scan_vcf(vcf, identity), "incomplete record")
})
