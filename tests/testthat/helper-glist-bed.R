# Small independent PLINK encoder. Never uses gsim's writer or decoder.
.stream_write_bed <- function(W, prefix, chromosome = "1") {
  con <- file(paste0(prefix, ".bed"), "wb")
  on.exit(close(con))
  writeBin(as.raw(c(108, 27, 1)), con)
  for (j in seq_len(ncol(W))) {
    x <- W[, j]
    code <- c(3L, 2L, 0L)[x + 1L]
    code[is.na(x)] <- 1L
    code <- c(code, rep(0L, (-length(code)) %% 4L))
    packed <- colSums(matrix(code * rep(c(1, 4, 16, 64), length.out = length(code)), 4L))
    writeBin(as.raw(packed), con)
  }
  write.table(data.frame(chromosome, colnames(W), 0, seq_len(ncol(W)), "G", "A"),
              paste0(prefix, ".bim"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(data.frame("fam", rownames(W), 0, 0, 0, -9),
              paste0(prefix, ".fam"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

.stream_fixture <- function(root) {
  dir.create(root, showWarnings = FALSE)
  set.seed(8401)
  W <- matrix(rbinom(32L * 24L, 2, rep(seq(.1, .8, length.out = 24L), each = 32L)), 32L)
  dimnames(W) <- list(paste0("id", 1:32), paste0("m", 1:24))
  W[c(1, 7, 19), 3] <- NA
  W[, 12] <- 1
  W[, 24] <- NA
  prefixes <- file.path(root, c("chr1", "chr2"))
  .stream_write_bed(W[, 1:12], prefixes[1], "1")
  .stream_write_bed(W[, 13:24], prefixes[2], "2")
  gl <- list(ids = rownames(W), n = nrow(W), rsids = list(colnames(W)[1:12], colnames(W)[13:24]),
             mchr = c(12L, 12L), bedfiles = paste0(prefixes, ".bed"),
             bimfiles = paste0(prefixes, ".bim"), famfiles = paste0(prefixes, ".fam"),
             af = list(rep(.31,12),rep(.27,12)), maf = list(rep(.31,12),rep(.27,12)))
  gl$rsidsLD <- list(gl$rsids[[1]][c(9,3,1,7)], gl$rsids[[2]][c(8,2,5,1)])
  list(W = W, Glist = gl, ids = rownames(W)[c(19,7,1,30,5,3,28,4,9,11,13,15,17,20,22,24)])
}
