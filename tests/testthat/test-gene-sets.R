.gene_fixture <- function() data.frame(
  snp = c("s1","s1","s2","s3","s4","s5","s6","s1"),
  gene = c("a","b","a","c","d","e","f","a"))

testthat::test_that("paired scenarios have exact composition and unique members", {
  x <- gsim_gene_sets(.gene_fixture(),c("s1","s2"),3,0:2,n_sets=3,seed=51)
  testthat::expect_identical(x$genes$gene,letters[1:6])
  testthat::expect_identical(x$genes$n_snps,c(2L,1L,1L,1L,1L,1L))
  testthat::expect_identical(x$genes$n_causal_snps,c(2L,1L,0L,0L,0L,0L))
  testthat::expect_identical(x$truth$realized_size,rep(3L,9))
  testthat::expect_identical(x$truth$realized_n_causal,rep(0:2,each=3))
  testthat::expect_identical(x$truth$causal_fraction,rep((0:2)/3,each=3))
  testthat::expect_identical(names(x$sets),x$truth$set_id)
  testthat::expect_identical(x$truth$replicate,rep(1:3,3))
  for (i in seq_along(x$sets)) {
    testthat::expect_identical(anyDuplicated(x$sets[[i]]),0L)
    rows <- x$membership[x$membership$set_id==names(x$sets)[i],]
    testthat::expect_identical(rows$gene,x$sets[[i]])
    testthat::expect_identical(sum(rows$causal),x$truth$requested_n_causal[i])
  }
  broadcast <- gsim_gene_sets(.gene_fixture(),"s1",2:3,1,seed=4)
  testthat::expect_identical(broadcast$truth$requested_size,2:3)
  testthat::expect_identical(broadcast$truth$requested_n_causal,c(1L,1L))
  testthat::expect_identical(unserialize(serialize(x,NULL)),x)
})

testthat::test_that("all-null, unmapped, and reused pools retain structural truth", {
  map <- .gene_fixture()
  null <- gsim_gene_sets(map,character(0),6,0,n_sets=2,seed=1)
  testthat::expect_false(any(null$genes$causal))
  testthat::expect_identical(null$sets[[1]],null$sets[[2]])
  testthat::expect_identical(intersect(null$sets[[1]],null$sets[[2]]),letters[1:6])
  testthat::expect_warning(x <- gsim_gene_sets(map,c("z","s1","y"),3,1,seed=3),
                          "2 causal SNP\\(s\\) absent")
  testthat::expect_identical(x$unmapped_causal_snps,c("y","z"))
  testthat::expect_identical(x$genes$gene[x$genes$causal],c("a","b"))
  testthat::expect_warning(only_unknown <- gsim_gene_sets(map,"unknown",3,0,seed=3),"1 causal")
  testthat::expect_false(any(only_unknown$genes$causal))
})

testthat::test_that("canonical order and the global RNG convention are stable", {
  map <- .gene_fixture()
  x <- gsim_gene_sets(map,c("s1","s2"),3,1,n_sets=2,seed=44)
  after <- .Random.seed
  shuffled <- gsim_gene_sets(map[rev(seq_len(nrow(map))),],c("s2","s1"),3,1,n_sets=2,seed=44)
  testthat::expect_identical(x,shuffled)
  testthat::expect_identical(.Random.seed,after)
  dedup <- gsim_gene_sets(unique(map),c("s1","s2"),3,1,n_sets=2,seed=44)
  testthat::expect_identical(x,dedup)
  # Independent reference draws: causal then noncausal, without replacement.
  set.seed(44)
  manual <- lapply(1:2,function(i) sort(c(c("a","b")[sample.int(2,1)],
                                           c("c","d","e","f")[sample.int(4,2)]),method="radix"))
  testthat::expect_identical(.Random.seed,after)
  testthat::expect_identical(unname(x$sets),manual)
  set.seed(44)
  from_stream <- gsim_gene_sets(map,c("s1","s2"),3,1,n_sets=2)
  testthat::expect_identical(from_stream$sets,x$sets)
  testthat::expect_identical(.Random.seed,after)
  set.seed(44); initial <- .Random.seed
  testthat::expect_false(identical(initial,after))
  # Exact identifiers are not normalized by trimming or case folding.
  exact <- data.frame(snp=c("S","s"," s"),gene=c("G","g"," g"))
  z <- gsim_gene_sets(exact,"s",1,1,seed=1)
  testthat::expect_identical(z$sets[[1]],"g")
})

testthat::test_that("invalid identifiers and infeasible counts are rejected", {
  map <- .gene_fixture()
  testthat::expect_error(gsim_gene_sets(as.matrix(map),"s1",3,1),"data.frame")
  for (bad in list(c(NA_character_,"s1"),c("","s1"),1:2,factor("s1"))) {
    testthat::expect_error(gsim_gene_sets(map,bad,3,1),"character IDs")
  }
  testthat::expect_error(gsim_gene_sets(map,c("s1","s1"),3,1),"unique")
  for (field in c("snp","gene")) {
    for (value in list(NA_character_,"",1,factor("s1"))) {
      bad <- map; bad[[field]] <- rep(value,nrow(map))
      testthat::expect_error(gsim_gene_sets(bad,"s1",3,1),"character IDs")
    }
  }
  bad <- map; bad$gene[1] <- "intergenic"
  testthat::expect_error(gsim_gene_sets(bad,"s1",3,1),"absent mappings")
  for (value in list(NA_real_,Inf,1.5,-1,numeric(0),"3",TRUE)) {
    testthat::expect_error(gsim_gene_sets(map,"s1",value,1),"integer counts")
    testthat::expect_error(gsim_gene_sets(map,"s1",3,value),"integer counts")
    testthat::expect_error(gsim_gene_sets(map,"s1",3,1,n_sets=value),"integer counts")
  }
  testthat::expect_error(gsim_gene_sets(map,"s1",0,0),"integer counts")
  testthat::expect_error(gsim_gene_sets(map,"s1",3,1,n_sets=1:2),"one positive")
  testthat::expect_error(gsim_gene_sets(map,"s1",3,1,n_sets=0),"integer counts")
  testthat::expect_error(gsim_gene_sets(map,"s1",2:3,0:2),"equal lengths")
  testthat::expect_error(gsim_gene_sets(map,"s1",c(3,3),c(1,3)),"Scenario 2.*2 causal.*4 noncausal")
  testthat::expect_error(gsim_gene_sets(map,"s1",5,0),"Scenario 1.*2 causal.*4 noncausal")
  testthat::expect_error(gsim_gene_sets(map,"s1",1,2),"Scenario 1")
  empty <- data.frame(snp=character(),gene=character())
  testthat::expect_error(gsim_gene_sets(empty,character(),1,0),"0 causal and 0 noncausal")
})

testthat::test_that("integration preserves effects and storage scales with memberships", {
  sim <- gsim(n=64,m=24,n_causal=4,seed=812)
  effects <- sim$B
  map <- data.frame(snp=rownames(sim$B),gene=paste0("gene",seq_len(nrow(sim$B))))
  x <- gsim_gene_sets(map,sim$causal_rsids,12,c(0,2,4),n_sets=3,seed=813)
  testthat::expect_identical(sim$B,effects)
  testthat::expect_identical(sum(x$genes$causal),4L)
  testthat::expect_identical(x$truth$realized_n_causal,rep(c(0L,2L,4L),each=3))
  testthat::expect_gt(length(intersect(x$sets[[1]],x$sets[[2]])),0L)
  a <- gsim_gene_sets(map,sim$causal_rsids,4,1,n_sets=20,seed=9)
  b <- gsim_gene_sets(map,sim$causal_rsids,4,1,n_sets=40,seed=9)
  testthat::expect_identical(nrow(a$membership),80L)
  testthat::expect_identical(nrow(b$membership),160L)
  testthat::expect_identical(sum(lengths(b$sets)),160L)
  testthat::expect_false(any(vapply(b,is.matrix,logical(1))))
  testthat::expect_lt(as.numeric(object.size(b)),2.5*as.numeric(object.size(a)))
})
