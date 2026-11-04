.boa_fixture <- function(M=31L) {
  tab<-data.frame(animal=c("AS","AD","BS","BD","CS","CD","F1S","F1D","F2","BC","TC"),
    sire=c(rep(NA,6),"AS","AS","F1S","AS","CS"),
    dam=c(rep(NA,6),"BD","BD","F1D","F1D","F1D"))
  ped<-gsim_pedigree_from_table(tab)
  h1<-outer(1:6,seq_len(M),function(i,j)(i+j%/%3+j)%%2)
  h2<-outer(1:6,seq_len(M),function(i,j)(i*3+j%/%4)%%2)
  dimnames(h1)<-dimnames(h2)<-list(tab$animal[1:6],paste0("m",seq_len(M)))
  chr<-rep(c("1","2"),c(M%/%2,M-M%/%2));pos<-c(seq(0,3,length.out=M%/%2),seq(0,3,length.out=M-M%/%2))
  args<-list(pedigree=ped,founder_haplotypes=list(h1=h1,h2=h2),
    founder_breed=setNames(rep(c("A","B","C"),each=2),tab$animal[1:6]),
    chromosome=chr,genetic_position=pos,seed=811)
  list(args=args,panel=do.call(gsim_boa_genotypes,args))
}

testthat::test_that("origin phases replay exactly the allele meiosis plan", {
  x<-.boa_fixture();p<-x$panel;ids<-p$pedigree$canonical_order
  tab<-p$pedigree$pedigree[match(ids,p$pedigree$pedigree$animal),]
  set.seed(77);rng<-.Random.seed;again<-do.call(gsim_boa_genotypes,x$args)
  testthat::expect_identical(.Random.seed,rng)
  for(k in names(p$blocks)) {
    b<-p$blocks[[k]];h1<-gsim:::.gsim_packed_unpack(b$h1);h2<-gsim:::.gsim_packed_unpack(b$h2)
    testthat::expect_identical(b$ancestry_h1,again$blocks[[k]]$ancestry_h1)
    for(i in which(!is.na(tab$sire)))for(side in 1:2) {
      parent<-match(if(side==1L)tab$sire[i] else tab$dam[i],ids)
      plan<-.Call(gsim:::C_gsim_meiosis_plan,b$variant_map$genetic_position,811,ids[i],k,as.integer(side))
      phase<-(plan$starting_haplotype-1L+vapply(b$variant_map$genetic_position,function(z)sum(plan$crossovers<=z),integer(1)))%%2L
      a<-b$ancestry_h1[parent,];a[phase==1]<-b$ancestry_h2[parent,phase==1]
      h<-h1[parent,];h[phase==1]<-h2[parent,phase==1]
      testthat::expect_identical(if(side==1L)b$ancestry_h1[i,] else b$ancestry_h2[i,],a)
      testthat::expect_identical(if(side==1L)h1[i,] else h2[i,],h)
    }
  }
  testthat::expect_equal(unname(p$pedigree_fractions["TC",]),c(.25,.25,.5))
  bad<-x$args;bad$maximum_ancestry_bytes<-1
  testthat::expect_error(do.call(gsim_boa_genotypes,bad),"ancestry buffers")
})

testthat::test_that("packed BOA accumulation matches independent dense formulas", {
  x<-.boa_fixture();p<-x$panel;B<-3L;Q<-2L;M<-31L
  S<-diag(6)*.02+.005;Omega<-diag(4)*.01+.002
  state<-gsim_boa_trait_state(p,S,channels=c("level","slope"),lbp_covariance=Omega,seed=177)
  v<-gsim_boa_values(p,state);expected<-lapply(p$breeds,function(b)matrix(0,11,2));lbp<-matrix(0,11,2)
  for(block in p$blocks) {
    h1<-matrix(as.integer(gsim:::.gsim_packed_unpack(block$h1)),11)
    h2<-matrix(as.integer(gsim:::.gsim_packed_unpack(block$h2)),11)
    a1<-matrix(as.integer(block$ancestry_h1),11);a2<-matrix(as.integer(block$ancestry_h2),11)
    at<-block$variant_ids
    for(b in seq_len(B)) {
      z<-(a1==b-1)*(h1-rep(block$frequencies[,b],each=11))+(a2==b-1)*(h2-rep(block$frequencies[,b],each=11))
      cols<-(b-1)*Q+seq_len(Q)
      expected[[b]]<-expected[[b]]+z%*%state$effects[at,cols]
      zlbp<-(a1==b-1)+(a2==b-1)-2*p$pedigree_fractions[,b]
      lbp<-lbp+zlbp%*%state$lbp_effects[at,cols]
    }
  }
  for(b in 1:3)testthat::expect_equal(unname(v$breed[[b]]),unname(expected[[b]]),tolerance=1e-12)
  testthat::expect_equal(unname(v$lbp),unname(lbp),tolerance=1e-12)
  testthat::expect_equal(unname(v$lbp[c("AS","F1S","F1D"),]),matrix(0,3,2))
  testthat::expect_equal(v$genetic,v$snp+v$lbp)
  # Congruent change of LBP reference preserves every individual value.
  other<-state;theta<-array(state$lbp_effects,c(M,Q,B))
  theta<-sweep(theta,c(1,2),theta[,,1],"-")
  other$lbp_effects[,]<-matrix(theta,M,B*Q);other$reference_breed<-"A"
  testthat::expect_equal(gsim_boa_values(p,other)$genetic,v$genetic,tolerance=1e-12)
  testthat::expect_identical(gsim_boa_trait_state(p,S,channels=c("level","slope"),lbp_covariance=Omega,seed=177)$effects,state$effects)
  markers<-rownames(state$effects);components<-setNames(rep(c("one","two"),length.out=M),markers)
  split<-gsim_boa_trait_state(p,list(one=S,two=S),channels=c("level","slope"),marker_components=components,seed=177)
  testthat::expect_identical(split$effects,state$effects)
  no_lbp<-gsim_boa_trait_state(p,S,channels=c("level","slope"),seed=177)
  testthat::expect_equal(gsim_boa_values(p,no_lbp)$lbp,matrix(0,11,2,dimnames=dimnames(v$lbp)))
})
