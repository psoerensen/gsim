testthat::test_that("explicit dominance and environment slopes match independent formulas", {
  W<-matrix(c(0,1,2,1,2,0,1,2,0,2,1,0),4,3,
    dimnames=list(paste0("i",1:4),paste0("m",1:3)))
  A<-matrix(c(1,-2,.5),3,1,dimnames=list(colnames(W),"D1"))
  D<-A*.3;I<-A*-.2;env<-setNames(c(-1,0,2,1),rownames(W))
  fit<-gsim(W=W,architecture="fixed",beta=A,dominance=D,interaction=I,
    environment=env[4:1],standardize_W=FALSE,scale_effects=FALSE,seed=7)
  testthat::expect_equal(fit$G_additive,W%*%A)
  testthat::expect_equal(fit$G_dominance,(W==1)%*%D)
  testthat::expect_equal(fit$G_slope,W%*%I)
  testthat::expect_equal(fit$G,W%*%A+(W==1)%*%D+(W%*%I)*as.double(env))
  testthat::expect_equal(fit$Y,fit$G+fit$E)
  frozen<-gsim_trait_state(fit)
  testthat::expect_equal(frozen$dominance,fit$dominance)
  testthat::expect_equal(frozen$interaction,fit$interaction)
  testthat::expect_equal(fit$G_interaction[2,],0)
  reverse<-gsim(W=2-W,architecture="fixed",beta=A,dominance=D,
    standardize_W=FALSE,scale_effects=FALSE,seed=7)
  testthat::expect_equal(reverse$G_dominance,fit$G_dominance)
  pure<-gsim(W=W,architecture="fixed",beta=A*0,dominance=D,
    standardize_W=FALSE,scale_effects=FALSE,seed=7)
  testthat::expect_equal(pure$G,(W==1)%*%D)
  testthat::expect_true(all(pure$causal$component==1L))
  calibrated<-gsim(W=W,architecture="fixed",beta=A,dominance=D,interaction=I,
    environment=env,standardize_W=FALSE,vg=2,seed=7)
  testthat::expect_equal(var(as.double(calibrated$G)),2)
  testthat::expect_equal(calibrated$G,(W%*%calibrated$B)+(W==1)%*%calibrated$dominance+
    (W%*%calibrated$interaction)*as.double(env))
  testthat::expect_error(gsim(W=W,beta=A,architecture="fixed",interaction=I),"together")
  testthat::expect_error(gsim(W=W,beta=A,architecture="fixed",environment=1),"together")
  testthat::expect_error(gsim(W=W,beta=A,architecture="fixed",interaction=I,environment=1:4),"name")
  bad<-W;bad[1,1]<-.5
  testthat::expect_error(gsim(W=bad,beta=A,architecture="fixed",dominance=D),"hard calls")
})

testthat::test_that("bounded BED extensions agree with dense truth including missing calls", {
  root<-tempfile("extensions-");dir.create(root);on.exit(unlink(root,recursive=TRUE))
  fixture<-.stream_fixture(root)
  # Persistent metadata readers use canonical tab-delimited triplets.
  for(path in c(fixture$Glist$famfiles,fixture$Glist$bimfiles)) {
    metadata<-read.table(path,stringsAsFactors=FALSE)
    con<-file(path,"wb")
    write.table(metadata,con,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    close(con)
  }
  W<-fixture$W[,1:11,drop=FALSE]
  gl<-fixture$Glist;gl$rsids<-gl$rsids[1];gl$rsids[[1]]<-colnames(W)
  gl$bedfiles<-gl$bedfiles[1];gl$bimfiles<-gl$bimfiles[1];gl$famfiles<-gl$famfiles[1];gl$mchr<-11L
  gl$rsidsLD<-NULL;gl$maf<-NULL;gl$af<-NULL
  ids<-fixture$ids
  A<-matrix(seq_len(22)/17,11,2,dimnames=list(colnames(W),c("D1","D2")))
  D<-A*.13;I<-A*-.09;env<-setNames(seq(-1,1,length.out=length(ids)),ids)
  args<-list(architecture="fixed",beta=A,dominance=D,interaction=I,environment=env,
    nt=2,rg=.2,standardize_W=FALSE,scale_effects=FALSE,seed=27)
  dense<-do.call(gsim,c(list(W=W[ids,,drop=FALSE]),args))
  bed<-do.call(gsim,c(list(Glist=gl,ids=ids),args))
  for(name in c("G","G_additive","G_dominance","G_slope","G_interaction","Y"))
    testthat::expect_equal(bed[[name]],dense[[name]],tolerance=1e-12)
  H<-1*(W[ids,,drop=FALSE]==1)
  for(j in seq_len(ncol(H)))H[is.na(H[,j]),j]<-mean(H[,j],na.rm=TRUE)
  testthat::expect_equal(bed$G_dominance,H%*%D)
  state<-gsim_trait_state(A,center=1,dominance=D,interaction=I)
  values<-gsim_genetic_values(sub("[.]bed$","",gl$bedfiles),state,ids,env)
  raw<-W[ids,,drop=FALSE];raw[is.na(raw)]<-1
  Hd<-1*(W[ids,,drop=FALSE]==1);Hd[is.na(Hd)]<-.5
  testthat::expect_equal(matrix(as.double(values),nrow(raw)),unname((raw-1)%*%A+Hd%*%D+(raw%*%I)*as.double(env)))
  parts<-attr(values,"genetic_components")
  testthat::expect_equal(parts$dominance,Hd%*%D)
  testthat::expect_identical(state$checksum,gsim_trait_state(A,1,D,I)$checksum)
  testthat::expect_error(gsim_genetic_values(sub("[.]bed$","",gl$bedfiles),state),"environment")
})

testthat::test_that("record-specific environments preserve independent noise and support repeated records", {
  G<-matrix(c(1,3,2,5),2,2,dimnames=list(c("a","b"),c("t1","t2")))
  S<-G*.2
  records<-data.frame(record=paste0("r",1:4),animal=c("a","a","b","b"),trait=c("t1","t1","t2","t2"))
  env<-setNames(c(-2,0,1,3),records$record)
  base<-gsim_records(G,records,seed=12)
  fit<-gsim_records(G,records,seed=12,interaction_values=S,environment=env[4:1])
  expected<-c(S[1,1],S[1,1],S[2,2],S[2,2])*as.double(env)
  testthat::expect_equal(fit$truth$interaction,expected)
  testthat::expect_equal(fit$truth$genetic,base$truth$genetic+expected)
  testthat::expect_identical(fit$truth$residual,base$truth$residual)
  testthat::expect_equal(fit$records$value,base$records$value+expected)
  testthat::expect_error(gsim_records(G,records,interaction_values=S),"together")
})

.mutation_fixture<-function(root,markers=1024L) {
  ids<-c("s","d");vcf<-file.path(root,"input.vcf")
  header<-paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",ids),collapse="\t")
  lines<-vapply(seq_len(markers),function(j)paste(c("1",j,paste0("m",j),"A","G",".",".",".","GT","0|0","0|0"),collapse="\t"),character(1))
  writeLines(c("##fileformat=VCFv4.2",header,lines),vcf)
  reference<-gsim_import_vcf(vcf,data.frame(chromosome="1",variant_id=paste0("m",seq_len(markers)),
    genetic_position_cm=seq(0,100,length.out=markers)),file.path(root,"parents"))
  pedigree<-gsim_pedigree_from_table(data.frame(animal=c(ids,"c"),sire=c(NA,NA,"s"),dam=c(NA,NA,"d"),sex=c("M","F","U")))
  list(reference=reference,pedigree=pedigree)
}
.mutation_decode<-function(x) .Call(gsim:::C_gsim_bed_read_selected,x$paths[["bed"]],
  as.integer(x$individual_count),as.integer(x$variant_count),seq_len(x$individual_count),1:64)

testthat::test_that("offspring mutation preserves founders and crossover streams with exact endpoints", {
  root<-tempfile("mutation-");dir.create(root);on.exit(unlink(root,recursive=TRUE))
  f<-.mutation_fixture(root)
  run<-function(rate,name)gsim:::.gsim_simulate_pedigree_impl(f$reference,f$pedigree,71,
    file.path(root,name),"bed",return_crossovers=TRUE,mutation_rate=rate)
  set.seed(37);before<-.Random.seed
  zero<-run(0,"zero");one<-run(1,"one");mid<-run(.25,"mid");again<-run(.25,"again")
  testthat::expect_identical(.Random.seed,before)
  testthat::expect_equal(.mutation_decode(zero),matrix(0,3,64))
  X<-.mutation_decode(one)
  testthat::expect_equal(X[match("c",one$sample_ids),],rep(2,64))
  testthat::expect_equal(X[match(c("s","d"),one$sample_ids),],matrix(0,2,64))
  testthat::expect_identical(zero$crossover_audit,one$crossover_audit)
  testthat::expect_identical(.mutation_decode(mid),.mutation_decode(again))
  testthat::expect_identical(mid$simulation$mutation_rate,.25)
  # Independent Binomial(2*1024,.25) expectation on homozygous-zero parents.
  total<-0
  for(start in seq.int(1L,1024L,64L)) {
    values<-.Call(gsim:::C_gsim_bed_read_selected,mid$paths[["bed"]],3L,1024L,
      as.integer(match("c",mid$sample_ids)),seq.int(start,start+63L))
    total<-total+sum(values)
  }
  testthat::expect_lt(abs(total-512),6*sqrt(2048*.25*.75))
  cohort<-gsim_simulate_cohort(f$reference,data.frame(animal="c",sire="s",dam="d"),71,file.path(root,"cohort"),mutation_rate=.25)
  testthat::expect_identical(as.double(.mutation_decode(cohort$bed)),as.double(.mutation_decode(mid)[match("c",mid$sample_ids),]))
  testthat::expect_error(run(-.1,"bad"),"mutation_rate")
  testthat::expect_error(run(NA_real_,"bad"),"mutation_rate")
})

testthat::test_that("mutations transmit to later descendants without changing breed-origin truth", {
  ped<-gsim_pedigree_from_table(data.frame(animal=c("s","d","c","g"),
    sire=c(NA,NA,"s","c"),dam=c(NA,NA,"d","d")))
  H<-matrix(0,2,16,dimnames=list(c("s","d"),paste0("m",1:16)))
  args<-list(pedigree=ped,founder_haplotypes=list(h1=H,h2=H),
    founder_breed=c(s="A",d="B"),chromosome=rep("1",16),genetic_position=seq(0,2,length.out=16),seed=71)
  base<-do.call(gsim_boa_genotypes,args)
  mutant<-do.call(gsim_boa_genotypes,c(args,list(mutation_rate=1)))
  b<-mutant$blocks[[1]]
  h1<-gsim:::.gsim_packed_unpack(b$h1);h2<-gsim:::.gsim_packed_unpack(b$h2)
  testthat::expect_true(all(h1["c",]==as.raw(1)))
  testthat::expect_true(all(h2["c",]==as.raw(1)))
  testthat::expect_true(all(h1["g",]==as.raw(0)))
  testthat::expect_true(all(h2["g",]==as.raw(1)))
  testthat::expect_identical(b$ancestry_h1,base$blocks[[1]]$ancestry_h1)
  testthat::expect_identical(b$ancestry_h2,base$blocks[[1]]$ancestry_h2)
  testthat::expect_identical(mutant$settings$mutation_rate,1)
})
