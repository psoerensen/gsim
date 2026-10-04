.breeding_fixture <- function(root) {
  ids<-c("s1","s2","d1","d2")
  vcf<-file.path(root,"input.vcf")
  header<-paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",ids),collapse="\t")
  lines<-c("##fileformat=VCFv4.2",header,
    "1\t10\tm1\tA\tG\t.\t.\t.\tGT\t0|0\t1|1\t0|1\t1|0",
    "1\t20\tm2\tC\tT\t.\t.\t.\tGT\t0|1\t1|0\t0|0\t1|1",
    "2\t10\tm3\tG\tA\t.\t.\t.\tGT\t1|1\t0|0\t1|0\t0|1")
  writeLines(lines,vcf,useBytes=TRUE)
  reference<-gsim_import_vcf(vcf,data.frame(chromosome=c("1","1","2"),
    variant_id=paste0("m",1:3),genetic_position_cm=c(0,100,0)),file.path(root,"parents"))
  ped<-gsim_pedigree_from_table(data.frame(animal=ids,sire=NA_character_,dam=NA_character_,sex=c("M","M","F","F"),herd=c("a","b","a","b")))
  list(ids=ids,reference=reference,pedigree=ped,
    W=matrix(c(0,2,1,1,1,1,0,2,2,0,1,1),4,3,dimnames=list(ids,paste0("m",1:3))))
}

testthat::test_that("supplied parentage and selection rebuild auditable contracts", {
  tab<-data.frame(animal=c("c","s","d"),sire=c("s",NA,NA),dam=c("d",NA,NA),sex=c("U","M","F"),herd="H")
  ped<-gsim_pedigree_from_table(tab)
  testthat::expect_identical(ped$canonical_order,c("s","d","c"))
  testthat::expect_identical(ped$pedigree$generation,c(2L,1L,1L))
  testthat::expect_identical(ped$pedigree$herd,tab$herd)
  testthat::expect_identical(gsim:::.gsim_block_checksum(list(tab,matrix(1:6,2))),gsim:::.gsim_checksum(list(tab,matrix(1:6,2))))
  bad<-tab;bad$sire[2]<-"c"
  testthat::expect_error(gsim_pedigree_from_table(bad),"cycle")
  bad<-tab;bad$sire[1]<-"absent"
  testthat::expect_error(gsim_pedigree_from_table(bad),"absent")
  bad<-tab;bad$sire[1]<-"c"
  testthat::expect_error(gsim_pedigree_from_table(bad),"own parent")
  bad<-tab;bad$sex[2]<-"F"
  testthat::expect_error(gsim_pedigree_from_table(bad),"sex")
  ids<-paste0("a",1:12);strata<-rep(c("H1","H2"),each=6)
  set.seed(123);before<-.Random.seed
  a<-gsim_select(ids,c(H1=2,H2=3),strata=strata,seed=52)
  b<-gsim_select(rev(ids),c(H2=3,H1=2),strata=rev(strata),seed=52)
  testthat::expect_identical(a$selected,b$selected[match(ids,b$animal)])
  testthat::expect_equal(a$inclusion_probability,c(rep(1/3,6),rep(.5,6)))
  testthat::expect_identical(.Random.seed,before)
  scored<-gsim_select(ids,2,scores=setNames(c(1,4,4,rep(0,9)),ids),method="score")
  testthat::expect_identical(scored$animal[scored$selected],ids[2:3])
  testthat::expect_identical(scored$inclusion_probability,as.double(scored$selected))
  testthat::expect_error(gsim_select(ids,13),"quota")
  testthat::expect_error(gsim_select(ids,c(unknown=2)),"stratum")
  testthat::expect_error(gsim_select(ids,2,scores=1:12,method="score"),"named")
  retained<-gsim_sample(ids,3,retain=ids[1:2],seed=52)
  testthat::expect_equal(sum(retained$selected),5)
  testthat::expect_equal(retained$inclusion_probability,c(1,1,rep(.3,10)))
  testthat::expect_identical(retained$retained,c(TRUE,TRUE,rep(FALSE,10)))
  testthat::expect_error(gsim_sample(ids,NA_real_,retain=ids),"integers")
})

testthat::test_that("incremental packed cohorts preserve transmission and marker truth", {
  root<-tempfile("breeding-");dir.create(root);on.exit(unlink(root,recursive=TRUE),add=TRUE)
  f<-.breeding_fixture(root)
  children<-paste0("c",1:8)
  mating<-gsim_mate(f$pedigree,c("s1","s2"),c("d1","d2"),children,seed=12,
    sex=setNames(rep(c("M","F"),4),children))
  reversed<-gsim_mate(f$pedigree,c("s2","s1"),c("d2","d1"),rev(children),seed=12)
  testthat::expect_identical(mating$offspring$sire,reversed$offspring$sire[match(children,reversed$offspring$animal)])
  testthat::expect_setequal(mating$parents_used,unique(c(mating$offspring$sire,mating$offspring$dam)))
  testthat::expect_identical(mating$pedigree$pedigree$herd[1:4],f$pedigree$pedigree$herd)
  set.seed(77);before<-.Random.seed
  cohort<-gsim_simulate_cohort(f$reference,mating,717,file.path(root,"cohort"))
  testthat::expect_s3_class(cohort,"gsim_cohort")
  testthat::expect_identical(cohort$bed$sample_ids,children)
  testthat::expect_false(cohort$simulation$dense_genotype_matrix_allocated)
  testthat::expect_identical(.Random.seed,before)
  full<-gsim_simulate_pedigree(f$reference,mating$pedigree,717,file.path(root,"full"),"bed")
  decode<-function(manifest,rows) .Call(gsim:::C_gsim_bed_read_selected,manifest$paths[["bed"]],
    as.integer(manifest$individual_count),as.integer(manifest$variant_count),as.integer(rows),1:3)
  expected<-decode(full,match(children,full$sample_ids))
  actual<-decode(cohort$bed,seq_along(children))
  testthat::expect_identical(actual,expected)
  fam<-.Call(gsim:::C_gsim_metadata_read_fam_external,cohort$bed$paths[["fam"]])
  testthat::expect_error(.Call(gsim:::C_gsim_metadata_read_fam,cohort$bed$paths[["fam"]]),"known parents")
  testthat::expect_identical(fam$paternal_id,mating$offspring$sire)
  testthat::expect_identical(fam$sex,rep(c(1L,2L),4))
  beta<-matrix(c(.2,-.1,.3,.1,.4,-.2),3,2,dimnames=list(paste0("m",1:3),c("milk","fat")))
  state<-gsim_trait_state(beta,center=setNames(c(.7,.5,.9),paste0("m",1:3)))
  values<-gsim_genetic_values(cohort,state)
  testthat::expect_equal(unname(values),unname(sweep(actual,2,state$center,"-")%*%beta),ignore_attr=TRUE)
  next_mating<-gsim_mate(mating$pedigree,children[c(1,3)],children[c(2,4)],c("g1","g2"),seed=7)
  second<-gsim_simulate_cohort(cohort,next_mating,717,file.path(root,"next"))
  second_full<-gsim_simulate_pedigree(f$reference,next_mating$pedigree,717,file.path(root,"full-next"),"bed")
  testthat::expect_identical(decode(second$bed,1:2),decode(second_full,match(c("g1","g2"),second_full$sample_ids)))
  testthat::expect_identical(attr(gsim_genetic_values(second,state),"trait_state_checksum"),state$checksum)
  architecture<-gsim(W=f$W,architecture="fixed",beta=setNames(c(.2,-.1,.3),colnames(f$W)),
    standardize_W=FALSE,scale_effects=FALSE,seed=3)
  founder_bed<-gsim_simulate_pedigree(f$reference,f$pedigree,717,file.path(root,"founder-bed"),"bed")
  testthat::expect_equal(as.numeric(gsim_genetic_values(founder_bed,gsim_trait_state(architecture))),as.numeric(architecture$G))
  standardized<-gsim(W=f$W,architecture="fixed",beta=setNames(c(.2,-.1,.3),colnames(f$W)),seed=3)
  testthat::expect_error(gsim_trait_state(standardized),"standardize_W=FALSE")
  testthat::expect_error(gsim_simulate_cohort(cohort,next_mating,717,file.path(root,"next")),"exists")
})

testthat::test_that("genomic records support multiple traits and appended visits", {
  G<-matrix(c(1,2,.2,.3),2,2,dimnames=list(c("a","b"),c("fat","milk")))
  R<-matrix(c(1,.2,.2,2),2,2,dimnames=list(colnames(G),colnames(G)))
  P<-matrix(c(.4,.1,.1,.8),2,2,dimnames=dimnames(R))
  rec<-data.frame(record=paste0("r",1:4),animal=rep(c("a","b"),each=2),
    trait=rep(colnames(G),2),observation_unit=rep(c("a-v1","b-v1"),each=2))
  set.seed(21);before<-.Random.seed
  a<-gsim_records(G,rec,residual_covariance=R,permanent_covariance=P,seed=8)
  testthat::expect_equal(a$truth$genetic,c(1,.2,2,.3))
  testthat::expect_equal(a$records$value,a$truth$genetic+a$truth$residual+a$truth$permanent)
  appended<-rec;appended$record<-paste0("later",1:4);appended$observation_unit<-paste0(appended$observation_unit,"-later")
  b<-gsim_records(G,rbind(rec,appended),residual_covariance=R,permanent_covariance=P,seed=8)
  testthat::expect_identical(a$records$value,b$records$value[1:4])
  testthat::expect_identical(b$truth$permanent[1:4],b$truth$permanent[5:8])
  reverse<-gsim_records(G,rec[4:1,],residual_covariance=R,permanent_covariance=P,seed=8)
  testthat::expect_identical(a$records$value,rev(reverse$records$value))
  testthat::expect_identical(.Random.seed,before)
  basis<-matrix(c(1,2,.5,.7),2,2,dimnames=list(c("a","b"),c("intercept","slope")))
  visits<-data.frame(record=c("v1","v2","v3"),animal=c("a","a","b"),trait="milk")
  design<-cbind(intercept=1,slope=c(0,1,2));rownames(design)<-visits$record
  longitudinal<-gsim_records(basis,visits,genetic_design=design[3:1,],residual_covariance=0,seed=8)
  testthat::expect_equal(longitudinal$truth$genetic,c(1,1.5,3.4))
  testthat::expect_identical(longitudinal$truth$residual,rep(0,3))
  testthat::expect_error(gsim_records(G,rec,residual_covariance=matrix(c(1,2,2,1),2)),"positive definite")
})

testthat::test_that("pool summaries and overlap covariance match a dense oracle", {
  root<-tempfile("pools-");dir.create(root);on.exit(unlink(root,recursive=TRUE),add=TRUE)
  f<-.breeding_fixture(root)
  founders<-gsim_simulate_pedigree(f$reference,f$pedigree,1,file.path(root,"bed"),"bed")
  G<-matrix(c(1,2,3,4),4,1,dimnames=list(f$ids,"milk"))
  rec<-data.frame(record=paste0("r",1:6),animal=c("s1","s1","s2","d1","d2","d1"),trait="milk")
  records<-gsim_records(G,rec,residual_covariance=2,permanent_covariance=.5,seed=6)
  tab<-data.frame(pool=c("P1","P1","P2","P2","P3","P3"),animal=c("s1","s2","s1","d1","d2","d1"),
    record=c("r1","r3","r2","r4","r5","r6"),dna_weight=c(1,3,2,2,1,1),phenotype_weight=c(1,3,1,3,1,1))
  pool<-gsim_pool(founders,tab,records,chunk_size=1)
  C<-matrix(0,3,6,dimnames=list(pool$pool_ids,rec$record));C[cbind(match(tab$pool,pool$pool_ids),match(tab$record,rec$record))]<-pool$contributions$phenotype_weight
  D<-matrix(0,3,4);D[cbind(match(tab$pool,pool$pool_ids),match(tab$animal,f$ids))]<-pool$contributions$dna_weight
  testthat::expect_equal(unname(pool$frequencies),unname(D%*%f$W/2))
  testthat::expect_equal(pool$summaries$value,as.double(C%*%records$records$value))
  Sigma<-diag(2,6)+.5*outer(rec$animal,rec$animal,"==")
  testthat::expect_equal(unname(gsim_pool_covariance(pool)),unname(C%*%Sigma%*%t(C)),tolerance=1e-12)
  testthat::expect_equal(pool$summaries$noise_variance,unname(diag(C%*%Sigma%*%t(C))))
  testthat::expect_identical(pool$summaries$matching_animal_weights,c(TRUE,FALSE,TRUE))
  chunks<-list();streamed<-gsim_pool(founders,tab,records,allele_consumer=function(x)chunks[[length(chunks)+1L]]<<-x,chunk_size=2)
  testthat::expect_null(streamed$frequencies)
  testthat::expect_equal(do.call(cbind,chunks),pool$frequencies)
  noisy<-gsim_pool(founders,tab,records,assay_depth=20,seed=44,chunk_size=1)
  testthat::expect_equal(noisy$assay_variance,pool$frequencies*(1-pool$frequencies)/20)
  changed<-gsim_pool(founders,tab[6:1,],records,markers=rev(pool$marker_ids),assay_depth=20,seed=44,chunk_size=3)
  testthat::expect_equal(noisy$frequencies,changed$frequencies[,pool$marker_ids])
  testthat::expect_true(all(noisy$frequencies>=0 & noisy$frequencies<=1))
  bad<-tab;bad$phenotype_weight[1]<--1
  testthat::expect_error(gsim_pool(founders,bad,records),"weight")
  bad<-tab;bad$animal[1]<-"d2"
  testthat::expect_error(gsim_pool(founders,bad,records),"identities")
  # Shared visits carry residual cross-trait covariance, and repeated animals
  # share permanent environment across visits, even without record overlap.
  multi<-matrix(1,4,2,dimnames=list(f$ids,c("fat","milk")))
  tr<-data.frame(record=c("f1","m1","f2","m2"),animal="s1",trait=c("fat","milk","fat","milk"),observation_unit=c("v1","v1","v2","v2"))
  R<-matrix(c(1,.2,.2,2),2);P<-matrix(c(.4,.1,.1,.8),2)
  multirecords<-gsim_records(multi,tr,residual_covariance=R,permanent_covariance=P)
  mt<-data.frame(pool=paste0("Q",1:4),animal="s1",record=tr$record)
  mp<-gsim_pool(founders,mt,multirecords)
  expected<-matrix(0,4,4)
  for(i in 1:4)for(j in 1:4)expected[i,j]<-P[(i-1)%%2+1,(j-1)%%2+1]+if(tr$observation_unit[i]==tr$observation_unit[j])R[(i-1)%%2+1,(j-1)%%2+1] else 0
  testthat::expect_equal(unname(gsim_pool_covariance(mp)),expected,tolerance=1e-12)
  testthat::expect_error(gsim_pool_covariance(mp,left=rep("Q1",257)),"duplicates|256")
})
