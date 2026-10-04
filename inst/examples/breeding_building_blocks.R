# Download-free illustration, not a validation of prediction accuracy.
# Prediction fitting belongs to an external package. Supply its named scores
# to gsim_select(method="score") to replace the phenotype-selection policy.
breeding_building_blocks <- function(output, generations=3L) {
  stopifnot(!file.exists(output),length(generations)==1L,generations>=1L)
  dir.create(output,recursive=TRUE)
  ids<-paste0("base",1:16)
  vcf<-file.path(output,"reference.vcf")
  header<-paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",ids),collapse="\t")
  body<-vapply(1:6,function(marker)paste(c("1",10*marker,paste0("m",marker),"A","G",".",".",".","GT",
    rep(c("0|0","0|1","1|0","1|1"),length.out=16)[(seq_len(16)+marker-2L)%%16L+1L]),collapse="\t"),character(1))
  writeLines(c("##fileformat=VCFv4.2",header,body),vcf,useBytes=TRUE)
  reference<-gsim::gsim_import_vcf(vcf,data.frame(chromosome="1",variant_id=paste0("m",1:6),
    genetic_position_cm=seq(0,100,length.out=6)),file.path(output,"reference"))
  sex<-rep(c("M","F"),8)
  pedigree<-gsim::gsim_pedigree_from_table(data.frame(animal=ids,sire=NA_character_,dam=NA_character_,sex=sex))
  current<-gsim::gsim_simulate_pedigree(reference,pedigree,101,file.path(output,"base-phase"),"hap")
  bed<-gsim::gsim_simulate_pedigree(reference,pedigree,101,file.path(output,"base-dosage"),"bed")
  effects<-matrix(c(.2,-.1,.15,0,.05,.1,.04,.1,-.03,.07,.02,-.08),6,2,
    dimnames=list(paste0("m",1:6),c("fat","milk")))
  state<-gsim::gsim_trait_state(effects,center=1)
  residual<-matrix(c(.4,.1,.1,1),2,2,dimnames=list(c("fat","milk"),c("fat","milk")))
  history<-list()
  for(generation in seq_len(generations)) {
    G<-gsim::gsim_genetic_values(bed,state)
    records<-gsim::gsim_records(G,residual_covariance=residual,permanent_covariance=.2,seed=200+generation)
    milk<-records$records[records$records$trait=="milk",]
    strata<-sex
    scores<-setNames(milk$value,milk$animal)
    selected<-gsim::gsim_select(ids,c(M=3,F=5),scores=scores,method="score",strata=strata)
    children<-paste0("g",generation,"-",1:32)
    mating<-gsim::gsim_mate(pedigree,selected$animal[selected$selected & selected$stratum=="M"],
      selected$animal[selected$selected & selected$stratum=="F"],children,seed=300+generation,
      sex=setNames(rep(c("M","F"),16),children))
    sampling<-gsim::gsim_sample(ids,4,retain=mating$parents_used,seed=400+generation)
    pool_rows<-data.frame(pool=rep(c("tank-a","tank-b"),length.out=length(ids)),animal=ids,
      record=milk$record[match(ids,milk$animal)],dna_weight=1,phenotype_weight=1)
    matched<-gsim::gsim_pool(bed,pool_rows,records,assay_depth=100,seed=500+generation)
    pool_rows$dna_weight<-seq_along(ids)
    mismatch<-gsim::gsim_pool(bed,pool_rows,records)
    next_cohort<-gsim::gsim_simulate_cohort(current,mating,101,file.path(output,paste0("cohort",generation)))
    history[[generation]]<-list(records=records,parent_selection=selected,parents_used=mating$parents_used,
      sampling=sampling,pool=matched,mismatched_pool=mismatch,offspring=next_cohort)
    pedigree<-mating$pedigree;current<-next_cohort;bed<-next_cohort
    ids<-children;sex<-mating$offspring$sex
  }
  list(state=state,pedigree=pedigree,cohorts=history,output=normalizePath(output,winslash="/"),
    interpretation="Small illustrative workflow; no predictions fitted or unbiasedness established.")
}
