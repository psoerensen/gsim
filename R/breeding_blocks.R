# Composable breeding simulation. R owns contracts and labels; native kernels
# own ordering, selection, mating, stochastic records and packed data products.
.gsim_block_checksum <- function(x) {
  .Call(C_gsim_breeding_checksum,serialize(x,NULL,ascii=FALSE,version=2L))
}
.gsim_block_ids <- function(x, name, unique = TRUE) {
  x <- enc2utf8(as.character(x))
  if (!length(x) || length(x) > .Machine$integer.max || anyNA(x) || any(!nzchar(x)) ||
      (unique && anyDuplicated(x))) .gsim_stop(name, " requires nonempty, nonmissing identities", if(unique)" without duplicates" else "", ".")
  x
}
.gsim_block_seed <- function(seed) {
  if (!is.numeric(seed) || length(seed)!=1L || !is.finite(seed) || seed<0 || seed>2^53-1 || seed!=floor(seed))
    .gsim_stop("seed must be an integer-valued scalar in [0, 2^53-1].")
  as.double(seed)
}
.gsim_block_integer <- function(x,name,minimum=0L) {
  if(!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x!=floor(x)) || any(x<minimum) || any(x>.Machine$integer.max))
    .gsim_stop(name," must contain supported integers.")
  as.integer(x)
}
.gsim_block_matrix <- function(x,name) {
  if(!is.matrix(x))x<-as.matrix(x)
  if(!is.numeric(x) || !all(dim(x)>0L) || any(!is.finite(x))) .gsim_stop(name," must be a nonempty finite numeric matrix.")
  storage.mode(x)<-"double";x
}

#' Construct a validated pedigree from supplied parentage
#'
#' Rebuilds topological order, mappings, diagnostics and checksums. Input rows
#' may be arbitrarily ordered; existing identifiers and labels are preserved.
#' @param table Data frame with `animal`, `sire`, `dam`; optional `sex` (`M`,
#'   `F`, `U`), `generation`, `cohort`, and logical `phenotyped`.
#'   Missing parents use `NA`, never a sentinel ID such as `"0"`.
#'   Animal ID `"0"` is reserved for FAM missing-parent encoding. One animal
#'   cannot occupy both sire and dam roles.
#' @return A `gsim_pedigree` object. Partial parentage is allowed in the table;
#'   biological meiosis still requires both parents for every offspring.
#' @export
gsim_pedigree_from_table <- function(table) {
  if(!is.data.frame(table) || !all(c("animal","sire","dam")%in%names(table)))
    .gsim_stop("table requires animal, sire and dam columns.")
  tab<-table;tab$animal<-.gsim_block_ids(tab$animal,"animal")
  if(any(tab$animal=="0")) .gsim_stop("Animal ID '0' is reserved for missing parents in FAM files.")
  for(side in c("sire","dam")) {
    tab[[side]]<-enc2utf8(as.character(tab[[side]]))
    if(any(!is.na(tab[[side]]) & !nzchar(tab[[side]]))) .gsim_stop("Parent IDs cannot be empty.")
  }
  if(any(!is.na(tab$sire) & !is.na(tab$dam) & tab$sire==tab$dam)) .gsim_stop("Sire and dam must be different animals.")
  if(length(intersect(tab$sire[!is.na(tab$sire)],tab$dam[!is.na(tab$dam)]))) .gsim_stop("An animal cannot occupy both sire and dam roles.")
  native<-.Call(C_gsim_breeding_order,tab$animal,tab$sire,tab$dam)
  if(!"sex"%in%names(tab))tab$sex<-"U"
  tab$sex<-as.character(tab$sex)
  if(anyNA(tab$sex) || any(!tab$sex%in%c("M","F","U"))) .gsim_stop("sex must be M, F or U.")
  sire<-match(tab$sire,tab$animal);dam<-match(tab$dam,tab$animal)
  if(any(tab$sex[sire[!is.na(sire)]]=="F") || any(tab$sex[dam[!is.na(dam)]]=="M")) .gsim_stop("Recorded parent sex contradicts sire/dam role.")
  if(!"generation"%in%names(tab))tab$generation<-native$generation
  tab$generation<-.gsim_block_integer(tab$generation,"generation",1L)
  for(parent in list(sire,dam)) {
    at<-which(!is.na(parent))
    if(any(tab$generation[parent[at]]>=tab$generation[at])) .gsim_stop("Parents must be in an earlier generation.")
  }
  if(!"cohort"%in%names(tab))tab$cohort<-tab$generation
  tab$cohort<-.gsim_block_integer(tab$cohort,"cohort",1L)
  if(!"phenotyped"%in%names(tab))tab$phenotyped<-TRUE
  if(!is.logical(tab$phenotyped) || anyNA(tab$phenotyped)) .gsim_stop("phenotyped must be logical without missing values.")
  canonical<-tab$animal[native$order];founder<-is.na(tab$sire)&is.na(tab$dam)
  pairs<-paste(sire,dam,sep="|");pairs<-pairs[!is.na(sire)&!is.na(dam)]
  diagnostics<-list(animals=nrow(tab),founders=sum(founder),generations=length(unique(tab$generation)),
    generation_sizes=table(tab$generation),sires_used=length(unique(tab$sire[!is.na(sire)])),
    dams_used=length(unique(tab$dam[!is.na(dam)])),unknown_sires=sum(!founder & is.na(sire)),
    unknown_dams=sum(!founder & is.na(dam)),unphenotyped=sum(!tab$phenotyped),
    paternal_half_sib_sires=sum(table(tab$sire)>=2L),maternal_half_sib_dams=sum(table(tab$dam)>=2L),
    full_sib_families=sum(table(pairs)>=2L),topological=TRUE,storage_bytes=as.numeric(object.size(tab)))
  ordered<-tab[native$order,c("animal","sire","dam"),drop=FALSE];rownames(ordered)<-NULL
  structure(list(pedigree=tab,canonical_order=canonical,external_order=tab$animal,
    mapping=data.frame(animal=canonical,canonical_index=seq_along(canonical),external_index=native$order),
    settings=list(source="validated supplied parentage",generation_sizes=table(tab$generation)),
    diagnostics=diagnostics,checksums=list(identifiers_and_parents=.gsim_block_checksum(ordered),
      canonical_order=.gsim_block_checksum(canonical),external_order=.gsim_block_checksum(tab$animal))),class="gsim_pedigree")
}

#' Select animals by random sampling or supplied scores
#'
#' Selection never fits a prediction model. Equal scores are resolved by UTF-8
#' animal identity. Random ranks are keyed by identity and seed, independently
#' of input order; sampling without replacement has quota/eligible probabilities.
#' @param animals Unique candidate IDs.
#' @param n One quota or a named quota vector exactly covering `strata`.
#' @param scores Named finite scores for all eligible animals when `method`
#'   is `"score"`; predictions are supplied by an external caller.
#' @param method `"random"` or `"score"`.
#' @param strata Optional stratum labels aligned to `animals`.
#' @param eligible Logical eligibility mask aligned to `animals`.
#' @param seed Nonnegative identity-stream seed.
#' @param oracle Logical provenance flag explicitly marking simulation-truth selection.
#' @return Selection table for eligible animals, including scores, selected
#'   mask and inclusion probabilities. Score selection has conditional 0/1
#'   probabilities and does not establish positive population sampling support.
#' @export
gsim_select <- function(animals,n,scores=NULL,method=c("random","score"),strata=NULL,
  eligible=rep(TRUE,length(animals)),seed=1,oracle=FALSE) {
  animals<-.gsim_block_ids(animals,"animals");method<-match.arg(method);seed<-.gsim_block_seed(seed)
  if(!is.logical(eligible) || length(eligible)!=length(animals) || anyNA(eligible)) .gsim_stop("eligible must align to animals.")
  oracle<-.gsim_meiosis_flag(oracle,"oracle")
  if(is.null(strata))strata<-rep("all",length(animals))
  strata<-.gsim_block_ids(strata,"strata",unique=FALSE)
  if(length(strata)!=length(animals)) .gsim_stop("strata must align to animals.")
  ids<-animals[eligible];h<-strata[eligible];if(!length(ids)) .gsim_stop("At least one animal must be eligible.")
  labels<-sort(unique(h),method="radix")
  if(length(labels)==1L && length(n)==1L && is.null(names(n)))n<-setNames(n,labels)
  if(is.null(names(n)) || anyDuplicated(names(n)) || !setequal(names(n),labels)) .gsim_stop("n must name each eligible stratum exactly once.")
  quotas<-.gsim_block_integer(n[match(labels,names(n))],"n")
  values<-rep(0,length(ids))
  if(method=="score") {
    if(!is.numeric(scores) || is.null(names(scores)) || anyDuplicated(names(scores)) ||
       anyNA(match(ids,names(scores)))) .gsim_stop("scores must be uniquely named for eligible animals.")
    values<-as.double(scores[match(ids,names(scores))]);if(any(!is.finite(values))) .gsim_stop("Eligible scores must be finite.")
  }
  group<-match(h,labels);chosen<-.Call(C_gsim_breeding_select,ids,as.integer(group),quotas,values,as.integer(method=="random"),seed)
  selected<-seq_along(ids)%in%chosen
  probability<-if(method=="random")quotas[group]/tabulate(group,nbins=length(labels))[group] else as.numeric(selected)
  result<-data.frame(animal=ids,stratum=h,score=if(method=="score")values else NA_real_,
    selected=selected,inclusion_probability=probability,stringsAsFactors=FALSE)
  attr(result,"settings")<-list(method=method,seed=seed,oracle=oracle,rng="identity-keyed SplitMix64; R RNG unused")
  class(result)<-c("gsim_selection","data.frame");result
}

#' Mate selected sires and dams to construct a new cohort
#' @param pedigree Validated `gsim_pedigree` containing all candidate parents.
#' @param sires,dams Selected parent ID vectors. Known sex must match its role.
#' @param offspring New unique offspring IDs, disjoint from the pedigree.
#' @param seed Nonnegative identity-stream seed.
#' @param generation Optional offspring generation; defaults to one after the latest parent.
#' @param sex Optional named offspring sex vector (`M` or `F`); otherwise drawn natively.
#' @return List with `offspring`, freshly validated combined `pedigree`, actual
#'   `parents_used` and settings. Each child draws both parents independently,
#'   uniformly with replacement, using a stream keyed by child identity.
#' @export
gsim_mate <- function(pedigree,sires,dams,offspring,seed=1,generation=NULL,sex=NULL) {
  if(!inherits(pedigree,"gsim_pedigree")) .gsim_stop("pedigree must be a gsim_pedigree object.")
  pedigree<-gsim_pedigree_from_table(pedigree$pedigree);tab<-pedigree$pedigree
  sires<-.gsim_block_ids(sires,"sires");dams<-.gsim_block_ids(dams,"dams");offspring<-.gsim_block_ids(offspring,"offspring")
  if(anyNA(match(c(sires,dams),tab$animal)) || length(intersect(sires,dams)) || length(intersect(offspring,tab$animal)))
    .gsim_stop("Parents must be present, parent pools disjoint, and offspring IDs new.")
  if(any(tab$sex[match(sires,tab$animal)]=="F") || any(tab$sex[match(dams,tab$animal)]=="M")) .gsim_stop("Parent sex contradicts mating role.")
  seed<-.gsim_block_seed(seed);draw<-.Call(C_gsim_breeding_mate,offspring,sires,dams,seed)
  if(is.null(generation))generation<-max(tab$generation[match(c(sires,dams),tab$animal)])+1L
  generation<-.gsim_block_integer(generation,"generation",1L);if(length(generation)!=1L) .gsim_stop("generation must be scalar.")
  if(!is.null(sex)) {
    if(is.null(names(sex)) || anyDuplicated(names(sex)) || !setequal(names(sex),offspring)) .gsim_stop("sex must name all offspring exactly once.")
    draw$sex<-as.character(sex[match(offspring,names(sex))]);if(anyNA(draw$sex) || any(!draw$sex%in%c("M","F"))) .gsim_stop("Offspring sex must be M or F.")
  }
  children<-data.frame(animal=offspring,sire=draw$sire,dam=draw$dam,sex=draw$sex,
    generation=generation,cohort=generation,phenotyped=TRUE,stringsAsFactors=FALSE)
  # Preserve caller-owned labels on existing animals; new labels are unknown.
  for (label in setdiff(names(tab),names(children))) children[[label]]<-tab[[label]][rep(NA_integer_,nrow(children))]
  children<-children[,names(tab),drop=FALSE]
  list(offspring=children,pedigree=gsim_pedigree_from_table(rbind(tab,children)),
    parents_used=unique(c(children$sire,children$dam)),settings=list(seed=seed,mating="uniform independent parental draws"))
}

#' Create persistent marker-effect and reporting-base state
#' @param effects Finite marker-by-trait matrix with marker row names, or a
#'   `gsim` result created with `standardize_W=FALSE`, whose `B` matrix supplies
#'   raw-dosage effects. Standardized results are rejected because their dosage
#'   scaling is not retained. A supplied matrix must already use raw-dosage units.
#' @param center Scalar or marker-named dosage center. It remains fixed across cohorts.
#' @return Compact `gsim_trait_state`. Effects are never rescaled or redrawn
#'   when the state is applied to later cohorts. Missing BED calls are imputed
#'   with this fixed center, not with selected-cohort allele frequencies.
#' @export
#' @importFrom stats setNames
gsim_trait_state <- function(effects,center=0) {
  if(inherits(effects,"gsim")) {
    if(!identical(effects$settings$standardize_W,FALSE))
      .gsim_stop("A gsim result must use standardize_W=FALSE; standardized dosage scaling is not retained.")
    effects<-effects$B
  }
  effects<-.gsim_block_matrix(effects,"effects");markers<-.gsim_block_ids(rownames(effects),"effect marker names")
  if(length(markers)!=nrow(effects)) .gsim_stop("Effects require marker row names.")
  if(is.null(colnames(effects)))colnames(effects)<-paste0("D",seq_len(ncol(effects)))
  .gsim_block_ids(colnames(effects),"trait names")
  if(length(center)==1L)center<-setNames(rep(center,length(markers)),markers)
  if(!is.numeric(center) || is.null(names(center)) || anyDuplicated(names(center)) || !setequal(names(center),markers))
    .gsim_stop("center must be scalar or name all effect markers exactly once.")
  center<-as.double(center[match(markers,names(center))]);if(any(!is.finite(center)) || any(center<0|center>2)) .gsim_stop("Dosage centers must be in [0,2].")
  structure(list(effects=effects,center=setNames(center,markers),
    checksum=.gsim_block_checksum(list(effects,center))),class="gsim_trait_state")
}

.gsim_block_bed <- function(genotypes) {
  if(inherits(genotypes,"gsim_cohort"))genotypes<-genotypes$bed$paths[["bed"]]
  if(is.list(genotypes) && !is.null(genotypes$paths))genotypes<-genotypes$paths[["bed"]]
  if(!is.character(genotypes) || length(genotypes)!=1L || is.na(genotypes)) .gsim_stop("genotypes must identify a BED dataset.")
  prefix<-sub("[.]bed$","",genotypes);prefix<-normalizePath(prefix,winslash="/",mustWork=FALSE)
  paths<-setNames(paste0(prefix,c(".bed",".bim",".fam")),c("bed","bim","fam"))
  if(!all(file.exists(paths))) .gsim_stop("Complete BED/BIM/FAM input is required.")
  fam<-.Call(C_gsim_metadata_read_fam_external,enc2utf8(paths[["fam"]]));bim<-.Call(C_gsim_metadata_read_bim,enc2utf8(paths[["bim"]]))
  .gsim_block_ids(fam$individual_id,"FAM IDs");.gsim_block_ids(bim$variant_id,"BIM IDs")
  list(paths=paths,samples=fam$individual_id,markers=bim$variant_id)
}

#' Apply persistent genomic truth to a BED cohort
#' @param genotypes BED prefix, BED manifest, or `gsim_cohort`.
#' @param state A persistent `gsim_trait_state`.
#' @param ids Optional unique sample IDs, in desired result order.
#' @return Animal-by-trait genomic values accumulated natively from packed BED.
#'   No phenotype generation, realized-variance normalization, cohort recentering,
#'   qgg dependency, or full decoded genotype matrix is involved.
#' @export
gsim_genetic_values <- function(genotypes,state,ids=NULL) {
  if(!inherits(state,"gsim_trait_state")) .gsim_stop("state must be a gsim_trait_state.")
  state<-gsim_trait_state(state$effects,state$center);bed<-.gsim_block_bed(genotypes)
  if(is.null(ids))ids<-bed$samples else ids<-.gsim_block_ids(ids,"ids")
  rows<-match(ids,bed$samples);columns<-match(rownames(state$effects),bed$markers)
  if(anyNA(rows) || anyNA(columns)) .gsim_stop("Selected animal/effect marker identities must be present in BED.")
  descriptor<-list(enc2utf8(bed$paths[["bed"]]),as.integer(length(bed$samples)),as.integer(length(bed$markers)),as.integer(rows))
  result<-.Call(C_gsim_bed_accumulate,list(descriptor),rep.int(1L,length(columns)),as.integer(columns),
    as.double(state$center),as.double(state$center),rep(1,length(columns)),state$effects)
  dimnames(result)<-list(ids,colnames(state$effects));attr(result,"trait_state_checksum")<-state$checksum;result
}

.gsim_block_factor <- function(covariance,traits,name) {
  t<-length(traits)
  if(is.null(covariance))return(matrix(0,t,t,dimnames=list(traits,traits)))
  if(length(covariance)==1L)covariance<-diag(as.double(covariance),t)
  covariance<-.gsim_block_matrix(covariance,name)
  if(!identical(dim(covariance),c(t,t))) .gsim_stop(name," dimensions must match observed traits.")
  if(!is.null(rownames(covariance)) || !is.null(colnames(covariance))) {
    if(is.null(rownames(covariance)) || is.null(colnames(covariance)) || anyDuplicated(rownames(covariance)) ||
      anyDuplicated(colnames(covariance)) || !setequal(rownames(covariance),traits) || !setequal(colnames(covariance),traits))
      .gsim_stop(name," names must match observed traits.")
    covariance<-covariance[traits,traits,drop=FALSE]
  }
  if(max(abs(covariance-t(covariance)))>1e-12) .gsim_stop(name," must be symmetric.")
  if(all(covariance==0))return(matrix(0,t,t,dimnames=list(traits,traits)))
  factor<-tryCatch(t(chol(covariance)),error=function(e)NULL)
  if(is.null(factor)) .gsim_stop(name," must be positive definite or entirely zero.")
  dimnames(factor)<-list(traits,traits);factor
}

#' Generate single, multiple or longitudinal records from genomic values
#' @param genetic_values Finite animal-by-trait or animal-by-basis-coefficient
#'   matrix with animal row names and column names.
#' @param records Data frame with unique `record`, `animal`, and `trait` IDs.
#'   Optional `observation_unit` groups correlated residuals across traits.
#'   Default creates one record per animal and trait.
#' @param genetic_design Optional record-by-genetic-column design, for example
#'   longitudinal basis rows. Default uses the record's trait column.
#' @param fixed Scalar or one fixed-effect contribution per record.
#' @param residual_covariance,permanent_covariance Trait covariance matrices
#'   (or scalar diagonal values). Entirely zero matrices are supported.
#'   Residuals are shared within observation units; permanent effects within animals.
#' @param seed Nonnegative seed for separate identity-keyed residual and permanent streams.
#' @return List with observable `records`, separate simulation `truth`, and a
#'   covariance `noise` descriptor for subsequent pooled observations. The
#'   supplied genomic values, rather than parent-average workload latents,
#'   generate the genetic record contribution. At most 64 observed traits.
#' @export
gsim_records <- function(genetic_values,records=NULL,genetic_design=NULL,fixed=0,
  residual_covariance=1,permanent_covariance=NULL,seed=1) {
  G<-.gsim_block_matrix(genetic_values,"genetic_values")
  animals<-.gsim_block_ids(rownames(G),"genomic animal names");coefficients<-.gsim_block_ids(colnames(G),"genomic column names")
  if(length(animals)!=nrow(G) || length(coefficients)!=ncol(G)) .gsim_stop("Genomic values require row/column identities.")
  if(is.null(records)) {
    records<-expand.grid(animal=animals,trait=coefficients,stringsAsFactors=FALSE)
    records$record<-paste0(nchar(records$animal),":",records$animal,":",records$trait)
    records$observation_unit<-records$animal
  }
  if(!is.data.frame(records) || !all(c("record","animal","trait")%in%names(records))) .gsim_stop("records requires record, animal, trait.")
  records$record<-.gsim_block_ids(records$record,"record")
  records$animal<-.gsim_block_ids(records$animal,"record animal",FALSE)
  records$trait<-.gsim_block_ids(records$trait,"record trait",FALSE)
  if(!"observation_unit"%in%names(records))records$observation_unit<-records$record
  records$observation_unit<-.gsim_block_ids(records$observation_unit,"observation unit",FALSE)
  if(anyDuplicated(records[c("observation_unit","trait")])) .gsim_stop("An observation unit cannot repeat the same trait.")
  if(anyNA(match(records$animal,animals))) .gsim_stop("Record animal is absent from genomic values.")
  # Observation units represent one animal's simultaneous trait measurements.
  if(any(vapply(split(records$animal,records$observation_unit),function(x)length(unique(x))!=1L,logical(1))))
    .gsim_stop("An observation unit cannot span different animals.")
  traits<-sort(unique(records$trait),method="radix");if(length(traits)>64L) .gsim_stop("At most 64 observed traits are supported.")
  if(is.null(genetic_design)) {
    column<-match(records$trait,coefficients);if(anyNA(column)) .gsim_stop("Trait absent from genomic columns; supply genetic_design for basis coefficients.")
    genetic_design<-matrix(0,nrow(records),ncol(G));genetic_design[cbind(seq_len(nrow(records)),column)]<-1
  }
  design<-.gsim_block_matrix(genetic_design,"genetic_design")
  if(!identical(dim(design),c(nrow(records),ncol(G)))) .gsim_stop("genetic_design must align to records and genomic columns.")
  if(!is.null(rownames(design))) {
    if(anyDuplicated(rownames(design)) || !setequal(rownames(design),records$record)) .gsim_stop("genetic_design row names must match record IDs.")
    design<-design[records$record,,drop=FALSE]
  }
  if(!is.null(colnames(design))) {
    if(anyDuplicated(colnames(design)) || !setequal(colnames(design),coefficients)) .gsim_stop("genetic_design names must match genomic columns.")
    design<-design[,coefficients,drop=FALSE]
  }
  if(!is.numeric(fixed)) .gsim_stop("fixed must be numeric.")
  if(length(fixed)>1L && !is.null(names(fixed))) {
    if(anyDuplicated(names(fixed)) || !setequal(names(fixed),records$record)) .gsim_stop("fixed names must match record IDs.")
    fixed<-fixed[records$record]
  }
  fixed<-as.double(fixed);if(length(fixed)==1L)fixed<-rep(fixed,nrow(records))
  if(length(fixed)!=nrow(records) || any(!is.finite(fixed))) .gsim_stop("fixed must align to records.")
  R<-.gsim_block_factor(residual_covariance,traits,"residual_covariance")
  P<-.gsim_block_factor(permanent_covariance,traits,"permanent_covariance");seed<-.gsim_block_seed(seed)
  value<-.Call(C_gsim_breeding_records,G,as.integer(match(records$animal,animals)),design,fixed,
    records$animal,records$observation_unit,as.integer(match(records$trait,traits)),traits,R,P,seed)
  records$value<-value[,4]
  structure(list(records=records,truth=list(genetic=value[,1],fixed=fixed,permanent=value[,2],residual=value[,3]),
    noise=list(records=records[,c("record","animal","trait","observation_unit")],traits=traits,
      residual_factor=R,permanent_factor=P),settings=list(seed=seed,
      rng="identity-keyed SplitMix64/Box-Muller; R RNG unused",genetic_source="supplied genomic values")),class="gsim_records")
}

#' Sample individual records and retain explicitly selected animals
#' @param animals Eligible animal IDs.
#' @param n Quotas as in [gsim_select()].
#' @param strata Optional sampling strata.
#' @param retain Animal IDs included with probability one, for example selected parents.
#'   Quotas apply to the remaining animals, in addition to retention.
#' @param seed Nonnegative identity-stream seed.
#' @return Eligible-animal table with sampled/retained masks and conditional
#'   inclusion probabilities. Retained animals need not be representative.
#' @export
gsim_sample <- function(animals,n,strata=NULL,retain=character(),seed=1) {
  .gsim_block_integer(n,"n")
  animals<-.gsim_block_ids(animals,"animals");if(length(retain))retain<-.gsim_block_ids(retain,"retain")
  if(anyNA(match(retain,animals))) .gsim_stop("Retained animals must be eligible.")
  kept<-animals%in%retain
  if(is.null(strata))strata<-rep("all",length(animals))
  strata<-.gsim_block_ids(strata,"strata",FALSE);if(length(strata)!=length(animals)) .gsim_stop("strata must align to animals.")
  if(all(kept)) {
    if(any(n!=0)) .gsim_stop("No animals remain available for sampling.")
    result<-data.frame(animal=animals,stratum=strata,selected=TRUE,retained=TRUE,inclusion_probability=1)
  } else {
    sample<-gsim_select(animals,n,method="random",strata=strata,eligible=!kept,seed=seed)
    at<-match(animals,sample$animal);result<-data.frame(animal=animals,stratum=strata,
      selected=kept,retained=kept,inclusion_probability=as.numeric(kept))
    result$selected[!kept]<-sample$selected[at[!kept]]
    result$inclusion_probability[!kept]<-sample$inclusion_probability[at[!kept]]
  }
  attr(result,"settings")<-list(seed=.gsim_block_seed(seed),quota="additional random sample among non-retained animals")
  class(result)<-c("gsim_sampling","data.frame");result
}
