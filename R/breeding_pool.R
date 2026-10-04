#' Form paired pooled allele frequencies and phenotype summaries
#'
#' Reads packed BED one marker at a time and returns or streams pool-by-marker
#' frequencies. DNA and phenotype weights are normalized separately. This is
#' a simulation observation model, not a fitted bias correction or predictor.
#' @param genotypes BED prefix, BED manifest or `gsim_cohort`.
#' @param contributions Data frame with `pool`, `animal`, `record`, and optional
#'   nonnegative `dna_weight`, `phenotype_weight` (default one). Each pool/record
#'   pair must be unique. Repeated animals across records are allowed; DNA
#'   weights then add across those contributions.
#' @param records A [gsim_records()] object, supplying observable records and
#'   its noise descriptor. Genetic truth is never used to form pooled summaries.
#' @param markers Optional unique BIM marker IDs, in desired order.
#' @param allele_consumer Optional function called with each frequency matrix
#'   chunk. With a consumer, the full frequency matrix is not retained.
#' @param chunk_size Marker count per chunk, between 1 and 64.
#' @param assay_depth Zero for exact allele frequencies; positive integer for
#'   independent binomial read-count measurement per pool and marker. This
#'   model does not include extraction bias or read correlation.
#' @param seed Nonnegative integer-valued assay seed.
#' @return `gsim_pools` containing frequencies (or NULL), weighted phenotypes,
#'   measurement-noise variances, normalized contributions and a covariance
#'   query descriptor. Matching weights are diagnosed over animals within each
#'   pool. Matching weights alone cannot establish a common genetic design
#'   across different traits, visits or fixed effects. Assay error has variance
#'   p(1-p)/depth conditional on the true frequency; it is separate from the
#'   phenotype covariance. Callback side effects are caller-owned.
#'   `assay_variance` reports conditional frequency measurement variance using
#'   simulated true frequencies; streamed chunks carry it as an attribute.
#' @export
gsim_pool <- function(genotypes, contributions, records, markers=NULL,
                      allele_consumer=NULL, chunk_size=64L, assay_depth=0L, seed=1) {
  bed<-.gsim_block_bed(genotypes)
  if(!inherits(records,"gsim_records")) .gsim_stop("records must be a gsim_records object with a noise descriptor.")
  if(!is.data.frame(contributions) || !all(c("pool","animal","record")%in%names(contributions)))
    .gsim_stop("contributions requires pool, animal and record columns.")
  tab<-contributions
  for(field in c("pool","animal","record")) tab[[field]]<-.gsim_block_ids(tab[[field]],field,FALSE)
  if(anyDuplicated(tab[c("pool","record")])) .gsim_stop("A pool cannot repeat the same record.")
  pool_ids<-sort(unique(tab$pool),method="radix");group<-as.integer(match(tab$pool,pool_ids));count<-length(pool_ids)
  observation<-records$records;noise<-records$noise
  .gsim_block_ids(observation$record,"observable record IDs")
  .gsim_block_ids(noise$records$record,"noise record IDs")
  at<-match(tab$record,observation$record);noise_at<-match(tab$record,noise$records$record)
  rows<-match(tab$animal,bed$samples)
  if(anyNA(at) || anyNA(noise_at) || anyNA(rows) || any(tab$animal!=observation$animal[at]) ||
     any(tab$animal!=noise$records$animal[noise_at])) .gsim_stop("Pool record/animal identities must match observed records, noise and BED.")
  for(field in c("dna_weight","phenotype_weight")) {
    if(!field%in%names(tab))tab[[field]]<-1
    if(!is.numeric(tab[[field]])) .gsim_stop(field," must be numeric.")
    tab[[field]]<-.Call(C_gsim_breeding_normalize,group,as.double(tab[[field]]),as.integer(count))
  }
  # Summaries exclude any truth fields and use only declared observations.
  value<-as.double(observation$value[at]);if(any(!is.finite(value))) .gsim_stop("Observable pool phenotypes must be finite.")
  pooled<-.Call(C_gsim_breeding_pool_values,group,tab$phenotype_weight,value,as.integer(count))
  descriptor<-list(pool=group,weight=tab$phenotype_weight,
    animal=tab$animal,unit=as.character(noise$records$observation_unit[noise_at]),
    trait=as.integer(match(noise$records$trait[noise_at],noise$traits)),
    residual_factor=noise$residual_factor,permanent_factor=noise$permanent_factor,groups=as.integer(count))
  if(anyNA(descriptor$trait)) .gsim_stop("Noise trait identities do not align.")
  variance<-do.call(.Call,c(list(C_gsim_breeding_pool_covariance),descriptor,
    list(as.integer(seq_len(count)),integer(),1L)))
  compatible<-vapply(pool_ids,function(id){part<-tab[tab$pool==id,,drop=FALSE]
    summed<-rowsum(cbind(part$dna_weight,part$phenotype_weight),part$animal,reorder=FALSE)
    all(abs(summed[,1]-summed[,2])<=1e-12)},logical(1))
  chunk_size<-.gsim_block_integer(chunk_size,"chunk_size",1L)
  if(length(chunk_size)!=1L || chunk_size>64L) .gsim_stop("chunk_size must be between 1 and 64.")
  assay_depth<-.gsim_block_integer(assay_depth,"assay_depth")
  if(length(assay_depth)!=1L) .gsim_stop("assay_depth must be scalar.")
  seed<-.gsim_block_seed(seed)
  if(is.null(markers))markers<-bed$markers else markers<-.gsim_block_ids(markers,"markers")
  columns<-match(markers,bed$markers);if(anyNA(columns)) .gsim_stop("Pooled markers are absent from BIM.")
  if(!is.null(allele_consumer) && !is.function(allele_consumer)) .gsim_stop("allele_consumer must be a function.")
  frequencies<-if(is.null(allele_consumer))matrix(NA_real_,count,length(markers),dimnames=list(pool_ids,markers)) else NULL
  assay_variance<-if(is.null(allele_consumer))matrix(0,count,length(markers),dimnames=list(pool_ids,markers)) else NULL
  for(start in seq.int(1L,length(markers),by=chunk_size)) {
    chunk<-seq.int(start,min(length(markers),start+chunk_size-1L))
    result<-.Call(C_gsim_breeding_pool_bed,enc2utf8(bed$paths[["bed"]]),as.integer(length(bed$samples)),
      as.integer(length(bed$markers)),as.integer(rows),group,tab$dna_weight,pool_ids,
      as.integer(columns[chunk]),markers[chunk],assay_depth,seed)
    dimnames(result)<-list(pool_ids,markers[chunk])
    noise_variance<-attr(result,"assay_variance");dimnames(noise_variance)<-dimnames(result)
    attr(result,"assay_variance")<-noise_variance
    if(is.null(allele_consumer)) {
      frequencies[,chunk]<-result;assay_variance[,chunk]<-noise_variance
    } else allele_consumer(result)
  }
  structure(list(frequencies=frequencies,assay_variance=assay_variance,summaries=data.frame(pool=pool_ids,value=pooled,
    noise_variance=variance,matching_animal_weights=compatible),contributions=tab,
    noise=descriptor,pool_ids=pool_ids,marker_ids=markers,settings=list(seed=seed,
      assay_depth=assay_depth,chunk_size=chunk_size,streamed=!is.null(allele_consumer),
      covariance="residual and permanent environment; excludes genetic covariance",
      allele="BIM bit1 allele frequency; BED A1 dosage divided by two",
      rng="identity-keyed assay streams; binomial reproducibility scoped to the compiled implementation")),class="gsim_pools")
}

#' Query pooled phenotype measurement-noise covariance
#' @param object A `gsim_pools` object.
#' @param left,right Unique pool IDs; maximum 256 on each side. Defaults query
#'   all pools, and require an explicit bounded subset when there are over 256.
#' @return Covariance matrix including shared observation-unit residuals and
#'   shared-animal permanent effects. Genetic covariance and genotype-assay
#'   error are excluded. No complete pool-by-pool matrix is stored by gsim_pool.
#' @export
gsim_pool_covariance <- function(object,left=object$pool_ids,right=left) {
  if(!inherits(object,"gsim_pools")) .gsim_stop("object must be gsim_pools.")
  left<-.gsim_block_ids(left,"left");right<-.gsim_block_ids(right,"right")
  if(length(left)>256L || length(right)>256L) .gsim_stop("Query at most 256 pools per side.")
  a<-match(left,object$pool_ids);b<-match(right,object$pool_ids)
  if(anyNA(a) || anyNA(b)) .gsim_stop("Covariance pool identities are absent.")
  value<-do.call(.Call,c(list(C_gsim_breeding_pool_covariance),object$noise,
    list(as.integer(a),as.integer(b),0L)))
  dimnames(value)<-list(left,right);value
}
