#' Generate one cohort from already simulated phased parents
#'
#' Uses existing packed chromosome meiosis and writes only new offspring.
#' Parent genotypes are reused, with no ancestor regeneration or dense genotype
#' matrix. HAP and BED outputs have separate metadata prefixes.
#' @param parents Phased dataset prefix, phased manifest, or `gsim_cohort`.
#'   All used parents must be present in this one dataset.
#' @param offspring Data frame containing `animal`, `sire`, `dam`, optionally
#'   sex and cohort labels; alternatively the result of [gsim_mate()].
#' @param seed Nonnegative integer-valued meiosis seed.
#' @param output New directory for the cohort's `phase` and `dosage` triplets.
#'   Existing directories are refused, including after a partial publication.
#' @return `gsim_cohort` with HAP/BED manifests, offspring metadata, actual
#'   parent IDs and simulation settings. Meiosis is single threaded. Memory
#'   holds one packed source chromosome, used parents, and the new cohort;
#'   previous cohorts are not accumulated. Parent datasets can be consolidated
#'   externally when mating across several historical cohorts.
#' @export
gsim_simulate_cohort <- function(parents, offspring, seed, output) {
  if (inherits(parents,"gsim_cohort")) parents<-parents$hap
  if (is.list(parents) && !inherits(parents,"gsim_reference") && !is.null(parents$paths[["hap"]]))
    parents<-sub("[.]hap$","",parents$paths[["hap"]])
  parents<-.gsim_public_phased_input(parents,"parents")
  if (is.list(offspring) && !is.data.frame(offspring) && !is.null(offspring$offspring)) offspring<-offspring$offspring
  if (!is.data.frame(offspring) || !all(c("animal","sire","dam")%in%names(offspring)))
    .gsim_stop("offspring requires animal, sire and dam columns.")
  offspring$animal<-.gsim_block_ids(offspring$animal,"offspring IDs")
  offspring$sire<-.gsim_block_ids(offspring$sire,"offspring sires",FALSE)
  offspring$dam<-.gsim_block_ids(offspring$dam,"offspring dams",FALSE)
  if (any(offspring$sire==offspring$dam)) .gsim_stop("Sire and dam must differ.")
  if (length(intersect(offspring$sire,offspring$dam))) .gsim_stop("Parent sire and dam roles must be disjoint.")
  seed<-.gsim_block_seed(seed)
  reader<-.gsim_hap_dataset_open(parents$prefix)
  on.exit(try(.gsim_hap_dataset_close(reader),silent=TRUE),add=TRUE)
  used<-sort(unique(c(offspring$sire,offspring$dam)),method="radix")
  source_rows<-match(used,reader$samples$individual_id)
  if (anyNA(source_rows) || length(intersect(offspring$animal,reader$samples$individual_id)))
    .gsim_stop("Used parents must be present and offspring IDs must be new.")
  source_sex<-reader$samples$sex[source_rows]
  if (any(source_sex[used%in%offspring$sire]==2L) || any(source_sex[used%in%offspring$dam]==1L))
    .gsim_stop("Recorded parent sex contradicts mating role.")
  local_parents<-data.frame(animal=used,sire=NA_character_,dam=NA_character_,
    sex=ifelse(used%in%offspring$sire,"M","F"),generation=1L)
  local_children<-offspring[,c("animal","sire","dam"),drop=FALSE]
  local_children$sex<-if("sex"%in%names(offspring))offspring$sex else "U"
  local_children$generation<-2L
  local<-gsim_pedigree_from_table(rbind(local_parents,local_children))
  # Publication metadata retains global generation/cohort labels, while the
  # two-generation local pedigree is only an execution view of frozen parents.
  if (!"sex"%in%names(offspring))offspring$sex<-"U"
  local_sex<-local$pedigree$sex[match(local$canonical_order,local$pedigree$animal)]
  child_metadata<-.gsim_plink_pedigree_metadata(local,sex=as.integer(match(local_sex,c("U","M","F"))-1L))
  child_metadata<-child_metadata[match(offspring$animal,child_metadata$individual_id),,drop=FALSE]
  if (!is.character(output) || length(output)!=1L || is.na(output) || !nzchar(output)) .gsim_stop("output must be one directory path.")
  output<-normalizePath(output,winslash="/",mustWork=FALSE)
  if (file.exists(output)) .gsim_stop("Cohort output already exists.")
  if (!dir.create(output,recursive=TRUE)) .gsim_stop("Cannot create cohort output directory.")
  provenance<-list(operation="incremental packed cohort meiosis",parent_prefix=parents$prefix,seed=seed)
  hap<-bed<-NULL;complete<-FALSE
  on.exit({if(!complete) {
    if(!is.null(hap))try(.gsim_hap_dataset_cancel(hap),silent=TRUE)
    if(!is.null(bed))try(.gsim_plink_dataset_cancel(bed),silent=TRUE)
  }},add=TRUE)
  hap<-.gsim_hap_dataset_create(file.path(output,"phase"),child_metadata,provenance=provenance,allow_external_parents=TRUE)
  bed<-.gsim_plink_dataset_create(file.path(output,"dosage"),child_metadata,provenance=provenance,allow_external_parents=TRUE)
  process<-function(chromosome) {
    handles<-list()
    close_all<-function()for(handle in handles)try(.gsim_packed_close(handle),silent=TRUE)
    on.exit(close_all(),add=TRUE)
    keep<-function(handle){handles[[length(handles)+1L]]<<-handle;handle}
    base<-.gsim_hap_dataset_load_chromosome(reader,chromosome)
    keep(base$h1);keep(base$h2)
    index<-match(chromosome,reader$chromosome)
    start<-sum(head(reader$marker_count,index-1L))+1L
    rows<-seq.int(start,length.out=reader$marker_count[[index]])
    variants<-reader$variants[rows,,drop=FALSE];markers<-variants$variant_id
    selected<-list(h1=keep(.gsim_packed_zero(length(used),length(markers),used,markers)),
      h2=keep(.gsim_packed_zero(length(used),length(markers),used,markers)))
    for(side in c("h1","h2"))for(i in seq_along(used))
      .gsim_packed_copy_interval(selected[[side]],i,base[[side]],source_rows[[i]],1L,length(markers))
    generated<-.gsim_pedigree_genotypes_packed_chromosome(local,selected,
      rep.int(chromosome,length(markers)),variants$genetic_position_cm/100,seed,
      return_haplotypes=TRUE,return_genotypes=FALSE,return_crossovers=FALSE)
    keep(generated$h1);keep(generated$h2)
    child<-list(h1=keep(.gsim_packed_zero(nrow(offspring),length(markers),offspring$animal,markers)),
      h2=keep(.gsim_packed_zero(nrow(offspring),length(markers),offspring$animal,markers)))
    child_rows<-match(offspring$animal,local$canonical_order)
    for(side in c("h1","h2"))for(i in seq_len(nrow(offspring)))
      .gsim_packed_copy_interval(child[[side]],i,generated[[side]],child_rows[[i]],1L,length(markers))
    .gsim_hap_dataset_append(hap,chromosome,child$h1,child$h2,variants)
    .gsim_plink_dataset_append(bed,chromosome,child$h1,child$h2,variants)
    sum(vapply(handles,function(handle)as.double(.gsim_packed_info(handle)[[4L]]),numeric(1)))
  }
  payload<-vapply(reader$chromosome,process,numeric(1))
  hap_manifest<-.gsim_hap_dataset_finalize(hap)
  bed_manifest<-.gsim_plink_dataset_finalize(bed);complete<-TRUE
  structure(list(hap=hap_manifest,bed=bed_manifest,offspring=offspring,parent_ids=used,
    simulation=list(seed=seed,threads=1L,dense_genotype_matrix_allocated=FALSE,
      peak_packed_chromosome_bytes=max(payload),source_parent_count=nrow(reader$samples),
      used_parent_count=length(used),offspring_count=nrow(offspring),
      meiosis="identity-keyed chromosome-wise Poisson no-interference",
      publication="each triplet atomic; two-triplet publication is not atomic")),class="gsim_cohort")
}
