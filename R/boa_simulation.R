#' Transmit phased genotypes and breed origins through a pedigree
#'
#' Replays the existing identity-keyed meiosis plan for breed labels, preserving
#' the allele/origin phase at every marker. Founders must be purebred and every
#' offspring must have two known parents. Origin buffers use one byte per copy;
#' genotype buffers remain packed. No fitting or pedigree covariance draws occur.
#' @param pedigree A [gsim_pedigree_from_table()] pedigree.
#' @param founder_haplotypes List with binary founder-by-marker `h1`, `h2`
#'   matrices with identical animal and marker names.
#' @param founder_breed Named breed labels, exactly one per founder.
#' @param chromosome Chromosome label per marker; chromosomes must be contiguous.
#' @param genetic_position Nondecreasing Morgan positions within chromosome.
#' @param seed Identity-keyed meiosis seed; R's RNG is unused.
#' @param mutation_rate Per-marker, per-transmitted-copy allele-flip probability.
#'   Breed origins retain the inherited segment label after mutation.
#' @param maximum_ancestry_bytes Admission limit for the two byte origin buffers.
#' @return BOA panel with packed chromosome blocks, zero-based phased origins,
#'   pedigree breed fractions, frozen founder frequencies and crossover audits.
#' @export
gsim_boa_genotypes <- function(pedigree, founder_haplotypes, founder_breed,
    chromosome, genetic_position, seed=1, maximum_ancestry_bytes=256*1024^2, mutation_rate=0) {
  mutation_rate<-.gsim_mutation_rate(mutation_rate)
  if (!inherits(pedigree,"gsim_pedigree")) .gsim_stop("pedigree must be a gsim_pedigree.")
  ids <- pedigree$canonical_order
  tab <- pedigree$pedigree[match(ids,pedigree$pedigree$animal),,drop=FALSE]
  founders <- ids[is.na(tab$sire)&is.na(tab$dam)]
  if (is.null(names(founder_breed)) || anyDuplicated(names(founder_breed)) ||
      !setequal(names(founder_breed),founders)) .gsim_stop("founder_breed names must match founders.")
  founder_breed <- as.character(founder_breed[founders]); names(founder_breed)<-founders
  if(anyNA(founder_breed)||any(!nzchar(founder_breed))) .gsim_stop("Missing founder breed.")
  breeds<-sort(unique(founder_breed),method="radix");B<-length(breeds)
  if(B<2L||B>255L) .gsim_stop("BOA requires 2 to 255 founder breeds.")
  h1<-.gsim_hapnest_raw_matrix(founder_haplotypes$h1,"founder H1")
  h2<-.gsim_hapnest_raw_matrix(founder_haplotypes$h2,"founder H2")
  if(!identical(dimnames(h1),dimnames(h2))||!identical(dim(h1),dim(h2))||
      is.null(rownames(h1))||!setequal(rownames(h1),founders)||anyDuplicated(rownames(h1))||
      is.null(colnames(h1))||anyDuplicated(colnames(h1))) .gsim_stop("Founder phased dimensions and identities must agree.")
  M<-ncol(h1);N<-length(ids);chr<-as.character(chromosome);pos<-as.double(genetic_position)
  if(length(chr)!=M||anyNA(chr)||any(!nzchar(chr))||length(pos)!=M||any(!is.finite(pos))||
      anyDuplicated(rle(chr)$values)) .gsim_stop("Invalid contiguous chromosome map.")
  if(!is.numeric(maximum_ancestry_bytes)||length(maximum_ancestry_bytes)!=1L||
      !is.finite(maximum_ancestry_bytes)||maximum_ancestry_bytes<2*as.double(N)*M)
    .gsim_stop("BOA ancestry buffers exceed maximum_ancestry_bytes.")
  seed<-.gsim_block_seed(seed)
  fractions<-matrix(0,N,B,dimnames=list(ids,breeds))
  for(i in seq_len(N)) {
    if(ids[i]%in%founders) fractions[i,match(founder_breed[ids[i]],breeds)]<-1 else {
      parents<-match(c(tab$sire[i],tab$dam[i]),ids)
      if(anyNA(parents)||any(parents>=i)) .gsim_stop("BOA requires two earlier known parents.")
      fractions[i,]<-colMeans(fractions[parents,,drop=FALSE])
    }
  }
  blocks<-lapply(unique(chr),function(label) {
    at<-which(chr==label)
    block<-.gsim_pedigree_genotypes_packed_chromosome(pedigree,
      list(h1=h1[,at,drop=FALSE],h2=h2[,at,drop=FALSE]),chr[at],pos[at],seed,
      return_haplotypes=TRUE,return_genotypes=FALSE,return_crossovers=TRUE,mutation_rate=mutation_rate)
    a1<-a2<-matrix(as.raw(0),N,length(at),dimnames=list(ids,colnames(h1)[at]))
    for(i in seq_len(N)) {
      if(ids[i]%in%founders) {
        a1[i,]<-a2[i,]<-as.raw(match(founder_breed[ids[i]],breeds)-1L)
      } else for(side in 1:2) {
        parent<-match(if(side==1L)tab$sire[i] else tab$dam[i],ids)
        plan<-.Call(C_gsim_meiosis_plan,pos[at],seed,ids[i],label,as.integer(side))
        # Native boundaries are zero-based marker offsets, switching before that marker.
        phase<-(plan$starting_haplotype-1L+findInterval(seq_along(at)-1L,plan$boundaries))%%2L
        value<-a1[parent,];value[phase==1L]<-a2[parent,phase==1L]
        if(side==1L)a1[i,]<-value else a2[i,]<-value
      }
    }
    freq<-vapply(breeds,function(b) {
      rows<-match(founders[founder_breed==b],rownames(h1))
      (colMeans(matrix(as.integer(h1[rows,at,drop=FALSE]),length(rows)))+
       colMeans(matrix(as.integer(h2[rows,at,drop=FALSE]),length(rows))))/2
    },numeric(length(at)))
    freq<-matrix(freq,length(at),B,dimnames=list(colnames(h1)[at],breeds))
    block$ancestry_h1<-a1;block$ancestry_h2<-a2;block$frequencies<-freq;block
  });names(blocks)<-unique(chr)
  structure(list(blocks=blocks,breeds=breeds,pedigree=pedigree,
    pedigree_fractions=fractions,founder_breed=founder_breed,
    settings=list(seed=seed,mutation_rate=mutation_rate,ancestry_bytes=2*as.double(N)*M,
      model="supplied founder haplotypes with no-interference Mendelian meiosis",
      centering="purebred founder frequencies, frozen across generations",
      qualification="structured meiosis truth; not an exact hybrid-H Gaussian draw")),class="gsim_boa_panel")
}

.gsim_boa_effects <- function(covariance,markers,axes,seed,prefix) {
  factor<-.gsim_block_factor(covariance,axes,"BOA effect covariance")
  covariance<-factor%*%t(factor)
  G<-matrix(0,length(markers),length(axes),dimnames=list(paste0(prefix,markers),axes))
  draws<-gsim_records(G,residual_covariance=covariance,seed=seed)
  out<-matrix(0,length(markers),length(axes),dimnames=list(markers,axes))
  out[cbind(match(draws$records$animal,rownames(G)),match(draws$records$trait,axes))]<-draws$truth$residual
  out
}

#' Freeze breed-specific SNP and reference-breed LBP effects
#' @param panel A [gsim_boa_genotypes()] panel.
#' @param snp_covariance Physical per-marker covariance over breed-major channels.
#'   A matrix or named list for disjoint marker components. Block diagonal breed
#'   covariance gives BOA_UNCOR; full covariance gives BOA_COR. No rescaling occurs.
#' @param channels Names of genetic traits or random-regression coefficients.
#' @param marker_components Named component label per marker, required for a list.
#' @param lbp_covariance Optional covariance over non-reference breed-major channels.
#' @param reference_breed Reference breed for the LBP contrast, default last breed.
#' @param seed Identity-keyed marker-effect seed, distinct from meiosis/noise streams.
#' @return Frozen effects, frequencies, component assignments and covariance inputs.
#'   Use with [gsim_boa_values()] then [gsim_records()] for repeated/missing records.
#'   Up to 64 breed/channel axes are supported by the current normal-draw kernel.
#' @export
gsim_boa_trait_state <- function(panel,snp_covariance,channels="trait1",
    marker_components=NULL,lbp_covariance=NULL,reference_breed=utils::tail(panel$breeds,1),seed=1) {
  if(!inherits(panel,"gsim_boa_panel")) .gsim_stop("panel must be a BOA panel.")
  channels<-.gsim_block_ids(channels,"channels");breeds<-panel$breeds
  axes<-as.vector(t(outer(breeds,channels,paste,sep=":")))
  if(anyDuplicated(axes)) .gsim_stop("Breed/channel axis labels must be unique.")
  if(length(axes)>64L) .gsim_stop("At most 64 breed/channel axes are supported.")
  markers<-unlist(lapply(panel$blocks,function(b)b$variant_ids),use.names=FALSE)
  if(anyDuplicated(markers)) .gsim_stop("Marker IDs must be unique across chromosomes.")
  if(is.matrix(snp_covariance)||is.numeric(snp_covariance))snp_covariance<-list(snp=snp_covariance)
  if(!is.list(snp_covariance)||is.null(names(snp_covariance))||anyDuplicated(names(snp_covariance))||any(!nzchar(names(snp_covariance))))
    .gsim_stop("SNP covariance components require unique names.")
  if(is.null(marker_components)&&length(snp_covariance)==1L)marker_components<-setNames(rep(names(snp_covariance),length(markers)),markers)
  if(is.null(names(marker_components))||anyDuplicated(names(marker_components))||!setequal(names(marker_components),markers))
    .gsim_stop("marker_components must name every marker exactly once.")
  component<-as.character(marker_components[markers])
  if(anyNA(component)||!setequal(unique(component),names(snp_covariance))) .gsim_stop("Marker components must match nonempty covariance blocks.")
  seed<-.gsim_block_seed(seed);effects<-matrix(0,length(markers),length(axes),dimnames=list(markers,axes))
  for(name in names(snp_covariance)) {
    at<-which(component==name)
    effects[at,]<-.gsim_boa_effects(snp_covariance[[name]],markers[at],axes,seed,"snp:")
  }
  if(length(reference_breed)!=1L||!reference_breed%in%breeds) .gsim_stop("Unknown reference breed.")
  nonref<-setdiff(breeds,reference_breed);laxes<-as.vector(t(outer(nonref,channels,paste,sep=":")))
  lbp<-matrix(0,length(markers),length(axes),dimnames=list(markers,axes))
  if(!is.null(lbp_covariance))lbp[,laxes]<-.gsim_boa_effects(lbp_covariance,markers,laxes,seed,"lbp:")
  structure(list(effects=effects,lbp_effects=lbp,channels=channels,breeds=breeds,
    reference_breed=reference_breed,marker_components=setNames(component,markers),
    frequencies=do.call(rbind,lapply(panel$blocks,function(b)b$frequencies)),
    covariances=list(snp=snp_covariance,lbp=lbp_covariance),seed=seed,
    centering="frozen founder frequencies; no cohort variance normalization"),class="gsim_boa_trait_state")
}

#' Accumulate persistent breed-origin genetic truth from packed alleles
#' @param panel A [gsim_boa_genotypes()] panel.
#' @param state Frozen [gsim_boa_trait_state()] effects.
#' @return Separate breed values, total SNP values, total LBP and total genetic
#'   values. Covariances of these realized values are descriptive, not forced to
#'   equal nominal effect covariances. No pedigree approximation is added.
#' @export
gsim_boa_values <- function(panel,state) {
  if(!inherits(panel,"gsim_boa_panel")||!inherits(state,"gsim_boa_trait_state")||
      !identical(panel$breeds,state$breeds)) .gsim_stop("BOA panel/state breeds disagree.")
  B<-length(state$breeds);Q<-length(state$channels);N<-nrow(panel$pedigree_fractions)
  out<-matrix(0,N,(B+1L)*Q)
  for(block in panel$blocks) {
    ids<-block$variant_ids
    if(anyNA(match(ids,rownames(state$effects)))) .gsim_stop("Panel markers missing from frozen state.")
    out<-out+.Call(C_gsim_boa_values,block$h1,block$h2,block$ancestry_h1,block$ancestry_h2,
      state$frequencies[ids,,drop=FALSE],state$effects[ids,,drop=FALSE],
      state$lbp_effects[ids,,drop=FALSE],panel$pedigree_fractions)
  }
  values<-lapply(seq_len(B+1L),function(b) {
    x<-out[,(b-1L)*Q+seq_len(Q),drop=FALSE]
    dimnames(x)<-list(rownames(panel$pedigree_fractions),state$channels);x
  });names(values)<-c(state$breeds,"lbp")
  snp<-Reduce(`+`,values[seq_len(B)]);total<-snp+values[[B+1L]]
  list(breed=values[seq_len(B)],snp=snp,lbp=values[[B+1L]],genetic=total,
    realized_covariance=if(N>1L)stats::cov(total) else NULL)
}
