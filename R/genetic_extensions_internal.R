# Explicit diploid effects; no extra draws in the default additive model.
.gsim_mutation_rate <- function(x) {
  if (!is.numeric(x) || length(x)!=1L || !is.finite(x) || x<0 || x>1)
    .gsim_stop("mutation_rate must be one finite probability in [0,1].")
  as.double(x)
}
.gsim_environment <- function(x, ids, name="environment") {
  if (!is.numeric(x)) .gsim_stop(name," must be numeric.")
  if(length(x)==1L) x<-rep(x,length(ids)) else {
    if(is.null(names(x)) || anyNA(names(x)) || anyDuplicated(names(x)) || !setequal(names(x),ids))
      .gsim_stop(name," must be scalar or name every identity exactly once.")
    x<-x[match(ids,names(x))]
  }
  if(length(x)!=length(ids) || any(!is.finite(x))) .gsim_stop(name," must be finite and aligned.")
  setNames(as.double(x),ids)
}
.gsim_extension_effects <- function(x,markers,traits,name) {
  if(is.null(x))return(NULL)
  if(!is.numeric(x)) .gsim_stop(name," must be numeric effects.")
  if(is.null(rownames(x)) && is.null(names(x))) .gsim_stop(name," effects require marker names.")
  if(is.matrix(x) && !is.null(colnames(x))) {
    if(anyDuplicated(colnames(x)) || !setequal(colnames(x),traits)) .gsim_stop(name," trait names must match.")
    x<-x[,traits,drop=FALSE]
  }
  x<-.gsim_align_fixed_beta(x,markers,length(traits));colnames(x)<-traits;x
}
# Reads no more than 64 decoded columns. Native BED products retain one record.
.gsim_extension_values <- function(effects,heterozygote,markers,ids,load,plan=NULL,chunk_size=64L) {
  result<-matrix(0,length(ids),ncol(effects),dimnames=list(ids,colnames(effects)))
  active<-which(rowSums(abs(effects))>0)
  if(!length(active))return(result)
  for(start in seq.int(1L,length(active),by=chunk_size)) {
    at<-active[start:min(start+chunk_size-1L,length(active))];selected<-markers[at]
    W<-load(selected)
    if(any(!is.na(W) & (!is.finite(W) | !(W %in% 0:2))))
      .gsim_stop("Dominance and interaction models require diploid hard calls 0, 1, 2 or NA.")
    X<-if(heterozygote)1.0*(W==1) else W
    replacement<-colMeans(X,na.rm=TRUE)
    if(any(!is.finite(replacement))) .gsim_stop("Non-additive effect marker has no observed calls.")
    if(is.null(plan)) {
      for(j in seq_len(ncol(X)))X[is.na(X[,j]),j]<-replacement[j]
      result<-result+X%*%effects[at,,drop=FALSE]
    } else {
      pos<-match(selected,plan$marker_ids)
      native<-if(heterozygote)C_gsim_bed_dominance else C_gsim_bed_accumulate
      result<-result+.Call(native,plan$files,plan$file_index[pos],plan$variant_index[pos],
        as.double(replacement),rep(0,length(at)),rep(1,length(at)),effects[at,,drop=FALSE])
    }
  }
  result
}
