args<-commandArgs(TRUE)
if(length(args)!=4L)stop('usage: verify_equivalence.R isolated-library simulation.yml standalone-data fresh-R-output')
.libPaths(c(normalizePath(args[1],mustWork=TRUE),.libPaths()))
if(!requireNamespace('yaml',quietly=TRUE))stop('R equivalence check requires the optional yaml package')
library(gsim,lib.loc=normalizePath(args[1],mustWork=TRUE))
recipe<-normalizePath(args[2],mustWork=TRUE)
s<-yaml::read_yaml(recipe)
if(!is.null(s$residual_covariance))s$residual_covariance<-matrix(unlist(s$residual_covariance,use.names=FALSE),nrow=length(s$traits),byrow=TRUE)
set.seed(7721);before<-.Random.seed
gsim_simulate(s,args[4],base_directory=dirname(recipe))
stopifnot(identical(before,.Random.seed))
for(name in list.files(args[4])) {
 a<-readBin(file.path(args[4],name),'raw',n=file.info(file.path(args[4],name))$size)
 b<-readBin(file.path(args[3],name),'raw',n=file.info(file.path(args[3],name))$size)
 if(!identical(a,b))stop('R/standalone byte mismatch: ',name)
}
cat('R/standalone files and truth are byte-identical; caller RNG unchanged.\n')

# Cross-check the shared workflow against the pre-existing public R semantics.
base<-dirname(recipe)
resolve<-function(path)if(grepl("^(/|[A-Za-z]:[/\\\\])",path))path else file.path(base,path)
effects<-read.table(resolve(s$inputs$effects),colClasses=c('character',rep('numeric',length(s$traits))))
B<-as.matrix(effects[,-1,drop=FALSE]);rownames(B)<-effects[,1];colnames(B)<-s$traits
ph<-read.table(resolve(s$inputs$founders),colClasses='character')
W<-t(vapply(seq_len(nrow(ph)),function(i)as.numeric(strsplit(ph[i,2],'')[[1]])+as.numeric(strsplit(ph[i,3],'')[[1]]),numeric(nrow(B))))
state<-gsim_trait_state(B,center=setNames(colMeans(W),rownames(B)))
actual<-gsim_genetic_values(file.path(args[3],'genotypes'),state)
truth<-read.csv(file.path(args[3],'true_breeding_values.csv'))
G<-matrix(NA_real_,length(unique(truth$animal)),length(s$traits),dimnames=list(unique(truth$animal),s$traits))
G[cbind(match(truth$animal,rownames(G)),match(truth$trait,colnames(G)))]<-truth$true_breeding_value
stopifnot(max(abs(actual-G[rownames(actual),,drop=FALSE]))<1e-12)
rt<-read.csv(file.path(args[3],'record_truth.csv'))
rec<-data.frame(record=paste(rt$record_id,rt$trait,sep=':'),animal=rt$animal,trait=rt$trait,observation_unit=rt$record_id)
R<-if(is.null(s$residual_covariance))diag(length(s$traits)) else s$residual_covariance
seed<-if(is.null(s$seed))1 else s$seed
expected<-gsim_records(G,records=rec,residual_covariance=R,seed=seed+3)
stopifnot(max(abs(expected$records$value-rt$phenotype))<1e-12)
cat('Existing genetic-value and genomic-record R semantics also agree.\n')
