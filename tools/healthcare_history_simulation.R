# Experimental development helpers, not exported package APIs.
# Continuous event times; intensities held constant until an event/grid boundary.
.hc_check <- function(ok, message) if (!isTRUE(ok)) stop(message, call.=FALSE)
.hc_rates <- function(age, first, burden, global, module, sex, book, scenario) {
  seen <- is.finite(first); recent <- numeric(length(first))
  recent[seen] <- exp(-(age-first[seen])/8)
  lp <- log(book$base_rate)+book$age_effect*(age-50)+book$sex_effect*sex+
    .4*global+.6*module[book$module]-.26+.18*burden[book$module]
  target <- nrow(book)
  lp[target] <- lp[target]+.3*recent[1]-.2*recent[3]+.15*sum(recent[book$module==book$module[target]])
  if(scenario!='additive') {
    lp[target] <- lp[target]+.8*recent[1]*recent[2]-.4*recent[3]*recent[4]
    lp[6] <- lp[6]+.6*recent[1]*recent[3]
    lp[7] <- lp[7]+.5*recent[2]*recent[4]
  }
  if(scenario=='order' && seen[1] && seen[2]) {
    if(first[1]<first[2] && first[2]-first[1]<5) lp[target] <- lp[target]+recent[2]
    if(first[2]<first[1] && first[1]-first[2]<5) lp[target] <- lp[target]-.35*recent[1]
  }
  rate <- exp(lp); rate[seen] <- 0
  death <- .0012*exp(.075*(age-50)+.04*sum(seen)+.3*global+.8*seen[target])
  ans <- c(rate,death)
  .hc_check(all(is.finite(ans)) && all(ans>=0) && is.finite(sum(ans)) && sum(ans)>0,'Invalid process intensities')
  ans
}
.hc_draw <- function(rate) {
  total <- sum(rate)
  .hc_check(length(rate)>0 && all(is.finite(rate)) && all(rate>=0) && is.finite(total) && total>0,'Invalid event rates')
  c(wait=stats::rexp(1,rate=total),event=which(cumsum(rate)>stats::runif(1)*total)[1])
}
hc_simulate_histories <- function(n=10000L, codes=128L, modules=8L, seed=20261103L,
                                 scenario=c('interaction','additive','order'),
                                 dropout_rate=.02, delay_mean=0, grid_step=1,
                                 progress=500L) {
  scenario <- match.arg(scenario)
  .hc_check(length(n)==1 && is.finite(n) && n==as.integer(n) && n>0,'n must be a positive integer')
  .hc_check(codes==as.integer(codes) && codes>=8 && modules==as.integer(modules) && modules>=2 && modules<=codes,'Invalid code/module count')
  .hc_check(is.finite(seed) && seed==floor(seed) && seed>=0 && seed<2e9 && is.finite(dropout_rate) && dropout_rate>=0 && is.finite(delay_mean) && delay_mean>=0 && is.finite(grid_step) && grid_step>0,'Invalid process controls')
  old_kind <- RNGkind(); had_seed <- exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) old_seed <- get('.Random.seed',envir=.GlobalEnv)
  on.exit({do.call(RNGkind,as.list(old_kind));if(had_seed) assign('.Random.seed',old_seed,envir=.GlobalEnv) else if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv)},add=TRUE)
  RNGkind('Mersenne-Twister','Inversion','Rejection')
  j <- seq_len(codes)
  book <- data.frame(code=sprintf('D%03d',j-1L),module=(j-1L)%%modules+1L,
    base_rate=.035/j^.9,age_effect=.025+.01*sin(j),sex_effect=.2*cos(j*pi/4),
    target=j==codes,stringsAsFactors=FALSE)
  book$base_rate[codes] <- .005; book$age_effect[codes] <- .045
  # The ontology is observable; generating coefficients/frailties are truth-only.
  people <- vector('list',n); observed <- vector('list',n); clinical <- vector('list',n)
  latent <- matrix(NA_real_,n,modules+1L)
  for(i in seq_len(n)) {
    patient_seed <- as.integer((as.double(seed)+as.double(i)*104729)%%2147483646)+1L
    set.seed(patient_seed)
    split_u <- runif(1); split <- if(split_u<.7) 'train' else if(split_u<.85) 'validation' else 'test'
    landmark <- runif(1,50,70); sex <- rbinom(1,1,.5); global <- rnorm(1); module <- rnorm(modules)
    latent[i,] <- c(global,module)
    observation_start <- max(30,landmark-runif(1,10,25))
    # Dropout is drawn BEFORE the path, independent of its event count/outcomes.
    censor_age <- if(dropout_rate>0) landmark+rexp(1,dropout_rate) else Inf
    observation_end <- min(85,censor_age)
    age <- 30; next_grid <- min(85,30+grid_step); first <- rep(Inf,codes); burden <- numeric(modules)
    event_age <- numeric(codes+1L); event_code <- integer(codes+1L); count <- 0L; death_age <- Inf
    while(age<85) {
      rate <- .hc_rates(age,first,burden,global,module,sex,book,scenario)
      proposal <- .hc_draw(rate)
      if(age+proposal[1]>=next_grid) {
        burden <- burden*exp(-(next_grid-age)/3); age <- next_grid;next_grid <- min(85,next_grid+grid_step)
      } else {
        next_age <- age+proposal[1]
        .hc_check(is.finite(next_age) && next_age>age,'Event-time arithmetic failure')
        burden <- burden*exp(-(next_age-age)/3);age <- next_age
        event <- as.integer(proposal[2]);count <- count+1L;event_age[count] <- age;event_code[count] <- event
        if(event==codes+1L) {death_age <- age;break}
        .hc_check(!is.finite(first[event]),'Repeated first-diagnosis event')
        first[event] <- age;burden[book$module[event]] <- burden[book$module[event]]+1
      }
    }
    keep <- seq_len(count); event_age <- event_age[keep];event_code <- event_code[keep]
    set.seed(as.integer((as.double(patient_seed)+999983)%%2147483646)+1L)
    delay <- if(delay_mean>0) rexp(count,rate=1/delay_mean) else numeric(count)
    delay[event_code==codes+1L] <- 0
    available <- event_age+delay
    target_age <- if(is.finite(first[codes])) available[which(event_code==codes)[1]] else Inf
    eligible <- death_age>landmark && target_age>landmark
    id <- sprintf('P%06d',i)
    people[[i]] <- data.frame(person_id=id,split=split,landmark_age=landmark,sex_code=sex,
      observation_start_age=observation_start,eligible=eligible,
      eligibility_reason=if(death_age<=landmark) 'death_before_landmark' else if(target_age<=landmark) 'known_target_before_landmark' else 'eligible',
      stringsAsFactors=FALSE)
    clinical[[i]] <- data.frame(person_id=rep(id,count),code_index=event_code,
      event_age=event_age,available_age=available,stringsAsFactors=FALSE)
    o <- event_age>=observation_start & available<=min(observation_end,death_age)
    observed[[i]] <- clinical[[i]][o,,drop=FALSE]
    # Future information stays in truth, not the baseline patient table.
    attr(people[[i]],'followup') <- c(target=unname(target_age),death=unname(death_age),censor=unname(observation_end))
    if(progress>0 && i%%progress==0) {cat('Simulated ',i,'/',n,' patients\n',sep='');flush.console()}
  }
  followup <- lapply(seq_len(n),function(i) {
    p <- people[[i]];if(!p$eligible) return(NULL)
    f <- attr(p,'followup');stop_age <- min(p$landmark_age+5,f['censor'])
    terminal <- min(f['target'],f['death'])
    .hc_check(!(is.finite(terminal) && f['target']==f['death']),'Unresolved target/death tie')
    cause <- if(terminal<=stop_age) if(f['target']<f['death']) 1L else 2L else 0L
    time <- (if(cause>0L) terminal else stop_age)-p$landmark_age
    out <- data.frame(person_id=p$person_id,split=p$split,time=time,cause=cause,
      censor_kind=if(cause>0L) 'event' else if(f['censor']<p$landmark_age+5) 'dropout' else 'administrative',stringsAsFactors=FALSE)
    for(h in c(1,3,5)) out[[paste0('target_by_',h)]] <- if(cause>0L && time<=h) as.integer(cause==1L) else if(time>=h) 0L else NA_integer_
    out
  })
  followup <- do.call(rbind,followup)
  people <- do.call(rbind,people);rownames(people) <- NULL
  attr(people,'followup') <- NULL
  # Future follow-up must never survive as a patient-table attribute.
  observed <- do.call(rbind,observed);clinical <- do.call(rbind,clinical)
  names <- c(book$code,'DEATH');observed$event_code <- names[observed$code_index];clinical$event_code <- names[clinical$code_index]
  ix <- match(observed$person_id,people$person_id)
  history <- observed[people$eligible[ix] & observed$event_age<=people$landmark_age[ix] & observed$available_age<=people$landmark_age[ix] & observed$code_index<=codes,,drop=FALSE]
  history <- history[order(history$person_id,history$event_age,history$code_index),];rownames(history) <- NULL
  latent <- data.frame(person_id=people$person_id,latent);names(latent) <- c('person_id','global',paste0('module_',seq_len(modules)))
  list(patients=people,history=history,observed=observed,followup=followup,
    vocabulary=book[,c('code','module','target')],
    truth=list(coefficients=book,latent=latent,clinical=clinical),
    config=list(n=n,codes=codes,modules=modules,seed=seed,scenario=scenario,
      dropout_rate=dropout_rate,delay_mean=delay_mean,grid_step=grid_step,
      time_unit='years',start_age=30,end_age=85,landmark_range=c(50,70),
      endpoint='first available target diagnosis versus death',horizons=c(1,3,5)))
}
hc_prepare_views <- function(sim,min_pair_count=20L,max_pairs=256L) {
  .hc_check(requireNamespace('Matrix',quietly=TRUE),'Matrix is required for sparse feature preparation')
  .hc_check(length(min_pair_count)==1 && is.finite(min_pair_count) && min_pair_count==as.integer(min_pair_count) && min_pair_count>0 && length(max_pairs)==1 && is.finite(max_pairs) && max_pairs==as.integer(max_pairs) && max_pairs>0,'Invalid pair catalog controls')
  p <- sim$patients[sim$patients$eligible,,drop=FALSE];h <- sim$history;k <- nrow(sim$vocabulary);n <- nrow(p)
  .hc_check(n>0 && any(p$split=='train'),'No eligible training cohort')
  i <- match(h$person_id,p$person_id);code <- h$code_index
  .hc_check(all(!is.na(i)) && all(code %in% seq_len(k)) && all(is.finite(h$event_age)) && all(is.finite(h$available_age)) && all(h$event_age<=p$landmark_age[i]) && all(h$available_age<=p$landmark_age[i]) && all(h$event_age>=p$observation_start_age[i]) && !anyDuplicated(paste(h$person_id,code,sep=':')),'History is not a unique as-of first-occurrence prefix')
  lag <- p$landmark_age[i]-h$event_age
  .hc_check(all(is.finite(lag)) && all(lag>=0 & lag<=40),'History outside frozen basis support')
  tr <- p$split=='train';ntr <- sum(tr)
  norm <- function(x,width,group_names) {
    energy <- Matrix::colSums(x[tr,,drop=FALSE]^2)
    scale <- sqrt(rowSums(matrix(energy,nrow=length(group_names),byrow=TRUE))/ntr)
    scale[scale==0] <- 1
    labels <- dimnames(x)
    x <- x %*% Matrix::Diagonal(x=1/rep(scale,each=width))
    dimnames(x) <- labels
    list(x=x,normalizers=data.frame(group=group_names,scale=scale,stringsAsFactors=FALSE))
  }
  window <- findInterval(lag,c(0,1,3,5,Inf))
  xw <- Matrix::sparseMatrix(i=i,j=(code-1L)*4L+window,x=1,dims=c(n,k*4L))
  namesw <- as.vector(t(outer(sim$vocabulary$code,c('0_1','1_3','3_5','5_plus'),paste,sep='_')))
  colnames(xw) <- namesw
  w <- norm(xw,4,sim$vocabulary$code)
  basis <- splines::bs(lag,knots=c(1,3,5,10,20),Boundary.knots=c(0,40),degree=3,intercept=TRUE)
  width <- ncol(basis);ii <- rep(i,width);jj <- rep((code-1L)*width,width)+rep(seq_len(width),each=length(i));value <- as.numeric(basis);nz <- value!=0
  xs <- Matrix::sparseMatrix(i=ii[nz],j=jj[nz],x=value[nz],dims=c(n,k*width))
  colnames(xs) <- as.vector(t(outer(sim$vocabulary$code,paste0('b',seq_len(width)-1L),paste,sep='_')))
  s <- norm(xs,width,sim$vocabulary$code)
  presence <- Matrix::sparseMatrix(i=i,j=code,x=1,dims=c(n,k))
  cooccur <- as.matrix(Matrix::crossprod(presence[tr,,drop=FALSE]))
  candidates <- which(upper.tri(cooccur) & cooccur>=min_pair_count,arr.ind=TRUE)
  pairs <- data.frame(a=candidates[,1],b=candidates[,2],training_cooccurrence=cooccur[candidates],stringsAsFactors=FALSE)
  pairs <- pairs[order(-pairs$training_cooccurrence,pairs$a,pairs$b),,drop=FALSE]
  pairs <- head(pairs,max_pairs);rownames(pairs) <- NULL
  .hc_check(nrow(pairs)>0,'No training-supported interaction pairs')
  first <- matrix(NA_real_,n,k);first[cbind(i,code)] <- h$event_age
  pi <- pj <- integer();pv <- numeric()
  for(j in seq_len(nrow(pairs))) {
    a <- first[,pairs$a[j]];b <- first[,pairs$b[j]];both <- which(!is.na(a)&!is.na(b))
    ab <- both[a[both]<b[both]];ba <- both[b[both]<a[both]]
    pi <- c(pi,both,ab,ba);pj <- c(pj,rep(3*j-2L,length(both)),rep(3*j-1L,length(ab)),rep(3*j,length(ba)))
  }
  xp <- Matrix::sparseMatrix(i=pi,j=pj,x=1,dims=c(n,3*nrow(pairs)))
  pairnames <- paste(sim$vocabulary$code[pairs$a],sim$vocabulary$code[pairs$b],sep=':')
  colnames(xp) <- as.vector(t(outer(pairnames,c('both','a_before_b','b_before_a'),paste,sep='_')))
  pair <- norm(xp,3,pairnames)
  nuisance <- Matrix::Matrix(cbind(intercept=1,age_centered_decades=(p$landmark_age-60)/10,sex_code=p$sex_code,observed_history_decades=(p$landmark_age-p$observation_start_age)/10),sparse=TRUE)
  views <- list(baseline=nuisance,window=cbind(nuisance,w$x),spline=cbind(nuisance,s$x),window_interactions=cbind(nuisance,w$x,pair$x))
  for(name in names(views)) rownames(views[[name]]) <- p$person_id
  pairs$code_a <- sim$vocabulary$code[pairs$a];pairs$code_b <- sim$vocabulary$code[pairs$b]
  list(patients=p,views=views,pairs=pairs,
    normalizers=list(window=w$normalizers,spline=s$normalizers,pair=pair$normalizers),
    design=list(windows=c(0,1,3,5,Inf),cubic_knots=c(1,3,5,10,20),boundary=c(0,40),
      training_only=TRUE,min_pair_count=min_pair_count,max_pairs=max_pairs))
}
hc_complete_followup <- function(sim) {
  # Realized future outcomes with dropout removed, NOT conditional risk probabilities.
  p <- sim$patients;clinical <- sim$truth$clinical;k <- nrow(sim$vocabulary)
  row <- match(clinical$person_id,p$person_id)
  by_person <- split(clinical,factor(row,levels=seq_len(nrow(p))))
  out <- lapply(which(p$eligible),function(i) {
    d <- by_person[[i]];target <- min(c(d$available_age[d$code_index==k],Inf));death <- min(c(d$available_age[d$code_index==k+1L],Inf))
    terminal <- min(target,death);landmark <- p$landmark_age[i]
    .hc_check(terminal>landmark,'Complete follow-up disagrees with baseline eligibility')
    cause <- if(terminal<=landmark+5) if(target<death) 1L else 2L else 0L
    time <- min(terminal-landmark,5)
    z <- data.frame(person_id=p$person_id[i],split=p$split[i],time=time,cause=cause,stringsAsFactors=FALSE)
    for(h in c(1,3,5)) {z[[paste0('target_by_',h)]] <- as.integer(cause==1L && time<=h);z[[paste0('death_by_',h)]] <- as.integer(cause==2L && time<=h)}
    z
  })
  do.call(rbind,out)
}
hc_export_benchmark <- function(sim,prepared,directory) {
  .hc_check(!file.exists(directory) && dir.exists(dirname(directory)),'Output must be a new directory under an existing parent')
  dir.create(directory);dir.create(file.path(directory,'model_inputs'));dir.create(file.path(directory,'evaluation'));dir.create(file.path(directory,'truth'))
  write <- function(x,path) write.csv(x,file.path(directory,path),row.names=FALSE,na='NA')
  write(prepared$patients,'model_inputs/patients.csv');write(sim$vocabulary,'model_inputs/vocabulary.csv')
  write(sim$history,'model_inputs/landmark_history.csv');write(prepared$pairs,'model_inputs/candidate_pairs.csv')
  write(sim$followup[sim$followup$split=='train',],'model_inputs/train_followup.csv')
  write(sim$followup[sim$followup$split=='validation',],'model_inputs/validation_followup.csv')
  write(sim$followup[sim$followup$split=='test',],'evaluation/test_followup.csv')
  e <- sim$observed;ix <- match(e$person_id,sim$patients$person_id)
  selected <- sim$patients$eligible[ix] & e$available_age<=sim$patients$landmark_age[ix]+5
  e <- e[selected,];ix <- match(e$person_id,sim$patients$person_id)
  write(e[sim$patients$split[ix]=='train',],'model_inputs/train_trajectories.csv')
  write(e[sim$patients$split[ix]=='validation',],'model_inputs/validation_trajectories.csv')
  write(e[sim$patients$split[ix]=='test',],'evaluation/test_trajectories.csv')
  for(name in names(prepared$views)) {
    Matrix::writeMM(prepared$views[[name]],file.path(directory,'model_inputs',paste0('X_',name,'.mtx')))
    write(data.frame(column=seq_len(ncol(prepared$views[[name]])),feature_id=colnames(prepared$views[[name]])),paste0('model_inputs/features_',name,'.csv'))
  }
  saveRDS(list(design=prepared$design,normalizers=prepared$normalizers,rows=prepared$patients[,c('person_id','split')]),file.path(directory,'model_inputs/preprocessing.rds'))
  saveRDS(sim$config,file.path(directory,'simulation_config.rds'))
  saveRDS(sim$truth,file.path(directory,'truth/generating_state.rds'))
  write(hc_complete_followup(sim),'truth/complete_followup.csv')
  write(sim$patients,'evaluation/cohort_eligibility.csv')
  writeLines(c('SIMULATION ONLY; not fitted to Danish healthcare data.',
    'X views contain the same as-of observed histories. Latent states and true interactions are not predictors.',
    'Test futures are evaluation-only. Pair candidates/normalizers use eligible TRAINING histories only.',
    'Follow-up is censored: NA horizon labels are unknown, not negative. Do not score complete cases naively.',
    'Truth/complete_followup.csv contains realized outcomes without dropout, for simulation test scoring only; never fill training labels from it.',
    'No conditional-history oracle risks are claimed; latent-state risks differ from accessible prediction targets.',
    'Delphi vocabulary/training/loss and BLR inference qualification remain separate work.'),file.path(directory,'README.txt'))
  paths <- list.files(directory,recursive=TRUE,full.names=TRUE)
  write.csv(data.frame(file=substring(paths,nchar(directory)+2),md5=unname(tools::md5sum(paths))),file.path(directory,'file_hashes.csv'),row.names=FALSE)
  invisible(directory)
}
