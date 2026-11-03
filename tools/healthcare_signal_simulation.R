# Development-only positive control; source healthcare_history_simulation.R first.
hc_simulate_strong_additive <- function(n=10000L,codes=32L,seed=20261104L,dropout_rate=.02) {
  .hc_check(n>=100 && n==as.integer(n) && codes>=8 && codes==as.integer(codes), 'Invalid control dimensions')
  .hc_check(is.finite(dropout_rate) && dropout_rate>=0,'Invalid dropout rate')
  had <- exists('.Random.seed',.GlobalEnv,inherits=FALSE);if(had) old <- .Random.seed
  on.exit(if(had) assign('.Random.seed',old,.GlobalEnv) else if(exists('.Random.seed',.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv))
  set.seed(seed)
  id <- sprintf('P%06d',seq_len(n));split <- sample(c('train','validation','test'),n,TRUE,c(.7,.15,.15))
  age <- runif(n,50,70);sex <- rbinom(n,1,.5)
  p <- data.frame(person_id=id,split=split,landmark_age=age,sex_code=sex,observation_start_age=age-15,eligible=TRUE,eligibility_reason='eligible')
  book <- data.frame(code=sprintf('D%03d',0:(codes-1L)),module=(0:(codes-1L))%%4+1L,target=seq_len(codes)==codes)
  present <- matrix(rbinom(n*(codes-1L),1,.25),n,codes-1L)
  # Exposures are observed, common, and independent of subsequent outcome noise.
  beta <- c(1.6,1.2,-1.4,.9,rep(0,codes-5L))
  eta <- -4+.2*(age-60)/10+.15*sex+as.numeric(present %*% beta)
  eta_death <- -5+.3*(age-60)/10
  denom <- 1+exp(eta)+exp(eta_death);pc <- exp(eta)/denom;pd <- exp(eta_death)/denom
  event_time <- rep(5,n);cause <- integer(n)
  for(year in 1:5) {
    u <- runif(n);hit <- cause==0 & u<pc+pd
    cause[hit] <- ifelse(u[hit]<pc[hit],1L,2L)
    event_time[hit] <- year-1+runif(sum(hit),.001,.999)
  }
  censor <- if(dropout_rate>0) pmin(5,rexp(n,dropout_rate)) else rep(5,n)
  observed_cause <- ifelse(event_time<=censor,cause,0L);time <- pmin(event_time,censor)
  f <- data.frame(person_id=id,split=split,time=time,cause=observed_cause,censor_kind=ifelse(observed_cause>0,'event',ifelse(censor<5,'dropout','administrative')))
  for(H in c(1,3,5)) f[[paste0('target_by_',H)]] <- ifelse(f$cause>0 & time<=H,as.integer(f$cause==1),ifelse(time>=H,0L,NA_integer_))
  ij <- which(present==1,arr.ind=TRUE);t <- age[ij[,1]]-runif(nrow(ij),.05,14.95)
  h <- data.frame(person_id=id[ij[,1]],code_index=ij[,2],event_age=t,available_age=t,event_code=book$code[ij[,2]])
  h <- h[order(h$person_id,h$event_age),];rownames(h) <- NULL
  terminal <- which(cause>0);tc <- ifelse(cause[terminal]==1,codes,codes+1L)
  future <- data.frame(person_id=id[terminal],code_index=tc,event_age=age[terminal]+event_time[terminal],available_age=age[terminal]+event_time[terminal],event_code=ifelse(tc==codes,book$code[codes],'DEATH'))
  clinical <- rbind(h,future);obs <- clinical[clinical$available_age<=age[match(clinical$person_id,id)]+censor[match(clinical$person_id,id)],]
  risk <- data.frame(person_id=id,split=split,annual_cancer=pc,annual_death=pd)
  for(H in c(1,3,5)) risk[[paste0('risk_',H)]] <- pc/(pc+pd)*(-expm1(H*log1p(-pc-pd)))
  list(patients=p,history=h,observed=obs,followup=f,vocabulary=book,
       truth=list(coefficients=data.frame(code=book$code[-codes],log_odds=beta),latent=NULL,clinical=clinical,exact_risk=risk),
       config=list(n=n,codes=codes,modules=4L,seed=seed,scenario='strong_additive',dropout_rate=dropout_rate,delay_mean=0,time_unit='years',horizons=c(1,3,5),
       endpoint='annual multinomial cancer versus death',cancer_intercept=-4,death_intercept=-5,age_effect=.2,sex_effect=.15,
       note='Positive control, not clinical realism: fixed observed histories, no intermediate future diagnoses; uniform event time inside annual category.'))
}

# Exact full-state conditional forward simulation for the original continuous process.
# Uses truth-only history and frailties; never supplies these quantities to fitters.
hc_latent_forward_risk <- function(sim,ids,paths=128L,seed=4107L) {
  .hc_check(paths>=2 && paths==as.integer(paths),'Invalid path count')
  .hc_check(identical(as.numeric(sim$config$grid_step),1) && identical(as.numeric(sim$config$delay_mean),0),'Forward diagnostic requires unit annual grid and zero report delay')
  .hc_check(!anyDuplicated(ids) && all(ids %in% sim$patients$person_id),'Invalid diagnostic patient IDs')
  had <- exists('.Random.seed',.GlobalEnv,inherits=FALSE);if(had) old <- .Random.seed
  on.exit(if(had) assign('.Random.seed',old,.GlobalEnv) else rm('.Random.seed',envir=.GlobalEnv))
  set.seed(seed);book <- sim$truth$coefficients;k <- nrow(book);out <- vector('list',length(ids))
  by <- split(sim$truth$clinical,sim$truth$clinical$person_id)
  for(i in seq_along(ids)) {
    p <- sim$patients[match(ids[i],sim$patients$person_id),];L <- p$landmark_age
    d <- by[[ids[i]]];if(is.null(d)) d <- sim$truth$clinical[FALSE,,drop=FALSE]
    d <- d[d$event_age<=L & d$code_index<=k,,drop=FALSE]
    first0 <- rep(Inf,k);first0[d$code_index] <- d$event_age
    boundary <- max(floor(L),d$event_age)
    burden0 <- numeric(sim$config$modules)
    if(nrow(d)) for(j in seq_len(nrow(d))) burden0[book$module[d$code_index[j]]] <- burden0[book$module[d$code_index[j]]]+exp(-(boundary-d$event_age[j])/3)
    z <- sim$truth$latent[match(ids[i],sim$truth$latent$person_id),]
    global <- z$global;module <- as.numeric(z[paste0('module_',seq_len(sim$config$modules))])
    rate0 <- .hc_rates(boundary,first0,burden0,global,module,p$sex_code,book,sim$config$scenario)
    y <- matrix(0,paths,3)
    for(b in seq_len(paths)) {
      age <- L;first <- first0;burden <- burden0*exp(-(L-boundary)/3);rate <- rate0;grid <- floor(L)+1
      repeat {
        prop <- .hc_draw(rate);next_age <- age+prop[1]
        if(min(next_age,grid)>=L+5) break
        if(next_age>=grid) {burden <- burden*exp(-(grid-age)/3);age <- grid;grid <- grid+1}
        else {
          burden <- burden*exp(-(next_age-age)/3);age <- next_age;event <- as.integer(prop[2])
          if(event==k+1L) break
          if(event==k) {y[b,] <- as.numeric(age-L<=c(1,3,5));break}
          first[event] <- age;burden[book$module[event]] <- burden[book$module[event]]+1
        }
        rate <- .hc_rates(age,first,burden,global,module,p$sex_code,book,sim$config$scenario)
      }
    }
    out[[i]] <- data.frame(person_id=ids[i],horizon=c(1,3,5),risk=colMeans(y),mcse=sqrt(colMeans(y)*(1-colMeans(y))/paths),paths=paths)
    if(i%%200==0) {cat('Latent-risk patients ',i,'/',length(ids),'\n',sep='');flush.console()}
  }
  do.call(rbind,out)
}
