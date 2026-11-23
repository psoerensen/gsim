.gsim_gaussian_simulate <- function(simulation, output, base_directory) {
  allowed<-c('schema','model','design_inputs','inputs','seed','maximum_memory_bytes','maximum_dimension')
  if(is.null(names(simulation))||anyDuplicated(names(simulation))||any(!names(simulation)%in%allowed))
    stop('Invalid Gaussian simulation keys.',call.=FALSE)
  quote<-function(x)encodeString(enc2utf8(x),quote='"')
  resolve<-function(x){
    if(!is.character(x)||length(x)!=1L||is.na(x)||!nzchar(x))stop('Gaussian input paths must be nonempty strings.')
    absolute<-grepl('^(/|[A-Za-z]:[/\\\\])',x)
    normalizePath(if(absolute)x else file.path(base_directory,x),winslash='/',mustWork=TRUE)
  }
  if(!is.character(output)||length(output)!=1L||is.na(output)||!nzchar(output))stop('output must be a path')
  target<-file.path(normalizePath(dirname(output),winslash='/',mustWork=TRUE),basename(output))
  if(file.exists(target))stop('Preserve existing Gaussian simulation output.')
  for(key in c('seed','maximum_memory_bytes','maximum_dimension'))if(!is.null(simulation[[key]])) {
    v<-simulation[[key]]
    if(!is.numeric(v)||length(v)!=1L||!is.finite(v)||v!=floor(v)||v<0||v>2^53-3)
      stop(paste('Invalid Gaussian integer control:',key))
  }
  full<-!is.null(simulation$model)
  run<-function(command,args) {
    result<-system2(command[1],args=shQuote(c(command[-1],args)),stdout=TRUE,stderr=TRUE)
    status<-attr(result,'status');if(!is.null(status)&&status!=0L)stop(paste(c('Gaussian model preparation/rendering failed:',result),collapse='\n'),call.=FALSE)
  }
  if(full) {
    if(!is.null(simulation$inputs)||!is.list(simulation$design_inputs)||is.null(names(simulation$design_inputs))||anyDuplicated(names(simulation$design_inputs)))stop('Gaussian model simulation requires named design_inputs, without inputs.')
    model<-resolve(simulation$model);design<-character()
    for(key in names(simulation$design_inputs))design<-c(design,'--design-input',paste0(key,'=',resolve(simulation$design_inputs[[key]])))
    preparer<-getOption('gsim.model_preparer')
    if(is.null(preparer)){
      executable<-Sys.which('gsuite');if(!nzchar(executable))executable<-Sys.which('gsuite.cmd')
      if(!nzchar(executable))stop('Gaussian model preparation requires the matching gsuite command on PATH, or option gsim.model_preparer.')
      prepare_command<-c(executable,'prepare-simulation');collect_command<-c(executable,'collect-simulation')
    } else {
      if(!is.character(preparer)||!length(preparer)||anyNA(preparer)||any(!nzchar(preparer)))stop('gsim.model_preparer must be an executable/argument vector.')
      prepare_command<-c(preparer,'prepare');collect_command<-c(preparer,'collect')
    }
    prepared<-paste0(target,'.preparation')
    args<-c('--model',model,design,'--output',prepared)
    driver<-getOption('gsim.model_driver');if(!is.null(driver))args<-c(args,'--exe',normalizePath(driver,winslash='/',mustWork=TRUE))
    if(!is.null(simulation$maximum_dimension))args<-c(args,'--maximum-dimension',as.character(simulation$maximum_dimension))
    run(prepare_command,args);specification<-normalizePath(file.path(prepared,'specification.txt'),winslash='/',mustWork=TRUE)
  } else {
    if(!is.null(simulation$design_inputs)||!is.list(simulation$inputs)||!identical(names(simulation$inputs),'specification')||!is.null(simulation$maximum_dimension))stop('Prepared Gaussian mode requires only inputs.specification and optional seed/memory controls.')
    specification<-resolve(simulation$inputs$specification)
  }
  lines<-c('GSIM_GAUSSIAN 1',paste('specification',quote(specification)))
  for(key in c('seed','maximum_memory_bytes'))if(!is.null(simulation[[key]]))lines<-c(lines,paste(key,quote(sprintf('%.0f',simulation[[key]]))))
  .Call(C_gsim_simulation_workflow,paste0(paste(lines,collapse='\n'),'\n'),target)
  if(full)run(collect_command,c('--model',model,design,'--prepared',prepared,'--output',target))
  invisible(target)
}
